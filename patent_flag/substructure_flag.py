import argparse
import time
from pathlib import Path
from typing import List, Optional, Dict, Any
import pandas as pd
import requests

PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

def _collect_toc_headings(section: Dict[str, Any]) -> List[str]:
    out = []
    if toc := section.get("TOCHeading"): out.append(toc)
    for sub in section.get("Section", []):
        if isinstance(sub, dict): out.extend(_collect_toc_headings(sub))
    return out

def check_patent_in_toc(session: requests.Session, cid: int) -> bool:
    url = f"{PUBCHEM_VIEW}/index/compound/{cid}/JSON"
    try:
        r = session.get(url, timeout=10)
        if r.status_code != 200: return False
        sections = r.json().get("Record", {}).get("Section", [])
        headings = []
        for sec in sections: headings.extend(_collect_toc_headings(sec))
        return any("patent" in h.lower() for h in headings)
    except: return False

def get_substructure_cids(session: requests.Session, smiles: str) -> List[int]:
    """Handles async substructure search using ListKey polling."""
    init_url = f"{PUBCHEM_REST}/compound/substructure/smiles/{smiles}/JSON"
    try:
        req = session.get(init_url, timeout=20)
        
        # Check if we got a direct response (immediate results)
        if req.status_code == 200:
            data = req.json()
            if 'IdentifierList' in data:
                return data['IdentifierList'].get('CID', [])
        
        # Handle async response (status 202 or Waiting key present)
        if req.status_code == 202 or 'Waiting' in req.json():
            data = req.json()
            if 'Waiting' not in data:
                return []
            
            list_key = data['Waiting']['ListKey']
            print(f"  Substructure search initiated, ListKey: {list_key}")
            list_url = f"{PUBCHEM_REST}/compound/listkey/{list_key}/cids/JSON"
            
            # Poll until results are ready
            for attempt in range(15):  # Increased to 15 attempts
                time.sleep(2)
                r = session.get(list_url, timeout=10)
                result_data = r.json()
                
                if 'IdentifierList' in result_data:
                    cids = result_data['IdentifierList'].get('CID', [])
                    print(f"  Retrieved {len(cids)} CIDs after {attempt + 1} polling attempts")
                    return cids
                
                if 'Fault' in result_data:
                    print(f"  Error in polling: {result_data.get('Fault')}")
                    break
                
                # Still waiting, continue polling
                if 'Waiting' in result_data:
                    continue
            
            print(f"  Polling timed out after 15 attempts")
            
    except requests.exceptions.RequestException as e:
        print(f"  Network error: {e}")
    except Exception as e:
        print(f"  Unexpected error: {e}")
    
    return []

def process_molecules(input_path: str, output_path: str):
    p = Path(input_path)
    if not p.exists(): 
        print(f"Error: Input file {input_path} not found")
        return

    raw_lines = p.read_text().splitlines()
    session, results = requests.Session(), []

    for i, line in enumerate(raw_lines, 1):
        clean_line = line.strip()
        if not clean_line or clean_line.startswith(('#', '//')):
            results.append({"input_smile": "", "is_patented": False})
            continue

        input_smi = clean_line.split()[0]
        print(f"\nRow {i}: Substructure search for {input_smi}")
        
        cids = get_substructure_cids(session, input_smi)
        is_patented = False

        if not cids:
            print(f"  No substructure matches found")
        else:
            print(f"  Checking {min(len(cids), 50)} compounds for patent information...")
            
            # Check top 50 results for patent info to save time
            for cid in cids[:50]:
                if check_patent_in_toc(session, cid):
                    print(f"  Found patented substructure CID: {cid}")
                    is_patented = True
                    break
                time.sleep(0.2) 

        results.append({
            "input_smile": input_smi,
            "is_patented": is_patented
        })

    pd.DataFrame(results).to_csv(output_path, index=False)
    print(f"\nResults saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", default="substructure_results.csv")
    args = parser.parse_args()
    process_molecules(args.input, args.output)