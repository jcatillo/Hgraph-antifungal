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

def process_molecules(input_path: str, output_path: str, threshold: int):
    p = Path(input_path)
    if not p.exists(): return

    # Removed the filter to preserve row count consistency
    raw_lines = p.read_text().splitlines()
    session, results = requests.Session(), []

    for i, line in enumerate(raw_lines, 1):
        clean_line = line.strip()
        
        # Handle empty lines or comments by appending an empty result row
        if not clean_line or clean_line.startswith(('#', '//')):
            results.append({"input_smile": "", "is_patented": False})
            continue

        input_smi = clean_line.split()[0]
        print(f"Row {i}: Processing {input_smi}")
        
        search_url = f"{PUBCHEM_REST}/compound/fastsimilarity_2d/smiles/cids/JSON"
        try:
            r = session.post(search_url, data={'smiles': input_smi, 'Threshold': threshold, 'MaxRecords': 50}, timeout=15)
            cids = r.json().get('IdentifierList', {}).get('CID', [])
        except: cids = []

        is_patented = False

        for cid in cids:
            if check_patent_in_toc(session, cid):
                print(f"  Found patented compound CID: {cid}")
                is_patented = True
                break
            time.sleep(0.2) 

        results.append({
            "input_smile": input_smi,
            "is_patented": is_patented
        })

    pd.DataFrame(results).to_csv(output_path, index=False)
    print(f"Results saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", default="patent_results.csv")
    parser.add_argument("-t", "--threshold", type=int, default=95)
    args = parser.parse_args()
    process_molecules(args.input, args.output, args.threshold)