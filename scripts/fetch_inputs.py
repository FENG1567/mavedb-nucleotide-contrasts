#!/usr/bin/env python3
"""Fetch explicitly selected public inputs and enforce frozen input hashes."""
import argparse
import csv
import hashlib
from pathlib import Path
import urllib.request
import time

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset',action='append',required=True,
                   help='source id in data/external/sources.csv; repeat to fetch multiple inputs')
    a=p.parse_args()
    with (ROOT/'data/external/sources.csv').open(encoding='utf-8-sig',newline='') as f:
        sources={r['id']:r for r in csv.DictReader(f)}
    for key in a.dataset:
        if key not in sources: p.error(f'unknown source id: {key}')
        r=sources[key]
        if not r['sha256'] or not r['destination'] or r['destination'].endswith('/'):
            p.error(f'{key}: manual repository/table retrieval is required; see {r["url"]}')
        dest=(ROOT/r['destination']).resolve()
        if not dest.is_relative_to(ROOT/'data/external'): p.error('destination escapes data/external')
        if dest.exists():
            if hashlib.sha256(dest.read_bytes()).hexdigest()==r['sha256']:
                print(f'{key}: frozen input already verified');continue
            raise SystemExit(f'{key}: existing input hash mismatch; preserved for inspection: {dest}')
        if 'raw.githubusercontent.com' not in r['url'] and 'ftp.ncbi.nlm.nih.gov' not in r['url']:
            p.error(f'{key}: URL is a source landing page, not a download; bundled file should already exist')
        dest.parent.mkdir(parents=True,exist_ok=True)
        tmp=dest.with_suffix(dest.suffix+'.partial')
        for attempt in range(3):
            try:
                request=urllib.request.Request(r['url'],headers={'User-Agent':'MAVE-reproducibility/2026.11.24'})
                with urllib.request.urlopen(request,timeout=90) as response, tmp.open('wb') as output:
                    while True:
                        chunk=response.read(1024*1024)
                        if not chunk: break
                        output.write(chunk)
                break
            except Exception:
                if attempt==2: raise
                time.sleep(2*(attempt+1))
        digest=hashlib.sha256(tmp.read_bytes()).hexdigest()
        if digest!=r['sha256']:
            raise SystemExit(f'{key}: HASH MISMATCH {digest}; retained as {tmp}; do not use as frozen input')
        tmp.replace(dest)
        print(f'{key}: downloaded and SHA-256 verified')

if __name__=='__main__':main()
