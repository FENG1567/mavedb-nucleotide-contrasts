#!/usr/bin/env python3
"""Rebuild and compare BRCA2 exonic-SNV mappings without network access."""
import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from itertools import product

ROOT=Path(__file__).resolve().parents[1]
# Standard genetic code in T,C,A,G order; equivalent to the project's mapper.
AA='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE=dict(zip((''.join(c) for c in product('TCAG',repeat=3)),AA))
SNV=re.compile(r'^(ENST00000380152(?:\.\d+)?):c\.(\d+)([ACGTN])>([ACGTN])$')
PIN='c01ae0690f9267aabaaa39bb29e1f5c5cf0dd05354b9ddf236b72a8cdc810b9a'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=ROOT/'data/external/brca2_cc0')
    p.add_argument('--output-dir',type=Path,default=ROOT/'results/brca2_mapping')
    args=p.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    cdsfile=args.input_dir/'ENST00000380152.txt'
    if hashlib.sha256(cdsfile.read_bytes()).hexdigest()!=PIN:
        raise ValueError('Pinned CDS file hash mismatch')
    cds=cdsfile.read_text().strip().upper()
    checks=[]
    for urn in ['00001225-a-1','00001242-a-1']:
        records=[];excluded=0
        with (args.input_dir/f'brca2_{urn}.csv').open(encoding='utf-8-sig',newline='') as f:
            for r in csv.DictReader(f):
                try:score=float(r['score'])
                except (ValueError,TypeError):continue
                if not math.isfinite(score):continue
                m=SNV.fullmatch(r['hgvs_nt'])
                if not m:excluded+=1;continue
                pos=int(m[2]);ref,alt=m[3],m[4]
                if not 1<=pos<=len(cds):excluded+=1;continue
                if cds[pos-1]!=ref:raise ValueError(f'Reference mismatch: {r["hgvs_nt"]}')
                start=(pos-1)//3*3;rc=cds[start:start+3]
                ac=rc[:(pos-1)%3]+alt+rc[(pos-1)%3+1:]
                ra,aa=CODE.get(rc,'?'),CODE.get(ac,'?')
                cls='synonymous' if ra==aa else 'nonsense' if aa=='*' else 'stop_lost' if ra=='*' else 'missense'
                records.append(dict(gene='BRCA2',variant_id=r['accession'],hgvs_nt=r['hgvs_nt'],score=score,
                                    cds_pos=pos,codon_start=start+1,ref_codon=rc,alt_codon=ac,
                                    ref_aa=ra,alt_aa=aa,protein_pos=start//3+1,**{'class':cls}))
        with (args.input_dir/f'mapped_BRCA2_{urn}.csv').open(encoding='utf-8-sig',newline='') as f:
            frozen=list(csv.DictReader(f))
        assert len(records)==len(frozen),(urn,len(records),len(frozen))
        for fresh,old in zip(records,frozen):
            for key,value in fresh.items():
                if key=='score':assert math.isclose(value,float(old[key]),rel_tol=1e-12,abs_tol=1e-12)
                else:assert str(value)==old[key],(urn,key,value,old[key])
        with (args.output_dir/f'mapped_BRCA2_{urn}.csv').open('w',encoding='utf-8',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
        checks.append({'score_set':urn,'mapped_rows':len(records),'excluded_non_exonic_or_out_of_range':excluded,
                       'all_mapping_fields_match_frozen':True})
    (args.output_dir/'mapping_comparison.json').write_text(json.dumps(checks,indent=2),encoding='utf-8')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
