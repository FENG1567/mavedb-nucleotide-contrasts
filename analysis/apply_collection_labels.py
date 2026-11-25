#!/usr/bin/env python3
"""Apply the verified MaveDB target label boundary to new set-level outputs."""
import csv,json
from pathlib import Path
import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output-dir',type=Path,required=True)
OUT=parser.parse_args().output_dir.resolve()
URN='urn:mavedb:00001254-a-1'
def fix_csv(name,gene_col='gene'):
    p=OUT/name
    with open(p,newline='',encoding='utf-8-sig') as f:
        r=csv.DictReader(f);fields=list(r.fieldnames or []);data=list(r)
    hist='historical_manifest_gene_label'
    if hist not in fields:fields.append(hist)
    changed=0
    for row in data:
        if row.get('score_set')==URN:
            row[hist]=row.get(hist) or row.get(gene_col,'')
            row[gene_col]='LARGE1';changed+=1
        elif hist not in row or not row[hist]:
            row[hist]=row.get(gene_col,'')
    with open(p,'w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(data)
    print(f'{name}: LARGE1 resource-label rows corrected={changed}',flush=True)
for name in ('group_level_empirical_testability.csv','per_score_set_reconciliation.csv','source_provenance_per_score_set.csv'):
    fix_csv(name)
p=OUT/'collection_group_reconciliation.json'
report=json.loads(p.read_text(encoding='utf-8'))
report['resource_label_corrections']=[
 {'score_set':URN,'frozen_manifest_gene_label':'FKRP','frozen_manifest_label_prefix':'00001197',
  'verified_MaveDB_target_gene_label':'LARGE1',
  'evidence':'revision_20261122/dna_audit/mavedb_bare_c_construct_by_set.csv and archived MaveDB target metadata; target-sequence support remains conditional under the documented c.1 alignment assumption.'},
 {'score_set':'urn:mavedb:00001197-a-5','gene_label':'FKRP',
  'pairing_boundary':'The separate historical replication pairing is unresolved and is not used as pair-level evidence in this collection audit.'}
]
report['observed']['group_tables_label_corrected_for_resource_target']=True
p.write_text(json.dumps(report,indent=2),encoding='utf-8')
print('JSON resource label boundary recorded',flush=True)

