"""Verify SHA-256 inventory after download or before publication."""
from pathlib import Path
import hashlib
import sys

root=Path(__file__).resolve().parents[1]
fail=[]
count=0
for line in (root/'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
    digest,name=line.split('  ',1)
    p=root/name
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
        fail.append(name)
    count+=1
print(f'{count} files checked; {len(fail)} failures')
if fail:
    print('\n'.join(fail));sys.exit(1)
