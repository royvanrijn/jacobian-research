#!/usr/bin/env python3
"""Archive terminal302 arithmetic and bisection constructor checkpoints."""
from pathlib import Path
import zipfile
import retrospective as r
NAMES=['curve302-arithmetic-unit-class-v1','curve302-conic-triple-bisection-v1']
ARCHIVE=r.OUT/'rank_jump_curve302_transfer_attempts_evidence_v1.zip'
MANIFEST=r.OUT/'rank_jump_curve302_transfer_attempts_evidence_v1.json'

if __name__=='__main__':
    files=[]
    for name in NAMES:
        d=r.ROOT/('artifacts/local/rank-jump-'+name)
        files += [p for p in sorted(d.glob('*')) if p.is_file() and p.suffix in ['.json','.log']]
    entries=[]
    with zipfile.ZipFile(ARCHIVE,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for p in files:
            name=str(p.relative_to(r.ROOT));raw=p.read_bytes();info=zipfile.ZipInfo(name,date_time=(2026,9,7,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            archive.writestr(info,raw,compresslevel=9);entries.append({'path':name,'sha256':r.digest(raw),'bytes':len(raw)})
    r.write_new(MANIFEST,{'schema':'rank-jump.curve302-transfer-attempts-evidence.v1','files':entries,
        'archive_sha256':r.digest(ARCHIVE.read_bytes()),'producer_sha256':r.digest(Path(__file__).read_bytes()),
        'restore':'Extract into an empty replay checkout, never over live checkpoints.',
        'boundary':'Terminal equation/generic-only constructor evidence. These attempts produced no additional302 strict class.'})
    with zipfile.ZipFile(ARCHIVE) as archive:assert all(r.digest(archive.read(row['path']))==row['sha256'] for row in entries)
    print('PASS packaged',len(entries),'files',ARCHIVE.stat().st_size,'bytes')
