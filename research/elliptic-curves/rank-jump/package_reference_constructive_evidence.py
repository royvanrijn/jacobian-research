#!/usr/bin/env python3
"""Preserve raw bounded-run checkpoints without binaries or upstream source copies."""
from pathlib import Path
import zipfile
import retrospective as r

NAMES=['reference-strict-class-construction','seeded-reference-class','seeded-reference-relations',
       'seeded-reference-relations-memory','reference-partial-class-extraction',
       'reference-class-targeted-relations','reference-targeted-class-extraction']
ARCHIVE=r.OUT/'rank_jump_reference_constructive_evidence_v1.zip'
MANIFEST=r.OUT/'rank_jump_reference_constructive_evidence_v1.json'

if __name__=='__main__':
    files=[]
    for name in NAMES:
        d=r.ROOT/('artifacts/local/rank-jump-'+name+'-v1')
        files += [p for p in sorted(d.glob('*')) if p.is_file() and p.suffix in ['.json','.log','.gp']]
    assert len(files)>512
    entries=[]
    with zipfile.ZipFile(ARCHIVE,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for p in files:
            name=str(p.relative_to(r.ROOT));raw=p.read_bytes();z.writestr(name,raw)
            entries.append({'path':name,'sha256':r.digest(raw),'bytes':len(raw)})
    r.write_new(MANIFEST,{'schema':'rank-jump.reference-constructive-evidence.v1','files':entries,
        'archive_sha256':r.digest(ARCHIVE.read_bytes()),'producer_sha256':r.digest(Path(__file__).read_bytes()),
        'restore':'Extract into an empty replay checkout; do not overwrite live checkpoints. No search is required to restore these inputs.',
        'boundary':'Raw inputs, logs and terminal records for completed reference construction attempts. No additional class is certified by packaging.'})
    with zipfile.ZipFile(ARCHIVE) as z:
        assert all(r.digest(z.read(row['path']))==row['sha256'] for row in entries)
    print('PASS packaged',len(entries),'checkpoint files',flush=True)
