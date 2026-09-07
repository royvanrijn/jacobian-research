#!/usr/bin/env python3
"""Archive terminal positive-construction checkpoints, including preserved failures."""
from pathlib import Path
import zipfile
import retrospective as r

NAMES=['reference-adaptive-class-completion-v1','reference-small-representative-class-v1',
       'reference-large-support-extraction-v1','reference-additional-strict-verification-v1',
       'generic-only-class-anchors-v1','curve302-strict-constructor-arithmetic-v1',
       'curve302-section-star-v1','curve302-section-star-v2']
ARCHIVE=r.OUT/'rank_jump_positive_class_evidence_v1.zip'
MANIFEST=r.OUT/'rank_jump_positive_class_evidence_v1.json'

if __name__=='__main__':
    files=[]
    for name in NAMES:
        d=r.ROOT/('artifacts/local/rank-jump-'+name)
        files += [p for p in sorted(d.glob('*')) if p.is_file() and p.suffix in ['.json','.log','.gp']]
    assert r.read(r.OUT/'rank_jump_reference_additional_strict_verification_v1.json')['status']=='PASS'
    entries=[]
    with zipfile.ZipFile(ARCHIVE,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for path in files:
            name=str(path.relative_to(r.ROOT));raw=path.read_bytes();info=zipfile.ZipInfo(name,date_time=(2026,9,7,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            archive.writestr(info,raw,compresslevel=9);entries.append({'path':name,'sha256':r.digest(raw),'bytes':len(raw)})
    r.write_new(MANIFEST,{'schema':'rank-jump.positive-class-evidence.v1','files':entries,
        'archive_sha256':r.digest(ARCHIVE.read_bytes()),'producer_sha256':r.digest(Path(__file__).read_bytes()),
        'restore':'Extract into an empty replay checkout. Supplement with the previous reference constructive archive. Never overwrite live checkpoints.',
        'boundary':'Raw bounded-construction, arithmetic and geometric checkpoints, including failures. The separate independent verifier certifies two additional reference strict classes; packaging makes no further mathematical claim.'})
    with zipfile.ZipFile(ARCHIVE) as archive:assert all(r.digest(archive.read(row['path']))==row['sha256'] for row in entries)
    print('PASS packaged',len(entries),'files',ARCHIVE.stat().st_size,'bytes')
