#!/usr/bin/env python3
"""Preserve fixed class-block arithmetic and oracle-labelled solubility checkpoints."""
from pathlib import Path
import zipfile
import retrospective as r
NAMES=['constructed-class-half-ideal-v1','constructed-class-direct-artin-v1',
       'curve302-maximal-norm-form-v1','constructed-class-compaction-v1',
       'constructed-class-oracle-solubility-v1']
ARCHIVE=r.OUT/'rank_jump_class_block_lift_evidence_v1.zip'
MANIFEST=r.OUT/'rank_jump_class_block_lift_evidence_v1.json'

if __name__=='__main__':
    files=[]
    for name in NAMES:
        d=r.ROOT/('artifacts/local/rank-jump-'+name)
        files += [p for p in sorted(d.glob('*')) if p.is_file() and p.suffix in ['.json','.log']]
    assert r.read(r.OUT/'rank_jump_constructed_cover_points_verification_v1.json')['status']=='PASS'
    entries=[]
    with zipfile.ZipFile(ARCHIVE,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for p in files:
            name=str(p.relative_to(r.ROOT));raw=p.read_bytes();info=zipfile.ZipInfo(name,date_time=(2026,9,7,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
            archive.writestr(info,raw,compresslevel=9);entries.append({'path':name,'sha256':r.digest(raw),'bytes':len(raw),
                'oracle_evaluation': 'oracle-solubility' in name})
    r.write_new(MANIFEST,{'schema':'rank-jump.class-block-lift-evidence.v1','files':entries,
        'archive_sha256':r.digest(ARCHIVE.read_bytes()),'producer_sha256':r.digest(Path(__file__).read_bytes()),
        'restore':'Extract into an empty replay checkout; do not overwrite live data. Earlier constructive archives remain separate.',
        'boundary':'Includes both point-free arithmetic and explicitly oracle-labelled retrospective solubility outputs. Oracle files must not feed prospective selectors.'})
    with zipfile.ZipFile(ARCHIVE) as archive:assert all(r.digest(archive.read(row['path']))==row['sha256'] for row in entries)
    print('PASS packaged',len(entries),'files',ARCHIVE.stat().st_size,'bytes')
