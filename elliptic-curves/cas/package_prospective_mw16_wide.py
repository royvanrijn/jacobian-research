#!/usr/bin/env python3
"""Portable source, population, census and result evidence for the wider MW16 pilot."""
import argparse
from hashlib import sha256
from pathlib import Path
import shutil
import zipfile
import certify_compact_r17_candidates as cert
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'

def hashed(path):
    h=sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()

def package(directory,manifest):
    archive=manifest.with_suffix('.zip')
    if manifest.exists() or archive.exists():raise FileExistsError('preserve evidence bundle')
    ledger=cert.read(directory/'point-ledger.json')
    if ledger['status']!='COMPLETE_FIXED_BATCH_ATTEMPTS':raise ArithmeticError('fixed batch unfinished')
    paths=set()
    for pattern in ('*.json','*/generic-census.json','*/population.json','*/scan-*.json','*/preparation-supervision.json',
                    '*/*.log','*/*.supervisor.json','*/candidate-*/result.json','verification/*','previous/*/*.json'):
        paths.update(p for p in directory.glob(pattern) if p.is_file())
    for name in ('protocol.json','point-protocol.json'):
        paths.update(ROOT/p for p in cert.read(directory/name)['sources'])
    paths.update((Path(__file__).resolve(),ART/'prospective_mw16_census_audit_v1.json',ART/'prospective_mw16_wide_results_v1.json',
                  ROOT/'elliptic-curves/cas/replay_prospective_mw16_search.py',ROOT/'elliptic-curves/cas/certify_prospective_mw16_wide_results.py',
                  ROOT/'elliptic-curves/cas/run_prospective_mw16_wide_batch.py',ROOT/'elliptic-curves/cas/audit_prospective_mw16_census.py'))
    paths.update(ROOT/r['result_path'] for r in ledger['rows'])
    parent=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h1024-v1'
    paths.update(parent/name for name in ('point-protocol.json','point-launch.json','point-ledger.json','verification/results.json'))
    paths.update(ROOT/name for name in cert.read(parent/'point-protocol.json')['sources'])
    rows=[]
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(paths):
            name=str(p.relative_to(ROOT));info=zipfile.ZipInfo(name,date_time=(2026,9,5,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o644<<16
            with p.open('rb') as source,z.open(info,'w') as sink:shutil.copyfileobj(source,sink,1<<20)
            rows.append({'path':name,'bytes':p.stat().st_size,'sha256':hashed(p)})
    with zipfile.ZipFile(archive) as z:
        if z.testzip() is not None:raise ArithmeticError('zip integrity failed')
        for r in rows:
            with z.open(r['path']) as stream:
                h=sha256()
                for chunk in iter(lambda:stream.read(1<<20),b''):h.update(chunk)
            if h.hexdigest()!=r['sha256']:raise ArithmeticError('member hash mismatch')
    cert.write(manifest,{'schema':'elliptic-curves.prospective-mw16-wide-evidence.v1','builder_sha256':hashed(Path(__file__).resolve()),
        'archive':str(archive.relative_to(ROOT)),'archive_sha256':hashed(archive),'archive_bytes':archive.stat().st_size,'files':rows,
        'claim_boundary':'Retained fixed population, all generic parity witnesses, point transcripts and exact rank certificates. Prime trace arrays regenerate from frozen source; signed scanner table hashes are retained. No universal novelty or rank upper bounds.'})
    print('PACKAGED WIDER MW16 PILOT',len(rows),'files',archive.stat().st_size,'bytes',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--directory',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);a=p.parse_args();package(a.directory.resolve(),a.manifest.resolve())
