#!/usr/bin/env python3
"""Package the smaller generating-base certificate and exact relation intersection."""
from pathlib import Path
import sys
import zipfile
import extend_small_conductor_norm_batch as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits,run

ROOT, ART, cert = batch.ROOT,batch.ART,batch.cert
D = ROOT/'artifacts/local/elliptic-curves/small-conductor-smaller-base-v2'
PREFIX = 'small_conductor_smaller_base'
ZIP = ART/(PREFIX+'_evidence_v2.zip')
MANIFEST = ART/(PREFIX+'_evidence_v2.json')
REPORT = ART/(PREFIX+'_portable_replay_v2.json')


def package():
    if ZIP.exists() or MANIFEST.exists():
        raise FileExistsError('preserve evidence')
    result = ART/'small_conductor_smaller_base_v2.json'
    if cert.read(result)['status']!='PASS':
        raise ArithmeticError('matrix certificate required')
    base = ART/'small_conductor_special_primes_evidence_v1.zip'
    base_manifest = cert.read(base.with_suffix('.json'))
    if cert.hashed(base)!=base_manifest['zip_sha256']:
        raise ArithmeticError('previous evidence hash differs')
    files = {Path(__file__).resolve(),Path(batch.__file__).resolve(),result,
             ROOT/'elliptic-curves/cas/certify_small_conductor_smaller_base.sage',
             ROOT/'elliptic-curves/cas/certify_small_conductor_smaller_base_v2.sage',
             ROOT/'elliptic-curves/cas/target_small_conductor_prime_ideals.sage',
             ROOT/'elliptic-curves/cas/audit_small_conductor_norm_batch.sage',
             ROOT/'elliptic-curves/cas/audit_bnf_free_s_class_quotient.py',
             ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_DESCENT_SHORTCUT_2026-09-06.md'}
    files.update(p for p in D.rglob('*') if p.is_file())
    files.update(p for p in D.with_name('small-conductor-smaller-base-v1').rglob('*') if p.is_file())
    overrides = {str(p.relative_to(ROOT)):p for p in files}
    entries = []
    with zipfile.ZipFile(ZIP,'x',zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        with zipfile.ZipFile(base) as old:
            for row in base_manifest['files']:
                if row['path'] not in overrides:
                    archive.writestr(row['path'],old.read(row['path']))
                    entries.append(row)
        for name,path in sorted(overrides.items()):
            archive.write(path,name)
            entries.append({'path':name,'bytes':path.stat().st_size,'sha256':cert.hashed(path)})
    checkpoint(MANIFEST,{'schema':'elliptic-curves.small-conductor-smaller-base-evidence.v2',
                         'zip_sha256':cert.hashed(ZIP),'zip_bytes':ZIP.stat().st_size,
                         'files':sorted(entries,key=lambda r:r['path']),
                         'requirements':'Sage with PARI2.17.3, MPFI and Arb. Isolated replay verifies the complete inherited principal-relation matrix, interval explicit-formula generation tests and exact outside-column elimination. Historical scalar searches and the failed V1 API call are retained without rerunning them.'})
    print('PACKAGED',len(entries),'FILES',ZIP.stat().st_size,'BYTES',flush=True)


def replay():
    if REPORT.exists():
        raise FileExistsError('preserve isolated report')
    manifest=cert.read(MANIFEST)
    if cert.hashed(ZIP)!=manifest['zip_sha256']:
        raise ArithmeticError('archive differs')
    destination=D.with_name('small-conductor-smaller-base-isolated-v2')
    destination.mkdir(exist_ok=False)
    with zipfile.ZipFile(ZIP) as archive:
        if any(Path(n).is_absolute() or '..' in Path(n).parts for n in archive.namelist()):
            raise ArithmeticError('unsafe archive path')
        archive.extractall(destination)
    for row in manifest['files']:
        if cert.hashed(destination/row['path'])!=row['sha256']:
            raise ArithmeticError('extracted file differs')
    command=['/home/royvanrijn/.local/bin/sage','-python','elliptic-curves/cas/certify_small_conductor_smaller_base_v2.sage','--check']
    outcome=run(command,limits=Limits(180,1610612736),cwd=destination,
                log_path=destination/'replay.log',checkpoint_path=destination/'replay.supervisor.json')
    passed=outcome['outcome']=='completed' and outcome['returncode']==0
    checkpoint(REPORT,{'schema':'elliptic-curves.small-conductor-smaller-base-isolated.v2',
                       'status':'PASS' if passed else 'FAILED_OR_CENSORED',
                       'zip_sha256':manifest['zip_sha256'],'supervision':outcome})
    print('ISOLATED SMALLER-BASE REPLAY','PASS' if passed else 'FAILED_OR_CENSORED',flush=True)
    if not passed: raise SystemExit(1)


if __name__=='__main__':
    {'package':package,'replay':replay}[sys.argv[1]]()
