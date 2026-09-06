#!/usr/bin/env python3
"""Package audited class-target waves with the complete inherited evidence."""
from pathlib import Path
import sys
import zipfile
import extend_small_conductor_norm_batch as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits,run

ROOT, ART, cert = batch.ROOT,batch.ART,batch.cert
D = ROOT/'artifacts/local/elliptic-curves/small-conductor-class-target-v1'
PREFIX = 'small_conductor_class_target'
ZIP = ART/(PREFIX+'_evidence_v1.zip')
MANIFEST = ART/(PREFIX+'_evidence_v1.json')
REPORT = ART/(PREFIX+'_portable_replay_v1.json')


def package():
    if ZIP.exists() or MANIFEST.exists():
        raise FileExistsError('preserve evidence')
    audits = sorted(ART.glob('small_conductor_class_target_wave_*_v1.json'))
    if not audits or any(cert.read(p)['status']!='PASS' for p in audits):
        raise ArithmeticError('audited waves required')
    last = cert.read(audits[-1])['wave']
    base = ART/'small_conductor_small_base_targets_evidence_v1.zip'
    baseline = cert.read(base.with_suffix('.json'))
    if cert.hashed(base)!=baseline['zip_sha256']:
        raise ArithmeticError('baseline archive differs')
    files = {Path(__file__).resolve(),ROOT/'elliptic-curves/cas/pursue_small_conductor_class_target.sage',
             ROOT/'elliptic-curves/cas/record_small_conductor_class_target.py',
             ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_CLASS_TARGET_2026-09-06.md',*audits}
    for wave in range(1,last+1):
        files.update(p for p in (D/('wave_%03d'%wave)).rglob('*') if p.is_file())
    overrides = {str(p.relative_to(ROOT)):p for p in files}
    entries = []
    with zipfile.ZipFile(ZIP,'x',zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        with zipfile.ZipFile(base) as old:
            for row in baseline['files']:
                if row['path'] not in overrides:
                    archive.writestr(row['path'],old.read(row['path']))
                    entries.append(row)
        for name,path in sorted(overrides.items()):
            archive.write(path,name)
            entries.append({'path':name,'bytes':path.stat().st_size,'sha256':cert.hashed(path)})
    checkpoint(MANIFEST,{'schema':'elliptic-curves.small-conductor-class-target-evidence.v1',
        'zip_sha256':cert.hashed(ZIP),'zip_bytes':ZIP.stat().st_size,'last_wave':last,
        'files':sorted(entries,key=lambda r:r['path']),
        'requirements':'Sage with PARI2.17.3, MPFI and Arb. Checks the complete inherited proof chain and every new principal-ideal witness, selection and exact matrix transition. Rejected-candidate enumeration is not part of the rank certificate.'})
    print('PACKAGED',len(entries),'FILES',ZIP.stat().st_size,'BYTES',flush=True)


def replay():
    if REPORT.exists():
        raise FileExistsError('preserve isolated report')
    manifest = cert.read(MANIFEST)
    if cert.hashed(ZIP)!=manifest['zip_sha256']:
        raise ArithmeticError('archive differs')
    destination = D.with_name('small-conductor-class-target-isolated-v1')
    destination.mkdir(exist_ok=False)
    with zipfile.ZipFile(ZIP) as archive:
        if any(Path(n).is_absolute() or '..' in Path(n).parts for n in archive.namelist()):
            raise ArithmeticError('unsafe archive path')
        archive.extractall(destination)
    for row in manifest['files']:
        if cert.hashed(destination/row['path'])!=row['sha256']:
            raise ArithmeticError('extracted file differs')
    command = ['/home/royvanrijn/.local/bin/sage','-python','elliptic-curves/cas/pursue_small_conductor_class_target.sage','check','--wave',str(manifest['last_wave'])]
    outcome = run(command,limits=Limits(600,1610612736),cwd=destination,
                  log_path=destination/'replay.log',checkpoint_path=destination/'replay.supervisor.json')
    passed = outcome['outcome']=='completed' and outcome['returncode']==0
    checkpoint(REPORT,{'schema':'elliptic-curves.small-conductor-class-target-isolated.v1',
        'status':'PASS' if passed else 'FAILED_OR_CENSORED','last_wave':manifest['last_wave'],
        'zip_sha256':manifest['zip_sha256'],'supervision':outcome})
    print('ISOLATED CLASS-TARGET REPLAY','PASS' if passed else 'FAILED_OR_CENSORED',flush=True)
    if not passed:
        raise SystemExit(1)


if __name__=='__main__':
    {'package':package,'replay':replay}[sys.argv[1]]()
