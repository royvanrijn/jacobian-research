#!/usr/bin/env python3
"""Package and replay the general verifier without the large MW16 dependency."""
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys
import zipfile
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT/'artifacts/local/elliptic-curves'
OUT = ART/'class_span_machinery_v1.json'
ZIP = ART/'class_span_machinery_evidence_v1.zip'
MANIFEST = ZIP.with_suffix('.json')
REPORT = ART/'class_span_machinery_portable_replay_v1.json'


def read(p):return json.loads(p.read_text())
def hashed(p):return sha256(p.read_bytes()).hexdigest()


def package():
    if ZIP.exists() or MANIFEST.exists():raise FileExistsError('preserve archive')
    result = read(OUT)
    if result['status']!='PASS':raise ArithmeticError('validation required')
    files = {OUT, Path(__file__).resolve()}
    for name,h in result['sources'].items():
        p = ROOT/name
        if hashed(p)!=h:raise ArithmeticError('frozen source differs')
        files.add(p)
    files.update((ART/'class-span-v1').glob('*.json'))
    for name in ['__init__.py','store.py','supervisor.py']:
        files.add(ROOT/'elliptic-curves/cas/research_runtime'/name)
    entries = []
    with zipfile.ZipFile(ZIP,'x',zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for p in sorted(files):
            name = str(p.relative_to(ROOT))
            archive.write(p,name)
            entries.append({'path':name,'sha256':hashed(p),'bytes':p.stat().st_size})
    checkpoint(MANIFEST, {'schema':'number-fields.class-span-evidence.v1','files':entries,
        'zip_sha256':hashed(ZIP),'zip_bytes':ZIP.stat().st_size,
        'scope':'Self-contained general engine, seven small-field fixtures, exact arithmetic controls, tests, CLI and protocol. MW16 integration has a separate regression and depends on the original curve evidence.',
        'replay':'sage -python elliptic-curves/cas/certify_class_span_machinery.sage check'})
    print('PACKAGED',len(files),'FILES',ZIP.stat().st_size,'BYTES',flush=True)


def replay():
    if REPORT.exists():raise FileExistsError('preserve replay')
    manifest=read(MANIFEST)
    if hashed(ZIP)!=manifest['zip_sha256']:raise ArithmeticError('archive differs')
    destination=LOCAL/'class-span-machinery-isolated-v1'
    destination.mkdir(exist_ok=False)
    with zipfile.ZipFile(ZIP) as archive:
        if any(Path(n).is_absolute() or '..' in Path(n).parts for n in archive.namelist()):
            raise ArithmeticError('unsafe archive path')
        archive.extractall(destination)
    for row in manifest['files']:
        if hashed(destination/row['path'])!=row['sha256']:raise ArithmeticError('extracted file differs')
    sage=shutil.which('sage')
    if not sage:raise RuntimeError('Sage required')
    stages=[]
    commands=[('core',[sage,'-python','elliptic-curves/cas/certify_class_span_machinery.sage','check'])]
    args=['elliptic-curves/cas/verify_class_span_grh.py',
          'artifacts/generated-results/elliptic-curves/class-span-v1/imaginary_c2.input.json',
          '--output','example.certificate.json']
    commands += [('cli-build',[sys.executable,*args]),('cli-replay',[sys.executable,*args,'--check'])]
    for name,command in commands:
        result=run(command,limits=Limits(120,1610612736),cwd=destination,
            log_path=destination/(name+'.log'),checkpoint_path=destination/(name+'.supervisor.json'))
        stages.append({'stage':name,'supervision':result})
        if result['outcome']!='completed' or result['returncode']!=0:break
    passed=len(stages)==3 and all(s['supervision']['outcome']=='completed' and s['supervision']['returncode']==0 for s in stages)
    if passed:
        example=read(destination/'example.certificate.json')
        passed=example['status']=='CERTIFIED_UNDER_GRH' and example['class_two_rank_upper_bound_under_grh']==1
    checkpoint(REPORT,{'schema':'number-fields.class-span-portable-replay.v1',
        'status':'PASS' if passed else 'FAILED_OR_CENSORED','zip_sha256':manifest['zip_sha256'],
        'stages':stages,'scope':manifest['scope']})
    print('GENERAL MACHINERY PORTABLE REPLAY','PASS' if passed else 'FAILED_OR_CENSORED',flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':{'package':package,'replay':replay}[sys.argv[1]]()
