#!/usr/bin/env python3
"""Package every audited class-target stage and the character completion proof."""
from pathlib import Path
import sys
import zipfile
import extend_small_conductor_norm_batch as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits,run

ROOT,ART,cert=batch.ROOT,batch.ART,batch.cert
PREFIX='small_conductor_class_target'
ZIP=ART/(PREFIX+'_evidence_v1.zip')
MANIFEST=ART/(PREFIX+'_evidence_v1.json')
REPORT=ART/(PREFIX+'_portable_replay_v1.json')
COMPLETION=ART/'small_conductor_class_completion_v1.json'
KINDS=[('box','','small-conductor-class-target-v1',''),
       ('strip','strip_','small-conductor-class-target-strips-v1','_strips'),
       ('protected','protected_','small-conductor-class-target-protected-v1','_protected'),
       ('capped','capped_','small-conductor-class-target-capped-v1','_capped'),
       ('residual','residual_','small-conductor-class-target-residual-v1','_residual')]


def package():
    if ZIP.exists() or MANIFEST.exists():raise FileExistsError('preserve archive')
    files={Path(__file__).resolve(),ROOT/'elliptic-curves/cas/record_small_conductor_class_target.py',
           ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_CLASS_TARGET_2026-09-06.md'}
    last=None
    for family,prefix,directory,suffix in KINDS:
        audits=sorted(ART.glob('small_conductor_class_target_'+prefix+'wave_*_v1.json'))
        if not audits:continue
        files.add(ROOT/('elliptic-curves/cas/pursue_small_conductor_class_target'+suffix+'.sage'))
        for index,p in enumerate(audits,1):
            a=cert.read(p)
            if a['status']!='PASS' or a['wave']!=index:raise ArithmeticError('incomplete audit sequence')
            files.add(p)
            d=ROOT/'artifacts/local/elliptic-curves'/directory/('wave_%03d'%index)
            files.update(p for p in d.rglob('*') if p.is_file())
        last=(family,len(audits))
    if not last:raise ArithmeticError('no audited waves')
    family,wave=last
    combined=ART/('small_conductor_class_target_combined_%s_%03d_v1.json'%(family,wave))
    if family not in ['capped','residual'] and cert.read(combined)['status']!='PASS':raise ArithmeticError('combined proof required')
    files.update(ART.glob('small_conductor_class_target_combined_*_v1.json'))
    for name in ['certify_small_conductor_class_target_combined.sage','certify_small_conductor_class_lower16.sage','certify_small_conductor_class_characters.sage','capped_primorial_remainders.py','benchmark_capped_primorial_remainders.py']:
        files.add(ROOT/'elliptic-curves/cas'/name)
    files.update([ART/'small_conductor_capped_remainders_benchmark_v1.json',ART/'small_conductor_class_lower16_v1.json',ART/'small_conductor_class_characters_v1.json'])
    completion=cert.read(COMPLETION)
    if completion['status']!='PASS' or family!='residual' or completion['residual_wave']!=wave:
        raise ArithmeticError('matching completion certificate required')
    for name,h in completion['sources'].items():
        path=ROOT/name
        if cert.hashed(path)!=h:raise ArithmeticError('completion source differs')
        files.add(path)
    files.add(COMPLETION)
    files.update(p for p in (ROOT/'artifacts/local/elliptic-curves/small-conductor-class-completion-v1').rglob('*') if p.is_file())
    base=ART/'small_conductor_small_base_targets_evidence_v1.zip';baseline=cert.read(base.with_suffix('.json'))
    if cert.hashed(base)!=baseline['zip_sha256']:raise ArithmeticError('baseline archive differs')
    overrides={str(p.relative_to(ROOT)):p for p in files};entries=[]
    with zipfile.ZipFile(ZIP,'x',zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        with zipfile.ZipFile(base) as old:
            for row in baseline['files']:
                if row['path'] not in overrides:
                    archive.writestr(row['path'],old.read(row['path']));entries.append(row)
        for name,path in sorted(overrides.items()):
            archive.write(path,name)
            entries.append({'path':name,'bytes':path.stat().st_size,'sha256':cert.hashed(path)})
    checkpoint(MANIFEST,{'schema':'elliptic-curves.small-conductor-class-target-evidence.v1',
        'family':family,'last_wave':wave,'zip_sha256':cert.hashed(ZIP),'zip_bytes':ZIP.stat().st_size,
        'completion_certificate_sha256':cert.hashed(COMPLETION),
        'files':sorted(entries,key=lambda r:r['path']),
        'requirements':'Sage with PARI2.17.3, MPFI and Arb. Replay reconstructs all inherited proofs and norm waves, the lower-bound and anchor-character proofs, adaptive target skips, point-derived principal parity relations, exact prime membership classifications and both interval-certified character exclusion tests. No exhaustive rejected-candidate replay is claimed.'})
    print('PACKAGED',len(entries),'FILES',ZIP.stat().st_size,'BYTES',flush=True)


def replay():
    if REPORT.exists():raise FileExistsError('preserve replay')
    manifest=cert.read(MANIFEST)
    if cert.hashed(ZIP)!=manifest['zip_sha256']:raise ArithmeticError('archive differs')
    destination=ROOT/'artifacts/local/elliptic-curves/small-conductor-class-target-isolated-v1'
    destination.mkdir(exist_ok=False)
    with zipfile.ZipFile(ZIP) as archive:
        if any(Path(n).is_absolute() or '..' in Path(n).parts for n in archive.namelist()):raise ArithmeticError('unsafe archive path')
        archive.extractall(destination)
    for row in manifest['files']:
        if cert.hashed(destination/row['path'])!=row['sha256']:raise ArithmeticError('extracted file differs')
    if manifest.get('completion_certificate_sha256'):
        if cert.hashed(destination/COMPLETION.relative_to(ROOT))!=manifest['completion_certificate_sha256']:
            raise ArithmeticError('completion certificate differs')
        arguments=['elliptic-curves/cas/certify_small_conductor_class_completion.sage','check']
    elif manifest['family'] in ['capped','residual']:
        arguments=['elliptic-curves/cas/pursue_small_conductor_class_target_'+manifest['family']+'.sage','check','--wave',str(manifest['last_wave'])]
    else:
        arguments=['elliptic-curves/cas/certify_small_conductor_class_target_combined.sage','--family',manifest['family'],'--wave',str(manifest['last_wave']),'--check']
    command=['/home/royvanrijn/.local/bin/sage','-python',*arguments]
    outcome=run(command,limits=Limits(600,1610612736),cwd=destination,
                log_path=destination/'replay.log',checkpoint_path=destination/'replay.supervisor.json')
    passed=outcome['outcome']=='completed' and outcome['returncode']==0
    benchmark=run(['/home/royvanrijn/.local/bin/sage','-python','elliptic-curves/cas/benchmark_capped_primorial_remainders.py','--check'],
                  limits=Limits(60,1610612736),cwd=destination,
                  log_path=destination/'benchmark.log',checkpoint_path=destination/'benchmark.supervisor.json')
    passed=passed and benchmark['outcome']=='completed' and benchmark['returncode']==0
    checkpoint(REPORT,{'schema':'elliptic-curves.small-conductor-class-target-isolated.v1',
        'status':'PASS' if passed else 'FAILED_OR_CENSORED','family':manifest['family'],'last_wave':manifest['last_wave'],
        'completion_certificate_sha256':manifest.get('completion_certificate_sha256'),
        'zip_sha256':manifest['zip_sha256'],'supervision':outcome,'benchmark_supervision':benchmark})
    print('ISOLATED CLASS-TARGET REPLAY','PASS' if passed else 'FAILED_OR_CENSORED',flush=True)
    if not passed:raise SystemExit(1)


if __name__=='__main__':{'package':package,'replay':replay}[sys.argv[1]]()
