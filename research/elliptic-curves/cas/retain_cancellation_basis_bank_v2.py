#!/usr/bin/env sage -python
"""Complete the import archive; retain the first portability failure unchanged."""
import argparse
import gzip
import io
import json
from pathlib import Path
import resource
import subprocess
import tarfile
import tempfile
import time

import cancellation_basis_early as early
from cancellation_basis_bank_cost_v2 import guard, OUT, RAW
from cancellation_basis_epoch import read, sha, need
from cancellation_scheduler_fresh import source_closure
from cancellation_scheduler_prepare import new_write
from finite_cancellation_corpus import ROOT, LOCAL, canonical, digest

CAS = Path(__file__).resolve().parent
PRE = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v1'
PRE_RAW = LOCAL/'cancellation-basis-bank-certificate-v1'


def pack():
    plan,_ = guard(); need(read(OUT/'audit.json')['status'] == 'PASS_RETAINED_BANK_COST_ACCOUNTING','accounting absent')
    need(not (OUT/'replay-bundle-v2.tar.gz').exists(),'preserve replay bundle')
    import sys
    import run_complement_seed_v3
    loaded = [Path(m.__file__).resolve() for m in list(sys.modules.values()) if getattr(m,'__file__',None) and Path(m.__file__).resolve().is_relative_to(ROOT)]
    dependencies = source_closure([Path(__file__),CAS/'verify_cancellation_basis_bank_cost.py']+loaded)
    need(read(OUT/'portable-check.json')['status'] == 'RETAINED_PORTABLE_FAILURE','first portability failure missing')
    new_write(OUT/'dependency-audit-v2.json',{'status':'BOUND_POST_PROTOCOL_ACCOUNTING_AND_RETENTION',
        'primary_protocol_sha256':sha(OUT/'protocol.json'),
        'extra_source_sha256':{k:v for k,v in dependencies.items() if k not in plan['source_sha256']},
        'failed_portable_check_sha256':sha(OUT/'portable-check.json'),
        'loaded_local_modules_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(set(loaded))},
        'boundary':'Post-protocol import closure repair. The original archive omitted the ecsearch package reached through certify_compact_r17_candidates. Existing current files are added and hashed now; they are not retrospectively called part of the121-source primary seal. Primary arithmetic and timing outputs are unchanged; the first archive and failed portable check remain retained.'})
    paths = {ROOT/p for p in set(plan['source_sha256']) | set(dependencies) | set(plan['retained_input_sha256'])}
    for source in (PRE/'protocol.json',early.OUT/'protocol.json'):
        old = read(source); paths.add(source)
        paths.update(ROOT/p for p in old['source_sha256'])
        paths.update(source.parent/p for p in old['input_sha256'])
        paths.update(ROOT/p for p in old.get('retained_input_sha256',{}))
    for folder in (PRE,OUT,PRE_RAW,RAW):
        paths.update(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts and
                     p.suffix not in ('.pyc','.lock') and not p.name.endswith('.tar.gz') and p.name not in ('replay-bundle-v2.tar.gz','replay-manifest-v2.json','pack-v2.json','portable-check-v2.json','completion.json'))
    files = []
    for path in sorted(paths):
        need(path.is_file(),'missing retained file: '+str(path))
        files.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),'bytes':path.stat().st_size})
    manifest = {'status':'COMPLETE_BANK_COST_REPLAY_INPUTS','protocol_sha256':sha(OUT/'protocol.json'),
        'files':files,'file_count':len(files),
        'boundary':'Both timing versions, commissioning source snapshots and payloads, all retained bank inputs, all source closures, and every paired output. Includes early-protocol guard inputs, not the entire earlier point-search campaign. No missing file is reconstructed.'}
    new_write(OUT/'replay-manifest-v2.json',manifest)
    with (OUT/'replay-bundle-v2.tar.gz').open('xb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,compresslevel=6) as compressed:
            with tarfile.open(fileobj=compressed,mode='w|') as archive:
                for row in files+[{'path':str((OUT/'replay-manifest-v2.json').relative_to(ROOT))}]:
                    data = (ROOT/row['path']).read_bytes(); info = tarfile.TarInfo('research/'+row['path'])
                    info.size = len(data); info.mode = 0o644; info.mtime = 0
                    archive.addfile(info,io.BytesIO(data))
    result = {'status':'PASS_REPLAY_PACKAGING','protocol_sha256':sha(OUT/'protocol.json'),
        'manifest_sha256':sha(OUT/'replay-manifest-v2.json'),'archive_sha256':sha(OUT/'replay-bundle-v2.tar.gz'),
        'files':len(files),'compressed_bytes':(OUT/'replay-bundle-v2.tar.gz').stat().st_size}
    new_write(OUT/'pack-v2.json',result); print(json.dumps(result),flush=True)


def portable():
    packaging = read(OUT/'pack-v2.json'); manifest = read(OUT/'replay-manifest-v2.json')
    need(sha(OUT/'replay-bundle-v2.tar.gz') == packaging['archive_sha256'] and
         sha(OUT/'replay-manifest-v2.json') == packaging['manifest_sha256'],'archive changed')
    before = resource.getrusage(resource.RUSAGE_CHILDREN); start = time.process_time()
    with tempfile.TemporaryDirectory(prefix='cancellation-bank-portable-') as temporary:
        target = Path(temporary)
        with tarfile.open(OUT/'replay-bundle-v2.tar.gz','r:gz') as archive:
            for entry in archive:
                relative = Path(entry.name)
                need(entry.isfile() and not relative.is_absolute() and '..' not in relative.parts and
                     relative.parts[0] == 'research','invalid archive entry')
                path = target/relative; path.parent.mkdir(parents=True,exist_ok=True)
                path.write_bytes(archive.extractfile(entry).read())
        extracted = target/'research'
        for row in manifest['files']:
            path = extracted/row['path']; need(path.stat().st_size == row['bytes'] and sha(path) == row['sha256'],'empty-root bytes differ')
        code = ('import sys, resource; resource.setrlimit(resource.RLIMIT_CPU,(25,30)); '
                'sys.path.insert(0,'+repr(str(extracted/'elliptic-curves/cas'))+'); '
                'import cancellation_basis_bank_cost_v2 as cost; cost.guard(); '
                'import verify_cancellation_basis_bank_cost as audit; audit.verify(write_result=False)')
        run = subprocess.run(['sage','-python','-c',code],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=60,cwd=target)
        (RAW/'portable-check-v2.log').write_text(run.stdout)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        result = {'status':'PASS_EMPTY_ROOT_BANK_COST_CHECK' if run.returncode == 0 else 'RETAINED_PORTABLE_FAILURE',
            'returncode':run.returncode,'files':len(manifest['files']),'archive_sha256':packaging['archive_sha256'],
            'manifest_sha256':packaging['manifest_sha256'],'log_sha256':sha(RAW/'portable-check-v2.log'),
            'checker_sha256':sha(CAS/'verify_cancellation_basis_bank_cost.py'),
            'full_child_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            'packaging_check_cpu_seconds':time.process_time()-start,
            'boundary':'All bytes checked in an empty root, followed by source/input guards and independent retained-output integer-norm, parity, timing and bootstrap accounting. No new CVP, native rank, map or point-search replay.'}
        new_write(OUT/'portable-check-v2.json',result)
        need(run.returncode == 0,'portable check failed')
    completion = {'status':'RETAINED_COMPLETE_BANK_COST_COMPARISON','engineering_gate':read(OUT/'summary.json')['engineering_gate'],
        'protocol_sha256':sha(OUT/'protocol.json'),'result_sha256':{name:sha(OUT/name) for name in
            ('summary.json','audit.json','supervision.json','dependency-audit-v2.json','pack-v2.json','replay-manifest-v2.json','portable-check-v2.json')},
        'archive_sha256':packaging['archive_sha256'],'point_search_calls':0,
        'boundary':'Finite implementation-cost result retained; sustained later-gain amplification and external comparison remain open.'}
    new_write(OUT/'completion.json',completion); print(json.dumps(result),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command',choices=['pack','portable']); args = parser.parse_args()
    globals()[args.command]()
