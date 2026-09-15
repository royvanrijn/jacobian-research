#!/usr/bin/env sage -python
"""Retain both cost versions and verify the successor in an empty root."""
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
    need(not (OUT/'replay-bundle.tar.gz').exists(),'preserve replay bundle')
    dependencies = source_closure([Path(__file__),CAS/'verify_cancellation_basis_bank_cost.py'])
    new_write(OUT/'dependency-audit.json',{'status':'BOUND_POST_PROTOCOL_ACCOUNTING_AND_RETENTION',
        'primary_protocol_sha256':sha(OUT/'protocol.json'),
        'extra_source_sha256':{k:v for k,v in dependencies.items() if k not in plan['source_sha256']},
        'boundary':'Additional accounting and packaging only; no primary arithmetic, output, gate or source seal changes.'})
    paths = {ROOT/p for p in set(plan['source_sha256']) | set(dependencies) | set(plan['retained_input_sha256'])}
    for source in (PRE/'protocol.json',early.OUT/'protocol.json'):
        old = read(source); paths.add(source)
        paths.update(ROOT/p for p in old['source_sha256'])
        paths.update(source.parent/p for p in old['input_sha256'])
        paths.update(ROOT/p for p in old.get('retained_input_sha256',{}))
    for folder in (PRE,OUT,PRE_RAW,RAW):
        paths.update(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts and
                     p.suffix not in ('.pyc','.lock') and p.name not in ('replay-bundle.tar.gz','replay-manifest.json','pack.json','portable-check.json','completion.json'))
    files = []
    for path in sorted(paths):
        need(path.is_file(),'missing retained file: '+str(path))
        files.append({'path':str(path.relative_to(ROOT)),'sha256':sha(path),'bytes':path.stat().st_size})
    manifest = {'status':'COMPLETE_BANK_COST_REPLAY_INPUTS','protocol_sha256':sha(OUT/'protocol.json'),
        'files':files,'file_count':len(files),
        'boundary':'Both timing versions, commissioning source snapshots and payloads, all retained bank inputs, all source closures, and every paired output. Includes early-protocol guard inputs, not the entire earlier point-search campaign. No missing file is reconstructed.'}
    new_write(OUT/'replay-manifest.json',manifest)
    with (OUT/'replay-bundle.tar.gz').open('xb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,compresslevel=6) as compressed:
            with tarfile.open(fileobj=compressed,mode='w|') as archive:
                for row in files+[{'path':str((OUT/'replay-manifest.json').relative_to(ROOT))}]:
                    data = (ROOT/row['path']).read_bytes(); info = tarfile.TarInfo('research/'+row['path'])
                    info.size = len(data); info.mode = 0o644; info.mtime = 0
                    archive.addfile(info,io.BytesIO(data))
    result = {'status':'PASS_REPLAY_PACKAGING','protocol_sha256':sha(OUT/'protocol.json'),
        'manifest_sha256':sha(OUT/'replay-manifest.json'),'archive_sha256':sha(OUT/'replay-bundle.tar.gz'),
        'files':len(files),'compressed_bytes':(OUT/'replay-bundle.tar.gz').stat().st_size}
    new_write(OUT/'pack.json',result); print(json.dumps(result),flush=True)


def portable():
    packaging = read(OUT/'pack.json'); manifest = read(OUT/'replay-manifest.json')
    need(sha(OUT/'replay-bundle.tar.gz') == packaging['archive_sha256'] and
         sha(OUT/'replay-manifest.json') == packaging['manifest_sha256'],'archive changed')
    before = resource.getrusage(resource.RUSAGE_CHILDREN); start = time.process_time()
    with tempfile.TemporaryDirectory(prefix='cancellation-bank-portable-') as temporary:
        target = Path(temporary)
        with tarfile.open(OUT/'replay-bundle.tar.gz','r:gz') as archive:
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
        (RAW/'portable-check.log').write_text(run.stdout)
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        result = {'status':'PASS_EMPTY_ROOT_BANK_COST_CHECK' if run.returncode == 0 else 'RETAINED_PORTABLE_FAILURE',
            'returncode':run.returncode,'files':len(manifest['files']),'archive_sha256':packaging['archive_sha256'],
            'manifest_sha256':packaging['manifest_sha256'],'log_sha256':sha(RAW/'portable-check.log'),
            'checker_sha256':sha(CAS/'verify_cancellation_basis_bank_cost.py'),
            'full_child_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            'packaging_check_cpu_seconds':time.process_time()-start,
            'boundary':'All bytes checked in an empty root, followed by source/input guards and independent retained-output integer-norm, parity, timing and bootstrap accounting. No new CVP, native rank, map or point-search replay.'}
        new_write(OUT/'portable-check.json',result)
        need(run.returncode == 0,'portable check failed')
    completion = {'status':'RETAINED_COMPLETE_BANK_COST_COMPARISON','engineering_gate':read(OUT/'summary.json')['engineering_gate'],
        'protocol_sha256':sha(OUT/'protocol.json'),'result_sha256':{name:sha(OUT/name) for name in
            ('summary.json','audit.json','supervision.json','dependency-audit.json','pack.json','replay-manifest.json','portable-check.json')},
        'archive_sha256':packaging['archive_sha256'],'point_search_calls':0,
        'boundary':'Finite implementation-cost result retained; sustained later-gain amplification and external comparison remain open.'}
    new_write(OUT/'completion.json',completion); print(json.dumps(result),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command',choices=['pack','portable']); args = parser.parse_args()
    globals()[args.command]()
