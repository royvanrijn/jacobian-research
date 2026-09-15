#!/usr/bin/env sage -python
"""Package the completed transcript integration without repeating arithmetic."""
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

from check_cancellation_basis_certificate import OUT, RAW, COST, guard
from cancellation_basis_epoch import read, sha, need
from cancellation_scheduler_prepare import new_write
from cancellation_scheduler_fresh import source_closure
from finite_cancellation_corpus import ROOT, digest


def inspect():
    plan=guard(); summary=read(OUT/'summary.json'); supervision=read(OUT/'supervision.json')
    need(summary['status']=='PASS_COMPLETE_CERTIFICATE_ENGINE_INTEGRATION' and supervision['status']=='COMPLETE','integration incomplete')
    need(summary['protocol_sha256']==supervision['protocol_sha256']==sha(OUT/'protocol.json'),'protocol differs')
    need(sha(OUT/'negative-checks.json')==summary['negative_checks_sha256'],'negative receipt changed')
    rows=[]
    for job in plan['jobs']:
        folder=RAW/'jobs'/job['id']; result=read(folder/'result.json'); proof=read(folder/'independent-verification.json')
        receipt=read(folder/'integration.json'); row=next(r for r in summary['jobs'] if r['job']==job['id'])
        need(receipt==row,'integration row changed')
        need(sha(folder/'result.json')==row['result_sha256'] and sha(folder/'independent-verification.json')==row['verification_sha256'],'integration output changed')
        need(proof['status']=='PASS' and result['independent_verification_sha256']==row['verification_sha256'],'final proof differs')
        need(sha(folder/'events.json')==result['events_sha256'] and sha(folder/'final-rank.json')==result['final_rank_sha256'],'events or final subgroup changed')
        for call in result['calls']:need(sha(folder/call['file'])==call['sha256'],'retained point call changed')
        rows.append({'job':job['id'],'calls':len(result['calls']),'directions':result['new_directions']})
    need(all(r['status']=='REJECTED_AS_REQUIRED' for r in read(OUT/'negative-checks.json')['rows']),'corruption check missing')
    return {'status':'PASS_RETAINED_INTEGRATION_BINDINGS','jobs':rows,'point_search_calls':0}


def pack():
    plan=guard(); audit=inspect()
    need(not (OUT/'replay-bundle.tar.gz').exists(),'preserve bundle')
    dependencies=source_closure([Path(__file__)])
    new_write(OUT/'dependency-audit.json',{'primary_protocol_sha256':sha(OUT/'protocol.json'),
        'extra_source_sha256':{n:h for n,h in dependencies.items() if n not in plan['source_sha256']},
        'boundary':'Post-integration byte-retention entry point only; no change to the155-source execution seal.'})
    new_write(OUT/'retention-audit.json',audit)
    paths={ROOT/n for n in set(plan['source_sha256'])|set(plan['retained_input_sha256'])|set(dependencies)}
    paths.add(COST/'protocol.json')
    for folder in (OUT,RAW):
        paths.update(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ('.pyc','.lock'))
    files={str(p.relative_to(ROOT)):sha(p) for p in sorted(paths)}
    with (OUT/'replay-bundle.tar.gz').open('xb') as raw:
        with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as zipped:
            with tarfile.open(fileobj=zipped,mode='w|') as archive:
                for name in files:
                    data=(ROOT/name).read_bytes(); item=tarfile.TarInfo(name); item.size=len(data); item.mode=0o644; item.mtime=0
                    archive.addfile(item,io.BytesIO(data))
    manifest={'status':'RETAINED_COMPLETE_INTEGRATION_INPUTS','files_sha256':files,
        'archive_sha256':sha(OUT/'replay-bundle.tar.gz'),'protocol_sha256':sha(OUT/'protocol.json'),
        'boundary':'All four integrations, three updated-hash corruptions, original retained point transcripts, input guards and source snapshots. No new point search or arithmetic replay.'}
    new_write(OUT/'replay-manifest.json',manifest)
    print(json.dumps({'status':manifest['status'],'files':len(files),'bytes':(OUT/'replay-bundle.tar.gz').stat().st_size}),flush=True)


def portable():
    manifest=read(OUT/'replay-manifest.json'); archive=OUT/'replay-bundle.tar.gz'
    need(sha(archive)==manifest['archive_sha256'],'bundle changed')
    before=resource.getrusage(resource.RUSAGE_CHILDREN); start=time.process_time()
    with tempfile.TemporaryDirectory(prefix='certificate-integration-empty-') as temporary:
        root=Path(temporary); seen=set()
        with tarfile.open(archive,'r:gz') as bundle:
            for member in bundle:
                name=member.name
                need(member.isfile() and name not in seen and not Path(name).is_absolute() and '..' not in Path(name).parts,'unsafe member')
                data=bundle.extractfile(member).read(); need(digest(data)==manifest['files_sha256'].get(name),'member changed')
                p=root/name; p.parent.mkdir(parents=True,exist_ok=True); p.write_bytes(data); seen.add(name)
        need(seen==set(manifest['files_sha256']),'missing member')
        code="import sys;sys.path[:0]=['elliptic-curves/cas','elliptic-curves'];import ecsearch.q12o5867_specialization;import retain_cancellation_basis_certificate as r;print(r.inspect())"
        run=subprocess.run(['sage','-python','-c',code],cwd=root,text=True,capture_output=True,timeout=60)
        (RAW/'portable-check.log').write_text(run.stdout+run.stderr)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result={'status':'PASS_EMPTY_ROOT_INTEGRATION_BINDINGS' if run.returncode==0 else 'RETAINED_PORTABLE_FAILURE',
        'returncode':run.returncode,'files':len(seen),'archive_sha256':manifest['archive_sha256'],
        'manifest_sha256':sha(OUT/'replay-manifest.json'),'log_sha256':sha(RAW/'portable-check.log'),
        'full_child_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        'packaging_check_component_cpu_seconds':time.process_time()-start,'point_search_calls':0,
        'boundary':'Empty-root archive bytes, indirect runtime imports, original transcript/source guards and all integration output bindings. Exact arithmetic is evidenced by the original metered integrations, not repeated here.'}
    new_write(OUT/'portable-check.json',result); need(run.returncode==0,'portable check failed; preserve log')
    new_write(OUT/'completion.json',{'status':'RETAINED_COMPLETE_CERTIFICATE_ENGINE_INTEGRATION',
        'protocol_sha256':sha(OUT/'protocol.json'),'archive_sha256':manifest['archive_sha256'],
        'result_sha256':{n:sha(OUT/n) for n in ('summary.json','supervision.json','negative-checks.json','retention-audit.json','dependency-audit.json','portable-check.json','replay-manifest.json')},
        'point_search_calls':0,'boundary':read(OUT/'summary.json')['boundary']})
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('command',choices=['pack','portable']); args=parser.parse_args()
    globals()[args.command]()
