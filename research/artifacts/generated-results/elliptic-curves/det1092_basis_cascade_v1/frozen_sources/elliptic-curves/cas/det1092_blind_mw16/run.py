#!/usr/bin/env python3
"""Dispatch the frozen17 once. Does not open the parent or any evaluation input."""
import hashlib,json,os,signal,subprocess,tempfile,time,datetime,resource
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4];CAS=Path(__file__).resolve().parent
PKG=ROOT/'research/artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2,sort_keys=True);f.write('\n')
 p.chmod(0o444)
protocol=json.loads((PKG/'protocol.json').read_text());commit=protocol['starting_commit']
assert sha(PKG/'roster.json')==protocol['roster_sha256']
for name in ['worker.py','isolate.py']:
 assert sha(CAS/name)==protocol['source_bindings'][str((CAS/name).relative_to(ROOT))]
roster=json.loads((PKG/'roster.json').read_text())['arms'];assert len(roster)==17
for arm in roster:assert sha(PKG/arm['fixture'])==arm['fixture_sha256']
write(PKG/'dispatch.json',{'starting_commit':commit,'protocol_sha256':sha(PKG/'protocol.json'),'runner_sha256':sha(Path(__file__)),'launch_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'arms':[r['arm_id'] for r in roster],'sequential':True,'reruns':0})
(PKG/'runs').mkdir();runs=[]
for arm in roster:
 aid=arm['arm_id'];output=PKG/'runs'/aid;output.mkdir();before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic();timeout=False
 with tempfile.TemporaryDirectory(prefix='det1092-blind-') as scratch:
  cmd=['unshare','-Unpf','python3',str(CAS/'isolate.py'),str(CAS/'worker.py'),str(PKG/arm['fixture']),str(output),scratch]
  with (output/'stdout.txt').open('x') as log:
   log.write('Starting commit: '+commit+'\n');log.flush()
   proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
   try:code=proc.wait(timeout=240)
   except subprocess.TimeoutExpired:
    timeout=True;os.killpg(proc.pid,signal.SIGKILL);code=proc.wait()
 after=resource.getrusage(resource.RUSAGE_CHILDREN)
 accounting={'starting_commit':commit,'arm_id':aid,'returncode':code,'wall_cap_reached':timeout,'wall_seconds':time.monotonic()-start,'user_seconds':after.ru_utime-before.ru_utime,'system_seconds':after.ru_stime-before.ru_stime,'cumulative_child_max_rss_KiB':after.ru_maxrss,'command':cmd[:-1]+['<ephemeral-private-scratch>']}
 if not (output/'result.json').exists():
  write(output/'result.json',{'starting_commit':commit,'arm_id':aid,'fixture_sha256':arm['fixture_sha256'],'terminal_status':'CENSORED_RESOURCE_CAP' if timeout or code in [-9,-24,137,152] else 'EXPERIMENT_FAILURE','reason':'Worker exited without final certificate; all partial events and stdout retained. No retry.','resource':accounting})
 result=json.loads((output/'result.json').read_text());accounting['terminal_status']=result['terminal_status'];write(output/'accounting.json',accounting)
 for p in output.iterdir():p.chmod(0o444)
 runs.append({'arm_id':aid,'terminal_status':result['terminal_status'],'files':{str(p.relative_to(PKG)):sha(p) for p in output.iterdir()}})
 print(aid,result['terminal_status'],round(accounting['wall_seconds'],3),flush=True)
write(PKG/'primary_complete.json',{'starting_commit':commit,'protocol_sha256':sha(PKG/'protocol.json'),'roster_sha256':sha(PKG/'roster.json'),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'evaluation_has_not_started':True,'runs':runs})
print('ALL17 PRIMARY ARMS TERMINAL; outputs frozen before evaluation',flush=True)
