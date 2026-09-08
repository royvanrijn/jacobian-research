import os,json,subprocess,time
from pathlib import Path
ROOT=Path('/home/royvanrijn/src/jacobian-research/research');D=ROOT/'artifacts/local/elliptic-curves/det1092-basis-cascade-v1';PKG=ROOT/'artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1'
roster=json.loads((PKG/'roster.json').read_text());waiting={r['arm_id'] for r in roster['arms']};jobs=[]
while waiting or jobs:
 for aid in list(waiting):
  path=D/aid/'supervisor-result.json'
  if not path.exists():continue
  result=json.loads(path.read_text());waiting.remove(aid)
  if result['exit_code']!=0:
   print('PRIMARY FAILED; REPLAY NEEDS PREFIX REVIEW',aid,flush=True);continue
  log=(PKG/f'replay-{aid}.log').open('x');log.write('Starting commit: '+roster['starting_commit']+'\n');log.flush()
  env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
  command=['/home/royvanrijn/.local/bin/sage','-python',str(ROOT/'elliptic-curves/cas/det1092_basis_cascade_replay.sage'),aid]
  proc=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT);jobs.append((aid,proc,log,time.monotonic()));print('REPLAY LAUNCHED',aid,flush=True)
 for job in list(jobs):
  aid,proc,log,start=job;code=proc.poll()
  if code is None:continue
  log.close();(PKG/f'replay-process-{aid}.json').write_text(json.dumps({'starting_commit':roster['starting_commit'],'arm_id':aid,'exit_code':code,'wall_seconds':time.monotonic()-start},indent=2)+'\n');jobs.remove(job);print('REPLAY EXIT',aid,code,flush=True)
 if waiting or jobs:time.sleep(2)
