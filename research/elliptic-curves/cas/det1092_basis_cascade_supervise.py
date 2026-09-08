#!/usr/bin/env python3
"""Run each frozen arm once, preserve timeout/failure exit accounting."""
import hashlib,json,os,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'artifacts/local/elliptic-curves/det1092-basis-cascade-v1';PKG=ROOT/'artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1'
roster=json.loads((PKG/'roster.json').read_text());jobs=[]
for row in roster['arms']:
 aid=row['arm_id'];arm=D/aid
 if (arm/'launch.json').exists():raise FileExistsError('arm already launched; no silent rerun')
 env={**os.environ,'CASCADE_ARM':aid,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1','PYTHONDONTWRITEBYTECODE':'1'}
 command=['/usr/bin/timeout','--signal=TERM','--kill-after=15','14400','/home/royvanrijn/.local/bin/sage','-python',str(ROOT/'elliptic-curves/cas/det1092_basis_cascade_v1.sage')]
 start=time.time();record={'starting_commit':roster['starting_commit'],'command':command,'arm_id':aid,'start_unix':start,'supervisor_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (arm/'launch.json').write_text(json.dumps(record,indent=2)+'\n');log=(arm/'run.log').open('x')
 proc=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
 jobs.append((aid,proc,log,start));print('LAUNCHED',aid,proc.pid,flush=True)
while jobs:
 for item in list(jobs):
  aid,proc,log,start=item;code=proc.poll()
  if code is None:continue
  log.close();arm=D/aid
  record={'starting_commit':roster['starting_commit'],'arm_id':aid,'exit_code':code,'wall_seconds':time.time()-start,'terminal_file_exists':(arm/'run/terminal.json').exists(),'classification':'TERMINAL' if code==0 and (arm/'run/terminal.json').exists() else 'CENSORED_OR_FAILED'}
  (arm/'supervisor-result.json').write_text(json.dumps(record,indent=2)+'\n');print('EXIT',record,flush=True);jobs.remove(item)
 if jobs:time.sleep(2)
