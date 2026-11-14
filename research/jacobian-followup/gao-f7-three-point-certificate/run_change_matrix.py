"""Single bounded Groebner.jl change-matrix attempt, preserving every log."""
from pathlib import Path
import hashlib,json,os,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent
field=sys.argv[1];assert field=='QQ' or field.isdigit()
name='cubic-six-'+field
log=HERE/(name+'.log');out=HERE/(name+'-multipliers.txt')
command=['julia','--startup-file=no','--heap-size-hint=1G','--project='+str(HERE/'julia'),str(HERE/'cubic-six-change-matrix.jl'),field,str(out)]
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
 resource.setrlimit(resource.RLIMIT_FSIZE,(20*1024**2,20*1024**2))
env=dict(os.environ,JULIA_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
started=time.monotonic()
with log.open('wb') as fd:
 try:
  run=subprocess.run(command,stdout=fd,stderr=subprocess.STDOUT,timeout=90,env=env,preexec_fn=limits);status='completed';code=run.returncode
 except subprocess.TimeoutExpired:status='timeout';code=None
receipt={'field':field,'status':status,'code':code,'seconds':time.monotonic()-started,'time_cap':90,'virtual_memory_cap_GB':8,'heap_hint_GB':1,'output_cap_MB':20,'command':command,'script_sha256':hashlib.sha256((HERE/'cubic-six-change-matrix.jl').read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'certificate_exported':out.exists(),'boundary':'Only a literal rational identity replay proves a characteristic-zero exclusion; finite-field output and process success alone are insufficient.'}
(HERE/(name+'-run.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt));print(log.read_text()[-1500:])
