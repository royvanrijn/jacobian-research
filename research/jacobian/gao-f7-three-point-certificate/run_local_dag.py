from pathlib import Path
import json,os,resource,subprocess,sys,time,hashlib
HERE=Path(__file__).resolve().parent
primes=sys.argv[1:];assert primes and all(x.isdigit() for x in primes)
name='dag-'+primes[0]+'-'+primes[-1];log=HERE/(name+'.log')
def caps():
 resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
 resource.setrlimit(resource.RLIMIT_FSIZE,(100*1024**2,100*1024**2))
cmd=['julia','--startup-file=no','--heap-size-hint=1G','--project='+str(HERE/'julia'),str(HERE/'generate_local_dag.jl')]+primes
t=time.monotonic()
with log.open('wb') as fd:
 try:
  run=subprocess.run(cmd,stdout=fd,stderr=subprocess.STDOUT,timeout=120,preexec_fn=caps,env=dict(os.environ,JULIA_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'));status='completed';code=run.returncode
 except subprocess.TimeoutExpired:status='timeout';code=None
receipt={'command':cmd,'status':status,'code':code,'seconds':time.monotonic()-t,'time_cap':120,'virtual_memory_cap_GB':8,'output_cap_MB':100,'source_sha256':hashlib.sha256((HERE/'generate_local_dag.jl').read_bytes()).hexdigest(),'patch_sha256':hashlib.sha256((HERE/'local_dag_patch.jl').read_bytes()).hexdigest(),'boundary':'Modular local identities are not a characteristic-zero certificate.'}
(HERE/(name+'-run.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt));print(log.read_text()[-2000:])
