"""Bounded single Singular lift experiment, not a proof decision."""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent
name=sys.argv[1]
assert name.replace('-','').replace('_','').isalnum()
source=HERE/(name+'.sing');log=HERE/(name+'.log')
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
 resource.setrlimit(resource.RLIMIT_FSIZE,(5*1024**2,5*1024**2))
t=time.monotonic()
with log.open('wb') as out:
 try:
  proc=subprocess.run(['Singular','-q',str(source)],stdout=out,stderr=subprocess.STDOUT,timeout=30,preexec_fn=limits)
  outcome='completed';returncode=proc.returncode
 except subprocess.TimeoutExpired:outcome='timeout';returncode=None
receipt={'name':name,'outcome':outcome,'returncode':returncode,'seconds':time.monotonic()-t,'cap_seconds':30,'memory_cap_MB':1536,'output_cap_MB':5,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'log_bytes':log.stat().st_size,'boundary':'A modular lift is support information only; no characteristic-zero inference without literal rational identity replay.'}
(HERE/(name+'-run.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
print(log.read_text()[:500])
