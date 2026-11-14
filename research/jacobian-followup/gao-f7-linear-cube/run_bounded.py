"""Bounded algebra triage for the new four-variable seven-two reduction."""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent
name=sys.argv[1]
assert name.replace('-','').isalnum()
cap=20
source=HERE/(name+'.sing');log=HERE/(name+'.log')
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
 resource.setrlimit(resource.RLIMIT_FSIZE,(10*1024**2,10*1024**2))
start=time.monotonic()
with log.open('wb') as out:
 try:
  proc=subprocess.run(['/usr/bin/Singular','-q',str(source)],stdout=out,stderr=subprocess.STDOUT,timeout=cap,preexec_fn=limits)
  outcome='completed';code=proc.returncode
 except subprocess.TimeoutExpired:outcome='timeout';code=None
receipt={'outcome':outcome,'returncode':code,'seconds':time.monotonic()-start,'cap_seconds':cap,'memory_cap_MiB':1536,'output_cap_MiB':10,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),'log_bytes':log.stat().st_size,'theorem_exclusion':False,'boundary':'Finite-field unit output is triage only. Even characteristic-zero output requires an exported identity and separate literal replay.'}
(HERE/(name+'-run.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt));print(log.read_text()[:3000])
