"""One bounded characteristic-zero identity lift on the new delta-zero branch."""
from pathlib import Path
import hashlib,json,resource,subprocess,time
HERE=Path(__file__).resolve().parent;START=time.monotonic();script=HERE/'delta-zero-rational-lift.sing'
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(20*1024**2,20*1024**2))
with (HERE/'delta-zero-rational-lift.log').open('wb') as fd:
 try:
  proc=subprocess.run(['/usr/bin/Singular','-q',str(script)],stdout=fd,stderr=subprocess.STDOUT,timeout=30,preexec_fn=limits);outcome='completed';code=proc.returncode
 except subprocess.TimeoutExpired:outcome='timeout';code=None
artifact=HERE/'delta-zero-rational-identity.txt';rec={'outcome':outcome,'returncode':code,'seconds':time.monotonic()-START,'cap_seconds':30,'address_cap_GiB':2,'per_file_cap_MiB':20,'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((HERE/'delta-zero-fivevar-input.json').read_bytes()).hexdigest(),'theorem_exclusion':False,'boundary':'Only a literal complete characteristic-zero identity replay can prove exclusion. This run targets the newly reduced five-variable delta-zero branch; no unchanged six-variable restart or cap increase follows.'}
if artifact.exists():rec['artifact']={'bytes':artifact.stat().st_size,'sha256':hashlib.sha256(artifact.read_bytes()).hexdigest()}
(HERE/'delta-zero-rational-lift-run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec));print((HERE/'delta-zero-rational-lift.log').read_text()[-1500:])
