"""Bounded lift of x^3 on the new four-variable normalized H0 chart."""
from pathlib import Path
import hashlib,json,resource,subprocess,time
HERE=Path(__file__).resolve().parent;START=time.monotonic();script=HERE/'two-scale-pole-rational-lift.sing'
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(20*1024**2,20*1024**2))
with (HERE/'two-scale-pole-rational-lift.log').open('wb') as fd:
 try:
  proc=subprocess.run(['/usr/bin/Singular','-q',str(script)],stdout=fd,stderr=subprocess.STDOUT,timeout=30,preexec_fn=limits);outcome='completed';code=proc.returncode
 except subprocess.TimeoutExpired:outcome='timeout';code=None
rec={'outcome':outcome,'returncode':code,'seconds':time.monotonic()-START,'cap_seconds':30,'address_cap_GiB':2,'per_file_cap_MiB':20,'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((HERE/'two-scale-pole-power-input.json').read_bytes()).hexdigest(),'motivation':'New four-variable finite-scale x^3 target, supported by a0.115second modular membership computation; the previous five-variable rational unit lift is not restarted.','theorem_exclusion':False,'boundary':'Accept only after literal complete characteristic-zero identity and full source/chart transport replay.'}
f=HERE/'two-scale-pole-rational-identity.txt'
if f.exists():rec['artifact']={'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()}
(HERE/'two-scale-pole-rational-lift-run.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec));print((HERE/'two-scale-pole-rational-lift.log').read_text()[-1500:])
