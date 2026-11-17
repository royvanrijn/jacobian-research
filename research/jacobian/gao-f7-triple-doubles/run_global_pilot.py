"""One capped modular pilot for each new five-variable coefficient branch."""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent;branch=sys.argv[1];assert branch in ['delta-zero','delta-nonzero'];START=time.monotonic();script=HERE/(branch+'-p32003.sing')
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(20*1024**2,20*1024**2))
with (HERE/(branch+'-p32003.log')).open('wb') as fd:
 try:
  proc=subprocess.run(['/usr/bin/Singular','-q',str(script)],stdout=fd,stderr=subprocess.STDOUT,timeout=20,preexec_fn=limits);outcome='completed';code=proc.returncode
 except subprocess.TimeoutExpired:outcome='timeout';code=None
receipt={'branch':branch,'outcome':outcome,'returncode':code,'seconds':time.monotonic()-START,'cap_seconds':20,'address_cap_GiB':2,'per_file_cap_MiB':20,'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((HERE/(branch+'-fivevar-input.json')).read_bytes()).hexdigest(),'theorem_exclusion':False,'boundary':'Finite-field triage only. Both source branches retained. No rational identity or global theorem follows without literal characteristic-zero replay. No automatic prime or cap escalation follows this pilot.'}
(HERE/(branch+'-p32003-run.json')).write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt));print((HERE/(branch+'-p32003.log')).read_text()[-1500:])
