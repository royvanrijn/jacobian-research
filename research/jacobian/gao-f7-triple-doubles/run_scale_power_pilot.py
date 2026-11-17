"""One bounded pilot for each two-scale carrier-power certificate target."""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent;branch=sys.argv[1];assert branch in ['pole','generic'];START=time.monotonic();script=HERE/('two-scale-'+branch+'-power-p32003.sing')
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(10*1024**2,10*1024**2))
with (HERE/('two-scale-'+branch+'-power-p32003.log')).open('wb') as fd:
 try:
  proc=subprocess.run(['/usr/bin/Singular','-q',str(script)],stdout=fd,stderr=subprocess.STDOUT,timeout=20,preexec_fn=limits);outcome='completed';code=proc.returncode
 except subprocess.TimeoutExpired:outcome='timeout';code=None
r={'outcome':outcome,'returncode':code,'seconds':time.monotonic()-START,'cap_seconds':20,'address_cap_GiB':2,'per_file_cap_MiB':10,'branch':branch,'script_sha256':hashlib.sha256(script.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((HERE/('two-scale-'+branch+'-power-input.json')).read_bytes()).hexdigest(),'theorem_exclusion':False,'boundary':'New four-variable two-scale carrier-power representation. No rational exclusion follows from modular reductions. No automatic prime or cap increase.'}
(HERE/('two-scale-'+branch+'-power-p32003-run.json')).write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));print((HERE/('two-scale-'+branch+'-power-p32003.log')).read_text()[-2000:])
