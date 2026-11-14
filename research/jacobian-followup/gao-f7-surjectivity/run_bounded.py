"""Run one retained exact coefficient system with a hard 30-second cap."""
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
name=sys.argv[1]
assert name in ['quartic_cube','quadratic_cube_cubic_square','centered-quadratic','quartic-lift','quartic-readable','centered-lift']
source=HERE/(name+'-input.sing')
log=HERE/(name+'-groebner.log')


def limits():
    resource.setrlimit(resource.RLIMIT_AS,(1536*1024**2,1536*1024**2))
    resource.setrlimit(resource.RLIMIT_FSIZE,(1024**2,1024**2))


start=time.monotonic()
with log.open('wb') as out:
    try:
        proc=subprocess.run(['Singular','-q',str(source)],stdout=out,stderr=subprocess.STDOUT,
                            timeout=30,preexec_fn=limits)
        outcome='completed';returncode=proc.returncode
    except subprocess.TimeoutExpired:
        outcome='timeout';returncode=None
receipt={'case':name,'outcome':outcome,'returncode':returncode,
         'seconds':time.monotonic()-start,'cap_seconds':30,'memory_cap_MB':1536,
         'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest(),
         'log_bytes':log.stat().st_size,
         'boundary':'No theorem or unit certificate inferred solely from a completed process; inspect and replay an explicit ideal identity.'}
(HERE/(name+'-run.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
if log.stat().st_size<1000:
    print(log.read_text())
else:
    print('Polynomial lift retained in '+log.name)
