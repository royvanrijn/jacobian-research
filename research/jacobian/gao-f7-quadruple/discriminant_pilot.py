"""Bounded exact discriminant pilot on the verified quadruple plane chart.

Exploratory output only; source transport and independent replay remain required.
The new representation uses an octic in d, a-degree 7, rather than a global
four-variable Groebner identity. No automatic restart or prime campaign.
"""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent
START=time.monotonic()
def emit(x):
 print(json.dumps(x,default=int),flush=True)
def digest(p):
 return hashlib.sha256(p.read_bytes()).hexdigest()
def worker():
 import sympy as S
 from flint import fmpz_mpoly_ctx
 d,a,s=S.symbols('d a s')
 src=json.loads((HERE/'generic-triple-input.json').read_text())
 ctx=fmpz_mpoly_ctx.get(['d','a','s'])
 C=ctx.from_dict({ex:int(c) for ex,c in S.Poly(S.sympify(src['equations'][1]['polynomial']),d,a,s).terms()})
 emit({'stage':'input','degrees':list(map(int,C.degrees())),'terms':len(C.to_dict()),'seconds':time.monotonic()-START})
 delta=C.discriminant('d')
 terms=delta.to_dict()
 out={'status':'UNVERIFIED_EXACT_CANDIDATE','source_sha256':digest(HERE/'generic-triple-input.json'),'variables':['d','a','s'],'degrees':list(map(int,delta.degrees())),'terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(terms.items())],'boundary':'Candidate exact discriminant only. Independent source relation and determinant replay required.'}
 (HERE/'residual-discriminant-candidate.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
 emit({'stage':'discriminant_saved','degrees':list(map(int,delta.degrees())),'terms':len(terms),'coefficient_bits':max(abs(int(c)).bit_length() for c in terms.values()),'seconds':time.monotonic()-START})
if '--worker' in sys.argv:
 worker()
else:
 def limits():
  resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
  resource.setrlimit(resource.RLIMIT_FSIZE,(100*1024**2,100*1024**2))
 receipt={'name':'primitive-discriminant-pilot','cap_seconds':120,'address_cap_GiB':4,'per_file_cap_MiB':100,'source_sha256':digest(HERE/'generic-triple-input.json'),'script_sha256':digest(Path(__file__)),'theorem_exclusion':False,'motivation':'Compute degree8 residual discriminant before primitive a-remainder reduction on the degree8 plane curve; a different finite-norm representation from retained global lifts.'}
 (HERE/'discriminant-pilot.json').write_text(json.dumps({**receipt,'outcome':'RUNNING'},indent=2)+'\n')
 with (HERE/'discriminant-pilot.log').open('wb') as fd:
  try:
   proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=120,preexec_fn=limits)
   receipt.update(outcome='completed',returncode=proc.returncode)
  except subprocess.TimeoutExpired:
   receipt.update(outcome='timeout',returncode=None)
 receipt['seconds']=time.monotonic()-START
 artifact=HERE/'residual-discriminant-candidate.json'
 if artifact.exists():receipt.update(artifact_bytes=artifact.stat().st_size,artifact_sha256=digest(artifact))
 (HERE/'discriminant-pilot.json').write_text(json.dumps(receipt,indent=2)+'\n')
 emit(receipt)
 print((HERE/'discriminant-pilot.log').read_text()[-4000:])
