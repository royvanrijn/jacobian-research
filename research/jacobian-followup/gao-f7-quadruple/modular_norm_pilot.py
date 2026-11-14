"""One modular norm pilot to measure residual degree and repeated factors.
No characteristic-zero exclusion is inferred from modular ideals or gcds.
"""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent;START=time.monotonic();PRIME=1000000007
def emit(x):print(json.dumps(x,default=int),flush=True)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def worker():
 import sympy as S
 from flint import nmod_mpoly_ctx,nmod_poly
 a,s=S.symbols('a s');ctx=nmod_mpoly_ctx.get(['a','s'],PRIME)
 data=json.loads((HERE/'primitive-remainder-candidate.json').read_text());R=ctx.from_dict({tuple(ex):int(c)%PRIME for ex,c in data['terms']})
 src=json.loads((HERE/'generic-triple-input.json').read_text());G=ctx.from_dict({ex:int(c)%PRIME for ex,c in S.Poly(S.sympify(src['equations'][0]['polynomial']),a,s).terms()})
 emit({'stage':'input','seconds':time.monotonic()-START})
 N=G.resultant(R,'a');td=N.to_dict();assert all(ex[0]==0 for ex in td);co=[0]*(int(N.degrees()[1])+1)
 for (i,j),c in td.items():co[int(j)]=int(c)
 P=nmod_poly(co,PRIME);factors=[]
 for coeff in [[0,1],[-2,0,0,3],[-1,0,0,3],[1,0,0,3]]:
  f=nmod_poly(coeff,PRIME);power=0
  while P:
   Q,rem=divmod(P,f)
   if rem:break
   P=Q;power+=1
  factors.append({'coefficients':coeff,'power':power})
 gcd=P.gcd(P.derivative())
 record={'status':'MODULAR_TRIAGE_ONLY','prime':PRIME,'input_sha256':digest(HERE/'primitive-remainder-candidate.json'),'norm_degree':len(co)-1,'norm_coefficients':co,'pole_factors':factors,'residual_degree':int(P.degree()),'residual_coefficients':list(map(int,P.coeffs())),'gcd_degree':int(gcd.degree()),'gcd_coefficients':list(map(int,gcd.coeffs())),'boundary':'No characteristic-zero assertion. Exact degree preservation and all removed factors must be established before finite-free reasoning.'}
 (HERE/'modular-norm-candidate.json').write_text(json.dumps(record,separators=(',',':'))+'\n');emit({k:v for k,v in record.items() if not k.endswith('coefficients')})
if '--worker' in sys.argv:worker()
else:
 def limits():
  resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(10*1024**2,10*1024**2))
 receipt={'name':'modular-norm-pilot','prime':PRIME,'cap_seconds':30,'address_cap_GiB':2,'per_file_cap_MiB':10,'script_sha256':digest(Path(__file__)),'theorem_exclusion':False}
 with (HERE/'modular-norm-pilot.log').open('wb') as fd:
  try:
   proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=30,preexec_fn=limits);receipt.update(outcome='completed',returncode=proc.returncode)
  except subprocess.TimeoutExpired:receipt.update(outcome='timeout',returncode=None)
 receipt['seconds']=time.monotonic()-START
 artifact=HERE/'modular-norm-candidate.json'
 if artifact.exists():receipt.update(artifact_bytes=artifact.stat().st_size,artifact_sha256=digest(artifact))
 (HERE/'modular-norm-pilot.json').write_text(json.dumps(receipt,indent=2)+'\n');emit(receipt);print((HERE/'modular-norm-pilot.log').read_text()[-3000:])
