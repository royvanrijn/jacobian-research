"""Bounded integer norm projector; preserve all factors and failed attempts."""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def emit(x):print(json.dumps(x,default=int),flush=True)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def worker():
 import sympy as S
 from flint import fmpz_mpoly_ctx,fmpz_poly
 a,s=S.symbols('a s');ctx=fmpz_mpoly_ctx.get(['a','s'])
 data=json.loads((HERE/'primitive-remainder-candidate.json').read_text());R=ctx.from_dict({tuple(ex):int(c) for ex,c in data['terms']})
 src=json.loads((HERE/'generic-triple-input.json').read_text());G=ctx.from_dict({ex:int(c) for ex,c in S.Poly(S.sympify(src['equations'][0]['polynomial']),a,s).terms()})
 emit({'stage':'input','G_degrees':list(G.degrees()),'R_degrees':list(R.degrees()),'s_degree_upper_bound':8*int(R.degrees()[1])+7*27,'seconds':time.monotonic()-START})
 N=G.resultant(R,'a');td=N.to_dict();assert all(ex[0]==0 for ex in td)
 coeff=[0]*(int(N.degrees()[1])+1)
 for (i,j),c in td.items():coeff[int(j)]=int(c)
 out={'status':'UNVERIFIED_EXACT_CANDIDATE','input_sha256':digest(HERE/'primitive-remainder-candidate.json'),'curve_sha256':digest(HERE/'generic-triple-input.json'),'variable':'s','degree':len(coeff)-1,'coefficients':list(map(str,coeff)),'boundary':'Unverified resultant of primitive discriminant remainder. Includes presentation poles with high multiplicity; no exclusions follow without source transport, degree verification and finite-free tests.'}
 (HERE/'norm-projector-candidate.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
 emit({'stage':'saved','degree':len(coeff)-1,'max_coefficient_bits':max(abs(c).bit_length() for c in coeff),'seconds':time.monotonic()-START})
 # Known pole factors are separated only after the entire exact polynomial is saved.
 P=fmpz_poly(coeff);small=[fmpz_poly([0,1]),fmpz_poly([-2,0,0,3]),fmpz_poly([-1,0,0,3]),fmpz_poly([1,0,0,3])];factors=[]
 for f in small:
  power=0
  while P:
   Q,rem=divmod(P,f)
   if rem:break
   P=Q;power+=1
  factors.append({'coefficients':list(map(int,f.coeffs())),'power':power})
 content=int(P.content());P=P//content
 rest={'status':'UNVERIFIED_EXACT_CANDIDATE','norm_sha256':digest(HERE/'norm-projector-candidate.json'),'pole_factors':factors,'content':str(content),'residual_degree':int(P.degree()),'residual_coefficients':list(map(str,P.coeffs())),'squarefree_gcd_degree':int(P.gcd(P.derivative()).degree())}
 (HERE/'norm-projector-factors-candidate.json').write_text(json.dumps(rest,separators=(',',':'))+'\n')
 emit({'stage':'separated','pole_factors':factors,'residual_degree':P.degree(),'gcd_degree':rest['squarefree_gcd_degree'],'seconds':time.monotonic()-START})
if '--worker' in sys.argv:worker()
else:
 def limits():
  resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(100*1024**2,100*1024**2))
 receipt={'name':'exact-norm-projector-pilot','cap_seconds':120,'address_cap_GiB':4,'per_file_cap_MiB':100,'script_sha256':digest(Path(__file__)),'input_sha256':digest(HERE/'primitive-remainder-candidate.json'),'theorem_exclusion':False}
 (HERE/'norm-projector-pilot.json').write_text(json.dumps({**receipt,'outcome':'RUNNING'},indent=2)+'\n')
 with (HERE/'norm-projector-pilot.log').open('wb') as fd:
  try:
   proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=120,preexec_fn=limits);receipt.update(outcome='completed',returncode=proc.returncode)
  except subprocess.TimeoutExpired:receipt.update(outcome='timeout',returncode=None)
 receipt['seconds']=time.monotonic()-START
 for name in ['norm-projector-candidate.json','norm-projector-factors-candidate.json']:
  artifact=HERE/name
  if artifact.exists():receipt[name]={'bytes':artifact.stat().st_size,'sha256':digest(artifact)}
 (HERE/'norm-projector-pilot.json').write_text(json.dumps(receipt,indent=2)+'\n');emit(receipt);print((HERE/'norm-projector-pilot.log').read_text()[-3000:])
