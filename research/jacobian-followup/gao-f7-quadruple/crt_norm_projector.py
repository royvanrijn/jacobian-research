"""Bounded exact CRT reconstruction of a 15x15 integer polynomial determinant.

For every coefficient, |N_j| <= ||G||_1^7 ||R||_1^8: expand the
Sylvester determinant, bound by the product of row coefficient l1 norms.
CRT modulus M>2B proves signed reconstruction of ALL coefficients, including
zeros up to the proved degree bound. This generator is untrusted; separate
modular determinant interpolation, bound and source replay are required.
"""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time,math
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def emit(x):print(json.dumps(x,default=int),flush=True)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def worker():
 import sympy as S
 from flint import nmod_mpoly_ctx,fmpz_poly,fmpz
 a,s=S.symbols('a s');data=json.loads((HERE/'primitive-remainder-candidate.json').read_text());rt={tuple(ex):int(c) for ex,c in data['terms']}
 src=json.loads((HERE/'generic-triple-input.json').read_text());gt={ex:int(c) for ex,c in S.Poly(S.sympify(src['equations'][0]['polynomial']),a,s).terms()}
 bound=sum(abs(c) for c in gt.values())**7*sum(abs(c) for c in rt.values())**8;degree_bound=8*1435+7*27;acc=[0]*(degree_bound+1);M=1;primes=[];p=1000000007
 emit({'stage':'bound','coefficient_bound_bits':bound.bit_length(),'degree_bound':degree_bound,'max_primes':512,'seconds':time.monotonic()-START})
 while M<=2*bound:
  assert len(primes)<512
  while not fmpz(p).is_prime():p+=2
  ctx=nmod_mpoly_ctx.get(['a','s'],p);G=ctx.from_dict({ex:c%p for ex,c in gt.items()});R=ctx.from_dict({ex:c%p for ex,c in rt.items()});N=G.resultant(R,'a');residue=[0]*(degree_bound+1)
  for (i,j),c in N.to_dict().items():assert i==0 and j<=degree_bound;residue[int(j)]=int(c)
  inv=pow(M%p,-1,p)
  for j,r in enumerate(residue):acc[j]+=M*((r-acc[j])%p*inv%p)
  M*=p;primes.append(p);p+=2
  if len(primes)%16==0:
   checkpoint={'status':'INCOMPLETE_CRT','input_sha256':digest(HERE/'primitive-remainder-candidate.json'),'degree_bound':degree_bound,'coefficient_bound_bits':bound.bit_length(),'primes':primes,'modulus':str(M),'coefficients_modulus':list(map(str,acc))}
   (HERE/'crt-norm-checkpoint.json').write_text(json.dumps(checkpoint,separators=(',',':'))+'\n');emit({'stage':'checkpoint','primes':len(primes),'modulus_bits':M.bit_length(),'seconds':time.monotonic()-START})
 coeff=[r if 2*r<=M else r-M for r in acc]
 assert all(abs(c)<=bound for c in coeff)
 while coeff and coeff[-1]==0:coeff.pop()
 out={'status':'EXACT_CRT_CANDIDATE_REQUIRES_REPLAY','input_sha256':digest(HERE/'primitive-remainder-candidate.json'),'curve_sha256':digest(HERE/'generic-triple-input.json'),'variable':'s','degree':len(coeff)-1,'degree_bound':degree_bound,'coefficient_bound_formula':'l1(G)^7*l1(R)^8','coefficient_bound_bits':bound.bit_length(),'crt_modulus_bits':M.bit_length(),'primes':primes,'coefficients':list(map(str,coeff)),'boundary':'CRT generator only. Independent coefficient-bound, input, degree and modular determinant-interpolation replay required.'}
 (HERE/'crt-norm-candidate.json').write_text(json.dumps(out,separators=(',',':'))+'\n');emit({'stage':'saved','degree':len(coeff)-1,'primes':len(primes),'coefficient_bits':max(abs(c).bit_length() for c in coeff),'seconds':time.monotonic()-START})
 # Exact division preserves every pole factor in this separate record.
 P=fmpz_poly(coeff);factors=[]
 for co in [[0,1],[-2,0,0,3],[-1,0,0,3],[1,0,0,3]]:
  f=fmpz_poly(co);power=0
  while P:
   Q,rem=divmod(P,f)
   if rem:break
   P=Q;power+=1
  factors.append({'coefficients':co,'power':power})
 pole=json.loads((HERE/'coordinate-pole-resultant-candidate.json').read_text());B=fmpz_poly(list(map(int,next(r['coefficients'] for r in pole['factors'] if len(r['coefficients'])==88))));power=0
 while P:
  Q,rem=divmod(P,B)
  if rem:break
  P=Q;power+=1
 content=int(P.content());H=P//content
 rec={'status':'EXACT_FACTORIZATION_CANDIDATE_REQUIRES_REPLAY','norm_sha256':digest(HERE/'crt-norm-candidate.json'),'pole_factors':factors,'coordinate_pole_factor':{'coefficients':list(map(str,B.coeffs())),'power':power},'content':str(content),'residual_degree':int(H.degree()),'residual_coefficients':list(map(str,H.coeffs()))}
 (HERE/'crt-norm-factors-candidate.json').write_text(json.dumps(rec,separators=(',',':'))+'\n');emit({'stage':'factors_saved','pole_factors':factors,'coordinate_pole_power':power,'residual_degree':H.degree(),'seconds':time.monotonic()-START})
if '--worker' in sys.argv:worker()
else:
 def limits():
  resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(100*1024**2,100*1024**2))
 receipt={'name':'exact-crt-norm-projector','cap_seconds':300,'max_primes':512,'address_cap_GiB':4,'per_file_cap_MiB':100,'script_sha256':digest(Path(__file__)),'theorem_exclusion':False,'justification':'Direct integer resultant timed out; one modular norm revealed a degree912 squarefree residual and degree87 coordinate-pole factor. A proved determinant coefficient height bounds a finite exact CRT reconstruction, rather than an open-ended prime search.'}
 (HERE/'crt-norm-run.json').write_text(json.dumps({**receipt,'outcome':'RUNNING'},indent=2)+'\n')
 with (HERE/'crt-norm.log').open('wb') as fd:
  try:
   proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=300,preexec_fn=limits);receipt.update(outcome='completed',returncode=proc.returncode)
  except subprocess.TimeoutExpired:receipt.update(outcome='timeout',returncode=None)
 receipt['seconds']=time.monotonic()-START
 for name in ['crt-norm-candidate.json','crt-norm-factors-candidate.json','crt-norm-checkpoint.json']:
  artifact=HERE/name
  if artifact.exists():receipt[name]={'bytes':artifact.stat().st_size,'sha256':digest(artifact)}
 (HERE/'crt-norm-run.json').write_text(json.dumps(receipt,indent=2)+'\n');emit(receipt);print((HERE/'crt-norm.log').read_text()[-4000:])
