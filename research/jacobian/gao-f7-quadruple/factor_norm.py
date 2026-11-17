"""Exact product decomposition of a saved CRT resultant, no refactorization search."""
from pathlib import Path
import hashlib,json,time
from flint import fmpz_poly,nmod_poly
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
data=json.loads((HERE/'crt-norm-candidate.json').read_text());P=fmpz_poly(list(map(int,data['coefficients'])));factors=[]
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
prime=1000000007;hm=nmod_poly([int(c)%prime for c in H.coeffs()],prime);bm=nmod_poly([int(c)%prime for c in B.coeffs()],prime)
rec={'status':'EXACT_FACTORIZATION_CANDIDATE_REQUIRES_REPLAY','norm_sha256':digest(HERE/'crt-norm-candidate.json'),'pole_factors':factors,'coordinate_pole_factor':{'coefficients':list(map(str,B.coeffs())),'power':power},'content':str(content),'residual_degree':int(H.degree()),'residual_coefficients':list(map(str,H.coeffs())),'prime':prime,'residual_leading_mod':int(hm.leading_coefficient()),'squarefree_gcd_degree':int(hm.gcd(hm.derivative()).degree()),'coordinate_pole_gcd_degree':int(hm.gcd(bm).degree())}
(HERE/'crt-norm-factors-candidate.json').write_text(json.dumps(rec,separators=(',',':'))+'\n')
print(json.dumps({k:v for k,v in rec.items() if k not in ['residual_coefficients','content','coordinate_pole_factor']}));print('seconds',time.monotonic()-START)
