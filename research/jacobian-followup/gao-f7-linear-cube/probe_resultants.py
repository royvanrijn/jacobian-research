"""Bounded construction of univariate resultant obstructions, not an acceptance checker."""
from pathlib import Path
import hashlib,json,sys,time,resource
import sympy as S
from flint import fmpz_mpoly_ctx,fmpz_poly
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
HERE=Path(__file__).resolve().parent;START=time.monotonic()
source=json.loads((HERE/'high-root-bivariate-input.json').read_text())
a,s=S.symbols('a s');ctx=fmpz_mpoly_ctx.get(['a','s']);polys=[]
for row in source['equations']:
 polys.append(ctx.from_dict({e:int(c) for e,c in S.Poly(S.sympify(row['polynomial']),a,s).terms()}))
records=[];common=None
for i,j in [(0,1),(0,3),(1,3),(0,2),(1,2),(2,3)]:
 start=time.monotonic();result=polys[i].resultant(polys[j],'a');terms=result.to_dict()
 assert all(e[0]==0 for e in terms)
 if not terms:
  print(json.dumps({'pair':[i,j],'zero_resultant':True}),flush=True);continue
 coefficients=[0]*(1+max(e[1] for e in terms))
 for e,c in terms.items():coefficients[e[1]]=int(c)
 q=fmpz_poly(coefficients);original_degree=q.degree();q=q//q.content()
 exponents=[]
 for divisor in [fmpz_poly([0,1]),fmpz_poly([-2,0,0,3])]:
  count=0
  while True:
   quotient,remainder=divmod(q,divisor)
   if remainder:break
   q=quotient;count+=1
  exponents.append(count)
 q=q//q.content()
 common=q if common is None else common.gcd(q)
 record={'pair':[i,j],'resultant_coefficients':coefficients,'removed_factors':{'s':exponents[0],'3s^3-2':exponents[1]},'primitive_reduced_coefficients':list(map(int,q.coeffs())),'original_degree':original_degree,'reduced_degree':q.degree(),'seconds':time.monotonic()-start}
 records.append(record)
 (HERE/'resultants-candidate.json').write_text(json.dumps({'source':'high-root-bivariate-input.json','variables':['s'],'records':records,'gcd_degree':common.degree(),'status':'CANDIDATE_ONLY; separate replay required','seconds':time.monotonic()-START},indent=2)+'\n')
 print(json.dumps({'pair':[i,j],'resultant_degree':original_degree,'reduced_degree':q.degree(),'gcd_degree':common.degree(),'seconds':time.monotonic()-start}),flush=True)
 if common.degree()==0:break
