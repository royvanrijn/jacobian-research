"""Construct necessary linear-a cofactor relations; no theorem decision."""
from pathlib import Path
import hashlib,json,time,resource
import sympy as S
from flint import fmpz_poly
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
HERE=Path(__file__).resolve().parent;START=time.monotonic()
source=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text());a,s=S.symbols('a s')
def univariate_coefficients(text):
 p=S.Poly(S.sympify(text),a,s);out=[[] for _ in range(p.degree(a)+1)]
 for i in range(len(out)):
  coeff=S.Poly(p.as_expr(),a).nth(i);q=S.Poly(coeff,s);out[i]=fmpz_poly([int(q.nth(j)) for j in range(q.degree()+1)]) if not q.is_zero else fmpz_poly([])
 return out
f=univariate_coefficients(source['equations'][0]['polynomial']);g=univariate_coefficients(source['equations'][3]['polynomial']);m,n=len(f)-1,len(g)-1;N=m+n
matrix=[[fmpz_poly([]) for _ in range(N)] for _ in range(N)]
for i in range(n):
 for j,c in enumerate(reversed(f)):matrix[i][i+j]=c
for i in range(m):
 for j,c in enumerate(reversed(g)):matrix[n+i][i+j]=c
def bareiss(A):
 A=[row[:] for row in A];N=len(A);previous=fmpz_poly([1]);sign=1
 for k in range(N-1):
  pivot=next((i for i in range(k,N) if A[i][k]),None)
  if pivot is None:return fmpz_poly([])
  if pivot!=k:A[k],A[pivot]=A[pivot],A[k];sign=-sign
  value=A[k][k]
  for i in range(k+1,N):
   left=A[i][k]
   for j in range(k+1,N):
    numerator=A[i][j]*value-left*A[k][j]
    q,r=divmod(numerator,previous);assert not r
    A[i][j]=q
   A[i][k]=fmpz_poly([])
  previous=value
 return sign*A[-1][-1]
cert=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-certificate.json').read_text())
p=fmpz_poly(cert['resultants'][1]['primitive_reduced_coefficients']);squarefree=p//p.gcd(p.derivative());bad=squarefree
records=[]
for row in [0,n,N-1]:
 co=[];start=time.monotonic()
 for col in [N-2,N-1]:
  minor=[[c for j,c in enumerate(line) if j!=col] for i,line in enumerate(matrix) if i!=row]
  value=bareiss(minor)*((-1)**(row+col));co.append(value)
 common=co[0].gcd(co[1]);co=[q//common for q in co]
 # Removing a common factor may lose its zero locus; retain it explicitly.
 const=co[1];bad=bad.gcd(const)
 record={'row':row,'a_cofactor':list(map(int,co[0].coeffs())),'constant_cofactor':list(map(int,const.coeffs())),'removed_common_factor':list(map(int,common.coeffs())),'original_degrees':[q.degree()+common.degree() for q in co],'normalized_degrees':[q.degree() for q in co],'bad_degree_on_squarefree_resultant':squarefree.gcd(const).degree(),'common_bad_degree':bad.degree()}
 records.append(record)
 (HERE/'cofactor-candidate.json').write_text(json.dumps({'source':'../gao-f7-linear-cube/high-root-bivariate-input.json','pair':[0,3],'sylvester_size':N,'resultant_primitive_degree':p.degree(),'resultant_squarefree_degree':squarefree.degree(),'records':records,'scope':'Candidate cofactor relations only. Removed common-factor zeros and all denominator zeros still need coverage. No exclusion.'},indent=2)+'\n')
 print(json.dumps({'row':row,'normalized_degrees':record['normalized_degrees'],'removed_factor_degree':common.degree(),'bad_degree':record['bad_degree_on_squarefree_resultant'],'seconds':time.monotonic()-start}),flush=True)
 if bad.degree()==0:break
print(json.dumps({'seconds':time.monotonic()-START,'candidate_bytes':(HERE/'cofactor-candidate.json').stat().st_size}),flush=True)
