"""Construct the D(s)=0 quadruple-root finite witness; not an acceptance checker."""
from pathlib import Path
import json,time,math
import sympy as S
from flint import fmpz_poly,nmod_poly,nmod_mat,nmod_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic();data=json.loads((HERE/'exceptional-source.json').read_text());s,a,b,z,t=S.symbols('s a b z t')
carrier=S.Poly(S.sympify(data['carrier']),s);H=S.Poly(S.expand(S.sympify(data['linear_resultant']).subs(a,z-s)),s,z)
def s_coefficients(poly):
 q=S.Poly(poly.as_expr(),s);out=[]
 for i in range(q.degree()+1):
  r=S.Poly(q.nth(i),z);out.append(fmpz_poly([int(r.nth(j)) for j in range(r.degree()+1)]) if not r.is_zero else fmpz_poly([]))
 return out
f=s_coefficients(carrier);g=s_coefficients(H);m,n=len(f)-1,len(g)-1;N=m+n;M=[[fmpz_poly([]) for _ in range(N)] for _ in range(N)]
for i in range(n):
 for j,c in enumerate(reversed(f)):M[i][i+j]=c
for i in range(m):
 for j,c in enumerate(reversed(g)):M[n+i][i+j]=c
def determinant(matrix):
 A=[row[:] for row in matrix];N=len(A);previous=fmpz_poly([1]);sign=1
 for k in range(N-1):
  pivot=next((i for i in range(k,N) if A[i][k]),None)
  if pivot is None:return fmpz_poly([])
  if pivot!=k:A[k],A[pivot]=A[pivot],A[k];sign=-sign
  value=A[k][k]
  for i in range(k+1,N):
   left=A[i][k]
   for j in range(k+1,N):
    q,r=divmod(A[i][j]*value-left*A[k][j],previous);assert not r;A[i][j]=q
   A[i][k]=fmpz_poly([])
  previous=value
 return sign*A[-1][-1]
U=determinant(M);P=fmpz_poly(data['projection_coefficients']);assert U.degree()==P.degree()==24
row=n;co=[]
for col in [N-2,N-1]:
 minor=[[c for j,c in enumerate(line) if j!=col] for i,line in enumerate(M) if i!=row]
 co.append(determinant(minor)*((-1)**(row+col)))
num,den=co;prime=1000000007
np=lambda q:nmod_poly([int(c)%prime for c in q],prime)
h=np(P.coeffs());X=np([0,1])
def norm(modulus,value):
 d=modulus.degree();matrix=nmod_mat(d,d,prime);current=value%modulus
 for j in range(d):
  for i in range(d):matrix[i,j]=int(current[i])
  current=current*X%modulus
 return int(matrix.det())
def inverse(value):
 g,u,v=value.xgcd(h);assert g.degree()==0 and int(g[0])!=0
 return u*pow(int(g[0]),-1,prime)%h
norms=[]
for name,value in [('s_denominator',np(den.coeffs())),('projection_squarefree',h.derivative())]:
 det=norm(h,value);assert det;norms.append({'name':name,'dimension':24,'det_mod':det})
S0=np(num.coeffs())*inverse(np(den.coeffs()))%h;A0=(X-S0)%h
F2=S.Poly(S.sympify(data['equations'][0]),b);assert F2.degree()==1

def evaluate(expr,variables,values):
 terms={ex:int(c)%prime for ex,c in S.Poly(expr,*variables).terms()}
 def recurse(data,level):
  if not data:return np([])
  answer=np([])
  for k in range(max(ex[level] for ex in data),-1,-1):
   sub={ex[:level]:c for ex,c in data.items() if ex[level]==k}
   q=np([sum(sub.values())]) if level==0 else recurse(sub,level-1)
   answer=(answer*values[level]+q)%h
  return answer
 return recurse(terms,len(values)-1)
Blead=evaluate(F2.nth(1),[s,a],[S0,A0]);Bconst=evaluate(F2.nth(0),[s,a],[S0,A0]);det=norm(h,Blead);assert det;norms.append({'name':'b_denominator','dimension':24,'det_mod':det});B0=-Bconst*inverse(Blead)%h
source=json.loads((HERE/'quadruple-input.json').read_text());q=[evaluate(S.sympify(row['polynomial']),[s,a,b],[S0,A0,B0]) for row in source['quotient_coefficients']];assert q[-1]==1
ctx=nmod_mpoly_ctx.get(['t','z'],prime);F=ctx.from_dict({(j,i):int(c) for j,row in enumerate(q) for i,c in enumerate(row.coeffs()) if c});result=F.resultant(F.derivative('t'),'t');terms=result.to_dict();assert all(ex[0]==0 for ex in terms)
coeff=[0]*(max(ex[1] for ex in terms)+1)
for (_,j),c in terms.items():coeff[j]=int(c)
det=norm(h,np(coeff));assert det;norms.append({'name':'residual_squarefree','dimension':24,'det_mod':det})
marked=np([])
for row in reversed(q):marked=(marked*S0+row)%h
P23,rem=divmod(P,fmpz_poly([0,1]));assert not rem and P23.degree()==23 and P23[0]!=0
h23=np(P23.coeffs());det=norm(h23,marked);assert det;norms.append({'name':'marked_root_nonorigin','dimension':23,'det_mod':det})
# Exact positive source identities after clearing the s-recovery denominator.
forward=[]
for label,polynomial in [('carrier',S.Poly(carrier.as_expr(),s,z)),('linear_resultant',H)]:
 degree=polynomial.degree(s);Npow=[fmpz_poly([1])];Dpow=[fmpz_poly([1])]
 for _ in range(degree):Npow.append(Npow[-1]*num);Dpow.append(Dpow[-1]*den)
 E=fmpz_poly([])
 for j in range(degree+1):
  r=S.Poly(S.Poly(polynomial.as_expr(),s).nth(j),z);r0=fmpz_poly([int(r.nth(k)) for k in range(r.degree()+1)]) if not r.is_zero else fmpz_poly([])
  E+=r0*Npow[j]*Dpow[degree-j]
 quotient,remainder=divmod(E,P);assert not remainder
 forward.append({'name':label,'homogeneous_coefficients':list(map(int,E.coeffs())),'quotient_coefficients':list(map(int,quotient.coeffs()))})
certificate={'status':'CANDIDATE_ONLY; separate exact replay required','prime':prime,'sylvester_size':N,'cofactor_row':n,'resultant_coefficients':list(map(int,U.coeffs())),'projection_coefficients':list(map(int,P.coeffs())),'nonorigin_projection':list(map(int,P23.coeffs())),'s_numerator':list(map(int,num.coeffs())),'s_denominator':list(map(int,den.coeffs())),'forward_identities':forward,'norms':norms,'quotient_coefficients_mod':[list(map(int,r.coeffs())) for r in q],'quotient_resultant_mod':coeff}
(HERE/'exceptional-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
print(json.dumps({'status':certificate['status'],'projection_degree':P.degree(),'sylvester_size':N,'cofactor_degrees':[num.degree(),den.degree()],'norms':norms,'finite_resultant_degree':len(coeff)-1,'bytes':(HERE/'exceptional-certificate.json').stat().st_size,'seconds':time.monotonic()-START}),flush=True)
