"""Construct a good specialized degree8 fibre; no global exclusion inference."""
from pathlib import Path
import json,time
import sympy as S
from flint import fmpz_poly,nmod_poly,nmod_mat,nmod_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic();prime=1000000007
old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text());source=json.loads((HERE/'quadruple-input.json').read_text());s,a,b=S.symbols('s a b');curve=S.Poly(S.sympify(old['equations'][3]['polynomial']),a);L=S.Poly(S.sympify(old['linear_equations'][0]),b)
np=lambda coeff:nmod_poly([int(c)%prime for c in coeff],prime)
co=[]
for j in range(9):
 q=S.Poly(curve.nth(j),s);co.append(fmpz_poly([int(q.nth(k)) for k in range(q.degree()+1)]) if not q.is_zero else fmpz_poly([]))
first=co[0];second=sum(co[1:],fmpz_poly([]));f=list(map(int,first.coeffs()));g=list(map(int,second.coeffs()));m,n=len(f)-1,len(g)-1;N=m+n;matrix=[[0]*N for _ in range(N)]
for i in range(n):
 for j,c in enumerate(reversed(f)):matrix[i][i+j]=c%prime
for i in range(m):
 for j,c in enumerate(reversed(g)):matrix[n+i][i+j]=c%prime
vertical_det=int(nmod_mat(matrix,prime).det());assert vertical_det!=0
point=1;sample_poly=S.Poly(curve.as_expr().subs(s,point),a);hcoeff=[int(sample_poly.nth(j)) for j in range(sample_poly.degree()+1)];assert len(hcoeff)==9;h=np(hcoeff);A=np([0,1]);S0=np([point])
def inverse(q):
 g,u,v=q.xgcd(h);assert g.degree()==0 and int(g[0])!=0
 return u*pow(int(g[0]),-1,prime)%h
def norm(q):
 matrix=nmod_mat(8,8,prime);current=q%h
 for j in range(8):
  for i in range(8):matrix[i,j]=int(current[i])
  current=current*A%h
 return int(matrix.det())
def evaluate(expression,variables,values):
 terms={ex:int(c)%prime for ex,c in S.Poly(expression,*variables).terms()}
 def recur(data,level):
  if not data:return np([])
  answer=np([])
  for j in range(max(ex[level] for ex in data),-1,-1):
   sub={ex[:level]:c for ex,c in data.items() if ex[level]==j};q=np([sum(sub.values())]) if level==0 else recur(sub,level-1)
   answer=(answer*values[level]+q)%h
  return answer
 return recur(terms,len(values)-1)
P=evaluate(L.nth(1),[s,a],[S0,A]);Q=evaluate(L.nth(0),[s,a],[S0,A]);B=-Q*inverse(P)%h
quotient=[evaluate(S.sympify(r['polynomial']),[s,a,b],[S0,A,B]) for r in source['quotient_coefficients']];assert quotient[-1]==1
ctx=nmod_mpoly_ctx.get(['t','a'],prime);F=ctx.from_dict({(j,i):int(c) for j,q in enumerate(quotient) for i,c in enumerate(q.coeffs()) if c});D=F.resultant(F.derivative('t'),'t');terms=D.to_dict();delta=[0]*(max(ex[1] for ex in terms)+1);assert all(ex[0]==0 for ex in terms)
for (_,j),c in terms.items():delta[j]=int(c)
marked=np([])
for q in reversed(quotient):marked=(marked*S0+q)%h
norms=[]
for name,value in [('sample_curve_squarefree',h.derivative()),('sample_b_denominator',P),('sample_residual_squarefree',np(delta)),('sample_marked_root',marked)]:
 det=norm(value);assert det!=0,name;norms.append({'name':name,'dimension':8,'det_mod':det})
certificate={'status':'CANDIDATE_ONLY; finite good fibre proves generic units only after source transport and degree preservation','prime':prime,'sample_s':point,'sample_curve_coefficients':hcoeff,'coefficient_first':list(map(int,first.coeffs())),'coefficient_second':list(map(int,second.coeffs())),'coefficient_combination':'sum of all nonconstant-a coefficients','coefficient_sylvester_shape':[N,N],'coefficient_det_mod':vertical_det,'norms':norms,'quotient_coefficients_mod':[list(map(int,q.coeffs())) for q in quotient],'quotient_resultant_mod':delta}
(HERE/'generic-certificate.json').write_text(json.dumps(certificate,indent=2)+'\n')
print(json.dumps({'status':certificate['status'],'coefficient_shape':[N,N],'coefficient_det_mod':vertical_det,'norms':norms,'sample_resultant_degree':len(delta)-1,'seconds':time.monotonic()-START}),flush=True)
