"""Good-reduction finite-algebra probe for the generic fivefold piece.

Candidate only. Every characteristic-zero presentation and coverage step
requires a separate checker and written proof before registration.
"""
from pathlib import Path
import json,time,resource
import sympy as S
from flint import nmod_poly,nmod_mpoly_ctx,nmod_mat,fmpz_poly
resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
HERE=Path(__file__).resolve().parent;START=time.monotonic();prime=1000000007
rur=json.loads((HERE/'rur-candidate.json').read_text());cofactor=json.loads((HERE/'cofactor-candidate.json').read_text())['records'][0];source=json.loads((HERE/'quintuple-input.json').read_text())
def poly(rows):return nmod_poly([c%prime for c in rows],prime)
def norm(h,f):
 n=h.degree();columns=[];current=f%h;x=nmod_poly([0,1],prime)
 for _ in range(n):columns.append(current);current=current*x%h
 matrix=nmod_mat(n,n,prime)
 for j,c in enumerate(columns):
  for i in range(n):matrix[i,j]=int(c[i])
 return int(matrix.det())
checks={}
good=poly(rur['generic_initial']);retained=poly(rur['generic_retained'])
out_integer=fmpz_poly(rur['generic_initial'])//fmpz_poly(rur['generic_retained']);out=poly(list(map(int,out_integer.coeffs())))
for name,h,f in [('cofactor_common',good,poly(cofactor['removed_common_factor'])),('a_denominator_good',good,poly(rur['a_denominator'])),('excluded_projection',out,poly(rur['extra_equation_substitution'])),('b_denominator_retained',retained,poly(rur['b_denominator']))]:
 det=norm(h,f);checks[name]={'degree':h.degree(),'det_mod':det};print(json.dumps({'stage':name,**checks[name]}),flush=True)
 assert det!=0,name
h=retained
def inverse(f):
 g,u,v=f.xgcd(h);assert g.degree()==0 and int(g[0])!=0
 return u*pow(int(g[0]),-1,prime)%h
ss,aa,bb=S.symbols('s a b');s=nmod_poly([0,1],prime)
a=poly(rur['a_numerator'])*inverse(poly(rur['a_denominator']))%h
b=poly(rur['b_numerator'])*inverse(poly(rur['b_denominator']))%h
powers=[]
for value,cap in [(s,30),(a,12),(b,8)]:
 rows=[nmod_poly([1],prime)]
 for _ in range(cap):rows.append(rows[-1]*value%h)
 powers.append(rows)
q=[]
for row in source['quotient_coefficients']:
 value=nmod_poly([],prime)
 for (i,j,k),c in S.Poly(S.sympify(row['polynomial']),ss,aa,bb).terms():value=(value+int(c)*powers[0][i]*powers[1][j]*powers[2][k])%h
 q.append(value)
assert q[7]==1
ctx=nmod_mpoly_ctx.get(['t','s'],prime);T,Z=ctx.gens()
Q=ctx.from_dict({(j,i):int(c) for j,row in enumerate(q) for i,c in enumerate(row.coeffs()) if c})
D=Q.resultant(Q.derivative('t'),'t');terms=D.to_dict();assert all(e[0]==0 for e in terms)
coefficients=[0]*(1+max(e[1] for e in terms))
for (_,i),c in terms.items():coefficients[i]=int(c)
delta=nmod_poly(coefficients,prime)%h;det=norm(h,delta);checks['residual_squarefree']={'degree':h.degree(),'det_mod':det,'resultant_degree':len(coefficients)-1}
print(json.dumps({'stage':'residual_squarefree',**checks['residual_squarefree']}),flush=True)
data={'status':'CANDIDATE_ONLY; characteristic-zero coverage and replay pending','prime':prime,'norms':checks,'residual_quotient_coefficients_mod':[list(map(int,c.coeffs())) for c in q],'residual_resultant_mod':coefficients,'generic_out':list(map(int,out_integer.coeffs())),'seconds':time.monotonic()-START}
(HERE/'finite-generic-candidate.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'status':data['status'],'seconds':time.monotonic()-START,'bytes':(HERE/'finite-generic-candidate.json').stat().st_size}),flush=True)
