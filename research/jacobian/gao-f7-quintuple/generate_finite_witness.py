"""Finite free algebra witness construction; not a proof acceptance step."""
from pathlib import Path
import json,time
import sympy as S
from flint import nmod_poly,nmod_mat,nmod_mpoly_ctx
HERE=Path(__file__).resolve().parent;START=time.monotonic();cert=json.loads((HERE/'certificate.json').read_text());prime=cert['prime'];s,a,b=S.symbols('s a b')
old=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text());source=json.loads((HERE/'quintuple-input.json').read_text());co=cert['cofactor']
def poly(coeffs):return nmod_poly([c%prime for c in coeffs],prime)
def inverse(f,h):
 g,u,v=f.xgcd(h);assert g.degree()==0 and int(g[0])!=0
 return u*pow(int(g[0]),-1,prime)%h
def norm(h,f):
 d=h.degree();M=nmod_mat(d,d,prime);value=f%h;x=poly([0,1])
 for j in range(d):
  for i in range(d):M[i,j]=int(value[i])
  value=value*x%h
 return int(M.det())
def evaluate(expression,h,values):
 V=[s,a] if len(values)==2 else [s,a,b];q=S.Poly(expression,*V);powers=[]
 for var,value in zip(V,values):
  rows=[poly([1])]
  for _ in range(q.degree(var)):rows.append(rows[-1]*value%h)
  powers.append(rows)
 result=poly([])
 for ex,c in q.terms():
  term=poly([int(c)%prime])
  for j,e in enumerate(ex):term=term*powers[j][e]%h
  result=(result+term)%h
 return result
norms=[]
for name,field_key,expression in [('cofactor_common','generic_initial',poly(co['removed_common_factor'])),('a_denominator','generic_initial',poly(co['constant_cofactor']))]:
 h=poly(cert[field_key]);det=norm(h,expression);assert det
 norms.append({'name':name,'dimension':h.degree(),'det_mod':det})
h=poly(cert['generic_out']);A=poly(co['a_cofactor'])*inverse(poly(co['constant_cofactor']),h)%h
extra=evaluate(S.sympify(old['equations'][4]['polynomial']),h,[poly([0,1]),A]);det=norm(h,extra);assert det
norms.append({'name':'excluded_projection','dimension':h.degree(),'det_mod':det})
h=poly(cert['generic_retained']);A=poly(co['a_cofactor'])*inverse(poly(co['constant_cofactor']),h)%h;ss=poly([0,1]);L=S.Poly(S.sympify(old['linear_equations'][0]),b)
P=evaluate(L.nth(1),h,[ss,A]);Q=evaluate(L.nth(0),h,[ss,A]);det=norm(h,P);assert det
norms.append({'name':'b_denominator','dimension':h.degree(),'det_mod':det});B=-Q*inverse(P,h)%h
quotient=[evaluate(S.sympify(row['polynomial']),h,[ss,A,B]) for row in source['quotient_coefficients']];assert quotient[7]==1
ctx=nmod_mpoly_ctx.get(['t','s'],prime);F=ctx.from_dict({(j,i):int(c) for j,row in enumerate(quotient) for i,c in enumerate(row.coeffs()) if c});D=F.resultant(F.derivative('t'),'t');terms=D.to_dict();assert all(ex[0]==0 for ex in terms)
coeff=[0]*(1+max(ex[1] for ex in terms))
for (_,i),c in terms.items():coeff[i]=int(c)
det=norm(h,poly(coeff));assert det;norms.append({'name':'residual_squarefree','dimension':h.degree(),'det_mod':det})
marked=nmod_poly([],prime)
for q in reversed(quotient):marked=(marked*ss+q)%h
det=norm(h,marked);assert det;norms.append({'name':'marked_root','dimension':h.degree(),'det_mod':det})
for name,value in [('generic_chart',nmod_poly([0,-2,0,0,3],prime)),('projection_squarefree',h.derivative())]:
 det=norm(h,value);assert det;norms.append({'name':name,'dimension':h.degree(),'det_mod':det})
result={'status':'CANDIDATE_ONLY','prime':prime,'norms':norms,'quotient_coefficients_mod':[list(map(int,row.coeffs())) for row in quotient],'quotient_resultant_mod':coeff,'seconds':time.monotonic()-START}
(HERE/'finite-witness.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'norms':norms,'resultant_degree':len(coeff)-1,'seconds':time.monotonic()-START}),flush=True)
