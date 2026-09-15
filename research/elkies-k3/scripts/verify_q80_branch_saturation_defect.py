#!/usr/bin/env python3
"""Exact finite modulo3 quotients and a quadratic-field half identity."""
import hashlib,json,resource
from pathlib import Path
from fractions import Fraction
import sympy as S
from verify_q80_branch_trace_specialization import add
resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
root=Path(__file__).resolve().parents[2]
packet=root/'artifacts/generated-results/elkies-k3-q80-branch-saturation-defect-v1'
inp=json.loads((packet/'input.json').read_text())
for d in inp['sources'].values():assert hashlib.sha256((root/d['path']).read_bytes()).hexdigest()==d['sha256']
source=json.loads((root/inp['sources']['parent']['path']).read_text());bank=json.loads((root/inp['sources']['bank']['path']).read_text());b=bank['bisections'][0]
q=list(map(int,b['branch']['numerator_coefficients']));assert q==inp['branch']
def ev(cs,t,p):
 v=0
 for c in reversed(cs):
  z=Fraction(c);v=(v*t+z.numerator*pow(z.denominator,-1,p))%p
 return v
def group_quotient(a,b,p):
 squares={}
 for y in range(p):squares.setdefault(y*y%p,[]).append(y)
 points=[None]+[(x,y) for x in range(p) for y in squares.get((x**3+a*x+b)%p,[])]
 triples={add(add(P,P,a,p),P,a,p) for P in points};labels={P:0 for P in triples};dim=0
 while len(labels)<len(points):
  P=next(P for P in points if P not in labels);PP=add(P,P,a,p);old=list(labels.items())
  for Q,v in old:
   for k,R in [(1,P),(2,PP)]:
    z=add(Q,R,a,p);assert z not in labels;labels[z]=v+k*3**dim
  dim+=1
 assert len(labels)==len(points)==len(triples)*3**dim and dim<=2
 return labels,dim,len(points),len(triples)
def fibre(p,t):
 assert S.isprime(p) and p>3 and q[2]%p and (q[1]**2-4*q[0]*q[2])%p
 assert ev(q,t,p)==0
 model=source['weierstrass_model'];a=ev(model['A_coefficients_low_to_high'],t,p);b=ev(model['B_coefficients_low_to_high'],t,p)
 assert (4*a**3+27*b*b)%p
 return a,b
piv={};checked=[]
for row in inp['records']:
 p,t=row['prime'],row['t'];a,c=fibre(p,t);pts=[]
 for sec in source['sections']['records']:
  pts.append(tuple(ev(sec[k]['numerator_coefficients_low_to_high'],t,p)*pow(ev(sec[k]['denominator_coefficients_low_to_high'],t,p),-1,p)%p for k in ['X','Y']))
 assert all((y*y-x*x*x-a*x-c)%p==0 for x,y in pts)
 labels,dim,n,nt=group_quotient(a,c,p);rr=[[labels[P]//3**j%3 for P in pts] for j in range(dim)]
 assert [list(P) for P in pts]==row['points'] and rr==row['rows']
 assert (a,c,n,nt)==(row['A'],row['B'],row['order'],row['triple_subgroup_order'])
 for v in rr:
  z=v[:]
  for i in range(17):
   if z[i]==0:continue
   if i in piv:
    h=z[i];z=[(a-h*b)%3 for a,b in zip(z,piv[i])]
   else:
    h=pow(z[i],-1,3);piv[i]=[h*a%3 for a in z];break
 assert len(piv)==row['rank'];checked.append({'prime':p,'t':t,'rank':len(piv)})
assert len(piv)==17
w=inp['torsion_witness'];a,c=fibre(w['prime'],w['t']);_,_,n,_=group_quotient(a,c,w['prime']);assert n==w['order']==44 and n%3
# Exact quadratic residue-field arithmetic: the branch half has twice P16-P1.
t=S.symbols('t')
def pol(cs):return S.Poly(sum(S.Rational(c)*t**i for i,c in enumerate(cs)),t)
Q=pol(q).monic();A=pol(source['weierstrass_model']['A_coefficients_low_to_high'])%Q;B=pol(source['weierstrass_model']['B_coefficients_low_to_high'])%Q;zero=S.Poly(0,t)
def field_add(P,R):
 if P is None:return R
 if R is None:return P
 x,y=P;u,v=R
 if x==u and ((y+v)%Q).is_zero:return None
 m=((3*x*x+A)*S.invert(2*y,Q) if P==R else (v-y)*S.invert(u-x,Q))%Q
 z=(m*m-x-u)%Q;return z,(m*(x-z)-y)%Q
pts=[]
for i in [1,16]:
 r=source['sections']['records'][i]
 pts.append(tuple((pol(r[k]['numerator_coefficients_low_to_high'])*S.invert(pol(r[k]['denominator_coefficients_low_to_high']),Q))%Q for k in ['X','Y']))
half=tuple(pol(b['lifted_section'][k+'_coefficients'])%Q for k in ['x0','y0'])
assert ((half[1]**2-half[0]**3-A*half[0]-B)%Q).is_zero
trace=field_add(pts[1],(pts[0][0],-pts[0][1]));assert trace is not None and field_add(half,half)==trace
assert b['section_basis_w']==[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]
out={'schema':'q80-branch-saturation-defect-v1','status':'PASS','input_sha256':hashlib.sha256((packet/'input.json').read_bytes()).hexdigest(),
 'modulo3_rank':17,'places':checked,'no_3_torsion':w,
 'half_coordinates':[[str(c) for c in reversed(v.all_coeffs())] for v in half],
 'double_half':'P16-P1','specialization_kernel_rank':0,
 'first_halving_defect_dimension':1,
 'boundary':'Integral specialization injectivity uses modulo3 rank17 and no3-torsion. The exact half lies outside the inherited image but inside its rational span. Dimension1 uses the prior modulo2 kernel certificate and absence of branch-field2-torsion; this is not a new rank18 assertion for the number-field fibre, a full2-saturation index or a two-gain construction.'}
path=packet/'result.json'
if '--write' in __import__('sys').argv:path.write_text(json.dumps(out,indent=2)+'\n')
else:assert out==json.loads(path.read_text())
print('PASS: rank17 specialization; no3-torsion; exact noninherited half and one-dimensional first halving defect')
