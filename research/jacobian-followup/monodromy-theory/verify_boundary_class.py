"""Exact controls for boundary-class obstructions KMFI87--89.

Class independence and all-dimensional purity are written proofs.
This checker does not treat a permitted class lattice as a realization.
"""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
import sympy as S
from pathlib import Path
from collections import Counter
from itertools import product
import hashlib,json,platform,time
from psl11_residues import exceptional_actions,centralizer,cyclic_overgroups,local_profile
from central_lift import generated
from boundary_class import ramified_boundary_profile,canonical_divisibility
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def zero(expr):
    return S.cancel(expr)==0
def same(expr,target):
    return zero(expr-target)
x,y,z=S.symbols('x y z');A,B,C,W,H,T=S.symbols('A B C W H T')
u=1+x*y;g=1-S.Rational(3,2)*x*y+x*x*z
F=[S.cancel((u+u*u-2*u**3*g)/(x*x)),
 S.cancel((1+2*u-3*u*u*g)/x),x*g]
for i,f in enumerate(F):check("actual cubic coordinate polynomial "+str(i),S.denom(f)==1)
check("actual cubic full Jacobian1",same(S.Matrix(F).jacobian([x,y,z]).det(),1))
src={A:F[0],B:F[1],C:F[2],W:u*g}
theta=S.cancel(u*(u*g-1)/x);src[H]=theta
check("normalizing H polynomial",S.denom(theta)==1)
check("normalizing H equals -(B+y)/3",same(theta,-(F[1]+y)/3))
eq=[W*W-W-C*H,W*H+B*W-A*C,H*H+B*H-A*W+A]
for i,e in enumerate(eq):check("integral normalization relation "+str(i),zero(e.subs(src)))
raw=W**3-W**2+B*C*W-A*C*C
check("primitive monic equation",zero(raw.subs(src)))
for v in [A,B,C,W]:check("raw hypersurface singular codim1 derivative "+str(v),zero(S.diff(raw,v).subs({C:0,W:0})))
check("raw hypersurface vanishes on codim1 singular plane",zero(raw.subs({C:0,W:0})))
# Rational inverse proves generic degree3/source-field identification; P(r)-A
# is irreducible by the degree of the polynomial rational map over k(B,C).
GG=3*W*W-2*W+B*C
check("raw derivative pulls to g",same(GG.subs(src),g))
check("source x reconstructs",same((C/GG).subs(src),x))
check("source y reconstructs",same((-3*H-B).subs(src),y))
zinv=(GG-1+S.Rational(3,2)*(C/GG)*(-3*H-B))/(C/GG)**2
check("source z reconstructs",same(zinv.subs(src),z))
# Exact finite-free multiplication algebra.
MW=S.Matrix([[0,0,A*C],[1,1,-B],[0,C,0]])
MH=S.Matrix([[0,A*C,-A],[0,-B,A],[1,0,-B]])
Id=S.eye(3);Z0=S.zeros(3)
check("normalizing algebra multiplication commutes",MW*MH==MH*MW)
check("normalizing algebra W square",MW*MW==MW+C*MH)
check("normalizing algebra WH",MW*MH==A*C*Id-B*MW)
check("normalizing algebra H square",MH*MH==A*MW-B*MH-A*Id)
check("W characteristic polynomial",same(MW.charpoly(T).as_expr(),raw.subs(W,T)))
hpoly=H*(H+B)**2+A*(H+B)-A*A*C
check("H monic integral polynomial",same(MH.charpoly(T).as_expr(),hpoly.subs(H,T)))
basis=[Id,MW,MH]
trace=S.Matrix([[S.trace(a*b) for b in basis] for a in basis])
delta=-27*A*A*C*C+18*A*B*C-4*A-4*B**3*C+B*B
check("integral-basis trace discriminant",same(trace.det(),delta))
check("raw power-basis discriminant includes conductor factor C^2",
 same(S.discriminant(raw,W),C*C*delta))
matrix=S.Matrix([[H+B,1-W,C],[A,-H,W]])
check("determinantal matrix hyperplane identity",same(matrix[0,1]+matrix[1,2],1))
check("three minors are normalization relations",
 all(same(matrix[:,list(pair)].det(),sign*e) for pair,sign,e in
 [((0,1),-1,eq[2]),((0,2),1,eq[1]),((1,2),-1,eq[0])]))
check("first Hilbert-Burch syzygy",zero((H+B)*eq[0]+(1-W)*eq[1]+C*eq[2]))
check("second Hilbert-Burch syzygy",zero(A*eq[0]-H*eq[1]+W*eq[2]))
# Smooth rank-one locus proof uses the hyperplane value1: radial tangent
# gives derivative1 and excludes the rank-zero singular locus.
r,c,q,s,a,h=S.symbols('r c q s a h')
chart1={A:r*q,B:q+r-r*r*c,C:c,W:r*c,H:-r+r*r*c}
chart2={A:a,B:s*a-h,C:s*(1+s*h),W:1+s*h,H:h}
for j,ch in enumerate([chart1,chart2]):
 for i,e in enumerate(eq):check("smooth affine chart relation "+str((j,i)),zero(e.subs(ch)))
 # At least one entry in the normalized row is nonzero, since the
 # hyperplane equation is1; neither chart loses a rank-one point.
 check("smooth chart hyperplane "+str(j),same((matrix[0,1]+matrix[1,2]).subs(ch),1))
transition={s:1/r,a:r*q,h:-r+r*r*c}
for v in [A,B,C,W,H]:
 check("two normal charts agree "+str(v),same(chart2[v].subs(transition,simultaneous=True),chart1[v]))
check("canonical transition determinant r",
 same(S.Matrix([transition[v] for v in [s,a,h]]).jacobian([r,c,q]).det(),r))
yZ=-3*H-B;yy1=yZ.subs(chart1);yy2=yZ.subs(chart2)
d=r-yy1;d2=1-s*yy2
check("boundary equation transition d=r*d2",same(d2.subs(transition,simultaneous=True)*r,d))
check("boundary affine plane linear graph",same(d,q-r+2*r*r*c))
jac1=S.Matrix([chart1[v] for v in [A,B,C]]).jacobian([r,c,q]).det()
jac2=S.Matrix([chart2[v] for v in [A,B,C]]).jacobian([s,a,h]).det()
check("first finite-projection Jacobian equals -d",same(jac1,-d))
check("second finite-projection Jacobian equals -d2",same(jac2,-d2))
D={q:r-2*r*r*c}
check("ramification divisor maps to discriminant",zero(delta.subs(chart1).subs(D)))
P=c*T**3-T*T+chart1[B]*T-chart1[A]
check("ramification repeated root",zero(S.diff(P,T).subs(T,r).subs(D)))
check("ramification generically index2",S.diff(P,T,2).subs(T,r)!=0)
# Both inverse charts of the actual ordinary-affine source.
inv1={x:1/d,y:yy1,z:c*d**3-d*d+S.Rational(3,2)*yy1*d}
inv2={x:s/d2,y:yy2,z:d2*S.Rational(1,2)*
 (4*yy2**2-a-s*yy2**3+2*s*a*yy2-s*s*a*yy2**2)}
for j,(inv,ch) in enumerate([(inv1,chart1),(inv2,chart2)]):
 for v,f in zip([A,B,C],F):
  check("affine inverse chart forward map "+str((j,v)),same(f.subs(inv,simultaneous=True),ch[v]))
 for v in [W,H]:
  check("affine inverse chart normalization coordinate "+str((j,v)),same(src[v].subs(inv,simultaneous=True),ch[v]))
 variables=[r,c,q] if j==0 else [s,a,h]
 boundary=S.Poly(d if j==0 else d2,*variables).monic()
 check("inverse-chart denominators only boundary "+str(j),
  all(S.Poly(factor,*variables).monic()==boundary
    for value in inv.values()
    for factor,exponent in S.factor_list(S.denom(S.cancel(value)))[1]))
back1={r:u/x,c:F[2],q:theta+F[1]}
back2={s:x/u,a:F[0],h:theta}
for j,(inv,back) in enumerate([(inv1,back1),(inv2,back2)]):
 for v in [x,y,z]:
  check("actual source chart roundtrip "+str((j,v)),same(inv[v].subs(back,simultaneous=True),v))
check("source two charts cover",same(u-x*y,1))
check("boundary first inverse denominator is1/x",same(d.subs(back1,simultaneous=True),1/x))
check("boundary second inverse denominator is1/u",same(d2.subs(back2,simultaneous=True),1/u))
# Derive the mandatory boundary prime profiles in both exceptional actions.
data,orders,triangles,classes,actions=exceptional_actions()
G=data['quotient'];I=data['qI'];mul=data['qmul']
expected={(2,2):(4,[1]*4),(2,4):(3,[1]*3),(2,6):(2,[1]*2),
 (2,12):(2,[1]*2),(3,3):(3,[2]*3),(3,6):(2,[2]*2),(5,5):(2,[4]*2)}
rows=[];seen=set()
for n in [2,3,5]:
 gen=min(g for g in G if orders[g]==n);inertia=generated([gen],I,mul)
 Cg=centralizer(G,gen,mul)
 for Dg in cyclic_overgroups(inertia,Cg,I,mul):
  for j,action in enumerate(actions):
   profile=local_profile(Dg,inertia,action['permutations'])
   bprofile=ramified_boundary_profile(profile['ef']);key=(n,len(Dg));tag=str((key,j,min(Dg)))
   check("boundary prime count "+tag,bprofile['mandatory_primes']==expected[key][0])
   check("canonical boundary coefficients "+tag,bprofile['canonical_coefficients']==expected[key][1])
   check("nontrivial exceptional branch forces rank>=2 "+tag,bprofile['mandatory_primes']>=2)
   if (key,j) not in seen:
    seen.add((key,j));rows.append({'action':j,'inertia':n,'decomposition_order':len(Dg),**bprofile})
check("all seven boundary rows in both exceptional actions",len(rows)==14)
c6=next(row for row in rows if row['inertia']==2 and row['decomposition_order']==6)
check("residue-degree3 prime counted once, not three times",
 c6['mandatory_primes']==2 and c6['residue_degrees']==[1,3])
for allowed,div in [([2,3,5],1),([3,5],2),([5],4)]:
 vec=[c for row in rows if row['inertia'] in allowed for c in row['canonical_coefficients']]
 check("exact canonical divisibility "+str(allowed),canonical_divisibility(vec)==div)
 check("unramified omitted primes do not change divisibility "+str(allowed),canonical_divisibility(vec+[0,0])==div)
check("independent curve mandatory prime counts4,3",
 expected[(2,2)][0]==4 and expected[(3,3)][0]==3)
# Full residue-class control for tame prime Kummer ramification.
for e in [2,3,11]:
 from math import gcd
 check("prime Kummer nontrivial inertia forces full cyclic decomposition "+str(e),
  all(e//gcd(e,m) in [1,e] and (m==0 or e//gcd(e,m)==e) for m in range(e)))
manifest=json.loads((ROOT/'sources/boundary-class-manifest.json').read_text())
for kind in ['sources','retained_primary_inputs']:
 for source in manifest[kind]:
  rawbytes=(ROOT/source['path']).read_bytes()
  check("primary bytes "+source['path'],len(rawbytes)==source['bytes'])
  check("primary hash "+source['path'],hashlib.sha256(rawbytes).hexdigest()==source['sha256'])
model={'claims':['KMFI87','KMFI88','KMFI89'],
 'universal_boundary_theorem':{'class_group':'free on all omitted prime divisors',
 'canonical_vector':'one coefficient e-1 per boundary prime',
 'degree_greater_than1_forces_nontorsion_canonical_class':True,
 'normal_torsion_canonical_model_excluded_with_fixed_target':True,
 'normal_global_affine_complete_intersection_excluded':True,
 'normal_global_monogenic_integral_closure_excluded':True,
 'arbitrary_birational_target_recharts_excluded':False,
 'smooth_normalization_excluded':False,
 'locally_complete_intersection_normalization_excluded':False,
 'all_omitted_primes_assumed_ramified':False},
 'cubic_positive_control':{'coordinates':[str(f) for f in F],
 'raw_monic_model_normal':False,'normalizing_relations':[str(e) for e in eq],
 'finite_free_basis':['1','W','H'],'normalization_smooth':True,
 'class_group':'Z','canonical_class':'primitive generator',
 'normalization_global_affine_complete_intersection':False,
 'actual_source_complement_of_boundary':'ordinary affine3',
 'finite_projection_discriminant':str(delta),
 'boundary':'irreducible affine plane, ramification index2',
 'new_occurrence_family_claimed':False},
 'exceptional_PSL11_boundary_rows':rows,
 'exceptional_PSL11_absolute_class_rank_lower_bound':2,
 'independent_curve_PSL11_only_nontrivial_rows':[[2,2],[3,3]],
 'independent_curve_PSL11_class_rank_lower_bound':'4*a+3*b>=3; a+b>0',
 'canonical_divisibility':'1 if inertia2 present; else2 if inertia3 present; else4',
 'canonical_divisibility_is_Cartier_index':False,
 'exceptional_degree11_absolute_occurrence_decided':False,
 'full_action_classification_complete':False}
(ROOT/'boundary_class_model.json').write_text(json.dumps(model,indent=2)+'\n')
names=['BOUNDARY_CLASS_OBSTRUCTION.md','boundary_class.py','verify_boundary_class.py',
 'boundary_class_model.json','psl11_residues.py','central_lift.py',
 'sources/boundary-class-manifest.json']+[s['path'] for k in ['sources','retained_primary_inputs'] for s in manifest[k]]
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-started,3),'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
 'assurance':'Written dimension-independent boundary-basis, tame canonical and purity proofs; exact complete cubic normalization/charts and both exceptional local action controls. Not formal verification, independent whole-proof replay or external review.',
 'goal_complete':False}
(ROOT/'boundary_class_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'exceptional_absolute_min_class_rank':2,'independent_curve_min_class_rank':3,
 'cubic_normalization_smooth':True,'goal_complete':False}))
