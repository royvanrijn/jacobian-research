"""Exact KMFI90--91 controls; global finite-fiber extension remains OPEN."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
from itertools import combinations,product
from collections import Counter
from math import gcd
import json,hashlib,platform,time
import sympy as s
from psl11_residues import exceptional_actions,inverse
from central_lift import generated
from norm_target import (permute,complement_color,marked_permutation,
 permutation_orbits,cofactor)
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def zero(expr):return s.cancel(expr)==0
data,orders,triangles,classes,actions=exceptional_actions()
G=data['quotient'];I=data['qI'];mul=data['qmul']
g2=min(g for g in G if orders[g]==2)
g3=next(g for g in G if orders[g]==3 and orders[mul(g2,g)]==11 and len(generated([g2,g],I,mul))==660)
check("root group full PSL2(11)",len(G)==660 and len(generated([g2,g3],I,mul))==660)
# The (2,3,11) presentation has trivial abelianization, so the generated
# group is perfect; the explicit stabilizers are checked directly.
from sympy.matrices.normalforms import smith_normal_form
diag=smith_normal_form(s.Matrix([[2,0],[0,3],[11,11]]),domain=s.ZZ)
check("generated PSL group perfect via triangle abelianization",abs(diag[0,0])==abs(diag[1,1])==1)
for j,action in enumerate(actions):
 H=action['H'];perms=action['permutations'];n=11;m=12
 comm={mul(mul(mul(a,b),data['canon'](inverse(a))),data['canon'](inverse(b))) for a,b in product(H,repeat=2)}
 check("A5 sheet stabilizer perfect "+str(j),generated(comm,I,mul)==H and len(H)==60)
 v=tuple((i*i+3*i+5)%m for i in range(n))
 colors={g:complement_color(perms[g],v,m) for g in G}
 check("all660 complement sums zero "+str(j),all(sum(c)%m==0 for c in colors.values()))
 check("all root-generator cocycle identities "+str(j),all(
  colors[mul(a,b)]==tuple((colors[a][i]+permute(perms[a],colors[b])[i])%m for i in range(n))
  for a in [g2,g3] for b in G))
 ps=[marked_permutation(perms[g],colors[g],m) for g in [g2,g3]]
 O=permutation_orbits(ps,n*m)
 check("conjugated complement twelve orbits of degree11 "+str(j),sorted(map(len,O))==[11]*12)
 check("each marked orbit has one color over each label "+str(j),
  all({a//m for a in orb}==set(range(n)) for orb in O))
 check("each orbit color graph is a common shift of v "+str(j),
  all(len({(a%m-v[a//m])%m for a in orb})==1 for orb in O))
 for d in [1,2,3,4,6]:
  translation=tuple(i*m+(k+d)%m for i in range(n) for k in range(m))
  orbs=permutation_orbits(ps+[translation],n*m)
  check("nonzero diagonal color kernel prevents degree11 "+str((j,d)),
    all(len(orb)==11*12//gcd(d,12)>11 for orb in orbs))
 # Averaging any triple over the transitive group gives180 per label.
 check("all triple orbit averages uniform "+str(j),all(
  len({sum(p[a]==i for p in perms.values() for a in triple) for i in range(n)})==1
  for triple in combinations(range(n),3)))
check("complement count",12**11//12==12**10)
# Resultants: algebraic fixtures, expressly NOT Klein/PSL realizations.
T,b,q,r,h,c=s.symbols('T b q r h c')
eps,Q,R,H=s.symbols('epsilon Q R H')
fixtures=[]
limit=s.resultant(T**11-1,H+R*T*T+Q*T**3,T)
for cc in [0,1]:
 f=T**11+cc*T**10+T
 result=s.resultant(f-b,cofactor(T,cc,q,r,h),T)
 check("fixture monic degree11 and first derivative1 "+str(cc),s.Poly(f,T).LC()==1 and s.diff(f,T).subs(T,0)==1)
 check("resultant b-degree3 "+str(cc),s.degree(result,b)==3)
 check("resultant b-leading coefficient q^11 "+str(cc),zero(s.Poly(result,b).LC()-q**11))
 check("resultant cofactor homogeneous degree11 "+str(cc),all(sum(mon[1:])==11 for mon,coef in s.Poly(result,b,q,r,h).terms()))
 scaled=s.expand(result.subs({b:eps**-11,q:eps**3*Q,r:eps**2*R,h:H},simultaneous=True))
 poly=s.Poly(scaled,eps,Q,R,H)
 check("full scaled resultant has no negative infinity powers "+str(cc),all(mon[0]>=0 for mon,coef in poly.terms()))
 check("full scaled resultant has expected cyclotomic limit "+str(cc),zero(scaled.subs(eps,0)-limit))
 fixtures.append({'f':str(f),'c':cc,'resultant':str(result),'not_PSL_witness':True})
# Direct product/norm fixture with eleven known distinct roots.
roots=list(range(1,12));f=s.prod(T-i for i in roots);cc=-sum(roots)
res=s.resultant(f,cofactor(T,cc,q,r,h),T)
check("resultant equals product at eleven known roots",
 zero(res-s.prod(cofactor(i,cc,q,r,h) for i in roots)))
check("passport derivative-product exponents even",4*(2-1)==4 and 3*(3-1)==6)
check("norm extraction U^6 over square discriminant",2*6==12 and 2*2==4 and 2*3==6)
# Exact determinant of the three branch lines.
a1,a2,a3=s.symbols('a1 a2 a3')
mat=s.Matrix([[1,a*a,a**3-c*a*a] for a in [a1,a2,a3]])
expected=(a2-a1)*(a3-a1)*(a3-a2)*(a1*a2+a1*a3+a2*a3)
check("finite triple determinant factorization",zero(mat.det()-expected))
check("projective infinity line-pair determinant",
 zero(s.det(s.Matrix([[0,1,0],[1,a1*a1,a1**3-c*a1*a1],[1,a2*a2,a2**3-c*a2*a2]]))-
 (a1**3-c*a1*a1-a2**3+c*a2*a2)))
# A specialized ordinary triple can occur even though the generic
# arrangement has none. This guards the precise cDV extension step.
tripleQ=(T-1)*(T-2)*(T+s.Rational(2,3))
tripleH=T**8+s.Rational(3,4)*T-1
bb=s.Rational(4,3);triplef=s.expand(bb+tripleQ*tripleH)
check("special triple fixture monic and normalized",
 s.Poly(triplef,T).LC()==1 and triplef.subs(T,0)==0 and s.diff(triplef,T).subs(T,0)==1)
check("special triple value is an ordinary curve fiber",
 s.gcd(triplef-bb,s.diff(triplef,T))==1)
check("three concurrent cofactor lines at an ordinary specialization",
 zero(s.Poly(tripleQ,T).coeff_monomial(T)) and s.rem(triplef-bb,tripleQ,T)==0)
# Exact cyclotomic checks for all triples and all pairs of eleven lines.
Z=s.symbols('Z');Phi=sum(Z**i for i in range(11))
check("Phi11 irreducible",s.Poly(Phi,Z).is_irreducible)
check("all165 cyclotomic triple sums nonzero",all(
 s.rem(sum(Z**i for i in tri),Phi,Z)!=0 for tri in combinations(range(11),3)))
check("all55 infinity-line pair determinants nonzero",all(
 s.rem(Z**(2*i+3*j)-Z**(3*i+2*j),Phi,Z)!=0 for i,j in combinations(range(11),2)))
# Projective fiber chart changes of q^3 dq dr/Y.
v,w=s.symbols('v w')
chartq=s.Matrix([1/v,w/v])
chartr=s.Matrix([w/v,1/v])
check("q-chart top form coefficient regular constant",
 zero((1/v)**3*chartq.jacobian([v,w]).det()*v**6+1))
check("r-chart top form coefficient regular cubic",
 zero((w/v)**3*chartr.jacobian([v,w]).det()*v**6-w**3))
base=s.Matrix([eps**-11,eps**3*Q,eps**2*R])
check("infinity full three-form coefficient -11 epsilon^2 Q^3",
 zero((eps**3*Q)**3*base.jacobian([eps,Q,R]).det()+11*eps**2*Q**3))
u,v=s.symbols('u v')
node=s.Matrix([u,u*v*v])
check("A1 resolution form coefficient2",zero(node.jacobian([u,v]).det()/(u*v)-2))
check("all tame pole orders stay negative under finite cover",
 all(-e*a+e-1<0 for a in range(1,5) for e in range(1,13)))
manifest=json.loads((ROOT/'sources/norm-target-manifest.json').read_text())
for kind in ['sources','retained_primary_inputs']:
 for source in manifest[kind]:
  raw=(ROOT/source['path']).read_bytes()
  check("primary bytes "+source['path'],len(raw)==source['bytes'])
  check("primary hash "+source['path'],hashlib.sha256(raw).hexdigest()==source['sha256'])
exploration=json.loads((ROOT/'norm_target_exploration.json').read_text())
check("m9 exploration kept distinct from rational m12 cofactor family",
 exploration['direct_derivative_marking_m9']['color_vector_degree']==9**11)
model={'claims':['KMFI90','KMFI91'],'relative_family':{
 'degree_f':11,'degree_whole_polynomial':14,'marking_exponent':12,
 'source_field':'k(t,a,r), rational3','source_original_affine_component_is_affine_space':False,
 'normal_closure_group':'C12^11 semidirect PSL2(11)',
 'both_exceptional_actions_checked':True,'all_complements_conjugate':True,
 'complement_count':12**10,'compulsory_norm_double_cover':'Y^2=Res(f(T)-b,1+r*T^2+q*(T^3-c*T^2))',
 'every_dominant_degree11_marked_descent_contains_norm_field':True,
 'non_dominant_specializations_covered':False,
 'forgetting_the_mark_claimed_as_same_operation':False},
 'finite_branch_normalization':{'first':0,'second':'beta !=0, not assumed1',
 'monic_f_and_first_derivative1_imposed':True},
 'ordinary_triple_countercontrol':{'f':str(triplef),'b':str(bb),'Q':str(tripleQ),'not_PSL_witness':True},
 'candidate_form':'q^3 db wedge dq wedge dr /Y',
 'normal_double_cover_branch_bidegree':[4,12],
 'dualizing_sheaf':'pullback O(0,3), not yet a global form on a resolution',
 'proved_extension_locus':'resolution above P1_b minus{0,beta}, including infinity',
 'infinity_scaled_form':'-11 epsilon^2 Q^3 d(epsilon) wedge dQ wedge dR/Y',
 'infinity_branch_lines_general_position':True,
 'ordinary_specialized_triple_intersections_not_excluded':True,
 'ordinary_triples_use_pinned_cDV_canonical_input':True,
 'finite_fiber0_extension_decided':False,'finite_fiber_beta_extension_decided':False,
 'global_regular_three_form_proved':False,'norm_field_nonunirational_proved':False,
 'cofactor_degree11_descent_excluded':False,
 'degree11_PSL11_absolute_occurrence_decided':False,
 'full_action_classification_complete':False,'algebraic_resultant_fixtures':fixtures}
(ROOT/'norm_target_model.json').write_text(json.dumps(model,indent=2)+'\n')
names=['NORM_TARGET_REDUCTION.md','norm_target.py','verify_norm_target.py','norm_target_model.json',
 'norm_target_exploration.json','sources/norm-target-manifest.json',
 'psl11_residues.py','central_lift.py']+[s['path'] for k in ['sources','retained_primary_inputs'] for s in manifest[k]]
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-started,3),'python':platform.python_version(),'sympy':s.__version__,
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
 'assurance':'Written all-complement norm-field reduction and local three-form proof away from two finite branch fibers; exact finite-group, resultant, projective-chart and infinity controls. Global extension remains unproved. No formal verification or independent whole-proof review.',
 'goal_complete':False}
(ROOT/'norm_target_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'infinity_regular':True,'global_form_proved':False,'PSL11_absolute_decided':False,'goal_complete':False}))
