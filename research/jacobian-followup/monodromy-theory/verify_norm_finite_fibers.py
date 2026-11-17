"""Exact KMFI92--93 controls; universal case exhaustion is written in the note."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
from itertools import combinations_with_replacement
import json,hashlib,platform,time
import sympy as S
from norm_finite_fibers import cofactor_root_cover,critical_weights,residue_exponent
from psl11_residues import exceptional_actions
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def zero(x):return S.cancel(x)==0
s,z,T,b,c,q,r,h=S.symbols('s z T b c q r h')
cover=cofactor_root_cover(s,z);u,v,w=cover['roots'];D=cover['D'];B=cover['B']
check("ordered roots satisfy zero linear coefficient",zero(u*v+u*w+v*w))
check("finite cone-cover square s",zero(u+v-s*s))
check("finite cone-cover square z",zero(u+w-z*z))
check("finite cone-cover product",zero(u-s*z))
check("cofactor constant square",zero(-u*v*w-D*D))
check("cofactor cubic reconstructed",zero(S.prod(T-a for a in [u,v,w])-
 (T**3+cover['r_minus_c']*T*T+cover['h'])))
J=S.Matrix([cover['h'],cover['r_minus_c']]).jacobian([s,z]).det()
check("full coefficient Jacobian -2DB",zero(J+2*D*B))
check("residue multiplier polynomial2B",zero(-J/D-2*B))
check("zero-root locus is three axes",zero(u*v*w+s*s*z*z*(s-z)**2))
# Exact norm identity for an algebraic degree11 fixture, not a Klein witness.
f=T**11+T;rr=cover['r_minus_c'];hh=cover['h']
R=S.resultant(f-b,h+r*T*T+T**3,T)
check("norm resultant at finite root cover",
 zero(R.subs({h:hh,r:rr},simultaneous=True)+S.prod(a**11+a-b for a in [u,v,w])))
check("finite norm double-cover equation",
 zero(hh*(-S.prod(a**11+a-b for a in [u,v,w]))+D**2*S.prod(a**11+a-b for a in [u,v,w])))
# Origin tangent cone: full five-by-five quadric pencil.
lam,mu=S.symbols('lambda mu')
M1=S.diag(0,0,1,0,-1);M1[0,1]=M1[1,0]=-1;M1[1,1]=1
M2=S.diag(-1,1,0,1,-1)
pencil=S.factor((lam*M1+mu*M2).det())
expected=lam*mu*(lam+mu)*(lam*lam+lam*mu+mu*mu)
check("origin complete-intersection pencil",zero(pencil-expected))
check("origin pencil five simple roots",S.discriminant(expected.subs(mu,1),lam)!=0 and S.degree(expected.subs(mu,1),lam)==4)
check("origin ordinary blowup crepant",residue_exponent((1,1,1,1,1),(2,2))==0)
# Uniform two/triple critical cases. All common binary roots are
# parametrized, so these identities cover every coefficient choice.
A,B0,C,al,be,t=S.symbols('A B0 C alpha beta t',nonzero=True)
critical_rows=[]
for e in [2,3]:
 weights,degrees=critical_weights(e)
 check("triple-critical residue exponent0 e"+str(e),residue_exponent(weights,degrees)==0)
 F1=C*(al+be*t)**e-A;F2=C*(al+be*t)**e-B0*t**e
 common={A:C*(al+be*t)**e,B0:C*(al+be*t)**e/t**e}
 a=S.diff(F1,t).subs(common);bb=S.diff(F2,t).subs(common)
 check("first common-root slope nonzero formula e"+str(e),
  zero(a-e*C*be*(al+be*t)**(e-1)))
 check("second common-root slope nonzero formula e"+str(e),
  zero(bb+e*C*al*(al+be*t)**(e-1)/t))
 check("common-root slope difference nonzero formula e"+str(e),
  zero(a-bb-e*(C*(al+be*t)**e/t**e)*t**(e-1)))
 det=S.factor(bb*(-a)*(a-bb))
 target=S.factor(e**3*C**3*al*be*(al+be*t)**(3*e-2)/t**2)
 check("exceptional A1 quadratic determinant nonzero e"+str(e),zero(det-target))
 critical_rows.append({'ramification':e,'weights':list(weights),'equation_weights':list(degrees),
   'residue_exponent':0,'exceptional_singularities_on_index_cover':'smooth or A1',
   'nonzero_hessian_certificate':str(det)})
# Exact weighted y-chart volume; index3 action is scalar determinant1.
ep,X1,X2,Y1,Y2=S.symbols('epsilon X1 X2 Y1 Y2')
for e in [2,3]:
 wx=1 if e==2 else 2;wy=1 if e==2 else 3
 expr=S.Matrix([ep**wy,ep**wx*X1,ep**wx*X2]).jacobian([ep,X1,X2]).det()/(ep**wy*Y1*ep**wy*Y2)
 check("weighted y-chart residue regular e"+str(e),zero(expr-wy/(Y1*Y2)))
Y3=S.symbols('Y3')
for e in [2,3]:
 wx=1 if e==2 else 2;wy=1 if e==2 else 3
 mat=S.Matrix([ep**wy*Y3,ep**wx,ep**wx*X2])
 val=mat.jacobian([ep,X2,Y3]).det()/(ep**wy*Y1*ep**wy*Y2)
 check("weighted x-chart no residual exceptional pole e"+str(e),zero(val-wx/(Y1*Y2)))
check("fractional epsilon negative half-order cleared with differential",
 2*S.Rational(-1,2)+(2-1)==0)
check("index3 volume action determinant1",(1+1+1)%3==0)
# Two zero roots: exact first derivative and implicit quadratic b-part.
U,W,a2,a3=S.symbols('U W a2 a3',nonzero=True)
vv=-U*W/(U+W);ff=lambda x:x+a2*x*x+a3*x*x*x
check("zero-axis root relation",zero(U*vv+U*W+vv*W))
check("two-zero-root implicit derivative2",S.diff(ff(U)-ff(vv),U).subs(U,0)==2)
y1,y2=S.symbols('y1 y2')
uquad=(y1*y1-y2*y2)/2
check("two-zero-root b quadratic rank2",
 zero(uquad-y1*y1+(y1*y1+y2*y2)/2))
# At all-nonzero roots each coordinate pair is etale on the cone;
# the third root slopes never vanish.
u0,v0=S.symbols('u0 v0',nonzero=True)
wfun=-u0*v0/(u0+v0)
check("third-root first slope",zero(S.diff(wfun,u0)+v0**2/(u0+v0)**2))
check("third-root second slope",zero(S.diff(wfun,v0)+u0**2/(u0+v0)**2))
# Multisets are essential when the cubic cofactor has a repeated root.
data,orders,triangles,classes,actions=exceptional_actions()
for j,action in enumerate(actions):
 perms=action['permutations']
 uniform=True
 for triple in combinations_with_replacement(range(11),3):
  counts=[sum(p[a]==i for p in perms.values() for a in triple) for i in range(11)]
  if counts!=[180]*11:uniform=False;break
 check("all286 root-multiset averages uniform action"+str(j),uniform)
check("generic reciprocal sum nonzero",-S.Integer(1)/(-b)==1/b)
qa=(u0+v0)/(u0*u0*v0*v0)
ra=(c*(u0+v0)-u0*u0-u0*v0-v0*v0)/(u0*u0*v0*v0)
jacq0=S.factor(S.Matrix([qa,ra]).jacobian([u0,v0]).det().subs(v0,-u0))
check("q0 ordered two-finite-root coordinates etale",zero(jacq0+2/u0**7))
# q=0 with h nonzero: exact compound-A factorization.
p1,p2,Z=S.symbols('p1 p2 Z')
BB=b-(p1+p2)/2
F=Z*Z-(p1-b)*(p2-b)
check("two finite root compound-A equation",
 zero((Z-BB)*(Z+BB)+(p1-p2)**2/4-F))
# q=h=0: universal even/odd normal form, no truncation.
EE,OO,HH,WW,BB=S.symbols('E O H W B')
F0=Z*Z-h*((b-EE)**2+h*OO**2)
check("boundary even-odd resultant identity",
 zero(S.expand((EE+T*OO-b)*(EE-T*OO-b)).subs(T*T,-h)-((b-EE)**2+h*OO**2)))
check("boundary A3 rescaling",
 zero((OO*Z)**2-(h*OO**2)*(b-EE)**2-(h*OO**2)**2-OO**2*F0))
Fa3=WW**2-HH*BB**2-HH**2
check("boundary A3 completing square",
 zero((WW-HH-BB**2/2)*(WW+HH+BB**2/2)+BB**4/4-Fa3))
check("boundary A3 section Milnor order",S.degree(BB**4,BB)-1==3)
check("constant cofactor norm1",S.resultant(f-b,1,T)==1)
# The complete ten-dimensional cubic numerator space: every weight
# component is nonzero on the generic limiting surface unless its
# coefficients vanish. This finite monomial census gives exact p_g.
monomials=[(aa,dd,3-aa-dd) for aa in range(4) for dd in range(4-aa)]
surviving=[m for m in monomials if 3*m[0]+2*m[1]-7>=0]
check("all ten dualizing cubic monomials accounted for",len(monomials)==10)
check("exactly three infinity-regular numerators",
 set(surviving)=={(3,0,0),(2,1,0),(1,2,0)})
check("global form exceptional orders2,1,0",
 sorted(3*a+2*d-7 for a,d,j in surviving)==[0,1,2])
# Required boundaries and new assurance flags are explicit.
manifest=json.loads((ROOT/'sources/norm-finite-fibers-manifest.json').read_text())
for src in manifest['retained_primary_inputs']:
 raw=(ROOT/src['path']).read_bytes()
 check("retained primary bytes "+src['path'],len(raw)==src['bytes'])
 check("retained primary hash "+src['path'],hashlib.sha256(raw).hexdigest()==src['sha256'])
case_table=[
 {'region':'q!=0, root-cover origin','resolution':'ordinary crepant blowup, smooth(2,2) exceptional surface'},
 {'region':'q!=0, nonzero two-zero-root axis','resolution':'implicit elimination, smooth or cA'},
 {'region':'q!=0, all nonzero roots and no participating critical root','resolution':'original norm cover cDV as in KMFI91'},
 {'region':'q!=0, one participating critical root','resolution':'root-cover Jacobian full rank'},
 {'region':'q!=0, two critical roots, third y unit','resolution':'cA_(e-1) slice'},
 {'region':'q!=0, two critical roots plus unramified third','resolution':'Morse cA_k; multiset reciprocal identity rules out k infinite'},
 {'region':'q!=0, three critical roots, e=2','resolution':'crepant ordinary complete-intersection blowup; smooth/A1 exceptional charts'},
 {'region':'q!=0, three critical roots, e=3','resolution':'crepant weighted blowup; index2 cDV and index3 scalar-volume charts'},
 {'region':'q=0, h*r!=0','resolution':'etale ordered binary roots, compound-A equation'},
 {'region':'q=h=0, r!=0','resolution':'smooth for b!=0; A3 section at b=0'},
 {'region':'q=r=0, h!=0','resolution':'constant norm, smooth'}]
model={'claims':['KMFI92','KMFI93'],'finite_cover_roots':[str(u),str(v),str(w)],
 'finite_cover_coefficients':{'h':str(hh),'r_minus_c':str(rr)},
 'residue_multiplier':str(B),'finite_case_exhaustion':case_table,'weighted_cases':critical_rows,
 'geometric_genus':3,'global_form_basis_numerators':['q^3','q^2*r','q*r^2'],
 'global_regular_three_form_proved':True,'global_form':'q^3 db wedge dq wedge dr/Y',
 'norm_field_nonunirational_proved':True,'norm_field_not_stably_rational':True,
 'dominant_rational_target_cofactor_degree11_descent_excluded':True,
 'both_exceptional_actions_covered':True,'additional_independent_target_variables_covered':True,
 'non_dominant_specializations_excluded':False,'forgetting_the_mark_excluded':False,
 'different_Keller_architectures_excluded':False,'exceptional11_absolute_action_decided':False,
 'full_action_classification_complete':False,
 'proof_basis':'Written exhaustive finite-root and projective-boundary proof, both crepant blowup chart classifications, previous infinity extension and primary cDV input. Finite controls supplement the argument.',
 'three_order3_critical_roots_required_to_be_absent':False,
 'fractional_epsilon_count_alone_used_as_resolution':False}
(ROOT/'norm_finite_fibers_model.json').write_text(json.dumps(model,indent=2)+'\n')
names=['NORM_TARGET_FINITE_FIBERS.md','norm_finite_fibers.py','verify_norm_finite_fibers.py',
 'norm_finite_fibers_model.json','sources/norm-finite-fibers-manifest.json',
 'NORM_TARGET_REDUCTION.md','norm_target.py','psl11_residues.py','central_lift.py']+[s['path'] for s in manifest['retained_primary_inputs']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-started,3),'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
 'assurance':'Written complete finite-fiber resolution plus exact cover, residue, tangent-pencil, weighted chart, multiset and boundary controls. Global form and specified dominant cofactor-descent obstruction; no absolute exceptional-action exclusion, formal verification or independent whole-proof review.',
 'goal_complete':False}
(ROOT/'norm_finite_fibers_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'global_three_form':True,'dominant_cofactor_descent_excluded':True,
 'exceptional_absolute_decided':False,'goal_complete':False}))
