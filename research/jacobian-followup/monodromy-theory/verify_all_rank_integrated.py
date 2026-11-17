"""Exact controls for KMFI112--114; the unbounded proofs are in the note."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
import sympy as S
import json,hashlib,time,platform
from integrated_sweeps import family
ROOT=Path(__file__).resolve().parent;start=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
t,A,B=S.symbols('t A B')
def integral(r,d):
    return S.expand(sum((-1)**j*S.binomial(r,j)*A**(r-j)*t**(d*j+1)/S.Integer(d*j+1) for j in range(r+1)))
# Arbitrary h and all endpoints at every differentiation step.
for r in range(1,9):
    hf=S.Function('h')(t);H=[t]+[S.Function('H'+str(j))(t) for j in range(1,r+1)]
    I=sum((-1)**j*S.binomial(r,j)*A**(r-j)*H[j] for j in range(r+1))
    replacement={S.diff(H[j],t):hf**j for j in range(1,r+1)}
    previous=I.subs(A,hf)
    for j in range(1,r+1):
        result=S.diff(previous,t).xreplace(replacement)/S.diff(hf,t)
        wanted=S.diff(I,A,j).subs(A,hf)
        check('arbitrary-h rank'+str(r)+' critical derivative step'+str(j),zero(result-wanted))
        previous=wanted
    check('rank'+str(r)+' critical derivative recovers t',zero(previous-S.factorial(r)*t))
discriminants=[]
for r,d in [(2,2),(2,3),(4,2),(4,3),(6,2),(6,3)]:
    I=integral(r,d);P=I-B;N=r*d+1;delta=(-1)**(r*d//2)*N
    res=S.resultant(P,A-t**d,t);disc=S.discriminant(P,t)
    check('rank'+str(r)+' degreeh'+str(d)+' exact signed discriminant',zero(disc-delta*res**r))
    discriminants.append({'r':r,'d':d,'N':N,'arithmetic_square_class':delta})
# Uniform cofactor, branch centering and coefficient extension bound.
marked=[]
for r in range(1,9):
 for d in [2,3,4]:
    I=integral(r,d).subs(A,1);value=I.subs(t,1)
    beta=S.factorial(r)*d**r/S.prod(d*j+1 for j in range(r+1))
    Q,rem=S.div(I-value,(t-1)**(r+1),t)
    label='rank'+str(r)+' h=t^'+str(d)
    check(label+' exact rational critical integral',zero(value-beta) and value!=0)
    check(label+' marked cofactor exact division',rem==0)
    check(label+' exact cofactor degree bound',S.degree(Q,t)==r*(d-1))
    check(label+' all cofactor roots regular',S.gcd(S.Poly(Q,t),S.Poly(1-t**d,t)).degree()==0)
    check(label+' distinct root-of-unity critical values',S.rem(I-value,t**d-1,t)==value*(t-1))
    marked.append({'r':r,'d':d,'degree':r*d+1,'coefficient_field_degree_upper':r*(d-1),'Q':str(Q)})
# Constructor identities in ranks not covered by the preceding turn.
higher=[]
for r in [5,6]:
    f=family(r,t*t,t);ps=f['parameters'];g=f['g'];X=f['X'];ys=f['target_variables'];swept=f['S']
    check('rank'+str(r)+' actual unit-triangular source',all(zero(S.diff(X[i],ps[i])-1) and
       not any(X[i].has(p) for p in ps[i+1:]) for i in range(1,r)))
    check('rank'+str(r)+' exact entire integral target identity',
      zero(swept[-1]-f['target_polynomial'].subs(dict(zip(ys,list(swept)[:r])),simultaneous=True)
                      -(-1)**(r-1)*f['integral'].subs(ys[0],swept[0])))
    check('rank'+str(r)+' derivative giving normalized determinant',
          zero(S.diff(f['integral'],t)-(ys[0]-t*t)**r))
    higher.append({'r':r,'X':[str(e) for e in X],'degree':f['degree'],
      'whole_Jacobian_check_method':'exact full target polynomial identity plus unit-triangular inversion and derivative; not a dense expanded determinant'})
# Exact elementary A5 exception, no group-table inference.
def multiply(a,b):return tuple(a[b[i]] for i in range(len(a)))
identity=tuple(range(5));a=(1,2,0,3,4);b=(0,1,3,4,2)
elements={identity};todo=[identity]
while todo:
    q=todo.pop()
    for gen in [a,b]:
        v=multiply(q,gen)
        if v not in elements:elements.add(v);todo.append(v)
check('small primitive degree5 two3-cycle group has order60',len(elements)==60)
check('small degree5 group transitive',len({q[0] for q in elements})==5)
# Regular nilpotency must not be confused with Nijenhuis integrability.
f=family(4,t*t,t);ps=f['parameters'];X=f['X']
P=S.Matrix.hstack(-S.Matrix(ps),S.eye(4))*X.jacobian(ps)
torsion=S.expand(sum(P[l,0]*S.diff(P[0,1],ps[l])-P[l,1]*S.diff(P[0,0],ps[l])
    -P[0,l]*(S.diff(P[l,1],ps[0])-S.diff(P[l,0],ps[1])) for l in range(4)))
check('actual nilpotent P is not Nijenhuis',zero(torsion-(4*t-2)) and torsion!=0)
check('actual P has regular nilpotent rank',P.rank()==3 and all(zero(e) for e in P**4))
# All degree11 factorizations, including the cyclic linear case.
degree11=[]
for r in range(1,11):
 if 10%r==0:
    d=10//r;group='C11' if d==1 else ('S11' if r%2 else 'A11')
    degree11.append({'r':r,'d':d,'group':group})
check('all integrated degree11 ranks classified',degree11==[
 {'r':1,'d':10,'group':'S11'},{'r':2,'d':5,'group':'A11'},
 {'r':5,'d':2,'group':'S11'},{'r':10,'d':1,'group':'C11'}])
# An actual nontrivial affine equivalence with gamma scaling verifies all blocks.
f=family(2,t*t,t);g=f['g'];ps=f['parameters'];xx=f['X']
lam=S.Rational(3);bb=S.Matrix([[2,1],[1,1]]);vv=S.Matrix([1,-2]);cc=S.Matrix([2,-1,4])
M=S.Matrix([[lam,0,0],[vv[0],bb[0,0],bb[0,1]],[vv[1],bb[1,0],bb[1,1]]])
q=S.Matrix(S.symbols('q1:3'));newp=(vv+bb*S.Matrix(ps))/lam
oldp=bb.inv()*(lam*q-vv)
newX=(M*xx+cc).subs(dict(zip(ps,oldp)),simultaneous=True)
snew=newX+g*S.Matrix([1,*q])
check('nontrivial affine sweep equivalence entire polynomial identity',
 all(zero(e) for e in M*f['S']+cc-snew.subs({g:lam*g,**dict(zip(q,newp))},simultaneous=True)))
check('nontrivial affine equivalence retains normalized Jacobian',
 zero(snew.jacobian([g,*q]).det()-g**2))
# A gamma-translation would move the exact critical divisor.
shift=S.Symbol('b')
check('nonzero gamma translation changes normalized critical polynomial',
 S.Poly((g+shift)**4-g**4,g).coeff_monomial(g**3)==4*shift)
manifest=json.loads((ROOT/'sources/all-rank-integrated-manifest.json').read_text())
for src in manifest['retained_primary_inputs']:
    raw=(ROOT/src['path']).read_bytes()
    check('retained primary bytes '+src['path'],len(raw)==src['bytes'])
    check('retained primary hash '+src['path'],hashlib.sha256(raw).hexdigest()==src['sha256'])
model={'claims':['KMFI112','KMFI113','KMFI114'],
'all_rank_geometric':'d1 cyclic C_(r+1); d>=2 S_(rd+1) for oddr, A_(rd+1) for evenr',
'even_rank_arithmetic_constant_square_class':'(-1)^(rd/2)*(rd+1)',
'linear_case_constants':'k(zeta_(r+1)); regular actions absolutely excluded',
'composite_cycle_primary_locator':'Jones Corollary1.3, not just prime-cycle Theorem1.1',
'discriminant_controls':discriminants,'marked_cofactor_controls':marked,
'higher_rank_constructor_controls':higher,'small_A5_group_order':len(elements),
'non_Nijenhuis_countercontrol':str(torsion),'degree11_complete_integrated_actions':degree11,
'geometric_completed_normalization_r_ge2':'Cl=Z^(rd), canonical vector(r,0,...,0)',
'rank1_class_formula_claimed':False,'coefficient_extension_degree_upper':'r*(d-1)',
'coordinate_rigidity_scope':'source and inverse carry normalized gamma-lines onto gamma-lines',
'arbitrary_coordinate_equivalences_classified':False,
'actual_PSL11_field_supplied':False,'absolute_exceptional_action_decided':False,'goal_complete':False}
(ROOT/'all_rank_integrated_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['ALL_RANK_INTEGRATED_ACTIONS.md','verify_all_rank_integrated.py','all_rank_integrated_model.json',
 'sources/all-rank-integrated-manifest.json','integrated_sweeps.py','INTEGRATED_SWEEP_PENCILS.md',
 'HIGHER_RANK_COMPLETION.md','COMPLETION_EXCEPTIONAL_DIVISOR.md']
paths += [e['path'] for e in manifest['retained_primary_inputs']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,'runtime_seconds':round(time.monotonic()-start,3),
 'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 'assurance':'Written all-rank critical-value, indecomposability, geometric/arithmetic group, exact coefficient-degree and coordinate-rigidity proofs. Finite exact controls supplement them. No formal, independent whole-proof or external review claimed.',
 'goal_complete':False}
(ROOT/'all_rank_integrated_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'runtime_seconds':receipt['runtime_seconds'],
 'degree11_actions':degree11,'actual_PSL11_field_supplied':False,'goal_complete':False}))
