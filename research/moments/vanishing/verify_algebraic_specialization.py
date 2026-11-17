"""Exact small controls; not an implementation of the general family algorithm."""
from pathlib import Path
import hashlib, json
from itertools import product
import sympy as s

p=Path(__file__).resolve().parent
checks={}
def eq(a,b):
    assert s.expand(a-b)==0,(a,b)

def contraction(h,xis,us,xs):
    total=0
    for exps,c in s.Poly(s.expand(h),*(xis+us)).terms():
        aa,bb=exps[:len(xis)],exps[len(xis):]
        if any(a>b for a,b in zip(aa,bb)):continue
        v=c
        for a,b,x in zip(aa,bb,xs):
            v*=s.factorial(b)/s.factorial(b-a)*x**(b-a)
        total+=v
    return s.expand(total)

def wick(h,zs,ws):
    out=0
    for exps,c in s.Poly(s.expand(h),*(zs+ws)).terms():
        aa,bb=exps[:len(zs)],exps[len(zs):]
        if aa==bb:
            v=c
            for a in aa:v*=s.factorial(a)
            out+=v
    return s.expand(out)

# Polynomial-valued Ecal versus independent Gaussian translation.
xi1,xi2,U1,U2,x1,x2,Z1,Z2,W1,W2=s.symbols(
    'xi1 xi2 U1 U2 x1 x2 Z1 Z2 W1 W2')
families=[
 ([xi1],[U1],[x1],[Z1],[W1],
  xi1*U1+xi1**2+U1**3,[s.Integer(1),xi1+U1,xi1*U1**2]),
 ([xi1,xi2],[U1,U2],[x1,x2],[Z1,Z2],[W1,W2],
  xi1*U1+xi2*U2+xi1*U2+U1*U2,[s.Integer(1),xi1+U2,xi2*U1]),
 ([xi1,xi2],[U1,U2],[x1,x2],[Z1,Z2],[W1,W2],
  xi1*xi2+U1*U2+xi1*U2,[s.Integer(1),xi1*U2,xi2+U1])]
count=0
for xis,us,xs,zs,ws,f,gs in families:
    shift=dict(zip(xis,zs));shift.update({u:x+w for u,x,w in zip(us,xs,ws)})
    for g in gs:
        for m in range(4):
            h=s.expand(g*f**m)
            eq(contraction(h,xis,us,xs),wick(h.subs(shift,simultaneous=True),zs,ws))
            count+=1
checks['full_nonseparable_contraction_checks']=count
# Explicit coefficient-matrix rank excludes rank-one separability in the first case.
eq(s.Matrix([[0,0,0,1],[0,1,0,0],[1,0,0,0]]).rank(),3)
checks['nonseparable_one_pair_example_rank']=3

# Reverse translation has no degree increase.
count=0
for a in range(3):
    for b in range(4):
        h=sum(s.binomial(b,j)*(-x1)**(b-j)*Z1**a*(x1+W1)**j for j in range(b+1))
        eq(h,Z1**a*W1**b)
        assert all(a+j<=a+b for j in range(b+1))
        count+=1
checks['reverse_translation_checks']=count
for m in range(1,5):
    eq(wick(W1**m,[Z1],[W1]),0)
    eq(contraction(U1**m,[xi1],[U1],[x1]),x1**m)
checks['origin_specialization_shortcut_detected']=True

# Safe original SIC f=xi^2 U: a fixed g=xi^a U^b dies for m>b-a.
count=0
for a,b in product(range(4),repeat=2):
    for m in range(7):
        got=contraction(xi1**a*U1**b*(xi1**2*U1)**m,[xi1],[U1],[x1])
        expected=0 if a+2*m>b+m else s.factorial(b+m)/s.factorial(b-a-m)*x1**(b-a-m)
        eq(got,expected)
        if m>a+b:assert got==0
        count+=1
checks['safe_original_linear_cutoff_checks']=count

# Bounded rational polynomial spaces inject at an Eisenstein point.
T=s.symbols('T')
records=[]
for d,M in [(1,1),(1,3),(2,1),(2,2),(3,1)]:
    base=max(M+1,2);N=base**d
    exponents=[sum(e[i]*base**i for i in range(d)) for e in product(range(M+1),repeat=d)]
    assert len(set(exponents))==(M+1)**d and min(exponents)==0 and max(exponents)<N
    polynomial=s.Poly(T**N-2,T)
    assert polynomial.LC()%2!=0
    assert all(c%2==0 for c in polynomial.all_coeffs()[1:])
    assert polynomial.TC()%4!=0
    # The monomial images have independent coordinate columns in 1,T,...,T^(N-1).
    columns=s.zeros(N,len(exponents))
    for col,e in enumerate(exponents):columns[e,col]=1
    assert columns.rank()==len(exponents)
    records.append({'dimension':d,'individual_degree_bound':M,'base':base,
                    'number_field_degree':N,'preserved_space_dimension':len(exponents)})
checks['bounded_space_injection_controls']=records
a1,a2,z,x=s.symbols('a1 a2 z x')
eq((a2-a1**3).subs({a1:T,a2:T**3},simultaneous=True),0)
checks['degree_box_boundary_shortcut_detected']=True

# Infinite integer conditions are a uniformly bounded space, not a finite root sample.
# Norm_{Q(a)(y)/Q(a)}, y^2=a, of z-y is z^2-a.
a,y=s.symbols('a y')
matrix=s.Matrix([[z,-a],[-1,z]])
eq(matrix.det(),z**2-a)
eq(s.resultant(y**2-a,z-y,y),z**2-a)
# Specialize a=sqrt(2), y=2^(1/4): integer roots of z-y would satisfy z^4-2=0.
norm=s.resultant(T**2-2,z**2-T,T)
eq(norm,z**4-2)
assert s.polys.polytools.ground_roots(norm,z)=={}
checks['norm_extension_control']={'minimal_polynomial':'y^4-2','integer_roots':[],
                                 'all_conjugates_nonintegral':True}
# Denominator and one nonzero late coefficient must be preserved separately.
eq(s.Matrix([[a,0],[0,a]]).det(),a**2)
assert s.rem(T**2,T**2-2,T)==2
checks['finite_open_and_coefficient_norm_controls']=True

# The binomial family shows why a finite nonzero jet does not preserve nonpolynomiality.
t=s.symbols('t')
for n in [0,1,5,9]:
    f=(1+t)**n
    eq((1+t)*s.diff(f,t)-n*f,0)
    assert s.degree(f,t)==n
# L=(1+t)D-y has infinity polynomial z-y; a=y^2, a=sqrt(2) preserves no integer roots.
fjet=sum(s.prod(y-j for j in range(k))/s.factorial(k)*t**k for k in range(7))
ode=s.expand((1+t)*s.diff(fjet,t)-y*fjet)
for k in range(6):eq(ode.coeff(t,k),0)
eq(fjet.coeff(t,0),1)
checks['binomial_integer_specialization_shortcut_detected']=True
checks['binomial_nonpolynomial_algebraic_specialization']={
 'annihilator':'(1+t)D-y','infinity_indicial':'z-y',
 'generic_bound':-1,'preserved_coefficient_index':0,
 'specialized_minimal_polynomial':'y^4-2'}

# Tiny pure component and safe specialization control. No primary solver is invoked.
# For P=c Z, every pure positive moment vanishes and
# F_(Z^a W^b)=c^(b-a)t^(b-a)*b!/(b-a)! for b>=a.
c=s.symbols('c')
count=0
for aa,bb in product(range(4),repeat=2):
    expected=s.Integer(0) if bb<aa else c**(bb-aa)*t**(bb-aa)*s.factorial(bb)/s.factorial(bb-aa)
    got=sum(wick(Z1**aa*W1**bb*(c*Z1)**m,[Z1],[W1])*t**m/s.factorial(m) for m in range(7))
    eq(got,expected)
    for cc in [0,1,-2]:
        eq(got.subs(c,cc),expected.subs(c,cc))
        count+=1
checks['generic_safe_tail_specialization_controls']=count

receipt={'completed':True,'date':'2026-11-15','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'software':'SymPy '+s.__version__,**checks,
         'general_family_decision_implemented':False,'general_family_decision_executed':False,
         'general_Noether_normalization_executed':False,'general_primary_ideal_solver_executed':False,
         'universal_proofs_certified_by_checker':False,'new_SIC33_answer_computed':False,
         'new_Delta_minimum_computed':False,'numeric_cubic_frontier_changed':False,
         'safe_locus_constructibility_proved':False,'uniform_family_cutoff_proved':False,
         'independent_replay':False,'external_review':False,'formal_verification':False}
out=p/'algebraic_specialization_verification.json'
content=(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n').encode()
if out.exists():assert out.read_bytes()==content,'Existing receipt differs'
else:out.write_bytes(content)
print(json.dumps(receipt,indent=2))

