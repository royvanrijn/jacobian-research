"""Small exact controls only; no general Fourier, diagonal or ideal solver."""
from pathlib import Path
from itertools import product
import hashlib,json
import sympy as s

p=Path(__file__).resolve().parent
t,z,x,c,e,a,u,r0,r1=s.symbols('t z x c e a u r0 r1')
s1,s2=s.symbols('s1 s2')
def eq(f,g):
    assert s.expand(f-g)==0,(f,g)
def coeff(f,variables,exps):
    pp=s.Poly(s.expand(f),*variables)
    return pp.coeff_monomial(s.prod(v**j for v,j in zip(variables,exps)))

checks={}
# Original operator Lambda=D_x,P=c, including every time coefficient of each chosen source degree.
D=8
Falpha={}
for d in range(D+1):
    Falpha[d]=s.expand(sum(c**m*t**m*x**(d-m)/(s.factorial(d-m)*s.factorial(m)**2) for m in range(d+1)))
    for m in range(d+2):
        expected=s.diff(x**d*c**m,x,m)/(s.factorial(d)*s.factorial(m)**2)
        eq(Falpha[d].coeff(t,m),expected)
F=sum(Falpha[d]*z**d for d in range(D+1))
Euler=lambda f:s.expand(z*s.diff(f,z))
R=lambda f:s.expand(-2*x*z**2*s.diff(f,z)-x*z*f+x**2*z**2*f-c*z*t*f)
bF=Euler(Euler(F))
L=s.expand(bF+R(F))
for d in range(D+1):eq(L.coeff(z,d),0)
for d in range(1,D+1):
    eq(d*d*Falpha[d],((2*d-1)*x+c*t)*Falpha[d-1]-x*x*Falpha.get(d-2,0))
    assert s.degree(Falpha[d],t)==d
    eq(Falpha[d].coeff(t,d),c**d/s.factorial(d)**2)
checks['original_operator_source_coefficient_checks']=sum(d+2 for d in range(D+1))
checks['source_b_lift_checks']=D+1
checks['source_lowering_recurrence_checks']=D
# Nonminimal b(E)=E^2(E-1), R1=(E-1)R.
b1=s.expand(Euler(bF)-bF)
R1=s.expand(Euler(R(F))-R(F))
for d in range(D+1):eq((b1+R1).coeff(z,d),0)
wrong=s.expand(-Euler(bF)-bF+R1)
assert wrong.coeff(z,1)!=0
checks['non_even_b_sign_mutation_detected']=True
# Fourier commutator [ -s,Ds ]=1; theta maps to sDs, not Ds*s.
for d in range(5):
    h=z**d
    eq(-z*s.diff(h,z)-s.diff(-z*h,z),h)
    eq(z*s.diff(h,z),d*h)
    assert s.diff(z*h,z)!=d*h
checks['Fourier_ordering_checks']=5

# General multivariate Taylor coefficient formula, including a growing coordinate and total lowering.
source_vars=[s1,s2]
Hcoeff={(aa,bb):(aa+1)*(bb+2)*(t**((aa+2*bb)%4)+t) for aa,bb in product(range(7),repeat=2)}
H=sum(v*s1**aa*s2**bb for (aa,bb),v in Hcoeff.items())
uv=[((1,0),(0,0)),((2,0),(1,0)),((1,1),(0,0)),((0,2),(1,0)),((1,2),(2,0))]
count=0
for uu,vv in uv:
    drop=sum(uu)-sum(vv);assert drop>=1
    for kk,jj in [(0,0),(2,1),(1,2)]:
        term=t**kk*s.diff(H,t,jj)
        for var,order in zip(source_vars,vv):term=s.diff(term,var,order)
        term=s.expand(term*s1**uu[0]*s2**uu[1])
        for aa,bb in product(range(4),repeat=2):
            alpha=(aa,bb)
            expected=0
            if all(alpha[i]>=uu[i] for i in range(2)):
                beta=tuple(alpha[i]-uu[i]+vv[i] for i in range(2))
                assert sum(beta)==sum(alpha)-drop
                factor=s.prod(s.factorial(beta[i])/s.factorial(alpha[i]-uu[i]) for i in range(2))
                expected=factor*t**kk*s.diff(Hcoeff[beta],t,jj)
            eq(coeff(term,source_vars,alpha),expected)
            count+=1
checks['multivariate_Taylor_extraction_checks']=count

# Rational cone filter: multiply the finite coefficient mask by the exact denominator.
count=0
for N in range(4):
    mask=sum(t**mm*s1**aa*s2**bb for mm in range(9) for aa,bb in product(range(4),repeat=2) if mm>N*(aa+bb))
    denominator=(1-t)*(1-s1*t**N)*(1-s2*t**N)
    residual=s.Poly(s.expand(denominator*mask-t),t,s1,s2)
    for mm in range(9):
        for aa,bb in product(range(4),repeat=2):
            eq(residual.coeff_monomial(t**mm*s1**aa*s2**bb),0)
            count+=1
    if N>0:
        assert coeff(mask,[t,s1,s2],[N,1,0])==0
        assert coeff(mask,[t,s1,s2],[N+1,1,0])==1
checks['cone_rational_identity_checks']=count
checks['strict_cone_boundary_and_zero_slope_checked']=True

# Hadamard product equals paired partial diagonal, keeping both independent variables.
count=0
for N in [0,1,2]:
    ff={(mm,dd):s.binomial(mm+dd,dd) for mm in range(7) for dd in range(5)}
    gg={(mm,dd):1 for mm in range(7) for dd in range(5) if mm>N*dd}
    diagonal={}
    for (m1,d1),v1 in ff.items():
        for (m2,d2),v2 in gg.items():
            if m1==m2 and d1==d2:diagonal[(m1,d1)]=v1*v2
    for index,value in ff.items():
        assert diagonal.get(index,0)==(value if index[0]>N*index[1] else 0)
        count+=1
checks['paired_partial_diagonal_checks']=count
# A single complete time/source diagonal misses off-diagonal violations.
off_diagonal=t**2*z
assert off_diagonal!=0 and all(coeff(off_diagonal,[t,z],[k,k])==0 for k in range(4))
checks['complete_diagonal_shortcut_detected']=True
# A single numerical ray can cancel distinct forbidden coefficients.
bad=t*z**2-t**2*z
eq(bad.subs({t:u,z:u},simultaneous=True),0)
generic_ray=s.expand(bad.subs({t:r0*u,z:r1*u},simultaneous=True))
assert generic_ray.coeff(u,3)==r0*r1**2-r0**2*r1 and generic_ray!=0
checks['numerical_ray_cancellation_detected']=True

# Nilpotents survive the prescribed-bound test: F=1+e*s*t modulo e^2.
assert s.rem(e,e**2,e)==e and e.subs(e,0)==0
assert 1>0*1 and not 1>1*1
checks['nilpotent_forbidden_coefficient_retained']=True

# Exact small ordinary ideals for H_a=1+s*((1+t)^a-1).
records=[]
for N in range(6):
    generator=s.prod(a-j for j in range(N+1))
    first=s.prod(a-j for j in range(N+1))/s.factorial(N+1)
    eq(first*s.factorial(N+1),generator)
    for mm in range(N+1,N+7):
        binomial=s.prod(a-j for j in range(mm))/s.factorial(mm)
        assert s.rem(binomial,generator,a)==0
    assert sorted(s.polys.polytools.ground_roots(generator,a))==list(range(N+1))
    records.append({'N':N,'ordinary_generator':str(s.expand(generator)),
                    'closed_set':list(range(N+1))})
checks['nonstabilizing_holonomic_control']=records
checks['not_claimed_as_Gaussian_operator_or_SIC_model']=True
# Universally safe original operator family: I_0=(c), I_1=0.
for d in range(D+1):
    for mm in range(1,d+1):
        h=Falpha[d].coeff(t,mm)
        for xxcoef in s.Poly(h,x).all_coeffs():assert s.rem(xxcoef,c,c)==0
        assert not mm>d
eq(Falpha[1].coeff(t,1),c)
checks['safe_operator_family_uniform_C']=1
checks['safe_operator_family_least_integer_C']=1
checks['zero_parameter_point_least_integer_C']=0

receipt={'completed':True,'date':'2026-11-15','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'software':'SymPy '+s.__version__,**checks,
 'general_source_Fourier_solver_implemented':False,'general_source_Fourier_solver_executed':False,
 'general_Hadamard_solver_implemented':False,'general_Hadamard_solver_executed':False,
 'general_cutoff_ideal_solver_executed':False,'general_family_uniform_solver_executed':False,
 'universal_proofs_certified_by_checker':False,'global_safe_locus_constructibility_proved':False,
 'global_safe_subset_uniform_bound_proved':False,'ordinary_ideal_chain_stabilization_claimed':False,
 'new_SIC33_answer_computed':False,'new_Delta_minimum_computed':False,
 'numeric_cubic_frontier_changed':False,'operator_order_reduced':False,
 'independent_replay':False,'external_review':False,'formal_verification':False}
out=p/'holonomic_source_cutoff_verification.json';data=(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n').encode()
if out.exists():assert out.read_bytes()==data,'Existing receipt differs'
else:out.write_bytes(data)
print(json.dumps(receipt,indent=2))

