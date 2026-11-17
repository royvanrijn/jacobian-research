"""Selected exact orbit/weight controls; universal proofs are written separately."""
from pathlib import Path
from hashlib import sha256
from math import gcd, ceil
import json
import sympy as s

HERE=Path(__file__).resolve().parent
t,X,Y=s.symbols('t X Y')
smith=[]
for aa in (1,2,5):
    U=s.Matrix([[1,t],[0,1]])
    V=s.Matrix([[1,0],[t**2,1]])
    D=s.diag(t**aa,t**-aa)
    g=U*D*V
    assert s.simplify(g.det())==1
    assert s.simplify(U.inv()*g*V.inv()-D)==s.zeros(2)
    smith.append({'a':aa,'determinant_one':True,'regular_outer_factors':True})
moving=s.expand((X+t*Y)**2).subs({X:t*X,Y:Y/t},simultaneous=True).expand()
assert moving==t**2*X**2+2*t*X*Y+Y**2
assert moving.subs(t,0)==Y**2
assert (t**2*X**2).subs(t,0)==0
negative=s.expand((X**2+Y**2).subs({X:t*X,Y:Y/t},simultaneous=True))
assert s.expand(negative).coeff(Y,2)==t**-2


def ct(expr, variables):
    expr=s.expand(expr)
    out=0
    for term in s.Add.make_args(expr):
        powers=term.as_powers_dict()
        if all(powers.get(v,0)==0 for v in variables):
            out+=term
    return s.expand(out)


def multinomial(n, counts):
    return s.factorial(n)/s.prod(s.factorial(c) for c in counts)


weights_controls=[]
for weights, rel in [([(2,),(-4,)],(2,1)), ([(1,0),(0,1),(-1,-1)],(1,1,1)),
                     ([(0,0),(1,0),(0,1)],(1,0,0))]:
    rank=len(weights[0]);ss=len(weights)
    W=s.Matrix(rank,ss,lambda i,j:weights[j][i])
    augmented=W.col_join(s.ones(1,ss))
    assert augmented.rank()==ss and W*s.Matrix(rel)==s.zeros(rank,1)
    assert gcd(*rel)==1
    Z=s.symbols('Z0:'+str(rank));C=s.symbols('C0:'+str(ss))
    P=sum(C[j]*s.prod(Z[i]**weights[j][i] for i in range(rank)) for j in range(ss))
    M=sum(rel);B=s.prod(C[j]**rel[j] for j in range(ss))
    for m in range(1,7):
        expected=0 if m%M else multinomial(m,[a*(m//M) for a in rel])*B**(m//M)
        assert s.expand(ct(P**m,Z)-expected)==0
    weights_controls.append({'weights':weights,'primitive_relation':rel,'M':M,
                             'zero_barycentric_face':any(a==0 for a in rel)})

# A genuine two-dimensional shifted ray, with n0=(3,0,1).
Z0,Z1=s.symbols('Z0 Z1');C0,C1,C2=s.symbols('C0 C1 C2')
P=C0*Z0+C1*Z1+C2/(Z0*Z1)
for m in range(1,8):
    expected=0
    if m>=4 and (m-4)%3==0:
        ell=(m-4)//3
        expected=multinomial(m,[3+ell,ell,1+ell])*C0**(3+ell)*C1**ell*C2**(1+ell)
    assert s.expand(ct(Z0**-2*Z1*P**m,[Z0,Z1])-expected)==0

# Literal ordinary pure ideal can remain nonreduced.
A,B=s.symbols('A B')
ideal=s.groebner([A**2*B],A,B)
assert ideal.reduce(A**2*B)[1]==0
assert ideal.reduce(A*B)[1]==A*B
assert ideal.reduce((A*B)**2)[1]==0

x,y,z,w,c0,c1=s.symbols('x y z w c0 c1')
fplus=x**8+14*x**4*y**4+y**8
fminus=z**8+14*z**4*w**4+w**8


def haar(poly):
    ans=0
    for (a,b,c,d),coeff in s.Poly(s.expand(poly),x,y,z,w).terms():
        if c==b and d==a:
            ans+=coeff*(-1)**b*s.factorial(a)*s.factorial(b)/s.factorial(a+b+1)
    return s.expand(ans)


P=c0*fplus+c1*fminus;R=c0*c1*fplus*fminus
for m in range(1,5):
    expected=0 if m%2 else s.binomial(m,m//2)*haar(R**(m//2))
    assert s.expand(haar(P**m)-expected)==0
    expected_marked=0 if m%2==0 else s.binomial(m,(m+1)//2)*haar(c0*fplus*fminus*R**((m-1)//2))
    assert s.expand(haar(fminus*P**m)-expected_marked)==0
assert haar(fplus*fminus)==s.Rational(8,15)
ee,ff=s.symbols('ee ff');zero_product=s.groebner([ee*ff],ee,ff,domain=s.QQ)
for m in range(1,5):
    assert zero_product.reduce(haar((ee*fplus+ff*fminus)**m))[1]==0
assert zero_product.reduce(ee*haar(fplus*fminus))[1]==8*ee/15

nil_controls=[]
for aa,cc,nu,degree in [(2,4,2,5),(8,12,3,6),(12,20,2,0)]:
    gg=gcd(aa,cc);pp=cc//gg;qq=aa//gg;M=pp+qq
    cutoff=M*nu+ceil(degree/min(aa,cc))
    for m in range(cutoff,cutoff+3):
        for i in range(m+1):
            j=m-i;tau=cc*j-aa*i
            if abs(tau)<=degree:assert i>=pp*nu and j>=qq*nu
    nil_controls.append({'positive_weight':aa,'negative_magnitude':cc,'nilpotence_assumed':nu,
                         'marker_degree':degree,'sufficient_m_at_least':cutoff})

U,T=s.symbols('U T')
F=(1-1/T)*((1-U)+U*T)
interval=[]
for m in range(1,7):
    pure=s.integrate(ct(F**m,[T]),(U,0,1))
    marked=s.integrate(ct(F**m/T,[T]),(U,0,1))
    assert pure==0 and marked==s.Rational((-1)**(m+1),m+1)
    interval.append({'m':m,'pure':str(pure),'marked':str(marked)})

result={
 'completed':True,'sympy_version':s.__version__,
 'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
 'Smith_matrix_controls':smith,
 'moving_right_factor_arc_limit_differs_from_constant_factor_cocharacter_limit_exact':True,
 'negative_constant_weight_gives_unavoidable_pole_exact':True,
 'affinely_independent_pure_weight_controls':weights_controls,
 'two_dimensional_marker_shift_ray_exact':True,
 'ordinary_balanced_pure_ideal_can_be_nonreduced_exact':True,
 'actual_full_W8_opposite_charge_Haar_pure_and_marker_identities_exact':True,
 'actual_balanced_Haar_first_moment':'8/15',
 'zero_product_coefficient_ring_both_charge_blocks_nonzero_pure_control_exact':True,
 'nilpotent_balanced_product_cutoff_controls':nil_controls,
 'three_weight_interval_countercontrol':interval,
 'full_nonclosed_eventual_pure_sector':'exactly left-nullcone;ALL original markers safe with2m>degreeQ',
 'general_lifting_scope':'ALL affinely independent torus supports,ALL torus ranks,ALL original markers,ordinary pure ideals and Q-algebras,assuming MZ on weight-zero algebra',
 'full_two_weight_sector':'ALL degrees,ALL isotypic mixtures,ALL right multiplicities,ANY commutativeC-algebra',
 'full_failure_necessary_scope':'closed finite-cyclic left orbit afterEVERY algebraic automorphism;at leastTHREE weights for every maximallefttorus;prior mixed,image-interior,strict-leading gates retained',
 'automorphism_nullcone_saturation_proved_for_every_pure_phase':False,
 'closed_orbit_or_two_sided_support_implies_failure':False,
 'six_variable_prime_marker_imported_into_actual_exceptional_algebra':False,
 'zero_weight_MZ_sufficient_for_arbitrary_weight_supports':False,
 'universal_curve_Smith_orbit_and_weight_lifting_proofs_machine_verified':False,
 'new_degree_support_density_or_prefix_campaign':False,
 'external_review_next_task':False,
 'main_finish_line':'COMPLETE exceptionalMZ theorem;currently OPEN on BOTH quotients',
 'independent_replay':False,'formal_verification':False,'external_review':False}
(HERE/'orbit_weight_lifting_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
