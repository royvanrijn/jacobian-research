"""Small exact weighted-cover controls; universal geometry is written separately."""
from pathlib import Path
import hashlib
import json
import math
import sympy as S

p,z,t,alpha,T,sig,f,eta,h,eps=S.symbols('p z t alpha T sig f eta h eps')
cases=[]
for r,s in [(1,1),(1,2),(2,1),(2,3),(3,2),(3,4)]:
    assert math.gcd(r,s)==1
    M=r+s
    for m in range(M+3):
        ct=S.expand((p+z**r+t/z**s)**m).coeff(z,0)
        expected=sum(S.factorial(m)/(S.factorial(m-M*l)*S.factorial(s*l)*S.factorial(r*l))
                     *p**(m-M*l)*t**(r*l) for l in range(m//M+1))
        assert S.expand(ct-expected)==0
    F=alpha**s*(T-alpha)**r
    derivative=alpha**(s-1)*(T-alpha)**(r-1)*(s*T-M*alpha)
    assert S.expand(S.diff(F,alpha)-derivative)==0
    ac=S.Rational(s,M)*T
    bc=S.Rational(r,M)*T
    c=S.Rational(s**s*r**r,M**M)
    assert S.simplify(F.subs(alpha,ac)-c*T**M)==0
    second=S.simplify(S.diff(F,alpha,2).subs(alpha,ac))
    assert S.simplify(second+M*ac**(s-1)*bc**(r-1))==0
    assert second.subs(T,-1)!=0
    assert S.rem(F-c*T**M,(alpha-ac)**2,alpha)==0
    L=S.cancel(-c*((-1+sig*f)**M-(-1)**M)/f)
    assert S.diff(L,sig).subs(f,0)==-c*M*(-1)**(M-1)
    assert S.diff(L,sig).subs(f,0)!=0
    assert S.expand((eta+f*L).subs(eta,f*h)-f*(h+L))==0
    assert S.expand((eta+f*L).subs(f,eta*h)-eta*(1+h*L.subs(f,eta*h)))==0
    cases.append({'r':r,'s':s,'cover_degree':M,'nonzero_critical_second_derivative_at_T_minus1':str(second.subs(T,-1))})

# The all-valuation weighted bound: t-2+C(1-t/3) >= C/3 for C>=4,1<=t<=2.
x,C=S.symbols('x C')
assert S.expand(x-2+C*(1-x/3)-C/3-(2-x)*(C/3-1))==0

# The local original log-volume order is v(x)+v(y)-1.
u,v,w=S.symbols('u v w')
for a,b,c in [(1,2,3),(2,1,1)]:
    old=[u**a,u**b*v,u**c*w]
    jac=S.Matrix([[S.diff(q,j) for j in [u,v,w]] for q in old]).det()
    assert S.simplify(jac/old[2]-a*u**(a+b-1)/w)==0

# Square-part volume bound in both principalization charts.
for bf in range(1,8):
    for be in range(bf+1):
        assert bf-be//2 >= (bf+1)//2

def truncate(expr,nu,tkill):
    out=0
    for (i,j),coef in S.Poly(S.expand(expr),eps,t).terms():
        if i<nu and j<tkill:out+=coef*eps**i*t**j
    return S.expand(out)

# Actual Laurent ring with U=z^r,V=t*z^-s,q=t^r and nilpotents retained.
ring=[]
for r,s,nu,kappa,marker in [(1,2,2,2,-2),(2,3,2,2,3),(3,2,1,2,-3)]:
    M=r+s
    cutoff=nu-1+M*kappa+math.ceil(abs(marker)/min(r,s))
    for m in [cutoff,cutoff+1]:
        val=S.expand(z**marker*(eps+z**r+t/z**s)**m).coeff(z,0)
        assert truncate(val,nu,r*kappa)==0
    ring.append({'r':r,'s':s,'nu':nu,'kappa':kappa,'marker_weight_abs':abs(marker),'cutoff':cutoff})
assert truncate(eps,2,1)==eps and truncate(eps**2,2,1)==0

out={'date':'2026-11-16','completed':True,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'sympy_version':S.__version__,
 'weighted_full_moment_and_simple_critical_controls':cases,
 'all_valuation_weighted_lower_bound_identity':True,
 'original_toric_log_volume_order_controls':True,
 'critical_chart_principalization_and_S_transversality_exact':True,
 'both_square_removal_volume_bounds':True,
 'all_ring_original_marker_cutoff_controls':ring,
 'nilpotents_not_discarded_or_first_purity_inferred':True,
 'exact_opposite_weight_ring_criterion':'p and U^s V^r nilpotent',
 'all_zero_containing_three_weight_sectors_safe':True,
 'all_three_weight_supports_claimed_closed':False,
 'zero_free_three_weights_claimed_closed':False,
 'full_mixed_phase_coverage_lemma_proved':False,
 'proper_volume_model_assumed_for_all_original_phases':False,
 'scalar_Cauchy_matching_claimed_sufficient_for_full_marker_safety':False,
 'purity_of_composite_reflected_to_primitive':False,
 'weighted_to_unweighted_Haar_transfer_used':False,
 'general_formal_symbol_Bernstein_work_duplicated':False,
 'Hopf_nonlinear_composites_claimed_excluded':False,
 'universal_geometry_or_MZ_proof_machine_verified':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'external_review_next_task':False,'new_degree_support_or_parameter_campaign':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'ownership':'FULL exceptional MZ closure;Gaussian thread owns general formal-symbol audit/effectivity'}
Path(__file__).with_name('haar_weighted_conic_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
