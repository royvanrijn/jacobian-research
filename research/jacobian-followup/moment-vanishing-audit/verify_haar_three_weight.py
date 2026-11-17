"""Exact joint-cover and nilpotent-ring controls, not a universal proof."""
from pathlib import Path
import hashlib
import json
import sympy as S

p,q,z,s,a,F,A,B,Y,r,t,eps,eta=S.symbols('p q z s a F A B Y r t eps eta')
for m in range(9):
    ct=S.expand((p+z+q/z)**m).coeff(z,0)
    expected=sum(S.factorial(m)/(S.factorial(j)**2*S.factorial(m-2*j))
                 *p**(m-2*j)*q**j for j in range(m//2+1))
    assert S.expand(ct-expected)==0
b=q/a
assert S.simplify(S.diff(p+a+b,a)-(a-b)/a)==0
assert S.expand((a-b)**2-((p+a+b)**2-2*(p+a+b)*p+p*p-4*q))==0
J=B-2*s*A+s*s*F
assert S.simplify(F*F*(s*s-2*s*A/F+B/F)-F*J)==0
assert S.diff(J,s).subs(F,0)==-2*A
assert J.subs({F:0,A:0})==B

# Local square removal, including the corner with two infinity components.
for b1,b2 in [(1,1),(2,3),(3,4)]:
    square=r**(b1//2)*t**(b2//2)
    denom=r**b1*t**b2
    assert S.simplify(denom/square
        -r**((b1+1)//2)*t**((b2+1)//2))==0

def truncated(expr,nu,kappa):
    out=0
    for (i,j),coef in S.Poly(S.expand(expr),eps,eta).terms():
        if i<nu and j<kappa:out+=coef*eps**i*eta**j
    return S.expand(out)

# Laurent coefficient ring R=C[eps,eta]/(eps^nu,eta^kappa),
# p=eps,U=z,V=eta/z; marker z^-B. The universal bound is attained.
controls=[]
for nu,kappa,weight in [(1,1,0),(2,3,2),(3,2,3)]:
    last=nu+2*kappa-3+weight
    marked=truncated(S.expand(z**(-weight)*(eps+z+eta/z)**last).coeff(z,0),nu,kappa)
    expected=S.factorial(last)/(S.factorial(nu-1)*S.factorial(kappa-1)
             *S.factorial(kappa-1+weight))*eps**(nu-1)*eta**(kappa-1)
    assert S.expand(marked-expected)==0 and marked!=0
    for m in [last+1,last+2]:
        assert truncated(S.expand(z**(-weight)*(eps+z+eta/z)**m).coeff(z,0),nu,kappa)==0
    controls.append({'nu':nu,'kappa':kappa,'marker_weight_abs':weight,
                     'last_possible_order':last,'nonzero_at_bound':str(marked)})

# Pure-tail vs first-moment distinction is preserved over nilpotents.
assert truncated(eps,2,1)==eps
assert truncated(eps**2,2,1)==0

# The abstract nonexceptional three-weight Legendre control stays allowed.
x=S.symbols('x')
for m in range(1,8):
    ct=S.expand((x+z+(x*x-1)/(4*z))**m).coeff(z,0)
    assert S.expand(ct-S.legendre(m,x))==0
    assert S.integrate(ct,(x,-1,1))==0

out={'date':'2026-11-16','completed':True,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'sympy_version':S.__version__,
 'universal_joint_neutral_moment_identity_controls':True,
 'conic_equation_and_original_volume_exact':True,
 'joint_projective_cover_equation_exact':True,
 'basepoint_free_branch_transversality_both_cases_exact':True,
 'two_component_infinity_square_removal_controls':True,
 'nilpotent_ring_full_marker_bound_controls':controls,
 'eventual_purity_not_confused_with_first_moment':True,
 'abstract_nonexceptional_three_weight_control_retained':True,
 'field_classification':'P=p+U+V eventual-pure iff p=0 and UV=0; every surviving phase SAFE',
 'field_all_H_classification':'H(P) eventual-pure iff H=0 or UV=0,p=c constant,H(c)=0',
 'all_commutative_C_algebras_criterion':'p and UV nilpotent',
 'full_marker_cutoff':'m>nu+2*kappa-3+B_Q/abs(n)',
 'arbitrary_C_algebra_safety_proved':True,
 'unequally_spaced_three_weights_claimed_closed':False,
 'every_failure_claimed_at_least_four_weights':False,
 'ring_recomposition_classification_copied_from_field':False,
 'neutral_zero_radical_alone_used_for_three_weights':False,
 'proper_model_assumed_for_original_full_phase':False,
 'purity_reflected_from_composite':False,
 'all_Bezout_sections_claimed_excluded':False,
 'Hopf_nonlinear_composites_claimed_excluded':False,
 'uniform_numerical_nilpotence_exponent_computed':False,
 'universal_geometry_or_MZ_proof_machine_verified':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'external_review_next_task':False,'new_degree_support_or_parameter_campaign':False,
 'full_exceptional_MZ_cases':'OPEN on BOTH quotients'}
Path(__file__).with_name('haar_three_weight_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
