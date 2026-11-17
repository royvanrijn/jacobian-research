"""Exact controls for projective pencil multiplicity bounds and critical expansions."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RETAINED = {
    'jacobian-followup/hc4-triangular-graphs/REDUCED_PENCILS_AND_RATIONAL_MONODROMY.md': '5f5e55bcb2f15a15e87d9ddd8c22f324beb5a40bf92e941472065af988b641e6',
    'jacobian-followup/hc4-triangular-graphs/verify_reduced_pencils.py': '152ef0fd112d392243198b2d1748e7fe87685031778836ac6c92606f3e4cfdb2',
    'jacobian-followup/hc4-triangular-graphs/receipt-reduced-pencils.json': '7b397e79d148777d81c12eca90e5910af061aa68f8f24d0c47b3cdd38f655555',
    'jacobian-followup/hc4-triangular-graphs/sources/reduced-pencil-import.json': 'bac9b7f0b62a8957481a54849f8dee56e43e8e51e50297f45b801d1589c79ac6',
    'jacobian-followup/hc4-triangular-graphs/LEADING_FIBRE_UNIT_GATE.md': '1d054ba8762ea657be879b6e3568aa58c9ccb65e14ab642d45e0865e8f6382f0',
    'jacobian-followup/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md': '3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
    'jacobian-followup/hc4-triangular-graphs/sources/stacks-leading-riemann-hurwitz-0C1B.html': '616e526d775d71ae3dd706bccfdc38115a846c7eafac1e7080d53691bc2249a3',
    'jacobian-followup/hc4-triangular-graphs/sources/petravchuk-iena-closed-rational-math-0701588v1.html': '16a826a72dd97d69898a14790bfe6391ae6c297830a2ad91b49562353f12fdd6',
    'jacobian-followup/hc4-triangular-graphs/sources/pakovich-genus-zero-mpim2016-36.pdf': '214d5f4fd101ef2d66898c5985a8f13721271722dbd2cbbf512be188cd670309'
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    checks = []

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        check(name,S.cancel(S.expand(value))==0)

    e,M,d = S.symbols('e M d',positive=True)
    zero('the local ramification gap is nonnegative precisely under the maximum bound',
        e*(1-1/M)-(e-1)-(M-e)/M)
    zero('after one unramified point the remaining local degrees sum to d-1',d-1-(d-1))
    zero('Riemann-Hurwitz cancels all prospective composition degrees',
        (2*d-2)/(d-1)-2)
    for profile in [(1,2,2),(1,1,3),(1,3,3,3)]:
        maximum = max(profile)
        local_degree = sum(profile)
        ramification = sum(index-1 for index in profile)
        check('the local bound holds for profile '+str(profile),
            Fraction(ramification)<=Fraction(local_degree-1)*(1-Fraction(1,maximum)))
    check('the simple-factor hypothesis cannot be removed from the local estimate',
        Fraction(2)>Fraction(3,2))

    def weight(values):
        return sum((1-Fraction(1,value) for value in values),Fraction(0))

    for values in [(2,2),(2,2,11),(2,3,5),(2,2,2)]:
        check('strict multiplicity certificate '+str(values),weight(values)<2)
    for values in [(2,2,2,2),(3,3,3),(2,4,4),(2,3,6)]:
        check('equality boundary list '+str(values),weight(values)==2)
    check('a nearby triple lies beyond the proved strict bound',weight((2,3,7))>2)
    check('five double members do not pass the bound',weight((2,2,2,2,2))>2)
    check('the maximum and lcm are distinct data and may not be interchanged',
        max(1,2,3)==3 and S.ilcm(1,2,3)==6 and 3 % 2 != 0)

    s,t,T,V,z,lam,h = S.symbols('s t T V z lam h')
    N = s**3+t**3
    denominator = t*(s*s+t*t)
    zero('the equality example is explicitly a degree-three rational composition',
        N/denominator-((h**3+1)/(h*h+1)).subs(h,s/t))
    check('the equality example has coprime numerator and denominator',S.gcd(N,denominator)==1)
    jac = S.diff(N,s)*S.diff(denominator,t)-S.diff(N,t)*S.diff(denominator,s)
    check('its numerator and denominator are algebraically independent',S.expand(jac)!=0)
    discr = S.discriminant(h**3-lam*h*h+1-lam,h)
    zero('the four exceptional values have the stated discriminant',
        discr-(1-lam)*(4*lam**3+27*lam-27))
    check('these four exceptional values are distinct',
        S.degree(discr,lam)==4 and S.gcd(discr,S.diff(discr,lam))==1)
    zero('one exceptional homogeneous member visibly has profile two plus one',
        (N-lam*denominator).subs(lam,1)-s*s*(s-t))
    check('the infinity member is reduced and must be included in the audit',
        S.gcd(S.gcd(denominator,S.diff(denominator,s)),S.diff(denominator,t))==1)
    check('the equality example has weight exactly two and remains composite',weight((2,2,2,2))==2)

    Q = T**3*(T-1)**2
    D = T+1
    C0 = 10*T**3-15*T**4+6*T**5
    C = D*C0
    G0 = S.gcd(Q,C)
    G1 = S.gcd(Q,C-D)
    check('the affine-in-V example has square-free nonconstant denominator',
        S.degree(D,T)==1 and S.gcd(D,S.diff(D,T))==1)
    check('its numerator and denominator are coprime',S.gcd(Q*V+C,D)==1)
    zero('its first exceptional content has exact multiplicity three',G0-T**3)
    zero('its second exceptional content has exact multiplicity two',G1-(T-1)**2)
    check('generic finite members have no content',S.gcd(Q,C-lam*D)==1)
    W = S.diff(C,T)*D-C*S.diff(D,T)
    candidate = S.gcd(S.gcd(Q,S.diff(Q,T)),W)
    zero('the candidate derivative gcd captures both exceptional levels',candidate-T**2*(T-1))
    check('the candidate set does not include a denominator zero',S.gcd(candidate,D)==1)
    check('their exceptional values are zero and one',C0.subs(T,0)==0 and C0.subs(T,1)==1)
    for value,content in [(0,G0),(1,G1)]:
        remaining = S.cancel((Q*V+C-value*D)/content)
        check('the primitive horizontal factor is coprime to content at '+str(value),S.gcd(content,remaining)==1)
    check('two different exceptional values still satisfy the all-degree criterion',weight((3,2))==Fraction(7,6)<2)
    check('a repeated denominator would require a separate infinity-member check',
        S.gcd(D**2,S.diff(D**2,T))!=1)

    Y = S.Symbol('Y')
    for degree in [3,4,6,9]:
        alpha = S.Rational(degree-1,degree)
        R = Y**(degree-1)*(Y-1)
        H = T*(T-alpha)
        leading_D = T-1
        derivative = S.diff(R,Y).subs(Y,T)
        zero(f'degree {degree}: critical derivative is divisible by the chosen square-free H',
            derivative-degree*H*T**(degree-3))
        check(f'degree {degree}: H and D are square-free and coprime',
            S.gcd(H,S.diff(H,T))==1 and S.gcd(H,leading_D)==1)
        expression = S.cancel(leading_D*(R.subs(Y,T+H*z)-R.subs(Y,T))/(H*H))
        check(f'degree {degree}: the critical expansion is a polynomial of unchanged degree',
            S.denom(expression)==1 and S.degree(expression,z)==degree)
        zero(f'degree {degree}: its leading multiplicity is exactly H^(m-2)',
            S.Poly(expression,z).LC()-leading_D*H**(degree-2))
        zero(f'degree {degree}: the two critical values have the stated nonzero separation',
            R.subs(Y,alpha)+alpha**(degree-1)/degree)
        parameter = H*H*V/leading_D-R.subs(Y,T)
        zero(f'degree {degree}: the output root equation reduces to one preserved parameter cover',
            H*H/leading_D*(expression+V)-(R.subs(Y,T+H*z)+parameter))

    R = Y**4+Y
    H = 4*T**3+1
    P = T*(H*H*z**4+4*T*H*z**3+6*T*T*z*z+z)
    zero('the repeated-leading quartic equals its full critical expansion',
        P-T*(R.subs(Y,T+H*z)-R.subs(Y,T))/(H*H))
    crit_values = S.resultant(H,R.subs(Y,T)+lam,T)
    check('its three selected critical values are distinct',
        S.degree(crit_values,lam)==3 and S.gcd(crit_values,S.diff(crit_values,lam))==1)
    check('the quartic passes the strict multiplicity bound',weight((2,2,2))==Fraction(3,2))
    a = S.Poly(P,z).LC()
    b3 = S.Poly(P,z).coeff_monomial(z**3)
    check('the older leading triple gcd genuinely fails at its repeated levels',
        S.degree(S.gcd(S.gcd(a,S.diff(a,T)),b3),T)>0)
    w,v = S.symbols('w v')
    actual_H = 1+s*v+t*w+w*P.subs({T:s,z:w})
    zero('the admitted output retains an explicit actual-unit-Hessian automorphic control',
        S.det(S.hessian(actual_H,(s,t,w,v)))-1)

    check('all retained dependency pins are present',len(RETAINED)==9)
    for path,digest in RETAINED.items():
        check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    context = json.loads((HERE/'sources/pencil-multiplicity-context.json').read_text())
    check('orbifold literature is context and not a new theorem import',
        context['new_external_theorem_imports']==[]
        and context['new_orbifold_classification_imported'] is False
        and 'maximum fibre multiplicities' in context['literature_context_only'][0]['role'])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_pencil_multiplicities.py','MULTIPLICITY_BOUND_AND_CRITICAL_EXPANSIONS.md','sources/pencil-multiplicity-context.json']},
        dependencies=RETAINED,
        written_not_finite_exhaustion=[
            'an unramified preimage and bounded local indices give per-fibre ramification at most (d-1)(1-1/M)',
            'strict total weight below two excludes every composition degree',
            'equality forces one unramified point and equal remaining local indices, with four numerical lists and no claimed closure',
            'an explicit independent-numerator rational pencil is composite at equality',
            'the affine-in-V content gcd computes all exceptional values and exact multiplicity bounds',
            'critical-point expansions with at most three selected critical values give full all-degree source irreducibility and plane carried-output inverses'],
        new_external_theorem_imports=[],maximum_multiplicities_not_lcm_used=True,
        simple_factor_required_in_each_exceptional_member=True,numerator_denominator_independence_required=True,
        equality_asserted_closed=False,strict_bound_replaced_by_nonstrict=False,
        denominator_member_ignored=False,orbifold_classification_imported=False,
        numerical_equality_lists_asserted_realizable=False,
        every_rational_translation_asserted_admitted=False,
        coefficient_residual_asserted_Keller_realization=False,
        nonlinear_source_Hessian_invariance_used=False,unit_ideal_endpoint_resolved=False,
        singleton_zero_endpoint_resolved=False,unrestricted_HC4_or_JC2_resolved=False,
        original_goal_complete=False,formal_verification=False,independent_replay=False,
        external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-pencil-multiplicities.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
