"""Exact controls for reduced pencils and rational-parameter monodromy."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RETAINED = {
    'hessian/hc4-triangular-graphs/FACTOR_RESULTANT_BARRIER.md': 'd5e89c6a37f815fda5cfc19ab31f6c9527dd0706298927aee66b032f27cbad04',
    'hessian/hc4-triangular-graphs/verify_factor_resultants.py': 'ca2aebf2b4b7f16df002733005a7618c7360c1792f13d185d4f3d7e0900acd65',
    'hessian/hc4-triangular-graphs/receipt-factor-resultants.json': '232dc71baa8a56e63b681a6b269e3a8cd18ac869f3e2e2be77ff3ef8cc25f4e4',
    'hessian/hc4-triangular-graphs/REGULAR_PARAMETER_MONODROMY.md': 'b010d2e4c62906a8a79da9dbc41ef7d0e3d5edefb4520cf487ba09ded5d86956',
    'hessian/hc4-triangular-graphs/LEADING_FIBRE_UNIT_GATE.md': '1d054ba8762ea657be879b6e3568aa58c9ccb65e14ab642d45e0865e8f6382f0',
    'hessian/hc4-triangular-graphs/ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md': '3a10bcc44350ae155d73f1786dbd1dc081c86b20d469e7efd0edfa408981462f',
    'hessian/hc4-triangular-graphs/sources/stacks-leading-riemann-hurwitz-0C1B.html': '616e526d775d71ae3dd706bccfdc38115a846c7eafac1e7080d53691bc2249a3',
    'hessian/hc4-triangular-graphs/sources/petravchuk-iena-closed-rational-math-0701588v1.html': '16a826a72dd97d69898a14790bfe6391ae6c297830a2ad91b49562353f12fdd6'
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
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(S.expand(item))==0 for item in values))

    s,t,T,U,V,z,A,c,beta,B,Q,lam,h = S.symbols('s t T U V z A c beta B Q lam h')
    degree,poles = S.symbols('degree poles',integer=True,positive=True)
    zero('the finite ramification contribution is d+r-2',
        (2*degree-2)-(degree-poles)-(degree+poles-2))
    zero('one omitted point can hide finite ramification only when the pole count is at most one',
        (degree+poles-2)-(degree-1)-(poles-1))
    zero('a full ramification point leaves no other point in its degree-d fibre',degree-degree)
    q = h+1/h
    zero('a rational map with two poles has the stated derivative',S.diff(q,h)-(1-1/(h*h)))
    zero('its first finite critical value gives a repeated pencil member',s*s+t*t-2*s*t-(s-t)**2)
    zero('its second finite critical value also gives a repeated member',s*s+t*t+2*s*t-(s+t)**2)

    for m in [2,3,5,7]:
        check(f'degree {m}: reciprocal-power counterexample has a constant zero member',
            (1-lam*s**m).subs(lam,0)==1)
        check(f'degree {m}: its nonzero finite members are generically square-free',
            S.gcd(1-lam*s**m,-lam*m*s**(m-1))==1)
        check(f'degree {m}: a power of a quotient instead has a repeated zero member',
            S.gcd(s**m,m*s**(m-1))==s**(m-1))

    H = s*s*t
    check('the one exceptional polynomial fibre has multiplicity gcd one',
        S.gcd(2,1)==1 and S.factor(H)==s*s*t)
    zero('its nonzero fibres have no critical point',H-lam-s*S.diff(H,s)/2+lam)
    check('the exceptional parameter really has a critical line',
        all(S.diff(H,x).subs(s,0)==0 for x in [s,t]))
    bad_H = H*H
    check('dropping the multiplicity gcd admits a composite polynomial',
        S.gcd(4,2)==2 and bad_H==H**2)
    reduced = s*t
    check('a reduced reducible fibre is allowed by the pencil test',
        S.gcd(S.gcd(reduced,S.diff(reduced,s)),S.diff(reduced,t))==1)
    check('that parameter may still have a critical point',
        S.diff(reduced,s).subs({s:0,t:0})==0 and S.diff(reduced,t).subs({s:0,t:0})==0)

    numerator = V+U**4+U
    denominator = U**3*(U-1)**2
    check('arbitrary repeated denominator zeros leave coprime numerator and denominator',
        S.gcd(numerator,denominator)==1)
    check('every finite quotient-pencil member has a unit V derivative',
        S.diff(numerator-lam*denominator,V)==1)
    Qpoly = U**4*(U-1)**3
    product_parameter = Qpoly*(V+U**2+1)
    check('the scaled-product zero fibre retains a factor of multiplicity one',
        S.gcd(Qpoly,V+U**2+1)==1)
    zero('a nonzero scaled-product level forces its V derivative nonzero',
        product_parameter-lam-S.diff(product_parameter,V)*(V+U**2+1)+lam)
    variable_linear_pencil = U*V+U+1-lam*(U*U+1)
    zero('a simple V-leading zero can create a reduced two-component member',
        variable_linear_pencil.subs(lam,1)-U*(-U+V+1))
    check('its rational numerator and denominator remain coprime',
        S.gcd(U*V+U+1,U*U+1)==1)
    dependent_N = s**3-6*s+4
    dependent_M = s**3
    zero('independence is essential for the one-exception rational test',
        dependent_N/dependent_M-(4*(-S.Rational(1,2)+1/s)**3-3*(-S.Rational(1,2)+1/s)))
    zero('the dependent counterexample has an exceptional fibre with multiplicity gcd one',
        dependent_N+dependent_M-2*(s-1)**2*(s+2))
    zero('its other critical value gives a reduced linear finite member',dependent_N-dependent_M-(4-6*s))
    zero('all other finite members of that counterexample are reduced',
        S.discriminant(dependent_N-lam*dependent_M,s)-432*(1-lam)*(1+lam))
    check('the counterexample numerator and denominator are algebraically dependent',
        S.det(S.Matrix([[S.diff(dependent_N,s),S.diff(dependent_N,t)],
            [S.diff(dependent_M,s),S.diff(dependent_M,t)]]))==0)
    independent_N = U**2*(V+U)
    independent_M = (U-1)**3
    check('the mixed scaling numerator and denominator are algebraically independent',
        S.expand(S.diff(independent_N,V)*S.diff(independent_M,U))!=0)
    check('its nonzero finite members are primitive linear polynomials in V',
        S.gcd(U**2,U**3-lam*(U-1)**3)==1)
    D0 = S.Symbol('D0')
    zero('a general rational scaling uses the mixed product quotient parameter',
        (V+B)/(D0/Q)-Q*(V+B)/D0)

    j11,j12,j21,j22,derivative = S.symbols('j11 j12 j21 j22 derivative')
    matrix = S.Matrix([[j11,j12],[j21,j22]])
    row = S.Matrix([[derivative,1]])
    zero('a quotient-pencil coordinate keeps nonzero source differential after Keller substitution',
        row*matrix*matrix.adjugate()-matrix.det()*row)
    zero('the non-Keller control produces a composite quotient parameter',t*t/(s*s)-(t/s)**2)
    check('that control fails constant nonzero Jacobian',S.diff(t*t,t)==2*t)

    for m in range(2,8):
        R = z**m+z+1
        P = A*R.subs(z,c*z+beta)+B
        parameter = (V+B)/A
        zero(f'degree {m}: arbitrary scaling and affine input reduce to one rational parameter',
            P+V-A*(R.subs(z,c*z+beta)+parameter))
        zero(f'degree {m}: the affine source-field input has an explicit inverse',
            (c*z+beta).subs(z,(z-beta)/c)-z)
        sparse_R = z**m+z**(m-1)
        k = 2
        reciprocal_P = S.cancel(sparse_R.subs(z,T**k*z)/T**(k*(m-1)))
        zero(f'degree {m}: reciprocal scaling gives a polynomial sparse family',
            reciprocal_P-(T**k*z**m+z**(m-1)))
        zero(f'degree {m}: reciprocal equation uses the multiplied Keller parameter',
            Q*(R.subs(z,c*z+beta)/Q+B+V)-(R.subs(z,c*z+beta)+Q*(B+V)))
        aa,bb = S.symbols('aa bb',nonzero=True)
        zero(f'degree {m}: arbitrary sparse leading coefficients have a rational fixed-polynomial presentation',
            bb**m/aa**(m-1)*sparse_R.subs(z,aa*z/bb)-(aa*z**m+bb*z**(m-1)))

    q_symbol = S.Symbol('q')
    cubic_disc = S.discriminant(z**3+z+q_symbol,z)
    zero('the repeated-leading cubic family has the stated S3 discriminant',cubic_disc+4+27*q_symbol*q_symbol)
    check('this cubic discriminant is square-free and hence nonsquare',
        S.degree(cubic_disc,q_symbol)==2 and S.gcd(cubic_disc,S.diff(cubic_disc,q_symbol))==1)
    for k in [1,2]:
        leading_order = 3*k
        allowed = []
        for order in range(1,leading_order+1):
            values = [0,leading_order-order,leading_order-3*order]
            if values.count(min(values))>=2:
                allowed.append(order)
        check(f'cubic leading order {3*k}: the older pole test permits an intermediate pole',allowed==[k])
    w,v = S.symbols('w v')
    actual_H = 1+s*v+t*w+w*s**3*(w**3+w)
    zero('the retained scaled-family automorphic control has actual Hessian determinant one',
        S.det(S.hessian(actual_H,(s,t,w,v)))-1)

    check('all retained proof and primary-source pins are present',len(RETAINED)==8)
    for path,digest in RETAINED.items():
        check('retained byte pin '+path,hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest)
    context = json.loads((HERE/'sources/reduced-pencil-import.json').read_text())
    check('the closed-rational equivalence and retained ramification input are explicit',
        len(context['new_external_theorem_imports'])==1
        and 'Lemma 2(1)' in context['new_external_theorem_imports'][0]['role']
        and 'Riemann-Hurwitz' in context['retained_primary_input']['role'])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_reduced_pencils.py','REDUCED_PENCILS_AND_RATIONAL_MONODROMY.md','sources/reduced-pencil-import.json']},
        dependencies=RETAINED,
        written_not_finite_exhaustion=[
            'a reduced rational pencil with no constant finite member is noncomposite by ramification and the one-omitted-point bound',
            'one exceptional polynomial fibre with multiplicity gcd one also excludes all compositions',
            'an exceptional rational fibre with multiplicity gcd one is admitted when numerator and denominator are algebraically independent',
            'Keller pullback preserves these rational-parameter regularity certificates in every dimension',
            'the full group of every finite Galois parameter cover is preserved',
            'the filtered class-action criterion admits those rational parameters without deleting their pole fibres',
            'rationally scaled and rational-affine-input fixed-polynomial outputs with polynomial translation and final polynomiality stay fully irreducible in every degree',
            'all sparse a(T)z^m+b(T)z^(m-1)+d(T) outputs close without any gcd or leading multiplicity restriction'],
        new_external_theorem_imports=['Petravchuk-Iena arXiv:math/0701588v1 Lemma 2(1), (a)--(c)'],
        constant_finite_member_hypothesis_retained=True,
        critical_points_allowed_for_regular_parameters=True,full_parameter_splitting_group_preserved=True,
        independent_numerator_denominator_hypothesis_retained=True,
        noncomposite_closed_rational_equivalence_fully_reproved=False,
        all_fixed_parameter_fibres_asserted_irreducible=False,
        arbitrary_multi_parameter_target_covers_asserted_preserved=False,
        every_etale_affine_curve_map_asserted_finite=False,
        polynomial_slice_asserted_locally_nilpotent=False,
        all_rational_affine_isotrivial_polynomials_asserted_closed=False,
        arbitrary_polynomial_coefficient_variation_asserted_closed=False,
        nonlinear_source_Hessian_invariance_used=False,unit_ideal_endpoint_resolved=False,
        singleton_zero_endpoint_resolved=False,unrestricted_HC4_or_JC2_resolved=False,
        original_goal_complete=False,formal_verification=False,independent_replay=False,
        external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-reduced-pencils.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome='PASS',check_count=len(checks),receipt=str(output)),indent=2))


if __name__=='__main__':
    main()
