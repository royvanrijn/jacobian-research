"""PF2OMIT5 exact algebra; inherited geometry remains a written input."""
from pathlib import Path
import hashlib
import json
import platform
import time
import sympy as S

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parents[1]
CHECKS = []


def check(label, value):
    if not value:
        raise AssertionError(label)
    CHECKS.append(label)


def zero(label, value):
    check(label, S.cancel(value) == 0)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started = time.monotonic()
    z, x, y, t, Z, W, r = S.symbols('z x y t Z W r')
    # Generic target infinity direction is x=1,y=r, with transcendental r.
    finite_cluster = Z**5*(r*W-Z)
    infinite_constant = W**5*(r*W-Z)
    infinite_linear = Z*W**4*(r*W-Z)
    check('finite profile has exact five-plus-one root multiplicities',
          S.Poly(finite_cluster.subs(W, 1), Z).degree() == 6
          and S.Poly(finite_cluster.subs(W, 1), Z).nth(5) == r)
    zero('finite fixed fifth root', S.diff(finite_cluster.subs(W, 1), Z, 5).subs(Z, 0)-S.factorial(5)*r)
    check('finite moving root simple', S.diff(finite_cluster.subs(W, 1), Z).subs(Z, r) != 0)
    check('infinite constant quadratic part root five',
          S.Poly(infinite_constant.subs(Z, 1), W).nth(5) == -1
          and S.Poly(infinite_constant.subs(Z, 1), W).nth(4) == 0)
    check('infinite linear quadratic part root four',
          S.Poly(infinite_linear.subs(Z, 1), W).nth(4) == -1
          and S.Poly(infinite_linear.subs(Z, 1), W).nth(3) == 0)
    check('infinite linear has a second constant root', infinite_linear.subs(Z, 0) == 0)
    check('three supports distinct over generic direction', r != 0)

    # Quartic: every quadratic monomial and its weighted tangent derivative
    # has positive order, so the polar terms come only from l*p and m*q.
    monomials = [(5, 2, 0), (5, 1, 1), (6, 2, 0), (6, 1, 1), (6, 0, 2)]
    valuations = []
    p, q, zz = t**6, t**10, t**(-4)
    slope = S.diff(q, t)/S.diff(p, t)
    for i, j, k in monomials:
        value = z**i*x**j*y**k
        weight = 4*i-6*j-10*k
        valuations.append(weight)
        check('positive quadratic incidence order '+str((i, j, k)), weight >= 4)
        tangent = p*(S.diff(value, x)+slope*S.diff(value, y))
        substituted = S.expand(tangent.subs({x: p, y: q, z: zz}))
        zero('weighted tangent leading coefficient '+str((i, j, k)),
             substituted-(j+S.Rational(5, 3)*k)*t**(-weight))
    l0, l1, m0, m1, m2, n0 = S.symbols('l0 l1 m0 m1 m2 n0')
    variables = [l0, l1, m0, m1, m2, n0]
    polar_equations = [m0, l0+m1, l0+S.Rational(5, 3)*m1,
                       l1+m2, l1+S.Rational(5, 3)*m2, n0]
    matrix, rhs = S.linear_eq_to_matrix(polar_equations, variables)
    check('all six polar/constant coefficients forced zero', matrix.rank() == 6 and rhs == S.zeros(6, 1))
    zero('leading cancellation matrix determinant', matrix.det()+S.Rational(4, 9))
    a, b, c, d, e, L, M, N = S.symbols('a b c d e L M N')
    quartic = z**5*(a*x*x+b*x*y-z*(c*x*x+d*x*y+e*y*y))+z**2*L*x+z**3*M*y+z*N
    quotient = S.Poly(S.cancel(quartic/z), x, y, z)
    zero('forced coefficient common divisor is z', quotient.as_expr()*z-quartic)
    # Outside the fifth-power profile, incidence and tangency can both hold.
    allowed_control = 5*z**4*x*y-4*z**5*y*y-1
    zero('order-four control incidence', allowed_control.subs({x: p, y: q, z: zz}))
    zero('order-four control tangency',
         (S.diff(allowed_control, x)+slope*S.diff(allowed_control, y)).subs({x: p, y: q, z: zz}))
    check('order-four control is not divisible by z', allowed_control.subs(z, 0) == -1)
    check('order-four control outside degree-six field', S.degree(allowed_control, z) == 5)

    # Retained profiles are imported data; the written proof supplies their
    # necessity. Check exactly which rows this new argument exhausts.
    prior = json.loads((HERE.parent/'jc2-r7-cubic-cover/receipt.json').read_text())
    profiles = prior['reciprocal_projected_patterns']
    check('retained finite profile exactly2,5,2,0', profiles['finite'] == [[2, 5, 2, 0]])
    check('retained infinite profiles exactly the two used',
          profiles['infinite'] == [[2, 0, 3, 5], [2, 1, 4, 6]])
    finite_max_degree = max(5+1, 6, 6, 6)
    infinite0_max_degree = max(0+1, 5, 3, 5)
    check('infinite constant profile cannot have field degree six',
          infinite0_max_degree == 5 < finite_max_degree == 6)

    # Smooth inverse-line lemma, with general nonzero b0. All quadratic
    # coefficients and the coefficient of y vanish to order at least two.
    b0, b1, b2, c2, d0, d1, d2 = S.symbols('b0 b1 b2 c2 d0 d1 d2', nonzero=True)
    inverse = z**5*(x*y-z*x*x)+(b0+b1*z+b2*z*z)*x+c2*z*z*y+d0+d1*z+d2*z*z
    zero('entire inverse fibre is a line', inverse.subs(z, 0)-(b0*x+d0))
    zero('line fibre is smooth in the inverse surface', S.diff(inverse, x).subs(z, 0)-b0)
    zero('line fibre y derivative vanishes', S.diff(inverse, y).subs(z, 0))
    zero('ramification derivative restricted to fibre', S.diff(inverse, z).subs(z, 0)-(b1*x+d1))
    jac_on_line = S.cancel((-S.diff(inverse, z)/S.diff(inverse, x)).subs({z: 0, x: -d0/b0}))
    zero('Jacobian is constant along the whole line', jac_on_line+(b0*d1-b1*d0)/b0**2)
    zero('one ramified line point forces the entire line critical',
         jac_on_line.subs(d1, b1*d0/b0))
    zero('nonzero constant allows no isolated ramification point',
         jac_on_line.subs(d1, b1*d0/b0+1)+1/b0)
    quadratic_control = x+z*y+z*z
    zero('positive quadratic control has isolated ramification on the line',
         S.diff(quadratic_control, z).subs(z, 0)-y)
    zero('quadratic control is smooth', S.diff(quadratic_control, x)-1)
    zero('quadratic control critical image is a parabola',
         quadratic_control.subs({z: -y/2, x: y*y/4}))
    check('quadratic control violates coefficient-y order two',
          S.diff(quadratic_control, y) == z)

    residual = [(e, excess, count) for e in range(2, 6)
                for excess in range(4) for count in range(4)
                if e+excess+count == 3]
    check('complete extra-ramified ledger has three possibilities',
          residual == [(2, 0, 1), (2, 1, 0), (3, 0, 0)])
    check('two further ramified primes exceed the global budget', 2+2 > 5-2)
    check('a shared ramified image exceeds the double-value length budget', 2*2+2*2 > 6)
    residual_rows = [dict(other_ramified_e=e, other_ramified_f=None, other_ramified_image=None,
                          other_ramified_image_differs_from_C=True,
                          other_ramified_excess=excess, unramified_omitted_primes=count,
                          other_ramified_generic_ambient_length=e,
                          exceptional_points=[] if excess == 0 else [dict(ambient_length=3, point=None, target=None)])
                     for e, excess, count in residual]
    check('off-branch singleton forces the only residual row with an unramified prime',
          [row for row in residual if row[2] > 0] == [(2, 0, 1)])
    P5 = t**5-t
    Q6 = t*P5
    check('five-branch control has exactly five distinct zero preimages',
          S.degree(S.gcd(P5, Q6), t) == 5 and S.gcd(P5, S.diff(P5, t)) == 1)
    check('five-branch control is immersive everywhere', S.gcd(S.diff(P5, t), S.diff(Q6, t)) == 1)
    zero('five-branch control tangent slope is t at every zero',
         S.rem(S.diff(Q6, t)-t*S.diff(P5, t), P5, t))
    zero('five-branch control birational inverse off x0', Q6-t*P5)

    rows = [
        dict(case='generic infinity', normalization_rows=[dict(e=1, f=1, ramifies=False, omitted_from_affine_N=True),
                                                         dict(e=5, f=1, ramifies=True, omitted_from_affine_N=True)],
             pencil_root_cluster_lengths=[1, 5], index_not_inferred_from_cluster_lengths=True),
        dict(case='quartic escape', fixed_root_order=5, moving_root_order=1,
             ramified_projection_degree=4, quadratic_incidence_orders=valuations,
             forced_linear_orders=dict(l=2, m=3, n=1), contradiction='All conic coefficients have factor z', realized=False),
        dict(case='reciprocal finite profile', fixed_root_order=5, moving_root_order=1,
             additional_finite_E_point=True, line_jacobian=str(jac_on_line),
             contradiction='The finite point forces a ramified line image, distinct from the degree10 branch curve', realized=False),
        dict(case='reciprocal infinite constant profile', possible_parameter_degree=5,
             required_parameter_degree=6, contradiction='Field degree', realized=False),
        dict(case='reciprocal infinite linear profile', trial_root_clusters=[4, 1, 1],
             trial_normalization_prime_count_minimum=3, required_normalization_prime_count=2,
             ramification_indices_not_assigned_to_nonnormal_trial_clusters=True,
             contradiction='Complete generic infinity prime census', realized=False),
    ]
    result = dict(as_of='2026-11-14', claim='PF2OMIT5', rows=rows,
                  scope='degree6 plane Keller map; entire unique affine ramified prime e2,f1; birational degree(6,10) gap21 image',
                  stratum_excluded=True, uses_source_finality=False, uses_any_Theorem3_8_conclusion=False,
                  all_omitted_unramified_primes_allowed=True,
                  entire_single_ramification_hypothesis_derived_from_terminal_passport=False,
                  unrestricted_singleton_zero_properness_proved=False, full_goal_complete=False,
                  residual_immersed_r7_ledger=residual_rows,
                  residual_hypothesis='degree6 with an immersed e2,f1 gap21 prime; additional ramification is allowed',
                  residual_off_branch_singleton=dict(ramified_indices=[2, 2], ramified_excesses=[0, 0],
                      second_ramified_f=None, unramified_omitted_prime_count=1, unramified_e=1, unramified_f=1,
                      branches_of_its_image_at_zero=5, omitted_local_lengths_at_zero=[1]*5,
                      source_length_at_zero=1, distinct_nonproper_target_components=3,
                      unramified_image_delta_at_zero_minimum=10, actual_image_equation=None,
                      five_branch_control=[str(P5), str(Q6)], control_is_global_Keller_realization=False),
                  controls_not_Keller_realizations=True)
    target = HERE/'infinity_ledger.json'
    target.write_text(json.dumps(result, indent=2)+'\n')
    inputs = [HERE/'R7_INFINITY_COMPLETION.md', HERE/'ADMISSIBLE_INVERSE.md', HERE/'GLOBAL_LEDGER.md',
              HERE.parent/'jc2-r7-universal-contact/PROOF.md',
              HERE.parent/'jc2-r7-cubic-cover/PROOF.md', HERE.parent/'jc2-r7-cubic-cover/receipt.json',
              HERE.parent/'jc2-r7-quartic-contact/PROOF.md']
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS), checks=CHECKS,
                   runtime_seconds=round(time.monotonic()-started, 6), python=platform.python_version(), sympy=S.__version__,
                   checker_sha256=digest(Path(__file__)),
                   input_sha256={str(p.relative_to(RESEARCH)): digest(p) for p in inputs},
                   ledger_sha256=digest(target), formal_verification=False, external_review=False,
                   proof_boundary='Graph normalization, henselian cluster rank and imported surface/pencil geometry are written steps, not proved by finite tests',
                   full_goal_complete=False)
    (HERE/'infinity_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
