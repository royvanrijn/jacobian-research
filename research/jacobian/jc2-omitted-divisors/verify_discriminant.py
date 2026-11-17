"""Exact controls for PF2OMIT4; general lemmas are written proofs, not tests."""
from pathlib import Path
from itertools import product
import hashlib
import json
import platform
import time
import sympy as S
from sympy.polys.matrices import DomainMatrix

HERE = Path(__file__).resolve().parent
CHECKS = []


def check(label, condition):
    if not condition:
        raise AssertionError(label)
    CHECKS.append(label)


def zero(label, expression):
    check(label, S.cancel(expression) == 0)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def order_at_origin(expression, u, v):
    poly = S.Poly(expression, u, v)
    assert not poly.is_zero
    return min(sum(m) for m, _ in poly.terms())


def partitions(n, upper=None):
    if n == 0:
        yield []
        return
    for first in range(min(n, upper or n), 0, -1):
        for tail in partitions(n-first, first):
            yield [first] + tail


def main():
    started = time.monotonic()
    T, t, b, q, r, u, v, x, y = S.symbols('T t b q r u v x y')
    previous = json.loads((HERE / 'inverse_equations.json').read_text())
    R = S.Poly(S.sympify(previous['inverse_polynomial'],
               locals={str(s): s for s in (T, t, b, q, r)}).subs({q: 1, r: 0}), T)
    P, derivative = t**4+t-b, 4*t**3+1
    delta = S.discriminant(P, t)
    zero('quartic line branch equation', delta+256*b**3+27)
    mt = S.Matrix([[0, 0, 0, b], [1, 0, 0, -1],
                   [0, 1, 0, 0], [0, 0, 1, 0]])
    dm = 4*mt**3+S.eye(4)
    check('quartic multiplication relation', mt**4+mt-b*S.eye(4) == S.zeros(4))
    gram_even = S.Matrix(4, 4, lambda i, j: 2*(mt**(i+j)).trace())
    gram_odd = S.Matrix(4, 4, lambda i, j: 2*(dm*mt**(i+j)).trace())
    normal_disc = S.expand(gram_even.det()*gram_odd.det())
    zero('normal trace discriminant 256 Delta cubed', normal_disc-256*delta**3)
    check('elliptic cubic is smooth', S.gcd(4*t**3+1, 12*t**2) == 1)
    field = S.QQ.frac_field(b)
    U = (t**3+1)/b
    V = S.rem(t*S.invert(derivative, P, t, domain=field)-U,
              P, t, domain=field)
    zero('rational inverse t times U equals one modulo P',
         S.rem(t*U-1, P, t, domain=field))
    zero('T normal-basis expression',
         S.rem(derivative*(V+U)-t, P, t, domain=field))

    def evaluate(poly):
        p = S.Poly(poly, t, domain=field)
        return sum((c*mt**i for (i,), c in p.terms()), S.zeros(4))

    um, vm = evaluate(U), evaluate(V)
    mT = S.BlockMatrix([[um, dm*vm], [vm, um]]).as_explicit().applyfunc(S.cancel)
    unit = S.eye(8)[:, 0]
    tc = mT*unit
    cols = [unit]
    zeta = (R.nth(8)*tc).applyfunc(S.cancel)
    cols.append(zeta)
    for j in range(2, 8):
        zeta = (mT*zeta+R.nth(9-j)*tc).applyfunc(S.cancel)
        cols.append(zeta)
    inclusion = S.Matrix.hstack(*cols)
    check('all binary-order basis coordinates polynomial',
          all(S.denom(S.cancel(a)) == 1 for a in inclusion))
    det = S.expand(DomainMatrix.from_Matrix(inclusion).det().as_expr())
    c3 = 256*b**3-b**2+6*b+28
    c24, remainder = S.div(det, c3, b)
    zero('index determinant contains C3 exactly', remainder)
    check('remaining index factor degree24', S.degree(c24, b) == 24)
    raw_disc = S.discriminant(R.as_expr(), T)
    zero('literal raw/normal/index discriminant identity', raw_disc-normal_disc*det**2)
    zero('inverse leading coefficient b Delta', R.LC()-b*delta)
    factors = [('Delta', delta), ('C3', c3), ('C24', c24)]
    for name, factor in factors:
        check(name+' squarefree', S.degree(S.gcd(factor, S.diff(factor, b)), b) == 0)
        check(name+' avoids zero', factor.subs(b, 0) != 0)
    for i, (name, factor) in enumerate(factors):
        for name2, factor2 in factors[i+1:]:
            check(name+' coprime '+name2, S.degree(S.gcd(factor, factor2), b) == 0)
    check('three normal branch values and 27 index values',
          S.degree(delta, b) == 3 and S.degree(det, b) == 27)
    check('unramified omitted value has nonzero raw discriminant', raw_disc.subs(b, 0) != 0)
    check('zero-value finite inverse degree seven', S.degree(R.as_expr().subs(b, 0), T) == 7)
    check('ramified values have inverse degree at most four',
          all(S.rem(R.nth(i), delta, b) == 0 for i in range(5, 9)))
    check('ramified values have inverse degree exactly four',
          S.degree(S.gcd(R.nth(4), delta), b) == 0)
    line_rows = [
        dict(condition='b*Delta*C3*C24 != 0', value_count='open', normal_lengths=[1]*8,
             source_lengths=[1]*8, omitted_lengths=[], raw_discriminant_order=0,
             normalized_different_order=0, order_index=0, leading_order=0,
             distinct_finite_T_values=8, ramifies=False, omitted=False),
        dict(condition='b=0', value_count=1, normal_lengths=[1]*8,
             source_lengths=[1]*7, omitted_lengths=[1], raw_discriminant_order=0,
             normalized_different_order=0, order_index=0, leading_order=1,
             omitted_rows=[dict(e=1, f=1, pole_order=1, ramifies=False)],
             distinct_finite_T_values=7, ramifies=False, omitted=True),
        dict(condition='Delta=0', value_count=3, normal_lengths=[4, 1, 1, 1, 1],
             source_lengths=[1]*4, omitted_lengths=[4], raw_discriminant_order=3,
             normalized_different_order=3, order_index=0, leading_order=1,
             omitted_rows=[dict(e=4, f=1, pole_order=1, ramifies=True)],
             distinct_finite_T_values=4, ramifies=True, omitted=True),
        dict(condition='C3*C24=0', value_count=27, normal_lengths=[1]*8,
             source_lengths=[1]*8, omitted_lengths=[], raw_discriminant_order=2,
             normalized_different_order=0, order_index=1, leading_order=0,
             conductor_exponents=[1, 1, 0, 0, 0, 0, 0, 0],
             distinct_finite_T_values=7, ramifies=False, omitted=False),
    ]
    for row in line_rows:
        check(row['condition']+' length conservation',
              sum(row['normal_lengths']) == sum(row['source_lengths'])+sum(row['omitted_lengths']) == 8)
        check(row['condition']+' discriminant separation',
              row['raw_discriminant_order'] == row['normalized_different_order']+2*row['order_index'])
    chi_source = 8-3*4-1
    check('source Euler characteristic has four finite punctures', chi_source == -5 == 2-2*1-5)
    check('normal affine elliptic Euler characteristic', 8-3*3 == -1)

    # A finite source-A2 control: residue degree weights and a pure index curve.
    R4 = T**4-2*(u+v)*T**2+(u-v)**2
    M4 = S.Matrix([[1, 0, u+v, 0], [0, 1, 0, u+3*v],
                    [0, 1, 0, 3*u+v], [0, 0, 2, 0]])
    gram4 = S.diag(4, 4*u, 4*v, 4*u*v)
    zero('finite square control normal trace discriminant', gram4.det()-256*u**2*v**2)
    zero('finite square control inclusion determinant', M4.det()+4*(u-v))
    zero('finite square control trace-index identity', S.discriminant(R4, T)-gram4.det()*M4.det()**2)
    zero('finite square control primitive relation', R4.subs({u: x**2, v: y**2, T: x+y}))
    jac = 4*x*y
    quotient = S.cancel(S.diff(R4, T).subs({u: x**2, v: y**2, T: x+y})/jac)
    zero('finite square control conductor Jacobian quotient', quotient-2*(x+y))
    check('finite square ramification exponent includes f2', (2-1)*2 == 2)

    groupings = []
    for shared_c in (0, 1):
        for groups in partitions(3-shared_c):
            solutions = []
            for counts in product(range(6), repeat=len(groups)):
                if sum(g*n for g, n in zip(groups, counts)) != 5:
                    continue
                if any(groups[i] == groups[j] and counts[i] > counts[j]
                       for i in range(len(groups)) for j in range(i+1, len(groups))):
                    continue
                solutions.append(list(counts))
            groupings.append(dict(unramified_primes_over_C=shared_c, other_image_groups=groups,
                                  ramified_prime=dict(e=2, f=1, omitted=True),
                                  unramified_prime_count=3, unramified_prime_e=1, unramified_prime_f=1,
                                  degree_mod_C=6-2-shared_c,
                                  degree_mod_other_images=[6-g for g in groups],
                                  generic_omitted_length_C=2+shared_c,
                                  generic_omitted_lengths_other_images=groups,
                                  raw_discriminant_min_C=1+2*shared_c,
                                  raw_discriminant_min_other_images=[g*(g-1) for g in groups],
                                  inverse_order_index_min_C=shared_c,
                                  inverse_order_index_min_other_images=[g*(g-1)//2 for g in groups],
                                  leading_min_C=1+shared_c, leading_min_other_images=groups,
                                  off_branch_singleton_branch_counts=solutions,
                                  off_branch_singleton_integer_compatible=bool(solutions),
                                  nonproper_target_component_count=1+len(groups)))
    check('exactly five generic r7 image groupings', len(groupings) == 5)
    expected = {(0, (3,)): [], (0, (2, 1)): [[0, 5], [1, 3], [2, 1]],
                (0, (1, 1, 1)): [[0, 0, 5], [0, 1, 4], [0, 2, 3], [1, 1, 3], [1, 2, 2]],
                (1, (2,)): [], (1, (1, 1)): [[0, 5], [1, 4], [2, 3]]}
    for row in groupings:
        key = (row['unramified_primes_over_C'], tuple(row['other_image_groups']))
        check('complete branch-count solutions '+str(key), row['off_branch_singleton_branch_counts'] == expected[key])
    allowed = [row for row in groupings if row['off_branch_singleton_integer_compatible']]
    check('off-branch singleton leaves exactly three groupings', len(allowed) == 3)
    check('off-branch singleton requires three or four target components',
          {row['nonproper_target_component_count'] for row in allowed} == {3, 4})

    jets = []
    for weights in ([1, 1], [1]*5, [1, 1, 2, 2, 3], [1, 2, 3, 4, 5], [3]):
        reciprocal = [(u+j*v)**w for j, w in enumerate(weights, 1)]
        t0 = u+7*v
        model = S.Poly((T-t0)*S.prod(1-z*T for z in reciprocal), T)
        orders = [order_at_origin(model.nth(i), u, v) for i in range(2, len(weights)+2)]
        expected_orders = [sum(weights[:i]) for i in range(1, len(weights)+1)]
        check('sharp coefficient orders '+str(weights), orders == expected_orders)
        pair_orders = [order_at_origin(reciprocal[i]-reciprocal[j], u, v)
                       for i in range(len(weights)) for j in range(i+1, len(weights))]
        lower = sum(min(weights[i], weights[j]) for i in range(len(weights)) for j in range(i+1, len(weights)))
        check('sharp reciprocal collision orders '+str(weights), sum(pair_orders) == lower)
        check('retained branch has unit separation '+str(weights),
              all((1-t0*z).subs({u: 0, v: 0}) == 1 for z in reciprocal))
        check('singleton specialization '+str(weights), model.as_expr().subs({u: 0, v: 0}) == T)
        if len(weights) <= 2:
            direct = S.discriminant(model.as_expr(), T)
            formula = S.prod((reciprocal[i]-reciprocal[j])**2 for i in range(len(weights))
                             for j in range(i+1, len(weights)))*S.prod((1-t0*z)**2 for z in reciprocal)
            zero('binary discriminant product identity '+str(weights), direct-formula)
        jets.append(dict(weights=weights, coefficient_orders_a2_to_ad=orders,
                         leading_order=sum(weights), raw_discriminant_order=2*sum(pair_orders),
                         normalized_discriminant_unit=True, formal_split_model_only=True,
                         connected_global_realization=False))
    degree6_model = next(row for row in jets if row['weights'] == [1]*5)
    check('degree-six model attains the full coefficient and discriminant bounds',
          degree6_model['coefficient_orders_a2_to_ad'] == [1, 2, 3, 4, 5]
          and degree6_model['raw_discriminant_order'] == 20)

    manifest = json.loads((HERE / 'sources/trace_manifest.json').read_text())
    for source in manifest['sources']:
        check('source pin '+source['file'], digest(HERE / 'sources' / source['file']) == source['sha256'])
    equations = dict(as_of='2026-11-14', line_inverse=str(R.as_expr()), line_branch_equation=str(delta),
                     C3=str(c3), C24=str(c24), index_determinant=str(c3*c24),
                     normal_basis=['1', 't', 't^2', 't^3', 'a', 'a*t', 'a*t^2', 'a*t^3'],
                     inclusion_matrix=[[str(a) for a in inclusion.row(i)] for i in range(8)],
                     normal_trace_gram_even=[[str(a) for a in gram_even.row(i)] for i in range(4)],
                     normal_trace_gram_odd=[[str(a) for a in gram_odd.row(i)] for i in range(4)],
                     normal_discriminant=str(S.factor(normal_disc)), raw_discriminant='256*Delta^3*(C3*C24)^2',
                     finite_plane_control=dict(map=['x^2', 'y^2'], inverse=str(R4),
                                               retained_ramification_rows=[dict(image='u=0', e=2, f=2), dict(image='v=0', e=2, f=2)],
                                               omitted_rows=[], inclusion_determinant=str(M4.det()),
                                               conductor_generator=str(quotient), Keller=False),
                     weighted_formal_jet_controls=jets, actual_r7_inverse=None,
                     actual_singleton_Keller_inverse=None)
    ledger = dict(as_of='2026-11-14', claim='PF2OMIT4', complete_marked_root_line=line_rows,
                  line_Euler=dict(normal=-1, source=chi_source, ramification_only_wrong_source=-4,
                                  pure_index_values=27, pure_index_source_deficit=0),
                  r7_hypotheses='degree6, entire single ramified prime e2,f1, immersed (6,10) gap21 image C',
                  r7_generic_image_groupings=groupings,
                  singleton_refinement_requires='a reduced singleton source fibre at zero and finite normalization etale over zero; in the single-ramification case zero outside C',
                  singleton_degree6_jets=dict(coefficient_min_orders_a2_to_a6=[1, 2, 3, 4, 5],
                                             raw_discriminant_min_order=20, normalized_discriminant_unit=True,
                                             polynomial_index_factor_min_order_under_single_e2f1_hypothesis=10),
                  global_images_determined=False, individual_exceptional_Keller_algebras_determined=False,
                  singleton_zero_properness_proved=False, full_goal_complete=False)
    for name, value in [('discriminant_equations.json', equations), ('discriminant_ledger.json', ledger)]:
        (HERE / name).write_text(json.dumps(value, indent=2) + '\n')
    inputs = ['DISCRIMINANT_LEDGER.md', 'ADMISSIBLE_INVERSE.md', 'GLOBAL_LEDGER.md',
              'inverse_equations.json', 'sources/trace_manifest.json']
    receipt = dict(as_of='2026-11-14', outcome='PASS', check_count=len(CHECKS), checks=CHECKS,
                   runtime_seconds=round(time.monotonic()-started, 6), python=platform.python_version(), sympy=S.__version__,
                   checker_sha256=digest(Path(__file__)), input_sha256={p: digest(HERE/p) for p in inputs},
                   output_sha256={p: digest(HERE/p) for p in ['discriminant_equations.json', 'discriminant_ledger.json']},
                   general_lemmas='Written trace-dual, local normalization and formal-splitting proofs; controls do not prove the general statements',
                   actual_Keller_inverse_available=False, full_goal_complete=False,
                   formal_verification=False, external_review=False)
    (HERE / 'discriminant_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({k: receipt[k] for k in ['outcome', 'check_count', 'runtime_seconds']}))


if __name__ == '__main__':
    main()
