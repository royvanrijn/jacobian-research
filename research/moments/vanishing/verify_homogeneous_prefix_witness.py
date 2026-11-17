"""Certify an isolated algebraic cubic with first nonzero pure moment14.

Arb/Acb bounds have outward rounding. Rational input and quantitative
contraction inequalities, not floating point residuals, certify the root.
"""
from pathlib import Path
from math import factorial, prod
from fractions import Fraction
import hashlib
import json
import time
import flint
from flint import acb, acb_mat, arb, ctx, fmpq

P = Path(__file__).resolve().parent
BASE = [((3, 0, 0), 1), ((1, 1, 1), 1)]
BASIS = [[((0, 3, 0), 1)], [((2, 0, 1), 1)], [((0, 2, 1), 1)],
         [((2, 1, 0), 1), ((1, 0, 2), -2)],
         [((1, 2, 0), 1), ((0, 1, 2), -2)],
         [((0, 0, 3), 1), ((1, 1, 1), -3)]]


def double_factorial(n):
    return prod(range(n, 0, -2))


def polynomial(parameters):
    result = {power: acb(value) for power, value in BASE}
    for parameter, family in zip(parameters, BASIS):
        for power, factor in family:
            result[power] = result.get(power, acb(0))+factor*parameter
    return result


def normalized_wick(poly, total, shift=(0, 0, 0)):
    value = acb(0)
    for power, coefficient in poly.items():
        i, j, k = (power[t]+shift[t] for t in range(3))
        if i == j and k % 2 == 0:
            value += coefficient*(factorial(i)*double_factorial(k-1))
    return value/double_factorial(total+1)


def evaluate(parameters, maximum=14):
    terms = polynomial(parameters)
    power = {(0, 0, 0): acb(1)}
    values = []
    jacobian = []
    for n in range(1, maximum+1):
        if n % 2 == 0 and n <= 12:
            jacobian.append([n*sum((factor*normalized_wick(power, 3*n, shift)
                                    for shift, factor in family), acb(0))
                             for family in BASIS])
        result = {}
        for old, coefficient in power.items():
            for add, factor in terms.items():
                monomial = tuple(old[i]+add[i] for i in range(3))
                result[monomial] = result.get(monomial, acb(0))+coefficient*factor
        power = result
        if n % 2 == 0:
            values.append(normalized_wick(power, 3*n))
    return values, acb_mat(jacobian)


def gradient_matrix(parameters):
    quadratics = [(i, j, 2-i-j) for i in range(3) for j in range(3-i)]
    quartics = [(i, j, 4-i-j) for i in range(5) for j in range(5-i)]
    indices = {power: i for i, power in enumerate(quartics)}
    columns = []
    pol = polynomial(parameters)
    for variable in range(3):
        for multiplier in quadratics:
            column = [acb(0)]*15
            for power, coefficient in pol.items():
                if power[variable]:
                    exponent = tuple(power[i]+multiplier[i]-(i == variable) for i in range(3))
                    column[indices[exponent]] += power[variable]*coefficient
            columns.append(column)
    return acb_mat([[columns[j][i] for j in range(18)] for i in range(15)])


def rational(value):
    q = Fraction(value)
    return fmpq(q.numerator, q.denominator)


def rational_complex(pair, radius='0'):
    return acb(arb(rational(pair[0]), rational(radius)), arb(rational(pair[1]), rational(radius)))


def max_abs(values):
    return max(value.abs_upper().fmpq() for value in values)


def main():
    start = time.monotonic()
    ctx.prec = 320
    path = P/'homogeneous_prefix_witness_certificate.json'
    certificate = json.loads(path.read_text())
    assert certificate['equations'] == [2, 4, 6, 8, 10, 12]
    assert certificate['normalization'] == 'M_n/(3n+1)!!'
    assert certificate['coefficient_order'] == ['b', 'c', 'd', 'e', 'f', 'h']
    center = [rational_complex(pair) for pair in certificate['center']]
    radius = rational(certificate['complex_sup_radius'])
    assert radius > 0
    boxes = [rational_complex(pair, str(radius)) for pair in certificate['center']]
    B = acb_mat([[rational_complex(pair) for pair in row] for row in certificate['preconditioner']])
    identity = acb_mat([[int(i == j) for j in range(6)] for i in range(6)])
    center_values, center_jacobian = evaluate(center, maximum=12)
    box_values, box_jacobian = evaluate(boxes)
    eta = max_abs(B*acb_mat([[v] for v in center_values]))
    error = identity-B*box_jacobian
    q = max(sum((error[i, j].abs_upper().fmpq() for j in range(6)), fmpq(0)) for i in range(6))
    assert q < fmpq(1, 2)
    assert eta < radius/4
    assert eta+q*radius < radius
    assert eta < rational(certificate['claimed_eta_upper'])
    assert q < rational(certificate['claimed_q_upper'])
    m14 = box_values[-1]
    assert m14.abs_lower() > 0
    for coordinate, bounds in [(m14.real, certificate['normalized_M14_real_interval']),
                               (m14.imag, certificate['normalized_M14_imag_interval'])]:
        assert coordinate.lower() > arb(rational(bounds[0]))
        assert coordinate.upper() < arb(rational(bounds[1]))
    columns = certificate['smooth_gradient_minor_columns']
    assert len(columns) == len(set(columns)) == 15 and min(columns) >= 0 and max(columns) < 18
    gradient = gradient_matrix(boxes)
    minor = acb_mat([[gradient[i, j] for j in columns] for i in range(15)]).det()
    assert minor.abs_lower() > 0
    # Exact integer regressions use an independent symbolic polynomial power.
    import sympy as s
    Z, W, T = s.symbols('Z W T')
    symbolic = Z**3+Z*W*T
    for value, family in zip([1, -2, 3, 1, -1, 2], BASIS):
        symbolic += sum(value*factor*Z**i*W**j*T**k for (i, j, k), factor in family)
    assert s.expand(2*s.diff(symbolic, Z, W)+s.diff(symbolic, T, 2)) == 2*T
    integer_values, _ = evaluate([acb(v) for v in [1, -2, 3, 1, -1, 2]], maximum=6)
    regressions = []
    for n, actual in zip([2, 4, 6], integer_values):
        expected = 0
        for (i, j, k), coefficient in s.Poly(s.expand(symbolic**n), Z, W, T).terms():
            if i == j and k % 2 == 0:
                expected += coefficient*factorial(i)*double_factorial(k-1)
        expected = fmpq(int(expected), double_factorial(3*n+1))
        assert actual.contains(acb(expected))
        regressions.append({'n': n, 'normalized_value': str(expected)})
    report = {
        'completed': True, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'python_flint_version': flint.__version__, 'precision_bits': ctx.prec,
        'runtime_seconds': round(time.monotonic()-start, 3),
        'root_method': 'Banach contraction for x -> x-B*F(x) in a complex sup-norm polydisc',
        'complex_sup_radius': str(radius), 'eta_upper': str(eta), 'q_upper': str(q),
        'eta_lt_radius_over_4': True, 'q_lt_half': True, 'strict_self_mapping': True,
        'normalized_M14_enclosure': str(m14), 'normalized_M14_abs_lower': str(m14.abs_lower().fmpq()),
        'smooth_gradient_minor_enclosure': str(minor),
        'smooth_gradient_minor_abs_lower': str(minor.abs_lower().fmpq()),
        'first_nonzero_pure_moment': 14, 'initial_zero_pure_moments': 13,
        'homogeneous': True, 'projective_leading_curve_smooth': True,
        'coefficients_specified_by_unique_isolated_algebraic_root': True,
        'integer_regressions': regressions,
        'degree14_sufficiency_proved': False, 'odd19_zero_witness_provided': False,
        'boundary': 'The rational polynomial system and certified isolating polydisc specify algebraic coefficients. The exact zero moments follow from a written contraction argument with rigorous Acb inequalities. No claim of radical membership, global sufficiency, independent replay, external review or formal verification.'}
    (P/'homogeneous_prefix_witness_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
