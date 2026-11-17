"""Rigorous interval acceptance gate for a cubic with first nonzero moment20.

Exact rational isolating data plus a contraction proof certify zero moments.
The numerical exploration is not used by this checker as proof input.
"""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
from math import factorial, prod
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
         [((0, 0, 3), 1), ((1, 1, 1), -3)],
         [((1, 0, 0), 1)], [((0, 1, 0), 1)], [((0, 0, 1), 1)]]


def dfact(n):
    return prod(range(n, 0, -2))


def rational(value):
    q = Fraction(value)
    return fmpq(q.numerator, q.denominator)


def rational_complex(pair, radius='0'):
    return acb(arb(rational(pair[0]), rational(radius)), arb(rational(pair[1]), rational(radius)))


def polynomial(parameters):
    result = {power: acb(value) for power, value in BASE}
    for parameter, family in zip(parameters, BASIS):
        for power, factor in family:
            result[power] = result.get(power, acb(0))+factor*parameter
    return result


@lru_cache(None)
def wick_map(power, shift=(0, 0, 0)):
    entries = []
    for i in range(3*power+1):
        j = i+shift[0]-shift[1]
        if j < 0:
            continue
        for k in range(max(0, power-i-j), 3*power-i-j+1):
            if (i+j+k-power) % 2 or (k+shift[2]) % 2:
                continue
            entries.append(((i, j, k), factorial(i+shift[0])*dfact(k+shift[2]-1)))
    return entries


def wick(poly, power, moment, shift=(0, 0, 0)):
    out = acb(0)
    for monomial, weight in wick_map(power, shift):
        value = poly.get(monomial)
        if value is not None:
            out += weight*value
    return out/dfact(3*moment+1)


def evaluate(parameters, maximum=20):
    terms = polynomial(parameters)
    power = {(0, 0, 0): acb(1)}
    values = []
    jacobian = []
    for n in range(1, maximum+1):
        if n % 2 == 0 and n <= 18:
            jacobian.append([n*sum((factor*wick(power, n-1, n, shift)
                                    for shift, factor in family), acb(0)) for family in BASIS])
        result = {}
        for (i, j, k), coefficient in power.items():
            for (a, b, c), factor in terms.items():
                monomial = (i+a, j+b, k+c)
                result[monomial] = result.get(monomial, acb(0))+coefficient*factor
        power = result
        if n % 2 == 0:
            values.append(wick(power, n, n))
    return values, acb_mat(jacobian)


def gradient_matrix(parameters):
    quadratics = [(i, j, 2-i-j) for i in range(3) for j in range(3-i)]
    quartics = [(i, j, 4-i-j) for i in range(5) for j in range(5-i)]
    indices = {power: i for i, power in enumerate(quartics)}
    columns = []
    pol = {power: value for power, value in polynomial(parameters).items() if sum(power) == 3}
    for variable in range(3):
        for multiplier in quadratics:
            column = [acb(0)]*15
            for power, coefficient in pol.items():
                if power[variable]:
                    exponent = tuple(power[i]+multiplier[i]-(i == variable) for i in range(3))
                    column[indices[exponent]] += power[variable]*coefficient
            columns.append(column)
    return acb_mat([[columns[j][i] for j in range(18)] for i in range(15)])


def max_abs(values):
    return max(value.abs_upper().fmpq() for value in values)


def main():
    start = time.monotonic()
    ctx.prec = 384
    path = P/'odd_prefix_witness_certificate.json'
    certificate = json.loads(path.read_text())
    assert certificate['equations'] == list(range(2, 19, 2))
    assert certificate['normalization'] == 'M_n/(3n+1)!!'
    assert certificate['coefficient_order'] == ['b', 'c', 'd', 'e', 'f', 'h', 'u', 'v', 'w']
    center = [rational_complex(pair) for pair in certificate['center']]
    radius = rational(certificate['complex_sup_radius'])
    assert radius == rational('1e-40')
    boxes = [rational_complex(pair, str(radius)) for pair in certificate['center']]
    B = acb_mat([[rational_complex(pair) for pair in row] for row in certificate['preconditioner']])
    identity = acb_mat([[int(i == j) for j in range(9)] for i in range(9)])
    center_values, _ = evaluate(center, maximum=18)
    box_values, box_jacobian = evaluate(boxes)
    eta = max_abs(B*acb_mat([[value] for value in center_values]))
    error = identity-B*box_jacobian
    q = max(sum((error[i, j].abs_upper().fmpq() for j in range(9)), fmpq(0)) for i in range(9))
    assert q < fmpq(1, 2) and eta < radius/4 and eta+q*radius < radius
    assert eta < rational(certificate['claimed_eta_upper'])
    assert q < rational(certificate['claimed_q_upper'])
    m20 = box_values[-1]
    assert m20.abs_lower() > 0
    for coordinate, bounds in [(m20.real, certificate['normalized_M20_real_interval']),
                               (m20.imag, certificate['normalized_M20_imag_interval'])]:
        assert coordinate.lower() > arb(rational(bounds[0]))
        assert coordinate.upper() < arb(rational(bounds[1]))
    columns = certificate['smooth_gradient_minor_columns']
    assert len(columns) == len(set(columns)) == 15 and min(columns) >= 0 and max(columns) < 18
    gradient = gradient_matrix(boxes)
    minor = acb_mat([[gradient[i, j] for j in columns] for i in range(15)]).det()
    assert minor.abs_lower() > 0
    # Independent symbolic integer regression; retain lower-degree Wick terms.
    import sympy as s
    Z, W, T = s.symbols('Z W T')
    symbolic = Z**3+Z*W*T
    integers = [1, -2, 3, 1, -1, 2, 2, -1, 3]
    for value, family in zip(integers, BASIS):
        symbolic += sum(value*factor*Z**i*W**j*T**k for (i, j, k), factor in family)
    assert s.expand(2*s.diff(symbolic, Z, W)+s.diff(symbolic, T, 2)) == 2*T
    values, _ = evaluate([acb(v) for v in integers], maximum=6)
    regressions = []
    for n, actual in zip([2, 4, 6], values):
        expected = 0
        for (i, j, k), coefficient in s.Poly(s.expand(symbolic**n), Z, W, T).terms():
            if i == j and k % 2 == 0:
                expected += coefficient*factorial(i)*dfact(k-1)
        expected = fmpq(int(expected), dfact(3*n+1))
        assert actual.contains(acb(expected))
        regressions.append({'n': n, 'normalized_value': str(expected)})
    report = {
        'completed': True, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'python_flint_version': flint.__version__, 'FLINT_version': flint.__FLINT_VERSION__,
        'precision_bits': ctx.prec, 'runtime_seconds': round(time.monotonic()-start, 3),
        'root_method': 'Banach contraction for x -> x-B*F(x) in a nine-dimensional complex sup-norm polydisc',
        'complex_sup_radius': str(radius), 'eta_upper': str(eta), 'q_upper': str(q),
        'eta_lt_radius_over_4': True, 'q_lt_half': True, 'strict_self_mapping': True,
        'normalized_M20_enclosure': str(m20), 'normalized_M20_abs_lower': str(m20.abs_lower().fmpq()),
        'smooth_gradient_minor_enclosure': str(minor), 'smooth_gradient_minor_abs_lower': str(minor.abs_lower().fmpq()),
        'first_nonzero_pure_moment': 20, 'initial_zero_pure_moments': 19,
        'odd': True, 'homogeneous': False, 'projective_leading_curve_smooth': True,
        'coefficients_specified_by_unique_isolated_algebraic_root': True,
        'integer_regressions': regressions, 'prefix20_sufficiency_proved': False,
        'boundary': 'Rational equations and a certified isolating polydisc specify algebraic coefficients with exact19 initial zero moments and nonzero20th. The cubic leading part is smooth by a full gradient multiplication minor. The finite-prefix obstruction does not prove global/class sufficiency, radical membership, independent replay, external review or formal verification.'}
    (P/'odd_prefix_witness_verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
