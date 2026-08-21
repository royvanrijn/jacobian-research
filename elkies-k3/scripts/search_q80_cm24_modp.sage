#!/usr/bin/env sage
"""Search for a q=80 CM24 seed over a small finite field.

This is a bounded experiment, not a proof.  It uses the exact q=80 ambient
normal form and the transported generic basis at CM24.  The first basis
section has profile (A1,A3,D5,E6)=(1,1,1,1) and P.O=0, forcing

    x1 = T*(a + b*T + (d-a-b)*T^2).

The second section has trivial component profile and P.O=0.  We first filter
ambient parameters by the residual-double-root condition and P1, then test
all degree-four x-coordinates for P2.  The remaining simple-pole P3 search is
kept separate so that this inexpensive filter can be audited independently.
"""

from argparse import ArgumentParser
from itertools import product

from sage.all import GF, PolynomialRing, matrix, vector


parser = ArgumentParser()
parser.add_argument("--prime", type=int, default=5)
parser.add_argument("--max-survivors", type=int, default=0)
parser.add_argument("--search-p3", action="store_true")
args = parser.parse_args()

field = GF(args.prime)
assert field.characteristic() not in (2, 3)
ring = PolynomialRing(field, "T")
T = ring.gen()


def polynomial_square_root(value):
    """Return one polynomial square root, or None."""
    if not value:
        return ring.zero()
    degree = value.degree()
    if degree % 2:
        return None
    leading = value.leading_coefficient()
    if not leading.is_square():
        return None
    half = degree // 2
    coefficients = [field.zero()] * (half + 1)
    coefficients[half] = leading.sqrt()
    for index in range(half - 1, -1, -1):
        target_degree = half + index
        known = sum(
            coefficients[left] * coefficients[target_degree-left]
            for left in range(index + 1, half)
            if index < target_degree-left < half
        )
        coefficients[index] = (
            value[target_degree] - known
        ) / (2 * coefficients[half])
    root = ring(coefficients)
    return root if root**2 == value else None


def ambient_surface(d, p, q, e):
    r = -3 * d**2 + 3 - p - q
    A = T**2 * (-3 + p*T + q*T**2 + r*T**3)

    # Hermite interpolation of the multiplicative square-root branch at T=1.
    s_ring = PolynomialRing(field, "s")
    s = s_ring.gen()
    A_one = s_ring(A(T=1+s))
    u = (A_one + 3*d**2) / (-3*d**2)
    branch = 2*d**3 * (1 + field(3)/2*u + field(3)/8*u**2 - field(1)/16*u**3)
    target = vector(field, [branch[j] for j in range(4)])
    matrix_rows = []
    for jet in range(4):
        matrix_rows.append([s_ring((1+s)**degree)[jet] for degree in range(4, 8)])
    jet_matrix = matrix(field, matrix_rows)
    fixed = vector(
        field,
        [s_ring(2*(1+s)**3 + e*(1+s)**8)[jet] for jet in range(4)],
    )
    b1, b2, b3, b4 = jet_matrix.solve_right(target-fixed)
    B = T**3 * (2 + b1*T + b2*T**2 + b3*T**3 + b4*T**4 + e*T**5)
    discriminant = 4*A**3 + 27*B**2
    residual, remainder = discriminant.quo_rem(T**7*(T-1)**4)
    assert not remainder and residual.degree() <= 5
    return A, B, residual


def first_sections(A, B, d):
    answers = []
    for a, b in product(field, repeat=2):
        x = T * (a + b*T + (d-a-b)*T**2)
        y = polynomial_square_root(x**3 + A*x + B)
        if y is None:
            continue
        # These are the visible Weierstrass-model incidences at D5, A3, E6.
        if x(0) != 0 or y(0) != 0 or x(1) != d or y(1) != 0:
            continue
        if x.degree() > 3 or y.degree() > 5:
            continue
        answers.append((x, y))
    return answers


def second_sections(A, B, d):
    answers = []
    for coefficients in product(field, repeat=5):
        if coefficients[0] == 0 or coefficients[4] == 0:
            continue
        x = ring(coefficients)
        if x(1) == d:
            continue
        y = polynomial_square_root(x**3 + A*x + B)
        if y is None or y.degree() > 6:
            continue
        answers.append((x, y))
    return answers


def third_sections(A, B, d, repeated):
    """Search the transported generic P3 with P.O=1 at CM24."""
    rho = -repeated[0]
    cubic_ring = PolynomialRing(field, "x")
    x_variable = cubic_ring.gen()
    singular_cubic = x_variable**3 + A(rho)*x_variable + B(rho)
    node_factor = singular_cubic.gcd(singular_cubic.derivative())
    if node_factor.degree() != 1:
        return []
    node = -node_factor[0]
    answers = []
    for pole in field:
        if pole in (field.zero(), field.one(), rho):
            continue
        Z = T-pole
        interpolation = matrix(field, [[1, 1], [rho**5, rho**6]])
        if not interpolation.is_invertible():
            continue
        for low in product(field, repeat=5):
            partial = ring(low)
            targets = vector(
                field,
                [
                    d*Z(1)**2-partial(1),
                    node*Z(rho)**2-partial(rho),
                ],
            )
            c5, c6 = interpolation.solve_right(targets)
            if c6 == 0:
                continue
            X = partial + c5*T**5 + c6*T**6
            if X(0) == 0 or X.gcd(Z) != 1:
                continue
            Y = polynomial_square_root(X**3 + A*X*Z**4 + B*Z**6)
            if Y is None or Y.degree() > 9 or Y.gcd(Z) != 1:
                continue
            if X(1) != d*Z(1)**2 or Y(1) != 0:
                continue
            if X(rho) != node*Z(rho)**2 or Y(rho) != 0:
                continue
            answers.append((X, Y, Z))
    return answers


def affine_intersection(left, right):
    """Return the visible affine intersection divisor of two sections."""
    X_left, Y_left, Z_left = left
    X_right, Y_right, Z_right = right
    x_difference = X_left*Z_right**2 - X_right*Z_left**2
    y_difference = Y_left*Z_right**3 - Y_right*Z_left**3
    divisor = x_difference.gcd(y_difference)
    return divisor.monic() if divisor else divisor


def orientation_audit(p1s, p2s, p3s):
    """Compare all signs with the transported pair-intersection pattern."""
    rows = []
    for p1_index, (x1, y1) in enumerate(p1s, 1):
        for p2_index, (x2, y2) in enumerate(p2s, 1):
            for p3_index, (X3, Y3, Z3) in enumerate(p3s, 1):
                for sign1, sign2, sign3 in product((1, -1), repeat=3):
                    section1 = (x1, sign1*y1, ring.one())
                    section2 = (x2, sign2*y2, ring.one())
                    section3 = (X3, sign3*Y3, Z3)
                    divisors = (
                        affine_intersection(section1, section2),
                        affine_intersection(section1, section3),
                        affine_intersection(section2, section3),
                    )
                    degrees = tuple(
                        -1 if not divisor else divisor.degree()
                        for divisor in divisors
                    )
                    if degrees == (2, 3, 1):
                        rows.append(
                            (
                                p1_index, p2_index, p3_index,
                                sign1, sign2, sign3,
                                divisors,
                            )
                        )
    return rows


tested = 0
cm_candidates = 0
p1_candidates = 0
survivors = []
nonzero = [value for value in field if value]
for d in nonzero:
    for p, q, e in product(field, repeat=3):
        tested += 1
        A, B, residual = ambient_surface(d, p, q, e)
        if residual.degree() != 5:
            continue
        repeated = residual.gcd(residual.derivative())
        if repeated.degree() != 1:
            continue
        cm_candidates += 1
        p1s = first_sections(A, B, d)
        if not p1s:
            continue
        p1_candidates += 1
        p2s = second_sections(A, B, d)
        if not p2s:
            continue
        p3s = third_sections(A, B, d, repeated.monic()) if args.search_p3 else []
        if args.search_p3 and not p3s:
            continue
        orientations = orientation_audit(p1s, p2s, p3s) if args.search_p3 else []
        record = (d, p, q, e, repeated.monic(), p1s, p2s, p3s)
        survivors.append(record)
        print(
            "Q80MODP|SURVIVOR|prime={}|d,p,q,e={},{},{},{}|double={}|P1={}|P2={}|P3={}".format(
                args.prime,
                d, p, q, e,
                repeated.monic(),
                ";".join(f"x:{x},y:{y}" for x, y in p1s),
                ";".join(f"x:{x},y:{y}" for x, y in p2s),
                ";".join(f"X:{X},Y:{Y},Z:{Z}" for X, Y, Z in p3s),
            ),
            flush=True,
        )
        print(f"Q80MODP|MODEL|A={A}|B={B}|residual={residual}", flush=True)
        for orientation in orientations:
            p1_index, p2_index, p3_index, sign1, sign2, sign3, divisors = orientation
            print(
                "Q80MODP|ORIENTATION|P1={}|P2={}|P3={}|signs={},{},{}|"
                "affine_pair_divisors={};{};{}".format(
                    p1_index, p2_index, p3_index, sign1, sign2, sign3,
                    *divisors,
                ),
                flush=True,
            )
        if args.max_survivors and len(survivors) >= args.max_survivors:
            break
    if args.max_survivors and len(survivors) >= args.max_survivors:
        break

print(
    "Q80MODP|SUMMARY|prime={}|tested={}|cm={}|with_P1={}|with_P1_P2={}|"
    "status=BOUNDED_EXPERIMENT".format(
        args.prime, tested, cm_candidates, p1_candidates, len(survivors)
    ),
    flush=True,
)
