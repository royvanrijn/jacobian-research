#!/usr/bin/env python3
"""Effective tail certificate for the cancellation contact resultant at r=6.

This checker closes the gap between the asymptotic 29-branch theorem and the
exact finite certificates through m=40.  Put

    t = 1/m,  y = 1+t*x,
    P6(t,x) = t^61 H6(1/t,1+t*x).

For a partition of 0 <= t <= 1/41, numerical continuation supplies candidate
centres for the 29 roots of P6.  The candidates are not trusted: on every
parameter cell, Arb ball arithmetic verifies 29 pairwise-disjoint Rouche
disks, each containing exactly one root for every t in the cell.

At a common endpoint root, the binomial sixth-moment equation implies

    10 = Q(t,x) * (1+t*x)^(6/t+1),

where Q is an explicit polynomial of bidegree (6,6).  On every certified
root tube, interval logarithms separate the two sides either in modulus or
in phase.  The positive real branch is handled by a real logarithmic bound.
Consequently no common root exists for m >= 41.  The separate modular-grid
checker certifies 1 <= m <= 40.

Requirements beyond the base repository environment: python-flint >= 0.9.
All Arb enclosures are rigorous; floating-point root finding only proposes
the disks subsequently certified by Arb.
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

import mpmath as mp
import sympy as sp
from flint import acb, arb, ctx, fmpq
from sympy.polys.domains import QQ


ctx.prec = 320
mp.mp.dps = 90

m, y, z, t, x = sp.symbols("m y z t x")
R = 6
D = 1 - y
TAIL_START = 41
CELL_COUNT = 256
LOG_TERMS = 100


def beta(k: int) -> sp.Expr:
    return sp.factorial(k) / sp.prod(m * k + j for j in range(1, k + 2))


def endpoint_tail(k: int) -> sp.Expr:
    return sum(
        (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
        for j in range(k + 1)
    )


E6 = sum(
    (-1) ** k
    * sp.binomial(R, k)
    * z ** (R - k - 1)
    * (beta(k) - y * z**k * endpoint_tail(k))
    * D ** (R - k - 1)
    for k in range(R)
)
F6 = beta(R) - y * z**R * endpoint_tail(R)


def integral_polynomial(expression: sp.Expr) -> sp.Expr:
    parameter_field = QQ.frac_field(m)
    over_parameter_ring = sp.Poly(
        expression, z, y, domain=parameter_field
    ).clear_denoms(convert=True)[1]
    return sp.Poly(
        over_parameter_ring.as_expr(), z, y, m, domain=QQ
    ).clear_denoms(convert=True)[1].as_expr()


def singular_expression(expression: sp.Expr) -> str:
    return str(expression).replace("**", "^")


def build_chart_polynomial() -> sp.Poly:
    """Construct P6 exactly, coefficientwise, using Singular for H6."""

    singular = shutil.which("Singular")
    assert singular is not None, "the effective r=6 certificate requires Singular"
    integral_E6 = sp.expand(integral_polynomial(E6))
    integral_F6 = sp.expand(integral_polynomial(F6))
    program = f"""
ring rr=0,(z,y,m),dp;
poly E={singular_expression(integral_E6)};
poly F={singular_expression(integral_F6)};
poly H=resultant(E,F,z)/(y-1)^7;
matrix Cm=coeffs(H,m);
"BEGIN";
for (int j=90;j>=0;j--) {{
  if (Cm[j+1,1]<>0) {{ j; Cm[j+1,1]; }}
}}
"END";
"""
    completed = subprocess.run(
        [singular, "-q"], input=program, text=True, capture_output=True, check=True
    )
    assert completed.stderr == "", completed.stderr
    match = re.search(r"BEGIN\n(.*?)END", completed.stdout, flags=re.DOTALL)
    assert match is not None
    lines = [line.strip() for line in match.group(1).splitlines() if line.strip()]
    coefficients: dict[int, sp.Poly] = {}
    for index in range(0, len(lines), 2):
        degree = int(lines[index])
        serialized = re.sub(
            r"(?<=\d)y(\d*)",
            lambda found: "*y"
            + (f"**{found.group(1)}" if found.group(1) else ""),
            lines[index + 1],
        )
        coefficients[degree] = sp.Poly(
            sp.sympify(serialized.replace("^", "**")), y, domain=QQ
        )
    assert sorted(coefficients) == list(range(91))

    expression = 0
    for m_degree, coefficient in coefficients.items():
        transformed = sp.Poly(coefficient.as_expr().subs(y, 1 + t * x), t, x)
        for (t_degree, x_degree), value in transformed.terms():
            total_t_degree = 61 - m_degree + t_degree
            assert total_t_degree >= 0
            expression += value * t**total_t_degree * x**x_degree
    polynomial = sp.Poly(expression, x, t, domain=QQ)
    assert (polynomial.degree(x), polynomial.degree(t)) == (29, 89)
    return polynomial


P6 = build_chart_polynomial()
P0 = sp.Poly(P6.as_expr().subs(t, 0), x, domain=QQ)
assert P0.degree() == 29 and sp.gcd(P0, P0.diff()).degree() == 0

# The common denominator in the normalized F6 coefficients cancels.  This
# produces the regular polynomial Q in 10 = Q*y^(6/t+1).
F_chart = sp.cancel(F6.subs({m: 1 / t, y: 1 + t * x}) / t**7)
F_poly = sp.Poly(F_chart, z)
B = sp.cancel(F_poly.coeff_monomial(1))
C = sp.cancel(-F_poly.coeff_monomial(z**6))
Q_expression = sp.factor(10 * C / (B * (1 + t * x)))
assert sp.denom(Q_expression) == 1
Q_polynomial = sp.Poly(Q_expression, x, t, domain=QQ)
assert (Q_polynomial.degree(x), Q_polynomial.degree(t)) == (6, 6)


def to_fmpq(value: sp.Rational | int) -> fmpq:
    rational = sp.Rational(value)
    return fmpq(int(rational.p), int(rational.q))


def exact_arb(value: sp.Rational | int) -> arb:
    return arb(to_fmpq(value))


def rational_complex(value: mp.mpc, digits: int = 45) -> tuple[sp.Rational, sp.Rational]:
    return (
        sp.Rational(mp.nstr(value.real, digits)),
        sp.Rational(mp.nstr(value.imag, digits)),
    )


def exact_acb(real: sp.Rational, imaginary: sp.Rational = sp.Rational(0)) -> acb:
    return acb(exact_arb(real), exact_arb(imaginary))


def interval_arb(left: sp.Rational, right: sp.Rational) -> arb:
    assert left <= right
    midpoint = (left + right) / 2
    radius = (right - left) / 2
    return arb(to_fmpq(midpoint), to_fmpq(radius))


def coefficient_polynomials(polynomial: sp.Poly, outer: sp.Symbol) -> list[list[fmpq]]:
    """Ascending outer coefficients, each an ascending polynomial in t."""

    result: list[list[fmpq]] = []
    for degree in range(polynomial.degree(outer) + 1):
        member = sp.Poly(polynomial.as_expr().coeff(outer, degree), t, domain=QQ)
        result.append([to_fmpq(member.nth(j)) for j in range(member.degree() + 1)])
    return result


P_COEFFICIENTS = coefficient_polynomials(P6, x)
Q_COEFFICIENTS = coefficient_polynomials(Q_polynomial, x)


def evaluate_univariate(coefficients: list[fmpq], argument: arb | acb) -> arb | acb:
    value: arb | acb = arb(0)
    for coefficient in reversed(coefficients):
        value = value * argument + coefficient
    return value


def evaluate_bivariate(
    coefficients: list[list[fmpq]], t_ball: arb, x_ball: arb | acb
) -> arb | acb:
    value: arb | acb = arb(0)
    for member in reversed(coefficients):
        value = value * x_ball + evaluate_univariate(member, t_ball)
    return value


def numerical_coefficients(parameter: mp.mpf) -> list[mp.mpc]:
    values: list[mp.mpc] = []
    for member in P_COEFFICIENTS:
        value = mp.mpf("0")
        for coefficient in reversed(member):
            value = value * parameter + mp.mpf(int(coefficient.p)) / int(coefficient.q)
        values.append(mp.mpc(value))
    return values


def polynomial_value(coefficients: list[mp.mpc], argument: mp.mpc) -> mp.mpc:
    value = mp.mpc(0)
    for coefficient in reversed(coefficients):
        value = value * argument + coefficient
    return value


def refine_root(coefficients: list[mp.mpc], seed: mp.mpc) -> mp.mpc:
    root = seed
    derivative = [j * coefficients[j] for j in range(1, len(coefficients))]
    final_correction = mp.inf
    for _ in range(30):
        correction = polynomial_value(coefficients, root) / polynomial_value(
            derivative, root
        )
        root -= correction
        final_correction = abs(correction)
        if abs(correction) < mp.mpf("1e-75"):
            break
    assert final_correction < mp.mpf("1e-60")
    return root


initial_roots = [
    mp.mpc(str(sp.re(root)), str(sp.im(root)))
    for root in sp.nroots(P0.as_expr(), n=80)
]
initial_roots.sort(key=lambda root: (float(root.real), float(root.imag)))


@dataclass(frozen=True)
class Tube:
    center_real: sp.Rational
    center_imaginary: sp.Rational
    radius: sp.Rational
    is_positive_real: bool

    @property
    def center(self) -> acb:
        return exact_acb(self.center_real, self.center_imaginary)


def derivative_balls(p_values: list[arb], center: acb) -> list[acb]:
    """Taylor coefficients P^(j)(center)/j!, with t still a ball."""

    degree = len(p_values) - 1
    powers = [acb(1)]
    for _ in range(degree):
        powers.append(powers[-1] * center)
    result: list[acb] = []
    for j in range(degree + 1):
        value = acb(0)
        for i in range(j, degree + 1):
            value += math.comb(i, j) * p_values[i] * powers[i - j]
        result.append(value)
    return result


def exact_low_derivative_ball(order: int, center: acb, t_ball: arb) -> acb:
    """Evaluate the order-th x Taylor coefficient before intervalizing t.

    Forming the exact coefficient polynomial first preserves the substantial
    cancellation in P(t,center) that a rectangular evaluation of the x
    coefficients would otherwise lose.
    """

    t_degree = P6.degree(t)
    coefficients = [acb(0) for _ in range(t_degree + 1)]
    for i in range(order, len(P_COEFFICIENTS)):
        multiplier = math.comb(i, order) * center ** (i - order)
        for j, coefficient in enumerate(P_COEFFICIENTS[i]):
            coefficients[j] += multiplier * coefficient
    return evaluate_univariate(coefficients, t_ball)


def certify_rouche(taylor: list[acb], radius: sp.Rational) -> bool:
    r_ball = exact_arb(radius)
    left = taylor[1].abs_lower() * r_ball
    right = taylor[0].abs_upper()
    power = r_ball * r_ball
    for coefficient in taylor[2:]:
        right += coefficient.abs_upper() * power
        power *= r_ball
    return bool(left > right)


def choose_radius(taylor: list[acb], separation: mp.mpf) -> sp.Rational:
    base = sp.Rational(1, 10**18)
    upper = sp.Rational(mp.nstr(separation * mp.mpf("0.42"), 35))
    radius = base
    while radius < upper:
        if certify_rouche(taylor, radius):
            return radius
        radius *= 2
    raise AssertionError("no Rouche radius found; increase CELL_COUNT")


def regular_log_ratio(t_ball: arb, x_ball: acb, terms: int = LOG_TERMS) -> acb:
    """Enclose log(1+t*x)/t by its regular power series at t=0."""

    u = t_ball * x_ball
    u_upper = u.abs_upper()
    assert bool(u_upper < 1)
    total = acb(0)
    power = acb(1)
    for order in range(terms):
        total += power / (order + 1)
        power *= -u
    result = x_ball * total
    x_upper = x_ball.abs_upper()
    tail = x_upper * (u_upper**terms) / ((terms + 1) * (1 - u_upper))
    error = acb(arb(0, tail), arb(0, tail))
    return result + error


def excludes_zero(value: arb) -> bool:
    return bool(value.lower() > 0 or value.upper() < 0)


def phase_excludes_multiples_of_two_pi(value: arb) -> bool:
    lower = float(value.lower())
    upper = float(value.upper())
    first = math.floor(lower / (2 * math.pi)) - 1
    last = math.ceil(upper / (2 * math.pi)) + 1
    two_pi = 2 * arb.pi()
    return all(excludes_zero(value - k * two_pi) for k in range(first, last + 1))


def certify_separation(t_ball: arb, tube: Tube) -> str:
    radius = exact_arb(tube.radius)
    if tube.is_positive_real:
        # The disk is invariant under conjugation and contains one root, so
        # that root is real.  Restrict to its real interval here.
        x_ball = interval_arb(
            tube.center_real - tube.radius, tube.center_real + tube.radius
        )
        q_ball = evaluate_bivariate(Q_COEFFICIENTS, t_ball, x_ball)
        assert isinstance(q_ball, arb) and bool(q_ball.lower() > 0)
        y_ball = 1 + t_ball * x_ball
        assert bool(y_ball.lower() > 0)
        if bool(t_ball.lower() > 0):
            difference = (q_ball / 10).log() + (6 / t_ball + 1) * y_ball.log()
        else:
            logarithm_ratio = regular_log_ratio(t_ball, acb(x_ball))
            difference = (q_ball / 10).log() + (6 + t_ball) * logarithm_ratio.real
        assert bool(difference.lower() > 0)
        return "modulus"

    x_error = acb(arb(0, radius), arb(0, radius))
    x_ball = tube.center + x_error
    q_ball = evaluate_bivariate(Q_COEFFICIENTS, t_ball, x_ball)
    assert isinstance(q_ball, acb) and not q_ball.contains(0)
    logarithm_ratio = regular_log_ratio(t_ball, x_ball)
    difference = (q_ball / 10).log() + (6 + t_ball) * logarithm_ratio
    if excludes_zero(difference.real):
        return "modulus"
    assert phase_excludes_multiples_of_two_pi(difference.imag)
    return "phase"


previous_roots = initial_roots
modulus_count = 0
phase_count = 0
tube_count = 0
max_radius = sp.Rational(0)

for cell in range(CELL_COUNT):
    left = sp.Rational(cell, TAIL_START * CELL_COUNT)
    right = sp.Rational(cell + 1, TAIL_START * CELL_COUNT)
    midpoint = (left + right) / 2
    numerical_t = mp.mpf(int(midpoint.p)) / int(midpoint.q)
    numerical = numerical_coefficients(numerical_t)
    roots = [refine_root(numerical, seed) for seed in previous_roots]
    roots.sort(key=lambda root: (float(root.real), float(root.imag)))
    previous_roots = roots

    t_ball = interval_arb(left, right)
    p_values = [evaluate_univariate(member, t_ball) for member in P_COEFFICIENTS]
    tubes: list[Tube] = []
    for index, root in enumerate(roots):
        real, imaginary = rational_complex(root)
        is_positive_real = abs(root.imag) < mp.mpf("1e-60") and root.real > 0
        if is_positive_real:
            imaginary = sp.Rational(0)
        center = exact_acb(real, imaginary)
        taylor = derivative_balls(p_values, center)
        taylor[0] = exact_low_derivative_ball(0, center, t_ball)
        taylor[1] = exact_low_derivative_ball(1, center, t_ball)
        separation = min(
            abs(root - other) for other_index, other in enumerate(roots) if other_index != index
        )
        try:
            radius = choose_radius(taylor, separation)
        except AssertionError as error:
            raise AssertionError(
                f"Rouche failure in cell {cell}, root {index}, center {root}"
            ) from error
        tubes.append(Tube(real, imaginary, radius, is_positive_real))
        max_radius = max(max_radius, radius)

    assert sum(tube.is_positive_real for tube in tubes) == 1
    for left_index, first in enumerate(tubes):
        for second in tubes[left_index + 1 :]:
            distance = (first.center - second.center).abs_lower()
            assert bool(distance > exact_arb(first.radius + second.radius))

    for tube in tubes:
        method = certify_separation(t_ball, tube)
        modulus_count += method == "modulus"
        phase_count += method == "phase"
        tube_count += 1

assert tube_count == 29 * CELL_COUNT
print(
    "PASS contact resultant r=6 effective tail: "
    f"{tube_count} rigorous Rouche tubes on 0<=t<=1/{TAIL_START}"
)
print(
    "PASS contact resultant r=6 effective tail: "
    f"{modulus_count} modulus and {phase_count} phase separations"
)
print(
    "PASS contact resultant r=6 effective tail: "
    f"nonvanishing for every integer m>={TAIL_START}"
)

# Join the effective tail to an exact finite certificate.  Degree preservation
# and a modular gcd equal to one prove coprimality over QQ.
PRIME = 1_000_003


def inverse_mod(value: int) -> int:
    residue = value % PRIME
    assert residue != 0
    return pow(residue, -1, PRIME)


def endpoint_moment_coefficients(parameter: int, k: int) -> list[int]:
    denominator = math.prod(parameter * k + j for j in range(1, k + 2))
    normalized_beta = math.factorial(k) * inverse_mod(denominator) % PRIME
    return [
        normalized_beta * math.comb(j + k, k) % PRIME
        for j in range(parameter * k + 1)
    ]


def finite_endpoint_pair(parameter: int) -> tuple[sp.Poly, sp.Poly]:
    degree = parameter * R
    moments = [endpoint_moment_coefficients(parameter, k) for k in range(R + 1)]
    k_coefficients = moments[R]
    l_coefficients = [0] * (degree + 1)
    for k, moment in enumerate(moments):
        shift = parameter * (R - k)
        multiplier = (-1) ** k * math.comb(R, k)
        for index, coefficient in enumerate(moment):
            l_coefficients[index + shift] = (
                l_coefficients[index + shift] + multiplier * coefficient
            ) % PRIME
    K = sp.Poly.from_list(k_coefficients[::-1], gens=y, modulus=PRIME)
    L = sp.Poly.from_list(l_coefficients[::-1], gens=y, modulus=PRIME)
    assert K.degree() == L.degree() == degree
    assert K.LC() % PRIME != 0 and L.LC() % PRIME != 0
    return K, L


for finite_m in range(1, TAIL_START):
    finite_K, finite_L = finite_endpoint_pair(finite_m)
    assert sp.gcd(finite_K, finite_L).degree() == 0

print(
    "PASS contact resultant r=6 finite range: "
    f"40 modular gcd certificates modulo {PRIME}"
)
print("PASS contact resultant: uniform nonvanishing for every m>=1 and r=6")
