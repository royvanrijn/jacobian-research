#!/usr/bin/env python3
"""Explore the branches of the fixed-r=6 endpoint eliminant at m=infinity.

This starts from the same exact integral endpoint equations as the reduction
certificate.  Singular forms the residual eliminant H_6 and returns its
coefficients in descending asymptotic orders of m.  Numerical roots only seed
the atlas: disjoint rational Rouche disks certify the root locations, and
rational Taylor bounds certify every limiting modulus separation.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from math import factorial, isqrt

import sympy as sp
from sympy.polys.domains import QQ


m, y, z = sp.symbols("m y z")
R = 6
D = 1 - y


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


integral_E6 = sp.expand(integral_polynomial(E6))
integral_F6 = sp.expand(integral_polynomial(F6))
SINGULAR = shutil.which("Singular")
assert SINGULAR is not None

program = f"""
ring rr=0,(z,y,m),dp;
poly E={singular_expression(integral_E6)};
poly F={singular_expression(integral_F6)};
poly H=resultant(E,F,z)/(y-1)^7;
matrix Cm=coeffs(H,m);
matrix Cy=coeffs(H,y);
"BEGIN_M_COEFFICIENTS";
for (int j=90;j>=0;j--) {{
  if (Cm[j+1,1]<>0) {{ j; Cm[j+1,1]; }}
}}
"END_M_COEFFICIENTS";
"BEGIN_Y_DEGREES";
for (int j=0;j<=29;j++) {{
  deg(poly(Cy[j+1,1]),intvec(0,0,1));
}}
"END_Y_DEGREES";
"""

completed = subprocess.run(
    [SINGULAR, "-q"], input=program, text=True, capture_output=True, check=True
)
assert completed.stderr == "", completed.stderr

match = re.search(
    r"BEGIN_M_COEFFICIENTS\n(.*?)END_M_COEFFICIENTS",
    completed.stdout,
    flags=re.DOTALL,
)
assert match is not None
lines = [line.strip() for line in match.group(1).splitlines() if line.strip()]
m_coefficients: dict[int, sp.Poly] = {}
for index in range(0, len(lines), 2):
    degree = int(lines[index])
    serialized = re.sub(
        r"(?<=\d)y(\d*)",
        lambda match: "*y" + (f"**{match.group(1)}" if match.group(1) else ""),
        lines[index + 1],
    )
    expression = sp.sympify(serialized.replace("^", "**"))
    m_coefficients[degree] = sp.Poly(expression, y, domain=QQ)

degree_match = re.search(
    r"BEGIN_Y_DEGREES\n(.*?)END_Y_DEGREES", completed.stdout, flags=re.DOTALL
)
assert degree_match is not None
y_coefficient_m_degrees = [
    int(line.strip())
    for line in degree_match.group(1).splitlines()
    if re.fullmatch(r"-?\d+", line.strip())
]

assert sorted(m_coefficients) == list(range(91))
assert y_coefficient_m_degrees == [90] * 30
for m_degree in range(61, 91):
    forced_power = m_degree - 61
    divisor = sp.Poly((y - 1) ** forced_power, y, domain=QQ)
    assert sp.rem(m_coefficients[m_degree], divisor) == 0

leading = m_coefficients[max(m_coefficients)]
assert sp.factor(leading.as_expr() / (y - 1) ** 29).is_Rational
print("[m^90]H6 is a nonzero scalar times (y-1)^29")

# The exact staircase [m^(90-k)]H divisible by (y-1)^(29-k) puts every
# branch in the chart y=1+t*x, t=1/m.  The normalized transform is
# P_6(t,x)=t^61 H(1/t,1+t*x).  Build it coefficientwise, avoiding a second
# parse of the full 2724-term eliminant.
t, x = sp.symbols("t x")
P_expression = 0
for m_degree, coefficient in m_coefficients.items():
    transformed = sp.Poly(coefficient.as_expr().subs(y, 1 + t * x), t, x)
    for (t_degree, x_degree), value in transformed.terms():
        total_t_degree = 61 - m_degree + t_degree
        assert total_t_degree >= 0
        P_expression += value * t**total_t_degree * x**x_degree
P6 = sp.Poly(P_expression, x, t, domain=QQ)
P0_raw = sp.Poly(P6.as_expr().subs(t, 0), x, domain=QQ)
assert P6.degree(x) == P0_raw.degree() == 29
P0_content, P0 = P0_raw.primitive()
assert sp.gcd(P0, P0.diff()).degree() == 0
print("P6 bidegree (x,t):", (P6.degree(x), P6.degree(t)))
P0_factors = sp.factor_list(P0.as_expr())[1]
assert len(P0_factors) == 1 and P0_factors[0][1] == 1
print("P0: degree 29, squarefree and irreducible over Q")

# Limiting endpoint equations in this chart.  E6 has t-valuation 6 and F6
# has t-valuation 7.  Their common-root z is generically reconstructed by
# the penultimate (linear) subresultant.
E_chart = sp.cancel(E6.subs({m: 1 / t, y: 1 + t * x}) / t**6)
F_chart = sp.cancel(F6.subs({m: 1 / t, y: 1 + t * x}) / t**7)
E0 = sp.Poly(sp.cancel(E_chart).subs(t, 0), z, x, domain=QQ)
F0 = sp.Poly(sp.cancel(F_chart).subs(t, 0), z, x, domain=QQ)
assert E0.degree(z) == 5 and F0.degree(z) == 6
limit_subresultants = sp.subresultants(E0.as_expr(), F0.as_expr(), z)
limit_linear = sp.Poly(limit_subresultants[-2], z, domain=QQ.poly_ring(x))
assert limit_linear.degree() == 1
limit_A = sp.Poly(limit_linear.coeff_monomial(z), x, domain=QQ)
limit_B = sp.Poly(limit_linear.coeff_monomial(1), x, domain=QQ)
z_limit = sp.cancel(-limit_B.as_expr() / limit_A.as_expr())
limit_resultant = sp.Poly(limit_subresultants[-1], x, domain=QQ)
assert sp.rem(limit_resultant, P0) == 0
assert sp.Poly(limit_resultant.exquo(P0), x).monic().as_expr() == x**7

roots = sp.nroots(P0.as_expr(), n=40, maxsteps=1000)
assert P0.count_roots(-sp.oo, 0) == 0
assert P0.count_roots(0, sp.oo) == 1

# First Puiseux corrections.  If x(t)=x0+x1*t+..., then
# y=1+x0/m+x1/m^2+... .  The t-linear coefficient of P gives x1, and the
# t-linear coefficient of F gives z1 for z=z0+z1*t+... .
P1 = sp.Poly(P6.as_expr().coeff(t, 1), x, domain=QQ)
x1_function = sp.cancel(-P1.as_expr() / P0_raw.diff().as_expr())
F_t0 = sp.cancel(sp.diff(F_chart, t).subs(t, 0))
F0_expression = F0.as_expr()
z1_function = sp.cancel(
    -(
        sp.diff(F0_expression, x) * x1_function + F_t0
    )
    / sp.diff(F0_expression, z)
)

atlas_rows: list[tuple[object, ...]] = []
for index, x_root in enumerate(roots, start=1):
    z_root = sp.N(z_limit.subs(x, x_root), 30)
    x1_root = sp.N(x1_function.subs(x, x_root), 24)
    z1_root = sp.N(z1_function.subs({x: x_root, z: z_root}), 24)
    exponential = sp.N(sp.exp(x_root), 30)
    ratio = sp.N(z_root / exponential, 30)
    atlas_rows.append(
        (
            index,
            x_root,
            x1_root,
            z_root,
            z1_root,
            sp.N(abs(ratio), 18),
            sp.N(sp.arg(ratio), 18),
        )
    )


def rational_complex(value: sp.Expr, digits: int = 30) -> sp.Expr:
    """Round a numerical complex value to a Gaussian rational center."""

    real, imaginary = value.as_real_imag()
    return sp.Rational(str(sp.N(real, digits))) + sp.I * sp.Rational(
        str(sp.N(imaginary, digits))
    )


def rational_abs_bounds(value: sp.Expr, digits: int = 25) -> tuple[sp.Rational, sp.Rational]:
    """Exact decimal lower/upper bounds for the modulus of a Gaussian rational."""

    real, imaginary = map(sp.Rational, value.as_real_imag())
    square = sp.Rational(real**2 + imaginary**2)
    scale = 10**digits
    quotient = (int(square.p) * scale**2) // int(square.q)
    lower_integer = isqrt(quotient)
    lower = sp.Rational(lower_integer, scale)
    upper = sp.Rational(lower_integer + 1, scale)
    assert lower**2 <= square <= upper**2
    return lower, upper


def taylor_coefficients(polynomial: sp.Poly, center: sp.Expr) -> list[sp.Expr]:
    coefficients: list[sp.Expr] = []
    derivative = polynomial
    for order in range(polynomial.degree() + 1):
        coefficients.append(sp.cancel(derivative.eval(center) / factorial(order)))
        derivative = derivative.diff()
    return coefficients


def l1_modulus(value: sp.Expr) -> sp.Rational:
    real, imaginary = map(sp.Rational, value.as_real_imag())
    return abs(real) + abs(imaginary)


def polynomial_modulus_bounds(
    polynomial: sp.Poly, center: sp.Expr, radius: sp.Rational
) -> tuple[sp.Rational, sp.Rational]:
    coefficients = taylor_coefficients(polynomial, center)
    center_lower, center_upper = rational_abs_bounds(coefficients[0])
    variation = sum(
        l1_modulus(coefficient) * radius**order
        for order, coefficient in enumerate(coefficients[1:], start=1)
    )
    assert center_lower > variation
    return center_lower - variation, center_upper + variation


def exp_positive_bounds(q: sp.Rational, terms: int = 220) -> tuple[sp.Rational, sp.Rational]:
    """Exact Taylor/geometric-tail enclosure of exp(q), for q>=0."""

    assert q >= 0 and terms + 2 > q
    term = sp.Rational(1)
    partial = term
    for order in range(1, terms + 1):
        term *= q / order
        partial += term
    first_omitted = term * q / (terms + 1)
    tail_upper = first_omitted / (1 - q / (terms + 2))
    return partial, partial + tail_upper


def exp_bounds(q: sp.Rational) -> tuple[sp.Rational, sp.Rational]:
    if q >= 0:
        return exp_positive_bounds(q)
    positive_lower, positive_upper = exp_positive_bounds(-q)
    return 1 / positive_upper, 1 / positive_lower


# Rigorous disks for all 29 simple x0-roots.  Rouché compares the linear
# Taylor term with the constant and higher-order terms on each circle.  The
# disks are disjoint, so their root counts exhaust deg(P0)=29.
root_radius = sp.Rational(1, 10**24)
centers = [rational_complex(root, 30) for root in roots]
assert all(
    max(
        abs(sp.re(left - right)),
        abs(sp.im(left - right)),
    )
    > 2 * root_radius
    for left_index, left in enumerate(centers)
    for right in centers[left_index + 1 :]
)

z_numerator_expression, z_denominator_expression = sp.fraction(z_limit)
z_numerator = sp.Poly(z_numerator_expression, x, domain=QQ)
z_denominator = sp.Poly(z_denominator_expression, x, domain=QQ)
assert sp.gcd(P0, z_denominator).degree() == 0
modulus_certificates: list[str] = []
for center in centers:
    taylor = taylor_coefficients(P0, center)
    linear_lower = rational_abs_bounds(taylor[1])[0] * root_radius
    remainder_upper = l1_modulus(taylor[0]) + sum(
        l1_modulus(coefficient) * root_radius**order
        for order, coefficient in enumerate(taylor[2:], start=2)
    )
    assert linear_lower > remainder_upper

    numerator_lower, numerator_upper = polynomial_modulus_bounds(
        z_numerator, center, root_radius
    )
    denominator_lower, denominator_upper = polynomial_modulus_bounds(
        z_denominator, center, root_radius
    )
    z_lower = numerator_lower / denominator_upper
    z_upper = numerator_upper / denominator_lower
    exp_lower, _ = exp_bounds(sp.Rational(sp.re(center) - root_radius))
    _, exp_upper = exp_bounds(sp.Rational(sp.re(center) + root_radius))
    if z_upper < exp_lower:
        modulus_certificates.append("<")
    elif z_lower > exp_upper:
        modulus_certificates.append(">")
    else:
        raise AssertionError((center, z_lower, z_upper, exp_lower, exp_upper))

assert len(modulus_certificates) == 29


def short_complex(value: object) -> str:
    numerical = sp.N(value, 12)
    real, imaginary = numerical.as_real_imag()
    if abs(float(imaginary)) < 1e-11:
        return f"{float(real):.9g}"
    return f"{float(real):.9g}{float(imaginary):+.9g}i"


print("class | x0 | x1 | z0 | z1 | |z0/e^x0| | arg | gap")
print("---|---|---|---|---|---:|---:|---")
class_index = 0
for row, sign in zip(atlas_rows, modulus_certificates):
    index, x_root, x1_root, z_root, z1_root, ratio_modulus, ratio_argument = row
    if float(sp.im(x_root)) > 1e-20:
        continue
    class_index += 1
    print(
        f"{class_index} | {short_complex(x_root)} | {short_complex(x1_root)} | "
        f"{short_complex(z_root)} | {short_complex(z1_root)} | "
        f"{float(ratio_modulus):.10g} | {float(ratio_argument):.9g} | {sign}"
    )
assert class_index == 15
print("PASS: 29 disjoint Rouche disks and 29 strict limiting-modulus gaps")
