#!/usr/bin/env python3
"""Exact fixed-r=5 reduction for the cancellation contact resultant.

This checker constructs the bounded-degree endpoint equations without forming
the degree-5*m polynomials in w.  It uses the sparse quartic--binomial-quintic
resultant and subresultant formulas to obtain the primitive endpoint eliminant
H_5(m,y) and the rational reconstruction z=-B/A.

The final uniform separation certificate is developed below this algebraic
core.  No floating-point computation is used for the asserted identities or
the exact small-m gcd checks.
"""

from __future__ import annotations

import math
import re
import shutil
import subprocess

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.rings import ring


m, y, z = sp.symbols("m y z")
D = 1 - y


def beta(k: int) -> sp.Expr:
    return sp.factorial(k) / sp.prod(m * k + j for j in range(1, k + 2))


def endpoint_tail(k: int) -> sp.Expr:
    return sum(
        (-1) ** j * sp.binomial(k, j) * y**j / (m * k + j + 1)
        for j in range(k + 1)
    )


# At a common K_{m,5},L_{m,5} root, M_5=0 and z != 0.  Dividing the
# triangular identity by z gives a quartic E_5; M_5 gives the binomial
# quintic F_5.  Writing the coefficients directly preserves their sparsity.
a = (
    D**5
    + 5 * y * D**3 * endpoint_tail(1)
    - 10 * y * D**2 * endpoint_tail(2)
    + 10 * y * D * endpoint_tail(3)
    - 5 * y * endpoint_tail(4)
)
b = -5 * beta(1) * D**3
c = 10 * beta(2) * D**2
d = -10 * beta(3) * D
e = 5 * beta(4)
p = -y * endpoint_tail(5)
q = beta(5)

E5 = a * z**4 + b * z**3 + c * z**2 + d * z + e
F5 = p * z**5 + q

assert sp.degree(E5, z) == 4
assert sp.degree(F5, z) == 5


# Clear parameter denominators once.  Scaling either endpoint equation by a
# nonzero function of m changes neither its roots nor the reconstructed z.
parameter_field = QQ.frac_field(m)
_, integral_E5 = sp.Poly(E5, z, y, domain=parameter_field).clear_denoms(
    convert=True
)
_, integral_F5 = sp.Poly(F5, z, y, domain=parameter_field).clear_denoms(
    convert=True
)

aa, bb, cc, dd, ee = sp.Poly(integral_E5.as_expr(), z).all_coeffs()
pp = sp.Poly(integral_F5.as_expr(), z).coeff_monomial(z**5)
qq = sp.Poly(integral_F5.as_expr(), z).coeff_monomial(1)


# Universal sparse resultant for a quartic and p*z^5+q.  It is also the
# constant member of their subresultant recurrence.  Keeping this 26-term
# norm formula explicit avoids a prohibitively expensive generic PRS over
# QQ(m,y).
def sparse_resultant(
    a0: object,
    b0: object,
    c0: object,
    d0: object,
    e0: object,
    p0: object,
    q0: object,
) -> object:
    return (
        a0**5 * q0**4
        + 5 * a0**3 * b0 * e0 * p0 * q0**3
        + 5 * a0**3 * c0 * d0 * p0 * q0**3
        - 5 * a0**2 * b0**2 * d0 * p0 * q0**3
        - 5 * a0**2 * b0 * c0**2 * p0 * q0**3
        + 5 * a0**2 * c0 * e0**2 * p0**2 * q0**2
        + 5 * a0**2 * d0**2 * e0 * p0**2 * q0**2
        + 5 * a0 * b0**3 * c0 * p0 * q0**3
        + 5 * a0 * b0**2 * e0**2 * p0**2 * q0**2
        - 5 * a0 * b0 * c0 * d0 * e0 * p0**2 * q0**2
        - 5 * a0 * b0 * d0**3 * p0**2 * q0**2
        - 5 * a0 * c0**3 * e0 * p0**2 * q0**2
        + 5 * a0 * c0**2 * d0**2 * p0**2 * q0**2
        + 5 * a0 * d0 * e0**3 * p0**3 * q0
        - b0**5 * p0 * q0**3
        - 5 * b0**3 * d0 * e0 * p0**2 * q0**2
        + 5 * b0**2 * c0**2 * e0 * p0**2 * q0**2
        + 5 * b0**2 * c0 * d0**2 * p0**2 * q0**2
        - 5 * b0 * c0**3 * d0 * p0**2 * q0**2
        + 5 * b0 * c0 * e0**3 * p0**3 * q0
        - 5 * b0 * d0**2 * e0**2 * p0**3 * q0
        + c0**5 * p0**2 * q0**2
        - 5 * c0**2 * d0 * e0**2 * p0**3 * q0
        + 5 * c0 * d0**3 * e0 * p0**3 * q0
        - d0**5 * p0**3 * q0
        + e0**5 * p0**4
    )


# Verify the compact formula once over the universal coefficient ring.
A0, B0, C0, D0, E0, P0, Q0 = sp.symbols("A0 B0 C0 D0 E0 P0 Q0")
assert sp.expand(
    sp.resultant(
        A0 * z**4 + B0 * z**3 + C0 * z**2 + D0 * z + E0,
        P0 * z**5 + Q0,
        z,
    )
    - sparse_resultant(A0, B0, C0, D0, E0, P0, Q0)
) == 0


# Polynomial-ring arithmetic is much faster here than expanding a SymPy
# expression with 49 parameter degrees.  The endpoint factor is exactly
# (y-1)^5 and the remaining primitive eliminant has degree 20 in y.
coefficient_ring, ring_m, ring_y = ring("m,y", QQ)
ring_coefficients = [
    coefficient_ring.from_expr(value) for value in (aa, bb, cc, dd, ee, pp, qq)
]
endpoint_resultant5 = sparse_resultant(*ring_coefficients)
for _ in range(5):
    endpoint_resultant5 = endpoint_resultant5.exquo(ring_y - 1)
assert endpoint_resultant5.rem(ring_y - 1) != 0

eliminant_with_content = sp.Poly(
    endpoint_resultant5.as_expr(), y, domain=QQ.poly_ring(m)
)
eliminant_content, H5 = eliminant_with_content.primitive()
assert sp.factor(eliminant_content - (m + 1)) == 0
assert H5.degree() == 20
assert sp.degree(H5.as_expr(), m) == 49
assert H5.eval(1) != 0


# The penultimate subresultant is A_5(m,y) z+B_5(m,y).  These are the two
# coefficient formulas supplied by the same sparse PRS.
linear_A5 = (
    -aa**3 * cc * pp * qq**2
    + aa**2 * bb**2 * pp * qq**2
    - 2 * aa**2 * dd * ee * pp**2 * qq
    + aa * bb * cc * ee * pp**2 * qq
    + 3 * aa * bb * dd**2 * pp**2 * qq
    - 2 * aa * cc**2 * dd * pp**2 * qq
    - aa * ee**3 * pp**3
    + bb**3 * ee * pp**2 * qq
    - 2 * bb**2 * cc * dd * pp**2 * qq
    + bb * cc**3 * pp**2 * qq
    + 2 * bb * dd * ee**2 * pp**3
    + cc**2 * ee**2 * pp**3
    - 3 * cc * dd**2 * ee * pp**3
    + dd**4 * pp**3
)
linear_B5 = (
    aa**3 * dd * pp * qq**2
    - 2 * aa**2 * bb * cc * pp * qq**2
    + aa**2 * ee**2 * pp**2 * qq
    + aa * bb**3 * pp * qq**2
    - aa * bb * dd * ee * pp**2 * qq
    - 3 * aa * cc**2 * ee * pp**2 * qq
    + 2 * aa * cc * dd**2 * pp**2 * qq
    + 2 * bb**2 * cc * ee * pp**2 * qq
    + bb**2 * dd**2 * pp**2 * qq
    - 3 * bb * cc**2 * dd * pp**2 * qq
    + bb * ee**3 * pp**3
    + cc**4 * pp**2 * qq
    - 2 * cc * dd * ee**2 * pp**3
    + dd**3 * ee * pp**3
)

universal_subresultants = sp.subresultants(
    A0 * z**4 + B0 * z**3 + C0 * z**2 + D0 * z + E0,
    P0 * z**5 + Q0,
    z,
)
assert [sp.degree(member, z) for member in universal_subresultants] == [
    5,
    4,
    3,
    2,
    1,
    0,
]
assert sp.expand(
    universal_subresultants[-2].subs(
        {A0: aa, B0: bb, C0: cc, D0: dd, E0: ee, P0: pp, Q0: qq}
    )
    - (linear_A5 * z + linear_B5)
) == 0


# Exact finite front end.  Directly imposing z=y^m leaves only the excluded
# clearing-denominator root y=1.  The threshold is chosen so that the uniform
# tail certificate can use t=1/m <= 1/20.
excluded_endpoint = sp.Poly((y - 1) ** 5, y).monic()
for m_value in range(1, 20):
    substitutions = {m: m_value, z: y**m_value}
    e_value = sp.Poly(integral_E5.as_expr().subs(substitutions), y)
    f_value = sp.Poly(integral_F5.as_expr().subs(substitutions), y)
    assert sp.gcd(e_value, f_value).monic() == excluded_endpoint


# Uniform tail.  The scale x=m(y-1), t=1/m resolves the twenty H_5 roots at
# t=0.  Since H_5 has m-degree 49 and a t-valuation 20 after substitution,
# P=t^29 H_5(1/t,1+tx) is a genuine polynomial of bidegree (20,49).
t, x, s, v = sp.symbols("t x s v", real=True)
P5 = sp.Poly(
    sp.expand(t**29 * H5.as_expr().subs({m: 1 / t, y: 1 + t * x})),
    x,
    t,
    domain=QQ,
)
assert P5.degree_list() == (20, 49)
assert len(P5.terms()) == 840


def exact_routh_signs(polynomial: sp.Poly) -> list[int]:
    """First-column signs in the exact Routh table."""

    coefficients = polynomial.all_coeffs()
    columns = (len(coefficients) + 1) // 2
    rows = [
        coefficients[0::2] + [0] * (columns - len(coefficients[0::2])),
        coefficients[1::2] + [0] * (columns - len(coefficients[1::2])),
    ]
    for _ in range(2, len(coefficients)):
        row_above_two, row_above = rows[-2], rows[-1]
        assert row_above[0] != 0
        row = [
            sp.cancel(
                (
                    row_above[0] * row_above_two[column + 1]
                    - row_above_two[0] * row_above[column + 1]
                )
                / row_above[0]
            )
            if column + 1 < columns
            else sp.Rational(0)
            for column in range(columns)
        ]
        rows.append(row)
    signs = [int(sp.sign(row[0])) for row in rows]
    assert all(sign != 0 for sign in signs)
    return signs


# At t=0 the exceptional pair is the only pair to the left of Re(x)=-4/5.
# Replacing x by -4/5-s turns that into a right-half-plane count; the exact
# Routh table has two sign changes.
limit_left_test = sp.Poly(
    P5.as_expr().subs({t: 0, x: -sp.Rational(4, 5) - s}), s, domain=QQ
)
routh_signs = exact_routh_signs(limit_left_test)
assert sum(
    routh_signs[index] != routh_signs[index - 1]
    for index in range(1, len(routh_signs))
) == 2


def split_on_line(
    polynomial: sp.Poly, *, axis: str, value: sp.Rational
) -> tuple[sp.Poly, sp.Poly]:
    """Real and imaginary parts on x=value+i*v or x=v+i*value."""

    real_part = sp.Rational(0)
    imaginary_part = sp.Rational(0)
    for (x_power, t_power), coefficient in polynomial.terms():
        for imaginary_power in range(x_power + 1):
            if axis == "vertical":
                term = (
                    coefficient
                    * sp.binomial(x_power, imaginary_power)
                    * value ** (x_power - imaginary_power)
                    * v**imaginary_power
                    * t**t_power
                )
            else:
                term = (
                    coefficient
                    * sp.binomial(x_power, imaginary_power)
                    * v ** (x_power - imaginary_power)
                    * value**imaginary_power
                    * t**t_power
                )
            if imaginary_power % 2 == 0:
                real_part += (-1) ** (imaginary_power // 2) * term
            else:
                imaginary_part += (
                    (-1) ** ((imaginary_power - 1) // 2) * term
                )
    real_poly = sp.Poly(real_part, v, t, domain=QQ).clear_denoms(convert=True)[1]
    imaginary_poly = sp.Poly(
        imaginary_part, v, t, domain=QQ
    ).clear_denoms(convert=True)[1]
    return real_poly, imaginary_poly


def singular_polynomial(polynomial: sp.Poly, variables: tuple[sp.Symbol, ...]) -> str:
    """Serialize an integral polynomial without Singular's a^2/5 ambiguity."""

    pieces: list[str] = []
    for monomial, coefficient in polynomial.terms():
        piece = str(int(coefficient))
        for variable, power in zip(variables, monomial):
            if power:
                piece += f"*{variable}^{power}"
        pieces.append(piece)
    return "+".join(pieces).replace("+-", "-") or "0"


SINGULAR = shutil.which("Singular")
assert SINGULAR is not None, "the exact r=5 certificate requires Singular"


def boundary_resultant_certificate(
    real_part: sp.Poly, imaginary_part: sp.Poly
) -> tuple[int, int, int, int]:
    """Resultant degree and Bernstein signs on 0<=t<=1/20.

    For R(t), the coefficients of
        (1+T)^deg(R) R(T/(20(1+T)))
    are positive binomial multiples of its Bernstein coefficients.  A strict
    common sign proves that the boundary resultant never vanishes.
    """

    program = f"""
ring certificate_ring=0,(v,t),dp;
poly real_part={singular_polynomial(real_part, (v, t))};
poly imaginary_part={singular_polynomial(imaginary_part, (v, t))};
poly boundary_resultant=resultant(real_part,imaginary_part,v);
int resultant_degree=deg(boundary_resultant);
matrix power_coefficients=coeffs(boundary_resultant,t);
poly bernstein_transform=0;
poly remaining_power=(1+t)^resultant_degree;
number scale=1;
for (int j=0;j<=resultant_degree;j++) {{
  bernstein_transform=bernstein_transform
    +number(power_coefficients[j+1,1])*t^j*remaining_power/scale;
  remaining_power=remaining_power/(1+t);
  scale=scale*20;
}}
matrix bernstein_coefficients=coeffs(bernstein_transform,t);
int positive_count=0;
int negative_count=0;
int zero_count=0;
for (int k=1;k<=resultant_degree+1;k++) {{
  if (number(bernstein_coefficients[k,1])>0) {{ positive_count++; }}
  else {{
    if (number(bernstein_coefficients[k,1])<0) {{ negative_count++; }}
    else {{ zero_count++; }}
  }}
}}
resultant_degree; positive_count; negative_count; zero_count;
"""
    completed = subprocess.run(
        [SINGULAR, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert completed.stderr == ""
    values = [int(value) for value in re.findall(r"-?\d+", completed.stdout)]
    assert len(values) == 4
    return tuple(values)  # type: ignore[return-value]


def univariate_bernstein_signs(polynomial: sp.Poly) -> tuple[int, int, int]:
    """Bernstein sign counts for a QQ[t] polynomial on 0<=t<=1/20."""

    degree = polynomial.degree()
    power_coefficients = [sp.Rational(0)] * (degree + 1)
    for (power,), coefficient in polynomial.terms():
        power_coefficients[power] = coefficient
    coefficients = [
        sum(
            sp.Rational(math.comb(index, power), math.comb(degree, power))
            * power_coefficients[power]
            / 20**power
            for power in range(index + 1)
        )
        for index in range(degree + 1)
    ]
    return (
        sum(1 for coefficient in coefficients if coefficient > 0),
        sum(1 for coefficient in coefficients if coefficient < 0),
        sum(1 for coefficient in coefficients if coefficient == 0),
    )


def compress_vertical_frequency(
    real_part: sp.Poly, imaginary_part: sp.Poly
) -> tuple[sp.Poly, sp.Poly]:
    """Replace even powers v^(2j) and odd powers v^(2j+1) by v^j."""

    compressed_real = sp.Rational(0)
    compressed_imaginary = sp.Rational(0)
    for (frequency_power, t_power), coefficient in real_part.terms():
        assert frequency_power % 2 == 0
        compressed_real += coefficient * v ** (frequency_power // 2) * t**t_power
    for (frequency_power, t_power), coefficient in imaginary_part.terms():
        assert frequency_power % 2 == 1
        compressed_imaginary += (
            coefficient * v ** ((frequency_power - 1) // 2) * t**t_power
        )
    return (
        sp.Poly(compressed_real, v, t),
        sp.Poly(compressed_imaginary, v, t),
    )


# No root crosses Re(x)=-4/5.  The degree-931 crossing resultant has 932
# strictly negative Bernstein coefficients, so the Routh count remains two
# for every 0<=t<=1/20.
half_plane_parts = split_on_line(
    P5, axis="vertical", value=-sp.Rational(4, 5)
)
half_plane_zero_frequency = sp.Poly(
    half_plane_parts[0].as_expr().subs(v, 0), t, domain=QQ
)
assert univariate_bernstein_signs(half_plane_zero_frequency) in {
    (half_plane_zero_frequency.degree() + 1, 0, 0),
    (0, half_plane_zero_frequency.degree() + 1, 0),
}
assert boundary_resultant_certificate(
    *compress_vertical_frequency(*half_plane_parts)
) == (931, 0, 932, 0)


# The K-root disk is contained in this left half-plane for m>=20.  Indeed,
# if Re(x)>=-4/5 then |m+x|>=m-4/5>m^2/(m+1).
assert sp.factor(
    m - sp.Rational(4, 5) - m**2 / (m + 1)
    - (m - 4) / (5 * (m + 1))
) == 0


# Trap the upper exceptional root in a fixed rectangle.  At t=0 the exact
# complex root count in the rectangle is one.  A root could leave only by
# meeting one of its four lines, but all four degree-1911 boundary resultants
# have a strict one-sign Bernstein certificate.
rectangle_left = -sp.Rational(13, 4)
rectangle_right = -sp.Rational(29, 10)
rectangle_bottom = sp.Rational(9, 5)
rectangle_top = sp.Rational(11, 5)
limit_polynomial = sp.Poly(P5.as_expr().subs(t, 0), x, domain=QQ)
assert limit_polynomial.count_roots(
    rectangle_left + sp.I * rectangle_bottom,
    rectangle_right + sp.I * rectangle_top,
) == 1

rectangle_certificates = [
    boundary_resultant_certificate(
        *split_on_line(P5, axis="vertical", value=rectangle_left)
    ),
    boundary_resultant_certificate(
        *split_on_line(P5, axis="vertical", value=rectangle_right)
    ),
    boundary_resultant_certificate(
        *split_on_line(P5, axis="horizontal", value=rectangle_bottom)
    ),
    boundary_resultant_certificate(
        *split_on_line(P5, axis="horizontal", value=rectangle_top)
    ),
]
assert rectangle_certificates == [
    (1911, 0, 1912, 0),
    (1911, 0, 1912, 0),
    (1911, 0, 1912, 0),
    (1911, 1912, 0, 0),
]

# The whole rectangle lies strictly inside the transformed K-root disk.  The
# worst corner has Re(x)=-29/10 and |Im(x)|=11/5; after m=n+20 the numerator
# of the squared-radius difference has strictly positive coefficients.
tail_parameter = sp.symbols("tail_parameter", nonnegative=True)
rectangle_disk_difference = sp.together(
    m**4 / (m + 1) ** 2
    - ((m - sp.Rational(29, 10)) ** 2 + sp.Rational(11, 5) ** 2)
).as_numer_denom()[0]
rectangle_disk_certificate = sp.Poly(
    sp.expand(rectangle_disk_difference.subs(m, tail_parameter + 20)),
    tail_parameter,
)
assert rectangle_disk_certificate.all_coeffs() == [
    380,
    22535,
    443330,
    2891275,
]


def split_real_imaginary(
    polynomial: sp.Expr,
) -> tuple[sp.Poly, sp.Poly]:
    """Evaluate x=a+i*b by exact even/odd binomial splitting."""

    real_a, real_b = sp.symbols("real_a real_b", real=True)
    real_part = sp.Rational(0)
    imaginary_part = sp.Rational(0)
    for (x_power, t_power), coefficient in sp.Poly(
        polynomial, x, t, domain=QQ
    ).terms():
        for imaginary_power in range(x_power + 1):
            term = (
                coefficient
                * sp.binomial(x_power, imaginary_power)
                * real_a ** (x_power - imaginary_power)
                * real_b**imaginary_power
                * t**t_power
            )
            if imaginary_power % 2 == 0:
                real_part += (-1) ** (imaginary_power // 2) * term
            else:
                imaginary_part += (
                    (-1) ** ((imaginary_power - 1) // 2) * term
                )
    return (
        sp.Poly(real_part, real_a, real_b, t, domain=QQ),
        sp.Poly(imaginary_part, real_a, real_b, t, domain=QQ),
    )


def tensor_bernstein_certificate(polynomial: sp.Poly) -> tuple[int, int, int]:
    """Count tensor Bernstein signs on the exceptional rectangle and t-box."""

    integral_polynomial = polynomial.clear_denoms(convert=True)[1]
    real_a, real_b = polynomial.gens[:2]
    degree_a, degree_b, degree_t = polynomial.degree_list()

    def scaled_endpoints(left: sp.Rational, right: sp.Rational) -> tuple[int, int, int]:
        denominator = math.lcm(int(left.q), int(right.q))
        return int(left * denominator), int(right * denominator), denominator

    left_a, right_a, denominator_a = scaled_endpoints(
        rectangle_left, rectangle_right
    )
    left_b, right_b, denominator_b = scaled_endpoints(
        rectangle_bottom, rectangle_top
    )
    left_t, right_t, denominator_t = scaled_endpoints(
        sp.Rational(0), sp.Rational(1, 20)
    )
    serialized = singular_polynomial(
        integral_polynomial, (real_a, real_b, t)
    )
    program = f"""
ring box_ring=0,(real_a,real_b,t,A,B,T),dp;
poly box_polynomial={serialized};
matrix axis_coefficients=coeffs(box_polynomial,real_a);
poly transformed=0;
for (int i=0;i<={degree_a};i++) {{
  transformed=transformed+axis_coefficients[i+1,1]
    *({left_a}+{right_a}*A)^i
    *({denominator_a}+{denominator_a}*A)^({degree_a}-i);
}}
box_polynomial=transformed;
axis_coefficients=coeffs(box_polynomial,real_b);
transformed=0;
for (i=0;i<={degree_b};i++) {{
  transformed=transformed+axis_coefficients[i+1,1]
    *({left_b}+{right_b}*B)^i
    *({denominator_b}+{denominator_b}*B)^({degree_b}-i);
}}
box_polynomial=transformed;
axis_coefficients=coeffs(box_polynomial,t);
transformed=0;
for (i=0;i<={degree_t};i++) {{
  transformed=transformed+axis_coefficients[i+1,1]
    *({left_t}+{right_t}*T)^i
    *({denominator_t}+{denominator_t}*T)^({degree_t}-i);
}}
box_polynomial=transformed;
int positive_count=0;
int negative_count=0;
int zero_count=0;
while (box_polynomial!=0) {{
  number leading_coefficient=leadcoef(box_polynomial);
  if (leading_coefficient>0) {{ positive_count++; }}
  else {{
    if (leading_coefficient<0) {{ negative_count++; }}
    else {{ zero_count++; }}
  }}
  box_polynomial=box_polynomial-lead(box_polynomial);
}}
positive_count; negative_count; zero_count;
"""
    completed = subprocess.run(
        [SINGULAR, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert completed.stderr == ""
    values = [int(value) for value in re.findall(r"-?\d+", completed.stdout)]
    assert len(values) == 3
    return tuple(values)  # type: ignore[return-value]


# Every F_5-root has the same modulus, whether or not the linear
# subresultant degenerates:
#
#       z^5 = beta_5/(y T_5(y)).
#
# On the exceptional box prove |z|>1/9 by checking
# 9^10 |Numerator|^2-|Denominator|^2>0.  This is substantially smaller than
# a Bernstein certificate for the reconstructed rational branch.
fifth_power_branch = sp.cancel(
    (beta(5) / (y * endpoint_tail(5))).subs({m: 1 / t, y: 1 + t * x})
)
fifth_power_numerator, fifth_power_denominator = sp.together(
    fifth_power_branch
).as_numer_denom()
numerator_real, numerator_imaginary = split_real_imaginary(
    fifth_power_numerator
)
denominator_real, denominator_imaginary = split_real_imaginary(
    fifth_power_denominator
)
modulus_difference = (
    9**10 * (numerator_real**2 + numerator_imaginary**2)
    - (denominator_real**2 + denominator_imaginary**2)
)
assert all(degree >= 0 for degree in modulus_difference.degree_list())
modulus_signs = tensor_bernstein_certificate(modulus_difference)
assert modulus_signs[0] > 0 and modulus_signs[1:] == (0, 0)


# The same rectangle makes y^m strictly smaller than 1/12.  From
# log(1+s)<=s,
#
# log |1+x/m|^m <= Re(x)+(|x|^2)/(2m) < -5/2.
#
# The displayed finite Taylor sum proves exp(5/2)>12 exactly.
log_modulus_bound = rectangle_right + sp.Rational(1, 40) * (
    rectangle_left**2 + rectangle_top**2
)
assert log_modulus_bound < -sp.Rational(5, 2)
exp_lower_bound = sum(
    sp.Rational(5, 2) ** power / sp.factorial(power) for power in range(8)
)
assert exp_lower_bound > 12


print("PASS contact resultant r=5: sparse eliminant has bidegree (49,20)")
print("PASS contact resultant r=5: explicit degree 5,4,3,2,1,0 subresultant recurrence")
print("PASS contact resultant r=5: exact endpoint gcds for 1<=m<=19")
print("PASS contact resultant r=5: exactly two K-disk eliminant roots for m>=20")
print("PASS contact resultant r=5: exceptional-box modulus separation")
print("PASS contact resultant: uniform nonvanishing for every m and r=5")
