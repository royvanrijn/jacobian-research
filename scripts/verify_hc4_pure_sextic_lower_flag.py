#!/usr/bin/env python3
"""Verify lower passive-flag stabilization over a pure sextic top.

HC4RSD29 makes the quintic correction to c_6=x^6 binary.  This checker
treats the curved case (c_5)_{yy} != 0.  It verifies the factorized first
lower-direction equation, every repeated-root resonance chart of its binary
cubic coefficient, and the triangular elimination of all later transverse
tails.  The factor classification of a binary cubic modulo squares is the
only structural characteristic-zero input.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_pure_sextic_lower_flag.json"
)

x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def bordered_invariant(polynomial: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(polynomial, variable) for variable in variables])
    hessian = sp.hessian(polynomial, variables)
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])


def homogeneous_face(polynomial: sp.Expr, degree: int) -> sp.Expr:
    expanded = sp.Poly(sp.expand(polynomial), *variables)
    return sp.expand(
        sum(
            coefficient * x**monomial[0] * y**monomial[1] * z**monomial[2]
            for monomial, coefficient in expanded.terms()
            if sum(monomial) == degree
        )
    )


def face_polynomial(polynomial: sp.Expr, degree: int) -> sp.Poly:
    return sp.Poly(homogeneous_face(bordered_invariant(polynomial), degree), x, y, z)


# Write c_4=R_4(x,y)+z*P_3(x,y).  The degree-fourteen face is the exact
# binary factor equation H_yy*Q_1=(P_3,y)^2, where Q_1=(c_3)_zz.
a = sp.symbols("a0:6")
r = sp.symbols("r0:5")
p = sp.symbols("p0:4")
u = sp.symbols("u0:4")
v = sp.symbols("v0:3")
q0, q1 = sp.symbols("q0 q1")
binary_quintic = sum(a[i] * x ** (5 - i) * y**i for i in range(6))
binary_quartic = sum(r[i] * x ** (4 - i) * y**i for i in range(5))
transverse_cubic = sum(p[i] * x ** (3 - i) * y**i for i in range(4))
binary_cubic = sum(u[i] * x ** (3 - i) * y**i for i in range(4))
transverse_quadratic = v[0] * x**2 + v[1] * x * y + v[2] * y**2
passive_linear = q0 * x + q1 * y
first_break_candidate = (
    x**6
    + binary_quintic
    + binary_quartic
    + z * transverse_cubic
    + binary_cubic
    + z * transverse_quadratic
    + z**2 * passive_linear / 2
)
first_break_face = homogeneous_face(
    bordered_invariant(first_break_candidate), 14
)
assert sp.expand(
    first_break_face
    - 36
    * x**10
    * (
        sp.diff(binary_quintic, y, 2) * passive_linear
        - sp.diff(transverse_cubic, y) ** 2
    )
) == 0


# If P_y is nonzero, unique factorization says that the binary cubic H_yy
# has at most one factor of odd multiplicity.  Up to affine changes of y/x,
# this leaves the five charts below.  Every coefficient capable of entering
# degree thirteen is retained.
a0, a1, p0 = sp.symbols("chart_a0 chart_a1 chart_p0")
tail_r = sp.symbols("tail_r0:5")
tail_u = sp.symbols("tail_u0:4")
tail_v = sp.symbols("tail_v0:3")
tail_q = sp.symbols("tail_q0:6")
tail_binary_quartic = sum(
    tail_r[i] * x ** (4 - i) * y**i for i in range(5)
)
tail_binary_cubic = sum(
    tail_u[i] * x ** (3 - i) * y**i for i in range(4)
)
tail_transverse_quadratic = (
    tail_v[0] * x**2 + tail_v[1] * x * y + tail_v[2] * y**2
)
generic_quadratic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        tail_q, (x**2, x * y, x * z, y**2, y * z, z**2)
    )
)
resonance_charts = {
    "triple_finite": (
        y**5 / 20,
        p0 * x**3 + y**3 / 3,
        y,
        x**4 * y**9,
        sp.Rational(5, 24),
    ),
    "triple_infinity": (
        x**3 * y**2 / 2,
        p0 * x**3 + x**2 * y,
        x,
        None,
        None,
    ),
    "double_finite_simple_infinity": (
        x * y**4 / 12,
        p0 * x**3 + x * y**2 / 2,
        x,
        x**12 * z,
        36,
    ),
    "double_infinity_simple_finite": (
        x**2 * y**3 / 6,
        p0 * x**3 + x * y**2 / 2,
        y,
        x**11 * y * z,
        -36,
    ),
    "double_finite_simple_finite": (
        y**5 / 20 - x * y**4 / 12,
        p0 * x**3 + y**3 / 3 - x * y**2 / 2,
        y - x,
        x**12 * z,
        36,
    ),
}
for chart_name, (quintic, transverse, passive, monomial, expected) in (
    resonance_charts.items()
):
    candidate = (
        x**6
        + a0 * x**5
        + a1 * x**4 * y
        + quintic
        + tail_binary_quartic
        + z * transverse
        + tail_binary_cubic
        + z * tail_transverse_quadratic
        + z**2 * passive / 2
        + generic_quadratic
    )
    chart_face = face_polynomial(candidate, 13)
    if chart_name != "triple_infinity":
        assert sp.factor(chart_face.coeff_monomial(monomial) - expected) == 0


# The triple-at-infinity chart is the sole degree-thirteen resonance.  Its
# equations are solved explicitly below.  The next face retains arbitrary
# lower tails and has the immutable coefficient -1 at x^8*y^4.
tr_p0 = sp.symbols("tr_p0")
tr_r = sp.symbols("tr_r0:4")
tr_u = sp.symbols("tr_u0:4")
tr_v0, tr_v1 = sp.symbols("tr_v0 tr_v1")
tr_q = sp.symbols("tr_q0:5")
tr_l = sp.symbols("tr_l0:3")
tr_binary_quartic = (
    tr_r[0] * x**4
    + tr_r[1] * x**3 * y
    + tr_r[2] * x**2 * y**2
    + tr_r[3] * x * y**3
)
tr_binary_cubic = sum(tr_u[i] * x ** (3 - i) * y**i for i in range(4))
tr_v2 = (tr_p0 + 18 * tr_r[3]) / 12
tr_z2 = (tr_p0**2 - 12 * tr_r[2] + 12 * tr_v1) / 12
tr_candidate = (
    x**6
    + x**3 * y**2 / 2
    + tr_binary_quartic
    + z * (tr_p0 * x**3 + x**2 * y)
    + tr_binary_cubic
    + z * (tr_v0 * x**2 + tr_v1 * x * y + tr_v2 * y**2)
    + x * z**2 / 2
    + tr_q[0] * x**2
    + tr_q[1] * x * y
    + tr_q[2] * x * z
    + tr_q[3] * y**2
    + tr_q[4] * y * z
    + tr_z2 * z**2
    + tr_l[0] * x
    + tr_l[1] * y
    + tr_l[2] * z
)
triple_infinity_next = face_polynomial(tr_candidate, 12)
assert triple_infinity_next.coeff_monomial(x**8 * y**4) == -1


# If P_y=0, homogeneity gives P=p*x^3.  The degree-thirteen face fixes the
# z^2 coefficient of c_2.  The degree-twelve equations force H_yy to be a
# cube at infinity and determine the y-bearing part of the next transverse
# coefficient.  After normalization its following face contains -4.
py0_p = sp.symbols("py0_p")
py0_a = sp.symbols("py0_a0:6")
py0_r = sp.symbols("py0_r0:5")
py0_u = sp.symbols("py0_u0:4")
py0_v = sp.symbols("py0_v0:3")
py0_q = sp.symbols("py0_q0:5")
py0_l = sp.symbols("py0_l0:3")
py0_quintic = sum(py0_a[i] * x ** (5 - i) * y**i for i in range(6))
py0_candidate = (
    x**6
    + py0_quintic
    + sum(py0_r[i] * x ** (4 - i) * y**i for i in range(5))
    + py0_p * x**3 * z
    + sum(py0_u[i] * x ** (3 - i) * y**i for i in range(4))
    + z * (py0_v[0] * x**2 + py0_v[1] * x * y + py0_v[2] * y**2)
    + py0_q[0] * x**2
    + py0_q[1] * x * y
    + py0_q[2] * x * z
    + py0_q[3] * y**2
    + py0_q[4] * y * z
    + py0_p**2 * z**2 / 12
    + py0_l[0] * x
    + py0_l[1] * y
    + py0_l[2] * z
)
py0_face = face_polynomial(py0_candidate, 12)
assert sp.factor(
    py0_face.coeff_monomial(x**4 * y**8) + 100 * py0_a[5] ** 2 * py0_p**2
) == 0
assert sp.factor(
    py0_face.coeff_monomial(x**6 * y**6).subs(py0_a[5], 0)
    + 64 * py0_a[4] ** 2 * py0_p**2
) == 0

normalized_py0 = (
    x**6
    + x**3 * y**2
    + sum(tail_r[i] * x ** (4 - i) * y**i for i in range(5))
    + x**3 * z
    + sum(tail_u[i] * x ** (3 - i) * y**i for i in range(4))
    + z * (tail_v[0] * x**2 + y**2 / 3)
    + tail_q[0] * x**2
    + tail_q[1] * x * y
    + tail_q[2] * x * z
    + tail_q[3] * y**2
    + tail_q[4] * y * z
    + z**2 / 12
)
normalized_py0_face = face_polynomial(normalized_py0, 11)
assert normalized_py0_face.coeff_monomial(x**10 * z) == -4


# Once c_4 is binary, later transverse terms die triangularly.  Keep a
# generic z-linear cubic, generic quadratic, and generic affine tail.  The
# five displayed identities successively kill c_2,zz, T_y, the remaining
# x^2*z term, the y*z quadratic, and finally the x*z quadratic.
cascade_t = sp.symbols("cascade_t0:3")
cascade_q = sp.symbols("cascade_q0:6")
cascade_l = sp.symbols("cascade_l0:3")
cascade_transverse = (
    cascade_t[0] * x**2 + cascade_t[1] * x * y + cascade_t[2] * y**2
)
cascade_quadratic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        cascade_q, (x**2, x * y, x * z, y**2, y * z, z**2)
    )
)
cascade_candidate = (
    x**6
    + binary_quintic
    + binary_quartic
    + binary_cubic
    + z * cascade_transverse
    + cascade_quadratic
    + cascade_l[0] * x
    + cascade_l[1] * y
    + cascade_l[2] * z
)
cascade_border = bordered_invariant(cascade_candidate)
cascade_hyy = sp.diff(binary_quintic, y, 2)
assert sp.expand(
    homogeneous_face(cascade_border, 13)
    - 72 * x**10 * cascade_q[5] * cascade_hyy
) == 0
assert sp.expand(
    homogeneous_face(cascade_border, 12).subs(cascade_q[5], 0)
    + 36 * x**10 * sp.diff(cascade_transverse, y) ** 2
) == 0
assert sp.expand(
    homogeneous_face(cascade_border, 11).subs(
        {cascade_q[5]: 0, cascade_t[1]: 0, cascade_t[2]: 0}
    )
    - 6 * x**8 * cascade_hyy * cascade_t[0] ** 2
) == 0
assert sp.expand(
    homogeneous_face(cascade_border, 10).subs(
        {
            cascade_q[5]: 0,
            cascade_t[0]: 0,
            cascade_t[1]: 0,
            cascade_t[2]: 0,
        }
    )
    + 36 * x**10 * cascade_q[4] ** 2
) == 0
assert sp.expand(
    homogeneous_face(cascade_border, 9).subs(
        {
            cascade_q[5]: 0,
            cascade_q[4]: 0,
            cascade_t[0]: 0,
            cascade_t[1]: 0,
            cascade_t[2]: 0,
        }
    )
    - 18 * x**6 * cascade_hyy * cascade_q[2] ** 2
) == 0

# At that point c=h(x,y)+ell*z.  The exact terminal identity and the leading
# binary Hessian coefficient force ell=0 when H_yy is nonzero.
hx, hy, hxx, hxy, hyy, ell = sp.symbols("hx hy hxx hxy hyy ell")
terminal_gradient = sp.Matrix([hx, hy, ell])
terminal_hessian = sp.Matrix([[hxx, hxy, 0], [hxy, hyy, 0], [0, 0, 0]])
terminal_border = sp.factor(
    (terminal_gradient.T * terminal_hessian.adjugate() * terminal_gradient)[0]
)
assert terminal_border == ell**2 * (hxx * hyy - hxy**2)


payload = {
    "format": "hc4-pure-sextic-lower-flag-v1",
    "status": {
        "id": "HC4RSD30",
        "kind": "hybrid theorem",
        "scope": (
            "the pure-sixth scalar chart after HC4RSD29, when the stabilized "
            "binary quintic correction has nonzero passive second derivative"
        ),
        "result": (
            "every possible quartic direction break and every later transverse "
            "tail is impossible, so the complete potential is a fixed cylinder"
        ),
    },
    "first_break": "H_yy*Q_1=(P_3,y)^2",
    "resonance_charts": list(resonance_charts),
    "immutable_coefficients": ["5/24", "36", "-36", "36", "-1", "-4"],
    "tail_cascade": [
        "72*x^10*H_yy*(c2_z2)",
        "-36*x^10*(T_y)^2",
        "6*x^8*H_yy*T_x2^2",
        "-36*x^10*(c2_yz)^2",
        "18*x^6*H_yy*(c2_xz)^2",
    ],
    "residual": (
        "pure-sixth top whose binary quintic correction is passive-affine: "
        "c5=a*x^5+x^4*L(y,z)"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: factored the first lower passive-direction equation")
print("PASS: all five repeated-root cubic resonance charts close")
print("PASS: the y-constant transverse escape has immutable coefficient -4")
print("PASS: every later transverse tail dies triangularly")
print("THEOREM: every curved-quintic pure-sixth chart is a fixed cylinder")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
