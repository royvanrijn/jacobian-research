#!/usr/bin/env python3
"""Verify the first passive-flag collision over a pure sextic top.

For a ternary polynomial c put

    J(c) = grad(c)^T adj(Hess(c)) grad(c).

This checker treats the remaining scalar degree-six chart after HC4RSD28.
With c_6=x^6, the degree-sixteen face says that the passive Hessian of c_5
is singular.  The binary Hesse normal form is a structural characteristic-
zero input.  After that normal form, every identity and every terminal
collision used to show that c_5 has a fixed passive direction is checked
here exactly.
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
    / "hc4_pure_sextic_collision.json"
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


# The first nonzero face over x^6 is the passive Hessian determinant of an
# arbitrary homogeneous quintic correction.
quintic_monomials = [
    x**i * y**j * z ** (5 - i - j)
    for i in range(6)
    for j in range(6 - i)
]
quintic_coefficients = sp.symbols(f"u0:{len(quintic_monomials)}")
generic_quintic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(quintic_coefficients, quintic_monomials)
)
first_face = homogeneous_face(bordered_invariant(x**6 + generic_quintic), 16)
assert sp.expand(
    first_face
    - 36 * x**10 * sp.hessian(generic_quintic, (y, z)).det()
) == 0


# The passive binary Hesse theorem puts the quintic correction in the form
# H_5(x,y)+k*x^4*z.  If H_yy is nonzero, the next face fixes the complete
# z^2 part of the quartic correction.
a = sp.symbols("a0:6")
k = sp.symbols("k")
binary_quintic = sum(a[i] * x ** (5 - i) * y**i for i in range(6))
quartic_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]
quartic_coefficients = sp.symbols(f"b0:{len(quartic_monomials)}")
generic_quartic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(quartic_coefficients, quartic_monomials)
)
second_face = homogeneous_face(
    bordered_invariant(x**6 + binary_quintic + k * x**4 * z + generic_quartic),
    15,
)
assert sp.expand(
    second_face
    - 18
    * x**10
    * sp.diff(binary_quintic, y, 2)
    * (2 * sp.diff(generic_quartic, z, 2) - k**2 * x**2)
) == 0


# Hence c_4=R_4(x,y)+z*P_3(x,y)+k^2*x^2*z^2/4.  The following face first
# removes the y^5 and y^4 terms of H_5 and then quantizes the y^3 term of
# P_3.  It also fixes the z^3 coefficient of c_3 whenever H_yy is nonzero.
r = sp.symbols("r0:5")
p = sp.symbols("p0:4")
w = sp.symbols("w0:10")
binary_quartic = sum(r[i] * x ** (4 - i) * y**i for i in range(5))
transverse_cubic = sum(p[i] * x ** (3 - i) * y**i for i in range(4))
cubic_monomials = [
    x**i * y**j * z ** (3 - i - j)
    for i in range(4)
    for j in range(4 - i)
]
generic_cubic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(w, cubic_monomials)
)
collision_candidate = (
    x**6
    + binary_quintic
    + k * x**4 * z
    + binary_quartic
    + z * transverse_cubic
    + k**2 * x**2 * z**2 / 4
    + generic_cubic
)
collision_face = sp.Poly(
    homogeneous_face(bordered_invariant(collision_candidate), 14), x, y, z
)
assert sp.factor(
    collision_face.coeff_monomial(x**6 * y**8) + 25 * a[5] ** 2 * k**2
) == 0
assert sp.factor(
    collision_face.coeff_monomial(x**8 * y**6).subs(a[5], 0)
    + 24 * a[4] ** 2 * k**2
) == 0
top_cubic_collision = collision_face.coeff_monomial(x**10 * y**4).subs(
    {a[4]: 0, a[5]: 0}
)
assert sp.factor(
    top_cubic_collision
    + 3 * (18 * p[3] - 7 * a[3] * k) * (6 * p[3] - a[3] * k)
) == 0
assert sp.factor(
    collision_face.coeff_monomial(x**13 * z).subs({a[4]: 0, a[5]: 0})
    - 4 * a[2] * (-k**3 + 108 * w[0])
) == 0
assert sp.factor(
    collision_face.coeff_monomial(x**12 * y * z).subs({a[4]: 0, a[5]: 0})
    - 12 * a[3] * (-k**3 + 108 * w[0])
) == 0


# If k and H_yy are both nonzero, translations and the scaling torus leave
# three normalized charts.  Arbitrary terms capable of contributing to the
# next homogeneous face are retained.  Each chart has an immutable nonzero
# coefficient, so the moving collision is impossible.
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
tail_z_cubic = z * (
    tail_v[0] * x**2 + tail_v[1] * x * y + tail_v[2] * y**2
)
quadratic_monomials = (x**2, x * y, x * z, y**2, y * z, z**2)
generic_quadratic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(tail_q, quadratic_monomials)
)

q_p0, q_p1 = sp.symbols("q_p0 q_p1")
quadratic_passive_chart = (
    x**6
    + x**3 * y**2
    + x**4 * z
    + tail_binary_quartic
    + z * (q_p0 * x**3 + q_p1 * x**2 * y + x * y**2 / 3)
    + x**2 * z**2 / 4
    + tail_binary_cubic
    + tail_z_cubic
    + z**2
    * ((q_p0 / 3 + q_p1**2 / 4) * x + q_p1 * y / 6)
    + z**3 / 108
    + generic_quadratic
)
quadratic_passive_face = sp.Poly(
    homogeneous_face(bordered_invariant(quadratic_passive_chart), 13), x, y, z
)
assert quadratic_passive_face.coeff_monomial(x**9 * y**4) == 1

c_p0, c_p2 = sp.symbols("c_p0 c_p2")
cubic_passive_chart_one = (
    x**6
    + x**2 * y**3
    + x**4 * z
    + tail_binary_quartic
    + z * (c_p0 * x**3 + c_p2 * x * y**2 + y**3 / 6)
    + x**2 * z**2 / 4
    + tail_binary_cubic
    + tail_z_cubic
    + z**2 * (c_p0 * x + c_p2**2 * y) / 3
    + z**3 / 108
    + generic_quadratic
)
cubic_passive_face_one = sp.Poly(
    homogeneous_face(bordered_invariant(cubic_passive_chart_one), 13), x, y, z
)
assert cubic_passive_face_one.coeff_monomial(x**8 * y**4 * z) == 6

cubic_passive_chart_two = (
    x**6
    + x**2 * y**3
    + x**4 * z
    + tail_binary_quartic
    + z * (c_p0 * x**3 + sp.Rational(7, 18) * y**3)
    + x**2 * z**2 / 4
    + tail_binary_cubic
    + tail_z_cubic
    + c_p0 * x * z**2 / 3
    + z**3 / 108
    + generic_quadratic
)
cubic_passive_face_two = sp.Poly(
    homogeneous_face(bordered_invariant(cubic_passive_chart_two), 13), x, y, z
)
assert cubic_passive_face_two.coeff_monomial(x**8 * y**4 * z) == -sp.Rational(
    10, 9
)


payload = {
    "format": "hc4-pure-sextic-collision-v1",
    "status": {
        "id": "HC4RSD29",
        "kind": "hybrid theorem",
        "scope": "the pure-sixth leading chart in the scalar HC4RSD20 packet",
        "result": (
            "the degree-five correction has a fixed passive direction; a moving "
            "two-direction correction reduces to three charts with immutable "
            "nonzero next-face coefficients"
        ),
    },
    "faces": {
        "degree_16": "36*x^10*det Hess_(y,z)(c5)",
        "degree_15": "18*x^10*(c5)_yy*(2*(c4)_zz-k^2*x^2)",
        "degree_14": {
            "y5": "-25*a5^2*k^2",
            "y4_after_y5": "-24*a4^2*k^2",
            "cubic_split": "-3*(18*p3-7*a3*k)*(6*p3-a3*k)",
            "z3_alignment": "108*w9=k^3",
        },
        "degree_13_terminal_coefficients": ["1", "6", "-10/9"],
    },
    "residual": (
        "pure-sixth leading potential with its degree-five correction binary "
        "in one fixed passive direction; lower components may still break the flag"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: the pure-sixth first face is the passive quintic Hessian")
print("PASS: the next two faces reduce a moving correction to three charts")
print("PASS: the terminal chart coefficients are 1, 6, and -10/9")
print("THEOREM: the quintic correction has one fixed passive direction")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
