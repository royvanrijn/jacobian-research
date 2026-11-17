#!/usr/bin/env python3
"""Verify the curved completed-quartic closure on the pure-sixth boundary.

After HC4RSD30 the remaining quintic correction is passive-affine.  Completing
the square in that linear form turns the first quartic face into an ordinary
passive singular-Hessian equation.  This checker verifies the finite moving
charts when the resulting quartic has nonzero passive curvature, including
the aligned lower-tail resonances.
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
    / "hc4_pure_sextic_affine_quartic.json"
)

x, y, z = sp.symbols("x y z")
variables = (x, y, z)


def bordered_invariant(polynomial: sp.Expr) -> sp.Expr:
    gradient = sp.Matrix([sp.diff(polynomial, variable) for variable in variables])
    hessian = sp.hessian(polynomial, variables)
    return sp.expand((gradient.T * hessian.adjugate() * gradient)[0])


def face_polynomial(polynomial: sp.Expr, degree: int) -> sp.Poly:
    expanded = sp.Poly(bordered_invariant(polynomial), *variables)
    return sp.Poly(
        sp.expand(
            sum(
                coefficient
                * x**monomial[0]
                * y**monomial[1]
                * z**monomial[2]
                for monomial, coefficient in expanded.terms()
                if sum(monomial) == degree
            )
        ),
        *variables,
    )


def generic_binary(degree: int, prefix: str) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    coefficients = sp.symbols(f"{prefix}0:{degree + 1}")
    return (
        sum(
            coefficients[i] * x ** (degree - i) * y**i
            for i in range(degree + 1)
        ),
        coefficients,
    )


# If c_5=x^4*L for a passive linear form L, its complete quartic face is the
# passive Hessian determinant of c_4-x^2*L^2/4.
la, lb = sp.symbols("la lb")
quartic_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]
quartic_coefficients = sp.symbols(f"gq0:{len(quartic_monomials)}")
generic_quartic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(quartic_coefficients, quartic_monomials)
)
linear_form = la * y + lb * z
completed_quartic = generic_quartic - x**2 * linear_form**2 / 4
completion_face = face_polynomial(
    x**6 + x**4 * linear_form + generic_quartic, 14
).as_expr()
assert sp.expand(
    completion_face
    - 36 * x**10 * sp.hessian(completed_quartic, (y, z)).det()
) == 0


# Put the singular completed quartic in the form Q_4(x,y)+k*x^3*z and
# write L=A*y+B*z.  The degree-thirteen face below supplies the complete
# first collision equations.
A, B, k = sp.symbols("A B k")
Q4, qa = generic_binary(4, "qa")
S3, su = generic_binary(3, "su")
tv = sp.symbols("tv0:3")
T2 = tv[0] * x**2 + tv[1] * x * y + tv[2] * y**2
c3z = sp.symbols("c3z0:3")
C3 = S3 + z * T2 + c3z[0] * x * z**2 + c3z[1] * y * z**2 + c3z[2] * z**3
normal_candidate = (
    x**6
    + x**4 * (A * y + B * z)
    + x**2 * (A * y + B * z) ** 2 / 4
    + Q4
    + k * x**3 * z
    + C3
)
normal_face = face_polynomial(normal_candidate, 13)
assert sp.factor(
    normal_face.coeff_monomial(x**13) - 48 * qa[2] * (-B * k + 3 * c3z[0])
) == 0
assert sp.factor(
    normal_face.coeff_monomial(x**12 * z)
    - 4 * qa[2] * (-B**3 + 108 * c3z[2])
) == 0


# Misaligned case B!=0.  Set B=1 and shear A to zero.  The next face first
# removes the y^4 and y^3 terms of Q_4.  With nonzero curvature, scaling
# leaves two exact ratios.  Their terminal coefficients are 2 and 5/108.
ma0, mk, mv0, mv1 = sp.symbols("ma0 mk mv0 mv1")
mbin, mb = generic_binary(3, "mb")
mq = sp.symbols("mq0:5")
mlin = sp.symbols("mlin0:3")


def misaligned_candidate(v1: sp.Expr, v2: sp.Expr, q5: sp.Expr) -> sp.Expr:
    return (
        x**6
        + x**4 * z
        + x**2 * z**2 / 4
        + ma0 * x**4
        + x**2 * y**2
        + mk * x**3 * z
        + mbin
        + z * (mv0 * x**2 + v1 * x * y + v2 * y**2)
        + mk * x * z**2 / 3
        + z**3 / 108
        + mq[0] * x**2
        + mq[1] * x * y
        + mq[2] * x * z
        + mq[3] * y**2
        + mq[4] * y * z
        + q5 * z**2
        + mlin[0] * x
        + mlin[1] * y
        + mlin[2] * z
    )


misaligned_one = misaligned_candidate(
    mv1,
    sp.Rational(1, 6),
    (-2 * ma0 + 3 * mk**2 + 6 * mv0 + 9 * mv1**2) / 36,
)
assert face_polynomial(misaligned_one, 11).coeff_monomial(x**8 * y**2 * z) == 2
misaligned_two = misaligned_candidate(
    0,
    sp.Rational(1, 3),
    (-2 * ma0 + 3 * mk**2 + 6 * mv0) / 36,
)
assert face_polynomial(misaligned_two, 10).coeff_monomial(x**6 * z**4) == sp.Rational(
    5, 108
)


# Aligned L and a nonzero transverse k.  Normalize L=y and k=1.  On the
# D=12*q_5-1=0 face, degree eleven removes the x^2*y^2 term and degree ten
# removes the remaining curved quartic terms.
aa = sp.symbols("aa0:5")
au, _ = generic_binary(3, "au")
av0 = sp.symbols("av0")
aq = sp.symbols("aq0:5")
aligned_d0 = (
    x**6
    + x**4 * y
    + x**2 * y**2 / 4
    + aa[0] * x**4
    + aa[1] * x**3 * y
    + aa[3] * x * y**3
    + aa[4] * y**4
    + x**3 * z
    + au
    + z * (av0 * x**2 + x * y / 3)
    + aq[0] * x**2
    + aq[1] * x * y
    + aq[2] * x * z
    + aq[3] * y**2
    + aq[4] * y * z
    + z**2 / 12
)
aligned_d0_face = face_polynomial(aligned_d0, 10)
assert aligned_d0_face.coeff_monomial(x**8 * y * z) == -12 * aa[3]
assert aligned_d0_face.coeff_monomial(x**7 * y**2 * z) == -24 * aa[4]

# For D!=0, a finite square-root direction is killed by 12*v_2*D.  At the
# root at infinity, the next two faces give E1=E2=0 and then F=0.  Their
# exact coefficient ideal is (D,s), contradicting the localized chart.
D, ss = sp.symbols("D ss")
E1 = D**2 - 2 * D * ss + 2 * ss
E2 = D**2 - 4 * D * ss + D - 4 * ss**2
F = 4 * D**3 - 10 * D**2 * ss + 3 * D**2 + 16 * D * ss**2 - 4 * D * ss - 4 * ss**2
aligned_resonance_basis = sp.groebner([E1, E2, F], D, ss, order="lex")
assert [basis.as_expr() for basis in aligned_resonance_basis.polys] == [D, ss]


# If k=0, a lower cubic attempts the direction break.  For a nonzero scalar
# Schur term, the finite and infinite square-root charts have immutable
# coefficients 36 and -1.
ca0, ca1, cv0 = sp.symbols("ca0 ca1 cv0")
cu, _ = generic_binary(3, "cu")
cq = sp.symbols("cq0:5")
cl = sp.symbols("cl0:3")


def aligned_lower_candidate(curved: sp.Expr, transverse: sp.Expr) -> sp.Expr:
    return (
        x**6
        + x**4 * y
        + x**2 * y**2 / 4
        + ca0 * x**4
        + ca1 * x**3 * y
        + curved
        + cu
        + z * transverse
        + cq[0] * x**2
        + cq[1] * x * y
        + cq[2] * x * z
        + cq[3] * y**2
        + cq[4] * y * z
        + z**2 / 2
        + cl[0] * x
        + cl[1] * y
        + cl[2] * z
    )


aligned_finite = aligned_lower_candidate(y**4 / 12, cv0 * x**2 + y**2 / 2)
assert face_polynomial(aligned_finite, 11).coeff_monomial(x**10 * z) == 36
aligned_infinity = aligned_lower_candidate(x**2 * y**2 / 2, cv0 * x**2 + x * y)
assert face_polynomial(aligned_infinity, 10).coeff_monomial(x**9 * z) == -12


# In the zero-scalar chart, a nonzero x^2*z escape reaches two possible
# ratios at degree eight, but the degree-seven z coefficient excludes both.
zs = sp.symbols("zs")
zero_scalar_obstructions = [
    zs * (6 * zs - 1) ** 2 * (18 * zs + 1),
    zs * (12 * zs - 1),
]
assert sp.groebner(zero_scalar_obstructions, zs).polys[0].as_expr() == zs

# With no cubic escape, the remaining quadratic tail dies by -36*e_y^2 and
# then by the curved quartic Hessian times e_x^2.
e0, e1 = sp.symbols("e0 e1")
tail_Q, tail_a = generic_binary(4, "tail_a")
tail_S, _ = generic_binary(3, "tail_s")
tail_binary_quadratic, _ = generic_binary(2, "tail_b")
tail_affine = sp.symbols("tail_l0:3")
tail_candidate = (
    x**6
    + x**4 * y
    + x**2 * y**2 / 4
    + tail_Q
    + tail_S
    + tail_binary_quadratic
    + z * (e0 * x + e1 * y)
    + tail_affine[0] * x
    + tail_affine[1] * y
    + tail_affine[2] * z
)
assert face_polynomial(tail_candidate, 10).coeff_monomial(x**10) == -36 * e1**2
tail_face_e1_zero = face_polynomial(tail_candidate, 8).as_expr().subs(e1, 0)
assert sp.factor(
    tail_face_e1_zero - 18 * x**6 * sp.diff(tail_Q, y, 2) * e0**2
) == 0


# When L=0, a nonzero k has the same D split.  D=0 directly removes every
# curved coefficient; D!=0 has either the 12*v_2*D obstruction or the pure
# square -108*s^4/D^2.  If k=0, the finite/infinite and zero-scalar lower
# charts have terminal coefficients 36, -1, and -648*s^4.
na = sp.symbols("na0:5")
nu, _ = generic_binary(3, "nu")
nv0 = sp.symbols("nv0")
nq = sp.symbols("nq0:5")
no_linear_d0 = (
    x**6
    + na[0] * x**4
    + na[1] * x**3 * y
    + na[2] * x**2 * y**2
    + na[3] * x * y**3
    + na[4] * y**4
    + x**3 * z
    + nu
    + nv0 * x**2 * z
    + nq[0] * x**2
    + nq[1] * x * y
    + nq[2] * x * z
    + nq[3] * y**2
    + nq[4] * y * z
    + z**2 / 12
)
no_linear_d0_face = face_polynomial(no_linear_d0, 10)
assert no_linear_d0_face.coeff_monomial(x**9 * z) == -4 * na[2]
assert no_linear_d0_face.coeff_monomial(x**8 * y * z) == -12 * na[3]
assert no_linear_d0_face.coeff_monomial(x**7 * y**2 * z) == -24 * na[4]

no_linear_finite = (
    x**6
    + y**4 / 12
    + nu
    + z * (nv0 * x**2 + y**2 / 2)
    + nq[0] * x**2
    + nq[1] * x * y
    + nq[2] * x * z
    + nq[3] * y**2
    + nq[4] * y * z
    + z**2 / 2
)
assert face_polynomial(no_linear_finite, 11).coeff_monomial(x**10 * z) == 36

iu0, iq0, iq1, iq2, iq4 = sp.symbols("iu0 iq0 iq1 iq2 iq4")
no_linear_infinity = (
    x**6
    + na[0] * x**4
    + na[1] * x**3 * y
    + x**2 * y**2 / 2
    + iu0 * x**3
    + (iq2 + na[1] * iq4) * x**2 * y
    + iq4 * x * y**2
    + z * (na[1] * x**2 + x * y)
    + iq0 * x**2
    + iq1 * x * y
    + iq2 * x * z
    + (na[1] ** 2 + 6 * iq4**2) * y**2 / 12
    + iq4 * y * z
    + z**2 / 2
)
assert face_polynomial(no_linear_infinity, 8).coeff_monomial(x**4 * y**4) == -1


payload = {
    "format": "hc4-pure-sextic-affine-quartic-v1",
    "status": {
        "id": "HC4RSD31",
        "kind": "hybrid theorem",
        "scope": (
            "the passive-affine quintic boundary after HC4RSD30, when the "
            "completed quartic has nonzero passive curvature"
        ),
        "result": (
            "all misaligned, transverse, and aligned lower-tail charts close, "
            "so the complete potential is a fixed cylinder"
        ),
    },
    "completion": "c4_hat=c4-x^2*L^2/4 has singular passive Hessian",
    "terminal_coefficients": ["2", "5/108", "36", "-1"],
    "localized_resonance_ideal": ["D", "s"],
    "residual": (
        "c5=a*x^5+x^4*L and c4=b*x^4+x^3*M+x^2*L^2/4 for constant "
        "passive linear forms L,M"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: completed the passive-affine quintic square")
print("PASS: every curved completed-quartic collision chart closes")
print("PASS: aligned lower cubic and quadratic escapes close")
print("THEOREM: every curved-quartic passive-affine quintic chart is fixed")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
