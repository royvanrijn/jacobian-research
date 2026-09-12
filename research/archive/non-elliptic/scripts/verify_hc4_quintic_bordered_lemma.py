#!/usr/bin/env python3
"""Verify the pure-fifth-power branch of the quintic bordered lemma.

For a ternary polynomial c put

    J(c) = grad(c)^T adj(Hess(c)) grad(c).

HC4RSD23 closes every degree-five leading form except a fifth power.  This
checker treats that remaining chart.  After normalizing c_5=x^5/20 it
verifies the homogeneous faces, the two full lower-tail unit ideals, and
the exact radical alignments used to prove that c omits a constant
direction.  The binary singular-Hessian normal forms used between the
checked faces are structural characteristic-zero inputs.
"""

from __future__ import annotations

import re
import shutil
import subprocess

import sympy as sp


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


def coefficient_equations(polynomial: sp.Expr) -> list[sp.Expr]:
    equations: list[sp.Expr] = []
    for _, coefficient in sp.Poly(sp.expand(polynomial), *variables).terms():
        numerator = sp.expand(sp.together(coefficient).as_numer_denom()[0])
        if numerator != 0 and numerator not in equations:
            equations.append(numerator)
    return equations


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the exact radical checks")


def run_singular(program: str, timeout: int = 240) -> str:
    completed = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    if completed.stderr.strip():
        raise RuntimeError(completed.stderr)
    return completed.stdout


def singular_branch(
    polynomial: sp.Expr,
    parameters: tuple[sp.Symbol, ...],
    *,
    marker: str,
    unit: bool = False,
    radical_generators: tuple[sp.Expr, ...] = (),
    maximum_power: int = 12,
) -> str:
    equations = coefficient_equations(bordered_invariant(polynomial))
    program = f"""
ring rr=0,({','.join(map(str, parameters))}),dp;
option(redSB);
ideal I={','.join(map(singular_expression, equations))};
ideal G=slimgb(I);
print("{marker}_BASE "+string(size(G))+" "+string(reduce(1,G)));
"""
    for index, generator in enumerate(radical_generators, start=1):
        expression = singular_expression(generator)
        program += f"""
poly p{index}={expression};
int n{index}=1;
while (n{index}<={maximum_power} && reduce(p{index},G)!=0)
{{
  p{index}=p{index}*({expression});
  n{index}++;
}}
print(
  "{marker}_POWER {index} "+string(n{index})+" "
  +string(reduce(p{index},G)==0)
);
"""
    output = run_singular(program)
    base = re.search(rf"(?m)^{marker}_BASE (\d+) ([01])$", output)
    assert base is not None
    if unit:
        assert base.group(1) == "1" and base.group(2) == "0"
    else:
        assert base.group(2) == "1"
    powers = re.findall(rf"(?m)^{marker}_POWER (\d+) (\d+) ([01])$", output)
    assert len(powers) == len(radical_generators)
    assert all(success == "1" for _, _, success in powers)
    return output


# Normalize the leading term to x^5/20.  For a generic quartic correction
# c4, the first nonzero homogeneous face is exactly its passive binary
# Hessian determinant.
quartic_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]
quartic_coefficients = sp.symbols(f"u0:{len(quartic_monomials)}")
generic_quartic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(quartic_coefficients, quartic_monomials)
)
pure_fifth = x**5 / 20
first_face = homogeneous_face(
    bordered_invariant(pure_fifth + generic_quartic), 12
)
assert sp.expand(
    first_face
    - x**8 * sp.hessian(generic_quartic, (y, z)).det() / 16
) == 0


# The passive singular-Hessian normal form is
#
#   c4=Q4(x,y)+k*x^3*z.
#
# When Q_yy is nonzero, the next face first fixes the z^2 part of c3.
a0, a1, a2, a3, a4, k = sp.symbols("a0 a1 a2 a3 a4 k")
binary_quartic = (
    a0 * x**4
    + a1 * x**3 * y
    + a2 * x**2 * y**2
    + a3 * x * y**3
    + a4 * y**4
)
cubic_coefficients = sp.symbols("b0:10")
cubic_monomials = [
    x**i * y**j * z ** (3 - i - j)
    for i in range(4)
    for j in range(4 - i)
]
generic_cubic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(cubic_coefficients, cubic_monomials)
)
pure_chart = pure_fifth + binary_quartic + k * x**3 * z + generic_cubic
face_11 = homogeneous_face(bordered_invariant(pure_chart), 11)
b0, b1, b2, b3, b4, b5, b6, b7, b8, b9 = cubic_coefficients
assert sp.expand(
    face_11
    - x**8
    * sp.diff(binary_quartic, y, 2)
    * (3 * b0 * z + b1 * y + b4 * x - 4 * k**2 * x)
    / 8
) == 0


# On Q_yy!=0 and k!=0, the next face forces a3=a4=0 and then has exactly
# two branches.  The scaling torus normalizes a2=k=1.  Both branches have
# unit complete coefficient ideals even after arbitrary binary cubic,
# quadratic, and linear tails.
d = sp.symbols("d0:4")
e = sp.symbols("e0:6")
r = sp.symbols("r0:3")
binary_cubic_tail = sum(d[i] * x ** (3 - i) * y**i for i in range(4))
quadratic_tail = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        e,
        (z**2, y * z, y**2, x * z, x * y, x**2),
    )
)
linear_tail = r[0] * x + r[1] * y + r[2] * z

branch_parameters_common = (a0, a1, *d, *e, *r)
branch_8 = (
    pure_fifth
    + a0 * x**4
    + a1 * x**3 * y
    + x**2 * y**2
    + x**3 * z
    + binary_cubic_tail
    + z * (8 * y**2 + 8 * a1 * x * y + 8 * a0 * x**2)
    + 4 * x * z**2
    + quadratic_tail
    + linear_tail
)
singular_branch(
    branch_8,
    branch_parameters_common,
    marker="CURVED8",
    unit=True,
)

free_b5 = sp.symbols("free_b5")
forced_b7 = 8 * a0 - (8 * a1 - free_b5) ** 2 / 16
branch_4 = (
    pure_fifth
    + a0 * x**4
    + a1 * x**3 * y
    + x**2 * y**2
    + x**3 * z
    + binary_cubic_tail
    + z * (4 * y**2 + free_b5 * x * y + forced_b7 * x**2)
    + 4 * x * z**2
    + quadratic_tail
    + linear_tail
)
singular_branch(
    branch_4,
    (a0, a1, free_b5, *d, *e, *r),
    marker="CURVED4",
    unit=True,
)


# If c4 is passive-affine and nonzero, normalize it to a*x^4+x^3*z.
# Its first cubic face is the prime rank-one system displayed in the note.
# The finite nonzero ratio chart has q=1; the complete radical forces
# independence of y.
pa, pb0, pb4, pb7, pb8, pb9 = sp.symbols(
    "pa pb0 pb4 pb7 pb8 pb9"
)
passive_finite = (
    pure_fifth
    + pa * x**4
    + x**3 * z
    + pb0 * (z + y) ** 3
    + x * (pb4 * (z + y) ** 2 - 8 * (z + y) * y + 4 * y**2)
    + pb7 * x**2 * z
    + pb8 * x**2 * y
    + pb9 * x**3
    + quadratic_tail
    + linear_tail
)
finite_alignment = (
    pb0,
    pb4 - 4,
    pb8,
    e[1],
    e[2],
    e[4],
    r[1],
)
singular_branch(
    passive_finite,
    (pa, pb0, pb4, pb7, pb8, pb9, *e, *r),
    marker="PASSIVE_FINITE",
    radical_generators=finite_alignment,
)


# At ratio q=0, the radical again removes every y coefficient.  This is
# the large but still exact 325-element Groebner calculation.
passive_zero = (
    pure_fifth
    + pa * x**4
    + x**3 * z
    + pb0 * z**3
    + pb4 * x * z**2
    + pb7 * x**2 * z
    + pb8 * x**2 * y
    + pb9 * x**3
    + quadratic_tail
    + linear_tail
)
zero_alignment = (pb8, e[1], e[2], e[4], r[1])
zero_output = singular_branch(
    passive_zero,
    (pa, pb0, pb4, pb7, pb8, pb9, *e, *r),
    marker="PASSIVE_ZERO",
    radical_generators=zero_alignment,
)
zero_base = re.search(r"(?m)^PASSIVE_ZERO_BASE (\d+) 1$", zero_output)
assert zero_base is not None and zero_base.group(1) == "325"


# The ratio-infinity cube is inconsistent.  On the residual conic boundary
# b4=4, b5=0, the full ideal directly contains b6 and falls back to q=0.
pb6 = sp.symbols("pb6")
passive_infinity = (
    pure_fifth
    + pa * x**4
    + x**3 * z
    + y**3
    + x * (4 * z**2 + pb6 * y**2)
    + pb7 * x**2 * z
    + pb8 * x**2 * y
    + pb9 * x**3
    + quadratic_tail
    + linear_tail
)
singular_branch(
    passive_infinity,
    (pa, pb6, pb7, pb8, pb9, *e, *r),
    marker="PASSIVE_INFINITY",
    unit=True,
)

passive_boundary = passive_infinity - y**3
singular_branch(
    passive_boundary,
    (pa, pb6, pb7, pb8, pb9, *e, *r),
    marker="PASSIVE_BOUNDARY",
    radical_generators=(pb6, pb8, e[1], e[2], e[4], r[1]),
)


# If c4 is binary and the first transverse cubic coefficient is nonzero,
# normalize that coefficient to x^2*z.  The exact radical aligns every
# lower passive term with z+d1*y, hence supplies a fixed missing direction.
qa = sp.symbols("qa0:5")
qd = sp.symbols("qd0:4")
qe = sp.symbols("qe0:6")
qr = sp.symbols("qr0:3")
nonzero_transverse = (
    pure_fifth
    + sum(qa[i] * x ** (4 - i) * y**i for i in range(5))
    + sum(qd[i] * x ** (3 - i) * y**i for i in range(4))
    + x**2 * z
    + sum(
        coefficient * monomial
        for coefficient, monomial in zip(
            qe,
            (z**2, y * z, y**2, x * z, x * y, x**2),
        )
    )
    + qr[0] * x
    + qr[1] * y
    + qr[2] * z
)
transverse_alignment = (
    *qa[1:],
    qd[2],
    qd[3],
    qe[1] - 2 * qd[1] * qe[0],
    qe[2] - qd[1] ** 2 * qe[0],
    qe[4] - qd[1] * qe[3],
    qr[1] - qd[1] * qr[2],
)
singular_branch(
    nonzero_transverse,
    (*qa, *qd, *qe, *qr),
    marker="TRANSVERSE",
    radical_generators=transverse_alignment,
)


# In the zero-transverse branch c=h(x,y)+z*(M*x+N*y+R)+F*z^2/2.
# The z^2 and z coefficients of J(c) share the exact unbordered Schur
# factor K.  The remaining binary analysis in the note gives a fixed
# direction; the final linear-in-z normal form has the square below.
hx, hy, hxx, hxy, hyy = sp.symbols("hx hy hxx hxy hyy")
M, N, F, R = sp.symbols("M N F R")
formal_gradient = sp.Matrix(
    [hx + M * z, hy + N * z, M * x + N * y + F * z + R]
)
formal_hessian = sp.Matrix(
    [[hxx, hxy, M], [hxy, hyy, N], [M, N, F]]
)
formal_bordered = sp.expand(
    (formal_gradient.T * formal_hessian.adjugate() * formal_gradient)[0]
)
formal_in_z = sp.Poly(formal_bordered, z)
schur_factor = F * (hxx * hyy - hxy**2) - (
    M**2 * hyy - 2 * M * N * hxy + N**2 * hxx
)
assert sp.expand(formal_in_z.coeff_monomial(z**2) - F * schur_factor) == 0
assert sp.expand(
    formal_in_z.coeff_monomial(z) - 2 * (M * x + N * y + R) * schur_factor
) == 0

fp, fpp, g, gp, gpp = sp.symbols("fp fpp g gp gpp")
terminal_gradient = sp.Matrix([fp + y * gp + z, g, x])
terminal_hessian = sp.Matrix(
    [[fpp + y * gpp, gp, 1], [gp, 0, 0], [1, 0, 0]]
)
terminal_bordered = sp.factor(
    (terminal_gradient.T * terminal_hessian.adjugate() * terminal_gradient)[0]
)
assert terminal_bordered == -(x * gp - g) ** 2


print("PASS: the pure-fifth first face is the passive quartic Hessian")
print("PASS: both curved quartic-correction branches have unit full ideals")
print("PASS: every passive-affine quartic branch is fixed or inconsistent")
print("PASS: a nonzero transverse cubic aligns every lower passive term")
print("PASS: the zero-transverse branch reduces to a binary Schur square")
print("THEOREM: every pure-fifth-power degree-five border coefficient is fixed")
