#!/usr/bin/env python3
"""Verify the rank-one sextic faces in the complete coordinate Meng chart.

Let

    psi = q_2 + h_3 + h_4 + h_6

in four variables and suppose rank Hess(h_6)=1.  The small-rank Hessian
classification and homogeneity put h_6=c*x^6.  Write W=(u,v,w) for its
constant three-dimensional kernel.  This checker verifies the successive
determinant faces used in the proof:

* lambda^10 is c*det Hess_W(h_4);
* on the rank-two quartic stratum, lambda^9 forces the cubic to vanish
  twice in the constant quartic-kernel direction;
* on the rank-one quartic stratum, lambda^8 is the binary discriminant of
  Hess(h_3) on the constant quartic-kernel plane; the apparent competing
  det Hess(h_4) term vanishes;
* on the rank-zero quartic stratum, lambda^7 is c*det Hess_W(h_3).

The binary discriminant is the Sym^2 nullcone L^2 familiar from the SIC(2)
binary-root calculation.  Its UFD factorization gives a constant direction.
The conceptual inputs not checked here are the small-rank polynomial-Hessian
normal form, its homogeneous specialization h_6=c*x^6, the ternary
singular-Hessian theorem, and the already proved common-direction reduction
HC4T31.
"""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import runpy

import sympy as sp


# Replay the common-direction descent and its terminal HC(3)/JC(2)/HC(2)
# identities.
shared_path = Path(__file__).with_name(
    "verify_hc4_meng_triple_rank_three_reduction.py"
)
with contextlib.redirect_stdout(io.StringIO()):
    runpy.run_path(str(shared_path))


lam = sp.symbols("lam")


def generic_symmetric(prefix: str) -> sp.Matrix:
    entries = sp.symbols(f"{prefix}0:10")
    return sp.Matrix(
        [
            [entries[0], entries[1], entries[2], entries[3]],
            [entries[1], entries[4], entries[5], entries[6]],
            [entries[2], entries[5], entries[7], entries[8]],
            [entries[3], entries[6], entries[8], entries[9]],
        ]
    )


# Variables are ordered as W=(u,v,w), followed by the sextic coordinate x.
c = sp.symbols("c", nonzero=True)
C6 = sp.diag(0, 0, 0, c)
H0 = generic_symmetric("h")
A3 = generic_symmetric("a")
B4 = generic_symmetric("b")

spatial = sp.Poly(
    (H0 + lam * A3 + lam**2 * B4 + lam**4 * C6).det(
        method="berkowitz"
    ),
    lam,
)
assert sp.expand(
    spatial.coeff_monomial(lam**10) - c * B4[:3, :3].det()
) == 0


# Rank two of B_W.  After a congruence over the function field, take
# B_W=diag(beta_1,beta_2,0).  The lambda^9 face is then the nonzero
# quotient determinant times A_ww.
beta_1, beta_2 = sp.symbols("beta_1 beta_2", nonzero=True)
rank_two_substitutions = {
    B4[0, 0]: beta_1,
    B4[0, 1]: 0,
    B4[0, 2]: 0,
    B4[1, 1]: beta_2,
    B4[1, 2]: 0,
    B4[2, 2]: 0,
}
rank_two_spatial = sp.Poly(
    (
        H0
        + lam * A3
        + lam**2 * B4.subs(rank_two_substitutions)
        + lam**4 * C6
    ).det(method="berkowitz"),
    lam,
)
assert sp.expand(
    rank_two_spatial.coeff_monomial(lam**9)
    - c * beta_1 * beta_2 * A3[2, 2]
) == 0


# Rank one of B_W.  The same normalization makes its only nonzero entry
# beta.  At lambda^8 the C-selected term is beta times the determinant of
# A on the two-dimensional kernel.  The only no-C term of the same weight
# is det(B4), which vanishes identically: a symmetric 1+3 block whose
# 3x3 principal block has rank one has total rank at most three.
beta = sp.symbols("beta", nonzero=True)
rank_one_substitutions = {
    B4[0, 0]: beta,
    B4[0, 1]: 0,
    B4[0, 2]: 0,
    B4[1, 1]: 0,
    B4[1, 2]: 0,
    B4[2, 2]: 0,
}
rank_one_B = B4.subs(rank_one_substitutions)
assert sp.expand(rank_one_B.det(method="berkowitz")) == 0
rank_one_spatial = sp.Poly(
    (
        H0 + lam * A3 + lam**2 * rank_one_B + lam**4 * C6
    ).det(method="berkowitz"),
    lam,
)
binary_discriminant = A3[1, 1] * A3[2, 2] - A3[1, 2] ** 2
assert sp.expand(
    rank_one_spatial.coeff_monomial(lam**8)
    - c * beta * binary_discriminant
) == 0


# The SIC(2) Sym^2 nullcone is precisely r0*r2-r1^2=0.  Its standard
# repeated-root parametrization makes the corresponding symmetric matrix
# rank one with a constant kernel direction.
ell, root_0, root_1 = sp.symbols("ell root_0 root_1")
nullcone_matrix = ell * sp.Matrix(
    [
        [root_0**2, root_0 * root_1],
        [root_0 * root_1, root_1**2],
    ]
)
repeated_root_direction = sp.Matrix([-root_1, root_0])
assert sp.expand(nullcone_matrix.det()) == 0
assert nullcone_matrix * repeated_root_direction == sp.zeros(2, 1)


# Rank zero of B_W.  Such a quartic has the homogeneous form
# x^4*a+x^3*ell(W), hence its full Hessian has rank at most two.  The
# possible no-C contribution to lambda^7 uses three B columns and one A
# column and therefore vanishes.  Direct substitution verifies that the
# surviving face is c*det(A_W).
rank_zero_substitutions = {
    B4[0, 0]: 0,
    B4[0, 1]: 0,
    B4[0, 2]: 0,
    B4[1, 1]: 0,
    B4[1, 2]: 0,
    B4[2, 2]: 0,
}
rank_zero_B = B4.subs(rank_zero_substitutions)
rank_zero_spatial = sp.Poly(
    (
        H0 + lam * A3 + lam**2 * rank_zero_B + lam**4 * C6
    ).det(method="berkowitz"),
    lam,
)
assert sp.expand(
    rank_zero_spatial.coeff_monomial(lam**7)
    - c * A3[:3, :3].det()
) == 0


# Check the homogeneous normal forms used to dispose of the same-weight
# quartic terms.  The rank-one W-Hessian form depends nonlinearly on one
# constant W-linear form and linearly on at most one additional W-form.
u, v, w, x, scale = sp.symbols("u v w x scale")
a0, a1, a2, a3, a4, lv, lw = sp.symbols(
    "a0 a1 a2 a3 a4 lv lw"
)
rank_one_quartic = (
    a0 * x**4
    + a1 * x**3 * u
    + a2 * x**2 * u**2
    + a3 * x * u**3
    + a4 * u**4
    + x**3 * (lv * v + lw * w)
)
rank_one_quartic_hessian = sp.hessian(rank_one_quartic, (u, v, w, x))
assert sp.expand(rank_one_quartic_hessian.det()) == 0
assert sp.simplify(
    (
        sp.Matrix([0, lw, -lv, 0]).T
        * rank_one_quartic_hessian
    )
) == sp.zeros(1, 4)

lu = sp.symbols("lu")
rank_zero_quartic = x**4 * a0 + x**3 * (lu * u + lv * v + lw * w)
rank_zero_quartic_hessian = sp.hessian(rank_zero_quartic, (u, v, w, x))
assert rank_zero_quartic_hessian.rank(iszerofunc=lambda z: z == 0) <= 2


# Homogeneity makes every projective Hessian-kernel scheme constant over
# the one-dimensional sextic quotient: Hess_W(f)(scale*x,scale*W) has
# weight degree-2.  Check this on every quartic and cubic monomial.
variables = (x, u, v, w)


def homogeneous_monomials(total_degree: int) -> list[sp.Expr]:
    result: list[sp.Expr] = []
    for x_degree in range(total_degree + 1):
        for u_degree in range(total_degree + 1 - x_degree):
            for v_degree in range(
                total_degree + 1 - x_degree - u_degree
            ):
                w_degree = (
                    total_degree - x_degree - u_degree - v_degree
                )
                result.append(
                    x**x_degree
                    * u**u_degree
                    * v**v_degree
                    * w**w_degree
                )
    return result


for total_degree in (3, 4):
    for monomial in homogeneous_monomials(total_degree):
        hessian_w = sp.hessian(monomial, (u, v, w))
        scaled = hessian_w.subs(
            {
                x: scale * x,
                u: scale * u,
                v: scale * v,
                w: scale * w,
            }
        )
        assert (
            scaled - scale ** (total_degree - 2) * hessian_w
        ).applyfunc(sp.expand) == sp.zeros(3)


# The rank-one sextic normal form is the highest-weight binary-root point
# L^6 (a stronger condition than the SIC binary-sextic nullcone L^4 Q).
alpha = sp.symbols("alpha0:4")
linear_form = sum(
    coefficient * variable
    for coefficient, variable in zip(alpha, variables)
)
sextic_power_hessian = sp.hessian(linear_form**6, variables)
for row in range(1, 4):
    for column in range(4):
        assert sp.expand(
            alpha[0] * sextic_power_hessian[row, column]
            - alpha[row] * sextic_power_hessian[0, column]
        ) == 0


print("PASS: rank-one sextic normal form is the constant cone L^6")
print("PASS: lambda^10 makes the ternary quartic W-Hessian singular")
print("PASS: quartic rank two aligns the cubic at lambda^9")
print("PASS: quartic rank one gives the SIC Sym^2 nullcone at lambda^8")
print("PASS: quartic rank zero makes the cubic W-Hessian singular at lambda^7")
print("PASS: one-base homogeneity makes every residual cone constant")
print("SCOPE: this closes the Hess(h_6)-rank-one triple-layer chart")
