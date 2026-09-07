#!/usr/bin/env python3
"""Verify the exact algebra in the rank-two triple-layer HC(4) reduction.

For psi=q_2+h_3+h_4+h_6 with rank Hess(h_6)=2, choose the constant
two-plane U=(t,m) in the sextic kernel and quotient variables X=(x,y).
The top determinant faces first make Hess_U(h_4) singular.

The quartic cone-degree lemma leaves a constant kernel or the moving
square c*(X^T M U)^2.  In the moving case the next odd face forces the
dual-degree-at-least-two cubic part to be

    (X^T M U)*(alpha*t+beta*m).

This checker verifies that the later dual-degree-four determinant face is
still 48*c^4*det(M)^2*(X^T M U)^4 and cannot cancel.  If Hess_U(h_4)=0,
the next face makes Hess_U(h_3) singular; cubic homogeneity admits no
moving projective cone, so a constant common direction results.
"""

from __future__ import annotations

import contextlib
import io
from itertools import permutations
from pathlib import Path
import runpy

import sympy as sp


# Replay the quartic cone-degree lemma and the common-direction reduction.
for shared_name in (
    "verify_hc4_even_quartic_sextic_closure.py",
    "verify_hc4_meng_triple_rank_three_reduction.py",
):
    shared_path = Path(__file__).with_name(shared_name)
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


# 1. The top rank-two faces.  Variables are ordered as U=(t,m), X=(x,y).
c1, c2 = sp.symbols("c1 c2")
C6 = sp.diag(0, 0, c1, c2)
H0 = generic_symmetric("h")
A3 = generic_symmetric("a")
B4 = generic_symmetric("b")
spatial = sp.Poly(
    (H0 + lam * A3 + lam**2 * B4 + lam**4 * C6).det(
        method="berkowitz"
    ),
    lam,
)
quotient_determinant = c1 * c2
quartic_kernel_determinant = B4[:2, :2].det()
mixed_binary_determinant = (
    B4[0, 0] * A3[1, 1]
    + B4[1, 1] * A3[0, 0]
    - 2 * B4[0, 1] * A3[0, 1]
)
assert sp.expand(
    spatial.coeff_monomial(lam**12)
    - quotient_determinant * quartic_kernel_determinant
) == 0
assert sp.expand(
    spatial.coeff_monomial(lam**11)
    - quotient_determinant * mixed_binary_determinant
) == 0

# If the quartic UU block is zero, the lambda^10 face is clean even with
# arbitrary quartic cross and quotient blocks.
zero_quartic_uu = B4.subs(
    {B4[0, 0]: 0, B4[0, 1]: 0, B4[1, 1]: 0}
)
zero_uu_spatial = sp.Poly(
    (
        H0
        + lam * A3
        + lam**2 * zero_quartic_uu
        + lam**4 * C6
    ).det(method="berkowitz"),
    lam,
)
assert sp.expand(
    zero_uu_spatial.coeff_monomial(lam**10)
    - quotient_determinant * A3[:2, :2].det()
) == 0


# 2. Normalize an invertible M to the identity, so L=x*t+y*m and the
# moving quartic UU kernel is k=(-y,x).  Solve the degree-eleven condition
# k^T Hess_U(h_3) k=0 on all twenty cubic monomials.
x, y, t, m, epsilon, scalar = sp.symbols(
    "x y t m epsilon scalar"
)
variables = (t, m, x, y)
L = x * t + y * m
kernel = sp.Matrix([-y, x])

cubic_monomials: list[sp.Expr] = []
for x_degree in range(4):
    for y_degree in range(4 - x_degree):
        for t_degree in range(4 - x_degree - y_degree):
            m_degree = 3 - x_degree - y_degree - t_degree
            cubic_monomials.append(
                x**x_degree
                * y**y_degree
                * t**t_degree
                * m**m_degree
            )
cubic_coefficients = sp.symbols(f"u0:{len(cubic_monomials)}")
generic_cubic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        cubic_coefficients, cubic_monomials
    )
)
cubic_uu = sp.hessian(generic_cubic, (t, m))
odd_face = sp.Poly(
    sp.expand((kernel.T * cubic_uu * kernel)[0]),
    x,
    y,
    t,
    m,
)
odd_equations = [coefficient for _, coefficient in odd_face.terms()]
odd_matrix = sp.zeros(len(odd_equations), len(cubic_coefficients))
for row, equation in enumerate(odd_equations):
    for column, coefficient in enumerate(cubic_coefficients):
        odd_matrix[row, column] = sp.diff(equation, coefficient)
assert odd_matrix.rank() == 8

solution = next(iter(sp.linsolve(odd_equations, cubic_coefficients)))
restricted_cubic = sp.expand(
    generic_cubic.subs(dict(zip(cubic_coefficients, solution)))
)


high_dual_cubic = sum(
    coefficient
    * x ** exponents[0]
    * y ** exponents[1]
    * t ** exponents[2]
    * m ** exponents[3]
    for exponents, coefficient in sp.Poly(
        restricted_cubic, x, y, t, m
    ).terms()
    if exponents[2] + exponents[3] >= 2
)
alpha, beta = sp.symbols("alpha beta")
aligned_template = sp.expand(L * (alpha * t + beta * m))
aligned_coefficients = sp.solve(
    sp.Poly(
        sp.expand(high_dual_cubic - aligned_template),
        x,
        y,
        t,
        m,
    ).coeffs(),
    (alpha, beta),
    dict=True,
)
assert aligned_coefficients


# 3. The lambda^8 / dual-degree-four face after cubic alignment.  Besides
# the moving square, allow every Hessian block compatible with a term of
# dual degree at most one in both h_4 and h_3.
p, q = sp.symbols("p q")
aligned_cubic = L * (p * t + q * m)
moving_quartic = scalar * L**2
aligned_hessian = sp.hessian(aligned_cubic, variables)
moving_hessian = sp.hessian(moving_quartic, variables)


def symmetric_two(entries: tuple[sp.Expr, sp.Expr, sp.Expr]) -> sp.Matrix:
    return sp.Matrix(
        [[entries[0], entries[1]], [entries[1], entries[2]]]
    )


def lower_dual_hessian(prefix: str) -> sp.Matrix:
    ux_entries = sp.symbols(f"{prefix}u0:4")
    c0_entries = sp.symbols(f"{prefix}c0:3")
    ct_entries = sp.symbols(f"{prefix}t0:3")
    cm_entries = sp.symbols(f"{prefix}m0:3")
    result = sp.zeros(4)
    result[:2, 2:] = sp.Matrix(2, 2, ux_entries)
    result[2:, :2] = result[:2, 2:].T
    result[2:, 2:] = (
        symmetric_two(c0_entries)
        + t * symmetric_two(ct_entries)
        + m * symmetric_two(cm_entries)
    )
    return result


quartic_hessian = moving_hessian + lower_dual_hessian("b")
cubic_hessian = aligned_hessian + lower_dual_hessian("a")
sextic_entries = sp.symbols("s0:3")
sextic_hessian = sp.zeros(4)
sextic_hessian[2:, 2:] = symmetric_two(sextic_entries)
quadratic_hessian = generic_symmetric("q")

full_pencil = (
    quadratic_hessian
    + lam * cubic_hessian
    + lam**2 * quartic_hessian
    + lam**4 * sextic_hessian
).subs({t: epsilon * t, m: epsilon * m})


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


target_face = 0
for permutation in permutations(range(4)):
    term: sp.Expr = permutation_sign(permutation)
    for row, column in enumerate(permutation):
        term *= full_pencil[row, column]
    target_face += sp.Poly(
        sp.expand(term), lam, epsilon
    ).coeff_monomial(lam**8 * epsilon**4)

expected_face = 48 * scalar**4 * L**4
assert sp.expand(target_face - expected_face) == 0


# 4. If Hess_U(h_4)=0, the binary singular-Hessian form of h_3 cannot
# have a moving projective cone.  If d>=2 is its highest nonlinear dual
# degree and delta>=1 is the cone degree, the scalar source degree would
# be 3-d-d*delta, which is always negative.
moving_cubic_degrees = [
    (dual, cone, 3 - dual - dual * cone)
    for dual in range(2, 4)
    for cone in range(1, 4)
    if 3 - dual - dual * cone >= 0
]
assert moving_cubic_degrees == []


print("PASS: rank-two top faces give det(B_U)=0 and its A_U polarization")
print("PASS: when B_U=0, the next face is det(A_U)")
print("PASS: the moving quartic cone forces an aligned high-dual cubic")
print("PASS: the aligned cubic cannot change the uncancellable lambda^8 face")
print("PASS: a cubic binary cone cannot move by total-degree counting")
print("SCOPE: this excludes the full Hess(h_6)-rank-two triple-layer chart")
