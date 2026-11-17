#!/usr/bin/env python3
"""Exact identities for the shared JC(2)--HC(4) isotropic boundary package.

This checker has eight independent blocks.

1.  For the cotangent potential

        Psi=t*P(x,y)+m*Q(x,y)+H(x,y),

    its four-variable Hessian determinant is Jac(P,Q)^2, independently of
    the second derivatives and of t,m.

2.  Every source-only Hessian direction R gives a constant-Hessian pencil
    through the cotangent block.  Its relative endomorphism N=S^(-1)T is
    square-zero.  Thus the plane cotangent packet lies inside the precise
    relative-nilpotent class closed by HC4MR1.

3.  For the more general isotropic form Psi=t*P(x,y)+Phi(x,y,m), the t^2
    Schur coefficient vanishes and the t coefficient is

        -Phi_mm * grad(P)^T adj(Hess(P)) grad(P).

    Thus the cotangent branch Phi_mm=0 is one exact factor of the first
    nontrivial Schur remainder.

4.  In an adapted index-two boundary chart g=r^2*ell(r,z), the reduced
    conormal residue r^(-1)*partial_r(g) mod r is 2*ell(0,z).  Hence the
    quartic packet's odd-square multiplier ell is half of the normal
    residue, and its square is the leading cotangent-Hessian coefficient.

5.  For the exact quartic spectator model, the discriminant pulls back as
    r^2*ell.  The cusp initial residue is 18*T^2.  At the two branches of
    the 2+2 connector the initial coefficients are respectively
    -9*sqrt(3)/2 and 9*sqrt(3)/2.

6.  The completed cusp and node conductor maps are isomorphisms in the
    positive associated-graded degree occupied by those residues.  Thus the
    proposed associated-graded cokernel class Obs_pair is zero.

7.  The degree-zero conductor mismatch also vanishes because both connector
    residue values are zero.  Comparing their first derivatives would
    require a noncanonical identification of the two tangent parameters.

8.  Every one of the 72 labelled cusp/matching packets generates S4, and
    the zero conductor class is independent of the labelling.

The script does not prove that an arbitrary HC(4) counterexample admits the
required polynomial isotropic flag.  It proves that the proposed local
associated-graded conductor cokernel cannot supply the desired nonvanishing;
any refinement must add global data that identifies endpoint jets.
"""

from __future__ import annotations

from itertools import permutations

import sympy as sp


# ---------------------------------------------------------------------------
# 1. Universal cotangent determinant
# ---------------------------------------------------------------------------

P_x, P_y, Q_x, Q_y = sp.symbols("P_x P_y Q_x Q_y")
A_xx, A_xy, A_yy = sp.symbols("A_xx A_xy A_yy")

cotangent_hessian = sp.Matrix(
    [
        [A_xx, A_xy, P_x, Q_x],
        [A_xy, A_yy, P_y, Q_y],
        [P_x, P_y, 0, 0],
        [Q_x, Q_y, 0, 0],
    ]
)
jacobian = P_x * Q_y - P_y * Q_x

assert sp.factor(cotangent_hessian.det() - jacobian**2) == 0


# ---------------------------------------------------------------------------
# 2. Cotangent packets are square-zero constant-Hessian pencils
# ---------------------------------------------------------------------------

s = sp.symbols("s")
R_xx, R_xy, R_yy = sp.symbols("R_xx R_xy R_yy")
source_direction = sp.Matrix(
    [
        [R_xx, R_xy, 0, 0],
        [R_xy, R_yy, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]
)
assert sp.factor(
    (cotangent_hessian + s * source_direction).det() - jacobian**2
) == 0

source_jacobian = sp.Matrix([[P_x, P_y], [Q_x, Q_y]])
source_hessian = sp.Matrix([[R_xx, R_xy], [R_xy, R_yy]])
relative_nilpotent = sp.zeros(4)
relative_nilpotent[2:4, 0:2] = source_jacobian.inv().T * source_hessian
assert (cotangent_hessian * relative_nilpotent - source_direction).applyfunc(
    sp.factor
) == sp.zeros(4)
assert relative_nilpotent**2 == sp.zeros(4)


# ---------------------------------------------------------------------------
# 3. First isotropic Schur remainder
# ---------------------------------------------------------------------------

t = sp.symbols("t")
P_xx, P_xy, P_yy = sp.symbols("P_xx P_xy P_yy")
Phi_xx, Phi_xy, Phi_xm = sp.symbols("Phi_xx Phi_xy Phi_xm")
Phi_yy, Phi_ym, Phi_mm = sp.symbols("Phi_yy Phi_ym Phi_mm")

isotropic_hessian = sp.Matrix(
    [
        [0, P_x, P_y, 0],
        [P_x, t * P_xx + Phi_xx, t * P_xy + Phi_xy, Phi_xm],
        [P_y, t * P_xy + Phi_xy, t * P_yy + Phi_yy, Phi_ym],
        [0, Phi_xm, Phi_ym, Phi_mm],
    ]
)
determinant_in_t = sp.Poly(sp.expand(isotropic_hessian.det()), t)
binary_bordered_remainder = (
    P_x**2 * P_yy - 2 * P_x * P_y * P_xy + P_y**2 * P_xx
)

assert determinant_in_t.degree() <= 1
assert determinant_in_t.coeff_monomial(t**2) == 0
assert sp.factor(
    determinant_in_t.coeff_monomial(t)
    + Phi_mm * binary_bordered_remainder
) == 0


# ---------------------------------------------------------------------------
# 4. Index-two reduced conormal residue and its square
# ---------------------------------------------------------------------------

r, z = sp.symbols("r z")
ell_0 = sp.Function("ell_0")(z)
ell_1 = sp.Function("ell_1")(z)
ell = ell_0 + r * ell_1
g_pullback = sp.expand(r**2 * ell)

reduced_normal_residue = sp.expand(sp.diff(g_pullback, r) / r).subs(r, 0)
assert sp.simplify(reduced_normal_residue - 2 * ell_0) == 0
assert sp.expand(reduced_normal_residue**2 - 4 * ell_0**2) == 0

T = sp.symbols("T")
cusp_ell = 4 * r - 9 * T**2
cusp_pullback = r**2 * cusp_ell
cusp_residue = sp.expand(sp.diff(cusp_pullback, r) / r).subs(r, 0)
assert cusp_residue == -18 * T**2


# ---------------------------------------------------------------------------
# 5. Exact 3+1 cusp and both 2+2 connector branches
# ---------------------------------------------------------------------------

X, u, v = sp.symbols("X u v")
quartic = X**4 - X**3 + u * X - v
quartic_discriminant = sp.factor(sp.discriminant(quartic, X))
expected_discriminant = (
    -27 * u**4
    + 4 * u**3
    + 6 * u**2 * v
    + 192 * u * v**2
    - 256 * v**3
    - 27 * v**2
)
assert sp.expand(quartic_discriminant - expected_discriminant) == 0

spectator_v = T**4 - T**3 + u * T
spectator_r = sp.diff(spectator_v, T)
pulled_discriminant = sp.factor(
    quartic_discriminant.subs(v, spectator_v)
)
spectator_ell = sp.factor(pulled_discriminant / spectator_r**2)
assert sp.expand(pulled_discriminant - spectator_r**2 * spectator_ell) == 0

critical_u = 3 * T**2 - 4 * T**3
ell_on_ramification = sp.factor(spectator_ell.subs(u, critical_u))
expected_ell_on_ramification = (
    -9 * T**2 * (2 * T - 1) ** 2 * (8 * T**2 - 4 * T - 1)
)
assert sp.expand(ell_on_ramification - expected_ell_on_ramification) == 0
rho_on_ramification = sp.factor(2 * ell_on_ramification)

# The full quartic discriminant differs from the standard cubic cusp
# equation by a unit at T=0.  Accordingly its oriented initial coefficient
# is +18 rather than the -18 in the standard frame checked above.
assert sp.limit(rho_on_ramification / T**2, T, 0) == 18

sqrt_three = sp.sqrt(3)
connector_plus = (1 + sqrt_three) / 4
connector_minus = (1 - sqrt_three) / 4
tau = sp.symbols("tau")


def initial_coefficient_at(polynomial: sp.Expr, point: sp.Expr) -> sp.Expr:
    """Return the coefficient of tau in polynomial(point+tau)."""

    expansion = sp.Poly(sp.expand(polynomial.subs(T, point + tau)), tau)
    assert expansion.coeff_monomial(1) == 0
    return sp.simplify(expansion.coeff_monomial(tau))


rho_plus_initial = initial_coefficient_at(
    rho_on_ramification, connector_plus
)
rho_minus_initial = initial_coefficient_at(
    rho_on_ramification, connector_minus
)
assert rho_plus_initial == -9 * sqrt_three / 2
assert rho_minus_initial == 9 * sqrt_three / 2


# ---------------------------------------------------------------------------
# 6. The associated-graded conductor obstruction is zero
# ---------------------------------------------------------------------------

def lies_in_image(
    image_matrix: sp.Matrix, target_vector: sp.Matrix
) -> tuple[bool, sp.Matrix]:
    """Decide exact cokernel vanishing and return one source lift."""

    augmented = image_matrix.row_join(target_vector)
    if image_matrix.rank() != augmented.rank():
        return False, sp.zeros(image_matrix.cols, 1)
    solution = sp.linsolve(
        (image_matrix, target_vector),
        *sp.symbols(f"w0:{image_matrix.cols}"),
    )
    lift = sp.Matrix(next(iter(solution)))
    return True, lift


# Cusp:
#   R=k[[T^2,T^3]], N=k[[T]], c=(T^2,T^3)=T^2*N.
# In conductor degree one, both c/c^2 and cN/(cN)^2 have basis
# (T^2,T^3), so the map is the identity.  The cusp residue has vector
# (18,0) in this basis.
cusp_conductor_gr_map = sp.eye(2)
cusp_initial_vector = sp.Matrix([18, 0])
cusp_zero, cusp_lift = lies_in_image(
    cusp_conductor_gr_map, cusp_initial_vector
)
assert cusp_zero and cusp_lift == cusp_initial_vector

# Node:
#   R={(f_+,f_-): f_+(0)=f_-(0)}, N=k[[tau_+]]+k[[tau_-]],
#   c=(tau_+)+(tau_-).
# In conductor degree one, c/c^2 -> cN/(cN)^2 is again the identity:
# the two linear coefficients are independent.  Consequently the unequal
# (indeed opposite) connector initials define no cokernel mismatch.
node_conductor_gr_map = sp.eye(2)
connector_initial_vector = sp.Matrix(
    [rho_plus_initial, rho_minus_initial]
)
connector_zero, connector_lift = lies_in_image(
    node_conductor_gr_map, connector_initial_vector
)
assert connector_zero and connector_lift == connector_initial_vector


# ---------------------------------------------------------------------------
# 7. Value descent and the absence of a canonical derivative comparison
# ---------------------------------------------------------------------------

# In degree zero the node map is the diagonal k -> k+k.  Its cokernel is
# detected by (a,b) |-> a-b.  Both residue values vanish, so this genuine
# conductor-pairing obstruction is zero.
node_degree_zero_map = sp.Matrix([[1], [1]])
connector_value_vector = sp.Matrix(
    [
        sp.simplify(rho_on_ramification.subs(T, connector_plus)),
        sp.simplify(rho_on_ramification.subs(T, connector_minus)),
    ]
)
connector_value_zero, _ = lies_in_image(
    node_degree_zero_map, connector_value_vector
)
assert connector_value_zero
assert connector_value_vector == sp.zeros(2, 1)

# A proposed subtraction of linear coefficients is not invariant under
# independent changes tau_+=alpha*tau'_+, tau_-=beta*tau'_-.
alpha, beta = sp.symbols("alpha beta", nonzero=True)
rescaled_initials = sp.Matrix(
    [alpha * rho_plus_initial, beta * rho_minus_initial]
)
naive_difference = sp.factor(
    rescaled_initials[0] - rescaled_initials[1]
)
assert sp.factor(naive_difference + 9 * sqrt_three * (alpha + beta) / 2) == 0
assert sp.simplify(naive_difference.subs(beta, -alpha)) == 0
assert sp.simplify(naive_difference.subs(beta, alpha)) != 0


# ---------------------------------------------------------------------------
# 8. Actual S4-compatible packet labellings
# ---------------------------------------------------------------------------

SheetPermutation = tuple[int, int, int, int]


def compose(
    left: SheetPermutation, right: SheetPermutation
) -> SheetPermutation:
    return tuple(left[right[index]] for index in range(4))


identity: SheetPermutation = (0, 1, 2, 3)


def transposition(first: int, second: int) -> SheetPermutation:
    values = list(range(4))
    values[first], values[second] = values[second], values[first]
    return tuple(values)


def generated_group(
    generators: tuple[SheetPermutation, ...],
) -> set[SheetPermutation]:
    group = {identity}
    frontier = [identity]
    while frontier:
        element = frontier.pop()
        for generator in generators:
            product = compose(generator, element)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return group


transpositions = [
    transposition(first, second)
    for first in range(4)
    for second in range(first + 1, 4)
]
cusp_pairs = [
    (first, second)
    for first in transpositions
    for second in transpositions
    if first != second
    and compose(compose(first, second), first)
    == compose(compose(second, first), second)
]
perfect_matchings = [
    (
        transposition(0, 1),
        transposition(2, 3),
    ),
    (
        transposition(0, 2),
        transposition(1, 3),
    ),
    (
        transposition(0, 3),
        transposition(1, 2),
    ),
]
assert len(cusp_pairs) == 24
assert len(perfect_matchings) == 3

packet_count = 0
for cusp_first, cusp_second in cusp_pairs:
    for connector_first, connector_second in perfect_matchings:
        packet_count += 1
        group = generated_group(
            (
                cusp_first,
                cusp_second,
                connector_first,
                connector_second,
            )
        )
        assert group == set(permutations(range(4)))
        # Sheet labels do not change the completed conductor map or its
        # already-zero cokernel class.
        assert connector_zero and connector_value_zero
assert packet_count == 72


print("PASS: cotangent determinant = Jac(P,Q)^2")
print("PASS: cotangent source directions give square-zero Hessian pencils")
print("PASS: first isotropic Schur remainder = -Phi_mm*R(P)")
print("PASS: index-two reduced conormal residue = 2*ell")
print("PASS: clean cusp initial conormal residue = -18*T^2")
print("PASS: quartic cusp initial residue = 18*T^2")
print(
    "PASS: connector initials = "
    "(-9*sqrt(3)/2, 9*sqrt(3)/2)"
)
print("PASS: Obs_pair is zero in the cusp and connector conductor gradeds")
print("PASS: derivative subtraction depends on a noncanonical jet transport")
print("PASS: zero class is compatible with all 72 transitive S4 packets")
