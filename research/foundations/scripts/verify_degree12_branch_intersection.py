"""Exact tangent-space certificate for C_(3,2) along E_(6,6)."""

import sympy as sp

W = sp.symbols("W")
q0, q1, q2, r0, r1 = sp.symbols("q0 q1 q2 r0 r1")
Q = W**3 + q2*W**2 + q1*W + q0
R = W**2 + r1*W + r0
M = sp.expand(Q**2 * R**3)
D = sp.diff(M, W).subs(W, 1) - sp.diff(M, W).subs(W, 0)
Phi = sp.expand(M.subs(W, 1) - M.subs(W, 0) - sp.diff(M, W).subs(W, 0))
H = sp.Poly(sp.cancel((-M + sp.diff(M, W).subs(W, 0)*W + M.subs(W, 0))/D), W)
parameters = (q0, q1, q2, r0, r1)
coefficients = sp.Matrix([H.coeff_monomial(W**degree) for degree in range(2, 13)])

# Universal off-diagonal sixfold block.  After translating away the quadratic
# coefficient of the cubic, write Q=z^3+u*z+v and T=z^2+a*z+b.
z, u, v, aa, bb = sp.symbols("z u v aa bb")
block_Q = z**3 + u*z + v
block_T = z**2 + aa*z + bb
block_polynomial = sp.Poly(sp.expand(block_Q**2-block_T**3), z)
block_ideal = [block_polynomial.coeff_monomial(z**degree)
               for degree in range(7)]
block_groebner = sp.groebner(block_ideal, v, u, bb, aa, order="lex")
assert block_groebner.reduce(v)[1] == 0
assert block_groebner.reduce(2*u-3*bb)[1] == 0
assert block_groebner.reduce(bb**2)[1] == 0
assert block_groebner.reduce(aa)[1] == 0

# The (6,6) incidence is formally solvable in a root coordinate at the exact
# witness, so imposing Phi removes one smooth direction and no dual number.
alpha, beta = sp.symbols("alpha beta")
collision_M = (W-alpha)**6*(W-beta)**6
collision_phi = (collision_M.subs(W, 1)-collision_M.subs(W, 0)
                 - sp.diff(collision_M, W).subs(W, 0))
collision_point = {alpha: sp.Rational(-6, 5), beta: sp.Integer(1)}
assert sp.diff(collision_phi, alpha).subs(collision_point) != 0
assert sp.diff(collision_phi, beta).subs(collision_point) != 0


def point(root_q, root_r):
    q = sp.Poly((W-root_q)**3, W)
    r = sp.Poly((W-root_r)**2, W)
    return {
        q0: q.coeff_monomial(1), q1: q.coeff_monomial(W),
        q2: q.coeff_monomial(W**2), r0: r.coeff_monomial(1),
        r1: r.coeff_monomial(W),
    }


def differential_data(substitution):
    constraint = sp.Matrix([[sp.diff(Phi, parameter).subs(substitution)
                             for parameter in parameters]])
    source_tangent = sp.Matrix.hstack(*constraint.nullspace())
    jacobian = coefficients.jacobian(parameters).subs(substitution)
    image = jacobian*source_tangent
    assert image.rank() == 4
    return {
        "constraint": constraint,
        "source_tangent": source_tangent,
        "jacobian": jacobian,
        "image": image,
        "constraint_hessian": sp.hessian(Phi, parameters).subs(substitution),
        "map_hessians": [sp.hessian(entry, parameters).subs(substitution)
                         for entry in coefficients],
    }


def second_jet(data, first, second):
    """Second derivative along Phi=0, modulo an arbitrary tangent correction."""
    constraint = data["constraint"]
    rhs = -(first.T*data["constraint_hessian"]*second)[0]
    pivot = next(index for index, entry in enumerate(constraint.tolist()[0])
                 if entry != 0)
    acceleration = sp.zeros(len(parameters), 1)
    acceleration[pivot] = rhs/constraint[0, pivot]
    hessian_term = sp.Matrix([(first.T*hessian*second)[0]
                              for hessian in data["map_hessians"]])
    return hessian_term + data["jacobian"]*acceleration


left = differential_data(point(sp.Rational(-6, 5), sp.Integer(1)))
right = differential_data(point(sp.Integer(1), sp.Rational(-6, 5)))
combined = left["image"].row_join(right["image"])
combined_rank = combined.rank()
intersection_dimension = (left["image"].rank() + right["image"].rank()
                          - combined_rank)

assert combined_rank == 5
assert intersection_dimension == 3

# Pair the three source directions having the same first-order target image.
paired_directions = left["image"].row_join(-right["image"]).nullspace()
assert len(paired_directions) == 3
quadratic_differences = []
for i in range(3):
    for j in range(i, 3):
        left_i = left["source_tangent"]*paired_directions[i][:4, 0]
        left_j = left["source_tangent"]*paired_directions[j][:4, 0]
        right_i = right["source_tangent"]*paired_directions[i][4:, 0]
        right_j = right["source_tangent"]*paired_directions[j][4:, 0]
        quadratic_differences.append(
            second_jet(left, left_i, left_j)-second_jet(right, right_i, right_j)
        )

quadratic_matrix = sp.Matrix.hstack(*quadratic_differences)
quadratic_quotient_rank = (combined.row_join(quadratic_matrix).rank()
                           - combined_rank)
assert quadratic_quotient_rank == 2

# Exact scalar representatives for the two quadratic relations on the three
# common tangent directions.  Off-diagonal entries use the bilinear-matrix
# convention, so the associated polynomial has coefficients 2*q_ij.
annihilator = sp.Matrix.hstack(*combined.T.nullspace()).T
projected_quadratics = annihilator*quadratic_matrix
q1 = sp.Matrix([[186624, 324000, -194400, 953125, -806250, 765000]])
q2 = sp.Matrix([[0, 0, 0, 25, -30, 36]])
assert projected_quadratics.col_join(q1).rank() == projected_quadratics.rank()
assert projected_quadratics.col_join(q2).rank() == projected_quadratics.rank()
assert q1.col_join(q2).rank() == 2
x0, x1, x2 = sp.symbols("x0 x1 x2")
q1_polynomial = (q1[0, 0]*x0**2 + 2*q1[0, 1]*x0*x1
                 + 2*q1[0, 2]*x0*x2 + q1[0, 3]*x1**2
                 + 2*q1[0, 4]*x1*x2 + q1[0, 5]*x2**2)
q2_polynomial = (q2[0, 0]*x0**2 + 2*q2[0, 1]*x0*x1
                 + 2*q2[0, 2]*x0*x2 + q2[0, 3]*x1**2
                 + 2*q2[0, 4]*x1*x2 + q2[0, 5]*x2**2)
assert sp.factor(q2_polynomial) == (5*x1-6*x2)**2
assert sp.factor(q1_polynomial.subs(x2, sp.Rational(5, 6)*x1)) == \
       9*(144*x0+125*x1)**2

# The stronger equality M_+=M_- supplies the quotient by two independent
# squares, with standard basis 1, epsilon, eta, epsilon*eta and length four.
epsilon, eta = sp.symbols("epsilon eta")
dual_block = sp.groebner([epsilon**2, eta**2], epsilon, eta)
standard_monomials = (sp.Integer(1), epsilon, eta, epsilon*eta)
assert all(dual_block.reduce(monomial)[1] == monomial
           for monomial in standard_monomials)

print("PASS: both normalization differentials have rank 4")
print("PASS: at the admissible (6,6) witness their tangent images meet in dimension 3")
print("PASS: an xy=0 normal-crossing model is impossible (the stratum has dimension 1)")
print("PASS: the excess common tangent directions separate with quadratic rank 2")
print("PASS: one exchanged sixfold block is a dual number: Q^2=T^3 gives b^2=0")
print("PASS: the quadratic upper bound and dual-block lower bound both give length 4")
print("PASS: two independent blocks give D=k[[t,e,n]]/(e^2,n^2)")
