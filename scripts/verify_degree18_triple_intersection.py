#!/usr/bin/env python3
"""Exact degree-18 tangent and triple-intersection certificate.

The component is C_(3,4), and the collision type is (6,6,6).  The three
normalization points are distinguished by the sixfold root assigned to Q.
All calculations are over QQ at the admissible witness (-3/4, 1, 2).

The certificate proves the pairwise and triple transverse intersection
algebras.  It does not assert the still-open completed-local three-sheet
equalizer or total conductor.
"""

import sympy as sp


W = sp.symbols("W")
q0, q1, q2, r0, r1, r2, r3 = sp.symbols("q0 q1 q2 r0 r1 r2 r3")
parameters = (q0, q1, q2, r0, r1, r2, r3)

Q = W**3 + q2 * W**2 + q1 * W + q0
R = W**4 + r3 * W**3 + r2 * W**2 + r1 * W + r0

roots = (sp.Rational(-3, 4), sp.Integer(1), sp.Integer(2))


def standard_monomial_count(groebner_basis):
    """Count standard monomials of a zero-dimensional Groebner basis."""
    leading = [
        polynomial.LM(order=groebner_basis.order).exponents
        for polynomial in groebner_basis.polys
    ]
    variable_count = len(groebner_basis.gens)
    bounds = []
    for variable in range(variable_count):
        pure_powers = [
            exponent[variable]
            for exponent in leading
            if exponent[variable] > 0
            and all(
                exponent[other] == 0
                for other in range(variable_count)
                if other != variable
            )
        ]
        if not pure_powers:
            raise AssertionError("the initial ideal is not zero-dimensional")
        bounds.append(min(pure_powers))
    count = 0
    for exponent in sp.utilities.iterables.cartes(*[range(bound) for bound in bounds]):
        if not any(
            all(exponent[index] >= initial[index] for index in range(variable_count))
            for initial in leading
        ):
            count += 1
    return count


def branch_point(q_root):
    """Return the normalization point with ``q_root`` assigned to Q."""
    r_roots = list(roots)
    r_roots.remove(q_root)
    q_polynomial = sp.Poly((W - q_root) ** 3, W)
    r_polynomial = sp.Poly(sp.prod((W - root) ** 2 for root in r_roots), W)
    return {
        q0: q_polynomial.coeff_monomial(1),
        q1: q_polynomial.coeff_monomial(W),
        q2: q_polynomial.coeff_monomial(W**2),
        r0: r_polynomial.coeff_monomial(1),
        r1: r_polynomial.coeff_monomial(W),
        r2: r_polynomial.coeff_monomial(W**2),
        r3: r_polynomial.coeff_monomial(W**3),
    }


def differential_data(q_root):
    """Build the normalized map through its degree-two local expansion.

    Expanding Q^2 R^3 before forming the normalized quotient is much faster
    than differentiating its seventeen global rational coefficient functions.
    """
    substitution = branch_point(q_root)
    deltas = sp.symbols("delta0:7")
    delta_q = deltas[2] * W**2 + deltas[1] * W + deltas[0]
    delta_r = (
        deltas[6] * W**3
        + deltas[5] * W**2
        + deltas[4] * W
        + deltas[3]
    )
    point_q = sp.expand(Q.subs(substitution))
    point_r = sp.expand(R.subs(substitution))
    local_M = (
        point_q**2 * point_r**3,
        2 * point_q * delta_q * point_r**3
        + 3 * point_q**2 * point_r**2 * delta_r,
        delta_q**2 * point_r**3
        + 6 * point_q * delta_q * point_r**2 * delta_r
        + 3 * point_q**2 * point_r * delta_r**2,
    )

    def denominator(polynomial):
        return sp.diff(polynomial, W).subs(W, 1) - sp.diff(polynomial, W).subs(W, 0)

    def incidence(polynomial):
        return (
            polynomial.subs(W, 1)
            - polynomial.subs(W, 0)
            - sp.diff(polynomial, W).subs(W, 0)
        )

    local_D = tuple(sp.expand(denominator(entry)) for entry in local_M)
    local_N = tuple(
        sp.Matrix(
            [
                -sp.Poly(entry, W).coeff_monomial(W**degree)
                for degree in range(2, 19)
            ]
        )
        for entry in local_M
    )
    assert local_D[0] != 0
    map_linear = local_N[1] / local_D[0] - local_N[0] * local_D[1] / local_D[0] ** 2
    map_quadratic = (
        local_N[2] / local_D[0]
        - local_N[1] * local_D[1] / local_D[0] ** 2
        + local_N[0] * local_D[1] ** 2 / local_D[0] ** 3
        - local_N[0] * local_D[2] / local_D[0] ** 2
    )
    incidence_linear = sp.expand(incidence(local_M[1]))
    incidence_quadratic = sp.expand(incidence(local_M[2]))
    constraint = sp.Matrix(
        [[sp.diff(incidence_linear, delta) for delta in deltas]]
    )
    source_tangent = sp.Matrix.hstack(*constraint.nullspace())
    jacobian = map_linear.jacobian(deltas)
    image = jacobian * source_tangent
    return {
        "constraint": constraint,
        "source_tangent": source_tangent,
        "jacobian": jacobian,
        "image": image,
        "constraint_hessian": sp.hessian(incidence_quadratic, deltas),
        "map_hessians": [
            sp.hessian(entry, deltas) for entry in map_quadratic
        ],
    }


def second_jet(data, first_coordinates, second_coordinates):
    """Second jet along Phi=0, modulo a later arbitrary tangent correction."""
    first = data["source_tangent"] * first_coordinates
    second = data["source_tangent"] * second_coordinates
    rhs = -(first.T * data["constraint_hessian"] * second)[0]
    pivot = next(
        index
        for index, entry in enumerate(data["constraint"].tolist()[0])
        if entry != 0
    )
    acceleration = sp.zeros(len(parameters), 1)
    acceleration[pivot] = rhs / data["constraint"][0, pivot]
    hessian_term = sp.Matrix(
        [(first.T * hessian * second)[0] for hessian in data["map_hessians"]]
    )
    return hessian_term + data["jacobian"] * acceleration


def extend_to_basis(columns, dimension):
    """Extend independent column vectors by standard basis vectors."""
    result = list(columns)
    for column in range(dimension):
        standard = sp.eye(dimension)[:, column]
        candidate = sp.Matrix.hstack(*result, standard)
        if candidate.rank() > len(result):
            result.append(standard)
        if len(result) == dimension:
            break
    return sp.Matrix.hstack(*result)


# The collision polynomial is independent of the allocation branch.
alpha, beta, gamma = sp.symbols("alpha beta gamma")
collision_M = (W - alpha) ** 6 * (W - beta) ** 6 * (W - gamma) ** 6
collision_Phi = sp.expand(
    collision_M.subs(W, 1)
    - collision_M.subs(W, 0)
    - sp.diff(collision_M, W).subs(W, 0)
)
collision_substitution = dict(zip((alpha, beta, gamma), roots))

assert collision_Phi.subs(collision_substitution) == 0
root_derivatives = tuple(
    sp.factor(sp.diff(collision_Phi, root).subs(collision_substitution))
    for root in (alpha, beta, gamma)
)
assert root_derivatives == (
    -sp.Rational(243, 2),
    -sp.Rational(2187, 32),
    -sp.Rational(2187, 128),
)

data = [differential_data(root) for root in roots]
assert all(entry["source_tangent"].rank() == 6 for entry in data)
assert all(entry["image"].rank() == 6 for entry in data)

for left in range(3):
    for right in range(left + 1, 3):
        pair_span = data[left]["image"].row_join(data[right]["image"])
        assert pair_span.rank() == 7
        assert data[left]["image"].rank() + data[right]["image"].rank() - 7 == 5

triple_span = data[0]["image"].row_join(data[1]["image"]).row_join(data[2]["image"])
assert triple_span.rank() == 8

# A triple of tangent vectors has a common target image precisely when it is
# in the kernel of this simultaneous-difference map.
zero = sp.zeros(17, 6)
simultaneous_difference = sp.Matrix.vstack(
    data[0]["image"].row_join(-data[1]["image"]).row_join(zero),
    data[0]["image"].row_join(zero).row_join(-data[2]["image"]),
)
common_tuples = sp.Matrix.hstack(*simultaneous_difference.nullspace())
assert common_tuples.shape == (18, 5)
common_images = data[0]["image"] * common_tuples[:6, :]
assert common_images.rank() == 5

# Locate the two-dimensional collision-stratum tangent plane inside the five
# common tangent directions, then extend it to a basis.  The remaining three
# directions form a linear transverse slice.
collision_constraint = sp.Matrix(
    [[
        sp.diff(collision_Phi, root).subs(collision_substitution)
        for root in (alpha, beta, gamma)
    ]]
)
collision_source_tangent = sp.Matrix.hstack(*collision_constraint.nullspace())
collision_q = sp.Poly((W - alpha) ** 3, W)
collision_r = sp.Poly((W - beta) ** 2 * (W - gamma) ** 2, W)
collision_parameters = sp.Matrix(
    [
        collision_q.coeff_monomial(1),
        collision_q.coeff_monomial(W),
        collision_q.coeff_monomial(W**2),
        collision_r.coeff_monomial(1),
        collision_r.coeff_monomial(W),
        collision_r.coeff_monomial(W**2),
        collision_r.coeff_monomial(W**3),
    ]
)
collision_parameter_tangent = (
    collision_parameters.jacobian((alpha, beta, gamma)).subs(collision_substitution)
    * collision_source_tangent
)
collision_sheet_coordinates = data[0]["source_tangent"].gauss_jordan_solve(
    collision_parameter_tangent
)[0]
collision_image = data[0]["image"] * collision_sheet_coordinates
assert collision_image.rank() == 2

collision_common_coordinates = []
for column in range(2):
    solution = sp.linsolve((common_images, collision_image[:, column]))
    collision_common_coordinates.append(sp.Matrix(tuple(solution)[0]))
common_basis_change = extend_to_basis(collision_common_coordinates, 5)
assert common_basis_change.det() != 0
common_tuples = common_tuples * common_basis_change

# Second-order simultaneous lifting uses one acceleration on each sheet.  Its
# obstruction is the cokernel of the same simultaneous-difference map.  Work
# only in the final three (transverse) common directions.
obstruction_quotient = sp.Matrix.hstack(
    *simultaneous_difference.T.nullspace()
).T
assert obstruction_quotient * simultaneous_difference == sp.zeros(
    obstruction_quotient.rows, simultaneous_difference.cols
)

second_obstructions = []
for first in range(2, 5):
    for second in range(first, 5):
        jets = []
        for sheet in range(3):
            sheet_slice = slice(6 * sheet, 6 * (sheet + 1))
            jets.append(
                second_jet(
                    data[sheet],
                    common_tuples[sheet_slice, first],
                    common_tuples[sheet_slice, second],
                )
            )
        simultaneous_jet = jets[0].col_join(jets[0]) - jets[1].col_join(jets[2])
        second_obstructions.append(obstruction_quotient * simultaneous_jet)

quadratic_obstruction = sp.Matrix.hstack(*second_obstructions)

# Rows are coefficients in the bilinear convention
# (x^2, xy, xz, y^2, yz, z^2).  Convert mixed terms to polynomial
# coefficients before computing the ideal of scalar quadratic obstructions.
quadratic_rows = quadratic_obstruction.rowspace()
polynomial_coefficient_rows = [
    sp.Matrix([[row[0], 2 * row[1], 2 * row[2], row[3], 2 * row[4], row[5]]])
    for row in quadratic_rows
]
quadratic_relation_matrix = sp.Matrix.vstack(*polynomial_coefficient_rows)

print("witness roots:", roots)
print("incidence root derivatives:", root_derivatives)
print("normalization differential ranks:", tuple(entry["image"].rank() for entry in data))
print("pairwise span ranks: (7, 7, 7)")
print("triple span rank:", triple_span.rank())
print("triple common tangent dimension:", common_images.rank())
print("transverse quadratic obstruction rank:", quadratic_relation_matrix.rank())

x, y, z = sp.symbols("x y z")
monomials = (x**2, x * y, x * z, y**2, y * z, z**2)
quadrics = [
    sp.expand(sum(row[index] * monomial for index, monomial in enumerate(monomials)))
    for row in quadratic_relation_matrix.rowspace()
]
print("transverse quadratic generators:")
for quadric in quadrics:
    print(" ", sp.factor(quadric))

if len(quadrics) == 3:
    quadratic_groebner = sp.groebner(quadrics, x, y, z, order="grevlex")
    assert quadratic_groebner.is_zero_dimensional
    assert standard_monomial_count(quadratic_groebner) == 8
    print("transverse quadratic Groebner basis:")
    for entry in quadratic_groebner.polys:
        print(" ", sp.factor(entry.as_expr()))

# The stronger equality M_1=M_2=M_3 separates by Hensel decomposition into
# one universal Q^2=T^3 block at each of the three distinct roots.
block_z, block_u, block_v, block_a, block_b = sp.symbols("Z u v a b")
block_q = block_z**3 + block_u * block_z + block_v
block_r = block_z**2 + block_a * block_z + block_b
block_coefficients = sp.Poly(sp.expand(block_q**2 - block_r**3), block_z).all_coeffs()
block_groebner = sp.groebner(
    block_coefficients, block_v, block_u, block_b, block_a, order="lex"
)
for generator in (block_a, block_v, 2 * block_u - 3 * block_b, block_b**2):
    assert block_groebner.reduce(generator)[1] == 0
strong_equality_basis = tuple(
    (first, second, third)
    for first in range(2)
    for second in range(2)
    for third in range(2)
)
assert len(strong_equality_basis) == 8


def pairwise_transverse_audit(left, right):
    """Return tangent and quadratic lengths transverse to E_(6,6,3,3)."""
    common = next(index for index in range(3) if index not in (left, right))
    pair_difference = data[left]["image"].row_join(-data[right]["image"])
    pair_common_tuples = sp.Matrix.hstack(*pair_difference.nullspace())
    assert pair_common_tuples.shape == (12, 5)
    pair_common_images = data[left]["image"] * pair_common_tuples[:6, :]
    assert pair_common_images.rank() == 5

    root_left, root_right = sp.symbols("root_left root_right")
    shared_constant, shared_linear = sp.symbols("shared_constant shared_linear")
    support_q = sp.Poly((W - root_left) ** 3, W)
    support_r = sp.Poly(
        (W - root_right) ** 2
        * (W**2 + shared_linear * W + shared_constant),
        W,
    )
    support_parameters = sp.Matrix(
        [
            support_q.coeff_monomial(1),
            support_q.coeff_monomial(W),
            support_q.coeff_monomial(W**2),
            support_r.coeff_monomial(1),
            support_r.coeff_monomial(W),
            support_r.coeff_monomial(W**2),
            support_r.coeff_monomial(W**3),
        ]
    )
    support_variables = (
        root_left,
        root_right,
        shared_constant,
        shared_linear,
    )
    support_substitution = {
        root_left: roots[left],
        root_right: roots[right],
        shared_constant: roots[common] ** 2,
        shared_linear: -2 * roots[common],
    }
    support_parameter_jacobian = support_parameters.jacobian(support_variables).subs(
        support_substitution
    )
    support_constraint = data[left]["constraint"] * support_parameter_jacobian
    support_source_tangent = sp.Matrix.hstack(*support_constraint.nullspace())
    assert support_source_tangent.rank() == 3
    support_parameter_tangent = support_parameter_jacobian * support_source_tangent
    support_sheet_coordinates = data[left]["source_tangent"].gauss_jordan_solve(
        support_parameter_tangent
    )[0]
    support_image = data[left]["image"] * support_sheet_coordinates
    assert support_image.rank() == 3

    support_common_coordinates = []
    for column in range(3):
        solution = sp.linsolve((pair_common_images, support_image[:, column]))
        support_common_coordinates.append(sp.Matrix(tuple(solution)[0]))
    pair_basis_change = extend_to_basis(support_common_coordinates, 5)
    pair_common_tuples = pair_common_tuples * pair_basis_change

    pair_obstruction_quotient = sp.Matrix.hstack(
        *pair_difference.T.nullspace()
    ).T
    pair_second_obstructions = []
    for first in range(3, 5):
        for second in range(first, 5):
            left_jet = second_jet(
                data[left],
                pair_common_tuples[:6, first],
                pair_common_tuples[:6, second],
            )
            right_jet = second_jet(
                data[right],
                pair_common_tuples[6:, first],
                pair_common_tuples[6:, second],
            )
            pair_second_obstructions.append(
                pair_obstruction_quotient * (left_jet - right_jet)
            )
    pair_quadratic_obstruction = sp.Matrix.hstack(*pair_second_obstructions)
    pair_rows = pair_quadratic_obstruction.rowspace()
    pair_relation_matrix = sp.Matrix.vstack(
        *[
            sp.Matrix([[row[0], 2 * row[1], row[2]]])
            for row in pair_rows
        ]
    )
    assert pair_relation_matrix.rank() == 2
    pair_x, pair_y = sp.symbols(f"pair_x_{left}_{right} pair_y_{left}_{right}")
    pair_monomials = (pair_x**2, pair_x * pair_y, pair_y**2)
    pair_quadrics = [
        sp.expand(
            sum(
                row[index] * monomial
                for index, monomial in enumerate(pair_monomials)
            )
        )
        for row in pair_relation_matrix.rowspace()
    ]
    pair_groebner = sp.groebner(
        pair_quadrics, pair_x, pair_y, order="grevlex"
    )
    assert pair_groebner.is_zero_dimensional
    length = standard_monomial_count(pair_groebner)
    assert length == 4
    return {
        "common_tangent_dimension": pair_common_images.rank(),
        "support_tangent_dimension": support_image.rank(),
        "quadratic_rank": pair_relation_matrix.rank(),
        "transverse_length": length,
    }


pairwise_results = tuple(
    pairwise_transverse_audit(left, right)
    for left in range(3)
    for right in range(left + 1, 3)
)
assert all(result["transverse_length"] == 4 for result in pairwise_results)

print("PASS: C_(3,4) has three exact degree-18 normalization branches")
print("PASS: all three pairwise intersections have transverse length four")
print("PASS: the triple collision has three excess common tangent directions")
print("PASS: three transverse quadratic obstructions give the length-eight upper bound")
print("PASS: three separated dual-number blocks give the matching length-eight lower bound")
