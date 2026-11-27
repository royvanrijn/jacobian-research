"""Exact checks for the smooth-member nonbinary axial exclusion."""

from itertools import product

from sympy import Matrix, expand, factor, hessian, symbols, trace


def check_lower_ties():
    # Only smaller z-degrees are free in the degree-at-most-five tail.
    expected = {
        (3, 4, 6): ((28, {(0, 1, 4)}), (25, {(1, 1, 3)}),
                    (22, {(2, 1, 2), (0, 4, 1)})),
        (2, 3, 4): ((19, {(0, 1, 4)}), (17, {(1, 1, 3), (0, 3, 2)}),
                    (15, {(2, 1, 2), (1, 3, 1), (0, 5, 0)})),
        (3, 4, 5): ((24, {(0, 1, 4)}), (22, {(1, 1, 3), (0, 3, 2)}),
                    (20, {(2, 1, 2), (1, 3, 1), (0, 5, 0)})),
    }
    for weights, rows in expected.items():
        for index, (degree, ties) in enumerate(rows):
            z_degree = 4 - index
            candidates = {
                (i, j, l) for i, j, l in product(range(6), repeat=3)
                if i + j + l <= 5 and l < z_degree
            }
            candidates.add((index, 1, z_degree))
            if index == 2 and weights != (3, 4, 6):
                candidates.remove((0, 4, 1))
            weighted = {mon: sum(a*b for a, b in zip(weights, mon)) for mon in candidates}
            assert max(weighted.values()) == degree
            assert {mon for mon in weighted if weighted[mon] == degree} == ties


def check_restricted_determinant():
    x, y, z, nu, beta = symbols("x y z nu beta")
    conic = x*z + nu*y**2
    triple = Matrix([
        conic**3 + y*z**4,
        conic**2*x**2 + x*y*z**3 + beta*y**3*z**2,
        conic*x**4 + x**2*y*z**2 + 2*beta*x*y**3*z + beta**2*y**5,
    ])
    determinant = triple.jacobian((x, y, z)).subs(z, 0).det()
    expected = 6*beta*nu**4*(5*beta - 4*nu)*x**2*y**12
    assert factor(determinant - expected) == 0
    assert factor(determinant.subs(beta, -nu/3)) == 34*nu**6*x**2*y**12/3


def check_residue_by_inverse():
    # Only these perturbations can contribute to the specified residue sector.
    A, B, z, nu, beta = symbols("A B z nu beta")
    variables = (A, B)
    first = Matrix([B**3*z**6/A**3, B**4*z**6/A**4 + beta*A**3/z**10])
    second = Matrix([3*nu*B**2/z**4, 2*nu*B**3/(A*z**4)])

    # Invert the two level equations by substitution in each order.
    d_first = first.jacobian(variables)
    d_second = second.jacobian(variables)
    inverse_one = -first
    inverse_two = -second - d_first*inverse_one
    quadratic = Matrix([
        (inverse_one.T*hessian(term, variables)*inverse_one)[0]/2
        for term in first
    ])
    inverse_three = -d_first*inverse_two - d_second*inverse_one - quadratic
    jac_one, jac_two, jac_three = (
        term.jacobian(variables) for term in (inverse_one, inverse_two, inverse_three)
    )

    det_one = trace(jac_one)
    det_two = trace(jac_two) + jac_one.det()
    det_three = trace(jac_three) + (jac_one + jac_two).det() - jac_one.det() - jac_two.det()

    density = -1/(A*z**3)
    r, s, t = inverse_one[0], inverse_two[0], inverse_three[0]
    density_one = -density*r/A
    density_two = density*(-s/A + r**2/A**2)
    density_three = density*(-t/A + 2*r*s/A**2 - r**3/A**3)
    coefficient = expand(
        density_three + density_two*det_one + density_one*det_two + density*det_three
    ).coeff(z, -1)
    assert factor(coefficient - 6*B**5*(nu + 3*beta)/A**6) == 0


if __name__ == "__main__":
    check_lower_ties()
    check_restricted_determinant()
    check_residue_by_inverse()
    print("Checked axial ties, the restricted determinant, and the third-order residue.")
