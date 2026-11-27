"""Exact determinant replay for the remaining normalized ternary quintic form."""

from itertools import combinations_with_replacement, permutations
from subprocess import run

from sympy import Matrix, Poly, Rational, groebner, symbols

import verify_quintic_ternary_pencil_first_equations as first
from verify_quintic_ternary_pencil_unequal_parameters import homogeneous_terms


x, y, z = first.x, first.y, first.z
variables = (x, y, z)
p, q = first.p, first.q
v, w, k, b, c, d, e = symbols("v w k b c d e")


def polynomial_column(expression, degree):
    poly = Poly(expression, *variables)
    return Matrix([poly.coeff_monomial(term) for term in homogeneous_terms(degree)])


def determinant_part(jacobians, degree):
    pieces = []
    for triple in combinations_with_replacement(sorted(jacobians), 3):
        if sum(triple) != degree+3:
            continue
        for orders in sorted(set(permutations(triple))):
            columns = [jacobians[order][:, j] for j, order in enumerate(orders)]
            pieces.append(Matrix.hstack(*columns).det())
    return sum(pieces)


def solve_linear(degree, rhs, prefix):
    """Parameterize all solutions of L(vector) = -rhs when compatible."""
    terms = homogeneous_terms(degree)
    columns = []
    for coefficient in first.normal:
        for term in terms:
            f = coefficient*term
            columns.append(polynomial_column(5*z**3*(-z*f.diff(x)+2*x*f.diff(y)), degree+7))
    operator = Matrix.hstack(*columns)
    values = polynomial_column(rhs, degree+7)
    pivots = operator.rref()[1]
    rows = operator[:, pivots].T.rref()[1]
    selected = -operator.extract(rows, pivots).inv()*values.extract(rows, [0])

    kernel = operator.nullspace()
    free = symbols(f"{prefix}:{len(kernel)}")
    coefficients = Matrix.hstack(*kernel)*Matrix(free)
    for pivot, value in zip(pivots, selected, strict=True):
        coefficients[pivot] += value
    vector = Matrix([
        sum(coefficients[i*len(terms)+j]*term for j, term in enumerate(terms))
        for i in range(3)
    ])
    compatibility = [row.dot(values).expand() for row in operator.T.nullspace()]
    return vector.applyfunc(lambda f: f.expand()), free, compatibility


def has_witness(equations, witness):
    return witness.expand() in equations


def determinant_setup(k_value, b_value, transverse=1):
    U = k_value*q + transverse*x*z
    V = v*x**2+w*y*z
    quartic = Matrix([
        2*p*U+b_value*p*q+c*p**2,
        q*U+p*V+d*p**2,
        2*q*V+e*p**2,
    ])
    jacobians = {5: first.J5, 4: quartic.jacobian(variables)}
    cubic, free, compatibility = solve_linear(3, determinant_part(jacobians, 10), "f")
    assert not any(compatibility)
    jacobians[3] = cubic.jacobian(variables)
    return jacobians, free


def specialize(jacobians, substitution):
    return {degree: matrix.subs(substitution, simultaneous=True).applyfunc(lambda f: f.expand())
            for degree, matrix in jacobians.items()}


def check_transverse():
    jacobians, _ = determinant_setup(k, b)
    _, _, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert has_witness(conditions, 2*k*(3*b+10*k)/5)
    jacobians, _ = determinant_setup(k, -10*k/3)
    _, _, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert has_witness(conditions, 4*k)

    jacobians, f = determinant_setup(0, b)
    _, _, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert has_witness(conditions, -6*b*(5*f[4]+1)/5)
    assert has_witness(conditions, 9*b*f[5]-3*b*v+3*f[1]-3)
    jacobians = specialize(jacobians, {f[4]: Rational(-1, 5), f[1]: -3*b*f[5]+b*v+1})
    quadratic, _, _ = solve_linear(2, determinant_part(jacobians, 9), "g")
    jacobians[2] = quadratic.jacobian(variables)
    _, _, conditions = solve_linear(1, determinant_part(jacobians, 8), "h")
    assert has_witness(conditions, -12*b/25)

    jacobians, f = determinant_setup(0, 0)
    _, _, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert has_witness(conditions, 3*f[1]-3)
    jacobians = specialize(jacobians, {f[1]: 1})
    quadratic, g, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert not any(conditions)
    jacobians[2] = quadratic.jacobian(variables)

    _, _, conditions = solve_linear(1, determinant_part(jacobians, 8), "h")
    assert has_witness(conditions, -20*f[4]**2)
    assert has_witness(conditions, 5*(f[8]+4*v**2-4*v*w)**2)
    jacobians = specialize(jacobians, {f[4]: 0, f[8]: 4*v*(w-v)})
    _, _, conditions = solve_linear(1, determinant_part(jacobians, 8), "h")
    witnesses = [c, d, e, g[1]-2*(f[5]-v)]
    parameters = sorted(set().union(*(eq.free_symbols for eq in conditions)), key=str)
    basis = groebner(conditions, *parameters, domain="QQ")
    assert all(basis.reduce(witness)[1] == 0 for witness in witnesses)

    jacobians = specialize(jacobians, {c: 0, d: 0, e: 0, g[1]: 2*(f[5]-v)})
    linear, _, conditions = solve_linear(1, determinant_part(jacobians, 8), "h")
    assert not any(conditions)
    jacobians[1] = linear.jacobian(variables)

    verify_no_keller_completion(jacobians, (7, 6, 5))


def check_equal_parameters():
    shift = w-v
    jacobians, f = determinant_setup(k, b, transverse=0)
    _, _, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert has_witness(conditions, -2*(3*b+10*k)*(f[4]-2*k*shift))
    assert has_witness(conditions, 3*f[5]*(3*b+10*k))
    common = {f[4]: 2*k*shift, f[5]: 0}
    jacobians = specialize(jacobians, common)
    quadratic, _, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert not any(conditions)
    jacobians[2] = quadratic.jacobian(variables)
    _, _, conditions = solve_linear(1, determinant_part(jacobians, 8), "h")
    assert has_witness(conditions, 5*(f[8]-4*v*shift)**2)
    assert has_witness(conditions, -f[9]*k*(3*b+10*k))

    exceptional, f = determinant_setup(k, -10*k/3, transverse=0)
    quadratic, _, conditions = solve_linear(2, determinant_part(exceptional, 9), "g")
    assert not any(conditions)
    exceptional[2] = quadratic.jacobian(variables)
    _, _, conditions = solve_linear(1, determinant_part(exceptional, 8), "h")
    m = f[4]-2*k*shift
    s = -48*c*k-15*f[1]+220*k*v
    assert has_witness(conditions, 2*m*s/15)
    assert has_witness(conditions, 2*m*(s-100*m)/5)
    assert has_witness(conditions, 5*(f[8]-4*v*shift)**2)
    exceptional = specialize(exceptional, {f[4]: 2*k*shift, f[8]: 4*v*shift})
    linear, _, _ = solve_linear(1, determinant_part(exceptional, 8), "h")
    exceptional[1] = linear.jacobian(variables)
    degree_seven = Poly(determinant_part(exceptional, 7), *variables)
    assert degree_seven.coeff_monomial(x**7).expand() == 4*f[5]**2*k
    assert degree_seven.coeff_monomial(x**3*z**4).subs(f[5], 0).expand() == 5*f[9]**2*k

    # Both cases with k nonzero now have the same parameterized coefficients.
    jacobians = specialize(jacobians, {f[8]: 4*v*shift, f[9]: 0})
    linear, _, conditions = solve_linear(1, determinant_part(jacobians, 8), "h")
    assert not any(conditions)
    assert linear.jacobian(variables).det().expand() == 0


def check_zero_k():
    jacobians, f = determinant_setup(0, b, transverse=0)
    jacobians = specialize(jacobians, {f[4]: 0})
    _, _, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert has_witness(conditions, 9*b*f[5])

    # For b=0, the next three equations force f5 to vanish as well.
    branch = specialize(jacobians, {b: 0})
    quadratic, g, conditions = solve_linear(2, determinant_part(branch, 9), "g")
    assert not any(conditions)
    branch[2] = quadratic.jacobian(variables)
    _, _, conditions = solve_linear(1, determinant_part(branch, 8), "h")
    assert has_witness(conditions, f[1]*f[5])
    assert has_witness(conditions, 5*(f[8]+4*v**2-4*v*w)**2)
    branch = specialize(branch, {f[8]: 4*v*(w-v)})
    linear, _, _ = solve_linear(1, determinant_part(branch, 8), "h")
    branch[1] = linear.jacobian(variables)
    degree_seven = Poly(determinant_part(branch, 7), *variables)
    degree_six = Poly(determinant_part(branch, 6), *variables)
    assert degree_seven.coeff_monomial(x**6*z).subs(f[1], 0).expand() == -f[5]*g[1]
    assert degree_six.coeff_monomial(x**6).subs({f[1]: 0, g[1]: 0}).expand() == 3*f[5]**3

    jacobians = specialize(jacobians, {f[5]: 0})
    quadratic, _, conditions = solve_linear(2, determinant_part(jacobians, 9), "g")
    assert not any(conditions)
    jacobians[2] = quadratic.jacobian(variables)
    _, _, conditions = solve_linear(1, determinant_part(jacobians, 8), "h")
    assert has_witness(conditions, 5*(f[8]+4*v**2-4*v*w)**2)
    jacobians = specialize(jacobians, {f[8]: 4*v*(w-v)})
    linear, _, conditions = solve_linear(1, determinant_part(jacobians, 8), "h")
    assert not any(conditions)
    jacobians[1] = linear.jacobian(variables)
    verify_no_keller_completion(jacobians, (7, 6))


def verify_no_keller_completion(jacobians, degrees):
    equations = []
    for degree in degrees:
        equations.extend(Poly(determinant_part(jacobians, degree), *variables).coeffs())
    invdet = symbols("invdet")
    equations.append(invdet*jacobians[1].det()-1)
    parameters = sorted(set().union(*(eq.free_symbols for eq in equations)), key=str)
    equations = [Poly(eq, *parameters).clear_denoms()[1].as_expr() for eq in equations]
    generators = ",\n".join(str(eq).replace("**", "^") for eq in equations)
    script = f"""ring R=0,({','.join(map(str, parameters))}),dp;
ideal I={generators};
ideal B=std(I);
reduce(1,B);
quit;
"""
    result = run(["Singular", "-q"], input=script, text=True, capture_output=True, check=True)
    assert result.stdout.strip() == "0" and not result.stderr, (result.stdout, result.stderr)


if __name__ == "__main__":
    check_transverse()
    check_equal_parameters()
    check_zero_k()
    print("Verified the exact exclusions in all equal-parameter ternary quintic sectors.")
