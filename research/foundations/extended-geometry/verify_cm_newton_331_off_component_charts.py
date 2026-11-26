"""Exact 331-adic chart exclusions at the two off-L3 conductor-18 points."""

from itertools import combinations

from sage.all import GF, QQ, PolynomialRing, binomial, matrix

from verify_cm_newton_full_width_projective_boundary import (
    substitute_truncated,
    truncate,
    weighted_degree,
)
from verify_cm_newton_full_width_reconstruction_certificates import reduced_residual_rows


FIELD = GF(331)
NORMAL_WEIGHTS = (2, 1, 2, 1, 2)


def unit_vector(dimension, index):
    return tuple(int(i == index) for i in range(dimension))


def implicit_jets(rows, weights, order, names):
    """Solve the first coordinate using the first row, rescaling the free ones."""
    parameter = rows[0].base_ring()
    free = PolynomialRing(parameter, names=names)
    derivative = rows[0][unit_vector(len(weights) + 1, 0)]
    assert derivative
    values = [free.zero(), *(derivative**2*x for x in free.gens())]
    for degree in range(1, order + 1):
        error = substitute_truncated(rows[0].dict(), values, weights, degree)
        terms = {}
        for exponent, coefficient in error.dict().items():
            assert weighted_degree(exponent, weights) == degree
            quotient, remainder = coefficient.quo_rem(derivative)
            assert not remainder
            terms[exponent] = quotient
        values[0] -= free(terms)
    jets = [substitute_truncated(row.dict(), values, weights, order) for row in rows]
    assert not jets[0]
    return derivative, jets[1:]


def first_jets(rows, weights, order):
    """Use shift 7 to solve for r in the chart ell=1."""
    parameter = PolynomialRing(QQ, "v")
    v = parameter.gen()
    chart = PolynomialRing(parameter, names=("r", "h", "d", "tau"))
    r, h, d, tau = chart.gens()
    values = (r, chart.one(), h, chart(v), (8*v - 1 + d)/17, tau)
    polynomials = []
    for row in rows:
        degree = max(weighted_degree(e, NORMAL_WEIGHTS) for e in row.dict())
        terms = {}
        for exponent, coefficient in row.dict().items():
            FIELD(coefficient)  # Check that each original row is integral at 331.
            power, parity = divmod(degree - weighted_degree(exponent, NORMAL_WEIGHTS), 2)
            assert not parity
            terms[tuple(exponent) + (power,)] = coefficient
        polynomial = substitute_truncated(terms, values, (1, *weights), order)
        assert not polynomial.constant_coefficient()
        polynomials.append(polynomial)
    return implicit_jets(polynomials, weights, order, ("x", "y", "z"))


def centered_terms(rows, center, order=None):
    """Expand v=center+epsilon, with an optional total-degree cap."""
    terms = []
    for row in rows:
        result = {}
        for exponent, coefficient in row.dict().items():
            for power, value in (coefficient(coefficient.parent().gen() + center)).dict().items():
                if order is None or sum(exponent) + power <= order:
                    assert sum(exponent) + power >= 2
                    result[tuple(exponent) + (power,)] = value
        terms.append(result)
    return terms


def epsilon_chart(terms):
    """Set (x,y,z,epsilon)=s(X,a,b,1) and divide by s^2."""
    parameter = PolynomialRing(QQ, "s")
    s = parameter.gen()
    chart = PolynomialRing(parameter, names=("X", "a", "b"))
    result = []
    for row in terms:
        coefficients = {}
        for exponent, value in row.items():
            e = tuple(exponent[:3])
            coefficients[e] = coefficients.get(e, parameter.zero()) + value*s**(sum(exponent)-2)
        polynomial = chart(coefficients)
        assert not polynomial.constant_coefficient()
        result.append(polynomial)
    return result


def verify_epsilon_support(rows):
    """Compute exact s-jets through order 3 with a and b still free."""
    free = PolynomialRing(QQ, names=("a", "b"))
    a, b = free.gens()
    first = rows[0]
    dx, da, db = [first[unit_vector(3, i)][0] for i in range(3)]
    assert FIELD(dx)
    center = -(da*a + db*b)/dx
    expansions = []
    for row in rows:
        taylor = {
            (power, derivative): free(sum(
                coefficient[power]*binomial(i, derivative)*center**(i-derivative)*a**j*b**k
                for (i, j, k), coefficient in row.dict().items() if i >= derivative
            ))
            for power in range(4) for derivative in range(4-power)
        }
        assert not taylor[0, 0]
        expansions.append(taylor)
    pivot = expansions[0]
    D = pivot[0, 1]
    assert FIELD(D[0, 0])

    # X=center+u1*s/D+u2*s^2/D^3+u3*s^3/D^5 modulo s^4.
    def numerators(taylor, u1, u2, u3):
        return (
            taylor[1, 0]*D + taylor[0, 1]*u1,
            taylor[2, 0]*D**3 + taylor[1, 1]*u1*D**2
            + taylor[0, 2]*u1**2*D + taylor[0, 1]*u2,
            taylor[3, 0]*D**5 + taylor[2, 1]*u1*D**4
            + taylor[1, 2]*u1**2*D**3 + taylor[0, 3]*u1**3*D**2
            + taylor[1, 1]*u2*D**2 + 2*taylor[0, 2]*u1*u2*D + taylor[0, 1]*u3,
        )

    u1 = -pivot[1, 0]
    u2, remainder = (-numerators(pivot, u1, 0, 0)[1]).quo_rem(D)
    assert not remainder
    u3, remainder = (-numerators(pivot, u1, u2, 0)[2]).quo_rem(D)
    assert not remainder
    assert not any(numerators(pivot, u1, u2, u3))

    monomials = ((5, 0, 3), (3, 1, 2), (1, 2, 1), (0, 3, 1))
    balanced = PolynomialRing(QQ, names=("A", "U"))
    A, U = balanced.gens()

    def b_coefficient(polynomial, power):
        return sum(c*A**i for (i, j), c in polynomial.dict().items() if j == power)

    coefficient_matrix = []
    balanced_numerators = []
    for taylor in expansions[1:]:
        numerators_at_row = numerators(taylor, u1, u2, u3)
        for s_degree, numerator in enumerate(numerators_at_row, 1):
            assert all(
                any(i >= m and j >= n and s_degree >= q for m, n, q in monomials)
                for i, j in numerator.dict()
            )
        coefficient_matrix.append([
            FIELD(numerators_at_row[k-1][i,j]/D[0,0]**(2*k-1)) for i,j,k in monomials
        ])
        n1, n2, n3 = numerators_at_row
        # Coefficient of s^3 after b=s*A^2*U, with denominator D(A,0)^5.
        D0 = b_coefficient(D, 0)
        balanced_numerators.append(
            D0**4*A**4*U**2*b_coefficient(n1, 2)
            + D0**2*A**2*U*b_coefficient(n2, 1)
            + b_coefficient(n3, 0)
        )

    coefficients = matrix(FIELD, coefficient_matrix)
    assert all(coefficients[0, i] for i in range(4))
    assert coefficients.matrix_from_rows_and_columns([0, 1], [2, 3]).det()
    first, second = balanced_numerators[:2]
    assert first[5, 0] and first[5, 0]*second == second[5, 0]*first


def verify_epsilon_balance(rows, expected):
    """Compute the balanced leading polynomials with exact normal-weight jets."""
    weights = (1, 2)
    delta, residuals = implicit_jets(rows, weights, 6, ("a", "b"))
    assert FIELD(delta[0])
    assert all(not truncate(row, weights, 4) for row in residuals)
    ninth, tenth = residuals[:2]
    combination = ninth[5, 0]*tenth - tenth[5, 0]*ninth
    assert not truncate(combination, weights, 5)
    R = PolynomialRing(FIELD, "U")
    U = R.gen()

    def balanced_leading(row, normal_weight, s_degree):
        return sum(
            FIELD(coefficient[s_degree-j])*U**j
            for (i, j), coefficient in row.dict().items() if i+2*j == normal_weight
        )

    Q = balanced_leading(ninth, 5, 3)
    H = balanced_leading(combination, 6, 7)
    assert (Q.degree(), H.degree(), Q.resultant(H)) == (2, 3, FIELD(expected))
    print(f"  The epsilon-dominant chart is excluded (resultant {expected}).")


def verify_dominant_d(terms, kappa, branch_slope, exception):
    """Use the d-dominant chart where the shift 8 pivot is a unit."""
    parameter = PolynomialRing(QQ, "E")
    E = parameter.gen()
    chart = PolynomialRing(parameter, names=("X", "Z", "s"))
    X, Z, s = chart.gens()
    chart_rows = [
        sum(chart(c)*(kappa+X)**i*Z**k*E**q*s**(i+j+k+q-2)
            for (i,j,k,q), c in row.items())
        for row in terms
    ]
    assert all(not row.constant_coefficient() for row in chart_rows)
    assert all(not row(-branch_slope*Z, Z, 0) for row in chart_rows)
    weights = (1, 1)
    delta, residuals = implicit_jets(chart_rows, weights, 4, ("a", "b"))
    assert all(j for row in residuals for i,j in row.dict())
    free = residuals[0].parent()
    divided = [free({(i,j-1):c for (i,j),c in row.dict().items()}) for row in residuals]
    assert all(not truncate(row, weights, 1) for row in divided)

    initials = [truncate(row, weights, 2) for row in divided]
    for f, g in combinations(divided, 2):
        for exponent in ((2, 0), (1, 1), (0, 2)):
            combination = f*g[exponent] - g*f[exponent]
            assert not truncate(combination, weights, 2)
            initials.append(truncate(combination, weights, 3))

    projective = PolynomialRing(FIELD, names=("T", "inverse", "E"), order="lex")
    T, inverse, E = projective.gens()

    def reduce(coefficient):
        return sum(FIELD(c)*E**i for i,c in coefficient.dict().items())

    determinant = reduce(delta)
    assert determinant.degree() == 1 and not determinant(0, 0, FIELD(exception))
    for coordinate in (0, 1):
        equations = {
            sum(reduce(c)*T**exponent[coordinate] for exponent,c in form.dict().items())
            for form in initials
        }
        equations.add(inverse*determinant - 1)
        assert projective.ideal(list(equations)).groebner_basis() == [1]
    print(f"  The d-dominant chart is excluded outside E={exception}.")


def verify_tangent(rows, kappa):
    exact = PolynomialRing(QQ, names=("x", "y", "z", "epsilon"))
    jets = [exact(row) for row in rows]
    quadratic = [truncate(row, (1, 1, 1, 1), 2) for row in jets]
    monomials = sorted({e for row in quadratic for e in row.dict()})
    M = matrix(QQ, [[row[e] for e in monomials] for row in quadratic])
    mod = M.change_ring(FIELD)
    pivot_rows = mod.transpose().pivots()
    pivot_columns = mod.matrix_from_rows(pivot_rows).pivots()
    assert len(pivot_rows) == M.rank()
    square = M.matrix_from_rows_and_columns(pivot_rows, pivot_columns)
    determinant = square.det()
    assert FIELD(determinant)
    initial = quadratic[:]
    basis = M.matrix_from_columns(pivot_columns)*square.adjugate()
    for j, row in enumerate(jets):
        combination = determinant*row - sum(
            basis[j,i]*jets[k] for i,k in enumerate(pivot_rows)
        )
        assert not truncate(combination, (1, 1, 1, 1), 2)
        initial.append(combination)
    residue = PolynomialRing(FIELD, names=exact.variable_names())
    x, y, z, _ = residue.gens()
    ideal = residue.ideal([residue({e:FIELD(c) for e,c in row.dict().items()}) for row in initial])
    assert ideal.radical() == residue.ideal(z, x-FIELD(kappa)*y)


def verify_off_component_charts(rows):
    J, ordinary = first_jets(rows, (1, 1, 1), 6)
    _, weighted = first_jets(rows, (1, 1, 2), 6)
    for center, kappa, branch_slope, resultant, exception in (
        (QQ(2)/17, QQ(17)/18, QQ(4454)/67, 163, 275),
        (QQ(5)/56, QQ(7)/9, QQ(61246)/1139, 104, 170),
    ):
        print(f"At v/ell={center} modulo 331:")
        assert FIELD(J(center))
        verify_tangent(centered_terms(ordinary, center, 3), kappa)

        verify_epsilon_support(epsilon_chart(centered_terms(ordinary, center, 5)))
        verify_epsilon_balance(epsilon_chart(centered_terms(weighted, center)), resultant)
        verify_dominant_d(centered_terms(ordinary, center, 6), kappa, branch_slope, exception)


if __name__ == "__main__":
    residuals = reduced_residual_rows(5, 16, even_linear=True, primitive=False)
    verify_off_component_charts(residuals)
