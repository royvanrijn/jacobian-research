"""Intrinsic Prony and Ritt matrices in optimal Gaussian moment coordinates.

The constructions in this module never start with equations in seed
coefficients.  They first recover the normalized seed from the optimal
factorial-normalized moments by series reversion, then build:

* the log-Prony Toeplitz--Hessenberg matrix whose leading determinants test
  whether the reversed seed has a polynomial ``m``-th root;
* saturated Christoffel--Hankel/Sylvester rank conditions for every mixed
  fixed-weight partition; and
* the Ritt--Krylov coefficient matrix whose maximal minors test whether a
  vertical translate belongs to ``k[B]``.

All formulas are over characteristic zero and are intended on the
exact-degree open where the leading seed coefficient is invertible.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import factorial

import sympy as sp


@dataclass(frozen=True)
class OptimalMomentSeed:
    """A normalized degree-``N`` seed recovered from its optimal moments."""

    degree: int
    variable: sp.Symbol
    series_variable: sp.Symbol
    moments: tuple[sp.Symbol, ...]
    generating_series: sp.Expr
    inverse_series: sp.Expr
    bridge_polynomial: sp.Expr
    seed: sp.Expr
    coefficients: tuple[sp.Expr, ...]

    def coefficient(self, degree: int) -> sp.Expr:
        """Return the coefficient of ``variable**degree`` in the seed."""
        return self.coefficients[degree]

    @property
    def leading_coefficient(self) -> sp.Expr:
        return self.coefficients[self.degree]


@dataclass(frozen=True)
class PowerPronyData:
    """Formal power-root data and its log-Prony Hessenberg matrix."""

    multiplicity: int
    root_degree: int
    reversed_coefficients: tuple[sp.Expr, ...]
    logarithmic_weights: tuple[sp.Expr, ...]
    hessenberg_matrix: sp.Matrix
    root_coefficients: tuple[sp.Expr, ...]
    equation_indices: tuple[int, ...]

    def leading_minor(self, size: int) -> sp.Expr:
        """Return the size-``size`` leading determinant."""
        return self.hessenberg_matrix[:size, :size].det(method="domain-ge")


@dataclass(frozen=True)
class RittKrylovData:
    """Approximate inner polynomial and coefficient Krylov matrix."""

    outer_degree: int
    inner_degree: int
    polynomial: sp.Expr
    approximate_inner: sp.Expr
    matrix: sp.Matrix
    pivot_rows: tuple[int, ...]
    selected_rows: tuple[int, ...]
    selected_minors: tuple[sp.Expr, ...]


@dataclass(frozen=True)
class FixedWeightPronyData:
    """Marked fixed-weight Prony incidence for weights ``(2,3,3)``.

    The distinguished node has weight two.  The four relations say that,
    after subtracting that node and dividing by three, the remaining power
    sums come from two unit-weight nodes.  The Fitting matrix presents the
    finite incidence algebra over the power-sum base; its maximal minors have
    the correct support but need not give the reduced scheme at collisions.
    """

    distinguished_node: sp.Symbol
    power_sums: tuple[sp.Symbol, ...]
    relations: tuple[sp.Expr, ...]
    fitting_matrix: sp.Matrix


@dataclass(frozen=True)
class WeightSubresultantCondition:
    """One fixed-weight gcd condition on a Sylvester matrix."""

    weight: int
    required_gcd_degree: int
    weight_polynomial: sp.Expr
    sylvester_matrix: sp.Matrix
    rank_bound: int

    def minors(self) -> Iterator[sp.Expr]:
        """Yield the minors cutting out the required Sylvester rank drop."""

        yield from rank_drop_minors(self.sylvester_matrix, self.rank_bound)


@dataclass(frozen=True)
class FixedWeightHankelData:
    """Unmarked fixed-weight Prony equations on the distinct-node closure."""

    multiplicities: tuple[int, ...]
    power_sums: tuple[sp.Expr, ...]
    variable: sp.Symbol
    hankel_matrix: sp.Matrix
    shifted_hankel_matrix: sp.Matrix
    hankel_determinant: sp.Expr
    support_polynomial: sp.Expr
    christoffel_numerator: sp.Expr
    weight_conditions: tuple[WeightSubresultantCondition, ...]
    extension_equations: tuple[sp.Expr, ...]


def _coefficient(expression: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expression), variable).nth(degree)


def optimal_moment_seed(
    degree: int,
    *,
    variable: sp.Symbol | None = None,
    series_variable: sp.Symbol | None = None,
    moment_prefix: str = "mu",
) -> OptimalMomentSeed:
    """Recover the normalized seed from ``mu_3,...,mu_(N-1)``.

    If ``g(u)=u+sum(mu_j*u**j)`` and ``eta=g^{-1}``, then the bridge
    polynomial is ``h=z/eta``.  Its coefficients through degree ``N-2`` are
    recovered by the reciprocal Toeplitz recurrence.  The two endpoint
    equations ``H(1)=0`` and ``H'(1)=-1`` complete degrees ``N-1,N``.
    """

    if degree < 4:
        raise ValueError("optimal moment coordinates require degree at least four")
    w = variable or sp.Symbol("W")
    z = series_variable or sp.Symbol("_moment_z")
    u = sp.Symbol("_moment_u")
    moments = tuple(sp.Symbol(f"{moment_prefix}{j}") for j in range(3, degree))
    moment_by_order = dict(zip(range(3, degree), moments))
    g = u + sum(moment_by_order[j] * u**j for j in range(3, degree))

    # Since the linear coefficient of g is one, the new coefficient of eta
    # enters g(eta) with coefficient one.  No symbolic solve is needed.
    eta = z
    for n in range(2, degree):
        known = _coefficient(
            sp.series(g.subs(u, eta), z, 0, n + 1).removeO(), z, n
        )
        eta = sp.expand(eta - known * z**n)

    e = sp.expand(eta / z)
    h_coefficients: list[sp.Expr] = [sp.Integer(1)]
    for n in range(1, degree - 1):
        value = -sum(
            _coefficient(e, z, j) * h_coefficients[n - j]
            for j in range(1, n + 1)
        )
        h_coefficients.append(sp.expand(value))
    h = sum(h_coefficients[j] * w**j for j in range(degree - 1))

    coefficients: list[sp.Expr] = [sp.Integer(0)] * (degree + 1)
    for j in range(2, degree - 1):
        coefficients[j] = h_coefficients[j]
    total = -sum(coefficients[j] for j in range(2, degree - 1))
    weighted = -1 - sum(
        j * coefficients[j] for j in range(2, degree - 1)
    )
    coefficients[degree - 1] = sp.expand(degree * total - weighted)
    coefficients[degree] = sp.expand(weighted - (degree - 1) * total)
    H = sp.expand(
        sum(coefficients[j] * w**j for j in range(2, degree + 1))
    )

    assert sp.expand(H.subs(w, 1)) == 0
    assert sp.expand(sp.diff(H, w).subs(w, 1)) == -1
    return OptimalMomentSeed(
        degree=degree,
        variable=w,
        series_variable=z,
        moments=moments,
        generating_series=g,
        inverse_series=eta,
        bridge_polynomial=sp.expand(1 + H),
        seed=H,
        coefficients=tuple(coefficients),
    )


@lru_cache(maxsize=None)
def _universal_power_root(
    last_degree: int, multiplicity: int
) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Expr, ...]]:
    """Return coefficients of ``(1+a_1 z+...)**(1/m)`` universally."""

    z = sp.Symbol("_root_z")
    a = sp.symbols(f"_root_a1:{last_degree + 1}")
    root = sp.Integer(1)
    coefficients: list[sp.Expr] = []
    for n in range(1, last_degree + 1):
        already = _coefficient(root**multiplicity, z, n)
        value = sp.expand((a[n - 1] - already) / multiplicity)
        coefficients.append(value)
        root += value * z**n
    return tuple(a), tuple(coefficients)


def power_prony_data(
    seed: OptimalMomentSeed, multiplicity: int
) -> PowerPronyData:
    """Build the pure-multiplicity Prony matrix in optimal moments.

    For ``N=m*r``, a vertical affine translate of ``H`` is an ``m``-th power
    iff the formal ``m``-th root of

        A(T)=T^N (H(1/T)-linear terms) / c_N

    truncates after degree ``r``.  Only coefficients through ``T^(N-2)`` are
    needed, so the test is intrinsic to the optimal moment vector.
    """

    N = seed.degree
    if multiplicity < 2 or N % multiplicity:
        raise ValueError("multiplicity must be a nontrivial divisor of the degree")
    last = N - 2
    leading = seed.leading_coefficient
    reversed_coefficients = tuple(
        sp.cancel(seed.coefficient(N - j) / leading)
        for j in range(1, last + 1)
    )

    universal_a, universal_root = _universal_power_root(last, multiplicity)
    substitution = dict(zip(universal_a, reversed_coefficients))
    root_coefficients = tuple(
        value.subs(substitution) for value in universal_root
    )

    # If L=A'/A=sum L_(j-1) T^(j-1), the root B=A^(1/m) satisfies
    # n*b_n=sum r_j*b_(n-j), r_j=L_(j-1)/m.  Its Bell determinant is the
    # leading minor of the following Toeplitz--Hessenberg matrix.
    log_coefficients: list[sp.Expr] = []
    a_with_constant = (sp.Integer(1),) + reversed_coefficients
    for n in range(1, last + 1):
        value = n * a_with_constant[n] - sum(
            a_with_constant[j] * log_coefficients[n - j - 1]
            for j in range(1, n)
        )
        log_coefficients.append(value)
    weights = tuple(value / multiplicity for value in log_coefficients)
    matrix = sp.Matrix(
        last,
        last,
        lambda row, column: (
            weights[row - column]
            if column <= row
            else -sp.Integer(row + 1)
            if column == row + 1
            else sp.Integer(0)
        ),
    )
    root_degree = N // multiplicity
    return PowerPronyData(
        multiplicity=multiplicity,
        root_degree=root_degree,
        reversed_coefficients=reversed_coefficients,
        logarithmic_weights=weights,
        hessenberg_matrix=matrix,
        root_coefficients=root_coefficients,
        equation_indices=tuple(range(root_degree + 1, last + 1)),
    )


def primitive_numerator(
    expression: sp.Expr, variables: tuple[sp.Symbol, ...]
) -> sp.Expr:
    """Clear denominators and scalar content without changing a principal ideal."""

    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.primitive(sp.Poly(numerator, *variables))[1].as_expr()


def _approximate_inner(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    outer_degree: int,
    inner_degree: int,
) -> sp.Expr:
    """Recover the monic approximate inner polynomial from top coefficients."""

    N = outer_degree * inner_degree
    B = variable**inner_degree
    for j in range(inner_degree - 1, 0, -1):
        target_degree = N - inner_degree + j
        target = _coefficient(polynomial, variable, target_degree)
        known = _coefficient(B**outer_degree, variable, target_degree)
        value = sp.cancel((target - known) / outer_degree)
        B = sp.expand(B + value * variable**j)
    return B


def ritt_krylov_data(
    polynomial: sp.Expr,
    variable: sp.Symbol,
    outer_degree: int,
    inner_degree: int,
) -> RittKrylovData:
    """Build the coefficient Krylov matrix for ``A o B`` modulo affine terms.

    The input must be monic of degree ``outer_degree*inner_degree``.  The
    columns are coefficient vectors of

        1, W, B, B^2, ..., B^a, f.

    Thus all maximal minors vanish exactly when ``f`` differs by an affine
    polynomial from an element of ``k[B]``.  Once ``B`` is the approximate
    root, the fixed-pivot minors in ``selected_minors`` generate the same
    determinantal ideal; the top ``b-1`` nonpivot minors vanish identically.
    """

    a, b = outer_degree, inner_degree
    N = a * b
    polynomial = sp.expand(polynomial)
    if _coefficient(polynomial, variable, N) != 1:
        raise ValueError("Ritt--Krylov input must be monic")
    B = _approximate_inner(polynomial, variable, a, b)
    columns = [sp.Integer(1), variable]
    columns.extend(sp.expand(B**j) for j in range(1, a + 1))
    columns.append(polynomial)
    matrix = sp.Matrix(
        N + 1,
        len(columns),
        lambda row, column: _coefficient(columns[column], variable, row),
    )
    pivot_rows = tuple([0, 1] + [j * b for j in range(1, a + 1)])
    selected_rows = tuple(
        row
        for row in range(2, N - b + 1)
        if row not in pivot_rows
    )
    selected_minors = tuple(
        sp.factor(
            matrix.extract(pivot_rows + (row,), range(matrix.cols)).det(
                method="domain-ge"
            )
        )
        for row in selected_rows
    )
    assert len(selected_minors) == N - a - b
    return RittKrylovData(
        outer_degree=a,
        inner_degree=b,
        polynomial=polynomial,
        approximate_inner=B,
        matrix=matrix,
        pivot_rows=pivot_rows,
        selected_rows=selected_rows,
        selected_minors=selected_minors,
    )


def moment_ritt_krylov_data(
    seed: OptimalMomentSeed, outer_degree: int, inner_degree: int
) -> RittKrylovData:
    """Build a Ritt--Krylov matrix directly from optimal moments."""

    if outer_degree * inner_degree != seed.degree:
        raise ValueError("factor degrees must multiply to the seed degree")
    w = seed.variable
    leading = seed.leading_coefficient
    generic_coefficients = tuple(
        sp.Symbol(f"_monic_c{j}") for j in range(2, seed.degree)
    )
    generic_polynomial = w**seed.degree + sum(
        generic_coefficients[j - 2] * w**j
        for j in range(2, seed.degree)
    )
    generic = ritt_krylov_data(
        generic_polynomial, w, outer_degree, inner_degree
    )
    substitution = {
        generic_coefficients[j - 2]: seed.coefficient(j) / leading
        for j in range(2, seed.degree)
    }
    return RittKrylovData(
        outer_degree=outer_degree,
        inner_degree=inner_degree,
        polynomial=generic.polynomial.subs(substitution),
        approximate_inner=generic.approximate_inner.subs(substitution),
        matrix=generic.matrix.subs(substitution),
        pivot_rows=generic.pivot_rows,
        selected_rows=generic.selected_rows,
        selected_minors=tuple(
            minor.subs(substitution) for minor in generic.selected_minors
        ),
    )


def fixed_weight_332_prony() -> FixedWeightPronyData:
    """Return the intrinsic marked Prony incidence for partition ``3+3+2``."""

    x = sp.Symbol("_weight_two_node")
    power_sums = sp.symbols("p1:7")
    elementary: list[sp.Expr] = [sp.Integer(1)]
    for n in range(1, 7):
        # After removing the distinguished weight-two node, divide by three.
        # The result must be the power-sum sequence of two unit-weight nodes.
        value = sum(
            (-1) ** (k - 1)
            * elementary[n - k]
            * (power_sums[k - 1] - 2 * x**k)
            / 3
            for k in range(1, n + 1)
        ) / n
        elementary.append(sp.factor(value))
    relations = tuple(
        primitive_numerator(elementary[n], (x,) + power_sums)
        for n in range(3, 7)
    )

    # Work in the free cubic algebra R[x]/(E3).  Multiplication by the
    # remainders of E4,E5,E6 gives a 3x9 presentation matrix.
    cubic = sp.Poly(relations[0], x)
    if cubic.LC() == 0:
        raise AssertionError("the marked Prony cubic must have constant leading term")
    basis = (sp.Integer(1), x, x**2)
    blocks: list[sp.Matrix] = []
    for relation in relations[1:]:
        remainder = sp.rem(relation, relations[0], x)
        block = sp.Matrix(
            3,
            3,
            lambda row, column: sp.Poly(
                sp.rem(remainder * basis[column], relations[0], x), x
            ).nth(row),
        )
        blocks.append(block)
    fitting = blocks[0].row_join(blocks[1]).row_join(blocks[2])
    return FixedWeightPronyData(
        distinguished_node=x,
        power_sums=power_sums,
        relations=relations,
        fitting_matrix=fitting,
    )


def reverse_power_sums(seed: OptimalMomentSeed) -> tuple[sp.Expr, ...]:
    """Return root power sums ``p_1,...,p_(N-2)`` from optimal moments.

    For ``A(T)=prod_i(1-r_i*T)=1+a_1*T+...``, Newton's recurrence is

        p_n+a_1*p_(n-1)+...+a_(n-1)*p_1+n*a_n=0.

    The last two coefficients of ``A`` correspond to the freely adjustable
    affine terms of an omitted polynomial, so the available range is exactly
    ``1,...,N-2``.
    """

    leading = seed.leading_coefficient
    a = (sp.Integer(1),) + tuple(
        seed.coefficient(seed.degree - j) / leading
        for j in range(1, seed.degree - 1)
    )
    power_sums: list[sp.Expr] = []
    for n in range(1, seed.degree - 1):
        value = -n * a[n] - sum(
            a[j] * power_sums[n - j - 1] for j in range(1, n)
        )
        power_sums.append(value)
    return tuple(power_sums)


def sylvester_matrix(
    left: sp.Expr, right: sp.Expr, variable: sp.Symbol
) -> sp.Matrix:
    """Return the coefficient Sylvester matrix of two univariate polynomials."""

    left_poly = sp.Poly(left, variable)
    right_poly = sp.Poly(right, variable)
    left_degree = left_poly.degree()
    right_degree = right_poly.degree()
    width = left_degree + right_degree
    rows: list[list[sp.Expr]] = []
    for shift in range(right_degree):
        multiple = sp.Poly(variable**shift * left, variable)
        rows.append([multiple.nth(j) for j in range(width)])
    for shift in range(left_degree):
        multiple = sp.Poly(variable**shift * right, variable)
        rows.append([multiple.nth(j) for j in range(width)])
    return sp.Matrix(rows)


def rank_drop_minors(
    matrix: sp.Matrix, rank_bound: int
) -> Iterator[sp.Expr]:
    """Yield all ``(rank_bound+1)`` minors of ``matrix`` lazily."""

    size = rank_bound + 1
    if size > min(matrix.rows, matrix.cols):
        return
    for rows in combinations(range(matrix.rows), size):
        for columns in combinations(range(matrix.cols), size):
            yield matrix.extract(rows, columns).det(method="domain-ge")


def fixed_weight_hankel_prony(
    multiplicities: tuple[int, ...],
    power_sums: tuple[sp.Expr, ...] | None = None,
    *,
    variable: sp.Symbol | None = None,
) -> FixedWeightHankelData:
    """Build unmarked Hankel/subresultant equations for a mixed partition.

    Let ``r=len(multiplicities)`` and ``p_0=sum(multiplicities)``.  At least
    ``p_1,...,p_(2r-1)`` are required.  The Hankel pencil recovers the support
    polynomial on ``det(H)!=0``.  For a candidate multiplicity ``m``, the
    Christoffel polynomial

        m*v(T)^T*adj(H)*v(T)-det(H)

    shares exactly the nodes of weight ``m`` with the support polynomial.
    Requiring the corresponding Sylvester rank drops, plus the recurrence for
    any later supplied moments, gives the fixed-weight Prony locus.  Its
    reduced collision closure is obtained by saturating these determinantal
    equations by ``det(H)``.
    """

    if not multiplicities or any(m < 2 for m in multiplicities):
        raise ValueError("omitted-partition multiplicities must all be at least two")
    r = len(multiplicities)
    minimum_last = 2 * r - 1
    if power_sums is None:
        power_sums = tuple(sp.symbols(f"p1:{minimum_last + 1}"))
    if len(power_sums) < minimum_last:
        raise ValueError(
            f"the Hankel pencil needs power sums through p_{minimum_last}"
        )
    T = variable or sp.Symbol("_prony_T")
    sequence = (sp.Integer(sum(multiplicities)),) + tuple(power_sums)
    H = sp.Matrix(r, r, lambda i, j: sequence[i + j])
    H1 = sp.Matrix(r, r, lambda i, j: sequence[i + j + 1])
    determinant = H.det(method="domain-ge")
    support = (T * H - H1).det(method="domain-ge")
    vector = sp.Matrix([T**j for j in range(r)])
    christoffel = (vector.T * H.adjugate() * vector)[0]

    conditions: list[WeightSubresultantCondition] = []
    for weight, count in sorted(Counter(multiplicities).items()):
        weight_polynomial = weight * christoffel - determinant
        if sp.expand(weight_polynomial) == 0:
            # For a one-node partition the sole weight is already fixed by
            # p_0=N, so only the support recurrence remains.
            continue
        matrix = sylvester_matrix(support, weight_polynomial, T)
        conditions.append(
            WeightSubresultantCondition(
                weight=weight,
                required_gcd_degree=count,
                weight_polynomial=weight_polynomial,
                sylvester_matrix=matrix,
                rank_bound=matrix.rows - count,
            )
        )

    support_poly = sp.Poly(support, T)
    extension_equations: list[sp.Expr] = []
    for n in range(2 * r, len(sequence)):
        recurrence = sum(
            support_poly.nth(j) * sequence[n - r + j]
            for j in range(r + 1)
        )
        extension_equations.append(recurrence)
    return FixedWeightHankelData(
        multiplicities=tuple(multiplicities),
        power_sums=tuple(power_sums),
        variable=T,
        hankel_matrix=H,
        shifted_hankel_matrix=H1,
        hankel_determinant=determinant,
        support_polynomial=support,
        christoffel_numerator=christoffel,
        weight_conditions=tuple(conditions),
        extension_equations=tuple(extension_equations),
    )


def moment_fixed_weight_hankel(
    seed: OptimalMomentSeed, multiplicities: tuple[int, ...]
) -> FixedWeightHankelData:
    """Evaluate the unmarked fixed-weight matrices in optimal moments."""

    if sum(multiplicities) != seed.degree:
        raise ValueError("multiplicities must sum to the seed degree")
    power_sums = reverse_power_sums(seed)
    generic_power_sums = tuple(
        sp.symbols(f"_generic_p1:{len(power_sums) + 1}")
    )
    generic = fixed_weight_hankel_prony(
        multiplicities,
        generic_power_sums,
        variable=sp.Symbol("_moment_prony_T"),
    )
    substitution = dict(zip(generic_power_sums, power_sums))
    conditions = tuple(
        WeightSubresultantCondition(
            weight=condition.weight,
            required_gcd_degree=condition.required_gcd_degree,
            weight_polynomial=condition.weight_polynomial.subs(substitution),
            sylvester_matrix=condition.sylvester_matrix.subs(substitution),
            rank_bound=condition.rank_bound,
        )
        for condition in generic.weight_conditions
    )
    return FixedWeightHankelData(
        multiplicities=tuple(multiplicities),
        power_sums=power_sums,
        variable=generic.variable,
        hankel_matrix=generic.hankel_matrix.subs(substitution),
        shifted_hankel_matrix=generic.shifted_hankel_matrix.subs(substitution),
        hankel_determinant=generic.hankel_determinant.subs(substitution),
        support_polynomial=generic.support_polynomial.subs(substitution),
        christoffel_numerator=generic.christoffel_numerator.subs(substitution),
        weight_conditions=conditions,
        extension_equations=tuple(
            equation.subs(substitution)
            for equation in generic.extension_equations
        ),
    )


def hessenberg_determinant_identity(size: int) -> sp.Expr:
    """Return ``det(K_n)-n!*b_n`` for the universal log-Prony recurrence."""

    r = sp.symbols(f"r1:{size + 1}")
    matrix = sp.Matrix(
        size,
        size,
        lambda row, column: (
            r[row - column]
            if column <= row
            else -sp.Integer(row + 1)
            if column == row + 1
            else sp.Integer(0)
        ),
    )
    b = [sp.Integer(1)]
    for n in range(1, size + 1):
        b.append(
            sp.expand(sum(r[j - 1] * b[n - j] for j in range(1, n + 1)) / n)
        )
    return sp.expand(matrix.det(method="domain-ge") - factorial(size) * b[size])
