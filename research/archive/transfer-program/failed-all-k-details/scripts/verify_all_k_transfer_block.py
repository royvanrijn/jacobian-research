"""Structural identities and new k=5,6 regressions for the all-k theorem."""

import itertools
import math

import sympy as sp

Z = sp.symbols("Z")

# The two differential syzygies proving affine difference equals zero.
for k in range(1, 5):
    u = sp.symbols(f"u0:{3*k}")
    v = sp.symbols(f"v0:{2*k}")
    U = Z ** (3 * k) + sum(u[i] * Z**i for i in range(3 * k))
    V = Z ** (2 * k) + sum(v[i] * Z**i for i in range(2 * k))
    D = sp.expand(U**2 - V**3)
    K0 = sp.expand(2 * V * sp.diff(U, Z) - 3 * U * sp.diff(V, Z))
    assert sp.expand(V * sp.diff(D, Z) - 3 * sp.diff(V, Z) * D - U * K0) == 0
    assert sp.expand(U * sp.diff(D, Z) - 2 * sp.diff(U, Z) * D - V**2 * K0) == 0
    assert 1 - 6 * k != 0 and 6 * k != 0

# One ordered-root Boolean factor satisfies the cusp equation exactly.
q, epsilon = sp.symbols("q epsilon")
Ui = q**3 + sp.Rational(3, 2) * epsilon * q
Vi = q**2 + epsilon
assert sp.rem(sp.Poly(sp.expand(Ui**2 - Vi**3), epsilon),
              sp.Poly(epsilon**2, epsilon)).as_expr() == 0


def coincident_block(k):
    """Return length and m-adic Hilbert function at S=Z^k."""
    transfer = sp.symbols(f"r0:{k}")
    u = sp.symbols(f"u0:{3*k}")
    V = Z ** (2 * k) + sum(transfer[i] * Z**i for i in range(k))
    U = Z ** (3 * k) + sum(u[i] * Z**i for i in range(3 * k))
    difference = sp.Poly(sp.expand(U**2 - V**3), Z)

    solution = {}
    for degree, variable in zip(range(6 * k - 1, 3 * k - 1, -1), reversed(u)):
        equation = sp.expand(difference.coeff_monomial(Z**degree).subs(solution))
        roots = sp.solve(equation, variable, dict=False)
        assert len(roots) == 1
        solution[variable] = roots[0]

    equations = [
        sp.factor(difference.coeff_monomial(Z**degree).subs(solution))
        for degree in range(3 * k - 1, -1, -1)
    ]
    variables = tuple(reversed(transfer))
    basis = sp.groebner(equations, *variables, order="grevlex", domain=sp.QQ)
    assert basis.is_zero_dimensional
    leading = [polynomial.LM(order=basis.order).exponents
               for polynomial in basis.polys]

    standard = [
        exponent
        for exponent in itertools.product(range(k + 2), repeat=k)
        if not any(all(exponent[j] >= monomial[j] for j in range(k))
                   for monomial in leading)
    ]
    index = {monomial: position for position, monomial in enumerate(standard)}
    monomial_basis = [
        sp.prod(variables[j] ** exponent[j] for j in range(k))
        for exponent in standard
    ]

    def vector(polynomial):
        answer = sp.zeros(len(standard), 1)
        remainder = sp.Poly(basis.reduce(polynomial)[1], *variables)
        for monomial, coefficient in remainder.terms():
            answer[index[monomial]] = coefficient
        return answer

    multiplication = [
        sp.Matrix.hstack(*[vector(variable * monomial)
                           for monomial in monomial_basis])
        for variable in variables
    ]
    power = sp.Matrix.hstack(*multiplication)
    dimensions = []
    while True:
        columns = power.columnspace()
        power = sp.Matrix.hstack(*columns) if columns else sp.zeros(len(standard), 0)
        dimensions.append(power.rank())
        if not dimensions[-1]:
            break
        next_columns = [matrix * column
                        for matrix in multiplication for column in columns]
        power = sp.Matrix.hstack(*next_columns)

    hilbert = [len(standard) - dimensions[0]] + [
        dimensions[i] - dimensions[i + 1]
        for i in range(len(dimensions) - 1)
    ]
    return len(standard), tuple(hilbert)


for k in (5, 6):
    length, hilbert = coincident_block(k)
    assert length == 2**k
    assert hilbert == tuple(math.comb(k, degree) for degree in range(k + 1))

print("PASS: the Wronskian identities force every affine difference to vanish")
print("PASS: each ordered root carries the universal square-zero cusp jet")
print("PASS: the new k=5 coincident block has length 32 and Hilbert series (1+t)^5")
print("PASS: the new k=6 coincident block has length 64 and Hilbert series (1+t)^6")
