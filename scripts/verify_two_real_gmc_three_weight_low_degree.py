#!/usr/bin/env python3
"""Exact low-degree closure of the {-1,0,1} GMC(2) support.

For U=ZW, write

    P = W*A(U) + C(U) + Z*B(U),       D(U) = U*A(U)*B(U).

Constant-term extraction in the rotational variable gives both

    E(exp(tP)) = L(exp(tC) I_0(2*t*sqrt(D)))

and the simpler ordinary moment series

    sum_m E(P**m)t**m
      = L(((1-t*C)**2 - 4*t**2*D)**(-1/2)),

where L(U**j)=j!.  In degree n, the invariant coordinates are

    C = sum_{j=1}^{floor(n/2)} c_j*(U**j-j!),
    D = sum_{j=1}^{1+2*floor((n-1)/2)} d_j*U**j.

Conversely, over C every such nonzero D factors as U*A*B with the required
degree bounds.  Thus C,D exactly parameterize the three-weight problem
after quotienting the relative A/B scaling.

For degrees four, five, and six, this verifier covers C != 0 and D != 0
by coefficient charts.  Weighted scalar multiplication normalizes one
nonzero c_j to one; a selected nonzero d_j is localized.  The second
moment eliminates the highest D coefficient.  Every resulting ideal has
the literal reduced characteristic-zero basis [1].
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch import msolve


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_real_gmc_three_weight_low_degree.json"
)
DEGREE_CUTOFFS = {4: 6, 5: 8, 6: 9}


def factorial_functional(polynomial: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    """Apply L(variable**j)=j! to a univariate polynomial."""

    return sp.expand(
        sum(
            coefficient * math.factorial(exponent[0])
            for exponent, coefficient in sp.Poly(
                sp.expand(polynomial), variable
            ).terms()
        )
    )


def bessel_moment(
    order: int,
    radial: sp.Symbol,
    centered: sp.Expr,
    circuit: sp.Expr,
) -> sp.Expr:
    """Return E(P**order) as a polynomial in the C,D coefficients."""

    return sp.expand(
        sum(
            math.factorial(order)
            // (
                math.factorial(pair_count) ** 2
                * math.factorial(order - 2 * pair_count)
            )
            * factorial_functional(
                centered ** (order - 2 * pair_count)
                * circuit**pair_count,
                radial,
            )
            for pair_count in range(order // 2 + 1)
        )
    )


def sparse_invariant_moments(
    c_coefficients: tuple[sp.Symbol, ...],
    d_coefficients: tuple[sp.Symbol, ...],
    cutoff: int,
) -> dict[int, sp.Expr]:
    """Generate invariant moments without expanding powers in SymPy."""

    coefficients = c_coefficients + d_coefficients
    coefficient_count = len(coefficients)
    zero_exponents = (0,) * coefficient_count
    c_terms: list[tuple[int, tuple[int, ...], int]] = []
    d_terms: list[tuple[int, tuple[int, ...], int]] = []

    for power in range(1, len(c_coefficients) + 1):
        exponents = [0] * coefficient_count
        exponents[power - 1] = 1
        exponent_tuple = tuple(exponents)
        c_terms.extend(
            (
                (power, exponent_tuple, 1),
                (0, exponent_tuple, -math.factorial(power)),
            )
        )
    offset = len(c_coefficients)
    for power in range(1, len(d_coefficients) + 1):
        exponents = [0] * coefficient_count
        exponents[offset + power - 1] = 1
        d_terms.append((power, tuple(exponents), 1))

    SparsePolynomial = dict[tuple[int, tuple[int, ...]], int]

    def multiply(
        left: SparsePolynomial, right: SparsePolynomial
    ) -> SparsePolynomial:
        product: SparsePolynomial = {}
        for (left_u, left_exponents), left_coefficient in left.items():
            for (
                right_u,
                right_exponents,
            ), right_coefficient in right.items():
                key = (
                    left_u + right_u,
                    tuple(
                        left_power + right_power
                        for left_power, right_power in zip(
                            left_exponents, right_exponents
                        )
                    ),
                )
                product[key] = (
                    product.get(key, 0)
                    + left_coefficient * right_coefficient
                )
        return {
            key: coefficient
            for key, coefficient in product.items()
            if coefficient
        }

    def powers(
        terms: list[tuple[int, tuple[int, ...], int]],
        maximum: int,
    ) -> list[SparsePolynomial]:
        base = {
            (u_degree, exponents): coefficient
            for u_degree, exponents, coefficient in terms
        }
        result: list[SparsePolynomial] = [
            {(0, zero_exponents): 1}
        ]
        for _ in range(maximum):
            result.append(multiply(result[-1], base))
        return result

    c_powers = powers(c_terms, cutoff)
    d_powers = powers(d_terms, cutoff // 2)
    moments: dict[int, sp.Expr] = {1: sp.Integer(0)}
    for order in range(2, cutoff + 1):
        accumulated: dict[tuple[int, ...], int] = {}
        for pair_count in range(order // 2 + 1):
            c_power = order - 2 * pair_count
            outer_coefficient = math.factorial(order) // (
                math.factorial(pair_count) ** 2
                * math.factorial(c_power)
            )
            product = multiply(
                c_powers[c_power], d_powers[pair_count]
            )
            for (
                u_degree,
                exponents,
            ), coefficient in product.items():
                accumulated[exponents] = (
                    accumulated.get(exponents, 0)
                    + outer_coefficient
                    * coefficient
                    * math.factorial(u_degree)
                )
        moments[order] = sp.Add(
            *(
                coefficient
                * sp.prod(
                    variable**exponent
                    for variable, exponent in zip(
                        coefficients, exponents
                    )
                )
                for exponents, coefficient in accumulated.items()
                if coefficient
            )
        )
    return moments


def direct_wick_moment(
    order: int,
    polynomial: sp.Expr,
    z: sp.Symbol,
    w: sp.Symbol,
) -> sp.Expr:
    """Compute E(P**order) directly from E(Z**a W**b)=delta(a,b)a!."""

    expanded = sp.Poly(sp.expand(polynomial**order), z, w)
    return sp.expand(
        sum(
            coefficient * math.factorial(z_degree)
            for (z_degree, w_degree), coefficient in expanded.terms()
            if z_degree == w_degree
        )
    )


def identity_regressions() -> None:
    """Check the two generating identities independently in finite order."""

    t, abstract_c, abstract_d = sp.symbols("t abstract_c abstract_d")
    ordinary = sp.series(
        ((1 - t * abstract_c) ** 2 - 4 * t**2 * abstract_d)
        ** sp.Rational(-1, 2),
        t,
        0,
        10,
    ).removeO()
    for order in range(10):
        expected = sum(
            sp.Rational(
                math.factorial(order),
                math.factorial(pair_count) ** 2
                * math.factorial(order - 2 * pair_count),
            )
            * abstract_c ** (order - 2 * pair_count)
            * abstract_d**pair_count
            for pair_count in range(order // 2 + 1)
        )
        assert sp.expand(ordinary.coeff(t, order) - expected) == 0

    # Direct Wick check with generic quartic A,B,C, independent of the
    # constant-term derivation above.
    u, z, w = sp.symbols("u z w")
    a0, a1, b0, b1, c1, c2 = sp.symbols(
        "a0 a1 b0 b1 c1 c2"
    )
    d1, d2, d3 = sp.symbols("d1 d2 d3")
    a_polynomial = a0 + a1 * u
    b_polynomial = b0 + b1 * u
    centered = c1 * (u - 1) + c2 * (u**2 - 2)
    circuit = sp.expand(u * a_polynomial * b_polynomial)
    generic_circuit = d1 * u + d2 * u**2 + d3 * u**3
    sparse_moments = sparse_invariant_moments(
        (c1, c2), (d1, d2, d3), 6
    )
    for order in range(1, 7):
        assert (
            sp.expand(
                sparse_moments[order]
                - bessel_moment(
                    order, u, centered, generic_circuit
                )
            )
            == 0
        )
    direct_polynomial = (
        w * (a0 + a1 * z * w)
        + c1 * (z * w - 1)
        + c2 * ((z * w) ** 2 - 2)
        + z * (b0 + b1 * z * w)
    )
    for order in range(1, 7):
        assert (
            sp.expand(
                direct_wick_moment(order, direct_polynomial, z, w)
                - bessel_moment(order, u, centered, circuit)
            )
            == 0
        )


def degree_systems(degree: int, cutoff: int) -> list[dict[str, object]]:
    """Build and solve every invariant coefficient chart in one degree."""

    radial_degree = degree // 2
    side_degree = (degree - 1) // 2
    circuit_degree = 1 + 2 * side_degree
    c_coefficients = sp.symbols(f"c1:{radial_degree + 1}")
    d_coefficients = sp.symbols(f"d1:{circuit_degree + 1}")
    rho = sp.Symbol("rho")

    moments = sparse_invariant_moments(
        c_coefficients, d_coefficients, cutoff
    )
    assert moments[1] == 0

    records: list[dict[str, object]] = []
    eliminated = d_coefficients[-1]
    for selected_c in c_coefficients:
        normalization = {selected_c: 1}
        solution = sp.solve(
            moments[2].subs(normalization), eliminated
        )[0]
        variables = (
            tuple(
                coefficient
                for coefficient in c_coefficients
                if coefficient != selected_c
            )
            + d_coefficients[:-1]
            + (rho,)
        )
        base_equations = [
            sp.together(
                moments[order]
                .subs(normalization)
                .subs(eliminated, solution)
            ).as_numer_denom()[0]
            for order in range(3, cutoff + 1)
        ]

        for selected_d in d_coefficients:
            localization = sp.factor(
                selected_d.subs(eliminated, solution)
            )
            equations = base_equations + [
                sp.expand(rho * localization - 1)
            ]
            input_text = msolve.input_text(
                equations, variables, characteristic=0
            )
            result = msolve.run(
                equations,
                variables,
                prime=0,
                threads=4,
                timeout=600,
            )
            assert result.returncode == 0
            assert result.contains_one
            records.append(
                {
                    "degree": degree,
                    "selected_c": str(selected_c),
                    "selected_d": str(selected_d),
                    "variables": list(map(str, variables)),
                    "moment_cutoff": cutoff,
                    "characteristic": 0,
                    "eliminated_coefficient": str(eliminated),
                    "elimination_solution": str(solution),
                    "reduced_basis": "[1]",
                    "input_sha256": hashlib.sha256(
                        input_text.encode()
                    ).hexdigest(),
                }
            )
    expected_charts = radial_degree * circuit_degree
    assert len(records) == expected_charts
    return records


def main() -> None:
    if not msolve.available():
        raise RuntimeError("msolve executable not found")

    identity_regressions()
    all_records: list[dict[str, object]] = []
    degree_summary: dict[str, object] = {}
    for degree, cutoff in DEGREE_CUTOFFS.items():
        records = degree_systems(degree, cutoff)
        all_records.extend(records)
        degree_summary[str(degree)] = {
            "moment_cutoff": cutoff,
            "chart_count": len(records),
            "all_reduced_bases": "[1]",
        }
        print(
            f"PASS degree {degree}: {len(records)} exact QQ "
            f"unit ideals through moment {cutoff}"
        )

    payload = {
        "format": "two-real-gmc-three-weight-low-degree-v1",
        "support": [-1, 0, 1],
        "ordinary_moment_series": (
            "L(((1-t*C)^2-4*t^2*D)^(-1/2))"
        ),
        "exponential_moment_series": (
            "L(exp(t*C)*I_0(2*t*sqrt(D)))"
        ),
        "degrees": degree_summary,
        "charts": all_records,
        "conclusion": (
            "no total-degree-at-most-six pure-moment-zero polynomial "
            "has all three weight components -1,0,1 nonzero"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")

    print("PASS {-1,0,1}: ordinary and exponential moment identities")
    print(f"PASS {{-1,0,1}}: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
