#!/usr/bin/env python3
"""Exact exclusion on a swap-antisymmetric rank-two SIC slice.

For a binary quartic P, consider

    F_P = xi1^4 P(z1,z2) - xi2^4 P(z2,z1).

Odd pure contractions vanish by the simultaneous swap involution.  This
checker proves that the five even contractions of orders 2,4,6,8,10 force
P to be a fourth power (z2+z1)^4 or (z2-z1)^4, up to scale.  At both
points the two coefficient rows are proportional, so no exact-rank-two
SIC counterexample occurs on this slice.

The characteristic-zero chart calculation uses a reconstructed lex ideal
whose exact containment and colength are checked over Q.  A complete
projective boundary exclusion is certified at one good prime.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_rank_two_swap_slice.json"
)
EXPLORER = (
    ROOT
    / "scripts"
    / "explore_two_pair_sic_bidegree44_rank_two_swap_slice.py"
)
SPEC = importlib.util.spec_from_file_location("swap_slice_explorer", EXPLORER)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load swap-slice explorer")
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

a = module.a
y = module.y
PRIMES = [32003, 32009, 32027, 32029, 32051, 32057, 32059]
GOOD_PRIME = PRIMES[0]


def exact_chart_moment(order: int) -> sp.Expr:
    expression = module.modular_moment_expression(order, None, 4)
    local_symbols = {str(variable): variable for variable in a[:4]}
    return sp.expand(
        sp.sympify(expression.replace("^", "**"), locals=local_symbols)
    )


def special_fiber_certificate() -> tuple[int, dict[str, bool]]:
    base_orders = (2, 4, 6, 8, 10)
    chart_polynomials = [
        module.modular_moment_expression(order, GOOD_PRIME, 4)
        for order in base_orders
    ]
    script = "\n".join(
        [
            f"ring r={GOOD_PRIME},(a0,a1,a2,a3),dp;",
            "option(redSB);",
            f"ideal I={','.join(chart_polynomials)};",
            "ideal G=slimgb(I);",
            'print("VDIM="+string(vdim(G)));',
            "quit;",
        ]
    )
    result = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        check=True,
    )
    marker = next(
        line for line in result.stdout.splitlines() if line.startswith("VDIM=")
    )
    colength = int(marker.split("=", 1)[1])

    boundary: dict[str, bool] = {}
    for chart in range(4):
        variables = [
            str(variable)
            for index, variable in enumerate(a)
            if index != chart
        ]
        polynomials = [
            module.modular_moment_expression(order, GOOD_PRIME, chart)
            for order in base_orders
        ]
        script = "\n".join(
            [
                f"ring r={GOOD_PRIME},({','.join(variables)}),dp;",
                "option(redSB);",
                f"ideal I={','.join(polynomials)},a4;",
                "ideal G=slimgb(I);",
                'if (size(G)==1 && G[1]==1) { print("UNIT"); }',
                "quit;",
            ]
        )
        result = subprocess.run(
            ["Singular", "-q"],
            input=script,
            text=True,
            capture_output=True,
            check=True,
        )
        boundary[f"a{chart}=1,a4=0"] = "UNIT" in result.stdout.splitlines()
    return colength, boundary


def main() -> None:
    if not all(sp.isprime(prime) for prime in PRIMES):
        raise AssertionError("reconstruction list contains a composite")

    reconstructed = module.reconstruct_lex_basis(PRIMES)
    basis = sp.groebner(reconstructed, *a[:4], order="lex")
    assert len(basis.polys) == 5

    # Exact characteristic-zero containment: every defining moment reduces
    # to zero modulo the reconstructed triangular ideal.
    reductions: dict[str, bool] = {}
    for order in (2, 4, 6, 8, 10, 12, 14, 16):
        remainder = basis.reduce(exact_chart_moment(order))[1]
        reductions[str(order)] = remainder == 0
    assert all(reductions.values())

    h = a[3] ** 2 - 16
    assert sp.expand(reconstructed[0] - h**5) == 0
    for sign in (-1, 1):
        point_a3 = 4 * sign
        assert sp.factor(
            reconstructed[2].subs(a[3], point_a3)
        ) == (a[2] - 6) ** 2
        assert sp.factor(
            reconstructed[3].subs({a[3]: point_a3, a[2]: 6})
        ) == a[1] - point_a3
        assert sp.factor(
            reconstructed[4].subs({a[3]: point_a3, a[2]: 6})
        ) == a[0] - 1

    # The leading monomials are a0, a1, a2^2, a2*a3^2, a3^10.
    # Hence the quotient basis is 1,a3,...,a3^9,a2,a2*a3.
    leading_monomials = [
        polynomial.LM(order=basis.order).exponents
        for polynomial in basis.polys
    ]
    expected_leading = {
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 2, 0),
        (0, 0, 1, 2),
        (0, 0, 0, 10),
    }
    assert set(leading_monomials) == expected_leading
    rational_colength = 12

    special_colength, boundary = special_fiber_certificate()
    assert special_colength == rational_colength
    assert all(boundary.values())

    # The reduced points give P=(z2 +/- z1)^4.  The coefficient rows P and
    # -reverse(P) are proportional, so the coefficient rank is one.
    reduced_points: list[list[int]] = []
    for sign in (-1, 1):
        coefficients = [1, 4 * sign, 6, 4 * sign, 1]
        matrix = sp.Matrix(
            [
                coefficients,
                [-value for value in reversed(coefficients)],
            ]
        )
        assert matrix.rank() == 1
        reduced_points.append(coefficients)

    # Sharper parity-even sub-slice: moments 2,4,6 have no nonzero
    # projective zero.  Three affine charts cover P=a0+a2*y^2+a4*y^4.
    even_first = a[0] + a[2] * y**2 + a[4] * y**4
    even_second = -(a[4] + a[2] * y**2 + a[0] * y**4)
    even_moments = [
        module.moment_from(even_first, even_second, order)
        for order in (2, 4, 6)
    ]
    parity_even_charts: dict[str, bool] = {}
    for chart in (0, 2, 4):
        variables = [a[index] for index in (0, 2, 4) if index != chart]
        chart_basis = sp.groebner(
            [value.subs(a[chart], 1) for value in even_moments],
            *variables,
            order="grevlex",
        )
        is_unit = list(chart_basis) == [1]
        parity_even_charts[f"a{chart}=1"] = is_unit
    assert all(parity_even_charts.values())

    artifact = {
        "format": "two-pair-sic-bidegree44-rank-two-swap-slice-v1",
        "field": "characteristic zero",
        "family": (
            "F_P=xi1^4*P(z1,z2)-xi2^4*P(z2,z1), "
            "P=sum(a_j*z1^j*z2^(4-j))"
        ),
        "odd_moments": (
            "identically zero by the simultaneous xi/z coordinate swap"
        ),
        "defining_even_moment_orders": [2, 4, 6, 8, 10],
        "a4_nonzero_chart": {
            "normalization": "a4=1",
            "rational_colength": rational_colength,
            "special_fiber_prime": GOOD_PRIME,
            "special_fiber_colength": special_colength,
            "reconstruction_primes": PRIMES,
            "tail_reductions_zero": reductions,
            "lex_basis": [str(value) for value in reconstructed],
        },
        "projective_boundary_a4_zero": boundary,
        "radical_points": reduced_points,
        "radical_interpretation": [
            "P=(z2-z1)^4",
            "P=(z2+z1)^4",
        ],
        "coefficient_rank_at_every_reduced_point": 1,
        "conclusion": (
            "no exact-rank-two SIC counterexample occurs on the "
            "swap-antisymmetric two-fourth-power-symbol slice"
        ),
        "parity_even_subslice": {
            "P": "a0*z2^4+a2*z1^2*z2^2+a4*z1^4",
            "orders": [2, 4, 6],
            "projective_chart_unit_bases": parity_even_charts,
            "conclusion": "only the zero quartic survives",
        },
        "scope": (
            "This is an exact orbit-slice exclusion, not a global "
            "rank-two theorem."
        ),
        "written_source": (
            "extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS swap slice: odd moments vanish by involution")
    print("PASS swap slice: even orders 2,4,6,8,10 have colength 12")
    print("PASS swap slice: radical consists of two rank-one fourth powers")
    print("PASS swap slice: complete projective a4=0 boundary is empty")
    print("PASS parity-even core: orders 2,4,6 have no projective survivor")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
