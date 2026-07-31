#!/usr/bin/env python3
"""Verify the six-support bidegree-(3,3) Rodrigues survivor.

The form checked here is

    F = y + y^3 - x - 2*x^2*y - x^2*y^3 - x^3*y^2

in the dehomogenized coefficient convention
``F_C(x,y) = sum(c_ij*x^i*y^j)``.  Its pure contractions vanish in every
positive order, although its coefficient matrix is invertible.  The same
Rodrigues identity proves the stronger and logically separate fact that
this particular moment--nullcone counterexample is SIC-safe.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from itertools import product
import json
from math import comb, factorial
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rodrigues_survivor.json"
)

X, Y, U, T, Z = sp.symbols("x y u t z")
A = 1 - T

COEFFICIENT_MATRIX = sp.Matrix(
    [
        [0, 1, 0, 1],
        [-1, 0, 0, 0],
        [0, -2, 0, -1],
        [0, 0, -1, 0],
    ]
)

F = sum(
    COEFFICIENT_MATRIX[i, j] * X**i * Y**j
    for i in range(4)
    for j in range(4)
)
P = sp.expand(A**3 * F.subs({X: U * T / A, Y: U**-1}))
P_FACTORED = U**-3 * (U**2 + A) * (A**2 - T * U**2)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--through", type=int, default=30)
    parser.add_argument("--mixed-degree", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


@lru_cache(maxsize=None)
def expanded_power(order: int) -> sp.Poly:
    return sp.Poly(sp.expand(F**order), X, Y, domain=sp.ZZ)


def dehomogeneous_moment(order: int) -> int:
    expanded = expanded_power(order)
    answer = 0
    for diagonal_degree in range(3 * order + 1):
        answer += (
            factorial(3 * order - diagonal_degree)
            * factorial(diagonal_degree)
            * int(expanded.coeff_monomial(X**diagonal_degree * Y**diagonal_degree))
        )
    return answer


def constant_term(expression: sp.Expr) -> sp.Expr:
    """Return the Laurent constant term in U."""

    expanded = sp.expand(expression)
    return sp.expand(
        sum(
            term
            for term in sp.Add.make_args(expanded)
            if term.as_powers_dict().get(U, 0) == 0
        )
    )


def direct_mixed_contraction(order: int, degree: int, i: int, j: int) -> int:
    """Contract M_ij of bidegree (degree,degree) against F**order."""

    expanded = expanded_power(order)
    total_degree = 3 * order + degree
    answer = 0
    for x_degree in range(3 * order + 1):
        y_degree = x_degree + i - j
        if not 0 <= y_degree <= 3 * order:
            continue
        coefficient = int(
            expanded.coeff_monomial(X**x_degree * Y**y_degree)
        )
        if coefficient == 0:
            continue
        answer += (
            coefficient
            * factorial(3 * order - y_degree + degree - j)
            * factorial(y_degree + j)
        )
    return answer


def beta_mixed_contraction(order: int, degree: int, i: int, j: int) -> int:
    """Evaluate the same contraction by exact beta/binomial extraction."""

    phase = i - j
    if (3 * order - phase) % 2:
        return 0
    target = (3 * order - phase) // 2
    normalized = Fraction(0)
    for second_z_degree in range(order + 1):
        first_z_degree = target - second_z_degree
        if not 0 <= first_z_degree <= order:
            continue
        t_degree = second_z_degree + i
        a_degree = (
            3 * order
            - first_z_degree
            - 2 * second_z_degree
            + degree
            - i
        )
        if a_degree < 0:
            continue
        normalized += Fraction(
            (-1) ** second_z_degree
            * comb(order, first_z_degree)
            * comb(order, second_z_degree)
            * factorial(t_degree)
            * factorial(a_degree),
            factorial(t_degree + a_degree + 1),
        )
    raw = normalized * factorial(3 * order + degree + 1)
    if raw.denominator != 1:
        raise AssertionError("beta expression did not recover an integer")
    return raw.numerator


def verify_rodrigues_identity(max_degree: int) -> list[dict[str, int]]:
    checks = []
    for phase in range(-max_degree, max_degree + 1):
        for order in range(max(1, abs(phase)), max_degree + 5):
            if (order - phase) % 2:
                continue
            r = (order - phase) // 2
            s = (order + phase) // 2
            coefficient = sp.Poly(
                (Z + A) ** order * (A**2 - T * Z) ** order,
                Z,
            ).coeff_monomial(Z ** (2 * r + s))
            rodrigues = (
                (-1) ** r
                * A**phase
                * sp.diff(T**order * A**order, T, s)
                / factorial(s)
            )
            if sp.cancel(coefficient - rodrigues) != 0:
                raise AssertionError(
                    f"Rodrigues mismatch for phase={phase}, order={order}"
                )
            checks.append({"phase": phase, "order": order})
    return checks


def verify_mixed_formula(max_degree: int) -> dict[str, object]:
    direct_checks = []
    cutoff_checks = []
    for degree in range(max_degree + 1):
        for i, j in product(range(degree + 1), repeat=2):
            phase = i - j
            for order in range(1, min(5, 3 * degree + 3) + 1):
                direct = direct_mixed_contraction(order, degree, i, j)
                beta = beta_mixed_contraction(order, degree, i, j)
                if direct != beta:
                    raise AssertionError(
                        "mixed beta mismatch at "
                        f"(m,e,i,j)=({order},{degree},{i},{j})"
                    )
                direct_checks.append(
                    {
                        "order": order,
                        "degree": degree,
                        "i": i,
                        "j": j,
                    }
                )
            for order in range(max(1, 3 * degree + 1), 3 * degree + 4):
                beta = beta_mixed_contraction(order, degree, i, j)
                if beta != 0:
                    raise AssertionError(
                        "uniform mixed cutoff failed at "
                        f"(m,e,i,j)=({order},{degree},{i},{j})"
                    )
                cutoff_checks.append(
                    {
                        "order": order,
                        "degree": degree,
                        "i": i,
                        "j": j,
                        "phase": phase,
                    }
                )
    return {
        "direct_beta_checks": len(direct_checks),
        "uniform_cutoff_checks": len(cutoff_checks),
        "uniform_balanced_bound": "m > 3e",
        "arbitrary_multiplier_bound": (
            "m > 3 deg_coordinate(Q), componentwise after bihomogeneous "
            "decomposition"
        ),
    }


def main() -> None:
    arguments = parse_arguments()
    if arguments.through < 1:
        raise ValueError("--through must be positive")
    if arguments.mixed_degree < 0:
        raise ValueError("--mixed-degree must be nonnegative")

    expected_factorization = -(
        X * Y + Y**2 + 1
    ) * (
        X**2 * Y + X - Y
    )
    if sp.expand(F - expected_factorization) != 0:
        raise AssertionError("dehomogenized factorization failed")
    if sp.cancel(P - P_FACTORED) != 0:
        raise AssertionError("Hopf/beta factorization failed")
    determinant = int(COEFFICIENT_MATRIX.det())
    if determinant != 1:
        raise AssertionError("unexpected coefficient determinant")

    moments = [dehomogeneous_moment(order) for order in range(1, arguments.through + 1)]
    if any(moments):
        raise AssertionError("a checked pure moment is nonzero")

    rodrigues_checks = verify_rodrigues_identity(arguments.mixed_degree)
    mixed = verify_mixed_formula(arguments.mixed_degree)

    payload = {
        "calculation": "two_pair_sic_bidegree33_rodrigues_survivor",
        "status": "proved",
        "field": "characteristic zero",
        "coefficient_matrix": [
            [int(COEFFICIENT_MATRIX[i, j]) for j in range(4)]
            for i in range(4)
        ],
        "coefficient_determinant": determinant,
        "coefficient_rank": int(COEFFICIENT_MATRIX.rank()),
        "support_size": sum(
            int(COEFFICIENT_MATRIX[i, j] != 0)
            for i in range(4)
            for j in range(4)
        ),
        "dehomogenized_form": str(F),
        "factorization": str(expected_factorization),
        "hopf_beta_factorization": (
            "p(t,u)=u^-3*(u^2+1-t)*((1-t)^2-t*u^2)"
        ),
        "pure_moments_checked": {
            "orders": [1, arguments.through],
            "values": moments,
        },
        "all_order_pure_identity": {
            "odd_order": "phase parity gives zero constant term",
            "even_order": (
                "for m=2n the beta sum is proportional to "
                "sum_{r=0}^n (-1)^r/((n-r)!r!)=0"
            ),
        },
        "rodrigues_identity": {
            "formula": (
                "CT_u[u^d p(t,u)^m] = "
                "(-1)^r*(1-t)^d/s!*D_t^s[t^m(1-t)^m], "
                "r=(m-d)/2, s=(m+d)/2"
            ),
            "symbolic_checks": len(rodrigues_checks),
            "phase_range": [
                -arguments.mixed_degree,
                arguments.mixed_degree,
            ],
        },
        "sic_safety": {
            "proved": True,
            "reason": (
                "after multiplying by the balanced monomial M_ij of "
                "bidegree (e,e), integration by parts s times differentiates "
                "t^i(1-t)^(e-j); it is zero when s>e+i-j.  Hence every "
                "balanced multiplier vanishes for m>3e.  An unbalanced "
                "multiplier is either degree-forced to zero or each output "
                "coefficient is detected by a balanced completion."
            ),
            **mixed,
        },
        "mathematical_consequence": {
            "moment_nullcone": (
                "MN_3 is false because all pure moments vanish while det(C)=1"
            ),
            "sic": (
                "this witness is not an SIC counterexample; it is SIC-safe "
                "with the explicit eventual bound above"
            ),
        },
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "status": payload["status"],
                "determinant": determinant,
                "pure_orders_checked": arguments.through,
                "rodrigues_checks": len(rodrigues_checks),
                "mixed": mixed,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
