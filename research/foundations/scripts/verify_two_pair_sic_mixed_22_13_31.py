#!/usr/bin/env python3
"""Exact SIC(2) theorem on the mixed (2,2)+(1,3)+(3,1) stratum.

Write

    F = A + B + C,

where A, B, C have bidegrees (2,2), (1,3), and (3,1).  If B is nonzero,
the dual-linear theorem and a contraction-preserving GL_2 change put

    B = xi2*z1^3.

This checker verifies the remaining exact calculation.  The full
polynomial-valued contractions through order four force the only
contraction-central coefficients to vanish.  Every surviving monomial
then lies in a two-step contraction cone, which gives eventual vanishing
for every fixed multiplier.

The B=0 branch is discharged in the written proof by the complete
bidegree-(2,2) theorem.  It is not re-proved here.
"""

from __future__ import annotations

import json
from math import comb, factorial, gcd
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_mixed_22_13_31.json"
)

XI1, XI2, Z1, Z2 = sp.symbols("xi1 xi2 z1 z2")
VARIABLES = (XI1, XI2, Z1, Z2)


def contraction(expression: sp.Expr) -> sp.Expr:
    """Apply E_2 exactly and return a polynomial in z1,z2."""

    result = sp.Integer(0)
    for (xi1, xi2, z1, z2), coefficient in sp.Poly(
        sp.expand(expression), *VARIABLES
    ).terms():
        if xi1 > z1 or xi2 > z2:
            continue
        result += (
            coefficient
            * sp.Rational(factorial(z1), factorial(z1 - xi1))
            * sp.Rational(factorial(z2), factorial(z2 - xi2))
            * Z1 ** (z1 - xi1)
            * Z2 ** (z2 - xi2)
        )
    return sp.expand(result)


def coefficient(expression: sp.Expr, z1: int, z2: int) -> sp.Expr:
    return sp.Poly(sp.expand(expression), Z1, Z2).coeff_monomial(
        Z1**z1 * Z2**z2
    )


def bidegree_ledger(degree_cap: int = 7) -> list[dict[str, int | str]]:
    ledger: list[dict[str, int | str]] = []
    for ordinary_degree in range(1, degree_cap + 1):
        for dual_degree in range(ordinary_degree + 1):
            coordinate_degree = ordinary_degree - dual_degree
            weight = coordinate_degree - dual_degree
            if weight > 0:
                sector = "positive"
            elif weight < 0:
                sector = "negative"
            else:
                sector = "balanced"
            ledger.append(
                {
                    "dual_degree": dual_degree,
                    "coordinate_degree": coordinate_degree,
                    "ordinary_degree": ordinary_degree,
                    "central_weight": weight,
                    "sector": sector,
                }
            )
    return ledger


def monomial_cone_data(
    exponent: tuple[int, int, int, int],
) -> dict[str, int | list[int]]:
    xi1, xi2, z1, z2 = exponent
    weight = z1 + z2 - xi1 - xi2
    delta = z1 - xi1
    epsilon = xi2 - z2
    assert epsilon == delta - weight
    return {
        "exponent": list(exponent),
        "weight": weight,
        "delta": delta,
        "epsilon": epsilon,
    }


def primitive_two_block_circuits(
    ledger: list[dict[str, int | str]],
) -> list[dict[str, object]]:
    positive = [item for item in ledger if item["sector"] == "positive"]
    negative = [item for item in ledger if item["sector"] == "negative"]
    circuits: list[dict[str, object]] = []
    for positive_block in positive:
        positive_weight = int(positive_block["central_weight"])
        for negative_block in negative:
            negative_weight = -int(negative_block["central_weight"])
            divisor = gcd(positive_weight, negative_weight)
            positive_multiplicity = negative_weight // divisor
            negative_multiplicity = positive_weight // divisor
            central_dual_degree = (
                positive_multiplicity
                * int(positive_block["dual_degree"])
                + negative_multiplicity
                * int(negative_block["dual_degree"])
            )
            central_coordinate_degree = (
                positive_multiplicity
                * int(positive_block["coordinate_degree"])
                + negative_multiplicity
                * int(negative_block["coordinate_degree"])
            )
            circuits.append(
                {
                    "positive_bidegree": [
                        positive_block["dual_degree"],
                        positive_block["coordinate_degree"],
                    ],
                    "negative_bidegree": [
                        negative_block["dual_degree"],
                        negative_block["coordinate_degree"],
                    ],
                    "positive_multiplicity": positive_multiplicity,
                    "negative_multiplicity": negative_multiplicity,
                    "central_bidegree": [
                        central_dual_degree,
                        central_coordinate_degree,
                    ],
                }
            )
    return circuits


def main() -> None:
    a = sp.symbols("a0:9")
    c = sp.symbols("c0:8")

    A = sp.Integer(0)
    for i in range(3):
        for j in range(3):
            A += (
                a[3 * i + j]
                * XI1**i
                * XI2 ** (2 - i)
                * Z1**j
                * Z2 ** (2 - j)
            )

    B = XI2 * Z1**3

    C = sp.Integer(0)
    for i in range(4):
        for j in range(2):
            C += (
                c[2 * i + j]
                * XI1**i
                * XI2 ** (3 - i)
                * Z1**j
                * Z2 ** (1 - j)
            )

    highest_weight = {
        "E(A)": contraction(A),
        "E(A*B)": contraction(A * B),
        "E(A*B^2)": contraction(A * B**2),
    }
    assert highest_weight["E(A)"] == 2 * a[0] + a[4] + 2 * a[8]
    assert sp.expand(
        highest_weight["E(A*B)"]
        - (6 * a[3] + 12 * a[7]) * Z1**2
        - 12 * a[6] * Z1 * Z2
    ) == 0
    assert highest_weight["E(A*B^2)"] == 60 * a[6] * Z1**4

    first_substitution = {
        a[6]: 0,
        a[3]: -2 * a[7],
        a[4]: -2 * a[0] - 2 * a[8],
    }
    first_reduced = sp.expand((A + B + C).subs(first_substitution))
    third = contraction(first_reduced**3)
    assert sp.expand(
        coefficient(third, 2, 0) + 144 * a[7] * (a[0] - 7 * a[8])
    ) == 0
    assert coefficient(third, 1, 1) == 720 * a[7] ** 2

    triangular_substitution = dict(first_substitution)
    triangular_substitution[a[7]] = 0
    triangular_A = sp.expand(A.subs(triangular_substitution))
    triangular_F = sp.expand(triangular_A + B + C)

    p, q, h = sp.symbols("p q h")
    core = (
        p * XI2**2 * Z2**2
        - 2 * (p + q) * XI1 * XI2 * Z1 * Z2
        + q * XI1**2 * Z1**2
        + B
        + h * XI1**3 * Z2
    )
    core_substitution = {p: a[0], q: a[8], h: c[6]}

    full_moments: dict[int, sp.Expr] = {}
    core_moments: dict[int, sp.Expr] = {}
    for order in range(1, 5):
        full_value = contraction(triangular_F**order)
        core_value = contraction(core**order)
        assert full_value == sp.expand(core_value.subs(core_substitution))
        assert not (full_value.free_symbols & {Z1, Z2})
        full_moments[order] = full_value
        core_moments[order] = core_value

    f2 = 3 * h + 4 * p**2 - 2 * p * q + 4 * q**2
    f3 = (
        -5 * h * p
        + 2 * h * q
        + 4 * p**3
        - 2 * p**2 * q
        - 2 * p * q**2
        + 4 * q**3
    )
    f4 = (
        5 * h**2
        + 9 * h * p**2
        - 6 * h * p * q
        + 15 * h * q**2
        + 12 * p**4
        - 6 * p**3 * q
        + 6 * p**2 * q**2
        - 6 * p * q**3
        + 12 * q**4
    )
    assert core_moments[1] == 0
    assert core_moments[2] == 4 * f2
    assert core_moments[3] == 72 * f3
    assert core_moments[4] == 1728 * f4

    basis = sp.groebner((f2, f3, f4), p, q, h, order="lex")
    expected_basis = (
        3 * h + 4 * p**2 - 2 * p * q + 4 * q**2,
        4 * h * p - h * q + 3 * p * q**2 - 2 * q**3,
        -29 * h**2 + 60 * h * p * q - 33 * h * q**2,
        1027 * h**2 * p + 176 * h**2 * q,
        79 * h**2 + 198 * h * q**2 + 120 * q**4,
        691 * h**2 * q + 1027 * h * q**3,
        h**2 * q**2,
        h**3,
    )
    assert tuple(basis.polys) == tuple(
        sp.Poly(item, p, q, h) for item in expected_basis
    )

    boundary_basis = sp.groebner((f2, f3, h), p, q, h, order="lex")
    assert boundary_basis.reduce(q**4)[1] == 0
    assert boundary_basis.reduce(p**4)[1] == 0
    assert basis.reduce(h**3)[1] == 0

    # After p=q=h=0, these are exactly the possible surviving monomials.
    residual_exponents = {
        "B": (0, 1, 3, 0),
        "a01": (0, 2, 1, 1),
        "a02": (0, 2, 2, 0),
        "a12": (1, 1, 2, 0),
        "c00": (0, 3, 0, 1),
        "c01": (0, 3, 1, 0),
        "c10": (1, 2, 0, 1),
        "c11": (1, 2, 1, 0),
        "c20": (2, 1, 0, 1),
        "c21": (2, 1, 1, 0),
        "c31": (3, 0, 1, 0),
    }
    cone = {
        name: monomial_cone_data(exponent)
        for name, exponent in residual_exponents.items()
    }
    for item in cone.values():
        epsilon = int(item["epsilon"])
        delta = int(item["delta"])
        assert epsilon >= 0
        if epsilon == 0:
            assert delta == -2
        else:
            assert delta <= 3

    ledger = bidegree_ledger()
    assert len(ledger) == 35
    assert sum(item["sector"] == "balanced" for item in ledger) == 3
    assert sum(item["sector"] == "positive" for item in ledger) == 16
    assert sum(item["sector"] == "negative" for item in ledger) == 16
    circuits = primitive_two_block_circuits(ledger)
    assert len(circuits) == 256
    for circuit in circuits:
        assert circuit["central_bidegree"][0] == circuit["central_bidegree"][1]

    opposite_pair_checks: list[dict[str, object]] = []
    h_opposite = sp.symbols("h_opposite")
    for coordinate_degree in range(2, 7):
        positive = XI2 * Z1**coordinate_degree
        corner = XI1**coordinate_degree * Z2
        pair = positive + h_opposite * corner
        checked_moments: dict[str, str] = {}
        for order in range(1, 4):
            odd_moment = contraction(pair ** (2 * order - 1))
            even_moment = contraction(pair ** (2 * order))
            expected = (
                comb(2 * order, order)
                * h_opposite**order
                * factorial(coordinate_degree * order)
                * factorial(order)
            )
            assert odd_moment == 0
            assert even_moment == expected
            checked_moments[str(2 * order)] = str(expected)

        residual_corner_data = []
        for i in range(coordinate_degree + 1):
            for j in range(2):
                if (i, j) == (coordinate_degree, 0):
                    continue
                exponent = (
                    i,
                    coordinate_degree - i,
                    j,
                    1 - j,
                )
                residual_corner_data.append(monomial_cone_data(exponent))
        for item in residual_corner_data:
            epsilon = int(item["epsilon"])
            delta = int(item["delta"])
            assert epsilon >= 0
            if epsilon == 0:
                assert delta == -(coordinate_degree - 1)
            else:
                assert delta <= coordinate_degree

        opposite_pair_checks.append(
            {
                "bidegrees": [
                    [1, coordinate_degree],
                    [coordinate_degree, 1],
                ],
                "ordinary_degree": coordinate_degree + 1,
                "checked_central_moments": checked_moments,
                "residual_multiplier_cutoff_rule": (
                    "m > (delta_Q+(2*d-1)*R_Q)/(d-1), "
                    "R_Q=max(0,-epsilon_Q)"
                ),
            }
        )

    payload = {
        "claim": (
            "SIC(2) holds for every F=A_(2,2)+B_(1,3)+C_(3,1); "
            "when B is nonzero, the first four pure contractions suffice"
        ),
        "degree_cap_ledger": {
            "ordinary_degree_cap": 7,
            "bidegree_count_excluding_constant": len(ledger),
            "balanced_count": 3,
            "positive_weight_count": 16,
            "negative_weight_count": 16,
            "entries": ledger,
            "balanced_nonempty_collection_count": 7,
            "two_sided_collection_count": (
                (2**16 - 1) * (2**16 - 1) * 2**3
            ),
            "primitive_two_block_circuit_count": len(circuits),
            "primitive_two_block_circuits": circuits,
        },
        "normalized_positive_piece": "B=xi2*z1^3",
        "opposite_dual_linear_pair_theorem": {
            "scope": (
                "V_(1,d)+V_(d,1) is SIC-safe for every d>=2; "
                "the requested ordinary degrees 3 through 7 are replayed"
            ),
            "degree_cap_checks": opposite_pair_checks,
        },
        "highest_weight_contractions": {
            key: str(sp.factor(value))
            for key, value in highest_weight.items()
        },
        "third_contraction_positive_coefficients": {
            "z1^2": str(sp.factor(coefficient(third, 2, 0))),
            "z1*z2": str(sp.factor(coefficient(third, 1, 1))),
        },
        "normalized_core_moments": {
            "mu2_over_4": str(f2),
            "mu3_over_72": str(f3),
            "mu4_over_1728": str(f4),
        },
        "core_groebner_basis": [str(item.as_expr()) for item in basis.polys],
        "radical_certificate": {
            "h_power": 3,
            "q_power_after_h": 4,
            "p_power_after_h": 4,
            "conclusion": "radical(f2,f3,f4)=(p,q,h)",
        },
        "residual_cone": cone,
        "multiplier_cutoff_rule": (
            "For a multiplier monomial Q with epsilon_Q=xi2-z2 and "
            "delta_Q=z1-xi1, put R=max(0,-epsilon_Q).  No monomial of "
            "E_2(Q*F^m) survives when m>(delta_Q+5*R)/2."
        ),
        "status": "exact characteristic-zero theorem",
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    print(
        "PASS mixed (2,2)+(1,3)+(3,1): first-four contraction reduction, "
        "origin radical, and eventual multiplier cone"
    )


if __name__ == "__main__":
    main()
