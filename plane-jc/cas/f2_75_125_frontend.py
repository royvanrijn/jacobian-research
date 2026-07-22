#!/usr/bin/env python3
"""Exact source-level skeleton for the F2 j=1 (75,125) front end.

This module records everything that follows mechanically from the published
complete-chain theorem and the type-I endpoint theorem.  It deliberately does
not promote those facts to an exhaustive Laurent normal form: the lower
boundary after the Puiseux translation is not classified in the cited
sources.

Run this file to emit the residual certificate as JSON.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import json

import sympy as sp


Point = tuple[int, int]


@dataclass(frozen=True)
class ForcedEdge:
    name: str
    P_start: Point
    P_end: Point
    Q_start: Point
    Q_end: Point
    weight: Point
    vertex_nonvanishing: tuple[str, ...]


@dataclass(frozen=True)
class ResidualObligation:
    id: str
    statement: str
    needed_for: tuple[str, ...]
    source_status: str


def _fraction(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


def chain_data() -> dict[str, object]:
    """Return the exact arithmetic attached to family F2 at j=1."""

    m, n = 3, 5
    a0 = (5, 20)
    a0_prime = (1, 0)
    a1 = (Fraction(7, 5), 2)
    rho0_sigma0 = (5, -1)
    l1 = 5
    k = 1
    assert (m + n) * a1[1] * k - n * (l1 * a1[1] - a1[0].numerator) == k
    return {
        "family": "F2",
        "parameter_j": 1,
        "degree_pair": [75, 125],
        "multiplicities": [m, n],
        "initial_corner": list(a0),
        "initial_start": list(a0_prime),
        "initial_direction": list(rho0_sigma0),
        "final_corner": [_fraction(a1[0]), a1[1]],
        "puiseux_denominator": l1,
        "final_type_I_k": k,
        "puiseux_translation": "X^5=x; y -> y+lambda/X; lambda != 0",
        "translation_root_multiplicities": {"P": 2 * m, "Q": 2 * n},
        "transformed_bracket": "[P,Q]_(X,y)=X^4 after scalar normalization",
    }


def forced_edges() -> tuple[ForcedEdge, ForcedEdge]:
    """Two consecutive forced edges in the integral X=x^(1/5) lattice."""

    initial = ForcedEdge(
        name="translated_type_II_edge",
        P_start=(75, 60),
        P_end=(21, 6),
        Q_start=(125, 100),
        Q_end=(35, 10),
        weight=(1, -1),
        vertex_nonvanishing=("p_75_60", "p_21_6", "q_125_100", "q_35_10"),
    )
    terminal = ForcedEdge(
        name="final_type_I_edge",
        P_start=(21, 6),
        P_end=(4, 1),
        Q_start=(35, 10),
        Q_end=(1, 0),
        weight=(5, -17),
        vertex_nonvanishing=("p_21_6", "p_4_1", "q_35_10", "q_1_0"),
    )
    for edge in (initial, terminal):
        assert edge.weight[0] * edge.P_start[0] + edge.weight[1] * edge.P_start[1] == (
            edge.weight[0] * edge.P_end[0] + edge.weight[1] * edge.P_end[1]
        )
        assert edge.weight[0] * edge.Q_start[0] + edge.weight[1] * edge.Q_start[1] == (
            edge.weight[0] * edge.Q_end[0] + edge.weight[1] * edge.Q_end[1]
        )
    return initial, terminal


def normalized_terminal_edge() -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Return the uniquely normalized terminal edge and its exact bracket."""

    X, y = sp.symbols("X y")
    s = X**17 * y**5
    P = X**4 * y * (1 + s)
    Q = -X * (1 + 3 * s + sp.Rational(9, 5) * s**2)
    bracket = sp.expand(sp.diff(P, X) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, X))
    assert bracket == X**4
    return P, Q, bracket


def common_power_top_band() -> dict[str, object]:
    """Describe the forced common-power band in t=X*y, z=y^-1."""

    t, z = sp.symbols("t z")
    h = sp.Function("H")
    C = t**7 * h(t) * z**5
    return {
        "band_chart": {"t": "X*y", "z": "y^-1", "jacobian": "[t,z]_(X,y)=-1"},
        "root_band": str(C),
        "H_degree": 18,
        "H_endpoint_nonvanishing": ["coeff(H,t^0) != 0", "coeff(H,t^18) != 0"],
        "P_top_band": "t^21*H(t)^3*z^15",
        "Q_top_band": "t^35*H(t)^5*z^25",
        "top_bracket": "0",
        "rhs_band": "X^4=t^4*z^4",
        "unresolved_layer_gap": 35,
    }


def residual_obligations() -> tuple[ResidualObligation, ...]:
    return (
        ResidualObligation(
            id="F2-NF-1",
            statement=(
                "Prove an F2 j=1 reduction theorem that controls every support point "
                "after y -> y+lambda*x^(-1/5), not only the two regular edges."
            ),
            needed_for=("exhaustive polygon branches", "complete Laurent supports"),
            source_status=(
                "The 2017 complete-chain theorem supplies corners and the translation; "
                "it does not supply the lower Newton boundary."
            ),
        ),
        ResidualObligation(
            id="F2-NF-2",
            statement=(
                "Determine all gamma normalizations for the (m,n)=(3,5) member, "
                "including whether the j=0 values 2 and 3 persist, and prove the "
                "resulting list exhaustive with all vertex nonvanishing."
            ),
            needed_for=("normalizing transformations", "branch exhaustiveness"),
            source_status=(
                "The 2014 Section 5 statement treats only (m,n)=(2,3) and explicitly "
                "omits the proof of its preliminary normalization."
            ),
        ),
        ResidualObligation(
            id="F2-NF-3",
            statement=(
                "For each resulting polygon, enumerate every z-band down from the "
                "common-power layer 39 to the bracket layer 4 and transcribe its "
                "coefficient normalizations."
            ),
            needed_for=("weighted-Wronskian IR", "de Rham obstruction", "residual ideal"),
            source_status=(
                "The forced terminal type-I block has zero residual obstruction; the "
                "first potentially new obstruction depends on the missing lower bands."
            ),
        ),
    )


def machine_certificate() -> dict[str, object]:
    P, Q, bracket = normalized_terminal_edge()
    return {
        "schema": "plane-jc.f2-75-125-residual.v1",
        "status": "partial-source-certificate-not-an-exhaustive-normal-form",
        "chain": chain_data(),
        "forced_edges": [asdict(edge) for edge in forced_edges()],
        "terminal_edge_normalization": {
            "complete_supports": {
                "P": [[4, 1], [21, 6]],
                "Q": [[1, 0], [18, 5], [35, 10]],
            },
            "coefficient_system": {
                "variables": ["a", "b", "c", "d", "e"],
                "nonzero": ["a", "b", "c", "d", "e"],
                "equations": [
                    {"terms": [[-1, ["a", "c"]], [-1, []]], "equals": 0},
                    {"terms": [[2, ["a", "d"]], [-6, ["b", "c"]]], "equals": 0},
                    {"terms": [[5, ["a", "e"]], [-3, ["b", "d"]]], "equals": 0},
                ],
                "torus_gauge": ["a=1", "b=1"],
                "unique_normalized_solution": {"a": "1", "b": "1", "c": "-1", "d": "-3", "e": "-9/5"},
            },
            "P": str(P),
            "Q": str(Q),
            "bracket": str(bracket),
            "coefficient_normalizations": {
                "p_4_1": "1",
                "p_21_6": "1",
                "q_1_0": "-1",
                "q_18_5": "-3",
                "q_35_10": "-9/5",
            },
            "de_rham_obstruction_rank": 0,
        },
        "common_power_top_band": common_power_top_band(),
        "laurent_polygon_branches": {
            "status": "not-derived",
            "known_branch_count": None,
            "reason": "the lower Newton boundary is not fixed by the complete-chain theorem",
        },
        "residual_obligations": [asdict(item) for item in residual_obligations()],
        "frontend_complete": False,
    }


if __name__ == "__main__":
    print(json.dumps(machine_certificate(), indent=2, sort_keys=True))
