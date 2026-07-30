#!/usr/bin/env python3
"""Exact small census for branching the DC_2 classical-symbol search.

This is a pre-search certificate.  It verifies four inequivalent incidence
packets and base-changes the bounded canonical fiber/connection complex to
each packet.  It does not claim a polynomial A_2 quantization or a global
boundary descent for the reciprocal packet.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.deformation_complex import (  # noqa: E402
    relative_fiber_connection_complex,
)
from master_cancellation import (  # noqa: E402
    hensel_jet,
    parameter_polynomial,
)


@dataclass(frozen=True)
class IncidencePacket:
    key: str
    boundary_type: str
    base_ring: str
    geometric_degree: int
    controlled_exponent: int
    normalized_unit_rank: int
    conductor_length: int
    composition_blocks: str | None
    rank_two_status: str

    @property
    def fingerprint(self) -> tuple[object, ...]:
        return (
            self.geometric_degree,
            self.boundary_type,
            self.normalized_unit_rank,
            self.conductor_length,
            self.composition_blocks,
        )


def verify_weighted_packet() -> IncidencePacket:
    """The full normalized degree-six weighted-seed family."""

    w = sp.symbols("w")
    h4, h5, h6 = sp.symbols("h4 h5 h6")
    h3 = -1 - 2 * h4 - 3 * h5 - 4 * h6
    seed = sp.expand(
        h3 * (w**3 - w**2)
        + h4 * (w**4 - w**2)
        + h5 * (w**5 - w**2)
        + h6 * (w**6 - w**2)
    )
    assert seed.subs(w, 0) == 0
    assert sp.diff(seed, w).subs(w, 0) == 0
    assert seed.subs(w, 1) == 0
    assert sp.diff(seed, w).subs(w, 1) == -1
    assert sp.Poly(seed, w).degree() == 6
    return IncidencePacket(
        key="weighted_n6",
        boundary_type="weighted A1 normalization",
        base_ring=(
            "Q[h4,h5,h6,1/(h6*(H''(1)+2))] with "
            "h3=-1-2h4-3h5-4h6"
        ),
        geometric_degree=6,
        controlled_exponent=1,
        normalized_unit_rank=0,
        conductor_length=0,
        composition_blocks=None,
        rank_two_status="polynomial exact-symplectic family proved",
    )


def verify_reciprocal_packet() -> IncidencePacket:
    """The type-(m,r)=(2,1) quartic cancellation incidence."""

    q, A = sp.symbols("q A")
    modulus = sp.expand(parameter_polynomial(2, 1, q))
    assert modulus == q**2 - 4 * q + 6
    jet = hensel_jet(2, 1, A, q)
    jet_expected = q + (4 * q - 6) * A
    assert sp.rem(
        sp.Poly(sp.expand(jet - jet_expected), q),
        sp.Poly(modulus, q),
    ).is_zero

    T, P, Q, R = sp.symbols("T P Q R")
    incidence = (
        T
        - Q**2 * T**2 / 2
        + 2 * P * Q * T**3 / 3
        - P**2 * T**4 / 4
        - R
    )
    assert sp.Poly(incidence, T).degree() == 4
    assert sp.factor(
        sp.diff(incidence, T) - (1 - T * (Q - P * T) ** 2)
    ) == 0
    return IncidencePacket(
        key="reciprocal_m2_r1",
        boundary_type="reciprocal Gm normalization",
        base_ring="Q[q]/(q^2-4q+6)",
        geometric_degree=4,
        controlled_exponent=1,
        normalized_unit_rank=1,
        conductor_length=0,
        composition_blocks=None,
        rank_two_status=(
            "localized canonical candidate; polynomial A4 descent open"
        ),
    )


def verify_conductor_packet() -> IncidencePacket:
    """The conductor-selected foundational cusp incidence."""

    w, s, t = sp.symbols("w s t")
    S, V = sp.symbols("S V")
    seed = w**2 * (1 - w)
    incidence = seed - s * w + t
    discriminant = sp.factor(sp.discriminant(incidence, w))
    shifted = sp.factor(
        discriminant.subs(
            {
                s: sp.Rational(1, 3) + S,
                t: sp.Rational(1, 27) + S / 3 + V,
            },
            simultaneous=True,
        )
    )
    assert shifted == -4 * S**3 - 27 * V**2

    u = sp.symbols("u")
    cusp_relation = sp.expand(
        4 * (-3 * u**2) ** 3 + 27 * (-2 * u**3) ** 2
    )
    assert cusp_relation == 0
    # Q[u]/Q[u^2,u^3] has the single missing monomial u; its conductor is
    # u^2 Q[u].
    conductor_length = 1
    return IncidencePacket(
        key="conductor_cusp_n3",
        boundary_type="cusp conductor Q[u^2,u^3] in Q[u]",
        base_ring="Q (with conductor square at (u^2,u^3))",
        geometric_degree=3,
        controlled_exponent=1,
        normalized_unit_rank=0,
        conductor_length=conductor_length,
        composition_blocks=None,
        rank_two_status="polynomial c=-9 exact-symplectic completion proved",
    )


def verify_composition_packet() -> IncidencePacket:
    """The imprimitive tower obtained from F_6 after F_3."""

    inner_degree = 3
    outer_degree = 6
    total_degree = inner_degree * outer_degree
    assert total_degree == 18
    return IncidencePacket(
        key="composition_f6_after_f3",
        boundary_type="two-stage weighted incidence tower",
        base_ring=(
            "Q[h4,h5,h6,1/(h6*(H''(1)+2))] for the "
            "outer degree-six seed"
        ),
        geometric_degree=total_degree,
        controlled_exponent=1,
        normalized_unit_rank=0,
        conductor_length=0,
        composition_blocks="6 blocks of size 3",
        rank_two_status=(
            "composition of polynomial exact-symplectic rank-two completions"
        ),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        help="optional path for the deterministic JSON certificate",
    )
    args = parser.parse_args()

    packets = (
        verify_weighted_packet(),
        verify_reciprocal_packet(),
        verify_conductor_packet(),
        verify_composition_packet(),
    )
    fingerprints = [packet.fingerprint for packet in packets]
    assert len(set(fingerprints)) == len(fingerprints)

    data = relative_fiber_connection_complex()
    complex_ = data.complex
    assert complex_.dimensions == (35, 60, 26)
    assert complex_.ranks == (34, 26)
    assert complex_.cohomology_dimensions == (1, 0, 0)
    assert complex_.dual_obstruction_cocycles() == ()

    # A nonzero maximal minor is a unit over every Q-algebra in the census.
    # Hence Fitt_0(coker d1)=(1), the obstruction module is zero, and at this
    # first bounded stage the strong dual-cocycle module is also zero.
    _, pivot_columns = complex_.d1.rref()
    assert len(pivot_columns) == complex_.d1.rows
    unit_minor = sp.factor(
        complex_.d1[:, list(pivot_columns)].det()
    )
    assert unit_minor != 0

    prime_profile = complex_.prime_rank_profile((31991, 32003, 65521))
    assert set(prime_profile.values()) == {complex_.ranks}

    certificate = {
        "scope": (
            "bounded canonical degrees 4->3->2; no Ore-boundary descent "
            "or polynomial A_2 quantization claim"
        ),
        "families": [asdict(packet) for packet in packets],
        "relative_complex": {
            "dimensions": list(complex_.dimensions),
            "ranks": list(complex_.ranks),
            "cohomology_dimensions": list(
                complex_.cohomology_dimensions
            ),
            "unit_maximal_minor": str(unit_minor),
            "fitting_0_obstruction_module": "unit ideal",
            "strong_dual_cocycle_module": "zero",
            "prime_rank_profile": {
                str(prime): list(ranks)
                for prime, ranks in prime_profile.items()
            },
        },
    }
    rendered = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if args.output:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)

    print("PASS: four incidence packets have distinct exact fingerprints")
    print("PASS: weighted, reciprocal, conductor, and composition types occur")
    print("PASS: relative complex dimensions are 35 -> 60 -> 26")
    print("PASS: Fitt_0(E)=(1) and the first strong-cocycle module is zero")
    print("PASS: ranks agree over Q and three good primes")
    print(
        "SCOPE: bounded canonical chart only; reciprocal polynomial A4 "
        "descent and all global A2 quantizations remain open"
    )


if __name__ == "__main__":
    main()
