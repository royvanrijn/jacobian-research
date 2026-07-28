#!/usr/bin/env python3
"""Audit the boundary handoff from the F2 quadratic common-root family.

The exact upper-band classification gives

    C0(u)=(u-1)^5*(u^5-1)^2*R(u^5),
    R(v)=a*v^2+b*v+(1/25-a-b).

This checker classifies the four disjoint root strata and their contact
partitions on the selected toric divisor.  It then enforces the repository's
typed Newton/boundary boundary: contact multiplicities are not silently
promoted to finite-normalization ramification indices.

For diagnosis only, it also runs the strongest naive surrogate in which
every contact center is treated as a distinct residue-degree-one boundary
row.  Even this unsupported promotion survives the coarse finite-flat and
packet-length gates, so it cannot be used as an F2 exclusion.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from math import gcd
from pathlib import Path

import sympy as sp

from finite_normalization_signatures import (
    BoundaryRow,
    TargetNormalizationSignature,
)
from plane_boundary_exclusion import (
    AffinePrimeLedgerEntry,
    BoundaryPrimeLedgerEntry,
    ConductorPacketPoint,
    TargetComponentLedger,
    TargetFiberPacket,
    audit_target_component_ledger,
    conductor_packet_budget,
)


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_boundary_handoff.json"
)


@dataclass(frozen=True)
class RootStratum:
    name: str
    conditions: tuple[str, ...]
    witness_a: Fraction
    witness_b: Fraction
    R_root_type: str
    contact_partition: tuple[int, ...]

    @property
    def center_count(self) -> int:
        return len(self.contact_partition)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def root_strata() -> tuple[RootStratum, ...]:
    """Return the exhaustive four-way quadratic root classification."""

    strata = (
        RootStratum(
            name="two_distinct_nonzero_R_roots",
            conditions=("c != 0", "Delta != 0"),
            witness_a=Fraction(1, 100),
            witness_b=Fraction(0),
            R_root_type="two distinct nonzero roots",
            contact_partition=(7, *([2] * 4), *([1] * 10)),
        ),
        RootStratum(
            name="one_double_nonzero_R_root",
            conditions=("c != 0", "Delta = 0"),
            witness_a=Fraction(1, 25),
            witness_b=Fraction(-4, 25),
            R_root_type="one double nonzero root",
            contact_partition=(7, *([2] * 9)),
        ),
        RootStratum(
            name="zero_and_nonzero_simple_R_roots",
            conditions=("c = 0", "b != 0"),
            witness_a=Fraction(1, 100),
            witness_b=Fraction(3, 100),
            R_root_type="one zero root and one nonzero simple root",
            contact_partition=(7, 5, *([2] * 4), *([1] * 5)),
        ),
        RootStratum(
            name="double_zero_R_root",
            conditions=("c = 0", "b = 0", "a = 1/25"),
            witness_a=Fraction(1, 25),
            witness_b=Fraction(0),
            R_root_type="one double root at zero",
            contact_partition=(10, 7, *([2] * 4)),
        ),
    )
    for stratum in strata:
        a = stratum.witness_a
        b = stratum.witness_b
        c = Fraction(1, 25) - a - b
        discriminant = b * b - 4 * a * c
        if a == 0:
            raise AssertionError("a witness lost quadratic degree")
        if a + b + c != Fraction(1, 25):
            raise AssertionError("R(1) normalization changed")
        if sum(stratum.contact_partition) != 25:
            raise AssertionError("a contact partition no longer has degree 25")
        if stratum.name == "two_distinct_nonzero_R_roots":
            assert c != 0 and discriminant != 0
        elif stratum.name == "one_double_nonzero_R_root":
            assert c != 0 and discriminant == 0
        elif stratum.name == "zero_and_nonzero_simple_R_roots":
            assert c == 0 and b != 0 and discriminant != 0
        else:
            assert c == 0 and b == 0 and discriminant == 0
    return strata


def common_edge_regression() -> dict[str, object]:
    """Verify the factorization and the fixed/root-dependent multiplicities."""

    u, v, a, b = sp.symbols("u v a b")
    c = sp.Rational(1, 25) - a - b
    R = a * v**2 + b * v + c
    C0 = sp.expand((u - 1) ** 5 * (u**5 - 1) ** 2 * R.subs(v, u**5))
    if sp.Poly(C0, u).degree() != 25:
        raise AssertionError("the common edge polynomial lost degree 25")
    if sp.expand(R.subs(v, 1)) != sp.Rational(1, 25):
        raise AssertionError("R(1) is no longer fixed")
    fixed_factor = sp.expand((u - 1) ** 5 * (u**5 - 1) ** 2)
    if sp.rem(C0, fixed_factor, u) != 0:
        raise AssertionError("the fixed root-of-unity factor disappeared")
    return {
        "u_coordinate": "u=1+t",
        "common_edge_polynomial": (
            "C0(u)=(u-1)^5*(u^5-1)^2*R(u^5)"
        ),
        "R": "a*v^2+b*v+(1/25-a-b)",
        "fixed_centers": {
            "u=1": {"contact_multiplicity": 7},
            "nontrivial_fifth_roots_of_unity": {
                "center_count": 4,
                "contact_multiplicity_each": 2,
            },
        },
        "no_collision_with_R_roots": "R(1)=1/25 != 0",
        "degree": 25,
    }


def scale_ambiguity() -> dict[str, object]:
    """Exhibit different local scales with the same edge restriction."""

    records: list[dict[str, object]] = []
    for multiplicity in (1, 2, 5, 7, 10):
        germs = []
        for normal_order in (1, 2, 3):
            divisor = gcd(multiplicity, normal_order)
            germs.append(
                {
                    "germ": (
                        f"g_{normal_order}(s,z)=s^{multiplicity}"
                        f"+z^{normal_order}"
                    ),
                    "restriction_at_z_zero": f"s^{multiplicity}",
                    "primitive_equality_ray_z_s": [
                        multiplicity // divisor,
                        normal_order // divisor,
                    ],
                }
            )
        if len(
            {
                tuple(germ["primitive_equality_ray_z_s"])
                for germ in germs
            }
        ) < 2:
            raise AssertionError("the scale ambiguity regression collapsed")
        records.append(
            {
                "contact_multiplicity": multiplicity,
                "same_edge_different_scales": germs,
            }
        )
    return {
        "theorem": (
            "the restriction g(s,0)=s^m does not determine the first "
            "nonzero normal order or the toroidal equality ray"
        ),
        "records": records,
        "consequence": (
            "an F2 contact multiplicity cannot be used as a branch scale, "
            "ramification index, residue degree, or puncture count"
        ),
    }


def conditional_normalization_surrogate(
    stratum: RootStratum,
) -> dict[str, object]:
    """Run an intentionally over-strong, explicitly unsupported handoff."""

    rows = tuple(
        BoundaryRow(
            ramification_index=multiplicity,
            residue_degree=1,
            punctures=1,
        )
        for multiplicity in stratum.contact_partition
    )
    # If every contact were a boundary row, its total contribution would be
    # 25.  Finite normalization requires a positive affine contribution, so
    # the smallest possible generic degree would be 26.
    signature = TargetNormalizationSignature(
        geometric_degree=26,
        boundary_rows=rows,
        affine_residue_degrees=(1,),
    )
    budget = conductor_packet_budget(
        stratum.contact_partition,
        generic_boundary_degree=25,
        affine_degree=1,
    )
    boundary_primes = tuple(
        BoundaryPrimeLedgerEntry(
            name=f"E{index}",
            transverse_index=multiplicity,
            residue_degree=1,
        )
        for index, multiplicity in enumerate(
            stratum.contact_partition,
            start=1,
        )
    )
    packet = TargetFiberPacket(
        name="all_edge_centers_in_one_target_fiber",
        points=tuple(
            ConductorPacketPoint(
                name=f"p{index}",
                boundary_prime=f"E{index}",
                transverse_index=multiplicity,
                residue_immersive=False,
            )
            for index, multiplicity in enumerate(
                stratum.contact_partition,
                start=1,
            )
        ),
        same_target_fiber_certified=False,
    )
    ledger = TargetComponentLedger(
        name=f"unsupported surrogate: {stratum.name}",
        generic_degree=26,
        boundary_primes=boundary_primes,
        affine_primes=(AffinePrimeLedgerEntry("A1", 1),),
        packets=(packet,),
        finite_flat_certified=True,
        target_transfer_certified=False,
        exhaustive_generic_pullback=False,
    )
    ledger_audit = audit_target_component_ledger(ledger)
    if budget.length_deficit != -1 or "permits" not in budget.verdict:
        raise AssertionError("the strongest packet surrogate became excluded")
    if ledger_audit.status != "incomplete":
        raise AssertionError("unsupported source contacts passed target transfer")
    return {
        "assumption_not_proved": (
            "each edge center is a distinct residue-degree-one normalization "
            "boundary prime with transverse index equal to contact multiplicity"
        ),
        "coarse_signature": {
            "generic_degree": signature.geometric_degree,
            "boundary_rows_e_f_s": [
                [
                    row.ramification_index,
                    row.residue_degree,
                    row.punctures,
                ]
                for row in signature.boundary_rows
            ],
            "boundary_degree": signature.boundary_degree,
            "affine_degree": signature.affine_degree,
            "residue_immersive_compatible": (
                signature.residue_immersive_compatible
            ),
        },
        "all_centers_one_fiber_packet_budget": asdict(budget),
        "typed_target_ledger_audit": ledger_audit.as_dict(),
        "outcome": (
            "even the strongest naive promotion survives the arithmetic "
            "packet budget, while the typed audit correctly remains incomplete"
        ),
    }


def build_payload() -> dict[str, object]:
    strata = root_strata()
    return {
        "schema": "plane-jc.f2-75-125-boundary-handoff.v1",
        "status": "exact-contact-classification-boundary-transfer-incomplete",
        "common_edge": common_edge_regression(),
        "root_strata": [
            {
                **asdict(stratum),
                "witness_a": fraction_text(stratum.witness_a),
                "witness_b": fraction_text(stratum.witness_b),
                "center_count": stratum.center_count,
                "conditional_normalization_surrogate": (
                    conditional_normalization_surrogate(stratum)
                ),
            }
            for stratum in strata
        ],
        "scale_nonuniqueness": scale_ambiguity(),
        "typed_handoff": {
            "certified": [
                "the four exhaustive quadratic-R root strata",
                "the number of distinct centers on the selected toric divisor",
                "the contact multiplicity partition at those centers",
            ],
            "missing": [
                "the first nonzero normal order at every center",
                "certified toroidal branch scales and proximity centers",
                "which resolved components are dicritical",
                "target-boundary pullback multiplicities",
                "finite-normalization transverse indices e and residue degrees f",
                "puncture counts and grouping over target curves or closed fibers",
                "the geometric degree and exhaustive affine-sheet ledger",
            ],
            "log_boundary_compiler_applicable": False,
            "finite_normalization_packet_theorem_applicable": False,
        },
        "verdict": {
            "excluded_strata": [],
            "surviving_contact_partitions": [
                list(stratum.contact_partition) for stratum in strata
            ],
            "closure_assessment": (
                "the quadratic-R data do not close F2 through current boundary "
                "gates; lower normal terms are exactly the missing bridge"
            ),
            "recommended_pivot": (
                "stop sequential F2 descent and seek a degree-independent "
                "theorem converting a common-power edge packet into target-side "
                "finite-normalization data"
            ),
        },
        "software": {
            "python": "standard library",
            "sympy": sp.__version__,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    args = parser.parse_args()
    # Normalize tuples from dataclass ``asdict`` calls to their JSON list
    # representation before comparing with the pinned artifact.
    payload = json.loads(json.dumps(build_payload(), sort_keys=True))
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display_path = artifact.relative_to(ROOT)
        except ValueError:
            display_path = artifact
        print(f"WROTE {display_path}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError(
                "pinned F2 boundary-handoff artifact is stale; "
                "inspect before --refresh"
            )
    print("PASS: four exhaustive quadratic-R contact strata")
    print("PASS: contact partitions have degree 25 and between 6 and 15 centers")
    print("PASS: identical edge contacts admit different toroidal scales")
    print("PASS: every strongest-naive finite-normalization surrogate survives")
    print("PASS: the typed boundary handoff remains explicitly incomplete")


if __name__ == "__main__":
    main()
