#!/usr/bin/env sage-python
"""Certify the Galois/component height gate for the five product survivors.

status: EXACT_CERTIFICATE
claim: rational height-eight sections on the five surviving product twists
       are necessarily disjoint P.O=0 sections in nonzero Tate class
inputs: direct alternate-Q80 model, full native bisections, two-prime product
        classification, and exact rank-one V4 base generators
outputs: artifacts/generated-results/
         elkies-k3-r17-product-survivor-galois-height-gate-v1.json

At a good branch place of a quadratic twist, the three nonidentity I0*
component classes are permuted as the three nonzero two-torsion points of the
original smooth fibre.  We factor that two-division cubic over each quadratic
residue field.  Irreducibility leaves no rational nonidentity component, so
the Shioda height formula has no local correction for a QQ(u)-section.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from hashlib import sha256
import json
from pathlib import Path
import shlex
import sys

from sage.all import EllipticCurve, NumberField, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
DIRECT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
BISECTIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json"
)
CLASSIFICATION = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-all17-product-toric-frobenius-campaign-v1.json"
)
BASES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-rank-one-bases-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-product-survivor-galois-height-gate-v1.json"
)
SCHEMA = "elkies-k3.r17-product-survivor-galois-height-gate.v1"


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def projective_bits(value) -> int:
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def evaluate(coefficients, value):
    answer = value.parent()(0) if hasattr(value, "parent") else QQ(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + value.parent()(QQ(coefficient))
    return answer


def multiple_u_bits(base, bound: int = 4):
    """Return exact u-height samples, using the stored exceptional n=1 inverse."""

    paired = base["paired_v4_base"]
    curve = EllipticCurve(
        QQ, [QQ(value) for value in paired["pointed_a1_a2_a3_a4_a6"]]
    )
    generator = curve(
        QQ(paired["primitive_generator_pointed"][0]),
        QQ(paired["primitive_generator_pointed"][1]),
    )
    data = base["map_to_v4_cover"]
    origin = data["origin"]
    u0 = QQ(origin["u0"])
    s0 = QQ(origin["s0"])
    q_left = [QQ(value) for value in data["left_q_coefficients_low_to_high"]]
    q_right = [QQ(value) for value in data["right_q_coefficients_low_to_high"]]
    c = QQ(data["pointed_inverse_constants"]["c"])
    d = QQ(data["pointed_inverse_constants"]["d"])
    v0 = QQ(data["pointed_inverse_constants"]["v0"])

    samples = [
        {
            "n": 1,
            "u": str(QQ(data["primitive_generator_image"]["u"])),
            "u_projective_bits": projective_bits(
                data["primitive_generator_image"]["u"]
            ),
            "source": "stored exact exceptional-safe pointed-quartic inverse",
        }
    ]
    for multiple in range(2, bound + 1):
        point = multiple * generator
        if point.is_zero() or point[1] == 0 or v0 == 0:
            raise ArithmeticError("unexpected exceptional multiple in height sample")
        slope = (4 * v0**2 * (point[0] + c) - d**2) / (2 * v0 * point[1])
        denominator = 1 - q_left[2] * slope**2
        if slope == 0 or denominator == 0:
            raise ArithmeticError("height sample meets a map exceptional locus")
        u_value = u0 + (
            (q_left[1] + 2 * q_left[2] * u0) * slope**2 - 2 * s0 * slope
        ) / denominator
        s_value = s0 + (u_value - u0) / slope
        paired_ordinate = (point[0] * slope**2 - d * slope) / (2 * v0) - v0
        t_value = paired_ordinate / denominator
        if s_value**2 != evaluate(q_left, u_value):
            raise ArithmeticError("sample misses the left quadratic cover")
        if t_value**2 != evaluate(q_right, u_value):
            raise ArithmeticError("sample misses the right quadratic cover")
        samples.append(
            {
                "n": multiple,
                "u": str(u_value),
                "u_projective_bits": projective_bits(u_value),
                "source": "exact multiplication and pointed-quartic inverse",
            }
        )
    return samples


def build_payload():
    direct = json.loads(DIRECT.read_text())
    bisections = json.loads(BISECTIONS.read_text())
    classification = json.loads(CLASSIFICATION.read_text())
    bases = json.loads(BASES.read_text())
    if direct.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ValueError("unexpected direct-model status")
    if bisections.get("schema") != "elkies-k3.bisection-extension-input.v1":
        raise ValueError("unexpected bisection schema")
    if classification.get("schema") != "elkies-k3.r17-all17-product-toric-frobenius-campaign.v1":
        raise ValueError("unexpected product-classification schema")
    if bases.get("schema") != "elkies-k3.r17-norm12-11952-v4-rank-one-bases.v1":
        raise ValueError("unexpected rank-one-base schema")

    survivors = classification["persistent_tate_survivor_pair_keys"]
    if len(survivors) != 5:
        raise ArithmeticError("the persistent product worklist no longer has five rows")
    by_label = {row["label"]: row for row in bisections["bisections"]}
    by_pair = {row["pair_key"]: row for row in bases["bases"]}
    classification_by_pair = {
        row["pair_key"]: row for row in classification["targets"]
    }
    ring_u = PolynomialRing(QQ, "u")
    u = ring_u.gen()
    model = direct["weierstrass_model"]
    coefficient_a = ring_u(
        [QQ(value) for value in model["A_coefficients_low_to_high"]]
    )
    coefficient_b = ring_u(
        [QQ(value) for value in model["B_coefficients_low_to_high"]]
    )

    records = []
    p137_certificates = []
    for pair_key in survivors:
        base = by_pair[pair_key]
        classification_row = classification_by_pair[pair_key]
        reduction_137 = next(
            row
            for row in classification_row["reductions"]
            if int(row["prime"]) == 137
        )
        certificate_137 = ROOT / reduction_137["certificate"]
        if digest(certificate_137) != reduction_137["certificate_sha256"]:
            raise ArithmeticError(f"{pair_key}: stale p=137 certificate hash")
        frobenius_137 = json.loads(certificate_137.read_text())
        if (
            frobenius_137.get("pair_key") != pair_key
            or int(frobenius_137.get("prime", 0)) != 137
            or frobenius_137.get("status") != "PASS_COMPLETE_FROBENIUS_PICARD_BOUND"
        ):
            raise ArithmeticError(f"{pair_key}: invalid p=137 certificate")
        cyclotomic_hits = frobenius_137["elliptic_L"][
            "cyclotomic_hits_after_T_equals_pZ"
        ]
        fixed_multiplicity = sum(
            int(hit["total_degree"])
            for hit in cyclotomic_hits
            if int(hit["order"]) == 1
        )
        anti_fixed_multiplicity = sum(
            int(hit["total_degree"])
            for hit in cyclotomic_hits
            if int(hit["order"]) == 2
        )
        if fixed_multiplicity != 1 or anti_fixed_multiplicity != 1:
            raise ArithmeticError(
                f"{pair_key}: expected normalized Tate factor (Z-1)(Z+1) at p=137"
            )
        p137_certificates.append(certificate_137)
        branch_records = []
        for label in pair_key.split(":"):
            branch = by_label[label]["branch"]
            quadratic = ring_u(
                [QQ(value) for value in branch["numerator_coefficients"]]
            )
            if quadratic.degree() != 2 or not quadratic.is_irreducible():
                raise ArithmeticError(f"{label}: branch place is not quadratic")
            short_label = label.rsplit("-", 1)[-1]
            residue_field = NumberField(quadratic, names=f"a_{short_label}")
            branch_root = residue_field.gen()
            ring_x = PolynomialRing(residue_field, f"x_{short_label}")
            x = ring_x.gen()
            two_division = (
                x**3
                + residue_field(coefficient_a(branch_root)) * x
                + residue_field(coefficient_b(branch_root))
            )
            factor_degrees = sorted(
                int(factor.degree())
                for factor, multiplicity in two_division.factor()
                for _unused in range(int(multiplicity))
            )
            if factor_degrees != [3] or not two_division.is_squarefree():
                raise ArithmeticError(
                    f"{label}: the two-division cubic is not irreducible etale"
                )
            branch_records.append(
                {
                    "label": label,
                    "branch_quadratic_coefficients_low_to_high": [
                        str(value) for value in quadratic
                    ],
                    "branch_quadratic_discriminant": str(quadratic.discriminant()),
                    "two_division_factor_degrees_over_residue_field": factor_degrees,
                    "two_division_squarefree": True,
                    "nonidentity_component_fixed_by_residue_galois": False,
                }
            )
        generator_height = base["paired_v4_base"][
            "primitive_generator_height_approx"
        ]
        records.append(
            {
                "shortlist_rank": int(base["shortlist_rank"]),
                "pair_key": pair_key,
                "branch_places": branch_records,
                "arithmetic_rank_gate": {
                    "reduction_prime": 137,
                    "normalized_tate_factor": "(Z-1)*(Z+1)",
                    "frobenius_fixed_multiplicity": fixed_multiplicity,
                    "frobenius_anti_fixed_multiplicity": anti_fixed_multiplicity,
                    "geometric_mw_rank_upper_bound": 2,
                    "QQ_u_mw_rank_upper_bound": 1,
                    "derivation": (
                        "A QQ(u)-rational section has a Frobenius-fixed divisor "
                        "class after good reduction. Injectivity of specialization "
                        "therefore bounds its rank by the normalized Z=1 "
                        "multiplicity in the elliptic quotient."
                    ),
                },
                "constructor_metrics": {
                    "base_jacobian_rank": 1,
                    "primitive_generator_canonical_height_approx": generator_height,
                    "map_degree_to_u": 4,
                    "primitive_generator_u_projective_bits": projective_bits(
                        base["map_to_v4_cover"]["primitive_generator_image"]["u"]
                    ),
                    "u_height_growth_samples": multiple_u_bits(base),
                    "qualification": (
                        "The bit counts are exact samples. The canonical height is the "
                        "stored Sage approximation; map degree four controls quadratic "
                        "asymptotic growth, while finite-n constants remain material."
                    ),
                },
            }
        )

    seed_order = sorted(
        records,
        key=lambda row: (
            row["constructor_metrics"]["primitive_generator_u_projective_bits"],
            row["pair_key"],
        ),
    )
    asymptotic_order = sorted(
        records,
        key=lambda row: (
            Decimal(
                row["constructor_metrics"][
                    "primitive_generator_canonical_height_approx"
                ]
            ),
            row["pair_key"],
        ),
    )
    inputs = (
        Path(__file__).resolve(),
        DIRECT,
        BISECTIONS,
        CLASSIFICATION,
        BASES,
        *p137_certificates,
    )
    return {
        "schema": SCHEMA,
        "status": "PASS_EXACT_PRODUCT_SURVIVOR_GALOIS_HEIGHT_GATE",
        "survivor_count": len(records),
        "records": records,
        "component_galois_gate": {
            "branch_geometry": (
                "At each of the four good twist branch fibres, the three "
                "nonidentity I0* component classes are indexed by the three "
                "nonzero points of the original fibre's two-torsion."
            ),
            "exact_input": (
                "The two-division cubic is irreducible over every one of the ten "
                "quadratic residue fields occurring in the five product targets."
            ),
            "conclusion": (
                "A QQ(u)-rational section meets the identity component at all four "
                "I0* fibres, so every local height correction is zero."
            ),
        },
        "arithmetic_rank_gate": {
            "reduction_prime": 137,
            "normalized_tate_factor_on_all_five": "(Z-1)*(Z+1)",
            "QQ_u_mw_rank_interval_on_each_survivor": [0, 1],
            "geometric_mw_rank_interval_on_each_survivor": [0, 2],
            "construction_consequence": (
                "One nonzero QQ(u)-rational product-character section on a "
                "survivor would prove its arithmetic product-twist rank is exactly one."
            ),
        },
        "height_gate": {
            "arithmetic_genus_chi": 4,
            "formula_for_rational_sections": "height(P)=8+2*(P.O)",
            "height_eight_equivalence": "height(P)=8 iff P.O=0",
            "direct_height_eight_degree_box": {"degree_X_at_most": 8, "degree_Y_at_most": 12},
            "first_higher_pole_box": {
                "P.O": 1,
                "height": 10,
                "denominator_degree": 1,
                "degree_X_numerator_at_most": 10,
                "degree_Y_numerator_at_most": 15,
            },
            "second_higher_pole_box": {
                "P.O": 2,
                "height": 12,
                "denominator_degree": 2,
                "degree_X_numerator_at_most": 12,
                "degree_Y_numerator_at_most": 18,
            },
            "interaction_with_prior_exhaustion": (
                "The zero-Tate carrier certificates exclude only the zero class in "
                "the height-eight P.O=0 box. Any height-eight section on these five "
                "targets would therefore be a nonzero Tate class. Higher-pole rational "
                "sections begin at height ten and were not scanned."
            ),
        },
        "constructor_rankings": {
            "smallest_primitive_u_first": [row["pair_key"] for row in seed_order],
            "smallest_generator_canonical_height_first": [
                row["pair_key"] for row in asymptotic_order
            ],
            "primary_target": "alternate-orbit-19bad:alternate-orbit-083ad",
            "primary_reason": (
                "Its primitive generator maps to a 124-bit u-value, uniquely below "
                "the 635--792-bit range of the other survivors."
            ),
            "secondary_target": "alternate-orbit-11ee2:alternate-orbit-0c36e",
            "secondary_reason": (
                "It has the smallest primitive canonical height and therefore the "
                "best asymptotic integer-parameter growth among the five degree-four maps."
            ),
        },
        "next_exact_gate": (
            "Compute the target-specific two-Selmer/Kummer quotient for the primary "
            "twist, then solve only its nonzero-class slice in the P.O=0 degree-(8,12) "
            "box. If that box is empty, move to the P.O=1 height-ten denominator slice; "
            "do not enlarge the integral coboundary atlas."
        ),
        "proof_boundary": (
            "The residue-field factorizations and resulting rational-component and "
            "height conclusions, and the arithmetic rank-at-most-one bounds are exact. "
            "They do not compute a two-Selmer group, construct a product-character "
            "section, turn a persistent Tate factor into a characteristic-zero "
            "divisor, or exclude height at least ten."
        ),
        "inputs": {display_path(path): digest(path) for path in inputs},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    payload["reproducing_command"] = shlex.join(sys.argv)
    if args.check:
        stored = json.loads(args.output.read_text())
        stored.pop("reproducing_command", None)
        payload.pop("reproducing_command", None)
        if stored != payload:
            raise ArithmeticError("stored Galois/height gate does not replay")
        print(
            "R17PRODUCTGALOISCHECK|survivors=5|branch_places=10"
            "|height8=po0|status=PASS",
            flush=True,
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"R17PRODUCTGALOIS|survivors=5|branch_places=10|height8=po0"
        f"|output={display_path(args.output)}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
