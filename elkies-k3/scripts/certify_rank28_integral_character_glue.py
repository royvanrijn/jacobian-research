#!/usr/bin/env python3
"""Extract the common integral character glue behind the rank-28 targets.

The input genus-one-bisection certificate constructs, for each public
exceptional direction Q_i, a double cover with a lifted section R_i and
common trace tau.  This checker turns

    R_i + sigma(R_i) = tau

into an exact rank-two lattice statement.  It is a derived certificate: the
input remains responsible for the equation identities, cover smoothness, and
height computation.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import reduce
from hashlib import sha256
import json
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PILOT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
)
R17_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-r17-rank28-integral-character-glue-v1.json"
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def matmul(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def transpose(value: list[list[Fraction]]) -> list[list[Fraction]]:
    return [list(row) for row in zip(*value)]


def gram_change(
    change: list[list[Fraction]], gram: list[list[Fraction]]
) -> list[list[Fraction]]:
    return matmul(matmul(change, gram), transpose(change))


def determinant_2(value: list[list[Fraction]]) -> Fraction:
    return value[0][0] * value[1][1] - value[0][1] * value[1][0]


def quadratic_norm(vector: list[int], gram: list[list[int]]) -> int:
    return sum(
        vector[i] * gram[i][j] * vector[j]
        for i in range(len(vector))
        for j in range(len(vector))
    )


def lcm(left: int, right: int) -> int:
    return abs(left * right) // gcd(left, right) if left and right else 0


def primitive_projective_key(coefficients: list[str]) -> list[int]:
    values = [Fraction(value) for value in coefficients]
    denominator = reduce(lcm, (value.denominator for value in values), 1)
    integers = [
        value.numerator * (denominator // value.denominator) for value in values
    ]
    content = reduce(gcd, (abs(value) for value in integers if value), 0)
    integers = [value // content for value in integers]
    first = next(value for value in integers if value)
    if first < 0:
        integers = [-value for value in integers]
    return integers


def rows(value: list[list[Fraction]]) -> list[list[str]]:
    return [[str(entry) for entry in row] for row in value]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    pilot = json.loads(PILOT.read_text())
    assert pilot["status"] == "PASS_EXACT_R17_RANK28_GENUS_ONE_BISECTION_PILOT"
    assert pilot["selected_trace_count"] == 1
    assert pilot["target_count"] == 11
    assert pilot["successful_trace_target_pairs"] == 11
    trace = pilot["traces"][0]
    targets = trace["targets"]
    assert [target["target_label"] for target in targets] == [
        f"Q{index}" for index in range(1, 12)
    ]

    r17_gram = [
        [int(entry) for entry in line.split()]
        for line in R17_GRAM.read_text().splitlines()
        if line.strip()
    ]
    tau = [int(entry) for entry in trace["pinned_rank17_w"]]
    tau_norm = quadratic_norm(tau, r17_gram)
    assert tau_norm == 8
    assert trace["published_basis_w"] == [0, -1, 0, 0, -1] + [0] * 12

    cover = pilot["generic_cover_consequence"]
    assert cover["base_change_degree"] == 2
    assert cover["anti_invariant_height"] == 16
    assert cover["independent_from_invariant_R17"] is True
    pulled_tau_norm = cover["base_change_degree"] * tau_norm
    anti_norm = cover["anti_invariant_height"]
    assert pulled_tau_norm == anti_norm == 16

    pure = [[Fraction(16), Fraction(0)], [Fraction(0), Fraction(16)]]
    # Rows express (tau, R_i) in the eigencoordinates (tau, T_i), where
    # R_i=(tau+T_i)/2.
    tau_branch_change = [
        [Fraction(1), Fraction(0)],
        [Fraction(1, 2), Fraction(1, 2)],
    ]
    carrier = gram_change(tau_branch_change, pure)
    assert carrier == [[16, 8], [8, 8]]
    branch_change = [
        [Fraction(1, 2), Fraction(1, 2)],
        [Fraction(1, 2), Fraction(-1, 2)],
    ]
    branch_gram = gram_change(branch_change, pure)
    assert branch_gram == [[8, 0], [0, 8]]
    assert abs(determinant_2(tau_branch_change)) == Fraction(1, 2)
    assert determinant_2(pure) == 256
    assert determinant_2(carrier) == 64

    squareclass_keys: list[list[int]] = []
    records = []
    required_flags = (
        "branch_polynomial_coprime_to_surface_discriminant",
        "branch_polynomial_coprime_to_trace_denominator",
        "branch_polynomial_irreducible_over_Q",
        "branch_polynomial_squarefree",
        "exact_target_x_and_y_verified",
    )
    for target in targets:
        assert target["branch_polynomial_degree"] == 4
        assert all(target[flag] is True for flag in required_flags)
        assert target["kummer_barcode"]["equal_modulo_squares"] is True
        assert target["lifted_section"]["trace"] == "P(t,s)+P(t,-s)=tau(t)"
        key = primitive_projective_key(
            target["branch_polynomial_q_coefficients_low_to_high"]
        )
        squareclass_keys.append(key)
        key_payload = json.dumps(key, separators=(",", ":"))
        records.append(
            {
                "direction": target["target_label"],
                "source_public_point_index_one_based": target[
                    "source_public_point_index_one_based"
                ],
                "pencil_parameter_lambda": target["pencil_parameter_lambda"],
                "deck_character_squareclass_primitive_polynomial_sha256": sha256(
                    key_payload.encode()
                ).hexdigest(),
                "kummer_generator_coefficients": target["kummer_barcode"][
                    "target_generator_coefficients"
                ],
                "integral_character_glue_signature": {
                    "coordinate_convention": (
                        "Smith coordinates induced by the eigenbasis (tau,T_i)"
                    ),
                    "invariant_trace": "tau=-P2-P5",
                    "invariant_height_after_pullback": 16,
                    "anti_invariant_generator": (
                        f"T_{target['target_label']}=R-sigma(R)"
                    ),
                    "anti_invariant_height": 16,
                    "half_sum_identity": "2R=tau+T",
                    "discriminant_group_graph_generator": [8, 8],
                    "discriminant_group_moduli": [16, 16],
                    "glue_order": 2,
                    "carrier_isometry_class": "<8> direct_sum <8>",
                },
            }
        )

    # Distinct squarefree irreducible quartics represent distinct classes in
    # QQ(t)^*/QQ(t)^*2 unless they are proportional.  The primitive keys make
    # proportionality literal equality.
    assert len({tuple(key) for key in squareclass_keys}) == 11

    payload = {
        "schema": "elkies-k3.r17-rank28-integral-character-glue.v1",
        "status": "PASS_EXACT_COMMON_RANK_TWO_CHARACTER_GLUE_FOR_ALL_ELEVEN",
        "parameter": pilot["parameter"],
        "inputs": {
            display(PILOT): digest(PILOT),
            display(R17_GRAM): digest(R17_GRAM),
            display(Path(__file__)): digest(Path(__file__)),
        },
        "common_trace": {
            "published_basis_expression": "-P2-P5",
            "pinned_R17_coordinates": tau,
            "height_on_R17": tau_norm,
            "height_after_degree_two_pullback": pulled_tau_norm,
        },
        "common_rank_two_carrier": {
            "eigenbasis": [
                "tau (invariant)",
                "T_i=R_i-sigma(R_i) (anti-invariant)",
            ],
            "pure_eigensum_gram": rows(pure),
            "pure_eigensum_determinant": 256,
            "glue_subgroup": {
                "ambient_discriminant_group": "Z/16 direct_sum Z/16",
                "generator_in_smith_coordinates": [8, 8],
                "order": 2,
                "quadratic_value_in_Q_mod_2Z": "0",
                "identity": "2R_i=tau+T_i",
            },
            "saturated_carrier_gram_in_basis_tau_Ri": rows(carrier),
            "saturated_carrier_gram_in_basis_Ri_sigmaRi": rows(branch_gram),
            "saturated_carrier_isometry_class": "<8> direct_sum <8>",
            "saturated_carrier_determinant": 64,
            "index_of_pure_eigensum": 2,
            "determinant_ratio": 4,
            "involution_in_branch_basis": [[0, 1], [1, 0]],
        },
        "directions": records,
        "aggregate": {
            "direction_count": len(records),
            "common_glue_signature_count": len(records),
            "distinct_deck_character_squareclasses": len(
                {tuple(key) for key in squareclass_keys}
            ),
            "shared_trace_classes": 1,
        },
        "mechanism_assessment": {
            "positive": (
                "All eleven exact lifts fitted from the chosen trace pencil have "
                "the same rank-two index-two character-glue pattern, invariant "
                "trace, and abstract carrier lattice."
            ),
            "essential_distinction": (
                "The eleven quartics define eleven distinct deck characters. This is "
                "one repeated local integral pattern, not one double cover carrying "
                "eleven anti-invariant directions, and not an intrinsic invariant of "
                "an isolated fibre point."
            ),
            "reglue_boundary": (
                "This is character saturation by adjoining a half-sum, not a root-"
                "annihilating bridge reglue: it enlarges the pure eigensum and cannot "
                "remove a vector already present there."
            ),
            "search_consequence": (
                "A mechanism-first search should enumerate rational members of a fixed "
                "low-complexity trace pencil and seek simultaneous rational fibres of "
                "their distinct genus-one deck characters, rather than add another "
                "Nagao feature."
            ),
        },
        "proof_boundary": {
            "proved": (
                "The rank-two Gram, index-two graph glue, common trace, and pairwise "
                "distinct deck squareclasses follow exactly from the pinned R17 Gram "
                "and the exact genus-one-bisection certificate, relative to its "
                "chosen trace pencil."
            ),
            "inherited": (
                "Equation identities, smooth irreducible covers, literal target "
                "specialization, and anti-invariant height 16 are inherited from the "
                "input theorem certificate."
            ),
            "not_proved": (
                "No one cover supplies more than one certified exceptional character, "
                "no new specialization is found, and no rank-32 lower or upper bound "
                "is claimed."
            ),
        },
        "reproducing_command": (
            "python3 elkies-k3/scripts/certify_rank28_integral_character_glue.py"
        ),
    }

    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists() or output.read_text() != rendered:
            raise SystemExit(f"stale or missing output: {display(output)}")
        print("R17R28GLUE|status=PASS_CHECK|directions=11|distinct_characters=11")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered)
    print(f"R17R28GLUE|status=PASS|output={display(output)}")


if __name__ == "__main__":
    main()
