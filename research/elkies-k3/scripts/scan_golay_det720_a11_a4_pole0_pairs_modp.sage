#!/usr/bin/env sage-python
"""Exhaust the two pole-zero sections for the G720-S0260 fibre chart.

Both basis sections of the exact ``A11+A4/MW2`` source meet the identity
components of ``I12`` and ``I5`` and have pole order zero.  In the physical
basis selected by the Smith-quotient audit their height Gram is
``[[4,2],[2,4]]``.  Thus two polynomial sections realize the full MW lattice
exactly when both avoid the reducible-fibre nodes and have no intersection on
smooth fibres; their height determinant is then twelve, equal to the source
regulator determinant.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIBRES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-a11-a4-source-ansatz-mod5-v1.json"
)
DEFAULT_POLES = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-source-poles-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-a11-a4-pole0-pairs-mod5-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def polynomial_roots(right):
    if not right.is_square():
        return []
    positive = right.sqrt()
    return [positive] if not positive else [positive, -positive]


def serialize(poly):
    return [int(value) for value in poly]


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--fibres", type=Path, default=DEFAULT_FIBRES)
parser.add_argument("--poles", type=Path, default=DEFAULT_POLES)
parser.add_argument("--source-id", default="G720-S0260")
parser.add_argument("--quadratic-twist", type=int, default=1)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

fibres_path = arguments.fibres.resolve()
poles_path = arguments.poles.resolve()
output_path = arguments.output.resolve()
fibres = json.loads(fibres_path.read_text())
poles = json.loads(poles_path.read_text())
if fibres.get("schema") != "elkies-k3.golay-det720-a11-a4-source-ansatz-modp.v1":
    raise ValueError("unexpected G720 A11+A4 fibre schema")
if not fibres["scan"]["exhausted"] or fibres["source"]["source_id"] != arguments.source_id:
    raise ValueError("section scan requires the exhaustive selected-source fibre chart")
pole_entry = next(
    row
    for row in poles["audits"]
    if row["source_id"] == arguments.source_id
    and row["source_gram_sha256"] == fibres["source"]["source_gram_sha256"]
)
audit = pole_entry["audit"]
if audit["basis_sorted_pole_profile"] != [0, 0] or audit["torsion_order"] != 1:
    raise ValueError("selected source lost its primitive pole-zero basis")
height = matrix(QQ, audit["height_gram"])
basis_coordinates = matrix(
    QQ, [section["free_mw_coordinates"] for section in audit["basis"]]
)
basis_height = basis_coordinates * height * basis_coordinates.transpose()
if basis_height != matrix(QQ, [[4, 2], [2, 4]]) or basis_height.det() != 12:
    raise ArithmeticError("unexpected physical basis height Gram")
if any(any(section["simple_root_pairings"]) for section in audit["basis"]):
    raise ArithmeticError("selected basis no longer meets identity components")

prime = int(fibres["prime"])
field = GF(prime)
twist = field(arguments.quadratic_twist)
if not twist:
    raise ValueError("quadratic twist must be nonzero")
ring = PolynomialRing(field, "t")
t = ring.gen()

models = []
x_candidates = 0
sections_total = 0
pairs_total = 0
pairs_meeting_singular_fibres = 0
for example_index, example in enumerate(fibres["examples"]):
    A = twist**2 * ring(example["A_coefficients_low_to_high"])
    B = twist**3 * ring(example["B_coefficients_low_to_high"])
    discriminant_core = 4 * A**3 + 27 * B**2
    node_zero = -field(3) * B[0] / (field(2) * A[0])
    node_infinity = -field(3) * B[12] / (field(2) * A[8])
    sections = []
    for coefficients in itertools.product(field, repeat=5):
        X = ring(list(coefficients))
        x_candidates += 1
        for Y in polynomial_roots(X**3 + A * X + B):
            identity_zero = not (X[0] == node_zero and Y[0] == 0)
            identity_infinity = not (
                X[4] == node_infinity and Y[6] == 0
            )
            if not identity_zero or not identity_infinity:
                continue
            sections.append({"X": X, "Y": Y})
    if not sections:
        continue
    sections_total += len(sections)
    pairs = []
    for left, right in itertools.combinations(range(len(sections)), 2):
        left_section = sections[left]
        right_section = sections[right]
        common = (left_section["X"] - right_section["X"]).gcd(
            left_section["Y"] - right_section["Y"]
        )
        if common.gcd(discriminant_core).degree() != 0:
            pairs_meeting_singular_fibres += 1
            continue
        intersection = int(common.degree())
        if intersection != 0:
            continue
        pairs_total += 1
        pairs.append(
            {
                "left_section_index": left,
                "right_section_index": right,
                "intersection_on_smooth_fibres": 0,
                "component_cross_correction": "0",
                "shioda_height_pairing": "2",
                "height_gram": [["4", "2"], ["2", "4"]],
                "height_determinant": "12",
                "mw_index_from_determinant": 1,
            }
        )
    models.append(
        {
            "example_index": example_index,
            "sections": [
                {
                    "X_coefficients_low_to_high": serialize(section["X"]),
                    "Y_coefficients_low_to_high": serialize(section["Y"]),
                    "pole_order": 0,
                    "component_depths_at_I12_I5": [0, 0],
                    "shioda_height": "4",
                }
                for section in sections
            ],
            "marked_mw2_pairs": pairs,
        }
    )

output = {
    "schema": "elkies-k3.golay-det720-a11-a4-pole0-pairs-modp.v1",
    "status": (
        "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW2_PAIRS"
        if pairs_total
        else "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW2_PAIR_LOCUS"
    ),
    "prime": prime,
    "quadratic_twist": int(twist),
    "quadratic_twist_square_class": "square" if twist.is_square() else "nonsquare",
    "inputs": {
        relative(fibres_path): digest(fibres_path),
        relative(poles_path): digest(poles_path),
    },
    "source": {
        "source_id": arguments.source_id,
        "source_gram_sha256": fibres["source"]["source_gram_sha256"],
        "root_type": "A11+A4",
        "mw_rank": 2,
        "basis_pole_profile": [0, 0],
        "basis_on_identity_components": True,
        "physical_basis_height_gram": [["4", "2"], ["2", "4"]],
        "physical_basis_height_determinant": "12",
        "required_smooth_pair_intersection": 0,
    },
    "scope": {
        "fibre_census_exhaustive": True,
        "stored_fibre_models": len(fibres["examples"]),
        "polynomial_X_candidates_per_model": prime**5,
        "all_polynomial_Y_square_roots_retained": True,
        "identity_component_gate_at_I12_and_I5": True,
        "all_unordered_section_pairs_tested": True,
    },
    "accounting": {
        "models_with_sections": len(models),
        "X_candidates_scanned": x_candidates,
        "pole_zero_sections": sections_total,
        "marked_mw2_pairs": pairs_total,
        "pairs_meeting_singular_fibres": pairs_meeting_singular_fibres,
    },
    "models": models,
    "proof_boundary": {
        "proved": (
            "For every model in the exhaustive normalized fibre chart, all "
            "degree-at-most-four X polynomials and every polynomial Y square "
            "root are tested. Every marked pair has the exact determinant-twelve "
            "height Gram and therefore index one in the declared MW2 lattice."
        ),
        "not_proved": (
            "A finite-field marked pair is not a characteristic-zero family, a "
            "Q-rational source marking, a rational parameterization, a target "
            "multisection spectrum, or a physical neighbour corridor."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/scan_golay_det720_a11_a4_pole0_pairs_modp.sage"
        + (f" --quadratic-twist {int(twist)}" if twist != 1 else "")
        + (
            f" --fibres {relative(fibres_path)}"
            if fibres_path != DEFAULT_FIBRES.resolve()
            else ""
        )
        + (
            f" --output {relative(output_path)}"
            if output_path != DEFAULT_OUTPUT.resolve()
            else ""
        )
    ),
}
serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("Golay-720 A11+A4 marked-pair artifact is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

print(
    "GOLAY720A11A4PAIR|"
    f"p={prime}|twist={int(twist)}|models={len(models)}|"
    f"sections={sections_total}|pairs={pairs_total}|"
    f"status={'PASS' if pairs_total else 'EMPTY'}",
    flush=True,
)
