#!/usr/bin/env sage-python
"""Certify the exact determinant-720 Golay Neron--Severi lattice.

This is the saturation gate for the abstract Golay-octad K3 class.  It keeps
three notions separate: primitive closure in ``N(24A1)``, primitive embedding
in the K3 lattice, and maximality inside the Neron--Severi lattice of a chosen
arithmetic specialization.  The first two pass here.  The last one fails for
the currently stored ``s6=10`` rational specialization and is not asserted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import Genus, QQ, ZZ, QuadraticForm, block_diagonal_matrix, matrix


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_TARGET = (
    ROOT / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-det720-ns-saturation-v1.json"
)
U = matrix(ZZ, ((0, 1), (1, 0)))


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primitive_closure_index(basis) -> int:
    original = matrix(ZZ, basis)
    saturated = original.row_module(ZZ).saturation().basis_matrix()
    coordinates = original * saturated.pseudoinverse()
    if not all(entry.denominator() == 1 for entry in coordinates.list()):
        raise ArithmeticError("primitive-closure coordinates are not integral")
    return abs(int(matrix(ZZ, coordinates).det()))


def discriminant_form_key(gram):
    normal = Genus(gram).discriminant_form().normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": [
            [str(entry) for entry in row]
            for row in normal.gram_matrix_quadratic().rows()
        ],
        "value_module": str(normal.value_module_qf()),
    }


def negate_form_key(key):
    return {
        "invariants": key["invariants"],
        "quadratic_gram": [
            [str(-QQ(entry)) for entry in row]
            for row in key["quadratic_gram"]
        ],
        "value_module": key["value_module"],
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

target_path = arguments.target.resolve()
output_path = arguments.output.resolve()
target = json.loads(target_path.read_text())
assert target["status"] == "PASS_EXACT_GOLAY_OCTAD_RANK17_DET720_LATTICE_DESIGN"

auxiliary = target["auxiliary"]
frame = target["frame"]
auxiliary_gram = matrix(ZZ, auxiliary["gram"])
frame_gram = matrix(ZZ, frame["gram"])
ternary_gram = matrix(ZZ, target["k3_realizability"]["transcendental_gram"])

auxiliary_index = primitive_closure_index(auxiliary["ambient_basis"])
frame_index = primitive_closure_index(frame["ambient_basis"])
assert auxiliary_index == auxiliary["primitive_closure_index"] == 1
assert frame_index == 1
assert auxiliary_gram.det() == frame_gram.det() == 720
assert frame_gram.is_positive_definite()

ns_gram = block_diagonal_matrix(U, -frame_gram)
assert ns_gram.det() == 720
assert QuadraticForm(QQ, ns_gram).signature() == -17
assert ternary_gram.det() == -720
assert QuadraticForm(QQ, ternary_gram).signature() == 1

frame_form = discriminant_form_key(frame_gram)
ternary_form = discriminant_form_key(ternary_gram)
assert frame_form == frame["discriminant_form_normal_key"]
assert ternary_form == frame_form
assert frame_form["invariants"] == [2, 6, 60]
ns_form = negate_form_key(frame_form)

payload = {
    "schema": "elkies-k3.golay-det720-ns-saturation.v1",
    "status": "PASS_EXACT_PRIMITIVE_K3_NS_DET720_AND_DISCRIMINANT_FORM",
    "reproduce": (
        "sage -python elkies-k3/scripts/certify_golay_det720_ns_saturation.sage"
    ),
    "inputs": {relative(target_path): digest(target_path)},
    "ambient_saturation": {
        "auxiliary_rank": 7,
        "auxiliary_primitive_closure_index_in_N24A1": auxiliary_index,
        "orthogonal_frame_rank": 17,
        "orthogonal_frame_primitive_closure_index_in_N24A1": frame_index,
        "common_determinant": 720,
        "argument": (
            "The auxiliary and complement bases are saturated in the integral "
            "Niemeier coordinate lattice. Independently, an orthogonal complement "
            "is the kernel of an integral pairing map and is therefore primitive."
        ),
    },
    "neron_severi": {
        "gram_construction": "U plus the negative Golay frame",
        "rank": 19,
        "signature": [1, 18],
        "literal_gram_determinant": int(ns_gram.det()),
        "signed_discriminant": -abs(int(ns_gram.det())),
        "signed_discriminant_convention": (
            "minus the absolute determinant, matching the elliptic-K3 source "
            "ledgers; the literal signature-(1,18) Gram determinant is +720"
        ),
        "discriminant_group_invariants": frame_form["invariants"],
        "discriminant_group_primary_decomposition": [2, 2, 4, 3, 3, 5],
        "discriminant_form_normal_key": ns_form,
        "discriminant_form_value_module": "Q/2Z",
        "isometry_class_unique_in_genus": True,
        "uniqueness_argument": (
            "The lattice is indefinite of rank 19 and its discriminant length is "
            "3, so the rank-at-least-length-plus-two uniqueness criterion applies."
        ),
    },
    "transcendental_complement": {
        "rank": 3,
        "signature": [2, 1],
        "gram": [[int(entry) for entry in row] for row in ternary_gram.rows()],
        "determinant": int(ternary_gram.det()),
        "discriminant_form_normal_key": ternary_form,
    },
    "k3_embedding": {
        "relation": "q_T = q_frame = -q_NS",
        "primitive": True,
        "argument": (
            "The displayed exact normal forms identify q_T with -q_NS. Gluing "
            "NS and T along this anti-isometry gives an even unimodular lattice "
            "of signature (3,19), hence the K3 lattice, with both summands primitive."
        ),
        "generic_period_has_exact_NS": True,
        "picard_rank": 19,
    },
    "proof_boundary": {
        "proved": (
            "The abstract Golay lattice is saturated in N(24A1), has signed "
            "discriminant -720 (literal determinant +720), has the displayed exact "
            "finite quadratic module, and embeds primitively in the K3 lattice. A "
            "generic period therefore gives a complex Picard-19 K3 with exact NS."
        ),
        "not_proved": (
            "This certificate does not descend that generic K3 to QQ. It does not "
            "assert that the current s6=10 rational 3I6 model has this NS; the "
            "separate saturation audit proves that model instead has determinant 20."
        ),
    },
}

encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not output_path.exists() or output_path.read_text() != encoded:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded)

print(
    "GOLAY720NS|aux_index=1|frame_index=1|det_gram=720|disc_signed=-720|"
    "A=2,6,60|primitive_K3=1|rho=19|status=PASS"
)
