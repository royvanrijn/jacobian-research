#!/usr/bin/env sage
"""Canonicalize the known alternate rootless frame in pinned H3/R17 NS coordinates.

This is not a complete Kneser--Nishiyama classification.  It turns the exact
q80 alternate-q4/q6 certificate into a first classification seed: a primitive
U embedding in the pinned NS lattice, a rootless rank-17 frame distinct from
the published R17 frame at the J2 (frame-isometry) level, and an equation-cost
handoff to the existing characteristic-zero compiler.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, matrix, pari


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "elkies-k3/data/fibrations"
GENERATED = ROOT / "artifacts/generated-results"

ALTERNATE = GENERATED / "q80-alternate-fifth-q6-rootless-transport.json"
PINNED_FRAME = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
Q80_FRAME = DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt"
PINNED_TO_Q80 = DATA / "kumar_q80_rootless_target_to_q80_ns_transport.txt"
OUTPUT = GENERATED / "elkies-k3-other-rank17-candidate.json"


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

alternate_data = json.loads(ALTERNATE.read_text())
alternate_frame = matrix(ZZ, alternate_data["rootless_frame"])
alternate_to_q80 = matrix(ZZ, alternate_data["rootless_to_q80_ns_transport"])
pinned_frame = load_matrix(PINNED_FRAME)
q80_frame = load_matrix(Q80_FRAME)
pinned_to_q80 = load_matrix(PINNED_TO_Q80)
U = matrix(ZZ, [[0, 1], [1, 0]])

alternate_ns = block_diagonal_matrix(U, -alternate_frame)
pinned_ns = block_diagonal_matrix(U, -pinned_frame)
q80_ns = block_diagonal_matrix(U, -q80_frame)

assert alternate_frame.is_positive_definite()
assert alternate_frame.det() == 948
assert ZZ(pari(alternate_frame).qfminim(2)[0]) == 0
assert pari(alternate_frame).qfisom(pari(pinned_frame)) == 0
assert alternate_to_q80 * q80_ns * alternate_to_q80.transpose() == alternate_ns
assert abs(alternate_to_q80.det()) == 1
assert pinned_to_q80 * q80_ns * pinned_to_q80.transpose() == pinned_ns
assert abs(pinned_to_q80.det()) == 1

alternate_to_pinned = (alternate_to_q80 * pinned_to_q80.inverse()).change_ring(ZZ)
assert alternate_to_pinned * pinned_ns * alternate_to_pinned.transpose() == alternate_ns
assert abs(alternate_to_pinned.det()) == 1

fibre = alternate_to_pinned.row(0)
mate = alternate_to_pinned.row(1)
zero = mate - fibre
assert fibre * pinned_ns * fibre == 0
assert zero * pinned_ns * zero == -2
assert fibre * pinned_ns * zero == 1

result = {
    "schema": "elkies-k3.other-rank17-candidate.v1",
    "status": "PASS_EXACT_DISTINCT_ROOTLESS_U_EMBEDDING",
    "classification_scope": {
        "proved": (
            "One primitive rootless U embedding is transported into the pinned "
            "H3/R17 Neron--Severi lattice. Its determinant-948 positive frame "
            "has no norm-two roots and is not integrally isometric to the "
            "published R17 frame, so it is a distinct J2 frame class."
        ),
        "not_proved": (
            "No completeness or J1 classification up to surface automorphisms "
            "is claimed. The rank-29 curve is not yet identified as a fibre."
        ),
    },
    "pinned_ns_coordinates": {
        "fibre": list(map(int, fibre)),
        "isotropic_mate": list(map(int, mate)),
        "zero": list(map(int, zero)),
        "alternate_basis_to_pinned_ns": rows(alternate_to_pinned),
        "transport_determinant": int(alternate_to_pinned.det()),
    },
    "frame": {
        "rank": 17,
        "determinant": int(alternate_frame.det()),
        "norm_two_root_count": 0,
        "mw_rank_by_shioda_tate_at_rho19": 17,
        "integrally_isometric_to_published_r17_frame": False,
        "sha256": alternate_data["rootless_frame_sha256"],
    },
    "equation_cost": {
        "last_neighbor": "alternate A1/MW16 --q6--> rootless/MW17",
        "old_fibre_degree_after_weyl_reduction": 2,
        "reduced_divisor_old_zero_pairing": 1,
        "horizontal_section_old_zero_pairing": 4,
        "horizontal_identity": "D=O+S-F",
        "nef_ambient_line_bundle_dimension": 5,
        "local_elementary_transforms": 3,
        "expected_pencil_dimension": 2,
        "generic_fibre_generators": ["1", "(y+y(S))/(x-x(S))"],
        "exact_cost_replay": "elkies-k3/scripts/verify_q80_alternate_final_q6_nef.sage",
        "remaining_gate": (
            "Lift the generic characteristic-zero parent and its three saturated "
            "local transforms, compile the rootless Weierstrass equation, then "
            "solve j_candidate(t)=j(rank29)."
        ),
    },
    "kneser_nishiyama_next_gate": {
        "objective": (
            "Classify all rootless rank-17 J2 frame classes by primitive "
            "embeddings into Niemeier lattices, then refine to marked J1 "
            "classes only where surface automorphisms are controlled."
        ),
        "deduplication_key": (
            "integral frame isometry plus discriminant-form/glue data; retain "
            "the full pinned U marking for equation construction"
        ),
        "known_rootless_frame_classes_lower_bound": 2,
    },
    "inputs": {
        str(path.relative_to(ROOT)): file_sha256(path)
        for path in (ALTERNATE, PINNED_FRAME, Q80_FRAME, PINNED_TO_Q80)
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python -c "
        "'from sage.all import *; globals()[\"__file__\"]="
        "\"/home/royvanrijn/src/jacobian-research/elkies-k3/scripts/"
        "canonicalize_other_rank17_candidate.sage\"; load(__file__)'"
    ),
}

serialized = json.dumps(result, indent=2, sort_keys=True, default=int) + "\n"
if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != serialized:
        raise SystemExit("pinned other-rank17 candidate artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(serialized)

print(
    "OTHERR17|rootless=1|rank=17|det=948|published_frame_isometric=0|"
    f"transport_det={alternate_to_pinned.det()}|"
    "known_rootless_J2_classes_at_least=2|status=PASS",
    flush=True,
)
