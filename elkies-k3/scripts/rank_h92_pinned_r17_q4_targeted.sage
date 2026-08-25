#!/usr/bin/env sage -python
"""Exactly rank every nef pinned-R17 q4 bisection presentation by marked targets."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "elkies-k3/data/lattice"
GENERATED = ROOT / "artifacts/generated-results"
PINNED = DATA / "rank17_gram.txt"
SHORT_GRAM = DATA / "short_vector_basis_gram.txt"
SHORT_CHANGE = DATA / "short_vector_basis_coords.txt"
MARKING = GENERATED / "elkies-k3-h3-pinned-r17-equation-marking.json"
OUTPUT = GENERATED / "elkies-k3-h3-pinned-r17-q4-degree2-targeted-ranking.json"
INPUTS = (PINNED, SHORT_GRAM, SHORT_CHANGE, MARKING)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


def parity_mask(value):
    answer = 0
    for index, entry in enumerate(value):
        if ZZ(entry) % 2:
            answer |= 1 << index
    return answer


def qform(gram):
    coefficients = []
    for row in range(gram.nrows()):
        for column in range(row, gram.ncols()):
            coefficients.append(
                gram[row, row] // 2 if row == column else gram[row, column]
            )
    return QuadraticForm(ZZ, gram.nrows(), coefficients)


pinned = load_matrix(PINNED)
short_gram = load_matrix(SHORT_GRAM)
change = load_matrix(SHORT_CHANGE)
marking = json.loads(MARKING.read_text())
assert marking["status"] == "PASS_EXACT_REVERSE_HUB_EQUATION_MARKING"
assert abs(change.det()) == 1 and short_gram == change * pinned * change.transpose()
g = block_diagonal_matrix(U2, -pinned)
targets = {
    name: vector(ZZ, value)
    for name, value in marking["target_fibres_in_root_adapted_hub"].items()
}

# Sage's list is indexed by q(v)=v.M.v/2.  Index four is the complete
# norm-eight shell, one representative from each +/- pair.
shells = qform(short_gram).short_vector_list_up_to_length(5, True)
assert len(shells[2]) == 1311 and len(shells[3]) == 26672
assert len(shells[4]) == 230040
excluded_masks = {
    parity_mask(vector(ZZ, value)) for value in tuple(shells[2]) + tuple(shells[3])
}

linear = {
    name: pinned * vector(ZZ, target[2:]) for name, target in targets.items()
}
constant = {name: 2 * (target[0] + target[1]) for name, target in targets.items()}
rankings = {name: [] for name in targets}


def ranking_key(item, name):
    return (
        item["marked_target_degrees"][name],
        item["coordinate_growth_max"],
        item["candidate_id"]["norm8_pair_index"],
        -item["candidate_id"]["sign"],
    )


def retain(item, name):
    bucket = rankings[name]
    bucket.append(item)
    if len(bucket) >= 400:
        bucket.sort(key=lambda candidate: ranking_key(candidate, name))
        del bucket[200:]


nef_count = 0
primitive_count = 0
for pair_index, short_values in enumerate(shells[4], start=1):
    short_w = vector(ZZ, short_values)
    if parity_mask(short_w) in excluded_masks:
        continue
    nef_count += 2
    pinned_w0 = short_w * change
    for sign in (1, -1):
        pinned_w = sign * pinned_w0
        fibre = vector(ZZ, [2, 2] + list(pinned_w))
        assert fibre * g * fibre == 0
        if gcd(tuple(g * fibre)) != 1:
            continue
        primitive_count += 1
        degrees = {
            name: int(constant[name] - pinned_w * linear[name])
            for name in targets
        }
        assert min(degrees.values()) >= 0
        record = {
            "candidate_id": {
                "q": 4,
                "old_fibre_degree": 2,
                "norm8_pair_index": pair_index,
                "sign": sign,
            },
            "fibre_in_pinned_R17": list(map(int, fibre)),
            "P_dot_O": 0,
            "minimum_section_intersection": "0",
            "marked_target_degrees": degrees,
            "coordinate_growth_max": int(max(abs(value) for value in fibre)),
        }
        for name in targets:
            retain(record, name)

for name in targets:
    rankings[name].sort(key=lambda item: ranking_key(item, name))
    rankings[name] = rankings[name][:200]

payload = {
    "schema": "elkies-k3.h3-pinned-r17-q4-targeted-ranking.v1",
    "status": "PASS_EXACT_PINNED_R17_Q4_TARGETED_RANKING",
    "q": 4,
    "old_fibre_degree": 2,
    "norm8_unoriented_pair_count": len(shells[4]),
    "excluded_shorter_parity_coset_count": len(excluded_masks),
    "nef_oriented_candidate_count": nef_count,
    "primitive_nef_oriented_candidate_count": primitive_count,
    "rankings_top_200": rankings,
    "proof_boundary": (
        "Exact complete norm-eight enumeration in a unimodular short basis. "
        "For D=(2,2,w), all-section nefness is equivalent to the parity coset "
        "of w containing no norm-four or norm-six vector; the surviving minimum "
        "is norm eight and gives minimum section intersection zero. Marked target "
        "degrees and primitivity are exact. Selected candidates still require "
        "child-root and full bidirectional transport certificates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for target in ("orbit12", "equation_A11", "mw2_a5_a4_2a3_semistable", "q25_mw7"):
    best = rankings[target][0]
    print(
        "R17Q4|target={}|degree={}|pair={}|sign={}|max={}|status=PASS".format(
            target, best["marked_target_degrees"][target],
            best["candidate_id"]["norm8_pair_index"], best["candidate_id"]["sign"],
            best["coordinate_growth_max"],
        ),
        flush=True,
    )
print(
    "R17Q4|pairs={}|nef_oriented={}|primitive={}|status={}|output={}".format(
        len(shells[4]), nef_count, primitive_count, payload["status"], OUTPUT,
    ),
    flush=True,
)
