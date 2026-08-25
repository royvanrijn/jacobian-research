#!/usr/bin/env sage -python
"""Score larger-q explicit-curve-nef exits from the orbit12 child."""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--gate", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
GATE = args.gate.resolve()
OUTPUT = args.output.resolve()
FRAME = GENERATED / "elkies-k3-h3-a11-q8-orbit12-explicit-zero-frame.txt"


def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(value) for value in line.split()] for line in path.read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")],
    )


frame = load_matrix(FRAME)
root_rank = 10
root = frame[:root_rank, :root_rank]
coupling = frame[:root_rank, root_rank:]
lattice = IntegralLattice(root)
root_edges = [
    (left, right)
    for left in range(root_rank)
    for right in range(left + 1, root_rank)
    if root[left, right] == -1
]


def vertical_layers(coefficients):
    magnitudes = [abs(ZZ(value)) for value in coefficients]
    previous = 0
    total = 0
    for level in sorted(set(value for value in magnitudes if value)):
        active = set(index for index, value in enumerate(magnitudes) if value >= level)
        count = 0
        while active:
            count += 1
            todo = [active.pop()]
            while todo:
                node = todo.pop()
                for left, right in root_edges:
                    other = right if left == node else left if right == node else None
                    if other in active:
                        active.remove(other)
                        todo.append(other)
        total += (level - previous) * count
        previous = level
    return int(total)


profile_cache = {}


def horizontal_profile(z):
    key = tuple(z)
    if key in profile_cache:
        return profile_cache[key]
    z = vector(ZZ, z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * frame[:, :root_rank]) * root.inverse()
    iterator = lattice.enumerate_close_vectors(-dual)
    minimum = None
    choices = []
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * (17 - root_rank))
        norm = QQ(lifted * frame * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        pole = (norm - 4) / 2
        if pole in ZZ and pole >= 0:
            choices.append((ZZ(pole), vector(ZZ, [ZZ(pole) + 1, 1] + list(lifted))))
    assert choices
    choices.sort(key=lambda item: (item[0], tuple(item[1])))
    profile_cache[key] = choices[0]
    return choices[0]


gate = json.loads(GATE.read_text())
assert gate["status"] == "PASS_EXACT_A5A5_EXPLICIT_ZERO_LARGE_Q_GATE"
records = []
for source in gate["survivors"]:
    fibre = vector(ZZ, source["fibre_in_parent"])
    pole, section = horizontal_profile(source["mw_projection"])
    residual = fibre - vector(ZZ, [-1, 1] + [0] * 17) - section
    vertical = vector(ZZ, residual[2:2 + root_rank])
    layers = vertical_layers(vertical)
    rr = 2 + 2 * int(pole) + layers
    zeros = source["explicit_degree_zero_curves"]
    ones = source["explicit_degree_one_curves"]
    terms = {
        "P_dot_O": 900 * int(pole),
        "horizontal_degree": 500,
        "RR_ambient": 120 * rr,
        "vertical_layers": 60 * layers,
        "vertical_support": 25 * sum(value != 0 for value in vertical),
        "coordinate_growth": source["coordinate_growth_max"],
        "no_explicit_degree_one_curve": 4000 if not ones else 0,
        "explicit_degree_one_credit": -500 * min(len(ones), 6),
        "explicit_degree_zero_credit": -100 * min(len(zeros), 12),
    }
    records.append({
        **source,
        "horizontal": {
            "P_dot_O": int(pole),
            "section": list(map(int, section)),
            "vertical": list(map(int, vertical)),
            "vertical_layers": layers,
            "vertical_support": sum(value != 0 for value in vertical),
        },
        "expected_RR_ambient": rr,
        "equation_cost_terms": terms,
        "equation_cost_score": int(sum(terms.values())),
    })

records.sort(key=lambda item: (
    item["equation_cost_score"],
    item["marked_target_degrees"]["pinned_R17"],
    item["candidate_id"]["orbit_index"],
))
pareto = []
for item in records:
    if not any(
        other["equation_cost_score"] <= item["equation_cost_score"]
        and other["marked_target_degrees"]["pinned_R17"] <= item["marked_target_degrees"]["pinned_R17"]
        and (
            other["equation_cost_score"] < item["equation_cost_score"]
            or other["marked_target_degrees"]["pinned_R17"] < item["marked_target_degrees"]["pinned_R17"]
        )
        for other in records
    ):
        pareto.append(item)
pareto.sort(key=lambda item: item["marked_target_degrees"]["pinned_R17"])
payload = {
    "schema": "elkies-k3.h3-a5a5-explicit-zero-large-q-cost.v1",
    "status": "PASS_EXACT_A5A5_EXPLICIT_ZERO_LARGE_Q_COST_SCORING",
    "candidate_count": len(records),
    "ranked_candidates": records,
    "pinned_degree_equation_cost_pareto_front": pareto,
    "proof_boundary": (
        "Exact P.O and vertical profiles plus an RR ambient estimate for fibres passing the explicit-curve gate. "
        "Child root data, all-section nefness, marked U, and full transports remain uncertified."
    ),
    "inputs": {
        "paths": [str(GATE.relative_to(ROOT)), str(FRAME.relative_to(ROOT))],
        "sha256": {
            str(GATE.relative_to(ROOT)): hashlib.sha256(GATE.read_bytes()).hexdigest(),
            str(FRAME.relative_to(ROOT)): hashlib.sha256(FRAME.read_bytes()).hexdigest(),
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = records[0]
print(
    "A5A5LARGECOST|candidates={}|best=q{}o{}|score={}|PO={}|RR={}|pinned={}|deg0={}|deg1={}|pareto={}|status={}".format(
        len(records), best["candidate_id"]["q"], best["candidate_id"]["orbit_index"],
        best["equation_cost_score"], best["horizontal"]["P_dot_O"], best["expected_RR_ambient"],
        best["marked_target_degrees"]["pinned_R17"], len(best["explicit_degree_zero_curves"]),
        len(best["explicit_degree_one_curves"]), len(pareto), payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{OUTPUT}", flush=True)
