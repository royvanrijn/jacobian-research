#!/usr/bin/env sage -python
"""Enumerate q4/o164 P.O=0 sections with at most one node hit.

The physical q4/o164 frame is U plus the negative definite rank-17 frame.
Thus every P.O=0 section is [1,1,v] with vHv=4.  Exact finite enumeration,
simple-component nef inequalities, and affine-component inequalities recover
the physical section shell.  Transport through the stored unimodular basis
then records each curve's degree on the q4/o1584 parent.
"""

import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
FRAME = GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-old_a11_component_8-frame.txt"
MARKING = GENERATED / "elkies-k3-h3-q4o208-q4o1584-q4o164-old_a11_component_8-marking.json"
OUTPUT = ROOT / "artifacts/local/elkies-k3/q4o164-zero-one-node-parent-degree-audit.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


H = matrix(ZZ, [list(map(int, line.split())) for line in FRAME.read_text().splitlines() if line and not line.startswith("#")])
marking = json.loads(MARKING.read_text())
assert H.nrows() == H.ncols() == 17 and H.is_positive_definite()
assert marking["status"] == "PASS_EXACT_Q4O164_PHYSICAL_EFFECTIVE_ZERO_MARKING"
B = matrix(ZZ, marking["basis_in_source"])
assert abs(B.det()) == 1

# Connected simple-root indices in the exact physical frame: A3+A1+A1+A3.
fibres = ((0, 3, 4), (1,), (2,), (5, 6, 7))
half = matrix(ZZ, pari(H).qfminim(4)[2]).transpose().rows()
vectors = [vector(ZZ, row) for row in half]
vectors += [-row for row in vectors]
vectors = [vector(ZZ, row) for row in sorted(set(map(tuple, vectors)))]

sections = []
for tail in vectors:
    if tail * H * tail != 4:
        continue
    simple_pairings = [int(tail * H.column(index)) for index in range(8)]
    if any(value < 0 for value in simple_pairings):
        continue
    fibre_hits = [sum(simple_pairings[index] for index in component) for component in fibres]
    if any(value not in (0, 1) for value in fibre_hits):
        continue
    child = vector(ZZ, [1, 1] + list(tail))
    parent = child * B
    sections.append({
        "class_in_q4o164_basis": list(map(int, child)),
        "class_in_q4o1584_basis": list(map(int, parent)),
        "simple_component_pairings": simple_pairings,
        "nonidentity_fibre_hits": fibre_hits,
        "node_hit_count": sum(fibre_hits),
        "q4o1584_parent_degree": int(parent[1]),
    })

assert len(sections) == 206
hit_histogram = Counter(row["node_hit_count"] for row in sections)
assert hit_histogram == Counter({0: 18, 1: 60, 2: 94, 3: 24, 4: 10})
zero_one = [row for row in sections if row["node_hit_count"] <= 1]
degree_histograms = {
    str(hits): dict(sorted(Counter(
        row["q4o1584_parent_degree"] for row in zero_one if row["node_hit_count"] == hits
    ).items()))
    for hits in (0, 1)
}
zero_node_degree_two = [
    row for row in zero_one
    if row["node_hit_count"] == 0 and row["q4o1584_parent_degree"] == 2
]
assert len(zero_node_degree_two) == 1

payload = {
    "schema": "elkies-k3.q4o164-zero-one-node-parent-degree-audit.v1",
    "status": "PASS_EXACT_Q4O164_ZERO_ONE_NODE_PARENT_DEGREE_AUDIT",
    "section_shell": {
        "signed_norm_four_vectors": len(vectors),
        "physical_P_dot_O_zero_sections": len(sections),
        "node_hit_histogram": dict(sorted(hit_histogram.items())),
        "zero_one_node_sections": len(zero_one),
        "q4o1584_parent_degree_histograms": degree_histograms,
    },
    "unique_zero_node_parent_degree_two": zero_node_degree_two[0],
    "proof_boundary": (
        "This is a complete exact finite lattice enumeration of P.O=0 sections in the "
        "stored physical frame. It identifies a unique zero-node class of parent degree "
        "two; constructing its equation remains separate."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (FRAME, MARKING)],
        "sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in (FRAME, MARKING)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
target = zero_node_degree_two[0]
print(
    "Q4O164LOWNODE|sections={}|hist={}|target_child={}|target_parent={}|status={}|output={}".format(
        len(sections), dict(sorted(hit_histogram.items())), target["class_in_q4o164_basis"],
        target["class_in_q4o1584_basis"], payload["status"], OUTPUT,
    ),
    flush=True,
)
