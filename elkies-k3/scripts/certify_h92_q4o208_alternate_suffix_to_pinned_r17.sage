#!/usr/bin/env sage -python
"""Close the physical q4o208 alternate suffix at pinned rootless R17.

This is a composition/endpoint checker.  Each edge is independently certified
by the generic marked-neighbour checker; here we pin their hashes, verify both
integral NS transports, and identify the final rootless positive frame with
``rank17_gram.txt`` by an explicitly checked unimodular isometry.
"""
import hashlib, json
from pathlib import Path
from sage.all import ZZ, identity_matrix, matrix, pari

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
OUTPUT = GEN / "elkies-k3-h3-q4o208-alternate-q4q8q12q4q8q8-pinned-r17-route-certificate.json"

EDGES = (
    GEN / "elkies-k3-h3-q4o208-physical-q4o323-corrected-a3-2a2-certificate.json",
    GEN / "elkies-k3-h3-q4o208-corrected-a3-2a2-to-5a1-physical-q8-c10-certificate.json",
    GEN / "elkies-k3-h3-physical-5a1-to-4a1-q12-certificate.json",
    GEN / "elkies-k3-h3-physical-4a1-q4o21633-a2-certificate.json",
    GEN / "elkies-k3-h3-physical-a2-q8o2102-a1-certificate.json",
    GEN / "elkies-k3-h3-physical-a2-q8o2102-a1-q8o5165-rootless-certificate.json",
)

def load(path):
    return matrix(ZZ, [[ZZ(x) for x in line.split()] for line in path.read_text().splitlines()
                       if line.strip() and not line.lstrip().startswith("#")])

def rows(value):
    return [[int(x) for x in row] for row in value.rows()]

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def qfisom_row(source, target):
    raw = pari(source).qfisom(pari(target))
    assert str(raw) != "0"
    value = matrix(ZZ, raw)
    candidates = [value, value.transpose()]
    inverse = value.inverse()
    if inverse.change_ring(ZZ) == inverse:
        inverse = inverse.change_ring(ZZ)
        candidates += [inverse, inverse.transpose()]
    for candidate in candidates:
        if candidate * source * candidate.transpose() == target:
            assert abs(candidate.det()) == 1
            return candidate
    raise ArithmeticError("PARI qfisom returned no verified row-convention isometry")

payloads = [json.loads(path.read_text()) for path in EDGES]
records = []
for path, payload in zip(EDGES, payloads):
    assert payload["status"] == "PASS_EXACT_MARKED_DEGREE_TWO_CANDIDATE_CERTIFICATE"
    forward = matrix(ZZ, payload["source_to_child_basis"])
    backward = matrix(ZZ, payload["child_to_source_basis"])
    equation_forward = matrix(ZZ, payload["equation_A11_to_child_basis"])
    equation_backward = matrix(ZZ, payload["child_to_equation_A11_basis"])
    assert forward * backward == backward * forward == identity_matrix(ZZ, 19)
    assert equation_forward * equation_backward == equation_backward * equation_forward == identity_matrix(ZZ, 19)
    assert abs(forward.det()) == abs(equation_forward.det()) == 1
    assert payload["first_edge_exact_horizontal_nef_gate"]
    assert not payload["first_edge_exact_negative_horizontal_walls"]
    records.append({
        "certificate": str(path.relative_to(ROOT)),
        "sha256": sha(path),
        "candidate_id": payload["candidate_id"],
        "source_hub": payload["source_hub"],
        "child_root_data": payload["child"]["root_data"],
        "child_mw_rank": payload["child"]["mw_rank"],
        "source_to_child_determinant": int(forward.det()),
        "child_to_source_determinant": int(backward.det()),
    })

final = payloads[-1]
assert final["child"]["root_data"] == [0, 0, 1] and final["child"]["mw_rank"] == 17
assert final["target_profiles"]["pinned_R17"]["nef_audit"]["nef_in_selected_component_chamber"]
child_path = ROOT / final["frame_output"]
child = load(child_path)
pinned = load(PINNED)
assert child.det() == pinned.det() == 948
isometry = qfisom_row(child, pinned)

result = {
    "schema": "elkies-k3.h3-q4o208-alternate-pinned-r17-route-certificate.v1",
    "status": "PASS_EXACT_Q4O208_ALTERNATE_SUFFIX_TO_PINNED_R17",
    "source": "physical q4/o208 3A3 with equation-effective C5 zero",
    "route": records,
    "q_sequence": [4, 8, 12, 4, 8, 8],
    "old_fibre_degrees": [2] * 6,
    "ade_mw_sequence": [
        "A3+2A2/MW10", "5A1/MW12", "4A1/MW13",
        "A2/MW15", "A1/MW16", "rootless/MW17",
    ],
    "final_child_frame": str(child_path.relative_to(ROOT)),
    "final_child_frame_sha256": sha(child_path),
    "pinned_rank17_frame": str(PINNED.relative_to(ROOT)),
    "pinned_rank17_frame_sha256": sha(PINNED),
    "rootless_child_to_pinned_r17_isometry": rows(isometry),
    "isometry_relation": "Q * child * Q^t = pinned rank17_gram.txt",
    "isometry_determinant": int(isometry.det()),
    "equation_A11_to_rootless_child_basis": final["equation_A11_to_child_basis"],
    "rootless_child_to_equation_A11_basis": final["child_to_equation_A11_basis"],
    "proof_boundary": (
        "Every listed edge has an exact marked U, primitive nef isotropic fibre, "
        "complete finite horizontal-wall audit, full roots, and mutually inverse "
        "unimodular NS transports.  The terminal child is rootless and is identified "
        "with pinned R17 by the displayed exact integral isometry.  This proves the "
        "lattice route; it does not assert that the three new late edges have equations."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in EDGES + (PINNED,)],
        "sha256": {str(path.relative_to(ROOT)): sha(path) for path in EDGES + (PINNED,)},
    },
}
OUTPUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print(f"H92ALTROOTLESS|edges={len(EDGES)}|q=4,8,12,4,8,8|det={isometry.det()}|status={result['status']}|output={OUTPUT}")
