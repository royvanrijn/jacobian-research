#!/usr/bin/env sage -python
"""Price exact crossovers from compiler candidates into every certified suffix stage."""

import hashlib
import json
from pathlib import Path

from sage.all import *


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
FINGERPRINT = LOCAL / "q24-a11-q8-construction-fingerprint.json"
O32 = GENERATED / "elkies-k3-h3-a5a5-q4-o32-explicit-zero-frames.json"
O3372 = GENERATED / "elkies-k3-h3-a5a5-q6-o3372-explicit-zero-frames.json"
ZERO_MISMATCH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
OUTPUT = GENERATED / "elkies-k3-h3-candidate-current-suffix-crossovers.json"
INPUTS = (MANIFEST, FINGERPRINT, O32, O3372, ZERO_MISMATCH)
U2 = matrix(ZZ, ((0, 1), (1, 0)))


manifest = json.loads(MANIFEST.read_text())
fingerprint = json.loads(FINGERPRINT.read_text())
o32 = json.loads(O32.read_text())
o3372 = json.loads(O3372.read_text())
zero_mismatch = json.loads(ZERO_MISMATCH.read_text())

historical_in_equation = block_diagonal_matrix(
    identity_matrix(ZZ, 2),
    matrix(ZZ, fingerprint["selected"]["frame_isometry_historical_basis_in_equation_coordinates"]),
)
cumulative = identity_matrix(ZZ, 19)
suffix_fibres = []
for index, step in enumerate(manifest["forward_steps"]):
    if index < 2:
        continue
    cumulative = matrix(ZZ, step["transition"]) * cumulative
    basis_in_equation = cumulative * historical_in_equation
    suffix_fibres.append({
        "stage_index": index,
        "stage": step["child"],
        "historical_incoming_q": int(step["q"]),
        "historical_incoming_orbit": int(step["orbit"]),
        "fibre_in_equation_A11": vector(ZZ, basis_in_equation.row(0)),
    })

a0 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["oldI9_A0"]["child_coordinates"])
p24 = vector(ZZ, zero_mismatch["correct_selected_R3_transport"]["close_P24"]["child_coordinates"])
old_simple = [vector(ZZ, [0, 0] + [-ZZ(index == node) for index in range(17)]) for node in range(11)]
old_affine = vector(ZZ, [1, 0] + [1] * 11 + [0] * 6)
curve_names = [f"old_A11_component_{node}" for node in range(11)] + ["old_A11_affine", "oldI9_A0", "close_P24"]
curves = old_simple + [old_affine, a0, p24]

candidates = {
    "q4_orbit32": o32["selected"],
    "q6_orbit3372": o3372["selected"],
}
records = []
for candidate_name, candidate in candidates.items():
    equation_to_candidate = matrix(ZZ, candidate["equation_A11_to_explicit_zero_basis"])
    candidate_to_equation = equation_to_candidate.inverse().change_ring(ZZ)
    frame = matrix(ZZ, candidate["frame"])
    g_candidate = block_diagonal_matrix(U2, -frame)
    g_equation = candidate_to_equation * g_candidate * candidate_to_equation.transpose()
    lattice = IntegralLattice(frame)
    for target in suffix_fibres:
        fibre_equation = target["fibre_in_equation_A11"]
        fibre = fibre_equation * candidate_to_equation
        assert fibre * g_candidate * fibre == 0
        degree = ZZ(fibre[1])
        q = ZZ(fibre[0] * fibre[1])
        center = vector(QQ, fibre[2:]) / degree
        closest = vector(ZZ, next(lattice.enumerate_close_vectors(center)))
        distance = (closest - center) * frame * (closest - center)
        minimum_section_intersection = QQ(degree) * distance / 2 - degree
        explicit_degrees = [int(curve * g_equation * fibre_equation) for curve in curves]
        records.append({
            "candidate": candidate_name,
            "target_stage_index": target["stage_index"],
            "target_stage": target["stage"],
            "historical_incoming_q": target["historical_incoming_q"],
            "historical_incoming_orbit": target["historical_incoming_orbit"],
            "target_fibre_in_candidate": [int(value) for value in fibre],
            "old_fibre_degree": int(degree),
            "q_in_candidate_zero": int(q),
            "target_zero_intersection": int(fibre[0] - fibre[1]),
            "minimum_section_intersection": str(minimum_section_intersection),
            "explicit_curve_degrees": dict(zip(curve_names, explicit_degrees)),
            "explicit_degree_zero_curves": [curve_names[index] for index, value in enumerate(explicit_degrees) if value == 0],
            "explicit_degree_one_curves": [curve_names[index] for index, value in enumerate(explicit_degrees) if value == 1],
            "negative_explicit_curves": [curve_names[index] for index, value in enumerate(explicit_degrees) if value < 0],
        })

payload = {
    "schema": "elkies-k3.h3-candidate-current-suffix-crossovers.v1",
    "status": "PASS_EXACT_CANDIDATE_CURRENT_SUFFIX_CROSSOVER_AUDIT",
    "inputs": {"paths": [str(path.relative_to(ROOT)) for path in INPUTS], "sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in INPUTS}},
    "records": records,
    "best_by_candidate": {
        name: sorted([item for item in records if item["candidate"] == name], key=lambda item: (item["old_fibre_degree"], item["q_in_candidate_zero"]))[:5]
        for name in candidates
    },
    "proof_boundary": "Exact marked class transports and direct crossover costs. A crossover is not promoted until its fibre is separately certified nef in the candidate chamber.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
for name in candidates:
    best = payload["best_by_candidate"][name][0]
    print("SUFFIXCROSS|candidate={}|target={}|degree={}|q={}|Pmin={}|deg0={}|deg1={}".format(
        name, best["target_stage"], best["old_fibre_degree"], best["q_in_candidate_zero"], best["minimum_section_intersection"],
        len(best["explicit_degree_zero_curves"]), len(best["explicit_degree_one_curves"])), flush=True)
print(f"SUFFIXCROSS|status={payload['status']}|output={OUTPUT.resolve()}", flush=True)
