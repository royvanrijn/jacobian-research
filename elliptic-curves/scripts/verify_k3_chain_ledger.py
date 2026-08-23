#!/usr/bin/env python3
"""Audit the canonical H3 q6/q8 identities and the later lattice chain."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def digest(relative: str) -> str:
    return hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()


EXPECTED_FILES = {
    "elkies-k3/scripts/certify_h3_q6_actual_neighbor_hop.sage":
        "0bb27f766dc0dcc35f73c7bed2cdfb887b4fe8c58b0c5258b64c1fd0229e37e2",
    "elkies-k3/scripts/derive_h92_q6_child_q8_marking.sage":
        "fe1cf5d14bc91c93d7c4d9fbb1c3126869abdb01e4aeab7b1f938df04ef7d174",
    "elkies-k3/scripts/derive_h92_q6_child_q8_corrected2cover_qq.sage":
        "01da6b84ee1fa1e51aab4b194d2fb9224aaf733e0f054841392256bb6f73a2b6",
    "elkies-k3/data/fibrations/h3_q8_component_nef_physical_root_target.json":
        "569667c6ad6dc6a8fb82bd22ac9444e208e5ee1bb50788a0917ef1c5a80fdc46",
    "elkies-k3/scripts/verify_h3_d13_to_mw17_path.sage":
        "3a3a1b51553c23373e6889af9e129f5d3667663e127ef03c9cedb026f7ee1bdb",
    "artifacts/generated-results/elkies-k3-h3-d13-to-mw17-path.json":
        "f6eac2339c86de84b79a0ddfec3229df9b9c1617110bdd9c474443e7e39fd484",
}

for path, expected in EXPECTED_FILES.items():
    actual = digest(path)
    assert actual == expected, f"stale K3 chain input {path}: {actual}"

status = json.loads((ROOT / "MATH_STATUS.json").read_text())
entries = {entry["id"]: entry for entry in status["entries"]}

q6 = entries["EC-K3-H3-Q6"]
assert q6["state"] == "proved"
assert q6["checker"] == "elkies-k3/scripts/certify_h3_q6_actual_neighbor_hop.sage"
assert "E8+E6" in q6["scope"] and "rank three" in q6["scope"]

marking = entries["EC-K3-H3-Q8-CHILD-MARKING"]
assert "degree-ten" in marking["scope"]
assert "height 24" in marking["scope"]
assert "degree-46" not in marking["scope"]

q8 = entries["EC-K3-H3-Q8-QQ-D13"]
assert q8["state"] == "proved"
assert q8["checker"] == "elkies-k3/scripts/derive_h92_q6_child_q8_corrected2cover_qq.sage"
assert "including the Dx factor" in q8["scope"]
assert "D13" in q8["scope"] and "rank is four" in q8["scope"]

historical_q8_descendants = {
    "CHILD-IVSTAR-VERTICAL-IDEAL-PAIR",
    "CHILD-IVSTAR-ORIENTATION",
    "CHILD-ADDITIVE-CHORD-BLOCKS",
    "CHILD-SATURATED-ANSATZ-PROBE",
    "CHILD-SATURATED-PENCIL-OBSTRUCTION",
    "CHILD-FINITE-MODULE-MODULAR",
    "CHILD-FINITE-Q-MODULE-QQ",
    "CHILD-Q-FRAME-POLE-PROFILE",
    "CHILD-Q-FRAME-NORMALIZATION-MODULAR",
    "CHILD-Q-FRAME-NORMALIZATION-CRT",
    "CHILD-Q-REGULAR-FRAME",
    "CHILD-Q-REGULAR-FINITE-MODULE",
    "CHILD-Q-REGULAR-SMOOTH-FRAME",
    "CHILD-Q-REGULAR-GENERATOR-OBSTRUCTION",
    "CHILD-DIAGONAL-PENCIL-OBSTRUCTION",
    "CHILD-SMOOTH-MODULE",
}
assert all(
    entries[f"EC-K3-H3-Q8-{entry_id}"]["state"] == "archived"
    for entry_id in historical_q8_descendants
)

chain = entries["EC-K3-H3-D13-MW17-LATTICE-CHAIN"]
assert chain["state"] == "proved"
assert chain["checker"] == "elkies-k3/scripts/verify_h3_d13_to_mw17_path.sage"
assert "eleven" in chain["scope"] and "equation-level" in chain["scope"]

replay = json.loads(
    (ROOT / "artifacts/generated-results/elkies-k3-h3-d13-to-mw17-path.json").read_text()
)
assert replay["status"] == "PASS_H3_D13_TO_MW17_LATTICE_PATH"
assert replay["source_ade"] == "D13" and replay["source_mw_rank"] == 4
assert replay["final_ade"] == "rootless" and replay["final_mw_rank"] == 17
assert len(replay["steps"]) == 11
assert [step["q"] for step in replay["steps"]] == [24, 6, 8, 4, 4, 4, 4, 4, 4, 4, 6]
assert all(step["factor_order"][1] == 2 for step in replay["steps"])
for step in replay["steps"]:
    assert digest(step["artifact"]) == step["artifact_sha256"]

corrected = (ROOT / q8["checker"]).read_text()
assert "R*h*Dy == Ny*Dx mod Nx" in corrected
assert "h3_q8_component_nef_physical_root_target.json" in corrected

print(
    "PASS K3 chain ledger: canonical q6, corrected q8, pinned Dx target, "
    "and exact eleven-step D13-to-rootless lattice chain"
)
