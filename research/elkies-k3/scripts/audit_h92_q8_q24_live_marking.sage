#!/usr/bin/env sage -python
"""
Minimal live-state audit for the H92 q8/q24 marking contradiction.

Executes the CURRENT local q8 equation certifier and CURRENT q24 effective-zero
producer in isolated namespaces in the SAME Sage process.  No transport is
attempted until the claimed common NS marking is verified entry-for-entry.

This is deliberately diagnostic and does not modify tracked files.
"""

import hashlib
import json
import sys
from pathlib import Path

from sage.all import QQ, ZZ, matrix, vector


def locate_repo():
    cwd = Path.cwd().resolve()
    candidates = [cwd, *cwd.parents]
    h = Path.home()
    candidates += [
        h / "Documents" / "jacobian-research",
        h / "jacobian-research",
        h / "src" / "jacobian-research",
        h / "git" / "jacobian-research",
        h / "projects" / "jacobian-research",
    ]
    seen = set()
    for c in candidates:
        try:
            c = c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (
            (c / "elkies-k3/scripts").is_dir()
            and (c / "artifacts/generated-results").is_dir()
        ):
            return c
    raise SystemExit("Could not locate jacobian-research")


ROOT = locate_repo()
SCRIPTS = ROOT / "elkies-k3/scripts"
LOCAL = ROOT / "artifacts/local/elkies-k3"
EQ = SCRIPTS / "certify_h92_q8_equation_ns_divisor.sage"
Q24 = SCRIPTS / "audit_h92_q8_q24_effective_zero_choices.sage"
OUT = LOCAL / "q8-q24-live-marking-audit.json"

for p in (EQ, Q24):
    if not p.exists():
        raise SystemExit(f"missing {p}")


def run_script(path, argv):
    saved = list(sys.argv)
    scope = {"__name__": "__embedded__"}
    try:
        sys.argv = [str(path)] + list(argv)
        exec(compile(path.read_text(), str(path), "exec"), scope)
    finally:
        sys.argv = saved
    return scope


print("Q8Q24LIVE|stage=equation_certifier", flush=True)
eq = run_script(EQ, [])

q24_tmp = LOCAL / "q8-q24-live-producer-output.json"
print("Q8Q24LIVE|stage=q24_producer", flush=True)
qp = run_script(Q24, ["--output", str(q24_tmp)])

required_eq = ("ns", "physical", "F8eq", "target")
required_qp = ("source_ns", "F_actual", "D_actual", "target_F", "target")
for label, scope, names in (
    ("eq", eq, required_eq),
    ("q24", qp, required_qp),
):
    missing = [n for n in names if n not in scope]
    if missing:
        raise SystemExit(f"{label} missing live variables: {missing}")

ns_eq = matrix(ZZ, eq["ns"])
ns_qp = matrix(ZZ, qp["source_ns"])
F_eq_physical = vector(ZZ, eq["physical"])
F_qp = vector(ZZ, qp["F_actual"])
F_qp_target = vector(ZZ, qp["target_F"])
D_qp = vector(ZZ, qp["D_actual"])

target_eq = vector(
    ZZ, eq["target"]["selected_q8"]["source_h3_ns_vector"]
)
target_qp = vector(
    ZZ, qp["target"]["selected_q8"]["source_h3_ns_vector"]
)

ns_equal = ns_eq == ns_qp
f_live_equal = F_eq_physical == F_qp
targets_equal = target_eq == target_qp
qp_target_equal = F_qp == F_qp_target

deg_qp = ZZ(D_qp * ns_qp * F_qp)
deg_cross = ZZ(D_qp * ns_eq * F_eq_physical)
sq_qp = ZZ(D_qp * ns_qp * D_qp)
sq_cross = ZZ(D_qp * ns_eq * D_qp)

print(
    "Q8Q24LIVE_MARKING|"
    f"ns_equal={int(ns_equal)}|"
    f"targets_equal={int(targets_equal)}|"
    f"F_live_equal={int(f_live_equal)}|"
    f"q24_F_equals_q24_target={int(qp_target_equal)}|"
    f"q24_degree_native={deg_qp}|q24_degree_cross={deg_cross}|"
    f"q24_square_native={sq_qp}|q24_square_cross={sq_cross}|"
    "status=PASS_DIAGNOSTIC",
    flush=True,
)

# Compare exact target source paths, if exposed.
eq_target_path = str(eq.get("TARGET", "NA"))
qp_target_path = str(qp.get("TARGET", "NA"))
eq_frame_path = str(eq.get("FRAME", "NA"))
qp_frame_path = str(qp.get("FRAME", "NA"))
print(
    "Q8Q24LIVE_PATHS|"
    f"eq_TARGET={eq_target_path}|q24_TARGET={qp_target_path}|"
    f"eq_FRAME={eq_frame_path}|q24_FRAME={qp_frame_path}",
    flush=True,
)

# Read back exactly what the live q24 producer just wrote.
written = json.loads(q24_tmp.read_text())
D_written = vector(
    ZZ, written["transport"]["q24_divisor_source_h3_ns"]
)
written_equal_live = D_written == D_qp
written_deg_qp = ZZ(D_written * ns_qp * F_qp)
written_deg_eq = ZZ(D_written * ns_eq * F_eq_physical)

print(
    "Q8Q24LIVE_SERIALIZATION|"
    f"D_written_equals_live={int(written_equal_live)}|"
    f"written_degree_native={written_deg_qp}|"
    f"written_degree_cross={written_deg_eq}|"
    "status=PASS_DIAGNOSTIC",
    flush=True,
)

def vecdiff(a, b):
    d = vector(ZZ, a) - vector(ZZ, b)
    return {
        "equal": bool(not d),
        "nonzero": [
            [int(i), int(v)]
            for i, v in enumerate(d)
            if v
        ],
        "square_eq_ns": int(d * ns_eq * d),
    }

ns_diff_entries = []
if not ns_equal:
    for i in range(ns_eq.nrows()):
        for j in range(ns_eq.ncols()):
            if ns_eq[i,j] != ns_qp[i,j]:
                ns_diff_entries.append(
                    [i, j, int(ns_eq[i,j]), int(ns_qp[i,j])]
                )

payload = {
    "schema": "elkies-k3.h92-q8-q24-live-marking-audit.v1",
    "status": "PASS_LIVE_MARKING_AUDIT",
    "paths": {
        "equation_script": str(EQ.relative_to(ROOT)),
        "q24_script": str(Q24.relative_to(ROOT)),
        "equation_TARGET": eq_target_path,
        "q24_TARGET": qp_target_path,
        "equation_FRAME": eq_frame_path,
        "q24_FRAME": qp_frame_path,
    },
    "matrix": {
        "ns_equal": ns_equal,
        "ns_diff_entries": ns_diff_entries,
    },
    "fibres": {
        "target_vectors_equal": targets_equal,
        "live_fibres_equal": f_live_equal,
        "q24_live_equals_q24_target": qp_target_equal,
        "eq_target_vs_q24_target": vecdiff(target_eq, target_qp),
        "eq_physical_vs_q24_live": vecdiff(F_eq_physical, F_qp),
    },
    "q24_divisor": {
        "live_vector": list(map(int, D_qp)),
        "degree_native": int(deg_qp),
        "degree_cross": int(deg_cross),
        "square_native": int(sq_qp),
        "square_cross": int(sq_cross),
        "written_equals_live": written_equal_live,
        "written_degree_native": int(written_deg_qp),
        "written_degree_cross": int(written_deg_eq),
    },
}

OUT.parent.mkdir(parents=True, exist_ok=True)
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)

if (
    ns_equal
    and targets_equal
    and f_live_equal
    and qp_target_equal
    and deg_qp == deg_cross == 2
    and sq_qp == sq_cross == 0
    and written_equal_live
    and written_deg_qp == written_deg_eq == 2
):
    print(
        "Q8Q24LIVE_RESULT|"
        "common_marking=1|live_and_serialized_degree=2|"
        "status=PASS_COMMON_SOURCE_H3_MARKING",
        flush=True,
    )
else:
    print(
        "Q8Q24LIVE_RESULT|"
        f"common_marking={int(ns_equal and targets_equal and f_live_equal)}|"
        f"native_degree={deg_qp}|cross_degree={deg_cross}|"
        f"serialized_equal={int(written_equal_live)}|"
        "status=MARKING_OR_RUNTIME_MISMATCH_LOCATED",
        flush=True,
    )
