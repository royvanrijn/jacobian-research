#!/usr/bin/env python3
"""Finish q8/orbit376 after an exact unpointed QQ reconstruction.

Runs the exact marked-class identification of the two transported q4/o164
one-node sections, the known-section quartic pointing pass, and finally the
preferred P1229-origin certificate.

Invoke with the repository's Sage Python.
"""

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
IDENTIFY = ROOT / "elkies-k3/scripts/identify_h92_q4o164_transported_degree1_marked_classes_qq.sage"
POINT = ROOT / "elkies-k3/scripts/point_h92_q4o164_q8o376_from_known_sections_qq.sage"
CERTIFY = ROOT / "elkies-k3/scripts/certify_h92_q4o164_q8o376_p1229_origin_qq.sage"
RR = LOCAL / "q4o164-q8o376-rr-p2-qq.json"
IDENTIFICATION = LOCAL / "q4o164-transported-degree1-marked-classes-qq.json"
POINTING = LOCAL / "q4o164-q8o376-known-section-pointing-qq.json"
P1229 = LOCAL / "q4o164-q8o376-p1229-origin-qq.json"

python = Path(sys.executable)
if not RR.exists():
    raise SystemExit(
        "missing exact q8 reconstruction; first run run_h92_q4o164_q8o376_rr_p2_stack.py"
    )

subprocess.run([str(python), str(IDENTIFY)], cwd=ROOT, check=True)
identified = json.loads(IDENTIFICATION.read_text())
print(
    "Q8O376FINISH|identify_status={}|p1229={}|output={}".format(
        identified["status"], identified.get("P1229_direct_match_count", 0), IDENTIFICATION
    ),
    flush=True,
)

subprocess.run([str(python), str(POINT)], cwd=ROOT, check=True)
pointing = json.loads(POINTING.read_text())
print(
    "Q8O376FINISH|pointing_status={}|degree1={}|output={}".format(
        pointing["status"], pointing.get("degree_one_count", 0), POINTING
    ),
    flush=True,
)

completed = subprocess.run([str(python), str(CERTIFY)], cwd=ROOT)
if not P1229.exists():
    raise RuntimeError("P1229 certification did not emit its diagnostic artifact")
p1229 = json.loads(P1229.read_text())
print(
    "Q8O376FINISH|p1229_status={}|exit={}|output={}".format(
        p1229["status"], completed.returncode, P1229
    ),
    flush=True,
)
if completed.returncode:
    raise SystemExit(completed.returncode)
