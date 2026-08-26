#!/usr/bin/env python3
"""Finish q8/orbit376 after exact unpointed QQ reconstruction."""

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
IDENTIFY = ROOT / "elkies-k3/scripts/identify_h92_q4o164_transported_degree1_marked_classes_qq.sage"
POINT = ROOT / "elkies-k3/scripts/point_h92_q4o164_q8o376_from_known_sections_qq.sage"
CERTIFY = ROOT / "elkies-k3/scripts/certify_h92_q4o164_q8o376_p1229_origin_qq_v2.sage"
RR = LOCAL / "q4o164-q8o376-rr-p2-qq.json"
IDENTIFICATION = LOCAL / "q4o164-transported-degree1-marked-classes-qq.json"
POINTING = LOCAL / "q4o164-q8o376-known-section-pointing-qq.json"
P1229 = LOCAL / "q4o164-q8o376-p1229-origin-qq.json"
python = Path(sys.executable)

if not RR.exists():
    raise SystemExit("missing q4o164-q8o376-rr-p2-qq.json; run the RR stack first")

for stage, script in (("identify", IDENTIFY), ("point", POINT)):
    subprocess.run([str(python), str(script)], cwd=ROOT, check=True)
    path = IDENTIFICATION if stage == "identify" else POINTING
    data = json.loads(path.read_text())
    print(f"Q8O376FINISH|stage={stage}|status={data['status']}|output={path}", flush=True)

completed = subprocess.run([str(python), str(CERTIFY)], cwd=ROOT)
if not P1229.exists():
    raise RuntimeError("preferred-origin certificate emitted no artifact")
data = json.loads(P1229.read_text())
print(
    f"Q8O376FINISH|stage=p1229|status={data['status']}|exit={completed.returncode}|output={P1229}",
    flush=True,
)
if completed.returncode:
    raise SystemExit(completed.returncode)
