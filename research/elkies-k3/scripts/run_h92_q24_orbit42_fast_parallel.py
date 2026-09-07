#!/usr/bin/env python3
"""Reproducible rejection gate for the superseded orbit42 fast path.

The old runner tried to materialize O12 and P42 as q6 rational points and
then launch two orientation certifiers. Exact coordinate conversion proves
their q6 degrees are 435 and 703. Run that audit and stop cleanly at the
valid resolved-RR frontier; never launch the invalid transport/orientations.
"""

import argparse
import json
import subprocess
from pathlib import Path


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if (candidate / "elkies-k3/scripts").is_dir():
            return candidate
    raise SystemExit("Could not locate jacobian-research; use --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
args = parser.parse_args()

root = locate_repo(args.repo)
scripts = root / "elkies-k3/scripts"
local = root / "artifacts/local/elkies-k3"
audit_script = scripts / "preflight_h92_q24_o12_p42_exact_q6_points.sage"
artifact = local / "q24-o12-p42-q6-preflight.json"

print("Q42PAR|stage=FAST_Q6_PREMISE_AUDIT|status=BEGIN", flush=True)
subprocess.run(
    ["sage", "-python", str(audit_script)],
    cwd=str(root),
    check=True,
)

data = json.loads(artifact.read_text())
expected = "PASS_Q42_FAST_Q6_PREMISE_REJECTION"
if data.get("status") != expected:
    raise SystemExit(
        "fast-q6 premise audit did not produce the pinned rejection status: "
        + str(data.get("status"))
    )

targets = data["targets"]
print(
    "Q42PAR_RESULT|"
    f"O12_q6_degree={targets['O12']['q6_degree']}|"
    f"P42_q6_degree={targets['P42']['q6_degree']}|"
    "transport=NOT_RUN|orientation0=NOT_RUN|orientation1=NOT_RUN|"
    "next=Q42_RESOLVED_RR_TRIVIALIZATION|"
    "status=STOPPED_INVALID_FAST_Q6_ROUTE",
    flush=True,
)
