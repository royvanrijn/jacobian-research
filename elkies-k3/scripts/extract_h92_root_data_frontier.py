#!/usr/bin/env python3
"""Stream a compact exact root-data subfrontier from a large neighbour JSON.

This is a mechanical jq-backed extractor: it preserves all top-level search
metadata and complete selected neighbour records, changing only the neighbour
array.  It is useful when an exhaustive shell artifact is too large for the
downstream Sage ranker to load economically.
"""

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--root-data", required=True, help="comma-separated rank,count,determinant")
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root_data = [int(value) for value in args.root_data.split(",")]
if len(root_data) != 3:
    parser.error("--root-data must contain exactly three integers")

expression = (
    ". as $x | $x + {neighbors: [$x.neighbors[] | "
    + "select(.child_root_data == " + json.dumps(root_data) + ")]}"
)
args.output.parent.mkdir(parents=True, exist_ok=True)
with tempfile.NamedTemporaryFile(
    mode="wb", dir=args.output.parent, prefix=args.output.name + ".", delete=False
) as stream:
    temporary = Path(stream.name)
    result = subprocess.run(
        ["jq", expression, str(args.input)], stdout=stream, stderr=subprocess.PIPE,
        check=False,
    )
if result.returncode:
    temporary.unlink(missing_ok=True)
    raise SystemExit(result.stderr.decode(errors="replace"))
temporary.chmod(0o644)
temporary.replace(args.output)
payload = json.loads(args.output.read_text())
with args.input.open("rb") as source_stream:
    source_sha256 = hashlib.file_digest(source_stream, "sha256").hexdigest()
payload["extraction"] = {
    "method": "exact jq selection by child_root_data",
    "source": str(args.input),
    "source_sha256": source_sha256,
    "selected_root_data": root_data,
    "selected_count": len(payload["neighbors"]),
}
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"ROOTDATAEXTRACT|root_data={','.join(map(str, root_data))}|"
    f"selected={len(payload['neighbors'])}|status=PASS|output={args.output}",
    flush=True,
)
