#!/usr/bin/env sage-python
"""Run the complete NS0024 edge-1 handoff from either supported source format.

status: ACTIVE_COMPILER
claim: deterministic dispatch from a certified source family or compact marked point

For a certified family, this invokes the edge compiler directly.  For a compact
quadratic-extension P4 point, it first runs the independent source-marking
adapter and then invokes the same compiler.  Every subprocess is fail-closed;
no output status is synthesized by this runner.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ADAPTER = HERE / "adapt_lattice_foundry_ns0024_mw4_point_for_edge1.sage"
COMPILER = HERE / "compile_lattice_foundry_ns0024_edge1_modp.sage"
FAMILY_SCHEMA = "elkies-k3.lattice-foundry-ns0024-mw4-family-modp.v1"
POINT_SCHEMA = "elkies-k3.lattice-foundry-ns0024-mw4-point-modp.v1"


def display_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def run(command):
    subprocess.run([str(item) for item in command], cwd=ROOT, check=True)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--seed", type=Path, help="optional MW3 seed override for a compact point")
parser.add_argument(
    "--source-output",
    type=Path,
    help="required adapted-source output when --input is a compact point",
)
parser.add_argument("--edge-output", type=Path, required=True)
parser.add_argument("--check", action="store_true")
args = parser.parse_args()

input_path = args.input.resolve()
edge_output = args.edge_output.resolve()
record = json.loads(input_path.read_text())
schema = record.get("schema")

if schema == FAMILY_SCHEMA:
    if args.seed is not None or args.source_output is not None:
        raise SystemExit("--seed/--source-output apply only to a compact point input")
    source_path = input_path
    mode = "family"
elif schema == POINT_SCHEMA:
    if args.source_output is None:
        raise SystemExit("compact point input requires --source-output")
    source_path = args.source_output.resolve()
    adapter_command = [
        sys.executable,
        ADAPTER,
        "--point",
        input_path,
        "--output",
        source_path,
    ]
    if args.seed is not None:
        adapter_command.extend(("--seed", args.seed.resolve()))
    if args.check:
        adapter_command.append("--check")
    run(adapter_command)
    mode = "point"
else:
    raise SystemExit("input schema is neither a certified NS0024 family nor compact marked point")

compiler_command = [
    sys.executable,
    COMPILER,
    "--input",
    source_path,
    "--output",
    edge_output,
]
if args.check:
    compiler_command.append("--check")
run(compiler_command)

print(
    "NS0024EDGE1HANDOFF|mode={}|source={}|edge={}|check={}|status=PASS".format(
        mode,
        display_path(source_path),
        display_path(edge_output),
        int(args.check),
    ),
    flush=True,
)
