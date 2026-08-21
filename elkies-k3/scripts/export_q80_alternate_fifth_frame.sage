#!/usr/bin/env sage
"""Export an exact alternate fifth A1/MW16 frame for bounded searches.

The default preserves the first degree-47 witness used by the existing q6
certificate.  ``--v`` and ``--output`` allow the same exact constructor to
export another retained q4 witness, notably the degree-43 orbit.
"""

import argparse
from pathlib import Path

from sage.all import ZZ, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--v",
    default=(
        "-9,8,-11,10,-4,0,5,1,-6,6,1,-2,-1,-1,1,2,0"
    ),
    help="comma-separated norm-eight vector in the generic fourth frame",
)
parser.add_argument(
    "--output",
    type=Path,
    default=(
        ROOT / "artifacts/local/elkies-k3/"
        "q80-alternate-fifth-a1-mw16-frame.txt"
    ),
)
arguments = parser.parse_args()
load(str(HERE / "analyze_q80_fifth_q4_chamber.sage"))

alternate_v = vector(ZZ, [ZZ(value) for value in arguments.v.split(",")])
assert len(alternate_v) == 17
alternate_child, _ = neighbor(
    fourth_child_frame, ZZ(4), ZZ(2), ZZ(2), alternate_v
)
simple, positive = deterministic_simple_roots(alternate_child)
assert simple.nrows() == 1 and len(positive) == 1
output = arguments.output
output.parent.mkdir(parents=True, exist_ok=True)
lines = [
    "# bounded-search input; exact q4 child of generic fourth q80 frame",
    "# q = 4",
    "# a = 2",
    "# b = 2",
    "# v = " + ",".join(map(str, alternate_v)),
]
lines.extend(" ".join(map(str, row)) for row in alternate_child.rows())
output.write_text("\n".join(lines)+"\n")
print(
    "Q80ALTFIFTHEXPORT|"
    f"det={alternate_child.det()}|output={output}|status=PASS",
    flush=True,
)
