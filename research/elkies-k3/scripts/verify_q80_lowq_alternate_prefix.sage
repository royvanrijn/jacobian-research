#!/usr/bin/env sage
"""Replay the exact generic Q80 low-q alternate prefix found on 2026-08-22.

The canonical q4,q4 prefix is supplied by analyze_q80_second_neighbor_chamber.
This verifier pins the three new retained moves

    D7+D5/MW5 --q6--> D7+D4/MW6
               --q4--> A6+A4/MW7
               --q4--> A6+A3/MW8.

It is a lattice verifier only; the companion CM24 equation verifier checks the
first two new equation-level specializations.
"""

from pathlib import Path

from sage.all import ZZ, vector

HERE = Path(__file__).resolve().parent
load(str(HERE / "analyze_q80_second_neighbor_chamber.sage"))

assert root_components(second_frame) == [(5, 40, 4), (7, 84, 4)]

escape_v = vector(
    ZZ,
    (-5, -3, 6, 6, -8, -4, 2, 4, -1, 8, -16, -1, 0, 3, 5, -2, -2),
)
q6_frame, q6_transport = neighbor(second_frame, ZZ(6), ZZ(2), ZZ(3), escape_v)
assert q6_frame.det() == second_frame.det() == 948
assert root_components(q6_frame) == [(4, 24, 4), (7, 84, 4)]

orbit424_v = vector(
    ZZ,
    (32, 48, -21, 28, 8, -52, -34, 0, 18, 5, -23, 43, 9, -18, 16, -6, -6),
)
a6a4_frame, orbit424_transport = neighbor(
    q6_frame, ZZ(4), ZZ(2), ZZ(2), orbit424_v
)
assert a6a4_frame.det() == 948
assert root_components(a6a4_frame) == [(4, 20, 5), (6, 42, 7)]

orbit1222_v = vector(
    ZZ,
    (10, 53, -192, -114, 29, -256, -170, -12, -14, 74, -32, -14, -6, -26, -58, 84, -28),
)
a6a3_frame, orbit1222_transport = neighbor(
    a6a4_frame, ZZ(4), ZZ(2), ZZ(2), orbit1222_v
)
assert a6a3_frame.det() == 948
assert root_components(a6a3_frame) == [(3, 12, 4), (6, 42, 7)]

print(
    "Q80LOWQ|prefix="
    "D7+D5/MW5-q6-D7+D4/MW6-q4-A6+A4/MW7-q4-A6+A3/MW8|"
    "det=948|status=PASS_LOWQ_ALTERNATE_PREFIX",
    flush=True,
)
