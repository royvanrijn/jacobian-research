#!/usr/bin/env sage -python
"""Run the immutable H3 quotient sieve on the Sage 10.9 source successor.

The original sieve is retained with its frozen source-family digest.  This
wrapper verifies that source, substitutes only the present Sage 10.9 successor
input and successor output identity, then executes the original exact code.
It remains a bounded cross-check, never a global rational-points proof.
"""

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ORIGINAL = ROOT / "elkies-k3/scripts/sieve_h3_level474_rational_points.sage"
ORIGINAL_SHA256 = "831e113eb407f00c9f0d214ec514df661889fdc4539a54f91d4913034c55c473"
SOURCE_SUCCESSOR_SHA256 = "4e4d53f4357f09a553e68926438e6a49ba61f753b1c28b2edccc94071ec16dd5"

source = ORIGINAL.read_text()
assert hashlib.sha256(source.encode()).hexdigest() == ORIGINAL_SHA256
replacements = {
    "elkies-k3-h3-level474-source-family.json": "elkies-k3-h3-level474-source-family-sage109-replay.json",
    "8f5afd11e1d8979d57cb1a569833309f9664c19cd47194af0581a5cbbf8f1d59": SOURCE_SUCCESSOR_SHA256,
    "PASS_EXACT_H3_SOURCE_FAMILY": "PASS_EXACT_H3_SOURCE_FAMILY_SAGE109_REPLAY",
    "elkies-k3-h3-level474-point-sieve.json": "elkies-k3-h3-level474-point-sieve-sage109-replay.json",
    "elkies-k3.h3-level474-point-sieve.v1": "elkies-k3.h3-level474-point-sieve-sage109-replay.v1",
    "PASS_BOUNDED_H3_SOURCE_POINT_SIEVE": "PASS_BOUNDED_H3_SOURCE_POINT_SIEVE_SAGE109_REPLAY",
}
for previous, successor in replacements.items():
    assert previous in source
    source = source.replace(previous, successor)
exec(compile(source, str(ORIGINAL), "exec"))
