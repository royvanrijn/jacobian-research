#!/usr/bin/env python3
"""Fast replay suite for the current binary GVC frontier.

The current proved results and the remaining promotion obstruction are
documented in
``extended-geometry/BINARY_GVC_UNIFORM_FACE_TERMINATION.md``.  This suite
replays five computational interfaces:

1. the uniform Hall/weighted-face reductions;
2. the torsion--torus regular-trace identities;
3. signed base-p digit separation and the Newton endpoint;
4. weighted-trace classification and the affine/factorial obstruction;
5. the bounded translation-twist rigidity search.

The suite is a regression harness.  It does not prove scale-compatible
Hall/carry promotion or unrestricted GVC(2).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
CHECKERS = (
    "verify_binary_gvc_uniform_face_termination.py",
    "research_binary_gvc_torsion_torus_trace.py",
    "verify_binary_gvc_torsion_torus_digit_separation.py",
    "verify_binary_gvc_weighted_trace_obstruction.py",
    "search_binary_gvc_translation_isoperiodic_twists.py",
)


def verify() -> None:
    for checker in CHECKERS:
        subprocess.run(
            [sys.executable, str(SCRIPT_DIRECTORY / checker)],
            check=True,
        )
    print(
        "PASS binary GVC frontier suite: Hall faces, finite trace, "
        "repeated digits, affine/factorial obstruction, and bounded "
        "translation rigidity"
    )
    print(
        "STATUS: scale-compatible Hall/carry promotion remains open"
    )


if __name__ == "__main__":
    verify()
