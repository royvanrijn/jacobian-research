#!/usr/bin/env python3
"""Fast replay suite for binary GVC and the parked Hall/carry route.

The unrestricted proof is in
``extended-geometry/BINARY_GVC_ENVELOPE_CLOSURE.md``.  Its new step is a
written finite-envelope argument.  This suite replays the computational
interfaces which led to it and the still-valid Hall/carry route:

1. the first exact Ferrers face exposed by the envelope argument;
2. the uniform Hall/weighted-face reductions;
3. the torsion--torus regular-trace identities;
4. signed base-p digit separation and the Newton endpoint;
5. weighted-trace classification and the affine/factorial obstruction;
6. the proved primitive translation-tangent and large-prime reduction;
7. the bounded translation-twist rigidity search;
8. the all-scale affine-ray Wronskian, four-state, common-base, and affine
   factorial boundary-transfer and periodic additive structural certificate;
9. the exact prime-specific Kummer obstruction to common Cobham promotion.

The suite is a regression harness.  The parked scale-compatible Hall/carry
promotion statement is not proved, but it is no longer needed for GVC(2).
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPT_DIRECTORY = Path(__file__).resolve().parent
CHECKERS = (
    "verify_binary_gvc_ferrers_regression.py",
    "verify_binary_gvc_uniform_face_termination.py",
    "research_binary_gvc_torsion_torus_trace.py",
    "verify_binary_gvc_torsion_torus_digit_separation.py",
    "verify_binary_gvc_weighted_trace_obstruction.py",
    "verify_binary_gvc_translation_tangent_rigidity.py",
    "search_binary_gvc_translation_isoperiodic_twists.py",
    "verify_binary_gvc_cobham_carry_obstruction.py",
)


def verify() -> None:
    for checker in CHECKERS:
        subprocess.run(
            [sys.executable, str(SCRIPT_DIRECTORY / checker)],
            check=True,
        )
    subprocess.run(
        [
            sys.executable,
            str(
                SCRIPT_DIRECTORY
                / "research_binary_gvc_all_scale_orbit_circuits.py"
            ),
            "--structural-certificate",
            "--maximum-wronskian-rank",
            "7",
        ],
        check=True,
    )
    print(
        "PASS binary GVC regression suite: Ferrers face, Hall faces, finite trace, "
        "repeated digits, affine/factorial obstruction, proved tangent "
        "rigidity, bounded translation search, affine-ray Wronskians, "
        "four-state pairing, affine-factorial boundary transfers, and "
        "periodic additive elimination, plus the Cobham carry obstruction"
    )
    print(
        "STATUS: unrestricted GVC(2) is proved by Hall-envelope separation; "
        "the prime-specific Hall/carry route is parked"
    )


if __name__ == "__main__":
    verify()
