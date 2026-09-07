#!/usr/bin/env python3
"""Run the complete backward cubic theorem and application suite."""

from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CHECKERS = (
    "scripts/audit_macfarlane_g20_dimension_reduction.py",
    "scripts/audit_macfarlane_f13_low_degree_invariants.py",
    "scripts/verify_backward_cubic_reduction.py",
    "scripts/verify_backward_cubic_current_applications.py",
    "scripts/verify_rank_compressed_bcw_24_route.py",
    "scripts/verify_restricted_minima_frontier.py",
    "scripts/verify_minimal_counterexample_scoreboard.py",
)
PINNED_FILES = {
    "jcsearch/backward_cubic.py": (
        "598535621f2c63c94338a878c95f94fb92f2a673b2d7ebb92e5c8b87efd8cfa5"
    ),
    "scripts/audit_macfarlane_g20_dimension_reduction.py": (
        "d29ad30d8f4d14890dc21f24c70ced9a4cee2735cb05d8902abe53370518a404"
    ),
    "scripts/audit_macfarlane_f13_low_degree_invariants.py": (
        "47bf01f682b402fbca4b0e3f45b8165f572bb921add26d9fdfaaf19f61ab0666"
    ),
    "scripts/verify_backward_cubic_reduction.py": (
        "12312bd7dd154117de7979d9655fd78d8e6f16a701d14dae93af9a5776052492"
    ),
    "scripts/verify_backward_cubic_current_applications.py": (
        "b040b43ae3e4fe8d4524c3d98efe088b6edbfc4fb2769ac8ae928d1d48379e35"
    ),
    "scripts/rank_compressed_bcw_homogenization.py": (
        "f6863bd99a54945a29ffe5a43d1a14fbf5fb2e977b1cdd11e0387ae243526473"
    ),
    "scripts/search_restricted_bcw_circuits.py": (
        "4cba27d88ce2fd4c1ecc064d91e33ccbfcaf3f4363484c423c82ca1521c02d30"
    ),
    "scripts/verify_rank_compressed_bcw_24_route.py": (
        "f8a69ecbe59d90058592c05d26521557603c23419e87f5450b60aaad9f34f9f6"
    ),
    "scripts/verify_restricted_minima_frontier.py": (
        "29c5e0dd00a7dba8b593c209e693a3af126d19c8dd3394970a892aa1cae8107f"
    ),
    "scripts/verify_minimal_counterexample_scoreboard.py": (
        "d0e2676e0760ba25c0c599b355a82730625a5228c4c33a7a0e8b204817ce7f46"
    ),
}


def main() -> None:
    for filename, expected in PINNED_FILES.items():
        actual = hashlib.sha256((ROOT / filename).read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(
                f"stale backward-cubic suite pin for {filename}: {actual}"
            )
    for checker in CHECKERS:
        subprocess.run(
            [sys.executable, checker],
            cwd=ROOT,
            check=True,
        )
    print(
        "PASS backward cubic suite: theorem, MacFarlane calibration, "
        "restricted-minima/scoreboard applications, and forward BCW regression"
    )


if __name__ == "__main__":
    main()
