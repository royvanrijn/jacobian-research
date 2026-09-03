#!/usr/bin/env sage-python
"""Replay the complete pinned-ICARM sweep and the curve-12 quotient."""

from __future__ import annotations

from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = (
    ROOT / "elkies-k3/scripts/certify_r17_norm12_icarm_database_sweep.sage",
    ROOT / "elkies-k3/scripts/certify_r17_norm12_curve12_alternate_q80_quotient.sage",
)


def main() -> None:
    for script in SCRIPTS:
        subprocess.run(
            ["sage", "-python", str(script), "--check"],
            cwd=ROOT,
            check=True,
        )
    print(
        "R17ICARMSWEEPMASTER|curves=474|j_classes=6|curve12_quotient=Z^12|"
        "status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
