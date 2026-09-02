#!/usr/bin/env python3
"""Run the all-Niemeier catalogue, T-arithmetic, and Pareto stages in order.

This is the supported factory entry point.  It intentionally stops before any
equation solver: equation targets are emitted separately, and their adapter
requires the hash-matched T-arithmetic ledger produced here.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CATALOGUE = ROOT / "elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage"
T_ARITHMETIC = ROOT / "elkies-k3/scripts/build_rank7_t_arithmetic.sage"
PARETO = ROOT / "elkies-k3/scripts/build_rank7_surface_pareto.py"


def run(command: list[str]) -> None:
    print("RANK7FACTORY|run=" + " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--cm-disc-max", type=int, default=200)
    arguments = parser.parse_args()

    sage = shutil.which("sage")
    if sage is None:
        raise SystemExit("sage executable is required for the rank-seven factory")
    check = ["--check"] if arguments.check else []
    run([sage, "-python", str(CATALOGUE), *check])
    run(
        [
            sage,
            "-python",
            str(T_ARITHMETIC),
            "--cm-disc-max",
            str(arguments.cm_disc_max),
            *check,
        ]
    )
    run(["python3", str(PARETO), *check])
    print(
        "RANK7FACTORY|catalogue=PASS|t_arithmetic=PASS|pareto=PASS|"
        "equation_solvers=NOT_LAUNCHED|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
