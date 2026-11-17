#!/usr/bin/env python3
"""Run a Singular program and reject interpreter errors and incomplete output.

Singular can print ``?`` diagnostics, continue to a PASS marker and exit zero.
An exit code or a later marker alone therefore cannot certify a replay. This
runner checks both output streams and requires its appended completion marker.
Programs passed here must return to the interpreter; do not terminate with quit.
This validates execution, not the mathematical sufficiency of a program's checks.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import subprocess
import sys
from uuid import uuid4


class SingularError(RuntimeError):
    """The interpreter failed or did not finish the submitted program."""


def run_singular(
    source: str,
    *,
    timeout: float = 120,
    executable: str = "Singular",
    required_markers: tuple[str, ...] = (),
) -> str:
    """Return stdout only after error and completion checks; preserve diagnostics."""
    marker = "SINGULAR_RUN_COMPLETE_" + uuid4().hex
    completed = subprocess.run(
        [executable, "-q"],
        input=source + f'\nprint("{marker}");\nquit;\n',
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    lines = completed.stdout.splitlines()
    nonempty = [line.strip() for line in lines if line.strip()]
    errors = [line for line in lines if line.lstrip().startswith("?")]
    reason = None
    if completed.returncode != 0:
        reason = f"exit status {completed.returncode}"
    elif errors or completed.stderr.strip():
        reason = "interpreter diagnostics (even if the exit status is zero)"
    elif nonempty.count(marker) != 1 or nonempty[-1:] != [marker]:
        reason = "missing or misplaced completion marker"
    else:
        missing = [item for item in required_markers if nonempty.count(item) != 1]
        if missing:
            reason = "missing or repeated result markers: " + ", ".join(missing)
    if reason:
        # Include the first errors as well as the tail: a long census may print
        # many successful records after the diagnostic that invalidated it.
        diagnostics = "\n".join(errors[:10])
        raise SingularError(
            f"Singular rejected: {reason}\n{diagnostics}\n"
            + completed.stdout[-4000:]
            + completed.stderr[-4000:]
        )
    return "".join(
        line for line in completed.stdout.splitlines(keepends=True)
        if line.strip() != marker
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Singular source file without a quit command")
    parser.add_argument("--timeout", type=float, default=600, help="wall-clock limit in seconds")
    parser.add_argument("--executable", default="Singular")
    parser.add_argument("--require-marker", action="append", default=[])
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    try:
        output = run_singular(
            args.source.read_text(), timeout=args.timeout, executable=args.executable,
            required_markers=tuple(args.require_marker),
        )
    except (SingularError, subprocess.TimeoutExpired, OSError) as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(output, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
