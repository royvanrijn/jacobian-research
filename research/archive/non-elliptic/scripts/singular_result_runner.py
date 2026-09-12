"""Guard legacy Singular RESULT blocks while preserving timeout records.

The historical producers hash a program ending in an explicit quit. Replace
only that terminal command with the shared completion guard. On interpreter
failure return empty stdout, so a caller cannot retain an apparent unit or
rank printed after an error. Timeouts still reach the caller's existing
timeout handler. No mathematical result is inferred here.
"""

import subprocess

from singular_runner import SingularError, run_singular


def run_singular_result(
    source: str, *, executable: str = "Singular", timeout: float = 120,
) -> subprocess.CompletedProcess:
    program = source.rstrip()
    if not program.endswith("quit;"):
        raise ValueError("legacy RESULT program must end with its explicit quit;")
    try:
        output = run_singular(
            program[:-5], executable=executable, timeout=timeout,
            required_markers=("RESULT_BEGIN", "RESULT_END"),
        )
    except SingularError as error:
        return subprocess.CompletedProcess([executable, "-q"], 1, "", str(error))
    return subprocess.CompletedProcess([executable, "-q"], 0, output, "")
