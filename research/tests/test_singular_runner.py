"""Interpreter success must not be inferred from a zero exit code or PASS text."""

from pathlib import Path
import re
import subprocess
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from singular_runner import SingularError, run_singular


class SingularRunnerTests(unittest.TestCase):
    def execute(self, stdout="PASS\n", stderr="", returncode=0, finish=True, **kwargs):
        def process(command, **options):
            marker = re.search(r'print\("(SINGULAR_RUN_COMPLETE_[a-f0-9]+)"\)', options["input"])[1]
            output = stdout + (marker + "\n" if finish else "")
            return subprocess.CompletedProcess(command, returncode, output, stderr)
        with patch("singular_runner.subprocess.run", side_effect=process):
            return run_singular('print("PASS");', **kwargs)

    def test_error_before_later_pass_with_zero_exit_is_rejected(self):
        with self.assertRaisesRegex(SingularError, "intentional error"):
            self.execute(stdout="   ? intentional error\nPASS\n")

    def test_diagnostics_on_stderr_are_rejected(self):
        with self.assertRaisesRegex(SingularError, "diagnostics"):
            self.execute(stderr="library loading failed\n")

    def test_clean_output_is_returned_without_internal_marker(self):
        self.assertEqual(self.execute(required_markers=("PASS",)), "PASS\n")

    def test_early_quit_or_truncated_chunk_is_rejected(self):
        with self.assertRaisesRegex(SingularError, "completion marker"):
            self.execute(finish=False)

    def test_missing_repeated_or_substring_result_marker_is_rejected(self):
        for output in ("", "PASS\nPASS\n", "NOT_PASS\n"):
            with self.subTest(output=output):
                with self.assertRaisesRegex(SingularError, "result markers"):
                    self.execute(stdout=output, required_markers=("PASS",))

    def test_nonzero_exit_and_timeout_are_not_accepted(self):
        with self.assertRaisesRegex(SingularError, "exit status 7"):
            self.execute(returncode=7)
        with patch("singular_runner.subprocess.run", side_effect=subprocess.TimeoutExpired("Singular", 1)):
            with self.assertRaises(subprocess.TimeoutExpired):
                run_singular('print("PASS");', timeout=1)


if __name__ == "__main__":
    unittest.main()
