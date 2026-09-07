from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elliptic-curves/ecsearch/fermigier_score_sweep.cpp"


class FermigierScoreSweepTests(unittest.TestCase):
    def test_fixed_family_bad_prime_is_skipped(self) -> None:
        compiler = shutil.which("g++")
        if compiler is None:
            self.skipTest("g++ is unavailable")
        with tempfile.TemporaryDirectory() as directory:
            executable = Path(directory) / "sweep"
            subprocess.run(
                [compiler, "-O2", "-fopenmp", "-std=c++20", str(SOURCE), "-o", str(executable)],
                check=True,
                capture_output=True,
                text=True,
            )
            environment = dict(os.environ)
            environment["OMP_NUM_THREADS"] = "1"
            completed = subprocess.run(
                [str(executable), "20", "3", "5"],
                check=True,
                capture_output=True,
                text=True,
                env=environment,
                timeout=30,
            )
        self.assertIn("skipping family-bad prime=5", completed.stderr)
        lines = completed.stdout.splitlines()
        self.assertEqual(len(lines), 6)
        self.assertTrue(lines[0].startswith("# rank\tnumerator\tdenominator"))


if __name__ == "__main__":
    unittest.main()

