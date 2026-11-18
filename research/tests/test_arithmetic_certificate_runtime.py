"""Certificate verifiers must reject assertion-disabled Python before output."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ARITHMETIC = Path(__file__).resolve().parents[1] / "arithmetic"
CHECKERS = (
    "verify_primitive_norm_class3_profiles.py",
    "verify_primitive_norm_class4_profiles.py",
    "verify_primitive_norm_coordinate_membership.py",
    "verify_primitive_norm_exact_sumsets.py",
    "verify_primitive_norm_genus_filters.py",
    "verify_primitive_norm_local_structure.py",
    "verify_primitive_norm_maximal_orders.py",
    "verify_primitive_norm_minus11.py",
    "verify_primitive_norm_minus20.py",
    "verify_primitive_norm_minus23.py",
    "verify_primitive_norm_minus23_sumset.py",
    "verify_primitive_norm_minus56.py",
    "verify_primitive_norm_minus56_sumset.py",
    "verify_primitive_norm_minus7.py",
    "verify_primitive_norm_minus8.py",
    "verify_primitive_norm_odd_class_absorbers.py",
    "verify_primitive_norm_quaternary.py",
)


class CertificateRuntimeTests(unittest.TestCase):
    def test_optimized_execution_cannot_generate_a_certificate(self):
        for name in CHECKERS:
            for mode in ("-O", "-OO", "environment"):
                with self.subTest(checker=name, mode=mode):
                    env = os.environ.copy()
                    env.pop("PYTHONOPTIMIZE", None)
                    command = [sys.executable]
                    if mode == "environment":
                        env["PYTHONOPTIMIZE"] = "1"
                    else:
                        command.append(mode)
                    with tempfile.TemporaryDirectory() as directory:
                        output = Path(directory) / "unexpected-certificate.json"
                        command.extend((str(ARITHMETIC / name), "--output", str(output)))
                        result = subprocess.run(
                            command, env=env, capture_output=True, text=True, timeout=5)
                        self.assertNotEqual(result.returncode, 0)
                        self.assertIn("requires active Python assertions", result.stderr)
                        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()

