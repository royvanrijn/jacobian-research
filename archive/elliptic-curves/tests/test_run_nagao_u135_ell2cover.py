from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/tools/run_nagao_u135_ell2cover.py"
SPEC = importlib.util.spec_from_file_location("run_nagao_u135_ell2cover", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class U135Ell2CoverRunnerTests(unittest.TestCase):
    def test_program_has_exact_curve_and_cover_markers(self) -> None:
        program = MODULE.gp_program(1_000_000_000)
        self.assertIn("ell2cover(E)", program)
        self.assertIn(str(MODULE.SHORT_A.numerator), program)
        self.assertIn(str(MODULE.SHORT_B.numerator), program)
        self.assertIn("COVER_COUNT=", program)
        self.assertIn("COVER_", program)

    def test_pinned_certificate_hash(self) -> None:
        self.assertEqual(
            MODULE.sha256_file(MODULE.CERTIFICATE), MODULE.CERTIFICATE_SHA256
        )


if __name__ == "__main__":
    unittest.main()
