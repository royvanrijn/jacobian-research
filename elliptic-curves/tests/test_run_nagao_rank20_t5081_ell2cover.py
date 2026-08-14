from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/tools/run_nagao_rank20_t5081_ell2cover.py"
SPEC = importlib.util.spec_from_file_location(
    "run_nagao_rank20_t5081_ell2cover", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class T5081Ell2CoverRunnerTests(unittest.TestCase):
    def test_program_has_exact_curve_and_cover_markers(self) -> None:
        program = MODULE.gp_program(1_000_000_000)
        self.assertIn("ell2cover(E)", program)
        for coefficient in MODULE.MINIMAL_MODEL:
            self.assertIn(str(coefficient), program)
        self.assertIn("CONDUCTOR=", program)
        self.assertIn("COVER_COUNT=", program)
        self.assertIn("COVER_", program)

    def test_pinned_certificate_and_candidate(self) -> None:
        data = MODULE.validate_certificate()
        self.assertEqual(
            MODULE.sha256_file(MODULE.CERTIFICATE), MODULE.CERTIFICATE_SHA256
        )
        self.assertEqual(data["candidate"]["constructor_parameter_T"], "5081/47")
        self.assertEqual(
            data["exact_rank_certificate"]["certified_algebraic_rank_lower_bound"],
            20,
        )

    def test_marker_parser(self) -> None:
        self.assertEqual(MODULE.parse_marker_integer("A=1\nCOVER_COUNT=4\n", "COVER_COUNT"), 4)
        self.assertIsNone(MODULE.parse_marker_integer("A=1\n", "COVER_COUNT"))


if __name__ == "__main__":
    unittest.main()
