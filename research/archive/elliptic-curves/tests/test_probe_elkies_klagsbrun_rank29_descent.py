#!/usr/bin/env python3

from pathlib import Path
import sys
import unittest


TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))

from probe_elkies_klagsbrun_rank29_descent import gp_program  # noqa: E402


class ElkiesKlagsbrunDescentProbeTests(unittest.TestCase):
    def test_gp_program_pins_curve_points_and_markers(self) -> None:
        program = gp_program(200_000_000)
        self.assertIn("default(parisizemax,200000000)", program)
        self.assertIn("ellrank(E,0,P)", program)
        self.assertIn("EXACT_POINTS_ON_CURVE=", program)
        self.assertIn("RANK_LOWER=", program)
        self.assertEqual(program.count("],["), 28)


if __name__ == "__main__":
    unittest.main()
