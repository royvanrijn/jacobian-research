#!/usr/bin/env python3
"""Tests for the strictly bounded q12o5867 PARI 2-cover runner."""

from __future__ import annotations

from fractions import Fraction
import importlib.util
from pathlib import Path
import shutil
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "elliptic-curves/scripts/probe_q12o5867_pari_two_cover.py"
SPEC = importlib.util.spec_from_file_location("q12_pari_two_cover", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Q12O5867PariTwoCoverTests(unittest.TestCase):
    def test_protocol_parser_preserves_exact_candidates(self) -> None:
        output = "\n".join(
            (
                "Q12P2|stage=input|status=complete|pari_version=[2, 15, 4]",
                "Q12P2|stage=ell2cover|status=complete|milliseconds=3|cover_count=2",
                "Q12P2|cover_index=1|quartic=x^4 - 8*x + 4",
                "Q12P2|cover_index=1|map_x=4/y^2*x^3 - 4/y^2*x^2 + 4/y^2",
                "Q12P2|cover_index=1|map_y=-1/y^3*x^6 + 8/y^3",
                "Q12P2|stage=cover_search|status=complete|cover_index=1|milliseconds=2|raw_point_count=1",
                "Q12P2|stage=candidate|cover_index=1|cover_x=3/2|cover_y=-7/4|curve_x=5/6|curve_y=11/9",
                "Q12P2|stage=cover|status=complete|cover_index=1",
                "Q12P2|stage=all|status=complete",
            )
        )
        parsed = MODULE.parse_gp_output(output)
        self.assertTrue(parsed["ell2cover_completed"])
        self.assertEqual(parsed["cover_count"], 2)
        self.assertEqual(parsed["completed_cover_indices"], [1])
        self.assertEqual(parsed["covers"][0]["quartic"], "x^4 - 8*x + 4")
        self.assertEqual(parsed["covers"][0]["raw_point_count"], 1)
        self.assertEqual(
            tuple(Fraction(value) for value in parsed["raw_candidate_images"][0]["curve_point"]),
            (Fraction(5, 6), Fraction(11, 9)),
        )
        self.assertTrue(parsed["all_completed"])

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is unavailable")
    def test_pari_cover_map_smoke(self) -> None:
        program = MODULE.gp_program(
            (0, 0, 0, -1, 1), stack_bytes=64_000_000, search_height=100
        )
        completed = subprocess.run(
            ["gp", "-f", "-q", "-s", "64000000"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=20,
        )
        parsed = MODULE.parse_gp_output(completed.stdout)
        self.assertTrue(parsed["ell2cover_completed"])
        self.assertTrue(parsed["all_completed"])
        self.assertGreaterEqual(parsed["cover_count"], 1)
        self.assertGreaterEqual(len(parsed["raw_candidate_images"]), 1)
        for record in parsed["raw_candidate_images"]:
            point = tuple(Fraction(value) for value in record["curve_point"])
            self.assertTrue(MODULE.is_on_weierstrass_curve((0, 0, 0, -1, 1), point))


if __name__ == "__main__":
    unittest.main()
