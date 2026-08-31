#!/usr/bin/env python3
"""End-to-end fail-closed tests for every active q12 point-search command."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
CAS = ELLIPTIC_ROOT / "cas"
sys.path[:0] = [str(ELLIPTIC_ROOT), str(CAS)]

from ecsearch.rank_certification import add_rational_points  # noqa: E402
from elkies_residual_selmer_gate import INCOMPLETE_STATUS, SCHEMA  # noqa: E402


class Q12O5867DescentFirstEntrypointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory()
        directory = Path(cls.temporary.name)
        model = tuple(map(Fraction, (0, 0, 0, -1, 1)))
        generator = (Fraction(0), Fraction(1))
        points = []
        point = None
        for _ in range(17):
            point = add_rational_points(model, point, generator)
            assert point is not None
            points.append(point)
        assert len(set(points)) == 17
        cls.specialization = directory / "specialization.json"
        cls.specialization.write_text(
            json.dumps(
                {
                    "status": "PASS_EXACT_Q12O5867_SPECIALIZED_GENERIC_RANK17_LOWER_BOUND",
                    "parameter": {"affine_value": "1/2"},
                    "global_minimal_specialization": {
                        "model": [str(value) for value in model],
                        "points": [[str(x), str(y)] for x, y in points],
                    },
                }
            )
        )
        cls.incomplete_gate = directory / "incomplete-gate.json"
        cls.incomplete_gate.write_text(
            json.dumps(
                {
                    "schema": SCHEMA,
                    "status": INCOMPLETE_STATUS,
                    "parameter": "1/2",
                    "global_minimal_model": [str(value) for value in model],
                    "descent_backend": {
                        "unconditional": False,
                        "class_group_completeness_completed": False,
                        "all_local_solubility_conditions_completed": False,
                    },
                    "gate": {
                        "residual_two_selmer_quotient_dimension": None,
                        "required_residual_dimension": 15,
                        "expensive_search_authorized": False,
                    },
                }
            )
        )
        cls.wrong_parameter_gate = directory / "wrong-parameter-gate.json"
        cls.wrong_parameter_gate.write_text(
            json.dumps(
                {
                    "schema": SCHEMA,
                    "status": "PASS_RANK32_RESIDUAL_2_SELMER_GATE",
                    "parameter": "3/4",
                    "global_minimal_model": [str(value) for value in model],
                    "descent_backend": {
                        "unconditional": True,
                        "class_group_completeness_completed": True,
                        "all_local_solubility_conditions_completed": True,
                    },
                    "gate": {
                        "residual_two_selmer_quotient_dimension": 15,
                        "required_residual_dimension": 15,
                        "expensive_search_authorized": True,
                    },
                }
            )
        )
        cls.directory = directory

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temporary.cleanup()

    def assert_gate_rejects(self, script: str, arguments: list[str], output: Path) -> None:
        completed = subprocess.run(
            [sys.executable, str(ELLIPTIC_ROOT / "scripts" / script), *arguments],
            text=True,
            capture_output=True,
            check=False,
            timeout=10,
        )
        self.assertNotEqual(completed.returncode, 0, completed.stdout)
        self.assertIn("expensive search is forbidden", completed.stderr)
        self.assertFalse(output.exists(), f"{script} created output before its gate")

    def test_direct_ratpoints_rejects_before_search(self) -> None:
        output = self.directory / "ratpoints.json"
        self.assert_gate_rejects(
            "probe_q12o5867_ratpoints.py",
            [
                "--input", str(self.specialization),
                "--residual-selmer-gate", str(self.incomplete_gate),
                "--height", "1",
                "--denominator-bound", "1",
                "--output", str(output),
            ],
            output,
        )

    def test_affine_charts_reject_before_directory_creation(self) -> None:
        output = self.directory / "charts.json"
        self.assert_gate_rejects(
            "probe_q12o5867_section_charts.py",
            [
                "--input", str(self.specialization),
                "--residual-selmer-gate", str(self.incomplete_gate),
                "--chart-limit", "1",
                "--height", "1",
                "--denominator-bound", "1",
                "--output", str(output),
            ],
            output,
        )
        self.assertFalse(output.with_suffix("").exists())

    def test_mwrank_parent_rejects_before_sage_start(self) -> None:
        output = self.directory / "mwrank.json"
        self.assert_gate_rejects(
            "probe_q12o5867_mwrank.py",
            [
                "--input", str(self.specialization),
                "--residual-selmer-gate", str(self.incomplete_gate),
                "--sage", str(self.directory / "must-not-run-sage"),
                "--output", str(output),
            ],
            output,
        )

    def test_normalized_slope_search_rejects_before_search(self) -> None:
        output = self.directory / "slopes.json"
        self.assert_gate_rejects(
            "search_q12o5867_section_slope_slices.py",
            [
                str(self.specialization),
                "--residual-selmer-gate", str(self.incomplete_gate),
                "--output", str(output),
            ],
            output,
        )

    def test_pari_two_cover_rejects_before_gp_start(self) -> None:
        output = self.directory / "pari-cover.json"
        self.assert_gate_rejects(
            "probe_q12o5867_pari_two_cover.py",
            [
                "--input", str(self.specialization),
                "--residual-selmer-gate", str(self.incomplete_gate),
                "--output", str(output),
                "--search-height", "1",
                "--timeout", "1",
                "--stack-bytes", "64000000",
                "--rss-limit-bytes", "64000000",
            ],
            output,
        )

    def test_pari_two_cover_also_binds_the_parameter(self) -> None:
        output = self.directory / "pari-cover-wrong-parameter.json"
        completed = subprocess.run(
            [
                sys.executable,
                str(ELLIPTIC_ROOT / "scripts/probe_q12o5867_pari_two_cover.py"),
                "--input", str(self.specialization),
                "--residual-selmer-gate", str(self.wrong_parameter_gate),
                "--output", str(output),
                "--search-height", "1",
                "--timeout", "1",
                "--stack-bytes", "64000000",
                "--rss-limit-bytes", "64000000",
            ],
            text=True,
            capture_output=True,
            check=False,
            timeout=10,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("different fibre", completed.stderr)
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
