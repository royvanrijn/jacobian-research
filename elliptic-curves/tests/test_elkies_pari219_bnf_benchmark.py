from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import run_elkies_2026_pari219_bnf_benchmark as benchmark  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_relative_2selmer_open_bottleneck_benchmarks_v2.json"
)
PROGRAM = CAS / "run_elkies_2026_pari219_bnf_benchmark.py"


class Pari219BnfBenchmarkTests(unittest.TestCase):
    def test_progress_parser_retains_factorbase_and_latest_deficit(self) -> None:
        log = """\
LIMC = 721, LIMC2 = 96259
KCZ = 89, KC = 163, n = 170
#### Look for 116 relations in 111 ideals (rnd_rel)
LIMC = 1442, LIMC2 = 96259
KCZ = 151, KC = 261, n = 268
#### Look for 191 relations in 191 ideals (rnd_rel)
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "run.log"
            path.write_text(log)
            parsed = benchmark.parse_progress(path)
        self.assertEqual(
            parsed["factorbase_stages"],
            [
                {
                    "limc": 721,
                    "limc2": 96259,
                    "kcz": 89,
                    "factorbase_ideals": 163,
                    "relation_target_with_units": 170,
                },
                {
                    "limc": 1442,
                    "limc2": 96259,
                    "kcz": 151,
                    "factorbase_ideals": 261,
                    "relation_target_with_units": 268,
                },
            ],
        )
        self.assertEqual(parsed["random_relation_round_count"], 2)
        self.assertEqual(
            parsed["latest_random_relation_deficit"],
            {"relations_requested": 191, "ideals_requested": 191},
        )

    def test_default_is_the_pinned_rank21_reduced_cubic(self) -> None:
        self.assertEqual(
            benchmark.RANK21_REDUCED_CUBIC,
            "x^3-x^2-774250153578278482962797863407542*x+"
            "4105678984643853583390832544029019669185034999158",
        )

    def test_pinned_long_run_remains_fail_closed(self) -> None:
        data = json.loads(ARTIFACT.read_text())["pari_219_threaded_relation_experiments"]
        self.assertEqual(
            data["program_sha256"], sha256(PROGRAM.read_bytes()).hexdigest()
        )
        self.assertEqual(
            data["build"]["source_commit"],
            "6af5b91cfaeb6939331945f301e65bd775f6cdef",
        )
        long_run = next(run for run in data["runs"] if run["timeout_seconds"] == 300)
        self.assertEqual(long_run["outcome"], "strict_wall_timeout_with_large_remaining_relation_deficit")
        self.assertEqual(long_run["last_factorbase_ideals"], 1996)
        self.assertEqual(long_run["latest_requested_relations"], 1635)
        self.assertIn("not class-group", data["warning"])


if __name__ == "__main__":
    unittest.main()
