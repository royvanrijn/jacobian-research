from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
STRATA = GENERATED / "elkies-k3-icarm-11952-norm8-low-root-strata-v1.json"
ATLAS = GENERATED / "elkies-k3-icarm-11952-norm8-low-root-atlas-v2.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


class IcarmNormEightLowRootAtlasTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.strata = json.loads(STRATA.read_text())
        cls.atlas = json.loads(ATLAS.read_text())

    def test_complete_exact_stratification(self) -> None:
        self.assertEqual(
            self.strata["status"],
            "PASS_EXACT_COMPLETE_NORM8_LOW_ROOT_STRATIFICATION",
        )
        rows = self.strata["strata"]
        self.assertEqual(
            [row["class_count"] for row in rows],
            [1266, 8410, 20348, 21405, 9861, 2280, 331, 16],
        )
        self.assertEqual(sum(row["class_count"] for row in rows), 63917)
        self.assertEqual([row["root_rank"] for row in rows], list(range(1, 9)))
        self.assertEqual(
            [row["geometric_mw_rank_at_rho_19"] for row in rows],
            list(range(16, 8, -1)),
        )
        self.assertTrue(self.strata["scope"]["all_classes_have_degree_one_sections"])

    def test_curve302_misses_every_complete_stratum(self) -> None:
        self.assertEqual(
            self.atlas["status"],
            "PASS_EXACT_COMPLETE_PRIORITY_ICARM_NORM8_LOW_ROOT_ATLAS",
        )
        row = next(row for row in self.atlas["targets"] if row["curve_id"] == 302)
        self.assertEqual(row["total_class_count"], 63917)
        self.assertEqual(row["modular_excluded_count"], 63917)
        self.assertEqual(row["exact_survivor_count"], 0)
        self.assertTrue(
            all(
                stratum["outcome"] == "MISS_EXACT_COMPLETE_STRATUM"
                and stratum["modular_excluded_count"] == stratum["class_count"]
                for stratum in row["strata"]
            )
        )

    def test_curve398_positive_control_is_true_a1(self) -> None:
        row = next(row for row in self.atlas["targets"] if row["curve_id"] == 398)
        self.assertEqual(row["qq_isomorphic_hit_priority_ranks"], [16875, 63669])
        hit_strata = [
            stratum for stratum in row["strata"]
            if stratum["qq_isomorphic_hit_count"]
        ]
        self.assertEqual(len(hit_strata), 1)
        self.assertEqual(hit_strata[0]["root_lattice"], "A1")
        self.assertEqual(hit_strata[0]["geometric_mw_rank_at_rho_19"], 16)

    def test_hash_chains(self) -> None:
        for document in (self.strata, self.atlas):
            for path, expected in document["inputs"].items():
                self.assertEqual(digest(ROOT / path), expected)


if __name__ == "__main__":
    unittest.main()
