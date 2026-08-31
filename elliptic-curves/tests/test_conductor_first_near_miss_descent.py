from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

import build_conductor_first_near_miss_targets as targets  # noqa: E402
from conductor_first_pipeline import (  # noqa: E402
    discriminant_sieve_record,
    pareto_frontier,
    rank_first_order,
    validate_candidate,
)


def load_generator():
    path = CAS / "build_conductor_first_near_miss_magma.py"
    spec = importlib.util.spec_from_file_location("near_miss_magma", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ConductorFirstNearMissTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest = targets.build_manifest()

    def test_four_targets_have_full_known_mod2_images(self) -> None:
        records = self.manifest["targets"]
        self.assertEqual(
            [record["id"] for record in records],
            [
                "icarm-245",
                "fermigier-u28917-20",
                "family2-u483",
                "family3-u660",
            ],
        )
        self.assertEqual(
            [record["certified_known_rank"] for record in records],
            [20, 20, 19, 19],
        )
        for record in records:
            self.assertEqual(
                record["known_kummer_image_dimension"],
                record["certified_known_rank"],
            )
            self.assertEqual(record["rational_two_torsion_dimension"], 0)
            self.assertEqual(
                len(record["known_basis"]), record["certified_known_rank"]
            )

    def test_mestre_replacement_bases_are_not_promoted_to_saturation_theorems(self) -> None:
        for record in self.manifest["targets"][2:]:
            provenance = record["basis_provenance"]
            self.assertEqual(provenance["saturation_prime_bound"], 3)
            self.assertIn("no global saturation claim", provenance["claim"])
            self.assertEqual(
                record["known_basis_mod2_certificate"]["combined_binary_rank"], 19
            )

    def test_magma_job_is_relative_and_unconditional(self) -> None:
        generator = load_generator()
        target = self.manifest["targets"][0]
        program = generator.build_magma(
            target, manifest_sha256="0" * 64, pairing_cover_cap=63
        )
        self.assertIn("TwoDescent(", program)
        self.assertIn("RemoveGens := SequenceToSet(known)", program)
        self.assertIn("CLOSED_EXACT_RANK_20", program)
        self.assertIn("HIGHER_DESCENT_REQUIRED", program)
        self.assertNotIn('SetClassGroupBounds("GRH")', program.splitlines()[3:])


class ConductorFirstPipelineTest(unittest.TestCase):
    @staticmethod
    def certified(candidate_id: str, rank: int, conductor: int) -> dict:
        return {
            "id": candidate_id,
            "cheap_sieve": discriminant_sieve_record(((2, 4), (5, 1)), complete=True),
            "tate": {
                "status": "complete",
                "global_minimal": True,
                "conductor": str(conductor),
                "local_reductions": [{"prime": 2, "conductor_exponent": 1}],
            },
            "residual_selmer": {"status": "complete", "residual_dimension": 1},
            "residual_covers": {
                "status": "complete",
                "locally_surviving_cover_ids": ["cover-1"],
            },
            "point_recovery": {
                "status": "certified",
                "certified_rank_lower_bound": rank,
            },
        }

    def test_point_recovery_cannot_precede_residual_selmer(self) -> None:
        record = {
            "id": "bad-order",
            "cheap_sieve": discriminant_sieve_record(((2, 2),), complete=True),
            "point_recovery": {"status": "bounded"},
        }
        with self.assertRaisesRegex(ValueError, "require exact Tate data"):
            validate_candidate(record)

    def test_zero_residual_dimension_closes_before_points(self) -> None:
        record = self.certified("closed", 21, 100)
        record["residual_selmer"]["residual_dimension"] = 0
        record.pop("residual_covers")
        record.pop("point_recovery")
        self.assertEqual(validate_candidate(record), "closed")

    def test_rank_first_pareto_is_not_a_blended_score(self) -> None:
        records = [
            self.certified("r21-small", 21, 100),
            self.certified("r21-large", 21, 200),
            self.certified("r22-large", 22, 1_000),
            self.certified("r20-tiny", 20, 10),
        ]
        self.assertEqual(
            [record["id"] for record in rank_first_order(records)],
            ["r22-large", "r21-small", "r21-large", "r20-tiny"],
        )
        self.assertEqual(
            [record["id"] for record in pareto_frontier(records)],
            ["r22-large", "r21-small", "r20-tiny"],
        )


if __name__ == "__main__":
    unittest.main()
