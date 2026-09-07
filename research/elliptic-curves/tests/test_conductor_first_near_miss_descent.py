from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
PILOT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "conductor_first_family_anchor_pilot_v1.json"
)
S_CLASS_ENVELOPES = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "conductor_first_s_class_envelopes_v1.json"
)
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


def load_fermigier_pari_runner():
    path = CAS / "run_fermigier_rank20_pari_descent.py"
    spec = importlib.util.spec_from_file_location("fermigier_pari_descent", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_conductor_pari_runner():
    path = CAS / "run_conductor_first_pari_diagnostic.py"
    spec = importlib.util.spec_from_file_location("conductor_pari_diagnostic", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_quadratic_specialq_collector():
    path = CAS / "run_fermigier_rank20_fixedfb_quadratic_specialq.py"
    spec = importlib.util.spec_from_file_location("quadratic_specialq_collector", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
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

    def test_magma_dimension_job_avoids_cover_construction(self) -> None:
        generator = load_generator()
        target = self.manifest["targets"][0]
        program = generator.build_selmer_dimension_magma(
            target, manifest_sha256="0" * 64
        )
        self.assertIn("TwoSelmerGroup(E)", program)
        self.assertIn("residual_dim := selmer_dim - 20", program)
        self.assertIn("RESIDUAL_COVERS_REQUIRED", program)
        self.assertNotIn("TwoDescent(", program)
        self.assertNotIn("known :=", program)
        self.assertNotIn('SetClassGroupBounds("GRH")', program.splitlines()[4:])

    def test_magma_three_selmer_job_is_independent_and_unconditional(self) -> None:
        generator = load_generator()
        target = self.manifest["targets"][2]
        program = generator.build_three_selmer_dimension_magma(
            target, manifest_sha256="0" * 64
        )
        self.assertIn("ThreeSelmerGroup(E)", program)
        self.assertIn("rank_upper := selmer_dim - t3dim", program)
        self.assertIn("CLOSED_EXACT_RANK_19_BY_3_SELMER", program)
        self.assertNotIn("ThreeDescent(", program)
        self.assertNotIn('SetClassGroupBounds("GRH")', program.splitlines()[4:])

    def test_magma_generator_accepts_certified_family_pilot_target(self) -> None:
        generator = load_generator()
        _, target = generator.load_target(PILOT, "family2-u481")
        self.assertEqual(target["certified_known_rank"], 14)
        self.assertEqual(target["known_kummer_image_dimension"], 14)
        self.assertEqual(len(target["known_basis"]), 14)
        program = generator.build_selmer_dimension_magma(
            target, manifest_sha256="0" * 64
        )
        self.assertIn("target=family2-u481", program)
        self.assertIn("known_kummer_dim=14", program)

    def test_pari_upper_endpoint_is_not_promoted_without_bnf_certificate(self) -> None:
        runner = load_fermigier_pari_runner()
        lower, classification = runner.classify_bounds(18, 20)
        self.assertEqual(lower, 20)
        self.assertEqual(classification, "P0_grh_conditional_rank20")

    def test_pari_diagnostic_binds_points_to_descent_model(self) -> None:
        runner = load_conductor_pari_runner()
        _, target, model, points = runner.load_exact_input(
            targets.DEFAULT_OUTPUT, "fermigier-u28917-20"
        )
        self.assertEqual(target["certified_known_rank"], 20)
        self.assertEqual(tuple(target["descent_model"][:3]), ("0", "0", "0"))
        self.assertEqual(len(model), 5)
        self.assertEqual(len(points), 20)

    def test_pari_diagnostic_protocol_parser(self) -> None:
        runner = load_conductor_pari_runner()
        parsed = runner.parse_output(
            "\n".join(
                [
                    "CFNMPARI|stage=input|status=complete|pari=[2, 19, 0]",
                    "CFNMPARI|stage=ellrank|status=complete|milliseconds=7|lower=20|upper=21|sha=0|point_count=1",
                    "CFNMPARI|stage=point|index=1|point=[3/2, 5/4]",
                    "CFNMPARI|stage=all|status=complete",
                ]
            )
        )
        self.assertTrue(parsed["completed"])
        self.assertEqual(parsed["rank_upper"], 21)
        self.assertEqual(parsed["points"], [(runner.Q(3, 2), runner.Q(5, 4))])

    def test_s_class_envelope_orders_exact_materialized_factor_bases(self) -> None:
        payload = json.loads(S_CLASS_ENVELOPES.read_text())
        digest_payload = dict(payload)
        declared = digest_payload.pop("result_sha256")
        actual = hashlib.sha256(
            json.dumps(
                digest_payload, sort_keys=True, separators=(",", ":")
            ).encode()
        ).hexdigest()
        self.assertEqual(declared, actual)
        self.assertEqual(
            payload["collector_priority"],
            [
                "icarm-245",
                "family2-u483",
                "fermigier-u28917-20",
                "family3-u660",
            ],
        )
        records = {record["id"]: record for record in payload["targets"]}
        self.assertEqual(
            [
                records[target_id]["bach_plus_selmer_prime_ideal_count"]
                for target_id in payload["collector_priority"]
            ],
            [39904, 40754, 42251, 43512],
        )
        for record in records.values():
            self.assertEqual(
                record["two_division_cubic_coefficients_ascending"][-1], "1"
            )
            self.assertIn("no principal-relation completeness", record["scope"])

    def test_sparse_large_prime_hyperedges_retain_exact_dependency(self) -> None:
        collector = load_quadratic_specialq_collector()
        eliminator = collector.SparseLargePrimeEliminator()
        self.assertEqual(eliminator.add(("a", "b"), 0b001, 0), (None, None))
        self.assertEqual(eliminator.add(("b", "c"), 0b010, 1), (None, None))
        row, provenance = eliminator.add(("a", "c"), 0b100, 2)
        self.assertEqual(row, 0b111)
        self.assertEqual(provenance, {0, 1, 2})
        self.assertEqual(eliminator.edge_count, 3)
        self.assertEqual(len(eliminator.pivots), 2)
        self.assertEqual(eliminator.dependency_count, 1)

        # Repeated vertices cancel in a GF(2) hyperedge.
        row, provenance = eliminator.add(("d", "d"), 0b1000, 3)
        self.assertEqual(row, 0b1000)
        self.assertEqual(provenance, {3})
        self.assertEqual(eliminator.edge_count - len(eliminator.pivots), 2)


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
            "known_subgroup": {
                "status": "complete",
                "certified_rank_lower_bound": rank,
                "kummer_image_dimension": rank,
                "exact_point_membership_checked": True,
            },
            "residual_selmer": {
                "status": "complete",
                "known_kummer_dimension": rank,
                "residual_dimension": 1,
            },
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

    def test_selmer_cannot_precede_full_known_kummer_certificate(self) -> None:
        record = self.certified("missing-subgroup", 19, 100)
        record.pop("known_subgroup")
        with self.assertRaisesRegex(ValueError, "certified known subgroup"):
            validate_candidate(record)

        record = self.certified("thin-subgroup", 19, 100)
        record["known_subgroup"]["kummer_image_dimension"] = 18
        with self.assertRaisesRegex(ValueError, "full-dimensional"):
            validate_candidate(record)

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
