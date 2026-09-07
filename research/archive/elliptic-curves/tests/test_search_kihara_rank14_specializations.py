from __future__ import annotations

from decimal import Decimal
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS_DIRECTORY = ROOT / "elliptic-curves" / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from kihara_discriminant_geometry import (  # noqa: E402
    FACTOR_SIGNATURE,
    derive_discriminant_geometry,
    frontier_value_t,
    hensel_lift_simple_t_root,
    polynomial_derivative_value,
    roots_mod_prime_t,
)
from search_kihara_rank14_specializations import (  # noqa: E402
    TARGET_LOG_CONDUCTOR,
    normalized_specialization,
    select_conductor_candidates,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_kihara_rank14_specializations.json"
)


class KiharaRank14SpecializationSearchTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        cls.geometry = derive_discriminant_geometry()

    def test_degree_398_frontier_is_derived_and_pinned(self) -> None:
        self.assertEqual(self.geometry.discriminant_degree, 1014)
        self.assertEqual(self.geometry.factor_signature, FACTOR_SIGNATURE)
        self.assertEqual(self.geometry.frontier_degree_z, 199)
        digest = hashlib.sha256(
            "\n".join(
                str(value) for value in self.geometry.frontier_coefficients_z
            ).encode()
        ).hexdigest()
        stored = self.data["discriminant_geometry"]
        self.assertEqual(stored["frontier_coefficients_sha256"], digest)
        self.assertEqual(
            stored["factor_signature_degree_exponent"],
            [list(pair) for pair in FACTOR_SIGNATURE],
        )

    def test_simple_roots_hensel_lift_exactly(self) -> None:
        coefficients = self.geometry.frontier_coefficients_z
        roots = roots_mod_prime_t(coefficients, 11)
        simple = tuple(
            root
            for root in roots
            if 2 * root * polynomial_derivative_value(coefficients, root * root) % 11
        )
        self.assertEqual(simple, (3, 8))
        lifted = tuple(
            hensel_lift_simple_t_root(coefficients, root, 11, 2)
            for root in simple
        )
        self.assertEqual(lifted, (36, 85))
        self.assertTrue(all(frontier_value_t(coefficients, root) % 121 == 0 for root in lifted))

    def test_geometry_selection_is_replayed_without_point_data(self) -> None:
        prefilter = self.data["geometry_prefilter"]
        self.assertFalse(prefilter["uses_point_or_rank_data"])
        selected = select_conductor_candidates(
            prefilter["records"], geometry_keep=12, crt_origin_keep=2
        )
        self.assertEqual(list(selected), prefilter["selected_parameters"])
        for record in prefilter["records"]:
            for condition in record["forced_frontier_valuations"]:
                self.assertGreaterEqual(
                    condition["actual_exponent"], condition["requested_exponent"]
                )

    def test_completed_conductors_and_point_stage_status_are_consistent(self) -> None:
        conductor_records = self.data["conductor_stage"]["records"]
        completed = {
            record["parameter_t"]: record["pari"]
            for record in conductor_records
            if record["pari"]["status"] == "completed"
        }
        self.assertEqual(set(completed), {"2", "4", "3/2"})
        self.assertEqual(
            completed["2"]["log_conductor"],
            "213.86643972249849941610804848114713141213263991221224650783370654200357120123670",
        )
        self.assertTrue(
            all(
                record["below_target"]
                == (Decimal(record["log_conductor"]) < TARGET_LOG_CONDUCTOR)
                for record in completed.values()
            )
        )
        self.assertEqual(self.data["conductor_stage"]["subthreshold_count"], 0)
        self.assertEqual(self.data["point_stage_selection"]["parameters"], ["2", "4"])
        for record in self.data["point_stage"]["records"]:
            self.assertEqual(record["initial_unexpected_abscissa_count"], 0)
            self.assertEqual(record["height_stability"]["stable_numerical_rank"], 14)
            self.assertEqual(record["stable_numerical_rank_gain_over_generic_14"], 0)
        self.assertFalse(self.data["target_reached"])

    def test_t2_normalized_model_matches_stored_exactly(self) -> None:
        specialization = normalized_specialization(Q(2))
        record = next(
            row
            for row in self.data["conductor_stage"]["records"]
            if row["parameter_t"] == "2"
        )["pari"]
        self.assertEqual(
            [str(value) for value in specialization.short_coefficients],
            record["normalized_short_coefficients"],
        )


if __name__ == "__main__":
    unittest.main()

