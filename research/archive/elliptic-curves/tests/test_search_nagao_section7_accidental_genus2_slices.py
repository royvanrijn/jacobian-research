from __future__ import annotations

import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_nagao_section7_accidental_genus2_slices.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_nagao_section7_accidental_genus2_slices.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "search_nagao_section7_accidental_genus2_slices", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)

from search_nagao_rank21_accidental_slices import rational_square_root  # noqa: E402


class AccidentalGenusTwoSliceTests(unittest.TestCase):
    def test_exact_slice_inventory_and_pinned_prior_union(self) -> None:
        slices, metadata = MODULE.load_genus_two_slices()
        self.assertEqual(len(slices), 240)
        self.assertEqual(metadata["slopes"], list(MODULE.GENUS_TWO_SLOPES))
        self.assertEqual(
            metadata["manifest_sha256"],
            "f79641ff400960ab2e4fd4310a8140619934f3ca795e721d340b67c47803ed7b",
        )
        self.assertEqual(
            {item.slope: sum(other.slope == item.slope for other in slices) for item in slices},
            {slope: 16 for slope in MODULE.GENUS_TWO_SLOPES},
        )
        self.assertTrue(
            all(
                item.normalized.raw_degree == 6
                and item.normalized.normalized_degree == 6
                and item.normalized.genus == 2
                and item.normalized.factor_degrees_and_exponents == ((6, 1),)
                and item.normalized.removed_square_coefficients == (Q(1),)
                for item in slices
            )
        )
        prior, records = MODULE.load_prior_parameters()
        self.assertEqual(len(records), 5)
        self.assertEqual(len(prior), 7_629)
        self.assertEqual(
            MODULE.parameter_stream_sha256(prior),
            "241f0be3acd384b031764f98a9ad2935c4b7c91cc1cf921943eaa911496a8690",
        )

    def test_exact_generic_prior_and_new_point_classification(self) -> None:
        slices, _ = MODULE.load_genus_two_slices()
        item = next(item for item in slices if item.identifier == "a10_sp08")
        prior, _ = MODULE.load_prior_parameters()

        forced_root = rational_square_root(item.normalized.normalized_value(Q(163)))
        self.assertIsNotNone(forced_root)
        forced = MODULE.classify_slice_point(
            item, (Q(163), forced_root), prior, tier="test"
        )
        self.assertTrue(forced["accepted_new_parameter"])
        self.assertEqual(forced["classification"], "new_forced_non_generic_parameter")
        self.assertEqual(forced["generic_quartic_x_labels"], [])
        self.assertEqual(forced["generic_jacobian_sign_pair_labels"], [])

        generic_root = rational_square_root(item.normalized.normalized_value(Q(69)))
        self.assertIsNotNone(generic_root)
        generic = MODULE.classify_slice_point(
            item, (Q(69), generic_root), prior, tier="test"
        )
        self.assertFalse(generic["accepted_new_parameter"])
        self.assertEqual(generic["classification"], "known_generic_section_intersection")
        self.assertIn("visible-00", generic["generic_quartic_x_labels"])

        base_root = rational_square_root(item.normalized.normalized_value(MODULE.T0))
        self.assertIsNotNone(base_root)
        base = MODULE.classify_slice_point(
            item, (MODULE.T0, base_root), prior, tier="test"
        )
        self.assertFalse(base["accepted_new_parameter"])
        self.assertEqual(base["classification"], "prior_parameter_population")

    def test_generated_artifact_replays_the_bounded_frontier(self) -> None:
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            data["reproduction"]["script_sha256"], MODULE.file_sha256(SCRIPT)
        )
        self.assertEqual(data["slice_population"]["slice_count"], 240)
        self.assertEqual(data["search_budget"]["shallow_height"], 5_000)
        self.assertEqual(
            data["search_budget"]["escalated_slice_ids"], ["a10_sp08"]
        )
        self.assertEqual(data["search_runs"]["declared_call_count"], 241)
        self.assertEqual(data["search_runs"]["status_counts"], {"completed": 241})
        candidates = data["new_candidate_population"]
        self.assertEqual(candidates["unique_parameter_count"], 2)
        self.assertEqual(
            candidates["parameter_stream_sha256"],
            "b5982527b7981075effe5d9a312518087e1d3051aa335d522b03ee818181a69f",
        )
        self.assertEqual(
            [record["constructor_parameter_T"] for record in candidates["records_sorted_by_radical_proxy"]],
            ["163", "1049/10"],
        )
        self.assertEqual(data["proxy_filter"]["below_threshold_count"], 2)
        conductors = data["exact_conductors"]
        self.assertEqual(conductors["attempted"], 2)
        self.assertEqual(conductors["completed"], 2)
        self.assertEqual(conductors["sub_182_72_count"], 2)
        self.assertEqual(
            [record["stable_numerical_rank"] for record in data["rank_triage"]["records"]],
            [12, 12],
        )
        self.assertFalse(data["outcome"]["rank21_certified"])
        self.assertFalse(data["outcome"]["breakthrough_curve_found"])


if __name__ == "__main__":
    unittest.main()
