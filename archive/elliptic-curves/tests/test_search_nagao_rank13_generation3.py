from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
sys.path.insert(0, str(CAS))

from nagao_1994 import (  # noqa: E402
    RANK13_CONSTRUCTION,
    primitive_quartic_coefficients,
    quartic_value,
)
from search_nagao_rank13_generation3 import (  # noqa: E402
    EXACT_RANK17_PARAMETERS,
    generate_generation3_population,
    label_crt_classes,
    load_exact_rank17_labels,
    load_prior_population,
    local_trace_tables,
    rank17_label_residues,
    shifted_polynomial_coefficients,
)


Q = Fraction


class Generation3UnitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.prior_paths = (
            GENERATED / "elliptic_nagao_rank13_rank_gain_search.json",
            GENERATED / "elliptic_nagao_rank13_rank_gain_mutations.json",
        )
        cls.trace_tables = local_trace_tables()
        cls.residue_groups = rank17_label_residues(cls.trace_tables)

    def test_exact_rank17_labels_and_local_groups_are_pinned(self) -> None:
        labels = load_exact_rank17_labels(
            GENERATED / "elliptic_nagao_rank17_frontier_certificate.json"
        )
        self.assertEqual(labels, EXACT_RANK17_PARAMETERS)
        self.assertEqual(
            self.residue_groups,
            {29: (13, 16, 24), 41: (1, 6), 43: (3, 42)},
        )
        classes = label_crt_classes(self.residue_groups)
        self.assertEqual(len(classes), 12)
        self.assertTrue(all(modulus == 29 * 41 * 43 for _, modulus, _ in classes))

    def test_population_is_new_and_high_denominator(self) -> None:
        prior = load_prior_population(self.prior_paths)
        population = generate_generation3_population(self.residue_groups, prior)
        self.assertEqual(len(population), 25_782)
        self.assertTrue(
            all(candidate.parameter_u.denominator > 64 for candidate in population)
        )
        self.assertFalse({candidate.parameter_u for candidate in population} & prior)
        origin_counts = {
            origin: sum(origin in candidate.origins for candidate in population)
            for origin in {
                origin for candidate in population for origin in candidate.origins
            }
        }
        self.assertEqual(origin_counts["pell-small-T-shell"], 9_968)
        self.assertEqual(origin_counts["rank17-local-crt"], 2_408)

    def test_shifted_polynomial_is_exact(self) -> None:
        parameter_t = Q(49, 1228)
        center = Q(703, 15) + parameter_t / 15
        coefficients = primitive_quartic_coefficients(
            RANK13_CONSTRUCTION, parameter_t
        )
        shifted = shifted_polynomial_coefficients(coefficients, center)
        for value in (Q(-7, 11), Q(0), Q(13, 17), Q(5)):
            self.assertEqual(
                quartic_value(shifted, value),
                quartic_value(coefficients, value + center),
            )


class Generation3ArtifactTests(unittest.TestCase):
    def test_artifact_if_present(self) -> None:
        path = GENERATED / "elliptic_nagao_rank13_generation3.json"
        if not path.exists():
            self.skipTest("generation-3 search artifact has not been generated")
        data = json.loads(path.read_text())
        self.assertEqual(data["population"]["count"], 25_782)
        self.assertTrue(data["population"]["disjoint_from_both_prior_populations"])
        self.assertGreaterEqual(data["population"]["minimum_reduced_denominator"], 65)
        self.assertEqual(
            data["script_sha256"],
            hashlib.sha256(
                (CAS / "search_nagao_rank13_generation3.py").read_bytes()
            ).hexdigest(),
        )
        summary = data["summary"]
        self.assertIn("maximum_stable_numerical_rank", summary)
        self.assertIn("maximum_exact_rank_lower_bound_newly_certified", summary)
        for record in summary["frontier"]:
            self.assertTrue(record["height_rank_stable_across_precisions"])
            self.assertTrue(record["all_returned_points_checked_exactly"])


if __name__ == "__main__":
    unittest.main()
