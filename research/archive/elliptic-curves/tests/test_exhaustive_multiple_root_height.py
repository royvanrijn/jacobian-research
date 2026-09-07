from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from exhaustive_multiple_root_height import (  # noqa: E402
    brute_force_projective_height,
    enumerate_projective_height,
    homogeneous_discriminant_factor,
    negation_orbits,
    parse_keep_counts,
    verify_local_constraints,
)
from fermigier_mestre import FermigierMestreFamily  # noqa: E402
from search_multiple_root_crt import (  # noqa: E402
    DEFAULT_GROUP_ORDER,
    choose_groups,
    crt_classes,
)


class ExhaustiveMultipleRootHeightTests(unittest.TestCase):
    def test_general_enumerator_matches_tiny_brute_force(self) -> None:
        modulus = 15
        residues = (1, 4, 11, 14)
        classes = tuple(
            {
                "class_index": index,
                "crt_residue": residue,
                "crt_modulus": modulus,
            }
            for index, residue in enumerate(residues, 1)
        )
        height = 37  # Deliberately exceeds twice the toy modulus.
        optimized = enumerate_projective_height(classes, height)
        brute_force = brute_force_projective_height(classes, height)
        self.assertEqual(
            tuple((item.numerator, item.denominator) for item in optimized),
            brute_force,
        )
        self.assertGreater(len(optimized), 0)
        self.assertTrue(all(item.numerator >= 0 for item in optimized))

    def test_negation_stability_is_required(self) -> None:
        self.assertEqual(negation_orbits((1, 4, 11, 14), 15), ((1, 14), (4, 11)))
        with self.assertRaises(ValueError):
            negation_orbits((1, 4, 14), 15)

    def test_homogenization_agrees_with_fraction_evaluation(self) -> None:
        numerator, denominator = 70, 223
        integer = homogeneous_discriminant_factor(numerator, denominator)
        rational = FermigierMestreFamily.discriminant_factor(
            Fraction(numerator, denominator)
        )
        self.assertEqual(integer, rational * denominator**20)

    def test_known_shortest_candidate_passes_all_five_exact_checks(self) -> None:
        groups = choose_groups(DEFAULT_GROUP_ORDER, ())
        classes = crt_classes(groups)
        candidates = enumerate_projective_height(classes, 223)
        candidate = next(item for item in candidates if item.identifier == "70/223")
        local = verify_local_constraints(candidate, groups)
        self.assertFalse(local["singular"])
        self.assertEqual(
            local["h_valuations"],
            {"7": 18, "11": 5, "17": 5, "19": 4, "37": 3},
        )

    def test_keep_schedule_validation(self) -> None:
        self.assertEqual(parse_keep_counts("256,32,12"), (256, 32, 12))
        with self.assertRaises(Exception):
            parse_keep_counts("12,32")


if __name__ == "__main__":
    unittest.main()
