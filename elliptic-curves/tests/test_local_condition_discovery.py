from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest

CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from fermigier_mestre import FermigierMestreFamily  # noqa: E402
from local_condition_discovery import (  # noqa: E402
    INVARIANT_POLYNOMIALS,
    classify_ball,
    discover_prime,
    polynomial_value,
    select_condition_groups,
    validate_invariant_polynomials,
)
from multiple_root_lifting import RootBall  # noqa: E402


class LocalConditionDiscoveryTests(unittest.TestCase):
    def test_constructed_invariant_polynomials_match_family(self) -> None:
        validate_invariant_polynomials()
        for parameter in (-101, 4, 88):
            invariant_i, invariant_j = FermigierMestreFamily.binary_invariants(
                Fraction(parameter)
            )
            self.assertEqual(
                polynomial_value(INVARIANT_POLYNOMIALS["I"], parameter),
                invariant_i,
            )
            self.assertEqual(
                polynomial_value(INVARIANT_POLYNOMIALS["J"], parameter),
                invariant_j,
            )

    def test_known_clean_and_additive_balls_are_discovered_exactly(self) -> None:
        split = classify_ball(RootBall(11, 1, 5))
        self.assertEqual(split.forced_h_valuation, 4)
        self.assertEqual(split.forced_minimal_discriminant_valuation, 4)
        self.assertEqual(split.reduction, "split multiplicative")

        # The coarse 0 mod 17 ball is correctly marked unresolved because
        # some next digits have extra invariant valuation.  Refinement gives
        # clean additive balls; 17 mod 17^2 is a representative one.
        additive = classify_ball(RootBall(17, 2, 17))
        self.assertEqual(additive.forced_h_valuation, 4)
        self.assertEqual(additive.reduction, "additive")

        seven_child = classify_ball(RootBall(7, 2, 28))
        self.assertEqual(seven_child.presented_model_scaling, 1)
        self.assertEqual(seven_child.forced_minimal_discriminant_valuation, 6)
        self.assertEqual(seven_child.reduction, "split multiplicative")

    def test_prime_13_is_found_without_a_pinned_condition_table(self) -> None:
        record = discover_prime(
            13,
            lift_exponent=4,
            classification_exponent=2,
            max_roots=200_000,
        )
        groups = record["_groups"]
        split_cubic = next(
            group
            for group in groups
            if group.reduction == "split multiplicative"
            and group.forced_minimal_discriminant_valuation == 3
        )
        self.assertEqual(split_cubic.exponent, 1)
        self.assertEqual(split_cubic.residues, (2, 11))
        self.assertEqual(str(split_cubic.reciprocal_density), "13/2")

    def test_default_objective_recovers_four_old_primes_and_new_13(self) -> None:
        groups = []
        for prime in (5, 7, 11, 13, 17, 19, 37, 43, 53):
            groups.extend(
                discover_prime(
                    prime,
                    lift_exponent=4,
                    classification_exponent=2,
                    max_roots=200_000,
                )["_groups"]
            )
        selected = select_condition_groups(
            groups,
            count=5,
            maximum_crt_classes=500,
        )
        self.assertEqual([group.prime for group in selected], [11, 7, 17, 19, 13])
        self.assertEqual(
            [len(group.residues) for group in selected], [3, 6, 2, 2, 2]
        )


if __name__ == "__main__":
    unittest.main()
