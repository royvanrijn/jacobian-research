from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import shutil
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from nagao_1994 import (  # noqa: E402
    RANK13_CONSTRUCTION,
    rank13_base_changed_short_jacobian_coefficients,
    rank13_base_parameter,
)
from nagao_rank13_local import (  # noqa: E402
    BASE_CHANGED_DISCRIMINANT,
    BASE_CHANGED_INVARIANTS,
    DEFAULT_CRT_PRIMES,
    default_crt_balls,
    discover_local_conditions,
    polynomial_value_mod,
    rational_discriminant_valuation,
)
from pari_bridge import minimal_curve_data  # noqa: E402
from search_nagao_rank13_local_crt import (  # noqa: E402
    enumerate_crt_candidates,
)


def evaluate(polynomial: tuple[int, ...], value: Q) -> Q:
    answer = Q(0)
    for coefficient in reversed(polynomial):
        answer = answer * value + coefficient
    return answer


class NagaoRank13LocalTests(unittest.TestCase):
    def test_cleared_invariants_match_exact_family_coefficients(self) -> None:
        for parameter_u in (Q(1), Q(7, 3), Q(118), Q(-5, 2)):
            coefficients = rank13_base_changed_short_jacobian_coefficients(parameter_u)
            self.assertEqual(
                evaluate(BASE_CHANGED_INVARIANTS["a4"], parameter_u)
                / (2 * parameter_u) ** 8,
                coefficients[3],
            )
            self.assertEqual(
                evaluate(BASE_CHANGED_INVARIANTS["a6"], parameter_u)
                / (2 * parameter_u) ** 12,
                coefficients[4],
            )

    def test_cleared_discriminant_matches_quartic_discriminant(self) -> None:
        for parameter_u in (Q(1), Q(7, 3), Q(118), Q(-5, 2)):
            parameter_t = rank13_base_parameter(parameter_u)
            self.assertEqual(
                evaluate(BASE_CHANGED_DISCRIMINANT, parameter_u)
                / (2 * parameter_u) ** 20,
                RANK13_CONSTRUCTION.primitive_quartic_discriminant(parameter_t),
            )

    def test_declared_root_lift_profiles_are_exact(self) -> None:
        expected = {
            7: ((7, 49, 308, 2072, 12362), 2),
            11: ((5, 55, 506), 0),
            13: ((12, 156, 2028), 0),
            17: ((12, 204, 1156), 0),
            19: ((14, 266, 1444), 0),
            23: ((10, 230, 1058), 0),
            31: ((12, 372, 3844), 0),
        }
        targets = {7: 5, 11: 3, 13: 3, 17: 3, 19: 3, 23: 3, 31: 3}
        for prime, (counts, fixed) in expected.items():
            discovery = discover_local_conditions(prime, targets[prime])
            self.assertEqual(discovery.level_counts, counts)
            self.assertEqual(discovery.fixed_divisor_valuation, fixed)

    def test_split_and_additive_unit_symbols_are_distinguished(self) -> None:
        groups = default_crt_balls()
        self.assertEqual(tuple(groups), DEFAULT_CRT_PRIMES)
        self.assertEqual(
            {prime: {ball.residue for ball in balls} for prime, balls in groups.items()},
            {
                7: {1, 2, 5, 6},
                11: {3, 4, 7, 8},
                13: {1, 4, 5, 6, 7, 8, 9, 12},
                19: {4, 7, 12, 15},
                31: {6, 12, 19, 25},
            },
        )
        for balls in groups.values():
            self.assertTrue(all(ball.reduction == "split multiplicative" for ball in balls))
            self.assertTrue(all(ball.conductor_exponent == 1 for ball in balls))

        at_13 = discover_local_conditions(13, 3)
        additive_13 = {
            ball.residue
            for ball in at_13.balls
            if ball.reduction.startswith("additive")
        }
        self.assertEqual(additive_13, {2, 3, 10, 11})
        self.assertTrue(
            all(
                ball.reduction.startswith("additive")
                for ball in discover_local_conditions(17, 3).balls
            )
        )
        self.assertTrue(
            all(
                ball.reduction.startswith("additive")
                for ball in discover_local_conditions(23, 3).balls
            )
        )

    def test_u118_replays_all_forced_valuations(self) -> None:
        self.assertEqual(
            {
                prime: rational_discriminant_valuation(118, 1, prime)
                for prime in DEFAULT_CRT_PRIMES
            },
            {7: 6, 11: 4, 13: 4, 19: 3, 31: 3},
        )
        # The mod-p c4 numerators are units at all five engineered primes.
        self.assertTrue(
            all(
                polynomial_value_mod(BASE_CHANGED_INVARIANTS["c4"], 118, prime)
                for prime in DEFAULT_CRT_PRIMES
            )
        )

    def test_crt_lattice_search_identifies_sign_duplicate_curves(self) -> None:
        candidates = enumerate_crt_candidates(
            coefficient_radius=6, representatives_per_class=1
        )
        self.assertEqual(len(candidates), 1024)
        self.assertEqual(
            min(
                candidates,
                key=lambda candidate: max(
                    abs(candidate.parameter_t.numerator),
                    candidate.parameter_t.denominator,
                ),
            ).parameter_u,
            118,
        )
        self.assertEqual(len({abs(candidate.parameter_t) for candidate in candidates}), 1024)

    @unittest.skipUnless(shutil.which("gp"), "PARI/GP is optional")
    def test_u118_conductor_frontier(self) -> None:
        curve = minimal_curve_data(
            rank13_base_changed_short_jacobian_coefficients(Q(118)), timeout=10
        )
        self.assertEqual(
            curve["conductor"],
            39951290420847784070371655775988962021455042106946985910,
        )
        self.assertTrue(curve["log_conductor"].startswith("128.027255994266404107"))
        self.assertEqual(curve["root_number"], -1)


if __name__ == "__main__":
    unittest.main()

