from __future__ import annotations

from decimal import Decimal
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from search_nagao_rank21_mutations import (  # noqa: E402
    CRT_PRIMES,
    DEFAULT_PROFILES,
    DEFAULT_REGIONS,
    LEAD_PARAMETER,
    LEAD_PROFILE,
    REQUESTED_PROFILE,
    TARGET_LOG_CONDUCTOR,
    conductor_radical_proxy,
    best_good_prime_residues,
    enumerate_mutations,
    enumerate_rank_friendly_lattice_mutations,
    nearest_congruent_numerator,
    profile_contains,
    rational_residue,
)


class NagaoRank21MutationTests(unittest.TestCase):
    def test_lead_profile_correction_and_two_crt_classes_are_exact(self) -> None:
        self.assertEqual(
            tuple(rational_residue(LEAD_PARAMETER, prime) for prime in CRT_PRIMES),
            (1, 4, 9, 3),
        )
        self.assertEqual(LEAD_PROFILE.crt_residue, 1901)
        self.assertEqual(REQUESTED_PROFILE.crt_residue, 361)
        self.assertTrue(profile_contains(LEAD_PROFILE, LEAD_PARAMETER))
        self.assertFalse(profile_contains(REQUESTED_PROFILE, LEAD_PARAMETER))

    def test_nearest_congruent_numerator_is_locally_minimal(self) -> None:
        for profile in DEFAULT_PROFILES:
            for denominator in (1, 32, 97, 199, 503):
                if any(denominator % prime == 0 for prime in CRT_PRIMES):
                    continue
                base, center = nearest_congruent_numerator(
                    LEAD_PARAMETER, denominator, profile.crt_residue
                )
                numerator = base + center * 5005
                distance = abs(Q(numerator, denominator) - LEAD_PARAMETER)
                self.assertLessEqual(
                    distance,
                    abs(Q(numerator - 5005, denominator) - LEAD_PARAMETER),
                )
                self.assertLessEqual(
                    distance,
                    abs(Q(numerator + 5005, denominator) - LEAD_PARAMETER),
                )

    def test_lead_radical_proxy_replays_its_known_log_conductor(self) -> None:
        proxy = conductor_radical_proxy(LEAD_PARAMETER)
        self.assertLess(
            abs(
                Decimal(str(proxy["log_radical_upper_proxy"]))
                - Decimal("138.23366243820961694")
            ),
            Decimal("1e-12"),
        )
        valuations = {
            prime: exponent for prime, exponent in proxy["small_prime_valuations"]
        }
        self.assertEqual(
            {prime: valuations[prime] for prime in (5, 7, 11, 13)},
            {5: 7, 7: 6, 11: 5, 13: 4},
        )

    def test_default_population_is_new_exact_and_pinned(self) -> None:
        candidates = enumerate_mutations()
        self.assertEqual(len(candidates), 874)
        self.assertEqual(len({candidate.parameter for candidate in candidates}), 874)
        self.assertTrue(
            all(
                Decimal(str(candidate.radical_proxy["log_radical_upper_proxy"]))
                < TARGET_LOG_CONDUCTOR
                for candidate in candidates
            )
        )
        for candidate in candidates:
            self.assertTrue(profile_contains(candidate.profile, candidate.parameter))
            if candidate.region_label == DEFAULT_REGIONS[0].label:
                self.assertLessEqual(candidate.parameter.denominator, 200)
                self.assertGreater(abs(candidate.parameter.numerator), 10_000)
            elif candidate.region_label == DEFAULT_REGIONS[1].label:
                self.assertGreaterEqual(candidate.parameter.denominator, 201)
                self.assertLessEqual(candidate.parameter.denominator, 1_000)
            else:
                self.fail("unexpected mutation region")

    def test_rank_friendly_lookup_and_lattice_population_are_exact(self) -> None:
        self.assertEqual(
            [record["residue"] for record in best_good_prime_residues(31)],
            [7, 14, 17, 24],
        )
        self.assertTrue(
            all(record["ellap"] == -11 for record in best_good_prime_residues(31))
        )
        self.assertEqual(
            [record["residue"] for record in best_good_prime_residues(43)],
            [3, 5, 6, 17, 26, 37, 38, 40],
        )
        candidates = enumerate_rank_friendly_lattice_mutations()
        self.assertEqual(len(candidates), 78)
        self.assertEqual(len(candidates), len({candidate.parameter for candidate in candidates}))
        for candidate in candidates:
            self.assertTrue(profile_contains(LEAD_PROFILE, candidate.parameter))
            self.assertEqual(len(candidate.rank_friendly_residues), 2)
            for prime, residue in candidate.rank_friendly_residues:
                self.assertEqual(rational_residue(candidate.parameter, prime), residue)


if __name__ == "__main__":
    unittest.main()
