import unittest
from fractions import Fraction
from unittest.mock import patch

import screen_mestre_two_section_transverse_band as band


class MestreTwoSectionTransverseBandTest(unittest.TestCase):
    def test_all_prime_union_retains_reconstructions_from_later_primes(self):
        calls = []

        def reconstructed_at_prime(moduli, primes, precision):
            calls.append(tuple(primes))
            if tuple(primes) == (7,):
                return ((Fraction(0), Fraction(1)),)
            if tuple(primes) == (11,):
                return ((Fraction(2), Fraction(3)),)
            raise ValueError("bad reduction")

        with patch.object(band, "sections_at_seed", side_effect=reconstructed_at_prime), patch.object(
            band, "DEFAULT_PRIMES", (7, 11, 13)
        ):
            sections = band.sections_across_primes((Fraction(1),) * 5, precision=8)

        self.assertEqual(
            sections,
            ((Fraction(0), Fraction(1)), (Fraction(2), Fraction(3))),
        )
        self.assertEqual(calls, [(7,), (11,), (13,)])


if __name__ == "__main__":
    unittest.main()
