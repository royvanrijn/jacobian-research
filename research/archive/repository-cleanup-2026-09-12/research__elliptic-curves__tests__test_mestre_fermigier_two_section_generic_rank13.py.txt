from __future__ import annotations

import unittest

from verify_mestre_fermigier_two_section_generic_rank13 import replay


class FermigierTwoSectionGenericRank13Test(unittest.TestCase):
    def test_finite_reduction_certificate(self) -> None:
        result = replay()
        self.assertEqual(result["combined_exact_rank_over_F3"], 13)
        self.assertEqual(result["torsion_exclusion"]["finite_group_order"], 28)
        self.assertEqual(len(result["quotient_probes"]), 13)


if __name__ == "__main__":
    unittest.main()
