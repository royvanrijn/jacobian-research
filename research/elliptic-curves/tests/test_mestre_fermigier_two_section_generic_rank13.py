from __future__ import annotations

import unittest
from verify_mestre_fermigier_two_section_generic_rank13 import replay


class RejectedFermigierRank13Tests(unittest.TestCase):
    def test_rejected_generic_claim_cannot_be_reemitted(self):
        with self.assertRaisesRegex(RuntimeError, "REJECTED rank-13"):
            replay()


if __name__ == "__main__":
    unittest.main()
