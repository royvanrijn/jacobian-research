import unittest
import strict_class_blocks as blocks
import remaining_bad_primes as rem


class StrictClassBlocks(unittest.TestCase):
    def test_separate_surjectivity_does_not_imply_joint_surjectivity(self):
        # Two generic directions, two independent global exceptional directions;
        # both exceptional local images agree, and their sum is locally zero.
        rows = [{"prime":3,"point_signature_rows":[[1],[0],[1],[1]]},
                {"prime":5,"point_signature_rows":[[1],[0],[0],[0]]}]
        kernel,columns,_ = blocks.kernel(rows,4)
        self.assertEqual(rem.r.rank(kernel),2)
        self.assertTrue(all(rem.r.reduce(x,rem.r.basis(kernel))==0 for x in (2,12)))
        self.assertTrue(all((c&12).bit_count()%2==0 for c in columns))

    def test_complete_r17_pair_has_distinct_generic_baselines(self):
        high,low = blocks.calculate(4),blocks.calculate(5)
        self.assertEqual((high["generic_unramified_character_dimension"],
                          high["relative_unramified_character_dimension"]),(2,6))
        self.assertEqual((low["generic_unramified_character_dimension"],
                          low["relative_unramified_character_dimension"]),(6,0))
        self.assertEqual(high["joint_quotient_checks"][0]["exceptional_value_mask"],34)

    def test_incomplete_factorizations_do_not_claim_class_group_bounds(self):
        for i in (1,2,3):
            row = blocks.calculate(i)
            self.assertFalse(row["all_bad_places_complete"])
            self.assertIsNone(row["ordinary_class_group_two_rank_lower_bound"])
            self.assertIsNone(row["relative_unramified_character_dimension"])
            self.assertIsNone(row["unramified_compositum_degree"])


if __name__ == "__main__":
    unittest.main()
