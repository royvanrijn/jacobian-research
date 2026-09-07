import unittest
import strict_artin as art
import strict_selmer_model as model


class StrictArtin(unittest.TestCase):
    def test_dual_words_detect_exactly_one_character(self):
        rows = art.r.read(art.OUTPUT)["rows"]
        for row in rows:
            if row["full_dual_basis"]:
                self.assertEqual(len(row["dual_ideal_words"]),len(row["character_masks"]))
                for i,word in enumerate(row["dual_ideal_words"]):
                    self.assertEqual(art.lc.lift(word,row["artin_columns"]),1<<i)

    def test_small_pool_does_not_erase_a_proven_global_class(self):
        row = model.calculate()["rows"][0]
        self.assertEqual(row["known_strict_dimension"],10)
        self.assertEqual(row["ideal_pool_rank"],9)
        self.assertEqual(row["characters_unseparated_by_ideal_pool"],[74])
        self.assertEqual(row["unseparated_character_point_masks"],[2441006])

    def test_remaining_incidence_is_not_a_numerical_upper_bound(self):
        for row in model.calculate()["rows"]:
            self.assertEqual(row["unknown_local_image_dimension_bound"],[0,1])
            self.assertEqual(row["witness_rank"],
                             row["known_strict_dimension"]+row["witness_local_image_dimension"])
            self.assertIn("epsilon",row["full_Selmer_dimension_formula"])
            self.assertIn("not computed",row["claim_boundary"])


if __name__ == "__main__":
    unittest.main()
