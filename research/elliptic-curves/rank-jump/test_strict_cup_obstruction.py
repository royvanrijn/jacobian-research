import unittest
from itertools import permutations
import retrospective as r
import affine_selmer as af
import strict_cup_obstruction as cup
import strict_Sha_Artin as art
import torsion_difference as td


class StrictCupObstruction(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = r.read(cup.OUTPUT)
        cls.cases = {row["u"]: row for row in cls.report["cases"]}

    def test_cup_action_matches_six_branch_subsets(self):
        states = {td.encode_state(k, v): (k, v)
                  for k in range(4) for v in range(4)}
        for perm in permutations(range(3)):
            g = (perm[0]+1, perm[1]+1)
            for signs in range(8):
                for subset, (k, v) in states.items():
                    moved = states[td.permute_subset(subset, perm, signs)]
                    gv = td.act(g, v)
                    self.assertEqual(moved[1], gv)
                    self.assertEqual(cup.embed(moved[0] ^ td.act(g, k)),
                                     cup.projection(signs & cup.embed(gv)))

    def test_complete_jacobian_support_includes_parameter_primes(self):
        inputs = r.read(af.INPUT)["cases"]
        for u, row in self.cases.items():
            raw = next(x for x in inputs if x["u"] == u)
            places = {x["place"] for x in raw["local"]}
            self.assertTrue({2, "infinity"}.issubset(places))
            if abs(u) > 1:
                self.assertIn(abs(u), places)
            norm = 1
            for p, e in row["norm_support"]:
                self.assertIn(p, places)
                norm *= p**e
            self.assertEqual(norm, abs(int(row["norm_gamma"])))

    def test_partial_detection_does_not_certify_zero(self):
        for u in (-3, -2):
            row = self.cases[u]
            self.assertEqual(row["detected_cup_image_dimension_lower_bound"], 0)
            self.assertGreater(row["full_CT_cross_rank"], 0)
        row = self.cases[-1]
        self.assertEqual(row["detected_cup_image_dimension_lower_bound"], 4)
        self.assertEqual(row["retained_character_annihilator_anchor_masks"], [17108])
        self.assertEqual(row["full_CT_cross_rank"], 5)
        self.assertIn("does not prove", self.report["boundary"])

    def test_cup_projection_differs_from_ordinary_half_ideal(self):
        row = self.report["u_minus_one_elementary_factor_projection"]
        M = r.read(art.OUTPUT)["result"]["Artin_matrix_rows"]
        selected = row["selected_character_indices"]
        words = row["projected_cup_half_ideal_words"]
        self.assertEqual(words, [0, 4, 16, 0, 16])
        for i, word in enumerate(words):
            detected = [sum(M[j][k] for k in range(5) if word >> k & 1) % 2
                        for j in selected]
            self.assertEqual(detected,
                             [self.cases[-1]["strict_CT_matrix"][i][j] for j in selected])
        # The first ordinary half ideal is detected; its cup projection is zero.
        self.assertTrue(any(M[j][0] for j in selected))
        self.assertEqual(words[0], 0)
        self.assertEqual(row["cup_image_in_complement_dimension_lower_bound"], 2)

    def test_inferred_values_remain_distinct_from_new_arithmetic(self):
        protocol = r.read(cup.PROTOCOL)
        for key in ("new_parameters", "point_searches", "class_groups",
                    "norm_equation_searches", "new_CT_calculations"):
            self.assertEqual(protocol["limits"][key], 0)
        projection = self.report["u_minus_one_elementary_factor_projection"]
        self.assertIn("No independent norm-witness evaluation", projection["coordinate_source"])
        self.assertIn("remain uncomputed", self.report["boundary"])


if __name__ == "__main__":
    unittest.main()
