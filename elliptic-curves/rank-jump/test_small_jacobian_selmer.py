import unittest
from pathlib import Path
import retrospective as r
import small_jacobian_selmer as ex
import verify_small_jacobian_selmer as verify
import replay_small_jacobian_selmer as replay


class SmallJacobianSelmer(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=r.read(ex.OUTPUT)

    def test_correct_labels_change_the_positive_common_space(self):
        positive,negative=self.data["rows"]
        self.assertEqual(len(positive["common_Selmer_basis"]),3)
        self.assertEqual(len(negative["common_Selmer_basis"]),2)
        self.assertEqual(positive["pullback_of_standard_second_Selmer_basis"],[1,2,4])
        self.assertEqual(negative["pullback_of_standard_second_Selmer_basis"],[1,2,8])

    def test_complete_selmer_lifts_can_still_be_sha(self):
        positive,negative=self.data["rows"]
        self.assertEqual(negative["strict_image_dimension_from_Jacobian_Selmer"],2)
        self.assertEqual(negative["strict_image_dimension_from_rational_Jacobian_points"],0)
        self.assertEqual(negative["Jacobian_Sha_2_dimension"],4)
        self.assertEqual(positive["Jacobian_Sha_2_dimension"],0)
        self.assertEqual(positive["strict_image_dimension_from_Jacobian_Selmer"],0)

    def test_real_connecting_correction_is_retained(self):
        for row,expected in zip(self.data["rows"],[(1,1),(2,0)]):
            real=next(v for v in row["local_conditions"] if v["place"]=="infinity")
            self.assertEqual((real["connecting_rank"],real["Jacobian_local_dimension"]),expected)

    def test_halving_changes_index_without_changing_rank(self):
        positive,negative=self.data["rows"]
        self.assertEqual([positive["Jacobian_exact_rank"],negative["Jacobian_exact_rank"]],[4,2])
        self.assertEqual([positive["index_of_product_rational_image"],negative["index_of_product_rational_image"]],[2,1])

    def test_independent_sha_quotient_and_source_certificate(self):
        replay.replay()
        certificate=r.read(verify.OUTPUT)
        self.assertEqual(certificate["status"],"PASS")
        self.assertEqual(self.data["bindings"],ex.bindings())
        self.assertEqual(certificate["analysis_sha256"],r.digest(ex.OUTPUT.read_bytes()))
        self.assertEqual(certificate["verifier_sha256"],r.digest(Path(verify.__file__).read_bytes()))
        self.assertEqual([(v["left_Sha_quotient_dimension"],v["right_Sha_quotient_dimension"]) for v in certificate["rows"]],[(0,0),(2,2)])


if __name__=="__main__":
    unittest.main()
