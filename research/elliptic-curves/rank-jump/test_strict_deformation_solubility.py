import unittest
import retrospective as r
import local_collision as lc
import strict_deformation_solubility as strict
import strict_Sha_Artin as art


class StrictSolubility(unittest.TestCase):
    def test_each_strict_block_injects_into_the_CT_obstruction(self):
        rows=r.read(strict.OUTPUT)["single_deformations"]
        self.assertEqual([x["strict_dimension"] for x in rows],[2,1,5,3,2,4])
        for row in rows:
            K=row["strict_anchor_basis"];cross=row["CT_cross_report"]
            self.assertEqual(r.rank(K),len(K))
            self.assertEqual(r.rank(cross["cross_pairing_rows"]),len(K))
            self.assertEqual(cross["necessary_soluble_anchor_basis"],[])
            for place in row["places"]:
                self.assertTrue(all((a&b).bit_count()%2==0 for a in K for b in place["constraint_rows"]))

    def test_one_dimensional_self_pairing_does_not_erase_cross_obstruction(self):
        row=next(x for x in r.read(strict.OUTPUT)["single_deformations"] if x["u"]==-2)
        old=next(x for x in r.read(lc.INPUT)["ct"] if x["u"]==-2)
        w=row["CT_cross_report"]["strict_basis_in_inherited_coordinates"][0]
        self.assertEqual(lc.pairing(w,w,old["matrix"]),0)
        self.assertEqual(row["CT_cross_report"]["cross_pairing_rank"],1)

    def test_zero_joint_endpoint_is_preserved(self):
        data=r.read(strict.OUTPUT)
        self.assertEqual(data["all_six_common_strict_basis"],[])
        self.assertEqual(sum(x["common_strict_dimension"]>0 for x in data["pairs"]),9)
        self.assertTrue(all(x["common_necessary_soluble_basis"]==[] for x in data["pairs"]))

    def test_three_dimensional_direct_factor_and_its_dual_words(self):
        row=r.read(art.OUTPUT)["result"];M=row["Artin_matrix_rows"]
        self.assertEqual(row["elementary_S_class_direct_factor_dimension"],3)
        rows=row["selected_character_indices"];cols=row["selected_half_ideal_indices"]
        minor=[r.pack(M[i][j] for j in cols) for i in rows]
        self.assertEqual(r.rank(minor),3)
        for k,mask in enumerate(row["dual_half_ideal_words"]):
            self.assertEqual([sum(M[i][j] for j in range(5) if mask>>j&1)%2 for i in rows],
                             [int(i==k) for i in range(3)])

    def test_all_seven_selected_nonzero_characters_are_obstructed(self):
        row=r.read(art.OUTPUT)["result"]
        sr=next(x for x in r.read(strict.OUTPUT)["single_deformations"] if x["u"]==-1)
        cross=sr["CT_cross_report"]["cross_pairing_rows"]
        selected=[cross[i] for i in row["selected_character_indices"]]
        self.assertTrue(all(lc.lift(mask,selected)!=0 for mask in range(1,8)))
        self.assertEqual(row["strict_Sha_dimension"],5)


if __name__=="__main__":
    unittest.main()
