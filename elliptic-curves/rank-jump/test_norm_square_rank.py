import unittest
import retrospective as r
import norm_square_rank as gate
import torsion_difference as td


class NormSquareRank(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=r.read(gate.OUTPUT)

    def test_norm_one_is_not_complete_branch_splitting(self):
        c=self.data["control"]
        self.assertEqual(c["D"],1)
        self.assertEqual(c["gamma_sign_counts"]["negative"],2)
        self.assertFalse(c["gamma_is_square_in_totally_real_splitting_field"])
        self.assertEqual(c["relative_Kummer_degree_over_splitting_field"],4)
        self.assertEqual(c["specialized_Mordell_Weil_rank"],"UNKNOWN")

    def test_rank_gate_counts_reducible_components(self):
        s=self.data["surface"]
        self.assertEqual(s["finite_fibres"],["I2*"]*3)
        self.assertEqual(s["infinity"],"I0")
        self.assertEqual(s["Euler_number"],3*(2+6))
        self.assertEqual(s["trivial_lattice_rank"],2+3*6)
        self.assertEqual(s["geometric_generic_Mordell_Weil_rank"],
                         s["geometric_Neron_Severi_rank"]-s["trivial_lattice_rank"])
        self.assertEqual(s["new_generic_rank_from_norm_base_change"],0)

    def test_even_sign_group_still_has_nonscalar_unipotents(self):
        m=self.data["even_sign_module"]
        self.assertEqual(len(m["actions"]),24)
        self.assertTrue(all(row["sign_mask"].bit_count()%2==0 for row in m["actions"]))
        action=next(row for row in m["actions"] if row["permutation"]==[0,1,2] and row["sign_mask"]==3)
        M=tuple(action["matrix_columns"])
        self.assertNotEqual(M,(1,2,4,8))
        self.assertEqual(td.compose(M,M),(1,2,4,8))
        self.assertEqual(m["module"],"INDECOMPOSABLE_NONSPLIT")
        self.assertEqual(len(m["idempotents_packed"]),2)

    def test_no_specialization_rank_veto(self):
        self.assertIn("not a rank upper bound on any specialization",self.data["boundary"])
        p=r.read(gate.PROTOCOL)
        self.assertEqual(p["limits"]["point_searches"],0)
        self.assertEqual(p["limits"]["rank_computations_on_new_fibres"],0)


if __name__=="__main__":
    unittest.main()
