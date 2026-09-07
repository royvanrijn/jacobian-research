import unittest
import derivative_local_duality as derivative
import retrospective as r
import local_collision as lc
import strict_class_blocks as strict


class DerivativeBoundary(unittest.TestCase):
    def test_complete_relaxed_boundary_has_half_ambient_dimension(self):
        for row in r.read(derivative.INPUT)["cases"]:
            old = r.read(derivative.rem.bad.INPUT)["cases"][row["case_index"]]
            point_rows = [0]*old["witness_dimension"]
            beta = 0
            offset = 0
            for local,source in zip(row["local"],strict.local_rows(row["case_index"])):
                self.assertEqual(local["place"],source["prime"])
                signatures = list(map(r.pack,source["point_signature_rows"]))
                basis_old = [signatures[i] for i in local["basis_generic_indices"]]
                for i,s in enumerate(signatures):
                    coords = lc.coordinates(s,basis_old)
                    point_rows[i] |= lc.lift(coords,local["basis_signatures"])<<offset
                beta |= local["derivative_signature"]<<offset
                offset += local["width"]
            ell = sum(x["point_dimension"] for x in row["local"])
            self.assertEqual(r.rank(point_rows),ell-1)
            self.assertEqual(r.rank(point_rows+[beta]),ell)

    def test_negative_discriminant_is_resolved_at_finite_places(self):
        row = next(x for x in r.read(derivative.INPUT)["cases"] if x["case_index"]==4)
        self.assertLess(int(row["delta"]),0)
        self.assertEqual(row["beta_real_sign_bits"],[0])
        self.assertIn(3,row["outside_places"])
        self.assertNotIn("infinity",row["outside_places"])

    def test_nonmembership_checks_every_local_correction(self):
        for row in r.read(derivative.INPUT)["cases"]:
            for local in row["local"]:
                self.assertEqual(len(local["local_square_comparison_results"]),1<<local["point_dimension"])
                self.assertEqual(any(local["local_square_comparison_results"]),local["is_in_point_image"])
                if not local["is_in_point_image"]:
                    a = local["separating_character"]
                    self.assertEqual((a&local["derivative_signature"]).bit_count()%2,1)
                    self.assertTrue(all((a&v).bit_count()%2==0 for v in local["basis_signatures"]))


if __name__=="__main__":
    unittest.main()
