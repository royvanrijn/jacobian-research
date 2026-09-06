import unittest
import production_minus_twist as local
import production_twist_blocks as blocks
import retrospective as r


class ProductionTwistTests(unittest.TestCase):
    def test_intersection_against_exhaustive_small_spaces(self):
        def span(xs):
            out={0}
            for x in xs:out|={v^x for v in out}
            return out
        # Includes dependent, empty and nontrivially intersecting input columns.
        spaces=[[],[0],[1,2],[3,5],[1,2,3],[7,3],[4,5,6]]
        for a in spaces:
            for b in spaces:self.assertEqual(span(local.intersection(a,b)),span(a)&span(b))

    def test_generic_hyperbolic_pair_is_preserved(self):
        M=[2,1,8,4,0]
        pairs,radical=blocks.decompose(M,(4,8))
        self.assertEqual(pairs[0],[4,8]);self.assertEqual(len(pairs),2)
        self.assertEqual(r.rank(radical),1)

    def test_zero_form_does_not_manufacture_obstructions(self):
        pairs,radical=blocks.decompose([0,0,0])
        self.assertEqual(pairs,[]);self.assertEqual(r.rank(radical),3)

    def test_degenerate_designated_pair_is_rejected(self):
        with self.assertRaises(AssertionError):blocks.decompose([2,1,0],(1,4))


if __name__=='__main__':unittest.main()
