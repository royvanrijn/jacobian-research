import unittest
import retrospective as r
import completed_ordinate_blocks as blocks
import horizontal_norm_gate as gate


class CompletedOrdinateTests(unittest.TestCase):
    def test_negative_scaling_preserves_groups_up_to_sign(self):
        model=[0,11,0,-14,1];points=[(0,1),(-1,5),(2,-5)]
        transformed,mapped=blocks.change(model,points,u=-2,shift=4,s=3,t=-7)
        self.assertEqual(blocks.grouping(model,points),[[1,2]])
        self.assertEqual(blocks.grouping(transformed,mapped),[[1,2]])

    def test_raw_y_equality_can_be_manufactured(self):
        model=[0,11,0,-14,1];points=[(0,1),(-1,5)]
        transformed,mapped=blocks.change(model,points,u=1,shift=0,s=-4,t=0)
        self.assertEqual(mapped[0][1],mapped[1][1])
        self.assertEqual(blocks.grouping(transformed,mapped),[])

    def test_three_points_are_not_three_directions(self):
        groups,_=blocks.horizontal([0,-1,0,-2,1],[(0,1),(-1,1),(2,1)])
        self.assertEqual(len(groups),1);self.assertEqual(groups[0]['distinct_x_count'],3)
        self.assertTrue(groups[0]['third_point_already_listed'])

    def test_odd_valuation_gate_and_rational_denominators(self):
        self.assertEqual(gate.norm_gate(-23,[23])['status'],'EXCLUDED')
        self.assertEqual(gate.norm_gate(-r.F(1,23),[23])['status'],'EXCLUDED')
        self.assertEqual(gate.norm_gate(-r.F(1,23**2),[23])['status'],'UNKNOWN')
        self.assertEqual(gate.norm_gate(-23**2,[23])['status'],'UNKNOWN')

    def test_real_and_degenerate_exclusions(self):
        self.assertEqual(gate.norm_gate(1,[23])['status'],'EXCLUDED')
        self.assertEqual(gate.norm_gate(0,[23])['status'],'EXCLUDED')


if __name__=='__main__':unittest.main()
