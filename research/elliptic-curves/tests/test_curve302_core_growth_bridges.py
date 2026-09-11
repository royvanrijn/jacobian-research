import sys
from pathlib import Path
import unittest

CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

import curve302_core_growth_bridges as b


class CoreGrowthBridgeTests(unittest.TestCase):
    def test_quotient_line_mod_prefix(self):
        # Mod <e1>, (3,2,0) and (7,-2,0) define the same unoriented quotient line.
        q1 = b.quotient_line_mod_span(((1,0,0),), (3,2,0), 3)
        q2 = b.quotient_line_mod_span(((1,0,0),), (7,-2,0), 3)
        self.assertEqual(q1["ambient_line"], q2["ambient_line"])
        self.assertEqual(q1["ambient_line"], (0,1,0))

    def test_newly_contained_mixed_axis(self):
        # Prefix contains e1+e2. Adding e1-e2 makes both axes rationally/integrally
        # available after saturation although neither candidate is a direct axis.
        sig = b.newly_contained_axes(((1,1,0),), (1,-1,0), ("A","B","C"), 3)
        self.assertEqual(sig, ("A","B"))

    def test_core_bridge_line(self):
        prefix = ((1,0,0),)
        core = ((1,0,0),(0,1,0))
        line = b.core_bridge_line(prefix, (3,2,1), core, 3)
        # This candidate does not add a pure core line: intersection remains rank1.
        self.assertIsNone(line)
        line = b.core_bridge_line(prefix, (3,2,0), core, 3)
        self.assertEqual(line, (0,1,0))

    def test_actual_tied_vs_unique(self):
        prefix=((1,0,0,0),)
        core=((1,0,0,0),(0,1,0,0))
        positives={
            (0,1,0,0): {"chart_count":1,"first_list_order":0,"chart_ids":["a"]},
            (1,1,0,0): {"chart_count":2,"first_list_order":1,"chart_ids":["b","c"]},
            (0,0,1,0): {"chart_count":1,"first_list_order":2,"chart_ids":["d"]},
        }
        result=b.analyze_stage(seed="s",epoch=0,prefix=prefix,
            gains=({"primitive":(0,1,0,0),"word":(0,1,0,0)},), positives=positives,
            next_core=core,direction_ids=("L1","L2","L3","L4"),n=4)
        self.assertEqual(result["stage"]["best_positive_single_delta"],1)
        self.assertEqual(result["stage"]["tied_best_alternative_count"],1)
        self.assertEqual(result["stage"]["classification"],"ACTUAL_TIED_WITH_ALTERNATIVES")

    def test_shape_ignores_values(self):
        a={"kind":"x","params":[1,2],"nested":{"v":3}}
        c={"kind":"y","params":[9,8],"nested":{"v":7}}
        self.assertEqual(b.object_shape(a), b.object_shape(c))


if __name__ == "__main__":
    unittest.main()
