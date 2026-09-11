import sys
from pathlib import Path
import unittest

CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

import curve302_closure_transition_graph as g


class TransitionGraphTests(unittest.TestCase):
    def test_rref_identifies_same_saturated_state(self):
        a = g.canonical_rref(((1,0,0),(0,1,0)), 3)
        b = g.canonical_rref(((1,1,0),(1,-1,0)), 3)
        self.assertEqual(a, b)

    def test_intersection(self):
        u = ((1,0,0),(0,1,0))
        v = ((1,1,0),(0,0,1))
        inter = g.intersection_rref(u, v, 3)
        self.assertEqual(len(inter), 1)
        self.assertEqual(g.matrix_rank(inter + ((1,1,0),), 3), 1)

    def test_stable_id_order_independent_for_payload_dict(self):
        self.assertEqual(g.stable_id("X", {"a":1,"b":2}), g.stable_id("X", {"b":2,"a":1}))

    def test_small_transition_graph_tied_and_unique(self):
        data = {
            "cores": {
                1: ((1,0,0),),
                2: ((1,0,0),(0,1,0)),
                3: ((1,0,0),(0,1,0),(0,0,1)),
            },
            "runs": [
                {"seed":"s0","seed_index":0,"final_dimension":3,"stages":[
                    {"epoch":0,"gains":[{"word":(0,1,0),"primitive":(0,1,0)}]},
                    {"epoch":1,"gains":[{"word":(0,0,1),"primitive":(0,0,1)}]},
                ]},
            ],
        }
        ledger = {"direction_ids":["A","B","C"],"runs":[]}
        bridge = {
            "status":"PASS_27_CORE_GROWTH_BRIDGE_ANALYSIS",
            "stages":[
                {"seed":"s0","epoch":0,"pre_dimension":1,"post_dimension":2,"classification":"ACTUAL_TIED_WITH_ALTERNATIVES","current_next_core_intersection_rank":1,"actual_batch_next_core_delta":1},
            ],
            "bridge_candidates":[
                {"seed":"s0","epoch":0,"candidate_word":[0,1,0],"actual_member":True,"tied_best_positive":True,"quotient_line_mod_prefix":[0,1,0],"core_bridge_line":[0,1,0],"new_axis_signature":["B"],"chart_count_lower_bound":1},
                {"seed":"s0","epoch":0,"candidate_word":[1,1,0],"actual_member":False,"tied_best_positive":True,"quotient_line_mod_prefix":[0,1,0],"core_bridge_line":[0,1,0],"new_axis_signature":["B"],"chart_count_lower_bound":1},
            ],
        }
        out = g.transition_graph(data, ledger, bridge, 3)
        self.assertEqual(out["summary"]["decisive_stages"], 1)
        # Two rational points but one quotient bridge line / one core transition.
        self.assertEqual(out["summary"]["common_core_transition_types"], 1)
        self.assertEqual(out["core_transition_types"][0]["distinct_bridge_lines_mod_prefix"], 1)
        self.assertEqual(out["core_transition_types"][0]["bottleneck_class"], "NARROW_OBSERVED_GATE")

    def test_rank29_style_control_detects_missing_late_gate(self):
        # 4D toy analogue: one run stalls saturated at rank2, another reaches rank4.
        data = {
            "cores": {
                1: ((1,0,0,0),),
                2: ((1,0,0,0),(0,1,0,0)),
                3: ((1,0,0,0),(0,1,0,0),(0,0,1,0)),
                4: ((1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)),
            },
            "runs": [
                {"seed":"stall","seed_index":0,"final_dimension":2,"stages":[
                    {"epoch":0,"gains":[{"word":(0,1,0,0),"primitive":(0,1,0,0)}]},
                    {"epoch":1,"gains":[],"terminal_no_gain":True},
                ]},
                {"seed":"ok","seed_index":0,"final_dimension":4,"stages":[
                    {"epoch":0,"gains":[{"word":(0,1,0,0),"primitive":(0,1,0,0)}]},
                    {"epoch":1,"gains":[{"word":(0,0,1,0),"primitive":(0,0,1,0)}]},
                    {"epoch":2,"gains":[{"word":(0,0,0,1),"primitive":(0,0,0,1)}]},
                ]},
            ],
        }
        ledger = {"direction_ids":["A","B","C","D"],"runs":[
            {"seed":"stall","stages":[
                {"epoch":0,"charts":[]}, {"epoch":1,"charts":[]},
            ]},
            {"seed":"ok","stages":[
                {"epoch":0,"charts":[]}, {"epoch":1,"charts":[]}, {"epoch":2,"charts":[]},
            ]},
        ]}
        graph={"bridge_edges":[
            {"seed":"ok","epoch":1,"pre_dimension":2,"candidate_word":[0,0,1,0],"actual_member":True},
            {"seed":"ok","epoch":2,"pre_dimension":3,"candidate_word":[0,0,0,1],"actual_member":True},
        ]}
        out=g.rank29_control(data,ledger,{},graph,4,stall_dimension=2,next_dimension=3,final_dimension=4)
        self.assertEqual(out["dimension13_deficit"],1)
        self.assertEqual(out["dimension14_deficit"],2)
        self.assertEqual(out["terminal_positive_new_directions"],0)
        self.assertEqual(out["first_missing_gate_status"],"NO_POSITIVE_TERMINAL_C13_BRIDGE")
        self.assertGreaterEqual(out["successful_late_bridge_line_families_mod_stall"],1)


if __name__ == "__main__":
    unittest.main()
