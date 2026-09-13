"""Small failure controls for the exact one-node proof obligations."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location("one_node_replay",ROOT/"elkies-k3/scripts/verify_r17_one_node_correlated.py")
V=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class OneNodeReplay(unittest.TestCase):
    def test_root_at_infinity_blocks_affine_miss(self):
        # 101*x-1 has no affine root mod101, but has the rational root1/101.
        self.assertTrue(V.has_projective_root([-1,101],101))

    def test_composite_witness_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"nonprime"):
            V.has_projective_root([1,0,1],9)

    def test_content_and_denominators_are_removed_before_reduction(self):
        self.assertEqual(V.primitive([V.Q(202,3),V.Q(0),V.Q(202,3)]),[1,0,1])
        self.assertFalse(V.has_projective_root([1,0,1],3))

    def test_nonexact_node_division_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"nonexact"):
            V.divide_exact([V.Q(1),V.Q(0),V.Q(1)],[V.Q(0),V.Q(0),V.Q(1)])

    def test_split_control_is_inherited_and_tamper_is_rejected(self):
        packet=json.loads((V.DEFAULT/"input.json").read_text())
        record=V.read_traces(V.DEFAULT)[54]
        carrier=next(r for r in record["records"] if r["status"]=="RATIONAL_HALF_FOUND")["carriers"][0]
        basis=[(V.V.poly(p["x"]),V.V.poly(p["y"])) for p in packet["generic_source"]["basis"]]
        self.assertEqual(V.identify_split_branches(packet["traces"][54],carrier,basis),["P1+P3","-P1-P16"])
        bad=copy.deepcopy(carrier);bad["q"][0]=str(V.Q(bad["q"][0])+1)
        with self.assertRaisesRegex(ValueError,"split sextic"):
            V.identify_split_branches(packet["traces"][54],bad,basis)


if __name__=="__main__":unittest.main()
