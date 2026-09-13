"""Adversarial controls for the complete projective-pencil exclusion."""
import copy
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("correlated_replay", ROOT / "elkies-k3/scripts/verify_r17_correlated_genus_one_pencils.py")
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)


class CorrelatedPencilReplay(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = V.DEFAULT
        cls.packet, cls.pencils, cls.result = [json.loads((path/name).read_text())
                                             for name in ("input.json", "pencils.json", "result.json")]

    def test_altered_trace_ordinate_is_rejected(self):
        pencils = copy.deepcopy(self.pencils)
        pencils[0]["Ny"][0] = str(V.Q(pencils[0]["Ny"][0])+1)
        with self.assertRaisesRegex(ValueError, "trace y group-law"):
            V.verify(self.packet,pencils,self.result,False)

    def test_altered_branch_is_rejected(self):
        pencils = copy.deepcopy(self.pencils)
        pencils[0]["branch_matrix"][0][0] = str(V.Q(pencils[0]["branch_matrix"][0][0])+1)
        with self.assertRaisesRegex(ValueError, "branch coefficient"):
            V.verify(self.packet,pencils,self.result,False)

    def test_missing_pair_is_rejected(self):
        result = copy.deepcopy(self.result)
        result["modular_witnesses"][0]["rejected_pairs"].pop()
        with self.assertRaisesRegex(ValueError, "exclusion witness"):
            V.verify(self.packet,self.pencils,result,False)

    def test_duplicate_pencil_is_rejected(self):
        pencils = copy.deepcopy(self.pencils)
        pencils[1] = pencils[0]
        with self.assertRaisesRegex(ValueError, "trace attachment"):
            V.verify(self.packet,pencils,self.result,False)

    def test_singular_reduction_cannot_exclude(self):
        a = [[V.Q(int(i==j)) for j in range(5)] for i in range(5)]
        a[4][4] = V.Q(101)
        self.assertIsNone(V.reduce_matrix(a,101))
        self.assertIsNotNone(V.reduce_matrix(a,103))

    def test_constant_scaling_keeps_same_projective_cover(self):
        a = [[int(i==j) for j in range(5)] for i in range(5)]
        b = [[3*x for x in row] for row in a]
        self.assertEqual(V.projective_image(a,101),V.projective_image(b,101))

    def test_projective_infinity_is_included(self):
        a = [[int(i==j) for j in range(5)] for i in range(5)]
        image = V.projective_image(a,101)
        self.assertIn((0,0,0,0,1),image)
        self.assertEqual(len(image),102)


if __name__ == "__main__":
    unittest.main()
