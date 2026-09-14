"""Exact allocation and frozen-selector regression; no point search."""
import json
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from height_model_benchmark import scaled_height, selection_from_portfolio


class Allocation(unittest.TestCase):
    def test_exact_ceiling(self):
        for ratio in (F(1),F(1,16),F(7,31),F(1,1000000)):
            H=scaled_height(125000,ratio)
            self.assertGreaterEqual(F(H**4),125000**4*ratio)
            self.assertLess(F((H-1)**4),125000**4*ratio)

    def test_reject_expansion(self):
        for ratio in (0,-1,2):
            with self.assertRaises(ArithmeticError):scaled_height(125000,F(ratio))

    def test_original_pair_regression(self):
        folder=ROOT/'artifacts/generated-results/elliptic-curves/pointed_height_bounds_first_chart_v1'
        paths=[folder/'preconditioned_full-bounds.json']+sorted((folder/'neighbours').glob('*-bounds.json'))
        packets=[json.loads(p.read_text()) for p in paths]
        result=selection_from_portfolio(json.loads((folder/'portfolio-verified.json').read_text()),packets)
        self.assertEqual(result['single'],0)
        self.assertEqual(set(result['pair']),{0,2})
        self.assertLess(F(result['pair_squared_address_ratio']),1)
        self.assertEqual(scaled_height(125000,F(result['D_pair'])/F(result['D_single'])),result['pair_height'])


if __name__=='__main__':unittest.main()
