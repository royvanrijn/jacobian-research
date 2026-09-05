"""Independent small-lattice checks and fail-closed selection regressions."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import json
import sys
from tempfile import TemporaryDirectory
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.cvp import Hole, VoronoiIterator
from research_runtime.deep_centres import norm, parity, exact_coset_minimum, diverse_deep


class DeepCentreTests(unittest.TestCase):
    def test_exact_parity_minima_against_independent_box(self):
        gram = [[4, 3, -1], [3, 6, 1], [-1, 1, 4]]
        for mask in range(8):
            vectors = [z for z in product(range(-5,6), repeat=3) if parity(z) == mask]
            expected = min(norm(gram,z) for z in vectors)
            seed = tuple(((mask >> i) & 1)+4 for i in range(3))
            actual, rep, _ = exact_coset_minimum(gram, mask, seed)
            self.assertEqual(actual, expected)
            self.assertEqual(norm(gram, rep), expected)

    def test_nearest_policy_does_not_select_deepest_classes(self):
        gram = [[4,0,0], [0,4,0], [0,0,4]]
        nearest = VoronoiIterator(gram).next_holes(2, diversity_window=2)
        self.assertTrue(all(r.squared_distance < 3 for r in nearest))
        deepest = exact_coset_minimum(gram,7,(1,1,1))
        self.assertEqual(deepest[0],12)

    def test_torus_diversity_is_sign_invariant(self):
        rows = [Hole(m,tuple((m >> i)&1 for i in range(3)),Fraction(m.bit_count())) for m in (3,5,6,7)]
        distance = lambda mask: Fraction(mask.bit_count())
        first = diverse_deep(rows,distance,3)
        second = diverse_deep([Hole(r.mask,tuple(-x for x in r.doubled_coordinates),r.squared_distance) for r in rows],distance,3)
        self.assertEqual([r.mask for r in first],[r.mask for r in second])
        self.assertEqual(first[0].mask,7)

    def test_bad_witness_and_budget_fail_closed(self):
        with self.assertRaises(ValueError): exact_coset_minimum([[4,1],[1,4]],3,(1,0))
        with self.assertRaises(TimeoutError): exact_coset_minimum([[4,1],[1,4]],3,(9,9),node_budget=1)

    def test_recovery_gate_rejects_missing_and_incomplete_cells(self):
        from run_mw18_centre_experiment import summarize
        protocol={'protocol_hash':'test','cases':{f'c{i}':{'anchor_id':f'curve{i}'} for i in range(5)},
            'policies':['nearest_first','deepest','diverse_deep'],
            'success_gate':{'minimum_total_gain':35,'minimum_each_anchor_gain':5}}
        with TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'cells').mkdir()
            for case in protocol['cases']:
                for policy in protocol['policies']:
                    row={'protocol_hash':'test','status':'COMPLETE','certified_gain':7 if policy!='nearest_first' else 0,
                         'charts':[{'search':{'status':'bounded_search_complete'}}]*40}
                    (root/'cells'/f'{case}--{policy}.json').write_text(json.dumps(row))
            self.assertTrue(summarize(protocol,root)['success_gate_passed'])
            path=root/'cells/c4--diverse_deep.json';saved=path.read_text();path.unlink()
            self.assertFalse(summarize(protocol,root)['success_gate_passed'])
            row=json.loads(saved);row['status']='INCOMPLETE_BOXES';path.write_text(json.dumps(row))
            self.assertFalse(summarize(protocol,root)['success_gate_passed'])


if __name__ == '__main__': unittest.main()
