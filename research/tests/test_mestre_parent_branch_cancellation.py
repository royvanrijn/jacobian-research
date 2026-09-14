"""Small failure controls and exact local cancellation examples."""
import copy
from fractions import Fraction as Q
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'parent_branch_replay', ROOT/'elkies-k3/scripts/verify_mestre_parent_branch_cancellation.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
V = C.V


def quotient(a, b):
    a = list(a)
    out = [Q(0)]*max(0, len(a)-len(b)+1)
    while a and len(a) >= len(b):
        i, c = len(a)-len(b), a[-1]/b[-1]
        out[i] = c
        a = V.sub(a, [Q(0)]*i + V.scale(b, c))
    if a:
        raise ValueError('not divisible')
    return V.trim(out)


def valuation(a, f):
    assert a, 'zero polynomial has no finite valuation'
    answer = 0
    while len(a) >= len(f) and not V.remainder(a, f):
        a = quotient(a, f)
        answer += 1
    return answer


def numerator(A, B, U, W):
    s = V.add(V.power(U, 2), V.power(W, 2))
    k = V.add(V.add(V.power(U, 4), V.mul(V.power(U, 2), V.power(W, 2))), V.power(W, 4))
    bracket = V.add(V.mul(V.power(B, 2), V.power(k, 3)),
                    V.mul(V.power(A, 3), V.mul(V.mul(V.power(U, 4), V.power(W, 4)), V.power(s, 2))))
    return V.scale(V.mul(V.mul(V.mul(A, B), s), bracket), -1)


class BranchGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.packet = json.loads((C.DEFAULT/'input.json').read_text())
        cls.result = json.loads((C.DEFAULT/'result.json').read_text())

    def test_complete_replay(self):
        answer = C.verify(C.DEFAULT)
        self.assertEqual(len(answer['parents']), 4)
        self.assertFalse(answer['positive_endpoint_complete'])

    def test_source_projection_tamper(self):
        p = copy.deepcopy(self.packet['parents'][2])
        p['A'][0] = str(Q(p['A'][0])+1)
        with self.assertRaisesRegex(ValueError, 'short-model projection'):
            C.verify_parent(p, self.result['parents'][2])

    def test_wrong_local_root(self):
        r = copy.deepcopy(self.result['parents'][1])
        r['fields']['A']['local_place']['root'] = 0
        with self.assertRaisesRegex(ValueError, 'simple degree-one place'):
            C.verify_parent(self.packet['parents'][1], r)

    def test_wrong_cyclotomic_congruence(self):
        row = copy.deepcopy(self.result['parents'][1]['fields']['A'])
        row['local_place'] = {'prime': 73, 'root': 1}
        with self.assertRaisesRegex(ValueError, 'cyclotomic congruence'):
            C.verify_field(V.poly(self.packet['parents'][1]['A']), row, 'A')

    def test_missing_parent_and_missing_final_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)
            (path/'input.json').write_bytes((C.DEFAULT/'input.json').read_bytes())
            with self.assertRaises(FileNotFoundError):
                C.verify(path)
            result = copy.deepcopy(self.result)
            result['parents'].pop()
            (path/'result.json').write_text(json.dumps(result))
            with self.assertRaisesRegex(ValueError, 'result parent coverage'):
                C.verify(path)

    def test_degree_drop_and_nonsimple_place(self):
        with self.assertRaisesRegex(ValueError, 'degree drop'):
            C.M.reduction([Q(1), Q(71)], 71)
        f = [Q(1), Q(-2), Q(1)]
        self.assertFalse(C.M.irreducible(C.M.reduction(f, 71), 71))

    def test_cyclotomic_residue_alone_does_not_cancel(self):
        one, t, pi = [Q(1)], [Q(0), Q(1)], [Q(1), Q(0), Q(1)]
        u1 = t
        u2 = V.mul(t, V.add(one, V.scale(pi, Q(1, 2))))
        u3 = V.mul(t, V.add(V.add(one, V.scale(pi, Q(1, 2))), V.scale(V.power(pi, 2), Q(3, 8))))
        # At A=0, ord(u^2+1)=1 cancels, but order 2 does not.
        self.assertEqual([valuation(numerator(pi, one, u, one), pi) for u in [u1, u2]], [2, 3])
        # At B=0, the three jets give valuations 4,5,6 respectively.
        self.assertEqual([valuation(numerator(one, pi, u, one), pi) for u in [u1, u2, u3]], [4, 5, 6])

    def test_primitive_cube_residue_cancellation(self):
        one, t, pi = [Q(1)], [Q(0), Q(1)], [Q(1), Q(1), Q(1)]
        u = V.add(t, V.mul([Q(1, 3), Q(2, 3)], pi))
        k = V.add(V.add(V.power(u, 4), V.power(u, 2)), one)
        self.assertGreaterEqual(valuation(k, pi), 2)
        self.assertEqual(valuation(numerator(pi, one, u, one), pi), 4)

    def test_auxiliary_zeros_and_poles_retain_odd_branch(self):
        one, t = [Q(1)], [Q(0), Q(1)]
        for n in range(-2, 3):
            U, W = (V.power(t, n), one) if n >= 0 else (one, V.power(t, -n))
            denom_order = 14*max(0, -n)
            at_A = valuation(numerator(t, one, U, W), t)-denom_order
            at_B = valuation(numerator(one, t, U, W), t)-denom_order
            self.assertEqual(at_A, 1+14*n if n < 0 else 1)
            self.assertEqual(at_B, 3+14*n if n < 0 else (1 if n == 0 else 3))


if __name__ == '__main__':
    unittest.main()
