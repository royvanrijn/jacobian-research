import itertools
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"cas"))
from fixed_cubic_geometry import (generic_geometry, evaluate_cubic, negative_square_section,
    alternating_rank_distribution, rational_square)

A = -5750886029903523759416717668139307
B = 167347710468055045100164888198438918505621536951206


def rank2(rows):
    basis = {}
    for word in rows:
        while word:
            pivot = word.bit_length()-1
            if pivot in basis:
                word ^= basis[pivot]
            else:
                basis[pivot] = word
                break
    return len(basis)


class GenericGeometryTests(unittest.TestCase):
    def test_pinned_anchor_has_generic_rank_zero(self):
        data = generic_geometry(A, B)
        self.assertTrue(all(data["checks"].values()))
        self.assertEqual(data["geometric_generic_rank"], 1)
        self.assertEqual(data["arithmetic_generic_rank"], 0)
        self.assertEqual(data["section_height"], "1/2")
        self.assertIn("not a bound", data["rank_scope"])

    def test_polynomial_identities_and_square_character(self):
        for a, b in itertools.product(range(-4, 5), range(-4, 6)):
            if not b or not 4*a**3+27*b*b:
                continue
            result = generic_geometry(a, b)
            self.assertEqual(result["arithmetic_generic_rank"], int(rational_square(b)))
        self.assertEqual(generic_geometry(Q(1, 2), Q(4, 9))["arithmetic_generic_rank"], 1)
        self.assertEqual(generic_geometry(Q(1, 2), Q(-4, 9))["arithmetic_generic_rank"], 0)

    def test_excluded_hypotheses(self):
        for a, b in ((-3, 2), (1, 0), (1.0, 2), (True, 2)):
            with self.assertRaises(ValueError):
                generic_geometry(a, b)

    def test_negative_square_witness_not_inherited_rank(self):
        p = negative_square_section(A, B, 1)
        self.assertEqual(p, {"u": "-1", "x": str(A+1), "y": str(A-B+1)})
        for v in (Q(2), Q(-3), Q(1, 5)):
            p = negative_square_section(A, B, v)
            self.assertEqual(evaluate_cubic(A, B, p["u"], p["x"]), Q(p["y"])**2)
        with self.assertRaises(ValueError):
            negative_square_section(A, B, 0)

    def test_parameter_is_covariant_not_model_independent(self):
        for scale in (Q(2), Q(-3), Q(1, 5)):
            for u in (Q(-2), Q(1, 3)):
                x = Q(7, 4)
                lhs = evaluate_cubic(A, B, u, scale*scale*x)
                rhs = scale**6*evaluate_cubic(A/scale**4, B/scale**6, u*scale**2, x)
                self.assertEqual(lhs, rhs)
                self.assertEqual(generic_geometry(A, B)["arithmetic_generic_rank"],
                                 generic_geometry(A/scale**4, B/scale**6)["arithmetic_generic_rank"])

    def test_alternating_distribution_by_exhaustion(self):
        for n in range(6):
            pairs = list(itertools.combinations(range(n), 2))
            counts = {}
            for mask in range(1 << len(pairs)):
                rows = [0]*n
                for bit, (i, j) in enumerate(pairs):
                    if mask >> bit & 1:
                        rows[i] ^= 1 << j
                        rows[j] ^= 1 << i
                rank = rank2(rows)
                counts[rank] = counts.get(rank, 0)+1
            self.assertEqual(alternating_rank_distribution(n),
                             {r: Q(c, 1 << len(pairs)) for r, c in counts.items()})

    def test_odd_radical_is_forced_not_positive_evidence(self):
        for n in (13, 15, 17):
            p = alternating_rank_distribution(n)
            self.assertEqual(sum(p.values()), 1)
            self.assertEqual(max(p), n-1)
            self.assertTrue(Q(83, 100) < p[n-1] < Q(85, 100))
        for n in (-1, 257, True, 3.0):
            with self.assertRaises(ValueError):
                alternating_rank_distribution(n)


if __name__ == "__main__":
    unittest.main()
