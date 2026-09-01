from fractions import Fraction
import unittest

from latent_lattice import (
    EllipticCurve,
    ShortVectorRecord,
    build_relation_complex,
    canonical_unoriented,
    finite_quotient_block,
    point_complexity,
    primitive_column_closure,
    rational_nullspace,
    rational_rank,
    recover_exact_embedding,
)
from latent_lattice.subspace import (
    exact_span_mask,
    exact_row_space_intersection,
    independent_row_basis,
    independent_relation_growth_proposals,
    modular_row_space_key,
    modular_right_kernel_basis,
    primitive_span_basis,
)


class LatentLatticeTests(unittest.TestCase):
    def test_canonical_unoriented_is_primitive(self) -> None:
        self.assertEqual(canonical_unoriented((-2, 4, 0)), (1, -2, 0))
        self.assertEqual(canonical_unoriented((2, -4, 0)), (1, -2, 0))
        with self.assertRaises(ValueError):
            canonical_unoriented((0, 0))

    def test_exact_rank(self) -> None:
        self.assertEqual(rational_rank(((1, 2, 3), (2, 4, 6), (0, 1, 0))), 2)
        kernel = rational_nullspace(((1, 2, 3), (0, 1, 0)))
        self.assertEqual(len(kernel), 1)
        self.assertEqual(sum(a * b for a, b in zip((1, 2, 3), kernel[0])), 0)
        self.assertEqual(sum(a * b for a, b in zip((0, 1, 0), kernel[0])), 0)

    def test_general_weierstrass_group_law(self) -> None:
        curve = EllipticCurve((0, 0, 0, -1, 0))
        point = (Fraction(0), Fraction(0))
        self.assertTrue(curve.is_on_curve(point))
        self.assertEqual(curve.multiply(point, 2), None)
        self.assertEqual(curve.linear_combination((point,), (3,)), point)
        self.assertTrue(point_complexity(point)["integral"])

    def test_relation_complex_is_order_and_sign_independent(self) -> None:
        first = build_relation_complex(((1, 0), (0, 1), (1, 1), (1, -1)))
        second = build_relation_complex(((-1, 1), (-1, -1), (0, -3), (-2, 0)))
        self.assertEqual(first.vertices, second.vertices)
        self.assertEqual(first.ternary_relations, second.ternary_relations)
        self.assertEqual(first.canonical_digest, second.canonical_digest)
        self.assertGreater(len(first.ternary_relations), 0)
        self.assertIn(2, {edge[3] for edge in first.scaled_relations})

        # A nontrivial unimodular change of ambient basis must leave the
        # coordinate-free invariant unchanged.
        transformed = build_relation_complex(
            tuple((x + 2 * y, x + y) for x, y in first.vertices)
        )
        self.assertEqual(first.canonical_digest, transformed.canonical_digest)
        self.assertEqual(first.wl_profile(), transformed.wl_profile())

    def test_general_finite_quotient_code_is_a_homomorphism(self) -> None:
        curve = EllipticCurve((0, 0, 0, -1, 1))
        point = (Fraction(0), Fraction(1))
        double = curve.multiply(point, 2)
        block = finite_quotient_block(curve, (point, double), 5, 2)
        self.assertIn(block.quotient_dimension, (0, 1, 2))
        point_class = block.vector_class((1, 0))
        double_class = block.vector_class((0, 1))
        self.assertEqual(double_class, tuple(0 for _ in point_class))

    def test_exact_span_and_primitive_basis(self) -> None:
        basis = ((2, 0, 0), (0, 3, 0))
        primitive = primitive_span_basis(basis)
        mask = exact_span_mask(((1, 1, 0), (0, 0, 1), (4, -3, 0)), primitive)
        self.assertEqual(rational_rank(primitive), 2)
        self.assertEqual(mask.tolist(), [True, False, True])
        self.assertEqual(
            independent_row_basis(((1, 0, 0), (2, 0, 0), (0, 3, 0))),
            ((1, 0, 0), (0, 3, 0)),
        )
        self.assertEqual(
            modular_row_space_key(((1, 2, 0), (0, 1, 0))),
            modular_row_space_key(((1, 0, 0), (3, 1, 0))),
        )
        self.assertEqual(
            modular_right_kernel_basis(((1, 0, 0), (0, 1, 0))),
            ((0, 0, 1),),
        )
        intersection = exact_row_space_intersection(
            ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)),
            ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1)),
        )
        self.assertEqual(rational_rank(intersection), 2)
        self.assertTrue(
            exact_span_mask(((1, 0, 0, 0), (0, 1, 0, 0)), intersection).all()
        )

    def test_growth_proposal_indices_use_input_record_order(self) -> None:
        coordinates = ((0, 0, 1), (1, 1, 0), (1, 0, 0), (0, 1, 0))
        records = tuple(
            ShortVectorRecord(
                vector,
                str(index + 1),
                None,
                {
                    "integral": index % 2 == 0,
                    "total_bits": index + 1,
                    "x_numerator_bits": index + 1,
                    "x_denominator_bits": 1,
                    "y_numerator_bits": index + 2,
                    "y_denominator_bits": 1,
                },
            )
            for index, vector in enumerate(coordinates)
        )
        complex_ = build_relation_complex(coordinates)
        proposals = independent_relation_growth_proposals(
            records, complex_, dimension=2, seed_edges=10
        )
        self.assertEqual(len(proposals), 1)
        self.assertEqual(proposals[0].inlier_indices, (1, 2, 3))
        relation_only = independent_relation_growth_proposals(
            records,
            complex_,
            dimension=2,
            seed_edges=10,
            priority_mode="relations",
            seed_strategy="stratified",
        )
        self.assertEqual(relation_only[0].inlier_indices, (1, 2, 3))

    def test_height_dual_embedding_has_exact_replay(self) -> None:
        curve = EllipticCurve((0, 0, 0, -1, 0))
        # This rank-one subgroup has no torsion ambiguity after choosing a
        # non-torsion point on y^2=x^3-x+1.
        curve = EllipticCurve((0, 0, 0, -1, 1))
        point = (Fraction(0), Fraction(1))
        target = curve.multiply(point, 3)
        embedding = recover_exact_embedding(curve, (point,), (target,), digits=60)
        self.assertEqual(embedding.columns, ((3,),))
        self.assertEqual(embedding.smith_invariant_factors(), (3,))
        self.assertFalse(embedding.is_primitive())
        self.assertEqual(primitive_column_closure(((3,), (6,))), ((1,), (2,)))


if __name__ == "__main__":
    unittest.main()
