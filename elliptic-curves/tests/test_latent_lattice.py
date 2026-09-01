from fractions import Fraction
import unittest

from latent_lattice import (
    ComponentBlock,
    EllipticCurve,
    FiniteQuotientBlock,
    ShortVectorRecord,
    build_relation_complex,
    candidate_finite_signature,
    candidate_finite_signature_from_record,
    cloud_height_profile_distance,
    cloud_height_signature,
    cross_bound_intersection_proposals,
    canonical_unoriented,
    finite_quotient_block,
    finite_rarity_scores,
    finite_rarity_weights,
    finite_signature_distance,
    hermite_signature,
    hermite_signature_distance,
    point_complexity,
    primitive_column_closure,
    rational_nullspace,
    rational_rank,
    recover_exact_embedding,
    theta_profile_distance,
    theta_signature,
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

    def test_theta_signature_is_unimodular_invariant(self) -> None:
        first = theta_signature(((2.0, 0.0), (0.0, 3.0)), minimum_vectors=8)
        # U G U^T for U=[[1,1],[0,1]].
        second = theta_signature(((5.0, 3.0), (3.0, 3.0)), minimum_vectors=8)
        self.assertEqual(first.vector_count, second.vector_count)
        self.assertLess(theta_profile_distance(first, second), 1e-12)

    def test_cloud_height_signature_is_scale_and_order_invariant(self) -> None:
        first = cloud_height_signature((2, 3, 5, 7, 11), (0, 1, 2, 3, 4), quantiles=5)
        second = cloud_height_signature((55, 35, 25, 15, 10), (4, 3, 2, 1, 0), quantiles=5)
        self.assertLess(cloud_height_profile_distance(first, second), 1e-12)

    def test_hermite_signature_is_scale_and_unimodular_invariant(self) -> None:
        first = hermite_signature(((2.0, 0.0), (0.0, 3.0)))
        second = hermite_signature(((35.0, 21.0), (21.0, 21.0)))
        self.assertLess(hermite_signature_distance(first, second), 1e-12)

    def test_finite_candidate_signature_forgets_all_basis_choices(self) -> None:
        vectors = (
            (1, 0, 0),
            (0, 1, 0),
            (1, 1, 0),
            (1, -1, 0),
            (0, 0, 1),
            (1, 0, 1),
        )
        complex_ = build_relation_complex(vectors)
        finite = FiniteQuotientBlock(
            reduction_prime=5,
            relation_prime=3,
            group_order=9,
            multiple_subgroup_order=1,
            quotient_dimension=2,
            rows=((1, 0, 1), (0, 1, 1)),
        )
        component = ComponentBlock(7, "I6", 6, (1, 2, 3), "synthetic")
        first = candidate_finite_signature(
            ((1, 0, 0), (0, 1, 0)),
            complex_,
            finite_blocks=(finite,),
            component_blocks=(component,),
        )

        # Rebase the candidate by GL(2,Z), and independently rebase the
        # quotient by GL(2,F_3).  Neither choice is part of the invariant.
        rebased_finite = FiniteQuotientBlock(
            reduction_prime=11,
            relation_prime=3,
            group_order=12,
            multiple_subgroup_order=4,
            quotient_dimension=2,
            rows=((1, 1, 2), (0, 1, 1)),
        )
        rebased = candidate_finite_signature(
            ((1, 1, 0), (0, 1, 0)),
            complex_,
            finite_blocks=(rebased_finite,),
            component_blocks=(ComponentBlock(13, "I6", 6, (5, 4, 3), "synthetic"),),
        )
        self.assertEqual(first.canonical_digest, rebased.canonical_digest)
        self.assertEqual(first.to_record(), rebased.to_record())
        self.assertEqual(finite_signature_distance(first, rebased), 0.0)
        restored = candidate_finite_signature_from_record(first.to_record())
        self.assertEqual(first, restored)

        # Change the displayed ambient basis by v -> (x,x+y,z), transforming
        # both code maps contragrediently.  This is independent of the two
        # quotient/candidate rebasings checked above.
        transformed_complex = build_relation_complex(
            tuple((x, x + y, z) for x, y, z in vectors)
        )
        ambient_changed = candidate_finite_signature(
            ((1, 1, 0), (0, 1, 0)),
            transformed_complex,
            finite_blocks=(
                FiniteQuotientBlock(17, 3, 18, 6, 2, ((1, 0, 1), (2, 1, 1))),
            ),
            component_blocks=(ComponentBlock(19, "I6", 6, (5, 2, 3), "synthetic"),),
        )
        self.assertEqual(first.canonical_digest, ambient_changed.canonical_digest)

    def test_finite_blocks_are_an_unordered_multiset_and_seed_keys_are_local(self) -> None:
        vectors = ((1, 0), (0, 1), (1, 1), (1, -1))
        complex_ = build_relation_complex(vectors)
        blocks = (
            FiniteQuotientBlock(5, 2, 4, 2, 1, ((1, 0),)),
            FiniteQuotientBlock(7, 2, 8, 4, 1, ((0, 1),)),
        )
        forward = candidate_finite_signature(
            ((1, 0), (0, 1)), complex_, finite_blocks=blocks
        )
        backward = candidate_finite_signature(
            ((1, 0), (0, 1)), complex_, finite_blocks=tuple(reversed(blocks))
        )
        self.assertEqual(forward.canonical_digest, backward.canonical_digest)
        weights = finite_rarity_weights(vectors, finite_blocks=blocks)
        self.assertEqual(len(weights), len(vectors))
        self.assertTrue(all(weight >= 1 for weight in weights))
        scores = finite_rarity_scores(vectors, finite_blocks=blocks)
        self.assertEqual(len(scores), len(vectors))
        self.assertTrue(all(score >= 0.0 for score in scores))

    def test_finite_relation_growth_is_a_proposal_mode(self) -> None:
        coordinates = ((1, 1, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        records = tuple(
            ShortVectorRecord(
                vector, "1", None, {"integral": False, "total_bits": 1}
            )
            for vector in coordinates
        )
        complex_ = build_relation_complex(coordinates)
        block = FiniteQuotientBlock(5, 2, 4, 2, 1, ((1, 0, 1),))
        proposals = independent_relation_growth_proposals(
            records,
            complex_,
            dimension=2,
            seed_edges=10,
            priority_mode="finite",
            finite_blocks=(block,),
        )
        self.assertEqual(len(proposals), 1)
        self.assertEqual(proposals[0].support, 3)

    def test_cross_bound_intersections_are_exact_and_primitive(self) -> None:
        records = tuple(
            ShortVectorRecord(vector, "1", None, {"integral": False, "total_bits": 1})
            for vector in ((1, 0, 0, 0), (0, 1, 0, 0), (1, 1, 0, 0), (0, 0, 1, 0))
        )
        from latent_lattice import GrowthProposal

        left = GrowthProposal(3, ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)), (), 0, 0, 0.0, 0.0)
        right = GrowthProposal(3, ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1)), (), 0, 0, 0.0, 0.0)
        proposals = cross_bound_intersection_proposals(
            records, (left,), (right,), target_dimension=2
        )
        self.assertEqual(len(proposals), 1)
        self.assertEqual(proposals[0].support, 3)
        self.assertEqual(proposals[0].basis_rows, ((1, 0, 0, 0), (0, 1, 0, 0)))


if __name__ == "__main__":
    unittest.main()
