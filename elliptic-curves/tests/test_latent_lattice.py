from fractions import Fraction
import unittest

from latent_lattice import (
    ComponentBlock,
    EllipticCurve,
    FiniteQuotientBlock,
    ShortVectorRecord,
    aggregate_repeated_intersection_ledgers,
    build_relation_complex,
    biconnected_components,
    bounded_dense_hyperplanes,
    bounded_metric_relation_search,
    bounded_metric_star_component_search,
    candidate_finite_signature,
    candidate_finite_signatures,
    candidate_relation_fingerprint,
    candidate_finite_signature_from_record,
    cloud_height_profile_distance,
    cloud_height_signature,
    cross_bound_intersection_proposals,
    exact_partial_relation_replay,
    exact_intersection_consensus,
    exact_graph_walk_consensus,
    exact_rational_ranks,
    canonical_unoriented,
    canonical_rational_unoriented,
    finite_quotient_block,
    finite_rarity_scores,
    finite_rarity_weights,
    finite_signature_distance,
    dense_two_core_component,
    exact_rational_space_key,
    hermite_signature,
    hermite_signature_distance,
    height_angle_profile,
    height_angle_profile_distance,
    intrinsic_shell_signature,
    intrinsic_shell_profile_distance,
    joint_nearest_candidate_scores,
    lift_relation_vertex_bijection,
    lift_relation_vertex_injection,
    maximal_star_components,
    merge_component_vertex_maps,
    modular_rank,
    overlapping_star_components,
    replay_and_deduplicate_components,
    relation_metric_profile_distance,
    relation_metric_signature,
    point_complexity,
    partial_replay_finite_signature,
    primitive_column_closure,
    primitive_hermite_signatures,
    rational_nullspace,
    rational_rank,
    row_basis_coordinates,
    row_embedding_is_primitive,
    row_embedding_saturation_indices,
    row_embedding_smith_invariant_factors,
    recover_exact_embedding,
    repeated_cross_bound_intersection_ledger,
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
        self.assertEqual(
            canonical_rational_unoriented((Fraction(1, 2), Fraction(-3, 4))),
            (2, -3),
        )

    def test_exact_rank(self) -> None:
        self.assertEqual(rational_rank(((1, 2, 3), (2, 4, 6), (0, 1, 0))), 2)
        kernel = rational_nullspace(((1, 2, 3), (0, 1, 0)))
        self.assertEqual(len(kernel), 1)
        self.assertEqual(sum(a * b for a, b in zip((1, 2, 3), kernel[0])), 0)
        self.assertEqual(sum(a * b for a, b in zip((0, 1, 0), kernel[0])), 0)
        coordinates = row_basis_coordinates(
            ((3, 5, 8), (-1, 2, 1)),
            ((1, 1, 2), (0, 1, 1)),
        )
        self.assertEqual(coordinates, ((3, 2), (-1, 3)))
        with self.assertRaises(ValueError):
            row_basis_coordinates(((1, 0),), ((2, 0),))
        self.assertEqual(modular_rank(((1, 0), (0, 2)), 2), 1)
        self.assertEqual(modular_rank(((1, 0), (0, 2)), 3), 2)

    def test_exact_component_global_replay_and_deduplication(self) -> None:
        vectors = ((1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1))
        complex_ = build_relation_complex(vectors)
        self.assertEqual(
            exact_rational_space_key(((1, 0, 0), (0, 1, 0))),
            exact_rational_space_key(((1, 1, 0), (1, -1, 0))),
        )
        replayed = replay_and_deduplicate_components(
            vectors,
            complex_,
            (
                ((1, 0, 0), (0, 1, 0)),
                ((1, 1, 0), (1, -1, 0)),
                ((1, 0, 0), (0, 0, 1)),
            ),
            development_indices=(0, 1),
            held_out_indices=(2, 3),
        )
        self.assertEqual(len(replayed), 2)
        plane = next(item for item in replayed if len(item.origin_indices) == 2)
        self.assertEqual(plane.origin_indices, (0, 1))
        self.assertEqual(plane.development_replayed_ray_indices, (0, 1))
        self.assertEqual(plane.held_out_replayed_ray_indices, (2,))
        self.assertEqual(len(plane.full_replayed_relation_indices), 1)
        self.assertEqual(dict(plane.modular_ranks), {2: 2, 3: 2})

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
        target_index = {vector: index for index, vector in enumerate(transformed.vertices)}
        vertex_map = tuple(
            target_index[canonical_unoriented((x + 2 * y, x + y))]
            for x, y in first.vertices
        )
        lifts = lift_relation_vertex_bijection(first, transformed, vertex_map)
        self.assertIn(((1, 1), (2, 1)), lifts)
        self.assertIn(((-1, -1), (-2, -1)), lifts)
        larger = build_relation_complex(transformed.vertices + ((2, 1),))
        larger_index = {vector: index for index, vector in enumerate(larger.vertices)}
        injection = tuple(
            larger_index[canonical_unoriented((x + 2 * y, x + y))]
            for x, y in first.vertices
        )
        injected_lifts = lift_relation_vertex_injection(first, larger, injection)
        self.assertIn(((1, 1), (2, 1)), injected_lifts)

        rectangular = build_relation_complex(
            tuple((x + 2 * y, x + y, x - y) for x, y in first.vertices)
            + ((0, 0, 1),)
        )
        rectangular_index = {
            vector: index for index, vector in enumerate(rectangular.vertices)
        }
        rectangular_map = tuple(
            rectangular_index[
                canonical_unoriented((x + 2 * y, x + y, x - y))
            ]
            for x, y in first.vertices
        )
        rectangular_lifts = lift_relation_vertex_injection(
            first, rectangular, rectangular_map
        )
        expected_rectangular = ((1, 1, 1), (2, 1, -1))
        self.assertIn(expected_rectangular, rectangular_lifts)
        self.assertEqual(
            row_embedding_smith_invariant_factors(expected_rectangular), (1, 1)
        )
        self.assertTrue(row_embedding_is_primitive(expected_rectangular))
        self.assertFalse(row_embedding_is_primitive(((2, 0, 0), (0, 1, 0))))
        self.assertEqual(
            row_embedding_saturation_indices(
                (((1, 0, 0), (0, 1, 0)), ((2, 0, 0), (0, 1, 0)))
            ),
            (1, 2),
        )
        self.assertEqual(
            exact_rational_ranks((((1, 0), (0, 1)), ((1, 2), (2, 4)))),
            (2, 1),
        )
        with self.assertRaises(ValueError):
            row_embedding_saturation_indices((((1, 0), (2, 0)),))

        # A rank-one mapped component is saturated and replayed without an
        # arbitrary completion to the rank-two source ambient lattice.
        partial_map = [-1] * len(first.vertices)
        for source_index, vector in enumerate(first.vertices):
            if vector[1] == 0:
                partial_map[source_index] = rectangular_index[
                    canonical_unoriented((vector[0], vector[0], vector[0]))
                ]
        partial = exact_partial_relation_replay(first, rectangular, partial_map)
        self.assertTrue(partial)
        self.assertEqual(partial[0].source_rank, 1)
        self.assertEqual(partial[0].source_subspace_ray_count, 1)
        self.assertEqual(len(partial[0].integral_matrix), 1)
        self.assertEqual(partial[0].replayed_source_rank, 1)
        self.assertTrue(partial[0].primitive_target_image)
        finite_signature = partial_replay_finite_signature(
            partial[0], rectangular
        )
        self.assertEqual(finite_signature.candidate_dimension, 1)
        self.assertEqual(finite_signature.retained_ray_count, 1)
        self.assertEqual(
            candidate_finite_signatures(
                (expected_rectangular,), rectangular
            ),
            (candidate_finite_signature(expected_rectangular, rectangular),),
        )

        metric_search = bounded_metric_relation_search(
            first,
            rectangular,
            ((1.0, 0.0), (0.0, 1.0)),
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 3.0)),
            source_center_limit=4,
            target_center_limit=7,
            seed_edges_per_center=4,
            beam_width=100,
            maximum_steps=10,
            norm_log_tolerance=3.0,
            angle_tolerance=1.0,
            angle_hard_tolerance=1.0,
        )
        self.assertTrue(metric_search.embeddings)
        self.assertGreaterEqual(
            metric_search.embeddings[0].global_replay_ray_count, 3
        )
        self.assertEqual(metric_search.embeddings[0].global_replay_source_rank, 2)
        self.assertTrue(
            any(
                row_embedding_is_primitive(matrix)
                for embedding in metric_search.embeddings
                for matrix in embedding.integral_matrices
            )
        )

        rank_three = build_relation_complex(
            ((1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1))
        )
        rank_four = build_relation_complex(
            tuple((*vector, 0) for vector in rank_three.vertices)
            + ((0, 0, 0, 1),)
        )
        source_center = rank_three.vertices.index((1, 0, 0))
        target_center = rank_four.vertices.index((1, 0, 0, 0))
        star = bounded_metric_star_component_search(
            rank_three,
            rank_four,
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            (
                (1.0, 0.0, 0.0, 0.0),
                (0.0, 1.0, 0.0, 0.0),
                (0.0, 0.0, 1.0, 0.0),
                (0.0, 0.0, 0.0, 2.0),
            ),
            source_center_indices=(source_center,),
            target_center_indices=(target_center,),
            source_edge_limit=2,
            beam_width=64,
            maximum_partial_replay_attempts=16,
            minimum_partial_replay_rank=2,
            partial_replay_candidate_limit=16,
            minimum_partial_replay_rays=3,
            norm_log_tolerance=1.0,
            angle_tolerance=1.0,
            angle_hard_tolerance=1.0,
        )
        self.assertTrue(star.candidates)
        self.assertGreaterEqual(star.maximum_partial_replay_ray_count, 3)
        self.assertTrue(star.candidates[0].replay.primitive_target_image)

        stars = maximal_star_components(first)
        self.assertTrue(stars)
        self.assertEqual(stars[0].mod2_rank, 2)
        self.assertTrue(overlapping_star_components(first))
        self.assertEqual(
            dense_two_core_component(first, minimum_relation_degree=1).rational_rank,
            2,
        )
        self.assertTrue(biconnected_components(first))

        plane_heavy = build_relation_complex(
            (
                (1, 0, 0),
                (0, 1, 0),
                (1, 1, 0),
                (1, -1, 0),
                (2, 1, 0),
                (0, 0, 1),
            )
        )
        hyperplanes = bounded_dense_hyperplanes(
            plane_heavy,
            sample_count=100,
            random_seed=7,
            maximum_components=8,
        )
        self.assertEqual(hyperplanes.ambient_rank, 3)
        self.assertEqual(len(hyperplanes.components[0].vertex_indices), 5)
        self.assertEqual(hyperplanes.components[0].rational_rank, 2)
        angle_profile = height_angle_profile(
            first.vertices, ((1.0, 0.0), (0.0, 1.0))
        )
        scaled_angle_profile = height_angle_profile(
            first.vertices, ((3.0, 0.0), (0.0, 3.0))
        )
        self.assertAlmostEqual(
            height_angle_profile_distance(angle_profile, scaled_angle_profile),
            0.0,
        )

        left_map = [-1] * len(first.vertices)
        right_map = [-1] * len(first.vertices)
        for source_index, target_index in enumerate(vertex_map):
            if source_index in (0, 1, 2):
                left_map[source_index] = target_index
            if source_index in (1, 2, 3):
                right_map[source_index] = target_index
        merges = merge_component_vertex_maps(
            first,
            transformed,
            left_map,
            right_map,
            held_out_source_indices=(3,),
            source_gram=((1.0, 0.0), (0.0, 1.0)),
            target_gram=((2.0, -1.0), (-1.0, 1.0)),
            maximum_height_angle_rms=1.0,
        )
        self.assertTrue(merges)
        self.assertEqual(merges[0].replay.replayed_source_rank, 2)
        self.assertEqual(merges[0].held_out_replayed_ray_count, 1)

        fingerprint_vectors = first.vertices
        fingerprint_heights = (1.0, 1.3, 2.1, 2.8)
        fingerprint_arithmetic = tuple(
            {"integral": index % 2 == 0, "total_bits": 10 + index}
            for index in range(4)
        )
        fingerprint = candidate_relation_fingerprint(
            fingerprint_vectors,
            fingerprint_heights,
            fingerprint_arithmetic,
            range(4),
            first,
            dimension=2,
            quantiles=4,
            projective_multiplicities=4,
        )
        transformed_vectors = tuple(
            canonical_unoriented((x + 2 * y, x + y))
            for x, y in fingerprint_vectors
        )
        transformed_fingerprint = candidate_relation_fingerprint(
            transformed_vectors,
            fingerprint_heights,
            fingerprint_arithmetic,
            range(4),
            build_relation_complex(transformed_vectors),
            dimension=2,
            quantiles=4,
            projective_multiplicities=4,
        )
        self.assertEqual(fingerprint, transformed_fingerprint)
        joint = joint_nearest_candidate_scores(
            ((fingerprint,), (transformed_fingerprint,))
        )
        self.assertEqual(float(joint[0][0].mean_distance), 0.0)
        self.assertEqual(joint[0][0].mutual_neighbour_count, 1)

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

    def test_primitive_hermite_signature_saturates_before_restriction(self) -> None:
        signatures = primitive_hermite_signatures(
            ((2.0, 0.0), (0.0, 3.0)),
            (((1, 0), (0, 1)), ((2, 0), (0, 1))),
        )
        self.assertEqual(tuple(item.saturation_index for item in signatures), (1, 2))
        self.assertLess(
            hermite_signature_distance(signatures[0].hermite, signatures[1].hermite),
            1e-12,
        )

    def test_intersection_consensus_uses_exact_codimension_one_support(self) -> None:
        ledger = exact_intersection_consensus(
            (
                ((1, 0, 0, 0), (0, 1, 0, 0)),
                ((1, 0, 0, 0), (0, 0, 1, 0)),
                ((1, 0, 0, 0), (0, 0, 0, 1)),
                ((0, 1, 0, 0), (0, 0, 1, 0)),
            ),
            (4.0, 3.0, 2.0, 1.0),
            pool_size=4,
        )
        self.assertEqual(ledger.exact_pair_count, 6)
        self.assertEqual(ledger.candidates[0].codimension_one_intersection_count, 3)
        self.assertEqual(ledger.selected.source_index, 0)

    def test_graph_walk_consensus_records_exact_walk_counts(self) -> None:
        ledger = exact_graph_walk_consensus(
            (
                ((1, 0, 0, 0), (0, 1, 0, 0)),
                ((1, 0, 0, 0), (0, 0, 1, 0)),
                ((1, 0, 0, 0), (0, 0, 0, 1)),
                ((0, 1, 0, 0), (0, 0, 1, 0)),
            ),
            (1.0, 0.99, 0.98, 0.97),
            pool_size=4,
            shape_gap_threshold=0.02,
        )
        self.assertEqual(ledger.exact_pair_count, 6)
        self.assertEqual(ledger.codimension_one_edge_count, 5)
        self.assertEqual(ledger.selector_mode, "exact_graph_walk_consensus")
        self.assertGreater(ledger.candidates[0].triangle_count, 0)
        self.assertEqual(ledger.selected.source_index, 0)

    def test_intrinsic_shell_relations_are_scale_and_basis_invariant(self) -> None:
        first = intrinsic_shell_signature(
            ((2.0, 0.0), (0.0, 3.0)), minimum_vectors=4, quantiles=4
        )
        second = intrinsic_shell_signature(
            ((35.0, 21.0), (21.0, 21.0)), minimum_vectors=4, quantiles=4
        )
        self.assertEqual(first.primitive_vector_count, second.primitive_vector_count)
        self.assertEqual(
            first.relation_complex.canonical_digest,
            second.relation_complex.canonical_digest,
        )
        self.assertLess(intrinsic_shell_profile_distance(first, second).total, 1e-12)

    def test_relation_metric_signature_is_scale_and_basis_invariant(self) -> None:
        vectors = ((1, 0), (0, 1), (1, 1), (1, -1))
        complex_ = build_relation_complex(vectors)
        first = relation_metric_signature(vectors, (2, 3, 5, 7), range(4), complex_, quantiles=5)
        transformed = tuple((x + 2 * y, x + y) for x, y in vectors)
        second = relation_metric_signature(
            transformed,
            (22, 33, 55, 77),
            (3, 2, 1, 0),
            build_relation_complex(transformed),
            quantiles=5,
        )
        self.assertLess(relation_metric_profile_distance(first, second), 1e-12)

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

    def test_repeated_cross_bound_ledger_counts_exact_spaces(self) -> None:
        coordinates = (
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (1, 1, 0, 0),
            (0, 0, 1, 0),
        )
        records = tuple(
            ShortVectorRecord(
                vector,
                "1",
                None,
                {"integral": index < 3, "total_bits": index + 1},
            )
            for index, vector in enumerate(coordinates)
        )
        complex_ = build_relation_complex(coordinates)
        from latent_lattice import GrowthProposal

        left = (
            GrowthProposal(3, ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)), (), 0, 0, 0.0, 0.0),
            GrowthProposal(3, ((1, 1, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)), (), 0, 0, 0.0, 0.0),
        )
        right = (
            GrowthProposal(3, ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1)), (), 0, 0, 0.0, 0.0),
            GrowthProposal(3, ((1, 1, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1)), (), 0, 0, 0.0, 0.0),
        )
        finite = FiniteQuotientBlock(5, 2, 4, 2, 1, ((1, 0, 0, 1),))
        ledger = repeated_cross_bound_intersection_ledger(
            records,
            complex_,
            left,
            right,
            target_dimension=2,
            finite_blocks=(finite,),
        )
        self.assertEqual(ledger.tested_pair_count, 4)
        self.assertEqual(ledger.modular_surviving_pair_count, 4)
        self.assertEqual(ledger.exact_candidate_count, 1)
        self.assertEqual(len(ledger.scored_candidates), 1)
        self.assertEqual(len(ledger.proposals), 1)
        proposal = ledger.proposals[0]
        self.assertEqual(proposal.occurrence_count, 4)
        self.assertEqual(proposal.best_pair, (0, 0))
        self.assertEqual(proposal.support, 3)
        self.assertGreater(proposal.induced_ternary_relation_count, 0)
        self.assertIsNotNone(proposal.finite_signature)
        self.assertEqual(
            proposal.primitive_basis_rows,
            ((1, 0, 0, 0), (0, 1, 0, 0)),
        )
        aggregate = aggregate_repeated_intersection_ledgers(
            (("h1-h2", ledger), ("h2-h3", ledger)), minimum_pair_support=2
        )
        self.assertEqual(len(aggregate), 1)
        self.assertEqual(aggregate[0].pair_support, 2)
        self.assertEqual(aggregate[0].total_occurrence_count, 8)


if __name__ == "__main__":
    unittest.main()
