"""Basis-independent finite fingerprints for latent Mordell--Weil lattices.

The package deliberately separates three evidence layers:

``exact``
    integer-coordinate, elliptic group-law, reduction-code, and relation data;
``numerical``
    high-precision canonical heights and derived scores;
``search``
    bounded enumeration and matching results.

Nothing in this package promotes a bounded match to a generic-family theorem.
"""

from .elliptic import EllipticCurve, Point, point_complexity
from .integer import (
    canonical_unoriented,
    primitive_coordinates,
    rational_nullspace,
    rational_rank,
    row_basis_coordinates,
)
from .relations import RelationComplex, build_relation_complex
from .finite import FiniteQuotientBlock, discriminant, finite_quotient_block
from .local import ComponentBlock, multiplicative_component_block
from .codes import (
    CandidateFiniteSignature,
    ComponentCodeSignature,
    PrimeCodeSignature,
    candidate_finite_signature,
    candidate_finite_signature_from_record,
    finite_joint_class_key,
    finite_rarity_scores,
    finite_rarity_weights,
    finite_signature_distance,
)
from .height import (
    CloudHeightSignature,
    HermiteSignature,
    IntrinsicShellSignature,
    IntrinsicShellDistance,
    RelationMetricSignature,
    ThetaSignature,
    cloud_height_profile_distance,
    cloud_height_signature,
    cloud_height_signature_from_record,
    hermite_signature,
    hermite_signature_distance,
    intrinsic_shell_signature,
    intrinsic_shell_profile_distance,
    relation_metric_profile_distance,
    relation_metric_signature,
    restricted_height_gram,
    theta_profile_distance,
    theta_signature,
)
from .pari import (
    ExactEmbedding,
    ShortVectorRecord,
    enumerate_short_vectors,
    height_gram,
    primitive_column_closure,
    recover_exact_embedding,
)
from .subspace import (
    SubspaceCandidate,
    GrowthProposal,
    RecombinedSearchLedger,
    RepeatedIntersectionLedger,
    RepeatedIntersectionProposal,
    RepeatedIntersectionScore,
    MultiBoundIntersectionScore,
    aggregate_repeated_intersection_ledgers,
    beam_subspace_scan,
    cross_bound_intersection_proposals,
    exact_span_mask,
    exact_row_space_intersection,
    extend_core_proposals,
    independent_row_basis,
    independent_relation_growth_proposals,
    modular_row_space_key,
    modular_right_kernel_basis,
    primitive_span_basis,
    repeated_cross_bound_intersection_ledger,
    recombined_core_extension_search,
    relation_seeded_subspace_scan,
)

__all__ = [
    "EllipticCurve",
    "ExactEmbedding",
    "FiniteQuotientBlock",
    "GrowthProposal",
    "RecombinedSearchLedger",
    "RepeatedIntersectionLedger",
    "RepeatedIntersectionProposal",
    "RepeatedIntersectionScore",
    "MultiBoundIntersectionScore",
    "ComponentBlock",
    "CandidateFiniteSignature",
    "ComponentCodeSignature",
    "CloudHeightSignature",
    "HermiteSignature",
    "IntrinsicShellSignature",
    "IntrinsicShellDistance",
    "RelationMetricSignature",
    "Point",
    "PrimeCodeSignature",
    "RelationComplex",
    "ShortVectorRecord",
    "SubspaceCandidate",
    "ThetaSignature",
    "beam_subspace_scan",
    "aggregate_repeated_intersection_ledgers",
    "build_relation_complex",
    "candidate_finite_signature",
    "candidate_finite_signature_from_record",
    "cloud_height_profile_distance",
    "cloud_height_signature",
    "cloud_height_signature_from_record",
    "cross_bound_intersection_proposals",
    "hermite_signature",
    "hermite_signature_distance",
    "intrinsic_shell_signature",
    "intrinsic_shell_profile_distance",
    "relation_metric_profile_distance",
    "relation_metric_signature",
    "canonical_unoriented",
    "enumerate_short_vectors",
    "finite_quotient_block",
    "finite_joint_class_key",
    "finite_rarity_scores",
    "finite_rarity_weights",
    "finite_signature_distance",
    "exact_span_mask",
    "exact_row_space_intersection",
    "extend_core_proposals",
    "height_gram",
    "independent_row_basis",
    "independent_relation_growth_proposals",
    "modular_row_space_key",
    "modular_right_kernel_basis",
    "multiplicative_component_block",
    "point_complexity",
    "primitive_column_closure",
    "primitive_span_basis",
    "repeated_cross_bound_intersection_ledger",
    "primitive_coordinates",
    "rational_nullspace",
    "recombined_core_extension_search",
    "relation_seeded_subspace_scan",
    "rational_rank",
    "row_basis_coordinates",
    "recover_exact_embedding",
    "restricted_height_gram",
    "theta_profile_distance",
    "theta_signature",
    "discriminant",
]
