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
)
from .relations import RelationComplex, build_relation_complex
from .finite import FiniteQuotientBlock, discriminant, finite_quotient_block
from .local import ComponentBlock, multiplicative_component_block
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
    beam_subspace_scan,
    exact_span_mask,
    exact_row_space_intersection,
    extend_core_proposals,
    independent_row_basis,
    independent_relation_growth_proposals,
    modular_row_space_key,
    modular_right_kernel_basis,
    primitive_span_basis,
    recombined_core_extension_search,
    relation_seeded_subspace_scan,
)

__all__ = [
    "EllipticCurve",
    "ExactEmbedding",
    "FiniteQuotientBlock",
    "GrowthProposal",
    "RecombinedSearchLedger",
    "ComponentBlock",
    "Point",
    "RelationComplex",
    "ShortVectorRecord",
    "SubspaceCandidate",
    "beam_subspace_scan",
    "build_relation_complex",
    "canonical_unoriented",
    "enumerate_short_vectors",
    "finite_quotient_block",
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
    "primitive_coordinates",
    "rational_nullspace",
    "recombined_core_extension_search",
    "relation_seeded_subspace_scan",
    "rational_rank",
    "recover_exact_embedding",
    "discriminant",
]
