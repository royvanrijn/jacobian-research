#!/usr/bin/env python3
"""Regression tests for the generic Case-2 infinity resolution."""

from case2_infinity_resolution import audit


result = audit(singular=None)
assert result.coordinate_degrees == (8, 12)
assert result.leading_orders == (4, 12)
assert result.characteristic_orders == (4, 13)
assert result.characteristic_numerator_term_count == 19
assert result.characteristic_numerator_parameter_degree == 5
assert result.normalization_translation_invariant
assert result.exceptional_rays == (
    (1, 1),
    (1, 2),
    (1, 3),
    (4, 13),
    (3, 10),
    (2, 7),
    (1, 4),
)
assert result.exceptional_self_intersections == (-2, -2, -5, -1, -2, -2, -2)
assert result.strict_transform_meets_ray == (4, 13)
assert result.blowup_count == 7
assert result.compatibility_equation_count == 7
assert result.endpoint_certificate_indices == (0, 1, 2, 3, 4, 5, 6)
assert result.endpoint_certificate_constraint_count == 7
assert result.endpoint_certificate_term_counts == (8,) * 7
assert result.endpoint_certificate_parameter_degrees == (3,) * 7
assert result.degree_twelve_open_nonempty is None
assert result.case2_excluded_by_j1_endpoint is None

print("PASS: the first Case-2 infinity characteristic is exact and nonzero")
print("PASS: the generic infinity branch is the primitive (4,13) cusp")
print("PASS: the seven-ray toric resolution graph is regular")
