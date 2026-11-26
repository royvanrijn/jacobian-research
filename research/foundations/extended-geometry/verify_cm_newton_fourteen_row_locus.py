"""Check the optimal fourteen-row coherent window using SageMath."""
from sage.all import GF, QQ

from verify_cm_newton_331_l3_normal_charts import verify_l3_nonendpoint_chart
from verify_cm_newton_331_l3_vzero_endpoint import verify_vzero_endpoint
from verify_cm_newton_331_l3_wzero_endpoint import verify_wzero_endpoint
from verify_cm_newton_fourteen_row_integral_places import replay as verify_integral_places
from verify_cm_newton_full_width_projective_boundary import (
    verify_boundary,
    verify_first_weighted_component,
    verify_second_weighted_component,
    verify_weighted_component_intersection,
)
from verify_cm_newton_full_width_reconstruction_certificates import reduced_residual_rows
from verify_cm_newton_seventeen_row_locus import (
    exact_constant_span,
    specialize_rows,
    verify_angular_zero_slice,
)
from verify_cm_newton_slope_eliminant import certify as verify_full_slope_span
from verify_cm_newton_slope_residue_discs import (
    compressed_setup,
    exact_schur,
    local_pivots,
    series_mats,
)


def verify_prefix_slope_discs(rows):
    coefficients, _, _ = compressed_setup(rows)
    assert len(coefficients) == 428
    cases = (
        (QQ(17) / 2, [359, 2, 0, 0, 0, 2, 0, 1], 19, 680),
        (QQ(56) / 5, [361, 1, 0, 1, 0, 0, 1], 10, 682),
    )
    for alpha, expected, vanishing_order, affine_rank in cases:
        polynomials = specialize_rows(rows, "ell=alpha*v", alpha)
        exact_constant_span(polynomials, (1, 1, 1, 1), (), 9, (847, 715, affine_rank))
        print(f"Slope {alpha}: exact affine unit identity of rank {affine_rank}", flush=True)

        order = len(expected) - 1
        modular = series_mats(coefficients, alpha, order, GF(331))
        levels, row_indices, column_indices = local_pivots(modular, 364)
        assert levels == expected
        assert sum(i * count for i, count in enumerate(levels)) == vanishing_order
        complement = exact_schur(
            coefficients, row_indices, column_indices, levels[0], alpha, order,
        )
        exact_levels, _, _ = local_pivots(complement, 364 - levels[0])
        assert exact_levels == [0, *expected[1:]]
        print(f"Slope {alpha}: exact and residue orders both {vanishing_order}", flush=True)


def replay():
    rows = reduced_residual_rows(5, 13, even_linear=True, primitive=False)
    verify_integral_places((331,))
    verify_boundary(rows)
    verify_first_weighted_component(rows)
    verify_second_weighted_component(rows)
    verify_weighted_component_intersection(rows)
    verify_angular_zero_slice(rows)
    verify_l3_nonendpoint_chart(rows)
    verify_vzero_endpoint(rows)
    verify_wzero_endpoint(rows)
    # The shorter row span inherits its generic rank upper bound and
    # the constant direction from the exact full slope calculation.
    verify_full_slope_span(reduced_residual_rows(5, 16, even_linear=True, primitive=False))
    verify_prefix_slope_discs(rows)
    print("The optimal initial window has fourteen rows.", flush=True)


if __name__ == "__main__":
    replay()
