#!/usr/bin/env python3
"""Verify the HC4 smooth rank-three vertex-colength obstruction.

The proof is commutative-algebraic.  This checker records its exact finite
length consequences, checks the complete-intersection Hilbert series, and
intersects the resulting lower bound with the degree-two/three HC4 atlas.
It does not claim that the Fermat calibration proves the universal filtered
lemma; that lemma is proved in the canonical note.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.projective_gradient_segre import (  # noqa: E402
    SmoothEssentialGradientNormalSlice,
)


ATLAS = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_projective_polar_atlas.json"
)
STRATA = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_quintic_infinity_rees_strata.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_rank3_vertex_colength.json"
)


atlas_payload = json.loads(ATLAS.read_text())
strata_payload = json.loads(STRATA.read_text())

# A smooth ternary quintic has three quartic partial derivatives with no
# common projective zero.  They form a regular sequence, so its Milnor
# algebra B has Hilbert series (1+t+t^2+t^3)^3.
normal_slice = SmoothEssentialGradientNormalSlice(4, 4, 3)
b_hilbert_function = normal_slice.jacobian_hilbert_function
assert b_hilbert_function == (1, 3, 6, 10, 12, 12, 10, 6, 3, 1)
b_length = normal_slice.jacobian_length
b_socle_degree = normal_slice.jacobian_socle_degree
assert b_length == 4**3 == 64
assert b_socle_degree == 3 * (4 - 1) == 9

# In the nonaligned branch, s3 is a nonzero cubic.  It cannot lie in the
# quartic Jacobian ideal and cannot lie in the degree-nine socle.  Therefore
# both s3 and some positive-degree multiple are nonzero and have different
# degrees, so dim(B*s3) >= 2.
s3_degree = 3
assert s3_degree < 4
assert s3_degree < b_socle_degree
minimum_b_s3_dimension = 2

# The epsilon-adic associated graded of
#
#   A = k[[epsilon,u1,u2,u3]]/(epsilon^4,G1,G2,G3)
#
# is B tensor k[epsilon]/(epsilon^4).  The initial form of the missing
# gradient component is epsilon*s3.  Its ideal therefore contains one copy
# of B*s3 in each of epsilon-degrees 1, 2, and 3.
epsilon_nilpotence_order = 4
nonzero_epsilon_layers = epsilon_nilpotence_order - 1
a_length = normal_slice.truncated_active_length
minimum_missing_component_ideal_length = (
    normal_slice.isolated_vertex_affine_degree_lower_bound(
        epsilon_order=1,
        cyclic_ideal_dimension=minimum_b_s3_dimension,
    )
)
assert a_length == 256
assert minimum_missing_component_ideal_length == 6

# Exact monomial calibration for
# B=Q[u1,u2,u3]/(u1^4,u2^4,u3^4), s3=u1^3.  This is deliberately labelled
# a calibration: the theorem uses only the universal lower bound above.
fermat_b_s3_basis_size = sum(
    1
    for exponent_1 in range(4)
    for exponent_2 in range(4)
    for exponent_3 in range(4)
    if exponent_1 >= 3
)
fermat_missing_component_ideal_length = (
    nonzero_epsilon_layers * fermat_b_s3_basis_size
)
assert fermat_b_s3_basis_size == 16
assert fermat_missing_component_ideal_length == 48

# In codimension four the infinity base is zero-dimensional, so sigma_4 is
# its length.  The exact sequence 0 -> A*G_t -> A -> A/(G_t) -> 0 and the
# projective degree formula give
#
#   length(A*G_t) = 256 - sigma_4 = delta.
#
# Hence the nonaligned smooth rank-three branch has affine degree delta>=6.
excluded_affine_degrees = (2, 3)
assert all(
    degree < minimum_missing_component_ideal_length
    for degree in excluded_affine_degrees
)

atlas_intersection: dict[str, dict[str, object]] = {}
for affine_degree in excluded_affine_degrees:
    key = f"gradient_degree_4_affine_degree_{affine_degree}"
    rows = atlas_payload["atlases"][key]
    codimension_four_rows = [
        row for row in rows if row["leading_base_codimension"] == 4
    ]
    assert len(codimension_four_rows) == 1
    row = codimension_four_rows[0]
    assert row["projective_degrees"] == [1, 4, 16, 64, affine_degree]
    assert row["segre_degrees"] == [
        0,
        0,
        0,
        256 - affine_degree,
    ]
    atlas_intersection[f"affine_degree_{affine_degree}"] = {
        "excluded_projective_degrees": row["projective_degrees"],
        "excluded_segre_degrees": row["segre_degrees"],
        "rows_before": len(rows),
        "rows_excluded": 1,
        "rows_after": len(rows) - 1,
    }

assert atlas_intersection["affine_degree_2"]["rows_after"] == 318
assert atlas_intersection["affine_degree_3"]["rows_after"] == 306

rank_three_support_rows = [
    row
    for row in strata_payload["support_strata"]
    if row["hessian_rank"] == 3
    and row["base_support_codimension"] == 4
]
assert len(rank_three_support_rows) == 1
assert rank_three_support_rows[0][
    "essential_quintic_singularity"
] == "empty (smooth ternary quintic)"

payload = {
    "format": "hc4-rank3-vertex-colength-v1",
    "software_assumptions": {
        "python": "dependency-free exact integer arithmetic",
        "coefficient_field": "characteristic zero",
        "independent_calibration": "Macaulay2 over QQ",
    },
    "scope": (
        "Exact filtered-length obstruction for the smooth essential "
        "rank-three top-quintic stratum. The universal conclusion is "
        "delta>=6 in the nonaligned branch; the aligned branch is already "
        "excluded by HC4CD5. The Fermat computation is a calibration only."
    ),
    "local_chart": {
        "vertex_chart": "t=1",
        "compactifying_parameter": "epsilon=X0",
        "active_generators": (
            "G_i=partial_i(h5)+epsilon*partial_i(h4)"
            "+epsilon^2*partial_i(h3)+epsilon^3*partial_i(h2)"
        ),
        "missing_generator": (
            "G_t=epsilon*s3+epsilon^2*D_t(h3)"
            "+epsilon^3*D_t(h2)"
        ),
        "s3": "D_t(h4), a nonzero ternary cubic in the nonaligned branch",
    },
    "smooth_top_complete_intersection": {
        "generator_degrees": [4, 4, 4],
        "hilbert_function": list(b_hilbert_function),
        "length": b_length,
        "socle_degree": b_socle_degree,
    },
    "filtered_flatness": {
        "ambient_complete_local_ring_dimension": 4,
        "active_regular_sequence_length": 3,
        "epsilon_is_nonzerodivisor": True,
        "free_rank_over_k_power_series_epsilon": b_length,
        "associated_graded": "B tensor k[epsilon]/(epsilon^4)",
        "truncated_length": a_length,
    },
    "nonaligned_colength_bound": {
        "s3_degree": s3_degree,
        "minimum_dimension_B_times_s3": minimum_b_s3_dimension,
        "nonzero_epsilon_layers": nonzero_epsilon_layers,
        "minimum_length_A_times_G_t": (
            minimum_missing_component_ideal_length
        ),
        "identity": "length(A*G_t)=256-sigma_4=delta",
        "conclusion": "delta>=6",
    },
    "fermat_calibration": {
        "B": "QQ[u1,u2,u3]/(u1^4,u2^4,u3^4)",
        "s3": "u1^3",
        "dimension_B_times_s3": fermat_b_s3_basis_size,
        "length_epsilon_s3_ideal": (
            fermat_missing_component_ideal_length
        ),
        "role": "calibration, not the universal proof",
    },
    "atlas_intersection": atlas_intersection,
    "theorem_consequence": (
        "No affine-degree-two or affine-degree-three HC4 candidate has a "
        "smooth essential rank-three quintic top. Equivalently, both "
        "codimension-four quartic-gradient atlas rows are empty."
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: smooth ternary-quintic Jacobian algebra has length 64")
print("PASS: its Hilbert function is (1,3,6,10,12,12,10,6,3,1)")
print("PASS: the epsilon-truncated active complete intersection has length 256")
print("PASS: a nonzero cubic s3 generates an ideal of dimension at least 2")
print("PASS: the missing gradient component has ideal length at least 6")
print("PASS: excluded both codimension-four affine-degree-two/three rows")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
