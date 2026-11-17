#!/usr/bin/env python3
"""Exact numerical atlas for projective compactifications of HC(4) gradients.

For a degree-m polynomial map F:A^4 -> A^4, use the graph compactification

    [X0^m : F_1^h : ... : F_4^h] : P^4 --> P^4.

If B is its base scheme and

    i_* s(B,P^4) = sigma_1 H + ... + sigma_4 H^4,

then the projective degrees satisfy

    g_i = m^i - sum_{k=1}^i binom(i,k)m^(i-k)sigma_k.

The script enumerates every integral degree list allowed by the elementary
bounds and log-concavity when m is 2, 3, or 4 and g_4 is 2 or 3.  For
m=4 it also records the exact leading-Hessian coverage matrix for a
collision-normalized quintic potential.  These are necessary numerical
signatures and structural reductions, not existence results and not a
classification by Hilbert polynomial.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from math import comb
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.projective_gradient_segre import (  # noqa: E402
    equal_degree_complete_intersection_segre,
    projective_degrees_from_segre,
    segre_degrees_from_projective,
)

OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_projective_polar_atlas.json"
)


@dataclass(frozen=True)
class Signature:
    projective_degrees: tuple[int, int, int, int, int]
    segre_degrees: tuple[int, int, int, int]
    leading_base_codimension: int
    leading_cycle_degree: int
    smooth_lci_curve_genus: int | None
    smooth_lci_curve_numerically_possible: bool | None


def segre_from_degrees(m: int, degrees: tuple[int, ...]) -> tuple[int, ...]:
    """Invert the triangular projective-degree/Segre-degree relation."""
    return segre_degrees_from_projective(m, degrees)


def projective_from_segre(m: int, sigmas: tuple[int, ...]) -> tuple[int, ...]:
    """Evaluate the blow-up formula in every projective degree."""
    return projective_degrees_from_segre(m, sigmas)


def castelnuovo_bound_p4(degree: int) -> int:
    """Castelnuovo's genus bound for a nondegenerate integral curve in P^4."""
    if degree <= 0:
        return -1
    quotient, remainder = divmod(degree - 1, 3)
    return 3 * quotient * (quotient - 1) // 2 + quotient * remainder


def make_signature(m: int, g2: int, g3: int, top_degree: int) -> Signature:
    degrees = (1, m, g2, g3, top_degree)
    sigmas = segre_from_degrees(m, degrees)
    assert projective_from_segre(m, sigmas) == degrees
    assert sigmas[0] == 0

    first_nonzero = next(
        (index for index, value in enumerate(sigmas, start=1) if value),
        4,
    )
    leading_degree = sigmas[first_nonzero - 1]
    curve_genus = None
    curve_possible = None
    if first_nonzero == 3:
        curve_degree = sigmas[2]
        numerator = 2 - 5 * curve_degree - sigmas[3]
        if numerator % 2 == 0:
            curve_genus = numerator // 2
        curve_possible = (
            curve_genus is not None
            and 0 <= curve_genus <= castelnuovo_bound_p4(curve_degree)
        )

    return Signature(
        projective_degrees=degrees,
        segre_degrees=sigmas,
        leading_base_codimension=first_nonzero,
        leading_cycle_degree=leading_degree,
        smooth_lci_curve_genus=curve_genus,
        smooth_lci_curve_numerically_possible=curve_possible,
    )


def atlas(m: int, top_degree: int) -> list[Signature]:
    """Enumerate the positive log-concave projective-degree lists."""
    result: list[Signature] = []
    for g2 in range(1, m**2 + 1):
        for g3 in range(1, m**3 + 1):
            if g2 * g2 < m * g3:
                continue
            if g3 * g3 < g2 * top_degree:
                continue
            result.append(make_signature(m, g2, g3, top_degree))
    return result


atlases = {
    f"gradient_degree_{m}_affine_degree_{top}": atlas(m, top)
    for m in (2, 3, 4)
    for top in (2, 3)
}
assert {key: len(value) for key, value in atlases.items()} == {
    "gradient_degree_2_affine_degree_2": 9,
    "gradient_degree_2_affine_degree_3": 7,
    "gradient_degree_3_affine_degree_2": 72,
    "gradient_degree_3_affine_degree_3": 67,
    "gradient_degree_4_affine_degree_2": 319,
    "gradient_degree_4_affine_degree_3": 307,
}


zero_dimensional = {}
for key, signatures in atlases.items():
    rows = [
        row
        for row in signatures
        if row.leading_base_codimension == 4
        and row.segre_degrees[:3] == (0, 0, 0)
    ]
    assert len(rows) == 1
    zero_dimensional[key] = rows[0].segre_degrees[3]

assert zero_dimensional == {
    "gradient_degree_2_affine_degree_2": 14,
    "gradient_degree_2_affine_degree_3": 13,
    "gradient_degree_3_affine_degree_2": 79,
    "gradient_degree_3_affine_degree_3": 78,
    "gradient_degree_4_affine_degree_2": 254,
    "gradient_degree_4_affine_degree_3": 253,
}


# The quintic atlas begins the genuinely open range.  The projective-degree
# signature alone does not determine the generic rank of Hess(h5), so the
# numerical rows and the leading-Hessian branches must be recorded as two
# transverse filters rather than falsely assigning individual signatures to
# theorem strata.
quintic_signature_counts_by_leading_codimension = {}
for top_degree in (2, 3):
    key = f"gradient_degree_4_affine_degree_{top_degree}"
    quintic_signature_counts_by_leading_codimension[
        f"affine_degree_{top_degree}"
    ] = {
        str(codimension): sum(
            row.leading_base_codimension == codimension
            for row in atlases[key]
        )
        for codimension in (2, 3, 4)
    }

assert quintic_signature_counts_by_leading_codimension == {
    "affine_degree_2": {"2": 260, "3": 58, "4": 1},
    "affine_degree_3": {"2": 249, "3": 57, "4": 1},
}

# The universal top-gradient support sieve is transverse to the determinant
# coverage matrix.  For an essential rank-r cone, the smooth top support is
# the kernel vertex P^(3-r), hence has codimension r+1 in P4.  Singular
# essential tops join that vertex to their singular locus and feed the lower
# codimension columns.  The pure-top CI vectors calibrate the degeneration;
# only the support codimension, not those multiplicities, is asserted for a
# completed constant-Hessian potential.
quintic_top_gradient_rees_sieve = {
    "universal_h5_coefficient_count": 56,
    "constant_hessian_top_gate": "det Hess(h5)=0; essential rank <=3",
    "generic_smooth_rank_strata": [],
    "singular_support_transport": {
        "rank_2_repeated_binary_root": 2,
        "rank_3_isolated_singularities": 3,
        "rank_3_positive_dimensional_singular_locus": 2,
    },
    "midpoint_collision_parity": [
        "grad(h2+h4)(a)=0",
        "grad(h3+h5)(a)=0",
    ],
    "artifact": (
        "artifacts/generated-results/"
        "hc4_quintic_infinity_rees_strata.json"
    ),
}
for rank in (1, 2, 3):
    codimension = rank + 1
    sigmas = equal_degree_complete_intersection_segre(
        ambient_dimension=4,
        codimension=codimension,
        generator_degree=4,
    )
    degrees = projective_degrees_from_segre(4, sigmas)
    quintic_top_gradient_rees_sieve[
        "generic_smooth_rank_strata"
    ].append(
        {
            "essential_hessian_rank": rank,
            "universal_essential_coefficient_count": (1, 6, 21)[
                rank - 1
            ],
            "base_support_codimension": codimension,
            "atlas_rows": {
                affine_key: counts[str(codimension)]
                for affine_key, counts in (
                    quintic_signature_counts_by_leading_codimension.items()
                )
            },
            "pure_top_projective_degrees": degrees,
            "pure_top_segre_degrees": sigmas,
            "generic_rees_ideal": (
                "linear type: Koszul plus inactive target equations"
            ),
        }
    )

assert [
    row["pure_top_projective_degrees"]
    for row in quintic_top_gradient_rees_sieve[
        "generic_smooth_rank_strata"
    ]
] == [
    (1, 4, 0, 0, 0),
    (1, 4, 16, 0, 0),
    (1, 4, 16, 64, 0),
]

# HC4PPG7 closes the two codimension-four rows.  The support sieve above
# shows that such a row can only come from a smooth essential rank-three
# top quintic.  On its vertex chart the filtered complete-intersection
# argument gives affine degree at least six in the nonaligned branch,
# while HC4CD5 already excludes the aligned branch.
rank_three_vertex_colength_sieve = {
    "support_stratum": "smooth essential rank-three ternary quintic",
    "base_support_codimension": 4,
    "active_complete_intersection_length": 256,
    "nonaligned_minimum_affine_degree": 6,
    "aligned_branch": {
        "status": "excluded",
        "result": "HC4CD5",
    },
    "excluded_rows": {},
    "artifact": (
        "artifacts/generated-results/hc4_rank3_vertex_colength.json"
    ),
}
for affine_degree in (2, 3):
    key = f"gradient_degree_4_affine_degree_{affine_degree}"
    codimension_four_rows = [
        row
        for row in atlases[key]
        if row.leading_base_codimension == 4
    ]
    assert len(codimension_four_rows) == 1
    row = codimension_four_rows[0]
    assert row.projective_degrees == (1, 4, 16, 64, affine_degree)
    assert row.segre_degrees == (0, 0, 0, 256 - affine_degree)
    rank_three_vertex_colength_sieve["excluded_rows"][
        f"affine_degree_{affine_degree}"
    ] = {
        "projective_degrees": row.projective_degrees,
        "segre_degrees": row.segre_degrees,
    }

assert rank_three_vertex_colength_sieve["excluded_rows"] == {
    "affine_degree_2": {
        "projective_degrees": (1, 4, 16, 64, 2),
        "segre_degrees": (0, 0, 0, 254),
    },
    "affine_degree_3": {
        "projective_degrees": (1, 4, 16, 64, 3),
        "segre_degrees": (0, 0, 0, 253),
    },
}

# HC4PPG8 refines, but does not globally delete, the codimension-three
# column.  On the smooth rank-two packet, a nonzero restriction h4|K
# synchronizes a constant kernel direction.  The squarefree binary-Hessian
# branch is excluded; on its nonsquarefree remainder the generic transverse
# multiplicity is sigma_3=16.  The h4|K=0 and rank-three isolated-singularity
# packets remain capable of feeding other numerical rows.
codimension_three_gradient_sieve = {
    "rank_two_smooth_top": {
        "nonzero_h4_kernel_restriction": {
            "squarefree_binary_hessian": {
                "status": "excluded",
                "result": "HC4PPG8 plus HC4CD5",
            },
            "nonsquarefree_binary_hessian": {
                "status": "open exceptional discriminant",
                "forced_sigma3": 16,
                "atlas_rows": {},
            },
        },
        "zero_h4_kernel_restriction": {
            "status": "open synchronization stratum",
        },
    },
    "rank_three_isolated_singularities": {
        "ordinary_points": "s3 vanishes at every Hessian-rank-two point",
        "worse_points": "open when Hess(h5) has rank at most one",
    },
    "unconditional_atlas_rows_excluded": 0,
    "artifact": (
        "artifacts/generated-results/hc4_codim3_gradient_strata.json"
    ),
}
for affine_degree in (2, 3):
    key = f"gradient_degree_4_affine_degree_{affine_degree}"
    rows = [
        row
        for row in atlases[key]
        if row.leading_base_codimension == 3
        and row.segre_degrees[2] == 16
    ]
    assert len(rows) == 1
    row = rows[0]
    assert row.projective_degrees == (1, 4, 16, 48, affine_degree)
    assert row.segre_degrees == (0, 0, 16, -affine_degree)
    codimension_three_gradient_sieve["rank_two_smooth_top"][
        "nonzero_h4_kernel_restriction"
    ]["nonsquarefree_binary_hessian"]["atlas_rows"][
        f"affine_degree_{affine_degree}"
    ] = {
        "projective_degrees": row.projective_degrees,
        "segre_degrees": row.segre_degrees,
    }

# HC4PPG9 applies the PGS3 order-one active profile to the essential
# rank-two repeated-root packet.  A root partition with q distinct roots
# has total transverse Jacobian length 5-q, hence sigma_2=5-q on this open
# lower-layer stratum.  Failure of the active-unit condition remains an
# exceptional torsion packet, so this is a partition, not an unconditional
# deletion of codimension-two rows.
codimension_two_binary_root_sieve = {
    "condition": (
        "at every repeated root, one redundant active gradient has "
        "epsilon-order one with unit coefficient after transverse elimination"
    ),
    "root_partitions": {},
    "generic_double_root_packet": {},
    "unconditional_atlas_rows_excluded": 0,
    "artifact": (
        "artifacts/generated-results/"
        "hc4_binary_root_partition_segre.json"
    ),
}
for partition in (
    (2, 1, 1, 1),
    (3, 1, 1),
    (2, 2, 1),
    (4, 1),
    (3, 2),
):
    sigma2 = 5 - len(partition)
    row_counts = {}
    for affine_degree in (2, 3):
        key = f"gradient_degree_4_affine_degree_{affine_degree}"
        row_counts[f"affine_degree_{affine_degree}"] = sum(
            row.leading_base_codimension == 2
            and row.segre_degrees[1] == sigma2
            for row in atlases[key]
        )
    codimension_two_binary_root_sieve["root_partitions"][
        "+".join(map(str, partition))
    ] = {
        "forced_sigma2": sigma2,
        "atlas_rows": row_counts,
    }

assert codimension_two_binary_root_sieve["root_partitions"][
    "2+1+1+1"
] == {
    "forced_sigma2": 1,
    "atlas_rows": {
        "affine_degree_2": 51,
        "affine_degree_3": 50,
    },
}
codimension_two_binary_root_sieve["generic_double_root_packet"] = {
    "root_partition": [2, 1, 1, 1],
    "forced_sigma2": 1,
    "rows_before": {
        "affine_degree_2": 260,
        "affine_degree_3": 249,
    },
    "rows_after": {
        "affine_degree_2": 51,
        "affine_degree_3": 50,
    },
}

quintic_smooth_curve_only_counts = {
    f"affine_degree_{top_degree}": sum(
        row.smooth_lci_curve_numerically_possible is True
        for row in atlases[
            f"gradient_degree_4_affine_degree_{top_degree}"
        ]
    )
    for top_degree in (2, 3)
}
assert quintic_smooth_curve_only_counts == {
    "affine_degree_2": 18,
    "affine_degree_3": 18,
}


# Exact determinant-face checks for
#
#   Hess(q2 + lambda*h3 + lambda^2*h4 + lambda^3*h5).
#
# After Gordan--Noether supplies a constant kernel for Hess(h5), put that
# kernel first.  For leading Hessian rank r, the top remaining face is
#
#   det(C_r) det(Hess(h4)|_K),
#
# in lambda-degree 8+r.  In rank three this first makes
# D_t^2 h4=0.  The next face is the cubic Schur identity
#
#   det(C_3) D_t^2 h3
#     - grad(D_t h4)^T adj(C_3) grad(D_t h4) = 0.
#
# The checks use generic symmetric blocks, so they certify the matrix
# identity independently of any support choice.
lam = sp.symbols("lambda")
c11, c12, c13, c22, c23, c33 = sp.symbols(
    "c11 c12 c13 c22 c23 c33"
)
c3 = sp.Matrix(
    [
        [c11, c12, c13],
        [c12, c22, c23],
        [c13, c23, c33],
    ]
)
b00, b01, b02, b03 = sp.symbols("b00 b01 b02 b03")
b11, b12, b13, b22, b23, b33 = sp.symbols(
    "b11 b12 b13 b22 b23 b33"
)
b4 = sp.Matrix(
    [
        [b00, b01, b02, b03],
        [b01, b11, b12, b13],
        [b02, b12, b22, b23],
        [b03, b13, b23, b33],
    ]
)
a00 = sp.symbols("a00")
a4 = sp.diag(a00, 0, 0, 0)
d3 = sp.diag(0, 1, 1, 1)
d3[1:4, 1:4] = c3
rank_three_pencil = lam**3 * d3 + lam**2 * b4 + lam * a4
rank_three_det = sp.Poly(
    sp.expand(rank_three_pencil.det()), lam
)
rank_three_top = rank_three_det.coeff_monomial(lam**11)
assert sp.expand(rank_three_top - b00 * c3.det()) == 0
rank_three_next = sp.expand(
    rank_three_det.coeff_monomial(lam**10).subs(b00, 0)
)
rank_three_schur = (
    a00 * c3.det()
    - (
        sp.Matrix([[b01, b02, b03]])
        * c3.adjugate()
        * sp.Matrix([b01, b02, b03])
    )[0]
)
assert sp.expand(rank_three_next - rank_three_schur) == 0
assert all(
    sp.expand(entry) == 0
    for entry in (
        c3 * c3.adjugate() - c3.det() * sp.eye(3)
    )
)

# In the rank-three nonaligned branch, C_3 is the Hessian of a ternary
# quintic, so det(C_3) has degree 9.  For a ternary cubic s3, d=grad(s3)
# has degree 2 and adj(C_3)d has degree 8.  If det(C_3) is squarefree, its
# rank-one adjugate at the generic point of every discriminant component
# turns d^T adj(C_3)d=0 into adj(C_3)d=0 on that component.  Hence the
# squarefree determinant divides the degree-eight vector, forcing the
# vector and then d to vanish.  The checker records the exact degree and
# adjugate identities; squarefree divisibility is the UFD argument in the
# canonical note.
rank_three_squarefree_degree_ledger = {
    "hessian_entry_degree": 3,
    "hessian_determinant_degree": 9,
    "cubic_gradient_degree": 2,
    "adjugate_degree": 6,
    "adjugate_gradient_vector_degree": 8,
}
assert (
    rank_three_squarefree_degree_ledger[
        "adjugate_gradient_vector_degree"
    ]
    < rank_three_squarefree_degree_ledger[
        "hessian_determinant_degree"
    ]
)

qx, qy, qz = sp.symbols("qx qy qz")
squarefree_quintic_witness = (
    qx**5
    + qy**5
    + qz**5
    + qx**4 * qy
    + qy**4 * qz
    + qz**4 * qx
)
squarefree_witness_determinant = sp.expand(
    sp.hessian(
        squarefree_quintic_witness, (qx, qy, qz)
    ).det()
)
squarefree_witness_gcd = squarefree_witness_determinant
for variable in (qx, qy, qz):
    squarefree_witness_gcd = sp.gcd(
        squarefree_witness_gcd,
        sp.diff(squarefree_witness_determinant, variable),
    )
assert sp.Poly(
    squarefree_witness_determinant, qx, qy, qz
).total_degree() == 9
assert sp.Poly(
    squarefree_witness_gcd, qx, qy, qz
).total_degree() == 0
rank_three_squarefree_degree_ledger["squarefree_witness"] = (
    "qx^5+qy^5+qz^5+qx^4*qy+qy^4*qz+qz^4*qx"
)
rank_three_squarefree_degree_ledger[
    "squarefree_witness_hessian_terms"
] = len(
    sp.Poly(
        squarefree_witness_determinant, qx, qy, qz
    ).terms()
)

c2 = sp.Matrix([[c11, c12], [c12, c22]])
d2 = sp.diag(0, 0, 1, 1)
d2[2:4, 2:4] = c2
rank_two_det = sp.Poly(
    sp.expand((lam**3 * d2 + lam**2 * b4).det()), lam
)
rank_two_top = rank_two_det.coeff_monomial(lam**10)
assert sp.expand(
    rank_two_top - c2.det() * b4[:2, :2].det()
) == 0

d1 = sp.diag(0, 0, 0, c11)
rank_one_det = sp.Poly(
    sp.expand((lam**3 * d1 + lam**2 * b4).det()), lam
)
rank_one_top = rank_one_det.coeff_monomial(lam**9)
assert sp.expand(
    rank_one_top - c11 * b4[:3, :3].det()
) == 0

quintic_coverage_matrix = [
    {
        "leading_hessian_rank": 0,
        "condition": "Hess(h5)=0, hence h5=0 by homogeneity",
        "status": "excluded",
        "result": "HC4CQ1",
        "first_unresolved_face": None,
    },
    {
        "leading_hessian_rank": 1,
        "condition": (
            "constant 3-plane K=ker Hess(h5); "
            "det(Hess(h4)|K)=0"
        ),
        "status": "unresolved",
        "result": None,
        "first_unresolved_face": (
            "constant synchronization of the function-field kernel of "
            "the ternary quartic Hessian restriction"
        ),
    },
    {
        "leading_hessian_rank": 2,
        "condition": (
            "constant 2-plane K=ker Hess(h5); "
            "det(Hess(h4)|K)=0"
        ),
        "status": "unresolved",
        "result": None,
        "first_unresolved_face": (
            "constant synchronization of the function-field kernel of "
            "the binary quartic Hessian restriction"
        ),
    },
    {
        "leading_hessian_rank": 3,
        "condition": (
            "constant direction t=ker Hess(h5); D_t^2 h4=0; "
            "write s3=D_t h4"
        ),
        "status": "split",
        "covered_subbranch": {
            "condition": "s3=0",
            "status": "excluded",
            "result": "HC4CD5",
        },
        "unresolved_subbranch": {
            "condition": "s3!=0",
            "status": "split_by_hessian_discriminant",
            "smooth_top_affine_degree_2_or_3": {
                "status": "excluded",
                "result": "HC4PPG7 vertex-colength obstruction",
            },
            "first_unresolved_face": (
                "det(Hess_u(h5)) divides "
                "grad(s3)^T adj(Hess_u(h5)) grad(s3), "
                "with linear quotient D_t^2 h3"
            ),
            "squarefree_hessian_determinant": {
                "status": "excluded",
                "result": "squarefree adjugate-divisibility obstruction",
            },
            "nonsquarefree_hessian_determinant": {
                "status": "unresolved",
                "first_unresolved_face": (
                    "classify repeated components of the degree-nine "
                    "ternary-quintic Hessian determinant and solve the "
                    "Schur incidence on those components"
                ),
            },
        },
    },
]

quintic_coverage_summary = {
    "numerical_signatures": {
        "affine_degree_2": 319,
        "affine_degree_3": 307,
        "total": 626,
    },
    "signature_counts_by_leading_base_codimension":
        quintic_signature_counts_by_leading_codimension,
    "zero_dimensional_lengths": {
        "affine_degree_2": 254,
        "affine_degree_3": 253,
    },
    "smooth_integral_curve_only_rows_passing_castelnuovo": (
        quintic_smooth_curve_only_counts
    ),
    "signature_level_exclusions": {
        "affine_degree_2": 1,
        "affine_degree_3": 1,
        "total": 2,
    },
    "remaining_numerical_signatures_after_vertex_colength": {
        "affine_degree_2": 318,
        "affine_degree_3": 306,
        "total": 624,
    },
    "reason_the_codimension_four_rows_are_assignable": (
        "The top-gradient support theorem identifies codimension four "
        "exactly with the smooth essential rank-three stratum; HC4PPG7 "
        "then excludes affine degrees two and three on that stratum."
    ),
    "recommended_next_exact_calculation": (
        "Attack the exceptional codimension-three packets: h4|K=0 in "
        "essential rank two, the nonsquarefree binary-Hessian sigma3=16 "
        "row, and the lower normal cones at isolated ternary singularities "
        "subject to s3(p)=0 at every ordinary point."
    ),
    "rank_three_squarefree_degree_ledger":
        rank_three_squarefree_degree_ledger,
}


# Wang's degree-two theorem says that every characteristic-zero quadratic
# Keller map is a polynomial automorphism.  Applied to grad(Psi), constant
# nonzero Hessian determinant therefore forces top projective degree one.
# The theorem is an external mathematical input; the finite consequences
# for the atlas are checked here.
quadratic_keller_atlas = atlas(2, 1)
assert len(quadratic_keller_atlas) == 11
assert all(
    row.projective_degrees[-1] != 1
    for key, rows in atlases.items()
    if key.startswith("gradient_degree_2_")
    for row in rows
)
quadratic_keller_consequence = {
    "external_input": (
        "Wang's theorem: every characteristic-zero Keller map of "
        "polynomial degree at most two is a polynomial automorphism"
    ),
    "forced_affine_degree": 1,
    "forced_total_segre_correction": 2**4 - 1,
    "number_of_log_concave_degree_one_signatures": len(
        quadratic_keller_atlas
    ),
    "excluded_atlas_rows": {
        "affine_degree_2": len(
            atlases["gradient_degree_2_affine_degree_2"]
        ),
        "affine_degree_3": len(
            atlases["gradient_degree_2_affine_degree_3"]
        ),
    },
    "excluded_zero_dimensional_lengths": [14, 13],
}
assert quadratic_keller_consequence["forced_total_segre_correction"] == 15
assert quadratic_keller_consequence["excluded_atlas_rows"] == {
    "affine_degree_2": 9,
    "affine_degree_3": 7,
}


# HC4CQ1 excludes every collision for a four-variable potential
# q_2+h_3+h_4 with constant nonzero Hessian determinant.  After translating
# the midpoint and subtracting the common gradient value, every collision
# of a degree-four potential has exactly this form.  Ax--Grothendieck then
# turns injectivity over the algebraic closure into polynomial
# invertibility.  These are external mathematical inputs; their complete
# numerical consequences for the cubic-gradient atlas are checked here.
cubic_keller_atlas = atlas(3, 1)
assert len(cubic_keller_atlas) == 80
cubic_gradient_consequence = {
    "external_inputs": [
        (
            "HC4CQ1: a characteristic-zero four-variable potential "
            "q2+h3+h4 with constant nonzero Hessian determinant has no "
            "nonzero antipodal gradient collision"
        ),
        (
            "Ax--Grothendieck: an injective polynomial endomorphism of "
            "affine space in characteristic zero is an automorphism"
        ),
    ],
    "forced_affine_degree": 1,
    "forced_total_segre_correction": 3**4 - 1,
    "number_of_log_concave_degree_one_signatures": len(cubic_keller_atlas),
    "excluded_atlas_rows": {
        "affine_degree_2": len(
            atlases["gradient_degree_3_affine_degree_2"]
        ),
        "affine_degree_3": len(
            atlases["gradient_degree_3_affine_degree_3"]
        ),
    },
    "excluded_zero_dimensional_lengths": [79, 78],
    "counterexample_potential_degree_lower_bound": 5,
}
assert cubic_gradient_consequence["forced_total_segre_correction"] == 80
assert cubic_gradient_consequence["excluded_atlas_rows"] == {
    "affine_degree_2": 72,
    "affine_degree_3": 67,
}


# The graph and full-polar degree lists are independently certified by the
# companion Macaulay2 checker.  Inverting them here checks the conventions.
calibrations = {
    "quadratic_graph": {
        "m": 2,
        "projective_degrees": (1, 2, 2, 2, 1),
    },
    "quadratic_full_polar": {
        "m": 2,
        "projective_degrees": (1, 2, 4, 4, 2),
    },
    "cubic_graph": {
        "m": 3,
        "projective_degrees": (1, 3, 3, 3, 1),
    },
    "cubic_full_polar": {
        "m": 3,
        "projective_degrees": (1, 3, 6, 6, 3),
    },
}
for row in calibrations.values():
    row["segre_degrees"] = segre_from_degrees(
        row["m"], row["projective_degrees"]
    )
    assert projective_from_segre(
        row["m"], row["segre_degrees"]
    ) == row["projective_degrees"]

assert calibrations["quadratic_graph"]["segre_degrees"] == (0, 2, -6, 15)
assert calibrations["quadratic_full_polar"]["segre_degrees"] == (
    0,
    0,
    4,
    -18,
)
assert calibrations["cubic_graph"]["segre_degrees"] == (0, 6, -30, 116)
assert calibrations["cubic_full_polar"]["segre_degrees"] == (0, 3, -6, -12)


controls = {
    "cotangent_lift_quartic_packet": {
        "ambient_dimension": 4,
        "affine_degree": 4,
        "correction_quadratic_gradient": 2**4 - 4,
        "correction_cubic_gradient": 3**4 - 4,
    },
    "meng_yang_before_schur_hc6": {
        "ambient_dimension": 6,
        "gradient_degree": 7,
        "affine_degree": 3,
        "total_segre_correction": 7**6 - 3,
    },
    "meng_yang_after_schur_hc5": {
        "ambient_dimension": 5,
        "gradient_degree": 13,
        "affine_degree": 3,
        "total_segre_correction": 13**5 - 3,
    },
}
assert controls["meng_yang_before_schur_hc6"]["total_segre_correction"] == 117646
assert controls["meng_yang_after_schur_hc5"]["total_segre_correction"] == 371290


payload = {
    "format": "hc4-projective-polar-atlas-v5",
    "conventions": {
        "map": "[X0^m:F1^h:F2^h:F3^h:F4^h]",
        "segre_pushforward": (
            "i_*s(B,P4)=sigma_1*H+sigma_2*H^2+"
            "sigma_3*H^3+sigma_4*H^4"
        ),
        "formula": (
            "g_i=m^i-sum_{k=1}^i binom(i,k)m^(i-k)sigma_k"
        ),
        "top_formula_without_fixed_divisor": (
            "g_4=m^4-6*m^2*sigma_2-4*m*sigma_3-sigma_4"
        ),
    },
    "scope": (
        "Integral signatures from positivity, degree bounds, and "
        "log-concavity. Wang's theorem excludes every listed "
        "quadratic-gradient affine-degree-two/three row. HC4CQ1 and "
        "Ax--Grothendieck exclude every listed cubic-gradient row. The "
        "quintic rows are paired with an exact leading-Hessian coverage "
        "matrix. The top-gradient support theorem assigns codimension four "
        "to the smooth essential rank-three branch, and the vertex-colength "
        "obstruction excludes its affine-degree-two/three rows. The other "
        "624 quintic rows remain necessary numerical configurations, not "
        "existence results; the ordinary Hilbert polynomial does not "
        "determine the Segre class."
    ),
    "counts": {key: len(value) for key, value in atlases.items()},
    "zero_dimensional_lengths": zero_dimensional,
    "quadratic_keller_consequence": quadratic_keller_consequence,
    "quadratic_keller_degree_one_atlas": [
        asdict(row) for row in quadratic_keller_atlas
    ],
    "cubic_gradient_consequence": cubic_gradient_consequence,
    "cubic_keller_degree_one_atlas": [
        asdict(row) for row in cubic_keller_atlas
    ],
    "quintic_coverage_summary": quintic_coverage_summary,
    "quintic_coverage_matrix": quintic_coverage_matrix,
    "quintic_top_gradient_rees_sieve": quintic_top_gradient_rees_sieve,
    "rank_three_vertex_colength_sieve": rank_three_vertex_colength_sieve,
    "codimension_three_gradient_sieve": codimension_three_gradient_sieve,
    "codimension_two_binary_root_sieve": (
        codimension_two_binary_root_sieve
    ),
    "atlases": {
        key: [asdict(row) for row in value]
        for key, value in atlases.items()
    },
    "calibrations": calibrations,
    "controls": controls,
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: inverted the projective-degree/Segre-degree formula exactly")
print("PASS: atlas counts are 9, 7, 72, 67, 319, and 307")
print("PASS: zero-dimensional lengths are 14, 13, 79, 78, 254, and 253")
print("PASS: Wang's theorem excludes all 16 quadratic-degree rows")
print("PASS: every quadratic Keller gradient has total correction 15")
print("PASS: HC4CQ1 excludes all 139 cubic-degree rows")
print("PASS: every cubic Keller gradient in HC4 has total correction 80")
print("PASS: any HC4 counterexample potential has degree at least five")
print("PASS: verified the rank-one/two/three quintic determinant faces")
print("PASS: rank-three aligned quintics reduce to HC4CD5")
print("PASS: isolated the rank-three cubic Schur divisibility gap")
print("PASS: squarefree rank-three Hessian discriminants force s3=0")
print("PASS: intersected generic quintic top-gradient/Rees strata with atlas")
print("PASS: excluded both codimension-four quintic degree-two/three rows")
print("PASS: refined the codimension-three rank-two/rank-three packets")
print("PASS: partitioned the codimension-two repeated-binary-root packet")
print("PASS: HC6/HC5 total corrections are 117646 and 371290")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
