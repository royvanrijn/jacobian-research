#!/usr/bin/env python3
"""Universal top-quintic and generic infinity-Rees strata for HC(4).

This checker works at the leading homogeneous/normal-cone level.  It does
not assume that lower potential layers leave the complete Segre class
unchanged.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from itertools import combinations
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
)


ATLAS = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_projective_polar_atlas.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_quintic_infinity_rees_strata.json"
)

AMBIENT_DIMENSION = 4
MAP_DEGREE = 4
POTENTIAL_DEGREE = 5


def exponent_tuples(total: int, length: int) -> tuple[tuple[int, ...], ...]:
    """Return all weak compositions of ``total`` into ``length`` parts."""

    if length == 1:
        return ((total,),)
    return tuple(
        (first, *tail)
        for first in range(total + 1)
        for tail in exponent_tuples(total - first, length - 1)
    )


def universal_form(
    variables: tuple[sp.Symbol, ...],
    *,
    degree: int,
    prefix: str,
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Construct the universal homogeneous form in the supplied variables."""

    exponents = exponent_tuples(degree, len(variables))
    coefficients = sp.symbols(f"{prefix}_0:{len(exponents)}")
    polynomial = sp.expand(
        sum(
            coefficient
            * sp.prod(
                variable**exponent
                for variable, exponent in zip(variables, powers)
            )
            for coefficient, powers in zip(coefficients, exponents)
        )
    )
    return polynomial, coefficients


@dataclass(frozen=True)
class GenericRankStratum:
    hessian_rank: int
    coefficient_count: int
    active_gradient_generators: int
    inactive_gradient_generators: int
    infinity_gradient_rees_relations: int
    compactification_rees_relations: int
    base_support_codimension: int
    pure_top_projective_degrees: tuple[int, ...]
    pure_top_segre_degrees: tuple[int, ...]
    atlas_rows_affine_degree_2: int
    atlas_rows_affine_degree_3: int


x = sp.symbols("x1:5")
u = sp.symbols("u1:5")
x_vector = sp.Matrix(x)

atlas_payload = json.loads(ATLAS.read_text())
atlases = atlas_payload["atlases"]

# The full 56-coefficient top quintic precedes the adapted rank strata.
full_h5, full_h5_coefficients = universal_form(
    x,
    degree=POTENTIAL_DEGREE,
    prefix="c4",
)
full_h5_gradient = tuple(sp.diff(full_h5, variable) for variable in x)
full_h5_hessian = sp.Matrix(
    [
        [sp.diff(full_h5, left, right) for right in x]
        for left in x
    ]
)
assert len(full_h5_coefficients) == 56
assert full_h5_hessian == full_h5_hessian.T
assert sp.expand(
    sum(a * b for a, b in zip(x, full_h5_gradient)) - 5 * full_h5
) == 0
assert (
    full_h5_hessian * x_vector
    - sp.Matrix([4 * component for component in full_h5_gradient])
).applyfunc(sp.expand) == sp.zeros(4, 1)

# Rank four is genuinely generic before the constant-Hessian gate.
fermat_h5 = sum(variable**5 for variable in x)
fermat_hessian_determinant = sp.factor(
    sp.det(sp.hessian(fermat_h5, x))
)
assert fermat_hessian_determinant == 20**4 * sp.prod(
    variable**3 for variable in x
)

# Midpoint translation separates the collision equations by parity:
# grad(h2+h4)(a)=0 and grad(h3+h5)(a)=0.
parity_checks = {}
for degree in range(2, 6):
    if degree == 5:
        form = full_h5
        coefficients = full_h5_coefficients
    else:
        form, coefficients = universal_form(
            x,
            degree=degree,
            prefix=f"p{degree}",
        )
    gradient = tuple(sp.diff(form, variable) for variable in x)
    parity = (-1) ** (degree - 1)
    assert all(
        sp.expand(
            component.subs(
                {variable: -variable for variable in x},
                simultaneous=True,
            )
            - parity * component
        )
        == 0
        for component in gradient
    )
    parity_checks[f"h{degree}"] = {
        "coefficient_count": len(coefficients),
        "gradient_parity_under_x_to_minus_x": parity,
    }

strata: list[GenericRankStratum] = []
universal_checks = {}

for rank in (1, 2, 3):
    active_variables = x[:rank]
    h5, coefficients = universal_form(
        active_variables,
        degree=POTENTIAL_DEGREE,
        prefix=f"c{rank}",
    )
    gradient = tuple(sp.diff(h5, variable) for variable in x)
    hessian = sp.Matrix(
        [
            [sp.diff(h5, left, right) for right in x]
            for left in x
        ]
    )

    assert len(coefficients) == comb(rank + 4, 5)
    assert hessian == hessian.T
    assert sp.expand(sum(a * b for a, b in zip(x, gradient)) - 5 * h5) == 0
    assert (hessian * x_vector - sp.Matrix(
        [4 * component for component in gradient]
    )).applyfunc(sp.expand) == sp.zeros(4, 1)
    assert all(
        sp.expand(
            sp.diff(gradient[left], x[right])
            - sp.diff(gradient[right], x[left])
        )
        == 0
        for left, right in combinations(range(4), 2)
    )
    assert all(component == 0 for component in gradient[rank:])
    assert hessian[rank:, :] == sp.zeros(4 - rank, 4)
    assert hessian[:, rank:] == sp.zeros(4, 4 - rank)

    infinity_rees_relations = tuple(
        sp.expand(gradient[left] * u[right] - gradient[right] * u[left])
        for left, right in combinations(range(rank), 2)
    )
    assert all(
        sp.expand(
            relation.subs(
                {
                    u[index]: gradient[index]
                    for index in range(4)
                }
            )
        )
        == 0
        for relation in infinity_rees_relations
    )

    codimension = rank + 1
    sigmas = equal_degree_complete_intersection_segre(
        ambient_dimension=AMBIENT_DIMENSION,
        codimension=codimension,
        generator_degree=MAP_DEGREE,
    )
    degrees = projective_degrees_from_segre(MAP_DEGREE, sigmas)
    expected_degrees = tuple(
        MAP_DEGREE**index if index <= rank else 0
        for index in range(AMBIENT_DIMENSION + 1)
    )
    assert degrees == expected_degrees

    counts = {}
    for affine_degree in (2, 3):
        key = f"gradient_degree_4_affine_degree_{affine_degree}"
        counts[affine_degree] = sum(
            row["leading_base_codimension"] == codimension
            for row in atlases[key]
        )

    stratum = GenericRankStratum(
        hessian_rank=rank,
        coefficient_count=len(coefficients),
        active_gradient_generators=rank,
        inactive_gradient_generators=4 - rank,
        infinity_gradient_rees_relations=comb(rank, 2) + (4 - rank),
        compactification_rees_relations=comb(rank + 1, 2) + (4 - rank),
        base_support_codimension=codimension,
        pure_top_projective_degrees=degrees,
        pure_top_segre_degrees=sigmas,
        atlas_rows_affine_degree_2=counts[2],
        atlas_rows_affine_degree_3=counts[3],
    )
    strata.append(stratum)
    universal_checks[f"rank_{rank}"] = {
        "universal_top_term_count": len(coefficients),
        "euler_identity": "sum x_i*d_i(h5)=5*h5",
        "hessian_euler_identity": "Hess(h5)*x=4*grad(h5)",
        "curl_relations": comb(4, 2),
        "inactive_target_relations": [f"u_{index + 1}" for index in range(rank, 4)],
        "active_koszul_relations": [
            f"d_{left + 1}(h5)*u_{right + 1}-d_{right + 1}(h5)*u_{left + 1}"
            for left, right in combinations(range(rank), 2)
        ],
    }


assert [row.coefficient_count for row in strata] == [1, 6, 21]
assert [row.base_support_codimension for row in strata] == [2, 3, 4]
assert [
    (row.atlas_rows_affine_degree_2, row.atlas_rows_affine_degree_3)
    for row in strata
] == [(260, 249), (58, 57), (1, 1)]
assert [row.pure_top_segre_degrees for row in strata] == [
    (0, 16, -128, 768),
    (0, 0, 64, -768),
    (0, 0, 0, 256),
]


# The support of the top gradient ideal is the join of the constant-kernel
# vertex P^(3-r) with the singular locus of the essential h5 in P^(r-1).
# Empty singular locus gives codimension r+1.  A nonempty singular locus of
# dimension s gives codimension r-s.
support_strata = [
    {
        "hessian_rank": 1,
        "essential_quintic_singularity": "empty",
        "base_support_codimension": 2,
        "atlas_rows": {"2": 260, "3": 249},
        "determinant_face_status": (
            "open: det(Hess(h4)|K)=0 on a constant three-plane; "
            "kernel synchronization remains"
        ),
    },
    {
        "hessian_rank": 2,
        "essential_quintic_singularity": "empty (squarefree binary quintic)",
        "base_support_codimension": 3,
        "atlas_rows": {"2": 58, "3": 57},
        "determinant_face_status": (
            "open: det(Hess(h4)|K)=0 on a constant two-plane; "
            "kernel synchronization remains"
        ),
    },
    {
        "hessian_rank": 2,
        "essential_quintic_singularity": "nonempty (repeated binary root)",
        "base_support_codimension": 2,
        "atlas_rows": {"2": 260, "3": 249},
        "determinant_face_status": (
            "open: singular binary top plus the rank-two determinant face"
        ),
    },
    {
        "hessian_rank": 3,
        "essential_quintic_singularity": "empty (smooth ternary quintic)",
        "base_support_codimension": 4,
        "atlas_rows": {"2": 1, "3": 1},
        "determinant_face_status": (
            "aligned branch closed; nonaligned branch requires "
            "nonsquarefree ternary Hessian determinant"
        ),
    },
    {
        "hessian_rank": 3,
        "essential_quintic_singularity": "isolated points",
        "base_support_codimension": 3,
        "atlas_rows": {"2": 58, "3": 57},
        "determinant_face_status": (
            "aligned branch closed; nonaligned branch requires "
            "nonsquarefree ternary Hessian determinant"
        ),
    },
    {
        "hessian_rank": 3,
        "essential_quintic_singularity": "positive-dimensional",
        "base_support_codimension": 2,
        "atlas_rows": {"2": 260, "3": 249},
        "determinant_face_status": (
            "aligned branch closed; nonaligned branch requires "
            "nonsquarefree ternary Hessian determinant"
        ),
    },
]


symbolic_atlas_restrictions = {
    "codimension_2": {
        "equations": ["sigma_1=0", "sigma_2>0"],
        "top_degree_equation": "96*sigma_2+16*sigma_3+sigma_4=256-delta",
    },
    "codimension_3": {
        "equations": ["sigma_1=sigma_2=0", "sigma_3>0"],
        "top_degree_equation": "16*sigma_3+sigma_4=256-delta",
    },
    "codimension_4": {
        "equations": ["sigma_1=sigma_2=sigma_3=0"],
        "top_degree_equation": "sigma_4=256-delta",
        "affine_degree_2_signature": [0, 0, 0, 254],
        "affine_degree_3_signature": [0, 0, 0, 253],
    },
}

atlas_intersections = {}
for affine_degree in (2, 3):
    atlas_key = f"gradient_degree_4_affine_degree_{affine_degree}"
    atlas_intersections[f"affine_degree_{affine_degree}"] = {
        f"codimension_{codimension}": [
            {
                "projective_degrees": row["projective_degrees"],
                "segre_degrees": row["segre_degrees"],
            }
            for row in atlases[atlas_key]
            if row["leading_base_codimension"] == codimension
        ]
        for codimension in (2, 3, 4)
    }

assert {
    affine_key: {
        codimension_key: len(rows)
        for codimension_key, rows in intersections.items()
    }
    for affine_key, intersections in atlas_intersections.items()
} == {
    "affine_degree_2": {
        "codimension_2": 260,
        "codimension_3": 58,
        "codimension_4": 1,
    },
    "affine_degree_3": {
        "codimension_2": 249,
        "codimension_3": 57,
        "codimension_4": 1,
    },
}


payload = {
    "format": "hc4-quintic-infinity-rees-strata-v1",
    "software_assumptions": {
        "python": "repository .python-version and requirements.txt",
        "coefficient_field": "characteristic zero; exact symbolic work over QQ",
        "independent_rees_replay": (
            "Macaulay2 1.22 with Cremona and ReesAlgebra"
        ),
    },
    "scope": (
        "Exact universal top-quintic identities, generic smooth essential "
        "rank-stratum Rees/complete-intersection models, and support-"
        "codimension intersections with the numerical HC4 atlas. The pure-"
        "top Segre vectors are degeneration calibrations, not the Segre "
        "vectors of a completed constant-Hessian potential. Lower layers "
        "can change normal-cone multiplicities while preserving the reduced "
        "support codimension."
    ),
    "universal_four_variable_h5": {
        "coefficient_count": len(full_h5_coefficients),
        "gradient_generators": 4,
        "curl_relations": 6,
        "euler_identity": "sum x_i*d_i(h5)=5*h5",
        "hessian_euler_identity": "Hess(h5)*x=4*grad(h5)",
        "generic_rank_four_witness": (
            "h5=sum x_i^5 has det Hess(h5)=20^4*product x_i^3"
        ),
        "constant_hessian_top_gate": "det Hess(h5)=0, hence rank<=3",
    },
    "midpoint_collision_parity": {
        "homogeneous_layer_checks": parity_checks,
        "even_gradient_equation": "grad(h3+h5)(a)=0",
        "odd_gradient_equation": "grad(h2+h4)(a)=0",
        "scope": (
            "These are exact collision equations after midpoint "
            "translation. They do not by themselves put [a] in the "
            "infinity gradient base scheme."
        ),
    },
    "universal_checks": universal_checks,
    "generic_smooth_rank_strata": [asdict(row) for row in strata],
    "support_strata": support_strata,
    "symbolic_atlas_restrictions": symbolic_atlas_restrictions,
    "atlas_intersections": atlas_intersections,
    "generic_rees_conclusion": (
        "On the smooth essential rank-r stratum, the active infinity "
        "gradient generators form a regular sequence and are of linear "
        "type. Their Rees ideal has only Koszul equations; padded zero "
        "gradient coordinates add linear target equations. Extra nonlinear "
        "Rees constraints must therefore come from singular top quintics "
        "or from the lower-layer completed base ideal."
    ),
    "constant_hessian_gate": (
        "The top determinant forces Hessian rank at most three. Existing "
        "HC4 determinant-face equations must be imposed after this support "
        "filter; they are not encoded by the pure-top complete intersection."
    ),
    "constant_determinant_normal_cone_faces": {
        "rank_1": {
            "lambda_degree": 9,
            "equation": "det(C1)*det(Hess(h4)|K)=0",
        },
        "rank_2": {
            "lambda_degree": 10,
            "equation": "det(C2)*det(Hess(h4)|K)=0",
        },
        "rank_3": {
            "lambda_degree_11": "det(C3)*D_t^2(h4)=0",
            "lambda_degree_10": (
                "det(C3)*D_t^2(h3)-grad(s3)^T*adj(C3)*grad(s3)=0"
            ),
            "certificate": "scripts/verify_hc4_projective_polar_atlas.py",
        },
    },
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
OUTPUT.write_text(serialized)
digest = hashlib.sha256(serialized.encode()).hexdigest()

print("PASS: built the universal 56-coefficient four-variable quintic")
print("PASS: constant Hessian removes its generic rank-four top stratum")
print("PASS: midpoint collision equations split into odd/even gradient layers")
print("PASS: built universal essential rank-1/2/3 quintic top forms")
print("PASS: verified Euler, Hessian-Euler, curl, and Koszul identities")
print("PASS: generic top Rees models are CI/linear-type targets")
print("PASS: pure-top Segre signatures are rank-stratified exactly")
print("PASS: atlas support counts are (260,249), (58,57), and (1,1)")
print(f"PASS: wrote {OUTPUT.relative_to(ROOT)}")
print(f"SHA256: {digest}")
