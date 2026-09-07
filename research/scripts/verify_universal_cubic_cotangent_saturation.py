#!/usr/bin/env python3
"""Exact formal-gauge certificate for universal cotangent saturation.

For the smooth ternary-cubic symbol, let K be the graded module of
compatible tensor corrections

    (c_ijk) with z*c_(0ij)-y*c_(1ij)+x*c_(2ij)=0.

A homogeneous 3-by-3 matrix D defines a simultaneous infinitesimal change
of collision coordinates and coefficient-module generators.  If

    r=(z,-y,x)^T,  D*r=(V_z,-V_y,V_x)^T,

then x |-> x+V_x, y |-> y+V_y, z |-> z+V_z and e |-> (I+D)e preserve the
Koszul presentation to first order.  The tensor takes values in det(M), so
transport back to the original determinant line contributes
-tr(D)*phi.  Differentiating the smooth cubic tensor with this determinant
twist gives a graded A-linear gauge map G:A^9(-3)->K.

This checker proves the exact module identity

    K = im(G) + A*eta,       (x,y,z)*eta subset im(G),

where eta is the tensor of 3*X*Y*Z, the one cubic-modulus direction not
induced by an infinitesimal gauge.  Consequently K_d=im(G)_d for every
d>=4.  Successive homogeneous gauge
changes therefore remove every compatible term above the cubic symbol in
the (x,y,z)-adic completion.  The universal quartic family is formally
equivalent to its central homogeneous fiber.  The central cotangent
presentation is independently checked to be saturated, so faithful
detection of (x,y,z)-power torsion by completion proves the universal
claim.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import shutil
import subprocess
from pathlib import Path

import sympy as sp

import research_universal_cubic_quartic_kernel_saturation as frontier
import verify_cubic_symbol_double_saturation as cubic_audit


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "universal_cubic_cotangent_saturation.json"
)
TRIPLES = tuple(
    itertools.combinations_with_replacement(range(3), 3)
)
PAIRS = tuple(
    itertools.combinations_with_replacement(range(3), 2)
)


def central_tensor() -> sp.Matrix:
    """Return the ten components of the smooth cubic tensor."""

    cubic = cubic_audit.CUBIC_STRATA["smooth"]
    return sp.Matrix(
        [
            cubic_audit.polarized_value(
                cubic,
                *(
                    cubic_audit.RELATION.cross(
                        cubic_audit.STANDARD_BASIS[index]
                    )
                    for index in triple
                ),
            )
            for triple in TRIPLES
        ]
    )


def modulus_tensor() -> sp.Matrix:
    """Return an integral generator of the cubic gauge cokernel."""

    return 3 * sp.Matrix(
        [
            cubic_audit.polarized_value(
                cubic_audit.X * cubic_audit.Y * cubic_audit.Z,
                *(
                    cubic_audit.RELATION.cross(
                        cubic_audit.STANDARD_BASIS[index]
                    )
                    for index in triple
                ),
            )
            for triple in TRIPLES
        ]
    )


def compatibility_matrix() -> sp.Matrix:
    """Return C:A^10 -> A^6 whose kernel is the tensor module K."""

    columns: list[sp.Matrix] = []
    relation_terms = (
        (cubic_audit.z, 1),
        (cubic_audit.y, -1),
        (cubic_audit.x, 1),
    )
    for triple in TRIPLES:
        column = sp.zeros(6, 1)
        for pair_index, pair in enumerate(PAIRS):
            for first, (variable, sign) in enumerate(relation_terms):
                if tuple(sorted((first, *pair))) == triple:
                    column[pair_index] += sign * variable
        columns.append(column)
    return sp.Matrix.hstack(*columns)


def gauge_matrix(phi: sp.Matrix) -> sp.Matrix:
    """Return the determinant-twisted formal-gauge differential.

    The transported tensor is

        det(I+D)^(-1) * sigma(phi)((I+D)e)^3.

    Its linear term therefore contains -tr(D)*phi in addition to the
    coordinate and three input-slot variations.
    """

    variables = cubic_audit.BASE_VARIABLES
    columns: list[sp.Matrix] = []
    for row in range(3):
        for column in range(3):
            endomorphism = sp.zeros(3, 3)
            endomorphism[row, column] = 1
            relation_change = endomorphism * cubic_audit.RELATION
            coordinate_change = sp.Matrix(
                (
                    relation_change[2],
                    -relation_change[1],
                    relation_change[0],
                )
            )
            values: list[sp.Expr] = []
            for triple_index, triple in enumerate(TRIPLES):
                value = sum(
                    coordinate_change[index]
                    * sp.diff(phi[triple_index], variables[index])
                    for index in range(3)
                )
                for position in range(3):
                    for replacement in range(3):
                        changed = list(triple)
                        changed[position] = replacement
                        changed_index = TRIPLES.index(
                            tuple(sorted(changed))
                        )
                        value += (
                            endomorphism[
                                replacement, triple[position]
                            ]
                            * phi[changed_index]
                        )
                value -= sp.trace(endomorphism) * phi[triple_index]
                values.append(sp.expand(value))
            columns.append(sp.Matrix(values))
    return sp.Matrix.hstack(*columns)


def explicit_maximal_action_lift() -> sp.Matrix:
    """Return L with G*L=[x*eta,y*eta,z*eta]."""

    x, y, z = cubic_audit.BASE_VARIABLES
    return sp.Matrix(
        [
            [0, 0, 0],
            [-y, 0, 0],
            [0, -x, 0],
            [z, 0, 0],
            [0, 0, 0],
            [0, 0, x],
            [0, -z, 0],
            [0, 0, -y],
            [0, 0, 0],
        ]
    )


def verify_dual_number_gauge_action(
    phi: sp.Matrix,
    gauge: sp.Matrix,
) -> None:
    """Derive every column of G from the exact determinant-twisted action."""

    epsilon = sp.Symbol("gauge_epsilon")
    variables = cubic_audit.BASE_VARIABLES
    relation = cubic_audit.RELATION
    for gauge_column, (row, column) in enumerate(
        itertools.product(range(3), repeat=2)
    ):
        endomorphism = sp.zeros(3, 3)
        endomorphism[row, column] = 1
        relation_change = endomorphism * relation
        coordinate_change = sp.Matrix(
            (
                relation_change[2],
                -relation_change[1],
                relation_change[0],
            )
        )
        basis_change = sp.eye(3) + epsilon * endomorphism
        transformed_variables = (
            sp.Matrix(variables) + epsilon * coordinate_change
        )
        assert (
            basis_change * relation
            - sp.Matrix(
                (
                    transformed_variables[2],
                    -transformed_variables[1],
                    transformed_variables[0],
                )
            )
        ).applyfunc(sp.expand) == sp.zeros(3, 1)
        substitution = dict(zip(variables, transformed_variables))
        transformed_components: list[sp.Expr] = []
        for triple in TRIPLES:
            value = sp.Integer(0)
            for first, second, third in itertools.product(
                range(3), repeat=3
            ):
                changed_index = TRIPLES.index(
                    tuple(sorted((first, second, third)))
                )
                value += (
                    basis_change[first, triple[0]]
                    * basis_change[second, triple[1]]
                    * basis_change[third, triple[2]]
                    * phi[changed_index].subs(
                        substitution, simultaneous=True
                    )
                )
            # This is det(I+epsilon*D)^(-1) in Q[epsilon]/(epsilon^2).
            value *= 1 - epsilon * sp.trace(endomorphism)
            transformed_components.append(
                sp.expand(value).coeff(epsilon, 1)
            )
        assert (
            sp.Matrix(transformed_components) - gauge[:, gauge_column]
        ).applyfunc(sp.expand) == sp.zeros(10, 1)


def universal_quartic_gauge_lift(
    gauge: sp.Matrix,
    parameters: tuple[sp.Symbol, ...],
    universal: dict[tuple[int, int, int], sp.Expr],
) -> tuple[sp.Matrix, sp.Matrix, dict[str, int]]:
    """Construct Q with G*Q equal to all 24 universal quartic directions."""

    x, y, z = cubic_audit.BASE_VARIABLES
    variables = cubic_audit.BASE_VARIABLES
    monomials = cubic_audit.homogeneous_monomials(4)
    directions = sp.Matrix.hstack(
        *[
            sp.Matrix(
                [
                    sp.expand(universal[triple]).coeff(parameter)
                    for triple in TRIPLES
                ]
            )
            for parameter in parameters
        ]
    )
    for entry in directions:
        if entry == 0:
            continue
        polynomial = sp.Poly(sp.expand(entry), *variables)
        assert all(
            sum(monomial) == 4
            for monomial, coefficient in polynomial.terms()
            if coefficient
        )

    def coefficient_column(vector: sp.Matrix) -> sp.Matrix:
        return sp.Matrix(
            [
                sp.Poly(
                    sp.expand(vector[row]), *variables
                ).coeff_monomial(monomial)
                for row in range(vector.rows)
                for monomial in monomials
            ]
        )

    action_matrix = sp.Matrix.hstack(
        *[
            coefficient_column(gauge[:, gauge_column] * variable)
            for gauge_column in range(gauge.cols)
            for variable in variables
        ]
    )
    direction_matrix = sp.Matrix.hstack(
        *[
            coefficient_column(directions[:, column])
            for column in range(directions.cols)
        ]
    )
    action_rank = action_matrix.rank()
    direction_rank = direction_matrix.rank()
    assert action_rank == direction_rank == 24
    solution, _free_parameters = action_matrix.gauss_jordan_solve(
        direction_matrix
    )
    free_symbols = set().union(
        *(entry.free_symbols for entry in solution)
    )
    solution = solution.subs({symbol: 0 for symbol in free_symbols})
    assert action_matrix * solution == direction_matrix
    lift = sp.Matrix(
        gauge.cols,
        directions.cols,
        lambda row, column: sp.expand(
            solution[3 * row, column] * x
            + solution[3 * row + 1, column] * y
            + solution[3 * row + 2, column] * z
        ),
    )
    assert (gauge * lift - directions).applyfunc(
        sp.expand
    ) == sp.zeros(10, len(parameters))
    return directions, lift, {
        "quartic_compatible_space_dimension": direction_rank,
        "linear_gauge_action_rank": action_rank,
        "linear_gauge_kernel_dimension": action_matrix.cols - action_rank,
        "explicit_lift_nonzero_entries": sum(entry != 0 for entry in lift),
    }


def singular_vector(vector: sp.Matrix) -> str:
    return cubic_audit.singular_vector(list(vector))


def singular_module_certificate(
    compatibility: sp.Matrix,
    gauge: sp.Matrix,
    eta: sp.Matrix,
) -> dict[str, int]:
    """Prove ker(C)=im(G)+A*eta and K/im(G)=Q."""

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    compatibility_columns = [
        compatibility[:, index]
        for index in range(compatibility.cols)
    ]
    gauge_columns = [
        gauge[:, index] for index in range(gauge.cols)
    ]
    program = f"""
ring R=0,(x,y,z),dp;
module C={",".join(map(singular_vector, compatibility_columns))};
module G={",".join(map(singular_vector, gauge_columns))};
vector eta={singular_vector(eta)};
module K=syz(C);
module H=G,eta;
H=std(H);
module kernel_difference=simplify(reduce(K,H),2);
module quotient_presentation=std(modulo(K,G));
vector eta_remainder=reduce(eta,std(G));
module maximal_times_quotient=
  (x+y+z)*freemodule(size(K));
module maximal_action=simplify(
  reduce(maximal_times_quotient,quotient_presentation),2
);
print("KERNEL_GENERATORS="+string(size(K)));
print("KERNEL_DIFFERENCE="+string(size(kernel_difference)));
print("QUOTIENT_VECTOR_DIMENSION="+string(vdim(quotient_presentation)));
print("ETA_REMAINDER_SIZE="+string(size(eta_remainder)));
print("MAXIMAL_ACTION="+string(size(maximal_action)));
quit;
"""
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    wanted = {
        "KERNEL_GENERATORS",
        "KERNEL_DIFFERENCE",
        "QUOTIENT_VECTOR_DIMENSION",
        "ETA_REMAINDER_SIZE",
        "MAXIMAL_ACTION",
    }
    values: dict[str, int] = {}
    for line in result.stdout.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key in wanted:
            values[key] = int(value)
    assert set(values) == wanted, result.stdout + result.stderr
    assert values == {
        "KERNEL_GENERATORS": 10,
        "KERNEL_DIFFERENCE": 0,
        "QUOTIENT_VECTOR_DIMENSION": 1,
        "ETA_REMAINDER_SIZE": 6,
        "MAXIMAL_ACTION": 0,
    }, values
    return values


def matrix_record(matrix: sp.Matrix) -> list[list[str]]:
    return [
        [sp.sstr(matrix[row, column]) for column in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def main() -> None:
    cubic_audit.FACTOR_SINGULAR_EXPRESSIONS = False
    phi = central_tensor()
    eta = modulus_tensor()
    compatibility = compatibility_matrix()
    gauge = gauge_matrix(phi)
    maximal_lift = explicit_maximal_action_lift()
    verify_dual_number_gauge_action(phi, gauge)

    assert (
        compatibility * phi
    ).applyfunc(sp.expand) == sp.zeros(6, 1)
    assert (
        compatibility * eta
    ).applyfunc(sp.expand) == sp.zeros(6, 1)
    assert (
        compatibility * gauge
    ).applyfunc(sp.expand) == sp.zeros(6, 9)
    assert (
        gauge * maximal_lift
        - sp.Matrix.hstack(
            *[
                variable * eta
                for variable in cubic_audit.BASE_VARIABLES
            ]
        )
    ).applyfunc(sp.expand) == sp.zeros(10, 3)

    # The serialized universal quartic tensor is an element of K_4.
    parameters, universal = frontier.universal_tensor()
    universal_vector = sp.Matrix([universal[triple] for triple in TRIPLES])
    assert (
        compatibility * universal_vector
    ).applyfunc(sp.expand) == sp.zeros(6, 1)
    quartic_directions, quartic_lift, quartic_lift_data = (
        universal_quartic_gauge_lift(gauge, parameters, universal)
    )

    module_certificate = singular_module_certificate(
        compatibility, gauge, eta
    )

    # This is the independent exact central input transported by the
    # formal gauge equivalence.
    central_result = cubic_audit.run_singular(
        cubic_audit.CUBIC_STRATA["smooth"]
    )
    assert central_result[0] == 0

    exact_data = {
        "compatibility_matrix": matrix_record(compatibility),
        "gauge_matrix": matrix_record(gauge),
        "central_tensor": matrix_record(phi),
        "modulus_tensor": matrix_record(eta),
        "maximal_action_lift": matrix_record(maximal_lift),
        "universal_quartic_directions": matrix_record(
            quartic_directions
        ),
        "universal_quartic_gauge_lift": matrix_record(quartic_lift),
    }
    exact_sha256 = hashlib.sha256(
        json.dumps(
            exact_data, sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    artifact = {
        "schema": "universal-cubic-cotangent-saturation-v1",
        "mathematical_status": "exact formal-gauge theorem",
        "base_ring": "Q[u1,...,u24,x,y,z]",
        "basis_conventions": {
            "collision_variables": [
                str(variable)
                for variable in cubic_audit.BASE_VARIABLES
            ],
            "tensor_component_order": [list(triple) for triple in TRIPLES],
            "gauge_matrix_unit_order": [
                [row, column]
                for row, column in itertools.product(range(3), repeat=2)
            ],
            "quartic_parameter_order": [
                str(parameter) for parameter in parameters
            ],
            "coefficient_field": "Q",
            "singular_monomial_order": "dp",
        },
        "identities": {
            "gauge_convention": (
                "det(I+D)^(-1)*sigma_D(phi)((I+D)e)^3"
            ),
            "dual_number_gauge_action": (
                "the exact action modulo epsilon^2 equals G "
                "for all nine matrix units"
            ),
            "compatibility_times_gauge": "0",
            "compatibility_times_central_tensor": "0",
            "gauge_times_maximal_lift": "[x*eta,y*eta,z*eta]",
            "compatible_tensor_module": "ker(C)=im(G)+Q[x,y,z]*eta",
            "gauge_cokernel": "ker(C)/im(G)=Q concentrated in degree 3",
            "universal_quartic_gauge_lift": "G*Q=[psi_1,...,psi_24]",
        },
        "module_certificate": module_certificate,
        "exact_matrices": exact_data,
        "central_cotangent_result": {
            "saturation_generators": central_result[0],
            "support_dimension": central_result[1],
            "ext3_vector_dimension": central_result[2],
            "ext2_dimension": central_result[3],
            "ext2_vector_dimension": central_result[4],
            "ext2_top_dimension": central_result[5],
            "ext2_square_action_generators": central_result[6],
        },
        "universal_quartic_parameter_count": len(parameters),
        "universal_quartic_gauge_lift_data": quartic_lift_data,
        "exact_matrix_sha256": exact_sha256,
        "proved": [
            (
                "the determinant-twisted finite gauge action over dual "
                "numbers has differential G for all nine matrix units"
            ),
            (
                "the 24 universal quartic directions form the full "
                "compatible quartic space and admit the displayed "
                "linear-polynomial gauge lift"
            ),
            (
                "ker(C)/im(G) is one-dimensional in collision degree "
                "three and vanishes in every higher degree"
            ),
            (
                "formal gauge rigidity and torsion-faithful completion "
                "transfer central cotangent saturation to the universal "
                "smooth quartic family"
            ),
        ],
        "not_proved": [
            "the corresponding formal rigidity for singular cubic symbols",
            "normality of a nonhomogeneous lift",
            "existence of a compatible Keller open",
        ],
        "consequence": (
            "all compatible terms above degree three are removed "
            "successively in the collision-adic completion; universal "
            "relative cotangent saturation follows from the central fiber"
        ),
        "reproduce": (
            ".venv/bin/python "
            "scripts/verify_universal_cubic_cotangent_saturation.py"
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("PASS exact gauge identities C*G=0 and G*L=m*eta")
    print("PASS determinant-twisted dual-number action equals G")
    print("PASS explicit linear gauge lift removes all 24 quartic directions")
    print("PASS ker(C)=im(G)+A*eta and m*eta is in im(G)")
    print("PASS central cotangent presentation is saturated")
    print(
        "PASS universal cubic cotangent saturation by formal gauge "
        "rigidity and torsion-faithful completion"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
