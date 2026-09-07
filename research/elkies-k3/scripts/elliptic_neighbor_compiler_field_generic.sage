#!/usr/bin/env sage
"""Field-generic compatibility layer for elliptic_neighbor_compiler.sage.

Q80's CM24 third-q12 local gates naturally live over
QQ(sqrt(-6),sqrt(-3)).  The core compiler historically hard-coded QQ in
`quotient_condition`, `resolved_chart_quotient_condition`, and
`compile_resolved_conditions`, although `finite_ambient_image_condition` was
already field-generic.

Load this file instead of the core compiler when exact local modules use a
number field.  It loads the core first and then replaces only those three
entry points.  Historical QQ callers retain the same behavior.

This compatibility layer is intentionally small so it can later be folded
into the core file without changing callers.
"""

from pathlib import Path

from sage.all import QQ, matrix, vector

HERE = Path(__file__).resolve().parent
CORE = HERE / "elliptic_neighbor_compiler.sage"
exec(compile(CORE.read_text(), str(CORE), "exec"))


def quotient_condition(
    name,
    ambient_basis,
    evaluator,
    quotient_basis,
    provenance,
    coefficient_field=QQ,
):
    """Build one exact finite-quotient block over `coefficient_field`."""
    ambient_basis = tuple(ambient_basis)
    quotient_basis = tuple(quotient_basis)
    columns = []
    for basis_element in ambient_basis:
        residue = vector(coefficient_field, evaluator(basis_element))
        if len(residue) != len(quotient_basis):
            raise ValueError("{} returned a residue of the wrong length".format(name))
        columns.append(residue)
    return {
        "name": str(name),
        "matrix": matrix(coefficient_field, columns).transpose(),
        "quotient_basis": tuple(map(str, quotient_basis)),
        "provenance": str(provenance),
    }


def resolved_chart_quotient_condition(
    name,
    ambient_basis,
    local_ring,
    trivialized_pullback,
    quotient_ideal,
    quotient_basis,
    provenance,
):
    """Compile one actual resolved quotient over the local ring's base field."""
    coefficient_field = local_ring.base_ring()
    quotient_basis = tuple(local_ring(value) for value in quotient_basis)
    groebner_basis = quotient_ideal.groebner_basis()
    if not quotient_basis:
        raise ValueError("{} has an empty quotient basis".format(name))
    try:
        quotient_dimension = quotient_ideal.vector_space_dimension()
    except AttributeError:
        if local_ring.ngens() != 1 or len(groebner_basis) != 1:
            raise ValueError(
                "{} needs a finite-colength quotient with computable dimension".format(name)
            )
        quotient_dimension = local_ring(groebner_basis[0]).degree()
        if quotient_dimension <= 0:
            raise ValueError("{} quotient is not finite-dimensional".format(name))
    if quotient_dimension != len(quotient_basis):
        raise ValueError(
            "{} quotient basis has length {}, but quotient dimension is {}".format(
                name, len(quotient_basis), quotient_dimension
            )
        )

    def normal_form(value):
        value = local_ring(value)
        try:
            return value.reduce(groebner_basis)
        except AttributeError:
            if local_ring.ngens() != 1 or len(groebner_basis) != 1:
                raise
            return value.mod(local_ring(groebner_basis[0]))

    if any(normal_form(value) != local_ring(value) for value in quotient_basis):
        raise ValueError("{} quotient basis is not reduced".format(name))

    def evaluator(basis_element):
        remainder = normal_form(trivialized_pullback(basis_element))
        coordinates = vector(
            coefficient_field,
            tuple(
                remainder.monomial_coefficient(monomial)
                for monomial in quotient_basis
            ),
        )
        reconstructed = sum(
            coefficient*monomial
            for coefficient, monomial in zip(coordinates, quotient_basis)
        )
        if remainder != reconstructed:
            raise ValueError(
                "{} quotient basis does not span the pullback remainder".format(name)
            )
        return coordinates

    return quotient_condition(
        name,
        ambient_basis,
        evaluator,
        quotient_basis,
        provenance,
        coefficient_field=coefficient_field,
    )


def _common_exact_field(blocks, coefficient_field=None):
    if coefficient_field is not None:
        return coefficient_field
    fields = tuple(block["matrix"].base_ring() for block in blocks)
    if not fields:
        return QQ
    result = QQ
    for field in fields:
        if result.has_coerce_map_from(field):
            continue
        if field.has_coerce_map_from(result):
            result = field
            continue
        raise ValueError(
            "condition blocks do not have a common exact coefficient field: {} and {}".format(
                result, field
            )
        )
    return result


def compile_resolved_conditions(
    ambient_basis,
    condition_blocks,
    complete=False,
    compute_kernel=True,
    coefficient_field=None,
):
    """Stack exact local blocks over their inferred compatible exact field."""
    ambient_basis = tuple(ambient_basis)
    width = len(ambient_basis)
    blocks = tuple(condition_blocks)
    for block in blocks:
        if block["matrix"].ncols() != width:
            raise ValueError(
                "condition block {} has incompatible width".format(block["name"])
            )
    coefficient_field = _common_exact_field(blocks, coefficient_field)
    condition_matrix = matrix(coefficient_field, 0, width)
    for block in blocks:
        condition_matrix = condition_matrix.stack(
            matrix(coefficient_field, block["matrix"])
        )
    rank = condition_matrix.rank()
    nullity = width-rank
    kernel_basis = None
    kernel_materialization = "not_requested"
    if compute_kernel:
        zero_columns = tuple(
            column for column in range(width)
            if all(
                condition_matrix[row, column] == 0
                for row in range(condition_matrix.nrows())
            )
        )
        active_columns = tuple(
            column for column in range(width) if column not in zero_columns
        )
        coordinate_kernel = matrix(
            coefficient_field,
            len(zero_columns),
            width,
            lambda row, column: (
                coefficient_field(1)
                if column == zero_columns[row]
                else coefficient_field(0)
            ),
        )
        active_nullity = len(active_columns)-rank
        if active_nullity < 0:
            raise ArithmeticError("condition rank exceeds its active columns")
        if active_nullity == 0:
            kernel_basis = coordinate_kernel
            kernel_materialization = "zero_columns"
        else:
            active_matrix = condition_matrix.matrix_from_columns(active_columns)
            active_kernel = active_matrix.right_kernel()
            if active_kernel.dimension() != active_nullity:
                raise ArithmeticError("active rank and right-kernel dimensions disagree")
            active_basis = active_kernel.basis_matrix()
            active_index = {
                column: index for index, column in enumerate(active_columns)
            }
            lifted_kernel = matrix(
                coefficient_field,
                active_kernel.dimension(),
                width,
                lambda row, column: (
                    active_basis[row, active_index[column]]
                    if column in active_index
                    else coefficient_field(0)
                ),
            )
            kernel_basis = coordinate_kernel.stack(lifted_kernel)
            kernel_materialization = "zero_columns_plus_reduced_right_kernel"
        zero = matrix(coefficient_field, condition_matrix.nrows(), nullity)
        if (
            kernel_basis.nrows() != nullity
            or condition_matrix*kernel_basis.transpose() != zero
        ):
            raise ArithmeticError(
                "materialized kernel does not match the exact condition matrix"
            )
    return {
        "ambient_dimension": width,
        "condition_rows": condition_matrix.nrows(),
        "rank": rank,
        "codimension": rank,
        "kernel_dimension": nullity,
        "kernel_basis": kernel_basis,
        "kernel_materialization": kernel_materialization,
        "complete_resolved_chart_cover": bool(complete),
        "h0_certified": bool(complete and nullity == 2),
        "condition_matrix": condition_matrix,
        "coefficient_field": coefficient_field,
    }
