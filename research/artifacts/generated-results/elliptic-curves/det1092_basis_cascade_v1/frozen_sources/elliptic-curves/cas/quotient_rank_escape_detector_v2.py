#!/usr/bin/env python3
"""Exact GF(2) measurement layer for quotient-aware rank-escape detector v2.

This module does not perform a number-field descent.  It consumes the exact
global norm-squareclass basis, all-place local-condition rows, and point
Kummer images emitted by a completed descent.  It then produces the canonical
row-space, quotient dimensions, and leave-one-place-out measurements without
choosing or pairing an arbitrary quotient complement.

Rows are linear equations on a fixed, ordered basis of the global norm-square
subspace.  Point images are vectors in that same ambient space.  The only
canonical complement used below is a deterministic presentation device; no
bilinear form is asserted to descend to it.

status: ACTIVE_PROOF
claim: exact post-descent quotient and all-place matrix measurements
inputs: completed global norm basis, all local condition rows, and point images
outputs: canonical GF(2) detector record supplied by the caller
supersedes/superseded-by: none
"""

from __future__ import annotations

from hashlib import sha256
import json
from typing import Any, Iterable, Mapping, Sequence


class DetectorInputError(ValueError):
    """A completed-descent payload is inconsistent or incomplete."""


def _bits(row: Sequence[int], width: int, label: str = "row") -> list[int]:
    result = [int(value) for value in row]
    if len(result) != width or any(value not in (0, 1) for value in result):
        raise DetectorInputError(f"{label} is not a binary row of width {width}")
    return result


def canonical_rref(rows: Iterable[Sequence[int]], width: int) -> list[list[int]]:
    """Return the unique reduced row-echelon basis with leftmost pivots."""

    if width < 0:
        raise DetectorInputError("matrix width must be nonnegative")
    matrix = [_bits(row, width) for row in rows]
    pivot_row = 0
    for column in range(width):
        source = next(
            (index for index in range(pivot_row, len(matrix)) if matrix[index][column]),
            None,
        )
        if source is None:
            continue
        matrix[pivot_row], matrix[source] = matrix[source], matrix[pivot_row]
        for index in range(len(matrix)):
            if index != pivot_row and matrix[index][column]:
                matrix[index] = [
                    left ^ right
                    for left, right in zip(matrix[index], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return [row for row in matrix[:pivot_row] if any(row)]


def f2_rank(rows: Iterable[Sequence[int]], width: int) -> int:
    return len(canonical_rref(rows, width))


def pivot_columns(rref: Sequence[Sequence[int]]) -> list[int]:
    return [next(index for index, value in enumerate(row) if value) for row in rref]


def nullspace_basis(equations: Iterable[Sequence[int]], width: int) -> list[list[int]]:
    """Return a canonical basis of the common kernel of equation rows."""

    reduced = canonical_rref(equations, width)
    pivots = pivot_columns(reduced)
    free = [column for column in range(width) if column not in set(pivots)]
    basis = []
    for free_column in free:
        vector = [0] * width
        vector[free_column] = 1
        for row, pivot in zip(reduced, pivots):
            vector[pivot] = row[free_column]
        basis.append(vector)
    return basis


def dot(left: Sequence[int], right: Sequence[int]) -> int:
    if len(left) != len(right):
        raise DetectorInputError("dot-product rows have different widths")
    return sum(int(a) * int(b) for a, b in zip(left, right)) & 1


def reduce_mod_rowspace(row: Sequence[int], basis: Iterable[Sequence[int]]) -> list[int]:
    """Return the canonical coset representative modulo a row space."""

    row = [int(value) for value in row]
    reduced_basis = canonical_rref(basis, len(row))
    result = _bits(row, len(row))
    for pivot_row, pivot in zip(reduced_basis, pivot_columns(reduced_basis)):
        if result[pivot]:
            result = [left ^ right for left, right in zip(result, pivot_row)]
    return result


def canonical_quotient_presentation(
    subspace_rows: Iterable[Sequence[int]],
    quotient_by_rows: Iterable[Sequence[int]],
    width: int,
) -> list[list[int]]:
    """Choose deterministic representatives spanning ``subspace/quotient_by``.

    This is a presentation tied to the declared ambient coordinates, not an
    intrinsic splitting of the quotient.
    """

    subspace = canonical_rref(subspace_rows, width)
    quotient = canonical_rref(quotient_by_rows, width)
    if f2_rank([*subspace, *quotient], width) != len(subspace):
        raise DetectorInputError("the quotienting rows are not in the subspace")
    current = list(quotient)
    current_rank = len(current)
    representatives = []
    for row in subspace:
        if f2_rank([*current, row], width) > current_rank:
            representative = reduce_mod_rowspace(row, quotient)
            representatives.append(representative)
            current.append(row)
            current_rank += 1
    if current_rank != len(subspace):
        raise DetectorInputError("failed to present the complete quotient")
    return canonical_rref(representatives, width)


def _place_key(place: str) -> tuple[int, int | str]:
    if place == "infinity":
        return (1, "infinity")
    try:
        value = int(place)
    except (TypeError, ValueError) as error:
        raise DetectorInputError(f"invalid place label: {place!r}") from error
    if value < 2:
        raise DetectorInputError(f"invalid finite place: {place!r}")
    return (0, value)


def _matrix_hash(rows: Sequence[Sequence[int]]) -> str:
    canonical = json.dumps(rows, separators=(",", ":"))
    return sha256(canonical.encode()).hexdigest()


def _require_points_in_selmer(
    label: str,
    points: Sequence[Sequence[int]],
    equations: Sequence[Sequence[int]],
    width: int,
) -> list[list[int]]:
    checked = [_bits(row, width, label) for row in points]
    for point_index, point in enumerate(checked, start=1):
        if any(dot(equation, point) for equation in equations):
            raise DetectorInputError(
                f"{label} point {point_index} violates a global local condition"
            )
    return checked


def analyze_complete_descent(
    *,
    ambient_dimension: int,
    local_condition_rows: Mapping[str, Sequence[Sequence[int]]],
    required_places: Sequence[str],
    mw17_rows: Sequence[Sequence[int]],
    exceptional_rows: Sequence[Sequence[int]] | None = None,
    local_metadata: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Measure one completed all-place 2-descent exactly.

    ``exceptional_rows`` must be omitted during the blind phase.  Supplying it
    later labels how much of the already-frozen blind quotient is realized by
    the twelve record points.
    """

    if ambient_dimension < 0:
        raise DetectorInputError("ambient dimension must be nonnegative")
    if not local_condition_rows:
        raise DetectorInputError("all-place local-condition rows are required")
    metadata = {} if local_metadata is None else dict(local_metadata)
    places = sorted((str(place) for place in local_condition_rows), key=_place_key)
    required = [str(place) for place in required_places]
    if len(required) != len(set(required)):
        raise DetectorInputError("the required-place declaration contains duplicates")
    required = sorted(required, key=_place_key)
    if places != required:
        missing = sorted(set(required) - set(places), key=_place_key)
        extra = sorted(set(places) - set(required), key=_place_key)
        raise DetectorInputError(
            f"local-condition place set mismatch: missing={missing}, extra={extra}"
        )
    if "2" not in places or "infinity" not in places:
        raise DetectorInputError("places 2 and infinity must both be explicit")

    per_place_rows = {
        place: canonical_rref(local_condition_rows[place], ambient_dimension)
        for place in places
    }
    stacked = canonical_rref(
        (row for place in places for row in per_place_rows[place]),
        ambient_dimension,
    )
    selmer_basis = nullspace_basis(stacked, ambient_dimension)
    selmer_dimension = len(selmer_basis)

    known17 = _require_points_in_selmer(
        "MW17", mw17_rows, stacked, ambient_dimension
    )
    mw17_dimension = f2_rank(known17, ambient_dimension)
    if len(known17) != 17 or mw17_dimension != 17:
        raise DetectorInputError("the actual MW17 image must have dimension seventeen")
    if selmer_dimension < mw17_dimension:
        raise DetectorInputError("the Selmer dimension is smaller than the MW17 image")

    blind_quotient = canonical_quotient_presentation(
        selmer_basis, known17, ambient_dimension
    )
    residual_dimension = selmer_dimension - mw17_dimension
    if len(blind_quotient) != residual_dimension:
        raise DetectorInputError("the blind quotient presentation has wrong dimension")

    cumulative_rows: list[list[int]] = []
    cumulative_rank = 0
    place_records = []
    for place in places:
        rows = per_place_rows[place]
        place_rank = len(rows)
        new_rank = f2_rank([*cumulative_rows, *rows], ambient_dimension)
        without = canonical_rref(
            (
                row
                for retained in places
                if retained != place
                for row in per_place_rows[retained]
            ),
            ambient_dimension,
        )
        delete_dimension = ambient_dimension - len(without)
        if delete_dimension < mw17_dimension:
            raise DetectorInputError("deleting a place lost the known MW17 image")
        suppression = delete_dimension - selmer_dimension
        record = {
            "place": place,
            "canonical_condition_row_space": rows,
            "canonical_condition_row_space_sha256": _matrix_hash(rows),
            "local_condition_codimension_on_global_norm_space": place_rank,
            "incremental_independent_codimension": new_rank - cumulative_rank,
            "matrix_rank_after_deleting_this_place": len(without),
            "s_res_after_deleting_this_place": delete_dimension - mw17_dimension,
            "single_place_suppression_of_residual_intersection": suppression,
        }
        if place in metadata:
            record["local_metadata"] = metadata[place]
        place_records.append(record)
        cumulative_rows.extend(rows)
        cumulative_rank = new_rank

    result: dict[str, Any] = {
        "ambient": "ordered basis of the global norm-square S-squareclass subspace",
        "ambient_dimension": ambient_dimension,
        "two_selmer_dimension": selmer_dimension,
        "actual_mw17_image_dimension": mw17_dimension,
        "s_res": residual_dimension,
        "blind_quotient_presentation": {
            "basis_rows_in_global_norm_coordinates": blind_quotient,
            "dimension": residual_dimension,
            "warning": (
                "deterministic ambient-coordinate representatives only; no "
                "pairing or canonical splitting is inferred"
            ),
        },
        "global_local_condition_matrix": {
            "canonical_rref_rows": stacked,
            "canonical_rref_sha256": _matrix_hash(stacked),
            "rank": len(stacked),
            "rank_nullity_verified": len(stacked) + selmer_dimension
            == ambient_dimension,
            "summed_placewise_codimension": sum(len(per_place_rows[p]) for p in places),
            "independent_all_place_codimension": len(stacked),
        },
        "places": place_records,
        "maximum_single_place_suppression": max(
            row["single_place_suppression_of_residual_intersection"]
            for row in place_records
        ),
        "required_places_verified": required,
        "all_place_local_conditions_explicit": True,
        "pairing_claim": None,
    }

    if exceptional_rows is not None:
        known12 = _require_points_in_selmer(
            "exceptional", exceptional_rows, stacked, ambient_dimension
        )
        if len(known12) != 12:
            raise DetectorInputError("the record calibration requires twelve exceptional points")
        quotient_rows = [reduce_mod_rowspace(row, known17) for row in known12]
        exceptional_dimension = f2_rank(quotient_rows, ambient_dimension)
        mw29_dimension = f2_rank([*known17, *known12], ambient_dimension)
        if exceptional_dimension != 12 or mw29_dimension != 29:
            raise DetectorInputError(
                "the twelve exceptional images are not independent modulo MW17"
            )
        result["held_out_record_calibration"] = {
            "exceptional_generator_count": 12,
            "exceptional_rows_modulo_mw17": quotient_rows,
            "exceptional_image_dimension_modulo_mw17": exceptional_dimension,
            "actual_mw29_image_dimension": mw29_dimension,
            "selmer_modulo_mw29_dimension": selmer_dimension - mw29_dimension,
            "blind_quotient_contains_complete_exceptional_span": (
                f2_rank([*known17, *blind_quotient, *known12], ambient_dimension)
                == mw17_dimension + residual_dimension
            ),
            "exact_rank_29_if_selmer_dimension_29": selmer_dimension == 29,
        }
    return result


def checkpointed_simon_condition_rows(selmer: Mapping[str, Any]) -> dict[str, list[list[int]]]:
    """Convert a completed checkpointed-Simon output to v2 equation rows.

    The conversion works entirely over GF(2).  It takes the inverse image of
    each stored allowed S-squareclass subspace inside the global norm subspace,
    then returns that inverse image's annihilator in norm-basis coordinates.
    """

    matrix = selmer.get("local_condition_matrix")
    if not isinstance(matrix, Mapping):
        raise DetectorInputError("missing checkpointed local-condition matrix")
    norm_columns = matrix.get(
        "global_norm_square_subspace_basis_columns_in_s_squareclasses"
    )
    global_dimension = int(matrix.get("global_s_squareclass_dimension", -1))
    norm_dimension = int(matrix.get("global_norm_square_subspace_dimension", -1))
    if not isinstance(norm_columns, list) or len(norm_columns) != norm_dimension:
        raise DetectorInputError("invalid stored global norm-subspace basis")
    norm_columns = [
        _bits(column, global_dimension, "norm-subspace column")
        for column in norm_columns
    ]

    result: dict[str, list[list[int]]] = {}
    for place_record in matrix.get("places", []):
        place = str(place_record["place"])
        allowed_columns = place_record.get(
            "allowed_subspace_basis_columns_in_global_s_squareclasses"
        )
        if not isinstance(allowed_columns, list):
            raise DetectorInputError(f"missing allowed subspace at {place}")
        allowed_columns = [
            _bits(column, global_dimension, f"allowed column at {place}")
            for column in allowed_columns
        ]
        variable_width = norm_dimension + len(allowed_columns)
        equations = []
        for coordinate in range(global_dimension):
            equations.append(
                [column[coordinate] for column in norm_columns]
                + [column[coordinate] for column in allowed_columns]
            )
        joint_kernel = nullspace_basis(equations, variable_width)
        allowed_in_norm_coordinates = canonical_rref(
            (row[:norm_dimension] for row in joint_kernel), norm_dimension
        )
        result[place] = nullspace_basis(
            allowed_in_norm_coordinates, norm_dimension
        )
        expected = int(
            place_record["norm_subspace_intersection_dimension_for_this_place_alone"]
        )
        if norm_dimension - len(result[place]) != expected:
            raise DetectorInputError(
                f"checkpointed intersection dimension mismatch at {place}"
            )
    return result


def checkpointed_point_rows_in_normspace(
    selmer: Mapping[str, Any], embedding: Mapping[str, Any]
) -> list[list[int]]:
    """Transport point coordinates from a Selmer basis to norm coordinates."""

    matrix = selmer.get("local_condition_matrix")
    if not isinstance(matrix, Mapping):
        raise DetectorInputError("missing checkpointed local-condition matrix")
    columns = matrix.get("selmer_basis_columns_in_global_norm_square_subspace")
    norm_dimension = int(matrix.get("global_norm_square_subspace_dimension", -1))
    selmer_dimension = int(selmer.get("two_selmer_dimension", -1))
    if not isinstance(columns, list) or len(columns) != selmer_dimension:
        raise DetectorInputError(
            "the completed descent lacks Selmer-basis columns in norm coordinates"
        )
    columns = [
        _bits(column, norm_dimension, "Selmer basis column") for column in columns
    ]
    point_rows = embedding.get("point_selmer_rows")
    if not isinstance(point_rows, list):
        raise DetectorInputError("missing point-to-Selmer embedding rows")
    result = []
    for point_index, coordinates in enumerate(point_rows, start=1):
        coordinates = _bits(
            coordinates, selmer_dimension, f"point {point_index} Selmer coordinates"
        )
        result.append(
            [
                sum(bit * column[ambient] for bit, column in zip(coordinates, columns))
                & 1
                for ambient in range(norm_dimension)
            ]
        )
    return result


def checkpointed_local_metadata(selmer: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    """Extract exact local Kummer dimensions from a completed v2 worker."""

    matrix = selmer.get("local_condition_matrix")
    if not isinstance(matrix, Mapping):
        raise DetectorInputError("missing checkpointed local-condition matrix")
    result = {}
    for row in matrix.get("places", []):
        required = (
            "ambient_local_kummer_dimension",
            "computed_local_kummer_image_dimension",
            "localized_global_s_squareclass_image_dimension",
        )
        if any(name not in row for name in required):
            raise DetectorInputError(
                f"the completed worker lacks v2 local dimensions at {row.get('place')}"
            )
        metadata: dict[str, Any] = {
            "ambient_local_kummer_dimension": int(
                row["ambient_local_kummer_dimension"]
            ),
            "local_kummer_image_dimension": int(
                row["computed_local_kummer_image_dimension"]
            ),
            "localized_global_s_squareclass_image_dimension": int(
                row["localized_global_s_squareclass_image_dimension"]
            ),
        }
        for name in ("elliptic_bad_place", "auxiliary_descent_place"):
            if name in row:
                metadata[name] = bool(row[name])
        if row.get("component_group_data") is not None:
            metadata["component_group_data"] = dict(row["component_group_data"])
        result[str(row["place"])] = metadata
    return result


def analyze_checkpointed_simon(
    selmer: Mapping[str, Any],
    generic_embedding: Mapping[str, Any],
    exceptional_embedding: Mapping[str, Any] | None = None,
    component_metadata: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Run detector v2 on completed checkpointed-Simon worker outputs."""

    matrix = selmer.get("local_condition_matrix")
    if not isinstance(matrix, Mapping):
        raise DetectorInputError("missing checkpointed local-condition matrix")
    local = checkpointed_local_metadata(selmer)
    if component_metadata:
        for place, values in component_metadata.items():
            local.setdefault(str(place), {}).update(dict(values))
    return analyze_complete_descent(
        ambient_dimension=int(matrix["global_norm_square_subspace_dimension"]),
        local_condition_rows=checkpointed_simon_condition_rows(selmer),
        required_places=[
            *selmer.get("finite_local_condition_primes", []),
            "infinity",
        ],
        mw17_rows=checkpointed_point_rows_in_normspace(
            selmer, generic_embedding
        ),
        exceptional_rows=(
            None
            if exceptional_embedding is None
            else checkpointed_point_rows_in_normspace(
                selmer, exceptional_embedding
            )
        ),
        local_metadata=local,
    )


def matrix_rank_profile(measurement: Mapping[str, Any]) -> list[list[int]]:
    """Return the predeclared cross-curve comparable rank profile.

    Complete row spaces live in different cubic-field squareclass spaces, so a
    direct Grassmann distance between curves is undefined.  This profile keeps
    only basis-independent per-place ranks: local codimension, single-place
    suppression, and local Kummer-image dimension (or -1 if absent).
    """

    profile = []
    for place in measurement["places"]:
        local = place.get("local_metadata", {})
        profile.append(
            [
                int(place["local_condition_codimension_on_global_norm_space"]),
                int(place["single_place_suppression_of_residual_intersection"]),
                int(local.get("local_kummer_image_dimension", -1)),
            ]
        )
    return sorted(profile)


def matrix_rank_profile_distance(
    left: Mapping[str, Any], right: Mapping[str, Any]
) -> int:
    """L1 distance between sorted, zero-padded all-place rank profiles."""

    a = matrix_rank_profile(left)
    b = matrix_rank_profile(right)
    width = max(len(a), len(b))
    a.extend([[0, 0, 0]] * (width - len(a)))
    b.extend([[0, 0, 0]] * (width - len(b)))
    return sum(abs(x - y) for row_a, row_b in zip(a, b) for x, y in zip(row_a, row_b))
