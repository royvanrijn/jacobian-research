#!/usr/bin/env python3
"""Exact common primitives for versioned elliptic-curve candidate records.

The candidate key is a mathematical family identifier plus a canonical
rational parameter.  Display strings and legacy parameter conventions are
aliases, never keys.  This module also supplies exact Weierstrass changes and
small-prime finite-quotient certificates.  The latter are deliberately
independent of ``ecsearch.rank_certification``'s cyclic-discrete-log engine.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path
from typing import Any, Iterable, Sequence


Q = Fraction
RationalPoint = tuple[Fraction, Fraction]
FinitePoint = tuple[int, int] | None
WeierstrassModel = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
CANDIDATE_SCHEMA = "elliptic-curves.candidate-record.v1"


def parse_fraction(value: str | int | Fraction) -> Fraction:
    """Parse and canonically reduce an exact rational value."""

    if isinstance(value, bool):
        raise ValueError("booleans are not rational parameters")
    return Q(value)


def fraction_text(value: str | int | Fraction) -> str:
    value = parse_fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def point_record(point: RationalPoint) -> dict[str, str]:
    return {"x": fraction_text(point[0]), "y": fraction_text(point[1])}


def point_from_record(record: dict[str, Any]) -> RationalPoint:
    return parse_fraction(record["x"]), parse_fraction(record["y"])


def model_record(model: Sequence[Fraction | int]) -> list[str]:
    if len(model) != 5:
        raise ValueError("five Weierstrass coefficients are required")
    return [fraction_text(value) for value in model]


def model_from_record(record: Sequence[str | int]) -> WeierstrassModel:
    if len(record) != 5:
        raise ValueError("five Weierstrass coefficients are required")
    return tuple(parse_fraction(value) for value in record)  # type: ignore[return-value]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json_sha256(value: Any) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def point_sequence_sha256(points: Sequence[RationalPoint]) -> str:
    digest = hashlib.sha256()
    for x_coordinate, y_coordinate in points:
        digest.update(
            f"{fraction_text(x_coordinate)}\t{fraction_text(y_coordinate)}\n".encode()
        )
    return digest.hexdigest()


def canonical_candidate_identity(
    family_id: str,
    parameter_name: str,
    parameter: Fraction | int,
    *,
    sign_quotient: bool = False,
) -> dict[str, Any]:
    """Return the stable identity block used by common candidate records."""

    if not family_id or not parameter_name:
        raise ValueError("the family and parameter names must be nonempty")
    parameter = parse_fraction(parameter)
    if sign_quotient:
        parameter = abs(parameter)
    text = fraction_text(parameter)
    return {
        "family_id": family_id,
        "canonical_parameter": {"name": parameter_name, "value": text},
        "sign_quotient": sign_quotient,
        "candidate_key": f"{family_id}:{parameter_name}={text}",
        "raw_parameter_strings_are_aliases_only": True,
    }


def is_on_weierstrass_curve(
    model: Sequence[Fraction | int], point: RationalPoint
) -> bool:
    a1, a2, a3, a4, a6 = model_from_record(model)
    x_coordinate, y_coordinate = map(Q, point)
    return (
        y_coordinate**2 + a1 * x_coordinate * y_coordinate + a3 * y_coordinate
        == x_coordinate**3
        + a2 * x_coordinate**2
        + a4 * x_coordinate
        + a6
    )


def weierstrass_invariants(
    model: Sequence[Fraction | int],
) -> dict[str, Fraction]:
    a1, a2, a3, a4, a6 = model_from_record(model)
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    discriminant = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    return {
        "b2": b2,
        "b4": b4,
        "b6": b6,
        "b8": b8,
        "c4": c4,
        "c6": c6,
        "discriminant": discriminant,
    }


@dataclass(frozen=True)
class WeierstrassChange:
    """A change ``x=u^2*x'+r, y=u^3*y'+s*u^2*x'+t``.

    The unprimed model is the source and the primed model is the target, as
    in PARI's ``ellchangecurve(source,[u,r,s,t])`` convention.
    """

    u: Fraction
    r: Fraction
    s: Fraction
    t: Fraction

    @classmethod
    def from_values(
        cls, values: Sequence[Fraction | int | str]
    ) -> "WeierstrassChange":
        if len(values) != 4:
            raise ValueError("a Weierstrass change requires [u,r,s,t]")
        answer = cls(*(parse_fraction(value) for value in values))
        if answer.u == 0:
            raise ValueError("the scale u must be nonzero")
        return answer

    def to_record(self) -> list[str]:
        return [fraction_text(value) for value in (self.u, self.r, self.s, self.t)]


def change_weierstrass_model(
    source: Sequence[Fraction | int], change: WeierstrassChange
) -> WeierstrassModel:
    """Apply an exact admissible change to a generalized Weierstrass model."""

    a1, a2, a3, a4, a6 = model_from_record(source)
    u, r, s, t = change.u, change.r, change.s, change.t
    return (
        (a1 + 2 * s) / u,
        (a2 - s * a1 + 3 * r - s**2) / u**2,
        (a3 + r * a1 + 2 * t) / u**3,
        (
            a4
            - s * a3
            + 2 * r * a2
            - (t + r * s) * a1
            + 3 * r**2
            - 2 * s * t
        )
        / u**4,
        (
            a6
            + r * a4
            + r**2 * a2
            + r**3
            - t * a3
            - r * t * a1
            - t**2
        )
        / u**6,
    )


def source_point_to_target(
    point: RationalPoint, change: WeierstrassChange
) -> RationalPoint:
    """Transport an affine source point to the target of ``change``."""

    x_source, y_source = map(Q, point)
    u, r, s, t = change.u, change.r, change.s, change.t
    x_target = (x_source - r) / u**2
    y_target = (y_source - s * (x_source - r) - t) / u**3
    return x_target, y_target


def target_point_to_source(
    point: RationalPoint, change: WeierstrassChange
) -> RationalPoint:
    """Transport an affine target point back to the source of ``change``."""

    x_target, y_target = map(Q, point)
    u, r, s, t = change.u, change.r, change.s, change.t
    return (
        u**2 * x_target + r,
        u**3 * y_target + s * u**2 * x_target + t,
    )


def binary_quartic_invariants_low_to_high(
    coefficients: Sequence[Fraction | int],
) -> tuple[Fraction, Fraction]:
    """Return classical ``(I,J)`` for ``e+d*x+c*x^2+b*x^3+a*x^4``."""

    if len(coefficients) != 5:
        raise ValueError("five binary-quartic coefficients are required")
    e, d, c, b, a = map(Q, coefficients)
    invariant_i = 12 * a * e - 3 * b * d + c**2
    invariant_j = (
        72 * a * c * e
        + 9 * b * c * d
        - 27 * a * d**2
        - 27 * b**2 * e
        - 2 * c**3
    )
    return invariant_i, invariant_j


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def primes_up_to(bound: int) -> tuple[int, ...]:
    if bound < 2:
        return ()
    return tuple(value for value in range(2, bound + 1) if _is_prime(value))


def matrix_rank_and_pivots_mod_prime(
    rows: Iterable[Sequence[int]], column_count: int, prime: int
) -> tuple[int, tuple[int, ...]]:
    if not _is_prime(prime):
        raise ValueError("the matrix modulus must be prime")
    matrix = []
    for row in rows:
        normalized = [int(value) % prime for value in row]
        if len(normalized) != column_count:
            raise ValueError("a matrix row has the wrong width")
        if any(normalized):
            matrix.append(normalized)
    rank = 0
    pivots: list[int] = []
    for column in range(column_count):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(value * inverse) % prime for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or matrix[index][column] == 0:
                continue
            multiple = matrix[index][column]
            matrix[index] = [
                (left - multiple * right) % prime
                for left, right in zip(matrix[index], matrix[rank])
            ]
        pivots.append(column)
        rank += 1
        if rank == len(matrix):
            break
    return rank, tuple(pivots)


def _reduce_rational(value: Fraction | int, prime: int) -> int:
    value = Q(value)
    denominator = value.denominator % prime
    if denominator == 0:
        raise ValueError("a rational denominator is not invertible")
    return value.numerator * pow(denominator, -1, prime) % prime


def _finite_add(
    left: FinitePoint, right: FinitePoint, coefficient_a: int, prime: int
) -> FinitePoint:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if (y1 + y2) % prime == 0:
            return None
        slope = (3 * x1**2 + coefficient_a) * pow(2 * y1, -1, prime)
    else:
        slope = (y2 - y1) * pow((x2 - x1) % prime, -1, prime)
    slope %= prime
    x3 = (slope**2 - x1 - x2) % prime
    y3 = (-y1 + slope * (x1 - x3)) % prime
    return x3, y3


def _finite_multiply(
    point: FinitePoint, scalar: int, coefficient_a: int, prime: int
) -> FinitePoint:
    if scalar < 0:
        if point is not None:
            point = point[0], (-point[1]) % prime
        return _finite_multiply(point, -scalar, coefficient_a, prime)
    answer: FinitePoint = None
    addend = point
    while scalar:
        if scalar & 1:
            answer = _finite_add(answer, addend, coefficient_a, prime)
        addend = _finite_add(addend, addend, coefficient_a, prime)
        scalar >>= 1
    return answer


def _finite_subtract(
    left: FinitePoint, right: FinitePoint, coefficient_a: int, prime: int
) -> FinitePoint:
    if right is not None:
        right = right[0], (-right[1]) % prime
    return _finite_add(left, right, coefficient_a, prime)


def _finite_curve_points(
    coefficient_a: int, coefficient_b: int, prime: int
) -> tuple[FinitePoint, ...]:
    roots: dict[int, list[int]] = {}
    for ordinate in range(prime):
        roots.setdefault(ordinate**2 % prime, []).append(ordinate)
    points: list[FinitePoint] = [None]
    for abscissa in range(prime):
        rhs = (abscissa**3 + coefficient_a * abscissa + coefficient_b) % prime
        points.extend((abscissa, ordinate) for ordinate in roots.get(rhs, ()))
    return tuple(points)


def finite_quotient_signature(
    model: Sequence[Fraction | int],
    points: Sequence[RationalPoint],
    reduction_prime: int,
    relation_prime: int,
) -> dict[str, Any]:
    """Compute coordinates in ``E(F_p)/ell E(F_p)`` by enumeration."""

    if not _is_prime(reduction_prime) or reduction_prime == relation_prime:
        raise ValueError("the reduction prime must be prime and differ from ell")
    if not _is_prime(relation_prime):
        raise ValueError("the relation modulus ell must be prime")
    coefficients = model_from_record(model)
    if any(coefficients[:3]):
        raise ValueError("finite quotient certificates require a short model")
    if any(not is_on_weierstrass_curve(coefficients, point) for point in points):
        raise ValueError("a rational point is not on the short curve")
    coefficient_a = _reduce_rational(coefficients[3], reduction_prime)
    coefficient_b = _reduce_rational(coefficients[4], reduction_prime)
    discriminant = -16 * (4 * coefficient_a**3 + 27 * coefficient_b**2)
    if discriminant % reduction_prime == 0:
        raise ValueError("the curve has bad reduction")
    finite_points = _finite_curve_points(
        coefficient_a, coefficient_b, reduction_prime
    )
    multiples = {
        _finite_multiply(point, relation_prime, coefficient_a, reduction_prime)
        for point in finite_points
    }
    span: list[FinitePoint] = [None]
    coordinates: list[tuple[int, ...]] = [()]
    basis: list[FinitePoint] = []
    for point in finite_points:
        if any(
            _finite_subtract(point, representative, coefficient_a, reduction_prime)
            in multiples
            for representative in span
        ):
            continue
        basis.append(point)
        old_span = tuple(span)
        old_coordinates = tuple(coordinates)
        span = []
        coordinates = []
        for scalar in range(relation_prime):
            multiple = _finite_multiply(
                point, scalar, coefficient_a, reduction_prime
            )
            span.extend(
                _finite_add(
                    representative, multiple, coefficient_a, reduction_prime
                )
                for representative in old_span
            )
            coordinates.extend(
                coordinate + (scalar,) for coordinate in old_coordinates
            )
    if len(span) * len(multiples) != len(finite_points):
        raise AssertionError("the finite quotient representatives do not cover")
    quotient_order = len(span)
    quotient_dimension = 0
    while quotient_order > 1 and quotient_order % relation_prime == 0:
        quotient_order //= relation_prime
        quotient_dimension += 1
    if quotient_order != 1 or quotient_dimension != len(basis):
        raise AssertionError("the finite quotient is not an ell-vector space")
    rows = [[0] * len(points) for _ in basis]
    for point_index, rational_point in enumerate(points):
        reduced = (
            _reduce_rational(rational_point[0], reduction_prime),
            _reduce_rational(rational_point[1], reduction_prime),
        )
        coordinate_index = next(
            (
                index
                for index, representative in enumerate(span)
                if _finite_subtract(
                    reduced, representative, coefficient_a, reduction_prime
                )
                in multiples
            ),
            None,
        )
        if coordinate_index is None:
            raise AssertionError("a reduced rational point missed every quotient coset")
        for basis_index, value in enumerate(coordinates[coordinate_index]):
            rows[basis_index][point_index] = value
    return {
        "prime": reduction_prime,
        "group_order": len(finite_points),
        "multiple_subgroup_order": len(multiples),
        "quotient_dimension": quotient_dimension,
        "rows": rows,
    }


def build_finite_quotient_certificate(
    model: Sequence[Fraction | int],
    points: Sequence[RationalPoint],
    *,
    relation_prime: int,
    prime_bound: int,
) -> dict[str, Any]:
    """Greedily build an exact, deterministic mod-``ell`` certificate.

    A full combined rank plus a good-reduction group order prime to ``ell``
    proves rational independence by infinite descent.  A rank-deficient
    result is only a bounded negative result for this certificate method.
    """

    if not points:
        raise ValueError("at least one rational point is required")
    if not _is_prime(relation_prime):
        raise ValueError("the relation modulus must be prime")
    if prime_bound < 3:
        raise ValueError("the reduction-prime bound must be at least three")
    coefficients = model_from_record(model)
    rational_points = tuple((Q(x), Q(y)) for x, y in points)
    if any(not is_on_weierstrass_curve(coefficients, point) for point in rational_points):
        raise ValueError("a rational point is not on the supplied model")
    selected: list[dict[str, Any]] = []
    rows: list[list[int]] = []
    rank = 0
    torsion_witness: dict[str, int] | None = None
    for reduction_prime in primes_up_to(prime_bound):
        if reduction_prime in (2, relation_prime):
            continue
        try:
            signature = finite_quotient_signature(
                coefficients, rational_points, reduction_prime, relation_prime
            )
        except ValueError:
            continue
        if torsion_witness is None and signature["group_order"] % relation_prime:
            torsion_witness = {
                "prime": reduction_prime,
                "group_order": signature["group_order"],
            }
        candidate_rows = [*rows, *signature["rows"]]
        candidate_rank, _ = matrix_rank_and_pivots_mod_prime(
            candidate_rows, len(rational_points), relation_prime
        )
        if candidate_rank > rank:
            selected.append(signature)
            rows.extend(signature["rows"])
            rank = candidate_rank
        if rank == len(rational_points) and torsion_witness is not None:
            break
    rank, pivots = matrix_rank_and_pivots_mod_prime(
        rows, len(rational_points), relation_prime
    )
    certified = rank == len(rational_points) and torsion_witness is not None
    return {
        "certificate_type": "finite-quotient-infinite-descent",
        "implementation": "elliptic-curves/cas/elliptic_candidate_record.py",
        "relation_prime": relation_prime,
        "reduction_prime_bound": prime_bound,
        "point_count": len(rational_points),
        "point_sequence_sha256": point_sequence_sha256(rational_points),
        "combined_rank_over_relation_field": rank,
        "pivot_columns_zero_based": list(pivots),
        "torsion_witness": torsion_witness,
        "signatures": selected,
        "certificate_primes": [item["prime"] for item in selected],
        "certified_independent": certified,
        "certified_rank_lower_bound": len(rational_points) if certified else None,
        "rank_deficient_interpretation": (
            None
            if certified
            else "bounded certificate-method failure; not evidence of dependence"
        ),
    }


def verify_finite_quotient_certificate(
    model: Sequence[Fraction | int],
    points: Sequence[RationalPoint],
    certificate: dict[str, Any],
) -> None:
    """Replay every exact finite-field assertion in a common certificate."""

    relation_prime = int(certificate["relation_prime"])
    rational_points = tuple((Q(x), Q(y)) for x, y in points)
    if certificate["point_count"] != len(rational_points):
        raise AssertionError("the certificate point count changed")
    if certificate["point_sequence_sha256"] != point_sequence_sha256(rational_points):
        raise AssertionError("the certificate point sequence changed")
    rows: list[list[int]] = []
    for signature in certificate["signatures"]:
        replay = finite_quotient_signature(
            model, rational_points, int(signature["prime"]), relation_prime
        )
        if replay != signature:
            raise AssertionError("a finite quotient signature changed")
        rows.extend(signature["rows"])
    rank, pivots = matrix_rank_and_pivots_mod_prime(
        rows, len(rational_points), relation_prime
    )
    if rank != certificate["combined_rank_over_relation_field"]:
        raise AssertionError("the combined finite quotient rank changed")
    if list(pivots) != certificate["pivot_columns_zero_based"]:
        raise AssertionError("the finite quotient pivots changed")
    witness = certificate["torsion_witness"]
    if witness is not None:
        replay = finite_quotient_signature(
            model, rational_points, int(witness["prime"]), relation_prime
        )
        if replay["group_order"] != witness["group_order"]:
            raise AssertionError("the torsion-witness group order changed")
        if witness["group_order"] % relation_prime == 0:
            raise AssertionError("the witness does not exclude rational ell-torsion")
    certified = rank == len(rational_points) and witness is not None
    if bool(certificate["certified_independent"]) != certified:
        raise AssertionError("the certificate conclusion changed")


def validate_candidate_identity(record: dict[str, Any]) -> None:
    """Validate the common schema's stable identity invariants."""

    if record.get("schema") != CANDIDATE_SCHEMA:
        raise AssertionError("unknown common candidate schema")
    identity = record["identity"]
    parameter = identity["canonical_parameter"]
    expected = canonical_candidate_identity(
        identity["family_id"],
        parameter["name"],
        parse_fraction(parameter["value"]),
        sign_quotient=bool(identity["sign_quotient"]),
    )
    for key, value in expected.items():
        if identity.get(key) != value:
            raise AssertionError(f"candidate identity field {key!r} changed")
