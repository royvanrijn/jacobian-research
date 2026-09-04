"""Exact finite-reduction certificates for independence of rational points.

If the reductions of ``n`` rational points have a full-rank discrete-log
matrix modulo a prime ``ell``, every integral relation has all coefficients
divisible by ``ell``.  A separate good-reduction prime whose group order is
not divisible by ``ell`` rules out rational ``ell``-torsion.  Dividing a
hypothetical relation and repeating then proves independence by infinite
descent.  The certificate below records enough elementary finite-field data
to replay this argument without trusting a numerical height calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Mapping, Sequence

from .fermigier import weierstrass_discriminant


AffinePoint = tuple[Fraction, Fraction]
FinitePoint = tuple[int, int] | None
WeierstrassModel = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]


@dataclass(frozen=True)
class ReductionRow:
    prime: int
    group_order: int
    generator: tuple[int, int]
    logs: tuple[int, ...]

    def to_json_object(self) -> dict[str, object]:
        return {
            "prime": self.prime,
            "group_order": self.group_order,
            "generator": list(self.generator),
            "logs": list(self.logs),
        }

    @classmethod
    def from_json_object(cls, value: Mapping[str, object]) -> "ReductionRow":
        generator = value["generator"]
        logs = value["logs"]
        if not isinstance(generator, list) or not isinstance(logs, list):
            raise ValueError("malformed reduction row")
        return cls(
            prime=int(value["prime"]),
            group_order=int(value["group_order"]),
            generator=(int(generator[0]), int(generator[1])),
            logs=tuple(map(int, logs)),
        )


@dataclass(frozen=True)
class IndependenceCertificate:
    relation_prime: int
    torsion_witness_prime: int
    torsion_witness_group_order: int
    rows: tuple[ReductionRow, ...]

    def to_json_object(self) -> dict[str, object]:
        return {
            "relation_prime": self.relation_prime,
            "torsion_witness": {
                "prime": self.torsion_witness_prime,
                "group_order": self.torsion_witness_group_order,
            },
            "rows": [row.to_json_object() for row in self.rows],
        }

    @classmethod
    def from_json_object(
        cls, value: Mapping[str, object]
    ) -> "IndependenceCertificate":
        witness = value["torsion_witness"]
        rows = value["rows"]
        if not isinstance(witness, Mapping) or not isinstance(rows, list):
            raise ValueError("malformed independence certificate")
        return cls(
            relation_prime=int(value["relation_prime"]),
            torsion_witness_prime=int(witness["prime"]),
            torsion_witness_group_order=int(witness["group_order"]),
            rows=tuple(ReductionRow.from_json_object(row) for row in rows),
        )


def _model(coefficients: Sequence[Fraction | int]) -> WeierstrassModel:
    if len(coefficients) != 5:
        raise ValueError("five Weierstrass coefficients are required")
    return tuple(map(Fraction, coefficients))  # type: ignore[return-value]


def is_on_weierstrass_curve(
    coefficients: Sequence[Fraction | int],
    point: AffinePoint,
) -> bool:
    a1, a2, a3, a4, a6 = _model(coefficients)
    x_coordinate, y_coordinate = map(Fraction, point)
    return (
        y_coordinate * y_coordinate
        + a1 * x_coordinate * y_coordinate
        + a3 * y_coordinate
        == x_coordinate**3
        + a2 * x_coordinate**2
        + a4 * x_coordinate
        + a6
    )


def negate_rational_point(
    coefficients: Sequence[Fraction | int],
    point: AffinePoint | None,
) -> AffinePoint | None:
    if point is None:
        return None
    a1, _, a3, _, _ = _model(coefficients)
    x_coordinate, y_coordinate = point
    return x_coordinate, -y_coordinate - a1 * x_coordinate - a3


def add_rational_points(
    coefficients: Sequence[Fraction | int],
    left: AffinePoint | None,
    right: AffinePoint | None,
) -> AffinePoint | None:
    """Add points exactly on a generalized Weierstrass model."""

    model = _model(coefficients)
    a1, a2, a3, a4, _ = model
    if left is None:
        return right
    if right is None:
        return left
    if not is_on_weierstrass_curve(model, left) or not is_on_weierstrass_curve(
        model, right
    ):
        raise ValueError("both operands must lie on the curve")
    x1, y1 = map(Fraction, left)
    x2, y2 = map(Fraction, right)
    if x1 == x2 and y2 == -y1 - a1 * x1 - a3:
        return None
    if left == right:
        denominator = 2 * y1 + a1 * x1 + a3
        if denominator == 0:
            return None
        slope = (3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1) / denominator
    else:
        if x1 == x2:
            raise ArithmeticError("invalid same-abscissa point pair")
        slope = (y2 - y1) / (x2 - x1)
    intercept = y1 - slope * x1
    x3 = slope * slope + a1 * slope - a2 - x1 - x2
    y3 = -(slope + a1) * x3 - intercept - a3
    result = (x3, y3)
    assert is_on_weierstrass_curve(model, result)
    return result


def subtract_rational_points(
    coefficients: Sequence[Fraction | int],
    left: AffinePoint | None,
    right: AffinePoint | None,
) -> AffinePoint | None:
    return add_rational_points(coefficients, left, negate_rational_point(coefficients, right))


def matrix_rank_mod_prime(matrix: Sequence[Sequence[int]], prime: int) -> int:
    """Return the exact row rank over ``F_prime``."""

    if not _is_prime(prime):
        raise ValueError("the modulus must be prime")
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix rows must have equal length")
    work = [[entry % prime for entry in row] for row in matrix]
    rank = 0
    for column in range(width):
        pivot = next(
            (row for row in range(rank, len(work)) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        inverse = pow(work[rank][column], -1, prime)
        work[rank] = [(entry * inverse) % prime for entry in work[rank]]
        for row in range(len(work)):
            if row == rank or work[row][column] == 0:
                continue
            multiple = work[row][column]
            work[row] = [
                (entry - multiple * pivot_entry) % prime
                for entry, pivot_entry in zip(work[row], work[rank])
            ]
        rank += 1
        if rank == len(work):
            break
    return rank


def _is_prime(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def _primes_between(lower: int, upper: int) -> Iterable[int]:
    for candidate in range(max(2, lower), upper + 1):
        if _is_prime(candidate):
            yield candidate


def _prime_divisors(value: int) -> tuple[int, ...]:
    result: list[int] = []
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            result.append(divisor)
            while value % divisor == 0:
                value //= divisor
        divisor += 1 if divisor == 2 else 2
    if value > 1:
        result.append(value)
    return tuple(result)


def _fraction_mod_prime(value: Fraction | int, prime: int) -> int:
    value = Fraction(value)
    denominator = value.denominator % prime
    if denominator == 0:
        raise ZeroDivisionError("a rational denominator vanishes modulo the prime")
    return value.numerator * pow(denominator, -1, prime) % prime


def _finite_model(
    coefficients: Sequence[Fraction | int], prime: int
) -> tuple[int, int, int, int, int]:
    return tuple(_fraction_mod_prime(value, prime) for value in coefficients)  # type: ignore[return-value]


def _finite_negate(
    coefficients: tuple[int, int, int, int, int],
    point: FinitePoint,
    prime: int,
) -> FinitePoint:
    if point is None:
        return None
    a1, _, a3, _, _ = coefficients
    x_coordinate, y_coordinate = point
    return x_coordinate, (-y_coordinate - a1 * x_coordinate - a3) % prime


def _finite_add(
    coefficients: tuple[int, int, int, int, int],
    left: FinitePoint,
    right: FinitePoint,
    prime: int,
) -> FinitePoint:
    a1, a2, a3, a4, _ = coefficients
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2 and right == _finite_negate(coefficients, left, prime):
        return None
    if left == right:
        denominator = (2 * y1 + a1 * x1 + a3) % prime
        if denominator == 0:
            return None
        slope = (
            (3 * x1 * x1 + 2 * a2 * x1 + a4 - a1 * y1)
            * pow(denominator, -1, prime)
        ) % prime
    else:
        if x1 == x2:
            raise ArithmeticError("invalid finite-field same-abscissa pair")
        slope = ((y2 - y1) * pow((x2 - x1) % prime, -1, prime)) % prime
    intercept = (y1 - slope * x1) % prime
    x3 = (slope * slope + a1 * slope - a2 - x1 - x2) % prime
    y3 = (-(slope + a1) * x3 - intercept - a3) % prime
    return x3, y3


def _finite_multiply(
    coefficients: tuple[int, int, int, int, int],
    point: FinitePoint,
    scalar: int,
    prime: int,
) -> FinitePoint:
    if scalar < 0:
        return _finite_multiply(
            coefficients, _finite_negate(coefficients, point, prime), -scalar, prime
        )
    result: FinitePoint = None
    addend = point
    while scalar:
        if scalar & 1:
            result = _finite_add(coefficients, result, addend, prime)
        addend = _finite_add(coefficients, addend, addend, prime)
        scalar >>= 1
    return result


def _finite_points(
    coefficients: tuple[int, int, int, int, int], prime: int
) -> tuple[FinitePoint, ...]:
    a1, a2, a3, a4, a6 = coefficients
    square_roots: dict[int, list[int]] = {}
    for value in range(prime):
        square_roots.setdefault(value * value % prime, []).append(value)
    inverse_two = pow(2, -1, prime)
    points: list[FinitePoint] = [None]
    for x_coordinate in range(prime):
        linear_y = (a1 * x_coordinate + a3) % prime
        rhs = (
            x_coordinate**3
            + a2 * x_coordinate**2
            + a4 * x_coordinate
            + a6
        ) % prime
        discriminant = (linear_y * linear_y + 4 * rhs) % prime
        for square_root in square_roots.get(discriminant, ()):
            points.append(
                (x_coordinate, ((square_root - linear_y) * inverse_two) % prime)
            )
    return tuple(points)


def _reduce_point(point: AffinePoint, prime: int) -> tuple[int, int]:
    return (
        _fraction_mod_prime(point[0], prime),
        _fraction_mod_prime(point[1], prime),
    )


def _is_good_reduction(
    coefficients: Sequence[Fraction | int], prime: int
) -> bool:
    try:
        return _fraction_mod_prime(weierstrass_discriminant(coefficients), prime) != 0
    except ZeroDivisionError:
        return False


def _cyclic_generator_and_logs(
    coefficients: tuple[int, int, int, int, int],
    points: Sequence[FinitePoint],
    prime: int,
) -> tuple[tuple[int, int], dict[FinitePoint, int]] | None:
    order = len(points)
    prime_divisors = _prime_divisors(order)
    for point in points:
        if point is None:
            continue
        if all(
            _finite_multiply(coefficients, point, order // divisor, prime) is not None
            for divisor in prime_divisors
        ):
            logs: dict[FinitePoint, int] = {}
            current: FinitePoint = None
            for scalar in range(order):
                if current in logs:
                    raise ArithmeticError("generator orbit repeated prematurely")
                logs[current] = scalar
                current = _finite_add(coefficients, current, point, prime)
            if current is not None or len(logs) != order:
                raise ArithmeticError("generator did not have the asserted order")
            return point, logs
    return None


def _scan_independence_rows(
    model: WeierstrassModel,
    rational_points: tuple[AffinePoint, ...],
    relation_prime: int,
    maximum_reduction_prime: int,
    *,
    stop_at_full_rank: bool,
) -> tuple[tuple[int, int] | None, tuple[ReductionRow, ...], int]:
    rows: list[ReductionRow] = []
    torsion_witness: tuple[int, int] | None = None
    rank = 0
    for prime in _primes_between(5, maximum_reduction_prime):
        if prime == relation_prime or not _is_good_reduction(model, prime):
            continue
        try:
            finite_model = _finite_model(model, prime)
        except ZeroDivisionError:
            continue
        finite_points = _finite_points(finite_model, prime)
        group_order = len(finite_points)
        if torsion_witness is None and group_order % relation_prime:
            torsion_witness = prime, group_order
        if group_order % relation_prime:
            if (
                stop_at_full_rank
                and rank == len(rational_points)
                and torsion_witness is not None
            ):
                break
            continue
        try:
            reduced_points = tuple(
                _reduce_point(point, prime) for point in rational_points
            )
        except ZeroDivisionError:
            continue
        cyclic_data = _cyclic_generator_and_logs(finite_model, finite_points, prime)
        if cyclic_data is None:
            continue
        generator, log_table = cyclic_data
        logs = tuple(log_table[point] for point in reduced_points)
        candidate_rank = matrix_rank_mod_prime(
            [row.logs for row in rows] + [logs], relation_prime
        )
        if candidate_rank > rank:
            rows.append(ReductionRow(prime, group_order, generator, logs))
            rank = candidate_rank
        if (
            stop_at_full_rank
            and rank == len(rational_points)
            and torsion_witness is not None
        ):
            break
    return torsion_witness, tuple(rows), rank


def build_independence_certificate(
    coefficients: Sequence[Fraction | int],
    points: Sequence[AffinePoint],
    *,
    relation_prime: int = 2,
    maximum_reduction_prime: int = 2000,
) -> IndependenceCertificate:
    """Build a deterministic exact independence certificate.

    Failure to find a certificate within the prime bound says nothing about
    dependence: the chosen points may be independent but fail this particular
    finite-reduction test (for example because of index divisible by
    ``relation_prime``).
    """

    model = _model(coefficients)
    rational_points = tuple((Fraction(x), Fraction(y)) for x, y in points)
    if not rational_points:
        raise ValueError("at least one point is required")
    if not _is_prime(relation_prime):
        raise ValueError("relation_prime must be prime")
    if any(not is_on_weierstrass_curve(model, point) for point in rational_points):
        raise ValueError("every point must lie on the supplied curve")

    torsion_witness, rows, rank = _scan_independence_rows(
        model,
        rational_points,
        relation_prime,
        maximum_reduction_prime,
        stop_at_full_rank=True,
    )
    if torsion_witness is None:
        raise ArithmeticError(
            "no good-reduction prime excluding rational relation-prime torsion was found"
        )
    if rank != len(rational_points):
        raise ArithmeticError(
            f"only finite-reduction rank {rank}/{len(rational_points)} was found "
            f"through prime {maximum_reduction_prime}"
        )
    return IndependenceCertificate(
        relation_prime=relation_prime,
        torsion_witness_prime=torsion_witness[0],
        torsion_witness_group_order=torsion_witness[1],
        rows=rows,
    )


def select_independent_subset(
    coefficients: Sequence[Fraction | int],
    points: Sequence[AffinePoint],
    *,
    relation_prime: int = 2,
    maximum_reduction_prime: int = 2000,
) -> tuple[tuple[int, ...], IndependenceCertificate]:
    """Select and certify a largest subset seen by the bounded reduction scan.

    Columns are considered in input order, so callers can place a preferred
    known-independent baseline first.  The result is maximal only for the
    collected modular rows; it need not equal the rational rank of the full
    point cloud.
    """

    model = _model(coefficients)
    rational_points = tuple((Fraction(x), Fraction(y)) for x, y in points)
    if not rational_points:
        raise ValueError("at least one point is required")
    if not _is_prime(relation_prime):
        raise ValueError("relation_prime must be prime")
    if any(not is_on_weierstrass_curve(model, point) for point in rational_points):
        raise ValueError("every point must lie on the supplied curve")
    torsion_witness, rows, rank = _scan_independence_rows(
        model,
        rational_points,
        relation_prime,
        maximum_reduction_prime,
        stop_at_full_rank=False,
    )
    if torsion_witness is None:
        raise ArithmeticError(
            "no good-reduction prime excluding rational relation-prime torsion was found"
        )
    if rank == 0:
        raise ArithmeticError(
            f"no nonzero finite-reduction rank was found through prime "
            f"{maximum_reduction_prime}"
        )

    selected: list[int] = []
    selected_rank = 0
    for column in range(len(rational_points)):
        candidate_columns = selected + [column]
        candidate_matrix = [
            [row.logs[index] for index in candidate_columns] for row in rows
        ]
        candidate_rank = matrix_rank_mod_prime(candidate_matrix, relation_prime)
        if candidate_rank > selected_rank:
            selected.append(column)
            selected_rank = candidate_rank
        if selected_rank == rank:
            break
    assert len(selected) == rank
    restricted_rows = tuple(
        ReductionRow(
            prime=row.prime,
            group_order=row.group_order,
            generator=row.generator,
            logs=tuple(row.logs[index] for index in selected),
        )
        for row in rows
    )
    certificate = IndependenceCertificate(
        relation_prime=relation_prime,
        torsion_witness_prime=torsion_witness[0],
        torsion_witness_group_order=torsion_witness[1],
        rows=restricted_rows,
    )
    verify_independence_certificate(
        model, tuple(rational_points[index] for index in selected), certificate
    )
    return tuple(selected), certificate


def verify_independence_certificate(
    coefficients: Sequence[Fraction | int],
    points: Sequence[AffinePoint],
    certificate: IndependenceCertificate,
) -> None:
    """Replay all exact claims needed for the infinite-descent proof."""

    model = _model(coefficients)
    rational_points = tuple((Fraction(x), Fraction(y)) for x, y in points)
    relation_prime = certificate.relation_prime
    if not _is_prime(relation_prime):
        raise AssertionError("relation modulus is not prime")
    if any(not is_on_weierstrass_curve(model, point) for point in rational_points):
        raise AssertionError("a certified point is not on the supplied curve")

    witness_prime = certificate.torsion_witness_prime
    assert witness_prime != relation_prime and _is_prime(witness_prime)
    assert _is_good_reduction(model, witness_prime)
    witness_model = _finite_model(model, witness_prime)
    assert len(_finite_points(witness_model, witness_prime)) == (
        certificate.torsion_witness_group_order
    )
    assert certificate.torsion_witness_group_order % relation_prime != 0

    log_rows: list[tuple[int, ...]] = []
    for row in certificate.rows:
        assert _is_prime(row.prime) and row.prime != relation_prime
        assert _is_good_reduction(model, row.prime)
        finite_model = _finite_model(model, row.prime)
        finite_points = _finite_points(finite_model, row.prime)
        assert len(finite_points) == row.group_order
        assert row.group_order % relation_prime == 0
        generator: FinitePoint = row.generator
        assert generator in finite_points
        assert _finite_multiply(
            finite_model, generator, row.group_order, row.prime
        ) is None
        for divisor in _prime_divisors(row.group_order):
            assert _finite_multiply(
                finite_model, generator, row.group_order // divisor, row.prime
            ) is not None
        assert len(row.logs) == len(rational_points)
        for point, logarithm in zip(rational_points, row.logs):
            assert 0 <= logarithm < row.group_order
            assert _finite_multiply(
                finite_model, generator, logarithm, row.prime
            ) == _reduce_point(point, row.prime)
        log_rows.append(row.logs)
    assert matrix_rank_mod_prime(log_rows, relation_prime) == len(rational_points)
