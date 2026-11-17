"""Sparse Macaulay systems and modular block-Wiedemann reconstruction.

The polynomial representation is deliberately small and dependency free:
a polynomial is a dictionary from exponent tuples to rational coefficients.
The computational driver may use python-flint for dense projected solves, but
the certificate checker only needs the elementary arithmetic helpers here.

The block-Wiedemann routine uses two block Krylov spaces.  For a nonsingular
square matrix ``A`` it forms

    K = [V, AV, ...],        P = [U, A^T U, ...]

with exactly ``n`` columns, solves the projected system

    (P^T A K) y = P^T b,

and returns ``x = K y``.  Every answer is replayed against the original
sparse matrix, so an unlucky projection is rejected rather than trusted.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
from fractions import Fraction
from typing import Iterable, Sequence


Exponent = tuple[int, ...]
RationalPolynomial = dict[Exponent, Fraction]
IntegerPolynomial = dict[Exponent, int]
SparseColumn = dict[int, int]


def add_polynomial(
    destination: RationalPolynomial,
    source: RationalPolynomial,
    scale: Fraction = Fraction(1),
) -> None:
    """Add ``scale * source`` to ``destination`` in place."""

    for monomial, coefficient in source.items():
        value = destination.get(monomial, Fraction(0))
        value += scale * coefficient
        if value:
            destination[monomial] = value
        else:
            destination.pop(monomial, None)


def multiply_polynomial(
    first: RationalPolynomial,
    second: IntegerPolynomial,
    normal_variables: int | None = None,
    normal_cutoff: int | None = None,
) -> RationalPolynomial:
    """Multiply two sparse polynomials, optionally modulo a normal power."""

    result: RationalPolynomial = {}
    for first_monomial, first_coefficient in first.items():
        for second_monomial, second_coefficient in second.items():
            monomial = tuple(
                left + right
                for left, right in zip(
                    first_monomial,
                    second_monomial,
                    strict=True,
                )
            )
            if (
                normal_cutoff is not None
                and normal_variables is not None
                and sum(monomial[:normal_variables]) >= normal_cutoff
            ):
                continue
            value = result.get(monomial, Fraction(0))
            value += first_coefficient * second_coefficient
            if value:
                result[monomial] = value
            else:
                result.pop(monomial, None)
    return result


def monomials(
    variable_count: int,
    maximum_degree: int,
) -> tuple[Exponent, ...]:
    """Return all exponent tuples of total degree at most the bound."""

    if variable_count == 0:
        return ((),)
    output: list[Exponent] = []

    def visit(prefix: tuple[int, ...], remaining: int) -> None:
        if len(prefix) == variable_count - 1:
            for last in range(remaining + 1):
                output.append(prefix + (last,))
            return
        for exponent in range(remaining + 1):
            visit(prefix + (exponent,), remaining - exponent)

    visit((), maximum_degree)
    return tuple(sorted(output, key=lambda item: (sum(item), item)))


def primitive_integer_polynomial(
    polynomial: RationalPolynomial,
) -> tuple[IntegerPolynomial, Fraction]:
    """Clear denominators and content, returning a primitive polynomial."""

    if not polynomial:
        return {}, 1
    denominator = 1
    for coefficient in polynomial.values():
        denominator = math.lcm(denominator, coefficient.denominator)
    integers = {
        monomial: int(coefficient * denominator)
        for monomial, coefficient in polynomial.items()
        if coefficient
    }
    content = 0
    for coefficient in integers.values():
        content = math.gcd(content, abs(coefficient))
    integers = {
        monomial: coefficient // content
        for monomial, coefficient in integers.items()
    }
    leading = integers[min(integers, key=lambda item: (sum(item), item))]
    if leading < 0:
        integers = {
            monomial: -coefficient
            for monomial, coefficient in integers.items()
        }
        content = -content
    return integers, Fraction(content, denominator)


def serialize_polynomial(
    polynomial: IntegerPolynomial | RationalPolynomial,
) -> list[list[object]]:
    """Return a canonical JSON-compatible sparse polynomial."""

    output: list[list[object]] = []
    for monomial in sorted(polynomial, key=lambda item: (sum(item), item)):
        coefficient = polynomial[monomial]
        if isinstance(coefficient, Fraction):
            encoded: object = (
                coefficient.numerator
                if coefficient.denominator == 1
                else [coefficient.numerator, coefficient.denominator]
            )
        else:
            encoded = coefficient
        output.append([list(monomial), encoded])
    return output


def deserialize_polynomial(
    data: Sequence[Sequence[object]],
) -> RationalPolynomial:
    """Decode the canonical sparse-polynomial JSON representation."""

    result: RationalPolynomial = {}
    for raw_monomial, raw_coefficient in data:
        monomial = tuple(int(value) for value in raw_monomial)
        if isinstance(raw_coefficient, list):
            coefficient = Fraction(
                int(raw_coefficient[0]),
                int(raw_coefficient[1]),
            )
        else:
            coefficient = Fraction(int(raw_coefficient))
        if coefficient:
            result[monomial] = coefficient
    return result


def canonical_digest(value: object) -> str:
    """Hash a JSON object with stable separators and key ordering."""

    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def sparse_matvec(
    rows: Sequence[dict[int, int]],
    vector: Sequence[int],
    prime: int,
) -> list[int]:
    """Multiply a sparse row matrix by one vector modulo ``prime``."""

    return [
        sum(coefficient * vector[column] for column, coefficient in row.items())
        % prime
        for row in rows
    ]


def sparse_block_matvec(
    rows: Sequence[dict[int, int]],
    block: Sequence[Sequence[int]],
    prime: int,
) -> list[list[int]]:
    """Multiply a sparse row matrix by a narrow dense block."""

    if not block:
        return [[] for _row in rows]
    width = len(block[0])
    output: list[list[int]] = []
    for row in rows:
        values = [0] * width
        for column, coefficient in row.items():
            source = block[column]
            for index in range(width):
                values[index] += coefficient * source[index]
        output.append([value % prime for value in values])
    return output


def transpose_rows(
    rows: Sequence[dict[int, int]],
    column_count: int,
) -> list[dict[int, int]]:
    """Transpose sparse rows."""

    output = [{} for _index in range(column_count)]
    for row_index, row in enumerate(rows):
        for column_index, coefficient in row.items():
            output[column_index][row_index] = coefficient
    return output


def _random_block(
    dimension: int,
    width: int,
    prime: int,
    rng: random.Random,
) -> list[list[int]]:
    return [
        [rng.randrange(prime) for _index in range(width)]
        for _row in range(dimension)
    ]


def _krylov_basis(
    rows: Sequence[dict[int, int]],
    dimension: int,
    block_size: int,
    prime: int,
    rng: random.Random,
) -> list[list[int]]:
    """Build exactly ``dimension`` block-Krylov columns."""

    first_width = min(block_size, dimension)
    block = _random_block(dimension, first_width, prime, rng)
    columns: list[list[int]] = []
    while len(columns) < dimension:
        take = min(len(block[0]), dimension - len(columns))
        columns.extend(
            [
                [block[row][column] for row in range(dimension)]
                for column in range(take)
            ]
        )
        if len(columns) == dimension:
            break
        block = sparse_block_matvec(rows, block, prime)
    return columns


def block_wiedemann_solve(
    rows: Sequence[dict[int, int]],
    right_sides: Sequence[Sequence[int]],
    prime: int,
    block_size: int,
    seed: int,
) -> tuple[list[list[int]], dict[str, object]]:
    """Solve a nonsingular sparse square system by block Wiedemann.

    The returned solutions are row-major: one length-``n`` list per right
    side.  A failed random projection raises ``ZeroDivisionError`` so callers
    can retry with a different seed.
    """

    from flint import nmod_mat

    dimension = len(rows)
    if any(
        column < 0 or column >= dimension
        for row in rows
        for column in row
    ):
        raise ValueError("block-Wiedemann input must be square")
    rng = random.Random(seed)
    right_rows = transpose_rows(rows, dimension)
    krylov = _krylov_basis(
        rows,
        dimension,
        block_size,
        prime,
        rng,
    )
    test_krylov = _krylov_basis(
        right_rows,
        dimension,
        block_size,
        prime,
        rng,
    )
    # Matrices are stored as lists of columns above.  The projected matrix is
    # P^T A K.  Apply sparse A in narrow blocks; only the final dense products
    # and solve are delegated to FLINT.
    applied_rows = [[] for _index in range(dimension)]
    for start in range(0, dimension, block_size):
        narrow_columns = krylov[start : start + block_size]
        narrow_rows = [
            [column[row] for column in narrow_columns]
            for row in range(dimension)
        ]
        applied = sparse_block_matvec(rows, narrow_rows, prime)
        for row in range(dimension):
            applied_rows[row].extend(applied[row])
    test_transpose = nmod_mat(test_krylov, prime)
    projected_matrix = test_transpose * nmod_mat(applied_rows, prime)
    right_rows = [
        [right_side[row] for right_side in right_sides]
        for row in range(dimension)
    ]
    right_matrix = test_transpose * nmod_mat(right_rows, prime)
    try:
        coordinates = projected_matrix.solve(right_matrix)
    except ZeroDivisionError:
        raise
    krylov_rows = [
        [column[row] for column in krylov]
        for row in range(dimension)
    ]
    solution_matrix = nmod_mat(krylov_rows, prime) * coordinates
    solutions: list[list[int]] = []
    for right_index in range(len(right_sides)):
        solution = [
            int(solution_matrix[row, right_index])
            for row in range(dimension)
        ]
        if sparse_matvec(rows, solution, prime) != [
            value % prime for value in right_sides[right_index]
        ]:
            raise ZeroDivisionError("unlucky block-Wiedemann projection")
        solutions.append(solution)
    transcript = {
        "algorithm": "two-sided-block-krylov-wiedemann",
        "block_size": block_size,
        "dimension": dimension,
        "seed": seed,
    }
    return solutions, transcript


def crt_pair(
    current: int,
    modulus: int,
    residue: int,
    prime: int,
) -> int:
    """Combine one residue with a current Chinese remainder."""

    correction = (
        (residue - current) * pow(modulus % prime, -1, prime)
    ) % prime
    return current + modulus * correction


def rational_reconstruction(
    residue: int,
    modulus: int,
) -> Fraction | None:
    """Reconstruct a rational using the symmetric square-root bound."""

    residue %= modulus
    bound = math.isqrt(modulus // 2)
    old_r, r = modulus, residue
    old_t, t = 0, 1
    while abs(r) > bound:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_t, t = t, old_t - quotient * t
    if t == 0 or abs(t) > bound or math.gcd(r, t) != 1:
        return None
    if (residue * t - r) % modulus:
        return None
    if t < 0:
        r, t = -r, -t
    return Fraction(r, t)


def reconstruct_vector(
    residues: Sequence[int],
    modulus: int,
) -> list[Fraction] | None:
    """Rationally reconstruct a whole vector, or return ``None``."""

    output: list[Fraction] = []
    for residue in residues:
        reconstructed = rational_reconstruction(residue, modulus)
        if reconstructed is None:
            return None
        output.append(reconstructed)
    return output


def evaluate_functional(
    functional: RationalPolynomial,
    polynomial: RationalPolynomial,
) -> Fraction:
    """Pair a finite-support coefficient functional with a polynomial."""

    return sum(
        coefficient * polynomial.get(monomial, Fraction(0))
        for monomial, coefficient in functional.items()
    )


def annihilator_multipliers(
    support: Iterable[Exponent],
    generator: IntegerPolynomial,
) -> tuple[Exponent, ...]:
    """List every monomial multiplier that can meet functional support."""

    output: set[Exponent] = set()
    for support_monomial in support:
        for generator_monomial in generator:
            if all(
                left >= right
                for left, right in zip(
                    support_monomial,
                    generator_monomial,
                    strict=True,
                )
            ):
                output.add(
                    tuple(
                        left - right
                        for left, right in zip(
                            support_monomial,
                            generator_monomial,
                            strict=True,
                        )
                    )
                )
    return tuple(sorted(output, key=lambda item: (sum(item), item)))
