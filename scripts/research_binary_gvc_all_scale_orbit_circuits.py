#!/usr/bin/env python3
"""Exact all-scale-prefix elimination of small projected Hall shells.

The projected two-colour return configuration has states

    R_i -> (1, 0, i),    B_i -> (0, 1, i),    0 <= i <= span.

This script enumerates every fixed-cardinality subset of every fibre with
bounded positive colour counts.  For each sign pattern (modulo one common sign), it
forms the factorially weighted scaled rows

    F_N = sum_j phase_j(N) W_j(N) M_j^N,

where ``W_j(N)`` is the product of the two exact multinomial coefficients.
It saturates the complete one-colour Taylor derivative ideal by every
coefficient channel.  The scale-two test is run with both natural phase laws:

``fixed``
    ``phase_j(N) = phase_j(1)``;
``character_power``
    ``phase_j(N) = phase_j(1)**N``.

For three states the bounded Singular census is a regression for a stronger
elementary theorem.  Over any characteristic-zero field, allow arbitrary
nonzero scalar coefficients at scales one and two.  If three monomial orbit
functions obey both rows, absorb the scale-one coefficients and write

    A_1 + A_2 + A_3 = 0,
    lambda_1 A_1^2 + lambda_2 A_2^2 + lambda_3 A_3^2 = 0.

After setting ``z=A_2/A_1``, the second equation is a nonzero quadratic in
``z`` with constant coefficients (its linear coefficient is
``2*lambda_3``).  Hence ``z`` is constant.  All three orbit functions are
proportional, and their Taylor degrees agree.  Therefore a three-state fibre
containing two different marked-side levels has empty two-scale coefficient
torus.  A remaining equal-level triple is a zero-transfer product block;
each pair splits into one-colour relations and is outside the genuinely
mixed primitive-transfer obstruction.

The ``--structural-certificate`` mode checks five stronger reductions.  A
confluent-Vandermonde determinant certifies the leading term of the
Wronskian which separates every fixed finite affine-ray family
``h_j(t)*f_j(t)**N`` for all sufficiently large ``N``, with arbitrary scalar
rows.  A Singular saturation certifies that the smooth-conic branch of the
arbitrary-coefficient four-state scale-1/2/3 moment ideal is empty.  The
remaining rank-two conic branch is classified in the canonical note by an
elementary four-point Vandermonde argument and consists only of two pair
blocks.  Exact rational row reduction also replays the common-base theorem:
after factoring an invertible base `f**N`, the infinite Cartesian ideal
depends only on the finite span of the scale coefficient vectors.
Finally, an exact Normaliz Graver census audits every primitive relation in
a bounded integer-affine factorial universe against the unbounded
boundary-transfer presentation.  A rational-function regression also
checks the exact residue-wise elimination of eventually periodic additive
factorial coefficients.

For four states, ``--require-factorial-pairing`` implements the all-scale
filter supplied by factorial trace independence.  A fixed-sign all-scale row
first splits by equality of its factorial rays; power-sum partial fractions
then force opposite-sign equal-monomial pairs.  The script tests those pair
orbit ideals in either colour and under the two independent derivations.  At
span two it additionally compares every independent survivor, in both ideal
directions, with the product of the quadratic Veronese ideals.

The script is deliberately a projected, bounded experiment.  It does not
enumerate actual affine Hall carry alphabets and is not the proof of
unrestricted GVC(2).  Its exact output distinguishes a torus survivor from a
GVC counterexample; the unrestricted proof is the separate Hall-envelope
theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import shutil
import subprocess
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable

from research_binary_gvc_prime_power_tomography import (
    State,
    normaliz_version,
)
from research_binary_gvc_translation_observability import (
    Polynomial,
    derivative_closure,
    logical_hash,
    polynomial_text,
    singular_torus_test,
    singular_version,
    state_exponent,
    state_weight,
)


PHASE_LAWS = ("fixed", "character_power")


def compositions(total: int, length: int) -> Iterable[tuple[int, ...]]:
    """Yield the weak compositions of ``total`` into ``length`` parts."""

    if length == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, length - 1):
            yield (first,) + tail


def side_level(side: tuple[int, ...]) -> int:
    return sum(index * multiplicity for index, multiplicity in enumerate(side))


def scaled_state(state: State, scale: int) -> State:
    return State(
        tuple(scale * value for value in state.operator),
        tuple(scale * value for value in state.polynomial),
    )


def scaled_row(
    states: tuple[State, ...],
    signs: tuple[int, ...],
    scale: int,
    phase_law: str,
) -> Polynomial:
    """Return the exact scaled multinomial row as a sparse polynomial."""

    answer: Polynomial = {}
    for state, sign in zip(states, signs, strict=True):
        scaled = scaled_state(state, scale)
        phase = sign if phase_law == "fixed" else sign**scale
        exponent = state_exponent(scaled)
        answer[exponent] = answer.get(exponent, 0) + phase * state_weight(scaled)
    if len(answer) != len(states):
        raise AssertionError("the states must stay distinct after scaling")
    return answer


def candidate_key(
    states: tuple[State, ...],
    signs: tuple[int, ...],
) -> str:
    payload = [
        {
            "operator": list(state.operator),
            "polynomial": list(state.polynomial),
            "sign": sign,
        }
        for state, sign in zip(states, signs, strict=True)
    ]
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def transform_state(state: State, *, swap: bool, reverse: bool) -> State:
    operator = state.operator
    polynomial = state.polynomial
    if swap:
        operator, polynomial = polynomial, operator
    if reverse:
        operator = tuple(reversed(operator))
        polynomial = tuple(reversed(polynomial))
    return State(operator, polynomial)


def symmetry_key(
    states: tuple[State, ...],
    signs: tuple[int, ...],
) -> tuple[int, ...]:
    """Canonicalize side exchange, common reversal, and common row sign."""

    images = []
    for swap in (False, True):
        for reverse in (False, True):
            terms = sorted(
                (
                    transform_state(state, swap=swap, reverse=reverse),
                    sign,
                )
                for state, sign in zip(states, signs, strict=True)
            )
            common_sign = terms[0][1]
            normalized = []
            for state, sign in terms:
                normalized.extend(state.operator)
                normalized.extend(state.polynomial)
                normalized.append(sign * common_sign)
            images.append(tuple(normalized))
    return min(images)


def affine_solutions_mod_two(
    matrix: tuple[tuple[int, ...], ...],
    target: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    """Solve a small affine system over F_2, returning every solution."""

    if len(matrix) != len(target):
        raise ValueError("matrix and target height differ")
    variable_count = len(matrix[0]) if matrix else 0
    rows = [
        [value % 2 for value in row] + [right % 2]
        for row, right in zip(matrix, target, strict=True)
    ]
    pivot_columns = []
    pivot_row = 0
    for column in range(variable_count):
        source = next(
            (index for index in range(pivot_row, len(rows)) if rows[index][column]),
            None,
        )
        if source is None:
            continue
        rows[pivot_row], rows[source] = rows[source], rows[pivot_row]
        for index in range(len(rows)):
            if index != pivot_row and rows[index][column]:
                rows[index] = [
                    left ^ right
                    for left, right in zip(rows[index], rows[pivot_row], strict=True)
                ]
        pivot_columns.append(column)
        pivot_row += 1

    if any(not any(row[:variable_count]) and row[-1] for row in rows):
        return ()

    free_columns = [
        column for column in range(variable_count) if column not in pivot_columns
    ]
    solutions = []
    for free_values in itertools.product((0, 1), repeat=len(free_columns)):
        solution = [0] * variable_count
        for column, value in zip(free_columns, free_values, strict=True):
            solution[column] = value
        for index, column in reversed(tuple(enumerate(pivot_columns))):
            solution[column] = rows[index][-1] ^ sum(
                rows[index][free] * solution[free] for free in free_columns
            ) % 2
        solutions.append(tuple(solution))
    return tuple(solutions)


def character_feasible(
    states: tuple[State, ...],
    signs: tuple[int, ...],
    modulus: int,
) -> bool:
    """Whether one channel label realizes the prescribed real phases."""

    matrix = tuple(state_exponent(state) for state in states)
    if modulus == 2:
        target = tuple(0 if sign == 1 else 1 for sign in signs)
        return bool(affine_solutions_mod_two(matrix, target))
    if modulus != 4:
        raise ValueError("only C2 and real C4 feasibility are implemented")

    target = tuple(0 if sign == 1 else 2 for sign in signs)
    # Write x=x_0+2y.  First solve A*x_0=0 mod 2, then lift the
    # residual by solving one further affine system over F_2.
    zero = (0,) * len(states)
    for low_bits in affine_solutions_mod_two(matrix, zero):
        residual = tuple(
            ((right - sum(value * bit for value, bit in zip(row, low_bits))) % 4)
            // 2
            for row, right in zip(matrix, target, strict=True)
        )
        if affine_solutions_mod_two(matrix, residual):
            return True
    return False


def fibre_states(
    span: int,
    minimum_state_count: int,
    maximum_operator_count: int,
    maximum_polynomial_count: int,
) -> dict[tuple[int, int, int], tuple[State, ...]]:
    width = span + 1
    fibres: dict[tuple[int, int, int], list[State]] = defaultdict(list)
    for operator_count in range(1, maximum_operator_count + 1):
        operators = tuple(compositions(operator_count, width))
        for polynomial_count in range(1, maximum_polynomial_count + 1):
            polynomials = tuple(compositions(polynomial_count, width))
            for operator in operators:
                operator_level = side_level(operator)
                for polynomial in polynomials:
                    state = State(operator, polynomial)
                    key = (
                        operator_count,
                        polynomial_count,
                        operator_level + side_level(polynomial),
                    )
                    fibres[key].append(state)
    return {
        key: tuple(sorted(states))
        for key, states in sorted(fibres.items())
        if len(states) >= minimum_state_count
    }


def state_record(state: State) -> dict[str, Any]:
    return {
        "operator": list(state.operator),
        "polynomial": list(state.polynomial),
        "operator_level": side_level(state.operator),
        "polynomial_level": side_level(state.polynomial),
        "scale_1_weight": state_weight(state),
        "scale_2_weight": state_weight(scaled_state(state, 2)),
    }


def factorial_vector(state: State) -> tuple[int, ...]:
    """Canonical denominator slopes for the common-fibre factorial ray."""

    return tuple(
        sorted(
            multiplicity
            for multiplicity in state.operator + state.polynomial
            if multiplicity
        )
    )


def four_state_factorial_pairings(
    states: tuple[State, ...],
    signs: tuple[int, ...],
) -> tuple[tuple[tuple[int, int], tuple[int, int]], ...]:
    """Opposite-sign pairings with equal all-scale scalar factorial rays."""

    if len(states) != 4:
        return ()
    pairings = (
        ((0, 1), (2, 3)),
        ((0, 2), (1, 3)),
        ((0, 3), (1, 2)),
    )
    return tuple(
        pairing
        for pairing in pairings
        if all(
            signs[left] == -signs[right]
            and factorial_vector(states[left]) == factorial_vector(states[right])
            for left, right in pairing
        )
    )


def quadratic_veronese_pair_tests(
    requests: list[
        tuple[
            str,
            int,
            tuple[State, ...],
            tuple[tuple[int, int], tuple[int, int]],
        ]
    ],
    singular: str,
) -> dict[tuple[str, int], bool]:
    """Compare independent pair ideals with the span-two Veronese product."""

    if not requests:
        return {}
    variable_names = ("R0", "R1", "R2", "B0", "B1", "B2")
    lines = []
    for index, (digest, pairing_index, states, pairing) in enumerate(requests):
        generators = []
        for left, right in pairing:
            polynomial = {
                state_exponent(states[left]): 1,
                state_exponent(states[right]): -1,
            }
            generators.extend(derivative_closure(polynomial, "independent", 3))
        generator_texts = [
            polynomial_text(generator, variable_names)
            for generator in generators
        ]
        lines.extend(
            (
                f"ring q{index}=0,({','.join(variable_names)},T),dp;",
                f"ideal I={','.join(generator_texts)},"
                "1-T*R0*R1*R2*B0*B1*B2;",
                "ideal GI=std(I);",
                "ideal J=R1^2-4*R0*R2,B1^2-4*B0*B2,"
                "1-T*R0*R1*R2*B0*B1*B2;",
                "ideal GJ=std(J);",
                "int forward=1;",
                "if (reduce(R1^2-4*R0*R2,GI)!=0) {forward=0;}",
                "if (reduce(B1^2-4*B0*B2,GI)!=0) {forward=0;}",
                "int reverse=1;",
            )
        )
        for generator_text in generator_texts:
            lines.append(
                f"if (reduce({generator_text},GJ)!=0) {{reverse=0;}}"
            )
        marker = f"{digest}|{pairing_index}"
        lines.extend(
            (
                f'print("BEGIN|{marker}");',
                "print(forward);",
                "print(reverse);",
                f'print("END|{marker}");',
            )
        )

    completed = subprocess.run(
        [singular, "-q"],
        input="\n".join(lines) + "\n",
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode:
        raise RuntimeError(
            "Singular Veronese comparison failed:\n"
            + completed.stdout
            + completed.stderr
        )
    results = {}
    output = completed.stdout.splitlines()
    cursor = 0
    while cursor < len(output):
        line = output[cursor].strip()
        cursor += 1
        if not line.startswith("BEGIN|"):
            continue
        marker = line.removeprefix("BEGIN|")
        digest, pairing_text = marker.split("|", 1)
        values = []
        while cursor < len(output) and not output[cursor].startswith("END|"):
            stripped = output[cursor].strip()
            cursor += 1
            if stripped:
                values.append(int(stripped))
        if cursor >= len(output) or output[cursor].strip() != f"END|{marker}":
            raise RuntimeError(f"truncated Veronese block for {marker}")
        cursor += 1
        if len(values) != 2 or any(value not in (0, 1) for value in values):
            raise RuntimeError(f"unexpected Veronese block {marker}: {values}")
        results[(digest, int(pairing_text))] = values == [1, 1]
    if len(results) != len(requests):
        raise RuntimeError(
            f"parsed {len(results)} Veronese blocks, expected {len(requests)}"
        )
    return results


def integer_partitions(
    total: int,
    maximum_part: int | None = None,
) -> Iterable[tuple[int, ...]]:
    """Yield nonincreasing positive partitions of ``total``."""

    if total == 0:
        yield ()
        return
    maximum = min(total, maximum_part or total)
    for first in range(maximum, 0, -1):
        for tail in integer_partitions(total - first, first):
            yield (first,) + tail


def confluent_wronskian_regression(maximum_rank: int) -> dict[str, Any]:
    """Check the exact confluent-Vandermonde leading symbol through a rank."""

    try:
        import sympy
    except ImportError as error:
        raise RuntimeError(
            "the structural certificate requires SymPy; run it with "
            ".venv/bin/python"
        ) from error

    scale = sympy.symbols("N")
    records = []
    for rank in range(1, maximum_rank + 1):
        for block_sizes in integer_partitions(rank):
            logarithmic_slopes = tuple(range(2, 2 + len(block_sizes)))
            columns = []
            for slope, block_size in zip(
                logarithmic_slopes,
                block_sizes,
                strict=True,
            ):
                for correction_order in range(block_size):
                    columns.append(
                        tuple(
                            (
                                sympy.binomial(row, correction_order)
                                * (scale * slope)
                                ** (row - correction_order)
                            )
                            if row >= correction_order
                            else 0
                            for row in range(rank)
                        )
                    )
            matrix = sympy.Matrix(
                rank,
                rank,
                lambda row, column: columns[column][row],
            )
            determinant = sympy.factor(matrix.det())
            cross_degree = sum(
                block_sizes[left] * block_sizes[right]
                for left in range(len(block_sizes))
                for right in range(left + 1, len(block_sizes))
            )
            expected = scale**cross_degree
            for left in range(len(block_sizes)):
                for right in range(left + 1, len(block_sizes)):
                    expected *= (
                        logarithmic_slopes[right]
                        - logarithmic_slopes[left]
                    ) ** (block_sizes[left] * block_sizes[right])
            quotient = sympy.factor(determinant / expected)
            if quotient != 1:
                raise AssertionError(
                    "confluent Vandermonde mismatch for "
                    f"{block_sizes}: {quotient}"
                )
            records.append(
                {
                    "rank": rank,
                    "block_sizes": list(block_sizes),
                    "cross_degree": cross_degree,
                    "exact_quotient": 1,
                }
            )
    return {
        "maximum_rank": maximum_rank,
        "partitions_checked": len(records),
        "identity": (
            "det(binomial(k,a)*(N*u_C)^(k-a)) = "
            "N^sum_(C<D)(d_C*d_D) * "
            "product_(C<D)(u_D-u_C)^(d_C*d_D)"
        ),
        "records": records,
    }


def rational_row_basis(
    rows: tuple[tuple[int, ...], ...],
) -> tuple[int, tuple[int, ...]]:
    """Return the exact rational row rank and lexicographic basis indices."""

    if not rows:
        return 0, ()
    width = len(rows[0])
    echelon: list[list[Fraction]] = []
    pivots: list[int] = []
    basis_indices = []
    for row_index, row in enumerate(rows):
        if len(row) != width:
            raise ValueError("coefficient rows have different widths")
        reduced = [Fraction(value) for value in row]
        for pivot, basis_row in zip(pivots, echelon, strict=True):
            factor = reduced[pivot]
            if factor:
                reduced = [
                    left - factor * right
                    for left, right in zip(reduced, basis_row, strict=True)
                ]
        pivot = next(
            (index for index, value in enumerate(reduced) if value),
            None,
        )
        if pivot is None:
            continue
        factor = reduced[pivot]
        reduced = [value / factor for value in reduced]
        for index, basis_row in enumerate(echelon):
            factor = basis_row[pivot]
            if factor:
                echelon[index] = [
                    left - factor * right
                    for left, right in zip(
                        basis_row,
                        reduced,
                        strict=True,
                    )
                ]
        insertion = next(
            (
                index
                for index, old_pivot in enumerate(pivots)
                if pivot < old_pivot
            ),
            len(pivots),
        )
        pivots.insert(insertion, pivot)
        echelon.insert(insertion, reduced)
        basis_indices.append(row_index)
    return len(echelon), tuple(basis_indices)


def common_base_span_regression() -> dict[str, Any]:
    """Replay the finite coefficient-span alternatives of Theorem 9.1."""

    import sympy

    vandermonde_rows = tuple(
        tuple(scale**power for power in range(4))
        for scale in range(1, 9)
    )
    constant_circuit_rows = tuple(
        tuple((scale + 1) * value for value in (1, -4, -1))
        for scale in range(8)
    )
    full_rank, full_basis = rational_row_basis(vandermonde_rows)
    circuit_rank, circuit_basis = rational_row_basis(constant_circuit_rows)
    if full_rank != 4 or circuit_rank != 1:
        raise AssertionError(
            "unexpected common-base coefficient ranks: "
            f"{full_rank}, {circuit_rank}"
        )
    scale = sympy.symbols("k")
    circuit = (1, -2, 3)
    rational_classes = tuple(
        (
            value * (scale + 1),
            value * (2 * scale - 1) / (scale + 2),
        )
        for value in circuit
    )
    projective_minors = tuple(
        sympy.factor(
            circuit[left] * rational_classes[right][class_index]
            - circuit[right] * rational_classes[left][class_index]
        )
        for left in range(len(circuit))
        for right in range(left + 1, len(circuit))
        for class_index in range(2)
    )
    if any(projective_minors):
        raise AssertionError("periodic additive projective criterion failed")
    perturbed = sympy.factor(
        circuit[0]
        * (rational_classes[2][1] + 1 / (scale + 3))
        - circuit[2] * rational_classes[0][1]
    )
    if perturbed == 0:
        raise AssertionError("periodic additive perturbation was not detected")
    elementary_additive_identity = sympy.factor(
        scale + (scale + 1) - (2 * scale + 1)
    )
    if elementary_additive_identity != 0:
        raise AssertionError("affine factorial additive identity failed")
    return {
        "full_span_example": {
            "row_law": "(1,N,N^2,N^3), 1<=N<=8",
            "rank": full_rank,
            "basis_row_indices": list(full_basis),
            "conclusion": "all four Laurent correction monomials enter the ideal",
        },
        "projectively_constant_circuit_example": {
            "row_law": "(N+1)*(1,-4,-1), 0<=N<8",
            "rank": circuit_rank,
            "basis_row_indices": list(circuit_basis),
            "conclusion": (
                "all scales impose only one correction row; this is the "
                "quadratic discriminant calibration"
            ),
        },
        "periodic_additive_factorial_example": {
            "rational_slope_classes": [
                "(k+1)*u",
                "((2*k-1)/(k+2))*u",
            ],
            "circuit_vector": list(circuit),
            "projective_minors": [0] * len(projective_minors),
            "detected_perturbation": str(perturbed),
            "zero_vector_additive_identity": (
                "k + (k+1) - (2*k+1) = 0"
            ),
            "conclusion": (
                "after residue restriction and factorial-shift "
                "normalization, projective constancy is a finite list of "
                "rational-function minors"
            ),
        },
    }


def arbitrary_four_state_singular_certificate(
    singular: str,
) -> dict[str, Any]:
    """Certify the empty smooth-conic branch and the pair-block converse."""

    source = r"""
ring r=0,(l1,l2,l3,l4,m1,m2,m3,m4,a,b,c,T),dp;
poly detQ=l1*l2*l3+l4*(l1*l2+l1*l3+l2*l3);
ideal I=
 (l1+l4)*a-(m1-m4),
 (l2+l4)*b-(m2-m4),
 (l3+l4)*c-(m3-m4),
 (l1+l4)*b+2*l4*a+3*m4,
 (l1+l4)*c+2*l4*a+3*m4,
 (l2+l4)*a+2*l4*b+3*m4,
 (l2+l4)*c+2*l4*b+3*m4,
 (l3+l4)*a+2*l4*c+3*m4,
 (l3+l4)*b+2*l4*c+3*m4,
 l4*(a+b+c)+3*m4,
 1-T*l1*l2*l3*l4*m1*m2*m3*m4*detQ;
ideal G=std(I);
int smoothUnit=0;
if (reduce(1,G)==0) {smoothUnit=1;}
print("BEGIN_SMOOTH");
print(smoothUnit);
print("END_SMOOTH");

ring s=0,(l1,l2,l3,l4,m1,m2,m3,m4,x,y,z),dp;
poly Q=l1*x^2+l2*y^2+l3*z^2+l4*(x+y+z)^2;
poly C=m1*x^3+m2*y^3+m3*z^3-m4*(x+y+z)^3;
ideal J=l2+l1,l4+l3,m2-m1,m4-m3,y+x;
ideal H=std(J);
int pairConverse=1;
if (reduce(Q,H)!=0) {pairConverse=0;}
if (reduce(C,H)!=0) {pairConverse=0;}
print("BEGIN_PAIR");
print(pairConverse);
print("END_PAIR");
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=source,
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode:
        raise RuntimeError(
            "Singular four-state certificate failed:\n"
            + completed.stdout
            + completed.stderr
        )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]

    def marked_value(begin: str, end: str) -> int:
        try:
            start = lines.index(begin)
            stop = lines.index(end, start + 1)
        except ValueError as error:
            raise RuntimeError(
                f"missing Singular marker {begin}/{end}: {lines}"
            ) from error
        values = lines[start + 1 : stop]
        if len(values) != 1 or values[0] not in ("0", "1"):
            raise RuntimeError(f"unexpected Singular block {begin}: {values}")
        return int(values[0])

    smooth_unit = marked_value("BEGIN_SMOOTH", "END_SMOOTH")
    pair_converse = marked_value("BEGIN_PAIR", "END_PAIR")
    if smooth_unit != 1 or pair_converse != 1:
        raise AssertionError(
            "four-state certificate failed: "
            f"smooth={smooth_unit}, pair={pair_converse}"
        )
    return {
        "smooth_conic_divisor_torus_saturation_is_unit": True,
        "saturation_factors": (
            "all lambda_i, all mu_i, and the restricted-conic determinant"
        ),
        "representative_rank_two_pair_block_annihilates_Q_and_C": True,
        "singular_rank_two_classification": (
            "hand proof: the four cubic line functions "
            "(1+t_i*s)^3 form a Vandermonde family; nonzero weights force "
            "a 2+2 partition, opposite reciprocal lambdas, and equal mus"
        ),
    }


def affine_factorial_graver_certificate(
    *,
    maximum_slope: int,
    maximum_offset: int,
    normaliz: str,
) -> dict[str, Any]:
    """Audit primitive exact-divisor relations against boundary transfers."""

    import sympy

    from research_binary_gvc_nonfree_factorization import (
        factorial_graver_basis,
    )
    from verify_factorial_trace_independence import (
        affine_factorial_value,
        affine_slope_vector,
        affine_successor_divisor,
        boundary_transfer_constant,
        boundary_transfer_decomposition,
    )

    atoms = tuple(
        (slope, offset)
        for slope in range(1, maximum_slope + 1)
        for offset in range(-maximum_offset, maximum_offset + 1)
    )
    roots = tuple(
        sorted(
            {
                Fraction(offset + shift, slope)
                for slope, offset in atoms
                for shift in range(1, slope + 1)
            }
        )
    )
    matrix = sympy.Matrix(
        [
            [
                sum(
                    Fraction(offset + shift, slope) == root
                    for shift in range(1, slope + 1)
                )
                for slope, offset in atoms
            ]
            for root in roots
        ]
    )
    graver = factorial_graver_basis(matrix, normaliz)
    records = []
    degree_distribution: Counter[int] = Counter()
    for left, right in graver:
        profile = {
            atom: left[index] - right[index]
            for index, atom in enumerate(atoms)
            if left[index] != right[index]
        }
        if affine_successor_divisor(profile):
            raise AssertionError("a reported Graver move has nonzero divisor")
        if affine_slope_vector(profile):
            raise AssertionError("an exact-divisor Graver move changes slopes")
        transfers = boundary_transfer_decomposition(profile)
        constant = boundary_transfer_constant(transfers)
        if (
            affine_factorial_value(profile, maximum_offset + 5) != constant
            or affine_factorial_value(profile, maximum_offset + 6) != constant
        ):
            raise AssertionError("Graver boundary constant failed exact replay")
        left_degree = sum(left)
        right_degree = sum(right)
        if left_degree != right_degree:
            raise AssertionError("an exact-divisor relation changed factor count")
        degree_distribution[left_degree] += 1
        records.append(
            {
                "left": [
                    {
                        "slope": atoms[index][0],
                        "offset": atoms[index][1],
                        "multiplicity": multiplicity,
                    }
                    for index, multiplicity in enumerate(left)
                    if multiplicity
                ],
                "right": [
                    {
                        "slope": atoms[index][0],
                        "offset": atoms[index][1],
                        "multiplicity": multiplicity,
                    }
                    for index, multiplicity in enumerate(right)
                    if multiplicity
                ],
                "factorization_degree": left_degree,
                "boundary_transfer_count": sum(
                    abs(multiplicity)
                    for multiplicity, _factor, _reference in transfers
                ),
                "constant_ratio": [
                    constant.numerator,
                    constant.denominator,
                ],
                "boundary_transfers": [
                    {
                        "multiplicity": multiplicity,
                        "factor": list(factor),
                        "reference": list(reference),
                        "rational_boundary": [
                            Fraction(factor[1], factor[0]).numerator,
                            Fraction(factor[1], factor[0]).denominator,
                        ],
                    }
                    for multiplicity, factor, reference in transfers
                ],
            }
        )

    rank = int(matrix.rank())
    if maximum_slope == 5 and maximum_offset == 3:
        expected_distribution = {2: 16, 3: 25, 4: 38, 5: 24, 6: 3}
        if (
            len(atoms) != 35
            or len(roots) != 31
            or rank != 26
            or len(graver) != 106
            or dict(sorted(degree_distribution.items()))
            != expected_distribution
        ):
            raise AssertionError(
                "the default affine-factorial Graver census changed"
            )

    return {
        "atom_universe": {
            "maximum_slope": maximum_slope,
            "offset_interval": [-maximum_offset, maximum_offset],
            "atoms": len(atoms),
            "exact_rational_roots": len(roots),
        },
        "incidence_matrix": {
            "rows": matrix.rows,
            "columns": matrix.cols,
            "rank": rank,
            "integer_kernel_rank": matrix.cols - rank,
        },
        "graver_basis_size": len(graver),
        "factorization_degree_distribution": {
            str(degree): count
            for degree, count in sorted(degree_distribution.items())
        },
        "maximum_factorization_degree": max(
            degree_distribution,
            default=0,
        ),
        "all_primitive_relations_boundary_generated": True,
        "primitive_relations": records,
    }


def build_structural_certificate(
    *,
    maximum_wronskian_rank: int,
    maximum_affine_slope: int,
    maximum_affine_offset: int,
    singular: str,
    normaliz: str,
) -> dict[str, Any]:
    """Build exact regressions accompanying the structural theorems."""

    wronskian = confluent_wronskian_regression(maximum_wronskian_rank)
    common_base = common_base_span_regression()
    four_state = arbitrary_four_state_singular_certificate(singular)
    affine_factorial = affine_factorial_graver_certificate(
        maximum_slope=maximum_affine_slope,
        maximum_offset=maximum_affine_offset,
        normaliz=normaliz,
    )
    return {
        "model": (
            "Exact regressions for proved characteristic-zero structural "
            "theorems; this is not a Hall-promotion certificate or a "
            "GVC(2) counterexample."
        ),
        "parameters": {
            "maximum_wronskian_rank": maximum_wronskian_rank,
            "maximum_affine_slope": maximum_affine_slope,
            "maximum_affine_offset": maximum_affine_offset,
            "sympy": __import__("sympy").__version__,
            "singular": singular_version(singular),
            "normaliz": normaliz_version(normaliz),
        },
        "arbitrary_coefficient_affine_ray_splitting": {
            "statement": (
                "For a fixed finite family h_j(t)*f_j(t)^N, arbitrary "
                "scalar rows split for all sufficiently large N by "
                "proportionality classes of f_j."
            ),
            "proof": (
                "The normalized confluent Wronskian is a nonzero polynomial "
                "in N with the checked leading Vandermonde symbol."
            ),
            "regression": wronskian,
        },
        "arbitrary_coefficient_four_state_three_scale_pairing": {
            "statement": (
                "Every positive-dimensional ratio component of nonzero "
                "scale-1/2/3 diagonal moment rows on four states is a 2+2 "
                "pair block."
            ),
            "certificate": four_state,
        },
        "common_base_all_scale_ideal_collapse": {
            "statement": (
                "For F_N=f^N*sum_j c_j(N)h_j with f invertible on the "
                "coefficient torus, the complete Cartesian derivative ideal "
                "equals the derivative ideal generated by a basis of the "
                "finite coefficient span of the vectors c(N)."
            ),
            "proof": (
                "Leibniz is triangular in derivatives of the correction "
                "row after multiplying by f^(-N); constant-field linearity "
                "then replaces all scale rows by a coefficient-span basis."
            ),
            "regression": common_base,
        },
        "periodic_additive_coefficient_elimination": {
            "statement": (
                "On each residue class of an eventually periodic finite "
                "additive affine-factorial law, projective constancy is "
                "equivalent to rational-function proportionality inside "
                "every canonical signed slope-vector class."
            ),
            "proof": (
                "Substitute N=L*k+r, normalize integer factorial offsets "
                "to rational multiples of zero-offset rays, and apply "
                "factorial-trace independence over K(k)."
            ),
            "regression": common_base[
                "periodic_additive_factorial_example"
            ],
        },
        "integer_affine_factorial_boundary_presentation": {
            "statement": (
                "The kernel of the exact successor-divisor map on "
                "integer-affine factorial atoms is generated by elementary "
                "same-rational-boundary transfers; every proportional "
                "factorial ratio is constant."
            ),
            "proof": (
                "Telescope offsets to zero; zero-offset slope divisors are "
                "independent, and the residual translation edges lie on "
                "acyclic rational-orbit lines."
            ),
            "regression": affine_factorial,
        },
    }


def build_experiment(
    *,
    span: int,
    state_count: int,
    maximum_operator_count: int,
    maximum_polynomial_count: int,
    maximum_scale: int,
    require_factorial_pairing: bool,
    singular: str,
) -> dict[str, Any]:
    fibres = fibre_states(
        span,
        state_count,
        maximum_operator_count,
        maximum_polynomial_count,
    )
    candidates = []
    scale_one_closures = []
    for fibre, states in fibres.items():
        for selected_states in itertools.combinations(states, state_count):
            for tail_signs in itertools.product(
                (1, -1), repeat=state_count - 1
            ):
                signs = (1,) + tail_signs
                pairings = four_state_factorial_pairings(
                    selected_states,
                    signs,
                )
                if require_factorial_pairing and not pairings:
                    continue
                digest = candidate_key(selected_states, signs)
                polynomial = scaled_row(selected_states, signs, 1, "fixed")
                closure = derivative_closure(polynomial, "polynomial", span + 1)
                scale_one_closures.append((digest, "N1", span + 1, closure))
                candidates.append(
                    {
                        "digest": digest,
                        "fibre": fibre,
                        "states": selected_states,
                        "signs": signs,
                        "symmetry_key": symmetry_key(selected_states, signs),
                        "c2_character_feasible": character_feasible(
                            selected_states, signs, 2
                        ),
                        "c4_real_character_feasible": character_feasible(
                            selected_states, signs, 4
                        ),
                        "factorial_pairings": pairings,
                    }
                )

    scale_one = singular_torus_test(scale_one_closures, singular)
    survivors = []
    for candidate in candidates:
        digest = candidate["digest"]
        result = scale_one[(digest, "N1")]
        if result["torus_empty"]:
            continue
        candidate["scale_1_ideal"] = result
        candidate["scale_ideals"] = {
            law: {"1": result} for law in PHASE_LAWS
        }
        candidate["first_torus_empty_scale"] = {
            law: None for law in PHASE_LAWS
        }
        survivors.append(candidate)

    active = {law: list(survivors) for law in PHASE_LAWS}
    for maximum_imposed_scale in range(2, maximum_scale + 1):
        closures = []
        for phase_law in PHASE_LAWS:
            for candidate in active[phase_law]:
                generators = []
                for scale in range(1, maximum_imposed_scale + 1):
                    polynomial = scaled_row(
                        candidate["states"],
                        candidate["signs"],
                        scale,
                        phase_law,
                    )
                    generators.extend(
                        derivative_closure(
                            polynomial,
                            "polynomial",
                            span + 1,
                        )
                    )
                mode = f"{phase_law}_K{maximum_imposed_scale}"
                closures.append(
                    (candidate["digest"], mode, span + 1, tuple(generators))
                )

        scale_results = singular_torus_test(closures, singular)
        for phase_law in PHASE_LAWS:
            next_active = []
            mode = f"{phase_law}_K{maximum_imposed_scale}"
            for candidate in active[phase_law]:
                result = scale_results[(candidate["digest"], mode)]
                candidate["scale_ideals"][phase_law][
                    str(maximum_imposed_scale)
                ] = result
                if result["torus_empty"]:
                    candidate["first_torus_empty_scale"][phase_law] = (
                        maximum_imposed_scale
                    )
                else:
                    next_active.append(candidate)
            active[phase_law] = next_active

    # A fixed-sign all-scale row with four states can survive factorial-trace
    # independence only through such pairings.  Test the pair equalities
    # themselves.  A proper ideal is then an exact all-N pair cancellation,
    # not merely a finite-prefix survivor.
    pair_closures = []
    for candidate in active["fixed"]:
        for pairing_index, pairing in enumerate(candidate["factorial_pairings"]):
            for mode in ("polynomial", "operator", "independent"):
                generators = []
                for left, right in pairing:
                    polynomial = {
                        state_exponent(candidate["states"][left]): 1,
                        state_exponent(candidate["states"][right]): -1,
                    }
                    generators.extend(
                        derivative_closure(polynomial, mode, span + 1)
                    )
                pair_closures.append(
                    (
                        candidate["digest"],
                        f"pair_{pairing_index}_{mode}",
                        span + 1,
                        tuple(generators),
                    )
                )
    pair_results = singular_torus_test(pair_closures, singular)
    veronese_requests = []
    for candidate in active["fixed"]:
        pairing_records = []
        for pairing_index, pairing in enumerate(candidate["factorial_pairings"]):
            orbit_ideals = {
                mode: pair_results[
                    (
                        candidate["digest"],
                        f"pair_{pairing_index}_{mode}",
                    )
                ]
                for mode in ("polynomial", "operator", "independent")
            }
            pairing_records.append(
                {
                    "pairs": [list(pair) for pair in pairing],
                    "orbit_ideals": orbit_ideals,
                }
            )
            if span == 2 and not orbit_ideals["independent"]["torus_empty"]:
                veronese_requests.append(
                    (
                        candidate["digest"],
                        pairing_index,
                        candidate["states"],
                        pairing,
                    )
                )
        candidate["fixed_all_scale_factorial_pairings"] = pairing_records
    veronese_results = quadratic_veronese_pair_tests(
        veronese_requests,
        singular,
    )
    for candidate in active["fixed"]:
        for pairing_index, pairing_record in enumerate(
            candidate["fixed_all_scale_factorial_pairings"]
        ):
            pairing_record["quadratic_veronese_product_ideal"] = (
                veronese_results.get((candidate["digest"], pairing_index), False)
            )

    phase_counts = {}
    for phase_law in PHASE_LAWS:
        obstruction_scales = Counter(
            candidate["first_torus_empty_scale"][phase_law]
            for candidate in survivors
        )
        phase_counts[phase_law] = {
            "first_torus_empty_scale_distribution": {
                ("survives" if scale is None else str(scale)): count
                for scale, count in sorted(
                    obstruction_scales.items(),
                    key=lambda item: (
                        item[0] is None,
                        item[0] if item[0] is not None else maximum_scale + 1,
                    ),
                )
            },
            "survivors_through_maximum_scale": len(active[phase_law]),
        }

    records = []
    for candidate in survivors:
        polynomial_levels = tuple(
            side_level(state.polynomial) for state in candidate["states"]
        )
        records.append(
            {
                "digest": candidate["digest"],
                "fibre": {
                    "operator_count": candidate["fibre"][0],
                    "polynomial_count": candidate["fibre"][1],
                    "total_level": candidate["fibre"][2],
                },
                "states": [state_record(state) for state in candidate["states"]],
                "signs": list(candidate["signs"]),
                "distinct_polynomial_levels": sorted(set(polynomial_levels)),
                "transfers_level": len(set(polynomial_levels)) > 1,
                "c2_character_feasible": candidate["c2_character_feasible"],
                "c4_real_character_feasible": candidate[
                    "c4_real_character_feasible"
                ],
                "scale_ideals": candidate["scale_ideals"],
                "first_torus_empty_scale": candidate[
                    "first_torus_empty_scale"
                ],
                "fixed_all_scale_factorial_pairings": candidate.get(
                    "fixed_all_scale_factorial_pairings",
                    [],
                ),
            }
        )

    candidate_orbits = {candidate["symmetry_key"] for candidate in candidates}
    survivor_orbits = {candidate["symmetry_key"] for candidate in survivors}
    transferring_survivors = sum(record["transfers_level"] for record in records)
    if state_count == 3 and any(
        not record["scale_ideals"][law]["2"]["torus_empty"]
        for record in records
        for law in PHASE_LAWS
        if record["transfers_level"]
    ):
        raise AssertionError("a transferring three-state two-scale survivor exists")
    fixed_scale_survivors = {
        law: [
            record["digest"]
            for record in records
            if record["first_torus_empty_scale"][law] is None
        ]
        for law in PHASE_LAWS
    }
    pair_survivors_by_mode = {
        mode: [
            candidate["digest"]
            for candidate in active["fixed"]
            if any(
                not pairing["orbit_ideals"][mode]["torus_empty"]
                for pairing in candidate.get(
                    "fixed_all_scale_factorial_pairings",
                    [],
                )
            )
        ]
        for mode in ("polynomial", "operator", "independent")
    }

    sign_candidates = len(candidates)
    theorem_certificate: dict[str, Any]
    if state_count == 3:
        theorem_certificate = {
            "characteristic": 0,
            "number_of_states": 3,
            "scales": [1, 2],
            "scale_coefficients": "arbitrary nonzero scalars",
            "elimination": (
                "A3=-A1-A2 makes the scale-two row a nonzero homogeneous "
                "quadratic; after z=A2/A1 its linear coefficient is "
                "2*lambda3, so z is constant."
            ),
            "conclusion": (
                "All three marked-side orbit functions are proportional.  "
                "Their Taylor degrees, hence marked-side levels, agree."
            ),
            "remaining_case": (
                "Equal marked-side levels; every pair has zero colour-level "
                "transfer and splits into one-colour return relations."
            ),
        }
    else:
        theorem_certificate = {
            "characteristic": 0,
            "number_of_states": state_count,
            "projective_moment_reduction": (
                "After the scale-one row eliminates one orbit function, "
                f"the scale-N rows cut P^{state_count - 2}.  A transferring "
                "orbit can survive only on a positive-dimensional common "
                "component; a zero-dimensional ratio variety makes every "
                "orbit-function ratio constant and forces equal Taylor "
                "degrees."
            ),
            "four_state_specialization": (
                "For four states the scale-two row is a conic in P^2.  "
                "Scale three adds a cubic.  With arbitrary nonzero row "
                "coefficients, a common positive-dimensional component "
                "forces the conic to have rank two; its common line is "
                "exactly a 2+2 proportional pair block.  A smooth conic "
                "cannot divide the cubic."
                if state_count == 4
                else None
            ),
            "fixed_template_all_scale_pairing": (
                "Factorial trace independence splits unequal scalar "
                "factorial rays.  In each same-ray class, partial fractions "
                "of sum_j epsilon_j*M_j*z/(1-M_j*z) pair every positive "
                "orbit monomial with an equal negative one."
                if state_count == 4
                else None
            ),
        }
    summary = {
        "fibres_with_enough_states": len(fibres),
        "state_count": state_count,
        "state_subsets": len(
            {candidate["states"] for candidate in candidates}
        ),
        "signed_candidates_modulo_common_sign": sign_candidates,
        "candidate_symmetry_orbits": len(candidate_orbits),
        "c2_character_feasible_candidates": sum(
            candidate["c2_character_feasible"] for candidate in candidates
        ),
        "c4_real_character_feasible_candidates": sum(
            candidate["c4_real_character_feasible"] for candidate in candidates
        ),
        "scale_1_torus_survivors": len(survivors),
        "scale_1_survivor_symmetry_orbits": len(survivor_orbits),
        "scale_1_transferring_survivors": transferring_survivors,
        "scale_1_zero_transfer_survivors": len(survivors) - transferring_survivors,
        "all_scale_prefixes": phase_counts,
        "survivor_digests_through_maximum_scale": fixed_scale_survivors,
        "fixed_exact_all_scale_pair_survivors": {
            mode: {
                "count": len(digests),
                "digests": digests,
            }
            for mode, digests in pair_survivors_by_mode.items()
        },
        "fixed_independent_survivors_on_quadratic_veronese_product": sum(
            any(
                pairing["quadratic_veronese_product_ideal"]
                for pairing in candidate.get(
                    "fixed_all_scale_factorial_pairings",
                    [],
                )
            )
            for candidate in active["fixed"]
        ),
    }
    return {
        "model": (
            "Bounded projected two-colour small-shell census.  A proper "
            "ideal is only a projected orbit survivor, not an affine Hall "
            "packet or GVC(2) counterexample."
        ),
        "parameters": {
            "span": span,
            "state_count": state_count,
            "maximum_operator_count": maximum_operator_count,
            "maximum_polynomial_count": maximum_polynomial_count,
            "maximum_scale": maximum_scale,
            "require_factorial_pairing": require_factorial_pairing,
            "translation_mode": "polynomial colour only",
            "phase_laws": list(PHASE_LAWS),
            "singular": singular_version(singular),
            "saturation": "all R_i and B_i coefficient channels",
        },
        "theorem_certificate": theorem_certificate,
        "summary": summary,
        "scale_1_survivors": records,
    }


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--span", type=int, default=2)
    parser.add_argument("--state-count", type=int, default=3)
    parser.add_argument("--maximum-operator-count", type=int, default=2)
    parser.add_argument("--maximum-polynomial-count", type=int, default=2)
    parser.add_argument("--maximum-scale", type=int, default=3)
    parser.add_argument(
        "--require-factorial-pairing",
        action="store_true",
        help=(
            "for four states, retain only the opposite-sign equal-factorial "
            "pairings which can survive an actual fixed-sign all-scale row"
        ),
    )
    parser.add_argument(
        "--structural-certificate",
        action="store_true",
        help=(
            "skip the bounded census and verify the affine-ray Wronskian "
            "symbol plus the four-state conic/cubic certificate"
        ),
    )
    parser.add_argument(
        "--maximum-wronskian-rank",
        type=int,
        default=7,
        help="largest rank in the exact confluent-Vandermonde regression",
    )
    parser.add_argument(
        "--maximum-affine-slope",
        type=int,
        default=5,
        help="largest factorial slope in the exact-divisor Graver regression",
    )
    parser.add_argument(
        "--maximum-affine-offset",
        type=int,
        default=3,
        help="symmetric integer offset radius in the Graver regression",
    )
    parser.add_argument("--singular", default="Singular")
    parser.add_argument("--normaliz", default="normaliz")
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    if arguments.span < 1:
        parser.error("--span must be positive")
    if arguments.state_count < 3:
        parser.error("--state-count must be at least three")
    if arguments.maximum_operator_count < 1:
        parser.error("--maximum-operator-count must be positive")
    if arguments.maximum_polynomial_count < 1:
        parser.error("--maximum-polynomial-count must be positive")
    if arguments.maximum_scale < 2:
        parser.error("--maximum-scale must be at least two")
    if arguments.maximum_wronskian_rank < 1:
        parser.error("--maximum-wronskian-rank must be positive")
    if arguments.maximum_affine_slope < 1:
        parser.error("--maximum-affine-slope must be positive")
    if arguments.maximum_affine_offset < 0:
        parser.error("--maximum-affine-offset must be nonnegative")
    if arguments.require_factorial_pairing and arguments.state_count != 4:
        parser.error("--require-factorial-pairing requires --state-count 4")
    singular = shutil.which(arguments.singular)
    if singular is None:
        parser.error(f"Singular executable not found: {arguments.singular}")

    if arguments.structural_certificate:
        normaliz = shutil.which(arguments.normaliz)
        if normaliz is None:
            parser.error(f"Normaliz executable not found: {arguments.normaliz}")
        result = build_structural_certificate(
            maximum_wronskian_rank=arguments.maximum_wronskian_rank,
            maximum_affine_slope=arguments.maximum_affine_slope,
            maximum_affine_offset=arguments.maximum_affine_offset,
            singular=singular,
            normaliz=normaliz,
        )
    else:
        result = build_experiment(
            span=arguments.span,
            state_count=arguments.state_count,
            maximum_operator_count=arguments.maximum_operator_count,
            maximum_polynomial_count=arguments.maximum_polynomial_count,
            maximum_scale=arguments.maximum_scale,
            require_factorial_pairing=arguments.require_factorial_pairing,
            singular=singular,
        )
    if arguments.output is not None:
        write_json(arguments.output, result)

    print(f"Singular: {result['parameters']['singular']}")
    if arguments.structural_certificate:
        regression = result["arbitrary_coefficient_affine_ray_splitting"][
            "regression"
        ]
        print(
            "Structural certificate: "
            f"Wronskian ranks<= {regression['maximum_rank']}, "
            f"partitions={regression['partitions_checked']}, "
            "common-base ideal=coefficient span, "
            "smooth four-state branch=empty, rank-two branch=pair blocks, "
            "affine-factorial Graver moves=boundary transfers, periodic "
            "additive rows=rational elimination"
        )
        print(f"logical_sha256={logical_hash(result)}")
        return

    summary = result["summary"]
    print(
        f"{arguments.state_count}-state census: "
        f"fibres={summary['fibres_with_enough_states']}, "
        f"subsets={summary['state_subsets']}, "
        f"signed={summary['signed_candidates_modulo_common_sign']}, "
        f"scale-1 survivors={summary['scale_1_torus_survivors']}"
    )
    for law in PHASE_LAWS:
        counts = summary["all_scale_prefixes"][law]
        print(
            f"{law}, through scale {arguments.maximum_scale}: "
            "first obstruction="
            f"{counts['first_torus_empty_scale_distribution']}, "
            "survivors="
            f"{counts['survivors_through_maximum_scale']}"
        )
    print(f"logical_sha256={logical_hash(result)}")


if __name__ == "__main__":
    main()
