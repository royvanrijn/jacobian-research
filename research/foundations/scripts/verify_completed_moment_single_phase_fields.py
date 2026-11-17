#!/usr/bin/env python3
"""Exact single-phase moment-field theorems in degrees three and five.

For each nonzero phase h, this checker chooses matching tau-eigendirections
B of phase +h and C of phase -h.  The +1 eigenspace is used whenever it
is nonzero.  In odd degree the extreme phase h=d is one-dimensional and
tau-odd, so the forced -1 eigendirections are used there.  On

    diag(a_0,...,a_d) + b B + c C

the residual diagonal torus acts with opposite weights on b and c.
Moments therefore descend to the quotient coordinates

    a_0,...,a_d,z=bc,

and apolar adjunction reverses the a_i while fixing z.

For d=3 and d=5, and for every phase h=1,...,d, the checker proves:

* the first d+2 moments have full Jacobian rank, while their common zero
  fiber modulo 32003 is one-dimensional;
* adding mu_(d+3) makes the moment-origin fiber zero-dimensional;
* the first d+3 moment fiber through the recorded integral point is
  exactly the two reduced points related by coordinate reversal.

Origin finiteness and the two-point fiber are certified in characteristic
32003.  For the fixed-target fiber, weighted homogenization introduces a
variable t and equations mu_r-target_r*t^r.  Origin finiteness excludes
points at t=0.  Properness over Z_(32003) and Nakayama therefore bound
the characteristic-zero fiber length by its two-dimensional reduced
special fiber.  The two explicit rational reversal points give equality.
These are slice theorems; they do not determine the moment-field degree
on the full invariant quotient R_d.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from functools import reduce
from math import factorial
from operator import mul
from pathlib import Path

from research_completed_moment_algebra import invariant_values_exact
from verify_degree_four_diagonal_moment_field import run_singular


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "completed_moment_single_phase_fields.json"
)
PARAMETER_PRIME = 32003
POINTS = {
    3: (2, 3, 5, 7, 221),
    5: (2, 3, 5, 7, 11, 13, 221),
}
EXPECTED_ORIGIN_LENGTHS = {3: 54, 5: 1934}
Exponent = tuple[int, ...]
Polynomial = dict[Exponent, int]
Direction = dict[tuple[int, int], int]


def tau_position(
    d: int,
    position: tuple[int, int],
) -> tuple[int, int]:
    row, column = position
    return d - column, d - row


def tau_eigen_directions(
    d: int,
    phase: int,
    positive: bool,
    eigenvalue: int,
) -> list[Direction]:
    assert eigenvalue in (-1, 1)
    positions = [
        (
            (index + phase, index)
            if positive
            else (index, index + phase)
        )
        for index in range(d + 1 - phase)
    ]
    seen: set[int] = set()
    result = []
    for index, position in enumerate(positions):
        if index in seen:
            continue
        partner = d - phase - index
        seen.update((index, partner))
        direction = {position: 1}
        if partner != index:
            direction[positions[partner]] = (
                eigenvalue * (-1) ** phase
            )
        for target, coefficient in direction.items():
            source = tau_position(d, target)
            transformed = (
                (-1) ** sum(target) * direction.get(source, 0)
            )
            if transformed != eigenvalue * coefficient:
                break
        else:
            result.append(direction)
    return result


def chosen_tau_eigendirections(
    d: int,
    phase: int,
) -> tuple[int, Direction, Direction]:
    for eigenvalue in (1, -1):
        positive = tau_eigen_directions(
            d, phase, True, eigenvalue
        )
        negative = tau_eigen_directions(
            d, phase, False, eigenvalue
        )
        if positive and negative:
            # A cross pair detects apolar orientation whenever both
            # eigendirection spaces have dimension at least two.
            return eigenvalue, positive[0], negative[-1]
    raise AssertionError("every phase has a tau eigendirection")


def restricted_moments(
    d: int,
    positive: Direction,
    negative: Direction,
    cutoff: int,
) -> list[Polynomial]:
    """Return moments in quotient variables a_0,...,a_d,z=bc."""

    b_index = d + 1
    c_index = d + 2
    terms = [(index, index, index, 1) for index in range(d + 1)]
    terms.extend(
        (row, column, b_index, coefficient)
        for (row, column), coefficient in positive.items()
    )
    terms.extend(
        (row, column, c_index, coefficient)
        for (row, column), coefficient in negative.items()
    )
    state: dict[tuple[int, int, tuple[int, ...]], int] = {
        (0, 0, (0,) * (d + 3)): 1
    }
    moments = []
    for _order in range(1, cutoff + 1):
        updated: dict[tuple[int, int, tuple[int, ...]], int] = {}
        for (left, right, exponents), value in state.items():
            for delta_left, delta_right, variable, coefficient in terms:
                new_exponents = list(exponents)
                new_exponents[variable] += 1
                key = (
                    left + delta_left,
                    right + delta_right,
                    tuple(new_exponents),
                )
                updated[key] = updated.get(key, 0) + value * coefficient
        state = {key: value for key, value in updated.items() if value}

        moment: Polynomial = {}
        for (left, right, exponents), value in state.items():
            if left != right or exponents[b_index] != exponents[c_index]:
                continue
            quotient_exponents = (
                exponents[: d + 1] + (exponents[b_index],)
            )
            moment[quotient_exponents] = (
                moment.get(quotient_exponents, 0)
                + factorial(left)
                * factorial(d * (_order) - left)
                * value
            )
        moments.append(
            {
                exponents: value
                for exponents, value in moment.items()
                if value
            }
        )
    return moments


def variable_names(d: int) -> tuple[str, ...]:
    return tuple([f"a{index}" for index in range(d + 1)] + ["z"])


def polynomial_string(d: int, polynomial: Polynomial) -> str:
    terms = []
    for exponents, coefficient in polynomial.items():
        sign = "+" if coefficient >= 0 and terms else ""
        absolute = abs(coefficient)
        factors = []
        if absolute != 1 or not any(exponents):
            factors.append(str(absolute))
        for variable, exponent in zip(
            variable_names(d), exponents, strict=True
        ):
            if exponent == 1:
                factors.append(variable)
            elif exponent:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors) if factors else "0"
        terms.append(sign + ("" if coefficient >= 0 else "-") + monomial)
    return "".join(terms) if terms else "0"


def evaluate(
    polynomial: Polynomial,
    point: tuple[int, ...],
) -> int:
    return sum(
        coefficient
        * reduce(
            mul,
            (
                value**exponent
                for value, exponent in zip(point, exponents, strict=True)
            ),
            1,
        )
        for exponents, coefficient in polynomial.items()
    )


def direction_payload(direction: Direction) -> list[dict[str, int]]:
    return [
        {"row": row, "column": column, "coefficient": coefficient}
        for (row, column), coefficient in sorted(direction.items())
    ]


def odd_invariant_value(
    d: int,
    positive: Direction,
    negative: Direction,
    point: tuple[int, ...],
) -> str:
    matrix = [[0] * (d + 1) for _ in range(d + 1)]
    for index, value in enumerate(point[:-1]):
        matrix[index][index] = value
    for (row, column), value in positive.items():
        matrix[row][column] += value
    for (row, column), value in negative.items():
        matrix[row][column] += point[-1] * value
    odd = invariant_values_exact(matrix, d)[1]
    return str(odd)


def modular_determinant(matrix: list[list[int]], prime: int) -> int:
    work = [[entry % prime for entry in row] for row in matrix]
    determinant = 1
    for column in range(len(work)):
        pivot = next(
            (
                row
                for row in range(column, len(work))
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            determinant = -determinant
        pivot_value = work[column][column]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, len(work)):
            factor = work[row][column] * inverse % prime
            for index in range(column, len(work)):
                work[row][index] = (
                    work[row][index]
                    - factor * work[column][index]
                ) % prime
    return determinant % prime


def jacobian_determinant(
    moments: list[Polynomial],
    point: tuple[int, ...],
) -> int:
    matrix = []
    for polynomial in moments[: len(point)]:
        row = []
        for variable_index in range(len(point)):
            derivative: Polynomial = {}
            for exponents, coefficient in polynomial.items():
                exponent = exponents[variable_index]
                if not exponent:
                    continue
                derived_exponents = list(exponents)
                derived_exponents[variable_index] -= 1
                derivative[tuple(derived_exponents)] = (
                    coefficient * exponent
                )
            row.append(evaluate(derivative, point) % PARAMETER_PRIME)
        matrix.append(row)
    determinant = modular_determinant(matrix, PARAMETER_PRIME)
    assert determinant
    return determinant


def origin_finiteness_certificate(
    d: int,
    moment_strings: list[str],
) -> dict[str, int]:
    names = ",".join(variable_names(d))
    code = (
        f"ring r={PARAMETER_PRIME},({names}),dp;\n"
        + "ideal Ishort="
        + ",".join(moment_strings[: d + 2])
        + ";\n"
        + "ideal Gshort=std(Ishort);\n"
        + "dim(Gshort);\n"
        + "vdim(Gshort);\n"
        + "ideal Ifinite="
        + ",".join(moment_strings[: d + 3])
        + ";\n"
        + "ideal Gfinite=std(Ifinite);\n"
        + "dim(Gfinite);\n"
        + "vdim(Gfinite);\n"
        + "size(Gfinite);\n"
    )
    output = run_singular(code)
    short_dimension = int(output[0])
    short_vdim = int(output[1])
    finite_dimension = int(output[2])
    finite_length = int(output[3])
    finite_basis_size = int(output[4])
    assert short_dimension == 1
    assert short_vdim == -1
    assert finite_dimension == 0
    assert finite_length == EXPECTED_ORIGIN_LENGTHS[d]
    return {
        "special_first_d_plus_2_origin_fiber_dimension": short_dimension,
        "special_first_d_plus_3_origin_fiber_dimension": finite_dimension,
        "special_first_d_plus_3_origin_quotient_length": finite_length,
        "special_first_d_plus_3_origin_standard_basis_size": (
            finite_basis_size
        ),
    }


def adapted_fiber_certificate(
    d: int,
    moment_strings: list[str],
    targets: list[int],
    point: tuple[int, ...],
) -> dict[str, object]:
    reversed_point = tuple(reversed(point[:-1])) + (point[-1],)
    midpoint = [
        Fraction(left + right, 2)
        for left, right in zip(point, reversed_point, strict=True)
    ]
    direction = [
        Fraction(left - right, 2)
        for left, right in zip(point, reversed_point, strict=True)
    ]
    assert direction[d] != 0

    def rational_string(value: Fraction) -> str:
        if value.denominator == 1:
            return str(value.numerator)
        return f"{value.numerator}/{value.denominator}"

    y_names = [f"y{index}" for index in range(d)]
    target_names = ["s"] + y_names + ["w"]
    map_entries = [
        (
            f"({rational_string(midpoint[index])})"
            f"+({rational_string(direction[index])})*s"
            f"+{y_names[index]}"
        )
        for index in range(d)
    ]
    map_entries.append(
        f"({rational_string(midpoint[d])})"
        f"+({rational_string(direction[d])})*s"
    )
    map_entries.append(f"{point[-1]}+w")
    source_ideal = ",".join(
        f"({moment})-({target})"
        for moment, target in zip(
            moment_strings, targets, strict=True
        )
    )
    expected_generators = y_names + ["w", "s^2-1"]
    code = (
        f"ring ra={PARAMETER_PRIME},"
        f"({','.join(variable_names(d))}),dp;\n"
        + f"ideal IA={source_ideal};\n"
        + f"ring rb={PARAMETER_PRIME},"
        f"({','.join(target_names)}),dp;\n"
        + f"map phi=ra,{','.join(map_entries)};\n"
        + "ideal I=phi(IA);\n"
        + f"ideal J={','.join(expected_generators)};\n"
        + "ideal GI=std(I);\n"
        + "ideal GJ=std(J);\n"
        + "vdim(GI);\n"
        + "size(GI);\n"
        + "size(GJ);\n"
        + "reduce(I,GJ);\n"
        + "reduce(J,GI);\n"
    )
    output = run_singular(code)
    fiber_length = int(output[0])
    fiber_basis_size = int(output[1])
    expected_basis_size = int(output[2])
    reduction_lines = output[3:]
    assert fiber_length == 2
    assert expected_basis_size == d + 2
    assert len(reduction_lines) == (d + 3) + (d + 2)
    assert all(line.endswith("=0") for line in reduction_lines)
    return {
        "test_point": list(point),
        "reversed_point": list(reversed_point),
        "first_d_plus_3_moment_values": [str(value) for value in targets],
        "adapted_coordinates": {
            "midpoint": [rational_string(value) for value in midpoint],
            "direction": [rational_string(value) for value in direction],
            "map": map_entries,
            "expected_fiber_ideal": expected_generators,
        },
        "fiber_length": fiber_length,
        "fiber_standard_basis_size": fiber_basis_size,
        "expected_ideal_standard_basis_size": expected_basis_size,
        "mutual_special_fiber_ideal_reductions_zero": True,
        "characteristic_zero_fiber_length": fiber_length,
        "characteristic_zero_lift_argument": (
            "weighted homogenization has no points at infinity because "
            "the moment-origin fiber is finite; properness and Nakayama "
            "bound the characteristic-zero fiber length by the reduced "
            "special fiber length, while the two rational reversal "
            "points give equality"
        ),
    }


def verify_phase(d: int, phase: int) -> dict[str, object]:
    eigenvalue, positive, negative = chosen_tau_eigendirections(d, phase)
    moments = restricted_moments(d, positive, negative, d + 3)
    moment_strings = [
        polynomial_string(d, moment) for moment in moments
    ]
    point = POINTS[d]
    targets = [evaluate(moment, point) for moment in moments]
    reversed_point = tuple(reversed(point[:-1])) + (point[-1],)
    assert [
        evaluate(moment, reversed_point) for moment in moments
    ] == targets

    determinant = jacobian_determinant(moments, point)
    origin = origin_finiteness_certificate(d, moment_strings)
    fiber = adapted_fiber_certificate(
        d, moment_strings, targets, point
    )
    odd_value = odd_invariant_value(d, positive, negative, point)
    expected_odd = (
        "-273686400/7" if d == 5 and phase in (1, 2) else "0"
    )
    assert odd_value == expected_odd
    print(
        f"PASS d={d}, h={phase}: full Jacobian rank, finite "
        "moment origin, characteristic-zero fiber length 2"
    )
    return {
        "phase": phase,
        "tau_eigenvalue_of_directions": eigenvalue,
        "positive_direction": direction_payload(positive),
        "negative_direction": direction_payload(negative),
        "first_apolar_odd_invariant_at_lift_b_1_c_z": odd_value,
        "genuinely_apolar_moving_quotient_test": odd_value != "0",
        "moment_term_counts_orders_1_through_d_plus_3": [
            len(moment) for moment in moments
        ],
        "parameter_prime": PARAMETER_PRIME,
        "first_d_plus_2_jacobian_determinant_mod_prime": determinant,
        **origin,
        "characteristic_zero_origin_finiteness_argument": (
            "weighted-homogeneous finite special fiber makes the "
            "weighted-projective characteristic-zero fiber empty by "
            "properness, so the slice coordinate ring is integral over "
            "the algebra generated by the first d+3 moments"
        ),
        **fiber,
        "generic_full_moment_degree_on_slice": 2,
        "full_moment_field_on_slice": "apolar-reversal-fixed field",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--degrees",
        nargs="+",
        type=int,
        choices=(3, 5),
        default=(3, 5),
    )
    parser.add_argument(
        "--phases",
        nargs="+",
        type=int,
        help="optional phase subset, applied to every selected degree",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    payload: dict[str, object] = {
        "format": "completed-moment-single-phase-fields-v1",
        "degrees": {},
    }
    degrees_payload: dict[str, object] = payload["degrees"]  # type: ignore[assignment]
    for d in arguments.degrees:
        phases = arguments.phases or list(range(1, d + 1))
        if any(phase < 1 or phase > d for phase in phases):
            raise SystemExit(f"phases must lie in 1,...,{d}")
        phase_payloads = [verify_phase(d, phase) for phase in phases]
        degrees_payload[str(d)] = {
            "quotient_coordinates": [
                *[f"a{index}" for index in range(d + 1)],
                "z=bc",
            ],
            "apolar_involution": (
                f"(a0,...,a{d},z) -> "
                f"(a{d},...,a0,z)"
            ),
            "certified_phase_count": len(phase_payloads),
            "phases": phase_payloads,
            "open_family_consequence": (
                "in every certified phase, finiteness and reduced "
                "two-point fiber persist on a nonempty Zariski-open "
                "family of matching tau-eigendirection pairs"
            ),
            "scope_warning": (
                "single-phase quotient slices only; no determination "
                f"of the degree on Frac(R_{d})"
            ),
        }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
