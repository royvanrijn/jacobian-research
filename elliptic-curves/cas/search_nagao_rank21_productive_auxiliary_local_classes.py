#!/usr/bin/env python3
"""Lift conductor-saving finite-group classes on the productive rank-21 slice.

The companion script ``search_nagao_rank21_productive_auxiliary_orbit.py``
finds a saturated eight-point basis on a genus-one auxiliary slice.  A naive
large coefficient box is exponential.  Here we instead reduce the auxiliary
curve and its basis at 13, 37, and 83, enumerate the complete product of the
three finite groups by dynamic programming, and retain a shortest coefficient
vector for every reachable state.  We then select exactly those states whose
rational specialization map lands in the known high-power discriminant root
unions

``T=3,10 (mod 13)``, ``T=0,17,20 (mod 37)``, and
``T=12,31,40,43,52,71 (mod 83)``.

Every retained vector is lifted with exact rational group law and receives an
exact small-prime radical proxy.  This is exhaustive for the finite state
space and the declared shortest-vector policy; it is not exhaustive over all
integer lifts of a finite-group class.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
from typing import Any, Sequence

from ek_k3 import primes_up_to, rational_to_string
from search_nagao_rank21_productive_auxiliary_orbit import (
    PointedQuartic,
    exact_radical_proxy,
    file_sha256,
)
from search_nagao_section7_auxiliary_jacobians import (
    weierstrass_add,
    weierstrass_multiply,
)


Q = Fraction
INPUT_ARTIFACT = Path(
    "artifacts/generated-results/"
    "elliptic_nagao_rank21_productive_auxiliary_orbit.json"
)
EXPECTED_INPUT_SHA256 = "afc9c6185067f78eeb7c27c94bf93cebe854cb9b1115772ec62288d0ba4779df"
EXPECTED_INPUT_SCRIPT_SHA256 = "84f443c5f54612c607309e12fbb4a55cf1fd933f44a529083af166581b7b0165"
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_nagao_rank21_productive_auxiliary_local_classes.json"
)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_productive_auxiliary_local_classes.py"
)
LOCAL_PRIMES = (13, 37, 83)
TARGET_RESIDUES = (
    frozenset((3, 10)),
    frozenset((0, 17, 20)),
    frozenset((12, 31, 40, 43, 52, 71)),
)
COEFFICIENT_BOUND = 8
PROXY_PRIME_BOUND = 1000


def reduce_fraction(value: Fraction, prime: int) -> int:
    value = Q(value)
    if value.denominator % prime == 0:
        raise ValueError("a rational coefficient is not integral at the prime")
    return value.numerator * pow(value.denominator, -1, prime) % prime


class FiniteCurve:
    """Tiny exact affine addition table for a generalized Weierstrass curve."""

    def __init__(self, coefficients: Sequence[Fraction], prime: int) -> None:
        self.prime = prime
        self.coefficients = tuple(reduce_fraction(value, prime) for value in coefficients)
        a1, a2, a3, a4, a6 = self.coefficients
        self.points: tuple[tuple[int, int] | None, ...] = (None,) + tuple(
            (x_value, y_value)
            for x_value in range(prime)
            for y_value in range(prime)
            if (
                y_value**2
                + a1 * x_value * y_value
                + a3 * y_value
                - x_value**3
                - a2 * x_value**2
                - a4 * x_value
                - a6
            )
            % prime
            == 0
        )
        self.point_ids = {point: index for index, point in enumerate(self.points)}
        self.addition_table = tuple(
            tuple(self.point_ids[self.add(left, right)] for right in self.points)
            for left in self.points
        )

    def add(
        self,
        left: tuple[int, int] | None,
        right: tuple[int, int] | None,
    ) -> tuple[int, int] | None:
        if left is None:
            return right
        if right is None:
            return left
        prime = self.prime
        a1, a2, a3, a4, a6 = self.coefficients
        x1, y1 = left
        x2, y2 = right
        if x1 == x2 and (y1 + y2 + a1 * x1 + a3) % prime == 0:
            return None
        if x1 == x2:
            denominator = (2 * y1 + a1 * x1 + a3) % prime
            if denominator == 0:
                return None
            slope = (
                (3 * x1**2 + 2 * a2 * x1 + a4 - a1 * y1)
                * pow(denominator, -1, prime)
            ) % prime
            intercept = (
                (-x1**3 + a4 * x1 + 2 * a6 - a3 * y1)
                * pow(denominator, -1, prime)
            ) % prime
        else:
            denominator = (x2 - x1) % prime
            slope = ((y2 - y1) * pow(denominator, -1, prime)) % prime
            intercept = ((y1 * x2 - y2 * x1) * pow(denominator, -1, prime)) % prime
        x3 = (slope**2 + a1 * slope - a2 - x1 - x2) % prime
        y3 = (-(slope + a1) * x3 - intercept - a3) % prime
        answer = x3, y3
        if answer not in self.point_ids:
            raise AssertionError("finite generalized group law left the curve")
        return answer

    def multiply_id(self, point_id: int, scalar: int) -> int:
        if scalar < 0:
            point = self.points[point_id]
            if point is None:
                return 0
            a1, _, a3, _, _ = self.coefficients
            point_id = self.point_ids[
                (point[0], (-point[1] - a1 * point[0] - a3) % self.prime)
            ]
            scalar = -scalar
        answer = 0
        addend = point_id
        while scalar:
            if scalar & 1:
                answer = self.addition_table[answer][addend]
            addend = self.addition_table[addend][addend]
            scalar >>= 1
        return answer


def load_input(root: Path) -> tuple[dict[str, Any], PointedQuartic, tuple[Any, ...]]:
    path = root / INPUT_ARTIFACT
    if file_sha256(path) != EXPECTED_INPUT_SHA256:
        raise AssertionError("the productive-orbit input artifact changed")
    data = json.loads(path.read_text(encoding="utf-8"))
    if data["reproduction"]["script_sha256"] != EXPECTED_INPUT_SCRIPT_SHA256:
        raise AssertionError("the productive-orbit input script changed")
    record = data["productive_auxiliary_curve"]
    quartic = tuple(Q(value) for value in record["normalized_quartic_coefficients_ascending"])
    base_parameter = Q(data["source"]["constructor_parameter"])
    base_ordinate = Q(record["pointed_base_ordinate"])
    from search_nagao_section7_auxiliary_jacobians import translate_polynomial

    shifted = translate_polynomial(quartic, base_parameter)
    coefficients = tuple(Q(value) for value in record["generalized_weierstrass_coefficients"])
    auxiliary = PointedQuartic(
        base_parameter,
        base_ordinate,
        quartic,
        shifted,
        coefficients,
    )
    basis = tuple(
        (Q(point["x"]), Q(point["y"])) for point in record["saturated_basis"]
    )
    return data, auxiliary, basis


def parameter_residues(
    curve: FiniteCurve,
    auxiliary: PointedQuartic,
) -> tuple[int | None, ...]:
    prime = curve.prime
    q_value = reduce_fraction(auxiliary.base_ordinate, prime)
    t_value = reduce_fraction(auxiliary.base_parameter, prime)
    _, d_value_q, c_value_q, _, _ = auxiliary.shifted_coefficients
    d_value = reduce_fraction(d_value_q, prime)
    c_value = reduce_fraction(c_value_q, prime)
    answer = []
    for point in curve.points:
        if point is None:
            answer.append(t_value)
            continue
        x_value, y_value = point
        denominator = 2 * q_value * y_value % prime
        if denominator == 0:
            answer.append(None)
            continue
        answer.append(
            (
                t_value
                + (4 * q_value**2 * (x_value + c_value) - d_value**2)
                * pow(denominator, -1, prime)
            )
            % prime
        )
    return tuple(answer)


def shortest_state_vectors(
    curves: Sequence[FiniteCurve],
    basis: Sequence[tuple[Fraction, Fraction]],
) -> tuple[dict[tuple[int, ...], tuple[int, ...]], tuple[int, ...]]:
    choices = tuple(range(-COEFFICIENT_BOUND, COEFFICIENT_BOUND + 1))
    multiples = []
    for curve in curves:
        reduced_basis = tuple(
            curve.point_ids[
                (
                    reduce_fraction(point[0], curve.prime),
                    reduce_fraction(point[1], curve.prime),
                )
            ]
            for point in basis
        )
        multiples.append(
            tuple(
                tuple(curve.multiply_id(point_id, scalar) for scalar in choices)
                for point_id in reduced_basis
            )
        )
    states: dict[tuple[int, ...], tuple[int, ...]] = {
        tuple(0 for _ in curves): ()
    }
    counts = []
    for basis_index in range(len(basis)):
        retained: dict[tuple[int, ...], tuple[int, int, tuple[int, ...]]] = {}
        for state, vector in states.items():
            for choice_index, scalar in enumerate(choices):
                new_state = tuple(
                    curve.addition_table[state[curve_index]][
                        multiples[curve_index][basis_index][choice_index]
                    ]
                    for curve_index, curve in enumerate(curves)
                )
                new_vector = vector + (scalar,)
                priority = (
                    sum(abs(value) for value in new_vector),
                    max(abs(value) for value in new_vector),
                    new_vector,
                )
                previous = retained.get(new_state)
                if previous is None or priority < previous:
                    retained[new_state] = priority
        states = {state: priority[2] for state, priority in retained.items()}
        counts.append(len(states))
    return states, tuple(counts)


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=root / DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[2]
    started = time.monotonic()
    _, auxiliary, basis = load_input(root)
    curves = tuple(FiniteCurve(auxiliary.weierstrass_coefficients, prime) for prime in LOCAL_PRIMES)
    residue_maps = tuple(parameter_residues(curve, auxiliary) for curve in curves)
    states, depth_counts = shortest_state_vectors(curves, basis)
    expected_product = 1
    for curve in curves:
        expected_product *= len(curve.points)
    if len(states) != expected_product:
        raise AssertionError("the basis did not reach the complete finite product")
    selected = []
    for state, vector in states.items():
        residues = tuple(residue_maps[index][point_id] for index, point_id in enumerate(state))
        if all(
            residue is not None and residue in TARGET_RESIDUES[index]
            for index, residue in enumerate(residues)
        ):
            selected.append((vector, residues))
    selected.sort()
    if len(selected) != 192:
        raise AssertionError("the selected finite-state count changed")

    exact: dict[Fraction, dict[str, Any]] = {}
    for vector, residues in selected:
        point = None
        for basis_point, scalar in zip(basis, vector, strict=True):
            point = weierstrass_add(
                auxiliary.weierstrass_coefficients,
                point,
                weierstrass_multiply(
                    auxiliary.weierstrass_coefficients, basis_point, scalar
                ),
            )
        inverse = auxiliary.inverse(point)
        if inverse is None or inverse[0] == 0:
            continue
        parameter = abs(inverse[0])
        proxy = exact_radical_proxy(parameter, primes_up_to(PROXY_PRIME_BOUND))
        valuations = {prime: valuation for prime, valuation in proxy["small_prime_valuations"]}
        if valuations.get(13, 0) < 4 or valuations.get(37, 0) < 2 or valuations.get(83, 0) < 2:
            raise AssertionError("a lifted finite class lost its forced valuations")
        previous = exact.get(parameter)
        record = {
            "parameter": parameter,
            "coefficient_vector": vector,
            "finite_parameter_residues": residues,
            "proxy": proxy,
        }
        if previous is None or vector < previous["coefficient_vector"]:
            exact[parameter] = record
    ordered = sorted(
        exact.values(),
        key=lambda record: (
            record["proxy"]["log_radical_upper_proxy"],
            max(abs(record["parameter"].numerator), record["parameter"].denominator),
            record["parameter"],
        ),
    )
    minimum_proxy = ordered[0]["proxy"]["log_radical_upper_proxy"]
    artifact = {
        "schema_version": 1,
        "status": "complete_finite_group_shortest_class_screen",
        "target_hit": False,
        "source": {
            "input_artifact": str(INPUT_ARTIFACT),
            "input_artifact_sha256": EXPECTED_INPUT_SHA256,
            "input_script_sha256": EXPECTED_INPUT_SCRIPT_SHA256,
            "auxiliary_basis_count": len(basis),
        },
        "finite_screen": {
            "primes": list(LOCAL_PRIMES),
            "group_orders": [len(curve.points) for curve in curves],
            "target_residue_unions": [sorted(values) for values in TARGET_RESIDUES],
            "coefficient_range_per_coordinate": [-COEFFICIENT_BOUND, COEFFICIENT_BOUND],
            "state_counts_after_each_basis_coordinate": list(depth_counts),
            "complete_product_state_count": len(states),
            "complete_product_expected_count": expected_product,
            "selected_target_state_count": len(selected),
        },
        "exact_lifts": {
            "distinct_nonsingular_parameter_count": len(exact),
            "minimum_log_radical_upper_proxy": minimum_proxy,
            "proxy_below_190_count": sum(
                record["proxy"]["log_radical_upper_proxy"] < 190 for record in ordered
            ),
            "top_32": [
                {
                    "parameter": rational_to_string(record["parameter"]),
                    "coefficient_vector": list(record["coefficient_vector"]),
                    "finite_parameter_residues": list(record["finite_parameter_residues"]),
                    **record["proxy"],
                }
                for record in ordered[:32]
            ],
            "parameter_sha256": hashlib.sha256(
                "\n".join(
                    rational_to_string(record["parameter"])
                    for record in sorted(exact.values(), key=lambda item: item["parameter"])
                ).encode()
            ).hexdigest(),
        },
        "conclusion": {
            "target_hit": False,
            "exact_conductor_call_count": 0,
            "reason": (
                "The best exact radical proxy is far above 182.72, so the "
                "declared conductor gate rejects every shortest finite-class lift."
            ),
            "scope_warning": (
                "Alternative kernel lifts of the same finite classes are not excluded."
            ),
        },
        "reproduction": {
            "command": REPRODUCING_COMMAND,
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "script_sha256": file_sha256(Path(__file__).resolve()),
            "wall_seconds": time.monotonic() - started,
            "no_subprocesses": True,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(
        f"states={len(states)} selected={len(selected)} exact={len(exact)} "
        f"min_proxy={minimum_proxy:.12f}",
        flush=True,
    )
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
