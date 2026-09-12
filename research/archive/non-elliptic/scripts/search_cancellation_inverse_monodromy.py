#!/usr/bin/env python3
"""Numerical branch-cycle search for cancellation inverse polynomials.

For positive ``m,r`` we use the affine-normalized specialization

    F'(x) = (1 - lambda*x*(1-x)^m)^r,
    lambda = 4*(m+1).

It is obtained from the cancellation family by putting
``x=P*T/Q`` and ``lambda=Q^(m+1)/P``.  Hence it has the same generic
monodromy whenever the specialization avoids critical-point and
critical-value collisions.  Both avoidance conditions are checked exactly
over QQ before numerical continuation starts.

The continuation uses only SymPy and mpmath.  The resulting permutations are
analysed with SymPy's Schreier--Sims implementation and are also printed in
GAP syntax for an independent replay.

Examples
--------

    .venv/bin/python scripts/search_cancellation_inverse_monodromy.py 2 2
    .venv/bin/python scripts/search_cancellation_inverse_monodromy.py \
        --grid 1:5 1:4 --dps 70
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import mpmath as mp
import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup


X, V = sp.symbols("X V")


@dataclass(frozen=True)
class ExactModel:
    m: int
    r: int
    degree: int
    lam: int
    critical_polynomial: sp.Poly
    antiderivative: sp.Poly
    critical_value_polynomial: sp.Poly


@dataclass(frozen=True)
class SearchResult:
    m: int
    r: int
    degree: int
    generators: tuple[tuple[int, ...], ...]
    cycle_lengths: tuple[tuple[int, ...], ...]
    transitive: bool
    primitive: bool
    order: int
    classification: str


def exact_model(m: int, r: int) -> ExactModel:
    """Build a specialization and certify its branch values are distinct."""
    if m < 1 or r < 1:
        raise ValueError("m and r must be positive")
    lam = 4 * (m + 1)
    critical = sp.Poly(1 - lam * X * (1 - X) ** m, X, domain=sp.QQ)
    if sp.discriminant(critical.as_expr(), X) == 0:
        raise ArithmeticError("chosen specialization has colliding critical points")

    antiderivative_expr = sp.integrate(critical.as_expr() ** r, X)
    antiderivative = sp.Poly(antiderivative_expr, X, domain=sp.QQ)
    # Clearing coefficient denominators makes the resultant integral and
    # stable; ``fraction`` does not extract them from a QQ polynomial.
    denominator = sp.ilcm(
        *[term.q for term in antiderivative.all_coeffs()]
    )
    numerator = sp.expand(denominator * antiderivative.as_expr())
    branch_resultant = sp.resultant(
        critical.as_expr(), numerator - denominator * V, X
    )
    branch_values = sp.Poly(branch_resultant, V, domain=sp.QQ).sqf_part().monic()
    if branch_values.degree() != m + 1:
        raise ArithmeticError("chosen specialization has colliding critical values")
    if sp.discriminant(branch_values.as_expr(), V) == 0:
        raise ArithmeticError("critical-value polynomial is not separable")
    return ExactModel(
        m=m,
        r=r,
        degree=r * (m + 1) + 1,
        lam=lam,
        critical_polynomial=critical,
        antiderivative=antiderivative,
        critical_value_polynomial=branch_values,
    )


class Continuation:
    """Adaptive root continuation for one exact polynomial."""

    def __init__(self, model: ExactModel, dps: int) -> None:
        mp.mp.dps = dps
        self.model = model
        self.dps = dps
        self.coefficients = tuple(
            self._mpc(value) for value in model.antiderivative.all_coeffs()
        )
        n = model.degree
        self.derivative_coefficients = tuple(
            (n - index) * value
            for index, value in enumerate(self.coefficients[:-1])
        )
        self.residual_tolerance = mp.mpf(10) ** (-(dps * 2 // 3))
        self.collision_tolerance = mp.mpf(10) ** (-(dps // 3))

    @staticmethod
    def _mpc(value: object) -> mp.mpc:
        numeric = sp.N(value, mp.mp.dps)
        real, imag = numeric.as_real_imag()
        return mp.mpc(str(real), str(imag))

    @staticmethod
    def _horner(coefficients: Sequence[mp.mpc], value: mp.mpc) -> mp.mpc:
        answer = coefficients[0]
        for coefficient in coefficients[1:]:
            answer = answer * value + coefficient
        return answer

    def roots_at(self, target: mp.mpc) -> list[mp.mpc]:
        target_sympy = sp.Float(str(target.real), self.dps) + sp.I * sp.Float(
            str(target.imag), self.dps
        )
        roots = sp.nroots(
            self.model.antiderivative.as_expr() - target_sympy,
            n=self.dps,
            maxsteps=3000,
        )
        return [self._mpc(root) for root in roots]

    def _try_step(
        self, old_target: mp.mpc, new_target: mp.mpc, roots: Sequence[mp.mpc]
    ) -> list[mp.mpc] | None:
        delta = new_target - old_target
        output: list[mp.mpc] = []
        predictors: list[mp.mpc] = []
        scale = 1 + abs(new_target)
        for root in roots:
            derivative = self._horner(self.derivative_coefficients, root)
            if abs(derivative) <= self.residual_tolerance:
                return None
            candidate = root + delta / derivative
            predictor = candidate
            converged = False
            for _ in range(24):
                residual = self._horner(self.coefficients, candidate) - new_target
                if abs(residual) <= self.residual_tolerance * scale:
                    converged = True
                    break
                derivative = self._horner(
                    self.derivative_coefficients, candidate
                )
                if abs(derivative) <= self.residual_tolerance:
                    break
                candidate -= residual / derivative
            if not converged:
                return None
            predictors.append(predictor)
            output.append(candidate)

        if len(output) > 1:
            minimum_separation = min(
                abs(output[i] - output[j])
                for i in range(len(output))
                for j in range(i)
            )
            if minimum_separation <= self.collision_tolerance:
                return None
            # Residual convergence alone is insufficient: a large Newton
            # correction can land in a neighbouring root's basin and create
            # a spurious permutation.  First-order prediction has quadratic
            # error as the target step shrinks, so this test forces adaptive
            # subdivision until each lift is unambiguous relative to the
            # root separation at the new target.
            if any(
                abs(root - predictor) > minimum_separation / 4
                for root, predictor in zip(output, predictors)
            ):
                return None
        return output

    def step(
        self,
        old_target: mp.mpc,
        new_target: mp.mpc,
        roots: Sequence[mp.mpc],
        depth: int = 0,
    ) -> list[mp.mpc]:
        answer = self._try_step(old_target, new_target, roots)
        if answer is not None:
            return answer
        if depth >= 24:
            raise ArithmeticError(
                f"continuation failed near target {new_target} at depth {depth}"
            )
        midpoint = (old_target + new_target) / 2
        halfway = self.step(old_target, midpoint, roots, depth + 1)
        return self.step(midpoint, new_target, halfway, depth + 1)

    def polyline(
        self,
        points: Sequence[mp.mpc],
        roots: Sequence[mp.mpc],
    ) -> list[mp.mpc]:
        answer = list(roots)
        for old_target, new_target in zip(points, points[1:]):
            answer = self.step(old_target, new_target, answer)
        return answer


def _segment_points(start: mp.mpc, end: mp.mpc, count: int) -> list[mp.mpc]:
    return [start + (end - start) * index / count for index in range(count + 1)]


def _circle_points(
    center: mp.mpc,
    radius_vector: mp.mpc,
    turns: float,
    count: int,
) -> list[mp.mpc]:
    return [
        center
        + radius_vector * mp.exp(2j * mp.pi * turns * index / count)
        for index in range(count + 1)
    ]


def _distance_to_segment(point: mp.mpc, start: mp.mpc, end: mp.mpc) -> mp.mpf:
    vector = end - start
    parameter = mp.re((point - start) * mp.conj(vector)) / abs(vector) ** 2
    parameter = max(mp.mpf(0), min(mp.mpf(1), parameter))
    return abs(point - (start + parameter * vector))


def _connector(
    branch_values: Sequence[mp.mpc], index: int
) -> tuple[list[mp.mpc], mp.mpc, mp.mpf]:
    """Return a collision-free lollipop stem from a common base point."""
    center = sum(branch_values) / len(branch_values)
    spread = max(abs(value - center) for value in branch_values)
    global_separation = min(
        abs(branch_values[i] - branch_values[j])
        for i in range(len(branch_values))
        for j in range(i)
    )
    outer_radius = spread + max(mp.mpf(1), 10 * global_separation)
    base = center + outer_radius
    branch = branch_values[index]
    local_separation = min(
        abs(branch - other)
        for j, other in enumerate(branch_values)
        if j != index
    )
    loop_radius = local_separation * mp.mpf("0.16")
    required_clearance = min(
        global_separation * mp.mpf("0.06"), loop_radius * mp.mpf("0.3")
    )

    candidates: list[tuple[mp.mpf, mp.mpc, mp.mpc]] = []
    for sample in range(256):
        angle = 2 * mp.pi * sample / 256
        direction = mp.exp(1j * angle)
        displacement = branch - center
        projection = mp.re(displacement * mp.conj(direction))
        discriminant = (
            projection**2 + outer_radius**2 - abs(displacement) ** 2
        )
        if discriminant <= 0:
            continue
        distance = -projection + mp.sqrt(discriminant)
        endpoint = branch + loop_radius * direction
        outer = branch + distance * direction
        clearance = min(
            _distance_to_segment(other, endpoint, outer)
            for j, other in enumerate(branch_values)
            if j != index
        )
        if clearance > required_clearance:
            candidates.append((clearance, endpoint, outer))
    if not candidates:
        raise ArithmeticError(f"could not construct a safe loop around branch {index}")
    _, endpoint, outer = max(candidates, key=lambda item: item[0])

    outer_angle = mp.arg(outer - center)
    arc_count = max(12, math.ceil(abs(float(outer_angle)) * 24))
    arc = [
        center
        + outer_radius * mp.exp(1j * outer_angle * step / arc_count)
        for step in range(arc_count + 1)
    ]
    stem = arc + _segment_points(outer, endpoint, 64)[1:]
    assert abs(stem[0] - base) < mp.mpf(10) ** (-(mp.mp.dps // 2))
    return stem, endpoint - branch, loop_radius


def _permutation_from_endpoint(
    initial: Sequence[mp.mpc], final: Sequence[mp.mpc], dps: int
) -> tuple[int, ...]:
    tolerance = mp.mpf(10) ** (-(dps // 3))
    permutation: list[int] = []
    for root in final:
        distances = [abs(root - candidate) for candidate in initial]
        closest = min(range(len(initial)), key=distances.__getitem__)
        if distances[closest] > tolerance:
            raise ArithmeticError(
                f"endpoint root-matching error {mp.nstr(distances[closest], 8)}"
            )
        permutation.append(closest)
    if len(set(permutation)) != len(initial):
        raise ArithmeticError("endpoint matching did not produce a permutation")
    return tuple(permutation)


def _cycle_lengths(permutation: Sequence[int]) -> tuple[int, ...]:
    seen: set[int] = set()
    lengths: list[int] = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        current = start
        length = 0
        while current not in seen:
            seen.add(current)
            length += 1
            current = permutation[current]
        if length > 1:
            lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def search_case(m: int, r: int, dps: int = 70) -> SearchResult:
    model = exact_model(m, r)
    continuation = Continuation(model, dps)
    branch_values = sorted(
        (
            continuation._mpc(root)
            for root in sp.nroots(
                model.critical_value_polynomial.as_expr(),
                n=dps,
                maxsteps=3000,
            )
        ),
        key=lambda value: (float(value.real), float(value.imag)),
    )
    center = sum(branch_values) / len(branch_values)
    spread = max(abs(value - center) for value in branch_values)
    separation = min(
        abs(branch_values[i] - branch_values[j])
        for i in range(len(branch_values))
        for j in range(i)
    )
    base = center + spread + max(mp.mpf(1), 10 * separation)
    initial = continuation.roots_at(base)

    generators: list[tuple[int, ...]] = []
    loop_count = max(160, 48 * (r + 1))
    for index, branch in enumerate(branch_values):
        stem, radius_vector, _ = _connector(branch_values, index)
        if abs(stem[0] - base) > mp.mpf(10) ** (-(dps // 2)):
            raise AssertionError("connectors do not share their base point")
        roots = continuation.polyline(stem, initial)
        loop = _circle_points(branch, radius_vector, 1.0, loop_count)
        roots = continuation.polyline(loop, roots)
        roots = continuation.polyline(list(reversed(stem)), roots)
        generators.append(_permutation_from_endpoint(initial, roots, dps))

    cycle_lengths = tuple(_cycle_lengths(generator) for generator in generators)
    expected = (r + 1,)
    if any(lengths != expected for lengths in cycle_lengths):
        raise ArithmeticError(
            f"unexpected local cycle types {cycle_lengths}; expected {expected}"
        )

    group = PermutationGroup(*(Permutation(generator) for generator in generators))
    order = int(group.order())
    alternating_order = math.factorial(model.degree) // 2
    symmetric_order = math.factorial(model.degree)
    if order == symmetric_order:
        classification = f"S_{model.degree}"
    elif order == alternating_order:
        classification = f"A_{model.degree}"
    else:
        classification = "proper exceptional"
    return SearchResult(
        m=m,
        r=r,
        degree=model.degree,
        generators=tuple(generators),
        cycle_lengths=cycle_lengths,
        transitive=bool(group.is_transitive()),
        primitive=bool(group.is_primitive()),
        order=order,
        classification=classification,
    )


def _gap_permutation(permutation: Sequence[int]) -> str:
    cycles: list[str] = []
    seen: set[int] = set()
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle: list[int] = []
        current = start
        while current not in seen:
            seen.add(current)
            cycle.append(current + 1)
            current = permutation[current]
        if len(cycle) > 1:
            cycles.append("(" + ",".join(map(str, cycle)) + ")")
    return "".join(cycles) or "()"


def _parse_range(text: str) -> range:
    if ":" not in text:
        value = int(text)
        return range(value, value + 1)
    start_text, end_text = text.split(":", 1)
    start, end = int(start_text), int(end_text)
    if start < 1 or end < start:
        raise argparse.ArgumentTypeError("ranges must have form positive_start:end")
    return range(start, end + 1)


def _result_dict(result: SearchResult) -> dict[str, object]:
    return {
        "m": result.m,
        "r": result.r,
        "degree": result.degree,
        "cycle_lengths": result.cycle_lengths,
        "generators_zero_based": result.generators,
        "gap_generators": [_gap_permutation(p) for p in result.generators],
        "transitive": result.transitive,
        "primitive": result.primitive,
        "order": str(result.order),
        "classification": result.classification,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("m", nargs="?", type=int)
    parser.add_argument("r", nargs="?", type=int)
    parser.add_argument(
        "--grid",
        nargs=2,
        metavar=("M_RANGE", "R_RANGE"),
        help="inclusive ranges such as 1:5 1:4",
    )
    parser.add_argument("--dps", type=int, default=70)
    parser.add_argument("--json", type=Path, help="write all tuples to this file")
    args = parser.parse_args()
    if args.dps < 40:
        parser.error("--dps must be at least 40")
    if args.grid:
        cases: Iterable[tuple[int, int]] = (
            (m, r)
            for m in _parse_range(args.grid[0])
            for r in _parse_range(args.grid[1])
        )
    elif args.m is not None and args.r is not None:
        cases = [(args.m, args.r)]
    else:
        parser.error("supply m r or --grid M_RANGE R_RANGE")

    records: list[dict[str, object]] = []
    for m, r in cases:
        result = search_case(m, r, args.dps)
        record = _result_dict(result)
        records.append(record)
        parity_prediction = "A" if r % 2 == 0 else "S"
        print(
            f"(m,r)=({m},{r}) N={result.degree}: {result.classification}; "
            f"transitive={result.transitive}, primitive={result.primitive}, "
            f"predicted={parity_prediction}_{result.degree}"
        )
        print("  GAP generators := [" + ", ".join(record["gap_generators"]) + "];")
        if args.json:
            # Preserve completed cases if a later, harder continuation needs
            # more precision or exposes a genuine exceptional passport.
            args.json.write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
