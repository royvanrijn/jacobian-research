#!/usr/bin/env python3
"""Nagao sieve on bounded points of the rank-four base curve E0.

The four certified Mordell--Weil generators are reduced once at every sieve
prime.  Their finite subgroup image is enumerated in at most ``#E0(F_p)``
steps, the exact map ``E0 -> t`` is evaluated on that image, and Nagao
contributions are cached by reduced E0 point.  The rational lattice box then
uses only four tabled finite-group additions per prime and coefficient vector.

This is a bounded heuristic sieve.  It does not assert that the four displayed
generators saturate E0(Q), nor does a large score prove a rank jump.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import importlib.util
from itertools import product
import json
from pathlib import Path
import shlex
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC_ROOT), str(ELLIPTIC_ROOT / "cas")]

from ecsearch.q12o5867_specialization import short_certificate_model  # noqa: E402
from elliptic_candidate_record import (  # noqa: E402
    _finite_add,
    _finite_curve_points,
    _finite_multiply,
    _reduce_rational,
    source_point_to_target,
)


COMMON_PATH = ROOT / "elkies-k3/scripts/search_h92_q12o5867_rootless_nagao.py"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
PAIRED = ROOT / "elkies-k3/data/fibrations/elkies_2026_rank19_paired_cover.json"
Q = Fraction


def fraction_text(value):
    value = Q(value)
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def e0_add(left, right):
    """Exact group law on y^2=x^3+1029367969*x^2-42900734074705920*x."""

    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = left
    x2, y2 = right
    if x1 == x2:
        if y1 == -y2 or y1 == 0:
            return None
        slope = (3 * x1**2 + 2 * Q(1029367969) * x1 - Q(42900734074705920)) / (2 * y1)
    else:
        slope = (y2 - y1) / (x2 - x1)
    x3 = slope**2 - Q(1029367969) - x1 - x2
    y3 = -y1 + slope * (x1 - x3)
    return x3, y3


def e0_multiply(point, scalar):
    if scalar < 0:
        return e0_multiply((point[0], -point[1]), -scalar)
    answer = None
    addend = point
    while scalar:
        if scalar & 1:
            answer = e0_add(answer, addend)
        addend = e0_add(addend, addend)
        scalar >>= 1
    return answer


def exact_e0_record(vector, generators):
    point = None
    for coefficient, generator in zip(vector, generators):
        point = e0_add(point, e0_multiply(generator, coefficient))
    if point is None:
        return {"E0_point": "infinity", "base_t": "infinity"}
    x, y = point
    denominator_r = Q(52185783681) * x + Q(1059345) * y - Q(6263787913107172000)
    numerator_r = (
        Q(353115) * x**2
        - Q(11678422142720) * x
        + Q(5060365752) * y
        - Q(11169586606652709110400)
    )
    if denominator_r == 0:
        r = None
        t = None
    else:
        r = numerator_r / (6 * denominator_r)
        t_denominator = Q(130) * r - Q(38636)
        t = None if t_denominator == 0 else (Q(289444) - r**2) / t_denominator
    return {
        "E0_point": [fraction_text(x), fraction_text(y)],
        "quartic_r": "infinity" if r is None else fraction_text(r),
        "base_t": "infinity" if t is None else fraction_text(t),
    }


def load_common():
    spec = importlib.util.spec_from_file_location("elkies_r17_nagao_common", COMMON_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {COMMON_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def subgroup_image(generators, coefficient_a, prime):
    image = {None}
    for generator in generators:
        cyclic = []
        point = None
        while point not in cyclic:
            cyclic.append(point)
            point = _finite_add(point, generator, coefficient_a, prime)
        if point is not None:
            raise AssertionError("a finite cyclic subgroup failed to return to zero")
        image = {
            _finite_add(left, right, coefficient_a, prime)
            for left in image
            for right in cyclic
        }
    return image


def source_coordinates_mod(point, change, prime):
    if point is None:
        return None
    x_target, y_target = point
    u = _reduce_rational(change.u, prime)
    r = _reduce_rational(change.r, prime)
    s = _reduce_rational(change.s, prime)
    shift_t = _reduce_rational(change.t, prime)
    x_source = (u**2 * x_target + r) % prime
    y_source = (u**3 * y_target + s * u**2 * x_target + shift_t) % prime
    return x_source, y_source


def t_index_from_e0_point(point, change, prime):
    if point is None:
        return prime
    x, y = source_coordinates_mod(point, change, prime)
    denominator_r = (52185783681 * x + 1059345 * y - 6263787913107172000) % prime
    numerator_r = (
        353115 * x**2
        - 11678422142720 * x
        + 5060365752 * y
        - 11169586606652709110400
    ) % prime
    r_numerator = numerator_r
    r_denominator = 6 * denominator_r % prime
    t_numerator = (289444 * r_denominator**2 - r_numerator**2) % prime
    t_denominator = (
        130 * r_numerator * r_denominator - 38636 * r_denominator**2
    ) % prime
    if t_numerator == 0 and t_denominator == 0:
        return None
    if t_denominator == 0:
        return prime
    return t_numerator * pow(t_denominator, -1, prime) % prime


def main() -> None:
    common = load_common()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--paired-cover", type=Path, default=PAIRED)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--coefficient-bound", type=int, default=6)
    parser.add_argument("--primes", default="11-97")
    parser.add_argument("--finalists", type=int, default=1000)
    args = parser.parse_args()
    if args.coefficient_bound < 0 or args.finalists < 1:
        parser.error("bounds must be nonnegative and --finalists positive")

    started = perf_counter()
    paired = json.loads(args.paired_cover.read_text())
    if paired.get("status") != "PASS_EXACT_ELKIES_2026_PAIRED_COVER_DATA":
        raise ValueError("the paired-cover exact data status is missing")
    model = common.load_family_model(args.model)
    primes = common.parse_prime_blocks(args.primes)[0]
    e0_model = tuple(Q(value) for value in paired["E0"]["weierstrass_coefficients_a1_a2_a3_a4_a6"])
    e0_generators = tuple(
        (Q(point[0]), Q(point[1])) for point in paired["E0"]["generators"]
    )
    short_model, short_change = short_certificate_model(e0_model)
    short_generators = tuple(
        source_point_to_target(point, short_change) for point in e0_generators
    )

    local_records = []
    local_data = []
    rejected = []
    coefficient_values = tuple(
        range(-args.coefficient_bound, args.coefficient_bound + 1)
    )
    for prime in primes:
        try:
            coefficient_a = _reduce_rational(short_model[3], prime)
            coefficient_b = _reduce_rational(short_model[4], prime)
            if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
                raise ValueError("E0 has bad reduction")
            reduced_generators = tuple(
                (
                    _reduce_rational(point[0], prime),
                    _reduce_rational(point[1], prime),
                )
                for point in short_generators
            )
            if any(
                (point[1] ** 2 - point[0] ** 3 - coefficient_a * point[0] - coefficient_b)
                % prime
                for point in reduced_generators
            ):
                raise ValueError("an E0 generator failed reduction")
            base_table = common.residue_table(model, prime)
            image = subgroup_image(reduced_generators, coefficient_a, prime)
            point_symbols = {}
            undefined = 0
            for point in image:
                index = t_index_from_e0_point(point, short_change, prime)
                if index is None:
                    point_symbols[point] = None
                    undefined += 1
                else:
                    point_symbols[point] = base_table[index]
            multiples = tuple(
                {
                    coefficient: _finite_multiply(
                        generator, coefficient, coefficient_a, prime
                    )
                    for coefficient in coefficient_values
                }
                for generator in reduced_generators
            )
        except (ValueError, ZeroDivisionError) as error:
            rejected.append({"prime": prime, "reason": str(error)})
            continue
        group_order = len(_finite_curve_points(coefficient_a, coefficient_b, prime))
        local_records.append(
            {
                "prime": prime,
                "E0_group_order": group_order,
                "mw_reduction_image_size": len(image),
                "mw_reduction_index": group_order // len(image),
                "undefined_map_points": undefined,
                "reduced_generators": [
                    None if point is None else list(point) for point in reduced_generators
                ],
            }
        )
        local_data.append((prime, coefficient_a, multiples, point_symbols))
    if not local_data:
        raise ValueError("no usable E0 reduction primes")

    candidates = []
    for vector in product(coefficient_values, repeat=4):
        if vector == (0, 0, 0, 0):
            continue
        score_units = 0
        good = 0
        bad = 0
        undefined = 0
        for _prime, coefficient_a, multiples, point_symbols in local_data:
            point = None
            for index, coefficient in enumerate(vector):
                point = _finite_add(
                    point, multiples[index][coefficient], coefficient_a, _prime
                )
            symbol = point_symbols[point]
            if symbol is None:
                undefined += 1
            elif symbol.good_reduction:
                score_units += symbol.contribution_units
                good += 1
            else:
                bad += 1
        candidates.append(
            {
                "mw_coefficients": list(vector),
                "coefficient_linf": max(abs(value) for value in vector),
                "coefficient_l1": sum(abs(value) for value in vector),
                "score_units_1e12": score_units,
                "score": score_units / common.SCORE_SCALE,
                "good_prime_count": good,
                "bad_reduction_prime_count": bad,
                "undefined_map_prime_count": undefined,
            }
        )
    candidates.sort(
        key=lambda record: (
            -record["score_units_1e12"],
            record["coefficient_linf"],
            record["coefficient_l1"],
            record["mw_coefficients"],
        )
    )
    finalists = candidates[: args.finalists]
    for record in finalists:
        record.update(exact_e0_record(record["mw_coefficients"], e0_generators))
    payload = {
        "schema": "elkies-k3.elkies-2026-E0-mw-nagao-sieve.v1",
        "status": "PASS_BOUNDED_HEURISTIC_E0_MW_NAGAO_SIEVE",
        "E0": {
            "source": str(args.paired_cover.resolve()),
            "rank_lower_bound": 4,
            "generator_count": 4,
            "coefficient_box": [-args.coefficient_bound, args.coefficient_bound],
        },
        "search": {
            "requested_primes": list(primes),
            "usable_prime_count": len(local_data),
            "rejected_primes": rejected,
            "enumerated_nonzero_lattice_vectors": len(candidates),
            "method": (
                "enumerate each four-generator MW reduction image once, cache E0->t "
                "Nagao symbols, then scan the bounded coefficient lattice"
            ),
        },
        "local_mw_reductions": local_records,
        "finalists": finalists,
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "The E0 generators, their independence, and E0->t map are exact. The "
            "coefficient box and Nagao ranking are bounded heuristics; candidates need "
            "exact rational evaluation and specialization certificates."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"PASS E0_mw vectors={len(candidates)} primes={len(local_data)} "
        f"finalists={len(finalists)} seconds={payload['runtime_seconds']:.3f} "
        f"output={args.output}"
    )


if __name__ == "__main__":
    main()
