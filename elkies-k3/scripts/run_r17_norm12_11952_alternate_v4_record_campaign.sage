#!/usr/bin/env sage-python
"""Run the native alternate-Q80 record campaign on rank-one V4 bases.

status: ACTIVE_SEARCH
claim: bounded periodic base sieve, exact all-bisection split tests, and rank lower bounds
inputs: seventeen prepared rank-one V4 bases, direct alternate model/sections,
        and all 39,147 native alternate bisections
outputs: artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-record-campaign-v1.json

For each base and each good prime, the complete period of ``u(nP)`` is
precomputed.  Integers n are ranked by centered, population-standardized
surface-trace and bisection-splitting signals on three disjoint prime blocks.
Only the declared best n are formed over QQ.  At each such u, all native
bisections pass through a no-false-negative modular square sieve before exact
rational-square tests.  The seventeen generic sections, the two defining
character sections, and every additional split section are then materialized
and certified through exact finite quotients.

The scoring is heuristic.  Rank is promoted only from a stored independence
certificate.  Residual descent is required only for a rank-32 candidate; with
the current generic rank-19 subgroup, residual 2-Selmer dimension below 13 is
an exact rejection gate before any unrestricted point search.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import importlib.util
import json
from math import log, sqrt
from pathlib import Path
import shlex
import sys
from time import perf_counter

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
BASES = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-rank-one-bases-v1.json"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-full-v1.json"
FINITE_HELPER = ROOT / "elliptic-curves/cas/elliptic_candidate_record.py"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-record-campaign-v1.json"
SCHEMA = "elkies-k3.r17-norm12-11952-v4-record-campaign.v1"
DEFAULT_BLOCKS = ((131, 137, 151), (157, 167, 173), (181, 191, 193))
SCORE_SCALE = 10**6


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def load_module(name: str, path: Path):
    specification = importlib.util.spec_from_loader(name, SourceFileLoader(name, str(path)))
    if specification is None or specification.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def point_text(point) -> list[str]:
    if point.is_zero():
        return ["0"]
    return [rational_text(point[0]), rational_text(point[1])]


def projective_bits(value) -> int:
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


def evaluate(coefficients, value):
    answer = value.parent()(0) if hasattr(value, "parent") else QQ(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + value.parent()(QQ(coefficient)) if hasattr(value, "parent") else answer * value + QQ(coefficient)
    return answer


def evaluate_qq(coefficients, value):
    answer = QQ(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + QQ(coefficient)
    return answer


def normalized_homogeneous_quadratic(record):
    numerator = [QQ(value) for value in record["branch"]["numerator_coefficients"]]
    denominator_polynomial = [QQ(value) for value in record["branch"]["denominator_coefficients"]]
    if len(numerator) != 3 or len(denominator_polynomial) != 1 or denominator_polynomial[0] == 0:
        raise ArithmeticError(f"{record['label']}: branch is not a finite quadratic")
    coefficients = [value / denominator_polynomial[0] for value in numerator]
    denominator = ZZ(1)
    for value in coefficients:
        denominator = denominator.lcm(value.denominator())
    integers = tuple(ZZ(value * denominator) for value in coefficients)
    content = denominator
    for value in integers:
        content = content.gcd(value)
    if content > 1:
        integers = tuple(value // content for value in integers)
        denominator //= content
    return denominator, integers


def exact_square_value(quadratic, numerator, denominator):
    scale, (c0, c1, c2) = quadratic
    homogeneous = c0 * denominator**2 + c1 * numerator * denominator + c2 * numerator**2
    square_integer = scale * homogeneous
    if square_integer < 0 or not square_integer.is_square():
        return None
    root = QQ(square_integer.sqrt()) / (scale * denominator)
    value = QQ(homogeneous) / (scale * denominator**2)
    if root**2 != value:
        raise ArithmeticError("integer square reconstruction failed")
    return value, root


def quadratic_residue_flags(prime):
    flags = [False] * prime
    flags[0] = True
    for value in range(1, prime):
        flags[value * value % prime] = True
    return flags


def split_count_table(quadratics, prime):
    residues = quadratic_residue_flags(prime)
    counts = [0] * (prime + 1)
    usable = 0
    for scale, coefficients in quadratics:
        scale_mod = int(scale % prime)
        if scale_mod == 0:
            continue
        c0, c1, c2 = (int(value % prime) for value in coefficients)
        usable += 1
        for value in range(prime):
            if residues[scale_mod * (c0 + c1 * value + c2 * value * value) % prime]:
                counts[value] += 1
        if residues[scale_mod * c2 % prime]:
            counts[prime] += 1
    return counts, usable


def coefficient_mod(value, prime):
    value = QQ(value)
    denominator = int(value.denominator() % prime)
    if denominator == 0:
        raise ZeroDivisionError
    return int(value.numerator() % prime) * pow(denominator, -1, prime) % prime


def homogeneous_value(coefficients, degree, numerator, denominator, prime):
    return sum(
        coefficient * pow(numerator, index, prime) * pow(denominator, degree - index, prime)
        for index, coefficient in enumerate(coefficients)
    ) % prime


def surface_trace_table(direct, prime):
    model = direct["weierstrass_model"]
    a_coefficients = [coefficient_mod(value, prime) for value in model["A_coefficients_low_to_high"]]
    b_coefficients = [coefficient_mod(value, prime) for value in model["B_coefficients_low_to_high"]]
    characters = [-1] * prime
    characters[0] = 0
    for value in range(1, prime):
        characters[value * value % prime] = 1
    table = []
    for index in range(prime + 1):
        numerator, denominator = (index, 1) if index < prime else (1, 0)
        coefficient_a = homogeneous_value(a_coefficients, 8, numerator, denominator, prime)
        coefficient_b = homogeneous_value(b_coefficients, 12, numerator, denominator, prime)
        if (4 * coefficient_a**3 + 27 * coefficient_b**2) % prime == 0:
            table.append(None)
            continue
        trace = -sum(characters[(x**3 + coefficient_a * x + coefficient_b) % prime] for x in range(prime))
        point_count = prime + 1 - trace
        table.append((2.0 - trace) / point_count * log(float(prime)))
    return table


def reduce_rational(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


def periodic_u_values(base, prime):
    field = GF(prime)
    paired = base["paired_v4_base"]
    try:
        curve = EllipticCurve(field, [reduce_rational(value, field) for value in paired["pointed_a1_a2_a3_a4_a6"]])
        generator_data = paired["primitive_generator_pointed"]
        generator = curve(reduce_rational(generator_data[0], field), reduce_rational(generator_data[1], field))
        data = base["map_to_v4_cover"]
        origin = data["origin"]
        u0 = reduce_rational(origin["u0"], field)
        s0 = reduce_rational(origin["s0"], field)
        q_left = [reduce_rational(value, field) for value in data["left_q_coefficients_low_to_high"]]
        c = reduce_rational(data["pointed_inverse_constants"]["c"], field)
        d = reduce_rational(data["pointed_inverse_constants"]["d"], field)
        v0 = reduce_rational(data["pointed_inverse_constants"]["v0"], field)
    except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
        return None
    if curve.discriminant() == 0 or generator.is_zero():
        return None
    order = int(generator.order())
    values = [None] * order
    point = curve(0)
    for index in range(order):
        if not point.is_zero() and point[1] != 0 and v0 != 0:
            numerator = 4 * v0**2 * (point[0] + c) - d**2
            denominator = 2 * v0 * point[1]
            if denominator != 0:
                slope = numerator / denominator
                conic_denominator = 1 - q_left[2] * slope**2
                if conic_denominator != 0:
                    u_value = u0 + (
                        (q_left[1] + 2 * q_left[2] * u0) * slope**2 - 2 * s0 * slope
                    ) / conic_denominator
                    values[index] = int(u_value)
        point += generator
    return order, values


def standardize(values):
    finite = [float(value) for value in values if value is not None]
    if not finite:
        return [0] * len(values), None
    mean = sum(finite) / len(finite)
    variance = sum((value - mean) ** 2 for value in finite) / len(finite)
    deviation = sqrt(variance)
    if deviation == 0:
        standardized = [0 if value is not None else 0 for value in values]
    else:
        standardized = [
            0 if value is None else int(round((float(value) - mean) / deviation * SCORE_SCALE))
            for value in values
        ]
    return standardized, {"mean": mean, "population_standard_deviation": deviation, "valid_period_positions": len(finite)}


def build_base_periods(base, blocks, local_tables):
    periods = {}
    rejected = []
    for block_index, block in enumerate(blocks):
        for prime in block:
            reduced = periodic_u_values(base, prime)
            if reduced is None:
                rejected.append({"prime": prime, "block": block_index + 1, "reason": "bad_base_reduction_or_generator"})
                continue
            order, u_values = reduced
            traces = [None if value is None else local_tables[prime]["trace"][value] for value in u_values]
            splits = [None if value is None else local_tables[prime]["split"][value] for value in u_values]
            trace_z, trace_stats = standardize(traces)
            split_z, split_stats = standardize(splits)
            periods[prime] = {
                "block": block_index + 1,
                "generator_order": order,
                "u_projective_indices": u_values,
                "trace_z_units": trace_z,
                "split_z_units": split_z,
                "trace_statistics": trace_stats,
                "split_statistics": split_stats,
            }
    return periods, rejected


def rank_integers(periods, blocks, n_bound, minimum_abs_n):
    candidates = []
    for n_value in range(-n_bound, n_bound + 1):
        if abs(n_value) < minimum_abs_n:
            continue
        trace_blocks = []
        split_blocks = []
        usable_blocks = []
        for block in blocks:
            trace_score = 0
            split_score = 0
            usable = 0
            for prime in block:
                period = periods.get(prime)
                if period is None:
                    continue
                index = n_value % period["generator_order"]
                if period["u_projective_indices"][index] is None:
                    continue
                trace_score += period["trace_z_units"][index]
                split_score += period["split_z_units"][index]
                usable += 1
            trace_blocks.append(trace_score)
            split_blocks.append(split_score)
            usable_blocks.append(usable)
        if any(value == 0 for value in usable_blocks):
            continue
        combined = [left + right for left, right in zip(trace_blocks, split_blocks)]
        candidates.append(
            {
                "n": n_value,
                "trace_block_score_units_1e6": trace_blocks,
                "split_block_score_units_1e6": split_blocks,
                "combined_block_score_units_1e6": combined,
                "usable_primes_per_block": usable_blocks,
                "ranking_key": [min(combined), min(split_blocks), min(trace_blocks), sum(combined)],
            }
        )
    candidates.sort(key=lambda row: tuple(-int(value) for value in row["ranking_key"]) + (abs(row["n"]), row["n"]))
    for rank, row in enumerate(candidates, start=1):
        row["period_score_rank"] = rank
    return candidates


def map_point_exact(base, n_value):
    paired = base["paired_v4_base"]
    curve = EllipticCurve(QQ, [QQ(value) for value in paired["pointed_a1_a2_a3_a4_a6"]])
    generator = curve(QQ(paired["primitive_generator_pointed"][0]), QQ(paired["primitive_generator_pointed"][1]))
    point = ZZ(n_value) * generator
    if point.is_zero() or point[1] == 0:
        raise ArithmeticError("selected multiple lies on the pointed-map exceptional locus")
    data = base["map_to_v4_cover"]
    origin = data["origin"]
    u0, s0, t0 = (QQ(origin[key]) for key in ("u0", "s0", "t0"))
    q_left = [QQ(value) for value in data["left_q_coefficients_low_to_high"]]
    q_right = [QQ(value) for value in data["right_q_coefficients_low_to_high"]]
    c = QQ(data["pointed_inverse_constants"]["c"])
    d = QQ(data["pointed_inverse_constants"]["d"])
    v0 = QQ(data["pointed_inverse_constants"]["v0"])
    slope = (4 * v0**2 * (point[0] + c) - d**2) / (2 * v0 * point[1])
    if slope == 0:
        u_value, s_value, t_value = u0, s0, -t0
    else:
        paired_ordinate = (point[0] * slope**2 - d * slope) / (2 * v0) - v0
        denominator = 1 - q_left[2] * slope**2
        if denominator == 0:
            raise ArithmeticError("selected multiple maps to the conic infinity chart")
        derivative = q_left[1] + 2 * q_left[2] * u0
        u_value = u0 + (derivative * slope**2 - 2 * s0 * slope) / denominator
        s_value = s0 + (u_value - u0) / slope
        t_value = paired_ordinate / denominator
    if s_value**2 != evaluate_qq(q_left, u_value) or t_value**2 != evaluate_qq(q_right, u_value):
        raise ArithmeticError("exact multiple misses the V4 cover")
    return point, u_value, s_value, t_value


def modular_survivors(quadratics, numerator, denominator, primes):
    survivors = list(range(len(quadratics)))
    counts = []
    for prime in primes:
        residues = quadratic_residue_flags(prime)
        next_survivors = []
        if denominator % prime:
            value = int(numerator % prime) * pow(int(denominator % prime), -1, prime) % prime
            infinity = False
        else:
            value = prime
            infinity = True
        for index in survivors:
            scale, coefficients = quadratics[index]
            scale_mod = int(scale % prime)
            if scale_mod == 0:
                next_survivors.append(index)
                continue
            c0, c1, c2 = (int(coefficient % prime) for coefficient in coefficients)
            evaluation = c2 if infinity else c0 + c1 * value + c2 * value * value
            if residues[scale_mod * evaluation % prime]:
                next_survivors.append(index)
        survivors = next_survivors
        counts.append(len(survivors))
        if not survivors:
            break
    return survivors, counts


def rational_function_value(record, coordinate, u_value):
    data = record[coordinate]
    numerator = evaluate_qq(data["numerator_coefficients_low_to_high"], u_value)
    denominator = evaluate_qq(data["denominator_coefficients_low_to_high"], u_value)
    if denominator == 0:
        raise ZeroDivisionError("section specializes at a pole")
    return numerator / denominator


def generic_sections(direct, u_value, curve):
    points = []
    for expected, record in enumerate(direct["sections"]["records"]):
        if int(record["basis_index"]) != expected:
            raise ArithmeticError("alternate section order changed")
        points.append(curve(rational_function_value(record, "X", u_value), rational_function_value(record, "Y", u_value)))
    if len(points) != 17:
        raise ArithmeticError("expected seventeen generic alternate sections")
    return points


def lifted_point(record, u_value, square_root, curve):
    lifted = record["lifted_section"]
    x_value = evaluate_qq(lifted["x0_coefficients"], u_value) + evaluate_qq(lifted["x1_coefficients"], u_value) * square_root
    y_value = evaluate_qq(lifted["y0_coefficients"], u_value) + evaluate_qq(lifted["y1_coefficients"], u_value) * square_root
    return curve(x_value, y_value)


def to_fraction(value):
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def finite_certificate(helper, model, points, prime_bound):
    fraction_points = tuple((to_fraction(point[0]), to_fraction(point[1])) for point in points)
    if any(not helper.is_on_weierstrass_curve(model, point) for point in fraction_points):
        raise ArithmeticError("a displayed point misses the specialized alternate curve")
    original_check = helper.is_on_weierstrass_curve
    helper.is_on_weierstrass_curve = lambda _model, _point: True
    try:
        attempts = []
        for relation_prime in (2, 3):
            certificate = helper.build_finite_quotient_certificate(
                model, fraction_points, relation_prime=relation_prime, prime_bound=prime_bound
            )
            attempts.append(certificate)
            if certificate["certified_independent"]:
                helper.verify_finite_quotient_certificate(model, fraction_points, certificate)
                return {
                    "certified_rank_lower_bound": len(points),
                    "all_displayed_points_independent": True,
                    "successful_certificate": certificate,
                    "attempt_ranks": [int(row["combined_rank_over_relation_field"]) for row in attempts],
                }
        best = max(attempts, key=lambda row: int(row["combined_rank_over_relation_field"]))
        return {
            "certified_rank_lower_bound": int(best["combined_rank_over_relation_field"]),
            "all_displayed_points_independent": False,
            "successful_certificate": None,
            "attempt_ranks": [int(row["combined_rank_over_relation_field"]) for row in attempts],
        }
    finally:
        helper.is_on_weierstrass_curve = original_check


def promotion(rank_lower_bound):
    if rank_lower_bound >= 32:
        return "RECORD_CANDIDATE"
    if rank_lower_bound >= 28:
        return "SERIOUS_ALTERNATE_Q80_POSITIVE_CONTROL"
    if rank_lower_bound >= 24:
        return "USEFUL_SCORE_CALIBRATION"
    if rank_lower_bound >= 20:
        return "NATIVE_MECHANISM_CONTROL"
    return "BELOW_PROMOTION_THRESHOLD"


def parse_blocks(text):
    blocks = tuple(tuple(int(value) for value in block.split(",") if value) for block in text.split(";"))
    if len(blocks) < 3 or any(not block for block in blocks):
        raise ValueError("at least three nonempty disjoint prime blocks are required")
    flattened = [prime for block in blocks for prime in block]
    if len(set(flattened)) != len(flattened) or any(prime < 5 or not ZZ(prime).is_prime() for prime in flattened):
        raise ValueError("sieve moduli must be distinct primes at least five")
    return blocks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bases", type=Path, default=BASES)
    parser.add_argument("--direct", type=Path, default=DIRECT)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--n-bound", type=int, default=4)
    parser.add_argument("--minimum-abs-n", type=int, default=2)
    parser.add_argument("--exact-per-base", type=int, default=1)
    parser.add_argument("--base-limit", type=int, default=17)
    parser.add_argument("--prime-blocks", default=";".join(",".join(map(str, block)) for block in DEFAULT_BLOCKS))
    parser.add_argument("--certificate-prime-bound", type=int, default=300)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        document = json.loads(args.output.read_text())
        if document.get("schema") != SCHEMA:
            raise ValueError("unexpected alternate V4 campaign artifact")
        for name, expected in document["inputs"].items():
            if digest(ROOT / name) != expected:
                raise ArithmeticError(f"campaign input changed: {name}")
        print(f"ALTV4CAMPAIGNCHECK|fibres={len(document['fibres'])}|output={display_path(args.output)}")
        return
    if args.n_bound < 1 or args.minimum_abs_n < 1 or args.minimum_abs_n > args.n_bound or args.exact_per_base < 1 or not 1 <= args.base_limit <= 17:
        parser.error("invalid n bounds or exact-per-base")
    blocks = parse_blocks(args.prime_blocks)
    all_primes = tuple(prime for block in blocks for prime in block)
    started = perf_counter()

    bases_document = json.loads(args.bases.read_text())
    direct = json.loads(args.direct.read_text())
    bisection_document = json.loads(args.bisections.read_text())
    if bases_document.get("status") != "PASS_EXACT_RANK_ONE_V4_BASE_GENERATORS_AND_MAPS" or len(bases_document["bases"]) != 17:
        raise ValueError("the seventeen prepared rank-one bases are required")
    if direct.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ValueError("unexpected direct alternate model")
    records = sorted(bisection_document["bisections"], key=lambda row: int(row["priority_rank"]))
    if len(records) != 39147 or [int(row["priority_rank"]) for row in records] != list(range(1, 39148)):
        raise ValueError("the complete 39,147-class native bisection batch is required")
    quadratics = [normalized_homogeneous_quadratic(record) for record in records]
    by_label = {record["label"]: (index, record) for index, record in enumerate(records)}

    local_tables = {}
    for ordinal, prime in enumerate(all_primes, start=1):
        split_table, usable = split_count_table(quadratics, prime)
        trace_table = surface_trace_table(direct, prime)
        local_tables[prime] = {"split": split_table, "trace": trace_table, "usable_bisections": usable}
        print(f"ALTV4PERIODTABLE|prime={prime}|progress={ordinal}/{len(all_primes)}|usable_bisections={usable}", flush=True)

    helper = load_module("alternate_v4_finite_helper", FINITE_HELPER)
    a_coefficients = direct["weierstrass_model"]["A_coefficients_low_to_high"]
    b_coefficients = direct["weierstrass_model"]["B_coefficients_low_to_high"]
    base_rows = []
    fibres = []
    total_exact_tests = 0
    selected_bases = bases_document["bases"][: args.base_limit]
    for base_ordinal, base in enumerate(selected_bases, start=1):
        periods, rejected = build_base_periods(base, blocks, local_tables)
        candidates = rank_integers(periods, blocks, args.n_bound, args.minimum_abs_n)
        selected = candidates[: args.exact_per_base]
        base_rows.append(
            {
                "shortlist_rank": base["shortlist_rank"],
                "pair_key": base["pair_key"],
                "rejected_primes": rejected,
                "periods": {str(prime): data for prime, data in periods.items()},
                "ranked_integers": candidates,
                "selected_n": [row["n"] for row in selected],
            }
        )
        for selected_row in selected:
            n_value = int(selected_row["n"])
            base_point, u_value, s_value, t_value = map_point_exact(base, n_value)
            numerator = ZZ(u_value.numerator())
            denominator = ZZ(u_value.denominator())
            survivor_indices, survivor_counts = modular_survivors(quadratics, numerator, denominator, all_primes)
            split = []
            collisions = []
            for index in survivor_indices:
                total_exact_tests += 1
                exact = exact_square_value(quadratics[index], numerator, denominator)
                if exact is None:
                    continue
                q_value, square_root = exact
                entry = {
                    "priority_rank": int(records[index]["priority_rank"]),
                    "label": records[index]["label"],
                    "lattice_orbit_mask": int(records[index]["lattice_orbit_mask"]),
                    "q_value": rational_text(q_value),
                    "canonical_positive_square_root": rational_text(square_root),
                }
                (collisions if q_value == 0 else split).append((entry, records[index], square_root))
            split_labels = {entry["label"] for entry, _record, _root in split}
            if not set(base["labels"]) <= split_labels:
                raise ArithmeticError("a V4 base point failed one of its defining split tests")

            coefficient_a = evaluate_qq(a_coefficients, u_value)
            coefficient_b = evaluate_qq(b_coefficients, u_value)
            curve = EllipticCurve(QQ, [coefficient_a, coefficient_b])
            generic = generic_sections(direct, u_value, curve)
            defining_points = []
            for label in base["labels"]:
                _index, record = by_label[label]
                defining_root = next(root for entry, _record, root in split if entry["label"] == label)
                defining_points.append(lifted_point(record, u_value, defining_root, curve))
            known_points = generic + defining_points
            model = tuple(to_fraction(value) for value in (0, 0, 0, coefficient_a, coefficient_b))
            known_certificate = finite_certificate(helper, model, known_points, args.certificate_prime_bound)
            if known_certificate["certified_rank_lower_bound"] < 19:
                raise ArithmeticError("the known rank-19 subgroup failed finite-quotient certification")

            seen = {(point[0], point[1]) for point in known_points}
            extras = []
            split_records = []
            for entry, record, square_root in split:
                point = lifted_point(record, u_value, square_root, curve)
                split_records.append({**entry, "point": point_text(point), "defining_character": entry["label"] in base["labels"]})
                if entry["label"] not in base["labels"] and (point[0], point[1]) not in seen:
                    extras.append(point)
                    seen.add((point[0], point[1]))
            total_points = known_points + extras
            total_certificate = known_certificate if not extras else finite_certificate(
                helper, model, total_points, args.certificate_prime_bound
            )
            rank_lower_bound = int(total_certificate["certified_rank_lower_bound"])
            promotion_status = promotion(rank_lower_bound)
            residual_threshold = 13
            residual = {
                "known_generic_subgroup_rank": 19,
                "rank_32_required_residual_2_selmer_dimension": residual_threshold,
                "status": (
                    "REQUIRED_BEFORE_UNRESTRICTED_POINT_SEARCH"
                    if rank_lower_bound >= 32
                    else "NOT_TRIGGERED_BELOW_RECORD_THRESHOLD"
                ),
                "result": None,
                "rejection_rule": "reject rank-32 candidacy if residual 2-Selmer dimension is below 13",
            }
            fibres.append(
                {
                    "pair_key": base["pair_key"],
                    "shortlist_rank": base["shortlist_rank"],
                    "n": n_value,
                    "period_score_rank": selected_row["period_score_rank"],
                    "period_scores": {key: value for key, value in selected_row.items() if "score" in key or key == "usable_primes_per_block"},
                    "base_point_on_paired_curve": point_text(base_point),
                    "v4_point": {"u": rational_text(u_value), "s": rational_text(s_value), "t": rational_text(t_value)},
                    "u_projective_bits": projective_bits(u_value),
                    "modular_survivors_after_each_prime": [
                        {"prime": prime, "survivors": count} for prime, count in zip(all_primes, survivor_counts)
                    ],
                    "exact_square_test_count": len(survivor_indices),
                    "split_bisection_count": len(split),
                    "branch_collisions": [entry for entry, _record, _root in collisions],
                    "specialized_model_a1_a2_a3_a4_a6": ["0", "0", "0", rational_text(coefficient_a), rational_text(coefficient_b)],
                    "generic_section_points": [point_text(point) for point in generic],
                    "defining_character_points": [point_text(point) for point in defining_points],
                    "known_rank_19_certificate": known_certificate,
                    "split_bisections": split_records,
                    "additional_distinct_split_point_count": len(extras),
                    "additional_split_points": [point_text(point) for point in extras],
                    "total_rank_certificate": total_certificate,
                    "certified_rank_lower_bound": rank_lower_bound,
                    "certified_quotient_gain_beyond_19": max(0, rank_lower_bound - 19),
                    "promotion": promotion_status,
                    "residual_descent_gate": residual,
                }
            )
            print(
                f"ALTV4FIBRE|base={base_ordinal}/{len(selected_bases)}|pair={base['pair_key']}|n={n_value}|"
                f"splits={len(split)}|rank_lb={rank_lower_bound}|promotion={promotion_status}",
                flush=True,
            )

    promotion_counts = {}
    rank_counts = {}
    for fibre in fibres:
        promotion_counts[fibre["promotion"]] = promotion_counts.get(fibre["promotion"], 0) + 1
        key = str(fibre["certified_rank_lower_bound"])
        rank_counts[key] = rank_counts.get(key, 0) + 1
    result = {
        "schema": SCHEMA,
        "status": "PASS_BOUNDED_NATIVE_ALTERNATE_V4_RECORD_CAMPAIGN",
        "inputs": {
            display_path(path): digest(path)
            for path in (Path(__file__).resolve(), args.bases, args.direct, args.bisections, FINITE_HELPER)
        },
        "bounds": {
            "integer_multiple_range": [-args.n_bound, args.n_bound],
            "minimum_absolute_multiple": args.minimum_abs_n,
            "exact_rational_fibres_per_base": args.exact_per_base,
            "rank_one_base_prefix": args.base_limit,
            "prime_blocks": [list(block) for block in blocks],
            "certificate_prime_bound": args.certificate_prime_bound,
            "native_bisection_count": len(records),
        },
        "promotion_thresholds": {
            "native_mechanism_control": 20,
            "useful_score_calibration": 24,
            "serious_alternate_q80_positive_control": 28,
            "record_candidate": 32,
        },
        "period_method": {
            "sequence": "n -> u(nP) modulo p",
            "complete_periods_precomputed": True,
            "local_signals": ["alternate-surface Nagao trace contribution", "native-bisection quadratic-residue count"],
            "normalization": "centered and population-standardized separately on every good base period",
            "ranking": "maximize the weakest combined disjoint-block score, then weakest split and trace blocks",
            "heuristic_only": True,
        },
        "local_tables": {
            str(prime): {"usable_bisections": data["usable_bisections"]} for prime, data in local_tables.items()
        },
        "bases": base_rows,
        "fibres": fibres,
        "summary": {
            "rank_one_bases": len(base_rows),
            "exact_fibres": len(fibres),
            "naive_square_tests": len(fibres) * len(records),
            "exact_square_tests_after_modular_rejection": total_exact_tests,
            "rank_lower_bound_counts": dict(sorted(rank_counts.items(), key=lambda item: int(item[0]))),
            "promotion_counts": promotion_counts,
            "maximum_certified_rank_lower_bound": max(int(row["certified_rank_lower_bound"]) for row in fibres),
            "maximum_split_bisection_count": max(int(row["split_bisection_count"]) for row in fibres),
        },
        "residual_descent_policy": {
            "rank_19_generic_subgroup_rejection_threshold": 13,
            "rank_20_generic_subgroup_rejection_threshold": 12,
            "rule": (
                "A completed residual 2-Selmer dimension below the applicable threshold rejects rank 32. "
                "Timeouts and incomplete descents remain UNKNOWN and authorize no unrestricted point search."
            ),
        },
        "runtime_seconds": perf_counter() - started,
        "software_assumptions": {"sage": SAGE_VERSION},
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "The period scores only rank the declared bounded multiples. Every stored rational split, point, "
            "and finite-quotient rank lower bound is exact. Absence of an additional split is bounded-negative, "
            "not a rank upper bound. No rank-32 point search is authorized unless the stored residual gate passes."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        f"ALTV4CAMPAIGN|bases={len(base_rows)}|fibres={len(fibres)}|"
        f"max_rank_lb={result['summary']['maximum_certified_rank_lower_bound']}|"
        f"status={result['status']}|output={display_path(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
