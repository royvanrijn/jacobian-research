#!/usr/bin/env sage-python
"""Certify rigid MW18 covers anchored at the refreshed extreme R17 fibres.

The historical rank-18 scan used one fixed published-R17 bisection.  This
replay instead transports the complete 39,120 norm-ten trace-character
inventory to the native ``07ca9`` and ``08234`` equations, reconstructs the
quadratic bisections there, and finds those splitting at the predeclared
``+10/+11/+12`` refreshed fibres.

Most nonsquares are rejected by exact finite-field nonresidues.  Every
survivor is reconstructed over QQ.  For every nonzero split the script checks
the cover and lifted-section identities, the two specialized branches and
their trace, and a rational conic parameterization through the extreme fibre.
It then uses numerical heights only to propose relations with the independent
displayed public subgroup and verifies the retained relations by exact group
addition.  A cover is called ``extreme anchored`` only when the resulting
class modulo specialized generic MW17 is proved nonzero.

The local checkpoint is only a performance cache.  ``--check`` rebuilds the
certificate from the canonical inputs and compares it byte-for-byte.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import hashlib
from math import isqrt
import json
from pathlib import Path
import re
import runpy
import sys
import time

from sage.all import (
    EllipticCurve,
    GF,
    PolynomialRing,
    QQ,
    ZZ,
    matrix,
    prime_range,
    vector,
)
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
PRIORITY = ROOT / "artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv"
PINNED_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
CHORD_HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
QUOTIENT_HELPER = ROOT / "elliptic-curves/scripts/evaluate_elkies_2026_bisections_at_controls.py"
REFRESH = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-refresh-priority-quotients-v1.json"
SWEEP = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
OVERVIEW = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve_refresh_475_573_overview_v1.json"
HISTORICAL = ROOT / "artifacts/generated-results/elliptic-curves/elkies_2026_bisection_specialization_controls_v1.json"
PUBLISHED_BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json"
CHECKPOINT_DIRECTORY = ROOT / "artifacts/local/elkies-k3/r17-extreme-anchored-mw18"

EXPECTED_RECORD_COUNT = 39120
# Large primes keep accidental collisions with the 24 singular fibres rare.
# This is the same window used by the complete 074d9 carrier-transfer replay.
SIEVE_PRIMES = tuple(map(int, prime_range(1009, 1300)))
PROTOCOL = "R17EXTREMEANCHORMW18"
STATUS = "PASS_EXACT_EXTREME_ANCHORED_MW18_COVERS"

CHARTS = {
    "07ca9": {
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit07ca9-direct-fibration-v1.json",
        "target_ids": (543, 544, 545),
    },
    "08234": {
        "direct": ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit08234-direct-fibration-v1.json",
        "target_ids": (531, 534, 535, 536, 537),
    },
}


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_text(polynomial) -> list[str]:
    if not polynomial:
        return ["0"]
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def rational_function_text(value) -> dict[str, list[str]]:
    return {
        "numerator_coefficients_low_to_high": polynomial_text(value.numerator()),
        "denominator_coefficients_low_to_high": polynomial_text(value.denominator()),
    }


def fraction(value) -> Fraction:
    value = QQ(value)
    return Fraction(int(value.numerator()), int(value.denominator()))


def point_text(point) -> list[str]:
    return [rational_text(point[0]), rational_text(point[1])]


def exact_square_root(value):
    value = QQ(value)
    if value < 0:
        return None
    numerator = isqrt(int(value.numerator()))
    denominator = isqrt(int(value.denominator()))
    if numerator * numerator != value.numerator() or denominator * denominator != value.denominator():
        return None
    return QQ(numerator) / QQ(denominator)


def load_integer_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(entry) for entry in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def parse_priority_rows(path: Path):
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    if len(rows) != EXPECTED_RECORD_COUNT:
        raise ArithmeticError("the published-R17 priority table changed size")
    masks = [int(row["orbit_mask"], 0) for row in rows]
    if len(set(masks)) != len(masks):
        raise ArithmeticError("the priority table has duplicate orbit masks")
    return rows


def parse_rational_function(record, ring, field):
    numerator = ring([QQ(value) for value in record["numerator_coefficients_low_to_high"]])
    denominator = ring([QQ(value) for value in record["denominator_coefficients_low_to_high"]])
    return field(numerator / denominator)


def transported_words(direct, priority_rows):
    pinned = load_integer_matrix(PINNED_GRAM)
    frame = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    isometry = matrix(
        ZZ, direct["frame_certificate"]["integral_isometry_to_published_R17"]
    )
    coordinate_matrix = matrix(
        ZZ, direct["sections"]["coordinate_matrix_in_compiled_frame"]
    )
    basis_gram = matrix(ZZ, direct["sections"]["height_gram"])
    if isometry * pinned * isometry.transpose() != frame:
        raise ArithmeticError("the stored published-R17 isometry has the wrong orientation")
    if coordinate_matrix * frame * coordinate_matrix.transpose() != basis_gram:
        raise ArithmeticError("the compiled section coordinates have the wrong Gram matrix")
    if abs(isometry.det()) != 1 or abs(coordinate_matrix.det()) != 1:
        raise ArithmeticError("the extreme chart does not have an integral saturated basis")
    transform = isometry.inverse() * coordinate_matrix.inverse()
    if any(entry not in ZZ for entry in transform.list()):
        raise ArithmeticError("the pinned-to-equation basis transform is not integral")
    words = []
    for row in priority_rows:
        pinned_word = vector(ZZ, [ZZ(value) for value in row["pinned_rank17_w"].split()])
        word = pinned_word * transform
        if any(entry not in ZZ for entry in word):
            raise ArithmeticError("a transported trace word is not integral")
        word = vector(ZZ, word)
        if word * basis_gram * word != 10:
            raise ArithmeticError("a transported trace word lost norm ten")
        words.append(tuple(map(int, word)))
    if len(set(words)) != EXPECTED_RECORD_COUNT:
        raise ArithmeticError("the integral transport identified two trace characters")
    return pinned, frame, basis_gram, isometry, coordinate_matrix, transform, words


def short_invariants(ainvs):
    a1, a2, a3, a4, a6 = map(QQ, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    return a1, a3, b2, -c4 / 48, -c6 / 864


def build_exact_context(direct, words):
    ring = PolynomialRing(QQ, "u")
    field = ring.fraction_field()
    A = ring([QQ(value) for value in direct["weierstrass_model"]["A_coefficients_low_to_high"]])
    B = ring([QQ(value) for value in direct["weierstrass_model"]["B_coefficients_low_to_high"]])
    discriminant = ring(-16 * (4 * A**3 + 27 * B**2))
    stored_delta = ring(
        [QQ(value) for value in direct["weierstrass_model"]["discriminant_coefficients_low_to_high"]]
    )
    if stored_delta != discriminant:
        raise ArithmeticError("the direct model discriminant changed")
    curve = EllipticCurve(field, [A, B])
    basis = []
    for expected_index, record in enumerate(direct["sections"]["records"]):
        if record["basis_index"] != expected_index or record["equation_verified"] is not True:
            raise ArithmeticError("the direct section order or exact status changed")
        basis.append(
            curve(
                parse_rational_function(record["X"], ring, field),
                parse_rational_function(record["Y"], ring, field),
            )
        )
    coefficient_sets = [set(word[index] for word in words) for index in range(17)]
    multiples = [
        {coefficient: coefficient * point for coefficient in coefficients}
        for point, coefficients in zip(basis, coefficient_sets)
    ]
    return ring, field, A, B, discriminant, curve, basis, multiples


def reduce_rational(value, finite_field):
    value = QQ(value)
    return finite_field(int(value.numerator())) / finite_field(int(value.denominator()))


def build_modular_context(direct, words, prime):
    finite_field = GF(prime)
    ring = PolynomialRing(finite_field, "u")
    field = ring.fraction_field()
    convert = lambda value: reduce_rational(value, finite_field)
    A = ring([convert(value) for value in direct["weierstrass_model"]["A_coefficients_low_to_high"]])
    B = ring([convert(value) for value in direct["weierstrass_model"]["B_coefficients_low_to_high"]])
    discriminant = ring(-16 * (4 * A**3 + 27 * B**2))
    curve = EllipticCurve(field, [A, B])
    basis = []
    for record in direct["sections"]["records"]:
        numerator_x = ring([convert(value) for value in record["X"]["numerator_coefficients_low_to_high"]])
        denominator_x = ring([convert(value) for value in record["X"]["denominator_coefficients_low_to_high"]])
        numerator_y = ring([convert(value) for value in record["Y"]["numerator_coefficients_low_to_high"]])
        denominator_y = ring([convert(value) for value in record["Y"]["denominator_coefficients_low_to_high"]])
        basis.append(curve(field(numerator_x / denominator_x), field(numerator_y / denominator_y)))
    coefficient_sets = [set(word[index] for word in words) for index in range(17)]
    multiples = [
        {coefficient: coefficient * point for coefficient in coefficients}
        for point, coefficients in zip(basis, coefficient_sets)
    ]
    return finite_field, ring, field, A, B, discriminant, curve, multiples


def quotient_from_frame(frame, A, ring):
    h, Nx, Ny, M = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    numerator = M**4 - 6 * M**2 * Nx - 8 * M * Ny - 3 * Nx**2 - 4 * A * h**4
    quotient, remainder = numerator.quo_rem(h**6)
    if remainder:
        raise ArithmeticError("the residual-chord numerator is not divisible by h^6")
    return ring(quotient)


def raw_quadratic(trace, A, ring, field, helper):
    X, Y = field(trace[0]), field(trace[1])
    frame = helper["trace_chord_frame"](X, Y, ring)
    if frame["h"].degree() == 3:
        return quotient_from_frame(frame, A, ring)
    reciprocal = helper["reciprocal_with_bound"]
    invert = helper["invert_rational"]
    inverse_frame = helper["trace_chord_frame"](
        invert(X, 4, ring, field), invert(Y, 6, ring, field), ring
    )
    if inverse_frame["h"].degree() != 3:
        raise ArithmeticError("neither affine chart exposes the three trace poles")
    inverse_q = quotient_from_frame(inverse_frame, reciprocal(A, 8, ring), ring)
    return reciprocal(inverse_q, 2, ring)


def exact_chord_data(trace, A, B, discriminant, ring, field, helper):
    X, Y = field(trace[0]), field(trace[1])
    frame = helper["trace_chord_frame"](X, Y, ring)
    if frame["h"].degree() == 3:
        data = helper["local_chord_data"](X, Y, A, B, discriminant, ring, field)
        data["construction_chart"] = "finite"
        return data
    reciprocal = helper["reciprocal_with_bound"]
    invert = helper["invert_rational"]
    inverse = helper["local_chord_data"](
        invert(X, 4, ring, field),
        invert(Y, 6, ring, field),
        reciprocal(A, 8, ring),
        reciprocal(B, 12, ring),
        reciprocal(discriminant, 24, ring),
        ring,
        field,
    )
    data = {
        key: reciprocal(inverse[key], bound, ring)
        for key, bound in (
            ("h", 3), ("M", 5), ("q", 2), ("sum_x", 4),
            ("product_x", 8), ("x0", 4), ("x1", 3), ("y0", 6), ("y1", 5),
        )
    }
    data["Nx"] = ring(X * data["h"] ** 2)
    data["Ny"] = ring(Y * data["h"] ** 3)
    data["construction_chart"] = "inverted_at_infinity"
    if data["sum_x"] ** 2 - 4 * data["product_x"] != data["h"] ** 2 * data["q"]:
        raise ArithmeticError("the reciprocal quadratic identity failed")
    if data["y0"] ** 2 + data["y1"] ** 2 * data["q"] != (
        data["x0"] ** 3 + 3 * data["x0"] * data["x1"] ** 2 * data["q"]
        + A * data["x0"] + B
    ):
        raise ArithmeticError("the reciprocal lifted-section identity failed")
    if 2 * data["y0"] * data["y1"] != (
        3 * data["x0"] ** 2 * data["x1"] + data["x1"] ** 3 * data["q"]
        + A * data["x1"]
    ):
        raise ArithmeticError("the reciprocal lifted-section linear identity failed")
    return data


def trace_from_word(word, curve, multiples):
    trace = sum(
        (multiples[index][coefficient] for index, coefficient in enumerate(word)),
        curve(0),
    )
    if trace.is_zero():
        raise ArithmeticError("a norm-ten trace specialized identically to zero")
    return trace


def input_records(refresh, sweep, overview, chart_key):
    target_ids = CHARTS[chart_key]["target_ids"]
    refresh_by_id = {int(record["curve_id"]): record for record in refresh["fibres"]}
    hits_by_id = {int(record["curve_id"]): record for record in sweep["rational_j_hits_and_twists"]}
    public_by_id = {int(record["id"]): record for record in overview["snapshot"]["records"]}
    result = {}
    for curve_id in target_ids:
        record = refresh_by_id[curve_id]
        expected_chart = f"norm12-orbit-{chart_key}"
        if record["native_chart"] != expected_chart:
            raise ArithmeticError(f"curve {curve_id} moved out of {expected_chart}")
        hit = hits_by_id[curve_id]
        native = next(item for item in hit["native_chart_twists"] if item["chart"] == expected_chart)
        if native["twist"]["status"] != "QQ_ISOMORPHIC_UNTWISTED":
            raise ArithmeticError(f"curve {curve_id} is no longer an untwisted rational fibre")
        parameter = QQ(native["native_parameter"]["numerator"]) / QQ(
            native["native_parameter"]["denominator"]
        )
        if rational_text(parameter) != record["native_parameter"]:
            raise ArithmeticError(f"curve {curve_id} parameter sources disagree")
        result[curve_id] = {
            "refresh": record,
            "hit": hit,
            "native": native,
            "public": public_by_id[curve_id],
            "parameter": parameter,
        }
    return result


def checkpoint_payload(input_key, chart_key, processed_primes, obstruction, validation):
    return {
        "schema": "elkies-k3.r17-extreme-anchored-mw18-sieve-checkpoint.v1",
        "input_key": input_key,
        "chart": chart_key,
        "processed_primes": list(processed_primes),
        "target_ids": list(CHARTS[chart_key]["target_ids"]),
        "obstruction_primes": {str(key): value for key, value in obstruction.items()},
        "global_validation_primes": validation,
    }


def write_checkpoint(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")) + "\n")


def modular_sieve(*, direct, words, targets, helper, input_key, chart_key, checkpoint, resume):
    target_ids = CHARTS[chart_key]["target_ids"]
    obstruction = {curve_id: [0] * len(words) for curve_id in target_ids}
    validation = [0] * len(words)
    processed_primes = []
    if resume and checkpoint.is_file():
        saved = json.loads(checkpoint.read_text())
        if saved.get("input_key") != input_key or saved.get("chart") != chart_key:
            raise ArithmeticError("the local sieve checkpoint has stale inputs")
        if saved.get("target_ids") != list(target_ids):
            raise ArithmeticError("the local sieve checkpoint has stale targets")
        processed_primes = [int(value) for value in saved["processed_primes"]]
        obstruction = {
            curve_id: [int(value) for value in saved["obstruction_primes"][str(curve_id)]]
            for curve_id in target_ids
        }
        validation = [int(value) for value in saved["global_validation_primes"]]

    for prime in SIEVE_PRIMES:
        if prime in processed_primes:
            continue
        started = time.monotonic()
        try:
            finite_field, ring, field, A, _B, discriminant, curve, multiples = (
                build_modular_context(direct, words, prime)
            )
            reduced_parameters = {
                curve_id: reduce_rational(targets[curve_id]["parameter"], finite_field)
                for curve_id in target_ids
            }
        except (ArithmeticError, ZeroDivisionError, ValueError):
            continue
        active = [
            index for index in range(len(words))
            if not validation[index]
            or any(not obstruction[curve_id][index] for curve_id in target_ids)
        ]
        valid_count = 0
        new_obstructions = 0
        for index in active:
            try:
                trace = trace_from_word(words[index], curve, multiples)
                q = raw_quadratic(trace, A, ring, field, helper)
            except (ArithmeticError, AssertionError, ZeroDivisionError, ValueError):
                continue
            if q.degree() != 2 or q.gcd(q.derivative()).degree() or q.gcd(discriminant).degree():
                continue
            valid_count += 1
            if not validation[index]:
                validation[index] = prime
            for curve_id in target_ids:
                if obstruction[curve_id][index]:
                    continue
                value = q(reduced_parameters[curve_id])
                if value and not value.is_square():
                    obstruction[curve_id][index] = prime
                    new_obstructions += 1
        processed_primes.append(prime)
        write_checkpoint(
            checkpoint,
            checkpoint_payload(input_key, chart_key, processed_primes, obstruction, validation),
        )
        unresolved_pairs = sum(
            not value for values in obstruction.values() for value in values
        )
        unvalidated = sum(not value for value in validation)
        print(
            f"{PROTOCOL}|chart={chart_key}|stage=sieve|prime={prime}|active={len(active)}"
            f"|valid={valid_count}|new_obstructions={new_obstructions}"
            f"|unresolved_pairs={unresolved_pairs}|unvalidated={unvalidated}"
            f"|seconds={time.monotonic()-started:.3f}",
            flush=True,
        )
        if not unresolved_pairs and not unvalidated:
            break
    return processed_primes, obstruction, validation


def conic_parameterization(q, t0, u0):
    parameter_ring = PolynomialRing(QQ, "r")
    parameter_field = parameter_ring.fraction_field()
    r = parameter_field(parameter_ring.gen())
    a0, a1, a2 = map(QQ, (q[0], q[1], q[2]))
    derivative_at_anchor = a1 + 2 * a2 * t0
    delta = parameter_field((derivative_at_anchor - 2 * u0 * r) / (r**2 - a2))
    t_of_r = parameter_field(t0 + delta)
    u_of_r = parameter_field(u0 + r * delta)
    q_of_t = a0 + a1 * t_of_r + a2 * t_of_r**2
    if u_of_r**2 != q_of_t:
        raise ArithmeticError("the anchor-line conic parameterization failed")
    return parameter_ring, parameter_field, t_of_r, u_of_r


def exact_survivors(*, direct, words, priority_rows, targets, obstruction, validation, helper, chart_key):
    target_ids = CHARTS[chart_key]["target_ids"]
    ring, field, A, B, discriminant, curve, _basis, multiples = build_exact_context(direct, words)
    survivor_indices = sorted(
        {
            index for index in range(len(words))
            if not validation[index]
            or any(not obstruction[curve_id][index] for curve_id in target_ids)
        }
    )
    splits = {curve_id: [] for curve_id in target_ids}
    ramified = {curve_id: [] for curve_id in target_ids}
    exact_nonsquare = {curve_id: [] for curve_id in target_ids}
    exact_global_validations = []
    for position, index in enumerate(survivor_indices, start=1):
        trace = trace_from_word(words[index], curve, multiples)
        data = exact_chord_data(trace, A, B, discriminant, ring, field, helper)
        q = ring(data["q"])
        if q.degree() != 2 or q.gcd(q.derivative()).degree() or q.gcd(discriminant).degree():
            raise ArithmeticError("an exact survivor is not a smooth rigid bisection")
        if not validation[index]:
            exact_global_validations.append(index)
        for curve_id in target_ids:
            if obstruction[curve_id][index]:
                continue
            t0 = targets[curve_id]["parameter"]
            value = QQ(q(t0))
            square_root = exact_square_root(value)
            if square_root is None:
                exact_nonsquare[curve_id].append(index)
                continue
            if not square_root:
                ramified[curve_id].append(index)
                continue
            x0 = QQ(data["x0"](t0))
            x1 = QQ(data["x1"](t0))
            y0 = QQ(data["y0"](t0))
            y1 = QQ(data["y1"](t0))
            positive = (x0 + x1 * square_root, y0 + y1 * square_root)
            negative = (x0 - x1 * square_root, y0 - y1 * square_root)
            fibre_curve = EllipticCurve(QQ, [QQ(A(t0)), QQ(B(t0))])
            positive_point = fibre_curve(positive)
            negative_point = fibre_curve(negative)
            trace_at_anchor = fibre_curve(QQ(trace[0](t0)), QQ(trace[1](t0)))
            if positive_point + negative_point != trace_at_anchor:
                raise ArithmeticError("the split branches do not add to the trace")
            parameter_ring, parameter_field, t_of_r, u_of_r = conic_parameterization(
                q, t0, square_root
            )
            evaluate = lambda poly: sum(
                parameter_field(poly[i]) * t_of_r**i for i in range(poly.degree() + 1)
            )
            section_x = evaluate(data["x0"]) + evaluate(data["x1"]) * u_of_r
            section_y = evaluate(data["y0"]) + evaluate(data["y1"]) * u_of_r
            if section_y**2 != section_x**3 + evaluate(A) * section_x + evaluate(B):
                raise ArithmeticError("the parameterized eighteenth section is off the surface")
            if not data["x1"] and not data["y1"]:
                raise ArithmeticError("the purported eighteenth section is Galois invariant")
            row = priority_rows[index]
            coefficient_bits = [
                max(abs(int(value.numerator())).bit_length(), int(value.denominator()).bit_length())
                for value in q
            ]
            splits[curve_id].append(
                {
                    "priority_index_zero_based": index,
                    "published_priority_rank": int(row["priority_rank"]),
                    "label": f"{chart_key}-orbit-{int(row['orbit_mask'], 0):05x}",
                    "lattice_orbit_mask": int(row["orbit_mask"], 0),
                    "pinned_rank17_trace_word": [int(value) for value in row["pinned_rank17_w"].split()],
                    "equation_basis_trace_word": list(words[index]),
                    "trace_height": 10,
                    "construction_chart": data["construction_chart"],
                    "branch_quadratic_coefficients_low_to_high": polynomial_text(q),
                    "branch_quadratic_coefficient_bits": coefficient_bits,
                    "branch_quadratic_maximum_coefficient_bits": max(coefficient_bits),
                    "q_at_anchor": rational_text(value),
                    "canonical_positive_square_root": rational_text(square_root),
                    "positive_chart_point": point_text(positive),
                    "negative_chart_point": point_text(negative),
                    "anchor_line_parameterization": {
                        "parameter": "r=(u-u0)/(t-t0)",
                        "t_of_r": rational_function_text(t_of_r),
                        "u_of_r": rational_function_text(u_of_r),
                        "passes_through_anchor_as_r_infinity": True,
                    },
                    "eighteenth_section": {
                        "x0_coefficients_low_to_high": polynomial_text(data["x0"]),
                        "x1_coefficients_low_to_high": polynomial_text(data["x1"]),
                        "y0_coefficients_low_to_high": polynomial_text(data["y0"]),
                        "y1_coefficients_low_to_high": polynomial_text(data["y1"]),
                    },
                    "checks": {
                        "exact_cover_and_lifted_section_identities": True,
                        "generic_galois_anti_invariant_nonzero": True,
                        "specialized_branches_on_fibre": True,
                        "specialized_branch_sum_is_trace": True,
                        "rational_parameterization_through_anchor": True,
                        "parameterized_section_on_base_changed_surface": True,
                    },
                }
            )
        print(
            f"{PROTOCOL}|chart={chart_key}|stage=exact|completed={position}/{len(survivor_indices)}"
            f"|index={index}|splits={sum(len(values) for values in splits.values())}",
            flush=True,
        )
    return splits, ramified, exact_nonsquare, exact_global_validations, (ring, A, B)


def exact_linear_combination(curve, coefficients, points):
    return sum(
        (int(coefficient) * point for coefficient, point in zip(coefficients, points)),
        curve(0),
    )


def public_and_generic_points(target, direct, ring, A, B):
    public_record = target["public"]
    a1, a3, b2, target_a, target_b = short_invariants(public_record["ainvs"])
    target_curve = EllipticCurve(QQ, [target_a, target_b])
    public_points = [
        target_curve(
            QQ(x) + b2 / 12,
            QQ(y) + (a1 * QQ(x) + a3) / 2,
        )
        for x, y in public_record["points"]
    ]
    t0 = target["parameter"]
    fibre_a, fibre_b = QQ(A(t0)), QQ(B(t0))
    scale_q = target_b * fibre_a / (fibre_b * target_a)
    scale_s = exact_square_root(scale_q)
    if scale_s is None or not scale_s:
        raise ArithmeticError("the exact native fibre acquired a nontrivial twist")
    if target_a != scale_q**2 * fibre_a or target_b != scale_q**3 * fibre_b:
        raise ArithmeticError("the exact native-to-public short isomorphism failed")
    generic_points = []
    for record in direct["sections"]["records"]:
        X = parse_rational_function(record["X"], ring, ring.fraction_field())
        Y = parse_rational_function(record["Y"], ring, ring.fraction_field())
        generic_points.append(target_curve(scale_q * X(t0), scale_s**3 * Y(t0)))
    stored = matrix(
        ZZ,
        target["refresh"]["specialized_generic_subgroup"][
            "coordinate_matrix_rows_in_ordered_public_points"
        ],
    )
    if stored.nrows() != len(public_points) or stored.ncols() != 17:
        raise ArithmeticError("the stored generic/public coordinate matrix changed shape")
    for column, point in enumerate(generic_points):
        if exact_linear_combination(target_curve, stored.column(column), public_points) != point:
            raise ArithmeticError("the direct generic basis disagrees with the refreshed quotient audit")
    return target_curve, public_points, generic_points, scale_q, scale_s


def discover_exact_relations(curve, public_points, split_points, quotient_helper):
    all_points = public_points + split_points
    model = (QQ(0), QQ(0), QQ(0), curve.a4(), curve.a6())
    model_text = ",".join(quotient_helper["gp_rational"](fraction(value)) for value in model)
    point_text = ",".join(
        "[{},{}]".format(
            quotient_helper["gp_rational"](fraction(point[0])),
            quotient_helper["gp_rational"](fraction(point[1])),
        )
        for point in all_points
    )
    program = f"""default(realprecision,140);
E=ellinit([{model_text}]);
P=[{point_text}];
M=ellheightmatrix(E,P);
T=qflllgram(M);
print("HEIGHTRANK|",matrank(M));
for(j=1,#P,if(abs(T[,j]~*M*T[,j])<1e-70,print("REL|",Vec(T[,j]))));
quit
"""
    output, _elapsed = quotient_helper["run_gp"](
        program, timeout=120.0, stack_bytes=1_024_000_000
    )
    rank_match = re.search(r"^HEIGHTRANK\|(\d+)$", output, re.MULTILINE)
    if rank_match is None:
        raise ArithmeticError("PARI omitted the numerical height rank")
    candidates = []
    for payload in re.findall(r"^REL\|\[([^]]*)\]$", output, re.MULTILINE):
        row = [int(value.strip()) for value in payload.split(",")]
        if len(row) != len(all_points):
            raise ArithmeticError("PARI returned a relation of the wrong length")
        candidates.append(row)
    public_count = len(public_points)
    new_count = len(split_points)
    selected = []
    new_rank = 0
    for row in candidates:
        candidate = selected + [row]
        candidate_rank = matrix(QQ, [item[public_count:] for item in candidate]).rank()
        if candidate_rank > new_rank:
            selected.append(row)
            new_rank = candidate_rank
        if new_rank == new_count:
            break
    if new_rank != new_count:
        raise ArithmeticError("PARI did not expose a full relation block for the split points")
    for row in selected:
        if exact_linear_combination(curve, row, all_points) != curve(0):
            raise ArithmeticError("a proposed height relation failed exact group addition")
    old_block = matrix(ZZ, [row[:public_count] for row in selected])
    new_block = matrix(ZZ, [row[public_count:] for row in selected])
    if new_block.nrows() != new_count or new_block.det() == 0:
        raise ArithmeticError("the exact relation new-point block is singular")
    coordinates = -new_block.inverse() * old_block
    return {
        "discovery_engine": "PARI/GP ellheightmatrix plus qflllgram",
        "real_precision_digits": 140,
        "zero_height_threshold": "1e-70",
        "numerical_height_matrix_rank": int(rank_match.group(1)),
        "relations": selected,
        "all_relations_verified_by_exact_group_addition": True,
        "new_point_block_determinant": int(new_block.det()),
        "coordinates_in_independent_public_basis": [
            [rational_text(value) for value in row] for row in coordinates.rows()
        ],
        "numerical_heights_used_in_proof": False,
    }, coordinates


def resolve_anchor_quotients(*, chart_key, direct, targets, splits, exact_context, quotient_helper):
    ring, A, B = exact_context
    results = {}
    for curve_id in CHARTS[chart_key]["target_ids"]:
        records = splits[curve_id]
        quotient_rank = int(
            targets[curve_id]["refresh"]["displayed_exceptional_quotient"]["free_rank"]
        )
        target_curve, public_points, _generic_points, scale_q, scale_s = public_and_generic_points(
            targets[curve_id], direct, ring, A, B
        )
        split_points = []
        for record in records:
            x, y = map(QQ, record["positive_chart_point"])
            point = target_curve(scale_q * x, scale_s**3 * y)
            split_points.append(point)
            record["positive_public_short_point"] = point_text(point)
        if not split_points:
            results[curve_id] = {
                "displayed_exceptional_quotient_rank": quotient_rank,
                "split_count": 0,
                "anchored_nonzero_count": 0,
                "anchored_exceptional_span_rank": 0,
                "remaining_displayed_directions_after_one_generic_section": None,
                "relation_certificate": None,
            }
            continue
        relation_record, coordinates = discover_exact_relations(
            target_curve, public_points, split_points, quotient_helper
        )
        quotient_coordinates = coordinates[:, 17:]
        anchored = []
        for index, record in enumerate(records):
            coordinates_row = list(coordinates.row(index))
            quotient_row = list(quotient_coordinates.row(index))
            nonzero = any(quotient_row)
            record["exact_class_in_displayed_public_group_tensor_Q"] = {
                "ordered_public_point_coordinates": [rational_text(value) for value in coordinates_row],
                "exceptional_basis": targets[curve_id]["refresh"]["displayed_exceptional_quotient"][
                    "preferred_public_quotient_basis"
                ],
                "exceptional_coordinates": [rational_text(value) for value in quotient_row],
                "nonzero_modulo_specialized_generic_MW17": nonzero,
                "derived_from_exact_integer_relations_and_public_independence": True,
            }
            record["extreme_anchored"] = nonzero
            if nonzero:
                anchored.append(index)
        span_rank = int(quotient_coordinates.rank())
        if span_rank == 0:
            remaining = None
        else:
            remaining = quotient_rank - 1
        results[curve_id] = {
            "displayed_exceptional_quotient_rank": quotient_rank,
            "split_count": len(records),
            "anchored_nonzero_count": len(anchored),
            "anchored_record_indices_zero_based": anchored,
            "anchored_exceptional_span_rank": span_rank,
            "remaining_displayed_directions_after_one_generic_section": remaining,
            "rank32_additional_directions_needed_beyond_anchor": (
                14 - remaining if remaining is not None else None
            ),
            "relation_certificate": relation_record,
        }
        print(
            f"{PROTOCOL}|chart={chart_key}|stage=quotient|curve={curve_id}"
            f"|splits={len(records)}|anchored={len(anchored)}|span={span_rank}",
            flush=True,
        )
    return results


def historical_anchor(priority_rows):
    historical = json.loads(HISTORICAL.read_text())
    fibre = next(record for record in historical["fibres"] if record["parameter"] == "-9529/5471")
    if fibre["split_bisection_count"] != 1:
        raise ArithmeticError("the historical rank-28 anchor split count changed")
    hit = fibre["hits"][0]
    label = hit["label"]
    full = json.loads(PUBLISHED_BISECTIONS.read_text())
    cover = next(record for record in full["bisections"] if record["label"] == label)
    q = [QQ(value) for value in cover["residual_chord"]["q_coefficients"]]
    t0 = QQ(-9529) / QQ(5471)
    q0 = sum(coefficient * t0**index for index, coefficient in enumerate(q))
    u0 = exact_square_root(q0)
    if u0 is None or not u0:
        raise ArithmeticError("the historical rank-28 cover no longer splits nontrivially")
    return {
        "source": relative(HISTORICAL),
        "curve_label": "published rank-at-least-28 control",
        "parameter": "-9529/5471",
        "displayed_exceptional_quotient_rank": 11,
        "split_count": 1,
        "anchored_nonzero_count": 1,
        "cover_label": label,
        "branch_quadratic_coefficients_low_to_high": [rational_text(value) for value in q],
        "q_at_anchor": rational_text(q0),
        "canonical_positive_square_root": rational_text(u0),
        "exact_nonzero_class_modulo_generic_MW17": hit[
            "finite_quotient_class_modulo_generic_17"
        ]["coordinates_over_f2"],
        "remaining_displayed_directions_after_one_generic_section": 10,
        "rank32_additional_directions_needed_beyond_anchor": 4,
        "generic_rank_consequence": "at least 18",
    }


def build_chart(chart_key, priority_rows, refresh, sweep, overview, helper, quotient_helper, resume):
    direct_path = CHARTS[chart_key]["direct"]
    direct = json.loads(direct_path.read_text())
    if direct["status"] != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ArithmeticError(f"the {chart_key} direct fibration is not exact")
    if direct["sections"]["status"] != "PASS_EXACT_SATURATED_RANK17_BASIS":
        raise ArithmeticError(f"the {chart_key} direct basis is not saturated")
    targets = input_records(refresh, sweep, overview, chart_key)
    pinned, frame, basis_gram, isometry, coordinate_matrix, transform, words = (
        transported_words(direct, priority_rows)
    )
    input_hashes = {
        relative(path): digest(path)
        for path in (direct_path, PRIORITY, PINNED_GRAM, CHORD_HELPER, QUOTIENT_HELPER, REFRESH, SWEEP, OVERVIEW)
    }
    input_key = hashlib.sha256(
        json.dumps(input_hashes, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    checkpoint = CHECKPOINT_DIRECTORY / f"{chart_key}-sieve-v1.json"
    processed_primes, obstruction, validation = modular_sieve(
        direct=direct,
        words=words,
        targets=targets,
        helper=helper,
        input_key=input_key,
        chart_key=chart_key,
        checkpoint=checkpoint,
        resume=resume,
    )
    splits, ramified, exact_nonsquare, exact_global_validations, exact_model = exact_survivors(
        direct=direct,
        words=words,
        priority_rows=priority_rows,
        targets=targets,
        obstruction=obstruction,
        validation=validation,
        helper=helper,
        chart_key=chart_key,
    )
    unresolved = {
        curve_id: [index for index, prime in enumerate(obstruction[curve_id]) if not prime]
        for curve_id in CHARTS[chart_key]["target_ids"]
    }
    for curve_id in CHARTS[chart_key]["target_ids"]:
        decided = len(exact_nonsquare[curve_id]) + len(ramified[curve_id]) + len(splits[curve_id])
        if decided != len(unresolved[curve_id]):
            raise ArithmeticError("the exact survivor decisions are incomplete")
    quotient_results = resolve_anchor_quotients(
        chart_key=chart_key,
        direct=direct,
        targets=targets,
        splits=splits,
        exact_context=(exact_model[0], exact_model[1], exact_model[2]),
        quotient_helper=quotient_helper,
    )
    fibres = []
    for curve_id in CHARTS[chart_key]["target_ids"]:
        records = splits[curve_id]
        records.sort(
            key=lambda record: (
                not record.get("extreme_anchored", False),
                record["branch_quadratic_maximum_coefficient_bits"],
                record["published_priority_rank"],
            )
        )
        quotient = quotient_results[curve_id]
        fibres.append(
            {
                "curve_id": curve_id,
                "native_parameter": rational_text(targets[curve_id]["parameter"]),
                "snapshot_rank_lower_bound": int(targets[curve_id]["refresh"]["snapshot_rank_lower_bound"]),
                "displayed_jump_over_MW17": quotient["displayed_exceptional_quotient_rank"],
                "ramified_cover_count": len(ramified[curve_id]),
                **quotient,
                "covers": records,
            }
        )
    return {
        "chart": f"norm12-orbit-{chart_key}",
        "direct_model": relative(direct_path),
        "target_ids": list(CHARTS[chart_key]["target_ids"]),
        "lattice_transport": {
            "record_count": len(words),
            "all_words_integral_distinct_and_norm_ten": True,
            "matrix_convention": "S*G_published*S^T=G_frame; C*G_frame*C^T=G_equation; w_equation=w_published*S^-1*C^-1",
            "isometry_to_published_R17": [list(map(int, row)) for row in isometry.rows()],
            "equation_basis_coordinates_in_frame": [list(map(int, row)) for row in coordinate_matrix.rows()],
            "pinned_to_equation_word_transform": [list(map(int, row)) for row in transform.rows()],
            "determinants": {
                "published": int(pinned.det()),
                "frame": int(frame.det()),
                "equation_basis": int(basis_gram.det()),
            },
        },
        "sieve": {
            "declared_primes": list(SIEVE_PRIMES),
            "processed_primes": processed_primes,
            "global_validation_primes_in_priority_order": validation,
            "exact_global_validation_indices_zero_based": exact_global_validations,
            "obstruction_primes_in_priority_order": {
                str(curve_id): obstruction[curve_id]
                for curve_id in CHARTS[chart_key]["target_ids"]
            },
            "obstruction_prime_histograms": {
                str(curve_id): {
                    str(prime): count
                    for prime, count in sorted(Counter(obstruction[curve_id]).items())
                    if prime
                }
                for curve_id in CHARTS[chart_key]["target_ids"]
            },
            "exact_nonsquare_survivor_indices_zero_based": {
                str(curve_id): exact_nonsquare[curve_id]
                for curve_id in CHARTS[chart_key]["target_ids"]
            },
            "exact_ramified_indices_zero_based": {
                str(curve_id): ramified[curve_id]
                for curve_id in CHARTS[chart_key]["target_ids"]
            },
        },
        "fibres": fibres,
        "inputs": input_hashes,
    }


def build(selected_charts, resume):
    priority_rows = parse_priority_rows(PRIORITY)
    refresh = json.loads(REFRESH.read_text())
    sweep = json.loads(SWEEP.read_text())
    overview = json.loads(OVERVIEW.read_text())
    if refresh["status"] != "PASS_EXACT_REFRESH_ATLAS_HIT_SPECIALIZATION_AUDIT":
        raise ArithmeticError("the refreshed quotient certificate is not exact")
    if sweep["status"] != "PASS_EXACT_COMPLETE_PINNED_ICARM_J_PREIMAGE_AND_TWIST_SWEEP":
        raise ArithmeticError("the refreshed atlas sweep is not exact")
    if overview["status"] != "PASS_EXACT_OVERVIEW_OF_ICARM_CURVES_475_THROUGH_573":
        raise ArithmeticError("the refreshed ICARM projection is not exact")
    helper = runpy.run_path(str(CHORD_HELPER))
    quotient_helper = runpy.run_path(str(QUOTIENT_HELPER))
    chart_records = [
        build_chart(
            chart_key, priority_rows, refresh, sweep, overview, helper, quotient_helper, resume
        )
        for chart_key in selected_charts
    ]
    all_fibres = [fibre for chart in chart_records for fibre in chart["fibres"]]
    anchored_count = sum(fibre["anchored_nonzero_count"] for fibre in all_fibres)
    return {
        "schema": "elkies-k3.r17-extreme-anchored-mw18-covers.v1",
        "status": STATUS,
        "claim": (
            "Complete exact rigid-bisection split census on the selected native "
            "published-R17-frame charts at the refreshed +10/+11/+12 fibres, with "
            "cover, section, anchor-parameterization, and displayed-quotient checks."
        ),
        "design_change": {
            "retired_default": "continue the historical first-conic Nagao scan",
            "new_default": (
                "search only after selecting an exact rational bisection whose new "
                "section specializes nontrivially in a known extreme quotient"
            ),
            "historical_first_cover_is_retained_as_regression": True,
        },
        "charts": chart_records,
        "historical_rank28_anchor": historical_anchor(priority_rows),
        "summary": {
            "selected_charts": list(selected_charts),
            "refreshed_fibre_count": len(all_fibres),
            "refreshed_curve_ids": [fibre["curve_id"] for fibre in all_fibres],
            "refreshed_split_cover_count": sum(fibre["split_count"] for fibre in all_fibres),
            "refreshed_extreme_anchored_cover_count": anchored_count,
            "generic_rank_lower_bound_on_every_listed_cover": 18,
            "generic_rank_is_not_proved_exact": True,
        },
        "rank_accounting": {
            "rank32_requirement_over_a_generic_rank_at_least_18_cover": "at most 14 further independent specialization directions",
            "plus12_anchor": "one nonzero anchored direction leaves at least eleven displayed exceptional directions, so three further directions would reach 14",
            "plus11_anchor": "one nonzero anchored direction leaves at least ten displayed exceptional directions, so four further directions would reach 14",
            "plus10_anchor": "one nonzero anchored direction leaves at least nine displayed exceptional directions, so five further directions would reach 14",
        },
        "claim_boundary": [
            "Every listed extreme-anchored cover has generic Mordell-Weil rank at least 18, not a proved exact generic rank of 18.",
            "Displayed jumps and remaining directions refer to certified independent displayed subgroups, not full Mordell-Weil rank upper bounds.",
            "The finite-field sieve proves nonsquareness only when it stores a good-prime nonresidue; every survivor is decided over QQ.",
            "Canonical heights are used only to propose relations; every retained relation is verified by exact elliptic-curve addition.",
            "No Nagao score, residual-Selmer computation, or new rank-32 specialization search is performed by this certificate.",
        ],
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact QQ and rational-function elliptic-curve arithmetic",
                "exact finite-field polynomial and square tests",
                "PARI/GP height relations followed by exact group-law verification",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/certify_r17_extreme_anchored_mw18_covers.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chart", choices=("all", *CHARTS), default="all")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    selected = tuple(CHARTS) if args.chart == "all" else (args.chart,)
    payload = build(selected, not args.no_resume and not args.check)
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists() or args.output.read_text() != serialized:
            raise SystemExit("anchored MW18 certificate differs from the pinned artifact")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"{PROTOCOL}|status={payload['status']}"
        f"|fibres={payload['summary']['refreshed_fibre_count']}"
        f"|anchored={payload['summary']['refreshed_extreme_anchored_cover_count']}"
        f"|output={relative(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
