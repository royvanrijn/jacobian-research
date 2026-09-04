#!/usr/bin/env sage-python
"""Certify complete rigid-bisection transfer across the five 074d9 fibres.

The historical 39,120-row bisection file is written in the original
published-R17 equation coordinates.  The five wgxli fibres do not occur at
rational parameters on that particular equation.  They occur on the
``norm12-orbit-074d9`` chart, whose saturated Mordell--Weil lattice is
integrally isometric to the published R17 lattice.  This replay transports
all 39,120 norm-ten trace characters through one exact integral isometry,
reconstructs their Proposition-F1 residual quadratics on 074d9, and tests
them at curves 351, 356, 376, 377, and 385.

Most nonsquares are certified cheaply by a good modular nonresidue.  A
modular calculation is accepted only when the reconstructed quadratic has
degree two, is squarefree, and is coprime to the reduced 24I1 discriminant.
Every pair surviving the declared prime list is reconstructed over QQ and
decided by exact integer square tests.  Exact split points are constructed in
both branches and transported to the public curve.  Their classes are then
resolved in the displayed free quotient whenever exact group relations prove
that they belong to the displayed subgroup.

The local checkpoint is a performance aid, not a proof input.  ``--check``
always rebuilds the result from the canonical inputs.
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
    pari,
    prime_range,
    vector,
)
from sage.env import SAGE_VERSION


sys.set_int_max_str_digits(0)
ROOT = Path(__file__).resolve().parents[2]
LINEAGE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)
PUBLIC = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
)
PRIORITY = (
    ROOT
    / "artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv"
)
PUBLISHED_ATLAS = (
    ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
)
PINNED_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
CHORD_HELPER = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
QUOTIENT_HELPER = (
    ROOT / "elliptic-curves/scripts/evaluate_elkies_2026_bisections_at_controls.py"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json"
)
CHECKPOINT = (
    ROOT
    / "artifacts/local/elkies-k3/r17-074d9-cross-fibre-bisection-sieve-v1.json"
)

CHART = "norm12-orbit-074d9"
TARGET_IDS = (351, 356, 376, 377, 385)
SIEVE_PRIMES = tuple(map(int, prime_range(1009, 1300)))
EXPECTED_RECORD_COUNT = 39120
STATUS = "PASS_EXACT_COMPLETE_074D9_CROSS_FIBRE_BISECTION_TRANSFER"
PROTOCOL = "R17074D9CROSSBISECT"


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


def representative_words(lineage, priority_rows):
    representative_gram = matrix(ZZ, lineage["generic_basis"]["height_gram"])
    pinned_gram = load_integer_matrix(PINNED_GRAM)
    isometry = matrix(ZZ, pari(representative_gram).qfisom(pari(pinned_gram)))
    if not isometry or abs(isometry.det()) != 1:
        raise ArithmeticError("no integral published-R17 lattice isometry was recovered")
    if isometry.transpose() * pinned_gram * isometry != representative_gram:
        raise ArithmeticError("the recovered lattice isometry has the wrong orientation")
    inverse = isometry.inverse()
    words = []
    for row in priority_rows:
        pinned_word = vector(ZZ, [ZZ(entry) for entry in row["pinned_rank17_w"].split()])
        word = inverse * pinned_word
        if any(entry not in ZZ for entry in word):
            raise ArithmeticError("a transported trace word is not integral")
        word = vector(ZZ, word)
        if word * representative_gram * word != 10:
            raise ArithmeticError("a transported trace lost norm ten")
        words.append(tuple(map(int, word)))
    if len(set(words)) != EXPECTED_RECORD_COUNT:
        raise ArithmeticError("the integral transport identified two trace characters")
    return representative_gram, pinned_gram, isometry, words


def short_invariants(ainvs):
    a1, a2, a3, a4, a6 = map(QQ, ainvs)
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    return {
        "a1": a1,
        "a3": a3,
        "b2": b2,
        "A": -c4 / 48,
        "B": -c6 / 864,
    }


def build_exact_context(lineage):
    ring = PolynomialRing(QQ, "u")
    field = ring.fraction_field()
    representative = lineage["representative"]
    A = ring([QQ(value) for value in representative["A_coefficients_low_to_high"]])
    B = ring([QQ(value) for value in representative["B_coefficients_low_to_high"]])
    discriminant = ring(-16 * (4 * A**3 + 27 * B**2))
    curve = EllipticCurve(field, [A, B])
    basis = [
        curve(
            ring([QQ(value) for value in record["representative_x_coefficients_low_to_high"]]),
            ring([QQ(value) for value in record["representative_y_coefficients_low_to_high"]]),
        )
        for record in representative["sections"]
    ]
    if len(basis) != 17:
        raise ArithmeticError("the 074d9 representative lost its seventeen sections")
    multiples = [
        {coefficient: coefficient * point for coefficient in range(-7, 8)}
        for point in basis
    ]
    return ring, field, A, B, discriminant, curve, basis, multiples


def reduce_rational(value, finite_field):
    value = QQ(value)
    return finite_field(int(value.numerator())) / finite_field(int(value.denominator()))


def build_modular_context(lineage, prime):
    finite_field = GF(prime)
    ring = PolynomialRing(finite_field, "u")
    field = ring.fraction_field()
    representative = lineage["representative"]
    convert = lambda value: reduce_rational(value, finite_field)
    A = ring([convert(value) for value in representative["A_coefficients_low_to_high"]])
    B = ring([convert(value) for value in representative["B_coefficients_low_to_high"]])
    discriminant = ring(-16 * (4 * A**3 + 27 * B**2))
    curve = EllipticCurve(field, [A, B])
    basis = [
        curve(
            ring([convert(value) for value in record["representative_x_coefficients_low_to_high"]]),
            ring([convert(value) for value in record["representative_y_coefficients_low_to_high"]]),
        )
        for record in representative["sections"]
    ]
    multiples = [
        {coefficient: coefficient * point for coefficient in range(-7, 8)}
        for point in basis
    ]
    return finite_field, ring, field, A, B, discriminant, curve, multiples


def quotient_from_frame(frame, A, ring):
    h, Nx, Ny, M = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    numerator = M**4 - 6 * M**2 * Nx - 8 * M * Ny - 3 * Nx**2 - 4 * A * h**4
    quotient, remainder = numerator.quo_rem(h**6)
    if remainder:
        raise ArithmeticError("the residual-chord numerator is not divisible by h^6")
    return ring(quotient)


def raw_quadratic(trace, A, B, ring, field, chord_helper):
    X, Y = field(trace[0]), field(trace[1])
    frame = chord_helper["trace_chord_frame"](X, Y, ring)
    if frame["h"].degree() == 3:
        return quotient_from_frame(frame, A, ring), "finite"
    reciprocal = chord_helper["reciprocal_with_bound"]
    invert = chord_helper["invert_rational"]
    inverse_frame = chord_helper["trace_chord_frame"](
        invert(X, 4, ring, field), invert(Y, 6, ring, field), ring
    )
    if inverse_frame["h"].degree() != 3:
        raise ArithmeticError("neither affine chart exposes the three trace poles")
    inverse_q = quotient_from_frame(
        inverse_frame, reciprocal(A, 8, ring), ring
    )
    return reciprocal(inverse_q, 2, ring), "inverted_at_infinity"


def exact_chord_data(trace, A, B, discriminant, ring, field, chord_helper):
    X, Y = field(trace[0]), field(trace[1])
    frame = chord_helper["trace_chord_frame"](X, Y, ring)
    if frame["h"].degree() == 3:
        data = chord_helper["local_chord_data"](
            X, Y, A, B, discriminant, ring, field
        )
        data["construction_chart"] = "finite"
        return data
    reciprocal = chord_helper["reciprocal_with_bound"]
    invert = chord_helper["invert_rational"]
    inverse = chord_helper["local_chord_data"](
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
            ("h", 3),
            ("M", 5),
            ("q", 2),
            ("sum_x", 4),
            ("product_x", 8),
            ("x0", 4),
            ("x1", 3),
            ("y0", 6),
            ("y1", 5),
        )
    }
    data["Nx"] = ring(X * data["h"] ** 2)
    data["Ny"] = ring(Y * data["h"] ** 3)
    data["construction_chart"] = "inverted_at_infinity"
    if data["sum_x"] ** 2 - 4 * data["product_x"] != data["h"] ** 2 * data["q"]:
        raise ArithmeticError("the reciprocal quadratic identity failed")
    if data["y0"] ** 2 + data["y1"] ** 2 * data["q"] != (
        data["x0"] ** 3
        + 3 * data["x0"] * data["x1"] ** 2 * data["q"]
        + A * data["x0"]
        + B
    ):
        raise ArithmeticError("the reciprocal lifted-section identity failed")
    if 2 * data["y0"] * data["y1"] != (
        3 * data["x0"] ** 2 * data["x1"]
        + data["x1"] ** 3 * data["q"]
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


def target_parameters(lineage):
    records = {
        int(record["curve_id"]): record
        for record in lineage["target_isomorphisms"]
        if record["chart"] == CHART
    }
    if tuple(sorted(records)) != TARGET_IDS:
        raise ArithmeticError("the five 074d9 target isomorphisms changed")
    return {curve_id: QQ(records[curve_id]["parameter"]) for curve_id in TARGET_IDS}, records


def checkpoint_payload(input_key, processed_primes, obstruction, validation):
    return {
        "schema": "elkies-k3.r17-074d9-cross-fibre-bisection-sieve-checkpoint.v1",
        "input_key": input_key,
        "processed_primes": list(processed_primes),
        "target_ids": list(TARGET_IDS),
        "obstruction_primes": {str(key): value for key, value in obstruction.items()},
        "global_validation_primes": validation,
    }


def write_checkpoint(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, separators=(",", ":")) + "\n")


def modular_sieve(
    *, lineage, words, parameters, chord_helper, input_key, checkpoint, resume
):
    obstruction = {curve_id: [0] * len(words) for curve_id in TARGET_IDS}
    validation = [0] * len(words)
    processed_primes = []
    if resume and checkpoint.is_file():
        saved = json.loads(checkpoint.read_text())
        if saved.get("input_key") != input_key:
            raise ArithmeticError("the local sieve checkpoint has stale inputs")
        if saved.get("target_ids") != list(TARGET_IDS):
            raise ArithmeticError("the local sieve checkpoint has stale targets")
        processed_primes = [int(value) for value in saved["processed_primes"]]
        obstruction = {
            curve_id: [int(value) for value in saved["obstruction_primes"][str(curve_id)]]
            for curve_id in TARGET_IDS
        }
        validation = [int(value) for value in saved["global_validation_primes"]]
        if any(len(values) != len(words) for values in obstruction.values()) or len(validation) != len(words):
            raise ArithmeticError("the local sieve checkpoint has stale dimensions")

    for prime in SIEVE_PRIMES:
        if prime in processed_primes:
            continue
        started = time.monotonic()
        try:
            finite_field, ring, field, A, B, discriminant, curve, multiples = (
                build_modular_context(lineage, prime)
            )
            reduced_parameters = {
                curve_id: reduce_rational(parameters[curve_id], finite_field)
                for curve_id in TARGET_IDS
            }
        except (ArithmeticError, ZeroDivisionError, ValueError):
            continue

        active = [
            index
            for index in range(len(words))
            if not validation[index]
            or any(not obstruction[curve_id][index] for curve_id in TARGET_IDS)
        ]
        valid_count = 0
        new_obstructions = 0
        for index in active:
            try:
                trace = trace_from_word(words[index], curve, multiples)
                q, _chart = raw_quadratic(trace, A, B, ring, field, chord_helper)
            except (ArithmeticError, AssertionError, ZeroDivisionError, ValueError):
                continue
            globally_good = (
                q.degree() == 2
                and q.gcd(q.derivative()).degree() == 0
                and q.gcd(discriminant).degree() == 0
            )
            if not globally_good:
                continue
            valid_count += 1
            if not validation[index]:
                validation[index] = prime
            for curve_id in TARGET_IDS:
                if obstruction[curve_id][index]:
                    continue
                value = q(reduced_parameters[curve_id])
                if value and not value.is_square():
                    obstruction[curve_id][index] = prime
                    new_obstructions += 1
        processed_primes.append(prime)
        write_checkpoint(
            checkpoint,
            checkpoint_payload(input_key, processed_primes, obstruction, validation),
        )
        unresolved_pairs = sum(
            not value for values in obstruction.values() for value in values
        )
        unvalidated = sum(not value for value in validation)
        print(
            f"{PROTOCOL}|stage=sieve|prime={prime}|active={len(active)}"
            f"|valid={valid_count}|new_obstructions={new_obstructions}"
            f"|unresolved_pairs={unresolved_pairs}|unvalidated={unvalidated}"
            f"|seconds={time.monotonic()-started:.3f}",
            flush=True,
        )
        if not unresolved_pairs and not unvalidated:
            break
    return processed_primes, obstruction, validation


def reduce_exact_polynomial(polynomial, prime):
    finite_field = GF(prime)
    ring = PolynomialRing(finite_field, "u")
    return ring([reduce_rational(value, finite_field) for value in polynomial])


def exact_survivors(
    *, lineage, words, priority_rows, parameters, obstruction, validation, chord_helper
):
    ring, field, A, B, discriminant, curve, _basis, multiples = build_exact_context(lineage)
    survivor_indices = sorted(
        {
            index
            for index in range(len(words))
            if not validation[index]
            or any(not obstruction[curve_id][index] for curve_id in TARGET_IDS)
        }
    )
    splits = {curve_id: [] for curve_id in TARGET_IDS}
    exact_nonsquare = {curve_id: [] for curve_id in TARGET_IDS}
    exact_global_validations = []
    for position, index in enumerate(survivor_indices, start=1):
        trace = trace_from_word(words[index], curve, multiples)
        data = exact_chord_data(trace, A, B, discriminant, ring, field, chord_helper)
        q = ring(data["q"])
        if q.degree() != 2 or q.gcd(q.derivative()).degree() or q.gcd(discriminant).degree():
            raise ArithmeticError("an exact survivor is not a smooth rigid bisection")
        if not validation[index]:
            exact_global_validations.append(index)
        else:
            prime = validation[index]
            modular_q = reduce_exact_polynomial(q, prime)
            if (
                modular_q.degree() != 2
                or modular_q.gcd(modular_q.derivative()).degree()
                or modular_q.gcd(reduce_exact_polynomial(discriminant, prime)).degree()
            ):
                raise ArithmeticError("an exact cover disagrees with its validation prime")

        for curve_id in TARGET_IDS:
            if obstruction[curve_id][index]:
                continue
            value = QQ(q(parameters[curve_id]))
            square_root = exact_square_root(value)
            if square_root is None:
                exact_nonsquare[curve_id].append(index)
                continue
            x0 = QQ(data["x0"](parameters[curve_id]))
            x1 = QQ(data["x1"](parameters[curve_id]))
            y0 = QQ(data["y0"](parameters[curve_id]))
            y1 = QQ(data["y1"](parameters[curve_id]))
            positive = (x0 + x1 * square_root, y0 + y1 * square_root)
            negative = (x0 - x1 * square_root, y0 - y1 * square_root)
            fibre_curve = EllipticCurve(
                QQ, [QQ(A(parameters[curve_id])), QQ(B(parameters[curve_id]))]
            )
            positive_point = fibre_curve(positive)
            negative_point = fibre_curve(negative)
            if positive_point + negative_point != fibre_curve(trace[0](parameters[curve_id]), trace[1](parameters[curve_id])):
                raise ArithmeticError("the split branches do not add to the frozen trace")
            row = priority_rows[index]
            splits[curve_id].append(
                {
                    "priority_index_zero_based": index,
                    "priority_rank": int(row["priority_rank"]),
                    "label": f"074d9-orbit-{int(row['orbit_mask'], 0):05x}",
                    "lattice_orbit_mask": int(row["orbit_mask"], 0),
                    "pinned_rank17_trace_word": [
                        int(value) for value in row["pinned_rank17_w"].split()
                    ],
                    "representative_rank17_trace_word": list(words[index]),
                    "construction_chart": data["construction_chart"],
                    "branch_quadratic_coefficients_low_to_high": polynomial_text(q),
                    "q_at_fibre": rational_text(value),
                    "canonical_square_root": rational_text(square_root),
                    "positive_chart_point": point_text(positive),
                    "negative_chart_point": point_text(negative),
                    "exact_cover_and_branch_identities_verified": True,
                }
            )
        print(
            f"{PROTOCOL}|stage=exact|completed={position}/{len(survivor_indices)}"
            f"|index={index}|splits={sum(len(values) for values in splits.values())}",
            flush=True,
        )
    return splits, exact_nonsquare, exact_global_validations


def public_short_points(public_record):
    invariants = short_invariants(public_record["ainvs"])
    curve = EllipticCurve(QQ, [invariants["A"], invariants["B"]])
    points = [
        curve(
            QQ(x) + invariants["b2"] / 12,
            QQ(y) + (invariants["a1"] * QQ(x) + invariants["a3"]) / 2,
        )
        for x, y in public_record["points"]
    ]
    return invariants, curve, points


def discover_relations(
    short_model, points, existing_count, quotient_helper
):
    """Discover a full height-relation block with PARI/GP."""

    model_text = ",".join(
        quotient_helper["gp_rational"](fraction(value)) for value in short_model
    )
    point_text = ",".join(
        "[{},{}]".format(
            quotient_helper["gp_rational"](fraction(point[0])),
            quotient_helper["gp_rational"](fraction(point[1])),
        )
        for point in points
    )
    program = f"""default(realprecision,100);
E=ellinit([{model_text}]);
P=[{point_text}];
M=ellheightmatrix(E,P);
T=qflllgram(M);
print("HEIGHTRANK|",matrank(M));
for(j=1,#P,if(abs(T[,j]~*M*T[,j])<1e-50,print("REL|",Vec(T[,j]))));
quit
"""
    output, _elapsed = quotient_helper["run_gp"](
        program, timeout=60.0, stack_bytes=512_000_000
    )
    rank_match = re.search(r"^HEIGHTRANK\|(\d+)$", output, re.MULTILINE)
    if rank_match is None:
        raise ArithmeticError("PARI omitted the numerical height rank")
    discovered = []
    for payload in re.findall(r"^REL\|\[([^]]*)\]$", output, re.MULTILINE):
        row = [int(value.strip()) for value in payload.split(",")]
        if len(row) != len(points):
            raise ArithmeticError("PARI returned a relation of the wrong length")
        discovered.append(row)
    new_count = len(points) - existing_count
    selected = []
    new_rank = 0
    for row in discovered:
        candidate = selected + [row]
        candidate_rank = matrix(QQ, [item[existing_count:] for item in candidate]).rank()
        if candidate_rank > new_rank:
            selected.append(row)
            new_rank = candidate_rank
        if new_rank == new_count:
            break
    if new_rank != new_count:
        raise ArithmeticError("PARI did not expose a full relation block")
    return {
        "discovery_engine": "PARI/GP ellheightmatrix plus qflllgram",
        "real_precision_digits": 100,
        "zero_height_threshold": "1e-50",
        "numerical_height_matrix_rank": int(rank_match.group(1)),
        "relation_count": len(selected),
        "new_coefficient_block_rank_over_q": new_rank,
        "candidate_relations_pending_exact_coordinate_verification": True,
        "relations": selected,
    }


def integer_coordinates_from_relations(base_points, new_points, relations):
    base_count = len(base_points)
    new_count = len(new_points)
    if len(relations) != new_count:
        return None
    old_block = matrix(ZZ, [row[:base_count] for row in relations])
    new_block = matrix(ZZ, [row[base_count:] for row in relations])
    if new_block.nrows() != new_count or new_block.ncols() != new_count:
        return None
    if abs(new_block.det()) != 1:
        return None
    coordinates = -new_block.inverse() * old_block
    if any(entry not in ZZ for entry in coordinates.list()):
        return None
    print(
        f"{PROTOCOL}|stage=coordinate_verify|new_points={new_count}"
        f"|relation_block_det={new_block.det()}"
        f"|max_relation_coefficient={max(map(abs, old_block.list()+new_block.list()))}"
        f"|max_solved_coefficient={max(map(abs, coordinates.list()))}",
        flush=True,
    )
    curve = base_points[0].curve()
    for row, point in zip(coordinates.rows(), new_points):
        reconstructed = sum(
            (int(coefficient) * base_point for coefficient, base_point in zip(row, base_points)),
            curve(0),
        )
        if reconstructed != point:
            raise ArithmeticError("an integral quotient coordinate failed exact addition")
    return [list(map(int, row)) for row in coordinates.rows()]


def resolve_quotients(lineage, public, splits, isomorphisms, quotient_helper):
    public_by_id = {int(record["id"]): record for record in public["records"]}
    representative = lineage["representative"]
    ring = PolynomialRing(QQ, "u")
    quotient_rank_by_id = {
        int(record["curve_id"]): int(record["free_rank"])
        for record in lineage["exceptional_quotients"]
    }
    results = {}
    for curve_id in TARGET_IDS:
        records = splits[curve_id]
        if not records:
            results[curve_id] = {
                "split_count": 0,
                "displayed_quotient_rank": quotient_rank_by_id[curve_id],
                "class_span_rank": 0,
                "records": [],
            }
            continue
        public_record = public_by_id[curve_id]
        invariants, target_curve, target_points = public_short_points(public_record)
        target_iso = isomorphisms[curve_id]
        parameter = QQ(target_iso["parameter"])
        scale_q = QQ(target_iso["q"])
        scale_s = QQ(target_iso["s_with_s_squared_q"])
        if scale_s**2 != scale_q:
            raise ArithmeticError("a target isomorphism has inconsistent scales")
        generic_points = []
        for section in representative["sections"]:
            x = ring([QQ(value) for value in section["representative_x_coefficients_low_to_high"]])(parameter)
            y = ring([QQ(value) for value in section["representative_y_coefficients_low_to_high"]])(parameter)
            generic_points.append(target_curve(scale_q * x, scale_s**3 * y))
        word_rows = lineage["generic_basis"]["word_matrix_rows"]
        for point, word in zip(generic_points, word_rows):
            expected = sum(
                (int(coefficient) * target_points[index] for index, coefficient in enumerate(word)),
                target_curve(0),
            )
            if point != expected:
                raise ArithmeticError("the generic section/public-point alignment changed")

        quotient_rank = quotient_rank_by_id[curve_id]
        quotient_points = target_points[17 : 17 + quotient_rank]
        base_points = generic_points + quotient_points
        new_points = []
        for record in records:
            x, y = map(QQ, record["positive_chart_point"])
            point = target_curve(scale_q * x, scale_s**3 * y)
            new_points.append(point)
            record["positive_target_short_point"] = point_text(point)
            negative_x, negative_y = map(QQ, record["negative_chart_point"])
            negative = target_curve(scale_q * negative_x, scale_s**3 * negative_y)
            record["negative_target_short_point"] = point_text(negative)

        short_model = tuple(
            map(fraction, (0, 0, 0, invariants["A"], invariants["B"]))
        )
        all_points = base_points + new_points
        base_count = len(base_points)
        exact_relation_record = None
        integral_coordinates = None
        relation_error = None
        try:
            exact_relation_record = discover_relations(
                short_model, all_points, base_count, quotient_helper
            )
            integral_coordinates = integer_coordinates_from_relations(
                base_points, new_points, exact_relation_record["relations"]
            )
            if integral_coordinates is not None:
                exact_relation_record[
                    "all_relations_verified_by_exact_group_addition"
                ] = True
                exact_relation_record[
                    "verification_method"
                ] = (
                    "solve the unimodular new-point block, then verify every solved "
                    "point as an exact integral combination of the independent base"
                )
        except (ArithmeticError, RuntimeError, TimeoutError) as error:
            relation_error = str(error)

        finite_record = None
        if exact_relation_record is None:
            all_pairs = [(fraction(point[0]), fraction(point[1])) for point in all_points]
            joint = quotient_helper["build_finite_quotient_certificate"](
                short_model, all_pairs, relation_prime=2, prime_bound=1000
            )
            rows = quotient_helper["signature_rows"](joint)
            base_rank, _pivots = quotient_helper["matrix_rank_and_pivots_mod_prime"](
                [row[:base_count] for row in rows], base_count, 2
            )
            if base_rank != base_count:
                raise ArithmeticError("the exact known displayed basis lost finite-quotient rank")
            joint_rank = int(joint["combined_rank_over_relation_field"])
            finite_record = {
                "base_rank": base_rank,
                "rank_after_adjoining_split_points": joint_rank,
                "gain": joint_rank - base_count,
                "certificate_primes": joint["certificate_primes"],
                "stacked_rows_sha256": quotient_helper["rows_sha256"](rows),
            }

        exact_quotient_columns = []
        for index, record in enumerate(records):
            if integral_coordinates is None:
                record["exact_displayed_free_quotient_class"] = {
                    "status": "UNKNOWN_NOT_PROVED_IN_DISPLAYED_SUBGROUP",
                    "basis": [f"P{entry}" for entry in range(18, 18 + quotient_rank)],
                    "coordinates": None,
                }
            else:
                coordinates = integral_coordinates[index]
                quotient_coordinates = coordinates[17:]
                exact_quotient_columns.append(quotient_coordinates)
                record["exact_displayed_free_quotient_class"] = {
                    "status": "PROVED_EXACT_IN_DISPLAYED_SUBGROUP",
                    "basis": [f"P{entry}" for entry in range(18, 18 + quotient_rank)],
                    "coordinates": quotient_coordinates,
                    "generic_correction_coordinates": coordinates[:17],
                }
                record["finite_quotient_class_modulo_generic_17"] = {
                    "generic_correction_over_f2": [value % 2 for value in coordinates[:17]],
                    "displayed_quotient_coordinates_over_f2": [
                        value % 2 for value in quotient_coordinates
                    ],
                    "derived_from_exact_integral_relation": True,
                }

        class_span = (
            matrix(QQ, exact_quotient_columns).rank() if exact_quotient_columns else None
        )
        results[curve_id] = {
            "split_count": len(records),
            "displayed_quotient_rank": quotient_rank,
            "finite_quotient_fallback": finite_record,
            "exact_relation_closeout": exact_relation_record,
            "relation_closeout_error": relation_error,
            "all_split_points_proved_in_displayed_subgroup": integral_coordinates is not None,
            "class_span_rank": int(class_span) if class_span is not None else None,
            "records": records,
        }
        print(
            f"{PROTOCOL}|stage=quotient|curve={curve_id}|splits={len(records)}"
            f"|relations={exact_relation_record is not None}"
            f"|exact_span={results[curve_id]['class_span_rank']}",
            flush=True,
        )
    return results


def build_payload(checkpoint: Path, resume: bool):
    lineage = json.loads(LINEAGE.read_text())
    public = json.loads(PUBLIC.read_text())
    if lineage["status"] != "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        raise ArithmeticError("the 074d9 lineage certificate is not exact")
    if lineage["representative"]["chart"] != CHART:
        raise ArithmeticError("the lineage representative changed")
    priority_rows = parse_priority_rows(PRIORITY)
    representative_gram, pinned_gram, isometry, words = representative_words(
        lineage, priority_rows
    )
    parameters, isomorphisms = target_parameters(lineage)
    chord_helper = runpy.run_path(str(CHORD_HELPER))
    quotient_helper = runpy.run_path(str(QUOTIENT_HELPER))
    input_hashes = {
        relative(path): digest(path)
        for path in (
            LINEAGE,
            PUBLIC,
            PRIORITY,
            PUBLISHED_ATLAS,
            PINNED_GRAM,
            CHORD_HELPER,
            QUOTIENT_HELPER,
        )
    }
    input_key = hashlib.sha256(
        json.dumps(input_hashes, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    processed_primes, obstruction, validation = modular_sieve(
        lineage=lineage,
        words=words,
        parameters=parameters,
        chord_helper=chord_helper,
        input_key=input_key,
        checkpoint=checkpoint,
        resume=resume,
    )
    splits, exact_nonsquare, exact_global_validations = exact_survivors(
        lineage=lineage,
        words=words,
        priority_rows=priority_rows,
        parameters=parameters,
        obstruction=obstruction,
        validation=validation,
        chord_helper=chord_helper,
    )
    unresolved = {
        curve_id: [
            index
            for index, prime in enumerate(obstruction[curve_id])
            if not prime
        ]
        for curve_id in TARGET_IDS
    }
    classified = {
        curve_id: len(exact_nonsquare[curve_id]) + len(splits[curve_id])
        for curve_id in TARGET_IDS
    }
    if any(classified[curve_id] != len(unresolved[curve_id]) for curve_id in TARGET_IDS):
        raise ArithmeticError("the exact survivor decisions are incomplete")
    if any(not value for value in validation) and len(exact_global_validations) != sum(
        not value for value in validation
    ):
        raise ArithmeticError("the global cover validations are incomplete")

    quotient_results = resolve_quotients(
        lineage, public, splits, isomorphisms, quotient_helper
    )
    split_sets = {
        curve_id: {
            int(record["priority_index_zero_based"])
            for record in quotient_results[curve_id]["records"]
        }
        for curve_id in TARGET_IDS
    }
    cross_record = sorted(split_sets[356] & split_sets[385])
    all_exact_classes = all(
        result["all_split_points_proved_in_displayed_subgroup"]
        for result in quotient_results.values()
    )
    histograms = {
        curve_id: {
            str(prime): count
            for prime, count in sorted(Counter(obstruction[curve_id]).items())
            if prime
        }
        for curve_id in TARGET_IDS
    }
    return {
        "schema": "elkies-k3.r17-074d9-cross-fibre-bisection-transfer.v1",
        "status": STATUS,
        "claim": (
            "Complete exact cross-fibre split census for the 39,120 rigid norm-ten "
            "published-R17 bisection characters reconstructed on the 074d9 chart at "
            "ICARM curves 351, 356, 376, 377, and 385."
        ),
        "claim_boundary": (
            "A modular nonresidue proves a rational value is not a square. All modular "
            "calculations used as witnesses also certify a degree-two squarefree branch "
            "polynomial coprime to the reduced surface discriminant. Every modular "
            "survivor is decided over QQ. Exact quotient coordinates refer only to the "
            "free quotient of the subgroup generated by the displayed public points; "
            "they are not a full Mordell--Weil quotient or a rank upper bound."
        ),
        "source_coordinate_correction": {
            "original_published_equation": (
                "The five targets have no rational j-preimage on the original published "
                "equation, so its stored q(t) polynomials cannot be directly specialized."
            ),
            "tested_chart": CHART,
            "transport": (
                "all norm-ten trace characters are transported through the displayed "
                "integral lattice isometry and their residual quadratics are reconstructed "
                "on 074d9"
            ),
        },
        "lattice_transport": {
            "record_count": len(words),
            "all_transported_words_integral_distinct_and_norm_ten": True,
            "matrix_convention": "M^T * G_pinned * M = G_074d9; w_074d9=M^-1*w_pinned",
            "isometry_matrix": [list(map(int, row)) for row in isometry.rows()],
            "pinned_gram_determinant": int(pinned_gram.det()),
            "representative_gram_determinant": int(representative_gram.det()),
        },
        "sieve": {
            "declared_primes": list(SIEVE_PRIMES),
            "processed_primes": processed_primes,
            "method": (
                "exact Proposition-F1 residual-chord arithmetic over GF(p); accept a "
                "nonresidue only on a degree-two squarefree branch coprime to Delta"
            ),
            "global_validation_primes_in_priority_order": validation,
            "exact_global_validation_indices_zero_based": exact_global_validations,
            "obstruction_primes_in_priority_order": {
                str(curve_id): obstruction[curve_id] for curve_id in TARGET_IDS
            },
            "obstruction_prime_histograms": {
                str(curve_id): histograms[curve_id] for curve_id in TARGET_IDS
            },
            "exact_nonsquare_survivor_indices_zero_based": {
                str(curve_id): exact_nonsquare[curve_id] for curve_id in TARGET_IDS
            },
        },
        "fibres": [
            {
                "curve_id": curve_id,
                "parameter": rational_text(parameters[curve_id]),
                **quotient_results[curve_id],
            }
            for curve_id in TARGET_IDS
        ],
        "cross_fibre_transfer": {
            "source_curve_356_split_count": len(split_sets[356]),
            "target_curve_385_split_count": len(split_sets[385]),
            "covers_splitting_at_both_records": len(cross_record),
            "shared_priority_indices_zero_based": cross_record,
            "reverse_direction_is_the_same_frozen_cover_intersection": True,
            "record_anchored_rank_at_least_18_cover_count": len(cross_record),
            "rank_at_least_19_two_character_compositum_available": len(cross_record) >= 2,
            "all_split_points_have_exact_displayed_quotient_coordinates": all_exact_classes,
        },
        "summary": {
            "target_ids": list(TARGET_IDS),
            "split_counts": [len(split_sets[curve_id]) for curve_id in TARGET_IDS],
            "split_class_span_ranks": [
                quotient_results[curve_id]["class_span_rank"] for curve_id in TARGET_IDS
            ],
            "cross_record_cover_count": len(cross_record),
        },
        "generation": {
            "command": (
                "sage -python elkies-k3/scripts/"
                "certify_r17_074d9_cross_fibre_bisection_transfer.sage"
            ),
            "checker_sha256": digest(Path(__file__)),
            "inputs": input_hashes,
            "sieve_checkpoint": relative(checkpoint),
            "sieve_checkpoint_is_proof_input": False,
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact rational and finite-field function fields",
                "exact elliptic-curve group law",
                "PARI qfisom, ellheightmatrix, and qflllgram",
            ],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--checkpoint", type=Path, default=CHECKPOINT)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload(arguments.checkpoint, arguments.resume and not arguments.check)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit("stale 074d9 cross-fibre bisection certificate")
        terminal = "PASS"
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        terminal = "WROTE"
    print(
        f"{PROTOCOL}|splits={','.join(map(str, payload['summary']['split_counts']))}"
        f"|spans={','.join(map(str, payload['summary']['split_class_span_ranks']))}"
        f"|cross_record={payload['summary']['cross_record_cover_count']}"
        f"|status={terminal}|output={arguments.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
