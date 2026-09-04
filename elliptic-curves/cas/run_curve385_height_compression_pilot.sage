#!/usr/bin/env sage -python
"""Build and run a small height-compression chart pilot on curve 385 at M29.

The build phase never calls ``hyperellratpoints``.  Starting from the frozen
M29 primary ledger, it samples the complete 29-bit parity space with
deterministic bit-flip ascent on the current canonical-height form.  It then
builds and reduces the quartics for a small stable pool, calibrates each
reduced coordinate on the signed M29 basis, and commits a diverse chart order.

The search phase is a separate invocation bound to the committed protocol.
It uses the historical height/time/stack budget and classifies every returned
point against M29 exactly.  A miss remains a bounded pilot outcome.
"""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import gzip
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
import math
from math import comb
from pathlib import Path
import statistics
import sys
from typing import Any, Iterable, Sequence

from sage.all import pari


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
ART = ROOT / "artifacts/generated-results/elliptic-curves"
SOURCE = ART / "curve385_sparse_quotient_rank32_primary_ledger_v1.json.gz"
PROTOCOL = ART / "curve385_height_compression_pilot_protocol_v1.json"
RESULT = ART / "curve385_height_compression_pilot_blind_v1.json"
LEGACY_SOURCE = CAS / "run_curve385_iterated_half_lattice_search.sage"
ENGINE_SOURCE = CAS / "half_lattice_fake_descent_replay.sage"

EXPECTED_SOURCE_SHA256 = "08a2e416255910f733ef98283332e3a60a947350646329e4e2045cbc08d802c0"
DIMENSION = 29
SEED_COUNT = 256
STABLE_REDUCED_POOL_SIZE = 32
COMPRESSION_SHORTLIST_SIZE = 24
SEARCH_CHART_COUNT = 16
AUDIT_SCALE = 100_000
OPERATIVE_SCALE = 1_000_000
HEIGHT_BOUND = 100_000
TIMEOUT_SECONDS = 15.0
STACK_BYTES = 1_000_000_000
RELATION_CHUNK_SIZE = 64
RELATION_TIMEOUT_SECONDS = 180.0
SEED_DOMAIN = "curve385-height-compression-pilot-v1/full-M29-parity-seed"

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]
legacy = SourceFileLoader("curve385_height_compression_legacy", str(LEGACY_SOURCE)).load_module()


Point = tuple[Fraction, Fraction]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def read_source() -> tuple[dict[str, Any], tuple[Fraction, ...], tuple[Point, ...]]:
    if digest(SOURCE) != EXPECTED_SOURCE_SHA256:
        raise ArithmeticError("the frozen curve-385 M29 primary ledger changed")
    with gzip.open(SOURCE, "rt") as handle:
        source = json.load(handle)
    if source.get("status") != "STOPPED_AFTER_DECLARED_SPARSE_STAGE_LIMIT":
        raise ArithmeticError("the source ledger is not the completed M29 primary campaign")
    basis = tuple(legacy.read_point(item) for item in source["current_basis"])
    model = tuple(Fraction(item) for item in source["curve"]["short_model"])
    if len(basis) != DIMENSION or int(source["stop"]["basis_rank"]) != DIMENSION:
        raise ArithmeticError("the source is not the expected rank-29 lattice state")
    return source, model, basis


def quadratic(gram: Sequence[Sequence[Decimal]], vector: Sequence[int]) -> Decimal:
    return sum(
        Decimal(vector[i]) * gram[i][j] * Decimal(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def residue(mask: int) -> tuple[int, ...]:
    return tuple((int(mask) >> index) & 1 for index in range(DIMENSION))


def primitive_matrix(matrix: Sequence[Sequence[int]]) -> list[list[int]]:
    entries = [int(item) for row in matrix for item in row]
    common = 0
    for item in entries:
        common = math.gcd(common, abs(item))
    if common:
        entries = [item // common for item in entries]
    first = next((item for item in entries if item), 1)
    if first < 0:
        entries = [-item for item in entries]
    return [entries[:2], entries[2:]]


def matrix_product(
    left: Sequence[Sequence[int]], right: Sequence[Sequence[int]]
) -> list[list[int]]:
    return [
        [sum(int(left[i][k]) * int(right[k][j]) for k in range(2)) for j in range(2)]
        for i in range(2)
    ]


def polynomial(coefficients: Sequence[int]) -> str:
    return "+".join(
        f"({int(value)})*x^{index}" for index, value in enumerate(coefficients)
    ) or "0"


def polynomial_coefficients(item: Any, length: int = 5) -> list[int]:
    values = [int(item.polcoef(index)) for index in range(length)]
    while values and values[-1] == 0:
        values.pop()
    return values


def horizontal_matrix(transformation: Any) -> list[list[int]]:
    matrix = transformation[1]
    return [[int(matrix[i, j]) for j in range(2)] for i in range(2)]


def reduce_chart(model: Sequence[Fraction], base_point: Point) -> dict[str, Any]:
    cover = legacy.engine.alternate_cover(model, base_point)
    denominator = 1
    for coefficient in cover.coefficients:
        denominator = math.lcm(denominator, Fraction(coefficient).denominator)
    integral = [int(Fraction(value) * denominator * denominator) for value in cover.coefficients]
    result = pari(
        "my(m1,m2,C0=["
        + polynomial(integral)
        + ",0],C1,C2);C1=hyperellminimalmodel(C0,&m1);"
        + "C2=hyperellred(C1,&m2);[C2,m1,m2]"
    )
    reduced, first, second = result[0], result[1], result[2]
    reduced_p = polynomial_coefficients(reduced[0])
    reduced_q = polynomial_coefficients(reduced[1])
    composed = primitive_matrix(
        matrix_product(horizontal_matrix(first), horizontal_matrix(second))
    )
    determinant = composed[0][0] * composed[1][1] - composed[0][1] * composed[1][0]
    if determinant == 0:
        raise ArithmeticError("a reduced horizontal coordinate map is singular")
    return {
        "raw_quartic_coefficients_ascending": [str(value) for value in cover.coefficients],
        "raw_rational_coefficient_bits": max(legacy.engine.bit_height(value) for value in cover.coefficients),
        "denominator_clearing_factor_bits": denominator.bit_length(),
        "integral_coefficient_bits": max(abs(value).bit_length() for value in integral),
        "reduced_model": {
            "P_coefficients_ascending": reduced_p,
            "Q_coefficients_ascending": reduced_q,
            "maximum_coefficient_bits": max(abs(value).bit_length() for value in reduced_p + reduced_q),
            "discriminant": str(pari(f"hyperelldisc([{polynomial(reduced_p)},{polynomial(reduced_q)}])")),
        },
        "horizontal_reduced_to_raw_matrix": composed,
        "horizontal_determinant": str(determinant),
        "horizontal_map_bits": max(abs(value).bit_length() for row in composed for value in row),
    }


def raw_coordinate_homogeneous(
    curve_a: Fraction, base_point: Point, target: Point
) -> tuple[int, int]:
    if target == base_point:
        return 1, 0
    x_value, y_value = target
    x_base, y_base = base_point
    if x_value == x_base:
        if y_value != -y_base or y_base == 0:
            raise ArithmeticError("unexpected exceptional point in the quartic chart")
        parameter = -(3 * x_base * x_base + curve_a) / (2 * y_base)
    else:
        parameter = (y_value + y_base) / (x_value - x_base)
    parameter = Fraction(parameter)
    return parameter.numerator, parameter.denominator


def reduced_coordinate_data(
    raw: tuple[int, int], matrix: Sequence[Sequence[int]]
) -> tuple[int, int, int, int]:
    t_numerator, t_denominator = raw
    a_value, b_value = map(int, matrix[0])
    c_value, d_value = map(int, matrix[1])
    numerator = d_value * t_numerator - b_value * t_denominator
    denominator = -c_value * t_numerator + a_value * t_denominator
    common = math.gcd(abs(numerator), abs(denominator))
    if not common:
        raise ArithmeticError("the inverse horizontal map vanished projectively")
    numerator //= common
    denominator //= common
    height = max(abs(numerator), abs(denominator))
    return numerator, denominator, height, common


def calibrate_chart(
    model: Sequence[Fraction],
    basis: Sequence[Point],
    gram: Sequence[Sequence[Decimal]],
    representative: Sequence[int],
    base_point: Point,
    matrix: Sequence[Sequence[int]],
) -> dict[str, Any]:
    log_heights = []
    distortions = []
    cancellation_bits = []
    visible = 0
    for index, basis_point in enumerate(basis):
        for sign in (1, -1):
            target = basis_point if sign == 1 else (basis_point[0], -basis_point[1])
            raw = raw_coordinate_homogeneous(Fraction(model[3]), base_point, target)
            unused_numerator, unused_denominator, height, cancellation = reduced_coordinate_data(
                raw, matrix
            )
            log_height = math.log(height)
            vector = [-int(value) for value in representative]
            vector[index] += 2 * sign
            centered_degree_two_height = float(quadratic(gram, vector) / Decimal(2))
            log_heights.append(log_height)
            distortions.append(log_height - centered_degree_two_height)
            cancellation_bits.append(cancellation.bit_length())
            visible += height <= HEIGHT_BOUND
    depth = float(quadratic(gram, representative) / Decimal(4))
    minimum_distortion = min(distortions)
    return {
        "calibration_set": "signed_M29_basis",
        "calibration_point_count": len(log_heights),
        "signed_basis_visible_at_search_bound": visible,
        "reduced_log_height": {
            "minimum": min(log_heights),
            "median": statistics.median(log_heights),
            "maximum": max(log_heights),
        },
        "height_distortion": {
            "definition": "log H(s(P))-hhat(2P-Q)/2",
            "minimum": minimum_distortion,
            "median": statistics.median(distortions),
            "maximum": max(distortions),
        },
        "horizontal_preimage_gcd_bits": {
            "minimum": min(cancellation_bits),
            "median": statistics.median(cancellation_bits),
            "maximum": max(cancellation_bits),
        },
        "estimated_old_point_exclusion_margin": (
            2 * depth + minimum_distortion - math.log(HEIGHT_BOUND)
        ),
        "margin_is_certified_global_bound": False,
    }


def deterministic_seed(index: int) -> int:
    encoded = f"{SEED_DOMAIN}/{index}".encode()
    return int.from_bytes(sha256(encoded).digest(), "big") & ((1 << DIMENSION) - 1)


def build_protocol() -> dict[str, Any]:
    source, model, basis = read_source()
    gram, asymmetry = legacy.canonical_height_gram(model, basis)
    oracles = {
        scale: legacy.CosetOracle(legacy.rounded_gram(gram, scale))
        for scale in (AUDIT_SCALE, OPERATIVE_SCALE)
    }
    cache: dict[tuple[int, int], tuple[Decimal, tuple[int, ...], float]] = {}

    def evaluate(mask: int, scale: int = OPERATIVE_SCALE):
        key = scale, int(mask)
        if key not in cache:
            unused_norm, representative, error = oracles[scale].solve(residue(mask))
            depth = quadratic(gram, representative) / Decimal(4)
            cache[key] = depth, representative, error
        return cache[key]

    def ascend(mask: int) -> tuple[int, int]:
        steps = 0
        while True:
            current_depth = evaluate(mask)[0]
            neighbors = [mask ^ (1 << index) for index in range(DIMENSION)]
            best = max(neighbors, key=lambda item: (evaluate(item)[0], -item))
            if evaluate(best)[0] <= current_depth:
                return mask, steps
            mask = best
            steps += 1

    maxima: dict[int, list[int]] = {}
    ascent_steps = []
    for index in range(SEED_COUNT):
        maximum, steps = ascend(deterministic_seed(index))
        maxima.setdefault(maximum, []).append(index)
        ascent_steps.append(steps)

    candidates = []
    for mask, source_seeds in maxima.items():
        operative_depth, operative_representative, operative_error = evaluate(mask)
        audit_depth, audit_representative, audit_error = evaluate(mask, AUDIT_SCALE)
        candidates.append(
            {
                "mask": mask,
                "hex": f"0x{mask:08x}",
                "source_seed_indices": source_seeds,
                "canonical_depth": str(operative_depth),
                "representative": list(operative_representative),
                "audit_representative": list(audit_representative),
                "representative_stable_between_scales": (
                    operative_representative == audit_representative
                ),
                "depth_difference_between_scales": str(abs(operative_depth - audit_depth)),
                "maximum_cvp_distance_error": max(operative_error, audit_error),
            }
        )
    candidates.sort(key=lambda row: (-Decimal(row["canonical_depth"]), row["mask"]))

    searched_keys = set(map(str, source["searched_base_point_keys"]))
    stable_pool = []
    for row in candidates:
        if not row["representative_stable_between_scales"]:
            continue
        representative = tuple(map(int, row["representative"]))
        base_point = legacy.exact_linear_combination(model[3], basis, representative)
        if base_point is None:
            continue
        base_key = legacy.point_identifier(base_point)
        if base_key in searched_keys:
            continue
        reduced = reduce_chart(model, base_point)
        calibration = calibrate_chart(
            model,
            basis,
            gram,
            representative,
            base_point,
            reduced["horizontal_reduced_to_raw_matrix"],
        )
        stable_pool.append(
            {
                **row,
                "base_point": legacy.point_record(base_point),
                "base_point_key": base_key,
                "quartic": reduced,
                "coordinate_calibration": calibration,
            }
        )
        if len(stable_pool) == STABLE_REDUCED_POOL_SIZE:
            break
    if len(stable_pool) != STABLE_REDUCED_POOL_SIZE:
        raise ArithmeticError("too few stable unsearched local maxima for the pilot pool")

    compression_order = sorted(
        stable_pool,
        key=lambda row: (
            row["coordinate_calibration"]["signed_basis_visible_at_search_bound"],
            -row["coordinate_calibration"]["estimated_old_point_exclusion_margin"],
            row["quartic"]["reduced_model"]["maximum_coefficient_bits"],
            -float(row["canonical_depth"]),
            row["mask"],
        ),
    )
    shortlist = compression_order[:COMPRESSION_SHORTLIST_SIZE]

    distance_cache: dict[int, float] = {0: 0.0}

    def torus_distance(left: int, right: int) -> float:
        difference = int(left) ^ int(right)
        if difference not in distance_cache:
            unused_norm, representative, unused_error = oracles[OPERATIVE_SCALE].solve(
                residue(difference)
            )
            distance_cache[difference] = float(quadratic(gram, representative) / Decimal(4))
        return distance_cache[difference]

    selected = [shortlist[0]]
    remaining = {int(row["mask"]): row for row in shortlist[1:]}
    while len(selected) < SEARCH_CHART_COUNT:
        chosen = max(
            remaining.values(),
            key=lambda row: (
                min(torus_distance(row["mask"], old["mask"]) for old in selected),
                row["coordinate_calibration"]["estimated_old_point_exclusion_margin"],
                float(row["canonical_depth"]),
                -row["mask"],
            ),
        )
        selected.append(chosen)
        del remaining[int(chosen["mask"])]

    selected_rows = []
    for priority, row in enumerate(selected, 1):
        selected_rows.append(
            {
                "priority": priority,
                "mask": row["mask"],
                "hex": row["hex"],
                "representative": row["representative"],
                "base_point": row["base_point"],
                "base_point_key": row["base_point_key"],
                "canonical_depth": row["canonical_depth"],
                "minimum_torus_distance_to_earlier_center": (
                    None
                    if priority == 1
                    else min(
                        torus_distance(row["mask"], earlier["mask"])
                        for earlier in selected[: priority - 1]
                    )
                ),
                "quartic": row["quartic"],
                "coordinate_calibration": row["coordinate_calibration"],
            }
        )

    return {
        "schema": "elliptic-curves.curve385-height-compression-pilot-protocol.v1",
        "status": "BUILT_OUTCOME_FREE_BOUNDED_PILOT_PROTOCOL",
        "input_hashes": {
            str(SOURCE.relative_to(ROOT)): digest(SOURCE),
            str(LEGACY_SOURCE.relative_to(ROOT)): digest(LEGACY_SOURCE),
            str(ENGINE_SOURCE.relative_to(ROOT)): digest(ENGINE_SOURCE),
            str(Path(__file__).resolve().relative_to(ROOT)): digest(Path(__file__).resolve()),
        },
        "source_state": {
            "curve": source["curve"],
            "basis_rank": len(basis),
            "basis_sha256": legacy.canonical_hash([legacy.point_record(item) for item in basis]),
            "prior_searched_base_point_count": len(searched_keys),
            "prior_campaign_status": source["status"],
        },
        "candidate_generation": {
            "universe": "all 2^29 parity classes of the frozen current M29 basis",
            "deterministic_seed_domain": SEED_DOMAIN,
            "seed_count": SEED_COUNT,
            "one_bit_strict_ascent": True,
            "distinct_local_maximum_count": len(maxima),
            "mean_ascent_steps": statistics.mean(ascent_steps),
            "maximum_ascent_steps": max(ascent_steps),
            "cvp_rounding_scales": [AUDIT_SCALE, OPERATIVE_SCALE],
            "canonical_height_maximum_asymmetry": str(asymmetry),
            "stable_reduced_pool_size": len(stable_pool),
            "stable_reduced_pool": stable_pool,
        },
        "selection": {
            "pre_diversity_shortlist_size": len(shortlist),
            "search_chart_count": len(selected_rows),
            "primary_score": (
                "signed-basis visibility count, then estimated old-point exclusion "
                "margin 2D+min(delta)-log(B)"
            ),
            "secondary_selection": (
                "maximin current-lattice torus distance inside the top compression shortlist"
            ),
            "selected_charts": selected_rows,
            "selected_order_sha256": canonical_hash(
                [(row["mask"], row["representative"]) for row in selected_rows]
            ),
        },
        "search_budget": {
            "reduced_coordinate_height_bound_each_quartic": HEIGHT_BOUND,
            "wall_timeout_seconds_each_quartic": TIMEOUT_SECONDS,
            "gp_stack_bytes_each_quartic": STACK_BYTES,
            "retries": 0,
            "checkpointing": "one immutable result record after every chart",
        },
        "claim_boundary": [
            "The candidate generator reads no new pilot-search outcome and calls no point search.",
            "The source M29 ledger is used only for the current curve/lattice and to exclude repeated base points.",
            "The signed-basis distortion minimum is a calibration statistic, not a certified global height bound.",
            "The 256-seed local-maxima sample is not a complete enumeration of M29/2M29.",
            "The selected charts form a bounded mechanism pilot, not a rank, saturation, or point-absence test.",
        ],
        "reproducing_commands": [
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/run_curve385_height_compression_pilot.sage --phase build",
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/run_curve385_height_compression_pilot.sage --phase search",
        ],
    }


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def search(protocol: dict[str, Any]) -> dict[str, Any]:
    source, model, basis = read_source()
    script_key = str(Path(__file__).resolve().relative_to(ROOT))
    if protocol.get("status") != "BUILT_OUTCOME_FREE_BOUNDED_PILOT_PROTOCOL":
        raise ArithmeticError("the pilot protocol is not frozen at its presearch boundary")
    if protocol["input_hashes"].get(script_key) != digest(Path(__file__).resolve()):
        raise ArithmeticError("the builder/runner changed after the protocol was frozen")
    if protocol["source_state"]["basis_sha256"] != legacy.canonical_hash(
        [legacy.point_record(item) for item in basis]
    ):
        raise ArithmeticError("the protocol is bound to a different M29 basis")

    records = []
    discoveries: dict[Point, set[str]] = {}
    partial = {
        "schema": "elliptic-curves.curve385-height-compression-pilot-blind.v1",
        "status": "PARTIAL_CHECKPOINT",
        "input_hashes": {
            str(SOURCE.relative_to(ROOT)): digest(SOURCE),
            str(PROTOCOL.relative_to(ROOT)): digest(PROTOCOL),
            script_key: digest(Path(__file__).resolve()),
        },
        "protocol_selected_order_sha256": protocol["selection"]["selected_order_sha256"],
        "curve": source["curve"],
        "basis_rank_before": len(basis),
        "cover_records": records,
        "discoveries": [],
        "claim_boundary": protocol["claim_boundary"],
    }
    for chart in protocol["selection"]["selected_charts"]:
        outcome = legacy.engine.run_quartic_search(
            mask=int(chart["mask"]),
            representative=tuple(map(int, chart["representative"])),
            short_model=model,
            generic_points=basis,
            height_bound=HEIGHT_BOUND,
            timeout_seconds=TIMEOUT_SECONDS,
            stack_bytes=STACK_BYTES,
        )
        if outcome.record.get("base_point") != chart["base_point"]:
            raise ArithmeticError("the searched chart base point changed after commitment")
        source_label = f"height-compression-pilot:priority:{chart['priority']}"
        for item in outcome.curve_points:
            item = legacy.canonical_point(item)
            discoveries.setdefault(item, set()).add(source_label)
        records.append(
            {
                "priority": chart["priority"],
                "mask": chart["mask"],
                "canonical_depth": chart["canonical_depth"],
                "coordinate_calibration": chart["coordinate_calibration"],
                "search": outcome.record,
            }
        )
        partial["discoveries"] = legacy.discovery_records(discoveries)
        write_json(RESULT, partial)
        print(
            f"C385COMPRESSION|priority={chart['priority']}/{SEARCH_CHART_COUNT}|"
            f"mask={chart['hex']}|status={outcome.record['status']}|"
            f"points={len(outcome.curve_points)}",
            flush=True,
        )

    final_basis, classification = legacy.classify_discovered_group(
        model=model,
        basis=basis,
        discoveries=discoveries,
        relation_chunk_size=RELATION_CHUNK_SIZE,
        relation_timeout_seconds=RELATION_TIMEOUT_SECONDS,
        stack_bytes=STACK_BYTES,
    )
    partial["classification"] = classification
    partial["basis_rank_after"] = len(final_basis)
    partial["finite_point_occurrence_count"] = sum(
        len(row["search"].get("finite_curve_points", [])) for row in records
    )
    partial["distinct_returned_point_count"] = len(discoveries)
    partial["timeout_count"] = sum(
        row["search"]["status"] == "bounded_search_timeout" for row in records
    )
    partial["pari_failure_count"] = sum(
        row["search"]["status"] == "pari_failure" for row in records
    )
    prior_records = [
        row
        for state in source["lattice_states"]
        for stage in state["stages"]
        for row in stage["cover_records"]
    ]
    prior_hits = sum(bool(row["search"].get("finite_curve_points")) for row in prior_records)
    pilot_hits = sum(bool(row["search"].get("finite_curve_points")) for row in records)
    total_charts = len(prior_records) + len(records)
    total_hits = prior_hits + pilot_hits
    lower_tail_numerator = sum(
        comb(total_hits, count) * comb(total_charts - total_hits, len(records) - count)
        for count in range(pilot_hits + 1)
        if count <= total_hits and len(records) - count <= total_charts - total_hits
    )
    lower_tail_denominator = comb(total_charts, len(records))
    partial["historical_old_point_hit_comparison"] = {
        "prior_natural_weight_one_two": {
            "chart_count": len(prior_records),
            "chart_with_finite_point_count": prior_hits,
        },
        "height_compression_pilot": {
            "chart_count": len(records),
            "chart_with_finite_point_count": pilot_hits,
        },
        "descriptive_fixed_margin_fisher_lower_tail": {
            "numerator": str(lower_tail_numerator),
            "denominator": str(lower_tail_denominator),
            "decimal": lower_tail_numerator / lower_tail_denominator,
        },
        "inferential_status": (
            "descriptive only: the pilot and historical charts were selected by "
            "different deterministic policies and were not randomized arms"
        ),
    }
    if partial["timeout_count"] or partial["pari_failure_count"]:
        partial["status"] = "INCOMPLETE_BOUNDED_PILOT"
    elif classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        partial["status"] = "UNKNOWN_UNCLASSIFIED_RETURNED_POINTS"
    elif len(final_basis) > len(basis):
        partial["status"] = "PASS_NEW_INDEPENDENT_DIRECTION"
    elif classification["events"]:
        partial["status"] = "PASS_FINITE_INDEX_ENLARGEMENT"
    else:
        partial["status"] = "COMPLETE_BOUNDED_NO_GROUP_GROWTH"
    partial["interpretation"] = (
        "This is an operational test of fresh current-lattice holes and reduced-coordinate "
        "calibration. A no-growth result is a bounded miss and does not refute the exact "
        "midpoint identity or imply saturation."
    )
    write_json(RESULT, partial)
    return partial


def verify_existing() -> None:
    if not PROTOCOL.is_file() or not RESULT.is_file():
        raise SystemExit("the pilot protocol and result must both exist")
    protocol = json.loads(PROTOCOL.read_text())
    result = json.loads(RESULT.read_text())
    if result["input_hashes"][str(PROTOCOL.relative_to(ROOT))] != digest(PROTOCOL):
        raise ArithmeticError("the pilot result is bound to a different protocol")
    if result["protocol_selected_order_sha256"] != protocol["selection"][
        "selected_order_sha256"
    ]:
        raise ArithmeticError("the searched pilot order differs from the committed order")
    if len(result["cover_records"]) != protocol["selection"]["search_chart_count"]:
        raise ArithmeticError("the pilot result is incomplete")
    print(
        "C385COMPRESSION|status=VERIFIED|"
        f"outcome={result['status']}|charts={len(result['cover_records'])}|"
        f"rank={result['basis_rank_before']}->{result['basis_rank_after']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("build", "search", "verify"), default="verify")
    args = parser.parse_args()
    if args.phase == "build":
        payload = build_protocol()
        write_json(PROTOCOL, payload)
        print(
            "C385COMPRESSION|status=BUILT|"
            f"pool={payload['candidate_generation']['stable_reduced_pool_size']}|"
            f"selected={payload['selection']['search_chart_count']}"
        )
    elif args.phase == "search":
        if not PROTOCOL.is_file():
            raise SystemExit("build the pilot protocol before searching")
        if RESULT.exists():
            raise SystemExit("refusing to overwrite the existing pilot result")
        payload = search(json.loads(PROTOCOL.read_text()))
        print(
            f"C385COMPRESSION|status={payload['status']}|"
            f"rank={payload['basis_rank_before']}->{payload['basis_rank_after']}"
        )
    else:
        verify_existing()


if __name__ == "__main__":
    main()
