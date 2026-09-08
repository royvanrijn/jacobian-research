#!/usr/bin/env sage-python
"""Run the checkpointed MW17-jump-v2 exact quotient-rank campaign.

Each worker specializes one frozen candidate, verifies the generic MW17, runs
the calibrated generic-depth 43-chart phase, and (after any certified initial
escape) the adaptive 301-chart phase.  Only exact independent directions are
scored.  Chunk files are resumable and source population files are read-only.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import lcm
from pathlib import Path
import platform
import resource
import runpy
import shutil
import subprocess
import sys
import time
from typing import Any

from sage.all import EllipticCurve, Matrix, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = ROOT / "artifacts/generated-results/elkies-k3-mw17-jump-v2-campaign-v1.json"
CHUNK_DIR = ROOT / "artifacts/local/elkies-k3/mw17-jump-v2-pointed-v1"
LEDGER = ROOT / "artifacts/generated-results/elkies-k3-mw17-jump-v2-pointed-ledger-v1.json"
STOP_SENTINEL = CHUNK_DIR / "STOP_GAIN15.json"
LEGACY = ROOT / "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage"
ENGINE_SOURCE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
POLICY_SOURCE = ROOT / "elliptic-curves/cas/half_lattice_chart_policy.py"
LADDER_V1 = ROOT / "elliptic-curves/cas/run_r17_refresh_jump_ladder_blind.sage"
FAMILY_074D9 = ROOT / "elkies-k3/scripts/audit_r17_prospective_crt_local_stability.sage"
LINEAGE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
ATLAS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-record-lineage-atlas-v1.json"
DIRECT = {
    family: ROOT / f"artifacts/generated-results/elkies-k3-r17-norm12-orbit{family}-direct-fibration-v1.json"
    for family in ("07ca9", "08234", "08f72", "11952", "103b2")
}

GENERIC_DIMENSION = 17
INITIAL_CHARTS = 43
ADAPTIVE_CHARTS = 301
HEIGHT_BOUND = 100_000
CHART_TIMEOUT_SECONDS = 15.0
STACK_BYTES = 1_000_000_000
RELATION_CHUNK_SIZE = 64
RELATION_TIMEOUT_SECONDS = 180.0
EXPECTED_STATUS = "FROZEN_SELECTED_IMMUTABLE_POPULATIONS_BEFORE_V2_EVALUATION"

GENERIC_WORDS_074D9 = (
    ((2, 1),), ((3, 1),), ((4, 1),), ((5, 1),), ((8, 1),),
    ((11, 1),), ((13, 1),), ((15, 1),), ((16, 1),), ((17, 1),),
    ((1, 1), (2, -1)), ((1, 1), (6, -1)), ((1, 1), (7, -1)),
    ((1, 1), (9, -1)), ((1, 1), (10, -1)), ((1, 1), (12, -1)),
    ((1, 1), (14, -1)),
)

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(ROOT / "elliptic-curves/cas")]


from pointed_quartic_migration import validate_frozen_sources, runtime_search, require_runtime
from research_runtime.supervisor import captured_run, Limits

def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def atomic_write(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def cpu_clock() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def as_fraction(value) -> Fraction:
    return Fraction(str(value))


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def polynomial_value(coefficients, parameter):
    parameter = QQ(parameter)
    result = QQ(0)
    for coefficient in reversed(coefficients):
        result = result * parameter + QQ(coefficient)
    return result


def rational_function_value(record, parameter):
    numerator = polynomial_value(record["numerator_coefficients_low_to_high"], parameter)
    denominator = polynomial_value(record["denominator_coefficients_low_to_high"], parameter)
    if denominator == 0:
        raise ArithmeticError("a selected parameter is a pole of a generic section")
    return numerator / denominator


def load_campaign() -> dict[str, Any]:
    campaign = json.loads(CAMPAIGN.read_text())
    if campaign.get("status") != EXPECTED_STATUS or campaign.get("candidate_count") != 2_239:
        raise ArithmeticError("the MW17-jump-v2 campaign is not the frozen 2,239-row manifest")
    if campaign["candidate_list_sha256"] != canonical_hash(campaign["rows"]):
        raise ArithmeticError("the campaign candidate-list commitment changed")
    definition = {key: value for key, value in campaign.items() if key != "protocol_definition_sha256"}
    if campaign["protocol_definition_sha256"] != canonical_hash(definition):
        raise ArithmeticError("the campaign protocol definition does not replay")
    for name, expected in campaign["immutability"]["source_file_sha256"].items():
        path = ROOT / name
        if digest(path) != expected:
            raise ArithmeticError(f"immutable campaign source changed: {name}")
    validate_frozen_sources(campaign["implementation_hashes"])
    detector = campaign["detector"]
    expected_budget = (
        detector["initial_chart_count"],
        detector["adaptive_chart_count_after_nonzero_certified_initial_gain"],
        detector["height_bound_each_chart"],
        detector["wall_timeout_seconds_each_chart"],
        detector["gp_stack_bytes_each_chart"],
        detector["relation_chunk_size"],
        detector["relation_timeout_seconds"],
    )
    if expected_budget != (
        INITIAL_CHARTS, ADAPTIVE_CHARTS, HEIGHT_BOUND, CHART_TIMEOUT_SECONDS,
        STACK_BYTES, RELATION_CHUNK_SIZE, RELATION_TIMEOUT_SECONDS,
    ):
        raise ArithmeticError("the executable detector budget differs from the frozen campaign")
    return campaign


def generic_gram_074d9():
    lineage = json.loads(LINEAGE.read_text())
    word_gram = Matrix(ZZ, lineage["generic_basis"]["height_gram"])
    word_matrix = Matrix(
        ZZ,
        17,
        17,
        lambda row, column: next(
            (coefficient for index, coefficient in GENERIC_WORDS_074D9[row] if index - 1 == column),
            0,
        ),
    )
    if abs(word_matrix.det()) != 1:
        raise ArithmeticError("the 074d9 generic basis change stopped being unimodular")
    inverse = word_matrix.inverse()
    gram = Matrix(ZZ, inverse * word_gram * inverse.transpose())
    if gram.det() != 948:
        raise ArithmeticError("the 074d9 generic determinant changed")
    return tuple(tuple(int(value) for value in row) for row in gram.rows())


class Families:
    def __init__(self):
        self.family_074d9 = runpy.run_path(str(FAMILY_074D9))["Family"]()
        self.direct = {family: json.loads(path.read_text()) for family, path in DIRECT.items()}
        for family, payload in self.direct.items():
            if payload.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
                raise ArithmeticError(f"direct family {family} is not exact")
            if payload["sections"]["rank"] != 17 or payload["sections"]["height_gram_determinant"] != 948:
                raise ArithmeticError(f"direct family {family} lost its saturated MW17 certificate")
        atlas = json.loads(ATLAS.read_text())
        member = next(
            member
            for equivalence_class in atlas["atlas"]["pgl2_equivalence_classes"]
            if equivalence_class["representative"] == "norm12-orbit-0e80b"
            for member in equivalence_class["members"]
            if member["label"] == "norm12-orbit-103b2"
        )
        if not member["solve_certificate"]["identity_verified"]:
            raise ArithmeticError("the 0e80b-to-103b2 PGL2 map is not exact")
        self.map_0e80b_to_103b2 = tuple(map(ZZ, member["representative_to_member_pgl2_matrix_a_b_c_d"]))
        self.gram_074d9 = generic_gram_074d9()

    def specialize(self, row: dict[str, Any]):
        source_family = row["family"]
        direct_family = row["direct_model_family"]
        source_parameter = QQ(row["parameter"])
        if source_family == "074d9":
            curve, points = self.family_074d9.specialize(source_parameter)
            return curve, tuple(points), self.gram_074d9, {
                "source_family": source_family,
                "direct_model_family": direct_family,
                "source_parameter": rational_text(source_parameter),
                "direct_parameter": rational_text(source_parameter),
                "base_transform": "identity",
            }
        direct_parameter = source_parameter
        transform = "identity"
        if source_family == "0e80b":
            a, b, c, d = self.map_0e80b_to_103b2
            denominator = c * source_parameter + d
            if denominator == 0:
                raise ArithmeticError("the selected 0e80b parameter maps to infinity on 103b2")
            direct_parameter = (a * source_parameter + b) / denominator
            transform = {
                "kind": "exact_representative_to_member_PGL2",
                "matrix_a_b_c_d": list(map(int, (a, b, c, d))),
                "identity_verified_in_atlas": True,
            }
        payload = self.direct[direct_family]
        model = payload["weierstrass_model"]
        coefficient_a = polynomial_value(model["A_coefficients_low_to_high"], direct_parameter)
        coefficient_b = polynomial_value(model["B_coefficients_low_to_high"], direct_parameter)
        curve = EllipticCurve(QQ, [coefficient_a, coefficient_b])
        if curve.discriminant() == 0:
            raise ArithmeticError("singular selected specialization")
        points = []
        for expected_index, section in enumerate(payload["sections"]["records"]):
            if int(section["basis_index"]) != expected_index:
                raise ArithmeticError("the direct generic basis order changed")
            x_coordinate = rational_function_value(section["X"], direct_parameter)
            y_coordinate = rational_function_value(section["Y"], direct_parameter)
            if y_coordinate**2 != x_coordinate**3 + coefficient_a * x_coordinate + coefficient_b:
                raise ArithmeticError(f"direct section {expected_index + 1} failed exact specialization")
            points.append(curve(x_coordinate, y_coordinate))
        gram = tuple(tuple(int(value) for value in line) for line in payload["sections"]["height_gram"])
        return curve, tuple(points), gram, {
            "source_family": source_family,
            "direct_model_family": direct_family,
            "source_parameter": rational_text(source_parameter),
            "direct_parameter": rational_text(direct_parameter),
            "base_transform": transform,
        }


def normalize_curve(curve, known, source_family: str):
    candidates = []

    # Direct family gauges can have enormous rational denominators.  Clearing
    # them on the already-short equation gives the smaller exact pointed
    # quartics used by the raw bounded search; choosing the apparently smaller
    # p=2 candidate by equation-coefficient bits instead exhausted GP memory on
    # the first four 07ca9 preflight fibres.  The established 074d9 path retains
    # its prior p=2 normalization.  Both exact candidates and sizes are recorded.
    integral = curve.integral_model()
    to_integral = curve.isomorphisms(integral)
    if not to_integral or any(integral.a_invariants()[index] != 0 for index in range(3)):
        raise ArithmeticError("the direct integral normalization is not an exact short model")
    integral_points = tuple(to_integral[0](point) for point in known)
    integral_model = tuple(as_fraction(value) for value in integral.a_invariants())
    candidates.append(("direct_short_integral_model", integral_model, integral_points))

    local_minimal = curve.local_data(2).minimal_model()
    isomorphisms = curve.isomorphisms(local_minimal)
    if not isomorphisms:
        raise ArithmeticError("no exact isomorphism to the deterministic p=2 minimal model")
    to_local = isomorphisms[0]
    local_points = tuple(to_local(point) for point in known)
    a1, a2, a3, a4, a6 = map(QQ, local_minimal.a_invariants())
    if any(value.denominator() != 1 for value in (a1, a2, a3, a4, a6)):
        raise ArithmeticError("the deterministic p=2 minimal model is not integral")
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    p2_model = tuple(as_fraction(value) for value in (0, 0, 0, -27 * c4, -54 * c6))
    p2_points = tuple(
        (
            as_fraction(36 * point[0] + 3 * b2),
            as_fraction(108 * (2 * point[1] + a1 * point[0] + a3)),
        )
        for point in local_points
    )
    candidates.append(("p2_minimal_then_short_invariants", p2_model, p2_points))

    def coefficient_bits(candidate_model):
        return max(
            max(abs(value.numerator).bit_length(), value.denominator.bit_length())
            for value in candidate_model
        )

    candidate_sizes = {name: coefficient_bits(candidate_model) for name, candidate_model, unused in candidates}
    selected_name = (
        "p2_minimal_then_short_invariants"
        if source_family == "074d9"
        else "direct_short_integral_model"
    )
    unused_name, search_model, search_points = next(item for item in candidates if item[0] == selected_name)
    certificate_model = p2_model
    certificate_points = p2_points
    search_curve = EllipticCurve(QQ, list(search_model))
    certificate_curve = EllipticCurve(QQ, list(certificate_model))
    certificate_point_objects = tuple(certificate_curve(point) for point in certificate_points)
    search_to_certificate_maps = search_curve.isomorphisms(certificate_curve)
    if not search_to_certificate_maps:
        raise ArithmeticError("the exact search and certificate short models are not isomorphic")
    search_to_certificate_map = next(
        (
            candidate
            for candidate in search_to_certificate_maps
            if tuple(candidate(point) for point in search_points) == certificate_point_objects
        ),
        None,
    )
    if search_to_certificate_map is None:
        raise ArithmeticError("no exact search-to-certificate isomorphism retains the ordered generic basis")
    certificate_to_search_map = ~search_to_certificate_map
    if len(search_points) != 17 or any(search_curve(point) == search_curve(0) for point in search_points):
        raise ArithmeticError("normalization did not retain seventeen nonzero generic points")
    def search_to_certificate(point):
        image = search_to_certificate_map(search_curve(QQ(point[0]), QQ(point[1])))
        return as_fraction(image[0]), as_fraction(image[1])

    def certificate_to_search(point):
        image = certificate_to_search_map(certificate_curve(QQ(point[0]), QQ(point[1])))
        return as_fraction(image[0]), as_fraction(image[1])

    return (
        certificate_model,
        tuple((as_fraction(point[0]), as_fraction(point[1])) for point in certificate_points),
        search_model,
        tuple((as_fraction(point[0]), as_fraction(point[1])) for point in search_points),
        search_to_certificate,
        certificate_to_search,
        {
        "source_model": [rational_text(value) for value in curve.a_invariants()],
        "p2_minimal_model": [rational_text(value) for value in local_minimal.a_invariants()],
        "short_integral_candidate_maximum_coefficient_bits": candidate_sizes["direct_short_integral_model"],
        "p2_minimal_then_short_candidate_maximum_coefficient_bits": candidate_sizes["p2_minimal_then_short_invariants"],
        "selected_search_normalization": selected_name,
        "certificate_normalization": "p2_minimal_then_short_invariants",
        "search_short_model": [str(value) for value in search_model],
        "search_short_model_sha256": canonical_hash([str(value) for value in search_model]),
        "certificate_short_model": [str(value) for value in certificate_model],
        "certificate_short_model_sha256": canonical_hash([str(value) for value in certificate_model]),
        "exact_search_certificate_isomorphism_verified_on_ordered_mw17": True,
        "all_seventeen_points_transported_exactly": True,
        },
    )


def signature_record(legacy, signatures, column_count: int):
    return {
        "column_count": column_count,
        "combined_rank": legacy.combined_mod2_rank(signatures, column_count),
        "prime_bound": legacy.CERTIFICATE_PRIME_BOUND,
        "signatures": [
            {
                "prime": row.prime,
                "group_order": row.group_order,
                "doubled_subgroup_order": row.doubled_subgroup_order,
                "quotient_dimension": row.quotient_dimension,
                "rows": [list(vector) for vector in row.rows],
            }
            for row in signatures
        ],
    }


def complete_generic_census(legacy, gram):
    oracle = legacy.CosetOracle(gram)
    rows = []
    maximum_error = 0.0
    for mask in range(1 << GENERIC_DIMENSION):
        residue = tuple((mask >> index) & 1 for index in range(GENERIC_DIMENSION))
        norm, representative, error = oracle.solve(residue)
        maximum_error = max(maximum_error, error)
        rows.append((norm, mask, representative))
    rows.sort(key=lambda row: (-row[0], row[1]))
    if len(rows) != 1 << GENERIC_DIMENSION:
        raise ArithmeticError("the complete generic half-lattice census is incomplete")
    return rows, maximum_error


def run_quartic_search_raw(engine, **kwargs):
    """Compatibility name; all active searches use PointedQuarticSearch."""
    from pointed_quartic_search import run_quartic_search
    return run_quartic_search(**kwargs)


class DetectorArgs:
    relation_chunk_size = RELATION_CHUNK_SIZE
    relation_timeout_seconds = RELATION_TIMEOUT_SECONDS
    stack_bytes = STACK_BYTES


def run_fibre(row: dict[str, Any], families: Families, modules) -> dict[str, Any]:
    ladder, legacy, engine, chart_policy = modules
    started_wall = time.monotonic()
    started_cpu = cpu_clock()
    curve, known, generic_gram, specialization = families.specialize(row)
    (
        model,
        generic,
        search_model,
        search_generic,
        search_to_certificate,
        certificate_to_search,
        normalization,
    ) = normalize_curve(curve, known, row["family"])
    signatures = legacy.find_mod2_reduction_certificate(
        model, generic, prime_bound=legacy.CERTIFICATE_PRIME_BOUND
    )
    if legacy.combined_mod2_rank(signatures, GENERIC_DIMENSION) != GENERIC_DIMENSION:
        raise ArithmeticError("the specialized generic MW17 lost exact independence")

    generic_rows, generic_cvp_error = complete_generic_census(legacy, generic_gram)
    initial_gram, initial_ranked, initial_ranking = ladder.rank_initial_charts(
        legacy, model, generic, generic_rows
    )
    universe = f"mw17-jump-v2:{row['sample_id']}:generic-depth43"
    chart_ids = [f"mask:{entry[1]:05x}" for entry in initial_ranked]
    initial_certificate = chart_policy.bind_ordering(
        basis_records=[legacy.point_record(point) for point in generic],
        height_gram_rows=initial_gram,
        generic_coordinate_rows=ladder.identity_rows(GENERIC_DIMENSION),
        quotient_coordinate_rows=[],
        chart_universe_id=universe,
        ordered_chart_ids=chart_ids,
        heuristics=["legacy_half_lattice_depth"],
    )
    chart_policy.validate_ordering(
        initial_certificate,
        basis_records=[legacy.point_record(point) for point in generic],
        height_gram_rows=initial_gram,
        generic_coordinate_rows=ladder.identity_rows(GENERIC_DIMENSION),
        quotient_coordinate_rows=[],
        chart_universe_id=universe,
        ordered_chart_ids=chart_ids,
    )

    discoveries = {}
    searched_keys = set()
    initial_covers = []
    for priority, (depth, mask, representative, generic_representative, generic_norm) in enumerate(initial_ranked, 1):
        base_point = legacy.exact_linear_combination(model[3], generic, representative)
        if base_point is None:
            raise ArithmeticError("an initial half-class produced infinity")
        base_key = legacy.point_identifier(base_point)
        if base_key in searched_keys:
            raise ArithmeticError("an initial pointed chart was duplicated")
        outcome = run_quartic_search_raw(
            engine,
            mask=mask,
            representative=representative,
            short_model=search_model,
            generic_points=search_generic,
            height_bound=HEIGHT_BOUND,
            timeout_seconds=CHART_TIMEOUT_SECONDS,
            stack_bytes=STACK_BYTES,
        )
        searched_keys.add(base_key)
        source = f"initial:priority:{priority}:mask:{mask:#07x}"
        for point in outcome.curve_points:
            certificate_point = search_to_certificate(point)
            discoveries.setdefault(legacy.canonical_point(certificate_point), set()).add(source)
        search_record = ladder.compact_search_record(outcome.record)
        search_record["error"] = outcome.record.get("error")
        initial_covers.append(
            {
                "priority": priority,
                "mask": mask,
                "exact_generic_norm": generic_norm,
                "generic_representative": list(generic_representative),
                "specialized_representative": list(representative),
                "specialized_depth": str(depth),
                "base_point_key": base_key,
                "search": search_record,
            }
        )
    basis, initial_classification = legacy.classify_discovered_group(
        model=model,
        basis=generic,
        discoveries=discoveries,
        relation_chunk_size=RELATION_CHUNK_SIZE,
        relation_timeout_seconds=RELATION_TIMEOUT_SECONDS,
        stack_bytes=STACK_BYTES,
    )
    if initial_classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        raise ArithmeticError("initial discoveries could not be exactly classified")
    initial_gain = len(basis) - GENERIC_DIMENSION

    adaptive = {
        "status": "NOT_APPLICABLE_ZERO_DISCOVERED_QUOTIENT",
        "cover_records": [],
        "exact_incremental_quotient_rank_recovered": 0,
    }
    if initial_gain:
        adaptive_gram, adaptive_ranked, generic_coordinates, complement, adaptive_ranking = ladder.rank_adaptive_pool(
            legacy, model, basis, generic, generic_rows, DetectorArgs()
        )
        adaptive_ids = [f"gmask:{entry[1]:05x}:qword:{entry[2]:x}" for entry in adaptive_ranked]
        adaptive_universe = f"mw17-jump-v2:{row['sample_id']}:adaptive301"
        adaptive_certificate = chart_policy.bind_ordering(
            basis_records=[legacy.point_record(point) for point in basis],
            height_gram_rows=adaptive_gram,
            generic_coordinate_rows=generic_coordinates,
            quotient_coordinate_rows=complement,
            chart_universe_id=adaptive_universe,
            ordered_chart_ids=adaptive_ids,
            heuristics=["legacy_half_lattice_depth", "quotient_hamming_weight"],
        )
        chart_policy.validate_ordering(
            adaptive_certificate,
            basis_records=[legacy.point_record(point) for point in basis],
            height_gram_rows=adaptive_gram,
            generic_coordinate_rows=generic_coordinates,
            quotient_coordinate_rows=complement,
            chart_universe_id=adaptive_universe,
            ordered_chart_ids=adaptive_ids,
        )
        adaptive_covers = []
        search_basis = tuple(certificate_to_search(point) for point in basis)
        for priority, (depth, generic_mask, quotient_word, residue, representative, generic_norm) in enumerate(adaptive_ranked, 1):
            base_point = legacy.exact_linear_combination(model[3], basis, representative)
            if base_point is None:
                raise ArithmeticError("an adaptive half-class produced infinity")
            base_key = legacy.point_identifier(base_point)
            if base_key in searched_keys:
                raise ArithmeticError("an adaptive pointed chart was already searched")
            mask = sum(int(bit) << index for index, bit in enumerate(residue))
            outcome = run_quartic_search_raw(
                engine,
                mask=mask,
                representative=representative,
                short_model=search_model,
                generic_points=search_basis,
                height_bound=HEIGHT_BOUND,
                timeout_seconds=CHART_TIMEOUT_SECONDS,
                stack_bytes=STACK_BYTES,
            )
            searched_keys.add(base_key)
            source = f"adaptive:priority:{priority}:gmask:{generic_mask:#07x}:qword:{quotient_word:#x}"
            for point in outcome.curve_points:
                certificate_point = search_to_certificate(point)
                discoveries.setdefault(legacy.canonical_point(certificate_point), set()).add(source)
            search_record = ladder.compact_search_record(outcome.record)
            search_record["error"] = outcome.record.get("error")
            adaptive_covers.append(
                {
                    "priority": priority,
                    "generic_mask": generic_mask,
                    "exact_generic_norm": generic_norm,
                    "quotient_word": quotient_word,
                    "current_basis_residue": list(residue),
                    "representative": list(representative),
                    "canonical_depth": str(depth),
                    "base_point_key": base_key,
                    "search": search_record,
                }
            )
        basis, final_classification = legacy.classify_discovered_group(
            model=model,
            basis=basis,
            discoveries=discoveries,
            relation_chunk_size=RELATION_CHUNK_SIZE,
            relation_timeout_seconds=RELATION_TIMEOUT_SECONDS,
            stack_bytes=STACK_BYTES,
        )
        if final_classification["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            raise ArithmeticError("adaptive discoveries could not be exactly classified")
        adaptive = {
            "status": "PASS_EXACTLY_CLASSIFIED",
            "ranking": adaptive_ranking,
            "ordering_certificate": adaptive_certificate,
            "cover_records": adaptive_covers,
            "discovered_group_classification": final_classification,
            "exact_incremental_quotient_rank_recovered": len(basis) - GENERIC_DIMENSION - initial_gain,
        }

    final_gain = len(basis) - GENERIC_DIMENSION
    all_covers = initial_covers + adaptive["cover_records"]
    backend_failures = sum(entry["search"]["status"] == "pari_failure" for entry in all_covers)
    timeouts = sum(entry["search"]["status"] == "bounded_search_timeout" for entry in all_covers)
    result = {
        "sample_id": row["sample_id"],
        "campaign_index": row["campaign_index"],
        "source_population": row["source_population"],
        "family": row["family"],
        "frame_class": row["frame_class"],
        "parameter": row["parameter"],
        "status": (
            "CENSORED_COVER_BACKEND_FAILURE"
            if backend_failures
            else "PASS_EXACT_CERTIFIED_QUOTIENT_GAIN"
        ),
        "specialization": specialization,
        "normalization": normalization,
        "generic_subgroup": {
            "rank": GENERIC_DIMENSION,
            "generic_height_gram": [list(line) for line in generic_gram],
            "generic_height_gram_determinant": int(Matrix(ZZ, generic_gram).det()),
            "finite_reduction_independence_certificate": signature_record(legacy, signatures, GENERIC_DIMENSION),
        },
        "generic_census": {
            "complete_class_count": 1 << GENERIC_DIMENSION,
            "maximum_exact_generic_norm": generic_rows[0][0],
            "deepest_class_count": sum(entry[0] == generic_rows[0][0] for entry in generic_rows),
            "maximum_cvp_distance_error": generic_cvp_error,
            "top301_masks_sha256": legacy.canonical_hash([entry[1] for entry in generic_rows[:ADAPTIVE_CHARTS]]),
        },
        "initial": {
            "status": "PASS_EXACTLY_CLASSIFIED",
            "ranking": initial_ranking,
            "ordering_certificate": initial_certificate,
            "cover_records": initial_covers,
            "discovered_group_classification": initial_classification,
            "initial_43_chart_certified_gain_diagnostic_only": initial_gain,
            "used_as_candidate_or_leaderboard_filter": False,
        },
        "adaptive": adaptive,
        "actual_certified_quotient_rank_gain": (
            final_gain if final_gain or not backend_failures else None
        ),
        "certified_rank_lower_bound": (
            GENERIC_DIMENSION + final_gain if final_gain or not backend_failures else None
        ),
        "attempted_chart_count": len(all_covers),
        "bounded_cover_timeout_count": timeouts,
        "cover_backend_failure_count": backend_failures,
        "final_basis": [legacy.point_record(point) for point in basis],
        "timing": {
            "worker_wall_seconds": time.monotonic() - started_wall,
            "worker_cpu_seconds_parent_plus_children": cpu_clock() - started_cpu,
        },
        "claim_boundary": (
            "The score is the exact independent rank gained in the discovered subgroup. "
            "It is a rank lower bound; misses are bounded by the frozen detector budget."
        ),
    }
    if final_gain >= 15:
        atomic_write(
            STOP_SENTINEL,
            {
                "schema": "elkies-k3.mw17-jump-v2-stop.v1",
                "status": "CERTIFIED_GAIN_AT_LEAST_15_GLOBAL_STOP",
                "sample_id": row["sample_id"],
                "campaign_index": row["campaign_index"],
                "actual_certified_quotient_rank_gain": final_gain,
                "certified_rank_lower_bound": 17 + final_gain,
                "campaign_sha256": digest(CAMPAIGN),
                "runtime_search": runtime_search(),
            },
        )
    return result


def load_modules():
    ladder = SourceFileLoader("mw17_jump_v2_ladder_helpers", str(LADDER_V1)).load_module()
    legacy = SourceFileLoader("mw17_jump_v2_legacy", str(LEGACY)).load_module()
    engine = SourceFileLoader("mw17_jump_v2_engine", str(ENGINE_SOURCE)).load_module()
    chart_policy = SourceFileLoader("mw17_jump_v2_policy", str(POLICY_SOURCE)).load_module()
    legacy.GENERIC_DIMENSION = GENERIC_DIMENSION
    legacy.OLD_CLASS_COUNT = INITIAL_CHARTS
    return ladder, legacy, engine, chart_policy


def failure_record(row, status: str, failure: Any, elapsed: float):
    return {
        "sample_id": row["sample_id"],
        "campaign_index": row["campaign_index"],
        "source_population": row["source_population"],
        "family": row["family"],
        "parameter": row["parameter"],
        "status": status,
        "failure": failure,
        "actual_certified_quotient_rank_gain": None,
        "certified_rank_lower_bound": None,
        "supervisor_wall_seconds": elapsed,
    }


def run_single(index: int) -> None:
    campaign = load_campaign()
    if not 0 <= index < len(campaign["rows"]):
        raise ValueError("single index is outside the frozen campaign")
    memory = campaign["resource_gate"]["worker_address_space_bytes"]
    if memory is not None:
        resource.setrlimit(resource.RLIMIT_AS, (memory, memory))
    result = run_fibre(campaign["rows"][index], Families(), load_modules())
    print("RESULT_JSON=" + canonical_text(result), flush=True)


def write_checkpoint(path: Path, campaign: dict[str, Any], chunk_index: int, chunk_count: int, indices, records):
    atomic_write(
        path,
        {
            "schema": "elkies-k3.mw17-jump-v2-chunk.v1",
            "status": "COMPLETE_CHUNK" if len(records) == len(indices) else "PARTIAL_CHECKPOINT",
            "campaign_sha256": digest(CAMPAIGN),
            "runtime_search": runtime_search(),
            "candidate_list_sha256": campaign["candidate_list_sha256"],
            "protocol_definition_sha256": campaign["protocol_definition_sha256"],
            "chunk_index": chunk_index,
            "chunk_count": chunk_count,
            "scheduled_indices": indices,
            "completed_record_count": len(records),
            "records": records,
        },
    )


def run_chunk(chunk_index: int, chunk_count: int, output: Path, max_new: int | None):
    campaign = load_campaign()
    if not 0 <= chunk_index < chunk_count:
        raise ValueError("chunk index must lie in [0, chunk count)")
    indices = [index for index in range(len(campaign["rows"])) if index % chunk_count == chunk_index]
    records = []
    if output.exists():
        old = json.loads(output.read_text())
        require_runtime(old)
        if (
            old.get("campaign_sha256") != digest(CAMPAIGN)
            or old.get("chunk_index") != chunk_index
            or old.get("chunk_count") != chunk_count
            or old.get("scheduled_indices") != indices
        ):
            raise ArithmeticError("an existing checkpoint belongs to another campaign or chunk layout")
        records = old["records"]
    completed_ids = {record["sample_id"] for record in records}
    new_count = 0
    timeout = campaign["resource_gate"]["worker_wall_timeout_seconds"]
    for position, index in enumerate(indices, 1):
        row = campaign["rows"][index]
        if row["sample_id"] in completed_ids:
            continue
        if STOP_SENTINEL.exists():
            print(f"MW17JUMPV2CHUNK|chunk={chunk_index}/{chunk_count}|status=STOP_GAIN15", flush=True)
            break
        if max_new is not None and new_count >= max_new:
            break
        command = [sys.executable, str(Path(__file__).resolve()), "--single-index", str(index), "--legacy-census-regression"]
        started = time.monotonic()
        try:
            completed = captured_run(
                command,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=False,
            )
            lines = [line for line in completed.stdout.splitlines() if line.strip()]
            result_line = next((line for line in reversed(lines) if line.startswith("RESULT_JSON=")), None)
            if completed.returncode == 0 and result_line is not None:
                result = json.loads(result_line[len("RESULT_JSON="):])
                result["supervisor_wall_seconds"] = time.monotonic() - started
            else:
                result = failure_record(
                    row,
                    "CENSORED_FIBRE_WORKER_FAILURE",
                    {"returncode": completed.returncode, "output_tail": lines[-40:]},
                    time.monotonic() - started,
                )
        except subprocess.TimeoutExpired as error:
            output_text = error.stdout or ""
            if isinstance(output_text, bytes):
                output_text = output_text.decode(errors="replace")
            result = failure_record(
                row,
                "CENSORED_FIBRE_WORKER_TIMEOUT",
                {"output_tail": output_text.splitlines()[-40:]},
                time.monotonic() - started,
            )
        records.append(result)
        new_count += 1
        write_checkpoint(output, campaign, chunk_index, chunk_count, indices, records)
        print(
            f"MW17JUMPV2CHUNK|chunk={chunk_index}/{chunk_count}|position={position}/{len(indices)}|"
            f"sample={row['sample_id']}|status={result['status']}|gain={result.get('actual_certified_quotient_rank_gain')}",
            flush=True,
        )
    write_checkpoint(output, campaign, chunk_index, chunk_count, indices, records)


def merge_chunks(chunk_dir: Path, chunk_count: int, output: Path) -> None:
    campaign = load_campaign()
    records_by_index = {}
    provenance = []
    for chunk_index in range(chunk_count):
        path = chunk_dir / f"chunk-{chunk_index:02d}-of-{chunk_count:02d}.json"
        if not path.exists():
            continue
        chunk = json.loads(path.read_text())
        require_runtime(chunk)
        if (
            chunk.get("campaign_sha256") != digest(CAMPAIGN)
            or chunk.get("candidate_list_sha256") != campaign["candidate_list_sha256"]
            or chunk.get("chunk_index") != chunk_index
            or chunk.get("chunk_count") != chunk_count
        ):
            raise ArithmeticError(f"chunk {chunk_index} is not bound to this campaign")
        for record in chunk["records"]:
            index = int(record["campaign_index"])
            if index in records_by_index:
                raise ArithmeticError("duplicate campaign index across chunks")
            records_by_index[index] = record
        provenance.append({"path": relative(path), "sha256": digest(path), "record_count": len(chunk["records"])})
    scored = [record for record in records_by_index.values() if record.get("actual_certified_quotient_rank_gain") is not None]
    leaderboard = sorted(scored, key=lambda record: (-record["actual_certified_quotient_rank_gain"], record["campaign_index"]))
    measurements = [
        {
            "leaderboard_position": position,
            "campaign_index": record["campaign_index"],
            "sample_id": record["sample_id"],
            "source_population": record["source_population"],
            "family": record["family"],
            "parameter": record["parameter"],
            "initial_43_chart_gain_diagnostic_only": record["initial"]["initial_43_chart_certified_gain_diagnostic_only"],
            "actual_certified_quotient_rank_gain": record["actual_certified_quotient_rank_gain"],
            "certified_rank_lower_bound": record["certified_rank_lower_bound"],
            "detailed_chunk_record_sha256": canonical_hash(record),
        }
        for position, record in enumerate(leaderboard, 1)
    ]
    complete = len(records_by_index) == campaign["candidate_count"]
    status_counts = Counter(record["status"] for record in records_by_index.values())
    document = {
        "schema": "elkies-k3.mw17-jump-v2-ledger.v1",
        "status": "COMPLETE_SCHEDULED_EVALUATION" if complete else "PARTIAL_CHECKPOINTED_EVALUATION",
        "campaign_sha256": digest(CAMPAIGN),
        "runtime_search": runtime_search(),
        "candidate_list_sha256": campaign["candidate_list_sha256"],
        "scheduled_candidate_count": campaign["candidate_count"],
        "completed_worker_count": len(records_by_index),
        "measured_certified_gain_count": len(measurements),
        "status_counts": dict(sorted(status_counts.items())),
        "measurements_ranked_only_by_actual_certified_quotient_rank_gain": measurements,
        "gain_at_least_15_found": any(row["actual_certified_quotient_rank_gain"] >= 15 for row in measurements),
        "chunk_provenance": provenance,
        "claim_boundary": campaign["claim_boundary"],
        "generation": {
            "script": relative(Path(__file__)),
            "python": platform.python_version(),
            "command": f"sage -python {relative(Path(__file__))} --merge --chunk-count {chunk_count}",
        },
    }
    atomic_write(output, document)
    print(
        f"MW17JUMPV2MERGE|status={document['status']}|completed={len(records_by_index)}/2239|"
        f"measured={len(measurements)}|output={relative(output)}",
        flush=True,
    )


def run_lazy_single(index, args):
    """Production entry: raw MWState and a bounded lazy CVP frontier."""
    from run_mw_search import search_request
    campaign=load_campaign()
    if not 0<=index<len(campaign['rows']):raise ValueError('candidate index outside the declared input')
    row=campaign['rows'][index]
    curve,points,gram,_=Families().specialize(row)
    request={'schema':'elliptic-curves.lazy-mw-search.v1',
             'curve':list(map(str,curve.a_invariants())),
             'points':[[str(p[0]),str(p[1])] for p in points],
             'metric_gram':[[str(v) for v in r] for r in gram],
             'policy':{'enumeration_backend':'gmp-pointed-sieve','chart_metric_weight':args.metric_weight,'diversity_window':2},
             'next_holes':args.next_holes,'height':args.height,'seconds_per_chart':args.seconds,
             'cvp_node_budget':args.cvp_node_budget,'extra_primes':[211,223,227,229,233,239,241,251]}
    output=args.output or ROOT/f'artifacts/local/runtime-mw17/index-{index}-lazy.json'
    if output.exists():raise FileExistsError('retain the prior search result and choose a new output')
    result=search_request(request,output)
    print(f"MW17LAZY|index={index}|charts={len(result['charts'])}|gain={result['certified_rank_gain']}|output={output}",flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single-index", type=int)
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--chunk-count", type=int, default=4)
    parser.add_argument("--chunk-dir", type=Path, default=CHUNK_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-new", type=int)
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument('--legacy-census-regression',action='store_true',help='explicitly regenerate the frozen 43/301-centre campaign')
    parser.add_argument('--next-holes',type=int,default=12)
    parser.add_argument('--height',type=int,default=10000)
    parser.add_argument('--seconds',default='1',help='exact decimal seconds per chart')
    parser.add_argument('--metric-weight',default='16')
    parser.add_argument('--cvp-node-budget',type=int,default=100000)
    parser.add_argument('--wall-seconds',type=int,default=300)
    parser.add_argument('--rss-bytes',type=int,default=1073741824)
    parser.add_argument('--lazy-worker',action='store_true',help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.single_index is not None and not args.legacy_census_regression:
        if args.lazy_worker:
            run_lazy_single(args.single_index,args)
        else:
            completed=captured_run([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--lazy-worker'],
                limits=Limits(args.wall_seconds,args.rss_bytes),text=True,check=True)
        return
    if args.chunk_index is not None and not args.legacy_census_regression:
        parser.error('production searches use --single-index with lazy budgets; --legacy-census-regression explicitly regenerates historical chunks')
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP is required")
    if args.check:
        campaign = load_campaign()
        Families()
        load_modules()
        print(f"MW17JUMPV2|status=PASS_PREFLIGHT|candidates={campaign['candidate_count']}")
        return
    if args.single_index is not None:
        run_single(args.single_index)
        return
    if args.merge:
        merge_chunks(args.chunk_dir, args.chunk_count, args.output or LEDGER)
        return
    if args.chunk_index is None:
        raise SystemExit("choose --single-index, --chunk-index, --merge, or --check")
    output = args.output or args.chunk_dir / f"chunk-{args.chunk_index:02d}-of-{args.chunk_count:02d}.json"
    run_chunk(args.chunk_index, args.chunk_count, output, args.max_new)


if __name__ == "__main__":
    main()
