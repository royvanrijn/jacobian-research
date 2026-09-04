#!/usr/bin/env sage-python
"""Compute the sealed pre-search class quotient for the frozen R17 cohort.

For each fibre this program certifies the maximal cubic order and the PARI BNF,
computes every S-prime and specialized generic-MW17 half-ideal class, and then
forms

    Q_t = Cl(K_t) / (2 Cl(K_t) + <S_t, c(MW17_t)>).

Only ``bnfcertify == 1`` rows receive a value of ``dim_Q``.  Chunk supervisors
turn timeouts and backend failures into explicit censored feature rows.  The
merge artifact is eligible to unlock point search only if all 100 rows are
exactly certified.  This executable never invokes a Mordell--Weil point search.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import subprocess
import sys
import time
from typing import Any

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, QQ, ZZ, pari
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
COHORT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-small-field-class-quotient-cohort-v1.json"
)
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
DEFAULT_CHUNK_DIR = ROOT / "artifacts/local/elkies-k3/r17-small-field-class-quotient-features-v1"
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-small-field-class-quotient-features-v1.json"
)

COHORT_STATUS = "FROZEN_RANK_BLIND_PRE_CLASS_GROUP_COHORT"
CHUNK_SCHEMA = "elkies-k3.r17-small-field-class-quotient-feature-chunk.v1"
SCHEMA = "elkies-k3.r17-small-field-class-quotient-features.v1"
GENERIC_RANK = 17
DEFAULT_ROW_TIMEOUT_SECONDS = 3600
PARI_STACK_BYTES = 4_000_000_000

sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def f2_rank(rows: list[list[int]], width: int) -> int:
    if width == 0 or not rows:
        return 0
    return int(Matrix(GF(2), len(rows), width, [bit for row in rows for bit in row]).rank())


def qtext(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def homogeneous_value(coefficients, numerator: int, denominator: int):
    degree = len(coefficients) - 1
    return sum(
        QQ(coefficient) * numerator**index * denominator ** (degree - index)
        for index, coefficient in enumerate(coefficients)
    )


def load_inputs():
    cohort = json.loads(COHORT.read_text())
    model = json.loads(MODEL.read_text())
    sections = json.loads(SECTIONS.read_text())
    if cohort.get("status") != COHORT_STATUS:
        raise ArithmeticError("the small-field cohort is not frozen")
    if len(cohort.get("rows", [])) != 100:
        raise ArithmeticError("the frozen small-field cohort stopped having 100 rows")
    if any(
        row.get("feature_status") != "NOT_OPENED"
        or row.get("outcome_status") != "SEALED_UNTIL_ALL_FEATURES_FREEZE"
        for row in cohort["rows"]
    ):
        raise ArithmeticError("the cohort was opened out of phase")
    committed = [
        {
            key: row[key]
            for key in (
                "sample_id",
                "family",
                "parameter",
                "projective_pair",
                "projective_height",
                "cubic_order_discriminant",
                "absolute_cubic_order_discriminant_bits",
                "j_invariant",
            )
        }
        for row in cohort["rows"]
    ]
    if canonical_hash(committed) != cohort["commitment"]["candidate_list_sha256"]:
        raise ArithmeticError("the frozen candidate-list hash does not replay")
    if model.get("status") != "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
        raise ArithmeticError("the exact R17 model input is not passing")
    if sections.get("status") != "PASS_TRANSCRIBED_PUBLISHED_R17_SECTIONS_AND_CHORDS":
        raise ArithmeticError("the exact R17 section input is not passing")
    return cohort, model, sections


def specialized_curve_and_points(row, model, section_input):
    numerator, denominator = map(int, row["projective_pair"])
    coefficient_a = homogeneous_value(
        model["A_coefficients_low_to_high"], numerator, denominator
    )
    coefficient_b = homogeneous_value(
        model["B_coefficients_low_to_high"], numerator, denominator
    )
    if coefficient_a not in ZZ or coefficient_b not in ZZ:
        raise ArithmeticError("the published homogeneous specialization is not integral")
    coefficient_a, coefficient_b = ZZ(coefficient_a), ZZ(coefficient_b)
    curve = EllipticCurve(QQ, [coefficient_a, coefficient_b])

    records = section_input["sections"]
    x_coordinates = [
        QQ(homogeneous_value(record["x_coefficients_low_to_high"], numerator, denominator))
        for record in records
    ]
    y_coordinates: list[QQ | None] = [None] * len(records)
    for index, record in enumerate(records):
        if "y_coefficients_low_to_high" in record:
            y_coordinates[index] = QQ(
                homogeneous_value(
                    record["y_coefficients_low_to_high"], numerator, denominator
                )
            )
            continue
        chord = record["chord"]
        reference = int(chord["reference_basis_index"])
        if y_coordinates[reference] is None:
            raise ArithmeticError("a chord references an unreconstructed section")
        slope = QQ(
            homogeneous_value(
                chord["slope_coefficients_low_to_high"], numerator, denominator
            )
        )
        y_coordinates[index] = (
            y_coordinates[reference]
            + slope * (x_coordinates[index] - x_coordinates[reference])
        )
    points = tuple(
        curve(x_coordinate, y_coordinate)
        for x_coordinate, y_coordinate in zip(x_coordinates, y_coordinates)
    )
    if len(points) != GENERIC_RANK or any(point not in curve for point in points):
        raise ArithmeticError("a generic R17 section missed the specialization")

    minimal = curve.global_minimal_model()
    isomorphism = curve.isomorphism_to(minimal)
    minimal_points = tuple(isomorphism(point) for point in points)
    if any(point not in minimal for point in minimal_points):
        raise ArithmeticError("minimal-model point transport failed")
    return curve, minimal, minimal_points, isomorphism


def finite_reduction_record(curve, points):
    model = tuple(Fraction(str(value)) for value in curve.a_invariants())
    affine = tuple(
        (Fraction(str(point[0])), Fraction(str(point[1]))) for point in points
    )
    signatures = find_mod2_reduction_certificate(model, affine, prime_bound=1000)
    rank = combined_mod2_rank(signatures, len(points))
    return {
        "combined_mod2_rank": rank,
        "column_count": len(points),
        "proves_generic_subgroup_primitive_and_independent": rank == len(points),
        "prime_bound": 1000,
        "signatures": [
            {
                "prime": record.prime,
                "group_order": record.group_order,
                "doubled_subgroup_order": record.doubled_subgroup_order,
                "quotient_dimension": record.quotient_dimension,
                "rows": [list(row) for row in record.rows],
            }
            for record in signatures
        ],
    }


def divide_supported_rational(value, primes):
    value = QQ(str(value))
    numerator = abs(ZZ(value.numerator()))
    denominator = ZZ(value.denominator())
    support = []
    for prime in primes:
        exponent = 0
        while numerator % prime == 0:
            numerator //= prime
            exponent += 1
        while denominator % prime == 0:
            denominator //= prime
            exponent -= 1
        if exponent:
            support.append({"rational_prime": int(prime), "exponent": exponent})
    if numerator != 1 or denominator != 1:
        raise ArithmeticError("a Kummer half-ideal correction has support outside S")
    return support


def class_coordinates(bnf, ideal, cyclic_invariants):
    coordinates = [int(value) for value in pari.bnfisprincipal(bnf, ideal, 0)]
    if len(coordinates) != len(cyclic_invariants):
        raise ArithmeticError("PARI returned the wrong class-coordinate width")
    mod2 = [
        coordinates[index] & 1
        for index, invariant in enumerate(cyclic_invariants)
        if invariant % 2 == 0
    ]
    return coordinates, mod2


def half_ideal_record(nf, bnf, theta, curve, point, bad_primes, cyclic_invariants):
    a1, _a2, a3, _a4, _a6 = map(QQ, curve.a_invariants())
    x_coordinate, y_coordinate = map(QQ, (point[0], point[1]))
    scaled_x = 4 * x_coordinate
    denominator_root = isqrt(int(scaled_x.denominator()))
    if denominator_root**2 != scaled_x.denominator():
        raise ArithmeticError("the denominator of 4*x is not a square")
    integral_alpha = (
        pari(str(denominator_root**2 * scaled_x))
        - denominator_root**2 * theta
    )
    norm_root = QQ(
        denominator_root**3
        * 4
        * (2 * y_coordinate + a1 * x_coordinate + a3)
    )
    if norm_root.denominator() != 1:
        raise ArithmeticError("the scaled Kummer norm root is not integral")
    if QQ(str(pari.nfeltnorm(nf, integral_alpha))) != norm_root**2:
        raise ArithmeticError("the Kummer norm-square identity failed")
    half_ideal = pari.idealadd(nf, integral_alpha, pari(str(norm_root)))
    correction = pari.idealdiv(
        nf, pari.idealpow(nf, half_ideal, 2), integral_alpha
    )
    correction_norm = pari.idealnorm(nf, correction)
    coordinates, mod2 = class_coordinates(
        bnf, half_ideal, cyclic_invariants
    )
    return {
        "denominator_root_for_4x": str(denominator_root),
        "half_ideal_hnf": str(half_ideal),
        "half_ideal_norm": str(pari.idealnorm(nf, half_ideal)),
        "localized_square_correction_ideal_hnf": str(correction),
        "localized_square_correction_norm": str(correction_norm),
        "localized_square_correction_support": divide_supported_rational(
            correction_norm, bad_primes
        ),
        "class_group_coordinates": coordinates,
        "class_group_mod_2_coordinates": mod2,
    }


def compute_feature(index: int):
    started = time.monotonic()
    cohort, model, sections = load_inputs()
    if not 0 <= index < len(cohort["rows"]):
        raise ValueError("feature index lies outside the frozen cohort")
    row = cohort["rows"][index]
    source_curve, curve, points, isomorphism = specialized_curve_and_points(
        row, model, sections
    )
    finite_reduction = finite_reduction_record(curve, points)

    discriminant = ZZ(curve.discriminant())
    factorization = list(abs(discriminant).factor(proof=True))
    bad_primes = sorted({ZZ(2), *(prime for prime, _exponent in factorization)})
    if any(not prime.is_prime(proof=True) for prime in bad_primes):
        raise ArithmeticError("the bad-prime factorization is not proved prime")

    a1, a2, a3, a4, a6 = map(ZZ, curve.a_invariants())
    b2 = a1**2 + 4 * a2
    b4 = a1 * a3 + 2 * a4
    b6 = a3**2 + 4 * a6
    polynomial_ring = PolynomialRing(ZZ, "z")
    z = polynomial_ring.gen()
    polynomial = z**3 + b2 * z**2 + 8 * b4 * z + 16 * b6
    if not polynomial.is_irreducible():
        raise ArithmeticError("the minimal completed-square cubic is reducible")
    if abs(ZZ(polynomial.discriminant())) != 256 * abs(discriminant):
        raise ArithmeticError("the completed-square discriminant identity failed")

    # BNF on the raw completed-square polynomial can overflow purely because
    # its power basis has a large index.  Move the distinguished root into a
    # canonical polredabs presentation first; all Kummer elements below use
    # this exact embedding, so the field change is algebraic, not heuristic.
    reduced_polynomial = pari.polredabs(pari(polynomial))
    theta = pari.nfisincl(pari(polynomial), reduced_polynomial, 1)
    if theta == 0:
        raise ArithmeticError("polredabs did not return an isomorphic cubic field")
    pari.addprimes(bad_primes)
    nf = pari.nfinit([reduced_polynomial, bad_primes])
    if list(pari.nfcertify(nf)):
        raise ArithmeticError("the maximal cubic order failed nfcertify")
    pari.allocatemem(PARI_STACK_BYTES)
    bnf_started = time.monotonic()
    bnf = pari.bnfinit(nf, 1)
    bnf_seconds = time.monotonic() - bnf_started
    certification_started = time.monotonic()
    bnf_certified = int(pari.bnfcertify(bnf)) == 1
    certification_seconds = time.monotonic() - certification_started
    if not bnf_certified:
        raise ArithmeticError("PARI did not unconditionally certify the BNF")

    cyclic_invariants = [int(value) for value in bnf.bnf_get_cyc()]
    class_number = 1
    for invariant in cyclic_invariants:
        class_number *= invariant
    even_indices = [
        index for index, invariant in enumerate(cyclic_invariants) if invariant % 2 == 0
    ]
    class_group_2_quotient_dimension = len(even_indices)
    signature = [int(value) for value in nf.nf_get_sign()]

    s_rows = []
    s_prime_records = []
    for rational_prime in bad_primes:
        for prime_index, ideal in enumerate(
            pari.idealprimedec(nf, rational_prime), start=1
        ):
            coordinates, mod2 = class_coordinates(bnf, ideal, cyclic_invariants)
            s_rows.append(mod2)
            s_prime_records.append(
                {
                    "rational_prime": int(rational_prime),
                    "prime_index_one_based": prime_index,
                    "ramification_index": int(ideal[2]),
                    "residue_degree": int(ideal[3]),
                    "ideal_hnf": str(pari.idealhnf(nf, ideal)),
                    "class_group_coordinates": coordinates,
                    "class_group_mod_2_coordinates": mod2,
                }
            )

    generic_records = []
    generic_rows = []
    for point_index, point in enumerate(points, start=1):
        record = half_ideal_record(
            nf,
            bnf,
            theta,
            curve,
            point,
            bad_primes,
            cyclic_invariants,
        )
        generic_rows.append(record["class_group_mod_2_coordinates"])
        generic_records.append({"label": f"G{point_index}", **record})

    killed_rows = s_rows + generic_rows
    killed_rank = f2_rank(killed_rows, class_group_2_quotient_dimension)
    dimension_q = class_group_2_quotient_dimension - killed_rank
    if dimension_q < 0:
        raise ArithmeticError("the computed quotient dimension is negative")

    fundamental_units = list(bnf.bnf_get_fu())
    unit_records = [
        {
            "unit": str(unit),
            "norm": qtext(pari.nfeltnorm(nf, unit)),
        }
        for unit in fundamental_units
    ]
    local_data = []
    for prime in bad_primes:
        data = curve.local_data(prime)
        local_data.append(
            {
                "prime": int(prime),
                "kodaira_symbol": str(data.kodaira_symbol()),
                "conductor_valuation": int(data.conductor_valuation()),
                "tamagawa_number": int(data.tamagawa_number()),
            }
        )

    feature = {
        "sample_id": row["sample_id"],
        "manifest_index": index,
        "family": row["family"],
        "parameter": row["parameter"],
        "status": "PASS_COMPLETE_UNCONDITIONAL_CLASS_QUOTIENT",
        "dim_Q": dimension_q,
        "definition": (
            "Cl(K)/(2Cl(K)+<prime ideals above 2*Delta_min, "
            "localized half-ideal classes of specialized generic MW17>)"
        ),
        "curve": {
            "homogeneous_source_model": [str(value) for value in source_curve.a_invariants()],
            "global_minimal_model": [str(value) for value in curve.a_invariants()],
            "global_minimal_discriminant": str(discriminant),
            "global_minimal_discriminant_factorization": [
                [str(prime), int(exponent)] for prime, exponent in factorization
            ],
            "root_number": int(curve.root_number()),
            "source_to_minimal_isomorphism_u_r_s_t": [
                qtext(value) for value in isomorphism.tuple()
            ],
        },
        "generic_mw17": {
            "point_count": len(points),
            "finite_reduction_certificate": finite_reduction,
            "localized_half_ideal_classes": generic_records,
            "image_rank_in_class_group_mod_2": f2_rank(
                generic_rows, class_group_2_quotient_dimension
            ),
        },
        "cubic_field": {
            "completed_square_polynomial": str(polynomial),
            "polredabs_polynomial": str(reduced_polynomial),
            "completed_square_root_in_polredabs_field": str(theta),
            "polynomial_discriminant": str(polynomial.discriminant()),
            "field_discriminant": str(nf.nf_get_disc()),
            "absolute_field_discriminant_bits": abs(int(nf.nf_get_disc())).bit_length(),
            "signature": signature,
            "nfcertify_passed": True,
        },
        "class_group": {
            "cyclic_invariants": cyclic_invariants,
            "class_number": str(class_number),
            "class_group_mod_2_coordinate_indices": even_indices,
            "dim_Cl_mod_2Cl": class_group_2_quotient_dimension,
            "bnfcertify_passed": True,
        },
        "S": {
            "definition": "all prime ideals above rational primes dividing 2*Delta_min",
            "bad_rational_primes": [int(prime) for prime in bad_primes],
            "prime_ideal_classes": s_prime_records,
            "image_rank_in_class_group_mod_2": f2_rank(
                s_rows, class_group_2_quotient_dimension
            ),
        },
        "killed_class_span_rank_mod_2": killed_rank,
        "units": {
            "signature": signature,
            "free_unit_rank": sum(signature) - 1,
            "norm_positive_unit_squareclass_dimension_bound": sum(signature) - 1,
            "fundamental_units": unit_records,
        },
        "local_data": local_data,
        "timing": {
            "bnfinit_seconds": bnf_seconds,
            "bnfcertify_seconds": certification_seconds,
            "total_seconds": time.monotonic() - started,
        },
        "proof_boundary": (
            "The class group, S-image, and generic Kummer image are exact.  This "
            "does not compute a Selmer group or any exceptional Mordell-Weil point."
        ),
    }
    feature["feature_sha256"] = canonical_hash(feature)
    return feature


def result_line(feature):
    print("RESULT_JSON=" + canonical_text(feature), flush=True)


def failure_record(row, index: int, status: str, detail: Any, elapsed: float):
    return {
        "sample_id": row["sample_id"],
        "manifest_index": index,
        "family": row["family"],
        "parameter": row["parameter"],
        "status": status,
        "dim_Q": None,
        "failure": detail,
        "elapsed_seconds": elapsed,
    }


def parse_worker(completed, row, index: int, elapsed: float):
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    result = next(
        (line for line in reversed(lines) if line.startswith("RESULT_JSON=")), None
    )
    if completed.returncode == 0 and result is not None:
        return json.loads(result[len("RESULT_JSON=") :])
    return failure_record(
        row,
        index,
        "CENSORED_FEATURE_BACKEND_FAILURE",
        {"returncode": completed.returncode, "output_tail": lines[-40:]},
        elapsed,
    )


def write_chunk(path, chunk_index, chunk_count, indices, records, candidate_hash):
    document = {
        "schema": CHUNK_SCHEMA,
        "status": (
            "COMPLETE_SCHEDULED_FEATURE_CHUNK"
            if len(records) == len(indices)
            else "PARTIAL_FEATURE_CHECKPOINT"
        ),
        "candidate_list_sha256": candidate_hash,
        "chunk_index": chunk_index,
        "chunk_count": chunk_count,
        "scheduled_indices": indices,
        "completed_record_count": len(records),
        "records": records,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")


def run_chunk(chunk_index, chunk_count, output, row_timeout_seconds, limit):
    cohort, _model, _sections = load_inputs()
    indices = [
        index for index in range(len(cohort["rows"]))
        if index % chunk_count == chunk_index
    ]
    if limit is not None:
        indices = indices[:limit]
    records = []
    if output.exists():
        old = json.loads(output.read_text())
        if (
            old.get("candidate_list_sha256")
            != cohort["commitment"]["candidate_list_sha256"]
            or old.get("chunk_index") != chunk_index
            or old.get("chunk_count") != chunk_count
            or old.get("scheduled_indices") != indices
        ):
            raise ArithmeticError("the existing feature checkpoint belongs to another schedule")
        records = old["records"]
    completed_ids = {record["sample_id"] for record in records}
    for position, index in enumerate(indices, start=1):
        row = cohort["rows"][index]
        if row["sample_id"] in completed_ids:
            continue
        command = [sys.executable, str(Path(__file__).resolve()), "--single-index", str(index)]
        started = time.monotonic()
        try:
            completed = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=row_timeout_seconds,
                check=False,
            )
            record = parse_worker(completed, row, index, time.monotonic() - started)
        except subprocess.TimeoutExpired as error:
            output_tail = error.stdout or ""
            if isinstance(output_tail, bytes):
                output_tail = output_tail.decode(errors="replace")
            record = failure_record(
                row,
                index,
                "CENSORED_FEATURE_TIMEOUT",
                {"output_tail": output_tail.splitlines()[-40:]},
                time.monotonic() - started,
            )
        records.append(record)
        write_chunk(
            output,
            chunk_index,
            chunk_count,
            indices,
            records,
            cohort["commitment"]["candidate_list_sha256"],
        )
        print(
            f"R17SMALLFIELDFEATURE|chunk={chunk_index}/{chunk_count}"
            f"|row={position}/{len(indices)}|sample={row['sample_id']}"
            f"|status={record['status']}|dimQ={record.get('dim_Q')}",
            flush=True,
        )


def merge_chunks(chunk_dir, chunk_count, output):
    cohort, _model, _sections = load_inputs()
    candidate_hash = cohort["commitment"]["candidate_list_sha256"]
    records_by_id = {}
    chunks = []
    for chunk_index in range(chunk_count):
        path = chunk_dir / f"chunk-{chunk_index:02d}-of-{chunk_count:02d}.json"
        chunk = json.loads(path.read_text())
        if chunk.get("status") != "COMPLETE_SCHEDULED_FEATURE_CHUNK":
            raise ArithmeticError(f"feature chunk {chunk_index} is incomplete")
        if chunk.get("candidate_list_sha256") != candidate_hash:
            raise ArithmeticError(f"feature chunk {chunk_index} names another cohort")
        for record in chunk["records"]:
            if record["sample_id"] in records_by_id:
                raise ArithmeticError("duplicate feature row across chunks")
            records_by_id[record["sample_id"]] = record
        chunks.append(
            {
                "path": relative(path),
                "sha256": digest(path),
                "record_count": len(chunk["records"]),
            }
        )
    expected_ids = [row["sample_id"] for row in cohort["rows"]]
    if set(records_by_id) != set(expected_ids):
        raise ArithmeticError("feature chunks do not cover the frozen cohort")
    records = [records_by_id[sample_id] for sample_id in expected_ids]
    exact = all(
        record.get("status") == "PASS_COMPLETE_UNCONDITIONAL_CLASS_QUOTIENT"
        and record.get("dim_Q") is not None
        and record.get("class_group", {}).get("bnfcertify_passed") is True
        for record in records
    )
    feature_commitment = canonical_hash(records) if exact else None
    document = {
        "schema": SCHEMA,
        "status": (
            "FROZEN_COMPLETE_UNCONDITIONAL_PRE_SEARCH_FEATURES"
            if exact
            else "INCOMPLETE_FEATURES_POINT_SEARCH_REMAINS_SEALED"
        ),
        "candidate_list_sha256": candidate_hash,
        "feature_commitment_sha256": feature_commitment,
        "point_search_unlocked": exact,
        "summary": {
            "scheduled_rows": len(records),
            "status_counts": dict(sorted(Counter(record["status"] for record in records).items())),
            "dim_Q_histogram": (
                dict(
                    sorted(
                        Counter(str(record["dim_Q"]) for record in records).items(),
                        key=lambda item: int(item[0]),
                    )
                )
                if exact
                else None
            ),
        },
        "records": records,
        "chunk_provenance": chunks,
        "inputs": {
            relative(path): digest(path) for path in (COHORT, MODEL, SECTIONS)
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
            "pari_stack_bytes_per_worker": PARI_STACK_BYTES,
            "recommended_maximum_concurrent_workers_on_32GiB_host": 4,
            "commands": [
                f"sage -python {relative(Path(__file__))} --chunk-index I --chunk-count {chunk_count}",
                f"sage -python {relative(Path(__file__))} --merge --chunk-count {chunk_count}",
            ],
        },
        "phase_boundary": (
            "No detector executable may run unless this status is the complete frozen "
            "status and a separately generated detector protocol pins this whole file."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(
        f"R17SMALLFIELDFEATUREMERGE|rows={len(records)}|complete={exact}"
        f"|hash={feature_commitment}|output={relative(output)}",
        flush=True,
    )


def audit_feature_artifact(path):
    cohort, _model, _sections = load_inputs()
    document = json.loads(path.read_text())
    records = document.get("records", [])
    if document.get("candidate_list_sha256") != cohort["commitment"]["candidate_list_sha256"]:
        raise ArithmeticError("the feature artifact names another cohort")
    complete = document.get("status") == "FROZEN_COMPLETE_UNCONDITIONAL_PRE_SEARCH_FEATURES"
    if complete:
        if len(records) != len(cohort["rows"]):
            raise ArithmeticError("the complete feature artifact has the wrong row count")
        if canonical_hash(records) != document.get("feature_commitment_sha256"):
            raise ArithmeticError("the frozen feature commitment does not replay")
        if not all(
            row.get("class_group", {}).get("bnfcertify_passed") is True
            and row.get("dim_Q") is not None
            for row in records
        ):
            raise ArithmeticError("an uncertified row entered the complete feature artifact")
    print(
        f"R17SMALLFIELDFEATUREAUDIT|status={document.get('status')}"
        f"|rows={len(records)}|complete={complete}",
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single-index", type=int)
    parser.add_argument("--chunk-index", type=int)
    parser.add_argument("--chunk-count", type=int, default=16)
    parser.add_argument("--chunk-dir", type=Path, default=DEFAULT_CHUNK_DIR)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--row-timeout-seconds", type=int, default=DEFAULT_ROW_TIMEOUT_SECONDS)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--audit", type=Path)
    args = parser.parse_args()
    modes = sum(
        value is not None for value in (args.single_index, args.chunk_index, args.audit)
    ) + int(args.merge)
    if modes != 1:
        raise SystemExit("choose exactly one of --single-index, --chunk-index, --merge, or --audit")
    if args.single_index is not None:
        result_line(compute_feature(args.single_index))
    elif args.chunk_index is not None:
        if not 0 <= args.chunk_index < args.chunk_count:
            raise SystemExit("chunk index is outside chunk count")
        output = args.output or (
            args.chunk_dir
            / f"chunk-{args.chunk_index:02d}-of-{args.chunk_count:02d}.json"
        )
        run_chunk(
            args.chunk_index,
            args.chunk_count,
            output.resolve(),
            args.row_timeout_seconds,
            args.limit,
        )
    elif args.merge:
        merge_chunks(
            args.chunk_dir.resolve(),
            args.chunk_count,
            (args.output or OUTPUT).resolve(),
        )
    else:
        audit_feature_artifact(args.audit.resolve())


if __name__ == "__main__":
    main()
