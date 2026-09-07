#!/usr/bin/env python3
"""Run the R17 relative 2-Selmer suite with open-source backends.

The primary backend is PARI/GP's ``ell2cover`` through Sage.  PARI returns a
basis of the everywhere locally soluble 2-covers as binary quartics together
with their maps to the elliptic curve.  The worker receives the curve only:
neither the specialized generic sections nor the held-out exceptional points
are supplied during descent or bounded quartic point search.

After a worker completes, the supervisor replays the exact input sources and
uses finite reductions modulo good primes to identify any rational points
found blindly on the returned covers.  If those recovered basis classes span
the known Mordell--Weil image, this also gives the point-to-Selmer embedding
and an explicit quotient basis.  A timeout, failed BNF certification, or an
incomplete recovered span is retained as missing evidence, never promoted to
a Selmer or Mordell--Weil conclusion.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import time
from typing import Any, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
INPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-suite-input.v1"
OUTPUT_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-open-run.v1"
CASE_SCHEMA = "elliptic-curves.elkies-2026-relative-2selmer-pari-case.v1"
PROTOCOL = "ELKIESR17OPENREL2"
GENERIC_RANK = 17
DEFAULT_SAGE_PYTHON = Path(
    "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python"
)


PARI_WORKER = r'''from __future__ import annotations
import json
from pathlib import Path
import time

from sage.all import EllipticCurve, QQ, pari
from sage.version import version as sage_version

payload = json.loads(Path(INPUT_PATH).read_text())
pari.allocatemem(int(payload["pari_stack_bytes"]))
factor_hints = [int(value) for value in payload.get("factor_hint_primes", [])]
if factor_hints:
    pari.addprimes(factor_hints)

def stage(name, status, **values):
    suffix = "".join(f"|{key}={value}" for key, value in values.items())
    print(f"ELKIESR17OPENREL2|case={payload['case_id']}|stage={name}|status={status}{suffix}", flush=True)

def rational_text(value):
    return str(value)

model = [QQ(value) for value in payload["global_minimal_model"]]
curve = EllipticCurve(QQ, model)
pari_curve = pari(curve)
two_torsion_dimension = int(curve.two_torsion_rank())
stage("ellrankinit", "start", factor_hints=len(factor_hints))
started = time.monotonic()
context = pari_curve.ellrankinit()
init_seconds = time.monotonic() - started
stage("ellrankinit", "complete", seconds=f"{init_seconds:.6f}")

stage("bnfcertify", "start")
started = time.monotonic()
certificates = []
try:
    field_data = context[2]
    for index in range(len(field_data)):
        certificates.append(int(pari.bnfcertify(field_data[index])))
except Exception as error:
    stage("bnfcertify", "error", error=type(error).__name__)
    raise
if not certificates or any(value != 1 for value in certificates):
    raise ArithmeticError(f"PARI did not certify every descent field: {certificates}")
certify_seconds = time.monotonic() - started
stage("bnfcertify", "complete", seconds=f"{certify_seconds:.6f}", fields=len(certificates))

stage("ell2cover", "start")
started = time.monotonic()
covers = context.ell2cover()
descent_seconds = time.monotonic() - started
stage("ell2cover", "complete", seconds=f"{descent_seconds:.6f}", dimension=len(covers))

cover_records = []
xy_variables = pari('[x,y]')
stage("blind_search", "start", covers=len(covers), bound=payload["search_bound"])
blind_started = time.monotonic()
for index in range(len(covers)):
    quartic = covers[index][0]
    cover_map = covers[index][1]
    search_started = time.monotonic()
    points = pari.hyperellratpoints(quartic, int(payload["search_bound"]))
    search_seconds = time.monotonic() - search_started
    record = {
        "basis_index": index + 1,
        "selmer_basis_bits": [1 if column == index else 0 for column in range(len(covers))],
        "quartic": str(quartic),
        "map_x": str(cover_map[0]),
        "map_y": str(cover_map[1]),
        "search_bound": int(payload["search_bound"]),
        "search_seconds": search_seconds,
        "signed_affine_cover_point_count": len(points),
        "search_status": "point_found" if len(points) else "no_point_within_bound",
    }
    if len(points):
        cover_point = points[0]
        image = pari.substvec(cover_map, xy_variables, cover_point)
        if int(pari.ellisoncurve(pari_curve, image)) != 1:
            raise ArithmeticError("a PARI cover map returned a point off the curve")
        record["first_cover_point"] = [rational_text(value) for value in cover_point]
        record["first_elliptic_image"] = [rational_text(image[0]), rational_text(image[1])]
    cover_records.append(record)
    stage(
        "blind_cover",
        "complete",
        index=index + 1,
        found=(1 if len(points) else 0),
        seconds=f"{search_seconds:.6f}",
    )
blind_seconds = time.monotonic() - blind_started
stage("blind_search", "complete", seconds=f"{blind_seconds:.6f}")

result = {
    "schema": "elliptic-curves.elkies-2026-relative-2selmer-pari-worker.v1",
    "case_id": payload["case_id"],
    "backend": "SageMath/PARI ellrankinit+bnfcertify+ell2cover",
    "sage_version": str(sage_version),
    "pari_version": str(pari.version()),
    "global_minimal_model": payload["global_minimal_model"],
    "two_torsion_dimension": two_torsion_dimension,
    "descent_field_bnf_certificates": certificates,
    "class_group_data_certified": True,
    "everywhere_locally_soluble_cover_basis_complete": True,
    "total_two_selmer_dimension": len(covers),
    "timings": {
        "ellrankinit_seconds": init_seconds,
        "bnfcertify_seconds": certify_seconds,
        "ell2cover_seconds": descent_seconds,
        "blind_cover_search_seconds": blind_seconds,
    },
    "factor_hint_primes": [str(value) for value in factor_hints],
    "blind_search": {
        "public_exceptional_points_supplied": False,
        "generic_points_supplied": False,
        "search_bound": int(payload["search_bound"]),
        "covers": cover_records,
    },
}
Path(OUTPUT_PATH).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
stage("worker", "complete", output=OUTPUT_PATH)
'''


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def query_open_source_versions(sage_python: str) -> dict[str, str]:
    program = (
        "import json; from sage.all import pari; "
        "from sage.version import version; "
        "print(json.dumps({'sage':str(version),'pari':str(pari.version())},sort_keys=True))"
    )
    completed = subprocess.run(
        [sage_python, "-c", program],
        check=True,
        capture_output=True,
        text=True,
        timeout=30,
    )
    return json.loads(completed.stdout.strip().splitlines()[-1])


def f2_rank(rows: Iterable[Sequence[int]]) -> int:
    basis: dict[int, int] = {}
    for row in rows:
        value = sum((int(bit) & 1) << index for index, bit in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def f2_linear_combination_coefficients(
    vectors: Sequence[Sequence[int]], target: Sequence[int]
) -> list[int] | None:
    """Return coefficients expressing ``target`` in the span of ``vectors``."""

    width = len(target)
    if any(len(vector) != width for vector in vectors):
        raise ValueError("inconsistent binary-vector widths")
    augmented = []
    for index, vector in enumerate(vectors):
        value = sum((int(bit) & 1) << column for column, bit in enumerate(vector))
        augmented.append((value, 1 << index))
    basis: dict[int, tuple[int, int]] = {}
    for value, combination in augmented:
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                old_value, old_combination = basis[pivot]
                value ^= old_value
                combination ^= old_combination
            else:
                basis[pivot] = (value, combination)
                break
    value = sum((int(bit) & 1) << column for column, bit in enumerate(target))
    combination = 0
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            return None
        old_value, old_combination = basis[pivot]
        value ^= old_value
        combination ^= old_combination
    return [(combination >> index) & 1 for index in range(len(vectors))]


def extend_rows_to_standard_basis(
    initial_rows: Sequence[Sequence[int]], width: int
) -> list[list[int]]:
    """Choose standard unit rows extending an independent row family to a basis."""

    rows = [list(map(int, row)) for row in initial_rows]
    if any(len(row) != width for row in rows):
        raise ValueError("an initial row has the wrong width")
    initial_rank = f2_rank(rows)
    if initial_rank != len(rows):
        raise ValueError("the initial rows are not independent")
    extension = []
    rank = initial_rank
    for index in range(width):
        unit = [1 if column == index else 0 for column in range(width)]
        candidate_rank = f2_rank([*rows, *extension, unit])
        if candidate_rank > rank:
            extension.append(unit)
            rank = candidate_rank
        if rank == width:
            break
    if rank != width:
        raise ArithmeticError("failed to extend rows to a binary basis")
    return extension


def flatten_certificate_rows(certificate: dict[str, Any]) -> list[list[int]]:
    return [
        [int(value) for value in row]
        for signature in certificate["signatures"]
        for row in signature["rows"]
    ]


def load_authoritative_cases() -> dict[str, Any]:
    import sys

    sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]
    from build_elkies_2026_relative_2selmer_suite import (
        load_high_rank_cases,
        load_nagao_cases,
        load_rank21_case,
        load_record_pair_cases,
    )

    cases = [
        load_rank21_case(),
        *load_high_rank_cases(),
        *load_record_pair_cases(),
        *load_nagao_cases(1000),
    ]
    return {case.case_id: case for case in cases}


def factor_hints(case_id: str) -> tuple[int, ...]:
    import sys

    sys.path.insert(0, str(CAS))
    if case_id == "control-r21-t3_8":
        from verify_icarm_curve394_rank21 import DISCRIMINANT_FACTORIZATION

        return tuple(prime for prime, _exponent in DISCRIMINANT_FACTORIZATION)
    if case_id == "control-r28":
        from build_elkies_2026_rank28_bad_place_ledger import (
            DISCRIMINANT_FACTORIZATION,
        )

        return tuple(prime for prime, _exponent in DISCRIMINANT_FACTORIZATION)
    if case_id in ("record-r29-356", "record-r29-385"):
        source = (
            ROOT
            / "artifacts/generated-results"
            / "elkies-k3-r17-norm12-icarm-public-fibres-v1.json"
        )
        curve_id = int(case_id.rsplit("-", 1)[1])
        document = json.loads(source.read_text())
        row = next(record for record in document["records"] if int(record["id"]) == curve_id)
        return tuple(int(prime) for prime in row["bad_primes"])
    return ()


def read_rss_bytes(pid: int) -> int:
    for line in Path(f"/proc/{pid}/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) * 1024
    return 0


def stop_process_group(process: subprocess.Popen[str], sig: signal.Signals) -> None:
    if process.poll() is None:
        os.killpg(process.pid, sig)


def supervise_source(
    sage_python: str,
    worker_source: str,
    payload: dict[str, Any],
    result_path: Path,
    log_path: Path,
    *,
    timeout: float,
    rss_limit_bytes: int,
) -> dict[str, Any]:
    result_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as input_handle:
        json.dump(payload, input_handle)
        input_path = Path(input_handle.name)
    worker_text = worker_source.replace("INPUT_PATH", repr(str(input_path))).replace(
        "OUTPUT_PATH", repr(str(result_path))
    )
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as worker_handle:
        worker_handle.write(worker_text)
        worker_path = Path(worker_handle.name)
    try:
        with log_path.open("w") as log:
            process = subprocess.Popen(
                [sage_python, str(worker_path)],
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                start_new_session=True,
            )
            started = time.monotonic()
            peak_rss = 0
            outcome = "running"
            try:
                while process.poll() is None:
                    elapsed = time.monotonic() - started
                    if elapsed >= timeout:
                        outcome = "strict_wall_timeout"
                        stop_process_group(process, signal.SIGTERM)
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            stop_process_group(process, signal.SIGKILL)
                            process.wait()
                        break
                    try:
                        peak_rss = max(peak_rss, read_rss_bytes(process.pid))
                    except (FileNotFoundError, ProcessLookupError):
                        pass
                    if peak_rss > rss_limit_bytes:
                        outcome = "strict_rss_limit"
                        stop_process_group(process, signal.SIGTERM)
                        try:
                            process.wait(timeout=5)
                        except subprocess.TimeoutExpired:
                            stop_process_group(process, signal.SIGKILL)
                            process.wait()
                        break
                    time.sleep(0.25)
            except BaseException:
                stop_process_group(process, signal.SIGTERM)
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    stop_process_group(process, signal.SIGKILL)
                    process.wait()
                raise
            wall_seconds = time.monotonic() - started
        if outcome == "running":
            outcome = (
                "completed"
                if process.returncode == 0 and result_path.exists()
                else "backend_failure"
            )
        return {
            "outcome": outcome,
            "returncode": process.returncode,
            "wall_seconds": wall_seconds,
            "peak_observed_rss_bytes": peak_rss,
            "timeout_seconds": timeout,
            "rss_limit_bytes": rss_limit_bytes,
            "log": str(log_path),
            "log_sha256": file_sha256(log_path),
            "worker_result": str(result_path) if result_path.exists() else None,
            "worker_result_sha256": (
                file_sha256(result_path) if result_path.exists() else None
            ),
        }
    finally:
        input_path.unlink(missing_ok=True)
        worker_path.unlink(missing_ok=True)


def supervise_worker(
    sage_python: str,
    payload: dict[str, Any],
    result_path: Path,
    log_path: Path,
    *,
    timeout: float,
    rss_limit_bytes: int,
) -> dict[str, Any]:
    return supervise_source(
        sage_python,
        PARI_WORKER,
        payload,
        result_path,
        log_path,
        timeout=timeout,
        rss_limit_bytes=rss_limit_bytes,
    )


def _fraction_point(values: Sequence[str]) -> tuple[Fraction, Fraction]:
    if len(values) != 2:
        raise ValueError("expected an affine elliptic point")
    return Fraction(values[0]), Fraction(values[1])


def identify_blind_recoveries(
    case: Any, worker: dict[str, Any], *, prime_bound: int
) -> dict[str, Any]:
    """Identify blindly found cover images in the known point basis when possible."""

    import sys

    sys.path[:0] = [str(ROOT / "elliptic-curves"), str(CAS)]
    from ecsearch.q12o5867_specialization import short_certificate_model
    from elliptic_candidate_record import (
        build_finite_quotient_certificate,
        source_point_to_target,
    )

    recovered_covers = [
        cover
        for cover in worker["blind_search"]["covers"]
        if "first_elliptic_image" in cover
    ]
    known_points = tuple(case.generic_points) + tuple(case.exceptional_points)
    recovered_points = tuple(
        _fraction_point(cover["first_elliptic_image"]) for cover in recovered_covers
    )
    short_model, change = short_certificate_model(case.model)
    all_points = known_points + recovered_points
    short_points = tuple(source_point_to_target(point, change) for point in all_points)
    certificate = build_finite_quotient_certificate(
        short_model,
        short_points,
        relation_prime=2,
        prime_bound=prime_bound,
    )
    rows = flatten_certificate_rows(certificate)
    known_count = len(known_points)
    known_vectors = [
        [row[column] for row in rows] for column in range(known_count)
    ]
    labels = []
    for offset, cover in enumerate(recovered_covers):
        target = [row[known_count + offset] for row in rows]
        coordinates = f2_linear_combination_coefficients(known_vectors, target)
        labels.append(
            {
                "basis_index": cover["basis_index"],
                "known_mw_coordinates": coordinates,
                "generic_coordinates": (
                    None if coordinates is None else coordinates[:GENERIC_RANK]
                ),
                "exceptional_quotient_coordinates": (
                    None if coordinates is None else coordinates[GENERIC_RANK:]
                ),
                "outside_known_mw_span_detected": coordinates is None,
            }
        )

    total_dimension = int(worker["total_two_selmer_dimension"])
    recovered_known = [row for row in labels if row["known_mw_coordinates"] is not None]
    known_coordinate_rows = [row["known_mw_coordinates"] for row in recovered_known]
    recovered_known_rank = f2_rank(known_coordinate_rows)
    selmer_embedding = None
    if recovered_known_rank == known_count:
        recovered_selmer_vectors = []
        for row in recovered_known:
            vector = [0] * total_dimension
            vector[int(row["basis_index"]) - 1] = 1
            recovered_selmer_vectors.append(vector)
        point_rows = []
        for point_index in range(known_count):
            target = [1 if index == point_index else 0 for index in range(known_count)]
            coefficients = f2_linear_combination_coefficients(
                known_coordinate_rows, target
            )
            assert coefficients is not None
            selmer_row = [0] * total_dimension
            for coefficient, vector in zip(coefficients, recovered_selmer_vectors):
                if coefficient:
                    selmer_row = [a ^ b for a, b in zip(selmer_row, vector)]
            point_rows.append(selmer_row)
        generic_selmer_rows = point_rows[:GENERIC_RANK]
        exceptional_selmer_rows = point_rows[GENERIC_RANK:]
        quotient_selmer_rows = extend_rows_to_standard_basis(
            generic_selmer_rows, total_dimension
        )
        aligned_basis = [*generic_selmer_rows, *quotient_selmer_rows]
        exceptional_quotient_rows = []
        for row in exceptional_selmer_rows:
            aligned = f2_linear_combination_coefficients(aligned_basis, row)
            assert aligned is not None
            exceptional_quotient_rows.append(aligned[GENERIC_RANK:])
        selmer_embedding = {
            "status": "complete_from_blindly_recovered_pari_basis_covers",
            "generic_selmer_rows": generic_selmer_rows,
            "exceptional_selmer_rows": exceptional_selmer_rows,
            "quotient_basis": [
                {
                    "quotient_basis_index": index + 1,
                    "selmer_basis_bits": row,
                    "pari_cover_basis_index": row.index(1) + 1,
                }
                for index, row in enumerate(quotient_selmer_rows)
            ],
            "exceptional_quotient_rows": exceptional_quotient_rows,
            "exceptional_quotient_rank": f2_rank(exceptional_quotient_rows),
        }

    quotient_rows = [
        row["exceptional_quotient_coordinates"]
        for row in labels
        if row["exceptional_quotient_coordinates"] is not None
    ]
    return {
        "finite_reduction_certificate": certificate,
        "known_mw_point_count": known_count,
        "known_mw_mod2_rank_in_certificate": f2_rank(known_vectors),
        "blind_basis_cover_point_count": len(recovered_covers),
        "blind_recovered_known_mw_rank": recovered_known_rank,
        "blind_recovered_exceptional_quotient_rank": f2_rank(quotient_rows),
        "cover_class_labels": labels,
        "point_to_selmer_embedding": selmer_embedding,
    }


def quotient_classification(case: Any, worker: dict[str, Any]) -> dict[str, Any]:
    total_dimension = int(worker["total_two_selmer_dimension"])
    residual_dimension = total_dimension - GENERIC_RANK
    exceptional_rank = len(case.exceptional_points)
    if residual_dimension < 0:
        raise ArithmeticError("the Selmer dimension is below the certified generic rank")
    if exceptional_rank > residual_dimension:
        raise ArithmeticError("the Selmer result contradicts the certified control subgroup")
    return {
        "total_two_selmer_dimension": total_dimension,
        "generic_mw17_image_dimension": GENERIC_RANK,
        "relative_quotient_dimension": residual_dimension,
        "known_exceptional_quotient_dimension": exceptional_rank,
        "known_realized_quotient_class_count_including_zero": 2**exceptional_rank,
        "unexplained_quotient_dimension": residual_dimension - exceptional_rank,
        "classes_not_realized_by_known_exceptional_subgroup": (
            2**residual_dimension - 2**exceptional_rank
        ),
        "rank_status": (
            "exact_rank_if_two_torsion_is_trivial_and_known_rank_matches_selmer_dimension"
            if case.certified_rank_lower_bound == total_dimension
            and int(worker["two_torsion_dimension"]) == 0
            else "selmer_upper_bound_and_existing_rank_lower_bound_do_not_match"
        ),
    }


def selected_manifest_cases(
    manifest: dict[str, Any], case_ids: set[str], controls_only: bool
) -> list[dict[str, Any]]:
    rows = []
    for row in manifest["cases"]:
        if case_ids and row["case_id"] not in case_ids:
            continue
        if controls_only and row["role"] != "held-out-positive-control":
            continue
        rows.append(row)
    if case_ids - {row["case_id"] for row in rows}:
        missing = sorted(case_ids - {row["case_id"] for row in rows})
        raise ValueError(f"case ids are absent from the manifest: {missing}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--controls-only", action="store_true")
    parser.add_argument("--sage-python", type=Path, default=DEFAULT_SAGE_PYTHON)
    parser.add_argument("--timeout-per-case", type=float, default=3600.0)
    parser.add_argument("--rss-limit-bytes", type=int, default=8_000_000_000)
    parser.add_argument("--pari-stack-bytes", type=int, default=4_000_000_000)
    parser.add_argument("--search-bound", type=int, default=1000)
    parser.add_argument("--certificate-prime-bound", type=int, default=1000)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if (
        args.timeout_per_case <= 0
        or args.rss_limit_bytes < 64_000_000
        or args.pari_stack_bytes < 64_000_000
        or args.search_bound < 0
        or args.certificate_prime_bound < 3
    ):
        parser.error("invalid time, memory, search, or certificate bound")
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema") != INPUT_SCHEMA:
        raise SystemExit("unexpected relative 2-Selmer input manifest")
    sage_python = shutil.which(str(args.sage_python))
    if sage_python is None:
        raise SystemExit(f"Sage Python is unavailable: {args.sage_python}")
    versions = query_open_source_versions(sage_python)
    authoritative = load_authoritative_cases()
    selected = selected_manifest_cases(manifest, set(args.case), args.controls_only)
    runs = []
    for manifest_case in selected:
        case_id = manifest_case["case_id"]
        case = authoritative.get(case_id)
        if case is None:
            raise SystemExit(f"authoritative input source lacks {case_id}")
        if [str(value) for value in case.model] != manifest_case["global_minimal_model"]:
            raise SystemExit(f"manifest model changed for {case_id}")
        result_path = args.output_dir / "cases" / f"{case_id}.json"
        log_path = args.output_dir / "logs" / f"{case_id}.log"
        if (result_path.exists() or log_path.exists()) and not args.overwrite:
            raise FileExistsError(result_path if result_path.exists() else log_path)
        result_path.unlink(missing_ok=True)
        payload = {
            "case_id": case_id,
            "global_minimal_model": manifest_case["global_minimal_model"],
            "search_bound": args.search_bound,
            "pari_stack_bytes": args.pari_stack_bytes,
            "factor_hint_primes": list(factor_hints(case_id)),
        }
        print(f"{PROTOCOL}|case={case_id}|stage=supervisor|status=start", flush=True)
        supervised = supervise_worker(
            sage_python,
            payload,
            result_path,
            log_path,
            timeout=args.timeout_per_case,
            rss_limit_bytes=args.rss_limit_bytes,
        )
        record: dict[str, Any] = {
            "case_id": case_id,
            "role": case.role,
            "parameter": case.parameter,
            "supervisor": supervised,
            "factor_hint_prime_count": len(payload["factor_hint_primes"]),
        }
        if supervised["outcome"] == "completed":
            worker = json.loads(result_path.read_text())
            if (
                worker.get("schema")
                != "elliptic-curves.elkies-2026-relative-2selmer-pari-worker.v1"
                or worker.get("case_id") != case_id
                or worker.get("global_minimal_model")
                != manifest_case["global_minimal_model"]
            ):
                raise ArithmeticError(f"worker result identity mismatch for {case_id}")
            record["worker"] = worker
            record["classification"] = quotient_classification(case, worker)
            record["blind_recovery"] = identify_blind_recoveries(
                case, worker, prime_bound=args.certificate_prime_bound
            )
            record["status"] = "COMPLETE_CERTIFIED_PARI_TWO_SELMER_BASIS"
        else:
            record["worker"] = None
            record["classification"] = None
            record["blind_recovery"] = None
            record["status"] = "INCOMPLETE_OPEN_SOURCE_DESCENT"
        runs.append(record)
        print(
            f"{PROTOCOL}|case={case_id}|stage=supervisor|status={record['status']}"
            f"|outcome={supervised['outcome']}|seconds={supervised['wall_seconds']:.6f}",
            flush=True,
        )
    complete = sum(run["status"].startswith("COMPLETE") for run in runs)
    output = {
        "schema": OUTPUT_SCHEMA,
        "status": (
            "COMPLETE_ALL_SELECTED_OPEN_SOURCE_DESCENTS"
            if complete == len(runs)
            else "INCOMPLETE_ONE_OR_MORE_OPEN_SOURCE_DESCENTS"
        ),
        "backend": {
            "name": "SageMath/PARI ellrankinit+bnfcertify+ell2cover",
            "license": "open_source",
            "full_selmer_basis_interface": "PARI ell2cover binary quartics",
            "bnf_certification_required": True,
            "sage_python": sage_python,
            "sage_version": versions["sage"],
            "pari_version": versions["pari"],
        },
        "run_parameters": {
            "selected_case_ids": [run["case_id"] for run in runs],
            "controls_only": args.controls_only,
            "timeout_per_case_seconds": args.timeout_per_case,
            "rss_limit_bytes": args.rss_limit_bytes,
            "pari_stack_bytes": args.pari_stack_bytes,
            "blind_search_bound": args.search_bound,
            "finite_reduction_certificate_prime_bound": args.certificate_prime_bound,
        },
        "input_manifest": {
            "path": str(args.manifest),
            "sha256": file_sha256(args.manifest),
        },
        "runner": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": file_sha256(Path(__file__).resolve()),
        },
        "selected_case_count": len(runs),
        "completed_case_count": complete,
        "runs": runs,
        "claim_boundary": [
            "Only a completed worker with successful bnfcertify and ell2cover is recorded as a full 2-Selmer basis.",
            "The blind worker receives neither generic nor public exceptional points.",
            "Finite-reduction labels are exact for the displayed point images but do not turn a bounded quartic search miss into non-solubility.",
            "PARI exposes basis quartics, not an addition operation on arbitrary cover classes; only returned basis classes have explicit quartics here.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(args.output)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|selected={len(runs)}|completed={complete}"
        f"|status={output['status']}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
