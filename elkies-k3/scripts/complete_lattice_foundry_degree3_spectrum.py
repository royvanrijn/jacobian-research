#!/usr/bin/env python3
"""Compute complete degree-three coset-minimum spectra for rootless targets.

For a rootless frame ``NS=U+M(-1)``, section translation identifies a
degree-three class with a residue in ``M/3M``.  Its all-section intersection
minimum is determined by the exact minimum norm of that residue class.  This
script visits every one of the ``3^17`` residues, using inversion to reduce the
CVP workload to ``(3^17+1)/2`` representatives, and emits the complete minimum
norm histogram.

The default batch is the five leading route-aware source classes from the
pre-prescribed-root ranking that also occurred in the 256-coset pilot, plus
published R17 as a control.  Each returned CVP vector has its norm recomputed
with the integral Gram matrix.  A deterministic subset is independently
repeated with MPFR GSO arithmetic.  Work is split into checkpointed chunks so
an interrupted multi-hour census can resume without repeating completed
chunks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing
import os
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATABASE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-top5-v1.json"
)
DEFAULT_FRAMES = [
    "NS0024-F005",
    "NS0005-F008",
    "NS0022-F011",
    "NS0033-F026",
    "NS0002-F007",
    "NS0001-F001",
]
PRE_PRESCRIBED_ROUTE_AWARE_TOP_FIVE = DEFAULT_FRAMES[:5]
CURRENT_MW2_SOURCE_RANKED_TOP_FIVE = [
    "NS0028-F005",
    "NS0011-F002",
    "NS0022-F011",
    "NS0005-F008",
    "NS0001-F001",
]
PUBLISHED_R17_CONTROL = DEFAULT_FRAMES[5]
DIMENSION = 17
DEGREE = 3
TOTAL_COSETS = DEGREE**DIMENSION
TOTAL_INVERSION_REPRESENTATIVES = (TOTAL_COSETS + 1) // 2


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gram_digest(gram: list[list[int]]) -> str:
    text = "\n".join(" ".join(map(str, row)) for row in gram) + "\n"
    return hashlib.sha256(text.encode()).hexdigest()


def exact_norm(vector: list[int], gram: list[list[int]]) -> int:
    total = 0
    for i, value in enumerate(vector):
        total += gram[i][i] * value * value
        total += 2 * value * sum(
            gram[i][j] * vector[j] for j in range(i + 1, DIMENSION)
        )
    return total


def target_coordinates(residue: list[int], mu: list[list[float]]) -> list[float]:
    # If y=-residue/3 in the original basis, its Gram--Schmidt coordinates are
    # y_i + sum_{j>i} y_j*mu[j,i].
    return [
        -(
            residue[i]
            + sum(residue[j] * mu[j][i] for j in range(i + 1, DIMENSION))
        )
        / DEGREE
        for i in range(DIMENSION)
    ]


def residue_from_group(first_nonzero: int, quotient: int) -> list[int]:
    residue = [0] * DIMENSION
    residue[first_nonzero] = 1
    for index in range(first_nonzero + 1, DIMENSION):
        residue[index] = quotient % DEGREE
        quotient //= DEGREE
    assert quotient == 0
    return residue


def gso_data(gram: list[list[int]], float_type: str, precision: int):
    from fpylll import FPLLL, GSO, IntegerMatrix

    if float_type == "mpfr":
        FPLLL.set_precision(precision)
    matrix = IntegerMatrix.from_matrix(gram)
    gso = GSO.Mat(matrix, gram=True, float_type=float_type, update=True)
    mu = [
        [gso.get_mu(i, j) if i > j else 0.0 for j in range(DIMENSION)]
        for i in range(DIMENSION)
    ]
    return gso, mu


def solve_residue(
    residue: list[int],
    gram: list[list[int]],
    gso,
    mu: list[list[float]],
    distance_bound: float,
) -> tuple[int, list[int], float]:
    from fpylll import Enumeration

    target = target_coordinates(residue, mu)
    solutions = Enumeration(gso).enumerate(
        0, DIMENSION, distance_bound, 0, target=target
    )
    if not solutions:
        raise RuntimeError("CVP enumeration returned no solution")
    reported_distance, coordinates = solutions[0]
    closest = [int(round(value)) for value in coordinates]
    if any(abs(value - integer) > 1e-7 for value, integer in zip(coordinates, closest)):
        raise RuntimeError("CVP coordinates are not numerically integral")
    representative = [
        residue[index] + DEGREE * closest[index] for index in range(DIMENSION)
    ]
    if any(
        (representative[index] - residue[index]) % DEGREE
        for index in range(DIMENSION)
    ):
        raise RuntimeError("CVP representative is in the wrong residue class")
    norm = exact_norm(representative, gram)
    rounding_error = abs(DEGREE * DEGREE * reported_distance - norm)
    if rounding_error > 1e-7:
        raise RuntimeError(
            f"CVP distance/integral-norm mismatch: {rounding_error}"
        )
    return norm, representative, rounding_error


def census_chunk(task: dict) -> dict:
    gram = task["gram"]
    primary_gso, primary_mu = gso_data(
        gram, task["float_type"], task["precision"]
    )
    audit_gso, audit_mu = gso_data(gram, "mpfr", task["audit_precision"])
    # Every residue vector with digits 0,1,2 is itself a valid representative.
    # This uniform absolute-value bound therefore guarantees an initial CVP ball.
    distance_bound = (
        4 * sum(abs(entry) for row in gram for entry in row) / 9.0 + 1.0
    )
    histogram = Counter()
    maximum_rounding_error = 0.0
    audited = 0
    checksum_mod_1 = 0
    checksum_mod_2 = 0
    modulus_1 = 1_000_000_007
    modulus_2 = 1_000_000_009
    power = DEGREE ** (task["first_nonzero"] + 1)
    base = DEGREE ** task["first_nonzero"]
    for quotient in range(task["quotient_start"], task["quotient_stop"]):
        residue = residue_from_group(task["first_nonzero"], quotient)
        norm, representative, rounding_error = solve_residue(
            residue, gram, primary_gso, primary_mu, distance_bound
        )
        if norm < 0 or norm % 2:
            raise RuntimeError(f"invalid even-lattice norm {norm}")
        maximum_rounding_error = max(maximum_rounding_error, rounding_error)
        histogram[norm] += 1
        residue_id = base + power * quotient
        checksum_mod_1 = (checksum_mod_1 + residue_id * (norm + 1)) % modulus_1
        checksum_mod_2 = (checksum_mod_2 + residue_id * (norm + 1)) % modulus_2
        if residue_id % task["audit_stride"] == task["audit_residue"]:
            audit_norm, audit_representative, audit_error = solve_residue(
                residue, gram, audit_gso, audit_mu, distance_bound
            )
            if audit_norm != norm:
                raise RuntimeError(
                    f"cross-precision minimum mismatch {norm} != {audit_norm}"
                )
            if exact_norm(audit_representative, gram) != norm:
                raise RuntimeError("MPFR audit representative norm mismatch")
            maximum_rounding_error = max(maximum_rounding_error, audit_error)
            audited += 1
    return {
        "task_id": task["task_id"],
        "frame_id": task["frame_id"],
        "representatives": task["quotient_stop"] - task["quotient_start"],
        "histogram": {str(key): value for key, value in sorted(histogram.items())},
        "maximum_rounding_error": maximum_rounding_error,
        "cross_precision_audits": audited,
        "checksum_mod_1000000007": checksum_mod_1,
        "checksum_mod_1000000009": checksum_mod_2,
    }


def task_specifications(
    frames: list[dict], arguments: argparse.Namespace
) -> list[dict]:
    tasks = []
    for frame in frames:
        for first_nonzero in range(DIMENSION):
            group_size = DEGREE ** (DIMENSION - first_nonzero - 1)
            for start in range(0, group_size, arguments.chunk_size):
                stop = min(group_size, start + arguments.chunk_size)
                task_id = (
                    f"{frame['frame_id']}:p{first_nonzero}:q{start}-{stop}"
                )
                tasks.append(
                    {
                        "task_id": task_id,
                        "frame_id": frame["frame_id"],
                        "gram": frame["gram"],
                        "first_nonzero": first_nonzero,
                        "quotient_start": start,
                        "quotient_stop": stop,
                        "float_type": arguments.float_type,
                        "precision": arguments.precision,
                        "audit_precision": arguments.audit_precision,
                        "audit_stride": arguments.audit_stride,
                        "audit_residue": int(
                            hashlib.sha256(task_id.encode()).hexdigest()[:8], 16
                        )
                        % arguments.audit_stride,
                    }
                )
    return tasks


def configuration_record(
    database_path: Path,
    frames: list[dict],
    arguments: argparse.Namespace,
    tasks: list[dict],
) -> dict:
    task_text = "\n".join(task["task_id"] for task in tasks) + "\n"
    return {
        "database": relative(database_path),
        "database_sha256": digest(database_path),
        "frame_ids": [frame["frame_id"] for frame in frames],
        "frame_gram_sha256": {
            frame["frame_id"]: gram_digest(frame["gram"]) for frame in frames
        },
        "dimension": DIMENSION,
        "degree": DEGREE,
        "total_cosets_per_frame": TOTAL_COSETS,
        "inversion_representatives_per_frame": TOTAL_INVERSION_REPRESENTATIVES,
        "chunk_size": arguments.chunk_size,
        "float_type": arguments.float_type,
        "precision": arguments.precision,
        "audit_precision": arguments.audit_precision,
        "audit_stride": arguments.audit_stride,
        "task_count": len(tasks),
        "task_list_sha256": hashlib.sha256(task_text.encode()).hexdigest(),
    }


def empty_aggregate(frames: list[dict]) -> dict:
    return {
        frame["frame_id"]: {
            "representatives": 1,  # the zero coset
            "histogram": {"0": 1},
            "maximum_rounding_error": 0.0,
            "cross_precision_audits": 0,
            "checksum_mod_1000000007": 0,
            "checksum_mod_1000000009": 0,
        }
        for frame in frames
    }


def merge_result(aggregate: dict, result: dict) -> None:
    row = aggregate[result["frame_id"]]
    row["representatives"] += result["representatives"]
    for norm, count in result["histogram"].items():
        row["histogram"][norm] = row["histogram"].get(norm, 0) + count
    row["maximum_rounding_error"] = max(
        row["maximum_rounding_error"], result["maximum_rounding_error"]
    )
    row["cross_precision_audits"] += result["cross_precision_audits"]
    for modulus in (1_000_000_007, 1_000_000_009):
        key = f"checksum_mod_{modulus}"
        row[key] = (row[key] + result[key]) % modulus


def write_checkpoint(
    path: Path,
    configuration: dict,
    completed: set[str],
    aggregate: dict,
) -> None:
    payload = {
        "schema": "elkies-k3.lattice-foundry-degree3-complete-checkpoint.v1",
        "configuration": configuration,
        "completed_task_ids": sorted(completed),
        "aggregate": aggregate,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def full_histogram(representative_histogram: dict[str, int]) -> dict[str, int]:
    result = {}
    for norm, count in representative_histogram.items():
        value = int(norm)
        result[norm] = count if value == 0 else 2 * count
    return dict(sorted(result.items(), key=lambda row: int(row[0])))


def spectrum_row(frame: dict, aggregate: dict) -> dict:
    if aggregate["representatives"] != TOTAL_INVERSION_REPRESENTATIVES:
        raise RuntimeError(
            f"incomplete representative coverage for {frame['frame_id']}"
        )
    histogram = full_histogram(aggregate["histogram"])
    if sum(histogram.values()) != TOTAL_COSETS:
        raise RuntimeError(f"full coset count mismatch for {frame['frame_id']}")
    maximum = max(map(int, histogram))
    rational = sum(
        count
        for norm, count in histogram.items()
        if int(norm) >= 20 and (int(norm) - 2) % 6 == 0
    )
    genus_one = sum(
        count
        for norm, count in histogram.items()
        if int(norm) >= 18 and int(norm) % 6 == 0
    )
    return {
        "frame_id": frame["frame_id"],
        "ns_id": frame["ns_id"],
        "determinant": int(frame["determinant"]),
        "published_R17_control": frame["frame_id"] == PUBLISHED_R17_CONTROL,
        "pre_prescribed_route_aware_top_five": (
            frame["frame_id"] in PRE_PRESCRIBED_ROUTE_AWARE_TOP_FIVE
        ),
        "current_mw2_source_ranked_top_five": (
            frame["frame_id"] in CURRENT_MW2_SOURCE_RANKED_TOP_FIVE
        ),
        "automorphism_group_order": int(
            frame["rootless_intrinsics"]["automorphism_group_order"]
        ),
        "translation_cosets": TOTAL_COSETS,
        "inversion_orbits": TOTAL_INVERSION_REPRESENTATIVES,
        "minimum_norm_histogram_all_translation_cosets": histogram,
        "maximum_coset_minimum_norm": maximum,
        "rational_trisection_translation_cosets": rational,
        "genus_one_trisection_translation_cosets": genus_one,
        "geometric_boundary": (
            "Counts impose integrality and nonnegative intersection with every "
            "section. Effectivity, global nefness, irreducibility, arithmetic "
            "descent, and specialization rank gain remain separate gates."
        ),
        "numerical_certificate": {
            "every_returned_candidate_norm_recomputed_integrally": True,
            "maximum_distance_to_integral_norm_error": aggregate[
                "maximum_rounding_error"
            ],
            "cross_precision_mpfr_audits": aggregate[
                "cross_precision_audits"
            ],
            "checksum_mod_1000000007": aggregate[
                "checksum_mod_1000000007"
            ],
            "checksum_mod_1000000009": aggregate[
                "checksum_mod_1000000009"
            ],
        },
    }


def selected_frames(database: dict, frame_ids: list[str]) -> list[dict]:
    by_id = {
        frame["frame_id"]: {"ns_id": ns["ns_id"], **frame}
        for ns in database["ns_classes"]
        for frame in ns["frames"]
    }
    missing = set(frame_ids) - set(by_id)
    if missing:
        raise ValueError(f"unknown frame ids: {sorted(missing)}")
    frames = [by_id[frame_id] for frame_id in frame_ids]
    for frame in frames:
        if int(frame["root_rank"]) != 0 or int(frame["mw_rank_for_rho_19"]) != 17:
            raise ValueError(f"frame is not rootless MW17: {frame['frame_id']}")
        if int(frame["rootless_intrinsics"]["minimum_squared_norm"]) < 4:
            raise ValueError(f"frame has unexpected minimum: {frame['frame_id']}")
    return frames


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DATABASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument(
        "--reuse-checkpoint",
        type=Path,
        action="append",
        default=[],
        help="reuse fully completed frame blocks from a compatible checkpoint",
    )
    parser.add_argument("--frame-id", action="append", default=[])
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--chunk-size", type=int, default=1_000_000)
    parser.add_argument("--float-type", choices=("double", "long double", "dd", "qd", "mpfr"), default="dd")
    parser.add_argument("--precision", type=int, default=160)
    parser.add_argument("--audit-precision", type=int, default=256)
    parser.add_argument("--audit-stride", type=int, default=4096)
    parser.add_argument("--smoke-count", type=int, default=0)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "verify the completed checkpoint, all coverage identities, and the "
            "stored final spectra without rerunning the full CVP census"
        ),
    )
    arguments = parser.parse_args()
    if arguments.workers <= 0 or arguments.chunk_size <= 0:
        parser.error("workers and chunk size must be positive")
    if arguments.audit_stride <= 0:
        parser.error("audit stride must be positive")

    database_path = arguments.database.resolve()
    database = json.loads(database_path.read_text())
    frame_ids = arguments.frame_id or DEFAULT_FRAMES
    frames = selected_frames(database, frame_ids)
    tasks = task_specifications(frames, arguments)
    if arguments.smoke_count:
        task = dict(tasks[0])
        task["quotient_stop"] = min(
            task["quotient_stop"], task["quotient_start"] + arguments.smoke_count
        )
        task["task_id"] += f":smoke{arguments.smoke_count}"
        started = time.time()
        result = census_chunk(task)
        print(
            "FOUNDRYDEG3COMPLETE|stage=smoke|"
            f"frame={result['frame_id']}|representatives={result['representatives']}|"
            f"seconds={time.time()-started:.3f}|histogram={result['histogram']}|status=PASS",
            flush=True,
        )
        return

    configuration = configuration_record(database_path, frames, arguments, tasks)
    output_path = arguments.output.resolve()
    checkpoint_path = (
        arguments.checkpoint.resolve()
        if arguments.checkpoint is not None
        else Path(str(output_path) + ".partial")
    )
    aggregate = empty_aggregate(frames)
    completed: set[str] = set()
    if checkpoint_path.is_file() and not arguments.rebuild:
        checkpoint = json.loads(checkpoint_path.read_text())
        if checkpoint.get("configuration") != configuration:
            raise SystemExit("degree-three checkpoint configuration mismatch")
        completed = set(checkpoint["completed_task_ids"])
        aggregate = checkpoint["aggregate"]
    elif arguments.reuse_checkpoint and not arguments.rebuild:
        compatibility_keys = (
            "database_sha256",
            "dimension",
            "degree",
            "total_cosets_per_frame",
            "inversion_representatives_per_frame",
            "chunk_size",
            "float_type",
            "precision",
            "audit_precision",
            "audit_stride",
        )
        tasks_by_frame = {
            frame["frame_id"]: {
                task["task_id"]
                for task in tasks
                if task["frame_id"] == frame["frame_id"]
            }
            for frame in frames
        }
        for reuse_path_argument in arguments.reuse_checkpoint:
            reuse_path = reuse_path_argument.resolve()
            reuse = json.loads(reuse_path.read_text())
            if reuse.get("schema") != (
                "elkies-k3.lattice-foundry-degree3-complete-checkpoint.v1"
            ):
                raise SystemExit(f"invalid reusable checkpoint: {reuse_path}")
            reuse_configuration = reuse["configuration"]
            if any(
                reuse_configuration.get(key) != configuration.get(key)
                for key in compatibility_keys
            ):
                raise SystemExit(f"incompatible reusable checkpoint: {reuse_path}")
            reuse_completed = set(reuse["completed_task_ids"])
            for frame_id, frame_task_ids in tasks_by_frame.items():
                if not frame_task_ids <= reuse_completed:
                    continue
                if (
                    reuse_configuration["frame_gram_sha256"].get(frame_id)
                    != configuration["frame_gram_sha256"][frame_id]
                ):
                    raise SystemExit(
                        f"reusable checkpoint Gram mismatch for {frame_id}"
                    )
                aggregate[frame_id] = reuse["aggregate"][frame_id]
                completed.update(frame_task_ids)
        if completed:
            write_checkpoint(checkpoint_path, configuration, completed, aggregate)

    if arguments.check:
        task_ids = {task["task_id"] for task in tasks}
        if completed != task_ids:
            raise SystemExit(
                f"incomplete degree-three checkpoint: {len(completed)}/{len(task_ids)} tasks"
            )
        expected_spectra = [
            spectrum_row(frame, aggregate[frame["frame_id"]]) for frame in frames
        ]
        stored = json.loads(output_path.read_text())
        if stored.get("schema") != "elkies-k3.lattice-foundry-degree3-complete-spectrum.v1":
            raise SystemExit("stored degree-three artifact has the wrong schema")
        if stored.get("spectra") != expected_spectra:
            raise SystemExit("stored degree-three spectra do not match the checkpoint")
        if stored.get("method", {}).get("task_list_sha256") != configuration["task_list_sha256"]:
            raise SystemExit("stored degree-three task digest is stale")
        if stored.get("inputs", {}).get(relative(database_path)) != digest(database_path):
            raise SystemExit("stored degree-three database digest is stale")
        print(
            "FOUNDRYDEG3COMPLETE|stage=check|"
            f"frames={len(frames)}|tasks={len(tasks)}|"
            f"cosets={len(frames)*TOTAL_COSETS}|status=PASS",
            flush=True,
        )
        return

    pending = [task for task in tasks if task["task_id"] not in completed]
    started = time.time()
    last_checkpoint = time.time()
    with ProcessPoolExecutor(
        max_workers=arguments.workers,
        mp_context=multiprocessing.get_context("fork"),
    ) as executor:
        future_to_task = {
            executor.submit(census_chunk, task): task for task in pending
        }
        for index, future in enumerate(as_completed(future_to_task), 1):
            result = future.result()
            if result["task_id"] in completed:
                raise RuntimeError("duplicate completed degree-three task")
            merge_result(aggregate, result)
            completed.add(result["task_id"])
            now = time.time()
            if now - last_checkpoint >= 30 or index == len(pending):
                write_checkpoint(
                    checkpoint_path, configuration, completed, aggregate
                )
                last_checkpoint = now
            if index % 10 == 0 or index == len(pending):
                representatives = sum(
                    row["representatives"] for row in aggregate.values()
                )
                total = len(frames) * TOTAL_INVERSION_REPRESENTATIVES
                print(
                    "FOUNDRYDEG3COMPLETE|stage=progress|"
                    f"tasks={len(completed)}/{len(tasks)}|"
                    f"representatives={representatives}/{total}|"
                    f"elapsed_seconds={now-started:.1f}|status=RUNNING",
                    flush=True,
                )

    spectra = [spectrum_row(frame, aggregate[frame["frame_id"]]) for frame in frames]
    output = {
        "schema": "elkies-k3.lattice-foundry-degree3-complete-spectrum.v1",
        "status": "PASS_COMPLETE_ALL_DEGREE3_COSETS_CROSS_PRECISION_AUDITED",
        "selection": {
            "pre_prescribed_route_aware_top_five_frame_ids": [
                frame_id
                for frame_id in PRE_PRESCRIBED_ROUTE_AWARE_TOP_FIVE
                if frame_id in frame_ids
            ],
            "current_mw2_source_ranked_top_five_frame_ids": [
                frame_id
                for frame_id in CURRENT_MW2_SOURCE_RANKED_TOP_FIVE
                if frame_id in frame_ids
            ],
            "published_R17_control_frame_id": (
                PUBLISHED_R17_CONTROL if PUBLISHED_R17_CONTROL in frame_ids else None
            ),
            "boundary": (
                "The default frames are the leading rootless-target surfaces in the "
                "pre-prescribed-root, route-aware source ledger and were already in "
                "the 256-coset pilot. The current MW2 source-ranked top five are "
                "identified separately when present. The formula used here requires "
                "NS=U+M(-1)."
            ),
        },
        "method": {
            **configuration,
            "inversion_representative_rule": (
                "zero plus ternary vectors whose first nonzero coordinate is one"
            ),
            "primary_cvp": (
                "complete fplll enumeration on the integral Gram GSO with exact "
                "integral recomputation of every returned candidate norm"
            ),
            "cross_precision_audit": (
                "deterministic residue subset independently repeated with MPFR GSO"
            ),
        },
        "proof_boundary": {
            "complete": (
                "Every residue in M/3M is represented exactly once after restoring "
                "inversion pairs; histogram totals equal 3^17 for every frame."
            ),
            "numerical": (
                "CVP branch decisions use the declared fplll floating backend. "
                "Every candidate norm is integral-exact and a deterministic subset "
                "is independently MPFR-audited; this is a complete computational "
                "census, not a formal proof of the floating enumerator."
            ),
            "geometric": (
                "Coset minima certify the lattice and all-section gates only. Curve "
                "effectivity, nefness, irreducibility, descent, and rank gain are open."
            ),
        },
        "spectra": spectra,
        "inputs": {relative(database_path): digest(database_path)},
        "reproduce": (
            f"{Path(os.sys.executable)} {relative(Path(__file__))} "
            + " ".join(f"--frame-id {frame_id}" for frame_id in frame_ids)
            + " "
            f"--workers {arguments.workers} --chunk-size {arguments.chunk_size} "
            f"--float-type '{arguments.float_type}' --audit-precision "
            f"{arguments.audit_precision} --audit-stride {arguments.audit_stride} "
            f"--output {relative(output_path)}"
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        "FOUNDRYDEG3COMPLETE|stage=done|"
        f"frames={len(frames)}|cosets_per_frame={TOTAL_COSETS}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
