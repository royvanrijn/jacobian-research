#!/usr/bin/env python3
"""Calibrate and run the bounded alternate-R17-in-rank29 embedding gate.

This replaces public-point subset selection by a search over integral vectors
in the complete PARI-enumerated short-vector ball of the displayed rank-29
height lattice.  A negative result is explicitly bounded by the norm cutoff
and trial count; a numerical fit is never promoted to an exact lattice
embedding or a family identification.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import random
import re
import subprocess
import sys
import time

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ELLIPTIC_ROOT = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC_ROOT), str(ELLIPTIC_ROOT / "cas")]

from compare_record_height_lattices import gp_vector  # noqa: E402
from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
)
from elliptic_candidate_record import source_point_to_target  # noqa: E402


INVARIANTS = ROOT / "artifacts/generated-results/elkies-k3-other-rank17-invariants.json"
HEIGHT_LATTICES = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_height_lattices_28_29_273_302_v1.json"
)
POSITIVE_CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_high_rank_positive_controls_v2.json"
)
PUBLISHED_TARGET = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
PINNED_GRAM = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-other-rank17-rank29-gate-a.json"

PARAMETERS = (
    ("rank25", -2, 377),
    ("rank26", -308, 251),
    ("rank27", 2456, 135),
    ("rank28", -9529, 5471),
)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gp_matrix_assignments(name: str, values) -> str:
    rows = list(values)
    result = [f"{name}=matrix({len(rows)},{len(rows[0])});"]
    for index, row in enumerate(rows, 1):
        result.append(f"{name}[{index},]=[" + ",".join(map(str, row)) + "];")
    return "\n".join(result)


def run_gp(program: str) -> list[str]:
    completed = subprocess.run(
        ["gp", "-q"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    stderr = re.sub(r"\x1b\[[0-9;]*m", "", completed.stderr)
    failures = [
        line
        for line in stderr.splitlines()
        if line.strip()
        and "Warning: new" not in line
        and "Warning: increasing stack" not in line
    ]
    if failures:
        raise RuntimeError("\n".join(failures))
    return [line.strip() for line in completed.stdout.splitlines() if line.strip()]


def parse_matrix(lines: list[str], begin: str, end: str) -> np.ndarray:
    start = lines.index(begin) + 1
    stop = lines.index(end, start)
    return np.array(
        [[float(value) for value in row.split("|")] for row in lines[start:stop]],
        dtype=np.float64,
    )


def published_section_gram() -> np.ndarray:
    pinned = np.loadtxt(PINNED_GRAM, dtype=np.int64)
    target = json.loads(PUBLISHED_TARGET.read_text())
    change = np.array(
        target["pinned_identification"]["basis_change_matrix"], dtype=np.int64
    )
    inverse = np.rint(np.linalg.inv(change)).astype(np.int64)
    if not np.array_equal(change @ inverse, np.eye(17, dtype=np.int64)):
        raise ArithmeticError("published/pinned basis change is not unimodular")
    gram = inverse.T @ pinned @ inverse
    if int(round(np.linalg.det(gram))) != 948 or not np.all(np.diag(gram) == 4):
        raise ArithmeticError("reconstructed published-section Gram changed")
    return gram


def specialized_generic_height_matrices() -> dict[str, np.ndarray]:
    data = load_q12o5867_data(MODEL, SECTIONS)
    result = {}
    for label, numerator, denominator in PARAMETERS:
        specialization = evaluate_projective_specialization(
            data, numerator, denominator
        )
        minimal_model, minimal_change, _ = global_minimal_model_with_change(
            specialization.model
        )
        points = tuple(
            source_point_to_target(point, minimal_change)
            for point in specialization.points
        )
        program = f"""
default(realprecision,100);
E=ellinit({gp_vector(minimal_model)});
P=[{','.join(gp_vector(point) for point in points)}];
H=ellheightmatrix(E,P);
print("BEGIN_HEIGHT");
for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(H[i,j]));print());
print("END_HEIGHT");
"""
        result[label] = parse_matrix(
            run_gp(program), "BEGIN_HEIGHT", "END_HEIGHT"
        )
    return result


def fit_metrics(observed: np.ndarray, target: np.ndarray) -> dict[str, float]:
    scale = float(np.sum(observed * target) / np.sum(target * target))
    residual = observed - scale * target
    normalized_entries = np.abs(residual).ravel() / scale
    return {
        "scale": scale,
        "relative_frobenius": float(
            np.linalg.norm(residual, "fro")
            / (scale * np.linalg.norm(target, "fro"))
        ),
        "maximum_absolute_residual_over_scale": float(normalized_entries.max()),
        "residual_over_scale_median": float(np.quantile(normalized_entries, 0.5)),
        "residual_over_scale_q90": float(np.quantile(normalized_entries, 0.9)),
        "residual_over_scale_q95": float(np.quantile(normalized_entries, 0.95)),
        "diagonal_ratio_minimum": float(
            np.min(np.diag(observed) / (scale * np.diag(target)))
        ),
        "diagonal_ratio_maximum": float(
            np.max(np.diag(observed) / (scale * np.diag(target)))
        ),
    }


def exact_integer_rank(values: np.ndarray) -> int:
    """Return the exact rational rank of a small integral matrix."""
    rows = [
        [Fraction(int(entry)) for entry in row]
        for row in values.tolist()
    ]
    rank = 0
    for column in range(len(rows[0])):
        pivot = next(
            (index for index in range(rank, len(rows)) if rows[index][column]),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [entry / pivot_value for entry in rows[rank]]
        for index in range(len(rows)):
            if index == rank or not rows[index][column]:
                continue
            multiplier = rows[index][column]
            rows[index] = [
                entry - multiplier * pivot_entry
                for entry, pivot_entry in zip(rows[index], rows[rank])
            ]
        rank += 1
        if rank == min(len(rows), len(rows[0])):
            break
    return rank


def enumerate_rank29_short_lines(
    record: dict[str, object], bound: float
) -> tuple[np.ndarray, np.ndarray, int]:
    reduced = record["lll_reduced_gram"]
    program = f"""
default(parisizemax,2000000000);
default(parisize,200000000);
default(realprecision,100);
{gp_matrix_assignments('R', reduced)}
Q=qfminim(R,{bound},1000000,2);
print("SIGNED|",Q[1]);
print("BEGIN_SHORT");
for(j=1,matsize(Q[3])[2],for(i=1,matsize(Q[3])[1],if(i>1,print1("|"));print1(Q[3][i,j]));print("|",Q[3][,j]~*R*Q[3][,j]));
print("END_SHORT");
"""
    lines = run_gp(program)
    signed = int(
        next(line.split("|", 1)[1] for line in lines if line.startswith("SIGNED|"))
    )
    start = lines.index("BEGIN_SHORT") + 1
    stop = lines.index("END_SHORT", start)
    reduced_vectors = []
    norms = []
    for line in lines[start:stop]:
        fields = line.split("|")
        reduced_vectors.append(list(map(int, fields[:-1])))
        norms.append(float(fields[-1]))
    if signed != 2 * len(reduced_vectors):
        raise ArithmeticError("PARI short-vector sign count changed")
    transform = np.array(record["lll_transform_columns"], dtype=np.int64)
    public_vectors = (transform @ np.array(reduced_vectors, dtype=np.int64).T).T
    height = np.array(record["height_gram"], dtype=np.float64)
    replayed_norms = np.einsum("ij,jk,ik->i", public_vectors, height, public_vectors)
    if np.max(np.abs(replayed_norms - np.array(norms))) > 1e-8:
        raise ArithmeticError("short-vector coordinates do not replay in public basis")
    return public_vectors, np.array(norms), signed


def target_order(gram: np.ndarray) -> list[int]:
    first = max(
        range(len(gram)),
        key=lambda index: (
            np.sum(np.abs(gram[index]) == 2),
            np.sum(np.abs(gram[index])),
        ),
    )
    order = [first]
    remaining = set(range(len(gram))) - {first}
    while remaining:
        following = max(
            remaining,
            key=lambda index: (
                sum(abs(int(gram[index, old])) == 2 for old in order),
                sum(abs(int(gram[index, old])) for old in order),
            ),
        )
        order.append(following)
        remaining.remove(following)
    return order


def embedding_search(
    height: np.ndarray,
    pool: np.ndarray,
    norms: np.ndarray,
    target: np.ndarray,
    *,
    trials: int,
    seed: int,
    top_choices: int,
    refinement_passes: int,
) -> tuple[dict[str, object], list[float]]:
    rng = random.Random(seed)
    order = target_order(target)
    pool_height = pool @ height
    trial_scores = []
    best = None

    # Cover the whole norm interval deterministically before random anchors.
    anchor_order = np.argsort(norms)
    anchors = [
        int(anchor_order[min(len(anchor_order) - 1, index * len(anchor_order) // trials)])
        for index in range(min(trials, len(anchor_order)))
    ]
    while len(anchors) < trials:
        anchors.append(rng.randrange(len(pool)))
    rng.shuffle(anchors)

    for trial, anchor in enumerate(anchors):
        selected = {order[0]: (anchor, 1)}
        used = {anchor}
        scale = norms[anchor] / target[order[0], order[0]]

        for target_index in order[1:]:
            prior = list(selected)
            oriented = np.array(
                [sign * pool[index] for index, sign in (selected[item] for item in prior)],
                dtype=np.int64,
            )
            pairings = pool_height @ oriented.T
            desired = scale * target[target_index, prior]
            plus = np.sum((pairings - desired) ** 2, axis=1)
            minus = np.sum((-pairings - desired) ** 2, axis=1)
            diagonal = (norms - scale * target[target_index, target_index]) ** 2 / 2
            scores = np.minimum(plus, minus) + diagonal
            signs = np.where(plus <= minus, 1, -1)
            if used:
                scores[list(used)] = np.inf
            candidates = np.argpartition(scores, min(top_choices, len(scores) - 1))[
                :top_choices
            ]
            candidates = candidates[np.argsort(scores[candidates])]
            # Mostly greedy, with bounded diversity among near-best candidates.
            rank = 0 if rng.random() < 0.7 else rng.randrange(len(candidates))
            choice = int(candidates[rank])
            selected[target_index] = (choice, int(signs[choice]))
            used.add(choice)

        columns = np.zeros((height.shape[0], target.shape[0]), dtype=np.int64)
        for target_index, (index, sign) in selected.items():
            columns[:, target_index] = sign * pool[index]

        def score_columns(value):
            observed = value.T @ height @ value
            metrics = fit_metrics(observed, target)
            return metrics["relative_frobenius"], observed, metrics

        current_score, observed, metrics = score_columns(columns)

        # Exact-coordinate local replacement in the enumerated pool.  Each
        # move is accepted only when the full 17-by-17 objective improves.
        for _ in range(refinement_passes):
            changed = False
            scale = metrics["scale"]
            for target_index in order:
                other = [index for index in range(17) if index != target_index]
                other_columns = columns[:, other].T
                pairings = pool_height @ other_columns.T
                desired = scale * target[target_index, other]
                plus = np.sum((pairings - desired) ** 2, axis=1)
                minus = np.sum((-pairings - desired) ** 2, axis=1)
                diagonal = (
                    norms - scale * target[target_index, target_index]
                ) ** 2 / 2
                local_scores = np.minimum(plus, minus) + diagonal
                signs = np.where(plus <= minus, 1, -1)
                occupied = {
                    selected[index][0]
                    for index in selected
                    if index != target_index
                }
                if occupied:
                    local_scores[list(occupied)] = np.inf
                candidates = np.argpartition(local_scores, min(8, len(local_scores) - 1))[:8]
                candidates = candidates[np.argsort(local_scores[candidates])]
                accepted = False
                for candidate in candidates:
                    proposal = columns.copy()
                    proposal[:, target_index] = int(signs[candidate]) * pool[candidate]
                    proposal_score, proposal_observed, proposal_metrics = score_columns(proposal)
                    if proposal_score + 1e-14 < current_score:
                        columns = proposal
                        current_score = proposal_score
                        observed = proposal_observed
                        metrics = proposal_metrics
                        selected[target_index] = (int(candidate), int(signs[candidate]))
                        changed = accepted = True
                        break
                if accepted:
                    continue
            if not changed:
                break

        exact_rank = exact_integer_rank(columns)
        if exact_rank < 17:
            continue
        trial_scores.append(current_score)
        record = {
            "trial": trial,
            "metrics": metrics,
            "embedding_columns_public_point_basis": columns.tolist(),
            "observed_gram": observed.tolist(),
            "residual": (observed - metrics["scale"] * target).tolist(),
            "maximum_embedding_coefficient": int(np.max(np.abs(columns))),
            "selected_vector_norms": [
                float(columns[:, index] @ height @ columns[:, index])
                for index in range(17)
            ],
            "rational_rank": exact_rank,
        }
        if best is None or current_score < best["metrics"]["relative_frobenius"]:
            best = record

    if best is None:
        raise ArithmeticError("all bounded embedding trials were rank-deficient")
    return best, trial_scores


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--short-bound", type=float, default=65.0)
parser.add_argument("--trials", type=int, default=1000)
parser.add_argument("--seed", type=int, default=291317)
parser.add_argument("--top-choices", type=int, default=8)
parser.add_argument("--refinement-passes", type=int, default=2)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

start_time = time.time()
invariants = json.loads(INVARIANTS.read_text())
alternate = np.array(invariants["frame"]["reduced_gram"], dtype=np.int64)
published = published_section_gram()

control_heights = specialized_generic_height_matrices()
control_metrics = {
    label: fit_metrics(height, published)
    for label, height in control_heights.items()
}
calibration = {
    "fibres": control_metrics,
    "relative_frobenius_envelope": {
        "minimum": min(item["relative_frobenius"] for item in control_metrics.values()),
        "maximum": max(item["relative_frobenius"] for item in control_metrics.values()),
    },
    "maximum_absolute_residual_over_scale_envelope": {
        "minimum": min(
            item["maximum_absolute_residual_over_scale"]
            for item in control_metrics.values()
        ),
        "maximum": max(
            item["maximum_absolute_residual_over_scale"]
            for item in control_metrics.values()
        ),
    },
    "interpretation": (
        "These four matrices use the known specialized generic seventeen in "
        "their exact published basis. The envelopes calibrate specialization "
        "distortion; they are not universal error bounds."
    ),
}

height_data = json.loads(HEIGHT_LATTICES.read_text())
rank29 = next(item for item in height_data["curves"] if item["label"] == "rank29")
height29 = np.array(rank29["height_gram"], dtype=np.float64)
pool, norms, signed_count = enumerate_rank29_short_lines(rank29, arguments.short_bound)
best, trial_scores = embedding_search(
    height29,
    pool,
    norms,
    alternate,
    trials=arguments.trials,
    seed=arguments.seed,
    top_choices=arguments.top_choices,
    refinement_passes=arguments.refinement_passes,
)

within_frobenius_envelope = (
    best["metrics"]["relative_frobenius"]
    <= calibration["relative_frobenius_envelope"]["maximum"]
)
within_entry_envelope = (
    best["metrics"]["maximum_absolute_residual_over_scale"]
    <= calibration["maximum_absolute_residual_over_scale_envelope"]["maximum"]
)
passes_calibrated_gate = within_frobenius_envelope and within_entry_envelope

payload = {
    "schema": "elkies-k3.other-rank17-rank29-gate-a.v1",
    "status": (
        "BOUNDED_NUMERICAL_CANDIDATE_PASSES_CONTROL_ENVELOPE"
        if passes_calibrated_gate
        else "NO_BOUNDED_NUMERICAL_EMBEDDING_IN_CONTROL_ENVELOPE"
    ),
    "calibration": calibration,
    "rank29_search": {
        "ambient_dimension": 29,
        "target_dimension": 17,
        "target": "alternate determinant-948 rootless frame",
        "short_vector_squared_height_bound": arguments.short_bound,
        "enumerated_unoriented_short_lines": len(pool),
        "enumerated_signed_short_vectors": signed_count,
        "trials": arguments.trials,
        "seed": arguments.seed,
        "top_choices": arguments.top_choices,
        "refinement_passes": arguments.refinement_passes,
        "best": best,
        "trial_relative_frobenius_quantiles": {
            str(quantile): float(np.quantile(trial_scores, quantile))
            for quantile in (0, 0.01, 0.1, 0.5, 0.9, 1)
        },
        "within_control_frobenius_envelope": within_frobenius_envelope,
        "within_control_maximum_entry_envelope": within_entry_envelope,
        "passes_both_calibrated_gates": passes_calibrated_gate,
    },
    "proof_boundary": {
        "proved": (
            "The four positive-control distortion measurements and every tested "
            "integer matrix are replayable from pinned height matrices."
        ),
        "not_proved": (
            "The search is bounded by its short-vector norm and trial count. A "
            "passing approximate fit would still not prove that rank29 is a "
            "specialization; a failing search is not a non-embedding theorem."
        ),
    },
    "inputs": {
        str(path.relative_to(ROOT)): file_sha256(path)
        for path in (
            INVARIANTS,
            HEIGHT_LATTICES,
            POSITIVE_CONTROLS,
            PUBLISHED_TARGET,
            PINNED_GRAM,
            MODEL,
            SECTIONS,
        )
    },
    "runtime_seconds": time.time() - start_time,
    "reproduce": (
        "python3 elkies-k3/scripts/gate_a_alternate_rank17_rank29.py "
        f"--short-bound {arguments.short_bound:g} --trials {arguments.trials} "
        f"--seed {arguments.seed} --top-choices {arguments.top_choices} "
        f"--refinement-passes {arguments.refinement_passes}"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not arguments.output.exists():
        raise SystemExit("Gate A artifact is missing")
    existing = json.loads(arguments.output.read_text())
    # Runtime is informational and intentionally excluded from replay equality.
    existing.pop("runtime_seconds", None)
    replay = dict(payload)
    replay.pop("runtime_seconds", None)
    if existing != replay:
        raise SystemExit("Gate A artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(serialized)

print(
    "OTHERR17GATEA|controls_relF={:.6f}-{:.6f}|pool={}|trials={}|"
    "best_relF={:.6f}|best_max_over_h={:.6f}|passes={}|status={}".format(
        calibration["relative_frobenius_envelope"]["minimum"],
        calibration["relative_frobenius_envelope"]["maximum"],
        len(pool),
        arguments.trials,
        best["metrics"]["relative_frobenius"],
        best["metrics"]["maximum_absolute_residual_over_scale"],
        int(passes_calibrated_gate),
        payload["status"],
    ),
    flush=True,
)
