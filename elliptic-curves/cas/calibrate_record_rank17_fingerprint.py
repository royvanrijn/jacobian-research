#!/usr/bin/env python3
"""Calibrate the numerical R17 fingerprint against a known negative control.

The four record-curve candidates were selected by a bounded rank-17
short-vector search.  That selection can overfit.  This script applies the
same construction to ICARM curve 245, whose exact Fermigier--Mestre rank-12
parent is independently reconstructed in this repository, and fits an exact
unimodular basis of the determinant-948 R17 lattice to every numerical Gram.

The output is a calibration experiment, not a K3-specialization certificate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import re
import subprocess

import numpy as np

from compare_record_height_lattices import CURVES, CurveData, gp_vector
import icarm_curve245
from elliptic_candidate_record import (
    WeierstrassChange,
    change_weierstrass_model,
    fraction_text,
    source_point_to_target,
)
from icarm_curve245_mestre import (
    ANCHOR_SHORT_TO_PUBLIC_CHANGE,
    CANONICAL_PARAMETER,
    CONSTRUCTION,
    PUBLIC_MODEL,
    extra_quartic_point,
    primitive_short_model,
)
from nagao_1994 import primitive_visible_points, quartic_point_to_short_jacobian
from search_record_rank17_core import (
    canonical_line,
    enumerate_vectors,
    exact_in_span_mask,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_rank17_core_candidates_v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_rank17_fingerprint_calibration_v1.json"
)
R17_PATH = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"

# Deterministic output of the declared additive-RANSAC search
#
#   search_record_rank17_core.py curve245 --bound 28 \
#     --additive-pair-limit 1800 --pool 1800 --trials 1000 \
#     --dimension 17 --seed 24517
#
# Rows are in the twenty public-point coordinates.  The subsequent double
# integer kernel saturates their rational span; the rows themselves are not
# asserted to be a distinguished subgroup.
CURVE245_SEED_BASIS = (
    (-6, 1, -7, 1, 2, -5, 6, 1, 3, 0, -5, 1, 1, -2, 2, -1, -1, 1, 5, 0),
    (-3, 1, -3, 1, 1, -2, 3, 0, 0, 0, -3, 0, 1, 0, 1, 0, 0, 1, 1, 0),
    (-2, 1, -2, 0, 0, -1, 2, 0, 0, 1, -2, 1, 1, 0, 0, 0, 0, 0, 0, 1),
    (-1, 1, -1, 1, 1, -1, 0, 0, 1, 0, -1, 0, 0, -1, 0, -1, 0, 1, 1, 0),
    (2, 0, 2, 0, -1, 1, -2, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, 0, -1, 0),
    (2, 0, 3, -1, -1, 2, -2, 0, -3, 0, 1, 0, 0, 2, -1, 1, 1, -1, -3, 1),
    (5, -1, 6, -1, -2, 5, -4, -1, -4, 0, 4, -1, -1, 2, -1, 1, 2, -2, -4, 0),
    (1, 0, 2, 0, 0, 1, -1, 0, -1, -1, 2, -1, 0, 0, 0, 0, 0, -1, 0, 0),
    (1, 0, 2, 0, -1, 2, 0, -1, -2, 0, 1, -1, 0, 1, 0, 0, 1, -1, -1, 0),
    (-2, 0, -1, 1, 0, -1, 2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0),
    (1, 0, 1, 0, 0, 0, -1, 1, 1, -1, 1, 0, 0, -1, 0, -1, -1, -1, 1, 0),
    (-2, 0, -2, 0, 0, -1, 2, 0, 1, 1, -2, 1, 0, 0, 0, 0, 0, 1, 0, 0),
    (-1, 0, -1, 1, 1, -1, 1, -1, 1, 0, 0, 0, 0, -1, 0, 0, 0, 1, 1, -1),
    (-3, 0, -3, 0, 1, -2, 3, 0, 1, 1, -3, 1, 0, 0, 0, 0, 0, 1, 1, 0),
    (-1, 0, -2, 0, 1, -2, 0, 1, 3, 0, -1, 1, 0, -2, 0, -1, -1, 1, 2, 0),
    (-1, 0, -1, 0, 0, -1, 1, 0, 1, 0, -1, 1, 0, 0, 0, 0, 0, 1, 0, 0),
    (-2, 0, -1, 1, 1, -1, 1, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 1, 1, 0),
)


def gp_matrix(rows: list[list[int]] | tuple[tuple[int, ...], ...]) -> str:
    return "[" + ";".join(",".join(map(str, row)) for row in rows) + "]"


def parse_matrix(lines: list[str], begin: str, end: str, cast) -> list[list]:
    start = lines.index(begin) + 1
    stop = lines.index(end, start)
    return [[cast(value) for value in line.split("|")] for line in lines[start:stop]]


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
    errors = [line for line in stderr.splitlines() if line.strip() and "Warning: new" not in line]
    if errors:
        raise RuntimeError("\n".join(errors))
    return [line.strip() for line in completed.stdout.splitlines() if line.strip()]


def exact_r17_reduced_form_and_shell() -> tuple[np.ndarray, np.ndarray]:
    rows = [list(map(int, line.split())) for line in R17_PATH.read_text().splitlines()]
    program = f"""
R={gp_matrix(rows)};
U=qflllgram(R);Q=U~*R*U;M=qfminim(Q,4,100000,2);V=M[3];
print("DET|",matdet(Q));print("SIGNED_MINIMAL_VECTORS|",M[1]);
print("BEGIN_Q");
for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(Q[i,j]));print());
print("END_Q");
print("BEGIN_SHELL");
for(j=1,matsize(V)[2],for(i=1,17,if(i>1,print1("|"));print1(V[i,j]));print());
print("END_SHELL");
"""
    lines = run_gp(program)
    if "DET|948" not in lines or "SIGNED_MINIMAL_VECTORS|2622" not in lines:
        raise RuntimeError("pinned R17 determinant or minimal shell changed")
    form = np.array(parse_matrix(lines, "BEGIN_Q", "END_Q", int), dtype=np.int64)
    unoriented = np.array(
        parse_matrix(lines, "BEGIN_SHELL", "END_SHELL", int), dtype=np.int64
    )
    return form, np.vstack((unoriented, -unoriented))


def curve245_data() -> CurveData:
    return CurveData(
        "curve245-negative-control",
        icarm_curve245.GENERAL_WEIERSTRASS_COEFFICIENTS,
        icarm_curve245.POINTS,
        "https://elliptic-rank.icarm.cloud/curve/245",
    )


def exact_rank(rows: np.ndarray) -> int:
    """Return the rational row rank of a small integral matrix."""

    matrix = [[Fraction(int(value)) for value in row] for row in rows.tolist()]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                value - multiplier * pivot_entry
                for value, pivot_entry in zip(matrix[index], matrix[rank])
            ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def curve245_true_generic_subgroup(
    digits: int,
    forced_basis_rows: np.ndarray,
) -> tuple[dict[str, object], np.ndarray]:
    """Recover the known generic rank-12 subgroup in public coordinates."""

    parameter = CANONICAL_PARAMETER
    change = WeierstrassChange.from_values(ANCHOR_SHORT_TO_PUBLIC_CHANGE)
    if change_weierstrass_model(primitive_short_model(parameter), change) != PUBLIC_MODEL:
        raise RuntimeError("the pinned curve-245 short/public change ceased to agree")
    quartic_points = list(primitive_visible_points(CONSTRUCTION, parameter))
    quartic_points.append(extra_quartic_point(parameter))
    short_points = [
        quartic_point_to_short_jacobian(CONSTRUCTION, parameter, point)
        for point in quartic_points
    ]
    public_points = [source_point_to_target(point, change) for point in short_points]

    curve = curve245_data()
    basis_text = "[" + ",".join(gp_vector(point) for point in curve.points) + "]"
    generic_text = "[" + ",".join(gp_vector(point) for point in public_points) + "]"
    program = f"""
default(realprecision,{max(digits, 150)});
E=ellinit({gp_vector(curve.coefficients)});P={basis_text};Q={generic_text};A=concat(P,Q);
H=ellheightmatrix(E,A);HP=H[1..#P,1..#P];C=H[1..#P,#P+1..#A];X=HP^-1*C;Z=matrix(#P,#Q,i,j,round(X[i,j]));K=matkerint(Z);
ok=1;for(j=1,#Q,T=[0];for(i=1,#P,if(Z[i,j],T=elladd(E,T,ellmul(E,P[i],Z[i,j]))));if(ellsub(E,T,Q[j])!=[0],ok=0));
print("COORDINATE_RANK|",matrank(Z));print("EXACT_REPLAY|",ok);print1("RELATION|");for(i=1,matsize(K)[1],if(i>1,print1(","));print1(K[i,1]));print();
print("BEGIN_COORDINATES");for(j=1,#Q,for(i=1,#P,if(i>1,print1("|"));print1(Z[i,j]));print());print("END_COORDINATES");
"""
    lines = run_gp(program)
    coordinates = np.array(
        parse_matrix(lines, "BEGIN_COORDINATES", "END_COORDINATES", int),
        dtype=np.int64,
    )
    coordinate_rank = int(
        next(line.split("|", 1)[1] for line in lines if line.startswith("COORDINATE_RANK|"))
    )
    replay = next(line.split("|", 1)[1] for line in lines if line.startswith("EXACT_REPLAY|"))
    relation = next(line.split("|", 1)[1] for line in lines if line.startswith("RELATION|"))
    expected_relation = ",".join(map(str, [1] * 12 + [0]))
    if coordinate_rank != 12 or replay != "1" or relation != expected_relation:
        raise RuntimeError("the exact curve-245 generic subgroup certificate changed")

    # The twelve visible points have one relation; retain eleven of them and
    # Fermigier's extra section as a basis of the known generic subgroup.
    true_basis_rows = coordinates[[*range(11), 12]]
    forced_rank = exact_rank(forced_basis_rows)
    true_rank = exact_rank(true_basis_rows)
    union_rank = exact_rank(np.vstack((forced_basis_rows, true_basis_rows)))
    intersection_rank = forced_rank + true_rank - union_rank
    fingerprint = integrality_fingerprint(
        curve,
        true_basis_rows,
        bound=28.0,
        digits=digits,
    )
    return (
        {
            "construction": "Fermigier--Mestre six-root family at canonical T=5801/10",
            "specialized_quartic_points": 13,
            "coordinate_rank": coordinate_rank,
            "unique_relation_in_quartic_point_order": [1] * 12 + [0],
            "generic_basis_point_indices_one_based": list(range(1, 12)) + [13],
            "exact_group_law_replay": True,
            "transported_public_points": [
                {"x": fraction_text(point[0]), "y": fraction_text(point[1])}
                for point in public_points
            ],
            "coordinates_by_generic_point_in_public_basis": coordinates.tolist(),
            "true_generic_rank": true_rank,
            "forced_candidate_rank": forced_rank,
            "rank_of_sum": union_rank,
            "intersection_rank_with_forced_rank17_candidate": intersection_rank,
            "short_vector_integrality_fingerprint": fingerprint,
        },
        true_basis_rows,
    )


def curve245_control(digits: int) -> tuple[dict[str, object], np.ndarray]:
    curve = curve245_data()
    points = "[" + ",".join(gp_vector(point) for point in curve.points) + "]"
    program = f"""
default(realprecision,{digits});
E=ellinit({gp_vector(curve.coefficients)});P={points};H=ellheightmatrix(E,P);
B={gp_matrix(CURVE245_SEED_BASIS)};K=matkerint(B);S=matkerint(K~);
G=S~*H*S;U=qflllgram(G);R=U~*G*U;
Q=qfminim(R,35,100000,2);
V=vecsort(vector(matsize(Q[3])[2],i,Q[3][,i]~*R*Q[3][,i]));
print("AMBIENT_DETERMINANT|",matdet(H));print("CORE_DETERMINANT|",matdet(G));
print("SHELL_LINES|",#V);
print("SHELL|",V[1],"|",V[132],"|",V[328],"|",V[656],"|",V[984],"|",V[1180],"|",V[1311],"|",V[1312]);
print("BEGIN_SATURATED_BASIS");
for(i=1,matsize(S)[1],for(j=1,matsize(S)[2],if(j>1,print1("|"));print1(S[i,j]));print());
print("END_SATURATED_BASIS");
print("BEGIN_REDUCED_GRAM");
for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(R[i,j]));print());
print("END_REDUCED_GRAM");
"""
    lines = run_gp(program)
    scalar = {
        line.split("|", 1)[0]: line.split("|", 1)[1]
        for line in lines
        if line.startswith(("AMBIENT_DETERMINANT|", "CORE_DETERMINANT|", "SHELL_LINES|"))
    }
    shell = next(line for line in lines if line.startswith("SHELL|")).split("|")[1:]
    reduced_text = parse_matrix(lines, "BEGIN_REDUCED_GRAM", "END_REDUCED_GRAM", str)
    reduced = np.array(reduced_text, dtype=float)
    saturated = parse_matrix(lines, "BEGIN_SATURATED_BASIS", "END_SATURATED_BASIS", int)

    vectors, _heights = enumerate_vectors(curve, digits, 28.0)
    mask = exact_in_span_mask(vectors, np.array(CURVE245_SEED_BASIS, dtype=np.int64))
    core_det = float(scalar["CORE_DETERMINANT"])
    determinant_scale = math.exp((math.log(core_det) - math.log(948.0)) / 17.0)
    record = {
        "label": curve.label,
        "role": "known non-R17 parent-family control",
        "exact_parent_reconstruction": "Fermigier--Mestre generic rank 12",
        "exact_parent_module": "elliptic-curves/cas/icarm_curve245_mestre.py",
        "search_bound": 28,
        "ambient_short_vector_lines_at_bound": len(vectors),
        "candidate_rank17_lines_at_bound": int(mask.sum()),
        "seed_basis_rows_in_public_point_coordinates": [list(row) for row in CURVE245_SEED_BASIS],
        "saturated_basis_columns_in_public_point_coordinates": saturated,
        "ambient_determinant": scalar["AMBIENT_DETERMINANT"],
        "candidate_core_determinant": scalar["CORE_DETERMINANT"],
        "determinant_forced_lambda": f"{determinant_scale:.17g}",
        "shell_enumerated_lines_to_height_35": int(scalar["SHELL_LINES"]),
        "first_1311_height_quantiles": shell,
        "first_1311_quantiles_over_determinant_lambda": [
            f"{float(value) / determinant_scale:.17g}" for value in shell[:7]
        ],
        "core_lll_gram": reduced_text,
    }
    return record, reduced


def integrality_fingerprint(
    curve: CurveData,
    candidate_basis_rows: np.ndarray,
    *,
    bound: float,
    digits: int,
) -> dict[str, object]:
    """Measure integral-point enrichment without using it in core selection."""

    points = "[" + ",".join(gp_vector(point) for point in curve.points) + "]"
    program = f"""
default(parisizemax,2000000000);
default(parisize,500000000);
default(realprecision,{digits});
E=ellinit({gp_vector(curve.coefficients)});P={points};H=ellheightmatrix(E,P);
U=qflllgram(H);R=U~*H*U;Q=qfminim(R,{bound},100000,2);V=U*Q[3];
for(j=1,matsize(V)[2],T=[0];for(i=1,#P,if(V[i,j],T=elladd(E,T,ellmul(E,P[i],V[i,j]))));for(i=1,matsize(V)[1],if(i>1,print1(","));print1(V[i,j]));print("|",denominator(T[1])));
"""
    lines = run_gp(program)
    vectors = []
    denominator_by_line = {}
    for line in lines:
        coordinates, denominator = line.rsplit("|", 1)
        vector = np.array(list(map(int, coordinates.split(","))), dtype=np.int64)
        vectors.append(vector)
        denominator_by_line[canonical_line(vector)] = int(denominator)
    vectors_array = np.array(vectors, dtype=np.int64)
    in_core = exact_in_span_mask(vectors_array, candidate_basis_rows)
    integral = np.array(
        [denominator_by_line[canonical_line(vector)] == 1 for vector in vectors_array]
    )
    core_integral = int(np.sum(integral & in_core))
    outside_integral = int(np.sum(integral & ~in_core))
    core_total = int(np.sum(in_core))
    outside_total = int(np.sum(~in_core))
    core_rate = core_integral / core_total
    outside_rate = outside_integral / outside_total
    odds_ratio = (
        core_integral
        * (outside_total - outside_integral)
        / ((core_total - core_integral) * outside_integral)
    )
    return {
        "label": curve.label,
        "height_bound": bound,
        "selection_used_integrality": False,
        "ambient_lines": len(vectors),
        "core_lines": core_total,
        "outside_lines": outside_total,
        "core_integral_lines": core_integral,
        "outside_integral_lines": outside_integral,
        "core_integral_rate": f"{core_rate:.17g}",
        "outside_integral_rate": f"{outside_rate:.17g}",
        "core_over_outside_rate_ratio": f"{core_rate / outside_rate:.17g}",
        "core_over_outside_odds_ratio": f"{odds_ratio:.17g}",
    }


def integral_inverse(matrix: np.ndarray) -> np.ndarray:
    inverse = np.rint(np.linalg.inv(matrix)).astype(np.int64)
    identity = np.eye(len(matrix), dtype=np.int64)
    if not np.array_equal(inverse @ matrix, identity):
        raise RuntimeError("basis ceased to be unimodular")
    return inverse


def fit_score(height_gram: np.ndarray, exact_form: np.ndarray) -> tuple[float, float, float]:
    scale = float(np.sum(height_gram * exact_form) / np.sum(exact_form * exact_form))
    error = height_gram / scale - exact_form
    return scale, float(np.sqrt(np.mean(error * error))), float(np.max(np.abs(error)))


def mapped_shell_dispersion(
    height_gram: np.ndarray,
    fitted_basis: np.ndarray,
    unoriented_r17_shell: np.ndarray,
) -> dict[str, object]:
    """Test all 1,311 exact minimal lines, not just the fitted Gram entries."""

    coordinates = unoriented_r17_shell @ integral_inverse(fitted_basis).T
    heights = np.einsum("ij,jk,ik->i", coordinates, height_gram, coordinates)
    mean = float(np.mean(heights))
    normalized = 4.0 * heights / mean
    quantiles = np.quantile(normalized, [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1])
    return {
        "unoriented_lines": len(unoriented_r17_shell),
        "normalization": "mean mapped height = exact R17 minimum 4",
        "raw_mean_height": f"{mean:.17g}",
        "normalized_quantiles_min_10_25_50_75_90_max": [
            f"{value:.17g}" for value in quantiles
        ],
        "normalized_rms_deviation_from_4": f"{float(np.sqrt(np.mean((normalized - 4.0) ** 2))):.17g}",
        "coefficient_of_variation": f"{float(np.std(heights) / mean):.17g}",
        "fraction_with_normalized_height_between_3_and_5": f"{float(np.mean((normalized >= 3.0) & (normalized <= 5.0))):.17g}",
    }


def fit_r17(
    height_gram: np.ndarray,
    r17_form: np.ndarray,
    signed_shell: np.ndarray,
    *,
    restarts: int,
    seed: int,
) -> dict[str, object]:
    rng = random.Random(seed)
    permutation_rng = np.random.default_rng(seed)

    def valid_replacements(basis: np.ndarray, column: int) -> np.ndarray:
        coordinates = signed_shell @ integral_inverse(basis).T
        return np.flatnonzero(np.abs(coordinates[:, column]) == 1)

    def objective(basis: np.ndarray) -> tuple[float, float, float]:
        exact_form = basis.T @ r17_form @ basis
        scale, rms, maximum = fit_score(height_gram, exact_form)
        return rms * rms, scale, maximum

    def descend(basis: np.ndarray) -> tuple[np.ndarray, tuple[float, float, float]]:
        exact_form = basis.T @ r17_form @ basis
        old = objective(basis)
        for _sweep in range(50):
            changed = False
            for column in permutation_rng.permutation(17):
                indices = valid_replacements(basis, int(column))
                pairings = signed_shell[indices] @ r17_form @ basis
                scale = fit_score(height_gram, exact_form)[0]
                screen = np.mean((pairings - height_gram[column] / scale) ** 2, axis=1)
                best = old
                best_vector = None
                for offset in np.argsort(screen, kind="stable")[:20]:
                    trial = basis.copy()
                    trial[:, column] = signed_shell[indices[offset]]
                    candidate = objective(trial)
                    if candidate[0] < best[0] - 1e-14:
                        best = candidate
                        best_vector = trial[:, column].copy()
                if best_vector is not None:
                    basis[:, column] = best_vector
                    exact_form = basis.T @ r17_form @ basis
                    old = best
                    changed = True
            if not changed:
                break
        return basis, objective(basis)

    best_basis = None
    best = (float("inf"), 0.0, 0.0)
    for _restart in range(restarts):
        basis = np.eye(17, dtype=np.int64)
        for _step in range(30 + rng.randrange(100)):
            column = rng.randrange(17)
            indices = valid_replacements(basis, column)
            basis[:, column] = signed_shell[rng.choice(indices)]
        basis, score = descend(basis)
        if score[0] < best[0]:
            best_basis, best = basis.copy(), score

    assert best_basis is not None
    exact_form = best_basis.T @ r17_form @ best_basis
    if round(np.linalg.det(best_basis)) not in (-1, 1):
        raise RuntimeError("reported R17 fit is not unimodular")
    return {
        "fit_scale": f"{best[1]:.17g}",
        "normalized_entrywise_rms": f"{math.sqrt(best[0]):.17g}",
        "normalized_maximum_absolute_entry_error": f"{best[2]:.17g}",
        "unimodular_r17_shell_basis_columns": best_basis.tolist(),
        "exact_isometric_r17_gram_in_fitted_basis": exact_form.tolist(),
        "mapped_full_minimal_shell": mapped_shell_dispersion(
            height_gram,
            best_basis,
            signed_shell[: len(signed_shell) // 2],
        ),
    }


def numerical_short_vectors(form: np.ndarray, bound: float) -> np.ndarray:
    matrix = "[" + ";".join(",".join(f"{value:.17g}" for value in row) for row in form) + "]"
    program = f"""
default(realprecision,50);A={matrix};Q=qfminim(A,{bound},200000,2);V=Q[3];
for(j=1,matsize(V)[2],for(i=1,17,if(i>1,print1("|"));print1(V[i,j]));print());
"""
    unoriented = np.array(
        [[int(value) for value in line.split("|")] for line in run_gp(program)],
        dtype=np.int64,
    )
    return np.vstack((unoriented, -unoriented))


def fit_numerical_core_pair(
    source: np.ndarray,
    target: np.ndarray,
    signed_source_vectors: np.ndarray,
    *,
    restarts: int,
    seed: int,
) -> dict[str, object]:
    """Fit two numerical cores up to GL(17,Z) and one scalar."""

    rng = random.Random(seed)
    permutation_rng = np.random.default_rng(seed)

    def valid_replacements(basis: np.ndarray, column: int) -> np.ndarray:
        coordinates = signed_source_vectors @ integral_inverse(basis).T
        return np.flatnonzero(np.abs(coordinates[:, column]) == 1)

    def objective(basis: np.ndarray) -> tuple[float, float, float]:
        transformed = basis.T @ source @ basis
        scale = float(np.sum(target * transformed) / np.sum(transformed * transformed))
        error = target / scale - transformed
        denominator = float(np.sqrt(np.mean(transformed * transformed)))
        relative_rms = float(np.sqrt(np.mean(error * error)) / denominator)
        relative_maximum = float(np.max(np.abs(error)) / denominator)
        return relative_rms * relative_rms, scale, relative_maximum

    def descend(basis: np.ndarray) -> tuple[np.ndarray, tuple[float, float, float]]:
        old = objective(basis)
        for _sweep in range(40):
            changed = False
            for column in permutation_rng.permutation(17):
                indices = valid_replacements(basis, int(column))
                pairings = signed_source_vectors[indices] @ source @ basis
                scale = old[1]
                screen = np.mean((pairings - target[column] / scale) ** 2, axis=1)
                best = old
                best_vector = None
                for offset in np.argsort(screen, kind="stable")[:20]:
                    trial = basis.copy()
                    trial[:, column] = signed_source_vectors[indices[offset]]
                    candidate = objective(trial)
                    if candidate[0] < best[0] - 1e-14:
                        best = candidate
                        best_vector = trial[:, column].copy()
                if best_vector is not None:
                    basis[:, column] = best_vector
                    old = best
                    changed = True
            if not changed:
                break
        return basis, objective(basis)

    best_basis = None
    best = (float("inf"), 0.0, 0.0)
    for _restart in range(restarts):
        basis = np.eye(17, dtype=np.int64)
        for _step in range(20 + rng.randrange(50)):
            column = rng.randrange(17)
            indices = valid_replacements(basis, column)
            basis[:, column] = signed_source_vectors[rng.choice(indices)]
        basis, score = descend(basis)
        if score[0] < best[0]:
            best_basis, best = basis.copy(), score

    assert best_basis is not None
    return {
        "fit_scale": f"{best[1]:.17g}",
        "relative_entrywise_rms": f"{math.sqrt(best[0]):.17g}",
        "relative_maximum_absolute_entry_error": f"{best[2]:.17g}",
        "unimodular_source_basis_columns": best_basis.tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--digits", type=int, default=100)
    parser.add_argument("--restarts", type=int, default=64)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    source = json.loads(args.input.read_text())
    r17_form, signed_shell = exact_r17_reduced_form_and_shell()
    control, control_gram = curve245_control(args.digits)
    true_generic, _true_basis = curve245_true_generic_subgroup(
        args.digits,
        np.array(CURVE245_SEED_BASIS, dtype=np.int64),
    )
    control["exact_true_generic_subgroup"] = true_generic

    grams = {
        record["label"]: np.array(record["core_lll_gram"], dtype=float)
        for record in source["curves"]
    }
    grams[control["label"]] = control_gram
    fits = []
    for index, (label, gram) in enumerate(grams.items()):
        fits.append(
            {
                "label": label,
                **fit_r17(
                    gram,
                    r17_form,
                    signed_shell,
                    restarts=args.restarts,
                    seed=948_000 + 101 * index,
                ),
            }
        )

    record_rms = [
        float(record["normalized_entrywise_rms"])
        for record in fits
        if record["label"] != control["label"]
    ]
    control_rms = float(fits[-1]["normalized_entrywise_rms"])
    ratio = control_rms / (sum(record_rms) / len(record_rms))
    source_by_label = {record["label"]: record for record in source["curves"]}
    search_bounds = {"rank28": 60.0, "rank29": 60.0, "curve273": 65.0, "curve302": 70.0}
    integrality = []
    for curve in CURVES:
        saturated_columns = np.array(
            source_by_label[curve.label]["saturated_basis_columns_in_public_point_coordinates"],
            dtype=np.int64,
        )
        integrality.append(
            integrality_fingerprint(
                curve,
                saturated_columns.T,
                bound=search_bounds[curve.label],
                digits=args.digits,
            )
        )
    integrality.append(
        integrality_fingerprint(
            curve245_data(),
            np.array(CURVE245_SEED_BASIS, dtype=np.int64),
            bound=28.0,
            digits=args.digits,
        )
    )
    pair_shell_bounds = {
        "rank28": 75.0,
        "rank29": 75.0,
        "curve273": 75.0,
        "curve302": 82.0,
        control["label"]: 35.0,
    }
    pair_shells = {
        label: numerical_short_vectors(gram, pair_shell_bounds[label])
        for label, gram in grams.items()
    }
    pair_labels = (
        ("curve273", "curve302"),
        ("rank29", "curve302"),
        ("rank28", "curve302"),
        (control["label"], "curve302"),
        ("rank29", "curve273"),
        (control["label"], "curve273"),
    )
    pairwise_fits = []
    for index, (source_label, target_label) in enumerate(pair_labels):
        pairwise_fits.append(
            {
                "source": source_label,
                "target": target_label,
                "source_unoriented_short_vectors": len(pair_shells[source_label]) // 2,
                **fit_numerical_core_pair(
                    grams[source_label],
                    grams[target_label],
                    pair_shells[source_label],
                    restarts=args.restarts,
                    seed=17_302 + 103 * index,
                ),
            }
        )
    pair_by_labels = {
        (record["source"], record["target"]): record for record in pairwise_fits
    }
    control_to_302 = float(
        pair_by_labels[(control["label"], "curve302")]["relative_entrywise_rms"]
    )
    record_to_302 = [
        float(pair_by_labels[(source_label, "curve302")]["relative_entrywise_rms"])
        for source_label in ("rank28", "rank29", "curve273")
    ]
    pairwise_control_ratio = control_to_302 / (sum(record_to_302) / len(record_to_302))
    payload = {
        "schema": "record-rank17-fingerprint-calibration-v1",
        "status": "bounded numerical negative-control experiment",
        "decimal_precision_digits": args.digits,
        "fit_random_restarts": args.restarts,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "target_lattice": {
            "rank": 17,
            "determinant": 948,
            "minimal_norm": 4,
            "unoriented_minimal_shell_size": 1311,
        },
        "method": (
            "Coordinate descent over unimodular bases made entirely from exact R17 "
            "minimal vectors; each numerical Gram is fitted by one scalar."
        ),
        "warning": (
            "The fit is deliberately optimized after a rank-17 subspace was selected. "
            "It measures basis-entry compatibility, not statistical significance or "
            "specialization; the separately reported mapped-shell diagnostic tests all "
            "1,311 exact R17 minimal lines."
        ),
        "negative_control": control,
        "exact_r17_fits": fits,
        "direct_pairwise_numerical_core_fits": {
            "method": (
                "Fit each target Gram by a scalar times a GL(17,Z) transform of a source "
                "Gram, using the source's declared bounded short-vector cloud."
            ),
            "source_shell_bounds": pair_shell_bounds,
            "fits": pairwise_fits,
            "curve245_to_302_rms_over_record_source_mean_rms": f"{pairwise_control_ratio:.17g}",
            "discrimination_result": (
                "FAILS_NEGATIVE_CONTROL"
                if pairwise_control_ratio < 1.20
                else "SEPARATES_CONTROL_AT_DECLARED_THRESHOLD"
            ),
        },
        "out_of_sample_integrality_fingerprint": integrality,
        "control_rms_over_record_mean_rms": f"{ratio:.17g}",
        "discrimination_result": "FAILS_NEGATIVE_CONTROL" if ratio < 1.20 else "SEPARATES_CONTROL_AT_DECLARED_THRESHOLD",
        "interpretation": (
            "At the declared 1.20 ratio threshold, the forced R17 basis-entry fit does "
            "not distinguish the record candidates from a known Fermigier--Mestre control, "
            "and its transported 1,311-vector shells are broadly dispersed rather than "
            "approximately minimal. The forced rank-17 control intersects the exactly "
            "recovered generic rank-12 subgroup in only dimension 9, the dimension-count "
            "minimum. "
            "Direct pairwise numerical-core fitting likewise fails its negative control. "
            "The separately evaluated integrality enrichment does distinguish the selected "
            "record cores, especially curves 273 and 302, as structured subspaces, but it "
            "does not identify their integral isometry class."
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != rendered:
            raise SystemExit(f"FAIL: {args.output} differs from recomputation")
        print(f"PASS|{args.output}|sha256={hashlib.sha256(rendered.encode()).hexdigest()}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(f"WROTE|{args.output}|sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
