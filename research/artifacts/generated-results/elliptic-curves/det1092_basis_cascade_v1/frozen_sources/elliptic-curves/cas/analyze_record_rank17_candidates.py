#!/usr/bin/env python3
"""Analyze the additive-RANSAC rank-17 candidate spaces in four record curves.

The inputs are local search outputs produced by ``search_record_rank17_core.py``.
Each rational 17-space is saturated exactly inside the public-point lattice
with PARI's integer-kernel routine before its canonical-height Gram is formed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from compare_record_height_lattices import CURVES, gp_vector
from search_record_rank17_core import enumerate_vectors, exact_in_span_mask


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_rank17_core_candidates_v1.json"
)
SEARCH_PARAMETERS = {
    "rank28": (60.0, 2423),
    "rank29": (60.0, 1507),
    "curve273": (65.0, 5936),
    "curve302": (70.0, 3458),
}
SHELL_ENUMERATION_BOUNDS = {
    "rank28": 75,
    "rank29": 75,
    "curve273": 75,
    "curve302": 82,
}


def gp_matrix(rows: list[list[int]]) -> str:
    return "[" + ";".join(",".join(str(x) for x in row) for row in rows) + "]"


def parse_matrix(lines: list[str], begin: str, end: str) -> list[list[str]]:
    i = lines.index(begin) + 1
    j = lines.index(end, i)
    return [line.split("|") for line in lines[i:j]]


def read_seed_basis(label: str) -> list[list[int]]:
    path = ROOT / f"artifacts/local/elliptic-curves/{label}-r17-additive-ransac.txt"
    lines = path.read_text().splitlines()
    start = lines.index("basis_rows_in_public_point_coordinates:") + 1
    return [[int(x) for x in line.split()] for line in lines[start : start + 17]]


def analyze_curve(curve, digits: int) -> dict[str, object]:
    seed_basis = read_seed_basis(curve.label)
    points = "[" + ",".join(gp_vector(point) for point in curve.points) + "]"
    shell_bound = SHELL_ENUMERATION_BOUNDS[curve.label]
    program = f"""
default(realprecision,{digits});
E=ellinit({gp_vector(curve.coefficients)});
P={points};
H=ellheightmatrix(E,P);
B={gp_matrix(seed_basis)};
K=matkerint(B);
S=matkerint(K~);
G=S~*H*S;
U=qflllgram(G);
R=U~*G*U;
Q=qfminim(R,{shell_bound},100000,2);
V=vecsort(vector(matsize(Q[3])[2],i,Q[3][,i]~*R*Q[3][,i]));
print("AMBIENT_DETERMINANT|",matdet(H));
print("CORE_DETERMINANT|",matdet(G));
print("SHELL_ENUMERATED_LINES|",#V);
print("SHELL_HEIGHTS|",V[1],"|",V[132],"|",V[328],"|",V[656],"|",V[984],"|",V[1180],"|",V[1311],"|",V[1312]);
print("BEGIN_SATURATED_BASIS");
for(i=1,matsize(S)[1],for(j=1,matsize(S)[2],if(j>1,print1("|"));print1(S[i,j]));print());
print("END_SATURATED_BASIS");
print("BEGIN_CORE_GRAM");
for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(G[i,j]));print());
print("END_CORE_GRAM");
print("BEGIN_CORE_LLL_TRANSFORM");
for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(U[i,j]));print());
print("END_CORE_LLL_TRANSFORM");
print("BEGIN_CORE_LLL_GRAM");
for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(R[i,j]));print());
print("END_CORE_LLL_GRAM");
"""
    completed = subprocess.run(
        ["gp", "-q"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    scalars = {}
    for line in lines:
        for key in (
            "AMBIENT_DETERMINANT",
            "CORE_DETERMINANT",
            "SHELL_ENUMERATED_LINES",
        ):
            if line.startswith(key + "|"):
                scalars[key.lower()] = line.split("|", 1)[1]
    saturated = [
        [int(x) for x in row]
        for row in parse_matrix(lines, "BEGIN_SATURATED_BASIS", "END_SATURATED_BASIS")
    ]
    gram = parse_matrix(lines, "BEGIN_CORE_GRAM", "END_CORE_GRAM")
    lll_transform = [
        [int(x) for x in row]
        for row in parse_matrix(
            lines, "BEGIN_CORE_LLL_TRANSFORM", "END_CORE_LLL_TRANSFORM"
        )
    ]
    reduced = parse_matrix(lines, "BEGIN_CORE_LLL_GRAM", "END_CORE_LLL_GRAM")
    shell_line = next(line for line in lines if line.startswith("SHELL_HEIGHTS|"))
    shell_heights = shell_line.split("|")[1:]

    ambient_det = float(scalars["ambient_determinant"])
    core_det = float(scalars["core_determinant"])
    determinant_scale = math.exp((math.log(core_det) - math.log(948.0)) / 17.0)
    diagonal = np.array([float(reduced[i][i]) for i in range(17)])

    bound, expected_lines = SEARCH_PARAMETERS[curve.label]
    vectors, heights = enumerate_vectors(curve, digits, bound)
    seed = np.array(seed_basis, dtype=np.int64)
    exact_mask = exact_in_span_mask(vectors, seed)
    if len(vectors) != expected_lines:
        raise RuntimeError(f"short-vector census changed for {curve.label}")

    return {
        "label": curve.label,
        "ambient_rank_lower_bound": len(curve.points),
        "exceptional_quotient_dimension": len(curve.points) - 17,
        "search_bound": bound,
        "ambient_short_vector_lines_at_bound": len(vectors),
        "candidate_core_lines_at_bound": int(exact_mask.sum()),
        "candidate_core_height_range_at_bound": [
            f"{heights[exact_mask].min():.17g}",
            f"{heights[exact_mask].max():.17g}",
        ],
        "approximate_1311_shell_enumeration_bound": shell_bound,
        "approximate_1311_shell_enumerated_lines": int(
            scalars.pop("shell_enumerated_lines")
        ),
        "approximate_1311_shell_height_quantiles": {
            "minimum": shell_heights[0],
            "p10_nearest_rank": shell_heights[1],
            "p25_nearest_rank": shell_heights[2],
            "median_nearest_rank": shell_heights[3],
            "p75_nearest_rank": shell_heights[4],
            "p90_nearest_rank": shell_heights[5],
            "pair_1311": shell_heights[6],
            "pair_1312": shell_heights[7],
        },
        "approximate_1311_shell_quantiles_over_lambda": [
            f"{float(value) / determinant_scale:.17g}" for value in shell_heights[:7]
        ],
        "seed_basis_rows_in_public_point_coordinates": seed_basis,
        "saturated_basis_columns_in_public_point_coordinates": saturated,
        "saturation_method": "S=matkerint(matkerint(B)^t)",
        **scalars,
        "quotient_schur_determinant": f"{ambient_det / core_det:.17g}",
        "determinant_forced_lambda_from_det_core_eq_948_lambda_pow_17": f"{determinant_scale:.17g}",
        "four_lambda": f"{4 * determinant_scale:.17g}",
        "core_gram": gram,
        "core_lll_transform_columns": lll_transform,
        "core_lll_gram": reduced,
        "core_lll_diagonal": [reduced[i][i] for i in range(17)],
        "core_lll_diagonal_over_lambda": [
            f"{value / determinant_scale:.17g}" for value in diagonal
        ],
        "core_lll_diagonal_log_rms_from_four_lambda": f"{np.sqrt(np.mean(np.log(diagonal / (4 * determinant_scale)) ** 2)):.17g}",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--digits", type=int, default=100)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = {
        "schema": "record-rank17-core-candidates-v1",
        "status": "bounded numerical provenance evidence",
        "target_lattice": {
            "rank": 17,
            "determinant": 948,
            "minimal_norm": 4,
            "unoriented_minimal_shell_size": 1311,
        },
        "warning": (
            "The dimension-17 spaces were selected from a bounded short-vector cloud. "
            "Even exact saturation and a close scaled determinant do not prove that they are transported generic sections."
        ),
        "decimal_precision_digits": args.digits,
        "curves": [analyze_curve(curve, args.digits) for curve in CURVES],
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
