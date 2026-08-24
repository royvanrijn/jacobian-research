#!/usr/bin/env python3
"""Compute high-precision height lattices for the 28/29/273/302 data set.

This is a numerical Archimedean comparison, not an independence certificate
or a K3-specialization certificate.  PARI/GP computes canonical heights; the
script records the public-point Gram, a unimodular LLL transform, and the
reduced Gram without rounding them to binary64.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
import platform
import subprocess
from typing import Iterable

import elkies_klagsbrun_rank29
import elkies_rank28
import icarm_curve273
import icarm_curve302


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_height_lattices_28_29_273_302_v1.json"
)


@dataclass(frozen=True)
class CurveData:
    label: str
    coefficients: tuple[Fraction, ...]
    points: tuple[tuple[Fraction, Fraction], ...]
    source: str


CURVES = (
    CurveData(
        "rank28",
        elkies_rank28.GENERAL_WEIERSTRASS_COEFFICIENTS,
        elkies_rank28.POINTS,
        "https://web.math.pmf.unizg.hr/~duje/tors/rk28.html",
    ),
    CurveData(
        "rank29",
        elkies_klagsbrun_rank29.GENERAL_WEIERSTRASS_COEFFICIENTS,
        elkies_klagsbrun_rank29.PUBLISHED_POINTS,
        "https://web.math.pmf.unizg.hr/~duje/tors/rk29.html",
    ),
    CurveData(
        "curve273",
        icarm_curve273.GENERAL_WEIERSTRASS_COEFFICIENTS,
        icarm_curve273.POINTS,
        "https://elliptic-rank.icarm.cloud/curve/273",
    ),
    CurveData(
        "curve302",
        icarm_curve302.GENERAL_WEIERSTRASS_COEFFICIENTS,
        icarm_curve302.POINTS,
        "https://elliptic-rank.icarm.cloud/curve/302",
    ),
)


def gp_rational(value: Fraction) -> str:
    value = Fraction(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"({value.numerator}/{value.denominator})"


def gp_vector(values: Iterable[Fraction]) -> str:
    return "[" + ",".join(gp_rational(value) for value in values) + "]"


def gp_program(curve: CurveData, digits: int) -> str:
    points = "[" + ",".join(gp_vector(point) for point in curve.points) + "]"
    return f"""
default(realprecision,{digits});
E=ellinit({gp_vector(curve.coefficients)});
P={points};
H=ellheightmatrix(E,P);
U=qflllgram(H);
R=U~*H*U;
Q=qfminim(R,,,2);
print("PARI_VERSION|",version());
print("DETERMINANT|",matdet(H));
print("MINIMUM_COUNT_SIGNED|",Q[1]);
print("MINIMUM|",Q[2]);
print("BEGIN_HEIGHT_GRAM");
for(i=1,matsize(H)[1],for(j=1,matsize(H)[2],if(j>1,print1("|"));print1(H[i,j]));print());
print("END_HEIGHT_GRAM");
print("BEGIN_LLL_TRANSFORM");
for(i=1,matsize(U)[1],for(j=1,matsize(U)[2],if(j>1,print1("|"));print1(U[i,j]));print());
print("END_LLL_TRANSFORM");
print("BEGIN_REDUCED_GRAM");
for(i=1,matsize(R)[1],for(j=1,matsize(R)[2],if(j>1,print1("|"));print1(R[i,j]));print());
print("END_REDUCED_GRAM");
"""


def parse_matrix(lines: list[str], begin: str, end: str) -> list[list[str]]:
    i = lines.index(begin) + 1
    j = lines.index(end, i)
    return [line.split("|") for line in lines[i:j]]


def run_curve(curve: CurveData, digits: int) -> dict[str, object]:
    completed = subprocess.run(
        ["gp", "-q"],
        input=gp_program(curve, digits),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    if completed.stderr.strip():
        raise RuntimeError(f"gp stderr for {curve.label}: {completed.stderr}")
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]

    scalar: dict[str, str] = {}
    for line in lines:
        for key in (
            "PARI_VERSION",
            "DETERMINANT",
            "MINIMUM_COUNT_SIGNED",
            "MINIMUM",
        ):
            prefix = key + "|"
            if line.startswith(prefix):
                scalar[key.lower()] = line[len(prefix) :]

    height_gram = parse_matrix(lines, "BEGIN_HEIGHT_GRAM", "END_HEIGHT_GRAM")
    lll_transform_text = parse_matrix(
        lines, "BEGIN_LLL_TRANSFORM", "END_LLL_TRANSFORM"
    )
    reduced_gram = parse_matrix(lines, "BEGIN_REDUCED_GRAM", "END_REDUCED_GRAM")
    lll_transform = [[int(value) for value in row] for row in lll_transform_text]

    n = len(curve.points)
    if not all(len(matrix) == n for matrix in (height_gram, lll_transform, reduced_gram)):
        raise RuntimeError(f"wrong matrix row count for {curve.label}")
    for matrix in (height_gram, lll_transform, reduced_gram):
        if not all(len(row) == n for row in matrix):
            raise RuntimeError(f"nonsquare matrix for {curve.label}")

    return {
        "label": curve.label,
        "rank_lower_bound": n,
        "public_source": curve.source,
        "weierstrass_coefficients": [gp_rational(value) for value in curve.coefficients],
        "point_count": n,
        **scalar,
        "height_gram": height_gram,
        "lll_transform_columns": lll_transform,
        "lll_reduced_gram": reduced_gram,
        "lll_reduced_diagonal": [reduced_gram[i][i] for i in range(n)],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--digits", type=int, default=100)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.digits < 50:
        raise SystemExit("--digits must be at least 50")

    records = [run_curve(curve, args.digits) for curve in CURVES]
    payload = {
        "schema": "record-height-lattices-28-29-273-302-v1",
        "status": "high-precision numerical computation",
        "interpretation_warning": (
            "These are Gram matrices of the displayed public-point subgroups. "
            "They do not identify a generic K3 subgroup or prove a family specialization."
        ),
        "decimal_precision_digits": args.digits,
        "python_version": platform.python_version(),
        "pari_version": records[0]["pari_version"],
        "curves": records,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"

    if args.check:
        current = args.output.read_text()
        if current != rendered:
            raise SystemExit(f"FAIL: {args.output} differs from recomputation")
        print(f"PASS|{args.output}|sha256={hashlib.sha256(current.encode()).hexdigest()}")
        return

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(f"WROTE|{args.output}|sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
