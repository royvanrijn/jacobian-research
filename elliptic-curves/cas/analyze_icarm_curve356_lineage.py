#!/usr/bin/env python3
"""Replay the bounded curve-351/356 common-section fingerprint.

This is numerical and source-comparative evidence, not a family-recognition
theorem.  The two ICARM JSON responses are hash-pinned as retrieved on
2026-08-25, and PARI/GP computes the canonical-height matrices of the first
seventeen displayed points in their displayed order.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from math import isqrt, sqrt
from pathlib import Path
import subprocess
import sys
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))

from icarm_curve356 import GENERAL_WEIERSTRASS_COEFFICIENTS, POINTS  # noqa: E402
from analyze_icarm_construction_fingerprints import (  # noqa: E402
    load_census_roots,
    recognize_families,
    weierstrass_invariants,
)


SOURCES = {
    351: (
        "https://elliptic-rank.icarm.cloud/curve/351.json",
        "02c0de1801d0c925dd6e42204f8461e99595926e95221e91da0c09466a6f67fd",
    ),
    356: (
        "https://elliptic-rank.icarm.cloud/curve/356.json",
        "58afbc62dbb6e01b47266c90edcf0e09bb003bb6a558333422b332e42546e89e",
    ),
}

EXPECTED_MATCHES = (
    (2, 1),
    (4, 1),
    (5, 71),
    (11, 5),
    (13, 679),
    (15, 1),
    (16, 7),
    (17, 41),
)
EXPECTED_SCALE = 1.4208782482875444
EXPECTED_RELATIVE_RESIDUAL = 0.11220111822209558
EXPECTED_CORRELATION = 0.9748839794656153


def rational_text(value: str) -> str:
    value_q = Fraction(value)
    if value_q.denominator == 1:
        return str(value_q.numerator)
    return f"({value_q.numerator}/{value_q.denominator})"


def fetch(curve_id: int) -> dict[str, object]:
    url, expected_hash = SOURCES[curve_id]
    with urlopen(url, timeout=30) as response:
        raw = response.read()
    observed_hash = hashlib.sha256(raw).hexdigest()
    if observed_hash != expected_hash:
        raise AssertionError(
            f"curve {curve_id} public JSON changed: {observed_hash} != {expected_hash}"
        )
    record = json.loads(raw)
    if int(record["id"]) != curve_id:
        raise AssertionError(f"wrong curve returned for {url}")
    return record


def denominator_roots(record: dict[str, object], count: int = 17) -> tuple[int, ...]:
    roots = []
    for x_value, _y_value in record["points"][:count]:
        denominator = Fraction(x_value).denominator
        root = isqrt(denominator)
        if root * root != denominator:
            raise AssertionError("an elliptic x-coordinate denominator is not a square")
        roots.append(root)
    return tuple(roots)


def height_matrix(record: dict[str, object], count: int = 17) -> tuple[list[float], str]:
    coefficients = ",".join(str(value) for value in record["ainvs"])
    points = ",".join(
        f"[{rational_text(x_value)},{rational_text(y_value)}]"
        for x_value, y_value in record["points"][:count]
    )
    program = f"""
default(realprecision,80);
E=ellinit([{coefficients}]);
P=[{points}];
H=ellheightmatrix(E,P);
print("PARI|",version());
print("MATRIX");
for(i=1,{count},for(j=1,{count},if(j>1,print1("|"));print1(H[i,j]));print());
"""
    completed = subprocess.run(
        ["gp", "-q"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=30,
        check=True,
    )
    if completed.stderr.strip():
        raise RuntimeError(f"PARI/GP stderr: {completed.stderr.strip()}")
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]
    pari_version = lines[0].removeprefix("PARI|")
    if lines[1] != "MATRIX":
        raise AssertionError("unexpected PARI/GP output")
    rows = [[float(value) for value in line.split("|")] for line in lines[2:]]
    if len(rows) != count or any(len(row) != count for row in rows):
        raise AssertionError("PARI/GP returned a nonsquare height matrix")
    return [value for row in rows for value in row], pari_version


def dot(left: list[float], right: list[float]) -> float:
    return sum(x_value * y_value for x_value, y_value in zip(left, right))


def norm(values: list[float]) -> float:
    return sqrt(dot(values, values))


def correlation(left: list[float], right: list[float]) -> float:
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    left_centered = [value - left_mean for value in left]
    right_centered = [value - right_mean for value in right]
    return dot(left_centered, right_centered) / (
        norm(left_centered) * norm(right_centered)
    )


def main() -> None:
    records = {curve_id: fetch(curve_id) for curve_id in SOURCES}
    curve356 = records[356]
    if tuple(Fraction(value) for value in curve356["ainvs"]) != GENERAL_WEIERSTRASS_COEFFICIENTS:
        raise AssertionError("the hard-coded curve-356 model differs from the public source")
    if tuple(
        (Fraction(x_value), Fraction(y_value)) for x_value, y_value in curve356["points"]
    ) != POINTS:
        raise AssertionError("the hard-coded curve-356 points differ from the public source")
    if (curve356["rank_lower_bound"], len(curve356["points"])) != (29, 29):
        raise AssertionError("curve 356 no longer has a 29-point witness")
    if (records[351]["rank_lower_bound"], len(records[351]["points"])) != (25, 25):
        raise AssertionError("curve 351 no longer has a 25-point witness")

    roots351 = denominator_roots(records[351])
    roots356 = denominator_roots(records[356])
    matches = tuple(
        (index, left)
        for index, (left, right) in enumerate(zip(roots351, roots356), 1)
        if left == right
    )
    if matches != EXPECTED_MATCHES:
        raise AssertionError(f"ordered denominator fingerprint changed: {matches}")

    height351, pari351 = height_matrix(records[351])
    height356, pari356 = height_matrix(records[356])
    if pari351 != pari356:
        raise AssertionError("PARI/GP version changed during the replay")
    scale = dot(height351, height356) / dot(height351, height351)
    residual = [right - scale * left for left, right in zip(height351, height356)]
    relative_residual = norm(residual) / norm(height356)
    pearson = correlation(height351, height356)
    for label, observed, expected in (
        ("scale", scale, EXPECTED_SCALE),
        ("relative residual", relative_residual, EXPECTED_RELATIVE_RESIDUAL),
        ("correlation", pearson, EXPECTED_CORRELATION),
    ):
        if abs(observed - expected) > 1e-12:
            raise AssertionError(f"{label} changed: {observed} != {expected}")

    c4, c6, discriminant, j_value = weierstrass_invariants(
        GENERAL_WEIERSTRASS_COEFFICIENTS
    )
    family_census = recognize_families(
        {
            356: {
                "ainvs": GENERAL_WEIERSTRASS_COEFFICIENTS,
                "points": POINTS,
                "c4": c4,
                "c6": c6,
                "discriminant": discriminant,
                "j": j_value,
            }
        },
        load_census_roots(),
    )
    census_target = family_census["targets"][0]
    if (
        family_census["census_family_count"] != 2330
        or census_target["exact_factorization_survivor_count"] != 111
        or census_target["survivors_with_rational_square_parameter"]
        or census_target["exact_j_matches"]
    ):
        raise AssertionError("the bounded six-root Mestre fingerprint changed")

    payload = {
        "status": "bounded numerical construction fingerprint; not family recognition",
        "public_json_sha256": {
            str(curve_id): expected_hash
            for curve_id, (_url, expected_hash) in SOURCES.items()
        },
        "ordered_first_17_x_denominator_root_matches": [list(item) for item in matches],
        "nontrivial_common_denominator_roots": sorted(
            {root for _index, root in matches if root != 1}
        ),
        "height_gram_fit": {
            "curve_351_to_curve_356_scale": scale,
            "relative_frobenius_residual": relative_residual,
            "pearson_correlation_all_289_entries": pearson,
            "pari_version": pari351,
        },
        "bounded_six_root_mestre_census": family_census,
        "interpretation": (
            "Strong evidence that the displayed order preserves a common seventeen-section "
            "template; it does not determine the generic rank or identify an H3/R17 family."
        ),
    }
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
