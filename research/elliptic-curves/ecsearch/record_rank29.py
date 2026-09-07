"""Exact local replay of the Elkies--Klagsbrun rank-at-least-29 curve."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Mapping

from .rank_certification import (
    AffinePoint,
    IndependenceCertificate,
    WeierstrassModel,
    build_independence_certificate,
    verify_independence_certificate,
)


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
BENCHMARKS = PROGRAM_ROOT / "data" / "benchmarks.json"
POINT_DATA = PROGRAM_ROOT / "data" / "elkies_klagsbrun_e29_points.json"


def load_rank29_baseline() -> tuple[WeierstrassModel, tuple[AffinePoint, ...]]:
    """Load the canonical local model and the 29 published points."""

    benchmarks = json.loads(BENCHMARKS.read_text(encoding="utf-8"))
    curve = benchmarks["curves"]["elkies_klagsbrun_e29"]
    point_data = json.loads(POINT_DATA.read_text(encoding="utf-8"))
    assert point_data["curve_id"] == "elkies_klagsbrun_e29"
    assert curve["published_rank_lower_bound"] == 29
    coefficients: WeierstrassModel = tuple(
        Fraction(value) for value in curve["weierstrass_coefficients"]
    )  # type: ignore[assignment]
    points = tuple(
        (Fraction(x_coordinate), Fraction(y_coordinate))
        for x_coordinate, y_coordinate in point_data["points"]
    )
    if len(points) != 29:
        raise AssertionError("the rank-29 baseline must contain exactly 29 points")
    return coefficients, points


def _arithmetic_digest(
    coefficients: WeierstrassModel,
    points: tuple[AffinePoint, ...],
) -> str:
    def text(value: Fraction) -> str:
        return (
            str(value.numerator)
            if value.denominator == 1
            else f"{value.numerator}/{value.denominator}"
        )

    payload = {
        "coefficients": [text(value) for value in coefficients],
        "points": [
            [text(x_coordinate), text(y_coordinate)]
            for x_coordinate, y_coordinate in points
        ],
    }
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def build_rank29_manifest(
    *, maximum_reduction_prime: int = 600
) -> dict[str, object]:
    """Discover and serialize a deterministic exact independence certificate."""

    coefficients, points = load_rank29_baseline()
    certificate = build_independence_certificate(
        coefficients,
        points,
        relation_prime=2,
        maximum_reduction_prime=maximum_reduction_prime,
    )
    return {
        "schema": "elliptic-curves.elkies-klagsbrun-rank29-certificate.v1",
        "claim": "the displayed curve has at least 29 independent rational points",
        "claim_status": "unconditional exact lower bound; no exact-rank claim",
        "target_status": (
            "public rank-29 calibration baseline only; no thirtieth point was found"
        ),
        "curve_id": "elkies_klagsbrun_e29",
        "points_source": "https://web.math.pmf.unizg.hr/~duje/tors/rk29.html",
        "announcement_source": (
            "https://listserv.nodak.edu/cgi-bin/wa.exe?"
            "A2=NMBRTHRY%3Bb9d018b1.2409&S=b"
        ),
        "arithmetic_sha256": _arithmetic_digest(coefficients, points),
        "independence_certificate": certificate.to_json_object(),
        "generation": {
            "command": (
                "python3 elliptic-curves/scripts/run_e29_independence.py --output "
                "artifacts/generated-results/elliptic-curves/"
                "elkies_klagsbrun_e29_independence_v1.json"
            ),
            "maximum_reduction_prime": maximum_reduction_prime,
            "arithmetic": "exact rational and finite-field group operations",
            "external_software": "none",
        },
        "limitations": {
            "conditional_upper_bound": (
                "the authors announce rank at most 29 under GRH for associated "
                "number-field zeta functions; this artifact does not replay it"
            ),
            "rank_30": "no thirtieth independent point is supplied or claimed",
        },
    }


def verify_rank29_manifest(manifest: Mapping[str, object]) -> None:
    """Replay pinned rows without searching for replacement reductions."""

    assert manifest["schema"] == (
        "elliptic-curves.elkies-klagsbrun-rank29-certificate.v1"
    )
    coefficients, points = load_rank29_baseline()
    assert manifest["arithmetic_sha256"] == _arithmetic_digest(coefficients, points)
    raw_certificate = manifest["independence_certificate"]
    if not isinstance(raw_certificate, Mapping):
        raise AssertionError("malformed rank-29 independence certificate")
    certificate = IndependenceCertificate.from_json_object(raw_certificate)
    verify_independence_certificate(coefficients, points, certificate)
