#!/usr/bin/env python3
"""Run a Magma-free relative rank test for the Fermigier rank-20 curve.

Run this file with Sage, not CPython::

    sage -python elliptic-curves/cas/run_fermigier_rank20_pari_descent.py

The exact rank-20 basis and its mod-2 certificate are loaded by the corrected
relative-descent builder.  PARI's ellrank is then called with those 20 known
points.  ellrank performs a 2-descent and Cassels-pairing restrictions, but
its cubic-field class-group data are GRH-conditional unless separately
certified with ``bnfcertify``.  Consequently this runner never promotes a
returned upper endpoint to an unconditional rank theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict
from fractions import Fraction
import json
from pathlib import Path
import sys
import time
from typing import Iterable, Sequence

CAS_ROOT = Path(__file__).resolve().parent
PROGRAM_ROOT = CAS_ROOT.parent
sys.path.insert(0, str(CAS_ROOT))
sys.path.insert(0, str(PROGRAM_ROOT))

PROTOCOL = "R20PARI"
PROTOCOL_VERSION = 1
KNOWN_RANK = 20

EXPECTED_MANIFEST_SHA256 = "8416e835887236e9e4eafcb01384a710ce4f1be0628701a97f4a7d7a07fe63b1"
EXPECTED_MINIMAL_BASIS_SHA256 = "6fbdc4367d52ca92cfdfef8b0cc71347b2943784df3780a7a72646b0caff898e"
EXPECTED_CANDIDATE_KEY = "fermigier-mestre-v1:u=28917/20"


@dataclass(frozen=True)
class DescentBasis:
    model: tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
    points: tuple[tuple[Fraction, Fraction], ...]
    basis_id: str
    basis_sha256: str
    candidate_record_sha256: str
    manifest_sha256: str
    mod2_rank: int
    mod2_certified: bool


def load_descent_basis(manifest_path: Path, candidate_record_path: Path) -> DescentBasis:
    """Replay the corrected bounded-saturation basis without requiring Magma."""
    import hashlib
    from elliptic_candidate_record import (
        WeierstrassChange,
        change_weierstrass_model,
        is_on_weierstrass_curve,
        model_from_record,
        point_from_record,
        point_sequence_sha256,
        source_point_to_target,
        target_point_to_source,
        validate_candidate_identity,
        verify_finite_quotient_certificate,
    )

    manifest_bytes = manifest_path.read_bytes()
    manifest_sha256 = hashlib.sha256(manifest_bytes).hexdigest()
    if manifest_sha256 != EXPECTED_MANIFEST_SHA256:
        raise ValueError("the pinned near-miss manifest changed")
    manifest = json.loads(manifest_bytes)
    if manifest.get("schema") != "elliptic-curves.fermigier-rank20-near-miss.v1":
        raise ValueError("unexpected near-miss schema")

    candidate_bytes = candidate_record_path.read_bytes()
    candidate_record_sha256 = hashlib.sha256(candidate_bytes).hexdigest()
    candidate = json.loads(candidate_bytes)
    if candidate.get("schema") != "elliptic-curves.candidate-record.v1":
        raise ValueError("unexpected candidate-record schema")
    validate_candidate_identity(candidate)
    if candidate["identity"]["candidate_key"] != EXPECTED_CANDIDATE_KEY:
        raise ValueError("unexpected Fermigier specialization")

    bridge = candidate["artifact_namespace_bridge"]["imported_ecsearch_namespace"]
    reference = next(
        (item for item in bridge if item["path"].endswith("fermigier_rank20_near_miss_v1.json")),
        None,
    )
    if reference is None or reference["sha256"] != manifest_sha256:
        raise ValueError("candidate record is not pinned to this near-miss manifest")

    models = candidate["models"]
    canonical = model_from_record(models["canonical_generalized"]["coefficients_a1_a2_a3_a4_a6"])
    legacy = model_from_record(models["legacy_normalized_short_jacobian"]["coefficients"])
    minimal = model_from_record(models["global_minimal"]["coefficients"])
    if minimal != model_from_record(manifest["global_curve"]["minimal_model"]):
        raise ValueError("candidate and manifest minimal models differ")

    transforms = candidate["exact_transformations"]
    c_to_l = WeierstrassChange.from_values(
        transforms["canonical_to_legacy_normalized_short"]["change_u_r_s_t"]
    )
    c_to_m = WeierstrassChange.from_values(
        transforms["canonical_to_global_minimal"]["change_u_r_s_t"]
    )
    if change_weierstrass_model(canonical, c_to_l) != legacy:
        raise ArithmeticError("canonical-to-legacy replay failed")
    if change_weierstrass_model(canonical, c_to_m) != minimal:
        raise ArithmeticError("canonical-to-minimal replay failed")

    saturation = candidate["bounded_saturation_status"]
    if saturation["returned_basis_count"] != 20:
        raise ValueError("bounded-saturation basis no longer has 20 points")
    legacy_points = tuple(point_from_record(item) for item in saturation["returned_legacy_basis"])
    if len(set(legacy_points)) != 20:
        raise ValueError("bounded-saturation basis is not distinct")
    if any(not is_on_weierstrass_curve(legacy, point) for point in legacy_points):
        raise ArithmeticError("bounded-saturation point is off the legacy model")

    certificate = candidate["independent_cas_cross_certificates"][
        "bounded_saturation_candidate_basis"
    ]["mod_2"]
    if (
        certificate["relation_prime"] != 2
        or certificate["point_count"] != 20
        or certificate["combined_rank_over_relation_field"] != 20
        or not certificate["certified_independent"]
    ):
        raise ValueError("full mod-2 certificate is missing")
    if certificate["point_sequence_sha256"] != point_sequence_sha256(legacy_points):
        raise ValueError("mod-2 certificate covers a different basis")
    verify_finite_quotient_certificate(legacy, legacy_points, certificate)

    canonical_points = tuple(target_point_to_source(point, c_to_l) for point in legacy_points)
    minimal_points = tuple(source_point_to_target(point, c_to_m) for point in canonical_points)
    if any(not is_on_weierstrass_curve(minimal, point) for point in minimal_points):
        raise ArithmeticError("transported point is off the global minimal model")
    basis_sha256 = point_sequence_sha256(minimal_points)
    if basis_sha256 != saturation["transported_basis_sha256"]["global_minimal"]:
        raise ValueError("transported basis hash differs from candidate record")
    if basis_sha256 != EXPECTED_MINIMAL_BASIS_SHA256:
        raise ValueError("pinned global-minimal basis changed")

    return DescentBasis(
        model=minimal,
        points=minimal_points,
        basis_id="bounded_saturation_candidate",
        basis_sha256=basis_sha256,
        candidate_record_sha256=candidate_record_sha256,
        manifest_sha256=manifest_sha256,
        mod2_rank=20,
        mod2_certified=True,
    )


@dataclass(frozen=True)
class PariPass:
    effort: int
    pari_lower: int
    pari_upper: int
    effective_lower: int
    sha_pairing_rank: int
    returned_points: int
    elapsed_seconds: float
    classification: str


def classify_bounds(pari_lower: int, pari_upper: int, known_rank: int = KNOWN_RANK) -> tuple[int, str]:
    """Combine PARI's interval with the repository's certified lower bound."""
    if pari_lower < 0 or pari_upper < pari_lower:
        raise ValueError("invalid PARI rank interval")
    if pari_upper < known_rank:
        raise ArithmeticError(
            f"PARI upper bound {pari_upper} contradicts certified rank >= {known_rank}"
        )
    lower = max(known_rank, pari_lower)
    if lower > pari_upper:
        raise ArithmeticError("combined rank interval is empty")
    if lower == pari_upper == known_rank:
        return lower, "P0_grh_conditional_rank20"
    if lower >= known_rank + 1:
        return lower, "P3_rank_at_least21_signal"
    return lower, "P2_grh_conditional_residual_rank_interval"


def parse_efforts(text: str) -> tuple[int, ...]:
    values = tuple(int(piece.strip()) for piece in text.split(",") if piece.strip())
    if not values or any(value < 0 for value in values):
        raise ValueError("--efforts must be a nonempty comma-separated list of nonnegative integers")
    # Preserve order while dropping duplicates.
    return tuple(dict.fromkeys(values))


def sage_q(value: Fraction | int):
    from sage.all import QQ

    value = Fraction(value)
    return QQ(value.numerator) / QQ(value.denominator)


def _pari_int(value) -> int:
    # cypari2 Gen -> Sage Integer/int conversion varies slightly across releases.
    try:
        return int(value)
    except TypeError:
        return int(value.sage())


def _pari_points(points) -> object:
    from sage.all import pari

    return pari([[point[0], point[1]] for point in points])


def _ellrankinit(pari_curve):
    # cypari2 exposes GP functions both as bound Gen methods and through pari.
    try:
        return pari_curve.ellrankinit()
    except (AttributeError, TypeError):
        from sage.all import pari

        return pari.ellrankinit(pari_curve)


def _ellrank(rank_context, effort: int, points):
    try:
        return rank_context.ellrank(effort, points)
    except (AttributeError, TypeError):
        from sage.all import pari

        return pari.ellrank(rank_context, effort, points)


def run_pari_pass(rank_context, known_points, effort: int) -> PariPass:
    started = time.monotonic()
    result = _ellrank(rank_context, effort, _pari_points(known_points))
    elapsed = time.monotonic() - started

    # PARI ellrank returns [r1, r2, s, L].  The documented s is the rank of
    # Sha[2]/2Sha[4] detected by Cassels-pairing restrictions.
    if len(result) != 4:
        raise RuntimeError(f"unexpected PARI ellrank result of length {len(result)}")
    pari_lower = _pari_int(result[0])
    pari_upper = _pari_int(result[1])
    sha_pairing_rank = _pari_int(result[2])
    returned_points = len(result[3])
    effective_lower, classification = classify_bounds(pari_lower, pari_upper)
    return PariPass(
        effort=effort,
        pari_lower=pari_lower,
        pari_upper=pari_upper,
        effective_lower=effective_lower,
        sha_pairing_rank=sha_pairing_rank,
        returned_points=returned_points,
        elapsed_seconds=elapsed,
        classification=classification,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(
            "artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json"
        ),
    )
    parser.add_argument(
        "--candidate-record",
        type=Path,
        default=Path(
            "artifacts/generated-results/elliptic-curves/elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
        ),
    )
    parser.add_argument(
        "--efforts",
        default="0,1,2",
        help="comma-separated PARI ellrank effort levels; runtime grows roughly cubically",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/local/elliptic-curves/fermigier_rank20_pari_descent.json"),
    )
    parser.add_argument(
        "--continue-after-new-rank",
        action="store_true",
        help="continue later effort passes even after PARI proves rank >= 21",
    )
    args = parser.parse_args()
    efforts = parse_efforts(args.efforts)

    try:
        from sage.all import EllipticCurve, QQ, version as sage_version, pari
    except ImportError as exc:
        raise SystemExit(
            "This program must run inside SageMath: sage -python " + str(Path(__file__))
        ) from exc

    basis = load_descent_basis(args.manifest, args.candidate_record)
    if basis.mod2_rank != 20 or not basis.mod2_certified:
        raise ArithmeticError("refusing descent without the exact full mod-2 basis certificate")

    ainvs = [sage_q(value) for value in basis.model]
    E = EllipticCurve(QQ, ainvs)
    known_points = [E(sage_q(x), sage_q(y)) for x, y in basis.points]
    if len(known_points) != 20 or len(set(known_points)) != 20:
        raise ArithmeticError("the pinned descent basis did not replay as 20 distinct Sage points")

    print(
        f"{PROTOCOL}|version={PROTOCOL_VERSION}|stage=input|sage={sage_version()}"
        f"|known=20|basis={basis.basis_id}|basis_sha256={basis.basis_sha256}",
        flush=True,
    )
    print(f"{PROTOCOL}|stage=ellrankinit|status=start", flush=True)
    init_started = time.monotonic()
    rank_context = _ellrankinit(E.pari_curve())
    init_elapsed = time.monotonic() - init_started
    print(f"{PROTOCOL}|stage=ellrankinit|status=complete|seconds={init_elapsed:.6f}")

    passes: list[PariPass] = []
    for effort in efforts:
        print(f"{PROTOCOL}|stage=ellrank|status=start|effort={effort}", flush=True)
        current = run_pari_pass(rank_context, known_points, effort)
        passes.append(current)
        print(
            f"{PROTOCOL}|stage=ellrank|status=complete|effort={effort}"
            f"|pari_lower={current.pari_lower}|pari_upper={current.pari_upper}"
            f"|effective_lower={current.effective_lower}"
            f"|sha_pairing_rank={current.sha_pairing_rank}"
            f"|returned_points={current.returned_points}"
            f"|seconds={current.elapsed_seconds:.6f}"
            f"|classification={current.classification}",
            flush=True,
        )
        if current.classification == "P0_grh_conditional_rank20":
            break
        if current.classification == "P3_rank_at_least21_signal" and not args.continue_after_new_rank:
            break

    best_lower = max(item.effective_lower for item in passes)
    best_upper = min(item.pari_upper for item in passes)
    if best_lower > best_upper:
        raise ArithmeticError("independent PARI passes produced contradictory rank intervals")
    if best_lower == best_upper == 20:
        final = "P0_grh_conditional_rank20"
    elif best_lower >= 21:
        final = "P3_rank_at_least21_signal"
    else:
        final = "P2_grh_conditional_residual_rank_interval"

    payload = {
        "schema": "elliptic-curves.fermigier-rank20-pari-descent.v1",
        "engine": "SageMath/PARI ellrank",
        "mathematical_status": "grh_conditional_external_cas_rank_interval",
        "upper_bound_hypothesis": "GRH unless the cubic-field BNF is separately certified",
        "basis": {
            "id": basis.basis_id,
            "count": 20,
            "sha256": basis.basis_sha256,
            "mod2_rank": basis.mod2_rank,
            "mod2_certified": basis.mod2_certified,
            "candidate_record_sha256": basis.candidate_record_sha256,
            "manifest_sha256": basis.manifest_sha256,
        },
        "rank_interval": [best_lower, best_upper],
        "classification": final,
        "passes": [asdict(item) for item in passes],
        "notes": [
            "The lower endpoint is combined with the repository's exact rank>=20 certificate.",
            "PARI ellrank's upper endpoint uses provisional cubic-field BNF data and is not unconditional without bnfcertify.",
            "A nonzero sha_pairing_rank records the obstruction rank detected by PARI; it is not a Mordell-Weil rank contribution.",
            "Returned rational points may still be checked independently to obtain unconditional lower bounds.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=summary|rank_lower={best_lower}|rank_upper={best_upper}"
        f"|classification={final}|output={args.output}",
        flush=True,
    )


if __name__ == "__main__":
    main()
