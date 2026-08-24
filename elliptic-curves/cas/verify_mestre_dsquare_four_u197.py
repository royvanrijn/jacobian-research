#!/usr/bin/env python3
"""Replay the pinned exact rank-17 promotion at family 2, ``u=197``.

The default input is the tracked, self-contained certificate in
``artifacts/generated-results``.  The verifier reconstructs the Mestre
quartic and Jacobian, checks all seventeen stored points, recomputes the full
mod-3 finite-reduction certificate, proves the strict conductor cutoff by
rational inequalities, and asks PARI/GP to reconstruct the global curve.
``--discovery-root`` additionally audits the ignored three-chart ratpoints
outputs and the original 28-column discovery certificate.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import factorial
from pathlib import Path
from typing import Any

from search_mestre_dsquare_four import (
    FAMILIES,
    POINT_PATTERN,
    RELATION_PRIME_BOUND,
    ROOT,
    base_parameter,
    canonical_digest,
    known_jacobian_points,
    rational_square_root,
    rational_text,
)
from search_mestre_root_tuple_scale import (
    capped_minimal_curve_data,
    point_digest,
    quartic_point_to_jacobian,
    quartic_value,
    sha256_file,
)
from search_mestre_root_tuple_scale_max200 import mod3_independence_certificate
from verify_mestre_02595143168205_rank13_section import (
    replay as replay_generic_rank13,
)
from verify_mestre_02595143168205_discriminants import (
    replay as replay_discriminants,
)


Q = Fraction
FAMILY_INDEX = 2
PARAMETER_U = Q(197)
PARAMETER_T = Q(337, 394)
PINNED_ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/elliptic_mestre_dsquare_four_u197_rank17.json"
)
EXPECTED_PINNED_ARTIFACT_SHA256 = (
    "f1235d845653219c53d906a06042d4904686feeb42c379ed7f3d83e01d7f0563"
)
EXPECTED_CERTIFICATE_FILE_SHA256 = (
    "5e72faf9211fc96fcb36be8802417a96a6306b6c2e365c643d5f7eae78f01796"
)
EXPECTED_CONDUCTOR_FILE_SHA256 = (
    "ecc92ae24aed9cb2ce09202d074d397ff2aa4dc837c8af6f94dcc7249b36cbd6"
)
EXPECTED_RAW_SHA256 = {
    "raw": "9dedf3c3a609b31341a451b809ed9e4a5c13c30f9c2fd26caa03ddaee89845aa",
    "plus-T": "3d526bc2074f9050c209f4478c8af450156a7ba71a42c407f656b68a49c36cd1",
    "minus-T": "d767044c47c5c15caabbb84c7757abb9e23ad35f9d6efe134f6379d067811b60",
}
EXPECTED_GLOBAL_CURVE = {
    "minimal_model": [
        1,
        1,
        1,
        -1163348683373499147707371416562962,
        15227131493689013260364706485730874765958430844575,
    ],
    "conductor": (
        "2462086522751621334987931952469307556796057284118717977320345864383117775914"
    ),
    "log_conductor": (
        "173.594891144976658977437999799761759808115912593376805006972"
    ),
    "minimal_discriminant": (
        "599083193576231086141477512884203049644313384166071298653060947888798901186742213814158847029018624"
    ),
    "root_number": -1,
}
EXPECTED_POINT_SHA256 = (
    "336945df96bd6d546035367429ab6161106335e7327bdd131e3480a3699a4258"
)
EXPECTED_SUBSET_SHA256 = (
    "a72bc351c1f8450e7fb5cf3032ae7302092b88d3a61f2bb55e87953ffa75f824"
)
EXPECTED_PIVOTS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 18, 21]


def assert_equal(actual: Any, expected: Any, label: str) -> None:
    if actual != expected:
        raise AssertionError(f"{label} changed: {actual!r} != {expected!r}")


def load_canonical_record(path: Path, expected_file_sha256: str) -> dict[str, Any]:
    assert_equal(sha256_file(path), expected_file_sha256, f"SHA256({path.name})")
    record = json.loads(path.read_text())
    declared_digest = record.pop("result_sha256")
    assert_equal(canonical_digest(record), declared_digest, f"digest({path.name})")
    record["result_sha256"] = declared_digest
    return record


def raw_abscissae(raw_root: Path) -> tuple[Fraction, ...]:
    offsets = {"raw": Q(0), "plus-T": PARAMETER_T, "minus-T": -PARAMETER_T}
    abscissae: set[Fraction] = set()
    for chart, offset in offsets.items():
        path = raw_root / f"{chart}.out"
        text = path.read_text()
        assert_equal(
            hashlib.sha256(text.encode()).hexdigest(),
            EXPECTED_RAW_SHA256[chart],
            f"raw SHA256({chart})",
        )
        for line in text.splitlines():
            match = POINT_PATTERN.fullmatch(line.strip())
            if match is None:
                raise AssertionError(f"malformed ratpoints line in {chart}: {line!r}")
            numerator, denominator = map(int, match.groups())
            if denominator:
                abscissae.add(Q(numerator, denominator) + offset)
    return tuple(sorted(abscissae))


def reconstruct_pool(raw_root: Path) -> tuple[
    tuple[tuple[Fraction, Fraction], ...], list[dict[str, str | bool]], int
]:
    family = FAMILIES[FAMILY_INDEX]
    construction = family.construction
    coefficients = construction.primitive_quartic_coefficients(PARAMETER_T)
    known = known_jacobian_points(family, PARAMETER_U)
    by_jacobian_x: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for point in known:
        by_jacobian_x.setdefault(point[0], point)
    known_column_count = len(by_jacobian_x)
    searched_records: list[dict[str, str | bool]] = []
    for x_value in raw_abscissae(raw_root):
        square_root = rational_square_root(quartic_value(coefficients, x_value))
        if square_root is None:
            raise AssertionError("ratpoints output failed the exact original quartic")
        if square_root == 0:
            continue
        point = quartic_point_to_jacobian(
            construction, PARAMETER_T, (x_value, square_root)
        )
        novel = point[0] not in by_jacobian_x
        by_jacobian_x.setdefault(point[0], point)
        searched_records.append(
            {
                "quartic_x": rational_text(x_value),
                "quartic_y_positive": rational_text(square_root),
                "jacobian_x": rational_text(point[0]),
                "novel_jacobian_x": novel,
            }
        )
    return tuple(by_jacobian_x.values()), searched_records, known_column_count


def verify_exact_log_bound(record: dict[str, Any], conductor: int) -> None:
    digits = len(str(conductor))
    assert conductor < 10**digits
    logarithm_bound = Q(231, 100)
    partial_sum = sum(logarithm_bound**degree / factorial(degree) for degree in range(8))
    assert partial_sum > 10
    upper = digits * logarithm_bound
    assert upper < Q(4568, 25)
    assert_equal(record["decimal_digit_count"], digits, "conductor digit count")
    assert_equal(record["conductor_less_than_power_of_ten"], f"10^{digits}", "power bound")
    assert_equal(record["exp_231_over_100_degree_7_partial_sum"], str(partial_sum), "exp bound")
    assert_equal(record["deduced_log_10_upper_bound"], str(logarithm_bound), "log(10) bound")
    assert_equal(record["deduced_log_conductor_upper_bound"], str(upper), "log(N) bound")
    assert_equal(record["strict_target_as_rational"], "4568/25", "strict target")
    assert record["strict_target_proved_exactly"] is True


def load_pinned_artifact(path: Path) -> dict[str, Any]:
    assert_equal(sha256_file(path), EXPECTED_PINNED_ARTIFACT_SHA256, "pinned artifact SHA256")
    artifact = json.loads(path.read_text())
    declared = artifact.pop("result_sha256")
    assert_equal(canonical_digest(artifact), declared, "pinned artifact result digest")
    artifact["result_sha256"] = declared
    return artifact


def replay(
    artifact_path: Path,
    *,
    verify_pari: bool,
    pari_timeout: float,
    discovery_root: Path | None = None,
) -> dict[str, Any]:
    generic_certificate = replay_generic_rank13()
    discriminant_certificate = replay_discriminants()
    assert_equal(
        generic_certificate["family_roots"],
        [0, 25, 95, 143, 168, 205],
        "generic-certificate family roots",
    )
    assert_equal(
        generic_certificate["generic_rank_lower_bound_after_base_change"],
        13,
        "generic rank lower bound",
    )
    assert_equal(
        generic_certificate["generic_companion_identity"][
            "identity_verified_over"
        ],
        "Q[T]",
        "generic companion identity base",
    )
    assert_equal(
        discriminant_certificate["direct"]["degree"],
        20,
        "direct discriminant degree",
    )
    assert_equal(
        discriminant_certificate["pullback"]["degree"],
        40,
        "base-changed discriminant degree",
    )
    assert discriminant_certificate["direct"]["irreducible_over_Q"] is True
    assert discriminant_certificate["direct"]["squarefree_over_Q"] is True
    assert discriminant_certificate["pullback"]["irreducible_over_Q"] is True
    assert discriminant_certificate["pullback"]["squarefree_over_Q"] is True

    family = FAMILIES[FAMILY_INDEX]
    assert_equal(family.roots, (0, 25, 95, 143, 168, 205), "family roots")
    assert_equal(base_parameter(family, PARAMETER_U), PARAMETER_T, "base parameter T")
    artifact = load_pinned_artifact(artifact_path)
    assert_equal(artifact["family"]["roots"], list(family.roots), "artifact roots")
    assert_equal(artifact["specialization"]["u"], "197", "artifact u")
    assert_equal(artifact["specialization"]["T"], "337/394", "artifact T")
    quartic = family.construction.primitive_quartic_coefficients(PARAMETER_T)
    coefficients = family.construction.primitive_jacobian_coefficients(PARAMETER_T)
    assert_equal(
        artifact["specialization"]["primitive_quartic_coefficients_ascending"],
        [rational_text(value) for value in quartic],
        "specialized quartic",
    )
    assert_equal(
        artifact["specialization"]["short_jacobian_coefficients"],
        [rational_text(value) for value in coefficients],
        "specialized Jacobian",
    )
    assert_equal(artifact["specialization"]["global_curve"], EXPECTED_GLOBAL_CURVE, "curve data")
    verify_exact_log_bound(
        artifact["specialization"]["exact_log_conductor_bound"],
        int(EXPECTED_GLOBAL_CURVE["conductor"]),
    )

    point_records = artifact["point_certificate"]["selected_points"]
    points = tuple((Q(record["x"]), Q(record["y"])) for record in point_records)
    assert_equal(len(points), 17, "pinned point count")
    coefficient_a, coefficient_b = coefficients[3], coefficients[4]
    for point in points:
        if point[1] ** 2 != point[0] ** 3 + coefficient_a * point[0] + coefficient_b:
            raise AssertionError("a pinned point missed the specialized Jacobian")
    replayed = mod3_independence_certificate(
        coefficients, points, prime_bound=RELATION_PRIME_BOUND
    )
    assert_equal(
        json.loads(json.dumps(replayed)),
        artifact["point_certificate"]["finite_reduction_certificate"],
        "finite-reduction certificate",
    )
    assert_equal(replayed["combined_exact_rank_over_F3"], 17, "exact rank lower bound")
    assert_equal(replayed["point_sha256"], EXPECTED_SUBSET_SHA256, "certificate points")
    assert_equal(replayed["independent_subset_sha256"], EXPECTED_SUBSET_SHA256, "subset digest")

    discovery_audited = False
    if discovery_root is not None:
        provenance = artifact["discovery_provenance"]
        assert_equal(
            sha256_file(discovery_root / "summary.json"),
            provenance["discovery_summary_sha256"],
            "discovery summary SHA256",
        )
        certificate_path = discovery_root / "point-certificates/f2_u197_1.json"
        conductor_path = discovery_root / "conductor-records/f2_u197_1.json"
        certificate_record = load_canonical_record(
            certificate_path, provenance["discovery_point_certificate_sha256"]
        )
        conductor_record = load_canonical_record(
            conductor_path, provenance["discovery_conductor_record_sha256"]
        )
        assert_equal(certificate_record["global_curve"], EXPECTED_GLOBAL_CURVE, "discovery curve")
        assert_equal(conductor_record["global_curve"], EXPECTED_GLOBAL_CURVE, "discovery conductor")
        pool, searched_records, known_column_count = reconstruct_pool(
            discovery_root / "ratpoints-raw/f2_u197_1"
        )
        assert_equal(known_column_count, 13, "known point columns modulo inverse")
        assert_equal(len(searched_records), 27, "searched finite quartic points")
        assert_equal(len(pool), 28, "discovery pool point count")
        assert_equal(point_digest(pool), EXPECTED_POINT_SHA256, "discovery pool digest")
        assert_equal(searched_records, certificate_record["point_search"]["searched_points"], "search inventory")
        full_replay = mod3_independence_certificate(
            coefficients, pool, prime_bound=RELATION_PRIME_BOUND
        )
        assert_equal(
            json.loads(json.dumps(full_replay)),
            certificate_record["finite_reduction_certificate"],
            "discovery finite-reduction certificate",
        )
        selected_from_pool = tuple(pool[index - 1] for index in EXPECTED_PIVOTS)
        assert_equal(selected_from_pool, points, "pinned subset from discovery pool")
        discovery_audited = True

    if verify_pari:
        recomputed_curve = capped_minimal_curve_data(
            coefficients, timeout=pari_timeout, stack_bytes=512_000_000
        )
        assert_equal(recomputed_curve, EXPECTED_GLOBAL_CURVE, "PARI global curve replay")

    return {
        "status": "verified exact rank lower bound and global curve data",
        "family_roots": list(family.roots),
        "u": rational_text(PARAMETER_U),
        "T": rational_text(PARAMETER_T),
        "minimal_model": EXPECTED_GLOBAL_CURVE["minimal_model"],
        "conductor": EXPECTED_GLOBAL_CURVE["conductor"],
        "log_conductor": EXPECTED_GLOBAL_CURVE["log_conductor"],
        "root_number": EXPECTED_GLOBAL_CURVE["root_number"],
        "pinned_point_count": len(points),
        "generic_rank_lower_bound_after_base_change": generic_certificate[
            "generic_rank_lower_bound_after_base_change"
        ],
        "generic_companion_identity_verified_over": generic_certificate[
            "generic_companion_identity"
        ]["identity_verified_over"],
        "primitive_discriminant_degrees": {
            "direct_T": discriminant_certificate["direct"]["degree"],
            "base_changed_u": discriminant_certificate["pullback"]["degree"],
        },
        "exact_rank_lower_bound": replayed["combined_exact_rank_over_F3"],
        "certificate_primes": replayed["certificate_primes"],
        "torsion_exclusion": replayed["rational_3_torsion_exclusion"],
        "point_sha256": EXPECTED_SUBSET_SHA256,
        "pinned_artifact_sha256": EXPECTED_PINNED_ARTIFACT_SHA256,
        "strict_log_conductor_cutoff_proved_exactly": True,
        "pari_recomputed": verify_pari,
        "discovery_raw_audited": discovery_audited,
        "claim_limit": "rank lower bound only; the bounded point search is not an upper bound",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--artifact",
        type=Path,
        default=PINNED_ARTIFACT,
    )
    parser.add_argument("--discovery-root", type=Path)
    parser.add_argument("--skip-pari", action="store_true")
    parser.add_argument("--pari-timeout", type=float, default=120.0)
    args = parser.parse_args()
    print(
        json.dumps(
            replay(
                args.artifact,
                verify_pari=not args.skip_pari,
                pari_timeout=args.pari_timeout,
                discovery_root=args.discovery_root,
            ),
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
