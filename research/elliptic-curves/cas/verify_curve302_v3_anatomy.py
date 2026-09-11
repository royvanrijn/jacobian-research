#!/usr/bin/env sage -python
"""Verify the locally checkable claims in the Curve302 V3 anatomy note.

The supplied verification JSON is an attestation to a larger external bundle.
This checker deliberately does not promote its numerical measurements to a
local replay: the files named by that hash manifest must be present before
their hashes can be checked.
"""

from __future__ import annotations

import hashlib
import json
import runpy
from pathlib import Path

import sympy as sp
from sage.all import EllipticCurve, GF, QQ, matrix, vector


ROOT = Path(__file__).resolve().parents[3]
EC = ROOT / "research/elliptic-curves"
DATA = EC / "data/curve302_v3_anatomy_verification_v1.json"
OUTPUT = (
    ROOT
    / "research/artifacts/generated-results/elliptic-curves"
    / "curve302_v3_anatomy_local_verification_v1.json"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def symbolic_chart_checks() -> dict[str, bool]:
    A, B, a, b, t, z = sp.symbols("A B a b t z")
    F = t**4 - 6 * a * t**2 - 8 * b * t - 3 * a**2 - 4 * A
    x = (t**2 - a + z) / 2
    y = t * (x - a) - b
    relations = sp.groebner(
        [z**2 - F, b**2 - a**3 - A * a - B],
        z,
        b,
        t,
        a,
        A,
        B,
        order="lex",
        domain=sp.QQ,
    )

    def numerator_remainder(expression):
        numerator = sp.together(expression).as_numer_denom()[0]
        return relations.reduce(sp.expand(numerator))[1]

    # Double P and then add -Q=(a,-b).
    doubling_slope = (3 * x**2 + A) / (2 * y)
    x2 = sp.cancel(doubling_slope**2 - 2 * x)
    y2 = sp.cancel(-y + doubling_slope * (x - x2))
    subtraction_slope = sp.cancel((-b - y2) / (a - x2))
    x_residual = sp.cancel(subtraction_slope**2 - x2 - a)
    N = (
        a * t**4
        + 4 * b * t**3
        + (6 * a**2 + 4 * A) * t**2
        + 4 * a * b * t
        + a**3
        + 4 * B
    )

    # The opposite quartic root maps to Q-P: the two images add to Q.
    x_other = (t**2 - a - z) / 2
    y_other = t * (x_other - a) - b
    deck_slope = sp.cancel((y_other - y) / (x_other - x))
    x_sum = sp.cancel(deck_slope**2 - x - x_other)
    y_sum = sp.cancel(-y + deck_slope * (x - x_sum))

    checks = {
        "quartic_point_maps_to_E": numerator_remainder(y**2 - x**3 - A * x - B) == 0,
        "inverse_ordinate_is_2x_plus_a_minus_t2": sp.expand(2 * x + a - t**2 - z) == 0,
        "deck_involution_sums_to_Q_x": numerator_remainder(x_sum - a) == 0,
        "deck_involution_sums_to_Q_y": numerator_remainder(y_sum - b) == 0,
        "covering_x_of_2P_minus_Q_is_N_over_F": numerator_remainder(x_residual - N / F) == 0,
    }
    if not all(checks.values()):
        raise AssertionError(f"symbolic chart verification failed: {checks}")
    return checks


def curve302_checks() -> dict:
    public = runpy.run_path(str(EC / "cas/icarm_curve302.py"))
    E = EllipticCurve(QQ, list(public["GENERAL_WEIERSTRASS_COEFFICIENTS"]))
    points = [E(QQ(x), QQ(y)) for x, y in public["POINTS"]]
    if len(points) != 31 or E.discriminant() == 0:
        raise AssertionError("public Curve302 model or point count changed")

    reduction_orders = {}
    for prime, expected in ((17, 26), (31, 43)):
        if E.discriminant().valuation(prime) != 0:
            raise AssertionError(f"bad reduction at {prime}")
        order = int(E.change_ring(GF(prime)).cardinality())
        if order != expected:
            raise AssertionError(f"unexpected reduction order at {prime}: {order}")
        reduction_orders[str(prime)] = order

    parent_path = (
        ROOT
        / "research/artifacts/generated-results/elliptic-curves"
        / "curve302_recovered_mw17_parent_v1.json"
    )
    parent = json.loads(parent_path.read_text())
    embedding = matrix(GF(2), parent["basis_embedding_in_public_D"])
    if embedding.dimensions() != (31, 17) or embedding.rank() != 17:
        raise AssertionError("generic basis embedding has unexpected dimensions or rank")
    columns = list(embedding.columns())
    current_rank = embedding.rank()
    complements = []
    for index in range(31):
        axis = vector(GF(2), [int(j == index) for j in range(31)])
        trial = matrix(GF(2), columns + [axis]).transpose()
        if trial.rank() > current_rank:
            columns.append(axis)
            complements.append(index)
            current_rank += 1
    expected_complements = [0, 2, 3, 5, 8, 12, 13, 16, 18, 21, 25, 28, 29, 30]
    if complements != expected_complements or current_rank != 31:
        raise AssertionError(f"unexpected complement indices: {complements}")

    focused_path = (
        ROOT
        / "research/artifacts/generated-results/elliptic-curves"
        / "curve302_focused_point_exposure_v2.json"
    )
    focused = json.loads(focused_path.read_text())
    if focused["parent_height_diagnostic"]["curve302_reduced_parameter"] != "-164518/143797":
        raise AssertionError("Curve302 reduced parent parameter changed")

    return {
        "public_model_nonsingular": True,
        "public_points_checked_exactly": len(points),
        "good_reduction_orders": reduction_orders,
        "rational_torsion_order_divides_gcd": 1,
        "generic_embedding_mod2_rank": int(embedding.rank()),
        "greedy_public_complement_indices_zero_based": complements,
        "reduced_parent_parameter_repository_crosscheck": "-164518/143797",
        "parent_artifact_sha256": sha256(parent_path),
        "focused_exposure_artifact_sha256": sha256(focused_path),
    }


def supplied_bundle_status(attestation: dict) -> dict:
    # The attachment contains only this manifest. Hashes are checkable if a
    # later complete bundle is placed beside it under this directory.
    payload_root = DATA.parent / "curve302_v3_anatomy_payload_v1"
    matched, missing, mismatched = [], [], []
    for relative, expected in sorted(attestation["hashes"].items()):
        path = payload_root / relative
        if not path.is_file():
            missing.append(relative)
        elif sha256(path) == expected:
            matched.append(relative)
        else:
            mismatched.append(relative)
    if mismatched:
        status = "FAIL_HASH_MISMATCH"
    elif missing:
        status = "UNVERIFIED_MISSING_HASHED_PAYLOADS"
    else:
        status = "PASS_ALL_SUPPLIED_PAYLOAD_HASHES"
    return {
        "status": status,
        "payload_root": str(payload_root.relative_to(ROOT)),
        "matched": matched,
        "missing": missing,
        "mismatched": mismatched,
    }


def build() -> dict:
    attestation = json.loads(DATA.read_text())
    result = {
        "schema": "elliptic-curves.curve302-v3-anatomy-local-verification.v1",
        "status": "PASS_LOCAL_CLAIMS_EXTERNAL_MEASUREMENTS_UNVERIFIED",
        "boundary": (
            "Exact chart identities, public Curve302 points, reduction orders, the mod-2 "
            "generic/complement split, and the repository's reduced parent parameter are "
            "checked locally. Numerical kappa tables, 12,082/1,734 pair counts, their "
            "byte-for-byte replay, and the plot remain an external attestation until all "
            "files named by the supplied hash manifest are present."
        ),
        "supplied_attestation": {
            "sha256": sha256(DATA),
            "declared_status": attestation.get("status"),
            "declared_tests_passed": attestation.get("tests_passed"),
            "declared_primary_point_anchor_pairs_checked": attestation.get(
                "primary_point_anchor_pairs_checked"
            ),
            "declared_same_class_point_anchor_pairs_checked": attestation.get(
                "same_class_point_anchor_pairs_checked"
            ),
        },
        "symbolic_chart_checks": symbolic_chart_checks(),
        "curve302_checks": curve302_checks(),
        "external_payload_check": supplied_bundle_status(attestation),
        "elementary_checks": {
            "deck_pair_generates_same_unsaturated_extension": (
                "For Q in S, Q-P lies in S+ZP and P=Q-(Q-P) lies in S+Z(Q-P)."
            ),
            "one_vector_fixed_core_intersection_gain_bound": 1,
            "reported_target31_height_factor_from_bits": 37386135.091249004,
        },
    }
    if result["external_payload_check"]["status"] == "FAIL_HASH_MISMATCH":
        raise AssertionError("one or more supplied payload hashes mismatch")
    return result


def canonical(value: dict) -> str:
    return json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True) + "\n"


def main() -> None:
    expected = canonical(build())
    import sys

    if len(sys.argv) == 2 and sys.argv[1] == "--check":
        if not OUTPUT.is_file() or OUTPUT.read_text() != expected:
            raise SystemExit("CURVE302_V3_ANATOMY_CHECK|status=FAIL")
        print("CURVE302_V3_ANATOMY_CHECK|status=PASS")
        return
    if len(sys.argv) != 1:
        raise SystemExit(f"usage: sage -python {Path(__file__).name} [--check]")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(expected)
    print(f"CURVE302_V3_ANATOMY_VERIFY|status=PASS|output={OUTPUT}")


if __name__ == "__main__":
    main()
