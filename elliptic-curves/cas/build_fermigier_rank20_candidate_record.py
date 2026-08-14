#!/usr/bin/env python3
"""Build the canonical cross-certified record for the Fermigier rank-20 anchor.

This is identity and certificate work only.  It replays a pinned bounded point
search; it never launches a point, score, parameter, or conductor search.  A
single capped PARI call independently checks the minimal model/global data and
asks for a small-prime saturation candidate for the already selected subgroup.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any, Sequence


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
sys.path.insert(0, str(PROGRAM_ROOT))
sys.path.insert(0, str(PROGRAM_ROOT / "cas"))

from ecsearch.fermigier import (  # noqa: E402
    evaluate_polynomial,
    fermigier_canonical_coefficients,
    quartic_point_to_canonical_point,
)
from ecsearch.fermigier_near_miss import canonical_ratpoints_output  # noqa: E402
from ecsearch.fermigier_rank import (  # noqa: E402
    parse_ratpoints_output,
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import (  # noqa: E402
    negate_rational_point,
    subtract_rational_points,
)
from elliptic_candidate_record import (  # noqa: E402
    CANDIDATE_SCHEMA,
    WeierstrassChange,
    binary_quartic_invariants_low_to_high,
    build_finite_quotient_certificate,
    canonical_candidate_identity,
    change_weierstrass_model,
    fraction_text,
    is_on_weierstrass_curve,
    model_record,
    point_record,
    point_sequence_sha256,
    sha256_file,
    source_point_to_target,
    stable_json_sha256,
    target_point_to_source,
    validate_candidate_identity,
    verify_finite_quotient_certificate,
    weierstrass_invariants,
)
from fermigier_mestre import FermigierMestreFamily  # noqa: E402


Q = Fraction
FAMILY_ID = "fermigier-mestre-v1"
ANCHOR_U = Q(28_917, 20)
LITERAL_SHIFT = 2 * ANCHOR_U
SATURATION_PRIME_BOUND = 20
CERTIFICATE_PRIME_BOUND = 2_000
GP_TIMEOUT_SECONDS = 30.0
GP_STACK_BYTES = 1_000_000_000

IMPORTED_NEAR_MISS = Path(
    "artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json"
)
IMPORTED_RANK_CERTIFICATES = Path(
    "artifacts/generated-results/elliptic-curves/fermigier_rank_certificates_v1.json"
)
GENERIC_RANK_THEOREM = Path(
    "artifacts/generated-results/elliptic_fermigier_generic_rank_exact.json"
)
EXCEPTIONAL_TRANSPORT = Path(
    "artifacts/generated-results/elliptic_fermigier_exceptional_transport.json"
)
EXCEPTIONAL_QUOTIENT_BALL = Path(
    "artifacts/generated-results/elliptic_fermigier_exceptional_quotient_ball.json"
)
BIDEGREE21_P13_R20E1_PILOT = Path(
    "artifacts/generated-results/"
    "elliptic_fermigier_bidegree21_p13_r20e1_pilot.json"
)
BIDEGREE21_ALL80 = Path(
    "artifacts/generated-results/elliptic_fermigier_bidegree21_all80.json"
)
EXCEPTIONAL_PAIR_SIMULTANEOUS_H200000 = Path(
    "artifacts/generated-results/"
    "elliptic_fermigier_exceptional_pair_simultaneous_h200000.json"
)
BIDEGREE21_P13_R20E1_NONLINEAR_POINTS_H1024 = Path(
    "artifacts/generated-results/"
    "elliptic_fermigier_bidegree21_p13_r20e1_nonlinear_points_h1024.json"
)
LEGACY_EXPLICIT_FORMULA = Path(
    "artifacts/generated-results/elliptic_fermigier_rank20_28917_20_explicit_formula_delta22.json"
)
LEGACY_NEIGHBORHOOD_AUDIT = Path(
    "artifacts/generated-results/elliptic_fermigier_rank20_adapter_neighborhood_audit.json"
)
LEGACY_HIGH_POWER_CRT = Path(
    "artifacts/generated-results/elliptic_fermigier_high_power_crt_gauss.json"
)
DEFAULT_OUTPUT = Path(
    "artifacts/generated-results/"
    "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
)

EXPECTED_HASHES = {
    IMPORTED_NEAR_MISS: "8416e835887236e9e4eafcb01384a710ce4f1be0628701a97f4a7d7a07fe63b1",
    IMPORTED_RANK_CERTIFICATES: "94fc64d7f1744f6a20a0396d32914cd36330107db2538e03ee95cc3e32927051",
    GENERIC_RANK_THEOREM: "61bf11ae14db1aedcf7809697c96e78f2c2978c22f4e7b8fd894de74628b3de7",
    EXCEPTIONAL_TRANSPORT: "a767e849119d4eb974eb8e85536031413c6d52a59151933239fa141235de5777",
    EXCEPTIONAL_QUOTIENT_BALL: "4e1f49f57fa9448b1172a3be4f16501138b21a3a7444ec39b27292f7430b1362",
    BIDEGREE21_P13_R20E1_PILOT: "423bec6bd9545783da0a550c1abb44bb6ac096c361011976ab3b209028341bae",
    BIDEGREE21_ALL80: "2c3aa7a8fc57ad7160397506e8db47bb07ea8c988bab87c9e51b1529000301f5",
    EXCEPTIONAL_PAIR_SIMULTANEOUS_H200000: "0a1a1ac50ac35689b4134106e4dd1469553363e15dd9c46a8c6f19358ec69394",
    BIDEGREE21_P13_R20E1_NONLINEAR_POINTS_H1024: "dd281569a1da8eb1c07a635faecb8b9f27269751c2639fb6b94f3a1bada46310",
    LEGACY_EXPLICIT_FORMULA: "9e8262d81026557999409860850cef55540a89a384d10723c4270e57261ebce5",
    LEGACY_NEIGHBORHOOD_AUDIT: "0eef1ad22211d9b8f6b8cdcec3e1c8829322f2889195a2f1527b03465e799615",
    LEGACY_HIGH_POWER_CRT: "2e42f162a7fffe11e43f5adf283900e3a6cee14a781f62c2ac252746eaae0c70",
}

EXPECTED_STRUCTURAL_SCRIPT_HASHES = {
    Path("elliptic-curves/cas/classify_fermigier_exceptional_quotient_ball.py"):
        "b5ed9a18d832289f06cf67c12329e64f82c8618ccdbd220fb10dd48e1d514f9c",
    Path("elliptic-curves/cas/analyze_fermigier_bidegree21_pilot.py"):
        "e35f3a78097a2ea9dbc43048d122b306cfcb261c716fed83ad60d734c44a9089",
    Path("elliptic-curves/cas/analyze_fermigier_bidegree21_all80.py"):
        "403238fc23d3a793906b37d75ce4d7eace5026e5ab41a632aa06e4e358ce859a",
    Path("elliptic-curves/cas/search_fermigier_exceptional_pair_simultaneous_h200000.py"):
        "cf1b6740e6127e7cbf92cf49dccd30c71c8ca6113c2d00fa1ec535a4fc0e0f01",
    Path("elliptic-curves/cas/search_fermigier_bidegree21_nonlinear_points.py"):
        "1dbe3cb7b95991d671f4df64109012afe36ac568556224298e8d8892941b3044",
}

EXPECTED_MINIMAL_CHANGE = WeierstrassChange.from_values(
    (
        "7/50",
        "207130170610471437/20000",
        "-43/100",
        "-5178254265262284553/1000000",
    )
)


def repository_path(relative: Path | str) -> Path:
    return REPOSITORY_ROOT / Path(relative)


def checked_json(relative: Path) -> dict[str, Any]:
    path = repository_path(relative)
    actual_hash = sha256_file(path)
    expected_hash = EXPECTED_HASHES[relative]
    if actual_hash != expected_hash:
        raise AssertionError(
            f"frozen input {relative} changed: {actual_hash} != {expected_hash}"
        )
    return json.loads(path.read_text())


def gp_rational(value: Fraction | int) -> str:
    value = Q(value)
    return f"({value.numerator}/{value.denominator})"


def gp_point(point: tuple[Fraction, Fraction]) -> str:
    return f"[{gp_rational(point[0])},{gp_rational(point[1])}]"


def run_capped_gp_cross_check(
    canonical_model: Sequence[Fraction],
    legacy_short_model: Sequence[Fraction],
    selected_legacy_points: Sequence[tuple[Fraction, Fraction]],
) -> dict[str, Any]:
    """Run one no-retry PARI model/global/small-prime-saturation replay."""

    gp = shutil.which("gp")
    if gp is None:
        raise FileNotFoundError("PARI/GP executable 'gp' is required")
    canonical = ",".join(gp_rational(value) for value in canonical_model)
    legacy = ",".join(gp_rational(value) for value in legacy_short_model)
    points = ",".join(gp_point(point) for point in selected_legacy_points)
    program = "\n".join(
        (
            "default(realprecision,100);",
            f"default(parisizemax,{GP_STACK_BYTES});",
            f"default(parisize,{min(GP_STACK_BYTES, 256_000_000)});",
            f"Ec=ellinit([{canonical}]);",
            "v=0;M=ellminimalmodel(Ec,&v);",
            'print("MIN_MODEL\\t",M.a1,"\\t",M.a2,"\\t",M.a3,"\\t",M.a4,"\\t",M.a6);',
            'print("MIN_CHANGE\\t",v[1],"\\t",v[2],"\\t",v[3],"\\t",v[4]);',
            'print("CONDUCTOR\\t",ellglobalred(M)[1]);',
            'print("MIN_DISC\\t",M.disc);',
            'print("ROOT_NUMBER\\t",ellrootno(M));',
            'print("PARI_VERSION\\t",version());',
            f"Es=ellinit([{legacy}]);",
            f"P=[{points}];",
            f"S=ellsaturation(Es,P,{SATURATION_PRIME_BOUND});",
            'print("SAT_COUNT\\t",#S);',
            'for(i=1,#S,print("SAT_POINT\\t",S[i][1],"\\t",S[i][2]));',
            "quit",
        )
    ) + "\n"
    try:
        completed = subprocess.run(
            [gp, "-q", "-f"],
            input=program,
            text=True,
            capture_output=True,
            check=True,
            timeout=GP_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"single PARI cross-check exceeded {GP_TIMEOUT_SECONDS}s; no retry"
        ) from error
    combined_output = completed.stdout + completed.stderr
    fatal_lines = [
        line
        for line in combined_output.splitlines()
        if "***" in line and "Warning:" not in line
    ]
    if fatal_lines:
        raise RuntimeError("\n".join(fatal_lines))
    values: dict[str, list[list[str]]] = {}
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if fields and fields[0] in {
            "MIN_MODEL",
            "MIN_CHANGE",
            "CONDUCTOR",
            "MIN_DISC",
            "ROOT_NUMBER",
            "PARI_VERSION",
            "SAT_COUNT",
            "SAT_POINT",
        }:
            values.setdefault(fields[0], []).append(fields[1:])

    def one(label: str) -> list[str]:
        records = values.get(label, [])
        if len(records) != 1:
            raise AssertionError(f"PARI emitted {len(records)} {label} records")
        return records[0]

    minimum = tuple(Q(value) for value in one("MIN_MODEL"))
    change = WeierstrassChange.from_values(one("MIN_CHANGE"))
    saturated = tuple(
        (Q(fields[0]), Q(fields[1])) for fields in values.get("SAT_POINT", [])
    )
    if int(one("SAT_COUNT")[0]) != len(saturated):
        raise AssertionError("PARI's saturation point count changed")
    return {
        "minimal_model": minimum,
        "minimal_change": change,
        "conductor": int(one("CONDUCTOR")[0]),
        "minimal_discriminant": int(one("MIN_DISC")[0]),
        "root_number": int(one("ROOT_NUMBER")[0]),
        "pari_version": one("PARI_VERSION")[0],
        "saturated_legacy_points": saturated,
        "program_sha256": hashlib.sha256(program.encode()).hexdigest(),
        "timeout_seconds": GP_TIMEOUT_SECONDS,
        "stack_bytes": GP_STACK_BYTES,
        "retried": False,
    }


def reconstruct_pool(
    specialization: Any,
    searched_quartic_points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[tuple[tuple[Fraction, Fraction], ...], list[dict[str, Any]]]:
    """Replay the 115 differences while retaining exact source provenance."""

    curve = specialization.canonical_model
    origin = specialization.canonical_points[0]
    differences = list(specialization.section_differences)
    provenance = [
        {
            "kind": "generic-section-difference",
            "definition": f"canonical_image(Q_{index + 1})-canonical_image(Q_0)",
            "generic_basis_index_zero_based": index,
        }
        for index in range(12)
    ]
    seen = set(differences)
    seen.update(negate_rational_point(curve, point) for point in differences)
    for signed_index, quartic_point in enumerate(searched_quartic_points):
        canonical_point = quartic_point_to_canonical_point(
            specialization.quartic_model, quartic_point
        )
        difference = subtract_rational_points(curve, canonical_point, origin)
        if difference is None or difference in seen:
            continue
        differences.append(difference)
        seen.add(difference)
        seen.add(negate_rational_point(curve, difference))
        provenance.append(
            {
                "kind": "bounded-quartic-point-difference",
                "definition": "canonical_image(R)-canonical_image(Q_0)",
                "signed_quartic_point_index_zero_based": signed_index,
                "abscissa_index_zero_based": signed_index // 2,
                "ordinate_sign": "positive" if signed_index % 2 == 0 else "negative",
                "quartic_x": fraction_text(quartic_point[0]),
                "quartic_y": fraction_text(quartic_point[1]),
            }
        )
    return tuple(differences), provenance


def experiment_ledger(
    audit: dict[str, Any],
    near_miss: dict[str, Any],
    explicit_formula: dict[str, Any],
    high_power: dict[str, Any],
    exceptional_transport: dict[str, Any],
    exceptional_quotient_ball: dict[str, Any],
    bidegree21_pilot: dict[str, Any],
    bidegree21_all80: dict[str, Any],
    simultaneous_h200000: dict[str, Any],
    nonlinear_points_h1024: dict[str, Any],
) -> dict[str, Any]:
    """Classify every audited Fermigier artifact without inferring absence."""

    paths = [
        item
        for item in audit["input_inventory"]["files"]
        if item["path"].startswith("artifacts/generated-results/")
        and "fermigier" in item["path"].lower()
    ]
    paths.extend(
        (
            {
                "path": str(LEGACY_NEIGHBORHOOD_AUDIT),
                "sha256": EXPECTED_HASHES[LEGACY_NEIGHBORHOOD_AUDIT],
                "size_bytes": repository_path(LEGACY_NEIGHBORHOOD_AUDIT).stat().st_size,
            },
            {
                "path": str(LEGACY_HIGH_POWER_CRT),
                "sha256": EXPECTED_HASHES[LEGACY_HIGH_POWER_CRT],
                "size_bytes": repository_path(LEGACY_HIGH_POWER_CRT).stat().st_size,
            },
            {
                "path": str(EXCEPTIONAL_TRANSPORT),
                "sha256": EXPECTED_HASHES[EXCEPTIONAL_TRANSPORT],
                "size_bytes": repository_path(EXCEPTIONAL_TRANSPORT).stat().st_size,
            },
            {
                "path": str(EXCEPTIONAL_QUOTIENT_BALL),
                "sha256": EXPECTED_HASHES[EXCEPTIONAL_QUOTIENT_BALL],
                "size_bytes": repository_path(EXCEPTIONAL_QUOTIENT_BALL).stat().st_size,
            },
            {
                "path": str(BIDEGREE21_P13_R20E1_PILOT),
                "sha256": EXPECTED_HASHES[BIDEGREE21_P13_R20E1_PILOT],
                "size_bytes": repository_path(BIDEGREE21_P13_R20E1_PILOT).stat().st_size,
            },
            {
                "path": str(BIDEGREE21_ALL80),
                "sha256": EXPECTED_HASHES[BIDEGREE21_ALL80],
                "size_bytes": repository_path(BIDEGREE21_ALL80).stat().st_size,
            },
            {
                "path": str(EXCEPTIONAL_PAIR_SIMULTANEOUS_H200000),
                "sha256": EXPECTED_HASHES[
                    EXCEPTIONAL_PAIR_SIMULTANEOUS_H200000
                ],
                "size_bytes": repository_path(
                    EXCEPTIONAL_PAIR_SIMULTANEOUS_H200000
                ).stat().st_size,
            },
            {
                "path": str(BIDEGREE21_P13_R20E1_NONLINEAR_POINTS_H1024),
                "sha256": EXPECTED_HASHES[
                    BIDEGREE21_P13_R20E1_NONLINEAR_POINTS_H1024
                ],
                "size_bytes": repository_path(
                    BIDEGREE21_P13_R20E1_NONLINEAR_POINTS_H1024
                ).stat().st_size,
            },
        )
    )
    unique: dict[str, dict[str, Any]] = {item["path"]: item for item in paths}
    fixed_fiber_markers = (
        "benchmark",
        "published_pair",
        "rank22_accidental",
        "rank22_auxiliary",
        "rank22_missing",
        "rank22_points",
        "rank22_record_group",
        "rank_certificates_v1",
    )
    entries = []
    for path, inventory in sorted(unique.items()):
        name = Path(path).name
        entry: dict[str, Any] = {
            "path": path,
            "sha256": inventory["sha256"],
            "identity": "fermigier-mestre-v1:u=28917/20",
        }
        if path == str(IMPORTED_NEAR_MISS):
            entry.update(
                {
                    "scope_membership": "exact identity",
                    "decision": "promoted-to-exact-rank-at-least-20-anchor",
                    "selection_leakage": (
                        "the 20-point subset was selected from the bounded point cloud "
                        "using mod-5 finite-reduction rank"
                    ),
                    "evidence": {
                        "selected_count": near_miss["point_cloud"]["selected_count"],
                        "log_conductor": near_miss["global_curve"]["log_conductor"],
                    },
                }
            )
        elif path == str(LEGACY_EXPLICIT_FORMULA):
            entry.update(
                {
                    "scope_membership": "exact identity",
                    "decision": "conditional-fixed-fiber-closure-not-identity-rejection",
                    "selection_leakage": "none; post-discovery analytic diagnostic",
                    "evidence": {
                        "delta": "11/5",
                        "conservative_upper_under_grh": explicit_formula[
                            "delta_11_over_5"
                        ]["conservative_explicit_formula_upper"],
                        "root_number": 1,
                        "interpretation": (
                            "GRH plus parity gives analytic rank at most 20; BSD+GRH "
                            "would make the exact rank 20"
                        ),
                    },
                }
            )
        elif path == str(LEGACY_NEIGHBORHOOD_AUDIT):
            entry.update(
                {
                    "scope_membership": "dense raw scope contains exact identity",
                    "decision": "excluded-before-selection-as-prior-anchor",
                    "selection_leakage": "none; exclusion preceded all scores and calls",
                    "evidence": {
                        "dense_excluded_parameters": audit["new_raw_scopes"][
                            "dense_window"
                        ]["excluded_parameters"],
                        "conductor_calls": audit["search_boundary"]["conductor_calls"],
                        "point_or_rank_calls": audit["search_boundary"][
                            "point_or_rank_calls"
                        ],
                    },
                }
            )
        elif path == str(LEGACY_HIGH_POWER_CRT):
            entry.update(
                {
                    "scope_membership": "excluded through the frozen prior set",
                    "decision": "not-searched-as-a-fresh-CRT-Gauss-candidate",
                    "selection_leakage": "none for the anchor; it was a prior exclusion",
                    "evidence": {
                        "fresh_population_count": high_power["population"]["fresh_count"],
                        "point_search_calls": high_power["outcome"]["point_search_calls"],
                        "rank_calls": high_power["outcome"]["rank_calls"],
                    },
                }
            )
        elif path == str(EXCEPTIONAL_TRANSPORT):
            entry.update(
                {
                    "scope_membership": "exact rank-20 anchor included",
                    "decision": "no-new-section-or-specialization-found",
                    "selection_leakage": (
                        "none; exhaustive exact structural classification in each "
                        "declared exceptional transport population"
                    ),
                    "evidence": {
                        "affine_pencils": sum(
                            item["count"]
                            for item in exceptional_transport["transport"]["affine"][
                                "histogram"
                            ]
                        ),
                        "quadratic_pencils": exceptional_transport["transport"][
                            "quadratic"
                        ]["histogram"][0]["count"],
                        "projective_mobius_pencils": exceptional_transport["transport"][
                            "mobius"
                        ]["pair_count"],
                        "true_fiber_products": exceptional_transport["fiber_products"][
                            "pair_count"
                        ],
                        "new_sections": exceptional_transport["outcome"]["new_sections"],
                        "new_specializations": exceptional_transport["outcome"][
                            "new_specializations"
                        ],
                    },
                }
            )
        elif path == str(EXCEPTIONAL_QUOTIENT_BALL):
            rank20_ball = exceptional_quotient_ball["direction_balls"]["rank20"]
            e22_ball = exceptional_quotient_ball["direction_balls"]["E22"]
            entry.update(
                {
                    "scope_membership": (
                        "exact rank-20 anchor included; complete signed exceptional-"
                        "quotient coefficient ball over {-1,0,1} with support weight "
                        "at most 2; no global-sign quotient"
                    ),
                    "decision": "no-new-low-genus-transport-in-support-at-most-2-ball",
                    "selection_leakage": (
                        "none inside the declared finite support<=2 populations; "
                        "directions of support>=3 are outside scope"
                    ),
                    "evidence": {
                        "coefficient_alphabet": rank20_ball["coefficient_alphabet"],
                        "maximum_support_weight": rank20_ball[
                            "maximum_support_weight"
                        ],
                        "global_sign_quotient": rank20_ball[
                            "global_sign_quotient"
                        ],
                        "signed_rank20_directions": len(rank20_ball["records"]),
                        "signed_E22_directions": len(e22_ball["records"]),
                        "cross_anchor_affine_interpolants": exceptional_quotient_ball[
                            "outcome"
                        ]["cross_anchor_affine_interpolants"],
                        "new_base_changes": exceptional_quotient_ball["outcome"][
                            "new_base_changes"
                        ],
                        "new_sections": exceptional_quotient_ball["outcome"][
                            "new_sections"
                        ],
                        "new_specializations": exceptional_quotient_ball["outcome"][
                            "new_specializations"
                        ],
                    },
                }
            )
        elif path == str(BIDEGREE21_P13_R20E1_PILOT):
            pilot_scope = bidegree21_pilot["scope"]
            entry.update(
                {
                    "scope_membership": (
                        "exact rank-20 anchor included in the single-pair "
                        "P13 x R20E1 bidegree-(2,1) finite-denominator-chart pilot"
                    ),
                    "decision": "no-new-section-or-specialization-in-single-pair-pilot",
                    "selection_leakage": (
                        "single representative pair only; no conclusion is transferred "
                        "to the other 79 possible independent pairs"
                    ),
                    "evidence": {
                        "completed_pairs": pilot_scope["completed_pairs"],
                        "completed_pair_count": pilot_scope["completed_pair_count"],
                        "possible_independent_pair_count": pilot_scope[
                            "possible_independent_pair_count"
                        ],
                        "all_pairs_classified": pilot_scope["all_pairs_classified"],
                        "nonlinear_rational_special_points_classified": False,
                        "genus_at_most_one_components": bidegree21_pilot["outcome"][
                            "genus_at_most_one_components"
                        ],
                        "new_sections": bidegree21_pilot["outcome"]["new_sections"],
                        "new_specializations": bidegree21_pilot["outcome"][
                            "new_specializations"
                        ],
                    },
                }
            )
        elif path == str(BIDEGREE21_ALL80):
            all80_scope = bidegree21_all80["scope"]
            entry.update(
                {
                    "scope_membership": (
                        "exact rank-20 anchor included in all 80 independent "
                        "E22 x rank20 bidegree-(2,1) pencils in the finite d=1 "
                        "denominator chart"
                    ),
                    "decision": (
                        "no-new-section-or-specialization-in-all-80-finite-"
                        "chart-classification"
                    ),
                    "selection_leakage": (
                        "none inside the declared 10 x 8 exact structural "
                        "population; the irreducible degree-32 components and "
                        "component intersections received no rational-point "
                        "search, and the infinity chart was not repeated"
                    ),
                    "evidence": {
                        "finite_denominator_chart": all80_scope[
                            "finite_denominator_chart"
                        ],
                        "all_80_independent_pairs_classified": all80_scope[
                            "all_80_independent_pairs_classified"
                        ],
                        "pair_count": bidegree21_all80["population"][
                            "pair_count"
                        ],
                        "irreducible_degree_32_pair_count": bidegree21_all80[
                            "residual_irreducibility"
                        ]["irreducible_pair_count"],
                        "unresolved_pair_count": len(
                            bidegree21_all80["residual_irreducibility"][
                                "unresolved_pairs"
                            ]
                        ),
                        "valid_genus_at_most_one_components": bidegree21_all80[
                            "outcome"
                        ]["valid_genus_at_most_one_components"],
                        "new_sections": bidegree21_all80["outcome"][
                            "new_sections"
                        ],
                        "new_specializations": bidegree21_all80["outcome"][
                            "new_specializations"
                        ],
                        "not_claimed": all80_scope["not_claimed"],
                    },
                }
            )
        elif path == str(EXCEPTIONAL_PAIR_SIMULTANEOUS_H200000):
            simultaneous_box = simultaneous_h200000["search_box"]
            simultaneous_outcome = simultaneous_h200000["outcome"]
            entry.update(
                {
                    "scope_membership": (
                        "exact rank-20 anchor included as a calibration parameter; "
                        "all 80 exceptional directions and all 3160 direction "
                        "pairs were checked in the declared projective-height "
                        "H<=200000 literal-shift box"
                    ),
                    "decision": (
                        "no-third-parameter-in-genuine-simultaneous-square-"
                        "H200000-box"
                    ),
                    "selection_leakage": (
                        "none inside the one-pass exact finite box; the genuine "
                        "simultaneous-square condition was used rather than a "
                        "product-square surrogate, and no absence is asserted "
                        "outside H=200000"
                    ),
                    "evidence": {
                        "coordinate": simultaneous_box["coordinate"],
                        "box_definition": simultaneous_box["definition"],
                        "projective_height_bound": simultaneous_box[
                            "projective_height_bound"
                        ],
                        "one_pass": simultaneous_box["one_pass"],
                        "retries": simultaneous_box["retries"],
                        "all_80_direction_searches_completed": simultaneous_outcome[
                            "all_80_direction_searches_completed"
                        ],
                        "direction_count": simultaneous_outcome[
                            "direction_count"
                        ],
                        "fiber_product_pair_count": simultaneous_outcome[
                            "fiber_product_pair_count"
                        ],
                        "product_square_surrogate_used": simultaneous_outcome[
                            "product_square_surrogate_used"
                        ],
                        "new_third_parameter_count": simultaneous_outcome[
                            "new_third_parameter_count"
                        ],
                        "rank_conductor_target_hits": simultaneous_outcome[
                            "rank_conductor_target_hits"
                        ],
                        "outside_box_absence_proved": False,
                    },
                }
            )
        elif path == str(BIDEGREE21_P13_R20E1_NONLINEAR_POINTS_H1024):
            nonlinear_scope = nonlinear_points_h1024["scope"]
            nonlinear_region = nonlinear_points_h1024["search_region"]
            entry.update(
                {
                    "scope_membership": (
                        "exact rank-20 anchor included in the single P13 x R20E1 "
                        "irreducible degree-32 component; complete primitive "
                        "projective-height H<=1024 affine box, with exact global "
                        "boundary and known-line intersection checks"
                    ),
                    "decision": (
                        "no-rational-point-on-single-degree32-component-in-"
                        "H1024-box"
                    ),
                    "selection_leakage": (
                        "single representative pair and one exact finite box only; "
                        "no conclusion is transferred outside H=1024 or to the "
                        "other 79 exceptional pairs"
                    ),
                    "evidence": {
                        "representative_pair": nonlinear_points_h1024[
                            "representative_pair"
                        ],
                        "box_definition": nonlinear_region["definition"],
                        "projective_height_bound": nonlinear_region["height"],
                        "affine_box_complete": nonlinear_scope[
                            "affine_box_complete"
                        ],
                        "C_D_pairs_scanned": nonlinear_region[
                            "C_D_pairs_scanned"
                        ],
                        "crt_root_combinations_tested": nonlinear_region[
                            "crt_root_combinations_tested"
                        ],
                        "primitive_exact_candidates": nonlinear_region[
                            "primitive_exact_candidates"
                        ],
                        "exact_homogeneous_evaluations": nonlinear_region[
                            "exact_homogeneous_evaluations"
                        ],
                        "exact_hits": nonlinear_region["exact_hits"],
                        "rational_points_at_projective_infinity": nonlinear_points_h1024[
                            "outcome"
                        ]["rational_points_at_projective_infinity"],
                        "distinct_rational_known_line_intersections": nonlinear_points_h1024[
                            "outcome"
                        ]["distinct_rational_known_line_intersections"],
                        "new_sections": nonlinear_points_h1024["outcome"][
                            "new_sections"
                        ],
                        "new_specializations": nonlinear_points_h1024["outcome"][
                            "new_specializations"
                        ],
                        "all_rational_points_on_degree32_component_classified": nonlinear_scope[
                            "all_rational_points_on_degree32_component_classified"
                        ],
                        "other_exceptional_pairs_classified": nonlinear_scope[
                            "other_exceptional_pairs_classified"
                        ],
                        "not_claimed": nonlinear_scope["not_claimed"],
                    },
                }
            )
        elif name == "elliptic_fermigier_global.json":
            global_data = json.loads(repository_path(path).read_text())
            entry.update(
                {
                    "scope_membership": (
                        "included exactly as primitive literal shift T=28917/10 "
                        "in the exhaustive 0<=a<=100000,1<=b<=1000 box"
                    ),
                    "decision": "enumerated-but-not-retained-in-discovery-union",
                    "selection_leakage": (
                        "the retained union used rank, power, composite, and per-"
                        "denominator frontier selection"
                    ),
                    "evidence": {
                        "primitive_pairs_enumerated": global_data["population"][
                            "primitive_pairs_enumerated"
                        ],
                        "retained_union_count": global_data["population"][
                            "retained_union_count"
                        ],
                        "anchor_absent_from_recorded_retained_parameters": True,
                    },
                }
            )
        elif any(marker in name for marker in fixed_fiber_markers):
            entry.update(
                {
                    "scope_membership": "different fixed identity u=19754/39",
                    "decision": "not-an-anchor-experiment",
                    "selection_leakage": "not applicable to this identity",
                    "evidence": (
                        "audit found no retained exact anchor parameter in this artifact"
                    ),
                }
            )
        elif name == "fermigier_crt_seed_v1.json":
            entry.update(
                {
                    "scope_membership": "different seed identity u=673709/29965",
                    "decision": "not-an-anchor-experiment",
                    "selection_leakage": "not applicable to this identity",
                    "evidence": "audit found no exact anchor parameter",
                }
            )
        else:
            entry.update(
                {
                    "scope_membership": "not established by recorded retained parameters",
                    "decision": "anchor-absent-from-recorded-parameter-set",
                    "selection_leakage": "unknown or artifact-specific",
                    "evidence": (
                        "the exact audit extracted no u=28917/20 source from this artifact; "
                        "this does not prove absence from an unstored raw population"
                    ),
                }
            )
        entries.append(entry)
    digest_input = [
        {"path": entry["path"], "sha256": entry["sha256"]} for entry in entries
    ]
    return {
        "coverage": (
            "every Fermigier JSON/JSONL experiment in the frozen 53-file audit, "
            "plus the audit itself, the later disjoint high-power CRT lane, "
            "the complete exceptional-transport classification, the exact signed "
            "support<=2 quotient balls, the single-pair P13 x R20E1 pilot, the "
            "all-80 finite-chart classification, the genuine simultaneous-square "
            "H<=200000 box, and the single-pair degree-32 H<=1024 point box"
        ),
        "entry_count": len(entries),
        "path_and_sha256_digest": stable_json_sha256(digest_input),
        "absence_rule": (
            "absence from retained exact parameter records is never upgraded to a "
            "claim of absence from an unstored raw population"
        ),
        "entries": entries,
    }


def source_reference(relative: Path | str) -> dict[str, Any]:
    path = repository_path(relative)
    return {
        "path": str(relative),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def replay_imported_verifier() -> dict[str, Any]:
    """Keep the imported cyclic-log verifier as a capped independent replay."""

    relative = Path("elliptic-curves/scripts/verify_fermigier_rank20_near_miss.py")
    path = repository_path(relative)
    timeout_seconds = 30.0
    try:
        completed = subprocess.run(
            [sys.executable, str(path)],
            cwd=REPOSITORY_ROOT,
            text=True,
            capture_output=True,
            check=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"the imported verifier exceeded {timeout_seconds}s; no retry"
        ) from error
    if completed.stderr.strip() or not completed.stdout.startswith(
        "PASS Fermigier near miss:"
    ):
        raise AssertionError("the imported independent verifier did not pass cleanly")
    return {
        **source_reference(relative),
        "status": "passed",
        "stdout": completed.stdout.strip(),
        "timeout_seconds": timeout_seconds,
        "retried": False,
        "implementation": "ecsearch cyclic finite-group discrete logarithms",
    }


def build_record() -> dict[str, Any]:
    near_miss = checked_json(IMPORTED_NEAR_MISS)
    rank_certificates = checked_json(IMPORTED_RANK_CERTIFICATES)
    generic_rank_theorem = checked_json(GENERIC_RANK_THEOREM)
    exceptional_transport = checked_json(EXCEPTIONAL_TRANSPORT)
    exceptional_quotient_ball = checked_json(EXCEPTIONAL_QUOTIENT_BALL)
    bidegree21_pilot = checked_json(BIDEGREE21_P13_R20E1_PILOT)
    bidegree21_all80 = checked_json(BIDEGREE21_ALL80)
    simultaneous_h200000 = checked_json(EXCEPTIONAL_PAIR_SIMULTANEOUS_H200000)
    nonlinear_points_h1024 = checked_json(
        BIDEGREE21_P13_R20E1_NONLINEAR_POINTS_H1024
    )
    explicit_formula = checked_json(LEGACY_EXPLICIT_FORMULA)
    audit = checked_json(LEGACY_NEIGHBORHOOD_AUDIT)
    high_power = checked_json(LEGACY_HIGH_POWER_CRT)
    for relative, expected_hash in EXPECTED_STRUCTURAL_SCRIPT_HASHES.items():
        actual_hash = sha256_file(repository_path(relative))
        if actual_hash != expected_hash:
            raise AssertionError(
                f"frozen structural script {relative} changed: "
                f"{actual_hash} != {expected_hash}"
            )
    imported_verifier_replay = replay_imported_verifier()
    if Q(near_miss["family"]["adapter_parameter"]) != ANCHOR_U:
        raise AssertionError("the imported near miss changed canonical identity")
    if Q(near_miss["family"]["literal_shift"]) != LITERAL_SHIFT:
        raise AssertionError("the imported literal-shift alias changed")

    specialization = specialize_fermigier_rank_sections(ANCHOR_U)
    canonical_model = specialization.canonical_model
    if canonical_model != fermigier_canonical_coefficients(ANCHOR_U):
        raise AssertionError("the canonical specialization changed")
    abscissas = tuple(Q(value) for value in near_miss["bounded_search"]["abscissas"])
    raw_ratpoints = canonical_ratpoints_output(abscissas)
    if hashlib.sha256(raw_ratpoints.encode()).hexdigest() != near_miss[
        "bounded_search"
    ]["canonical_output_sha256"]:
        raise AssertionError("the bounded-search abscissa stream changed")
    searched_quartic_points = parse_ratpoints_output(
        specialization.quartic_model, raw_ratpoints
    )
    cloud = section_and_point_cloud_differences(
        specialization, searched_quartic_points
    )
    replay_cloud, pool_provenance = reconstruct_pool(
        specialization, searched_quartic_points
    )
    if replay_cloud != cloud or len(cloud) != 115 or len(abscissas) != 58:
        raise AssertionError("the complete 58-abscissa/115-difference pool changed")

    selected_indices = tuple(near_miss["point_cloud"]["selected_indices"])
    selected_canonical = tuple(cloud[index] for index in selected_indices)
    if len(selected_canonical) != 20 or not set(range(12)).issubset(selected_indices):
        raise AssertionError("the imported selected basis changed")

    canonical_invariants = weierstrass_invariants(canonical_model)
    b2 = canonical_invariants["b2"]
    c4 = canonical_invariants["c4"]
    c6 = canonical_invariants["c6"]
    canonical_to_legacy = WeierstrassChange(
        Q(1, 6),
        -b2 / 12,
        -canonical_model[0] / 2,
        canonical_model[0] * b2 / 24 - canonical_model[2] / 2,
    )
    legacy_short_model = change_weierstrass_model(
        canonical_model, canonical_to_legacy
    )
    expected_legacy = FermigierMestreFamily.coefficients(LITERAL_SHIFT)
    if legacy_short_model != expected_legacy:
        raise AssertionError("the canonical/legacy short-model bridge changed")
    if legacy_short_model != (Q(0), Q(0), Q(0), -27 * c4, -54 * c6):
        raise AssertionError("the canonical standard-short invariant bridge failed")

    quartic = tuple(specialization.quartic_model.quartic)
    invariant_i, invariant_j = binary_quartic_invariants_low_to_high(quartic)
    raw_quartic_jacobian = (Q(0), Q(0), Q(0), -27 * invariant_i, -27 * invariant_j)
    scale = 101_232 * ANCHOR_U
    if invariant_i != scale**4 * c4 or invariant_j != 2 * scale**6 * c6:
        raise AssertionError("the binary-quartic/canonical invariant bridge failed")
    legacy_to_raw = WeierstrassChange(Q(1, 1) / scale, Q(0), Q(0), Q(0))
    if change_weierstrass_model(legacy_short_model, legacy_to_raw) != raw_quartic_jacobian:
        raise AssertionError("the normalized/raw binary-quartic scale changed")

    selected_legacy = tuple(
        source_point_to_target(point, canonical_to_legacy)
        for point in selected_canonical
    )
    gp = run_capped_gp_cross_check(
        canonical_model, legacy_short_model, selected_legacy
    )
    imported_minimal_model = tuple(Q(value) for value in near_miss["global_curve"]["minimal_model"])
    if gp["minimal_model"] != imported_minimal_model:
        raise AssertionError("the independently reduced minimal model changed")
    if gp["minimal_change"] != EXPECTED_MINIMAL_CHANGE:
        raise AssertionError("the exact canonical/minimal change changed")
    if change_weierstrass_model(canonical_model, EXPECTED_MINIMAL_CHANGE) != imported_minimal_model:
        raise AssertionError("the exact canonical/minimal transformation failed")
    for field, expected in (
        ("conductor", int(near_miss["global_curve"]["conductor"])),
        ("minimal_discriminant", int(near_miss["global_curve"]["minimal_discriminant"])),
        ("root_number", int(near_miss["global_curve"]["root_number"])),
    ):
        if gp[field] != expected:
            raise AssertionError(f"the independent PARI {field} replay changed")

    canonical_to_raw = WeierstrassChange(
        canonical_to_legacy.u * legacy_to_raw.u,
        canonical_to_legacy.r,
        canonical_to_legacy.s,
        canonical_to_legacy.t,
    )
    if change_weierstrass_model(canonical_model, canonical_to_raw) != raw_quartic_jacobian:
        raise AssertionError("the composed canonical/raw change failed")

    transported_cloud = {
        "canonical_generalized": cloud,
        "legacy_normalized_short": tuple(
            source_point_to_target(point, canonical_to_legacy) for point in cloud
        ),
        "raw_binary_quartic_jacobian": tuple(
            source_point_to_target(point, canonical_to_raw) for point in cloud
        ),
        "global_minimal": tuple(
            source_point_to_target(point, EXPECTED_MINIMAL_CHANGE) for point in cloud
        ),
    }
    models = {
        "canonical_generalized": canonical_model,
        "legacy_normalized_short": legacy_short_model,
        "raw_binary_quartic_jacobian": raw_quartic_jacobian,
        "global_minimal": imported_minimal_model,
    }
    for label, points in transported_cloud.items():
        if len(points) != 115 or any(
            not is_on_weierstrass_curve(models[label], point) for point in points
        ):
            raise AssertionError(f"the complete pool failed transport to {label}")
    selected_transports = {
        label: tuple(points[index] for index in selected_indices)
        for label, points in transported_cloud.items()
    }
    for label, points in selected_transports.items():
        if tuple(points) != tuple(
            transported_cloud[label][index] for index in selected_indices
        ):
            raise AssertionError("selected-point transport order changed")

    saturated_legacy = gp["saturated_legacy_points"]
    if len(saturated_legacy) != 20 or any(
        not is_on_weierstrass_curve(legacy_short_model, point)
        for point in saturated_legacy
    ):
        raise AssertionError("PARI returned an invalid bounded saturation candidate")
    saturated_canonical = tuple(
        target_point_to_source(point, canonical_to_legacy)
        for point in saturated_legacy
    )
    saturated_minimal = tuple(
        source_point_to_target(point, EXPECTED_MINIMAL_CHANGE)
        for point in saturated_canonical
    )
    saturated_raw = tuple(
        source_point_to_target(point, canonical_to_raw)
        for point in saturated_canonical
    )
    for label, points in (
        ("canonical_generalized", saturated_canonical),
        ("legacy_normalized_short", saturated_legacy),
        ("raw_binary_quartic_jacobian", saturated_raw),
        ("global_minimal", saturated_minimal),
    ):
        if any(not is_on_weierstrass_curve(models[label], point) for point in points):
            raise AssertionError(f"a saturation candidate failed transport to {label}")

    original_certificates = {
        f"mod_{relation_prime}": build_finite_quotient_certificate(
            legacy_short_model,
            selected_legacy,
            relation_prime=relation_prime,
            prime_bound=CERTIFICATE_PRIME_BOUND,
        )
        for relation_prime in (2, 3, 5)
    }
    saturated_certificates = {
        f"mod_{relation_prime}": build_finite_quotient_certificate(
            legacy_short_model,
            saturated_legacy,
            relation_prime=relation_prime,
            prime_bound=CERTIFICATE_PRIME_BOUND,
        )
        for relation_prime in (2, 3, 5)
    }
    for certificates, points in (
        (original_certificates, selected_legacy),
        (saturated_certificates, saturated_legacy),
    ):
        for certificate in certificates.values():
            verify_finite_quotient_certificate(
                legacy_short_model, points, certificate
            )
    expected_original_ranks = {"mod_2": 0, "mod_3": 19, "mod_5": 20}
    actual_original_ranks = {
        key: value["combined_rank_over_relation_field"]
        for key, value in original_certificates.items()
    }
    if actual_original_ranks != expected_original_ranks:
        raise AssertionError(
            "the imported subgroup's cross-modulus profile changed: "
            f"{actual_original_ranks!r}"
        )
    if any(
        certificate["combined_rank_over_relation_field"] != 20
        or not certificate["certified_independent"]
        for certificate in saturated_certificates.values()
    ):
        raise AssertionError("the saturated candidate lost a full cross-certificate")

    direct_canonical_images = tuple(
        quartic_point_to_canonical_point(specialization.quartic_model, point)
        for point in searched_quartic_points
    )
    if any(
        point[1] ** 2
        != evaluate_polynomial(specialization.quartic_model.quartic, point[0])
        for point in searched_quartic_points
    ):
        raise AssertionError("a reconstructed bounded-search quartic point failed")
    if any(
        not is_on_weierstrass_curve(canonical_model, point)
        for point in direct_canonical_images
    ):
        raise AssertionError("a bounded quartic point mapped off the canonical curve")

    pool_entries = []
    for index, (point, provenance) in enumerate(zip(cloud, pool_provenance)):
        pool_entries.append(
            {
                "index_zero_based": index,
                "canonical_point": point_record(point),
                "provenance": provenance,
                "selected_in_imported_basis": index in selected_indices,
            }
        )
    selected_entries = []
    for selected_position, pool_index in enumerate(selected_indices):
        selected_entries.append(
            {
                "basis_index_zero_based": selected_position,
                "pool_index_zero_based": pool_index,
                "points": {
                    label: point_record(points[selected_position])
                    for label, points in selected_transports.items()
                },
            }
        )

    generic_point_labels = []
    for root in (0, 55, 314, 378, 1007, 1036):
        for sign in (-1, 1):
            generic_point_labels.append(f"Q_root_{root}_{'minus' if sign < 0 else 'plus'}_s")
    generic_point_labels.append("Q_extra_1256_over_5_minus_17s_over_35")

    artifact: dict[str, Any] = {
        "schema": CANDIDATE_SCHEMA,
        "claim_level": "exact-identity-model-transforms-rank-lower-bound-and-global-data",
        "identity": {
            **canonical_candidate_identity(
                FAMILY_ID, "u", ANCHOR_U, sign_quotient=True
            ),
            "coordinate_definition": "canonical adapter u=s/2",
            "aliases": [
                {
                    "coordinate": "adapter_u",
                    "value": fraction_text(ANCHOR_U),
                    "maps_to_canonical_u": fraction_text(ANCHOR_U),
                    "role": "canonical",
                },
                {
                    "coordinate": "adapter_u",
                    "value": fraction_text(-ANCHOR_U),
                    "maps_to_canonical_u": fraction_text(ANCHOR_U),
                    "role": "sign-quotient alias",
                },
                {
                    "coordinate": "literal_symmetric_shift_s",
                    "value": fraction_text(LITERAL_SHIFT),
                    "map": "u=abs(s)/2",
                    "maps_to_canonical_u": fraction_text(ANCHOR_U),
                    "role": "legacy literal-shift alias only",
                },
                {
                    "coordinate": "literal_symmetric_shift_s",
                    "value": fraction_text(-LITERAL_SHIFT),
                    "map": "u=abs(s)/2",
                    "maps_to_canonical_u": fraction_text(ANCHOR_U),
                    "role": "legacy sign alias only",
                },
            ],
            "normalization_warning": (
                "the historical factor-two discrepancy remains unresolved; aliases "
                "do not create distinct candidate identities"
            ),
        },
        "artifact_namespace_bridge": {
            "imported_ecsearch_namespace": [
                source_reference(IMPORTED_NEAR_MISS),
                source_reference(IMPORTED_RANK_CERTIFICATES),
            ],
            "legacy_generated_results_namespace": [
                source_reference(GENERIC_RANK_THEOREM),
                source_reference(EXCEPTIONAL_TRANSPORT),
                source_reference(EXCEPTIONAL_QUOTIENT_BALL),
                source_reference(BIDEGREE21_P13_R20E1_PILOT),
                source_reference(BIDEGREE21_ALL80),
                source_reference(EXCEPTIONAL_PAIR_SIMULTANEOUS_H200000),
                source_reference(BIDEGREE21_P13_R20E1_NONLINEAR_POINTS_H1024),
                source_reference(LEGACY_EXPLICIT_FORMULA),
                source_reference(LEGACY_NEIGHBORHOOD_AUDIT),
                source_reference(LEGACY_HIGH_POWER_CRT),
            ],
            "new_structural_experiment_implementations": [
                source_reference(relative)
                for relative in EXPECTED_STRUCTURAL_SCRIPT_HASHES
            ],
            "files_moved_or_rewritten": [],
        },
        "models": {
            "canonical_generalized": {
                "coefficients_a1_a2_a3_a4_a6": model_record(canonical_model),
                "discriminant": fraction_text(canonical_invariants["discriminant"]),
                "c4": fraction_text(c4),
                "c6": fraction_text(c6),
            },
            "binary_quartic": {
                "equation": "z^2=e+d*x+c*x^2+b*x^3+a*x^4",
                "coefficients_low_to_high_e_d_c_b_a": [
                    fraction_text(value) for value in quartic
                ],
                "invariant_I": fraction_text(invariant_i),
                "invariant_J": fraction_text(invariant_j),
                "raw_jacobian_equation": "Y^2=X^3-27*I*X-27*J",
                "raw_jacobian_coefficients": model_record(raw_quartic_jacobian),
            },
            "legacy_normalized_short_jacobian": {
                "coefficients": model_record(legacy_short_model),
                "equation": "y^2=x^3-27*c4*x-54*c6",
                "raw_quartic_scale_101232u": fraction_text(scale),
            },
            "global_minimal": {
                "coefficients": model_record(imported_minimal_model),
                "minimal_discriminant": str(gp["minimal_discriminant"]),
                "conductor": str(gp["conductor"]),
                "log_conductor": near_miss["global_curve"]["log_conductor"],
                "root_number": gp["root_number"],
                "strict_log_target": "182.72",
                "below_strict_log_target": True,
            },
        },
        "exact_transformations": {
            "convention": (
                "[u,r,s,t] means x_source=u^2*x_target+r and "
                "y_source=u^3*y_target+s*u^2*x_target+t"
            ),
            "canonical_to_legacy_normalized_short": {
                "change_u_r_s_t": canonical_to_legacy.to_record(),
                "forward_coordinate_formula": {
                    "x_legacy": "36*x_canonical+3*b2",
                    "y_legacy": "108*(2*y_canonical+a1*x_canonical+a3)",
                },
                "exact_model_replay_checked": True,
            },
            "legacy_normalized_short_to_raw_binary_quartic_jacobian": {
                "change_u_r_s_t": legacy_to_raw.to_record(),
                "forward_coordinate_formula": {
                    "x_raw": "(101232*u)^2*x_legacy",
                    "y_raw": "(101232*u)^3*y_legacy",
                },
                "exact_model_replay_checked": True,
            },
            "canonical_to_raw_binary_quartic_jacobian": {
                "change_u_r_s_t": canonical_to_raw.to_record(),
                "forward_coordinate_formula": {
                    "x_raw": "(101232*u)^2*(36*x_canonical+3*b2)",
                    "y_raw": (
                        "108*(101232*u)^3*(2*y_canonical+"
                        "a1*x_canonical+a3)"
                    ),
                },
                "exact_model_replay_checked": True,
            },
            "canonical_to_global_minimal": {
                "change_u_r_s_t": EXPECTED_MINIMAL_CHANGE.to_record(),
                "exact_model_replay_checked": True,
                "independently_returned_by_pari": True,
            },
        },
        "generic_twelve_section_basis": {
            "generic_definition": "D_i=canonical_image(Q_i)-canonical_image(Q_0), i=1,...,12",
            "quartic_point_order": generic_point_labels,
            "specialized_at_canonical_u": fraction_text(ANCHOR_U),
            "specialized_canonical_basis": [
                point_record(point) for point in specialization.section_differences
            ],
            "specialized_basis_sha256": point_sequence_sha256(
                specialization.section_differences
            ),
            "generic_independence_reference": {
                "path": str(IMPORTED_RANK_CERTIFICATES),
                "sha256": EXPECTED_HASHES[IMPORTED_RANK_CERTIFICATES],
                "certificate_at_u": rank_certificates["generic_sections"][
                    "adapter_parameter"
                ],
                "argument": (
                    "one rank-12 specialization certificate proves generic "
                    "independence of these twelve section differences"
                ),
            },
            "exact_generic_rank_theorem": {
                "path": str(GENERIC_RANK_THEOREM),
                "sha256": EXPECTED_HASHES[GENERIC_RANK_THEOREM],
                "verifier": source_reference(
                    "elliptic-curves/cas/verify_fermigier_generic_rank_exact.py"
                ),
                "arithmetic_generic_rank_over_Q_of_u": generic_rank_theorem[
                    "theorem"
                ]["arithmetic_generic_Mordell_Weil_rank_over_Q_of_u"],
                "status": generic_rank_theorem["theorem"][
                    "arithmetic_rank_status"
                ],
                "geometric_generic_rank_interval_over_Qbar_of_u": generic_rank_theorem[
                    "theorem"
                ]["geometric_generic_Mordell_Weil_rank_interval_over_Qbar_of_u"],
            },
        },
        "complete_point_pool": {
            "bounded_search_is_not_complete": True,
            "height_bound": near_miss["bounded_search"]["height_bound"],
            "denominator_bound": near_miss["bounded_search"]["denominator_bound"],
            "abscissa_count": len(abscissas),
            "abscissas": [fraction_text(value) for value in abscissas],
            "signed_quartic_point_count": len(searched_quartic_points),
            "signed_quartic_points_sha256": point_sequence_sha256(
                searched_quartic_points
            ),
            "direct_canonical_images_sha256": point_sequence_sha256(
                direct_canonical_images
            ),
            "deduplicated_difference_count": len(cloud),
            "difference_pool": pool_entries,
            "transport_replay": {
                label: {
                    "point_count": len(points),
                    "point_sequence_sha256": point_sequence_sha256(points),
                    "all_exact_curve_memberships_checked": True,
                }
                for label, points in transported_cloud.items()
            },
        },
        "imported_selected_twenty_basis": {
            "selected_indices_zero_based": list(selected_indices),
            "basis_count": len(selected_entries),
            "basis": selected_entries,
            "canonical_point_sequence_sha256": point_sequence_sha256(
                selected_canonical
            ),
            "imported_ecsearch_cyclic_log_mod5_certificate": near_miss[
                "point_cloud"
            ]["certificate"],
            "imported_verifier": imported_verifier_replay,
            "claim": "the imported subgroup has exact rank 20",
        },
        "bounded_saturation_status": {
            "status": "not-globally-proved-saturated",
            "engine": "PARI/GP ellsaturation",
            "pari_version": gp["pari_version"],
            "prime_bound_strict_upper_limit": SATURATION_PRIME_BOUND,
            "input_basis_count": len(selected_legacy),
            "input_legacy_point_sequence_sha256": point_sequence_sha256(
                selected_legacy
            ),
            "returned_basis_count": len(saturated_legacy),
            "returned_basis_changed": saturated_legacy != selected_legacy,
            "returned_legacy_point_sequence_sha256": point_sequence_sha256(
                saturated_legacy
            ),
            "returned_legacy_basis": [point_record(point) for point in saturated_legacy],
            "returned_basis_all_exact_memberships_checked": True,
            "transported_basis_sha256": {
                "canonical_generalized": point_sequence_sha256(saturated_canonical),
                "raw_binary_quartic_jacobian": point_sequence_sha256(saturated_raw),
                "global_minimal": point_sequence_sha256(saturated_minimal),
            },
            "scope_warning": (
                "PARI documents ellsaturation under a finite-index hypothesis. "
                "No unconditional rank upper bound is known, so this is a bounded "
                "small-prime saturation candidate, not a saturation theorem."
            ),
        },
        "independent_cas_cross_certificates": {
            "implementation": source_reference(
                "elliptic-curves/cas/elliptic_candidate_record.py"
            ),
            "implementation_independence": (
                "finite quotients E(F_p)/ell E(F_p) are enumerated directly; "
                "the imported ecsearch certificate instead uses cyclic generators "
                "and discrete logarithms"
            ),
            "original_imported_basis": {
                "interpretation": (
                    "mod 2 and mod 3 expose subgroup-index obstructions; the new "
                    "general mod-l engine independently certifies the original basis at ell=5"
                ),
                **original_certificates,
            },
            "bounded_saturation_candidate_basis": {
                "interpretation": (
                    "the returned exact basis is independently certified at ell=2,3,5; "
                    "this proves its independence, not global saturation"
                ),
                **saturated_certificates,
            },
            "exact_rank_lower_bound": 20,
        },
        "global_arithmetic": {
            "conductor": str(gp["conductor"]),
            "log_conductor": near_miss["global_curve"]["log_conductor"],
            "root_number": gp["root_number"],
            "rank_lower_bound": 20,
            "unconditional_rank_upper_bound": None,
            "conditional_fixed_fiber_diagnostic": {
                "assumption": "GRH for the elliptic-curve L-function",
                "delta": "11/5",
                "conservative_upper": explicit_formula["delta_11_over_5"][
                    "conservative_explicit_formula_upper"
                ],
                "root_parity_used": "+1 implies even analytic rank",
                "conclusion_under_grh": "analytic rank at most 20",
                "conclusion_under_bsd_and_grh": "algebraic rank exactly 20",
            },
        },
        "promotion_and_rejection_ledger": experiment_ledger(
            audit,
            near_miss,
            explicit_formula,
            high_power,
            exceptional_transport,
            exceptional_quotient_ball,
            bidegree21_pilot,
            bidegree21_all80,
            simultaneous_h200000,
            nonlinear_points_h1024,
        ),
        "reproducibility": {
            "generator": "elliptic-curves/cas/build_fermigier_rank20_candidate_record.py",
            "generator_sha256": sha256_file(Path(__file__)),
            "canonical_command": (
                "PYTHONPATH=elliptic-curves/cas python3 "
                "elliptic-curves/cas/build_fermigier_rank20_candidate_record.py "
                "--output artifacts/generated-results/"
                "elliptic_curve_candidate_fermigier_mestre_v1_u28917_20.json"
            ),
            "pari_cross_check": {
                "program_sha256": gp["program_sha256"],
                "timeout_seconds": gp["timeout_seconds"],
                "stack_bytes": gp["stack_bytes"],
                "retried": gp["retried"],
                "point_search_calls": 0,
                "score_calls": 0,
                "parameter_search_calls": 0,
                "conductor_search_calls": 0,
            },
            "randomness": "none",
        },
        "limitations": {
            "rank": "no twenty-first independent point and no unconditional upper bound",
            "saturation": "no global saturation theorem",
            "bounded_point_pool": "absence of further points is not proved",
            "normalization": "historical factor-two discrepancy remains unresolved",
            "target_status": "rank-20 near miss, not a rank-21 breakthrough",
        },
    }
    validate_candidate_identity(artifact)
    artifact["result_sha256"] = stable_json_sha256(artifact)
    return artifact


def exclusive_write(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    with os.fdopen(descriptor, "w") as stream:
        json.dump(artifact, stream, indent=2, sort_keys=True)
        stream.write("\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    artifact = build_record()
    artifact["generated_at_utc"] = datetime.now(timezone.utc).isoformat()
    # Generation time is provenance, so keep the mathematical result digest
    # defined on the record before adding this volatile field.
    output = args.output
    if not output.is_absolute():
        output = REPOSITORY_ROOT / output
    exclusive_write(output, artifact)
    print(
        "PASS canonical Fermigier candidate "
        f"{artifact['identity']['candidate_key']}: rank>=20, "
        "all 115 pool points transported, independent mod-2/mod-3/mod-5 "
        "cross-certificates on the bounded saturation candidate"
    )


if __name__ == "__main__":
    main()
