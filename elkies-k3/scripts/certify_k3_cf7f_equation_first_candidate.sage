#!/usr/bin/env sage-python
"""Aggregate the determinant-714 MW1-to-MW16 equation-first promotion."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, matrix, pari


ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
ADAPTERS = [
    GEN / f"elkies-k3-k3-cf7f6c91a3a40d32-source-search-target-partner{i}-lattice-only-v1.json"
    for i in (1, 2)
]
SOURCES = [
    GEN / f"elkies-k3-k3-cf7f6c91a3a40d32-semistable-mw0-2-sources-large-a-partner{i}-v1.json"
    for i in (1, 2)
]
FIBRES = {
    5: GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-fibre-ansatz-mod5-v1.json",
    7: GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-fibre-ansatz-mod7-v1.json",
}
MARKINGS = {
    (5, "square"): GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-pole1-marking-mod5-square-v1.json",
    (5, "nonsquare"): GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-pole1-marking-mod5-nonsquare-v1.json",
    (7, "square"): GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-pole1-marking-mod7-square-v1.json",
    (7, "nonsquare"): GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-pole1-marking-mod7-nonsquare-v1.json",
}
HENSEL = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-marked-gf7-hensel-v1.json"
HENSEL_SECOND = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-marked-gf7-seed2-hensel-v1.json"
FORMAL = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-formal-smoothness-v1.json"
RATIONAL_SCAN = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-rational-parameter-scan-v1.json"
CORRIDOR = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-same-ns-compiler-routes-rankfirst-cap2000-v1.json"
DEFAULT_OUTPUT = GEN / "elkies-k3-k3-cf7f6c91a3a40d32-equation-first-candidate-v1.json"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

adapters = [load(path) for path in ADAPTERS]
source_payloads = [load(path) for path in SOURCES]
for index, adapter in enumerate(adapters, start=1):
    if (
        adapter["surface_id"] != "K3-cf7f6c91a3a40d32"
        or adapter["determinant"] != 714
        or adapter["frame"]["frame_id"] != "K3-cf7f6c91a3a40d32-F001"
        or adapter["frame"]["root_type"] != "A1"
        or adapter["equation_work_authorized"]
    ):
        raise ValueError(f"partner-{index} lattice-only adapter changed")

source_counts = []
for payload, expected_total, expected_ranks, expected_complete in zip(
    source_payloads,
    (463, 548),
    (Counter({2: 463}), Counter({2: 546, 1: 2})),
    (Counter({2: 186}), Counter({2: 216, 1: 2})),
):
    rows = payload["sources"]
    ranks = Counter(row["source"]["mw_rank_for_rho_19"] for row in rows)
    complete = Counter(
        row["source"]["mw_rank_for_rho_19"]
        for row in rows
        if bool(row["source"]["pole_audit"].get("basis_with_all_poles_at_most_two"))
    )
    if len(rows) != expected_total or ranks != expected_ranks or complete != expected_complete:
        raise ValueError("source census changed")
    source_counts.append(
        {
            "rows": len(rows),
            "mw_rank_histogram": {str(key): value for key, value in sorted(ranks.items())},
            "complete_basis_through_pole_two_by_mw_rank": {
                str(key): value for key, value in sorted(complete.items())
            },
        }
    )

mw1_rows = [
    row
    for row in source_payloads[1]["sources"]
    if row["source"]["mw_rank_for_rho_19"] == 1
]
if [row["source_id"] for row in mw1_rows] != [
    "K3-cf7f6c91a3a40d32-S0223",
    "K3-cf7f6c91a3a40d32-S0430",
]:
    raise ValueError("MW1 representatives changed")
for row in mw1_rows:
    source = row["source"]
    if (
        source["root_type"] != "A4+2A6"
        or source["torsion"] != 1
        or source["mw_height_gram"] != [["102/35"]]
        or [section["pole_order"] for section in source["pole_audit"]["basis"]] != [1]
    ):
        raise ValueError("MW1 source profile changed")
left, right = (matrix(ZZ, row["source"]["gram"]) for row in mw1_rows)
isometry = pari(left).qfisom(pari(right))
if isometry == 0:
    raise ArithmeticError("the two MW1 source Grams are no longer integrally isometric")
isometry_matrix = matrix(ZZ, isometry)
if abs(isometry_matrix.det()) != 1:
    raise ArithmeticError("MW1 isometry is not unimodular")

fibre_summary = {}
for prime, path in FIBRES.items():
    payload = load(path)
    expected_models = 6 if prime == 5 else 20
    if (
        payload["status"] != "PASS_EXACT_EXHAUSTIVE_MODULAR_SOURCE_FIBRE_ANSATZ"
        or not payload["scan"]["exhausted"]
        or payload["ansatz"]["normalized_reducible_supports"]
        != ["0:I5", "1:I7", "infinity:I7"]
        or len(payload["examples"]) != expected_models
    ):
        raise ValueError(f"GF({prime}) fibre census changed")
    fibre_summary[str(prime)] = {
        "squarefree_models": expected_models,
        "normalized_A_polynomials": payload["scan"]["normalized_A_polynomials"],
    }

marking_summary = []
for (prime, twist_class), path in MARKINGS.items():
    payload = load(path)
    sections = payload["accounting"]["marked_mw1_sections"]
    expected = 4 if (prime, twist_class) == (7, "nonsquare") else 0
    expected_status = (
        "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_WITH_MARKED_MW1_SECTION"
        if expected
        else "PASS_EXACT_EXHAUSTIVE_NORMALIZED_CHART_EMPTY_MARKED_MW1_LOCUS"
    )
    if (
        payload["status"] != expected_status
        or payload["source"]["component_depths_at_normalized_supports"] != [1, 2, 1]
        or sections != expected
    ):
        raise ValueError(f"marking result changed: {path}")
    marking_summary.append(
        {
            "prime": prime,
            "twist_square_class": twist_class,
            "marked_sections": sections,
            "models_with_marked_section": payload["accounting"]["models_with_marked_section"],
        }
    )

hensel = load(HENSEL)
hensel_second = load(HENSEL_SECOND)
formal = load(FORMAL)
if any(
    payload["status"]
    != "PASS_ONE_DIMENSIONAL_MARKED_TANGENT_AND_EXPLICIT_Z7_LIFT_TO_REQUESTED_PRECISION"
    for payload in (hensel, hensel_second)
):
    raise ValueError("one of the two marked seeds lost its smooth finite lift")
if [
    hensel["jacobian_certificate"]["pivot_minor_determinant_mod_prime"],
    hensel_second["jacobian_certificate"]["pivot_minor_determinant_mod_prime"],
] != [2, 6]:
    raise ValueError("marked-seed unit minors changed")
if formal["status"] != "PASS_ONE_DIMENSIONAL_FORMALLY_SMOOTH_Z7_MARKED_FAMILY":
    raise ValueError("formal smoothness certificate changed")

rational = load(RATIONAL_SCAN)
if (
    rational["search"]["candidate_count"] != 23
    or rational["search"]["status_counts"] != {"NO_FULL_RR": 23}
    or rational["exact_rational_points"]
):
    raise ValueError("bounded rational scan changed")

corridor = load(CORRIDOR)
route = corridor["results"][0]
depths = [
    {"depth": row["depth"], "best_root_rank": row["best_root_rank"]}
    for row in route["accounting"]
]
if corridor["status"] != "PASS_BOUNDED_SAME_NS_COMPILER_ROUTE_EMPTY" or depths != [
    {"depth": 1, "best_root_rank": 14},
    {"depth": 2, "best_root_rank": 10},
    {"depth": 3, "best_root_rank": 6},
    {"depth": 4, "best_root_rank": 3},
    {"depth": 5, "best_root_rank": 3},
    {"depth": 6, "best_root_rank": 2},
    {"depth": 7, "best_root_rank": 3},
    {"depth": 8, "best_root_rank": 3},
]:
    raise ValueError("bounded corridor accounting changed")

inputs = ADAPTERS + SOURCES + list(FIBRES.values()) + list(MARKINGS.values())
inputs += [HENSEL, HENSEL_SECOND, FORMAL, RATIONAL_SCAN, CORRIDOR]
payload = {
    "schema": "elkies-k3.k3-cf7f-equation-first-candidate.v1",
    "status": "PASS_FORMALLY_SMOOTH_MW1_SOURCE_PROMOTED_OPPOSITE_MW16_TARGET",
    "surface": {
        "surface_id": "K3-cf7f6c91a3a40d32",
        "determinant": 714,
        "target_frame_id": "K3-cf7f6c91a3a40d32-F001",
        "target_root_type": "A1",
        "target_mw_rank_for_rho_19": 16,
        "equation_work_authorized_by_t_arithmetic": False,
    },
    "source_censuses_by_auxiliary": source_counts,
    "promoted_source": {
        "representatives": [row["source_id"] for row in mw1_rows],
        "integral_isometry_class_count": 1,
        "representative_isometry_determinant": int(isometry_matrix.det()),
        "root_type": "A4+2A6",
        "reducible_fibres": ["I5", "I7", "I7"],
        "support_count": 3,
        "mw_rank": 1,
        "mw_height": "102/35",
        "torsion_order": 1,
        "basis_pole_profile": [1],
        "component_depths": [1, 2, 1],
    },
    "finite_field_equation_gates": {
        "fibre_censuses": fibre_summary,
        "marking_scans": marking_summary,
        "smooth_seed_count": 2,
        "seed_unit_minors_mod_7": [2, 6],
        "formal_dimension_over_Z7": 1,
    },
    "bounded_rational_scan": {
        "free_parameter": "m8",
        "integer_interval": [-40, 40],
        "candidate_count": 23,
        "exact_QQ_points": 0,
    },
    "same_ns_corridor": {
        "search": route["search"],
        "best_root_rank_by_depth": depths,
        "best_mw_rank_reached": 15,
        "target_hit": False,
    },
    "optimizer_decision": {
        "classification": "ACTIVE_SECOND_FORMALLY_SMOOTH_MW1_SOURCE",
        "comparison": (
            "This is the first formally smooth MW1 source found opposite an MW16 "
            "target and the strongest new rank-tradeoff candidate after determinant 500."
        ),
        "next_gates": [
            "rational algebraization or a rational point outside the bounded m8 scan",
            "torsion, divisibility, Picard-rank, and primitive-closure audit of any Q point",
            "targeted final corridor from the exact MW15 frontier or a mixed-degree search",
            "rootful-target D.F=2,3,4 multisection enumeration",
            "resolution of the T-arithmetic curve identification gate",
        ],
    },
    "inputs": {relative(path): digest(path) for path in inputs},
    "proof_boundary": {
        "proved": (
            "The two MW1 source Grams are integrally isometric; the finite-field fibre "
            "and marking censuses are exhaustive in their normalized charts; both positive "
            "seeds have one-dimensional smooth tangents and finite lifts; the exact identity "
            "proves a one-dimensional formally smooth Z_7 marked branch."
        ),
        "not_proved": (
            "No Q-rational source marking, primitive determinant-714 model, exact target "
            "corridor, target multisection spectrum, or specialization rank jump is proved."
        ),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/certify_k3_cf7f_equation_first_candidate.sage"
    ),
}
serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output_path = arguments.output.resolve()
if arguments.check:
    if not output_path.exists() or output_path.read_text() != serialized:
        raise SystemExit(f"stale artifact: {output_path}")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)
print("K3CF7FEQUATIONFIRST|mw1=2|formal=1|rational=0|corridor_mw=15|status=PASS")
