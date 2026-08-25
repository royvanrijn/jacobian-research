#!/usr/bin/env sage -python
"""Score exhaustive equation-side A11 neighbours by expected compiler cost.

This is a route-planning certificate, not an equation lift.  It combines the
exact orbit64 A11 marking, the exhaustive declared neighbour shells, and the
eighteen exact characteristic-zero identity-shell sections.  For every
primitive neighbour it computes:

* the cheapest section representative of the horizontal MW coset;
* P.O, vertical support/layers, and an RR-ambient planning estimate;
* the quotient obstruction modulo the exact-section rank-five lattice;
* how many already-explicit sections have degree zero or one over the new
  fibre; and
* a transparent weighted equation-cost score.

The score is deliberately compiler-facing.  It is not a theorem that the
lowest-scoring edge is cheapest to lift, and it does not replace nefness,
full transport, or pinned-endpoint certification for a retained route.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    IntegralLattice,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    matrix,
    vector,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GENERATED = ROOT / "artifacts/generated-results"
U2 = matrix(ZZ, ((0, 1), (1, 0)))

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--neighbors",
    type=Path,
    action="append",
    default=[],
    help="root-adapted-neighbour JSON; repeat to combine declared shells",
)
parser.add_argument(
    "--output",
    type=Path,
    default=GENERATED / "elkies-k3-h3-a11-equation-cost-neighbors.json",
)
parser.add_argument("--retain", type=int, default=100)
args = parser.parse_args()

DEFAULT_NEIGHBORS = LOCAL / "q24-a11-orbit64-q8-all.json"
NEIGHBOR_PATHS = [path.resolve() for path in args.neighbors] or [DEFAULT_NEIGHBORS]
PARENT_FRAME_PATH = LOCAL / "q24-downstream-lift/d12-c10a-zero-frame.txt"
Q6_PATH = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
IDENTITY_PATH = LOCAL / "q24-orbit42-identity-halving-audit.json"
MATCHING_PATH = LOCAL / "q24-orbit42-identity-halving-qq.json"
ZERO_PATH = LOCAL / "q24-orbit42-rational-zero-pole-sections-qq.json"
BRIDGE_PATH = LOCAL / "q24-a11-target-coset-bridge.json"
ZERO_MISMATCH_PATH = GENERATED / "elkies-k3-h3-a11-quintic-bridge-zero-mismatch.json"
INPUTS = (
    PARENT_FRAME_PATH,
    Q6_PATH,
    IDENTITY_PATH,
    MATCHING_PATH,
    ZERO_PATH,
    BRIDGE_PATH,
    ZERO_MISMATCH_PATH,
    *NEIGHBOR_PATHS,
)
for path in INPUTS:
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def connected_components(edges, active):
    active = set(active)
    count = 0
    while active:
        count += 1
        todo = [active.pop()]
        while todo:
            node = todo.pop()
            for left, right in edges:
                other = right if left == node else left if right == node else None
                if other in active:
                    active.remove(other)
                    todo.append(other)
    return count


def vertical_layers(coefficients, edges):
    magnitudes = [abs(ZZ(value)) for value in coefficients]
    total = ZZ(0)
    previous = ZZ(0)
    for level in sorted(set(value for value in magnitudes if value)):
        active = [index for index, value in enumerate(magnitudes) if value >= level]
        total += (level - previous) * connected_components(edges, active)
        previous = level
    return int(total)


parent = load_matrix(PARENT_FRAME_PATH)
q6 = json.loads(Q6_PATH.read_text())
identity = json.loads(IDENTITY_PATH.read_text())
matching = json.loads(MATCHING_PATH.read_text())
zero = json.loads(ZERO_PATH.read_text())
bridge = json.loads(BRIDGE_PATH.read_text())
zero_mismatch = json.loads(ZERO_MISMATCH_PATH.read_text())
assert q6["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert identity["status"] == "PASS_Q42_ORBIT42_IDENTITY_HALVING_LATTICE_GATE"
assert matching["status"] == "Q42_IDENTITY_HALVING_HAS_NO_A11_CHORD"
assert zero["status"] == "PASS_EXACT_Q42_RATIONAL_ZERO_POLE_SECTIONS_QQ"
assert bridge["status"] == "PASS_EXACT_A11_TARGET_COSET_BRIDGE"
assert zero_mismatch["status"] == "REJECT_A11_QUINTIC_BRIDGE_ZERO_MISMATCH"

selected_q6 = next(row for row in q6["neighbors"] if int(row["orbit_index"]) == 64)
transition = block_diagonal_matrix(
    identity_matrix(ZZ, 2), matrix(ZZ, selected_q6["child_root_adapted_basis"])
) * matrix(ZZ, selected_q6["neighbor_basis"])
assert abs(transition.det()) == 1
transition_inverse = transition.inverse().change_ring(ZZ)
a11 = matrix(ZZ, selected_q6["child_root_adapted_frame"])
g_a11 = block_diagonal_matrix(U2, -a11)
root_rank = 11
root = a11[:root_rank, :root_rank]
coupling = a11[:root_rank, root_rank:]
tail = a11[root_rank:, root_rank:]
height = tail - coupling.transpose() * root.inverse() * coupling
root_lattice = IntegralLattice(root)
root_edges = [
    (left, right)
    for left in range(root_rank)
    for right in range(left + 1, root_rank)
    if root[left, right] == -1
]

# Reconstruct the exact parent P.O=0 identity shell as full divisor classes,
# then transport those actual curves into the equation A11 marking.
parent_root_rank = 12
parent_root = parent[:parent_root_rank, :parent_root_rank]
parent_coupling = parent[:parent_root_rank, parent_root_rank:]
old_sections = []
for values in identity["exact_model_R3_zero"]["identity_vectors"]:
    z = vector(ZZ, values)
    root_coefficients = -(z * parent_coupling.transpose()) * parent_root.inverse()
    assert all(value in ZZ for value in root_coefficients)
    section = vector(
        ZZ,
        [1, 1] + list(map(ZZ, root_coefficients)) + list(z),
    )
    assert section * block_diagonal_matrix(U2, -parent) * section == -2
    old_sections.append(section)
explicit_sections = [section * transition_inverse for section in old_sections]
assert len(explicit_sections) == len(zero["sections"]) == 18
assert all(section * g_a11 * section == -2 for section in explicit_sections)

equation_mw = [vector(ZZ, values) for values in bridge["exact_identity_shell"]["MW_vectors_in_equation_order"]]
shell_mapping = matching["matching"]["mappings_abstract_to_equation"][7]
reordered_sections = [None] * len(explicit_sections)
for abstract_index, equation_index in enumerate(shell_mapping):
    reordered_sections[equation_index] = explicit_sections[abstract_index]
explicit_sections = reordered_sections
assert all(section is not None for section in explicit_sections)
assert [vector(ZZ, section[-6:]) for section in explicit_sections] == equation_mw

known_lattice = matrix(ZZ, [list(item[:5]) for item in equation_mw]).row_module()
assert known_lattice.rank() == 5
assert abs(known_lattice.basis_matrix().det()) == 5
a0_record = zero_mismatch["correct_selected_R3_transport"]["oldI9_A0"]
quintic_record = zero_mismatch["correct_selected_R3_transport"]["close_P24"]
a0_mw = vector(ZZ, a0_record["A11_MW_Abel_Jacobi"])
a0_curve = vector(ZZ, a0_record["child_coordinates"])
quintic_curve = vector(ZZ, quintic_record["child_coordinates"])
assert a0_mw[-1] == 0 and vector(ZZ, a0_mw[:5]) in known_lattice
saturated_known_lattice = matrix(
    ZZ, [list(item[:5]) for item in equation_mw] + [list(a0_mw[:5])]
).row_module()
assert saturated_known_lattice.rank() == 5
assert abs(saturated_known_lattice.basis_matrix().det()) == 5
assert a0_curve * g_a11 * a0_curve == quintic_curve * g_a11 * quintic_curve == -2
assert a0_curve[1] == 4 and quintic_curve[1] == 46


profile_cache = {}


def section_profile(z):
    """Return every closest root lift of one A11 MW vector."""
    z = vector(ZZ, z)
    key = tuple(z)
    if key in profile_cache:
        return profile_cache[key]
    horizontal_height = QQ(z * height * z)
    base = vector(ZZ, [0] * root_rank + list(z))
    dual = vector(QQ, base * a11[:, :root_rank]) * root.inverse()
    iterator = root_lattice.enumerate_close_vectors(-dual)
    minimum = None
    profiles = []
    for unused in range(100000):
        shift = vector(ZZ, next(iterator))
        lifted = base + vector(ZZ, list(shift) + [0] * 6)
        norm = QQ(lifted * a11 * lifted)
        if minimum is None:
            minimum = norm
        elif norm > minimum:
            break
        correction = norm - horizontal_height
        pole_order = (norm - 4) / 2
        if pole_order in ZZ and pole_order >= 0:
            profiles.append((lifted, horizontal_height, correction, ZZ(pole_order)))
    if not profiles:
        raise ArithmeticError(f"no effective section profile for MW vector {key}")
    profile_cache[key] = profiles
    return profiles


records = []
shell_summaries = []
seen_ids = set()
for neighbor_path in NEIGHBOR_PATHS:
    payload = json.loads(neighbor_path.read_text())
    assert payload["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    payload_frame = payload["frame"]
    payload_frame = (
        load_matrix(ROOT / payload_frame)
        if isinstance(payload_frame, str)
        else matrix(ZZ, payload_frame)
    )
    assert payload_frame == a11
    shell_summaries.extend(payload["summaries"])
    for raw in payload["neighbors"]:
        candidate_id = (int(raw["q"]), int(raw["old_fiber_degree"]), int(raw["orbit_index"]))
        if candidate_id in seen_ids:
            continue
        seen_ids.add(candidate_id)
        fibre = vector(ZZ, raw["fiber"])
        z = vector(ZZ, raw["mw_projection"])
        assert fibre * g_a11 * fibre == 0
        degree = ZZ(raw["old_fiber_degree"])
        old_zero = vector(ZZ, [-1, 1] + [0] * 17)
        old_fibre = vector(ZZ, [1, 0] + [0] * 17)

        profile_rows = []
        for lifted, horizontal_height, correction, pole_order in section_profile(z):
            section = vector(ZZ, [pole_order + 1, 1] + list(lifted))
            assert section * g_a11 * section == -2
            residual = fibre - (degree - 1) * old_zero - section
            assert residual[1] == 0
            assert not any(residual[2 + root_rank :])
            vertical = vector(ZZ, residual[2 : 2 + root_rank])
            layers = vertical_layers(vertical, root_edges)
            profile_rows.append(
                {
                    "height": str(horizontal_height),
                    "local_correction": str(correction),
                    "P_dot_O": int(pole_order),
                    "section": entries(section),
                    "vertical": entries(vertical),
                    "fibre_twist": int(residual[0]),
                    "vertical_support": int(sum(value != 0 for value in vertical)),
                    "vertical_L1": int(sum(abs(value) for value in vertical)),
                    "vertical_max": int(max([abs(value) for value in vertical] + [ZZ(0)])),
                    "vertical_layers": layers,
                }
            )
        profile_rows.sort(
            key=lambda row: (
                row["P_dot_O"],
                row["vertical_layers"],
                row["vertical_support"],
                row["vertical_L1"],
                tuple(row["section"]),
            )
        )
        horizontal = profile_rows[0]

        explicit_degrees = [int(section * g_a11 * fibre) for section in explicit_sections]
        a0_degree = int(a0_curve * g_a11 * fibre)
        quintic_degree = int(quintic_curve * g_a11 * fibre)
        dominant_labels = list(map(int, raw["dominant_labels"]))
        identity_component_degree = int(degree - sum(dominant_labels))
        physical_component_degrees = dominant_labels + [identity_component_degree]
        physical_degree_zero = physical_component_degrees.count(0)
        physical_degree_one = physical_component_degrees.count(1)
        negative = [index for index, value in enumerate(explicit_degrees) if value < 0]
        named_negative = [
            name
            for name, value in (("oldI9_A0", a0_degree), ("close_P24", quintic_degree))
            if value < 0
        ]
        if identity_component_degree < 0:
            named_negative.append("old_A11_identity_component")
        degree_zero = [index for index, value in enumerate(explicit_degrees) if value == 0]
        degree_one = [index for index, value in enumerate(explicit_degrees) if value == 1]
        named_degrees = (a0_degree, quintic_degree)
        total_explicit_degree_zero = len(degree_zero) + physical_degree_zero + named_degrees.count(0)
        total_explicit_degree_one = len(degree_one) + physical_degree_one + named_degrees.count(1)

        missing_coordinate = abs(int(z[-1]))
        in_shell_lattice = missing_coordinate == 0 and vector(ZZ, z[:5]) in known_lattice
        in_exact_hyperplane_lattice = in_shell_lattice
        finite_coset_gap = 0 if in_exact_hyperplane_lattice else 1 if missing_coordinate == 0 else None
        construction_tier = 0 if in_exact_hyperplane_lattice else 1 if missing_coordinate == 0 else 2
        rr_ambient = 2 + 2 * horizontal["P_dot_O"] + horizontal["vertical_layers"]
        coordinate_growth = max(abs(int(value)) for value in fibre)
        component_complexity = int(raw["child_root_data"][1])

        # Transparent planning weights.  Missing rank direction dominates a
        # finite index-five coset, then P.O/RR size, while explicit degree-one
        # curves receive a large credit.  Raw fields are retained so alternate
        # weightings can be replayed without rerunning Sage.
        cost_terms = {
            "declared_curve_non_nef_penalty": 1000000 if negative or named_negative else 0,
            "missing_rank_direction": 12000 * missing_coordinate,
            "finite_index_five_coset": 2500 * (finite_coset_gap or 0),
            "P_dot_O": 900 * horizontal["P_dot_O"],
            "horizontal_degree": 250 * int(degree),
            "RR_ambient": 120 * rr_ambient,
            "vertical_layers": 60 * horizontal["vertical_layers"],
            "vertical_support": 25 * horizontal["vertical_support"],
            "child_root_count": component_complexity,
            "coordinate_growth": coordinate_growth,
            "no_explicit_degree_one_curve": 3000 if not total_explicit_degree_one else 0,
            "explicit_degree_one_credit": -500 * min(total_explicit_degree_one, 4),
            "explicit_degree_zero_credit": -100 * min(total_explicit_degree_zero, 8),
        }
        score = sum(cost_terms.values())
        records.append(
            {
                "candidate_id": {
                    "q": candidate_id[0],
                    "old_fibre_degree": candidate_id[1],
                    "orbit_index": candidate_id[2],
                },
                "declared_curve_nef_gate": "PASS" if not negative and not named_negative else "REJECT",
                "child": {
                    "ade": raw["child_ade"],
                    "mw_rank": int(raw["child_mw_rank"]),
                    "root_data": raw["child_root_data"],
                },
                "fibre": entries(fibre),
                "mw_projection": entries(z),
                "horizontal": horizontal,
                "closest_horizontal_count": len(profile_rows),
                "expected_RR_ambient": rr_ambient,
                "explicit_curve_degrees": {
                    "identity_shell": explicit_degrees,
                    "negative_indices": negative,
                    "oldI9_A0": a0_degree,
                    "close_P24": quintic_degree,
                    "named_negative_curves": named_negative,
                    "physical_old_A11_fibre_components": physical_component_degrees,
                    "physical_component_degree_zero_count": physical_degree_zero,
                    "physical_component_degree_one_count": physical_degree_one,
                    "total_explicit_degree_zero_count": total_explicit_degree_zero,
                    "total_explicit_degree_one_count": total_explicit_degree_one,
                    "degree_zero_indices": degree_zero,
                    "degree_one_indices": degree_one,
                    "degree_zero_count": len(degree_zero),
                    "degree_one_count": len(degree_one),
                },
                "target_coset_mod_exact_sections": {
                    "exact_section_span_rank": 5,
                    "index_in_saturated_sixth_coordinate_zero_hyperplane": 5,
                    "missing_sixth_coordinate_absolute": missing_coordinate,
                    "in_exact_section_lattice": in_exact_hyperplane_lattice,
                    "in_identity_shell_lattice_before_A0_saturation": in_shell_lattice,
                    "finite_index_five_coset_gap": finite_coset_gap,
                    "construction_tier": construction_tier,
                },
                "coordinate_growth_max": coordinate_growth,
                "equation_cost_terms": cost_terms,
                "equation_cost_score": int(score),
            }
        )

records.sort(
    key=lambda row: (
        row["equation_cost_score"],
        row["target_coset_mod_exact_sections"]["construction_tier"],
        row["horizontal"]["P_dot_O"],
        -row["explicit_curve_degrees"]["degree_one_count"],
        row["expected_RR_ambient"],
        tuple(row["candidate_id"].values()),
    )
)
retained = records[: args.retain]

named = {}
for orbit in (12, 2162):
    hits = [
        row
        for row in records
        if row["candidate_id"]["q"] == 8
        and row["candidate_id"]["old_fibre_degree"] == 2
        and row["candidate_id"]["orbit_index"] == orbit
    ]
    if hits:
        assert len(hits) == 1
        named[f"orbit{orbit}"] = hits[0]

payload = {
    "schema": "elkies-k3.h3-a11-equation-cost-neighbors.v1",
    "status": "PASS_EXACT_A11_EQUATION_COST_SCORING",
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in INPUTS],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in INPUTS
        },
    },
    "declared_search_summaries": shell_summaries,
    "exact_section_subgroup": {
        "section_count": len(explicit_sections),
        "rank": 5,
        "saturated_hyperplane_index": 5,
        "index_after_adjoining_exact_oldI9_A0": 5,
        "oldI9_A0_MW": entries(a0_mw),
        "mw_vectors": [entries(value) for value in equation_mw],
    },
    "score_definition": {
        "ordering": "weighted integer score, then transparent deterministic tie-breaks",
        "weights_are_planning_convention": True,
        "formula": "sum(equation_cost_terms)",
    },
    "candidate_count": len(records),
    "retained_count": len(retained),
    "best_candidate": retained[0],
    "named_construction_candidates": named,
    "retained_candidates": retained,
    "proof_boundary": (
        "Exact lattice arithmetic in the selected equation-side A11 marking, "
        "using exhaustive declared neighbour artifacts and exact QQ section "
        "classes. RR ambient and the weighted total are route-planning estimates; "
        "they do not prove an equation lift, global route optimality, or a pinned "
        "R17 endpoint for any newly retained branch. A negative exact-section "
        "degree rejects nefness; nonnegative degrees certify only the declared curves."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
best = retained[0]
print(
    "A11EQCOST|candidates={}|best=q{}d{}o{}|child={}|PO={}|RR={}|"
    "deg0={}|deg1={}|tier={}|score={}|status={}".format(
        len(records),
        best["candidate_id"]["q"],
        best["candidate_id"]["old_fibre_degree"],
        best["candidate_id"]["orbit_index"],
        best["child"]["ade"],
        best["horizontal"]["P_dot_O"],
        best["expected_RR_ambient"],
        best["explicit_curve_degrees"]["degree_zero_count"],
        best["explicit_curve_degrees"]["degree_one_count"],
        best["target_coset_mod_exact_sections"]["construction_tier"],
        best["equation_cost_score"],
        payload["status"],
    ),
    flush=True,
)
print(f"OUTPUT|{args.output.resolve()}", flush=True)
