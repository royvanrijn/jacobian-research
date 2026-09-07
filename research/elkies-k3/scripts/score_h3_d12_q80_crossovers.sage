#!/usr/bin/env sage -python
"""Score retained generic Q80 fibre classes in the exact H3 D12 frame.

The comparison is made in one pinned Neron--Severi marking:

  Q80 stage -> initial Q80 -> pinned R17 -> current equation-side H3 D12.

The initial-Q80-to-pinned map is the already pinned canonical Q80 transport.
The pinned-to-D12 map is read from the current exact equation-D13-to-R17
certificate.  Thus no ADE-label identification or fresh qfisom is used.

For each transported fibre D, the script computes its exact D12 old-fibre
degree, the shortest marked section in its D12 Mordell--Weil coset, P.O,
height, pole-degree bounds, and a cheapest vertical decomposition

    D = (d-1) O + P + V + k F.

The reported RR ambient is a planning estimate, not an executed resolved-RR
calculation.  It uses the connected-layer rule established by the Q80
compiler:

    expected ambient = 2 + 2 P.O + vertical layer count,

where a vertical layer is a connected component of the support remaining at
each absolute coefficient level.  The exact orbit42 compiler is the pinned
regression: 2 + 2*3 + 1 = 9.
"""

from pathlib import Path
import csv
import hashlib
import json

from sage.all import (
    QQ, ZZ, QuadraticForm, block_diagonal_matrix, gcd, identity_matrix,
    matrix, pari, vector, xgcd,
)


def find_repo():
    candidates = (Path.cwd(), Path(__file__).resolve().parents[2])
    for candidate in candidates:
        if (candidate / "elkies-k3" / "AGENTS.md").exists():
            return candidate.resolve()
    raise SystemExit("run from the jacobian-research repository root")


ROOT = find_repo()
DATA = ROOT / "elkies-k3" / "data" / "fibrations"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"
GENERATED = ROOT / "artifacts" / "generated-results"

CURRENT = LOCAL / "q24-equation-d13-to-pinned-r17.json"
BRIDGE = LOCAL / "q24-orbit42-current-equation-bridge.json"
D12_PREFLIGHT = LOCAL / "q24-d12-to-a11-orbit42-divval-preflight.json"
Q80_START = DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt"
Q80_PREFIX = DATA / "kumar_q80_to_rootless_path.tsv"
Q80_SUFFIX = DATA / "kumar_q80_new_lowq_rootless_path.tsv"
PINNED_TO_Q80 = DATA / "kumar_q80_rootless_target_to_q80_ns_transport.txt"
CM24 = DATA / "kumar_q80_cm24_equation_progress.tsv"
OUTPUT = GENERATED / "h3-d12-q80-crossover-scores.json"

U2 = matrix(ZZ, ((0, 1), (1, 0)))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def ns(frame):
    return block_diagonal_matrix(U2, -matrix(ZZ, frame))


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def bezout_vector(pairings):
    current = ZZ(0)
    coefficients = [ZZ(0)] * len(pairings)
    for index, pairing in enumerate(pairings):
        if not pairing:
            continue
        new_gcd, left, right = xgcd(current, ZZ(pairing))
        coefficients = [left * value for value in coefficients]
        coefficients[index] += right
        current = new_gcd
    assert abs(current) == 1
    return vector(
        ZZ, coefficients if current == 1 else [-value for value in coefficients]
    )


def neighbor(parent, qnorm, a, b, witness):
    parent = matrix(ZZ, parent)
    witness = vector(ZZ, witness)
    gram = ns(parent)
    fibre = vector(ZZ, [a, b] + list(witness))
    assert a * b == qnorm
    assert witness * parent * witness == 2 * qnorm
    assert fibre * gram * fibre == 0
    assert gcd([abs(ZZ(value)) for value in gram * fibre]) == 1
    mate = bezout_vector(list(gram * fibre))
    mate -= ZZ(mate * gram * mate) // 2 * fibre
    complement = matrix(
        ZZ, [list(fibre * gram), list(mate * gram)]
    ).right_kernel_matrix()
    child = -(complement * gram * complement.transpose())
    transition = matrix(ZZ, [list(fibre), list(mate)] + complement.rows())
    assert abs(transition.det()) == 1
    assert transition * gram * transition.transpose() == ns(child)
    return child, transition


def connected_components(edges, active):
    active = set(active)
    components = 0
    while active:
        components += 1
        todo = [active.pop()]
        while todo:
            node = todo.pop()
            for left, right in edges:
                other = right if left == node else left if right == node else None
                if other in active:
                    active.remove(other)
                    todo.append(other)
    return components


def vertical_layers(coefficients, edges):
    magnitudes = [abs(ZZ(value)) for value in coefficients]
    total = 0
    previous = ZZ(0)
    # The coefficients transported from Q80 can be enormous.  The support is
    # constant between consecutive distinct magnitudes, so integrate the
    # connected-component count over those intervals instead of iterating once
    # per coefficient level.
    for level in sorted(set(value for value in magnitudes if value)):
        active = [index for index, value in enumerate(magnitudes) if value >= level]
        total += (level - previous) * connected_components(edges, active)
        previous = level
    return int(total)


def shortest_section_profiles(frame, divisor, root_rank, root_edges):
    """Return all closest D12-coset lifts and the cheapest vertical one."""
    frame = matrix(ZZ, frame)
    divisor = vector(ZZ, divisor)
    gram = ns(frame)
    root = frame[:root_rank, :root_rank]
    coupling = frame[:root_rank, root_rank:]
    tail = frame[root_rank:, root_rank:]
    height_gram = tail - coupling.transpose() * root.inverse() * coupling

    z = vector(ZZ, divisor[2 + root_rank:])
    height = QQ(z * height_gram * z)
    base = vector(ZZ, [0] * root_rank + list(z))
    pair = vector(ZZ, base * frame[:, :root_rank])
    dual = vector(QQ, pair) * root.inverse()

    # D12 has determinant four.  Enumerate the complete shortest shell in the
    # required dual coset via the integral adjugate form.
    det_root = ZZ(root.det())
    adjugate = (det_root * root.inverse()).change_ring(ZZ)
    candidates = []
    seen = set()
    for correction in (ZZ(0), ZZ(1), ZZ(3)):
        if correction == 0:
            lambdas = [vector(ZZ, [0] * root_rank)]
        else:
            qf = pari(adjugate).qfminim(correction * det_root)
            lambdas = []
            for column in matrix(ZZ, qf[2]).columns():
                lambdas.extend((vector(ZZ, column), -vector(ZZ, column)))
        for lam in lambdas:
            key = tuple(lam)
            if key in seen:
                continue
            seen.add(key)
            if lam * adjugate * lam != correction * det_root:
                continue
            shift_q = vector(QQ, lam - pair) * root.inverse()
            if not all(value in ZZ for value in shift_q):
                continue
            shift = vector(ZZ, [ZZ(value) for value in shift_q])
            pframe = base + vector(ZZ, list(shift) + [0] * (17 - root_rank))
            total_norm = QQ(pframe * frame * pframe)
            actual_correction = total_norm - height
            if actual_correction != correction:
                continue
            po = QQ((total_norm - 4) / 2)
            if po not in ZZ or po < 0:
                continue
            po = ZZ(po)
            section = vector(ZZ, [po + 1, 1] + list(pframe))
            assert section * gram * section == -2
            assert section * gram * vector(ZZ, [1, 0] + [0] * 17) == 1

            degree = ZZ(divisor[1])
            zero = vector(ZZ, [-1, 1] + [0] * 17)
            residual = divisor - (degree - 1) * zero - section
            assert residual[1] == 0
            assert not any(residual[2 + root_rank:])
            vertical = vector(ZZ, residual[2:2 + root_rank])
            layer_count = vertical_layers(vertical, root_edges)
            candidates.append({
                "height": height,
                "local_correction": actual_correction,
                "P_dot_O": po,
                "section": section,
                "vertical": vertical,
                "fibre_twist": ZZ(residual[0]),
                "vertical_support": sum(value != 0 for value in vertical),
                "vertical_L1": sum(abs(value) for value in vertical),
                "vertical_max": max([abs(value) for value in vertical] + [ZZ(0)]),
                "vertical_layers": layer_count,
            })
        if candidates:
            break

    if not candidates:
        raise ArithmeticError("no closest section lift in transported D12 coset")
    candidates.sort(key=lambda item: (
        item["vertical_layers"], item["vertical_support"],
        item["vertical_L1"], tuple(item["section"]),
    ))
    return candidates, candidates[0]


for path in (
    CURRENT, BRIDGE, D12_PREFLIGHT, Q80_START, Q80_PREFIX, Q80_SUFFIX,
    PINNED_TO_Q80, CM24
):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

current = json.loads(CURRENT.read_text())
bridge = json.loads(BRIDGE.read_text())
preflight = json.loads(D12_PREFLIGHT.read_text())
assert current["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert bridge["status"] == "PASS_Q24_ORBIT42_CURRENT_EQUATION_LATTICE_BRIDGE"
assert preflight["status"] == "PASS_Q42_DIVVAL_PREFLIGHT"

d12 = matrix(ZZ, current["q24"]["child_frame"])
assert d12 == matrix(ZZ, bridge["D12"]["frame"])
g_d12 = ns(d12)
pinned = load_matrix(ROOT / current["pinned_rank17_frame"])
g_pinned = ns(pinned)

# E24: D12 basis in equation-D13 coordinates.
# TOTAL: pinned basis in equation-D13 coordinates.
e24 = matrix(ZZ, current["q24"]["equation_d13_to_d12_transition"])
total = matrix(ZZ, current["equation_d13_to_pinned_r17_transition"])
pinned_basis_in_d12 = total * e24.inverse()
assert pinned_basis_in_d12.change_ring(ZZ) == pinned_basis_in_d12
pinned_basis_in_d12 = pinned_basis_in_d12.change_ring(ZZ)
assert pinned_basis_in_d12 * g_d12 * pinned_basis_in_d12.transpose() == g_pinned

q80 = load_matrix(Q80_START)
g_q80 = ns(q80)
pinned_basis_in_q80 = load_matrix(PINNED_TO_Q80)
assert pinned_basis_in_q80 * g_q80 * pinned_basis_in_q80.transpose() == g_pinned
q80_basis_in_pinned = pinned_basis_in_q80.inverse().change_ring(ZZ)

with Q80_PREFIX.open() as handle:
    prefix_rows = list(csv.DictReader(handle, delimiter="\t"))[:2]
with Q80_SUFFIX.open() as handle:
    suffix_rows = list(csv.DictReader(handle, delimiter="\t"))
assert len(prefix_rows) == 2 and len(suffix_rows) == 8

stage_specs = [{
    "stage_id": "Q80-START",
    "alias": "Q80 start",
    "ade": "E6+D5+A3",
    "mw_rank": 3,
    "basis_in_q80": identity_matrix(ZZ, 19),
}]

frame = q80
cumulative = identity_matrix(ZZ, 19)
for index, row in enumerate(prefix_rows, 1):
    child, transition = neighbor(
        frame, ZZ(row["q"]), ZZ(row["a"]), ZZ(row["b"]),
        vector(ZZ, [ZZ(value) for value in row["v"].split(",")]),
    )
    cumulative = transition * cumulative
    frame = child
    stage_specs.append({
        "stage_id": f"Q80-P{index}",
        "alias": "first q4" if index == 1 else "second q4",
        "ade": row["ADE"],
        "mw_rank": int(row["MW"]),
        "basis_in_q80": cumulative,
    })

suffix_names = (
    ("Q80-LQ1-D7D4", "escape"),
    ("Q80-LQ2-A6A4", "orbit424"),
    ("Q80-LQ3-A6A3", "orbit1222"),
    ("Q80-LQ4-A4A2A1", "q6_7774"),
    ("Q80-LQ5-A3A2", "q4_1938"),
    ("Q80-LQ6-4A1", "q4_6855"),
    ("Q80-LQ7-A1", "q4_a1_candidate1"),
    ("Q80-LQ8-ROOTLESS", "q6_rootless"),
)
for row, (stage_id, alias) in zip(suffix_rows, suffix_names):
    child, transition = neighbor(
        frame, ZZ(row["q"]), ZZ(row["a"]), ZZ(row["b"]),
        vector(ZZ, [ZZ(value) for value in row["v"].split(",")]),
    )
    cumulative = transition * cumulative
    frame = child
    stage_specs.append({
        "stage_id": stage_id,
        "alias": alias,
        "ade": row["ADE"],
        "mw_rank": int(row["MW"]),
        "basis_in_q80": cumulative,
    })

with CM24.open() as handle:
    cm24_by_stage = {
        row["stage"]: row for row in csv.DictReader(handle, delimiter="\t")
    }

cm24_labels = {
    "escape": "old q12 section",
    "orbit424": "rational 2-torsion chord",
    "orbit1222": "P.O=1 A7 saturated chord",
    "q6_7774": "P3",
    "q4_1938": "-P1+P2+2P3",
    "q4_6855": "2P1",
    "q4_a1_candidate1": "-P3",
    "q6_rootless": "P2-P3",
}
cm24_word_l1 = {
    "q6_7774": 1,
    "q4_1938": 4,
    "q4_6855": 2,
    "q4_a1_candidate1": 1,
    "q6_rootless": 2,
}

root_edges = [
    tuple(map(int, edge))
    for edge in preflight["abstract_D12_marking"]["root_edges"]
]

scores = []
for stage in stage_specs:
    fibre_q80 = vector(ZZ, stage["basis_in_q80"].row(0))
    fibre_pinned = vector(ZZ, fibre_q80 * q80_basis_in_pinned)
    fibre_d12 = vector(ZZ, fibre_pinned * pinned_basis_in_d12)
    assert fibre_d12 * g_d12 * fibre_d12 == 0
    if fibre_d12[1] < 0:
        fibre_d12 = -fibre_d12
    degree = ZZ(fibre_d12[1])
    assert degree > 0

    shell, best = shortest_section_profiles(d12, fibre_d12, 12, root_edges)
    po = ZZ(best["P_dot_O"])
    expected_ambient = ZZ(2 + 2 * po + best["vertical_layers"])
    alias = stage["alias"]
    cm = cm24_by_stage.get(alias)
    cm_record = None
    if cm is not None:
        cm_record = {
            "horizontal": cm24_labels[alias],
            "word_L1": cm24_word_l1.get(alias),
            "P_dot_O": int(cm["section_po"]),
            "height": cm["section_height"],
            "fibre_twist": int(cm["fiber_twist"]),
            "vertical_summary": cm["vertical_summary"],
            "equation_status": cm["equation_status"],
        }

    record = {
        "stage_id": stage["stage_id"],
        "alias": alias,
        "generic_frame": f"{stage['ade']}/MW{stage['mw_rank']}",
        "fibre_in_current_D12": [int(value) for value in fibre_d12],
        "d": int(degree),
        "P_dot_O": int(po),
        "MW_height": str(best["height"]),
        "local_correction": str(best["local_correction"]),
        "vertical_components": int(best["vertical_support"]),
        "vertical_L1": int(best["vertical_L1"]),
        "vertical_max": int(best["vertical_max"]),
        "vertical_layers": int(best["vertical_layers"]),
        "fibre_twist": int(best["fibre_twist"]),
        "closest_section_lifts": len(shell),
        "expected_RR_ambient": int(expected_ambient),
        "expected_RR_formula": "2 + 2*P.O + connected vertical layers",
        "pole_degrees": {
            "Z": int(po),
            "X_max": int(2 * po + 4),
            "Y_max": int(3 * po + 6),
        },
        "CM24_section_complexity": cm_record,
    }
    scores.append(record)
    print(
        "H3D12Q80|stage={}|frame={}|d={}|P.O={}|height={}|"
        "vcomponents={}|vlayers={}|twist={}|ambient={}|poles={},{},{}|"
        "cm24={}|status=PASS".format(
            record["stage_id"], record["generic_frame"], record["d"],
            record["P_dot_O"], record["MW_height"],
            record["vertical_components"], record["vertical_layers"],
            record["fibre_twist"], record["expected_RR_ambient"],
            record["pole_degrees"]["Z"], record["pole_degrees"]["X_max"],
            record["pole_degrees"]["Y_max"],
            "none" if cm_record is None else cm_record["horizontal"],
        ),
        flush=True,
    )

# Exact regression against the executed orbit42 D12 compiler.
orbit42 = bridge["D12"]
assert int(orbit42["orbit42_P_dot_O"]) == 3
selected = orbit42["selected_section_representative"]
assert int(selected["vertical_support"]) == 11
orbit42_layers = vertical_layers(
    selected["vertical_root_coefficients"], root_edges
)
assert orbit42_layers == 1
assert 2 + 2 * 3 + orbit42_layers == 9

ranked = sorted(scores, key=lambda row: (
    row["d"], row["expected_RR_ambient"], row["P_dot_O"],
    row["vertical_layers"], row["vertical_components"], row["stage_id"],
))
best = ranked[0]
candidate1 = next(row for row in scores if row["stage_id"] == "Q80-LQ7-A1")

payload = {
    "schema": "elkies-k3.h3-d12-q80-crossover-scores.v1",
    "status": "PASS_EXACT_H3_D12_Q80_CROSSOVER_AUDIT",
    "proof_boundary": (
        "Exact integral transport and D12 lattice scoring for every retained Q80 "
        "stage. Pole degrees are exact D12 section bounds. expected_RR_ambient "
        "is a connected-layer planning estimate calibrated by the executed "
        "orbit42 ambient 9; it is not a resolved-RR kernel or an equation lift. "
        "CM24 records describe specialized compiler complexity and do not define "
        "the generic transported classes."
    ),
    "inputs": {
        str(path.relative_to(ROOT)): digest(path)
        for path in (
            CURRENT, BRIDGE, D12_PREFLIGHT, Q80_START, Q80_PREFIX, Q80_SUFFIX,
            PINNED_TO_Q80, CM24,
        )
    },
    "transport": {
        "method": "Q80 stage -> initial Q80 -> pinned R17 -> current exact D12",
        "pinned_basis_in_D12": rows(pinned_basis_in_d12),
        "q80_basis_in_pinned": rows(q80_basis_in_pinned),
        "det_pinned_basis_in_D12": int(pinned_basis_in_d12.det()),
        "det_q80_basis_in_pinned": int(q80_basis_in_pinned.det()),
    },
    "orbit42_RR_estimate_regression": {
        "P_dot_O": 3,
        "vertical_layers": 1,
        "estimated_ambient": 9,
        "executed_exact_ambient": 9,
    },
    "scores": scores,
    "ranking": [row["stage_id"] for row in ranked],
    "best_stage": best["stage_id"],
    "candidate1": candidate1,
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(
    "H3D12Q80_RESULT|stages={}|best={}|best_d={}|candidate1_d={}|"
    "candidate1_P.O={}|candidate1_ambient={}|artifact={}|"
    "status=PASS_EXACT_H3_D12_Q80_CROSSOVER_AUDIT".format(
        len(scores), best["stage_id"], best["d"], candidate1["d"],
        candidate1["P_dot_O"], candidate1["expected_RR_ambient"],
        OUTPUT.relative_to(ROOT),
    ),
    flush=True,
)
