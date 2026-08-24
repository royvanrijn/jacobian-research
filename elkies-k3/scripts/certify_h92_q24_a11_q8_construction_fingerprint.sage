#!/usr/bin/env sage -python
"""Transport the historical A11 -> 2A5 construction to the equation A11.

The historical q8/orbit922 arrow has a much stronger fingerprint than its
root data: its new fibre is

    O + P - 2 F,

with no additional vertical-root divisor.  This script transports that exact
construction through every fibre-preserving integral isometry between the
historical A11 frame and the equation-side orbit64 A11 frame.

The rank-17 positive-definite frame isometry problem is split into the A11
root chain and its rank-6 Mordell--Weil quotient.  Only root-chain
orientations and MW isometries which satisfy the integral glue condition are
retained.  Their transported q8 rays are then matched literally against the
exhaustive equation-side q8 search.

This is a lattice/marking/construction certificate.  It selects a convenient
equation-side q8 target, but does not construct its section function or the
characteristic-zero 2-neighbour pencil.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    QQ,
    ZZ,
    block_diagonal_matrix,
    block_matrix,
    identity_matrix,
    lcm,
    matrix,
    pari,
    vector,
    zero_matrix,
)


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--output",
    type=Path,
    default=LOCAL / "q24-a11-q8-construction-fingerprint.json",
)
args = parser.parse_args()

MANIFEST = LOCAL / "h3-r17-backward-exact-lift-manifest.json"
CLOSEOUT = LOCAL / "q24-equation-d13-to-pinned-r17.json"
Q6_SEARCH = LOCAL / "q24-downstream-lift/d12-c10a-zero-q6-all.json"
Q8_SEARCH = LOCAL / "q24-a11-orbit64-q8-all.json"
MARKING = LOCAL / "q24-a11-equation-marking-orbit64-mod100003.json"

for path in (MANIFEST, CLOSEOUT, Q6_SEARCH, Q8_SEARCH, MARKING):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

manifest = json.loads(MANIFEST.read_text())
closeout = json.loads(CLOSEOUT.read_text())
q6_search = json.loads(Q6_SEARCH.read_text())
q8_search = json.loads(Q8_SEARCH.read_text())
marking = json.loads(MARKING.read_text())

assert marking["status"] == "PASS_Q42_A11_EQUATION_MARKING_ORBIT64_MOD100003"
assert marking["selected"]["orbit_index"] == 64
assert marking["selected"]["mapping_index"] == 7
assert q6_search["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
assert q8_search["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"

U2 = matrix(ZZ, ((0, 1), (1, 0)))


def ns(frame):
    return block_diagonal_matrix(U2, -matrix(ZZ, frame))


def rows(value):
    return [[int(entry) for entry in row] for row in matrix(ZZ, value).rows()]


def entries(value):
    return [int(entry) for entry in vector(ZZ, value)]


def height(frame, root_rank=11):
    frame = matrix(ZZ, frame)
    root = frame[:root_rank, :root_rank]
    coupling = frame[:root_rank, root_rank:]
    tail = frame[root_rank:, root_rank:]
    return tail - coupling.transpose() * root.inverse() * coupling


def manifest_step(parent, child):
    hits = [
        step
        for step in manifest["forward_steps"]
        if step["parent"] == parent and step["child"] == child
    ]
    assert len(hits) == 1
    return hits[0]


# Reconstruct the historical A11 frame from the stored D12 frame and the
# literal historical transition.  No missing generated neighbour file is
# needed for this calculation.
d12_to_a11 = manifest_step("D12/MW5", "A11/MW6")
a11_to_2a5 = manifest_step("A11/MW6", "2A5/MW7")
historical_d12 = matrix(ZZ, closeout["q24"]["child_frame"])
historical_transition = matrix(ZZ, d12_to_a11["transition"])
historical_a11 = -(
    historical_transition
    * ns(historical_d12)
    * historical_transition.transpose()
)[2:, 2:]

equation_q6 = [
    record
    for record in q6_search["neighbors"]
    if int(record["orbit_index"]) == 64
]
assert len(equation_q6) == 1
equation_a11 = matrix(ZZ, equation_q6[0]["child_root_adapted_frame"])
assert historical_a11.det() == equation_a11.det() == 948

root_rank = 11
Re = equation_a11[:root_rank, :root_rank]
Ce = equation_a11[:root_rank, root_rank:]
Te = equation_a11[root_rank:, root_rank:]
Rh = historical_a11[:root_rank, :root_rank]
Ch = historical_a11[:root_rank, root_rank:]
Th = historical_a11[root_rank:, root_rank:]
He = height(equation_a11)
Hh = height(historical_a11)

# PARI qfisom uses a column convention.  If Q=qfisom(E,H), then
# Q^t H Q=E; M=Q^(-t) is the desired row-convention isometry M E M^t=H.
scale = lcm(entry.denominator() for entry in list(He) + list(Hh))
Ei = (scale * He).change_ring(ZZ)
Hi = (scale * Hh).change_ring(ZZ)
qiso = matrix(ZZ, pari(Ei).qfisom(pari(Hi)))
M0 = qiso.inverse().transpose().change_ring(ZZ)
assert M0 * Ei * M0.transpose() == Hi

# qfauto generators also use the column convention.  Enumerate the full
# (order 16 here) row-convention MW automorphism group exactly.
automorphism_data = pari(Ei).qfauto()
claimed_automorphism_order = int(automorphism_data[0])
generators = [matrix(ZZ, item).transpose() for item in automorphism_data[1]]
identity = identity_matrix(ZZ, Ei.nrows())
automorphisms = {tuple(identity.list()): identity}
frontier = [identity]
while frontier:
    current = frontier.pop()
    for generator in generators:
        candidate = current * generator
        key = tuple(candidate.list())
        if key not in automorphisms:
            assert candidate * Ei * candidate.transpose() == Ei
            automorphisms[key] = candidate
            frontier.append(candidate)
assert len(automorphisms) == claimed_automorphism_order == 16


def endpoints(root):
    return [
        index
        for index in range(root_rank)
        if sum(root[index, other] == -1 for other in range(root_rank)) == 1
    ]


def chain_order(root, start):
    order = []
    previous = None
    current = start
    while True:
        order.append(current)
        following = [
            other
            for other in range(root_rank)
            if root[current, other] == -1 and other != previous
        ]
        if not following:
            return order
        assert len(following) == 1
        previous, current = current, following[0]


historical_chain = chain_order(Rh, endpoints(Rh)[0])
equation_chain = chain_order(Re, endpoints(Re)[0])
root_isometries = []
for equation_order in (equation_chain, list(reversed(equation_chain))):
    P = zero_matrix(ZZ, root_rank, root_rank)
    for position, historical_index in enumerate(historical_chain):
        P[historical_index, equation_order[position]] = 1
    assert P * Re * P.transpose() == Rh
    root_isometries.append(P)

# A fibre-preserving full-frame isometry has block form [[P,0],[K,M]].
# The root/MW coupling determines K; its integrality is exactly the glue gate.
glue_isometries = []
for root_orientation, P in enumerate(root_isometries):
    for M in automorphisms.values():
        quotient_isometry = M0 * M
        K_transpose = Re.inverse() * (
            P.inverse() * Ch - Ce * quotient_isometry.transpose()
        )
        if not all(entry in ZZ for entry in K_transpose.list()):
            continue
        K = K_transpose.transpose().change_ring(ZZ)
        full = block_matrix(
            [
                [P, zero_matrix(ZZ, root_rank, Ei.nrows())],
                [K, quotient_isometry],
            ]
        )
        assert abs(full.det()) == 1
        assert full * equation_a11 * full.transpose() == historical_a11
        glue_isometries.append((root_orientation, full))
assert len(glue_isometries) == 8

historical_fibre = vector(ZZ, a11_to_2a5["new_fibre_in_parent"])
historical_section = vector(ZZ, a11_to_2a5["horizontal"]["section_class"])
historical_horizontal = a11_to_2a5["horizontal"]
assert historical_horizontal["vertical_root_support"] == 0
assert historical_horizontal["vertical_root_L1"] == 0
assert historical_horizontal["vertical_root_coefficients"] == [0] * root_rank
assert historical_fibre == vector(ZZ, [4, 2] + list(historical_fibre[2:]))
assert historical_section == vector(ZZ, [7, 1] + list(historical_fibre[2:]))

q8_by_fibre = {
    tuple(map(int, record["fiber"])): record for record in q8_search["neighbors"]
}
transport_records = []
for root_orientation, full in glue_isometries:
    transported_witness = vector(ZZ, historical_fibre[2:]) * full
    transported_fibre = vector(ZZ, [4, 2] + list(transported_witness))
    transported_section = vector(ZZ, [7, 1] + list(transported_witness))
    assert transported_fibre * ns(equation_a11) * transported_fibre == 0
    assert transported_section * ns(equation_a11) * transported_section == -2
    old_fibre = vector(ZZ, [1, 0] + [0] * 17)
    old_zero = vector(ZZ, [-1, 1] + [0] * 17)
    assert transported_fibre * ns(equation_a11) * old_fibre == 2
    assert transported_section * ns(equation_a11) * old_fibre == 1
    assert transported_section * ns(equation_a11) * old_zero == 6
    assert transported_fibre == old_zero + transported_section - 2 * old_fibre

    match = q8_by_fibre.get(tuple(map(int, transported_fibre)))
    transport_records.append(
        {
            "root_chain_orientation": int(root_orientation),
            "frame_isometry_historical_basis_in_equation_coordinates": rows(full),
            "transported_fibre": entries(transported_fibre),
            "transported_section": entries(transported_section),
            "equation_nef_orbit": None if match is None else int(match["orbit_index"]),
            "equation_child_ade": None if match is None else match["child_ade"],
            "equation_mw_projection": (
                None if match is None else list(map(int, match["mw_projection"]))
            ),
        }
    )

nef_records = [record for record in transport_records if record["equation_nef_orbit"]]
compatible_orbits = sorted({record["equation_nef_orbit"] for record in nef_records})
assert compatible_orbits == [12, 2162]
for record in nef_records:
    assert record["equation_child_ade"] == "A5+A5"

# Both rays reproduce the historical construction exactly.  Choose the ray
# with the smallest MW-coordinate L1 norm; this is an explicit complexity
# convention, not a claim that the other ray is geometrically impossible.
orbit_profiles = []
for orbit in compatible_orbits:
    record = next(
        item for item in q8_search["neighbors"] if int(item["orbit_index"]) == orbit
    )
    orbit_profiles.append(
        {
            "orbit_index": orbit,
            "mw_projection": list(map(int, record["mw_projection"])),
            "mw_L1": sum(abs(int(entry)) for entry in record["mw_projection"]),
            "child_ade": record["child_ade"],
            "child_root_data": list(map(int, record["child_root_data"])),
            "child_mw_rank": int(record["child_mw_rank"]),
        }
    )
selected_profile = min(orbit_profiles, key=lambda item: (item["mw_L1"], item["orbit_index"]))
assert selected_profile["orbit_index"] == 12
selected_candidates = [
    record
    for record in nef_records
    if record["equation_nef_orbit"] == selected_profile["orbit_index"]
]
selected = min(
    selected_candidates,
    key=lambda item: tuple(
        entry
        for row in item["frame_isometry_historical_basis_in_equation_coordinates"]
        for entry in row
    ),
)

payload = {
    "schema": "elkies-k3.h3-q24-a11-q8-construction-fingerprint.v1",
    "status": "PASS_Q24_A11_Q8_CONSTRUCTION_FINGERPRINT",
    "inputs": {
        "paths": [
            str(path.relative_to(ROOT))
            for path in (MANIFEST, CLOSEOUT, Q6_SEARCH, Q8_SEARCH, MARKING)
        ],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (MANIFEST, CLOSEOUT, Q6_SEARCH, Q8_SEARCH, MARKING)
        },
    },
    "historical_construction": {
        "q": 8,
        "orbit_index": int(a11_to_2a5["orbit"]),
        "new_fibre_formula": "O + P - 2 F",
        "new_fibre_in_parent": entries(historical_fibre),
        "section_class": entries(historical_section),
        "P_dot_O": int(historical_horizontal["P_dot_O"]),
        "height": historical_horizontal["height"],
        "local_correction": historical_horizontal["local_correction"],
        "vertical_root_coefficients": historical_horizontal[
            "vertical_root_coefficients"
        ],
    },
    "isometry_decomposition": {
        "root_chain_orientations": len(root_isometries),
        "mw_automorphism_order": len(automorphisms),
        "integral_glue_isometries": len(glue_isometries),
        "nef_transports": len(nef_records),
        "distinct_nef_orbits": compatible_orbits,
    },
    "construction_compatible_orbits": orbit_profiles,
    "selection_rule": "minimum MW projection L1, then orbit index",
    "selected": selected,
    "all_integral_glue_transports": transport_records,
    "proof_boundary": (
        "Exact integral lattice/marking certificate. It transports the literal "
        "historical O+P-2F divisor through every root/MW/glue-compatible, "
        "fibre-preserving A11 frame isometry and matches the resulting nef rays "
        "against the exhaustive equation-side q8 search. Orbit12 is selected "
        "by an explicit coordinate-complexity convention. The section function, "
        "resolved RR space, quartic, and characteristic-zero Weierstrass child "
        "remain to be constructed."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q24A11Q8CONSTRUCTION|"
    f"root_orientations={len(root_isometries)}|mw_aut={len(automorphisms)}|"
    f"glue={len(glue_isometries)}|nef={len(nef_records)}|"
    f"orbits={','.join(map(str, compatible_orbits))}|selected=12|"
    "formula=O+P-2F|vertical_root=0|status=PASS",
    flush=True,
)
print(f"OUTPUT|{args.output}", flush=True)
