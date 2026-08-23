#!/usr/bin/env sage -python
"""
Close the R17-directed lattice/marking corridor from the CURRENT equation D13
frame through q24 orbit85 and the certified historical suffix to pinned R17.

Critical idea
-------------
The pinned/historical D13 frame and the current equation D13 frame are related
by a nontrivial exact NS marking change.  Therefore their q24 fibre classes
must not be compared as literal H3-coordinate rays.

We reconstruct the exact 19x19 coordinate map A such that

    A * G_equation_D13 * A^t = G_historical_D13,

using the current passing dominant->nef->physical->Weyl->translation transport.

Then the historical orbit85 q24 child basis H24 gives

    E24 = H24 * A,

and hence

    E24 * G_equation_D13 * E24^t
      = U + (- historical_D12_frame).

Thus the q24 child is deliberately expressed in the already-certified
R17-directed D12 coordinates.  From there every stored suffix transition
replays literally, including

    historical A1/MW16 --q6 orbit2247--> rootless/MW17,

followed by the pinned endpoint isometry.

This is an exact integral lattice/NS/marking certificate.  It does not execute
the downstream characteristic-zero Weierstrass pencils.
"""

import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import (
    ZZ, block_diagonal_matrix, identity_matrix, matrix, vector
)

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "elkies-k3" / "scripts"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"
GEN = ROOT / "artifacts" / "generated-results"

PROFILE = SCRIPTS / "compare_h92_three_q24_d12_current_equation_profiles.sage"
Q24_ART = GEN / "elkies-k3-h3-q6-q8-d13-q24-degree2.json"
PINNED = ROOT / "elkies-k3" / "data" / "lattice" / "rank17_gram.txt"
ENDPOINT_ISO = (
    ROOT / "elkies-k3" / "data" / "fibrations"
    / "h3_rootless_mw17_to_pinned_rank17_isometry.txt"
)
OUT = LOCAL / "q24-equation-d13-to-pinned-r17.json"

U2 = matrix(ZZ, ((0,1),(1,0)))

SUFFIX = [
    # artifact, orbit, q, root_data, ADE, MW
    ("elkies-k3-h3-d12-o85-q6-degree2.json",42,6,(11,132,12),"A11",6),
    ("elkies-k3-h3-a11-middle-q8-degree2.json",922,8,(10,60,36),"A5+A5",7),
    ("elkies-k3-h3-a5a5-c2-q4-degree2.json",472,4,(9,36,64),"A3+A3+A3",8),
    ("elkies-k3-h3-a3x3-q4-degree2.json",323,4,(7,24,36),"A2+A2+A3",10),
    ("elkies-k3-h3-mw10-a3a2a2-q4-degree2.json",207,4,(5,10,32),"A1+A1+A1+A1+A1",12),
    ("elkies-k3-h3-mw12-5a1-q4-degree2-first-hit.json",52,4,(4,8,16),"A1+A1+A1+A1",13),
    ("elkies-k3-h3-mw13-4a1-q4-degree2-first-hit.json",114,4,(3,6,8),"A1+A1+A1",14),
    ("elkies-k3-h3-mw14-3a1-q4-degree2-first-hit.json",498,4,(2,4,4),"A1+A1",15),
    ("elkies-k3-h3-mw15-2a1-q4-degree2-first-hit.json",981,4,(1,2,2),"A1",16),
    ("elkies-k3-h3-mw16-a1-q6-degree2-cap10000-stream-chunk001.json",2247,6,(0,0,1),"rootless",17),
]


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(x) for x in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


def ns(frame):
    return block_diagonal_matrix(U2, -frame)


def rows(M):
    return [[int(x) for x in row] for row in M.rows()]


def run_scope(path, argv=()):
    saved = list(sys.argv)
    scope = {"__name__":"__embedded_profile__","__file__":str(path)}
    buf = io.StringIO()
    try:
        sys.argv = [str(path)] + list(argv)
        with contextlib.redirect_stdout(buf):
            exec(compile(path.read_text(), str(path), "exec"), scope)
    finally:
        sys.argv = saved
    # Surface only the useful profile lines.
    for line in buf.getvalue().splitlines():
        if (
            line.startswith("Q24ALL3CURRENT_PROFILE|")
            or line.startswith("Q24ALL3CURRENT_REGRESSION|")
            or line.startswith("Q24ALL3CURRENT_RESULT|")
        ):
            print(line, flush=True)
    return scope


for p in [PROFILE,Q24_ART,PINNED,ENDPOINT_ISO] + [GEN/x[0] for x in SUFFIX]:
    if not p.exists():
        raise SystemExit(f"missing prerequisite: {p}")

# -------------------------------------------------------------------------
# 1. Current exact equation-D13 marking and historical->equation map.
# -------------------------------------------------------------------------
sc = run_scope(PROFILE)

need = (
    "Gpinned","ns","adapted","Badapt","coords","profiles","F8eq",
    "Bpinned_to_simple","q6_simple_ns","simple_root_classes","bridge",
    "Bsimple","source_ns","actual_roots","physical_reflections","reflect",
    "weyl_transport","tau",
)
missing = [k for k in need if k not in sc]
if missing:
    raise SystemExit("profile scope missing: " + ",".join(missing))

Ghist = matrix(ZZ, sc["Gpinned"])
Geq_frame = matrix(ZZ, sc["adapted"])
Geq = ns(Geq_frame)
Badapt = matrix(ZZ, sc["Badapt"])
coords = sc["coords"]
Bpinned_to_simple = matrix(ZZ, sc["Bpinned_to_simple"])
q6_simple_ns = matrix(ZZ, sc["q6_simple_ns"])
simple_root_classes = tuple(vector(ZZ,r) for r in sc["simple_root_classes"])
bridge = tuple(sc["bridge"])
Bsimple = matrix(ZZ, sc["Bsimple"])
source_ns = matrix(ZZ, sc["source_ns"])
actual_roots = tuple(vector(ZZ,r) for r in sc["actual_roots"])
physical_reflections = tuple(sc["physical_reflections"])
reflect = sc["reflect"]
weyl_transport = sc["weyl_transport"]
tau = sc["tau"]

assert Ghist.dimensions() == Geq.dimensions() == (19,19)

def historical_d13_to_equation_ambient(vh):
    # Unrestricted linear version of the current passing marking transport.
    # The older helpers deliberately assert q24-specific degree-two conditions.
    D = vector(ZZ, vh) * Bpinned_to_simple

    # dominant -> nef
    for i, unused_pairing in reversed(bridge):
        D = reflect(D, q6_simple_ns, simple_root_classes[i])

    # q6-simple -> source H3
    D = vector(ZZ, D * Bsimple)

    # exact physical component reflections
    for i, unused_pairing in physical_reflections:
        D = reflect(D, source_ns, actual_roots[i])

    # current q6 Weyl chamber -> equation zero normalization
    D = vector(ZZ, weyl_transport(D))
    D = vector(ZZ, tau(D))
    return D
assert abs(ZZ(Ghist.det())) == abs(ZZ(Geq.det()))

# Rows of A are historical-D13 basis vectors expressed in equation-D13 coords.
basis = identity_matrix(ZZ,19)
Arows = []
for i in range(19):
    vh = vector(ZZ, basis.row(i))
    ambient = historical_d13_to_equation_ambient(vh)
    ceq = vector(ZZ, coords(ambient))
    Arows.append(list(ceq))

A = matrix(ZZ, Arows)
assert abs(ZZ(A.det())) == 1
assert A * Geq * A.transpose() == Ghist

print(
    "Q24R17MAP_D13|"
    f"det={A.det()}|historical_basis_in_equation=1|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 2. Historical orbit85 q24 -> D12, transported into equation-D13 coords.
# -------------------------------------------------------------------------
q24data = json.loads(Q24_ART.read_text())
assert q24data["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
hits = [
    r for r in q24data["neighbors"]
    if int(r["orbit_index"]) == 85
    and int(r["q"]) == 24
    and tuple(r["child_root_data"]) == (12,264,4)
    and r["child_ade"] == "D12"
    and int(r["child_mw_rank"]) == 5
]
assert len(hits) == 1
r24 = hits[0]

H24 = (
    block_diagonal_matrix(
        identity_matrix(ZZ,2),
        matrix(ZZ,r24["child_root_adapted_basis"]),
    )
    * matrix(ZZ,r24["neighbor_basis"])
)
D12 = matrix(ZZ,r24["child_root_adapted_frame"])
assert H24 * Ghist * H24.transpose() == ns(D12)
assert abs(ZZ(H24.det())) == 1

# Equation-D13 -> historical/R17-directed D12 basis.
E24 = H24 * A
assert abs(ZZ(E24.det())) == 1
assert E24 * Geq * E24.transpose() == ns(D12)

f_hist = vector(ZZ,r24["fiber"])
f_eq = vector(ZZ, f_hist * A)
assert f_eq == vector(ZZ,E24.row(0))
assert f_eq * Geq * f_eq == 0
assert f_eq[1] == 2

# Independent equation-profile regression: orbit85 divisor from current
# transport stack must equal this mapped fibre.
p85 = next(x for x in sc["profiles"] if int(x["orbit"]) == 85)
Deq_ambient = vector(ZZ,p85["equation_divisor"])
Deq_coords = vector(ZZ,coords(Deq_ambient))
assert Deq_coords == f_eq

print(
    "Q24R17MAP_Q24|"
    f"orbit=85|q=24|degree={f_eq[1]}|"
    "equation_profile_match=1|child=D12/MW5|"
    f"transport_det={E24.det()}|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 3. Replay the certified R17-directed suffix literally from this D12 basis.
# -------------------------------------------------------------------------
current = D12
T = E24
steps = [{
    "stage":"D12",
    "q":24,
    "orbit":85,
    "root_data":[12,264,4],
    "mw_rank":5,
    "transition":rows(E24),
}]

a1_frame = None
a1_to_rootless = None

for index,(name,orbit,q,rdata,ade,mw) in enumerate(SUFFIX, start=1):
    path = GEN / name
    data = json.loads(path.read_text())
    assert data["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    assert load_gram(ROOT/data["frame"]) == current

    matches = [r for r in data["neighbors"] if int(r["orbit_index"]) == orbit]
    assert len(matches) == 1
    rec = matches[0]
    assert int(rec["q"]) == q
    assert int(rec["old_fiber_degree"]) == 2
    assert tuple(rec["child_root_data"]) == rdata
    assert rec["child_ade"] == ade
    assert int(rec["child_mw_rank"]) == mw

    child = matrix(ZZ,rec["child_root_adapted_frame"])
    stepT = (
        block_diagonal_matrix(
            identity_matrix(ZZ,2),
            matrix(ZZ,rec["child_root_adapted_basis"]),
        )
        * matrix(ZZ,rec["neighbor_basis"])
    )
    assert abs(ZZ(stepT.det())) == 1
    assert stepT * ns(current) * stepT.transpose() == ns(child)

    if ade == "A1":
        a1_frame = child

    if ade == "rootless":
        assert a1_frame is not None
        a1_to_rootless = stepT
        assert tuple(rec["child_root_data"]) == (0,0,1)

    T = stepT * T
    current = child

    steps.append({
        "stage":ade,
        "q":q,
        "orbit":orbit,
        "old_fibre_degree":2,
        "root_data":list(rdata),
        "mw_rank":mw,
        "witness":rec["witness"],
        "transition":rows(stepT),
    })

    print(
        "Q24R17MAP_SUFFIX|"
        f"step={index}|q={q}|orbit={orbit}|ADE={ade}|MW={mw}|status=PASS",
        flush=True,
    )

assert a1_frame is not None and a1_to_rootless is not None
rootless = current
assert T * Geq * T.transpose() == ns(rootless)

print(
    "Q24R17MAP_FINAL_Q6|"
    "parent=A1/MW16|q=6|orbit=2247|degree=2|"
    "child=rootless/MW17|status=PASS",
    flush=True,
)

# -------------------------------------------------------------------------
# 4. Pin historical rootless endpoint to recovered 17x17 R17.
# -------------------------------------------------------------------------
pinned = load_gram(PINNED)
C = load_gram(ENDPOINT_ISO)
assert C.transpose() * pinned * C == rootless
assert C.det() == 1

# pinned basis in historical rootless coords
P2R = block_diagonal_matrix(identity_matrix(ZZ,2), C.transpose())
assert P2R * ns(pinned) * P2R.transpose() == ns(rootless)

# rootless -> pinned
R2P = P2R.inverse()
assert R2P.change_ring(ZZ) == R2P
R2P = R2P.change_ring(ZZ)

TOTAL = R2P * T
assert abs(ZZ(TOTAL.det())) == 1
assert TOTAL * Geq * TOTAL.transpose() == ns(pinned)

# Isolate the final A1 -> pinned-R17 map as requested.
A1_FINAL = R2P * a1_to_rootless
assert abs(ZZ(A1_FINAL.det())) == 1
assert A1_FINAL * ns(a1_frame) * A1_FINAL.transpose() == ns(pinned)

print(
    "Q24R17MAP_PINNED|"
    f"final_transport_det={TOTAL.det()}|"
    f"a1_to_pinned_det={A1_FINAL.det()}|"
    "pinned_R17=1|status=PASS",
    flush=True,
)

payload = {
    "schema":"elkies-k3.h3-q24-equation-d13-to-pinned-r17.v1",
    "status":"PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH",
    "proof_boundary":(
        "Exact integral lattice/NS/marking certificate. The current equation "
        "D13 marking is related to the historical selected D13 marking by the "
        "fully reconstructed current transport stack. Orbit85 is transported "
        "into equation-D13 coordinates and its child basis is deliberately "
        "chosen as the certified R17-directed D12 basis. The stored suffix "
        "then replays literally through the final q6 orbit2247 and the exact "
        "rootless-to-pinned-R17 isometry. This does not execute the downstream "
        "characteristic-zero Weierstrass pencils."
    ),
    "equation_d13_frame":rows(Geq_frame),
    "historical_d13_basis_in_equation_d13":rows(A),
    "q24":{
        "orbit":85,
        "q":24,
        "old_fibre_degree":2,
        "equation_d13_fibre":list(map(int,f_eq)),
        "equation_profile_match":True,
        "child_frame":rows(D12),
        "equation_d13_to_d12_transition":rows(E24),
    },
    "steps":steps,
    "final_a1_frame":rows(a1_frame),
    "final_a1_to_rootless_q6_transition":rows(a1_to_rootless),
    "final_a1_to_pinned_r17_transition":rows(A1_FINAL),
    "rootless_frame":rows(rootless),
    "pinned_rank17_frame":str(PINNED.relative_to(ROOT)),
    "equation_d13_to_pinned_r17_transition":rows(TOTAL),
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24R17MAP_RESULT|"
    "q24_orbit85=1|D12_R17_directed=1|"
    "A1_q6_orbit2247=1|pinned_R17=1|"
    "status=PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH",
    flush=True,
)
