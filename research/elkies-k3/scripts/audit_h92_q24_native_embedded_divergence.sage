#!/usr/bin/env sage -python
"""
Locate the first divergence between the q24-native suffix and the selected
R17-directed H3 corridor by comparing the actual elliptic fibre classes in the
common H3 Neron--Severi lattice.

No qfisom is used.

Both constructions carry exact 19x19 basis maps T satisfying
    G_stage = T * G_H3 * T^t.
Therefore row 0 of T is the stage fibre class in H3 coordinates.  Equality
(up to simultaneous sign) is an exact test that the two stages define the
same elliptic fibration, much stronger and cheaper than abstract positive-frame
isometry.

At the first divergence, the script pulls the historical/R17-directed next
fibre into the last matching native parent coordinates and prints its exact
(a,b,w), q=a*b, and old-fibre degree b.  This is the correct repair witness.
"""

import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import ZZ, block_diagonal_matrix, identity_matrix, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "elkies-k3" / "scripts"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3" / "q24-native-suffix"
GEN = ROOT / "artifacts" / "generated-results"

ENGINE = SCRIPTS / "exact_neighbor_engine.sage"
CLOSE = SCRIPTS / "close_h92_q8_q24_by_q6_translation.sage"
SUFFIX = LOCAL / "q24-native-d12-to-a1.json"
REVERSE = GEN / "elkies-k3-rank17-to-h3-reverse-transport.json"
OUT = LOCAL / "q24-native-embedded-divergence.json"

U2 = matrix(ZZ, ((0,1),(1,0)))

STAGE_IDS = [
    ("D12", "H3-03-D12"),
    ("A11", "H3-04-A11"),
    ("A5+A5", "H3-05-2A5"),
    ("A3+A3+A3", "H3-06-3A3"),
    ("A2+A2+A3", "H3-07-A3-2A2"),
    ("A1+A1+A1+A1+A1", "H3-08-5A1"),
    ("A1+A1+A1+A1", "H3-09-4A1"),
    ("A1+A1+A1", "H3-10-3A1"),
    ("A1+A1", "H3-11-2A1"),
    ("A1", "H3-12-A1"),
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


def same_ray(a, b):
    a = vector(ZZ, a)
    b = vector(ZZ, b)
    if a == b:
        return 1
    if a == -b:
        return -1
    return 0


for p in (ENGINE, CLOSE, SUFFIX, REVERSE):
    if not p.exists():
        raise SystemExit(f"missing prerequisite: {p}")

# -------------------------------------------------------------------------
# Native H3 -> D12 from D24eq.
# -------------------------------------------------------------------------
eng = {"__name__":"__embedded_div_engine__","__file__":str(ENGINE)}
exec(compile(ENGINE.read_text(), str(ENGINE), "exec"), eng)
splitter = eng["primitive_hyperbolic_split"]
minimize = eng["minimize_child_frame"]

saved = list(sys.argv)
scope = {"__name__":"__embedded_div_close__","__file__":str(CLOSE)}
buf = io.StringIO()
try:
    sys.argv = [str(CLOSE)]
    with contextlib.redirect_stdout(buf):
        exec(compile(CLOSE.read_text(), str(CLOSE), "exec"), scope)
finally:
    sys.argv = saved

G_H3 = matrix(ZZ, scope["ns"])
D24 = vector(ZZ, scope["D24eq"])
assert D24 * G_H3 * D24 == 0

sp = splitter(G_H3, D24)
raw = matrix(ZZ, sp["child_frame"])
mi = minimize(raw)
D12 = matrix(ZZ, mi["frame"])
T = (
    block_diagonal_matrix(identity_matrix(ZZ,2), matrix(ZZ, mi["basis"]))
    * matrix(ZZ, sp["transport"])
)
assert T * G_H3 * T.transpose() == ns(D12)
assert vector(ZZ,T.row(0)) == D24

native = [{
    "name":"D12",
    "frame":D12,
    "T":T,
    "incoming": {
        "q":24,
        "native_orbit":"D24eq",
    },
}]

# -------------------------------------------------------------------------
# Native D12 -> A1 suffix, preserving every exact H3-coordinate transport.
# -------------------------------------------------------------------------
suffix = json.loads(SUFFIX.read_text())
assert suffix["status"] == "PASS_Q24_NATIVE_D12_TO_A1"
current = D12

for item in suffix["steps"]:
    artifact = ROOT / item["search_artifact"]
    data = json.loads(artifact.read_text())
    assert data["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
    assert load_gram(ROOT / data["frame"]) == current

    oi = int(item["discovered_orbit_index"])
    hits = [r for r in data["neighbors"] if int(r["orbit_index"]) == oi]
    assert len(hits) == 1
    rec = hits[0]

    raw = matrix(ZZ, rec["child_frame"])
    adapt = matrix(ZZ, rec["child_root_adapted_basis"])
    child = matrix(ZZ, rec["child_root_adapted_frame"])
    neighbor = matrix(ZZ, rec["neighbor_basis"])
    assert adapt * raw * adapt.transpose() == child

    stepT = block_diagonal_matrix(identity_matrix(ZZ,2), adapt) * neighbor
    assert stepT * ns(current) * stepT.transpose() == ns(child)

    T = stepT * T
    current = child
    assert T * G_H3 * T.transpose() == ns(current)

    native.append({
        "name": item["target_ade"],
        "frame": current,
        "T": T,
        "incoming": {
            "q": int(item["q"]),
            "orbit": oi,
            "old_fibre_degree": 2,
        },
    })

assert len(native) == 10

# -------------------------------------------------------------------------
# Historical selected corridor, already stored in common H3 coordinates.
# -------------------------------------------------------------------------
rev = json.loads(REVERSE.read_text())
assert rev["status"] == "PASS_EXACT_PINNED_R17_TO_H3_REVERSE_TRANSPORT"
if "stages" not in rev:
    raise SystemExit("reverse transport artifact has no 'stages' ledger")

hist_by_id = {s["stage_id"]: s for s in rev["stages"]}

records = []
first_div = None
last_match_index = None

for i, ((expected_name, stage_id), nstage) in enumerate(zip(STAGE_IDS, native)):
    assert nstage["name"] == expected_name, (i, nstage["name"], expected_name)
    if stage_id not in hist_by_id:
        raise SystemExit(f"reverse ledger missing stage {stage_id}")
    hstage = hist_by_id[stage_id]

    Tn = matrix(ZZ, nstage["T"])
    Th = matrix(ZZ, hstage["stage_basis_in_h3_ns"])
    assert Tn.dimensions() == Th.dimensions() == (19,19)

    Fn = vector(ZZ, Tn.row(0))
    Fh = vector(ZZ, Th.row(0))
    ray = same_ray(Fn, Fh)
    fibre_match = bool(ray)

    # Exact coordinate-change between the two stage bases.  This is always an
    # integral NS isometry because both are bases of the same H3 NS.
    X = Th * Tn.inverse()
    assert X.change_ring(ZZ) == X
    X = X.change_ring(ZZ)
    assert abs(ZZ(X.det())) == 1
    e0 = vector(ZZ, [1] + [0]*18)
    preserves_fibre = (
        vector(ZZ, X.row(0)) == e0
        or vector(ZZ, X.row(0)) == -e0
    )
    assert preserves_fibre == fibre_match

    same_standard_u = (
        X[:2,:] == identity_matrix(ZZ,19)[:2,:]
    )

    if fibre_match:
        last_match_index = i
    elif first_div is None:
        first_div = {
            "index": i,
            "stage": expected_name,
            "stage_id": stage_id,
        }

    records.append({
        "index":i,
        "stage":expected_name,
        "stage_id":stage_id,
        "native_incoming":nstage["incoming"],
        "same_fibre_class":fibre_match,
        "fibre_sign":ray,
        "same_standard_U":bool(same_standard_u),
        "native_fibre_in_h3":list(map(int,Fn)),
        "historical_fibre_in_h3":list(map(int,Fh)),
        "native_to_historical_stage_ns":rows(X),
    })

    print(
        "Q24EMBED|"
        f"stage={expected_name}|stage_id={stage_id}|"
        f"same_fibre={int(fibre_match)}|same_U={int(same_standard_u)}|"
        f"native_incoming_q={nstage['incoming'].get('q')}|"
        f"native_orbit={nstage['incoming'].get('orbit',nstage['incoming'].get('native_orbit'))}|"
        f"status={'MATCH' if fibre_match else 'DIVERGED'}",
        flush=True,
    )

# -------------------------------------------------------------------------
# Repair witness: historical next fibre in last matching native parent coords.
# -------------------------------------------------------------------------
repair = None

if first_div is not None:
    j = first_div["index"]

    if j == 0:
        # Divergence already at D12.  The repair target is historical D12's
        # q24 fibre expressed directly in H3 coordinates; D24eq is the native
        # q24 fibre, so report both.  A D13-equation-frame conversion would be
        # a separate presentation issue.
        hnext = hist_by_id["H3-03-D12"]
        Ftarget_h3 = vector(ZZ, matrix(ZZ,hnext["stage_basis_in_h3_ns"]).row(0))
        repair = {
            "parent":"H3 source / equation-side D13 presentation",
            "target_stage":"D12",
            "target_fibre_in_h3":list(map(int,Ftarget_h3)),
            "native_D24eq_in_h3":list(map(int,D24)),
            "same_ray_as_D24eq":bool(same_ray(Ftarget_h3,D24)),
            "note":(
                "The selected historical/R17-directed q24 fibre itself differs "
                "from D24eq in H3 coordinates; repair must occur at the q24 hop."
            ),
        }
        print(
            "Q24EMBED_REPAIR|parent=H3|target=D12|"
            f"same_as_D24eq={int(bool(same_ray(Ftarget_h3,D24)))}|"
            "status=REPAIR_Q24",
            flush=True,
        )
    else:
        parent_native = native[j-1]
        parent_name = STAGE_IDS[j-1][0]
        target_name = STAGE_IDS[j][0]
        target_id = STAGE_IDS[j][1]

        Tparent = matrix(ZZ, parent_native["T"])
        Thnext = matrix(ZZ, hist_by_id[target_id]["stage_basis_in_h3_ns"])
        Ftarget_h3 = vector(ZZ, Thnext.row(0))

        # Row coordinates in native parent basis.
        Ftarget_native = vector(ZZ, Ftarget_h3 * Tparent.inverse())
        if Ftarget_native[1] < 0:
            Ftarget_native = -Ftarget_native

        Gp = ns(parent_native["frame"])
        Fold = vector(ZZ, [1,0] + [0]*17)
        assert Ftarget_native * Gp * Ftarget_native == 0

        degree = ZZ(Ftarget_native * Gp * Fold)
        a = ZZ(Ftarget_native[0])
        b = ZZ(Ftarget_native[1])
        w = vector(ZZ, Ftarget_native[2:])
        q = a*b
        assert degree == b
        assert w * parent_native["frame"] * w == 2*q

        repair = {
            "parent_stage":parent_name,
            "target_stage":target_name,
            "historical_target_stage_id":target_id,
            "target_fibre_in_h3":list(map(int,Ftarget_h3)),
            "target_fibre_in_native_parent":list(map(int,Ftarget_native)),
            "a":int(a),
            "b":int(b),
            "q":int(q),
            "old_fibre_degree":int(degree),
            "witness":list(map(int,w)),
            "witness_norm":int(w*parent_native["frame"]*w),
            "primitive":bool(__import__("math").gcd(*[abs(int(x)) for x in Ftarget_native]) == 1),
        }

        print(
            "Q24EMBED_REPAIR|"
            f"parent={parent_name}|target={target_name}|"
            f"a={a}|b={b}|q={q}|degree={degree}|"
            f"witness_norm={w*parent_native['frame']*w}|"
            "status=PASS_EXACT_REPAIR_FIBRE",
            flush=True,
        )

last_match = None if last_match_index is None else STAGE_IDS[last_match_index][0]
first_name = None if first_div is None else first_div["stage"]

payload = {
    "schema":"elkies-k3.h3-q24-native-embedded-divergence.v1",
    "status":"PASS_Q24_NATIVE_EMBEDDED_DIVERGENCE_AUDIT",
    "comparison_basis":"common H3 Neron-Severi coordinates",
    "last_matching_stage":last_match,
    "first_divergent_stage":first_name,
    "stages":records,
    "repair":repair,
    "interpretation":(
        "same_fibre_class compares the actual primitive isotropic fibre rays "
        "inside the common H3 NS; it is not an abstract frame-isometry test. "
        "The repair record is the selected R17-directed next fibre pulled "
        "exactly into the last matching native parent coordinates."
    ),
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24EMBED_RESULT|"
    f"last_match={last_match}|first_divergence={first_name}|"
    "status=PASS_Q24_NATIVE_EMBEDDED_DIVERGENCE_AUDIT",
    flush=True,
)
