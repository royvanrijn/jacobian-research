#!/usr/bin/env sage -python
"""
Extract the selected historical orbit42 D12 -> A11 divisor/marked section
directly in the CURRENT equation-D13 coordinates using the exact
q24-equation-D13-to-pinned-R17 certificate.

No MW search, polynomial-section enumeration, q32 marking, or arbitrary
q24-native D12 split is used.

Outputs:
  * exact historical D12 zero O12 in equation-D13 coordinates;
  * exact orbit42 fibre D42 in D12 and equation-D13 coordinates;
  * exact orbit42 MW class z=(-1,0,-1,-1,0);
  * a convenient section representative P42 with correction 1, P.O=2;
  * exact decomposition D42 = O12 + P42 + V + F12;
  * q8/D13 degrees and MW projections of O12 and P42, to choose the
    equation-level recovery method.
"""
import argparse
import json
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

MAP = LOCAL / "q24-equation-d13-to-pinned-r17.json"
O42ART = GEN / "elkies-k3-h3-d12-o85-q6-degree2.json"
for path in (MAP, O42ART):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

mp = json.loads(MAP.read_text())
q6 = json.loads(O42ART.read_text())

assert mp["status"] == "PASS_Q24_EQUATION_D13_TO_PINNED_R17_LATTICE_PATH"
assert q6["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"

U2 = matrix(ZZ, ((0,1),(1,0)))

def ns(frame):
    return block_diagonal_matrix(U2, -matrix(ZZ, frame))

Geq_frame = matrix(ZZ, mp["equation_d13_frame"])
Geq = ns(Geq_frame)

q24 = mp["q24"]
D12 = matrix(ZZ, q24["child_frame"])
G12 = ns(D12)
E24 = matrix(ZZ, q24["equation_d13_to_d12_transition"])

assert E24.dimensions() == (19,19)
assert abs(ZZ(E24.det())) == 1
assert E24 * Geq * E24.transpose() == G12

# Historical selected orbit42.
hits = [
    r for r in q6["neighbors"]
    if int(r["orbit_index"]) == 42
    and int(r["q"]) == 6
    and int(r["old_fiber_degree"]) == 2
    and tuple(r["child_root_data"]) == (11,132,12)
    and r["child_ade"] == "A11"
    and int(r["child_mw_rank"]) == 6
]
assert len(hits) == 1
rec = hits[0]

# Cross-check the certificate's first suffix transition against the artifact.
step1 = mp["steps"][1]
assert step1["stage"] == "A11"
assert int(step1["q"]) == 6
assert int(step1["orbit"]) == 42
assert tuple(step1["root_data"]) == (11,132,12)
assert int(step1["mw_rank"]) == 6

stepT = matrix(ZZ, step1["transition"])
expected_stepT = (
    block_diagonal_matrix(
        identity_matrix(ZZ,2),
        matrix(ZZ, rec["child_root_adapted_basis"]),
    )
    * matrix(ZZ, rec["neighbor_basis"])
)
assert stepT == expected_stepT
assert abs(ZZ(stepT.det())) == 1

# Standard U objects in the historical/R17-directed D12 basis.
F12 = vector(ZZ, [1,0] + [0]*17)
O12 = vector(ZZ, [-1,1] + [0]*17)
D42 = vector(ZZ, rec["fiber"])
assert D42 == vector(ZZ, stepT.row(0))

assert F12 * G12 * F12 == 0
assert O12 * G12 * O12 == -2
assert O12 * G12 * F12 == 1
assert D42 * G12 * D42 == 0
assert D42 * G12 * F12 == 2

# Positive D12 frame split.
R = D12[:12,:12]
C = D12[:12,12:]
Tail = D12[12:,12:]
H = Tail - C.transpose()*R.inverse()*C

z = vector(ZZ, rec["mw_projection"])
w = vector(ZZ, rec["witness"])
assert vector(ZZ, w[12:]) == z
assert w * D12 * w == 12
height = QQ(z * H * z)
assert height == 7

candidates=[]

# Exact D12 spinor class: local correction 3.
#
# correction=3 is a discriminant-coset statement, not necessarily incidence
# with one simple component.  Put pair=(0,z)G_roots and dual=pair*R^{-1}.
# We need an integral root shift r with
#
#     (r+dual) R (r+dual)^t = 1.
#
# Equivalently lambda=(r+dual)R is an integral dual-lattice vector in the
# required coset and lambda R^{-1} lambda^t=3.
#
# For D12, det(R)=4, so enumerate norm-4 vectors of the integral adjugate
# form 4*R^{-1}, at norm 12; this is exact and tiny.
from sage.all import pari

pair = vector(ZZ, vector(ZZ,[0]*12 + list(z)) * D12[:,:12])
dual = vector(QQ, pair) * R.inverse()

detR = ZZ(R.det())
assert detR == 4
Radj = (detR * R.inverse()).change_ring(ZZ)

qf = pari(Radj).qfminim(ZZ(3)*detR)
Mmin = matrix(ZZ,qf[2])
seen_lam = set()

for col in Mmin.columns():
    for sg in (1,-1):
        lam = sg * vector(ZZ,col)
        key = tuple(map(int,lam))
        if key in seen_lam:
            continue
        seen_lam.add(key)

        if ZZ(lam * Radj * lam) != ZZ(3)*detR:
            continue

        rQ = vector(QQ, lam - pair) * R.inverse()
        if not all(x in ZZ for x in rQ):
            continue
        r = vector(ZZ,[ZZ(x) for x in rQ])

        s = vector(QQ,r) + dual
        corr = QQ(s * R * s)
        assert corr == 3

        pframe = vector(ZZ,list(r)+list(z))
        norm = QQ(pframe * D12 * pframe)
        assert norm == height + corr == 10

        P12 = vector(ZZ,[4,1]+list(pframe))
        assert P12 * G12 * P12 == -2
        assert P12 * G12 * F12 == 1
        assert P12 * G12 * O12 == 3

        residual = D42 - O12 - P12
        assert residual[1] == 0
        assert all(v == 0 for v in residual[14:])
        fibre_twist = ZZ(residual[0])
        root_coeff = vector(ZZ,residual[2:14])

        candidates.append({
            "dual_pairing":list(map(int,lam)),
            "root_coordinates":list(map(int,r)),
            "section_frame":list(map(int,pframe)),
            "section_class":list(map(int,P12)),
            "fibre_twist":int(fibre_twist),
            "vertical_root_coefficients":list(map(int,root_coeff)),
            "vertical_L1":sum(abs(int(v)) for v in root_coeff),
            "vertical_support":sum(bool(v) for v in root_coeff),
        })

print(
    "Q24O42EQ_DISCRIMINANT|"
    f"detR={detR}|correction=3|dual_min_vectors={len(seen_lam)}|"
    f"coset_section_candidates={len(candidates)}|status="
    + ("PASS" if candidates else "FAIL"),
    flush=True,
)

if not candidates:
    raise ArithmeticError(
        "orbit42 exact D12 coset requires correction 3, but no norm-3 representative "
        "was found in the target D12 discriminant coset"
    )


candidates.sort(key=lambda x:(x["vertical_L1"],x["vertical_support"],tuple(x["dual_pairing"])))
best = candidates[0]
P12 = vector(ZZ, best["section_class"])

# Map all relevant classes into current equation-D13 coordinates.
def to_eq(v12):
    return vector(ZZ, vector(ZZ,v12) * E24)

F12eq = to_eq(F12)
O12eq = to_eq(O12)
P12eq = to_eq(P12)
D42eq = to_eq(D42)

assert F12eq == vector(ZZ, q24["equation_d13_fibre"])
assert D42eq == vector(ZZ, (stepT*E24).row(0))

# Current q8/D13 coordinates are U + positive D13 frame.
F8 = vector(ZZ, [1,0] + [0]*17)
O8 = vector(ZZ, [-1,1] + [0]*17)

def q8_profile(name, v):
    v = vector(ZZ,v)
    degree = ZZ(v * Geq * F8)
    odot = ZZ(v * Geq * O8)
    square = ZZ(v * Geq * v)
    frame = vector(ZZ, v[2:])
    root = vector(ZZ, frame[:13])
    mw = vector(ZZ, frame[13:])
    return {
        "name":name,
        "square":int(square),
        "q8_degree":int(degree),
        "q8_O_intersection":int(odot),
        "equation_d13_coordinates":list(map(int,v)),
        "d13_root_coordinates":list(map(int,root)),
        "d13_mw_projection":list(map(int,mw)),
        "frame_L1":sum(abs(int(x)) for x in frame),
        "mw_L1":sum(abs(int(x)) for x in mw),
    }

profiles = [
    q8_profile("F12_q24",F12eq),
    q8_profile("O12_historical_zero",O12eq),
    q8_profile("P42_orbit42_marked",P12eq),
    q8_profile("D42_A11_fibre",D42eq),
]

for pr in profiles:
    print(
        "Q24O42EQ_CLASS|"
        f"name={pr['name']}|square={pr['square']}|"
        f"q8_degree={pr['q8_degree']}|q8_O={pr['q8_O_intersection']}|"
        f"d13_mw={','.join(map(str,pr['d13_mw_projection']))}|"
        f"frame_L1={pr['frame_L1']}|status=PASS",
        flush=True,
    )

print(
    "Q24O42EQ_D12|"
    f"mw={','.join(map(str,z))}|height={height}|corr=3|PdotO=3|"
    f"section_candidates={len(candidates)}|"
    f"chosen_dual_pairing={','.join(map(str,best['dual_pairing']))}|"
    f"vertical_F={best['fibre_twist']}|"
    f"vertical_L1={best['vertical_L1']}|"
    f"vertical_support={best['vertical_support']}|"
    "status=PASS_EXACT_ORBIT42_DECOMPOSITION",
    flush=True,
)

# Decide the cheapest likely recovery route from the q8 equation.
Opr = next(x for x in profiles if x["name"]=="O12_historical_zero")
Ppr = next(x for x in profiles if x["name"]=="P42_orbit42_marked")

if Opr["q8_degree"] == 0:
    zero_method = "OLD_D13_FIBRE_COMPONENT"
elif Opr["q8_degree"] == 1:
    zero_method = "Q8_SECTION"
else:
    zero_method = "Q8_MULTI-SECTION"

if Ppr["q8_degree"] == 0:
    marked_method = "OLD_D13_FIBRE_COMPONENT"
elif Ppr["q8_degree"] == 1:
    marked_method = "Q8_SECTION"
else:
    marked_method = "Q8_MULTI-SECTION"

print(
    "Q24O42EQ_FRONTIER|"
    f"zero_method={zero_method}|marked_method={marked_method}|"
    f"zero_q8_degree={Opr['q8_degree']}|marked_q8_degree={Ppr['q8_degree']}|"
    "next=RECOVER_MAPPED_O12_P42_ON_Q24_QUARTIC_THEN_COMPILE_Q6|status=ACTIONABLE",
    flush=True,
)

payload = {
    "schema":"elkies-k3.h3-q24-orbit42-current-equation-bridge.v1",
    "status":"PASS_Q24_ORBIT42_CURRENT_EQUATION_LATTICE_BRIDGE",
    "source_certificate":str(MAP.relative_to(ROOT)),
    "orbit42_artifact":str(O42ART.relative_to(ROOT)),
    "route_end":"pinned R17",
    "D12":{
        "frame":[[int(x) for x in row] for row in D12.rows()],
        "F":list(map(int,F12)),
        "O":list(map(int,O12)),
        "orbit42_fibre":list(map(int,D42)),
        "orbit42_mw_projection":list(map(int,z)),
        "orbit42_height":str(height),
        "orbit42_local_correction":"3",
        "orbit42_P_dot_O":3,
        "section_representatives":candidates,
        "selected_section_representative":best,
    },
    "current_equation_D13":{
        "q24_fibre":list(map(int,F12eq)),
        "historical_D12_zero":list(map(int,O12eq)),
        "orbit42_marked_section":list(map(int,P12eq)),
        "A11_fibre":list(map(int,D42eq)),
        "profiles":profiles,
    },
    "next":{
        "zero_recovery_method":zero_method,
        "marked_recovery_method":marked_method,
        "instruction":(
            "Recover the mapped historical D12 zero and orbit42 marked section "
            "on the already-certified q24 genus-one/quartic model, then compile "
            "the degree-two q6 pencil and require A11 root data (11,132,12)."
        ),
    },
    "proof_boundary":(
        "Exact integral lattice/NS/marking bridge only. It eliminates the D12 "
        "marking ambiguity and provides the selected orbit42 divisor and section "
        "classes in current equation-D13 coordinates. It does not yet supply "
        "their rational functions on the q24 D12 equation."
    ),
}

OUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-orbit42-current-equation-bridge.json"
)
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24O42EQ_RESULT|"
    "orbit85_current_equation=1|orbit42_selected=1|R17_directed=1|"
    "status=PASS_Q24_ORBIT42_CURRENT_EQUATION_LATTICE_BRIDGE",
    flush=True,
)
