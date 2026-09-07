#!/usr/bin/env sage
"""
Trace the retained Q80 low-q suffix as MARKED lattice objects.

No equation search is performed. This replays the exact retained neighbor
vectors, keeps every unimodular basis transport, reconstructs the final q6
horizontal section on candidate1 from the chamber-reduced divisor, and
transports that exact (-2)-class backwards to the original E6+D5+A3/MW3
frame.

Conventions:
  current full NS Gram = U + (-frame)
  U basis = (F, O+F), hence
      F = (1,0,0...)
      O = (-1,1,0...)
  neighbor(...) returns T whose rows are the new full NS basis expressed in
  the old full NS basis. Thus row coordinates transform by
      class_old = class_new * T.
"""

from pathlib import Path
import csv
import json

from sage.all import (
    ZZ, QuadraticForm, block_diagonal_matrix, gcd,
    identity_matrix, matrix, vector, xgcd
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "elkies-k3" / "data" / "fibrations"
OUT = ROOT / "artifacts" / "generated-results"
OUT.mkdir(parents=True, exist_ok=True)

U = matrix(ZZ, [[0,1],[1,0]])

def load_matrix(path):
    return matrix(
        ZZ,
        [[ZZ(x) for x in line.split()]
         for line in Path(path).read_text().splitlines()
         if line.strip() and not line.lstrip().startswith("#")]
    )

def bezout_vector(pairings):
    current = ZZ(0)
    coeffs = [ZZ(0)] * len(pairings)
    for i, pairing in enumerate(pairings):
        if not pairing:
            continue
        g, left, right = xgcd(current, ZZ(pairing))
        coeffs = [left*c for c in coeffs]
        coeffs[i] += right
        current = g
    assert abs(current) == 1
    return vector(ZZ, coeffs if current == 1 else [-c for c in coeffs])

def neighbor(parent, qnorm, a, b, coordinates):
    ns = block_diagonal_matrix(U, -parent)
    coordinates = vector(ZZ, coordinates)
    Fnew = vector(ZZ, [a,b] + list(coordinates))
    assert a*b == qnorm
    assert coordinates * parent * coordinates == 2*qnorm
    assert Fnew * ns * Fnew == 0
    assert gcd([abs(ZZ(x)) for x in ns*Fnew]) == 1

    mate = bezout_vector(list(ns*Fnew))
    mate -= ZZ(mate*ns*mate)//2 * Fnew

    complement = matrix(
        ZZ, [list(Fnew*ns), list(mate*ns)]
    ).right_kernel_matrix()

    child = -(complement * ns * complement.transpose())
    T = matrix(ZZ, [list(Fnew), list(mate)] + complement.rows())
    assert abs(T.det()) == 1

    child_ns = block_diagonal_matrix(U, -child)
    assert T * ns * T.transpose() == child_ns
    return child, T

def full_ns(frame):
    return block_diagonal_matrix(U, -frame)

def standard_F(frame):
    return vector(ZZ, [1,0] + [0]*frame.nrows())

def standard_O(frame):
    return vector(ZZ, [-1,1] + [0]*frame.nrows())

def intersection(x,y,frame):
    return ZZ(vector(ZZ,x) * full_ns(frame) * vector(ZZ,y))

def roots_norm_two(frame):
    half = QuadraticForm(ZZ, frame).short_vector_list_up_to_length(
        2, up_to_sign_flag=True
    )[1]
    half = [vector(ZZ,row) for row in half]
    return half + [-r for r in half]

def lex_positive(row):
    return next(v > 0 for v in row if v)

def a1_root(frame):
    roots = roots_norm_two(frame)
    positive = [r for r in roots if lex_positive(r)]
    assert len(positive) == 1, len(positive)
    r = positive[0]
    assert r * frame * r == 2
    return r

def chamber_reduce_a1(raw, frame):
    ns = full_ns(frame)
    F = standard_F(frame)
    O = standard_O(frame)
    r = a1_root(frame)
    R = vector(ZZ, [0,0] + list(r))
    Theta0 = F - R

    curves = (("O",O),("R1",R),("Theta0_1",Theta0))
    D = vector(ZZ, raw)
    sequence = []

    while True:
        changed = False
        for name,C in curves:
            pairing = ZZ(D * ns * C)
            if pairing < 0:
                D += pairing*C
                sequence.append((name,int(pairing)))
                assert D * ns * D == 0
                changed = True
                break
        if not changed:
            return D, tuple(sequence), r

def vec(v):
    return [int(x) for x in v]

def mat_rows(M):
    return [[int(x) for x in row] for row in M.rows()]

# 1. Original common prefix, retaining transports.
with (DATA / "kumar_q80_to_rootless_path.tsv").open() as handle:
    canonical = list(csv.DictReader(handle, delimiter="\t"))
assert len(canonical) >= 2

start = load_matrix(DATA / "kumar_q80_e6_d5_a3_mw3_frame.txt")
assert start.nrows() == 17 and start.ncols() == 17
assert start.det() == 948

first = canonical[0]
frame1, T_prefix1 = neighbor(
    start,
    ZZ(first["q"]), ZZ(first["a"]), ZZ(first["b"]),
    vector(ZZ, [ZZ(x) for x in first["v"].split(",")]),
)

second = canonical[1]
frame2, T_prefix2 = neighbor(
    frame1,
    ZZ(second["q"]), ZZ(second["a"]), ZZ(second["b"]),
    vector(ZZ, [ZZ(x) for x in second["v"].split(",")]),
)

print(
    "Q80MARKTRACEPREFIX|T1_det={}|T2_det={}|"
    "status=PASS_RETAIN_COMMON_PREFIX_TRANSPORTS".format(
        int(T_prefix1.det()), int(T_prefix2.det())
    ),
    flush=True,
)

# 2. Complete retained low-q path from D7+D5/MW5.
with (DATA / "kumar_q80_new_lowq_rootless_path.tsv").open() as handle:
    steps = list(csv.DictReader(handle, delimiter="\t"))
assert len(steps) == 8

expected_children = (
    "D7+D4/MW6",
    "A6+A4/MW7",
    "A6+A3/MW8",
    "A4+A2+A1/MW10",
    "A3+A2/MW12",
    "4A1/MW13",
    "A1/MW16",
    "rootless/MW17",
)

parents = []
children = []
transports = []
parent = frame2

for index, step in enumerate(steps):
    parents.append(parent)
    v = vector(ZZ, [ZZ(x) for x in step["v"].split(",")])
    q = ZZ(step["q"])
    a = ZZ(step["a"])
    b = ZZ(step["b"])

    child, T = neighbor(parent,q,a,b,v)
    assert child.det() == 948

    transports.append(T)
    children.append(child)

    print(
        "Q80MARKTRACESTEP|step={}|target={}|q={}|a={}|b={}|"
        "T_det={}|status=PASS_REPLAY_TRANSPORT".format(
            index+1,expected_children[index],q,a,b,int(T.det())
        ),
        flush=True,
    )
    parent = child

candidate1 = children[6]
assert candidate1 == parents[7]
assert len(roots_norm_two(candidate1)) == 2

# 3. Recover final-q6 horizontal P on candidate1.
final_step = steps[7]
v_final = vector(ZZ, [ZZ(x) for x in final_step["v"].split(",")])
raw_final = vector(
    ZZ,
    [ZZ(final_step["a"]),ZZ(final_step["b"])] + list(v_final)
)

Dred, reflection_sequence, a1r = chamber_reduce_a1(
    raw_final,candidate1
)

Fc = standard_F(candidate1)
Oc = standard_O(candidate1)
P_candidate1 = Dred - Oc + Fc

assert Dred * full_ns(candidate1) * Dred == 0
assert intersection(Dred,Fc,candidate1) == 2
assert intersection(Dred,Oc,candidate1) == 1
assert intersection(P_candidate1,Fc,candidate1) == 1
assert intersection(P_candidate1,P_candidate1,candidate1) == -2
assert intersection(P_candidate1,Oc,candidate1) == 4
assert Dred == Oc + P_candidate1 - Fc

print(
    "Q80MARKCAND1|raw={}|reduced={}|reflections={}|A1_root={}|"
    "P={}|P2={}|P_F={}|P_O={}|"
    "status=PASS_EXACT_CANDIDATE1_FINAL_HORIZONTAL_CLASS".format(
        tuple(vec(raw_final)),tuple(vec(Dred)),reflection_sequence,
        tuple(vec(a1r)),tuple(vec(P_candidate1)),
        int(intersection(P_candidate1,P_candidate1,candidate1)),
        int(intersection(P_candidate1,Fc,candidate1)),
        int(intersection(P_candidate1,Oc,candidate1)),
    ),
    flush=True,
)

# 4. Transport P backwards through selected child frames.
back_labels = (
    "A1/MW16 candidate1",
    "4A1/MW13 (6855 child)",
    "A3+A2/MW12 (1938 child)",
    "A4+A2+A1/MW10 (7774 child)",
    "A6+A3/MW8",
    "A6+A4/MW7",
    "D7+D4/MW6",
    "D7+D5/MW5",
    "D9+A4/MW4",
    "E6+D5+A3/MW3 start",
)

back_frames = [
    candidate1,
    parents[6],
    parents[5],
    parents[4],
    parents[3],
    parents[2],
    parents[1],
    parents[0],
    frame1,
    start,
]

snapshots = []

def record_snapshot(label,frame,P):
    F = standard_F(frame)
    O = standard_O(frame)
    row = {
        "frame":label,
        "class":vec(P),
        "square":int(intersection(P,P,frame)),
        "fiber_degree":int(intersection(P,F,frame)),
        "zero_intersection":int(intersection(P,O,frame)),
    }
    assert row["square"] == -2
    snapshots.append(row)
    print(
        "Q80MARKBACK|frame={}|class={}|square={}|Fdeg={}|Oint={}|"
        "status=PASS_TRANSPORTED_CLASS".format(
            label,tuple(row["class"]),row["square"],
            row["fiber_degree"],row["zero_intersection"]
        ),
        flush=True,
    )

P = vector(ZZ,P_candidate1)
record_snapshot(back_labels[0],back_frames[0],P)

# T7,T6,...,T1. T8 maps rootless -> candidate1 and is intentionally omitted.
for out_index, transport_index in enumerate(range(6,-1,-1), start=1):
    T = transports[transport_index]
    child_frame = children[transport_index]
    parent_frame = parents[transport_index]

    P_old = P*T
    assert (
        P * full_ns(child_frame) * P
        == P_old * full_ns(parent_frame) * P_old
    )
    P = vector(ZZ,P_old)
    record_snapshot(back_labels[out_index],back_frames[out_index],P)

# Common prefix D7+D5 -> D9+A4 -> E6+D5+A3.
for T,label,frame in (
    (T_prefix2,back_labels[8],back_frames[8]),
    (T_prefix1,back_labels[9],back_frames[9]),
):
    P = vector(ZZ,P*T)
    record_snapshot(label,frame,P)

P_start = vector(ZZ,P)

# 5. Explicit composed transports.
T_candidate1_to_start = identity_matrix(ZZ,19)
for transport_index in range(6,-1,-1):
    T_candidate1_to_start *= transports[transport_index]
T_candidate1_to_start *= T_prefix2
T_candidate1_to_start *= T_prefix1

assert abs(T_candidate1_to_start.det()) == 1
assert P_candidate1*T_candidate1_to_start == P_start
assert (
    T_candidate1_to_start
    * full_ns(start)
    * T_candidate1_to_start.transpose()
    == full_ns(candidate1)
)

# Candidate1 -> A6+A3 is T7*T6*T5*T4.
T_candidate1_to_a6a3 = identity_matrix(ZZ,19)
for transport_index in range(6,2,-1):
    T_candidate1_to_a6a3 *= transports[transport_index]

P_a6a3 = vector(ZZ,P_candidate1*T_candidate1_to_a6a3)
assert P_a6a3 == vector(ZZ,snapshots[4]["class"])

print(
    "Q80MARKCOMPOSE|candidate1_to_A6A3_det={}|candidate1_to_start_det={}|"
    "P_A6A3={}|P_start={}|"
    "status=PASS_COMPOSED_MARKED_TRANSPORT".format(
        int(T_candidate1_to_a6a3.det()),
        int(T_candidate1_to_start.det()),
        tuple(vec(P_a6a3)),tuple(vec(P_start)),
    ),
    flush=True,
)

payload = {
    "status":"PASS_Q80_CANDIDATE1_MARKED_TRANSPORT",
    "meaning":(
        "Exact final-q6 horizontal (-2)-class recovered on generic "
        "candidate1 and transported backwards through the retained greedy "
        "neighbor chain. No CM24 equation data used."
    ),
    "source_files":[
        "kumar_q80_to_rootless_path.tsv",
        "kumar_q80_new_lowq_rootless_path.tsv",
        "kumar_q80_e6_d5_a3_mw3_frame.txt",
    ],
    "candidate1":{
        "raw_final_q6_divisor":vec(raw_final),
        "reduced_final_q6_divisor":vec(Dred),
        "reflection_sequence":[list(x) for x in reflection_sequence],
        "a1_root":vec(a1r),
        "horizontal_P":vec(P_candidate1),
        "P_square":-2,
        "P_F":1,
        "P_O":4,
        "decomposition":"D_reduced = O + P - F",
    },
    "snapshots":snapshots,
    "P_A6A3":vec(P_a6a3),
    "P_start_E6_D5_A3":vec(P_start),
    "transport_candidate1_to_A6A3":mat_rows(T_candidate1_to_a6a3),
    "transport_candidate1_to_start":mat_rows(T_candidate1_to_start),
    "step_transport_determinants":[int(T.det()) for T in transports],
    "prefix_transport_determinants":[
        int(T_prefix1.det()),int(T_prefix2.det())
    ],
    "next":(
        "Express P_A6A3 or P_start in the corresponding explicit marked "
        "MW/section basis; then carry that marking forward through the exact "
        "equation-level neighbor maps and compare only at the end with the "
        "CM24 P2-P3 mod-73 marking."
    ),
}

artifact = OUT / "q80-candidate1-marked-transport.json"
artifact.write_text(json.dumps(payload,indent=2,default=int)+"\n")

print(
    "Q80MARKFINAL|artifact={}|P_candidate1={}|P_A6A3={}|P_start={}|"
    "status=PASS_Q80_CANDIDATE1_MARKED_TRANSPORT".format(
        artifact,tuple(vec(P_candidate1)),
        tuple(vec(P_a6a3)),tuple(vec(P_start))
    ),
    flush=True,
)
