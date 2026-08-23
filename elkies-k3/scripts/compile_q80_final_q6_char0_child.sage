#!/usr/bin/env sage
"""
Compile the certified exact Q80 final-q6 RR pencil.

Input:
  artifacts/generated-results/q80-final-q6-char0-rr.json
  artifacts/generated-results/q80-q4-candidate1-char0-family.json

The RR artifact already certifies:
  * exact H=P2-P3 on the selected characteristic-zero candidate1 parent;
  * whole-A4 quotient;
  * connected-A5 quotient;
  * ambient 4, rank 2, kernel dimension 2, h0(D)=2.

This script performs only the mechanical final hand-off:
  1. build the degree-two chord quartic from the exact 2D RR kernel;
  2. compute the Jacobian;
  3. finite-minimize and classify its fibres exactly;
  4. certify the CM24 characteristic-zero child
       4 I3 + I4 + I6 + 2 I1
     with smooth infinity, root rank 16 and root determinant 1944.

The generic lattice endpoint rootless/MW17 is a separate generic certificate;
this equation is the exact CM24 characteristic-zero specialization.
"""

from pathlib import Path
import json

from sage.all import (
    PolynomialRing, QuadraticField, vector, sage_eval
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
OUT = ROOT / "artifacts" / "generated-results"
DATA = ROOT / "elkies-k3" / "data" / "fibrations" / "q80-final-q6-char0"
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

RR_ART = OUT / "q80-final-q6-char0-rr.json"
PARENT_ART = OUT / "q80-q4-candidate1-char0-family.json"

if not RR_ART.exists():
    raise SystemExit(f"missing final q6 RR artifact: {RR_ART}")
if not PARENT_ART.exists():
    raise SystemExit(f"missing candidate1 parent artifact: {PARENT_ART}")

rr = json.loads(RR_ART.read_text())
parent_family = json.loads(PARENT_ART.read_text())

assert rr["status"] == "PASS_EXACT_FINAL_Q6_RR"
assert parent_family["status"] == "PASS_EXACT_Q4_CANDIDATE1_FAMILY"

selected = parent_family["selected"]
assert selected["target_candidate1"]
assert int(selected["index"]) == int(rr["source_candidate1_index"])

load(str(HERE / "elliptic_neighbor_compiler.sage"))
load(str(HERE / "elliptic_neighbor_compiler_field_generic.sage"))

K = QuadraticField(-3, "j")
j = K.gen()

# Old candidate1 base W.
RW0 = PolynomialRing(K, "W")
W0 = RW0.gen()
KW0 = RW0.fraction_field()

A = RW0(sage_eval(selected["minimal_A"], locals={"j":j, "W":W0}))
B = RW0(sage_eval(selected["minimal_B"], locals={"j":j, "W":W0}))
Hx = RW0(sage_eval(rr["horizontal"]["Hx"], locals={"j":j, "W":W0}))
Hy = RW0(sage_eval(rr["horizontal"]["Hy"], locals={"j":j, "W":W0}))

assert A.degree() == 8 and B.degree() == 12
assert Hx.degree() <= 4 and Hy.degree() <= 6
assert Hy^2 == Hx^3 + A*Hx + B
assert tuple(rr["horizontal"]["hits"]) == ("I3","I4","I6")

def parse_k(value):
    return K(sage_eval(str(value), locals={"j":j}))

def parse_row(values):
    return vector(K, [parse_k(value) for value in values])

k1 = parse_row(rr["rr"]["kernel"][0])
k2 = parse_row(rr["rr"]["kernel"][1])
assert len(k1) == 4 and len(k2) == 4

# New base Z.
RZ = PolynomialRing(K, "Z")
Z0 = RZ.gen()
KZ = RZ.fraction_field()
Z = KZ(Z0)

# Old W variable over the new-base function field.
RW = PolynomialRing(KZ, "W")
W = RW.gen()

def lift_w(poly):
    poly = RW0(poly)
    return RW([KZ(c) for c in poly.list()])

def chord_coeff(row):
    a = KZ(row[0]) + KZ(row[1])*W + KZ(row[2])*W^2
    b = KZ(row[3])
    return RW(a), b

a0,b0 = chord_coeff(k1)
a1,b1 = chord_coeff(k2)

print(
    "Q80FINALCOMPILE|phase=quartic|"
    "status=RUNNING_EXACT_FINAL_Q6_COMPILE",
    flush=True,
)

hop = compile_degree_two_chord_hop(
    RW,
    Z,
    a0,b0,a1,b1,
    lift_w(Hx),-lift_w(Hy),
    lift_w(A),lift_w(B),
)

quartic = hop["binary_quartic"]
if quartic.degree() != 4:
    raise ArithmeticError(
        f"final q6 produced quartic degree {quartic.degree()}, expected 4"
    )

Araw = KZ(hop["jacobian_a"])
Braw = KZ(hop["jacobian_b"])

print(
    "Q80FINALCOMPILE|quartic_degree=4|phase=exact_fibres|"
    "status=PASS_EXACT_FINAL_BINARY_QUARTIC",
    flush=True,
)

classification = classify_finite_short_weierstrass_fibres(
    RZ,Araw,Braw
)

finite_totals = {}
finite_dump = []
for item in classification["finite_fibres"]:
    symbol = str(item["kodaira"])
    degree = int(item["degree"])
    finite_totals[symbol] = finite_totals.get(symbol,0) + degree
    finite_dump.append({
        "factor":str(item["factor"]),
        "degree":degree,
        "kodaira":symbol,
        "raw_orders":[int(v) for v in item["raw_orders"]],
        "minimal_orders":[int(v) for v in item["minimal_orders"]],
        "scaling":int(item["scaling"]),
    })

infinity = classification["infinity_boundary"]
infinity_orders = tuple(int(v) for v in infinity["normalized_orders"])
infinity_smooth = (infinity_orders[2] == 0)

minim = classification["finite_minimization"]
Amin = RZ(minim["minimal_a"])
Bmin = RZ(minim["minimal_b"])
scale = KZ(minim["scaling_unit"])

root_rank = int(classification["finite_root_rank"])
root_det = int(classification["finite_root_determinant"])
euler = int(classification["finite_euler_number"])

expected_finite = {
    "I3":4,
    "I4":1,
    "I6":1,
    "I1":2,
}

is_target = (
    finite_totals == expected_finite
    and infinity_smooth
    and root_rank == 16
    and root_det == 1944
    and euler == 24
)

print(
    "Q80FINALCHILD|finite={}|infinity_orders={}|infinity_smooth={}|"
    "root_rank={}|root_det={}|euler={}|target={}|"
    "status=PASS_EXACT_FINAL_CHILD_CLASSIFICATION".format(
        finite_totals,infinity_orders,int(infinity_smooth),
        root_rank,root_det,euler,int(is_target)
    ),
    flush=True,
)

if not is_target:
    raise ArithmeticError(
        "exact final q6 child does not match "
        "4I3+I4+I6+2I1 / root rank16 det1944"
    )

jmap = KZ(
    1728 * 4*Araw^3 / (4*Araw^3 + 27*Braw^2)
)

payload = {
    "status":"PASS_EXACT_Q80_FINAL_Q6_CHILD",
    "field":"QQ(sqrt(-3))",
    "meaning":(
        "Exact CM24 characteristic-zero final q6 specialization compiled "
        "from the certified RR pencil. The separate generic lattice endpoint "
        "is rootless/MW17."
    ),
    "source_rr":"q80-final-q6-char0-rr.json",
    "source_candidate1_index":int(rr["source_candidate1_index"]),
    "horizontal":{
        "identity":"P2-P3 (historical marking; exact section recovered via MW basis difference)",
        "Hx":str(Hx),
        "Hy":str(Hy),
        "hits":["I3","I4","I6"],
    },
    "rr":{
        "kernel":[[str(v) for v in k1],[str(v) for v in k2]],
        "h0":2,
    },
    "quartic":str(quartic),
    "jacobian_A_raw":str(Araw),
    "jacobian_B_raw":str(Braw),
    "finite_scaling_unit":str(scale),
    "minimal_A":str(Amin),
    "minimal_B":str(Bmin),
    "finite_fibres":finite_dump,
    "finite_totals":finite_totals,
    "infinity_orders":[int(v) for v in infinity_orders],
    "infinity_smooth":bool(infinity_smooth),
    "root_rank":root_rank,
    "root_determinant":root_det,
    "finite_euler":euler,
    "specialized_root_lattice":"4A2+A3+A5",
    "specialized_MW_rank":2,
    "j_map":str(jmap),
    "generic_endpoint":{
        "root_lattice":"rootless",
        "MW_rank":17,
        "determinant":948,
        "status":"previously certified by generic lattice neighbor search",
    },
}

artifact = OUT / "q80-final-q6-char0-child.json"
artifact.write_text(json.dumps(payload,indent=2,default=int)+"\n")

model = DATA / "q80_char0_final_q6_child.sage"
model.write_text(
    "\n".join([
        "#!/usr/bin/env sage",
        "from sage.all import QuadraticField, PolynomialRing",
        'K = QuadraticField(-3, "j")',
        "j = K.gen()",
        'R = PolynomialRing(K, "Z")',
        "Z = R.gen()",
        f"A = {Amin}",
        f"B = {Bmin}",
        'print("Q80FINALMODEL|status=PASS_EXACT_Q80_FINAL_Q6_CHILD")',
    ]) + "\n"
)

note = DATA / "Q80_CHAR0_FINAL_Q6_CERTIFICATE.md"
note.write_text(
    "# Q80 final q6 — exact characteristic-zero CM24 certificate\n\n"
    "Status: **PASS_EXACT_Q80_FINAL_Q6_CHILD**\n\n"
    "- source: selected exact q4 candidate1 parent;\n"
    "- horizontal: exact `P2-P3`, recovered as a difference of easier exact MW sections;\n"
    "- resolved RR: ambient `4`, condition rank `2`, `h0(D)=2`;\n"
    "- exact quartic degree: `4`;\n"
    "- finite fibres: `4 I3 + I4 + I6 + 2 I1`;\n"
    "- infinity: smooth;\n"
    "- root lattice: `4A2+A3+A5`;\n"
    "- root rank: `16`;\n"
    "- root determinant: `1944`;\n"
    "- CM24 MW rank: `2`;\n"
    "- Euler number: `24`.\n\n"
    "This is the characteristic-zero CM24 specialization of the separately "
    "certified generic final neighbour `A1/MW16 -> rootless/MW17`.\n"
)

print(
    "Q80FINALCHILDFINAL|artifact={}|model={}|certificate={}|"
    "finite={}|root_rank={}|root_det={}|MW=2|"
    "status=PASS_EXACT_Q80_FINAL_Q6_CHILD".format(
        artifact,model,note,finite_totals,root_rank,root_det
    ),
    flush=True,
)
