#!/usr/bin/env sage
"""
Compile every exact q4-candidate1 RR survivor.

The previous probe proved a finite family of exact characteristic-zero
horizontals, and for each one a complete 4 -> 2 resolved RR pencil.
This script does not guess which representative is "-P3".  Instead it:

  * compiles the binary quartic/Jacobian for every exact survivor;
  * finite-minimizes and classifies every child exactly;
  * retains precisely the branches with the candidate1 root configuration
        I2 + I3 + I4 + I5 + I6 + 4 I1
    and smooth infinity;
  * records all surviving branches for the final q6 stage.

Thus any residual symmetry/marking ambiguity is propagated explicitly rather
than hidden behind a modular gauge that is bad at the supported D4 fibre.
"""

from pathlib import Path
import json

from sage.all import (
    QQ, PolynomialRing, QuadraticField, vector, sage_eval
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
PARENT_ART = ROOT / "artifacts" / "generated-results" / "q80-q4-6855-char0-resolved-rr.json"
FAMILY_ART = ROOT / "artifacts" / "generated-results" / "q80-q4-candidate1-char0-rr-family.json"
OUT = ROOT / "artifacts" / "generated-results"
DATA = ROOT / "elkies-k3" / "data" / "fibrations" / "q80-q4-candidate1-char0"
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

if not PARENT_ART.exists():
    raise SystemExit(f"missing q4_6855 artifact: {PARENT_ART}")
if not FAMILY_ART.exists():
    raise SystemExit(f"missing candidate1 RR family: {FAMILY_ART}")

parent_data = json.loads(PARENT_ART.read_text())
family = json.loads(FAMILY_ART.read_text())
assert parent_data["status"] == "PASS_EXACT_Q4_6855_RESOLVED_RR"
assert family["status"] == "PASS_EXACT_CANDIDATE1_RR_FAMILY"

load(str(HERE / "elliptic_neighbor_compiler.sage"))
load(str(HERE / "elliptic_neighbor_compiler_field_generic.sage"))

K = QuadraticField(-3, "j")
j = K.gen()

R = PolynomialRing(K, "U")
U0 = R.gen()
KF = R.fraction_field()
Ubase = KF(U0)

A = R(sage_eval(parent_data["child"]["minimal_A"], locals={"j":j, "U":U0}))
B = R(sage_eval(parent_data["child"]["minimal_B"], locals={"j":j, "U":U0}))
assert A.degree() == 8 and B.degree() == 12

WR = PolynomialRing(K, "W")
W0 = WR.gen()
KW = WR.fraction_field()
W = KW(W0)

RU = PolynomialRing(KW, "U")
U = RU.gen()

def lift_u(poly):
    poly = R(poly)
    return RU([KW(c) for c in poly.list()])

A_l = lift_u(A)
B_l = lift_u(B)

def parse_k(value):
    return K(sage_eval(str(value), locals={"j":j}))

def parse_poly(value):
    return R(sage_eval(str(value), locals={"j":j, "U":U0}))

def parse_row(row):
    return vector(K, [parse_k(v) for v in row])

def chord_coeff(row):
    aa = KW(row[0]) + KW(row[1])*U + KW(row[2])*U^2
    bb = KW(row[3])
    return RU(aa), bb

expected_finite = {
    "I2":1,
    "I3":1,
    "I4":1,
    "I5":1,
    "I6":1,
    "I1":4,
}

results = []
target_branches = []

for record in family["candidates"]:
    index = int(record["index"])
    Hx = parse_poly(record["Hx"])
    Hy = parse_poly(record["Hy"])
    k1 = parse_row(record["kernel"][0])
    k2 = parse_row(record["kernel"][1])

    assert Hy^2 == Hx^3 + A*Hx + B
    assert len(k1) == 4 and len(k2) == 4

    a0,b0 = chord_coeff(k1)
    a1,b1 = chord_coeff(k2)

    print(
        f"Q80CAND1COMPILE|index={index}|phase=quartic|status=RUNNING",
        flush=True,
    )

    hop = compile_degree_two_chord_hop(
        RU,
        W,
        a0,b0,a1,b1,
        lift_u(Hx),-lift_u(Hy),
        A_l,B_l,
    )
    quartic = hop["binary_quartic"]
    if quartic.degree() != 4:
        raise ArithmeticError(
            f"candidate {index} produced quartic degree {quartic.degree()}, expected 4"
        )

    Araw = KW(hop["jacobian_a"])
    Braw = KW(hop["jacobian_b"])

    print(
        f"Q80CAND1COMPILE|index={index}|quartic_degree=4|"
        "phase=exact_fibres|status=RUNNING",
        flush=True,
    )

    classification = classify_finite_short_weierstrass_fibres(WR,Araw,Braw)

    finite_totals = {}
    finite_dump = []
    for item in classification["finite_fibres"]:
        symbol = str(item["kodaira"])
        degree = int(item["degree"])
        finite_totals[symbol] = finite_totals.get(symbol,0)+degree
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
    Amin = WR(minim["minimal_a"])
    Bmin = WR(minim["minimal_b"])
    scale = KW(minim["scaling_unit"])

    # Intrinsic root-lattice invariants of the target configuration.
    root_rank = int(classification["finite_root_rank"])
    root_det = int(classification["finite_root_determinant"])
    euler = int(classification["finite_euler_number"])

    is_target = (
        finite_totals == expected_finite
        and infinity_smooth
        and root_rank == 15
        and root_det == 720
        and euler == 24
    )

    # j-map is useful to see whether surviving branches are literally the same
    # W-model or only related by a base transformation.
    jmap = KW(
        1728 * 4*Araw^3 / (4*Araw^3 + 27*Braw^2)
    )

    result = {
        "index":index,
        "horizontal":{
            "Hx":str(Hx),
            "Hy":str(Hy),
        },
        "rr":{
            "kernel":[[str(x) for x in k1],[str(x) for x in k2]],
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
        "root_rank":root_rank,
        "root_determinant":root_det,
        "finite_euler":euler,
        "j_map":str(jmap),
        "target_candidate1":bool(is_target),
    }
    results.append(result)

    print(
        "Q80CAND1CHILD|index={}|finite={}|infinity_smooth={}|"
        "root_rank={}|root_det={}|euler={}|target={}|"
        "status=PASS_EXACT_CHILD_CLASSIFICATION".format(
            index,finite_totals,int(infinity_smooth),
            root_rank,root_det,euler,int(is_target)
        ),
        flush=True,
    )

    if is_target:
        target_branches.append(result)

if not target_branches:
    raise ArithmeticError(
        "none of the exact candidate1 RR survivors produced the target "
        "A1+A2+A3+A4+A5 child"
    )

# Do not manufacture uniqueness.  Preserve every exact target branch.
# For convenience only, choose the smallest-index branch as a deterministic
# continuation representative; the full target list remains in the artifact.
selected = min(target_branches, key=lambda item: item["index"])

# Compare target j-maps literally in the chosen W coordinate.
same_j_as_selected = [
    item["index"]
    for item in target_branches
    if item["j_map"] == selected["j_map"]
]

payload = {
    "status":"PASS_EXACT_Q4_CANDIDATE1_FAMILY",
    "field":"QQ(sqrt(-3))",
    "parent":"q80-q4-6855-char0-resolved-rr.json",
    "rr_family":"q80-q4-candidate1-char0-rr-family.json",
    "candidate_count":len(results),
    "target_count":len(target_branches),
    "target_indices":[int(item["index"]) for item in target_branches],
    "selected_index":int(selected["index"]),
    "selected_same_j_indices":[int(v) for v in same_j_as_selected],
    "expected_child":{
        "root_lattice":"A1+A2+A3+A4+A5",
        "root_rank":15,
        "root_determinant":720,
        "finite_pattern":"I2 + I3 + I4 + I5 + I6 + 4 I1",
        "infinity":"smooth",
    },
    "candidates":results,
    "selected":selected,
    "next":"final q6 exact characteristic-zero branch propagation",
}

json_path = OUT / "q80-q4-candidate1-char0-family.json"
json_path.write_text(json.dumps(payload,indent=2,default=int)+"\n")

# Persist the selected exact minimal parent in a directly loadable Sage file.
parent_file = DATA / "q80_char0_candidate1_selected_parent.sage"
parent_file.write_text(
    "\n".join([
        "#!/usr/bin/env sage",
        "from sage.all import QuadraticField, PolynomialRing",
        'K = QuadraticField(-3, "j")',
        "j = K.gen()",
        'R = PolynomialRing(K, "W")',
        "W = R.gen()",
        f"A = {selected['minimal_A']}",
        f"B = {selected['minimal_B']}",
        f"source_candidate_index = {selected['index']}",
        'print("Q80CAND1PARENT|index={}|status=PASS_EXACT_SELECTED_PARENT".format(source_candidate_index))',
    ]) + "\n"
)

note = DATA / "Q80_CHAR0_Q4_CANDIDATE1_FAMILY_CERTIFICATE.md"
note.write_text(
    "# Q80 q4 candidate1 — characteristic-zero exact family certificate\n\n"
    "Status: **PASS_EXACT_Q4_CANDIDATE1_FAMILY**\n\n"
    f"- exact RR survivors compiled: `{len(results)}`;\n"
    f"- exact target children: `{len(target_branches)}` at indices "
    f"`{[item['index'] for item in target_branches]}`;\n"
    "- target fibre pattern: `I2 + I3 + I4 + I5 + I6 + 4 I1`;\n"
    "- root lattice: `A1+A2+A3+A4+A5`, rank `15`, determinant `720`;\n"
    "- infinity smooth;\n"
    f"- deterministic continuation representative: index `{selected['index']}`;\n"
    f"- target branches with the identical displayed `j(W)` map as the "
    f"selected branch: `{same_j_as_selected}`.\n\n"
    "All exact target branches are retained in the JSON artifact; selection "
    "of the smallest index is a continuation convention, not a uniqueness "
    "claim.\n"
)

print(
    "Q80CAND1FAMILYFINAL|candidate_count={}|target_count={}|"
    "target_indices={}|selected={}|same_j={}|json={}|parent={}|"
    "status=PASS_EXACT_Q4_CANDIDATE1_FAMILY".format(
        len(results),len(target_branches),
        tuple(item["index"] for item in target_branches),
        selected["index"],tuple(same_j_as_selected),
        json_path,parent_file,
    ),
    flush=True,
)
