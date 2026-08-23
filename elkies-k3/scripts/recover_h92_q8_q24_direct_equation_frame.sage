#!/usr/bin/env sage -python
"""
Recover the H92 q24 divisor directly in the CURRENT equation-level D13 frame.

We deliberately avoid transporting 19-coordinate vectors between the q24
producer and the equation certifier, because the live marking audit proved
those ambient coordinate realizations differ.

Instead:
  * obtain marking-INVARIANT q24 data from the native q24 producer at zero
    II*_E8_1;
  * compute the exact equation-frame AJ(S3) MW class;
  * require matching height / D13 correction / P.O;
  * enumerate integral lifts of that SAME MW class in the equation D13 frame,
    ordered by positive-frame norm;
  * use the native q24 vertical fibre coefficient to pin the required lift norm;
  * match vertical scalar invariants;
  * require the resulting degree-2 primitive isotropic divisor to have D12
    orthogonal child.

A unique hit closes the q24 divisor without any unsafe ambient-vector reuse.
"""

import json
import sys
from pathlib import Path

from sage.all import (
    IntegralLattice, QQ, ZZ, block_diagonal_matrix, gcd, identity_matrix,
    lcm, matrix, pari, vector
)


def locate_repo():
    cwd = Path.cwd().resolve()
    candidates = [cwd, *cwd.parents]
    h = Path.home()
    candidates += [
        h / "Documents" / "jacobian-research",
        h / "jacobian-research",
        h / "src" / "jacobian-research",
        h / "git" / "jacobian-research",
        h / "projects" / "jacobian-research",
    ]
    seen = set()
    for c in candidates:
        try:
            c = c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (
            (c / "elkies-k3/scripts").is_dir()
            and (c / "artifacts/generated-results").is_dir()
        ):
            return c
    raise SystemExit("Could not locate jacobian-research")


ROOT = locate_repo()
SCRIPTS = ROOT / "elkies-k3/scripts"
LOCAL = ROOT / "artifacts/local/elkies-k3"
EQ = SCRIPTS / "certify_h92_q8_equation_ns_divisor.sage"
Q24 = SCRIPTS / "audit_h92_q8_q24_effective_zero_choices.sage"
OUT = LOCAL / "q8-q24-equation-frame-direct-recovery.json"
TMP = LOCAL / "q8-q24-native-invariants-temp.json"

for path in (EQ, Q24):
    if not path.exists():
        raise SystemExit(f"missing {path}")


def run_script(path, argv):
    saved = list(sys.argv)
    scope = {"__name__": "__embedded__"}
    try:
        sys.argv = [str(path)] + list(argv)
        exec(compile(path.read_text(), str(path), "exec"), scope)
    finally:
        sys.argv = saved
    return scope


print("Q8Q24DIRECT|stage=equation_source", flush=True)
eq = run_script(EQ, [])

print("Q8Q24DIRECT|stage=native_q24_invariants", flush=True)
qp = run_script(Q24, ["--output", str(TMP)])

# ---------------------------------------------------------------------------
# Native q24 invariants at the SAME geometric zero II*_E8_1.
# ---------------------------------------------------------------------------

profiles = qp["profiles"]
native = next((p for p in profiles if p["zero"] == "II*_E8_1"), None)
if native is None:
    raise SystemExit("native q24 producer has no II*_E8_1 zero profile")

nsn = qp["source_ns"]
Dn = vector(ZZ, qp["D_actual"])
On = vector(ZZ, native["zero_source_h3_ns"])
Pn = vector(ZZ, native["section_source_h3_ns"])
Vn = Dn - On - Pn

assert Dn * nsn * Dn == 0
assert Dn * nsn * qp["F_actual"] == 2
assert Pn * nsn * Pn == -2
assert Pn * nsn * qp["F_actual"] == 1
assert Vn * nsn * qp["F_actual"] == 0

native_inv = {
    "height": str(native["height"]),
    "correction": str(native["D13_correction"]),
    "pole": ZZ(native["P_dot_O"]),
    "vertical_F": ZZ(native["vertical_fibre_coefficient"]),
    "V2": ZZ(Vn * nsn * Vn),
    "VO": ZZ(Vn * nsn * On),
    "VP": ZZ(Vn * nsn * Pn),
}

print(
    "Q8Q24DIRECT_NATIVE|"
    f"height={native_inv['height']}|corr={native_inv['correction']}|"
    f"PdotO={native_inv['pole']}|vertical_F={native_inv['vertical_F']}|"
    f"V2={native_inv['V2']}|VO={native_inv['VO']}|VP={native_inv['VP']}|"
    "status=PASS_NATIVE_INVARIANTS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Equation-frame D13 + AJ(S3).
# ---------------------------------------------------------------------------

ns = eq["ns"]
F8 = vector(ZZ, eq["F8eq"])
S3 = vector(ZZ, eq["S3std"])
selected = eq["selected"]
isotropic_mate = eq["isotropic_mate"]
roots_and_data = eq["roots_and_data"]

# Reuse the q24 producer's deterministic lattice utilities only on the
# equation-frame Gram matrix.  No ambient vectors cross the marking boundary.
child_frame_with_zero = qp["child_frame_with_zero"]
d13_root_adaptation = qp["d13_root_adaptation"]
section_for_mw = qp["section_for_mw"]

O = vector(ZZ, selected["E8"]["simple_root_vectors_in_source_h3_ns"][0])
assert O * ns * O == -2
assert O * ns * F8 == 1

child8, Bzero = child_frame_with_zero(ns, F8, O)
assert roots_and_data(child8)[2] == (13, 312, 4)
A13, adapted, H = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ, 2), A13) * Bzero
Gadapt = block_diagonal_matrix(matrix(ZZ, ((0,1),(1,0))), -adapted)
assert Badapt * ns * Badapt.transpose() == Gadapt


def coords(C):
    c = vector(QQ, C) * Badapt.inverse()
    assert all(x in ZZ for x in c)
    return vector(ZZ, c)


cS3 = coords(S3)
assert cS3[1] == 52
z = vector(ZZ, cS3[-4:])

Pcoords, h, class_ord, corr, pole = section_for_mw(adapted, z, 32768)
P = vector(ZZ, Pcoords) * Badapt

assert P * ns * P == -2
assert P * ns * F8 == 1
assert P * ns * O == pole

print(
    "Q8Q24DIRECT_AJ|"
    f"mw={','.join(map(str,z))}|height={h}|class_order={class_ord}|"
    f"corr={corr}|PdotO={pole}|status=PASS_EQUATION_AJ",
    flush=True,
)

horizontal_invariants_match = (
    str(h) == native_inv["height"]
    and str(corr) == native_inv["correction"]
    and ZZ(pole) == native_inv["pole"]
)

print(
    "Q8Q24DIRECT_HORIZONTAL|"
    f"invariants_match={int(horizontal_invariants_match)}|"
    f"native={native_inv['height']},{native_inv['correction']},{native_inv['pole']}|"
    f"equation={h},{corr},{pole}|"
    f"status={'PASS' if horizontal_invariants_match else 'MISMATCH'}",
    flush=True,
)

if not horizontal_invariants_match:
    print(
        "Q8Q24DIRECT_RESULT|status=AJ_S3_NOT_NATIVE_Q24_HORIZONTAL",
        flush=True,
    )
    raise SystemExit(2)

# ---------------------------------------------------------------------------
# Exact lift shell forced by the native vertical fibre coefficient.
#
# In adapted U+positive-frame coordinates:
#   O = [-1,1,0]
#   P = [aP,1,wP]
#   D = O + P + f F + Rroot
# so D's first U coefficient is aD = -1+aP+f.
# For degree two D=[aD,2,wD], isotropicity says ||wD|| = 4 aD.
# ---------------------------------------------------------------------------

aP = ZZ(Pcoords[0])
f_native = native_inv["vertical_F"]
target_norm = ZZ(4 * (-1 + aP + f_native))
minimum_norm = ZZ(vector(ZZ, Pcoords[2:]) * adapted * vector(ZZ, Pcoords[2:]))

print(
    "Q8Q24DIRECT_SHELL|"
    f"section_a={aP}|vertical_F={f_native}|"
    f"minimum_lift_norm={minimum_norm}|target_lift_norm={target_norm}|"
    f"extra_norm={target_norm-minimum_norm}|status=PASS",
    flush=True,
)

if target_norm < minimum_norm:
    print(
        "Q8Q24DIRECT_RESULT|status=NATIVE_FIBRE_COEFF_INCOMPATIBLE_WITH_AJ_CLASS",
        flush=True,
    )
    raise SystemExit(3)

root = adapted[:13, :13]
base = vector(ZZ, [0]*13 + list(z))
pairing = vector(QQ, base * adapted[:, :13])
dual = pairing * root.inverse()

L = IntegralLattice(root)
iterator = L.enumerate_close_vectors(-dual)

# Root basis in the equation ambient marking.
root_source = []
for i in range(13):
    e = vector(ZZ, [0,0] + [ZZ(j == i) for j in range(17)])
    r = e * Badapt
    assert r * ns * F8 == 0
    assert r * ns * O == 0
    root_source.append(r)
root_source = matrix(ZZ, [list(r) for r in root_source])

hits = []
seen_lifts = set()
examined = 0
target_shell = 0
invariant_matches = 0
classified = 0

for unused in range(250000):
    shift = vector(ZZ, next(iterator))
    lift = base + vector(ZZ, list(shift) + [0]*4)
    norm = ZZ(lift * adapted * lift)
    examined += 1

    if norm > target_norm:
        break
    if norm != target_norm:
        continue

    key = tuple(lift)
    if key in seen_lifts:
        continue
    seen_lifts.add(key)
    target_shell += 1

    # D = [aD,2,lift], with aD=norm/4.
    if norm % 4:
        continue
    aD = ZZ(norm // 4)
    Dcoords = vector(ZZ, [aD, 2] + list(lift))
    D = Dcoords * Badapt

    if D * ns * D != 0 or D * ns * F8 != 2:
        continue
    if gcd(tuple(D)) != 1:
        continue

    V = D - O - P
    inv = (
        ZZ(V * ns * V),
        ZZ(V * ns * O),
        ZZ(V * ns * P),
    )
    wanted = (
        native_inv["V2"],
        native_inv["VO"],
        native_inv["VP"],
    )
    if inv != wanted:
        continue

    invariant_matches += 1

    # Express the vertical correction exactly as f F + root.
    vertical_basis = matrix(
        QQ, [list(F8)] + [list(r) for r in root_source.rows()]
    )
    coeff = vertical_basis.solve_left(vector(QQ, V))
    assert all(x in ZZ for x in coeff)
    coeff = vector(ZZ, coeff)

    # The native fibre coefficient is itself part of the invariant package.
    if coeff[0] != f_native:
        continue

    mate = isotropic_mate(ns, D)
    orth = matrix(
        ZZ, [list(D*ns), list(mate*ns)]
    ).right_kernel_matrix()
    child = -(orth * ns * orth.transpose())
    rd = roots_and_data(child)[2]
    classified += 1

    print(
        "Q8Q24DIRECT_CANDIDATE|"
        f"index={classified}|lift_norm={norm}|"
        f"vertical_F={coeff[0]}|"
        f"vertical_roots={','.join(map(str,coeff[1:]))}|"
        f"root_data={rd[0]},{rd[1]},{rd[2]}|"
        f"status={'D12' if rd==(12,264,4) else 'OTHER'}",
        flush=True,
    )

    if rd == (12,264,4):
        hits.append({
            "D": D,
            "Dcoords": Dcoords,
            "lift": lift,
            "vertical": V,
            "vertical_coeffs": coeff,
            "root_data": rd,
        })

print(
    "Q8Q24DIRECT_SEARCH|"
    f"examined={examined}|target_shell={target_shell}|"
    f"invariant_matches={invariant_matches}|classified={classified}|"
    f"D12_hits={len(hits)}|status=PASS",
    flush=True,
)

if len(hits) != 1:
    status = "NO_D12_HIT" if not hits else "MULTIPLE_D12_HITS"
    payload = {
        "schema": "elkies-k3.h92-q8-q24-equation-frame-direct-recovery.v1",
        "status": status,
        "native_invariants": {
            k: int(v) if k in ("pole","vertical_F","V2","VO","VP") else v
            for k,v in native_inv.items()
        },
        "equation_AJ": {
            "mw": list(map(int,z)),
            "height": str(h),
            "correction": str(corr),
            "P_dot_O": int(pole),
        },
        "search": {
            "minimum_lift_norm": int(minimum_norm),
            "target_lift_norm": int(target_norm),
            "examined": examined,
            "target_shell": target_shell,
            "invariant_matches": invariant_matches,
            "classified": classified,
            "D12_hits": len(hits),
        },
    }
    OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(f"OUTPUT|{OUT}",flush=True)
    print(f"Q8Q24DIRECT_RESULT|status={status}",flush=True)
    raise SystemExit(4)

hit = hits[0]
D = hit["D"]
coeff = hit["vertical_coeffs"]

# Independent direct-profile regression.
DIRECT = LOCAL / "q8-s3-direct-x-mod-100003.json"
direct_match = None
if DIRECT.exists():
    dm = json.loads(DIRECT.read_text())
    direct_match = (
        ZZ(dm["weierstrass_structure"]["denominator_root_degree"]) == pole
        and (
            ZZ(dm["x"]["numerator_degree"]),
            ZZ(dm["x"]["denominator_degree"]),
        ) == (2*pole+4, 2*pole)
        and (
            ZZ(dm["weierstrass_structure"]["y_abs_numerator_degree"]),
            ZZ(dm["weierstrass_structure"]["y_denominator_degree"]),
        ) == (3*pole+6, 3*pole)
    )
    assert direct_match

exact_bound = False
EXACT = LOCAL / "q8-s3-direct-section-qq.json"
if EXACT.exists():
    exact = json.loads(EXACT.read_text())
    if exact.get("status") == "PASS_EXACT_Q8_S3_DIRECT_SECTION":
        assert exact["verification"]["exact_weierstrass_identity"] is True
        assert ZZ(exact["profile"]["P_dot_O_from_denominator"]) == pole
        exact_bound = True

payload = {
    "schema": "elkies-k3.h92-q8-q24-equation-frame-direct-recovery.v1",
    "status": "PASS_EXACT_Q24_EQUATION_FRAME_DIRECT_RECOVERY",
    "native_q24_invariants": {
        "zero": "II*_E8_1",
        "height": native_inv["height"],
        "D13_correction": native_inv["correction"],
        "P_dot_O": int(native_inv["pole"]),
        "vertical_fibre_coefficient": int(native_inv["vertical_F"]),
        "vertical_square": int(native_inv["V2"]),
        "vertical_dot_O": int(native_inv["VO"]),
        "vertical_dot_P": int(native_inv["VP"]),
    },
    "equation_AJ_S3": {
        "mw_coordinates": list(map(int,z)),
        "height": str(h),
        "D13_class_order": int(class_ord),
        "D13_correction": str(corr),
        "P_dot_O": int(pole),
        "section_source_h3_ns": list(map(int,P)),
        "direct_modular_profile_match": direct_match,
        "exact_characteristic_zero_section_bound": exact_bound,
    },
    "q24_equation_frame": {
        "divisor_source_h3_ns": list(map(int,D)),
        "adapted_D13_coordinates": list(map(int,hit["Dcoords"])),
        "lift_norm": int(target_norm),
        "vertical_fibre_coefficient": int(coeff[0]),
        "vertical_root_coefficients": list(map(int,coeff[1:])),
        "child_root_data": list(map(int,hit["root_data"])),
        "child_root_lattice": "D12",
        "MW_rank_if_rho19": 5,
    },
    "uniqueness": {
        "examined_close_vectors": examined,
        "target_shell_size": target_shell,
        "vertical_invariant_matches": invariant_matches,
        "classified_candidates": classified,
        "D12_hits": 1,
    },
    "boundary": (
        "This reconstructs the q24 divisor directly in the equation D13 frame "
        "without transporting ambient vectors across the mismatched historical "
        "source-H3 realizations. The remaining gap is equation-level "
        "algebraization of this certified D12-producing pencil."
    ),
}

OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q8Q24DIRECT_RESULT|"
    f"AJ_height={h}|AJ_corr={corr}|AJ_PdotO={pole}|"
    f"lift_norm={target_norm}|vertical_F={coeff[0]}|"
    f"vertical_support={sum(bool(x) for x in coeff[1:])}|"
    "child=D12|MW=5|unique=1|"
    f"char0_AJ_bound={int(exact_bound)}|"
    "status=PASS_EXACT_Q24_EQUATION_FRAME_DIRECT_RECOVERY",
    flush=True,
)
