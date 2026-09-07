#!/usr/bin/env sage -python
"""
Backtrack H92 q24 using only the canonical q8 D13 elliptic curve.

Trusted objects from the successful path:
  * direct AJ(S3) modular section A over GF(100003), including X/Z^2 and ±Y/Z^3;
  * exact canonical-D13 points G1 and G3 recovered from the IV* branch and
    old_E7_7 bisection;
  * native branch-zero MW lattice where q24=(2,-1,-1,1).

The exhaustive structural-isometry audit leaves two positive-orientation
possibilities for A in the native basis:
    A1=(0,-1, 1,1), giving R=G1-G3;
    A2=(0, 1,-1,1), giving R=G1-G2.

They are distinguished directly on the elliptic curve by the pole profile of
A+G3.  The predicted unordered pole pairs for the two global signs of A are:
    A1: {23,31}
    A2: {25,29}.

If A1 is selected, R=G1-G3 is already explicit over QQ(U), and
    q24 = A + 2R
is constructed modulo p immediately.  If the exact Hensel AJ artifact exists,
the exact QQ(U) q24 horizontal point is constructed as well.
"""

import json
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, sage_eval, vector
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
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"

AJMOD = LOCAL / "q8-s3-direct-x-mod-100003.json"
G3ART = LOCAL / "q8-d13-g3-from-e77-bisection.json"
OLDCURVES = LOCAL / "q8-explicit-old-curves.json"
EXACT_AJ = LOCAL / "q8-s3-direct-section-qq.json"
OUTPUT = LOCAL / "q8-q24-canonical-backtrack.json"
EXACT_Q24 = LOCAL / "q8-q24-horizontal-section-qq.json"

child_candidates = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
CHILD = next(
    (
        p for p in child_candidates
        if p.exists()
        and "minimal_A_coefficients_low_to_high"
            in json.loads(p.read_text()).get("child", {})
    ),
    None,
)

for p in (AJMOD, G3ART, OLDCURVES):
    if not p.exists():
        raise SystemExit(f"Missing prerequisite: {p}")
if CHILD is None:
    raise SystemExit("No complete q8 corrected D13 child artifact found")

aj = json.loads(AJMOD.read_text())
g3art = json.loads(G3ART.read_text())
old = json.loads(OLDCURVES.read_text())
child = json.loads(CHILD.read_text())

assert aj["status"] == "PASS_DIRECT_ANCHORED_Q8_S3_PROFILE_DISCOVERY"
assert g3art["status"] == "PASS_EXACT_D13_G3_FROM_E77_BISECTION"
assert old["status"] == "PASS_EXACT_Q8_EXPLICIT_OLD_CURVE_PROFILE"
assert child["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"

# ---------------------------------------------------------------------------
# Canonical D13 over QQ(U) and GF(p)(U).
# ---------------------------------------------------------------------------

RQ = PolynomialRing(QQ, "U")
UQ = RQ.gen()
KQ = RQ.fraction_field()

Aq = RQ([QQ(v) for v in child["child"]["minimal_A_coefficients_low_to_high"]])
Bq = RQ([QQ(v) for v in child["child"]["minimal_B_coefficients_low_to_high"]])
EQ = EllipticCurve(KQ, [0, 0, 0, KQ(Aq), KQ(Bq)])

p = ZZ(aj["prime"])
assert p == 100003
F = GF(p)
RF = PolynomialRing(F, "U")
UF = RF.gen()
KF = RF.fraction_field()

def modq(q):
    q = QQ(q)
    d = ZZ(q.denominator())
    if d % p == 0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {q}")
    return F(ZZ(q.numerator())) / F(d)

Af = RF([modq(v) for v in Aq.list()])
Bf = RF([modq(v) for v in Bq.list()])
EF = EllipticCurve(KF, [0, 0, 0, KF(Af), KF(Bf)])


def parse_exact_rf(text):
    return KQ(sage_eval(str(text), locals={"U": UQ}))


def reduce_poly(poly):
    poly = RQ(poly)
    return RF([modq(v) for v in poly.list()])


def reduce_rf(value):
    value = KQ(value)
    return KF(reduce_poly(value.numerator())) / KF(reduce_poly(value.denominator()))


def parse_point(record):
    assert not record.get("zero", False)
    x = parse_exact_rf(record["x"])
    y = parse_exact_rf(record["y"])
    P = EQ(x, y)
    assert P in EQ
    Pf = EF(reduce_rf(x), reduce_rf(y))
    return P, Pf


G1Q, G1F = parse_point(g3art["canonical_D13"]["G1"])
G3Q, G3F = parse_point(g3art["canonical_D13"]["G3"])
AJ77Q, AJ77F = parse_point(g3art["canonical_D13"]["AJ_old_E7_7"])
assert G1Q + G3Q == AJ77Q
assert G1F + G3F == AJ77F

# Direct AJ modular rational section.
Z = RF([F(v) for v in aj["weierstrass_structure"]["denominator_root_coefficients_low_to_high"]])
X = RF([F(v) for v in aj["x"]["numerator_coefficients_low_to_high"]])
Yabs = RF([F(v) for v in aj["weierstrass_structure"]["y_abs_numerator_coefficients_low_to_high"]])
assert Z.degree() == 24 and X.degree() == 52 and Yabs.degree() == 78

xA = KF(X) / KF(Z**2)
yAabs = KF(Yabs) / KF(Z**3)
Aplus = EF(xA, yAabs)
Aminus = -Aplus
assert Aplus in EF and Aminus in EF

# ---------------------------------------------------------------------------
# Native branch-zero height lattice: derive the two predicted A+G3 profiles.
# ---------------------------------------------------------------------------

HN = matrix(QQ, [[QQ(v) for v in row] for row in old["anchor"]["height_gram"]])
assert HN.det() == 237

G1mw = vector(ZZ, (1, 0, 0, 0))
G3mw = vector(ZZ, (0, 0, 1, 0))
q24mw = vector(ZZ, (2, -1, -1, 1))
assert QQ(G1mw * HN * G1mw) == QQ(3) / 4
assert QQ(G3mw * HN * G3mw) == QQ(11) / 4
assert QQ(q24mw * HN * q24mw) == 52

A1 = vector(ZZ, (0, -1, 1, 1))
A2 = vector(ZZ, (0, 1, -1, 1))
assert QQ(A1 * HN * A1) == 52
assert QQ(A2 * HN * A2) == 52

def unique_pole_for_height(h):
    hits = []
    for corr in (QQ(0), QQ(1), QQ(13)/4):
        pole = (QQ(h) + corr - 4) / 2
        if pole in ZZ and pole >= 0:
            hits.append((ZZ(pole), corr))
    assert len(hits) == 1, (h, hits)
    return hits[0]

def predicted_add_profiles(A):
    vals = {}
    for label, z in (
        ("A_plus_G3", A + G3mw),
        ("A_minus_G3", A - G3mw),
    ):
        h = QQ(z * HN * z)
        pole, corr = unique_pole_for_height(h)
        vals[label] = {
            "mw": list(map(int, z)),
            "height": str(h),
            "pole": int(pole),
            "correction": str(corr),
            "x_degrees": [int(2*pole + 4), int(2*pole)],
            "y_degrees": [int(3*pole + 6), int(3*pole)],
        }
    return vals

pred1 = predicted_add_profiles(A1)
pred2 = predicted_add_profiles(A2)
pair1 = sorted([pred1["A_plus_G3"]["pole"], pred1["A_minus_G3"]["pole"]])
pair2 = sorted([pred2["A_plus_G3"]["pole"], pred2["A_minus_G3"]["pole"]])
assert pair1 == [23, 31]
assert pair2 == [25, 29]

print(
    "Q8Q24BACKTRACK_PREDICT|"
    f"A1_pair={pair1[0]},{pair1[1]}|A2_pair={pair2[0]},{pair2[1]}|"
    "A1_R=G1-G3|A2_R=G1-G2|status=PASS",
    flush=True,
)


def rf_degrees(v):
    v = KF(v)
    return (
        int(RF(v.numerator()).degree()),
        int(RF(v.denominator()).degree()),
    )


def point_profile(P):
    if P.is_zero():
        return {"zero": True}
    x, y = P.xy()
    xd = rf_degrees(x)
    yd = rf_degrees(y)
    if xd[1] % 2 or yd[1] % 3:
        raise ArithmeticError(f"non-Weierstrass denominator profile x={xd} y={yd}")
    px = xd[1] // 2
    py = yd[1] // 3
    if px != py:
        raise ArithmeticError(f"x/y pole disagreement {px} vs {py}")
    return {
        "zero": False,
        "pole": px,
        "x_degrees": list(xd),
        "y_degrees": list(yd),
    }

# Both possible global orientations of direct A, added to the FIXED oriented G3.
tests = []
for sign, P in ((+1, Aplus), (-1, Aminus)):
    S = P + G3F
    prof = point_profile(S)
    tests.append((sign, P, S, prof))
    print(
        "Q8Q24BACKTRACK_G3|"
        f"A_sign={sign:+d}|PdotO={prof['pole']}|"
        f"x={prof['x_degrees'][0]}/{prof['x_degrees'][1]}|"
        f"y={prof['y_degrees'][0]}/{prof['y_degrees'][1]}|status=PASS",
        flush=True,
    )

actual_pair = sorted([rec[3]["pole"] for rec in tests])
if actual_pair == pair1:
    branch = "G1-G3"
    native_A_up_to_sign = "0,-1,1,1"
elif actual_pair == pair2:
    branch = "G1-G2"
    native_A_up_to_sign = "0,1,-1,1"
else:
    raise ArithmeticError(
        f"A± + G3 pole pair {actual_pair} matches neither predicted pair "
        f"{pair1} nor {pair2}"
    )

print(
    "Q8Q24BACKTRACK_BRANCH|"
    f"observed_pair={actual_pair[0]},{actual_pair[1]}|"
    f"selected={branch}|native_A_up_to_sign={native_A_up_to_sign}|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# If selected branch is G1-G3, orient A and construct q24 directly.
# ---------------------------------------------------------------------------

payload = {
    "schema": "elkies-k3.h92-q8-q24-canonical-backtrack.v1",
    "status": "PASS_CANONICAL_BACKTRACK_ORIENTATION",
    "prime": int(p),
    "predicted_pole_pairs": {
        "A_native_0_-1_1_1__R_G1_minus_G3": pair1,
        "A_native_0_1_-1_1__R_G1_minus_G2": pair2,
    },
    "observed_A_sign_plus_G3": [
        {
            "A_sign": int(sign),
            **prof,
        }
        for sign, unusedP, unusedS, prof in tests
    ],
    "selected_branch": branch,
    "native_A_up_to_sign": native_A_up_to_sign,
}

if branch == "G1-G3":
    # The orientation with A+G3 pole 31 is A1=(0,-1,1,1).
    chosen = next(rec for rec in tests if rec[3]["pole"] == pred1["A_plus_G3"]["pole"])
    asign, AF, unused, unused_prof = chosen

    RQpoint = G1Q - G3Q
    RFpoint = G1F - G3F
    assert reduce_rf(RQpoint[0]) == RFpoint[0]
    assert reduce_rf(RQpoint[1]) == RFpoint[1]

    rprof = point_profile(RFpoint)
    # h=4, corr=0 => P.O=0, so x,y are polynomial degrees 4,6.
    assert rprof["pole"] == 0
    assert rprof["x_degrees"][1] == 0
    assert rprof["y_degrees"][1] == 0
    assert rprof["x_degrees"][0] <= 4
    assert rprof["y_degrees"][0] <= 6

    QF = AF + 2*RFpoint
    qprof = point_profile(QF)

    # Native q24 has h=52, identity D13 component => P.O=24.
    assert qprof["pole"] == 24
    assert qprof["x_degrees"] == [52, 48]
    assert qprof["y_degrees"] == [78, 72]

    print(
        "Q8Q24BACKTRACK_R|"
        f"R=G1-G3|A_sign={asign:+d}|"
        f"R_x={rprof['x_degrees'][0]}/{rprof['x_degrees'][1]}|"
        f"R_y={rprof['y_degrees'][0]}/{rprof['y_degrees'][1]}|"
        "height=4|corr=0|PdotO=0|status=PASS_EXPLICIT_R",
        flush=True,
    )
    print(
        "Q8Q24BACKTRACK_Q24_MODP|"
        f"x={qprof['x_degrees'][0]}/{qprof['x_degrees'][1]}|"
        f"y={qprof['y_degrees'][0]}/{qprof['y_degrees'][1]}|PdotO={qprof['pole']}|"
        "formula=AJ(S3)+2*(G1-G3)|status=PASS_EXPLICIT_Q24_MODP",
        flush=True,
    )

    qx, qy = QF.xy()
    rx, ry = RQpoint.xy()
    payload.update({
        "status": "PASS_EXPLICIT_Q24_MODP_FROM_AJ_G1_G3",
        "A_orientation_sign_relative_to_stored_Yabs": int(asign),
        "R": {
            "formula": "G1-G3",
            "x_exact": str(rx),
            "y_exact": str(ry),
            "modp_profile": rprof,
        },
        "q24_modp": {
            "formula": "AJ(S3)+2*(G1-G3)",
            "profile": qprof,
            "x_numerator_coefficients_low_to_high":
                [int(v) for v in RF(qx.numerator()).list()],
            "x_denominator_coefficients_low_to_high":
                [int(v) for v in RF(qx.denominator()).list()],
            "y_numerator_coefficients_low_to_high":
                [int(v) for v in RF(qy.numerator()).list()],
            "y_denominator_coefficients_low_to_high":
                [int(v) for v in RF(qy.denominator()).list()],
        },
    })

    # ---------------------------------------------------------------
    # Optional characteristic-zero closure if the Hensel A section exists.
    # ---------------------------------------------------------------
    if EXACT_AJ.exists():
        ex = json.loads(EXACT_AJ.read_text())
        assert ex["status"] == "PASS_EXACT_Q8_S3_DIRECT_SECTION"
        assert ex["verification"]["exact_weierstrass_identity"] is True

        Ze = RQ([QQ(v) for v in ex["section"]["Z_coefficients_low_to_high"]])
        Xe = RQ([QQ(v) for v in ex["section"]["X_coefficients_low_to_high"]])
        Ye = RQ([QQ(v) for v in ex["section"]["Y_coefficients_low_to_high"]])
        Aex0 = EQ(KQ(Xe)/KQ(Ze**2), KQ(Ye)/KQ(Ze**3))

        # Orient exact A by reduction modulo p to the selected modular AF.
        candidates = [(+1, Aex0), (-1, -Aex0)]
        exact_choice = None
        for esign, Pex in candidates:
            px, py = Pex.xy()
            if reduce_rf(px) == AF[0] and reduce_rf(py) == AF[1]:
                exact_choice = (esign, Pex)
                break
        if exact_choice is None:
            raise ArithmeticError(
                "exact Hensel AJ reduces to neither selected modular orientation"
            )

        esign, Aex = exact_choice
        Qex = Aex + 2*RQpoint
        assert Qex in EQ
        qxe, qye = Qex.xy()

        def qq_degrees(v):
            v = KQ(v)
            return [
                int(RQ(v.numerator()).degree()),
                int(RQ(v.denominator()).degree()),
            ]

        qxdeg = qq_degrees(qxe)
        qydeg = qq_degrees(qye)
        assert qxdeg == [52, 48]
        assert qydeg == [78, 72]

        exact_payload = {
            "schema": "elkies-k3.h92-q8-q24-horizontal-section-qq.v1",
            "status": "PASS_EXACT_Q24_HORIZONTAL_SECTION",
            "zero": "II*_E8_1_branch_anchor",
            "formula": "Q24 = AJ(S3) + 2*(G1-G3)",
            "AJ_exact_orientation_sign": int(esign),
            "profile": {
                "P_dot_O": 24,
                "height": "52",
                "D13_local_correction": "0",
                "x_degrees": qxdeg,
                "y_degrees": qydeg,
            },
            "section": {
                "x_numerator_coefficients_low_to_high":
                    [str(v) for v in RQ(qxe.numerator()).list()],
                "x_denominator_coefficients_low_to_high":
                    [str(v) for v in RQ(qxe.denominator()).list()],
                "y_numerator_coefficients_low_to_high":
                    [str(v) for v in RQ(qye.numerator()).list()],
                "y_denominator_coefficients_low_to_high":
                    [str(v) for v in RQ(qye.denominator()).list()],
                "x": str(qxe),
                "y": str(qye),
            },
            "R": {
                "formula": "G1-G3",
                "x": str(rx),
                "y": str(ry),
            },
            "verification": {
                "exact_weierstrass_identity": True,
                "reduces_to_mod_100003_q24": (
                    reduce_rf(qxe) == qx and reduce_rf(qye) == qy
                ),
            },
            "boundary": (
                "This gives the exact q24 horizontal point on the canonical D13 "
                "q8 Jacobian. The next task is to build the degree-two neighbour "
                "pencil and derive the D12 Weierstrass equation."
            ),
        }
        assert exact_payload["verification"]["reduces_to_mod_100003_q24"]
        EXACT_Q24.write_text(json.dumps(exact_payload, indent=2, sort_keys=True) + "\n")

        payload["exact_q24_output"] = str(EXACT_Q24.relative_to(ROOT))
        payload["exact_q24_status"] = exact_payload["status"]

        print(f"EXACT_OUTPUT|{EXACT_Q24}", flush=True)
        print(
            "Q8Q24BACKTRACK_Q24_QQ|"
            f"x={qxdeg[0]}/{qxdeg[1]}|y={qydeg[0]}/{qydeg[1]}|"
            f"AJ_exact_sign={esign:+d}|"
            "status=PASS_EXACT_Q24_HORIZONTAL_SECTION",
            flush=True,
        )
    else:
        payload["exact_q24_status"] = "WAITING_FOR_EXACT_AJ_HENSEL_ARTIFACT"
        print(
            "Q8Q24BACKTRACK_Q24_QQ|"
            "status=WAITING_FOR_EXACT_AJ_HENSEL_ARTIFACT",
            flush=True,
        )

else:
    payload["status"] = "PASS_BRANCH_IS_G1_MINUS_G2_NEEDS_G2_OR_POLYNOMIAL_R_RECOVERY"
    print(
        "Q8Q24BACKTRACK_R|"
        "R=G1-G2|status=NEEDS_G2_OR_DIRECT_HEIGHT4_POLYNOMIAL_SECTION_RECOVERY",
        flush=True,
    )

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q8Q24BACKTRACK_RESULT|"
    f"branch={branch}|"
    f"exact_q24={int(payload.get('exact_q24_status') == 'PASS_EXACT_Q24_HORIZONTAL_SECTION')}|"
    f"status={payload['status']}",
    flush=True,
)
