#!/usr/bin/env sage -python
"""Exact rejection audit for the proposed orbit42 fast-q6 transport.

The bridge stores O12/P42 in the current equation-D13 basis.  Convert them
with the named unimodular equation-D13 -> ambient H3-NS matrix and certify
their actual q6/q8 degrees.  They have q6 degrees 435 and 703, not degree one,
so neither is a q6 rational section and the proposed point-transport route is
invalid.  A passing result here is a negative audit, not an equation lift.
"""

import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, vector, matrix,
    block_diagonal_matrix,
)

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "elkies-k3" / "scripts"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"
GEN = ROOT / "artifacts" / "generated-results"

BRIDGE_JSON = LOCAL / "q24-orbit42-current-equation-bridge.json"
Q6_WORDS = SCRIPTS / "search_h92_q24_bridge_equation_frame.sage"
Q8_CERT = SCRIPTS / "certify_h92_q8_equation_ns_divisor.sage"
Q24_TRANSLATION = LOCAL / "q8-q24-physical-to-equation-translation.json"
Q8_CHILD_CANDS = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
CHILD = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMP = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
S3BRIDGE = LOCAL / "q6-third-to-q8-bridge.json"
TRANSLATION = LOCAL / "q6-standard-zero-translation.json"
OUT = LOCAL / "q24-o12-p42-q6-preflight.json"

for path in (
    BRIDGE_JSON, Q6_WORDS, Q8_CERT, Q24_TRANSLATION,
    CHILD, ZERO, COMP, S3BRIDGE, TRANSLATION,
):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

Q8_CHILD = next(
    (
        p for p in Q8_CHILD_CANDS
        if p.exists()
        and json.loads(p.read_text()).get("status") == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
    ),
    None,
)
if Q8_CHILD is None:
    raise SystemExit("missing complete corrected q8 child artifact")

bridge = json.loads(BRIDGE_JSON.read_text())


def run_scope(path, argv=(), allow_assert=False, allow_clean_exit=False):
    saved = list(sys.argv)
    scope = {"__name__":"__embedded__", "__file__":str(path)}
    buf = io.StringIO()
    try:
        sys.argv = [str(path)] + list(argv)
        with contextlib.redirect_stdout(buf):
            try:
                exec(compile(path.read_text(), str(path), "exec"), scope)
            except AssertionError:
                if not allow_assert:
                    raise
            except SystemExit as exc:
                if not allow_clean_exit or exc.code not in (None, 0):
                    raise
    finally:
        sys.argv = saved
    return scope, buf.getvalue()


# -------------------------------------------------------------------------
# 1. Authoritative q6 MW-word machinery in common H3 NS.
# -------------------------------------------------------------------------
words, words_log = run_scope(
    Q6_WORDS,
    ("--setup-only",),
    allow_clean_exit=True,
)
for key in (
    "ns", "F6", "F8", "O6", "word_of", "section_from_word",
    "Badapt", "G13",
):
    if key not in words:
        raise SystemExit(f"q6 word scope missing {key}")

ns = words["ns"]
F6 = vector(ZZ, words["F6"])
F8 = vector(ZZ, words["F8"])
O6 = vector(ZZ, words["O6"])
word_of = words["word_of"]
section_from_word = words["section_from_word"]

# -------------------------------------------------------------------------
# 2. Exact NAMED current equation-D13 -> ambient H3-NS coordinate bridge.
#
# q24-orbit42-current-equation-bridge.json stores O12/P42 in the CURRENT
# equation-D13 basis. q6 word_of() uses the common ambient H3 NS.
#
# Exact conversion:
#
#     C_ambient = C_equation_D13 * Badapt
#
# A Gram-compatible D13 basis is not sufficient here: D13 has nontrivial
# marking automorphisms.  Use the exact Badapt constructed in Q6_WORDS, which
# is the named equation frame used by the bridge producer, and pin it by the
# q24-fibre round-trip to the independent physical->equation certificate.
# -------------------------------------------------------------------------
BADAPT_CACHE = LOCAL / "q24-equation-d13-badapt-cache.json"
Gambient = matrix(ZZ, ns)
Badapt = matrix(ZZ, words["Badapt"])
Geq_frame = matrix(ZZ, words["G13"])

if Badapt.dimensions() != (19,19) or abs(ZZ(Badapt.det())) != 1:
    raise ArithmeticError(
        f"Badapt dimensions/determinant invalid: "
        f"{Badapt.dimensions()}, det={Badapt.det()}"
    )

U2c = matrix(ZZ, ((0,1),(1,0)))
Geq = block_diagonal_matrix(U2c, -Geq_frame)
if Badapt * Gambient * Badapt.transpose() != Geq:
    raise ArithmeticError(
        "Badapt * Gambient * Badapt^t != named equation-D13 Gram"
    )

if vector(ZZ, Badapt.row(0)) != F8:
    raise ArithmeticError("named equation-D13 fibre row does not equal F8")

q24_translation = json.loads(Q24_TRANSLATION.read_text())
if q24_translation.get("status") != "PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION":
    raise SystemExit(
        "q24 physical->equation translation is not passing: "
        + str(q24_translation.get("status"))
    )

bridge_q24 = vector(ZZ, bridge["current_equation_D13"]["q24_fibre"])
expected_q24_ambient = vector(
    ZZ,
    q24_translation["q24_equation"]["equation_divisor_source_h3_ns"],
)
mapped_q24_ambient = vector(ZZ, bridge_q24 * Badapt)
if mapped_q24_ambient != expected_q24_ambient:
    raise ArithmeticError(
        "named equation-D13 Badapt fails the exact q24-fibre round-trip"
    )

BADAPT_CACHE.write_text(json.dumps({
    "schema": "elkies-k3.q24-equation-d13-badapt-cache.v2",
    "status": "PASS_Q24_NAMED_EQUATION_D13_BADAPT_CACHE",
    "source": str(Q6_WORDS.relative_to(ROOT)),
    "anchor": str(Q24_TRANSLATION.relative_to(ROOT)),
    "ambient_ns": [[int(x) for x in row] for row in Gambient.rows()],
    "Badapt": [[int(x) for x in row] for row in Badapt.rows()],
    "equation_d13_frame": [
        [int(x) for x in row] for row in Geq_frame.rows()
    ],
    "q24_fibre_equation_d13": list(map(int, bridge_q24)),
    "q24_fibre_ambient": list(map(int, mapped_q24_ambient)),
}, indent=2, sort_keys=True) + "\n")

print(
    "Q24O42MAP_COORDS|"
    "method=NAMED_Q6_EQUATION_FRAME|"
    f"det={Badapt.det()}|gram=PASS|q8_fibre=PASS|"
    "q24_fibre_roundtrip=PASS|cache=v2|status=PASS",
    flush=True,
)


if bridge.get("status") != "PASS_Q24_ORBIT42_CURRENT_EQUATION_LATTICE_BRIDGE":
    raise SystemExit(
        "q24 orbit42 current-equation bridge is not passing: "
        + str(bridge.get("status"))
    )

eq = bridge.get("current_equation_D13", {})
field_for = {
    "O12": "historical_D12_zero",
    "P42": "orbit42_marked_section",
}
missing = [field for field in field_for.values() if field not in eq]
if missing:
    raise SystemExit(
        "current-equation bridge missing direct class fields: "
        + ",".join(missing)
    )

direct_eq = {
    label: vector(ZZ, eq[field])
    for label, field in field_for.items()
}

direct = {
    label: vector(ZZ, row * Badapt)
    for label, row in direct_eq.items()
}

targets = {}
expected_q6_degree = {"O12": 435, "P42": 703}
expected_degree = {"O12": 30, "P42": 48}

for label, C in direct.items():
    if len(C) != 19:
        raise ArithmeticError(
            f"{label} converted class dimension={len(C)}, expected 19"
        )

    square = ZZ(C * ns * C)
    d6 = ZZ(C * ns * F6)
    d8 = ZZ(C * ns * F8)

    if square != -2:
        raise ArithmeticError(
            f"{label} converted class square={square}, expected -2"
        )
    if d6 != expected_q6_degree[label]:
        raise ArithmeticError(
            f"{label} converted q6 degree={d6}, "
            f"expected {expected_q6_degree[label]}"
        )
    if d8 != expected_degree[label]:
        raise ArithmeticError(
            f"{label} converted q8 degree={d8}, "
            f"expected {expected_degree[label]}"
        )

    U2c = matrix(ZZ, ((0,1),(1,0)))
    Geq = block_diagonal_matrix(U2c, -Geq_frame)
    eq_square = ZZ(direct_eq[label] * Geq * direct_eq[label])
    if eq_square != -2:
        raise ArithmeticError(
            f"{label} equation-D13 square={eq_square}, expected -2"
        )

    targets[label] = {
        "class": C,
        "q6_degree": int(d6),
        "q8_degree": int(d8),
        "paths": [
            "bridge/current_equation_D13/"
            + field_for[label]
            + " -> Badapt -> ambient"
        ],
    }

    print(
        "Q24O42MAP_CLASS|"
        f"curve={label}|square={square}|q6_degree={d6}|q8_degree={d8}|"
        "source=DIRECT_BRIDGE_VIA_BADAPT|status=PASS",
        flush=True,
    )

print(
    "Q24O42MAP_PRODUCER|"
    "script=SKIPPED|"
    "method=DIRECT_BRIDGE_VIA_BADAPT|"
    "status=PASS",
    flush=True,
)

# Terminal negative gate.  Everything below is retained only as the historical
# implementation of the disproved degree-one premise and is deliberately
# unreachable.
payload = {
    "schema": "elkies-k3.h3-q24-orbit42-fast-q6-premise-audit.v1",
    "status": "PASS_Q42_FAST_Q6_PREMISE_REJECTION",
    "bridge": str(BRIDGE_JSON.relative_to(ROOT)),
    "coordinate_conversion": {
        "method": "NAMED_EQUATION_D13_BADAPT",
        "determinant": int(Badapt.det()),
        "gram_verified": True,
        "q24_fibre_roundtrip_verified": True,
        "cache": str(BADAPT_CACHE.relative_to(ROOT)),
    },
    "targets": {
        label: {
            "square": -2,
            "q6_degree": rec["q6_degree"],
            "q8_degree": rec["q8_degree"],
            "equation_D13_class": list(map(int, direct_eq[label])),
            "ambient_H3_NS_class": list(map(int, rec["class"])),
        }
        for label, rec in targets.items()
    },
    "conclusion": (
        "O12 and P42 are high-degree q6 multisections, not q6 rational "
        "sections. Exact q6 MW-word extraction and rational-point transport "
        "are therefore inapplicable."
    ),
    "next": "Q42_RESOLVED_RR_TRIVIALIZATION",
    "proof_boundary": (
        "Exact coordinate and intersection-degree audit only. It disproves "
        "the fast q6-point premise and does not construct the orbit42 pencil."
    ),
}
OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUT}", flush=True)
print(
    "Q24O42MAP_RESULT|"
    "O12_q6_degree=435|P42_q6_degree=703|"
    "next=Q42_RESOLVED_RR_TRIVIALIZATION|"
    "status=PASS_Q42_FAST_Q6_PREMISE_REJECTION",
    flush=True,
)
raise SystemExit(0)

# Exact q6 words relative to O6 (= old zero), then standard-Weierstrass words.
z_old_std = vector(
    ZZ,
    json.loads(TRANSLATION.read_text())["standard_MW_coordinates"]["old_zero"],
)
assert z_old_std == vector(ZZ,(2,-1,0))

for label, rec in targets.items():
    C = rec["class"]
    w_old = vector(ZZ, word_of(C))
    roundtrip = section_from_word(w_old)[0]
    assert vector(ZZ, roundtrip) == C
    w_std = w_old + z_old_std
    rec["word_old"] = w_old
    rec["word_standard"] = w_std

    print(
        "Q24O42MAP_WORD|"
        f"curve={label}|q8_degree={rec['q8_degree']}|"
        f"old_mw={','.join(map(str,w_old))}|"
        f"standard_mw={','.join(map(str,w_std))}|"
        f"even12={int(w_std[0]%2==0 and w_std[1]%2==0)}|"
        "status=PASS_EXACT_Q6_WORD",
        flush=True,
    )

# -------------------------------------------------------------------------
# 3. Exact q6 rational point basis Pmap,Qmap,S3.
# -------------------------------------------------------------------------
child = json.loads(CHILD.read_text())
zero = json.loads(ZERO.read_text())
components = json.loads(COMP.read_text())
s3bridge = json.loads(S3BRIDGE.read_text())
q8 = json.loads(Q8_CHILD.read_text())

assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert components["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert s3bridge["status"] == "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"

R = PolynomialRing(QQ, "T")
T = R.gen()
K = R.fraction_field()


def poly(values):
    return R([QQ(v) for v in values])


def rational(data, nk, dk):
    return K(poly(data[nk])) / K(poly(data[dk]))


model = child["minimal_short_weierstrass"]
A = poly(model["A_coefficients_low_to_high"])
B = poly(model["B_coefficients_low_to_high"])
E = EllipticCurve(K,[0,0,0,K(A),K(B)])

zdata = zero["section"]
old_zero_point = E(
    rational(
        zdata,
        "x_numerator_coefficients_low_to_high",
        "x_denominator_coefficients_low_to_high",
    ),
    rational(
        zdata,
        "y_numerator_coefficients_low_to_high",
        "y_denominator_coefficients_low_to_high",
    ),
)

pts = {}
for entry in components["sections"]:
    pts[entry["sign"]] = E(
        rational(
            entry,
            "x_numerator_coefficients_low_to_high",
            "x_denominator_coefficients_low_to_high",
        ),
        rational(
            entry,
            "y_numerator_coefficients_low_to_high",
            "y_denominator_coefficients_low_to_high",
        ),
    )

e77 = pts[components["source"]["E7_7_sign"]]
affine = pts[components["source"]["affine_E7_sign"]]
Pmap = e77 - old_zero_point
Qmap = e77 - affine

s3data = s3bridge["third_section_canonical_q6"]
S3 = E(
    rational(
        s3data["x"],
        "numerator_coefficients_low_to_high",
        "denominator_coefficients_low_to_high",
    ),
    rational(
        s3data["y"],
        "numerator_coefficients_low_to_high",
        "denominator_coefficients_low_to_high",
    ),
)

# -------------------------------------------------------------------------
# 4. Reconstruct corrected q8 RR pencil, then restrict to O12/P42.
# -------------------------------------------------------------------------
mdata = q8["marking"]["section"]
sx = rational(
    mdata,
    "x_numerator_coefficients_low_to_high",
    "x_denominator_coefficients_low_to_high",
)
sy = rational(
    mdata,
    "y_numerator_coefficients_low_to_high",
    "y_denominator_coefficients_low_to_high",
)
Smark = E(sx,sy)
assert Smark == Pmap + Qmap


def monic_power_root(value, exponent):
    value = R(value)
    out = R.one()
    for fac, mult in value.factor():
        assert int(mult) % exponent == 0
        out *= fac.monic() ** (int(mult)//exponent)
    return out.monic()


nx,dx = R(sx.numerator()),R(sx.denominator())
ny,dy = R(sy.numerator()),R(sy.denominator())
h = monic_power_root(dx,2)
assert h == monic_power_root(dy,3)

ii = R(next(x for x in child["finite_fibres"] if x["kodaira"]=="II*")["factor"]).monic()
iv = R(next(x for x in child["finite_fibres"] if x["kodaira"]=="IV*")["factor"]).monic()
M = (ii**2 * iv**2).monic()

normalizer = (ny * dx * (h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0
p_fun = -sy/sx
rho = (normalizer * nx.inverse_mod(M)).mod(M)

pairs = []
for entry in q8["rr"]["kernel_polynomials"]:
    sp = R(entry["s"])
    tp = R(entry["t"])
    Bcoef = K(sp)/K(h)
    Acoef = (
        -K(sp)*p_fun/K(h)
        - K(sp)*K(normalizer)/K(nx)
        + K(sp*rho)
        + K(tp*M)
    )
    pairs.append((Acoef,Bcoef))
assert len(pairs)==2
(A0,B0),(A1,B1)=pairs


def rf_record(v):
    v = K(v)
    return {
        "num":[str(x) for x in R(v.numerator()).list()],
        "den":[str(x) for x in R(v.denominator()).list()],
        "num_degree":int(R(v.numerator()).degree()),
        "den_degree":int(R(v.denominator()).degree()),
    }


all_materialized = True
for label, rec in targets.items():
    a,b,c = map(ZZ,rec["word_standard"])
    if a % 2 or b % 2:
        rec["materialized"] = False
        all_materialized = False
        print(
            "Q24O42MAP_POINT|"
            f"curve={label}|standard_mw={a},{b},{c}|"
            "reason=primitive_halving_needed|status=NEEDS_NEW_Q6_SECTION",
            flush=True,
        )
        continue

    kp = -a//2
    kq = -b//2
    ks = c
    P = kp*Pmap + kq*Qmap + ks*S3
    if P.is_zero():
        raise ArithmeticError(f"{label} unexpectedly equals standard zero")
    px,py = P.xy()
    assert py**2 == px**3 + K(A)*px + K(B)

    # Degree against current equation q8 fibre via actual RR restriction.
    mP = (py+sy)/(px-sx)
    UP = K((A1+B1*mP)/(A0+B0*mP))
    un = R(UP.numerator())
    ud = R(UP.denominator())
    common = un.gcd(ud)
    if common.degree() > 0:
        un //= common
        ud //= common
    degree = max(un.degree(),ud.degree())
    assert degree == rec["q8_degree"], (label,degree,rec["q8_degree"])

    rec["materialized"] = True
    rec["point_coefficients"] = [int(kp),int(kq),int(ks)]
    rec["x"] = rf_record(px)
    rec["y"] = rf_record(py)
    rec["q8_parameter"] = rf_record(K(un)/K(ud))
    rec["q8_parameter_degree"] = int(degree)

    print(
        "Q24O42MAP_POINT|"
        f"curve={label}|formula={kp}Pmap+{kq}Qmap+{ks}S3|"
        f"q8_degree={degree}|"
        f"U_numdeg={un.degree()}|U_dendeg={ud.degree()}|"
        "status=PASS_EXACT_Q6_POINT_AND_Q8_RESTRICTION",
        flush=True,
    )

payload = {
    "schema":"elkies-k3.h3-q24-o12-p42-q6-preflight.v1",
    "status":(
        "PASS_Q24_O12_P42_EXACT_Q6_POINTS"
        if all_materialized
        else "Q24_O12_P42_NEEDS_PRIMITIVE_Q6_SECTION_RECOVERY"
    ),
    "bridge":str(BRIDGE_JSON.relative_to(ROOT)),
    "bridge_producer":"SKIPPED_DIRECT_BRIDGE_VIA_BADAPT",
    "targets":{
        label:{
            "q8_degree":rec["q8_degree"],
            "source_paths":rec["paths"],
            "class":[int(x) for x in rec["class"]],
            "q6_old_zero_mw":[int(x) for x in rec["word_old"]],
            "q6_standard_mw":[int(x) for x in rec["word_standard"]],
            "materialized":bool(rec.get("materialized",False)),
            **(
                {
                    "coefficients_on_Pmap_Qmap_S3":rec["point_coefficients"],
                    "x":rec["x"],
                    "y":rec["y"],
                    "q8_parameter":rec["q8_parameter"],
                    "q8_parameter_degree":rec["q8_parameter_degree"],
                }
                if rec.get("materialized") else {}
            ),
        }
        for label,rec in targets.items()
    },
    "next":(
        "Map the explicit q6 curves through the corrected q8 2-cover/covariant, "
        "then evaluate the q24 RR basis along them. Their q24 degree is one, "
        "so the resulting q24 base parameter must be Mobius in the q6 parameter."
        if all_materialized
        else
        "Recover the missing primitive q6 point(s) indicated above before "
        "attempting q8/q24 mapping."
    ),
}
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24O42MAP_RESULT|"
    f"O12_materialized={int(targets['O12'].get('materialized',False))}|"
    f"P42_materialized={int(targets['P42'].get('materialized',False))}|"
    f"status={payload['status']}",
    flush=True,
)
