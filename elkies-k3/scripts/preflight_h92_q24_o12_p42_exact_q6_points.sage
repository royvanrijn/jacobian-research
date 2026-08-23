#!/usr/bin/env sage -python
"""
Extract the current-equation q24 D12 zero O12 and marked orbit42 section P42
from the local bridge producer, identify their exact q6 MW words, materialize
them on the exact q6 Weierstrass model when possible, and independently replay
the q8 RR parameter on each curve.

This deliberately avoids the pointed-q24 zero-pole group law.  O12 and P42
are q6 sections / q8 multisections in the common H3 NS, and the bridge already
certifies q8 degrees 30 and 48 respectively.

PASS here gives explicit rational q6 points and exact U(T) maps for both
curves.  That is the correct input for mapping the curves through the q8
2-cover/covariant and then evaluating the q24 RR pencil.
"""

import contextlib
import io
import json
import sys
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, vector

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "elkies-k3" / "scripts"
LOCAL = ROOT / "artifacts" / "local" / "elkies-k3"
GEN = ROOT / "artifacts" / "generated-results"

BRIDGE_JSON = LOCAL / "q24-orbit42-current-equation-bridge.json"
Q6_WORDS = SCRIPTS / "search_h92_q24_bridge_equation_frame.sage"
Q8_CERT = SCRIPTS / "certify_h92_q8_equation_ns_divisor.sage"
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

for path in (BRIDGE_JSON, Q6_WORDS, Q8_CERT, CHILD, ZERO, COMP, S3BRIDGE, TRANSLATION):
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


def run_scope(path, argv=(), allow_assert=False):
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
    finally:
        sys.argv = saved
    return scope, buf.getvalue()


# -------------------------------------------------------------------------
# 1. Authoritative q6 MW-word machinery in common H3 NS.
# -------------------------------------------------------------------------
words, words_log = run_scope(Q6_WORDS)
for key in ("ns", "F6", "F8", "O6", "word_of", "section_from_word"):
    if key not in words:
        raise SystemExit(f"q6 word scope missing {key}")

ns = words["ns"]
F6 = vector(ZZ, words["F6"])
F8 = vector(ZZ, words["F8"])
O6 = vector(ZZ, words["O6"])
word_of = words["word_of"]
section_from_word = words["section_from_word"]

# -------------------------------------------------------------------------
# 2. Find the local producer that generated the passing current O42 bridge.
# -------------------------------------------------------------------------
producers = []
for path in SCRIPTS.glob("*.sage"):
    try:
        txt = path.read_text()
    except Exception:
        continue
    if "Q24O42EQ_RESULT|" in txt:
        producers.append(path)

if not producers:
    raise SystemExit(
        "Could not locate the local producer containing Q24O42EQ_RESULT|. "
        "Keep the script that generated q24-orbit42-current-equation-bridge.json "
        "in elkies-k3/scripts/ and rerun."
    )

# Prefer a producer mentioning the exact bridge filename/current-equation wording.
producers.sort(
    key=lambda p: (
        0 if "q24-orbit42-current-equation-bridge" in p.read_text() else 1,
        0 if "current" in p.name.lower() else 1,
        p.name,
    )
)
producer = producers[0]
ps, producer_log = run_scope(producer, allow_assert=True)

print(
    "Q24O42MAP_PRODUCER|"
    f"script={producer.relative_to(ROOT)}|candidates={len(producers)}|status=PASS",
    flush=True,
)


def iter_objects(obj, path="", seen=None, depth=0):
    if seen is None:
        seen = set()
    if depth > 8:
        return
    oid = id(obj)
    if oid in seen:
        return
    seen.add(oid)

    # Sage/integer row vector candidate.
    try:
        v = vector(ZZ, obj)
        if len(v) == 19:
            yield path, v
            return
    except Exception:
        pass

    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from iter_objects(v, f"{path}/{k}", seen, depth+1)
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            yield from iter_objects(v, f"{path}/{i}", seen, depth+1)


candidates = {}
for name, obj in ps.items():
    if name.startswith("__"):
        continue
    for subpath, v in iter_objects(obj, name):
        key = tuple(map(int, v))
        candidates.setdefault(key, {"vector":v, "paths":[]})["paths"].append(subpath)

# Also inspect serialized bridge content.
for subpath, v in iter_objects(bridge, "bridge"):
    key = tuple(map(int, v))
    candidates.setdefault(key, {"vector":v, "paths":[]})["paths"].append(subpath)

sections = []
for entry in candidates.values():
    C = entry["vector"]
    try:
        square = ZZ(C * ns * C)
        d6 = ZZ(C * ns * F6)
        d8 = ZZ(C * ns * F8)
    except Exception:
        continue
    if square == -2 and d6 == 1 and d8 in (30, 48):
        sections.append({
            "class":C,
            "q8_degree":int(d8),
            "paths":sorted(set(entry["paths"])),
        })

# Deduplicate by class.
uniq = {}
for rec in sections:
    uniq[tuple(map(int,rec["class"]))] = rec
sections = list(uniq.values())

for rec in sections:
    print(
        "Q24O42MAP_CLASS|"
        f"q8_degree={rec['q8_degree']}|"
        f"paths={';'.join(rec['paths'][:8])}|status=CANDIDATE",
        flush=True,
    )

by_degree = {30:[],48:[]}
for rec in sections:
    by_degree[rec["q8_degree"]].append(rec)

if len(by_degree[30]) != 1 or len(by_degree[48]) != 1:
    # Names/paths may allow an unambiguous choice if duplicate classes appeared
    # through different coordinate containers, but class dedup has already run.
    raise SystemExit(
        "Could not uniquely identify O12/P42 from the producer scope: "
        f"degree30={len(by_degree[30])}, degree48={len(by_degree[48])}. "
        "Paste the Q24O42EQ producer filename/output and we can pin its exact variables."
    )

targets = {
    "O12": by_degree[30][0],
    "P42": by_degree[48][0],
}

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
    "bridge_producer":str(producer.relative_to(ROOT)),
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
