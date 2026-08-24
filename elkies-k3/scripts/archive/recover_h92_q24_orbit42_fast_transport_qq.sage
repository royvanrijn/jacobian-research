#!/usr/bin/env sage -python
"""
status: HISTORICAL_DIAGNOSTIC
claim: rejected fast transport; O12/P42 have q6 degrees 435/703 and are not
       q6 rational sections.

Strategy
--------
Avoid a new RR search entirely.

1. Use the existing exact q6 preflight to materialize BOTH geometric q24
   degree-one curves O12 and P42 as exact q6 rational points and exact q8
   base maps.
2. Map each curve algebraically through the already-certified q8 branch
   quartic -> D13 birational map.
3. Evaluate the already-certified exact q24 2x56 pencil on the D13 point.
   Since O12/P42 have q24 degree one, require a Mobius q24 base map.
4. Map the resulting q24 quartic point to the exact D12 Jacobian, convert
   from the raw invariant Jacobian to the certified minimal D12 model, and
   recover each section globally over QQ(V).
5. Subtract O12 from P42 on the exact D12 model.  Require the corrected
   marked section profile deg Z=3 (P.O=3).
6. Compile the degree-two chord pencil directly.  The cleared discriminant
   has old-base degree <= 20, so use a tiny gcd-only Yun decomposition.
7. Compute the exact Jacobian and require A11/MW6.

No Hensel lifting, no 56D RR solve, no primary decomposition, no large
degree-160 GCD.

Output:
  artifacts/local/elkies-k3/q24-d12-orbit42-fast-transport-qq.json

Terminal:
  Q42FAST_RESULT|...|status=PASS_Q42_FAST_A11_CANDIDATES
"""

import argparse
import json
import subprocess
import time
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, matrix, sage_eval
)


# ---------------------------------------------------------------------------
# Paths.
# ---------------------------------------------------------------------------
def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents/jacobian-research",
        home / "jacobian-research",
        home / "src/jacobian-research",
        home / "git/jacobian-research",
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
            and (c / "artifacts/local/elkies-k3").is_dir()
        ):
            return c
    raise SystemExit("Could not locate jacobian-research")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--output", type=Path)
parser.add_argument(
    "--no-run-preflight",
    action="store_true",
    help="do not auto-run the O12/P42 exact-q6 preflight if its artifact is absent",
)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
LOCAL = ROOT / "artifacts/local/elkies-k3"
GEN = ROOT / "artifacts/generated-results"
SCRIPTS = ROOT / "elkies-k3/scripts"

PREFLIGHT = LOCAL / "q24-o12-p42-q6-preflight.json"
PREFLIGHT_SCRIPT = SCRIPTS / "preflight_h92_q24_o12_p42_exact_q6_points.sage"
Q8_CHILD_CANDIDATES = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8_ANCHOR = LOCAL / "q8-d13-branch-anchor.json"
Q24_PARENT = LOCAL / "q24-d13-to-d12-component-valuation-qq.json"
Q24_SECTION = LOCAL / "q8-q24-horizontal-section-qq.json"
CORE = SCRIPTS / "elliptic_neighbor_compiler.sage"

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / "q24-d12-orbit42-fast-transport-qq.json"
)

if not PREFLIGHT.exists() and not args.no_run_preflight:
    if not PREFLIGHT_SCRIPT.exists():
        raise SystemExit(f"missing preflight producer: {PREFLIGHT_SCRIPT}")
    print(
        "Q42FAST|stage=Q6_PREFLIGHT|status=BEGIN",
        flush=True,
    )
    subprocess.run(
        ["sage", "-python", str(PREFLIGHT_SCRIPT)],
        cwd=str(ROOT),
        check=True,
    )

for path in (PREFLIGHT, Q8_ANCHOR, Q24_PARENT, Q24_SECTION, CORE):
    if not path.exists():
        raise SystemExit(f"missing prerequisite: {path}")

Q8_CHILD = next(
    (
        p for p in Q8_CHILD_CANDIDATES
        if p.exists()
        and json.loads(p.read_text()).get("status")
        == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
        and "branch_quartic" in json.loads(p.read_text()).get("pencil", {})
    ),
    None,
)
if Q8_CHILD is None:
    raise SystemExit("missing complete exact q8/D13 child artifact")

pre = json.loads(PREFLIGHT.read_text())
q8 = json.loads(Q8_CHILD.read_text())
anchor = json.loads(Q8_ANCHOR.read_text())
q24 = json.loads(Q24_PARENT.read_text())
h24 = json.loads(Q24_SECTION.read_text())

if pre.get("status") != "PASS_Q24_O12_P42_EXACT_Q6_POINTS":
    print(
        "Q42FAST_BLOCKED|"
        f"preflight_status={pre.get('status')}|"
        "reason=O12_OR_P42_NOT_MATERIALIZED_AS_EXACT_Q6_POINT|"
        "status=NEEDS_PRIMITIVE_Q6_SECTION_RECOVERY",
        flush=True,
    )
    raise SystemExit(2)

for label in ("O12", "P42"):
    if not pre["targets"][label].get("materialized"):
        raise SystemExit(f"{label} is not materialized in exact q6 preflight")

assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert anchor["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"
assert q24["status"] == "PASS_EXACT_Q24_D13_TO_D12_COMPONENT_VALUATION_RR"
assert h24["status"] == "PASS_EXACT_Q24_HORIZONTAL_SECTION"

started = time.monotonic()


def log(stage, **fields):
    tail = "|".join(f"{k}={v}" for k, v in fields.items())
    print(
        f"Q42FAST|stage={stage}|elapsed={time.monotonic()-started:.3f}"
        + (f"|{tail}" if tail else ""),
        flush=True,
    )


# ---------------------------------------------------------------------------
# Exact arithmetic helpers.
# ---------------------------------------------------------------------------
def poly_from_list(ring, values):
    return ring([QQ(v) for v in values])


def normalize_rf(value, ring):
    value = value.parent()(value)
    n = ring(value.numerator())
    d = ring(value.denominator())
    g = n.gcd(d)
    if g.degree() > 0:
        n //= g
        d //= g
    if not d:
        raise ZeroDivisionError
    lc = QQ(d.leading_coefficient())
    n = ring(n / lc)
    d = ring(d / lc)
    return value.parent()(n) / value.parent()(d)


def eval_poly(poly, arg):
    out = arg.parent()(0)
    for c in reversed(list(poly)):
        out = out * arg + arg.parent()(c)
    return out


def eval_rf(value, arg, source_ring):
    value = value.parent()(value)
    n = source_ring(value.numerator())
    d = source_ring(value.denominator())
    return eval_poly(n, arg) / eval_poly(d, arg)


def monic_power_root(poly, exponent):
    """
    Exact coefficient-recursion root of a monic perfect power.
    No factorization.
    """
    ring = poly.parent()
    poly = ring(poly)
    if not poly:
        return ring.zero()
    if poly.leading_coefficient() != 1:
        raise ArithmeticError("monic_power_root requires monic input")
    if poly.degree() % exponent:
        raise ArithmeticError(
            f"degree {poly.degree()} not divisible by {exponent}"
        )
    d = poly.degree() // exponent
    coeffs = [QQ(0)] * d + [QQ(1)]
    for k in range(d - 1, -1, -1):
        current = ring(coeffs)
        target_degree = (exponent - 1) * d + k
        known = (current**exponent)[target_degree]
        coeffs[k] = QQ(poly[target_degree] - known) / QQ(exponent)
    root = ring(coeffs)
    if root**exponent != poly:
        raise ArithmeticError(
            f"coefficient-recursion {exponent}-th root failed"
        )
    return root


def fast_rational_square_root(value, ring):
    """
    Exact square root in Frac(ring), using coefficient recursion rather than
    factorization.  Returns None when the rational function is not a square.
    """
    field = ring.fraction_field()
    value = field(value)
    if not value:
        return field.zero()

    n = ring(value.numerator())
    d = ring(value.denominator())
    g = n.gcd(d)
    if g.degree() > 0:
        n //= g
        d //= g

    # Normalize denominator to monic.  If n/d is a square, both resulting
    # polynomials are squares over QQ.
    lc = QQ(d.leading_coefficient())
    n = ring(n / lc)
    d = ring(d / lc)

    if d.leading_coefficient() != 1:
        return None

    try:
        rd = monic_power_root(d, 2)
    except ArithmeticError:
        return None

    # Normalize numerator by its leading coefficient.
    nlc = QQ(n.leading_coefficient())
    if nlc < 0:
        return None
    nn = ZZ(nlc.numerator())
    nd = ZZ(nlc.denominator())
    if not nn.is_square() or not nd.is_square():
        return None
    root_lc = QQ(nn.sqrt()) / QQ(nd.sqrt())
    nmonic = ring(n / nlc)
    if nmonic.leading_coefficient() != 1:
        return None
    try:
        rn0 = monic_power_root(nmonic, 2)
    except ArithmeticError:
        return None
    rn = ring(root_lc * rn0)

    root = field(rn) / field(rd)
    if root**2 != value:
        return None
    return root


def rf_record(value, ring):
    value = normalize_rf(value, ring)
    return {
        "num": [str(v) for v in ring(value.numerator()).list()],
        "den": [str(v) for v in ring(value.denominator()).list()],
        "num_degree": int(ring(value.numerator()).degree()),
        "den_degree": int(ring(value.denominator()).degree()),
    }


def parse_record(rec, ring):
    field = ring.fraction_field()
    return (
        field(poly_from_list(ring, rec["num"]))
        / field(poly_from_list(ring, rec["den"]))
    )


def projectivize_section(xp, yp, A, B, ring, expected_z=None):
    field = ring.fraction_field()
    xp = normalize_rf(field(xp), ring)
    yp = normalize_rf(field(yp), ring)

    nx = ring(xp.numerator())
    dx = ring(xp.denominator())
    ny = ring(yp.numerator())
    dy = ring(yp.denominator())

    # normalize denominators monic
    lcx = QQ(dx.leading_coefficient())
    lcy = QQ(dy.leading_coefficient())
    nx = ring(nx / lcx)
    dx = ring(dx / lcx)
    ny = ring(ny / lcy)
    dy = ring(dy / lcy)

    if dx.leading_coefficient() != 1 or dy.leading_coefficient() != 1:
        raise ArithmeticError("section denominators did not normalize monic")

    Zx = monic_power_root(dx, 2)
    Zy = monic_power_root(dy, 3)
    if Zx != Zy:
        raise ArithmeticError(
            f"incompatible section denominators: Zx={Zx}, Zy={Zy}"
        )
    Z = Zx
    X = ring(nx)
    Y = ring(ny)

    if expected_z is not None and Z.degree() != expected_z:
        raise ArithmeticError(
            f"section Z degree={Z.degree()}, expected {expected_z}"
        )

    if Y**2 != X**3 + A*X*Z**4 + B*Z**6:
        raise ArithmeticError("projectivized section identity failed")
    return X, Y, Z


def homogeneous_compose(poly, num, den, weight):
    """
    den^weight * poly(num/den), assuming deg(poly)<=weight.
    """
    ring = poly.parent()
    poly = ring(poly)
    num = ring(num)
    den = ring(den)
    out = ring.zero()
    for i, c in enumerate(poly.list()):
        if c:
            out += QQ(c) * num**i * den**(weight - i)
    return ring(out)


# ---------------------------------------------------------------------------
# Rings and exact parent data.
# ---------------------------------------------------------------------------
RT = PolynomialRing(QQ, "T")
T = RT.gen()
KT = RT.fraction_field()

RU = PolynomialRing(QQ, "U")
U = RU.gen()
KU = RU.fraction_field()

RV = PolynomialRing(QQ, "V")
V = RV.gen()
KV = RV.fraction_field()

# D13 child = q24 parent.
A13 = poly_from_list(
    RU, q8["child"]["minimal_A_coefficients_low_to_high"]
)
B13 = poly_from_list(
    RU, q8["child"]["minimal_B_coefficients_low_to_high"]
)

# Certified minimal D12 parent.
A12 = poly_from_list(
    RV, q24["child"]["minimal_A_coefficients_low_to_high"]
)
B12 = poly_from_list(
    RV, q24["child"]["minimal_B_coefficients_low_to_high"]
)

# q24 marked section on D13.
s24 = h24["section"]
Z24 = poly_from_list(RU, s24["Z_coefficients_low_to_high"])
X24 = poly_from_list(RU, s24["X_coefficients_low_to_high"])
Y24 = poly_from_list(RU, s24["Y_coefficients_low_to_high"])
assert (Z24.degree(), X24.degree(), Y24.degree()) == (24, 52, 78)

# q24 exact 2x56 plane.
plane = q24["rr"]["plane_2x56"]
if len(plane) != 2 or any(len(row) != 56 for row in plane):
    raise ArithmeticError("q24 exact plane is not 2x56")

q24_pairs = []
for row in plane:
    AA = RU([QQ(v) for v in row[:41]])
    BB = RU([QQ(v) for v in row[41:]])
    q24_pairs.append((AA, BB))

# q24 quartic over QQ(V)[U].
q24_coeff = [
    KV(sage_eval(str(text), locals={"V": V}))
    for text in q24["quartic"]["coefficients_in_U_low_to_high"]
]
if len(q24_coeff) != 5:
    raise ArithmeticError("q24 quartic does not have five coefficients")

# old D13 I9* root, used as a known rational q24 quartic point.
i9 = next(
    item for item in q8["child"]["finite_fibres"]
    if item["kodaira"] == "I9*"
)
i9_factor = RU(sage_eval(str(i9["factor"]), locals={"U": U}))
alpha24 = QQ(-i9_factor[0] / i9_factor[1])

# q24 pointed quartic map, constructed ONCE over QQ(V).
qalpha = sum(q24_coeff[i] * KV(alpha24**i) for i in range(5))
q0 = fast_rational_square_root(qalpha, RV)
if q0 is None or not q0:
    raise ArithmeticError("exact old-I9* q24 quartic point is not a square")

qa0, qa1, qa2, qa3, qa4 = q24_coeff
a24 = qa4
b24 = qa3 + 4*alpha24*qa4
c24 = qa2 + 3*alpha24*qa3 + 6*alpha24**2*qa4
d24 = (
    qa1
    + 2*alpha24*qa2
    + 3*alpha24**2*qa3
    + 4*alpha24**3*qa4
)

a1w = d24 / q0
a2w = c24 - d24**2 / (4*q0**2)
a3w = 2*q0*b24
a4w = -4*q0**2*a24
a6w = a2w*a4w
b2w = a1w**2 + 4*a2w
b4w = 2*a4w + a1w*a3w
b6w = a3w**2 + 4*a6w
c4w = b2w**2 - 24*b4w
c6w = -b2w**3 + 36*b2w*b4w - 216*b6w
rawA_from_pointing = 81 * (-c4w / 48)
rawB_from_pointing = 729 * (-c6w / 864)

rawA = KV(sage_eval(str(q24["jacobian_raw"]["A"]), locals={"V": V}))
rawB = KV(sage_eval(str(q24["jacobian_raw"]["B"]), locals={"V": V}))

if rawA_from_pointing != rawA or rawB_from_pointing != rawB:
    raise ArithmeticError(
        "exact pointed-q24 quartic convention does not reproduce raw Jacobian"
    )

# Raw -> certified minimal D12 short-Weierstrass gauge.
cA = rawA / KV(A12)
cB = rawB / KV(B12)
scale2 = cB / cA
if scale2**2 != cA:
    raise ArithmeticError("raw/minimal D12 x-scaling relation failed")
scale3 = fast_rational_square_root(cB, RV)
if scale3 is None or scale3**2 != cB:
    raise ArithmeticError("raw/minimal D12 y-scaling is not an exact square")

log(
    "Q24_POINTED_MAP",
    alpha=alpha24,
    raw_match=1,
    scaling=1,
    status="PASS",
)


# ---------------------------------------------------------------------------
# q8 branch-anchor map helpers.
# ---------------------------------------------------------------------------
Tii = QQ(anchor["zero"]["old_base_T"])
q8_shift = anchor["quartic_to_anchor"]["shifted_coefficients"]
q8_urst = anchor["anchor_to_canonical"]["urst"]

def eval_U_text(text, Umap):
    val = KU(sage_eval(str(text), locals={"U": U}))
    return eval_rf(val, Umap, RU)


# q8 quartic string can be evaluated directly at U=Umap, T=current q6 base.
branch_quartic_text = str(q8["pencil"]["branch_quartic"])


# ---------------------------------------------------------------------------
# q24 base/pencil evaluation on an exact D13 point.
# ---------------------------------------------------------------------------
def q24_base_map(Umap, x13, y13):
    z = eval_poly(Z24, Umap)
    xx = eval_poly(X24, Umap)
    yy = eval_poly(Y24, Umap)
    xh = xx / z**2
    yh = yy / z**3

    if x13 == xh:
        raise ArithmeticError("D13 curve hits q24 marked section identically")

    m = (y13 + yh) / (x13 - xh)

    values = []
    for AA, BB in q24_pairs:
        aa = eval_poly(AA, Umap) / z**2
        bb = eval_poly(BB, Umap) / z
        values.append(aa + bb*m)

    if not values[0]:
        raise ArithmeticError("q24 first pencil section vanished identically")

    vmap = normalize_rf(values[1] / values[0], RT)
    n = RT(vmap.numerator())
    d = RT(vmap.denominator())
    degree = max(n.degree(), d.degree())
    return vmap, degree


def eval_KV_at_T(value, Vmap):
    return eval_rf(KV(value), Vmap, RV)


def q24_quartic_value(Umap, Vmap):
    coeff = [eval_KV_at_T(q, Vmap) for q in q24_coeff]
    return sum(coeff[i] * Umap**i for i in range(5))


def q24_point_to_minimal(Umap, Vmap, Wvalue):
    """
    Exact pointed quartic -> raw D12 -> certified minimal D12.
    """
    uu = Umap - KT(alpha24)
    if not uu:
        raise ArithmeticError("curve is the q24 pointing section itself")

    qt = eval_KV_at_T(q0, Vmap)
    at = eval_KV_at_T(a24, Vmap)
    bt = eval_KV_at_T(b24, Vmap)
    ct = eval_KV_at_T(c24, Vmap)
    dt = eval_KV_at_T(d24, Vmap)

    xg = (2*qt*(Wvalue + qt) + dt*uu) / uu**2
    yg = (
        4*qt**2*(Wvalue + qt)
        + 2*qt*(dt*uu + ct*uu**2)
        - dt**2*uu**2/(2*qt)
    ) / uu**3

    a1 = dt/qt
    a2 = ct - dt**2/(4*qt**2)
    a3 = 2*qt*bt
    a4 = -4*qt**2*at
    a6 = a2*a4
    b2 = a1**2 + 4*a2

    xs = xg + b2/12
    ys = yg + (a1*xg + a3)/2

    xraw = 9*xs
    yraw = 27*ys

    rawAt = eval_KV_at_T(rawA, Vmap)
    rawBt = eval_KV_at_T(rawB, Vmap)
    if yraw**2 != xraw**3 + rawAt*xraw + rawBt:
        raise ArithmeticError("q24 pointed map misses raw D12 Jacobian")

    s2 = eval_KV_at_T(scale2, Vmap)
    s3 = eval_KV_at_T(scale3, Vmap)
    xmin = normalize_rf(xraw/s2, RT)
    ymin = normalize_rf(yraw/s3, RT)

    Amin = eval_poly(A12, Vmap)
    Bmin = eval_poly(B12, Vmap)
    if ymin**2 != xmin**3 + Amin*xmin + Bmin:
        # harmless global y-scaling sign
        ymin = -ymin
        if ymin**2 != xmin**3 + Amin*xmin + Bmin:
            raise ArithmeticError("raw/minimal D12 point conversion failed")

    return xmin, ymin


def mobius_inverse(vmap):
    vmap = normalize_rf(vmap, RT)
    n = RT(vmap.numerator())
    d = RT(vmap.denominator())
    if max(n.degree(), d.degree()) != 1:
        raise ArithmeticError("base map is not Mobius")

    n0 = QQ(n[0])
    n1 = QQ(n[1] if n.degree() >= 1 else 0)
    d0 = QQ(d[0])
    d1 = QQ(d[1] if d.degree() >= 1 else 0)

    # V*(d0+d1*T)=n0+n1*T.
    num = RV(n0 - V*d0)
    den = RV(V*d1 - n1)
    if not den:
        raise ArithmeticError("degenerate Mobius inverse")
    return KV(num) / KV(den)


def globalize_section_from_T(xmin, ymin, vmap):
    """
    First pass through the PGL-transformed polynomial K3 gauge to expose and
    cancel the low-degree section, then map the low-degree section back to V.
    """
    vmap = normalize_rf(vmap, RT)
    n = RT(vmap.numerator())
    d = RT(vmap.denominator())

    At = homogeneous_compose(A12, n, d, 8)
    Bt = homogeneous_compose(B12, n, d, 12)

    xt = normalize_rf(KT(d**4) * xmin, RT)
    yt = normalize_rf(KT(d**6) * ymin, RT)

    XT, YT, ZT = projectivize_section(xt, yt, At, Bt, RT)

    # Map the already-cancelled low-degree section back to the original
    # q24 base coordinate.
    TinV = mobius_inverse(vmap)
    dV = eval_poly(d, TinV)

    xV = (
        eval_poly(XT, TinV) / eval_poly(ZT, TinV)**2
    ) / dV**4
    yV = (
        eval_poly(YT, TinV) / eval_poly(ZT, TinV)**3
    ) / dV**6

    xV = normalize_rf(KV(xV), RV)
    yV = normalize_rf(KV(yV), RV)

    XV, YV, ZV = projectivize_section(
        xV, yV, A12, B12, RV
    )
    return xV, yV, XV, YV, ZV


# ---------------------------------------------------------------------------
# Transport O12/P42.  q8 and q24 square-root signs are deliberately retained
# as small exact branches; no expensive global search.
# ---------------------------------------------------------------------------
curve_candidates = {"O12": [], "P42": []}

for label in ("O12", "P42"):
    rec = pre["targets"][label]
    Umap = parse_record(rec["q8_parameter"], RT)
    expected_q8_degree = int(rec["q8_parameter_degree"])

    q8_eval = KT(
        sage_eval(
            branch_quartic_text,
            locals={"U": Umap, "T": T},
        )
    )
    W8root = fast_rational_square_root(q8_eval, RT)
    if W8root is None:
        raise ArithmeticError(f"{label}: q8 quartic restriction is not square")

    for sign8 in (1, -1):
        W8 = sign8 * W8root
        r = KT(T - Tii)

        d8 = eval_U_text(q8_shift["d_r1"], Umap)
        uiso = eval_U_text(q8_urst[0], Umap)
        riso = eval_U_text(q8_urst[1], Umap)
        siso = eval_U_text(q8_urst[2], Umap)
        tiso = eval_U_text(q8_urst[3], Umap)

        Xa = d8 / r
        Ya = d8 * W8 / r**2
        x13 = normalize_rf((Xa - riso) / uiso**2, RT)
        y13 = normalize_rf(
            (Ya - siso*(Xa-riso) - tiso) / uiso**3,
            RT,
        )

        A13t = eval_poly(A13, Umap)
        B13t = eval_poly(B13, Umap)
        if y13**2 != x13**3 + A13t*x13 + B13t:
            raise ArithmeticError(
                f"{label}: q8 quartic->D13 map failed for sign {sign8}"
            )

        Vmap, q24_degree = q24_base_map(Umap, x13, y13)

        print(
            "Q42FAST_MAP|"
            f"curve={label}|q8_sign={sign8}|"
            f"q8_degree={expected_q8_degree}|q24_degree={q24_degree}|"
            f"status={'KEEP' if q24_degree==1 else 'REJECT_Q24_DEGREE'}",
            flush=True,
        )

        if q24_degree != 1:
            continue

        q24_eval = q24_quartic_value(Umap, Vmap)
        W24root = fast_rational_square_root(q24_eval, RT)
        if W24root is None:
            print(
                "Q42FAST_MAP|"
                f"curve={label}|q8_sign={sign8}|"
                "reason=Q24_QUARTIC_NOT_SQUARE|status=REJECT",
                flush=True,
            )
            continue

        for sign24 in (1, -1):
            xminT, yminT = q24_point_to_minimal(
                Umap, Vmap, sign24*W24root
            )

            try:
                xV, yV, XV, YV, ZV = globalize_section_from_T(
                    xminT, yminT, Vmap
                )
            except ArithmeticError as exc:
                print(
                    "Q42FAST_MAP|"
                    f"curve={label}|q8_sign={sign8}|q24_sign={sign24}|"
                    f"reason=GLOBALIZE:{type(exc).__name__}:{exc}|status=REJECT",
                    flush=True,
                )
                continue

            key = (str(xV), str(yV))
            if any(c["key"] == key for c in curve_candidates[label]):
                continue

            item = {
                "key": key,
                "q8_sign": sign8,
                "q24_sign": sign24,
                "q8_parameter": rf_record(Umap, RT),
                "q24_parameter": rf_record(Vmap, RT),
                "x": rf_record(xV, RV),
                "y": rf_record(yV, RV),
                "X": [str(v) for v in XV.list()],
                "Y": [str(v) for v in YV.list()],
                "Z": [str(v) for v in ZV.list()],
                "degrees": [
                    int(XV.degree()),
                    int(YV.degree()),
                    int(ZV.degree()),
                ],
            }
            curve_candidates[label].append(item)

            print(
                "Q42FAST_D12_SECTION|"
                f"curve={label}|q8_sign={sign8}|q24_sign={sign24}|"
                f"Xdeg={XV.degree()}|Ydeg={YV.degree()}|Zdeg={ZV.degree()}|"
                "status=PASS_EXACT_GLOBAL_D12_SECTION",
                flush=True,
            )

    if not curve_candidates[label]:
        raise ArithmeticError(f"{label}: no exact D12 section candidate survived")


# ---------------------------------------------------------------------------
# Historical P42 = geometric P42 - geometric O12 in the D12 group law.
# Keep every small sign pairing that yields the corrected P.O=3 profile.
# ---------------------------------------------------------------------------
E12 = EllipticCurve(KV, [0, 0, 0, KV(A12), KV(B12)])

marked_candidates = []
seen_marked = set()

for oi, oc in enumerate(curve_candidates["O12"]):
    Opoint = E12(
        parse_record(oc["x"], RV),
        parse_record(oc["y"], RV),
    )

    for pi, pc in enumerate(curve_candidates["P42"]):
        Pgeom = E12(
            parse_record(pc["x"], RV),
            parse_record(pc["y"], RV),
        )

        Pmarked = Pgeom - Opoint
        if Pmarked.is_zero():
            continue
        xp, yp = Pmarked.xy()

        try:
            XP, YP, ZP = projectivize_section(
                xp, yp, A12, B12, RV, expected_z=3
            )
        except ArithmeticError:
            continue

        if XP.degree() > 10 or YP.degree() > 15:
            continue

        key = (str(normalize_rf(xp, RV)), str(normalize_rf(yp, RV)))
        if key in seen_marked:
            continue
        seen_marked.add(key)

        marked_candidates.append({
            "O12_candidate": oi,
            "P42_candidate": pi,
            "x": rf_record(xp, RV),
            "y": rf_record(yp, RV),
            "X": [str(v) for v in XP.list()],
            "Y": [str(v) for v in YP.list()],
            "Z": [str(v) for v in ZP.list()],
            "degrees": [
                int(XP.degree()),
                int(YP.degree()),
                int(ZP.degree()),
            ],
        })

        print(
            "Q42FAST_MARKED|"
            f"O12={oi}|P42={pi}|"
            f"Xdeg={XP.degree()}|Ydeg={YP.degree()}|Zdeg={ZP.degree()}|"
            "PdotO=3|status=PASS_CORRECTED_PROFILE",
            flush=True,
        )

if not marked_candidates:
    raise ArithmeticError(
        "no O12/P42 sign pairing produced corrected P.O=3 D12 section"
    )


# ---------------------------------------------------------------------------
# Tiny exact chord compiler.
# ---------------------------------------------------------------------------
core = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), core)
binary_quartic_invariants = core["binary_quartic_invariants"]
kodaira_data_from_short_orders = core["kodaira_data_from_short_orders"]

RM = PolynomialRing(QQ, "M")
M = RM.gen()
KM = RM.fraction_field()

RVM = PolynomialRing(KM, "V")
Vo = RVM.gen()


def lift_V_poly(poly):
    poly = RV(poly)
    return RVM([KM(c) for c in poly.list()])


def exact_quo(poly, divisor, label):
    q, r = poly.quo_rem(divisor)
    if r:
        raise ArithmeticError(f"{label}: nonzero quotient remainder")
    return q


def tiny_yun_squareclass(F):
    F = RVM(F)
    lc = KM(F.leading_coefficient())
    f = RVM(F / lc)
    if not f.is_monic():
        raise ArithmeticError("tiny Yun monic normalization failed")

    g = f.gcd(f.derivative())
    w = exact_quo(f, g, "yun-initial")

    odd = RVM.one()
    square = RVM.one()
    multiplicities = []
    mult = 1

    while w.degree() > 0:
        y = w.gcd(g)
        z = exact_quo(w, y, f"yun-z-{mult}")
        if z.degree() > 0:
            multiplicities.append(mult)
            if mult % 2:
                odd *= z
            if mult // 2:
                square *= z**(mult // 2)
        w = y
        g = exact_quo(g, y, f"yun-g-{mult}")
        mult += 1

    if g.degree() > 0:
        raise ArithmeticError("tiny Yun left residual repeated factor")

    quartic = RVM(lc * odd)
    if F != RVM(quartic * square**2):
        raise ArithmeticError("tiny Yun exact reconstruction failed")

    return quartic, square, multiplicities


def squarefree_part(poly):
    poly = RM(poly)
    if poly.degree() <= 0:
        return RM.one()
    g = poly.gcd(poly.derivative())
    q, r = poly.quo_rem(g)
    if r:
        raise ArithmeticError("squarefree support quotient failed")
    return RM(q)


def poly_order(poly, factor):
    poly = RM(poly)
    if not poly:
        raise ArithmeticError("valuation of zero polynomial")
    count = 0
    while poly.degree() >= factor.degree():
        q, r = poly.quo_rem(factor)
        if r:
            break
        count += 1
        poly = RM(q)
    return count


def rf_order(value, factor):
    value = KM(value)
    return (
        poly_order(value.numerator(), factor)
        - poly_order(value.denominator(), factor)
    )


def fast_classify(coefficient_a, coefficient_b):
    coefficient_a = KM(coefficient_a)
    coefficient_b = KM(coefficient_b)
    delta = KM(-16) * (
        KM(4)*coefficient_a**3 + KM(27)*coefficient_b**2
    )
    if not delta:
        raise ArithmeticError("child Jacobian discriminant is zero")

    support = RM.one()
    for p in (
        RM(coefficient_a.denominator()),
        RM(coefficient_b.denominator()),
        RM(delta.numerator()),
        RM(delta.denominator()),
    ):
        sf = squarefree_part(p)
        common = support.gcd(sf)
        q, r = sf.quo_rem(common)
        if r:
            raise ArithmeticError("support union failed")
        support *= q

    if support.degree() > 0:
        support /= support.leading_coefficient()

    factors = tuple(f for f, unused in support.factor()) if support.degree() > 0 else ()

    scaling_unit = KM.one()
    raw_places = []
    for factor in sorted(factors, key=str):
        va = rf_order(coefficient_a, factor)
        vb = rf_order(coefficient_b, factor)
        vd = rf_order(delta, factor)
        scale = min(va // 4, vb // 6)
        scaling_unit *= KM(factor)**(-scale)
        raw_places.append((factor, va, vb, vd, scale))

    minA_rf = coefficient_a * scaling_unit**4
    minB_rf = coefficient_b * scaling_unit**6
    minD_rf = delta * scaling_unit**12

    if any(v.denominator() != 1 for v in (minA_rf, minB_rf, minD_rf)):
        raise ArithmeticError("finite minimization left denominators")

    minA = RM(minA_rf.numerator())
    minB = RM(minB_rf.numerator())
    minD = RM(minD_rf.numerator())

    root_rank = 0
    root_det = ZZ(1)
    euler = 0
    finite = []

    for factor, unused_va, unused_vb, unused_vd, unused_scale in raw_places:
        orders = (
            poly_order(minA, factor),
            poly_order(minB, factor),
            poly_order(minD, factor),
        )
        if orders[2] == 0:
            continue
        rr, ee, dd, symbol = kodaira_data_from_short_orders(*orders)
        degree = int(factor.degree())
        root_rank += degree * int(rr)
        euler += degree * int(ee)
        root_det *= ZZ(dd)**degree
        finite.append({
            "factor": str(factor),
            "degree": degree,
            "orders": list(map(int, orders)),
            "kodaira": symbol,
        })

    infinity_raw = (
        -minA.degree(),
        -minB.degree(),
        -minD.degree(),
    )
    infinity_scale = min(infinity_raw[0] // 4, infinity_raw[1] // 6)
    infinity_orders = tuple(
        infinity_raw[i] - (4, 6, 12)[i]*infinity_scale
        for i in range(3)
    )
    infinity_kind = "smooth"
    if infinity_orders[2] > 0:
        rr, ee, dd, infinity_kind = kodaira_data_from_short_orders(
            *infinity_orders
        )
        root_rank += int(rr)
        euler += int(ee)
        root_det *= ZZ(dd)

    return {
        "minimal_A": minA,
        "minimal_B": minB,
        "minimal_Delta": minD,
        "finite_fibres": finite,
        "infinity_orders": list(map(int, infinity_orders)),
        "infinity_kind": infinity_kind,
        "root_rank": int(root_rank),
        "root_det": int(root_det),
        "euler": int(euler),
    }


passing_children = []

for ci, cand in enumerate(marked_candidates):
    XP = RV([QQ(v) for v in cand["X"]])
    YP = RV([QQ(v) for v in cand["Y"]])
    ZP = RV([QQ(v) for v in cand["Z"]])

    Zm = lift_V_poly(ZP)
    Xm = lift_V_poly(XP)
    Ym = lift_V_poly(YP)
    Am = lift_V_poly(A12)
    mm = KM(M)

    # Marked quartic point is -P, because the chord function is
    #     m=(y+yP)/(x-xP).
    Fdisc = RVM(
        mm**4 * Zm**4
        - 6*Xm*mm**2*Zm**2
        - 8*Ym*mm*Zm
        - 3*Xm**2
        - 4*Am*Zm**4
    )

    quartic, square_factor, multiplicities = tiny_yun_squareclass(Fdisc)
    qdeg = int(quartic.degree())
    if qdeg not in (3, 4):
        print(
            "Q42FAST_CHILD|"
            f"candidate={ci}|quartic={qdeg}|status=REJECT_QUARTIC_DEGREE",
            flush=True,
        )
        continue

    I, J = binary_quartic_invariants(quartic)
    jacA = KM(-27) * KM(I)
    jacB = KM(-27) * KM(J)

    classification = fast_classify(jacA, jacB)

    is_a11 = (
        classification["root_rank"] == 11
        and classification["root_det"] == 12
        and classification["euler"] == 24
    )

    print(
        "Q42FAST_CHILD|"
        f"candidate={ci}|quartic={qdeg}|"
        f"root_rank={classification['root_rank']}|"
        f"root_det={classification['root_det']}|"
        f"euler={classification['euler']}|MW=6|"
        f"status={'PASS_A11' if is_a11 else 'REJECT_CHILD'}",
        flush=True,
    )

    child_payload = {
        **cand,
        "candidate_index": ci,
        "quartic": {
            "degree": qdeg,
            "coefficients_in_V_low_to_high": [
                str(v) for v in quartic.list()
            ],
            "yun_multiplicities": list(map(int, multiplicities)),
            "I": str(I),
            "J": str(J),
        },
        "jacobian_raw": {
            "A": str(jacA),
            "B": str(jacB),
        },
        "child": {
            "minimal_A_coefficients_low_to_high": [
                str(v) for v in classification["minimal_A"].list()
            ],
            "minimal_B_coefficients_low_to_high": [
                str(v) for v in classification["minimal_B"].list()
            ],
            "minimal_discriminant_coefficients_low_to_high": [
                str(v) for v in classification["minimal_Delta"].list()
            ],
            "finite_fibres": classification["finite_fibres"],
            "infinity_orders": classification["infinity_orders"],
            "infinity_kind": classification["infinity_kind"],
            "root_rank": classification["root_rank"],
            "root_det": classification["root_det"],
            "euler": classification["euler"],
            "MW_rank_if_rho19": 6,
        },
        "status": "PASS_A11" if is_a11 else "REJECT_CHILD",
    }

    if is_a11:
        passing_children.append(child_payload)


if not passing_children:
    raise ArithmeticError(
        "no corrected P.O=3 sign pairing compiled to exact A11/MW6"
    )

payload = {
    "schema": "elkies-k3.h3-q24-d12-orbit42-fast-transport-qq.v1",
    "status": "PASS_Q42_FAST_A11_CANDIDATES",
    "inputs": {
        "q6_preflight": str(PREFLIGHT.relative_to(ROOT)),
        "q8_child": str(Q8_CHILD.relative_to(ROOT)),
        "q8_anchor": str(Q8_ANCHOR.relative_to(ROOT)),
        "q24_parent": str(Q24_PARENT.relative_to(ROOT)),
        "q24_horizontal_section": str(Q24_SECTION.relative_to(ROOT)),
    },
    "method": {
        "q6_to_q8": "EXACT_Q6_POINT_PLUS_Q8_QUARTIC",
        "q8_to_d13": "EXACT_BRANCH_ANCHOR_BIRATIONAL_MAP",
        "q24_base": "EVALUATE_CERTIFIED_2X56_PENCIL",
        "q24_to_d12": "POINTED_QUARTIC_MAP_THEN_RAW_TO_MINIMAL_SCALING",
        "marked_section": "P42_GEOMETRIC_MINUS_O12_GEOMETRIC",
        "neighbor": "DIRECT_CHORD_SLOPE",
        "quartic_squareclass": "DEGREE_LE_20_GCD_ONLY_YUN",
        "rr_search_avoided": True,
        "primary_decomposition_avoided": True,
        "hensel_avoided": True,
    },
    "geometric_sections": {
        "O12_candidates": curve_candidates["O12"],
        "P42_candidates": curve_candidates["P42"],
    },
    "corrected_marked_candidates": marked_candidates,
    "passing_A11_candidates": passing_children,
    "next": (
        "Run the two physical I8* orientation certifiers in parallel. "
        "They only pull back the single chord slope and compare its exact "
        "exceptional valuations with the two C10/C11 spinor orientations."
    ),
}

OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42FAST_RESULT|"
    f"marked={len(marked_candidates)}|A11={len(passing_children)}|"
    "status=PASS_Q42_FAST_A11_CANDIDATES",
    flush=True,
)
