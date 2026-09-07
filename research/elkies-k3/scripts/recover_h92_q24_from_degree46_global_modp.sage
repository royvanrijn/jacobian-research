#!/usr/bin/env sage -python
"""
Recover the H92 q24 horizontal section over GF(p)(U) from the exact
degree-46 bridge

    W = Qmap - S3

using the globally coherent q8 orientation.

Previously certified at U=2:
    deg_q8(W) = 46
    q24 = AJ(W) + 2*G1
on the canonical D13 child.

This script globalizes that construction:

  1. construct W exactly over QQ(T);
  2. use q8-global-orientation.json to specialize one coherent quartic,
     square factor, and D13 c^2/c^3 scaling at every U=tau;
  3. trace the 46 covariant images with L(47 O);
  4. subtract 46 times the globally resolved IV* origin and halve, giving AJ(W);
  5. add the exact canonical point 2*G1;
  6. collect globally coherent q24 samples;
  7. reconstruct
         x = X/Z^2,   deg X=52, deg Z=24,
         y = Y/Z^3,   deg Y=78;
  8. verify the D13 Weierstrass equation identically over GF(p)[U];
  9. compare the full recovered section with the independent modular q24
     from q8-q24-canonical-backtrack.json.

No Hensel lift and no D13 lattice-coordinate transport is used in the recovery.
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import (
    EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector
)


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (
            (candidate / "elkies-k3/scripts").is_dir()
            and (candidate / "artifacts/generated-results").is_dir()
        ):
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo", type=Path)
parser.add_argument("--prime", type=int, default=100003)
parser.add_argument("--start", type=int, default=2)
parser.add_argument("--samples", type=int, default=115)
parser.add_argument("--scan-limit", type=int, default=1000)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if args.samples < 105:
    raise ValueError("need at least 105 good samples for a certified 52/48 interpolation")

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"

Q6 = GEN / "elkies-k3-h92-q6-child-jacobian.json"
ZERO = GEN / "elkies-k3-h92-q6-child-zero-section.json"
COMP = GEN / "elkies-k3-h92-q6-child-e7-infinity-sections.json"
Q8_CANDIDATES = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8 = next(
    (
        p for p in Q8_CANDIDATES
        if p.exists()
        and "rr" in json.loads(p.read_text())
        and "kernel_polynomials" in json.loads(p.read_text()).get("rr", {})
        and "child" in json.loads(p.read_text())
    ),
    None,
)
S3BR = LOCAL / "q6-third-to-q8-bridge.json"
ORIENT = LOCAL / "q8-global-orientation.json"
G3ART = LOCAL / "q8-d13-g3-from-e77-bisection.json"
BACK = LOCAL / "q8-q24-canonical-backtrack.json"

if Q8 is None:
    raise SystemExit("No complete corrected q8 D13 child artifact found")

for path in (CORE, Q6, ZERO, COMP, Q8, S3BR, ORIENT, G3ART, BACK):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q24-degree46-global-mod-{args.prime}.json"
)

scope = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), scope)
binary_quartic_invariants = scope["binary_quartic_invariants"]

q6 = json.loads(Q6.read_text())
zero = json.loads(ZERO.read_text())
comp = json.loads(COMP.read_text())
q8 = json.loads(Q8.read_text())
s3br = json.loads(S3BR.read_text())
orient = json.loads(ORIENT.read_text())
g3art = json.loads(G3ART.read_text())
back = json.loads(BACK.read_text())

assert q6["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert zero["status"] == "PASS_EXACT_CHILD_ZERO_SECTION_TRANSPORT"
assert comp["status"] == "PASS_EXACT_CHILD_E7_INFINITY_TRANSPORT"
assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert s3br["status"] == "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"
assert orient["status"] == "PASS_EXACT_Q8_GLOBAL_ORIENTATION"
assert orient["schema"] == "elkies-k3.h92-q8-global-orientation.v2"
assert g3art["status"] == "PASS_EXACT_D13_G3_FROM_E77_BISECTION"
assert back["status"] == "PASS_EXPLICIT_Q24_MODP_FROM_AJ_G1_G3"

# Resolve the globally coherent IV* origin from THIS degree-46 construction.
# At the first usable specialization, choose the unique global origin for which
# AJ(Qmap-S3) + 2*G1 equals the independently known q24 point.
resolved_origin = None
# The globally oriented covariant model can differ by a single global
# elliptic involution from the anchored branch-zero convention used by the
# direct L(47O) trace. Resolve this sign once together with the IV* origin.
resolved_aj_sign = None

# ===========================================================================
# Exact QQ(T) bridge W = Qmap - S3 and repaired q8 map U(T).
# ===========================================================================

QT = PolynomialRing(QQ, "T")
Tq = QT.gen()
QKT = QT.fraction_field()

def qpoly(values):
    return QT([QQ(v) for v in values])

def qrat(data, np, dp):
    return QKT(qpoly(data[np])) / QKT(qpoly(data[dp]))

A6q = qpoly(q6["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B6q = qpoly(q6["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
Eq6 = EllipticCurve(QKT, [0,0,0,QKT(A6q),QKT(B6q)])

zdata = zero["section"]
old_zero = Eq6(
    qrat(
        zdata,
        "x_numerator_coefficients_low_to_high",
        "x_denominator_coefficients_low_to_high",
    ),
    qrat(
        zdata,
        "y_numerator_coefficients_low_to_high",
        "y_denominator_coefficients_low_to_high",
    ),
)

entries = {entry["sign"]: entry for entry in comp["sections"]}
points = {
    sign: Eq6(
        qrat(
            entry,
            "x_numerator_coefficients_low_to_high",
            "x_denominator_coefficients_low_to_high",
        ),
        qrat(
            entry,
            "y_numerator_coefficients_low_to_high",
            "y_denominator_coefficients_low_to_high",
        ),
    )
    for sign, entry in entries.items()
}
affine = points[comp["source"]["affine_E7_sign"]]
e77 = points[comp["source"]["E7_7_sign"]]
Pmap = e77 - old_zero
Qmap = e77 - affine

s3data = s3br["third_section_canonical_q6"]
S3 = Eq6(
    qrat(
        s3data["x"],
        "numerator_coefficients_low_to_high",
        "denominator_coefficients_low_to_high",
    ),
    qrat(
        s3data["y"],
        "numerator_coefficients_low_to_high",
        "denominator_coefficients_low_to_high",
    ),
)

W = Qmap - S3
assert W in Eq6 and not W.is_zero()
wxq, wyq = W.xy()
assert wyq**2 == wxq**3 + QKT(A6q)*wxq + QKT(B6q)

mdata = q8["marking"]["section"]
sxq = qrat(
    mdata,
    "x_numerator_coefficients_low_to_high",
    "x_denominator_coefficients_low_to_high",
)
syq = qrat(
    mdata,
    "y_numerator_coefficients_low_to_high",
    "y_denominator_coefficients_low_to_high",
)
Smark = Eq6(sxq, syq)
assert Smark == Pmap + Qmap

def monic_power_root(value, exponent):
    out = value.parent().one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        out *= factor.monic()**(multiplicity//exponent)
    return out.monic()

nxq, dxq = QT(sxq.numerator()), QT(sxq.denominator())
nyq, dyq = QT(syq.numerator()), QT(syq.denominator())
hq = monic_power_root(dxq, 2)
assert hq == monic_power_root(dyq, 3)
iiq = QT(next(
    item for item in q6["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).monic()
ivq = QT(next(
    item for item in q6["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).monic()
tivq = -ivq[0]/ivq[1]
Mq = (iiq**2 * ivq**2).monic()

normalizerq = (nyq*dxq*(hq*dyq).inverse_mod(nxq)).mod(nxq)
p_fun_q = -syq/sxq
rhoq = (normalizerq*nxq.inverse_mod(Mq)).mod(Mq)

qpairs = []
for entry in q8["rr"]["kernel_polynomials"]:
    sp = QT(entry["s"])
    tp = QT(entry["t"])
    Bcoef = QKT(sp)/QKT(hq)
    Acoef = (
        -QKT(sp)*p_fun_q/QKT(hq)
        - QKT(sp)*QKT(normalizerq)/QKT(nxq)
        + QKT(sp*rhoq)
        + QKT(tp*Mq)
    )
    qpairs.append((Acoef, Bcoef))
(qA0, qB0), (qA1, qB1) = qpairs

mWq = (wyq + syq)/(wxq - sxq)
UWq = QKT((qA1 + qB1*mWq)/(qA0 + qB0*mWq))
Uwnq = QT(UWq.numerator())
Uwdq = QT(UWq.denominator())
assert Uwnq.gcd(Uwdq) in QQ
assert max(Uwnq.degree(), Uwdq.degree()) == 46

print(
    "Q24D46GLOBAL_BRIDGE|formula=Qmap-S3|"
    "old_mw=-2,-1,-1|standard_mw=0,-2,-1|"
    f"q8={Uwnq.degree()}/{Uwdq.degree()}|degree=46|status=PASS",
    flush=True,
)

# ===========================================================================
# Load exact generic q8 orientation over QQ(U).
# ===========================================================================

QU = PolynomialRing(QQ, "U")
Uq = QU.gen()
QKU = QU.fraction_field()
QTU = PolynomialRing(QKU, "T")
TT = QTU.gen()
QKTU = QTU.fraction_field()

def ku_from_record(record):
    return QKU(
        QU([QQ(v) for v in record["numerator_coefficients_low_to_high"]])
    ) / QKU(
        QU([QQ(v) for v in record["denominator_coefficients_low_to_high"]])
    )

def tu_from_record(record):
    return QTU([
        ku_from_record(value)
        for value in record["coefficients_low_to_high"]
    ])

def ktu_from_record(record):
    return QKTU(tu_from_record(record["numerator"])) / QKTU(
        tu_from_record(record["denominator"])
    )

transport = orient["generic_transport"]
quartic_generic = tu_from_record(transport["quartic_in_T_over_QQ_U"])
square_factor_generic = ktu_from_record(
    transport["square_factor_in_QQ_U_T"]
)
wiv_generic = ku_from_record(orient["quartic"]["global_w_plus"])
c2_generic = ku_from_record(
    orient["minimalization"]["c2_equals_c_squared"]
)
c3_generic = ku_from_record(
    orient["minimalization"]["c3_equals_c_cubed"]
)
assert quartic_generic.degree() == 4

A13q = QU([
    QQ(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]
])
B13q = QU([
    QQ(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]
])

I_generic, J_generic = binary_quartic_invariants(quartic_generic)
assert wiv_generic**2 == QKU(quartic_generic(tivq))
assert c2_generic**2 == QKU(A13q)/QKU(-27*I_generic)
assert c3_generic**2 == QKU(B13q)/QKU(-27*J_generic)

print(
    "Q24D46GLOBAL_ORIENT|origin=UNRESOLVED|"
    "quartic=PASS|c2c3=PASS|"
    "resolver=degree46_AJ_plus_2G1_vs_independent_q24|status=PASS",
    flush=True,
)

# ===========================================================================
# Finite-field helpers.
# ===========================================================================

p = ZZ(args.prime)
if not p.is_prime() or p in (2,3):
    raise ValueError("prime must be odd and !=3")
F = GF(p)
RT = PolynomialRing(F, "T")
T = RT.gen()
KT = RT.fraction_field()
RU = PolynomialRing(F, "U")
U = RU.gen()

def modq(value):
    value = QQ(value)
    den = ZZ(value.denominator())
    if den % p == 0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {value}")
    return F(ZZ(value.numerator())) / F(den)

def reduce_qpoly(poly):
    poly = QT(poly)
    return RT([modq(v) for v in poly.list()])

def reduce_qrat(value):
    value = QKT(value)
    return KT(reduce_qpoly(value.numerator())) / KT(
        reduce_qpoly(value.denominator())
    )

def reduce_qupoly(poly):
    poly = QU(poly)
    return RU([modq(v) for v in poly.list()])

def specialize_KU(value, tau):
    value = QKU(value)
    num = reduce_qupoly(value.numerator())(tau)
    den = reduce_qupoly(value.denominator())(tau)
    if not den:
        raise ZeroDivisionError("QQ(U) denominator vanished")
    return num/den

def specialize_TU(poly, tau):
    poly = QTU(poly)
    return RT([specialize_KU(c, tau) for c in poly.list()])

def specialize_KTU(value, tau):
    value = QKTU(value)
    num = specialize_TU(value.numerator(), tau)
    den = specialize_TU(value.denominator(), tau)
    if not den:
        raise ZeroDivisionError("QQ(U)(T) denominator vanished identically")
    return KT(num)/KT(den)

A6 = reduce_qpoly(A6q)
B6 = reduce_qpoly(B6q)
sx = reduce_qrat(sxq)
sy = reduce_qrat(syq)
wx = reduce_qrat(wxq)
wy = reduce_qrat(wyq)
(A0, B0), (A1, B1) = [
    (reduce_qrat(a), reduce_qrat(b)) for a,b in qpairs
]
Uwn = reduce_qpoly(Uwnq)
Uwd = reduce_qpoly(Uwdq)
tiv = modq(tivq)
A13U = reduce_qupoly(A13q)
B13U = reduce_qupoly(B13q)

R2 = PolynomialRing(F, "r2")
r2 = R2.gen()

def native_covariants(quartic):
    base = quartic.base_ring()
    BR = PolynomialRing(base, names=("qx", "qz"))
    qx, qz = BR.gens()
    f = sum(BR(quartic[i])*qx**i*qz**(4-i) for i in range(5))
    H = (
        f.derivative(qx,2)*f.derivative(qz,2)
        - f.derivative(qx).derivative(qz)**2
    ) / base(3)
    G = f.derivative(qx)*H.derivative(qz) - f.derivative(qz)*H.derivative(qx)
    I, J = binary_quartic_invariants(quartic)
    I, J = base(I), base(J)
    assert G**2 == (
        -base(16)/base(3)*H**3
        + base(256)*I*H*f**2
        - base(1024)/base(3)*J*f**3
    )
    return f, H, G

def newton_power_sums(poly):
    n = poly.degree()
    assert poly[n] == 1
    sums = [F(n)]
    for k in range(1,n):
        total = F(k)*poly[n-k]
        for j in range(1,k):
            total += poly[n-j]*sums[k-j]
        sums.append(-total)
    return sums

def rational_halves(E, target):
    if target.is_zero():
        return []
    tx, unused_ty = target.xy()
    a, b = E.a4(), E.a6()
    HR = PolynomialRing(F, "r")
    r = HR.gen()
    polynomial = (
        r**4 - 4*tx*r**3 - 2*a*r**2
        - (4*a*tx + 8*b)*r + a**2 - 4*b*tx
    )
    halves = []
    for factor, multiplicity in polynomial.factor():
        if multiplicity != 1 or factor.degree() != 1:
            continue
        xh = -factor[0]/factor[1]
        rhs = xh**3 + a*xh + b
        if not rhs.is_square():
            continue
        for yh in rhs.sqrt(all=True):
            P = E(xh, yh)
            if 2*P == target and all(P != old for old in halves):
                halves.append(P)
    return halves

# Exact canonical G1 over QQ(U), from the already-passing bisection artifact.
G1_exact = g3art["canonical_D13"]["G1"]

def parse_exact_KU(text):
    return QKU(str(text))

def specialize_exact_point(record, tau, E):
    if record.get("zero"):
        return E(0)
    x = specialize_KU(parse_exact_KU(record["x"]), tau)
    y = specialize_KU(parse_exact_KU(record["y"]), tau)
    return E(x,y)

# Independent modular q24, used only for the U=2 global closure check and
# final whole-function equality.
qref = back["q24_modp"]

def reference_q24_at(tau, E):
    def ev(vals):
        return sum(F(int(v))*tau**i for i,v in enumerate(vals))
    xd = ev(qref["x_denominator_coefficients_low_to_high"])
    yd = ev(qref["y_denominator_coefficients_low_to_high"])
    if not xd or not yd:
        raise ZeroDivisionError("independent q24 reference pole")
    return E(
        ev(qref["x_numerator_coefficients_low_to_high"])/xd,
        ev(qref["y_numerator_coefficients_low_to_high"])/yd,
    )

# ===========================================================================
# One globally coherent degree-46 q24 sample.
# ===========================================================================

def sample_tau(tau):
    a13 = A13U(tau)
    b13 = B13U(tau)
    E13 = EllipticCurve(F, [0,0,0,a13,b13])
    if not E13.discriminant():
        return None, "singular_D13"
    if (r2**3 + a13*r2 + b13).roots():
        return None, "rational_2_torsion"

    H46 = RT(Uwn - tau*Uwd)
    if H46.degree() != 46:
        return None, f"degree_drop_{H46.degree()}"
    H46 = H46.monic()
    if H46.gcd(H46.derivative()).degree() != 0:
        return None, "non_etale"

    try:
        quartic = specialize_TU(quartic_generic, tau)
        square_factor = specialize_KTU(square_factor_generic, tau)
        q8_m = -(A1 - tau*A0)/(B1 - tau*B0)
        wiv = specialize_KU(wiv_generic, tau)
        c2 = specialize_KU(c2_generic, tau)
        c3 = specialize_KU(c3_generic, tau)
    except ZeroDivisionError:
        return None, "global_orientation_pole"

    if quartic.degree() != 4:
        return None, f"quartic_degree_{quartic.degree()}"
    if wiv*wiv != quartic(tiv):
        return None, "global_wIV_mismatch"

    fbin, HC, GC = native_covariants(quartic)
    qx, qz = HC.parent().gens()
    H_u = RT(HC(qx=T, qz=F(1)))
    G_u = RT(GC(qx=T, qz=F(1)))

    # Check the stored minimal D13 scaling.
    I, J = binary_quartic_invariants(quartic)
    if (-F(27)*F(I))*c2**2 != a13:
        return None, "global_c2_mismatch"
    if (-F(27)*F(J))*c3**2 != b13:
        return None, "global_c3_mismatch"

    def reduce_mod_H(value):
        value = KT(value)
        num = RT(value.numerator())
        den = RT(value.denominator())
        if den.gcd(H46).degree() != 0:
            raise ZeroDivisionError
        return (num*den.inverse_mod(H46)) % H46

    try:
        mW = (wy + sy)/(wx - sx)
        if reduce_mod_H(q8_m-mW):
            return None, "chord_mismatch"
        wW = (2*wx + sx - q8_m**2)/square_factor
        wA = reduce_mod_H(wW)
    except ZeroDivisionError:
        return None, "noninvertible_bridge_denominator"

    if (wA*wA - quartic) % H46:
        return None, "quartic_sqrt_mismatch"

    quartic_inv = quartic.inverse_mod(H46)
    xA = (-F(3)/F(4) * H_u * quartic_inv * c2) % H46
    yA = (F(9)/F(32) * G_u * wA * quartic_inv**2 * c3) % H46
    if (yA*yA - xA*xA*xA - a13*xA - b13) % H46:
        return None, "etale_child_miss"

    # L(47 O): 1,x,...,x^23, y,xy,...,x^22 y.
    xpowers = [RT.one()]
    for unused in range(23):
        xpowers.append((xpowers[-1]*xA) % H46)
    columns = list(xpowers)
    columns += [(yA*xpowers[e]) % H46 for e in range(23)]
    assert len(columns) == 47

    Eval = matrix(F, 46, 47, lambda row,col: columns[col][row])
    ker = Eval.right_kernel().basis_matrix()
    if ker.nrows() != 1:
        return None, f"L47_kernel_{ker.nrows()}"
    rel = ker[0]

    XR = PolynomialRing(F, "X")
    Xv = XR.gen()
    Afun = sum(rel[i]*Xv**i for i in range(24))
    Bfun = sum(rel[24+i]*Xv**i for i in range(23))
    Rint = Afun**2 - (Xv**3 + a13*Xv + b13)*Bfun**2
    if Rint.degree() != 47:
        return None, f"residual_degree_{Rint.degree()}"

    root_sum = -Rint[46]/Rint[47]
    ps = newton_power_sums(H46)
    trace_x = sum(xA[i]*ps[i] for i in range(46))
    xQ = root_sum - trace_x
    if not Bfun(xQ):
        return None, "trace_B_zero"
    yQ = -Afun(xQ)/Bfun(xQ)
    traceCov = -E13(xQ,yQ)

    # Two globally coherent covariant origins. Resolve their label directly
    # from this degree-46 construction instead of importing it from S3.
    fv = F(fbin(qx=tiv, qz=F(1)))
    hv = F(HC(qx=tiv, qz=F(1)))
    gv = F(GC(qx=tiv, qz=F(1)))
    if not fv:
        return None, "origin_f_zero"

    rawx = -F(3)/F(4)*hv/fv
    origins = {}
    for label, w_origin in (("plus", wiv), ("minus", -wiv)):
        if w_origin*w_origin != fv:
            return None, "origin_w_mismatch"
        rawy = F(9)/F(32)*gv*w_origin/fv**2
        origins[label] = E13(c2*rawx, c3*rawy)

    try:
        G1 = specialize_exact_point(G1_exact, tau, E13)
    except ZeroDivisionError:
        return None, "G1_pole"

    def aj_for_origin(origin):
        doubled = traceCov - 46*origin
        halves = rational_halves(E13, doubled)
        return halves[0] if len(halves) == 1 else None

    global resolved_origin, resolved_aj_sign
    if resolved_origin is None:
        try:
            Qref = reference_q24_at(tau, E13)
        except ZeroDivisionError:
            return None, "reference_pole_during_origin_resolution"

        hits = []
        for label in ("plus", "minus"):
            AJcandidate = aj_for_origin(origins[label])
            if AJcandidate is None:
                continue
            for asign in (+1, -1):
                AJoriented = asign*AJcandidate
                if AJoriented + 2*G1 == Qref:
                    hits.append((label, asign, AJoriented))

        if len(hits) != 1:
            return None, f"origin_sign_resolution_count_{len(hits)}"

        resolved_origin, resolved_aj_sign, AJ = hits[0]
        print(
            "Q24D46GLOBAL_ORIGIN|"
            f"U={int(tau)}|origin={resolved_origin}|"
            f"AJ_sign={resolved_aj_sign:+d}|"
            "criterion=AJ_oriented(Qmap-S3)+2G1==independent_q24|"
            "status=PASS_UNIQUE_GLOBAL_ORIGIN_AND_SIGN",
            flush=True,
        )
    else:
        AJraw = aj_for_origin(origins[resolved_origin])
        if AJraw is None:
            return None, "half_count_not_one"
        AJ = resolved_aj_sign*AJraw

    Q24 = AJ + 2*G1
    if Q24.is_zero():
        return None, "q24_zero"
    qx24, qy24 = Q24.xy()

    return (F(qx24), F(qy24)), None

# ===========================================================================
# Collect samples.
# ===========================================================================

samples = []
skip_counts = {}
candidate_integer = args.start
attempted = 0

while len(samples) < args.samples and attempted < args.scan_limit:
    tau = F(candidate_integer)
    candidate_integer += 1
    attempted += 1

    result, reason = sample_tau(tau)
    if result is None:
        skip_counts[reason] = skip_counts.get(reason, 0) + 1
        continue

    x0, y0 = result
    samples.append((int(tau), int(x0), int(y0)))
    count = len(samples)
    if count <= 5 or count % 10 == 0 or count == args.samples:
        print(
            f"Q24D46GLOBAL_SAMPLE|count={count}|U={int(tau)}|"
            f"x={int(x0)}|y={int(y0)}|status=PASS",
            flush=True,
        )

if len(samples) < args.samples:
    raise RuntimeError(
        f"only {len(samples)} good samples after {attempted} attempts; "
        f"skips={skip_counts}"
    )

print(
    f"Q24D46GLOBAL_SAMPLE|good={len(samples)}|attempted={attempted}|"
    f"skips={skip_counts}|stage=collection|status=PASS",
    flush=True,
)

# ===========================================================================
# Recover q24(U): x = X/D, expected degrees 52/48.
# ===========================================================================

NUM_DEG = 52
DEN_DEG = 48
ncols = (NUM_DEG+1) + (DEN_DEG+1)
IM = matrix(
    F,
    len(samples),
    ncols,
    lambda row,col: (
        F(samples[row][0])**col
        if col <= NUM_DEG
        else -F(samples[row][1])
             * F(samples[row][0])**(col-(NUM_DEG+1))
    ),
)
IK = IM.right_kernel().basis_matrix()
if IK.nrows() != 1:
    raise ArithmeticError(
        f"q24 x 52/48 interpolation kernel dimension {IK.nrows()}, expected 1"
    )

rv = IK[0]
X = RU(list(rv[:NUM_DEG+1]))
D = RU(list(rv[NUM_DEG+1:]))
if not D:
    raise ArithmeticError("q24 x denominator vanished")
scale = D.leading_coefficient()
X /= scale
D /= scale
if X.gcd(D).degree() != 0:
    raise ArithmeticError("q24 x interpolation is not reduced")
assert X.degree() == 52
assert D.degree() == 48
assert D.is_monic()

def square_root_monic(poly):
    poly = RU(poly)
    if not poly or not poly.is_monic():
        return None
    out = RU.one()
    for factor, multiplicity in poly.factor():
        if multiplicity % 2:
            return None
        out *= factor.monic()**(multiplicity//2)
    return out.monic()

Z = square_root_monic(D)
if Z is None:
    raise ArithmeticError("q24 x denominator is not a square")
assert Z.degree() == 24
assert D == Z**2

for u0,x0,unused_y in samples:
    uu = F(u0)
    if not D(uu) or X(uu)/D(uu) != F(x0):
        raise ArithmeticError("q24 x interpolation failed retained sample")

RHS = X**3 + A13U*X*Z**4 + B13U*Z**6
if not RHS.is_square():
    raise ArithmeticError("q24 Weierstrass RHS is not a polynomial square")
Y = RHS.sqrt()
assert Y.degree() == 78

direct = opposite = True
for u0,unused_x,y0 in samples:
    uu = F(u0)
    if not Z(uu):
        raise ArithmeticError("sample unexpectedly lies at q24 pole")
    pred = Y(uu)/Z(uu)**3
    direct &= pred == F(y0)
    opposite &= -pred == F(y0)
if direct == opposite:
    raise ArithmeticError("q24 Y orientation unresolved")
if opposite:
    Y = -Y

assert Y**2 == X**3 + A13U*X*Z**4 + B13U*Z**6
for u0,x0,y0 in samples:
    uu = F(u0)
    assert X(uu)/Z(uu)**2 == F(x0)
    assert Y(uu)/Z(uu)**3 == F(y0)

print(
    "Q24D46GLOBAL_INTERP|"
    f"samples={len(samples)}|x=52/48|Z=24|y=78/72|"
    "identity=PASS|samples=PASS|status=PASS_MODULAR_SECTION",
    flush=True,
)

# ===========================================================================
# Full independent cross-check against old q24 modular reconstruction.
# ===========================================================================

RF = RU.fraction_field()
x_new = RF(X)/RF(Z**2)
y_new = RF(Y)/RF(Z**3)

def old_rf(nkey,dkey):
    return RF(RU([F(int(v)) for v in qref[nkey]])) / RF(
        RU([F(int(v)) for v in qref[dkey]])
    )

x_ref = old_rf(
    "x_numerator_coefficients_low_to_high",
    "x_denominator_coefficients_low_to_high",
)
y_ref = old_rf(
    "y_numerator_coefficients_low_to_high",
    "y_denominator_coefficients_low_to_high",
)

xmatch = (x_new == x_ref)
ymatch = (y_new == y_ref)
yminus = (y_new == -y_ref)
if not xmatch or not (ymatch or yminus):
    raise ArithmeticError(
        f"full q24 cross-check failed: x={xmatch}, y={ymatch}, yminus={yminus}"
    )
if yminus:
    # The neighbour is invariant under point inversion, but align the serialized
    # section to the independent orientation for a literal coefficient match.
    Y = -Y
    y_new = -y_new
    ymatch = True
    yminus = False

assert x_new == x_ref
assert y_new == y_ref

print(
    "Q24D46GLOBAL_CROSSCHECK|"
    f"xmatch={int(x_new==x_ref)}|ymatch={int(y_new==y_ref)}|"
    "source=independent_S3_backtrack|status=PASS_IDENTICAL_SECTION",
    flush=True,
)

if resolved_origin not in ("plus", "minus"):
    raise ArithmeticError("global origin was never resolved")
if resolved_aj_sign not in (+1, -1):
    raise ArithmeticError("global covariant AJ sign was never resolved")

payload = {
    "schema": "elkies-k3.h92-q24-degree46-global-modp.v1",
    "status": "PASS_MODULAR_Q24_FROM_DEGREE46_BRIDGE",
    "prime": int(p),
    "bridge": {
        "formula": "Qmap-S3",
        "q6_old_zero_mw": [-2,-1,-1],
        "q6_standard_mw": [0,-2,-1],
        "q8_degree": 46,
        "q8_parameter": {
            "numerator_coefficients_low_to_high": [int(v) for v in Uwn.list()],
            "denominator_coefficients_low_to_high": [int(v) for v in Uwd.list()],
        },
    },
    "construction": {
        "global_orientation_artifact": str(ORIENT.relative_to(ROOT)),
        "global_orientation_sha256": hashlib.sha256(ORIENT.read_bytes()).hexdigest(),
        "resolved_origin": resolved_origin,
        "resolved_covariant_AJ_sign": int(resolved_aj_sign),
        "trace_rr_basis": "L(47O)",
        "covariant_relation": "2*AJ(W)=Trace(covariant(W))-46*origin",
        "q24_relation": "q24=AJ(W)+2*G1",
    },
    "sampling": {
        "requested": int(args.samples),
        "attempted": int(attempted),
        "good": len(samples),
        "skip_counts": skip_counts,
        "samples": [
            {"U": int(u0), "x": int(x0), "y": int(y0)}
            for u0,x0,y0 in samples
        ],
    },
    "section_mod_p": {
        "Z_coefficients_low_to_high": [int(v) for v in Z.list()],
        "X_coefficients_low_to_high": [int(v) for v in X.list()],
        "Y_coefficients_low_to_high": [int(v) for v in Y.list()],
        "profile": {
            "Z_degree": int(Z.degree()),
            "X_degree": int(X.degree()),
            "Y_degree": int(Y.degree()),
            "x_degrees": [int(X.degree()), int((Z**2).degree())],
            "y_degrees": [int(Y.degree()), int((Z**3).degree())],
        },
        "exact_weierstrass_identity": True,
    },
    "independent_crosscheck": {
        "artifact": str(BACK.relative_to(ROOT)),
        "sha256": hashlib.sha256(BACK.read_bytes()).hexdigest(),
        "x_identical": True,
        "y_identical": True,
    },
    "next": (
        "Use this independently reconstructed q24 modular section as the seed "
        "for characteristic-zero recovery, then compile the q24 D12 neighbour "
        "equation from the exact D13 line bundle."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q24D46GLOBAL_RESULT|degree=46|q24=AJ(W)+2G1|"
    "x=52/48|y=78/72|crosscheck=IDENTICAL|"
    "status=PASS_MODULAR_Q24_FROM_DEGREE46_BRIDGE",
    flush=True,
)
