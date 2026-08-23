#!/usr/bin/env sage -python
"""
Recover AJ_q8(S3) modulo p with globally coherent q8 orientation.

Prerequisites already certified:
  * exact q6 third section S3 over QQ;
  * exact q8 degree-52 bridge;
  * exact corrected q8 -> D13 child;
  * exact global q8 orientation:
        w_IV(U) in QQ(U),
        c2(U)=c(U)^2 in QQ(U),
        c3(U)=c(U)^3 in QQ(U).

This script rebuilds the generic q8 quartic and its square factor ONCE over
QQ(U).  Every finite-field specialization is obtained by reducing those exact
generic objects.  It never chooses a square root or minimizing-unit sign
independently at a specialization.

It samples both globally coherent IV* origins, traces the degree-52 S3
multisection through the q8 binary-quartic 2-cover, halves exactly on the D13
child, and tests each global branch against the pinned target profile

    D13 MW = (0,-1,1,1)
    height = 47
    local D13 correction = 1
    P.O = 22
    x degrees = 48/44
    y degrees = 72/66.

A branch is accepted only if its x samples have a unique rational
interpolation of exact degree 48/44, the denominator is Z^2 with deg Z=22,
the curve RHS is an exact square Y^2 with deg Y=72, and all sampled x/y
coordinates replay.

Run:
  sage -python ~/Downloads/recover_h92_q8_s3_aj_samples_modp_v2.sage
"""

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix


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
parser.add_argument("--samples", type=int, default=100)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if args.samples < 94:
    raise ValueError("need at least 94 samples for degree-48/44 interpolation")

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
Q6 = GEN / "elkies-k3-h92-q6-child-jacobian.json"
Q8 = GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json"
BRIDGE = LOCAL / "q6-third-to-q8-bridge.json"
ORIENT = LOCAL / "q8-global-orientation.json"
OUTPUT = (
    args.output.resolve()
    if args.output and args.output.is_absolute()
    else ROOT / (
        args.output
        if args.output
        else Path("artifacts/local/elkies-k3/q8-s3-aj-oriented-mod-100003.json")
    )
)

for path in (CORE, Q6, Q8, BRIDGE, ORIENT):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

scope = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), scope)
squarefree_binary_quartic = scope["squarefree_binary_quartic"]
binary_quartic_invariants = scope["binary_quartic_invariants"]

q6 = json.loads(Q6.read_text())
q8 = json.loads(Q8.read_text())
bridge = json.loads(BRIDGE.read_text())
orient = json.loads(ORIENT.read_text())
assert q6["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert bridge["status"] == "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"
assert orient["status"] == "PASS_EXACT_Q8_GLOBAL_ORIENTATION"

# ===========================================================================
# Exact generic q8 objects over QQ(U), including the globally normalized
# quartic and square factor.
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
assert syq**2 == sxq**3 + QKT(A6q)*sxq + QKT(B6q)


def q_monic_power_root(value, exponent):
    result = QT.one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        result *= factor.monic()**(multiplicity//exponent)
    return result.monic()


nxq, dxq = QT(sxq.numerator()), QT(sxq.denominator())
nyq, dyq = QT(syq.numerator()), QT(syq.denominator())
hq = q_monic_power_root(dxq, 2)
assert hq == q_monic_power_root(dyq, 3)
assert hq.degree() == 10

iiq = QT(next(
    item for item in q6["finite_fibres"] if item["kodaira"] == "II*"
)["factor"]).monic()
ivq = QT(next(
    item for item in q6["finite_fibres"] if item["kodaira"] == "IV*"
)["factor"]).monic()
assert iiq.degree() == ivq.degree() == 1
tivq = -ivq[0]/ivq[1]
Mq = (iiq**2 * ivq**2).monic()

normalizerq = (nyq*dxq*(hq*dyq).inverse_mod(nxq)).mod(nxq)
assert (normalizerq*hq*dyq - nyq*dxq) % nxq == 0
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
assert len(qpairs) == 2
(qA0, qB0), (qA1, qB1) = qpairs

QU = PolynomialRing(QQ, "U")
Uq = QU.gen()
QKU = QU.fraction_field()
QTU = PolynomialRing(QKU, "T")
TT = QTU.gen()
QKTU = QTU.fraction_field()


def lift_qpoly(value):
    value = QT(value)
    return QTU([QKU(c) for c in value.list()])


def lift_qrat(value):
    value = QKT(value)
    return QKTU(lift_qpoly(value.numerator())) / QKTU(
        lift_qpoly(value.denominator())
    )


m_generic = -(
    lift_qrat(qA1) - QKU(Uq)*lift_qrat(qA0)
) / (
    lift_qrat(qB1) - QKU(Uq)*lift_qrat(qB0)
)
sx_generic = lift_qrat(sxq)
sy_generic = lift_qrat(syq)
A6_generic = lift_qpoly(A6q)

# Direct chord discriminant for the line through -S:
# m^4 - 6*xS*m^2 - 8*yS*m - 3*xS^2 - 4*A.
disc_generic = (
    m_generic**4
    - QKTU(6)*sx_generic*m_generic**2
    - QKTU(8)*sy_generic*m_generic
    - QKTU(3)*sx_generic**2
    - QKTU(4)*A6_generic
)
quartic_generic, square_factor_generic = squarefree_binary_quartic(
    disc_generic, QTU
)
assert quartic_generic.degree() == 4
I_generic, J_generic = binary_quartic_invariants(quartic_generic)
stdA_generic = QKU(-27*I_generic)
stdB_generic = QKU(-27*J_generic)

# Exact orientation functions from the separately passed certificate.
def ku_from_record(record):
    return QKU(QU([QQ(v) for v in record["numerator_coefficients_low_to_high"]])) / QKU(
        QU([QQ(v) for v in record["denominator_coefficients_low_to_high"]])
    )


wiv_generic = ku_from_record(orient["quartic"]["global_w_plus"])
c2_generic = ku_from_record(
    orient["minimalization"]["c2_equals_c_squared"]
)
c3_generic = ku_from_record(
    orient["minimalization"]["c3_equals_c_cubed"]
)

assert wiv_generic**2 == QKU(quartic_generic(tivq))

A13q = QU([
    QQ(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]
])
B13q = QU([
    QQ(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]
])
assert c2_generic**2 == QKU(A13q)/stdA_generic
assert c3_generic**2 == QKU(B13q)/stdB_generic

print(
    "Q8S3AJORIENT|generic_quartic=PASS|global_wIV=PASS|"
    "global_c2c3=PASS|status=PASS",
    flush=True,
)

# ===========================================================================
# Finite field and specialization helpers.
# ===========================================================================
p = ZZ(args.prime)
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
sx = KT(reduce_qpoly(sxq.numerator())) / KT(reduce_qpoly(sxq.denominator()))
sy = KT(reduce_qpoly(syq.numerator())) / KT(reduce_qpoly(syq.denominator()))
assert sy**2 == sx**3 + KT(A6)*sx + KT(B6)
tiv = modq(tivq)

section3 = bridge["third_section_canonical_q6"]


def rf_bridge(entry):
    return KT(RT([modq(v) for v in entry["numerator_coefficients_low_to_high"]])) / KT(
        RT([modq(v) for v in entry["denominator_coefficients_low_to_high"]])
    )


x3 = rf_bridge(section3["x"])
y3 = rf_bridge(section3["y"])
assert y3**2 == x3**3 + KT(A6)*x3 + KT(B6)

U3_num = RT([
    modq(v)
    for v in bridge["q8_parameter_on_third"][
        "numerator_coefficients_low_to_high"
    ]
])
U3_den = RT([
    modq(v)
    for v in bridge["q8_parameter_on_third"][
        "denominator_coefficients_low_to_high"
    ]
])
assert max(U3_num.degree(), U3_den.degree()) == 52

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
    return f, H, G, I, J


def newton_power_sums(poly):
    n = poly.degree()
    assert poly[n] == 1
    sums = [F(n)]
    for k in range(1, n):
        total = F(k)*poly[n-k]
        for j in range(1, k):
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


def sample_tau(tau):
    a13 = A13U(tau)
    b13 = B13U(tau)
    E13 = EllipticCurve(F, [0, 0, 0, a13, b13])
    if not E13.discriminant():
        return None, "singular_D13"
    if (r2**3 + a13*r2 + b13).roots():
        return None, "rational_2_torsion"

    H52 = RT(U3_num - tau*U3_den)
    if H52.degree() != 52:
        return None, f"degree_drop_{H52.degree()}"
    H52 = H52.monic()
    if H52.gcd(H52.derivative()).degree() != 0:
        return None, "non_etale"

    try:
        quartic = specialize_TU(quartic_generic, tau)
        square_factor = specialize_KTU(square_factor_generic, tau)
        q8_m = specialize_KTU(m_generic, tau)
        wiv = specialize_KU(wiv_generic, tau)
        c2 = specialize_KU(c2_generic, tau)
        c3 = specialize_KU(c3_generic, tau)
    except ZeroDivisionError:
        return None, "global_orientation_pole"

    if quartic.degree() != 4:
        return None, f"quartic_degree_{quartic.degree()}"
    if wiv*wiv != quartic(tiv):
        return None, "global_wIV_mismatch"

    # Verify global minimalization after specialization.
    fbin, HC, GC, I, J = native_covariants(quartic)
    if (-F(27)*I)*c2**2 != a13:
        return None, "global_c2_mismatch"
    if (-F(27)*J)*c3**2 != b13:
        return None, "global_c3_mismatch"

    def reduce_mod_H(value):
        value = KT(value)
        num = RT(value.numerator())
        den = RT(value.denominator())
        if den.gcd(H52).degree() != 0:
            raise ZeroDivisionError
        return (num*den.inverse_mod(H52)) % H52

    try:
        # Exact q8 fibre check against the already-certified bridge.
        m3 = (y3 + sy)/(x3 - sx)
        if reduce_mod_H(q8_m-m3):
            return None, "chord_mismatch"

        # IMPORTANT: square_factor is the specialization of the GENERIC
        # square factor, so w uses the same global quartic normalization.
        w3 = (2*x3 + sx - q8_m**2)/square_factor
        wA = reduce_mod_H(w3)
    except ZeroDivisionError:
        return None, "noninvertible_section_denominator"

    if (wA*wA - quartic) % H52:
        return None, "quartic_sqrt_mismatch"

    qx, qz = HC.parent().gens()
    H_u = RT(HC(qx=T, qz=F(1)))
    G_u = RT(GC(qx=T, qz=F(1)))
    quartic_inv = quartic.inverse_mod(H52)

    # Globally coherent scaling: c2=c^2 and c3=c^3.
    xA = (-F(3)/F(4) * H_u * quartic_inv * c2) % H52
    yA = (F(9)/F(32) * G_u * wA * quartic_inv**2 * c3) % H52
    if (yA*yA - xA*xA*xA - a13*xA - b13) % H52:
        return None, "etale_child_miss"

    # L(53 O) trace of the 52 covariant images.
    xpowers = [RT.one()]
    for unused in range(26):
        xpowers.append((xpowers[-1]*xA) % H52)
    columns = list(xpowers)
    columns += [(yA*xpowers[e]) % H52 for e in range(26)]
    assert len(columns) == 53

    Eval = matrix(F, 52, 53, lambda row, col: columns[col][row])
    ker = Eval.right_kernel().basis_matrix()
    if ker.nrows() != 1:
        return None, f"L53_kernel_{ker.nrows()}"
    rel = ker[0]

    XR = PolynomialRing(F, "X")
    Xv = XR.gen()
    Afun = sum(rel[i]*Xv**i for i in range(27))
    Bfun = sum(rel[27+i]*Xv**i for i in range(26))
    Rint = Afun**2 - (Xv**3 + a13*Xv + b13)*Bfun**2
    if Rint.degree() != 53:
        return None, f"residual_degree_{Rint.degree()}"

    root_sum = -Rint[52]/Rint[53]
    ps = newton_power_sums(H52)
    trace_x = sum(xA[i]*ps[i] for i in range(52))
    xQ = root_sum - trace_x
    if not Bfun(xQ):
        return None, "trace_B_zero"
    yQ = -Afun(xQ)/Bfun(xQ)
    traceCov = -E13(xQ, yQ)

    def map_origin(wvalue):
        fv = F(fbin(qx=tiv, qz=F(1)))
        hv = F(HC(qx=tiv, qz=F(1)))
        gv = F(GC(qx=tiv, qz=F(1)))
        if not fv or wvalue*wvalue != fv:
            raise ArithmeticError
        rawx = -F(3)/F(4)*hv/fv
        rawy = F(9)/F(32)*gv*wvalue/fv**2
        return E13(c2*rawx, c3*rawy)

    try:
        Qplus = map_origin(wiv)
        Qminus = map_origin(-wiv)
    except ArithmeticError:
        return None, "origin_map_failure"

    answers = []
    for origin in (Qplus, Qminus):
        doubled = traceCov - 52*origin
        halves = rational_halves(E13, doubled)
        if len(halves) != 1:
            return None, f"half_count_{len(halves)}"
        answers.append(halves[0])

    return tuple(answers), None


# ===========================================================================
# Collect globally coherent samples.
# ===========================================================================
samples = {"plus": [], "minus": []}
attempted = 0
candidate_integer = args.start
skip_counts = {}

while len(samples["plus"]) < args.samples:
    tau = F(candidate_integer)
    candidate_integer += 1
    attempted += 1
    result, reason = sample_tau(tau)
    if result is None:
        skip_counts[reason] = skip_counts.get(reason, 0) + 1
        continue

    plus, minus = result
    for label, point in (("plus", plus), ("minus", minus)):
        x, y = point.xy()
        samples[label].append((int(tau), int(x), int(y)))

    count = len(samples["plus"])
    if count <= 5 or count % 10 == 0 or count == args.samples:
        px, py = plus.xy()
        mx, my = minus.xy()
        print(
            f"Q8S3AJORIENT_SAMPLE|count={count}|U={int(tau)}|"
            f"plus={int(px)},{int(py)}|minus={int(mx)},{int(my)}|status=PASS",
            flush=True,
        )

print(
    f"Q8S3AJORIENT_SAMPLE|good={args.samples}|attempted={attempted}|"
    f"skips={skip_counts}|stage=collection|status=PASS",
    flush=True,
)

# ===========================================================================
# Target-profile interpolation.
# ===========================================================================
def interpolate_x(rows, num_degree=48, den_degree=44):
    ncols = (num_degree+1) + (den_degree+1)
    interp = matrix(
        F,
        len(rows),
        ncols,
        lambda row, col: (
            F(rows[row][0])**col
            if col <= num_degree
            else -F(rows[row][1])
            * F(rows[row][0])**(col-(num_degree+1))
        ),
    )
    ker = interp.right_kernel().basis_matrix()
    if ker.nrows() != 1:
        return None, int(ker.nrows())
    v = ker[0]
    P = RU(list(v[:num_degree+1]))
    D = RU(list(v[num_degree+1:]))
    if not D:
        return None, -1
    scale = D.leading_coefficient()
    return (RU(P/scale), RU(D/scale)), 1


def square_root_monic(poly):
    poly = RU(poly)
    if not poly or not poly.is_monic():
        return None
    result = RU.one()
    for factor, multiplicity in poly.factor():
        if multiplicity % 2:
            return None
        result *= factor.monic()**(multiplicity//2)
    return result.monic()


branch_results = {}
for label in ("plus", "minus"):
    rows = [(u0, x0) for u0, x0, unused_y in samples[label]]
    recovered, dim = interpolate_x(rows)

    if recovered is None:
        print(
            f"Q8S3AJORIENT_INTERP|branch={label}|kernel_dim={dim}|"
            "profile=48/44|status=NO_FIT",
            flush=True,
        )
        branch_results[label] = {"fit": False, "kernel_dim": dim}
        continue

    X, D = recovered
    if X.gcd(D).degree() != 0:
        print(
            f"Q8S3AJORIENT_INTERP|branch={label}|kernel_dim=1|"
            "status=NONREDUCED",
            flush=True,
        )
        branch_results[label] = {"fit": False}
        continue

    Z = square_root_monic(D)
    if (
        X.degree() != 48
        or D.degree() != 44
        or Z is None
        or Z.degree() != 22
    ):
        print(
            f"Q8S3AJORIENT_INTERP|branch={label}|kernel_dim=1|"
            f"x_degrees={X.degree()}/{D.degree()}|"
            f"Z_degree={-1 if Z is None else Z.degree()}|"
            "status=WRONG_PROFILE",
            flush=True,
        )
        branch_results[label] = {"fit": False}
        continue

    RHS = X**3 + A13U*X*Z**4 + B13U*Z**6
    if not RHS.is_square():
        print(
            f"Q8S3AJORIENT_INTERP|branch={label}|x=48/44|Z=22|"
            "rhs_square=0|status=WRONG_PROFILE",
            flush=True,
        )
        branch_results[label] = {"fit": False}
        continue

    Y = RHS.sqrt()
    if Y.degree() != 72:
        print(
            f"Q8S3AJORIENT_INTERP|branch={label}|x=48/44|Z=22|"
            f"Y_degree={Y.degree()}|status=WRONG_PROFILE",
            flush=True,
        )
        branch_results[label] = {"fit": False}
        continue

    direct = opposite = True
    for u0, unused_x0, y0 in samples[label]:
        uu = F(u0)
        if not Z(uu):
            raise ArithmeticError("sample unexpectedly lies at reconstructed pole")
        predicted = Y(uu)/Z(uu)**3
        direct &= predicted == F(y0)
        opposite &= -predicted == F(y0)

    if direct == opposite:
        raise ArithmeticError(f"branch {label}: Y orientation unresolved")
    if opposite:
        Y = -Y

    assert Y**2 == X**3 + A13U*X*Z**4 + B13U*Z**6
    for u0, x0, y0 in samples[label]:
        uu = F(u0)
        assert X(uu)/Z(uu)**2 == F(x0)
        assert Y(uu)/Z(uu)**3 == F(y0)

    print(
        f"Q8S3AJORIENT_INTERP|branch={label}|kernel_dim=1|"
        "x=48/44|Z=22|Y=72|identity=PASS|samples=PASS|"
        "status=TARGET_PROFILE",
        flush=True,
    )
    branch_results[label] = {"fit": True, "X": X, "Y": Y, "Z": Z}


winners = [label for label, data in branch_results.items() if data.get("fit")]
if len(winners) != 1:
    print(
        f"Q8S3AJORIENT_RESULT|winners={','.join(winners) if winners else 'none'}|"
        "status=ORIGIN_NOT_RESOLVED",
        flush=True,
    )
    raise SystemExit(0)

winner = winners[0]
X = branch_results[winner]["X"]
Y = branch_results[winner]["Y"]
Z = branch_results[winner]["Z"]

payload = {
    "schema": "elkies-k3.h92-q8-s3-aj-oriented-modp.v1",
    "status": "PASS_MODULAR_Q8_S3_AJ_TARGET_SECTION",
    "prime": int(p),
    "samples": int(args.samples),
    "attempted": int(attempted),
    "skip_counts": skip_counts,
    "global_orientation": {
        "artifact": str(ORIENT.relative_to(ROOT)),
        "sha256": hashlib.sha256(ORIENT.read_bytes()).hexdigest(),
        "resolved_origin": winner,
    },
    "target": {
        "D13_MW_coordinates": [0, -1, 1, 1],
        "height": "47",
        "local_D13_correction": "1",
        "O_intersection": 22,
        "coordinate_profile": {
            "Z_degree": 22,
            "X_degree": 48,
            "Y_degree": 72,
            "x_degrees": [48, 44],
            "y_degrees": [72, 66],
        },
    },
    "section_mod_p": {
        "Z_coefficients_low_to_high": [int(v) for v in Z.list()],
        "X_coefficients_low_to_high": [int(v) for v in X.list()],
        "Y_coefficients_low_to_high": [int(v) for v in Y.list()],
    },
    "boundary": (
        "This resolves the globally coherent q8 origin modulo p by the pinned "
        "D13 target profile and constructs the full modular target section. "
        "Characteristic-zero Hensel lifting remains separate."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    f"Q8S3AJORIENT_RESULT|origin={winner}|D13_AJ=0,-1,1,1|"
    "height=47|O=22|x=48/44|y=72/66|"
    "status=PASS_MODULAR_Q8_S3_AJ_TARGET_SECTION",
    flush=True,
)
