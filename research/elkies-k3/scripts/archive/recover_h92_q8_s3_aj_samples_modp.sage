#!/usr/bin/env sage -python
"""
Sample both q8 IV* origin branches and recover the modular D13 AJ section.

The single-specialization smoke test has already passed:
  degree-52 S3 fibre -> q8 covariant transport -> L(53O) trace -> halving.

This script repeats that exact computation over many good U=tau values.
It keeps both possible q8 origins ("plus" and "minus") and then tests each
branch against the pinned D13 target profile for MW vector (0,-1,1,1):

    height = 47
    local D13 correction = 1
    P.O = 22

hence, in the displayed short Weierstrass D13 model,

    x = X/Z^2,  deg(X,Z^2) = (48,44),
    y = Y/Z^3,  deg(Y,Z^3) = (72,66).

For each branch we rationally interpolate x(U) with bounds 48/44.  A branch
is accepted only if:
  * the interpolation kernel is exactly one-dimensional;
  * the denominator has exact degree 44 and is a square Z^2 with deg Z=22;
  * the numerator has exact degree 48 and is coprime to Z;
  * X^3 + A*X*Z^4 + B*Z^6 is an exact square Y^2 in GF(p)[U];
  * deg Y = 72;
  * one sign of Y/Z^3 agrees with every sampled y-coordinate.

This both resolves the IV* origin and constructs the complete modular D13
section corresponding to AJ_q8(S3).

Run:
  sage -python ~/Downloads/recover_h92_q8_s3_aj_samples_modp.sage

Optional:
  --prime 100003 --start 2 --samples 100
"""

import argparse
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
    raise ValueError("need at least 94 samples for a 48/44 rational interpolation")

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
Q6 = GEN / "elkies-k3-h92-q6-child-jacobian.json"
Q8 = GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json"
BRIDGE = LOCAL / "q6-third-to-q8-bridge.json"
OUTPUT = (
    args.output.resolve()
    if args.output and args.output.is_absolute()
    else ROOT / (
        args.output
        if args.output
        else Path("artifacts/local/elkies-k3/q8-s3-aj-samples-mod-100003.json")
    )
)

for path in (CORE, Q6, Q8, BRIDGE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

scope = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), scope)
squarefree_binary_quartic = scope["squarefree_binary_quartic"]
binary_quartic_invariants = scope["binary_quartic_invariants"]

q6 = json.loads(Q6.read_text())
q8 = json.loads(Q8.read_text())
bridge = json.loads(BRIDGE.read_text())
assert q6["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert bridge["status"] == "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"
assert bridge["q8_parameter_on_third"]["degree"] == 52

p = ZZ(args.prime)
F = GF(p)
RT = PolynomialRing(F, "T")
T = RT.gen()
KT = RT.fraction_field()
QQTR = PolynomialRing(QQ, "T")


def modq(value):
    value = QQ(value)
    den = ZZ(value.denominator())
    if den % p == 0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {value}")
    return F(ZZ(value.numerator())) / F(den)


def poly_strings(values):
    return RT([modq(value) for value in values])


def rf_from_bridge(entry):
    return KT(poly_strings(entry["numerator_coefficients_low_to_high"])) / KT(
        poly_strings(entry["denominator_coefficients_low_to_high"])
    )


def reduce_factor_string(text):
    source = QQTR(text)
    return RT([modq(value) for value in source.list()]).monic()


# ---------------------------------------------------------------------------
# Fixed q6 data, exact S3 and exact corrected q8 pencil.
# ---------------------------------------------------------------------------
A6 = poly_strings(q6["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B6 = poly_strings(q6["minimal_short_weierstrass"]["B_coefficients_low_to_high"])

section3 = bridge["third_section_canonical_q6"]
x3 = rf_from_bridge(section3["x"])
y3 = rf_from_bridge(section3["y"])
assert y3**2 == x3**3 + KT(A6)*x3 + KT(B6)

U3_num = poly_strings(
    bridge["q8_parameter_on_third"]["numerator_coefficients_low_to_high"]
)
U3_den = poly_strings(
    bridge["q8_parameter_on_third"]["denominator_coefficients_low_to_high"]
)
assert max(U3_num.degree(), U3_den.degree()) == 52

mdata = q8["marking"]["section"]
sx = KT(poly_strings(mdata["x_numerator_coefficients_low_to_high"])) / KT(
    poly_strings(mdata["x_denominator_coefficients_low_to_high"])
)
sy = KT(poly_strings(mdata["y_numerator_coefficients_low_to_high"])) / KT(
    poly_strings(mdata["y_denominator_coefficients_low_to_high"])
)
assert sy**2 == sx**3 + KT(A6)*sx + KT(B6)


def monic_power_root(value, exponent):
    result = RT.one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        result *= factor.monic()**(multiplicity//exponent)
    return result.monic()


nx, dx = RT(sx.numerator()), RT(sx.denominator())
ny, dy = RT(sy.numerator()), RT(sy.denominator())
h = monic_power_root(dx, 2)
assert h == monic_power_root(dy, 3)
assert h.degree() == 10

ii = reduce_factor_string(
    next(item for item in q6["finite_fibres"] if item["kodaira"] == "II*")["factor"]
)
iv = reduce_factor_string(
    next(item for item in q6["finite_fibres"] if item["kodaira"] == "IV*")["factor"]
)
assert ii.degree() == iv.degree() == 1
M = (ii**2 * iv**2).monic()
t_iv = -iv[0]/iv[1]

normalizer = (ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0
p_fun = -sy/sx
rho = (normalizer*nx.inverse_mod(M)).mod(M)

kernel_polys = q8["rr"]["kernel_polynomials"]
assert len(kernel_polys) == 2
pairs = []
for entry in kernel_polys:
    sp = RT([modq(v) for v in QQTR(entry["s"]).list()])
    tp = RT([modq(v) for v in QQTR(entry["t"]).list()])
    Bcoef = KT(sp)/KT(h)
    Acoef = (
        -KT(sp)*p_fun/KT(h)
        - KT(sp)*KT(normalizer)/KT(nx)
        + KT(sp*rho)
        + KT(tp*M)
    )
    pairs.append((Acoef, Bcoef))
(A0, B0), (A1, B1) = pairs

# D13 child polynomials in U.
RU = PolynomialRing(F, "U")
U = RU.gen()
KU = RU.fraction_field()
A13U = RU([modq(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]])
B13U = RU([modq(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]])

# Rational 2-torsion gate at finite specializations.
R2 = PolynomialRing(F, "r2")
r2 = R2.gen()


def native_covariants(quartic):
    base = quartic.base_ring()
    BR = PolynomialRing(base, names=("qx", "qz"))
    qx, qz = BR.gens()
    f = sum(
        BR(quartic[i])*qx**i*qz**(4-i)
        for i in range(5)
    )
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


def fourth_sixth_unit(std_a, std_b, target_a, target_b):
    CR = PolynomialRing(F, "c")
    c = CR.gen()
    roots = [
        root for root, mult in (c**4-target_a/std_a).roots()
        if mult == 1 and std_b*root**6 == target_b
    ]
    if not roots:
        raise ArithmeticError("no q8 fourth/sixth minimizing unit")
    return roots[0]


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
    """Return (plus_point, minus_point) or (None, reason)."""
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

    # Exact q8 chord on the U=tau fibre.
    q8_m = -(A1 - tau*A0)/(B1 - tau*B0)
    radicand = (
        q8_m**4
        - 6*sx*q8_m**2
        - 8*sy*q8_m
        - 3*sx**2
        - 4*KT(A6)
    )
    quartic, square_factor = squarefree_binary_quartic(radicand, RT)
    if quartic.degree() != 4:
        return None, f"quartic_degree_{quartic.degree()}"

    fbin, HC, GC, I, J = native_covariants(quartic)
    stdA, stdB = -F(27)*I, -F(27)*J
    try:
        unit = fourth_sixth_unit(stdA, stdB, a13, b13)
    except ArithmeticError:
        return None, "no_minimizing_unit"

    def reduce_mod_H(value):
        value = KT(value)
        num = RT(value.numerator())
        den = RT(value.denominator())
        if den.gcd(H52).degree() != 0:
            raise ZeroDivisionError
        return (num*den.inverse_mod(H52)) % H52

    try:
        m3 = (y3 + sy)/(x3 - sx)
        if reduce_mod_H(q8_m-m3):
            return None, "chord_mismatch"

        w3 = (2*x3 + sx - q8_m**2)/square_factor
        wA = reduce_mod_H(w3)
    except ZeroDivisionError:
        return None, "noninvertible_section_denominator"

    if (wA*wA - quartic) % H52:
        return None, "quartic_sqrt_mismatch"

    # Covariant images of all 52 conjugates in the etale algebra.
    qx, qz = HC.parent().gens()
    H_u = RT(HC(qx=T, qz=F(1)))
    G_u = RT(GC(qx=T, qz=F(1)))
    quartic_inv = quartic.inverse_mod(H52)
    xA = (-F(3)/F(4) * H_u * quartic_inv * unit**2) % H52
    yA = (F(9)/F(32) * G_u * wA * quartic_inv**2 * unit**3) % H52
    if (yA*yA - xA*xA*xA - a13*xA - b13) % H52:
        return None, "etale_child_miss"

    # L(53O) residual-intersection trace.
    one = RT.one()
    xpowers = [one]
    for unused in range(26):
        xpowers.append((xpowers[-1]*xA) % H52)
    columns = list(xpowers)
    for exponent in range(26):
        columns.append((yA*xpowers[exponent]) % H52)
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
    bQ = Bfun(xQ)
    if not bQ:
        return None, "trace_B_zero"
    yQ = -Afun(xQ)/bQ
    traceCov = -E13(xQ, yQ)

    # Two candidate q8 section origins at old IV*.
    qiv = quartic(t_iv)
    if not qiv or not qiv.is_square():
        return None, "IV_origin_not_rational"
    wroot = qiv.sqrt()

    def map_quartic_point(tvalue, wvalue):
        fv = F(fbin(qx=F(tvalue), qz=F(1)))
        hv = F(HC(qx=F(tvalue), qz=F(1)))
        gv = F(GC(qx=F(tvalue), qz=F(1)))
        if not fv or F(wvalue)**2 != fv:
            raise ArithmeticError
        rx = -F(3)/F(4)*hv/fv
        ry = F(9)/F(32)*gv*F(wvalue)/fv**2
        return E13(unit**2*rx, unit**3*ry)

    try:
        Qplus = map_quartic_point(t_iv, wroot)
        Qminus = map_quartic_point(t_iv, -wroot)
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


# ---------------------------------------------------------------------------
# Collect samples.
# ---------------------------------------------------------------------------
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
            f"Q8S3AJSAMPLE|count={count}|U={int(tau)}|"
            f"plus={int(px)},{int(py)}|minus={int(mx)},{int(my)}|status=PASS",
            flush=True,
        )

print(
    f"Q8S3AJSAMPLE|good={args.samples}|attempted={attempted}|"
    f"skips={skip_counts}|stage=collection|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Rational interpolation x=P/Q with bounds 48/44.
# ---------------------------------------------------------------------------
def interpolate_x(rows, num_degree=48, den_degree=44):
    ncols = (num_degree+1) + (den_degree+1)
    interp = matrix(
        F,
        len(rows),
        ncols,
        lambda row, col: (
            F(rows[row][0])**col
            if col <= num_degree
            else -F(rows[row][1]) * F(rows[row][0])**(col-(num_degree+1))
        ),
    )
    ker = interp.right_kernel().basis_matrix()
    if ker.nrows() != 1:
        return None, int(ker.nrows())

    v = ker[0]
    P = RU(list(v[:num_degree+1]))
    Q = RU(list(v[num_degree+1:]))
    if not Q:
        return None, -1

    scale = Q.leading_coefficient()
    P /= scale
    Q /= scale
    return (RU(P), RU(Q)), 1


def square_root_monic(poly):
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
    rows = [(u, x) for u, x, unused_y in samples[label]]
    recovered, dim = interpolate_x(rows)

    if recovered is None:
        print(
            f"Q8S3AJINTERP|branch={label}|kernel_dim={dim}|"
            "profile=48/44|status=NO_FIT",
            flush=True,
        )
        branch_results[label] = {"fit": False, "kernel_dim": dim}
        continue

    X, D = recovered
    if X.gcd(D).degree() != 0:
        print(
            f"Q8S3AJINTERP|branch={label}|kernel_dim=1|"
            "status=NONREDUCED",
            flush=True,
        )
        branch_results[label] = {"fit": False, "kernel_dim": 1}
        continue

    Z = square_root_monic(D)
    square_den = Z is not None
    degrees_ok = (
        X.degree() == 48
        and D.degree() == 44
        and square_den
        and Z.degree() == 22
    )

    if not degrees_ok:
        print(
            f"Q8S3AJINTERP|branch={label}|kernel_dim=1|"
            f"x_degrees={X.degree()}/{D.degree()}|"
            f"square_den={int(square_den)}|"
            f"Z_degree={-1 if Z is None else Z.degree()}|status=WRONG_PROFILE",
            flush=True,
        )
        branch_results[label] = {
            "fit": False,
            "kernel_dim": 1,
            "x_degrees": [int(X.degree()), int(D.degree())],
            "square_den": bool(square_den),
        }
        continue

    RHS = X**3 + A13U*X*Z**4 + B13U*Z**6
    if not RHS.is_square():
        print(
            f"Q8S3AJINTERP|branch={label}|kernel_dim=1|x=48/44|Z=22|"
            "rhs_square=0|status=WRONG_PROFILE",
            flush=True,
        )
        branch_results[label] = {"fit": False, "kernel_dim": 1}
        continue

    Y = RHS.sqrt()
    if Y.degree() != 72:
        print(
            f"Q8S3AJINTERP|branch={label}|kernel_dim=1|x=48/44|Z=22|"
            f"Y_degree={Y.degree()}|status=WRONG_PROFILE",
            flush=True,
        )
        branch_results[label] = {"fit": False, "kernel_dim": 1}
        continue

    # Determine the global Y sign from all samples.
    direct = True
    opposite = True
    for u0, unused_x0, y0 in samples[label]:
        uu = F(u0)
        if not Z(uu):
            raise ArithmeticError("sample unexpectedly lies at reconstructed pole")
        predicted = Y(uu)/Z(uu)**3
        direct &= predicted == F(y0)
        opposite &= -predicted == F(y0)

    if direct == opposite:
        raise ArithmeticError(
            f"branch {label}: global Y orientation was not uniquely determined"
        )
    if opposite:
        Y = -Y

    # Replay full curve identity and all x/y samples.
    assert Y**2 == X**3 + A13U*X*Z**4 + B13U*Z**6
    for u0, x0, y0 in samples[label]:
        uu = F(u0)
        assert X(uu)/Z(uu)**2 == F(x0)
        assert Y(uu)/Z(uu)**3 == F(y0)

    print(
        f"Q8S3AJINTERP|branch={label}|kernel_dim=1|"
        "x=48/44|Z=22|Y=72|identity=PASS|samples=PASS|"
        "status=TARGET_PROFILE",
        flush=True,
    )
    branch_results[label] = {
        "fit": True,
        "kernel_dim": 1,
        "X": X,
        "Z": Z,
        "Y": Y,
    }


winners = [label for label, data in branch_results.items() if data.get("fit")]
if len(winners) != 1:
    print(
        f"Q8S3AJSAMPLE_RESULT|winners={','.join(winners) if winners else 'none'}|"
        "status=ORIGIN_NOT_RESOLVED",
        flush=True,
    )
    raise SystemExit(0)

winner = winners[0]
result = branch_results[winner]
X, Z, Y = result["X"], result["Z"], result["Y"]

payload = {
    "schema": "elkies-k3.h92-q8-s3-aj-samples-modp.v1",
    "status": "PASS_MODULAR_Q8_S3_AJ_TARGET_SECTION",
    "prime": int(p),
    "samples": int(args.samples),
    "attempted": int(attempted),
    "skip_counts": skip_counts,
    "resolved_origin": winner,
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
    "sample_values": {
        label: [
            {"U": u0, "x": x0, "y": y0}
            for u0, x0, y0 in rows
        ]
        for label, rows in samples.items()
    },
    "boundary": (
        "This resolves the q8 IV* origin modulo p by the pinned D13 target "
        "pole profile and constructs the complete modular target section. "
        "Characteristic-zero lifting/certification remains separate."
    ),
}
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    f"Q8S3AJSAMPLE_RESULT|origin={winner}|"
    "D13_AJ=0,-1,1,1|height=47|O=22|x=48/44|y=72/66|"
    "status=PASS_MODULAR_Q8_S3_AJ_TARGET_SECTION",
    flush=True,
)
