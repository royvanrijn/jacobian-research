#!/usr/bin/env sage -python
"""
Recover the direct anchored q8 Abel-Jacobi x-coordinate of H92 S3 over GF(p).

Prerequisites:
  * PASS_EXACT_Q8_EQUATION_NS_DIVISOR
  * PASS_EXACT_D13_BRANCH_ANCHOR
  * PASS_EXACT_Q8_OLD_Q6_MW_IMAGES
  * PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52

The q8 zero is the actual II*_E8_1 branch section.  For many U=tau values,
map the degree-52 S3 divisor birationally through the anchored branch-point
quartic -> D13 map, sum the 52 points with L(53 O), and retain AJ_x(tau).

The anchored lattice profile determines P.O=22, so on the minimal K3 D13 model
the target x-coordinate has profile

    deg numerator = 48
    deg denominator = 44 = 2*22.

93 independent samples determine the rational function; extra samples certify it.
No covariant 2-cover, IV* origin subtraction, or halving is used.
"""

import argparse
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, sage_eval, vector


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
parser.add_argument("--samples", type=int, default=105)
parser.add_argument("--scan-limit", type=int, default=500)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if args.samples < 94:
    raise ValueError("need at least 94 good samples for a certified 48/44 interpolation")

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
Q6 = GEN / "elkies-k3-h92-q6-child-jacobian.json"

# Require the complete corrected-child schema; the generated-results copy can be stale.
q8_candidates = [
    LOCAL / "q8-corrected2cover-qq-child.json",
    GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8 = next(
    (
        path for path in q8_candidates
        if path.exists()
        and "branch_quartic" in json.loads(path.read_text()).get("pencil", {})
    ),
    None,
)
if Q8 is None:
    raise SystemExit("No complete corrected q8 child artifact with pencil.branch_quartic")

BRIDGE = LOCAL / "q6-third-to-q8-bridge.json"
ANCHOR = LOCAL / "q8-d13-branch-anchor.json"
EQNS = LOCAL / "q8-equation-ns-divisor.json"
OLDMW = LOCAL / "q8-old-q6-mw-images.json"
SEED = LOCAL / "q8-s3-direct-anchor-trace-mod-100003-tau-2.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q8-s3-direct-x-mod-{args.prime}.json"
)

for path in (CORE, Q6, Q8, BRIDGE, ANCHOR, EQNS, OLDMW):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

scope = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), scope)
squarefree_binary_quartic = scope["squarefree_binary_quartic"]

q6 = json.loads(Q6.read_text())
q8 = json.loads(Q8.read_text())
bridge = json.loads(BRIDGE.read_text())
anchor = json.loads(ANCHOR.read_text())
eqns = json.loads(EQNS.read_text())
oldmw = json.loads(OLDMW.read_text())

assert q6["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert bridge["status"] == "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"
assert anchor["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"
assert eqns["status"] == "PASS_EXACT_Q8_EQUATION_NS_DIVISOR"
assert oldmw["status"] == "PASS_EXACT_Q8_OLD_Q6_MW_IMAGES"
assert eqns["q8_equation_fibre"]["S3_degree"] == 52
assert eqns["q8_equation_fibre"]["root_data"] == [13, 312, 4]

# Generator 3 is the already-certified q6 third section S3.
s3prof = next(row for row in oldmw["old_q6_generators"] if int(row["index"]) == 3)
assert int(s3prof["q8_multisection_degree"]) == 52
print(
    "Q8S3MODP_PROFILE|"
    f"mw={','.join(map(str,s3prof['mw_coordinates']))}|"
    f"height={s3prof['height']}|"
    f"class_order={s3prof['D13_class_order']}|"
    f"correction={s3prof['D13_local_correction']}|"
    f"PdotO={s3prof['jacobian_section_P_dot_O']}|"
    f"P_q6_degree={s3prof['jacobian_section_q6_degree']}|"
    f"P_h3_degree={s3prof['jacobian_section_source_H3_degree']}|"
    "status=PASS_ANCHORED_PROFILE",
    flush=True,
)

# This is the profile expected from the historical pinned target, now tested in
# the actual II*_E8_1 anchored frame.
expected_profile = {
    "mw_coordinates": [0, -1, 1, 1],
    "height": "47",
    "D13_local_correction": "1",
    "jacobian_section_P_dot_O": 22,
}
for key, expected in expected_profile.items():
    if s3prof[key] != expected:
        raise ArithmeticError(
            f"anchored S3 profile mismatch for {key}: {s3prof[key]} != {expected}"
        )

p = ZZ(args.prime)
if not p.is_prime() or p in (2, 3):
    raise ValueError("prime must be odd and != 3")
F = GF(p)
R = PolynomialRing(F, "T")
T = R.gen()
K = R.fraction_field()


def modq(value):
    value = QQ(value)
    den = ZZ(value.denominator())
    if den % p == 0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {value}")
    return F(ZZ(value.numerator())) / F(den)


def poly_strings(values):
    return R([modq(v) for v in values])


def rf_from_bridge(entry):
    return K(poly_strings(entry["numerator_coefficients_low_to_high"])) / K(
        poly_strings(entry["denominator_coefficients_low_to_high"])
    )


# q6 curve and exact S3.
A6 = poly_strings(q6["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
B6 = poly_strings(q6["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
section3 = bridge["third_section_canonical_q6"]
x3 = rf_from_bridge(section3["x"])
y3 = rf_from_bridge(section3["y"])
assert y3**2 == x3**3 + K(A6)*x3 + K(B6)

U3_num = poly_strings(
    bridge["q8_parameter_on_third"]["numerator_coefficients_low_to_high"]
)
U3_den = poly_strings(
    bridge["q8_parameter_on_third"]["denominator_coefficients_low_to_high"]
)
assert max(U3_num.degree(), U3_den.degree()) == 52

# Corrected marked q8 section and pencil frame.
mdata = q8["marking"]["section"]
sx = K(poly_strings(mdata["x_numerator_coefficients_low_to_high"])) / K(
    poly_strings(mdata["x_denominator_coefficients_low_to_high"])
)
sy = K(poly_strings(mdata["y_numerator_coefficients_low_to_high"])) / K(
    poly_strings(mdata["y_denominator_coefficients_low_to_high"])
)
assert sy**2 == sx**3 + K(A6)*sx + K(B6)


def monic_power_root(value, exponent):
    out = R.one()
    for factor, mult in value.factor():
        assert mult % exponent == 0
        out *= factor.monic()**(mult//exponent)
    return out.monic()


nx, dx = R(sx.numerator()), R(sx.denominator())
ny, dy = R(sy.numerator()), R(sy.denominator())
h = monic_power_root(dx, 2)
assert h == monic_power_root(dy, 3)
assert h.degree() == 10

QQTR = PolynomialRing(QQ, "T")


def reduce_factor_string(text):
    src = QQTR(text)
    return R([modq(v) for v in src.list()]).monic()


ii = reduce_factor_string(
    next(item for item in q6["finite_fibres"] if item["kodaira"] == "II*")["factor"]
)
iv = reduce_factor_string(
    next(item for item in q6["finite_fibres"] if item["kodaira"] == "IV*")["factor"]
)
M = (ii**2 * iv**2).monic()
normalizer = (ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
p_fun = -sy/sx
rho = (normalizer*nx.inverse_mod(M)).mod(M)

pairs = []
for entry in q8["rr"]["kernel_polynomials"]:
    sp = R([modq(v) for v in QQTR(entry["s"]).list()])
    tp = R([modq(v) for v in QQTR(entry["t"]).list()])
    Bcoef = K(sp)/K(h)
    Acoef = (
        -K(sp)*p_fun/K(h)
        -K(sp)*K(normalizer)/K(nx)
        +K(sp*rho)
        +K(tp*M)
    )
    pairs.append((Acoef, Bcoef))
(A0, B0), (A1, B1) = pairs

# D13 child polynomials over F[U].
UR = PolynomialRing(F, "U")
U = UR.gen()
A13_poly = UR([modq(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]])
B13_poly = UR([modq(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]])

# Parse the exact branch-anchor rational functions over QQ(U), once.
QUQ = PolynomialRing(QQ, "U")
UQ = QUQ.gen()
KUQ = QUQ.fraction_field()


def parse_u(text):
    return KUQ(sage_eval(str(text), locals={"U": UQ}))


def spec_u(text, tau):
    value = parse_u(text)
    num = QUQ(value.numerator())
    den = QUQ(value.denominator())
    nv = sum(modq(c)*tau**i for i, c in enumerate(num.list()))
    dv = sum(modq(c)*tau**i for i, c in enumerate(den.list()))
    if not dv:
        raise ZeroDivisionError("branch-anchor U denominator vanished")
    return nv/dv


tii = modq(QQ(anchor["zero"]["old_base_T"]))
coef = anchor["quartic_to_anchor"]["shifted_coefficients"]
urst = anchor["anchor_to_canonical"]["urst"]


def reduce_mod_H(value, H):
    value = K(value)
    num = R(value.numerator())
    den = R(value.denominator())
    if den.gcd(H).degree() != 0:
        raise ZeroDivisionError("denominator not invertible modulo degree-52 fibre")
    return (num*den.inverse_mod(H)) % H


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


def direct_x_at(tau_int):
    tau = F(tau_int)

    A13 = A13_poly(tau)
    B13 = B13_poly(tau)
    E13 = EllipticCurve(F, [0, 0, 0, A13, B13])
    if not E13.discriminant():
        raise ArithmeticError("singular D13 child fibre")

    H = R(U3_num - tau*U3_den)
    if H.degree() != 52:
        raise ArithmeticError("degree-52 fibre dropped")
    H = H.monic()
    if H.gcd(H.derivative()).degree() != 0:
        raise ArithmeticError("S3 degree-52 fibre is not etale")

    q8_m = -(A1 - tau*A0)/(B1 - tau*B0)
    radicand = q8_m**4 - 6*sx*q8_m**2 - 8*sy*q8_m - 3*sx**2 - 4*K(A6)
    quartic, square_factor = squarefree_binary_quartic(radicand, R)
    if quartic.degree() != 4:
        raise ArithmeticError("q8 quartic degree dropped")

    m3 = (y3 + sy)/(x3 - sx)
    if reduce_mod_H(q8_m - m3, H):
        raise ArithmeticError("q8 chord mismatch on S3 fibre")

    w3 = (2*x3 + sx - q8_m**2)/square_factor
    wA = reduce_mod_H(w3, H)
    if (wA*wA - quartic) % H:
        raise ArithmeticError("S3 quartic square-root mismatch")

    a = spec_u(coef["a_r4"], tau)
    b = spec_u(coef["b_r3"], tau)
    c = spec_u(coef["c_r2"], tau)
    d = spec_u(coef["d_r1"], tau)

    rpoly = T - tii
    branch_poly = d*rpoly + c*rpoly**2 + b*rpoly**3 + a*rpoly**4
    if branch_poly.degree() != 4:
        raise ArithmeticError("anchored branch quartic degree dropped")

    scale = quartic[4]/branch_poly[4]
    if quartic != scale*branch_poly:
        raise ArithmeticError("anchored and direct quartics are not scalar-equivalent")
    if not scale.is_square():
        raise ArithmeticError("quartic scalar is nonsquare at this tau")

    # The choice of square root can flip every point, but AJ_x is invariant
    # under global negation, so independent finite-field sqrt signs are harmless.
    scale_root = scale.sqrt()
    wbA = wA/scale_root
    if (wbA*wbA - branch_poly) % H:
        raise ArithmeticError("branch W conversion failed")

    rA = rpoly % H
    if rA.gcd(H).degree() != 0:
        raise ArithmeticError("S3 divisor meets branch zero")
    rInv = rA.inverse_mod(H)
    Xa = (d*rInv) % H
    Ya = (d*wbA*rInv**2) % H

    u = spec_u(urst[0], tau)
    rr = spec_u(urst[1], tau)
    ss = spec_u(urst[2], tau)
    tt = spec_u(urst[3], tau)
    if not u:
        raise ArithmeticError("anchor isomorphism u vanished")

    xA = ((Xa - rr)/(u**2)) % H
    yA = ((Ya - ss*(Xa - rr) - tt)/(u**3)) % H

    if (yA*yA - xA*xA*xA - A13*xA - B13) % H:
        Ya = (-d*wbA*rInv**2) % H
        yA = ((Ya - ss*(Xa - rr) - tt)/(u**3)) % H
        if (yA*yA - xA*xA*xA - A13*xA - B13) % H:
            raise ArithmeticError("direct images miss D13 child")

    # L(53 O): 1,x,...,x^26,y,xy,...,x^25y.
    one = R.one()
    xp = [one]
    for unused in range(26):
        xp.append((xp[-1]*xA) % H)
    columns = list(xp) + [(yA*xp[e]) % H for e in range(26)]
    assert len(columns) == 53

    Eval = matrix(F, 52, 53, lambda row, col: columns[col][row])
    ker = Eval.right_kernel().basis_matrix()
    if ker.nrows() != 1:
        raise ArithmeticError(f"L(53O) trace kernel dimension {ker.nrows()}")

    rel = ker[0]
    XR = PolynomialRing(F, "X")
    X = XR.gen()
    Afun = sum(rel[i]*X**i for i in range(27))
    Bfun = sum(rel[27+i]*X**i for i in range(26))
    Rint = Afun**2 - (X**3 + A13*X + B13)*Bfun**2
    if Rint.degree() != 53:
        raise ArithmeticError("residual intersection degree is not 53")

    root_sum = -Rint[52]/Rint[53]
    ps = newton_power_sums(H)
    trace_x = sum(xA[i]*ps[i] for i in range(52))
    xQ = root_sum - trace_x

    bQ = Bfun(xQ)
    if not bQ:
        raise ArithmeticError("trace residual has B(x_Q)=0")
    yQ = -Afun(xQ)/bQ
    Qres = E13(xQ, yQ)
    AJ = -Qres
    if AJ.is_zero():
        raise ArithmeticError("AJ unexpectedly zero")
    ax, ay = AJ.xy()
    return F(ax), F(ay)


samples = []
skips = {}
tau = int(args.start)
stop = int(args.start) + int(args.scan_limit)

while tau < stop and len(samples) < args.samples:
    try:
        ax, ay = direct_x_at(tau)
        samples.append((F(tau), ax, ay))
        if tau == 2 or len(samples) % 10 == 0:
            print(
                f"Q8S3MODP_SAMPLE|good={len(samples)}|tau={tau}|"
                f"x={int(ax)}|y_local={int(ay)}|status=PASS",
                flush=True,
            )
    except Exception as exc:
        key = type(exc).__name__
        skips[key] = skips.get(key, 0) + 1
    tau += 1

if len(samples) < args.samples:
    raise ArithmeticError(
        f"only {len(samples)} good samples in scan window; skips={skips}"
    )

# Regression against the first direct anchored point already recovered.
if p == 100003:
    seed2 = next((row for row in samples if int(row[0]) == 2), None)
    if seed2 is None:
        raise ArithmeticError("tau=2 was not a good sample")
    assert int(seed2[1]) == 11524
    # y_local can differ by sign because local sqrt choices are irrelevant to x.

# Solve x(U) = P(U)/Q(U), deg P=48, Q monic deg 44.
NP = 49
NQFREE = 44
NUNKNOWN = NP + NQFREE  # 93

rows = []
rhs = []
for u0, x0, unused_y in samples:
    up = [F.one()]
    for unused in range(48):
        up.append(up[-1]*u0)
    rows.append(
        up[:49] + [-(x0*up[j]) for j in range(44)]
    )
    rhs.append(x0*(u0**44))

Mall = matrix(F, rows)
if Mall.rank() != NUNKNOWN:
    raise ArithmeticError(
        f"rational interpolation rank {Mall.rank()} != {NUNKNOWN}"
    )

# Pivot columns of transpose identify 93 independent sample rows.
independent_rows = list(Mall.transpose().pivots())
assert len(independent_rows) == NUNKNOWN
Msq = matrix(F, [rows[i] for i in independent_rows])
bsq = vector(F, [rhs[i] for i in independent_rows])
solution = Msq.solve_right(bsq)

P = UR(list(solution[:49]))
Q = UR(list(solution[49:]) + [F.one()])
assert P.degree() == 48
assert Q.degree() == 44
assert P.gcd(Q).degree() == 0

for u0, x0, unused_y in samples:
    assert Q(u0)
    assert P(u0)/Q(u0) == x0

# Denominator must be Z^2 with deg Z=P.O=22.
def polynomial_square_root(poly):
    poly = UR(poly)
    fac = poly.factor()
    unit = F(fac.unit())
    if not unit.is_square():
        raise ArithmeticError("polynomial square has nonsquare unit")
    root = UR(unit.sqrt())
    for factor, exponent in fac:
        if exponent % 2:
            raise ArithmeticError(
                f"odd factor multiplicity {exponent} in expected square"
            )
        root *= factor**(exponent//2)
    assert root*root == poly
    return root


Z = polynomial_square_root(Q)
assert Z.degree() == 22

# Curve equation predicts y^2 with numerator degree 144 and denominator Z^6.
N_y2 = P**3 + A13_poly*P*Q**2 + B13_poly*Q**3
Yabs = polynomial_square_root(N_y2)
assert Yabs.degree() == 72
assert (Z**3).degree() == 66

print(
    "Q8S3MODP_INTERP|"
    f"prime={p}|samples={len(samples)}|"
    "x_degree=48/44|Z_degree=22|"
    "y_square_degree=72/66|"
    f"skips={json.dumps(skips,sort_keys=True,separators=(',',':'))}|"
    "status=PASS_DIRECT_X_INTERPOLATION",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q8-s3-direct-x-modp.v1",
    "status": "PASS_DIRECT_ANCHORED_Q8_S3_X_MODP",
    "prime": int(p),
    "anchored_profile": s3prof,
    "method": (
        "direct II*_E8_1 branch-point birational trace; "
        "x is invariant under local W sqrt sign; no 2-cover halving"
    ),
    "sample_count": len(samples),
    "samples": [
        {"U": int(u0), "x": int(x0), "y_local_up_to_global_sign": int(y0)}
        for u0, x0, y0 in samples
    ],
    "x": {
        "numerator_coefficients_low_to_high": [int(v) for v in P.list()],
        "denominator_coefficients_low_to_high": [int(v) for v in Q.list()],
        "numerator_degree": int(P.degree()),
        "denominator_degree": int(Q.degree()),
        "denominator_root_coefficients_low_to_high": [int(v) for v in Z.list()],
        "denominator_root_degree": int(Z.degree()),
    },
    "y_square": {
        "numerator_square_root_coefficients_low_to_high":
            [int(v) for v in Yabs.list()],
        "numerator_degree": int(Yabs.degree()),
        "denominator_root_cubed_degree": int((Z**3).degree()),
        "global_sign": "unresolved_up_to_negation",
    },
    "skip_counts": skips,
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q8S3MODP_RESULT|x=48/44|Z=22|y_abs=72/66|"
    "status=PASS_DIRECT_ANCHORED_Q8_S3_X_MODP",
    flush=True,
)
