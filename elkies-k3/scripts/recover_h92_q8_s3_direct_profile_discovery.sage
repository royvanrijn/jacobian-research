#!/usr/bin/env sage -python
"""
Recover the direct anchored q8 Abel-Jacobi x-coordinate of H92 S3 over GF(p).

Prerequisites:
  * PASS_EXACT_Q8_EQUATION_NS_DIVISOR
  * PASS_EXACT_D13_BRANCH_ANCHOR
  * PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52

The q8 zero is the actual II*_E8_1 branch section.  For many U=tau values,
map the degree-52 S3 divisor birationally through the anchored branch-point
quartic -> D13 map, sum the 52 points with L(53 O), and retain AJ_x(tau).

No degree profile is assumed.  A training set determines the minimal rational
function x(U)=N(U)/D(U) by polynomial rational reconstruction; held-out direct
traces certify it.  The result is then required to have D=Z^2 and to make the
Weierstrass right-hand side an exact square Y^2.

No covariant 2-cover, IV* origin subtraction, halving, or historical q8-zero
profile is used.
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
parser.add_argument("--samples", type=int, default=135)
parser.add_argument("--scan-limit", type=int, default=800)
parser.add_argument("--holdout", type=int, default=20)
parser.add_argument("--output", type=Path)
args = parser.parse_args()
if args.samples < 60:
    raise ValueError("need at least 60 direct samples for profile discovery")
if args.holdout < 8 or args.holdout >= args.samples // 2:
    raise ValueError("holdout must be at least 8 and less than half the samples")

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
SEED = LOCAL / "q8-s3-direct-anchor-trace-mod-100003-tau-2.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q8-s3-direct-x-mod-{args.prime}.json"
)
SAMPLE_OUTPUT = LOCAL / f"q8-s3-direct-samples-mod-{args.prime}.json"

for path in (CORE, Q6, Q8, BRIDGE, ANCHOR, EQNS):
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

assert q6["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"
assert q8["status"] == "PASS_EXACT_CORRECTED_Q8_D13_CHILD"
assert bridge["status"] == "PASS_EXACT_Q6_THIRD_TO_Q8_DEGREE52"
assert anchor["status"] == "PASS_EXACT_D13_BRANCH_ANCHOR"
assert eqns["status"] == "PASS_EXACT_Q8_EQUATION_NS_DIVISOR"
assert eqns["q8_equation_fibre"]["S3_degree"] == 52
assert eqns["q8_equation_fibre"]["root_data"] == [13, 312, 4]

# Deliberately do not import a historical D13 MW/profile assertion here.
# The direct II*_E8_1 branch-point construction will discover its own global
# rational section profile from its specialization values.
print(
    "Q8S3DIRECT_SETUP|degree=52|zero=II*_E8_1|"
    "profile=DISCOVER_FROM_DIRECT_SAMPLES|status=PASS",
    flush=True,
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


# ---------------------------------------------------------------------------
# Direct sample collection, resumable.
# ---------------------------------------------------------------------------

samples = []
skips = {}
next_tau = int(args.start)

if SAMPLE_OUTPUT.exists():
    try:
        saved = json.loads(SAMPLE_OUTPUT.read_text())
        if (
            saved.get("schema") == "elkies-k3.h92-q8-s3-direct-samples-modp.v1"
            and int(saved.get("prime", -1)) == int(p)
            and int(saved.get("start", -1)) == int(args.start)
        ):
            samples = [
                (F(row["U"]), F(row["x"]), F(row["y_local_up_to_sign"]))
                for row in saved.get("samples", [])
            ]
            skips = dict(saved.get("skip_counts", {}))
            next_tau = int(saved.get("next_tau", args.start))
            print(
                f"Q8S3DIRECT_RESUME|good={len(samples)}|next_tau={next_tau}|"
                "status=PASS",
                flush=True,
            )
    except Exception as exc:
        print(
            f"Q8S3DIRECT_RESUME|status=IGNORE_INVALID|reason={type(exc).__name__}",
            flush=True,
        )
        samples = []
        skips = {}
        next_tau = int(args.start)


def save_checkpoint():
    payload = {
        "schema": "elkies-k3.h92-q8-s3-direct-samples-modp.v1",
        "status": "PARTIAL_DIRECT_ANCHORED_Q8_S3_SAMPLES",
        "prime": int(p),
        "start": int(args.start),
        "next_tau": int(next_tau),
        "samples": [
            {
                "U": int(u0),
                "x": int(x0),
                "y_local_up_to_sign": int(y0),
            }
            for u0, x0, y0 in samples
        ],
        "skip_counts": skips,
    }
    SAMPLE_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    SAMPLE_OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


stop = int(args.start) + int(args.scan_limit)
tau = max(next_tau, int(args.start))

while tau < stop and len(samples) < args.samples:
    try:
        ax, ay = direct_x_at(tau)
        samples.append((F(tau), ax, ay))
        if tau == 2 or len(samples) <= 5 or len(samples) % 10 == 0:
            print(
                f"Q8S3DIRECT_SAMPLE|good={len(samples)}|tau={tau}|"
                f"x={int(ax)}|y_local={int(ay)}|status=PASS",
                flush=True,
            )
    except Exception as exc:
        key = type(exc).__name__
        skips[key] = skips.get(key, 0) + 1
    tau += 1
    next_tau = tau
    if len(samples) % 10 == 0:
        save_checkpoint()

save_checkpoint()

if len(samples) < args.samples:
    raise ArithmeticError(
        f"only {len(samples)} good direct samples in scan window; "
        f"next_tau={next_tau}; skips={skips}"
    )

# Regression against the first direct anchored point already recovered.
if p == 100003:
    seed2 = next((row for row in samples if int(row[0]) == 2), None)
    if seed2 is None:
        raise ArithmeticError("tau=2 was not retained as a good direct sample")
    assert int(seed2[1]) == 11524

print(
    f"Q8S3DIRECT_COLLECTION|good={len(samples)}|holdout={args.holdout}|"
    f"skips={json.dumps(skips,sort_keys=True,separators=(',',':'))}|status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Profile-free rational reconstruction of x(U).
#
# Let M(U)=prod(U-u_i) on the training set and let f(U) interpolate the x_i.
# We seek N/D with N == D*f mod M.  Extended Euclid on (M,f) enumerates the
# rational-reconstruction convergents.  No numerator/denominator degree is
# supplied in advance.  Held-out direct traces select the genuine candidate.
# ---------------------------------------------------------------------------

train_count = len(samples) - int(args.holdout)
train = samples[:train_count]
holdout = samples[train_count:]
assert len(holdout) == int(args.holdout)

Mtrain = UR.one()
for u0, unused_x0, unused_y0 in train:
    Mtrain *= (U - u0)
Mtrain = UR(Mtrain.monic())

# Direct Lagrange interpolation, avoiding dependence on a Sage helper API.
Finterp = UR.zero()
for u0, x0, unused_y0 in train:
    quotient, remainder = Mtrain.quo_rem(U - u0)
    assert not remainder
    denom = quotient(u0)
    assert denom
    Finterp += x0 * quotient / denom
Finterp = UR(Finterp)
assert Finterp.degree() < Mtrain.degree()
for u0, x0, unused_y0 in train:
    assert Finterp(u0) == x0


def normalize_candidate(N, D):
    N = UR(N)
    D = UR(D)
    if not D:
        return None
    g = N.gcd(D)
    if g:
        N //= g
        D //= g
    lc = D.leading_coefficient()
    N /= lc
    D /= lc
    N = UR(N)
    D = UR(D)
    assert D.is_monic()
    return N, D


def replays(candidate, rows):
    N, D = candidate
    for u0, x0, unused_y0 in rows:
        dv = D(u0)
        if not dv or N(u0) != x0 * dv:
            return False
    return True


# r_i = s_i*M + t_i*f, hence r_i == t_i*f mod M.
r0, r1 = Mtrain, Finterp
t0, t1 = UR.zero(), UR.one()
valid = {}
step = 0

while r1:
    candidate = normalize_candidate(r1, t1)
    if candidate is not None:
        N, D = candidate
        # Only relations whose total degree is below the interpolation count
        # can be uniquely meaningful.
        if N.degree() + D.degree() + 1 < train_count:
            if replays(candidate, train) and replays(candidate, holdout):
                key = (
                    tuple(int(v) for v in N.list()),
                    tuple(int(v) for v in D.list()),
                )
                valid[key] = (N, D, step)
                print(
                    "Q8S3DIRECT_RR_CANDIDATE|"
                    f"step={step}|x_degree={N.degree()}/{D.degree()}|"
                    f"total={N.degree()+D.degree()}|holdout=PASS",
                    flush=True,
                )

    q, r2 = r0.quo_rem(r1)
    t2 = t0 - q*t1
    r0, r1 = r1, r2
    t0, t1 = t1, t2
    step += 1

if not valid:
    raise ArithmeticError(
        "no rational-reconstruction convergent predicts the held-out direct samples; "
        "collect more samples or change prime"
    )

candidates = list(valid.values())
candidates.sort(
    key=lambda item: (
        item[0].degree() + item[1].degree(),
        max(item[0].degree(), item[1].degree()),
        item[1].degree(),
        item[0].degree(),
    )
)
P, Q, rr_step = candidates[0]

# The minimal candidate must replay every direct specialization.
assert replays((P, Q), samples)

# Reject ambiguity at the same minimal complexity.
best_key = (
    P.degree() + Q.degree(),
    max(P.degree(), Q.degree()),
    Q.degree(),
    P.degree(),
)
ties = [
    (N, D, st)
    for N, D, st in candidates
    if (
        N.degree() + D.degree(),
        max(N.degree(), D.degree()),
        D.degree(),
        N.degree(),
    ) == best_key
]
if len(ties) != 1:
    raise ArithmeticError(
        f"ambiguous minimal rational reconstruction: {len(ties)} candidates at {best_key}"
    )

print(
    "Q8S3DIRECT_PROFILE_DISCOVERY|"
    f"train={train_count}|holdout={len(holdout)}|"
    f"rr_step={rr_step}|x={P.degree()}/{Q.degree()}|"
    f"candidate_count={len(candidates)}|status=PASS_UNIQUE_HELDOUT_PROFILE",
    flush=True,
)

# ---------------------------------------------------------------------------
# Intrinsic Weierstrass structure: denominator square and y-square.
# ---------------------------------------------------------------------------

def polynomial_square_root(poly):
    poly = UR(poly)
    if not poly:
        return UR.zero()
    fac = poly.factor()
    unit = F(fac.unit())
    if not unit.is_square():
        return None
    root = UR(unit.sqrt())
    for factor, exponent in fac:
        if exponent % 2:
            return None
        root *= factor**(exponent//2)
    if root*root != poly:
        return None
    return root


Z = polynomial_square_root(Q)
if Z is None:
    raise ArithmeticError(
        f"discovered x denominator degree {Q.degree()} is not a polynomial square"
    )

N_y2 = P**3 + A13_poly*P*Q**2 + B13_poly*Q**3
Yabs = polynomial_square_root(N_y2)
if Yabs is None:
    raise ArithmeticError(
        "discovered x(U) replays all direct samples but curve RHS numerator "
        "is not a polynomial square"
    )

# Every local direct y may differ by sign because scale.sqrt() is chosen
# independently at each specialization.  Its absolute section must agree.
for u0, unused_x0, y0 in samples:
    zv = Z(u0)
    assert zv
    predicted = Yabs(u0) / (zv**3)
    if predicted != y0 and -predicted != y0:
        raise ArithmeticError(
            f"global |y| fails direct specialization U={int(u0)}"
        )

x_num_deg = int(P.degree())
x_den_deg = int(Q.degree())
z_deg = int(Z.degree())
y_num_deg = int(Yabs.degree())
y_den_deg = int((Z**3).degree())

print(
    "Q8S3DIRECT_STRUCTURE|"
    f"x={x_num_deg}/{x_den_deg}|Z={z_deg}|"
    f"y_abs={y_num_deg}/{y_den_deg}|"
    "denominator_square=1|curve_rhs_square=1|"
    "local_y_up_to_sign=PASS|status=PASS",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q8-s3-direct-profile-discovery-modp.v1",
    "status": "PASS_DIRECT_ANCHORED_Q8_S3_PROFILE_DISCOVERY",
    "prime": int(p),
    "method": (
        "II*_E8_1 branch-point birational trace + L(53O); "
        "profile-free polynomial rational reconstruction with held-out replay; "
        "no covariant 2-cover, origin subtraction, halving, or assumed MW profile"
    ),
    "sample_count": len(samples),
    "training_count": train_count,
    "holdout_count": len(holdout),
    "skip_counts": skips,
    "known_direct_regression": {"U": 2, "x": 11524},
    "x": {
        "numerator_coefficients_low_to_high": [int(v) for v in P.list()],
        "denominator_coefficients_low_to_high": [int(v) for v in Q.list()],
        "numerator_degree": x_num_deg,
        "denominator_degree": x_den_deg,
        "rational_reconstruction_step": int(rr_step),
    },
    "weierstrass_structure": {
        "denominator_is_square": True,
        "denominator_root_coefficients_low_to_high": [int(v) for v in Z.list()],
        "denominator_root_degree": z_deg,
        "curve_rhs_numerator_is_square": True,
        "y_abs_numerator_coefficients_low_to_high": [int(v) for v in Yabs.list()],
        "y_abs_numerator_degree": y_num_deg,
        "y_denominator_degree": y_den_deg,
        "local_y_replay": "up_to_independent_specialization_sign",
    },
    "samples": [
        {"U": int(u0), "x": int(x0), "y_local_up_to_sign": int(y0)}
        for u0, x0, y0 in samples
    ],
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

# Mark checkpoint complete too.
checkpoint = json.loads(SAMPLE_OUTPUT.read_text())
checkpoint["status"] = "PASS_DIRECT_ANCHORED_Q8_S3_SAMPLES"
checkpoint["final_profile_output"] = str(OUTPUT.relative_to(ROOT))
SAMPLE_OUTPUT.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")

print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q8S3DIRECT_RESULT|"
    f"x={x_num_deg}/{x_den_deg}|Z={z_deg}|"
    f"y_abs={y_num_deg}/{y_den_deg}|"
    "status=PASS_DIRECT_ANCHORED_Q8_S3_PROFILE_DISCOVERY",
    flush=True,
)
