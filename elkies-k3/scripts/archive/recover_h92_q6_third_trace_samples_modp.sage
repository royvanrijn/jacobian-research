#!/usr/bin/env sage -python
"""
Recover the missing q6 third section modulo 100003 by fast trace sampling.

This avoids the expensive generic 50x51 kernel over GF(p)(T).

For each ordinary child-base value tau in GF(p):

  1. H_tau(u)=0 is the degree-50 fibre of the explicit old section C=-P2.
  2. Work directly in A_tau=GF(p)[u]/H_tau.
  3. Transport C to the certified q6 child in that algebra.
  4. Evaluate the 51 basis functions of L(51 O):
       1,x,...,x^25, y,x*y,...,x^24*y.
     A 50x51 matrix over GF(p) has a one-dimensional kernel.
  5. Recover the one residual point Q from the kernel relation using
     Newton sums of H_tau, hence Trace(C)=-Q.
  6. Orient the two chord-pole points intrinsically and form
       S3 = Trace(C) - 25*QO - 24*QA.

After collecting good samples, recover
    x(T)=X(T)/Z(T)^2,   deg X<=46, deg Z<=21,
then interpolate
    y(T)=Y(T)/Z(T)^3,   deg Y<=69.

The final polynomial Weierstrass identity is verified identically over GF(p).

Run:
  sage -python ~/Downloads/recover_h92_q6_third_trace_samples_modp.sage

Useful options:
  --samples 100
  --start 1
  --output /path/to/result.json
"""

import argparse
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector


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
parser.add_argument("--samples", type=int, default=100)
parser.add_argument("--start", type=int, default=1)
parser.add_argument("--scan-limit", type=int, default=1000)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
P1FILE = GEN / "elkies-k3-h92-p1-lift.json"
P2FILE = GEN / "elkies-k3-h92-p2-hensel-100003-p1024.json"
RRFILE = GEN / "elkies-k3-h92-q6-global-rr.json"
CHILDFILE = GEN / "elkies-k3-h92-q6-child-jacobian.json"
OUTPUT = (
    args.output.resolve()
    if args.output
    else ROOT / "artifacts/local/elkies-k3/q6-third-trace-samples-mod-100003.json"
)

for path in (CORE, ANCHOR, H92, P1FILE, P2FILE, RRFILE, CHILDFILE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

scope = {}
exec(compile(CORE.read_text(), str(CORE), "exec"), scope)
squarefree_binary_quartic = scope["squarefree_binary_quartic"]
binary_quartic_invariants = scope["binary_quartic_invariants"]

p = ZZ(args.prime)
if p != 100003:
    raise ValueError("the pinned Hensel -P2 input is tied to prime 100003")
F = GF(p)

R = PolynomialRing(F, "u")
u = R.gen()
K = R.fraction_field()

p1data = json.loads(P1FILE.read_text())
p2data = json.loads(P2FILE.read_text())
rr = json.loads(RRFILE.read_text())
child = json.loads(CHILDFILE.read_text())

assert p1data["status"] == "PASS_EXACT_H92_P1"
assert p2data["schema"] == "elkies-k3.h92-p2-hensel-lift.v1"
assert p2data["complete"]
assert rr["status"] == "PASS_EXACT_GLOBAL_RR_KERNEL"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"


def modq(value):
    value = QQ(value)
    if ZZ(value.denominator()) % p == 0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {value}")
    return F(ZZ(value.numerator())) / F(ZZ(value.denominator()))


def poly(values):
    return R([modq(value) for value in values])


xp = K(poly(p1data["x_entrance_base"]["numerator_coefficients"])) / K(
    poly(p1data["x_entrance_base"]["denominator_coefficients"])
)
yp = K(poly(p1data["y_entrance_base"]["numerator_coefficients"])) / K(
    poly(p1data["y_entrance_base"]["denominator_coefficients"])
)


def reciprocal(values):
    answer = K(0)
    for index, value in enumerate(values):
        answer += K(modq(value)) / K(u**index)
    return answer


z2 = reciprocal(p2data["Z"])
xC = reciprocal(p2data["X"]) / z2**2
yC = reciprocal(p2data["Y"]) / z2**3

anchor = SourceFileLoader("q6_trace_samples_anchor", str(ANCHOR)).load_module()
r0, s0 = anchor.EXPECTED_H92
unused_ring, formulas = anchor.parse_h92(H92)
A1q, Aq, B1q, Bq, B2q = tuple(QQ(value(r0, s0)) for value in formulas)
A1, A, B1old, Bold, B2old = map(modq, (A1q, Aq, B1q, Bq, B2q))
old_a = K(A1) / u**3 + K(A) / u**4
old_b = K(B1old) / u**5 + K(Bold) / u**6 + K(B2old) / u**7
assert yp**2 == xp**3 + old_a*xp + old_b
assert yC**2 == xC**3 + old_a*xC + old_b

h = poly(p1data["structured_denominator"]["Z4_coefficients"])


def coefficient_pair(entry):
    Ap = poly(entry["A_coefficients_low_to_high"])
    Bp = poly(entry["B_coefficients_low_to_high"])
    return K(Ap) / K(h**2), K(Bp) / K(h)


(a0, b0), (a1, b1) = tuple(
    coefficient_pair(entry) for entry in rr["kernel"]["sections"]
)

mC = (yC - yp) / (xC - xp)
TC = (a1 + b1*mC) / (a0 + b0*mC)
assert max(TC.numerator().degree(), TC.denominator().degree()) == 50

childA_coeff = [
    modq(value)
    for value in child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
]
childB_coeff = [
    modq(value)
    for value in child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
]


def child_coeff(coeffs, tau):
    return sum(value*tau**index for index, value in enumerate(coeffs))


def native_covariants(quartic):
    base = quartic.base_ring()
    BR = PolynomialRing(base, names=("qx", "qz"))
    qx, qz = BR.gens()
    f = sum(
        BR(quartic[index]) * qx**index * qz**(4-index)
        for index in range(5)
    )
    Hcov = (
        f.derivative(qx, 2)*f.derivative(qz, 2)
        - f.derivative(qx).derivative(qz)**2
    ) / base(3)
    Gcov = f.derivative(qx)*Hcov.derivative(qz) - f.derivative(qz)*Hcov.derivative(qx)
    I, J = binary_quartic_invariants(quartic)
    I, J = base(I), base(J)
    assert Gcov**2 == (
        -base(16)/base(3)*Hcov**3
        + base(256)*I*Hcov*f**2
        - base(1024)/base(3)*J*f**3
    )
    return f, Hcov, Gcov, I, J


def fourth_sixth_unit(std_a, std_b, target_a, target_b):
    if not std_a or not std_b or not target_a or not target_b:
        raise ArithmeticError("degenerate Jacobian coefficient")
    CR = PolynomialRing(F, "c")
    c = CR.gen()
    roots = [
        root for root, multiplicity in (c**4-target_a/std_a).roots()
        if multiplicity == 1 and std_b*root**6 == target_b
    ]
    if not roots:
        raise ArithmeticError("no fourth/sixth minimizing unit")
    return roots[0]


def reduce_mod_H(value, H):
    value = K(value)
    num = R(value.numerator())
    den = R(value.denominator())
    if den.gcd(H).degree() != 0:
        raise ZeroDivisionError("denominator not invertible modulo H")
    return (num*den.inverse_mod(H)) % H


def removable_eval(value, at):
    value = K(value)
    num = R(value.numerator())
    den = R(value.denominator())
    linear = u - F(at)
    while num(F(at)) == 0 and den(F(at)) == 0:
        num, rn = num.quo_rem(linear)
        den, rd = den.quo_rem(linear)
        if rn or rd:
            raise ArithmeticError("common-factor cancellation failed")
    if den(F(at)) == 0:
        raise ZeroDivisionError("genuine pole")
    return num(F(at))/den(F(at))


def oriented_pole_points(tau, quartic, covariants, square_factor, q6_m,
                         unit, child_a, child_b):
    B0 = poly(rr["kernel"]["sections"][0]["B_coefficients_low_to_high"])
    B1p = poly(rr["kernel"]["sections"][1]["B_coefficients_low_to_high"])
    line = B1p - tau*B0
    quotient, remainder = line.quo_rem(u**3)
    if remainder or quotient.degree() != 1:
        raise ArithmeticError("pole branch not linear")
    uO = -quotient[0]/quotient[1]
    qvalue = quartic(uO)
    if not qvalue:
        raise ArithmeticError("pole branch is a quartic branch point")
    wO = removable_eval(q6_m**2/square_factor, uO)
    if wO**2 != qvalue:
        raise ArithmeticError("old-O sign failed")

    f, Hcov, Gcov, I, J = covariants
    E = EllipticCurve(F, [0, 0, 0, child_a, child_b])

    def map_point(w):
        fv = F(f(qx=F(uO), qz=F(1)))
        hv = F(Hcov(qx=F(uO), qz=F(1)))
        gv = F(Gcov(qx=F(uO), qz=F(1)))
        rx = -F(3)/F(4)*hv/fv
        ry = F(9)/F(32)*gv*F(w)/fv**2
        return E(unit**2*rx, unit**3*ry)

    QO = map_point(wO)
    QA = map_point(-wO)
    assert QO != QA
    return QO, QA


def newton_power_sums(H):
    n = H.degree()
    assert H[n] == 1
    sums = [F(n)]
    for k in range(1, n):
        total = F(k)*H[n-k]
        for j in range(1, k):
            total += H[n-j]*sums[k-j]
        sums.append(-total)
    return sums


def section_at_tau(tau):
    child_a = child_coeff(childA_coeff, tau)
    child_b = child_coeff(childB_coeff, tau)
    E = EllipticCurve(F, [0, 0, 0, child_a, child_b])
    if not E.discriminant():
        raise ArithmeticError("singular child fibre")

    H = R(TC.numerator() - tau*TC.denominator())
    if H.degree() != 50:
        raise ArithmeticError("degree-50 fibre dropped")
    H = H.monic()
    if H.gcd(H.derivative()).degree() != 0:
        raise ArithmeticError("degree-50 fibre not etale")

    q6_m = (a1 - tau*a0) / (tau*b0 - b1)
    radicand = q6_m**4 - 6*xp*q6_m**2 + 8*yp*q6_m - 3*xp**2 - 4*old_a
    quartic, square_factor = squarefree_binary_quartic(radicand, R)
    if quartic.degree() != 4:
        raise ArithmeticError("quartic degree dropped")

    covariants = native_covariants(quartic)
    f, Hcov, Gcov, I, J = covariants
    std_a, std_b = -F(27)*I, -F(27)*J
    unit = fourth_sixth_unit(std_a, std_b, child_a, child_b)

    # C in the etale algebra.
    if reduce_mod_H(q6_m-mC, H):
        raise ArithmeticError("q6 chord mismatch in etale algebra")
    wC = (2*xC + xp - mC**2)/square_factor
    wA = reduce_mod_H(wC, H)
    if (wA*wA - quartic) % H:
        raise ArithmeticError("quartic square-root mismatch")

    H_u = R(Hcov(qx=u, qz=F(1)))
    G_u = R(Gcov(qx=u, qz=F(1)))
    quartic_inv = quartic.inverse_mod(H)

    xA = (
        -F(3)/F(4) * H_u * quartic_inv * unit**2
    ) % H
    yA = (
        F(9)/F(32) * G_u * wA * quartic_inv**2 * unit**3
    ) % H
    if (yA*yA - xA*xA*xA - child_a*xA - child_b) % H:
        raise ArithmeticError("child C point missed curve in etale algebra")

    # Unique f=A(X)+YB(X) in L(51O) vanishing at all 50 C-points.
    one = R(1)
    xpowers = [one]
    for _ in range(25):
        xpowers.append((xpowers[-1]*xA) % H)
    columns = list(xpowers)
    for exponent in range(25):
        columns.append((yA*xpowers[exponent]) % H)

    M = matrix(F, 50, 51, lambda row, col: columns[col][row])
    ker = M.right_kernel().basis_matrix()
    if ker.nrows() != 1:
        raise ArithmeticError(f"trace RR kernel dimension {ker.nrows()}")
    rel = ker[0]

    XR = PolynomialRing(F, "X")
    Xv = XR.gen()
    Afun = sum(rel[index]*Xv**index for index in range(26))
    Bfun = sum(rel[26+index]*Xv**index for index in range(25))
    Rint = Afun**2 - (Xv**3 + child_a*Xv + child_b)*Bfun**2
    if Rint.degree() != 51:
        raise ArithmeticError("residual intersection degree is not 51")

    root_sum = -Rint[50]/Rint[51]
    power_sums = newton_power_sums(H)
    trace_x = sum(xA[index]*power_sums[index] for index in range(50))
    xQ = root_sum - trace_x
    bq = Bfun(xQ)
    if not bq:
        raise ArithmeticError("residual point B(X)=0")
    yQ = -Afun(xQ)/bq
    Q = E(xQ, yQ)
    traceC = -Q

    QO, QA = oriented_pole_points(
        tau, quartic, covariants, square_factor, q6_m, unit, child_a, child_b
    )
    S3 = traceC - 25*QO - 24*QA
    if S3.is_zero():
        raise ArithmeticError("S3 specialized to zero")
    sx, sy = S3.xy()
    return F(sx), F(sy)


samples = []
attempted = 0
for integer in range(args.start, args.start + args.scan_limit):
    if len(samples) >= args.samples:
        break
    attempted += 1
    tau = F(integer)
    try:
        sx, sy = section_at_tau(tau)
    except Exception as error:
        print(
            f"Q6TRACESAMPLE_SKIP|T={integer}|"
            f"reason={type(error).__name__}:{error}",
            flush=True,
        )
        continue

    if integer == 2:
        assert (sx, sy) == (F(63169), F(42300))

    samples.append((tau, sx, sy))
    if len(samples) <= 5 or len(samples) % 10 == 0:
        print(
            f"Q6TRACESAMPLE|count={len(samples)}|T={integer}|"
            f"x={int(sx)}|y={int(sy)}|status=PASS",
            flush=True,
        )

if len(samples) < args.samples:
    raise RuntimeError(
        f"only {len(samples)} good samples after {attempted} attempts"
    )

print(
    f"Q6TRACESAMPLE|good={len(samples)}|attempted={attempted}|stage=collection|status=PASS",
    flush=True,
)

# Recover x=N/D with deg N<=46, deg D<=42.
num_deg, den_deg = 46, 42
interp_rows = []
for tau, sx, unused_sy in samples:
    interp_rows.append(
        [tau**j for j in range(num_deg+1)]
        + [-sx*tau**j for j in range(den_deg+1)]
    )
IM = matrix(F, interp_rows)
IK = IM.right_kernel().basis_matrix()
if IK.nrows() != 1:
    raise ArithmeticError(f"x rational-interpolation kernel dimension {IK.nrows()}")
rv = IK[0]

TR = PolynomialRing(F, "T")
T = TR.gen()
Nx = TR(list(rv[:num_deg+1]))
Dx = TR(list(rv[num_deg+1:]))

if not Dx:
    raise ArithmeticError("interpolated x denominator vanished")
scale = Dx.leading_coefficient()
Nx /= scale
Dx /= scale

# Trim is automatic in polynomial construction; verify every sample.
for tau, sx, unused_sy in samples:
    if not Dx(tau) or Nx(tau)/Dx(tau) != sx:
        raise ArithmeticError("x interpolation failed a retained sample")

if not Dx.is_square():
    raise ArithmeticError(
        f"x denominator degree {Dx.degree()} is not a square"
    )
Z = Dx.sqrt()
if Z.leading_coefficient() != 1:
    zscale = Z.leading_coefficient()
    Z /= zscale
    Nx /= zscale**2
    Dx = Z**2
assert Dx == Z**2

# Interpolate Y=y*Z^3, degree <=69.
usable_y = [
    (tau, sy*Z(tau)**3)
    for tau, unused_sx, sy in samples
    if Z(tau) != 0
]
if len(usable_y) < 70:
    raise ArithmeticError("fewer than 70 non-pole y samples")

Ydeg = 69
VM = matrix(F, [[tau**j for j in range(Ydeg+1)] for tau, value in usable_y[:70]])
rhs = vector(F, [value for tau, value in usable_y[:70]])
ycoeff = VM.solve_right(rhs)
Y = TR(list(ycoeff))

for tau, value in usable_y:
    if Y(tau) != value:
        raise ArithmeticError("Y interpolation failed a retained sample")

# Global polynomial child identity.
Apoly = TR(childA_coeff)
Bpoly = TR(childB_coeff)
assert Y**2 == Nx**3 + Apoly*Nx*Z**4 + Bpoly*Z**6

print(
    "Q6TRACESAMPLE_PROFILE|"
    f"Z_degree={Z.degree()}|X_degree={Nx.degree()}|Y_degree={Y.degree()}|"
    f"x={Nx.degree()}/{Dx.degree()}|"
    f"good_samples={len(samples)}|status=PASS",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q6-third-trace-samples-modp.v1",
    "status": "PASS_MODULAR_Q6_THIRD_TRACE_SECTION",
    "prime": int(p),
    "method": {
        "trace_multisection_degree": 50,
        "trace_rr_basis": "L(51*O)",
        "per_specialization_matrix": [50, 51],
        "orientation": "local m->infinity old-O asymptotic",
        "formula": "S3=Trace(C)-25*QO-24*QA",
        "x_interpolation_bounds": [num_deg, den_deg],
        "y_recovery": "interpolate polynomial Y=y*Z^3 after recovering Z",
    },
    "sample_count": len(samples),
    "attempted": attempted,
    "check_T2": {"x": 63169, "y": 42300},
    "profile": {
        "Z_degree": int(Z.degree()),
        "X_degree": int(Nx.degree()),
        "Y_degree": int(Y.degree()),
        "x_denominator_degree": int(Dx.degree()),
    },
    "Z": [int(value) for value in Z.list()],
    "X": [int(value) for value in Nx.list()],
    "Y": [int(value) for value in Y.list()],
    "samples": [
        {"T": int(tau), "x": int(sx), "y": int(sy)}
        for tau, sx, sy in samples
    ],
    "boundary": (
        "Exact over GF(100003)(T) by interpolation plus identity verification. "
        "Characteristic-zero p-adic lifting and exact Shioda/component audits "
        "remain downstream."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}")
print("Q6TRACESAMPLE_RESULT|status=PASS_MODULAR_Q6_THIRD_TRACE_SECTION", flush=True)
