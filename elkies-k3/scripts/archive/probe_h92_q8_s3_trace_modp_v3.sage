#!/usr/bin/env sage -python
"""
Smoke-test the degree-52 Abel-Jacobi trace of the exact q6 third section
through the corrected q8 neighbour onto the D13 child, modulo p.

This is the q8 analogue of the passing q6 trace/halving method.

For one good D13 base value U=tau:
  * H(T)=0 is the degree-52 intersection of exact S3 with the q8 fibre;
  * work in A=GF(p)[T]/H, no degree-52 splitting field;
  * map S3 through the q8 binary-quartic covariant 2-cover;
  * use the 53 functions of L(53 O_D13)
        1,x,...,x^26, y,xy,...,x^25 y
    to recover the sum of all 52 covariant image points;
  * at the q6 IV* old-base value, use the two quartic signs as two candidate
    q8 section origins R+ and R-;
  * for each origin form
        W = TraceCov(S3) - 52*phi(R) = epsilon*2*AJ_R(S3)
    and halve W exactly on the D13 child.

The origin/lattice alignment is intentionally left unresolved in this smoke
test.  A later sampler distinguishes the branch by the pinned D13 target
profile: MW height 47, D13 local correction 1, hence P.O=22.

Run:
  sage -python ~/Downloads/probe_h92_q8_s3_trace_modp.sage --tau 2
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
parser.add_argument("--tau", type=int, default=2)
args = parser.parse_args()

ROOT = locate_repo(args.repo)
GEN = ROOT / "artifacts/generated-results"
LOCAL = ROOT / "artifacts/local/elkies-k3"
CORE = ROOT / "elkies-k3/scripts/elliptic_neighbor_compiler.sage"
Q6 = GEN / "elkies-k3-h92-q6-child-jacobian.json"
Q8 = GEN / "elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json"
BRIDGE = LOCAL / "q6-third-to-q8-bridge.json"

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
R = PolynomialRing(F, "T")
T = R.gen()
K = R.fraction_field()


def modq(value):
    value = QQ(value)
    if ZZ(value.denominator()) % p == 0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {value}")
    return F(ZZ(value.numerator())) / F(ZZ(value.denominator()))


def poly_strings(values):
    return R([modq(value) for value in values])


def rf_from_bridge(entry):
    return K(poly_strings(entry["numerator_coefficients_low_to_high"])) / K(
        poly_strings(entry["denominator_coefficients_low_to_high"])
    )


# q6 child and exact S3 in canonical q6 T-coordinate.
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

# Corrected q8 marked section S on q6 child.
mdata = q8["marking"]["section"]
sx = K(poly_strings(mdata["x_numerator_coefficients_low_to_high"])) / K(
    poly_strings(mdata["x_denominator_coefficients_low_to_high"])
)
sy = K(poly_strings(mdata["y_numerator_coefficients_low_to_high"])) / K(
    poly_strings(mdata["y_denominator_coefficients_low_to_high"])
)
assert sy**2 == sx**3 + K(A6)*sx + K(B6)


def monic_power_root(value, exponent):
    result = R.one()
    for factor, multiplicity in value.factor():
        assert multiplicity % exponent == 0
        result *= factor.monic()**(multiplicity//exponent)
    return result.monic()


nx, dx = R(sx.numerator()), R(sx.denominator())
ny, dy = R(sy.numerator()), R(sy.denominator())
h = monic_power_root(dx, 2)
assert h == monic_power_root(dy, 3)
assert h.degree() == 10

# Fibre factors are serialized as QQ[T] strings.  Parse over QQ first,
# then reduce coefficient-wise; direct parsing inside GF(p)[T] leaves literal
# Rational Field coefficients in the expression tree and triggers coercion
# failures.
QQTR = PolynomialRing(QQ, "T")
def reduce_factor_string(text):
    source = QQTR(text)
    return R([modq(value) for value in source.list()]).monic()

ii = reduce_factor_string(
    next(item for item in q6["finite_fibres"] if item["kodaira"] == "II*")["factor"]
)
iv = reduce_factor_string(
    next(item for item in q6["finite_fibres"] if item["kodaira"] == "IV*")["factor"]
)
assert ii.degree() == iv.degree() == 1
M = (ii**2 * iv**2).monic()

normalizer = (ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
assert (normalizer*h*dy - ny*dx) % nx == 0
p_fun = -sy/sx
alpha = -p_fun/K(h) - K(normalizer)/K(nx)
beta = K(T**2)/K(h)
rho = (normalizer*nx.inverse_mod(M)).mod(M)

kernel_polys = q8["rr"]["kernel_polynomials"]
assert len(kernel_polys) == 2
pairs = []
for entry in kernel_polys:
    sp = R([modq(v) for v in QQTR(entry["s"]).list()])
    tp = R([modq(v) for v in QQTR(entry["t"]).list()])
    Bcoef = K(sp)/K(h)
    Acoef = (
        -K(sp)*p_fun/K(h)
        - K(sp)*K(normalizer)/K(nx)
        + K(sp*rho)
        + K(tp*M)
    )
    pairs.append((Acoef, Bcoef))
(A0, B0), (A1, B1) = pairs

# D13 child coefficients.  Finite-field specializations can acquire
# accidental rational 2-torsion even though the generic D13 MW group has no
# 2-torsion.  Since rational halves then differ by E[2](F_p), preselect a
# specialization with trivial rational 2-torsion before doing the expensive
# degree-52 trace.
A13_poly = [modq(v) for v in q8["child"]["minimal_A_coefficients_low_to_high"]]
B13_poly = [modq(v) for v in q8["child"]["minimal_B_coefficients_low_to_high"]]

def eval_coeff(values, at):
    return sum(value*at**i for i, value in enumerate(values))

XR2 = PolynomialRing(F, "r2")
r2 = XR2.gen()

tau = None
A13 = B13 = E13 = H = None
for offset in range(100):
    candidate = F(args.tau + offset)
    a13 = eval_coeff(A13_poly, candidate)
    b13 = eval_coeff(B13_poly, candidate)
    curve13 = EllipticCurve(F, [0, 0, 0, a13, b13])

    if not curve13.discriminant():
        print(
            f"Q8S3TRACE_SKIP|tau={int(candidate)}|reason=singular_D13_fibre",
            flush=True,
        )
        continue

    two_torsion_roots = (r2**3 + a13*r2 + b13).roots()
    if two_torsion_roots:
        print(
            f"Q8S3TRACE_SKIP|tau={int(candidate)}|"
            f"reason=rational_2_torsion|roots={len(two_torsion_roots)}",
            flush=True,
        )
        continue

    candidate_H = R(U3_num - candidate*U3_den)
    if candidate_H.degree() != 52:
        print(
            f"Q8S3TRACE_SKIP|tau={int(candidate)}|"
            f"reason=degree_drop|degree={candidate_H.degree()}",
            flush=True,
        )
        continue
    candidate_H = candidate_H.monic()

    if candidate_H.gcd(candidate_H.derivative()).degree() != 0:
        print(
            f"Q8S3TRACE_SKIP|tau={int(candidate)}|"
            "reason=non_etale_degree52_fibre",
            flush=True,
        )
        continue

    tau, A13, B13, E13, H = (
        candidate, a13, b13, curve13, candidate_H
    )
    break
else:
    raise ArithmeticError(
        "no good tau with trivial rational 2-torsion in 100 consecutive values"
    )

print(
    f"Q8S3TRACE|prime={p}|tau={int(tau)}|degree=52|"
    "rational_2_torsion=0|stage=setup|status=PASS",
    flush=True,
)

# q8 chord on this fibre.  The marked point is -S=(sx,-sy).
q8_m = -(A1 - tau*A0)/(B1 - tau*B0)
radicand = (
    q8_m**4
    - 6*sx*q8_m**2
    - 8*sy*q8_m
    - 3*sx**2
    - 4*K(A6)
)
quartic, square_factor = squarefree_binary_quartic(radicand, R)
if quartic.degree() != 4:
    raise ArithmeticError(f"q8 quartic degree dropped to {quartic.degree()}")


def native_covariants(quartic):
    base = quartic.base_ring()
    BR = PolynomialRing(base, names=("qx", "qz"))
    qx, qz = BR.gens()
    f = sum(
        BR(quartic[i])*qx**i*qz**(4-i)
        for i in range(5)
    )
    HC = (
        f.derivative(qx,2)*f.derivative(qz,2)
        - f.derivative(qx).derivative(qz)**2
    ) / base(3)
    GC = f.derivative(qx)*HC.derivative(qz) - f.derivative(qz)*HC.derivative(qx)
    I, J = binary_quartic_invariants(quartic)
    I, J = base(I), base(J)
    assert GC**2 == (
        -base(16)/base(3)*HC**3
        + base(256)*I*HC*f**2
        - base(1024)/base(3)*J*f**3
    )
    return f, HC, GC, I, J

covariants = native_covariants(quartic)
fbin, HC, GC, I, J = covariants
stdA, stdB = -F(27)*I, -F(27)*J


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

unit = fourth_sixth_unit(stdA, stdB, A13, B13)


def reduce_mod_H(value):
    value = K(value)
    num = R(value.numerator())
    den = R(value.denominator())
    if den.gcd(H).degree() != 0:
        raise ZeroDivisionError("denominator not invertible modulo degree-52 fibre")
    return (num*den.inverse_mod(H)) % H

# Check that the exact S3 chord agrees with q8_m in the etale algebra.
m3 = (y3 + sy)/(x3 - sx)
if reduce_mod_H(q8_m-m3):
    raise ArithmeticError("q8 chord mismatch on S3 degree-52 fibre")

# Quartic square root at S3: residual quadratic discriminant sqrt.
w3 = (2*x3 + sx - q8_m**2)/square_factor
wA = reduce_mod_H(w3)
if (wA*wA - quartic) % H:
    raise ArithmeticError("S3 q8 quartic square-root mismatch")

# Transport all 52 conjugate points at once into A=F[T]/H.
H_u = R(HC(qx=T, qz=F(1)))
G_u = R(GC(qx=T, qz=F(1)))
quartic_inv = quartic.inverse_mod(H)
xA = (-F(3)/F(4) * H_u * quartic_inv * unit**2) % H
yA = (F(9)/F(32) * G_u * wA * quartic_inv**2 * unit**3) % H
if (yA*yA - xA*xA*xA - A13*xA - B13) % H:
    raise ArithmeticError("S3 covariant images miss D13 child in etale algebra")
print("Q8S3TRACE|stage=etale_transport|status=PASS", flush=True)

# Trace the 52 covariant image points using L(53 O_D13).
one = R.one()
xpowers = [one]
for _ in range(26):
    xpowers.append((xpowers[-1]*xA) % H)
columns = list(xpowers)  # x^0..x^26, 27 columns
for exponent in range(26):
    columns.append((yA*xpowers[exponent]) % H)  # y*x^0..25, 26
assert len(columns) == 53

Eval = matrix(F, 52, 53, lambda row, col: columns[col][row])
ker = Eval.right_kernel().basis_matrix()
if ker.nrows() != 1:
    raise ArithmeticError(f"L(53O) trace kernel dimension {ker.nrows()}")
rel = ker[0]

XR = PolynomialRing(F, "X")
Xv = XR.gen()
Afun = sum(rel[i]*Xv**i for i in range(27))
Bfun = sum(rel[27+i]*Xv**i for i in range(26))
Rint = Afun**2 - (Xv**3 + A13*Xv + B13)*Bfun**2
if Rint.degree() != 53:
    raise ArithmeticError(f"residual intersection degree {Rint.degree()}, expected 53")
root_sum = -Rint[52]/Rint[53]

# Newton sums for Trace_A/F(xA).
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

ps = newton_power_sums(H)
trace_x = sum(xA[i]*ps[i] for i in range(52))
xQ = root_sum - trace_x
bQ = Bfun(xQ)
if not bQ:
    raise ArithmeticError("trace residual has B(x_Q)=0")
yQ = -Afun(xQ)/bQ
Qres = E13(xQ, yQ)
traceCov = -Qres
print("Q8S3TRACE|stage=covariant_trace|status=PASS", flush=True)

# Two candidate q8 section origins from the two IV* component branches.
t_iv = -iv[0]/iv[1]
qiv = quartic(t_iv)
if not qiv:
    raise ArithmeticError("IV* old-base point is a quartic branch point; need local limit origin")
if not qiv.is_square():
    raise ArithmeticError("IV* quartic value is not square over base field")
wroot = qiv.sqrt()


def map_quartic_point(tvalue, wvalue):
    fv = F(fbin(qx=F(tvalue), qz=F(1)))
    hv = F(HC(qx=F(tvalue), qz=F(1)))
    gv = F(GC(qx=F(tvalue), qz=F(1)))
    if not fv or F(wvalue)**2 != fv:
        raise ArithmeticError("candidate q8 origin is not a nonbranch quartic point")
    rx = -F(3)/F(4)*hv/fv
    ry = F(9)/F(32)*gv*F(wvalue)/fv**2
    return E13(unit**2*rx, unit**3*ry)

Qplus = map_quartic_point(t_iv, wroot)
Qminus = map_quartic_point(t_iv, -wroot)
if Qplus == Qminus:
    raise ArithmeticError("IV* origin candidates collapsed")


def rational_halves(E, target):
    if target.is_zero():
        raise ArithmeticError("unexpected zero doubled AJ target")
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


def ptxt(P):
    if P.is_zero():
        return "O"
    x, y = P.xy()
    return f"{int(x)},{int(y)}"

answers = []
for label, origin in (("plus", Qplus), ("minus", Qminus)):
    doubled = traceCov - 52*origin
    halves = rational_halves(E13, doubled)
    if len(halves) != 1:
        raise ArithmeticError(
            f"origin {label}: doubled AJ has {len(halves)} rational halves "
            "despite E[2](F_p)=0"
        )
    AJ = halves[0]
    assert 2*AJ == doubled
    answers.append((label, origin, AJ))

print(
    "Q8S3TRACE_SPECIALIZATION|"
    f"prime={p}|tau={int(tau)}|tIV={int(t_iv)}|"
    f"traceCov={ptxt(traceCov)}|"
    f"Qplus={ptxt(Qplus)}|AJplus={ptxt(answers[0][2])}|"
    f"Qminus={ptxt(Qminus)}|AJminus={ptxt(answers[1][2])}|"
    "status=PASS",
    flush=True,
)
print(
    "Q8S3TRACE_RESULT|degree52=PASS|etale_transport=PASS|"
    "L53_trace=PASS|rational_2_torsion=0|two_cover_halving=PASS|"
    "origin=UNRESOLVED_IVSTAR_PLUS_MINUS|status=PASS_SMOKE_TEST",
    flush=True,
)
