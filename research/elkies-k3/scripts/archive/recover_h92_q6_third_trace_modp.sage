#!/usr/bin/env sage -python
"""
Recover the missing q6 third section directly over GF(p)(T) by divisor trace.

This is the generic version of the passing specialization smoke test.

Input geometry:
    C = reconstructed old section -P2,
    deg(C -> q6 base) = 50,
    S3 = Trace_q6(C) - 25*QO - 24*QA

where QO and QA are the binary-quartic images of old O and A=-P1,
oriented intrinsically by the m->infinity local asymptotic.

Instead of sampling ~90 T-values, work over K=GF(p)(T).

Let H(T,u)=0 be the degree-50 equation of C above the q6 base.  In the
etale K-algebra A=K[u]/H, transport C to the certified q6 child and evaluate
the 51 standard basis functions of L(51 O_child):

    1,x,...,x^25, y,x*y,...,x^24*y.

Their 50 coefficient rows have a one-dimensional kernel.  The corresponding
function f=A(x)+y*B(x) has

    div(f) = D_C + Q - 51 O,

hence Trace(C)=-Q.

To recover x(Q) without a degree-50 resultant:
  * R(x)=A(x)^2-(x^3+a*x+b)B(x)^2 has the 50 x(C_i) plus x(Q) as roots;
  * the sum of all roots of R is read from its top two coefficients;
  * Trace_A/K(x_C) is computed from H by Newton power sums.

Thus x(Q)=sum_roots(R)-Trace(x_C), and y(Q) follows from f(Q)=0.

Finally subtract the oriented QO/QA correction and export a modular
x=X/Z^2, y=Y/Z^3 seed if the denominator profile is square/cubic.

Run:
  sage -python ~/Downloads/recover_h92_q6_third_trace_modp.sage
"""

import argparse
import json
from importlib.machinery import SourceFileLoader
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
parser.add_argument("--check-tau", type=int, default=2)
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
    else ROOT / "artifacts/local/elkies-k3/q6-third-trace-mod-100003.json"
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
    raise ValueError(
        "the default Hensel P2 input is tied to p=100003; use that prime"
    )
F = GF(p)

TR = PolynomialRing(F, "T")
T = TR.gen()
K = TR.fraction_field()

UR = PolynomialRing(K, "u")
u = UR.gen()
L = UR.fraction_field()


def modq(value):
    value = QQ(value)
    if ZZ(value.denominator()) % p == 0:
        raise ZeroDivisionError(f"denominator divisible by {p}: {value}")
    return F(ZZ(value.numerator())) / F(ZZ(value.denominator()))


def fpoly(values):
    return PolynomialRing(F, "v")([modq(value) for value in values])


# Use dedicated F[u] polynomials for source data, then lift coefficients to K.
FuR = PolynomialRing(F, "u0")
u0f = FuR.gen()


def polyF(values):
    return FuR([modq(value) for value in values])


def liftF(poly):
    return UR([K(value) for value in poly.list()])


def asL_from_Fpoly(poly):
    return L(liftF(poly))


p1data = json.loads(P1FILE.read_text())
p2data = json.loads(P2FILE.read_text())
rr = json.loads(RRFILE.read_text())
child = json.loads(CHILDFILE.read_text())

assert p1data["status"] == "PASS_EXACT_H92_P1"
assert p2data["schema"] == "elkies-k3.h92-p2-hensel-lift.v1"
assert p2data["complete"]
assert rr["status"] == "PASS_EXACT_GLOBAL_RR_KERNEL"
assert child["status"] == "PASS_EXACT_E8_E6_CHILD_JACOBIAN"

xp = asL_from_Fpoly(polyF(p1data["x_entrance_base"]["numerator_coefficients"])) / (
    asL_from_Fpoly(polyF(p1data["x_entrance_base"]["denominator_coefficients"]))
)
yp = asL_from_Fpoly(polyF(p1data["y_entrance_base"]["numerator_coefficients"])) / (
    asL_from_Fpoly(polyF(p1data["y_entrance_base"]["denominator_coefficients"]))
)


def reciprocal(values):
    answer = L(0)
    for index, value in enumerate(values):
        answer += L(K(modq(value))) / L(u**index)
    return answer


z2 = reciprocal(p2data["Z"])
xC = reciprocal(p2data["X"]) / z2**2
yC = reciprocal(p2data["Y"]) / z2**3

anchor = SourceFileLoader("q6_trace_generic_anchor", str(ANCHOR)).load_module()
r0, s0 = anchor.EXPECTED_H92
unused_ring, formulas = anchor.parse_h92(H92)
A1q, Aq, B1q, Bq, B2q = tuple(QQ(value(r0, s0)) for value in formulas)
A1, A, B1old, Bold, B2old = map(modq, (A1q, Aq, B1q, Bq, B2q))

old_a = L(K(A1)) / u**3 + L(K(A)) / u**4
old_b = (
    L(K(B1old)) / u**5
    + L(K(Bold)) / u**6
    + L(K(B2old)) / u**7
)
assert yp**2 == xp**3 + old_a*xp + old_b
assert yC**2 == xC**3 + old_a*xC + old_b

hF = polyF(p1data["structured_denominator"]["Z4_coefficients"])
h = liftF(hF)


def coefficient_pair(entry):
    Ap = liftF(polyF(entry["A_coefficients_low_to_high"]))
    Bp = liftF(polyF(entry["B_coefficients_low_to_high"]))
    return L(Ap) / L(h**2), L(Bp) / L(h)


(a0, b0), (a1, b1) = tuple(
    coefficient_pair(entry) for entry in rr["kernel"]["sections"]
)

mC = (yC - yp) / (xC - xp)
TC = (a1 + b1*mC) / (a0 + b0*mC)
assert max(TC.numerator().degree(), TC.denominator().degree()) == 50

# Degree-50 generic fibre of C over the new base T.
Hraw = UR((TC - L(K(T))).numerator())
if Hraw.degree() != 50:
    raise ArithmeticError(f"generic C fibre has degree {Hraw.degree()}, expected 50")
H = Hraw.monic()
assert H.gcd(H.derivative()).degree() == 0
print("Q6TRACEGEN|stage=degree50|status=PASS", flush=True)

# Reconstruct the q6 quartic over K=F_p(T).
q6_m = (a1 - L(K(T))*a0) / (L(K(T))*b0 - b1)
radicand = (
    q6_m**4
    - 6*xp*q6_m**2
    + 8*yp*q6_m
    - 3*xp**2
    - 4*old_a
)
quartic, square_factor = squarefree_binary_quartic(radicand, UR)
assert quartic.degree() == 4
print("Q6TRACEGEN|stage=quartic|degree=4|status=PASS", flush=True)


def native_covariants(quartic):
    base = quartic.base_ring()
    BR = PolynomialRing(base, names=("qx", "qz"))
    qx, qz = BR.gens()
    f = sum(
        BR(quartic[index]) * qx**index * qz**(4-index)
        for index in range(5)
    )
    HC = (
        f.derivative(qx, 2)*f.derivative(qz, 2)
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


fbin, HC, GC, I, J = native_covariants(quartic)
std_a = K(-27)*I
std_b = K(-27)*J

minA = TR([
    modq(value)
    for value in child["minimal_short_weierstrass"]["A_coefficients_low_to_high"]
])
minB = TR([
    modq(value)
    for value in child["minimal_short_weierstrass"]["B_coefficients_low_to_high"]
])
child_a = K(minA)
child_b = K(minB)


def rational_fourth_power_units(ratio):
    ratio = K(ratio)
    num = TR(ratio.numerator())
    den = TR(ratio.denominator())
    fn = num.factor()
    fd = den.factor()
    scalar = F(fn.unit()) / F(fd.unit())
    base_part = K(1)
    for factor, exponent in fn:
        if exponent % 4:
            raise ArithmeticError("numerator fourth-power exponent failure")
        base_part *= K(factor)**(exponent//4)
    for factor, exponent in fd:
        if exponent % 4:
            raise ArithmeticError("denominator fourth-power exponent failure")
        base_part /= K(factor)**(exponent//4)
    CR = PolynomialRing(F, "c")
    c = CR.gen()
    roots = [root for root, multiplicity in (c**4-scalar).roots() if multiplicity == 1]
    return [K(root)*base_part for root in roots]


units = [
    candidate
    for candidate in rational_fourth_power_units(child_a/std_a)
    if candidate**6*std_b == child_b
]
if not units:
    raise ArithmeticError("could not recover global minimalizing unit")
minimalizing_unit = units[0]
assert minimalizing_unit**4*std_a == child_a
assert minimalizing_unit**6*std_b == child_b
print(
    f"Q6TRACEGEN|stage=minimalizing_unit|candidates={len(units)}|status=PASS",
    flush=True,
)

# Verify q6_m and mC agree in A=K[u]/H.
def reduce_A(value):
    value = L(value)
    num = UR(value.numerator())
    den = UR(value.denominator())
    if den.gcd(H).degree() != 0:
        raise ZeroDivisionError("denominator is not a unit modulo H")
    return (num * den.inverse_mod(H)) % H


assert reduce_A(q6_m - mC) == 0

# C's quartic square root.
wC = (2*xC + xp - mC**2) / square_factor
assert reduce_A(wC**2 - quartic) == 0

# Evaluate the covariant formulas at z=1.
H_u = UR(HC(qx=u, qz=K(1)))
G_u = UR(GC(qx=u, qz=K(1)))
raw_xC = -L(K(3))/L(K(4)) * L(H_u)/L(quartic)
raw_yC = (
    L(K(9))/L(K(32))
    * L(G_u)*wC/L(quartic)**2
)
child_xC = minimalizing_unit**2 * raw_xC
child_yC = minimalizing_unit**3 * raw_yC

xA = reduce_A(child_xC)
yA = reduce_A(child_yC)
assert (
    (yA*yA - xA*xA*xA - UR(child_a)*xA - UR(child_b)) % H
) == 0
print("Q6TRACEGEN|stage=child_point_in_etale_algebra|status=PASS", flush=True)

# Build the 51 evaluations of L(51 O).
one = UR(1)
x_powers = [one]
for _ in range(25):
    x_powers.append((x_powers[-1]*xA) % H)
columns = list(x_powers)
for exponent in range(25):
    columns.append((yA*x_powers[exponent]) % H)
assert len(columns) == 51

evaluation = matrix(
    K,
    50,
    51,
    lambda row, col: K(columns[col][row]),
)
print(
    f"Q6TRACEGEN|stage=rr_matrix|dimensions={evaluation.dimensions()}",
    flush=True,
)
kernel = evaluation.right_kernel().basis_matrix()
if kernel.nrows() != 1:
    raise ArithmeticError(f"L(51O) evaluation kernel has dimension {kernel.nrows()}")
relation = kernel[0]
print("Q6TRACEGEN|stage=rr_kernel|dimension=1|status=PASS", flush=True)

XR = PolynomialRing(K, "X")
Xv = XR.gen()
Afun = sum(relation[index]*Xv**index for index in range(26))
Bfun = sum(relation[26+index]*Xv**index for index in range(25))
intersection = Afun**2 - (Xv**3 + child_a*Xv + child_b)*Bfun**2
if intersection.degree() != 51:
    raise ArithmeticError(
        f"residual intersection polynomial has degree {intersection.degree()}, expected 51"
    )
sum_intersection_roots = -intersection[50]/intersection[51]

# Newton power sums for roots of monic H.
n = H.degree()
assert n == 50 and H[n] == 1
power_sums = [K(n)]
for k in range(1, n):
    total = K(k)*K(H[n-k])
    for j in range(1, k):
        total += K(H[n-j])*power_sums[k-j]
    power_sums.append(-total)

trace_x = sum(K(xA[index])*power_sums[index] for index in range(n))
xQ = K(sum_intersection_roots - trace_x)
b_at_q = K(Bfun(xQ))
if not b_at_q:
    raise ArithmeticError("residual point has B(x_Q)=0; alternate recovery needed")
yQ = -K(Afun(xQ))/b_at_q

E = EllipticCurve(K, [0, 0, 0, child_a, child_b])
Q = E(xQ, yQ)
traceC = -Q
print("Q6TRACEGEN|stage=trace_point|status=PASS", flush=True)

# Intrinsic old-O / A=-P1 orientation at the common chord-pole branch.
B0F = polyF(rr["kernel"]["sections"][0]["B_coefficients_low_to_high"])
B1F = polyF(rr["kernel"]["sections"][1]["B_coefficients_low_to_high"])
B0 = liftF(B0F)
B1 = liftF(B1F)
line, rem = (B1 - UR(K(T))*B0).quo_rem(u**3)
if rem or line.degree() != 1:
    raise ArithmeticError("generic pole branch is not linear after u^3 removal")
uO = -K(line[0])/K(line[1])


def removable_eval(value, at):
    value = L(value)
    num = UR(value.numerator())
    den = UR(value.denominator())
    linear = u - K(at)
    while num(K(at)) == 0 and den(K(at)) == 0:
        num, rn = num.quo_rem(linear)
        den, rd = den.quo_rem(linear)
        if rn or rd:
            raise ArithmeticError("local common-factor cancellation failed")
    if den(K(at)) == 0:
        raise ZeroDivisionError("genuine pole at oriented branch")
    return K(num(K(at))/den(K(at)))


wO = removable_eval(q6_m**2/square_factor, uO)
if wO**2 != K(quartic(uO)):
    raise ArithmeticError("oriented old-O quartic sign failed")

def map_quartic_K(xvalue, zvalue, wvalue):
    fv = K(fbin(qx=K(xvalue), qz=K(zvalue)))
    hv = K(HC(qx=K(xvalue), qz=K(zvalue)))
    gv = K(GC(qx=K(xvalue), qz=K(zvalue)))
    if not fv or K(wvalue)**2 != fv:
        raise ArithmeticError("quartic K-point failed")
    rx = -K(3)/K(4)*hv/fv
    ry = K(9)/K(32)*gv*K(wvalue)/fv**2
    cx = minimalizing_unit**2*rx
    cy = minimalizing_unit**3*ry
    P = E(cx, cy)
    return P

QO = map_quartic_K(uO, 1, wO)
QA = map_quartic_K(uO, 1, -wO)
assert QO != QA

S3 = traceC - 25*QO - 24*QA

# Choose the global y-orientation by the already-certified tau=2 smoke value.
tau = F(args.check_tau)
if minA(tau) == 0 and minB(tau) == 0:
    raise ArithmeticError("check specialization is degenerate")

def specialize_K(value, at):
    value = K(value)
    top = TR(value.numerator())
    bottom = TR(value.denominator())
    if bottom(at) == 0:
        raise ZeroDivisionError("check specialization hits a pole")
    return F(top(at)/bottom(at))

sx, sy = S3.xy()
sx2, sy2 = specialize_K(sx, tau), specialize_K(sy, tau)
target_x, target_y = F(63169), F(42300)
if (sx2, sy2) == (target_x, target_y):
    orientation = "direct"
elif (sx2, sy2) == (target_x, -target_y):
    S3 = -S3
    sx, sy = S3.xy()
    sx2, sy2 = specialize_K(sx, tau), specialize_K(sy, tau)
    orientation = "global_y_negated"
else:
    raise ArithmeticError(
        f"generic trace specialization {(int(sx2), int(sy2))} "
        "does not match oriented tau=2 smoke value"
    )
assert (sx2, sy2) == (target_x, target_y)
print(
    f"Q6TRACEGEN|stage=orientation|tau={int(tau)}|"
    f"point={int(sx2)},{int(sy2)}|mode={orientation}|status=PASS",
    flush=True,
)

sx = K(S3[0])
sy = K(S3[1])
xnum, xden = TR(sx.numerator()), TR(sx.denominator())
ynum, yden = TR(sy.numerator()), TR(sy.denominator())

# Normalize through the common pole polynomial Z whenever possible.
if not xden.is_square():
    raise ArithmeticError(
        f"x denominator degree {xden.degree()} is not a square"
    )
Z = xden.sqrt()
if Z.leading_coefficient() != 1:
    Z /= Z.leading_coefficient()
Xcoord = K(sx) * K(Z**2)
Ycoord = K(sy) * K(Z**3)
if Xcoord.denominator() != 1:
    raise ArithmeticError("x*Z^2 is not polynomial")
if Ycoord.denominator() != 1:
    raise ArithmeticError("y*Z^3 is not polynomial")
Xcoord = TR(Xcoord.numerator())
Ycoord = TR(Ycoord.numerator())

# Direct child identity.
assert Ycoord**2 == Xcoord**3 + minA*Xcoord*Z**4 + minB*Z**6

print(
    "Q6TRACEGEN_PROFILE|"
    f"Z_degree={Z.degree()}|X_degree={Xcoord.degree()}|Y_degree={Ycoord.degree()}|"
    f"x_raw={xnum.degree()}/{xden.degree()}|"
    f"y_raw={ynum.degree()}/{yden.degree()}|status=PASS",
    flush=True,
)

payload = {
    "schema": "elkies-k3.h92-q6-third-trace-modp.v1",
    "status": "PASS_MODULAR_Q6_THIRD_TRACE_SECTION",
    "prime": int(p),
    "method": {
        "multisection_degree": 50,
        "trace_rr_space": "L(51*O)",
        "evaluation_matrix": [50, 51],
        "kernel_dimension": 1,
        "residual_point_recovery": (
            "x(Q)=sum_roots(A^2-E_rhs*B^2)-Trace_A/K(x_C), "
            "with Trace from Newton sums of H"
        ),
        "literal_support_correction": "S3=Trace(C)-25*QO-24*QA",
        "orientation": orientation,
    },
    "check_specialization": {
        "T": int(tau),
        "x": int(sx2),
        "y": int(sy2),
    },
    "profile": {
        "Z_degree": int(Z.degree()),
        "X_degree": int(Xcoord.degree()),
        "Y_degree": int(Ycoord.degree()),
        "x_numerator_degree": int(xnum.degree()),
        "x_denominator_degree": int(xden.degree()),
        "y_numerator_degree": int(ynum.degree()),
        "y_denominator_degree": int(yden.degree()),
    },
    "Z": [int(value) for value in Z.list()],
    "X": [int(value) for value in Xcoord.list()],
    "Y": [int(value) for value in Ycoord.list()],
    "boundary": (
        "This is an exact section over GF(100003)(T) on the certified q6 child. "
        "Characteristic-zero reconstruction and exact Shioda/component checks "
        "remain downstream."
    ),
}

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}")
print("Q6TRACEGEN_RESULT|status=PASS_MODULAR_Q6_THIRD_TRACE_SECTION", flush=True)
