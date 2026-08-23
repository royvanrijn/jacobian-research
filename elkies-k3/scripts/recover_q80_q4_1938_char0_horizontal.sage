#!/usr/bin/env sage
"""
Recover the exact q4_1938 horizontal section on the exact q6_7774 child.

Uses the pinned modular section
    X = Nx/(T+26)^2
    Y = Ny/(T+26)^3
and the exact I2 factor of the characteristic-zero 7774 child reducing to T=47.

The 14 section coefficients are lifted by a full-rank Q_73 Newton system and
recognized in QQ(sqrt(-3)), then verified by an exact Weierstrass identity.
"""

from pathlib import Path
import json

from sage.all import (
    QQ, ZZ, GF, Qp, PolynomialRing, QuadraticField,
    Matrix, vector, sage_eval
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ART = ROOT / "artifacts" / "generated-results" / "q80-q6-7774-char0-resolved-rr.json"
OUTDIR = ROOT / "elkies-k3" / "data" / "fibrations" / "q80-q6-7774-char0"
OUTDIR.mkdir(parents=True, exist_ok=True)

if not ART.exists():
    raise SystemExit(
        "Missing resolved 7774 artifact. First rerun certify_q80_q6_7774_char0_resolved_rr.sage"
    )

data = json.loads(ART.read_text())
if data.get("status") != "PASS_EXACT_Q6_7774_RESOLVED_RR":
    raise SystemExit("7774 artifact is not certified")

K = QuadraticField(-3, "j")
j = K.gen()
R = PolynomialRing(K, "T")
T = R.gen()

A = R(sage_eval(data["child"]["jacobian_A"], locals={"j":j, "T":T}))
B = R(sage_eval(data["child"]["jacobian_B"], locals={"j":j, "T":T}))
Delta = -16*(4*A^3+27*B^2)
assert (A.degree(), B.degree(), Delta.degree()) == (8,12,18)  # affine Delta; I6 at infinity supplies remaining 6

p = 73
F = GF(p)
JMOD = F(17)
assert JMOD^2 == F(-3)

def red_q(q):
    q = QQ(q)
    return F(q.numerator())/F(q.denominator())

def red_k(z):
    z = K(z)
    cc = list(z)+[QQ(0),QQ(0)]
    return red_q(cc[0]) + JMOD*red_q(cc[1])

Rp = PolynomialRing(F, "T")
Tp = Rp.gen()

def red_poly(poly):
    return Rp([red_k(c) for c in R(poly).list()])

# Pinned GF(73) 7774 model.
Apin = (
    46*Tp^8 + 5*Tp^7 + 16*Tp^6 + 44*Tp^5 + 6*Tp^4
    + 13*Tp^3 + Tp^2 + Tp
)
Bpin = (
    54*Tp^12 + 58*Tp^11 + 48*Tp^10 + 16*Tp^9 + 42*Tp^8
    + 67*Tp^7 + 25*Tp^6 + 19*Tp^5 + 27*Tp^4 + 45*Tp^3
    + 61*Tp^2 + 44*Tp + 49
)

scale_marker = int(data["child"]["gf73_scale_marker"])
if scale_marker <= 0:
    raise ArithmeticError(
        "This recovery currently expects the direct constant scaling orientation; "
        f"got scale marker {scale_marker}"
    )
u = F(scale_marker)
assert red_poly(A) == u^4*Apin
assert red_poly(B) == u^6*Bpin

# Identify exact I2 factor reducing to T=47.  Use already-classified finite
# factors from the 7774 artifact instead of refactoring Delta.
i2_candidates = []
for item in data["child"]["finite_fibres"]:
    if item["kodaira"] != "I2" or int(item["degree"]) != 1:
        continue
    f = R(sage_eval(item["factor"], locals={"j":j, "T":T}))
    root = K(-f[0]/f[1])
    i2_candidates.append((root, int(red_k(root))))

matches = [root for root, residue in i2_candidates if residue == 47]
if len(matches) != 1:
    raise ArithmeticError(f"expected one I2 root reducing to 47, got {i2_candidates}")
pole = matches[0]
q = R(T-pole)
assert red_poly(q) == Tp+26

print(
    f"Q801938H|pole={pole}|pole_mod73=47|q={q}|"
    "status=PASS_EXACT_POLE_FACTOR",
    flush=True,
)

# Pinned q4_1938 horizontal over the target modular 7774 model.
#
# IMPORTANT: the historical reconstruction was NOT a free 14-coefficient
# rational-section lift. The known component profile is
#
#     (I2,I2,I5,I5,I6) = (0,1,4,4,3),
#
# so H is nonidentity at both finite I5 fibres. On the singular Weierstrass
# model that means the section passes through both I5 nodes. We therefore
# impose the exact node constraints before lifting:
#
#     N = N0 + P5*Q,   deg Q <= 4
#     M = P5*U,        deg U <= 4
#
# where X=N/q^2, Y=M/q^3 and P5 is the product of the two exact I5 factors.
# This is the characteristic-zero port of the repo's original
# "node-constrained polynomial search".
Xnum_pin = Rp([60,40,34,58,48,30,3])
Ynum_pin = Rp([24,15,0,57,28,65,15])
qpin = Tp+26

# Convert the pinned section to the actual constant scaling of our exact
# q6_7774 Jacobian.
Nseed = u^2*Xnum_pin
Mseed = u^3*Ynum_pin
Ared, Bred = red_poly(A), red_poly(B)

if Mseed^2 != Nseed^3 + Ared*Nseed*qpin^4 + Bred*qpin^6:
    Mseed = -Mseed
assert Mseed^2 == Nseed^3 + Ared*Nseed*qpin^4 + Bred*qpin^6

def coeffp(poly, i):
    return poly[i] if i <= poly.degree() else F(0)

# Exact two I5 fibres, selected by their pinned reductions 14 and 42.
i5_candidates = []
for item in data["child"]["finite_fibres"]:
    if item["kodaira"] != "I5" or int(item["degree"]) != 1:
        continue
    f = R(sage_eval(item["factor"], locals={"j":j, "T":T}))
    root = K(-f[0]/f[1])
    i5_candidates.append((root, int(red_k(root))))

by_residue = {residue: root for root, residue in i5_candidates}
if set(by_residue) != {14,42}:
    raise ArithmeticError(f"expected I5 roots reducing to 14,42; got {i5_candidates}")
r14 = by_residue[14]
r42 = by_residue[42]
P5 = R((T-r14)*(T-r42))
assert red_poly(P5) == (Tp-14)*(Tp-42)

def node_x(root):
    ar, br = K(A(root)), K(B(root))
    if ar == 0:
        raise ArithmeticError("unexpected c4=0 at multiplicative fibre")
    x0 = K(-3*br/(2*ar))
    assert x0^3 + ar*x0 + br == 0
    assert 3*x0^2 + ar == 0
    return x0

# N(root)=x_node*q(root)^2 at each I5 node.
Ninterp = R.lagrange_polynomial([
    (r14, node_x(r14)*q(r14)^2),
    (r42, node_x(r42)*q(r42)^2),
])
assert Ninterp.degree() <= 1

P5p = red_poly(P5)
Ninterpp = red_poly(Ninterp)

Qseed, rem = (Nseed-Ninterpp).quo_rem(P5p)
assert rem == 0 and Qseed.degree() <= 4
Useed, rem = Mseed.quo_rem(P5p)
assert rem == 0 and Useed.degree() <= 4

print(
    f"Q801938HSTRUCT|I5_roots=({r14},{r42})|"
    f"P5={P5}|Ninterp={Ninterp}|"
    "status=PASS_EXACT_NODE_CONSTRAINTS",
    flush=True,
)

# Modular Jacobian of the structured ten-variable system.
baseN = -(3*Nseed^2 + Ared*qpin^4)
baseM = 2*Mseed
cols = (
    [baseN*P5p*Tp^i for i in range(5)]
    + [baseM*P5p*Tp^i for i in range(5)]
)
equation_rows = max(col.degree() for col in cols)+1
Jp = Matrix(
    F, equation_rows, 10,
    lambda row,col: coeffp(cols[col], row)
)
rank = Jp.rank()
print(
    f"Q801938HSTRUCT|modular_jacobian_rank={rank}|unknowns=10|"
    f"equations={equation_rows}|status=PASS_STRUCTURED_LINEARIZATION",
    flush=True,
)
if rank != 10:
    raise ArithmeticError(
        f"node-constrained section Jacobian rank {rank}, expected 10"
    )

pivot_rows = tuple(Jp.transpose().pivots())
assert len(pivot_rows) == 10

# Q_73 Newton lift of Q and U, keeping the exact node incidence built in.
PREC = 440
Q73 = Qp(73, PREC, type="capped-rel")
j73 = next(x for x in Q73(-3).sqrt(all=True) if int(x.residue()) == 17)

def K_to_Q73(c):
    c = K(c)
    cc = list(c)+[QQ(0),QQ(0)]
    return Q73(QQ(cc[0])) + j73*Q73(QQ(cc[1]))

R73 = PolynomialRing(Q73, "T")
T73 = R73.gen()

def poly73(poly):
    return R73([K_to_Q73(c) for c in R(poly).list()])

A73 = poly73(A)
B73 = poly73(B)
q73 = poly73(q)
P573 = poly73(P5)
Ninterp73 = poly73(Ninterp)

def coeff73(poly, i):
    return poly[i] if i <= poly.degree() else Q73(0)

def state(z):
    Q = R73(list(z[:5]))
    U = R73(list(z[5:10]))
    N = Ninterp73 + P573*Q
    M = P573*U
    E = M^2 - N^3 - A73*N*q73^4 - B73*q73^6
    return Q,U,N,M,E

def jacobian(N,M):
    baseN = -(3*N^2 + A73*q73^4)
    baseM = 2*M
    cc = (
        [baseN*P573*T73^i for i in range(5)]
        + [baseM*P573*T73^i for i in range(5)]
    )
    return Matrix(
        Q73, 10, 10,
        lambda rr,col: coeff73(cc[col], pivot_rows[rr])
    )

seed = (
    [int(coeffp(Qseed,i)) for i in range(5)]
    + [int(coeffp(Useed,i)) for i in range(5)]
)
z = vector(Q73, [Q73(v) for v in seed])

for iteration in range(12):
    Qcur,Ucur,Ncur,Mcur,Ecur = state(z)
    residual = vector(Q73, [coeff73(Ecur,row) for row in pivot_rows])
    J = jacobian(Ncur,Mcur)
    z -= J.solve_right(residual)

    _,_,_,_,Echeck = state(z)
    vals = [c.valuation() for c in Echeck.list() if c != 0]
    minval = min(vals) if vals else PREC
    print(
        f"Q801938HSTRUCT|newton={iteration+1}|min_v73={minval}|status=NEWTON",
        flush=True,
    )
    if minval >= 330:
        break

_,_,_,_,Epad = state(z)
vals = [c.valuation() for c in Epad.list() if c != 0]
assert not vals or min(vals) >= 260

def qq_sqrt(x):
    x = QQ(x)
    if x < 0:
        return None
    num, den = x.numerator(), x.denominator()
    if not num.is_square() or not den.is_square():
        return None
    return QQ(num.sqrt())/QQ(den.sqrt())

def recognize_K(alpha, residue):
    dep = alpha.algebraic_dependency(2)
    if dep.degree() == 1:
        val = K(-QQ(dep[0])/QQ(dep[1]))
        if red_k(val) == F(residue):
            return val
    if dep.degree() != 2:
        raise ArithmeticError(f"unexpected algdep degree: {dep}")
    c0,c1,c2 = map(QQ,(dep[0],dep[1],dep[2]))
    aa = -c1/(2*c2)
    norm = c0/c2
    bb2 = (norm-aa^2)/3
    bb = qq_sqrt(bb2)
    if bb is None:
        raise ArithmeticError(f"algdep is not QQ(sqrt(-3))-shaped: {dep}")
    for sign in (1,-1):
        val = K(aa+sign*bb*j)
        if red_k(val) == F(residue):
            return val
    raise ArithmeticError(
        f"could not choose QQ(sqrt(-3)) root matching residue {residue}: {dep}"
    )

# Only Q needs p-adic recognition. Once N is exact, the curve equation
# determines U^2 exactly, so reconstruct U algebraically instead of running
# fragile algdep on another five huge coefficients.
qexact_coeffs = []
for idx,(alpha,residue) in enumerate(zip(z[:5],seed[:5])):
    val = recognize_K(alpha,residue)
    qexact_coeffs.append(val)
    print(
        f"Q801938HSTRUCT|recognized={idx}|value={val}|status=ALGDEP",
        flush=True,
    )

Qexact = R(qexact_coeffs)
N = R(Ninterp + P5*Qexact)

rhs = R(N^3 + A*N*q^4 + B*q^6)
U2, rem = rhs.quo_rem(P5^2)
assert rem == 0
U2 = R(U2)

# Exact polynomial square-root certificate.
fac = U2.factor()
unit = K(fac.unit())
assert unit.is_square(), ("U^2 scalar unit is nonsquare", unit)
Uroot = R(unit.sqrt())
for f,e in fac:
    assert int(e) % 2 == 0, ("U^2 has odd factor", f, e)
    Uroot *= f^(int(e)//2)

Ucandidates = (R(Uroot), R(-Uroot))
Uexact = next((uu for uu in Ucandidates if red_poly(uu) == Useed), None)
if Uexact is None:
    raise ArithmeticError(
        "exact square root exists but neither sign matches modular U seed"
    )

M = R(P5*Uexact)
assert Uexact^2 == U2
print(
    f"Q801938HSTRUCT|U={Uexact}|"
    "status=PASS_EXACT_U_FROM_CURVE_SQUARE",
    flush=True,
)

assert M^2 == N^3 + A*N*q^4 + B*q^6
Hx = N/q^2
Hy = M/q^3
assert Hy^2 == Hx^3 + A*Hx + B

# Exact node-incidence audit on both I5 fibres.
for root in (r14,r42):
    assert K(Hx(root)) == node_x(root)
    assert M(root) == 0

assert red_poly(N) == Nseed
assert red_poly(M) == Mseed

print(
    f"Q801938H|N={N}|M={M}|q={q}|P5={P5}|"
    "status=PASS_EXACT_Q4_1938_HORIZONTAL",
    flush=True,
)

section_file = OUTDIR / "q80_char0_q6_7774_H1938.sage"
section_file.write_text(
    "\n".join([
        "#!/usr/bin/env sage",
        "from sage.all import PolynomialRing, QuadraticField",
        'K = QuadraticField(-3, "j")',
        "j = K.gen()",
        'R = PolynomialRing(K, "T")',
        "T = R.gen()",
        f"A = {A}",
        f"B = {B}",
        f"q = {q}",
        f"N_H = {N}",
        f"M_H = {M}",
        "H_x = N_H/q^2",
        "H_y = M_H/q^3",
        "assert H_y^2 == H_x^3 + A*H_x + B",
        'print("Q801938H|status=PASS_EXACT_Q4_1938_HORIZONTAL")',
    ]) + "\n"
)

note = OUTDIR / "Q80_CHAR0_Q4_1938_HORIZONTAL.md"
note.write_text(
    "# Q80 q4_1938 — exact horizontal on q6_7774\n\n"
    "Status: **PASS_EXACT_Q4_1938_HORIZONTAL**\n\n"
    f"- exact pole factor: `q={q}` (root reduces to `T=47`)\n"
    "- modular denominator regression: `q -> T+26`\n"
    f"- exact x numerator: `{N}`\n"
    f"- exact y numerator: `{M}`\n"
    "- exact Weierstrass section identity verified\n"
    "- reduction matches the pinned q4_1938 horizontal after the certified "
    f"7774 constant scaling `u={scale_marker}`\n\n"
    "Next: compile the smooth saturation and connected-A4 quotient to the "
    "two-dimensional q4_1938 Riemann--Roch pencil.\n"
)

print(
    f"Q801938HFINAL|section={section_file}|note={note}|"
    "status=PASS_EXACT_Q4_1938_HORIZONTAL",
    flush=True,
)
