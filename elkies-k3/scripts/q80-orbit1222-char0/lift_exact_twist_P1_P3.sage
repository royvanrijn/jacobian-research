#!/usr/bin/env sage
"""
Q80 orbit 1222: exact twist from P1/P3 with the FULL reducible-fiber
component multiplicities imposed.

Critical correction:
For split multiplicative I_n, a section on component k has y-vanishing order
min(k,n-k) at the Weierstrass node.

Pinned CM24 marking:
  P1:
    I7 labels (5,2) -> y-order 2 at each I7
    two nonidentity I2 fibers -> y-order 1 at each selected I2
    therefore W1 = P7^2 * Q2 * u1, with u1 constant.

  P3:
    I7 labels (5,5) -> y-order 2 at each I7
    no I2 support
    therefore W3 = P7^2 * U3, deg(U3)<=2.

Base exact j/fiber model:
    E0: y^2 = x^3 + A*x + B.

Twist representation:
    Z^3 + A Z + B = d W^2.
Then x=d Z, y=d^2 W lies on
    y^2 = x^3 + d^2 A x + d^3 B.

Z constraints:
  P1:
    M1 = P7*Q2
    Z1 = interp1 + lambda*M1
    F(Z1)/M1^2 = d * P7^2 * u1^2

  P3:
    Z3 = node7 + P7*(q0+q1 V+q2 V^2)
    F(Z3)/P7^2 = d * P7^2 * U3^2

Global square gauge:
    d -> d*s^2, u1,U3 -> u1/s,U3/s.
Fix lead(U3)=1.

Unknowns (8):
    d, lambda, u1, q0,q1,q2, U30,U31.

Equations:
    5 P1 quotient coefficients + 9 P3 quotient coefficients = 14.

The p=73 seed has full column Jacobian rank 8.
"""

import contextlib
import io
import itertools
import time
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, QQ, Qp, vector

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
ORBIT_DATA = (
    REPO_ROOT
    / "elkies-k3"
    / "data"
    / "fibrations"
    / "q80-orbit1222-char0"
)
BASE_MODEL = ORBIT_DATA / "q80_char0_orbit1222_weierstrass.sage"

OUT_MODEL = ORBIT_DATA / "q80_char0_orbit1222_jacobian.sage"
OUT_SECTIONS = ORBIT_DATA / "q80_char0_orbit1222_P1_P3.sage"
OUT_NOTE = ORBIT_DATA / "Q80_CHAR0_ORBIT1222_JACOBIAN_TWIST.md"

if not BASE_MODEL.exists():
    raise SystemExit(f"missing exact base model: {BASE_MODEL}")

with contextlib.redirect_stdout(io.StringIO()):
    load(str(BASE_MODEL))

K0, j0, R0, V0 = K, j, R, V
A0, B0, Delta0 = A, B, Delta

# ---------------------------------------------------------------------------
# Exact reducible-fiber support.
# ---------------------------------------------------------------------------

G1 = Delta0.gcd(Delta0.derivative()).monic()
G2 = G1.gcd(G1.derivative()).monic()
G3 = G2.gcd(G2.derivative()).monic()
assert (G1.degree(), G2.degree(), G3.degree()) == (15, 10, 8)

P7, rem = G2.quo_rem(G3)
assert rem == 0
P7 = P7.monic()

P7P2, rem = G1.quo_rem(G2)
assert rem == 0
P7P2 = P7P2.monic()
P2, rem = P7P2.quo_rem(P7)
assert rem == 0
P2 = P2.monic()

assert (P7.degree(), P2.degree()) == (2, 3)


def node_residue(modulus):
    den = (2*A0).mod(modulus)
    return (-3*B0*den.inverse_mod(modulus)).mod(modulus)


node7 = node_residue(P7)
node2 = node_residue(P2)

# ---------------------------------------------------------------------------
# p=73 pinned model and deterministic P1/P3.
# ---------------------------------------------------------------------------

p = 73
F = GF(p)
Rp = PolynomialRing(F, "V")
Vp = Rp.gen()
Xring = PolynomialRing(F, "x")
xx = Xring.gen()


def red_Q(q):
    q = QQ(q)
    return F(q.numerator()) / F(q.denominator())


def red_K(c):
    cc = K0(c)
    return red_Q(cc[0]) + F(17)*red_Q(cc[1])


def red_poly(poly):
    return Rp([red_K(c) for c in poly.list()])


Abase = red_poly(A0)
Bbase = red_poly(B0)
P7p = red_poly(P7)
P2p = red_poly(P2)
node7p = red_poly(node7)
node2p = red_poly(node2)

dlocal = F(17)
Aold = dlocal**2 * Abase
Bold = dlocal**3 * Bbase

Aold_expected = (
    6*Vp**8 + 16*Vp**7 + 47*Vp**6 + 33*Vp**5 + 58*Vp**4
    + 2*Vp**3 + 63*Vp**2 + 17*Vp + 23
)
Bold_expected = (
    33*Vp**12 + 64*Vp**10 + 61*Vp**9 + 45*Vp**8 + 14*Vp**7
    + 20*Vp**6 + 54*Vp**5 + 8*Vp**4 + 50*Vp**3
    + 57*Vp**2 + 47*Vp + 43
)
assert Aold == Aold_expected
assert Bold == Bold_expected

i7_roots = (F(-20), F(-67))
i2_roots = (F(-17), F(-30), F(-68))


def node_x_old(root):
    cubic = xx**3 + Aold(root)*xx + Bold(root)
    common = cubic.gcd(cubic.derivative())
    assert common.degree() == 1
    return -common[0]/common[1]


old_nodes = {r: node_x_old(r) for r in i7_roots+i2_roots}
old_node7 = Rp.lagrange_polynomial([(r, old_nodes[r]) for r in i7_roots])

old_P7 = Rp.one()
for r in i7_roots:
    old_P7 *= Vp-r
old_P7 = old_P7.monic()


def polynomial_square_roots(poly):
    assert poly.degree() <= 12
    if poly == 0:
        return (Rp.zero(),)

    shift = next(v for v in F if poly(v) != 0)
    shifted = poly(Vp+shift)
    constant = shifted[0]
    if not constant.is_square():
        return ()

    roots = []
    for first in constant.sqrt(all=True):
        coeffs = [first]
        for degree in range(1, 7):
            known = sum(
                coeffs[left]*coeffs[degree-left]
                for left in range(1, degree)
            )
            coeffs.append((shifted[degree]-known)/(2*first))
        candidate = Rp(coeffs)
        if candidate**2 == shifted:
            roots.append(candidate(Vp-shift))
    return tuple(roots)


started = time.monotonic()
candidates = []
for coeffs in itertools.product(F, repeat=3):
    X = old_node7 + Rp(coeffs)*old_P7
    for Y in polynomial_square_roots(X**3+Aold*X+Bold):
        candidates.append((X, Y))

assert len(candidates) == 30
P1oldx, P1oldy = candidates[20]
P3oldx, P3oldy = candidates[2]

P1_i2_support = tuple(
    r for r in i2_roots
    if P1oldx(r) == old_nodes[r] and P1oldy(r) == 0
)
P3_i2_support = tuple(
    r for r in i2_roots
    if P3oldx(r) == old_nodes[r] and P3oldy(r) == 0
)

assert len(P1_i2_support) == 2
assert len(P3_i2_support) == 0

print(
    "Q80TWCOMP|"
    f"P1_i2_support={tuple(int(r) for r in P1_i2_support)}|"
    f"search_seconds={time.monotonic()-started:.3f}|"
    "stage=PINNED_SUPPORT",
    flush=True,
)

# ---------------------------------------------------------------------------
# Exact degree-2 I2 support Q2 selected by P1.
# ---------------------------------------------------------------------------

target_Q2p = Rp.one()
for r in P1_i2_support:
    target_Q2p *= Vp-r
target_Q2p = target_Q2p.monic()

fac = tuple(P2.factor())
divisors = []


def enum_divisors(index, current):
    if index == len(fac):
        if current.degree() == 2:
            divisors.append(current.monic())
        return

    g, mult = fac[index]
    acc = R0.one()
    for _ in range(mult+1):
        nxt = current*acc
        if nxt.degree() <= 2:
            enum_divisors(index+1, nxt)
        acc *= g


enum_divisors(0, R0.one())

Q2matches = [
    q for q in divisors
    if red_poly(q).degree() == 2
    and red_poly(q).monic() == target_Q2p
]
assert len(Q2matches) == 1
Q2 = Q2matches[0]

# Exact P1 x-interpolation.
node2_Q2 = node2.mod(Q2)
rhs = (node2_Q2-node7.mod(Q2)).mod(Q2)
t = (rhs * P7.mod(Q2).inverse_mod(Q2)).mod(Q2)
interp1 = node7 + P7*t
M1 = (P7*Q2).monic()

assert interp1.mod(P7) == node7.mod(P7)
assert interp1.mod(Q2) == node2_Q2
assert M1.degree() == 4

# ---------------------------------------------------------------------------
# Convert pinned sections to base twist representation.
# ---------------------------------------------------------------------------


def to_ZW(X, Y):
    Z = X/dlocal
    W = Y/(dlocal**2)
    assert Z**3 + Abase*Z + Bbase == dlocal*W**2
    return Z, W


Z1p, W1p = to_ZW(P1oldx, P1oldy)
Z3p, W3p = to_ZW(P3oldx, P3oldy)

interp1p = red_poly(interp1)
M1p = red_poly(M1)

lambda_poly, rem = (Z1p-interp1p).quo_rem(M1p)
assert rem == 0 and lambda_poly.degree() <= 0
lambda_seed = lambda_poly[0]

q3p, rem = (Z3p-node7p).quo_rem(P7p)
assert rem == 0 and q3p.degree() <= 2

# FULL component multiplicity factors.
D1p = P7p**2 * target_Q2p
D3p = P7p**2

u1p, rem = W1p.quo_rem(D1p)
assert rem == 0 and u1p.degree() <= 0

U3p, rem = W3p.quo_rem(D3p)
assert rem == 0 and U3p.degree() <= 2
assert U3p.degree() == 2 and U3p[2] != 0

# Regression: exact modular vanishing orders encoded by factors.
assert P1oldy.valuation(Vp-i7_roots[0]) == 2
assert P1oldy.valuation(Vp-i7_roots[1]) == 2
for r in P1_i2_support:
    assert P1oldy.valuation(Vp-r) == 1

assert P3oldy.valuation(Vp-i7_roots[0]) == 2
assert P3oldy.valuation(Vp-i7_roots[1]) == 2

print(
    "Q80TWCOMP|"
    "P1_orders=I7:2,2;I2:1,1|"
    "P3_orders=I7:2,2|"
    "stage=PASS_COMPONENT_VANISHING_ORDERS",
    flush=True,
)

# One global square gauge: lead(U3)=1.
s = U3p[2]
dseed = dlocal*s**2
u1p /= s
U3p /= s
assert U3p[2] == 1

assert Z1p**3+Abase*Z1p+Bbase == dseed*(D1p*u1p)**2
assert Z3p**3+Abase*Z3p+Bbase == dseed*(D3p*U3p)**2

seed = [
    int(dseed),
    int(lambda_seed),
    int(u1p[0]),
    int(q3p[0]), int(q3p[1]), int(q3p[2]),
    int(U3p[0]), int(U3p[1]),
]
assert len(seed) == 8

print(
    "Q80TWCOMP|"
    f"seed_d={seed[0]}|lambda={seed[1]}|u1={seed[2]}|"
    f"q3={tuple(seed[3:6])}|U3={tuple(seed[6:])}+(1,)|"
    "stage=COMPONENT_SEED",
    flush=True,
)

# ---------------------------------------------------------------------------
# Quotient equations and 14x8 modular Jacobian.
# ---------------------------------------------------------------------------


def quotient_mod_p(Z, divisor_sq):
    num = Z**3 + Abase*Z + Bbase
    q, r = num.quo_rem(divisor_sq)
    assert r == 0
    return q


H1p = quotient_mod_p(Z1p, M1p**2)
H3p = quotient_mod_p(Z3p, P7p**2)

assert H1p == dseed * P7p**2 * u1p**2
assert H3p == dseed * P7p**2 * U3p**2
assert H1p.degree() <= 4
assert H3p.degree() <= 8


def coeffp(poly, i):
    return poly[i] if i <= poly.degree() else F(0)


T1, rem = (3*Z1p**2+Abase).quo_rem(M1p)
assert rem == 0

T3, rem = (3*Z3p**2+Abase).quo_rem(P7p)
assert rem == 0

P7sq = P7p**2

cols = []
# common d
cols.append((-P7sq*u1p**2, -P7sq*U3p**2))
# P1 lambda
cols.append((T1, Rp.zero()))
# P1 scalar u1
cols.append((-2*dseed*P7sq*u1p, Rp.zero()))
# P3 q0,q1,q2
for i in range(3):
    cols.append((Rp.zero(), T3*Vp**i))
# P3 U30,U31
for i in range(2):
    cols.append((Rp.zero(), -2*dseed*P7sq*U3p*Vp**i))

assert len(cols) == 8

Jp = Matrix(
    F,
    14,
    8,
    lambda row, col: (
        coeffp(cols[col][0], row)
        if row < 5
        else coeffp(cols[col][1], row-5)
    ),
)

rank = Jp.rank()

print(
    "Q80TWCOMP|"
    f"jacobian_rows=14|jacobian_cols=8|jacobian_rank={rank}|"
    f"right_kernel_dim={Jp.right_kernel().dimension()}|"
    f"left_kernel_dim={Jp.left_kernel().dimension()}|"
    "stage=COMPONENT_JACOBIAN",
    flush=True,
)

if rank != 8:
    raise ArithmeticError(
        f"component-resolved Jacobian rank={rank}, expected 8"
    )

pivot_rows = tuple(Jp.transpose().pivots())
assert len(pivot_rows) == 8

print(
    f"Q80TWCOMP|pivot_rows={pivot_rows}|stage=PASS_FULL_RANK_8",
    flush=True,
)

# ---------------------------------------------------------------------------
# Q_73 Newton lift of the correctly resolved system.
# ---------------------------------------------------------------------------

PREC = 500
Q73 = Qp(73, PREC, type="capped-rel")
jroots = Q73(-3).sqrt(all=True)
j73 = next(root for root in jroots if int(root.residue()) == 17)


def K_to_Q73(c):
    cc = K0(c)
    return Q73(QQ(cc[0])) + j73*Q73(QQ(cc[1]))


R73 = PolynomialRing(Q73, "V")
V73 = R73.gen()


def poly73(poly):
    return R73([K_to_Q73(c) for c in poly.list()])


A73 = poly73(A0)
B73 = poly73(B0)
P773 = poly73(P7)
node773 = poly73(node7)
M173 = poly73(M1)
interp173 = poly73(interp1)


def coeff73(poly, i):
    return poly[i] if i <= poly.degree() else Q73(0)


def exact_div(num, den):
    q, r = num.quo_rem(den)
    if r != 0:
        vals = [c.valuation() for c in r.list() if c != 0]
        if vals and min(vals) < PREC-30:
            raise ArithmeticError(
                f"unexpected quotient remainder valuation {min(vals)}"
            )
    return q


def state73(z):
    d = z[0]
    lam = z[1]
    u1 = z[2]
    q3 = R73(list(z[3:6]))
    U3 = R73(list(z[6:8])+[Q73(1)])

    Z1 = interp173 + lam*M173
    Z3 = node773 + P773*q3

    H1 = exact_div(Z1**3+A73*Z1+B73, M173**2)
    H3 = exact_div(Z3**3+A73*Z3+B73, P773**2)

    E1 = H1 - d*P773**2*u1**2
    E3 = H3 - d*P773**2*U3**2

    return d, Z1, u1, Z3, U3, E1, E3


def jac73(d, Z1, u1, Z3, U3):
    T1 = exact_div(3*Z1**2+A73, M173)
    T3 = exact_div(3*Z3**2+A73, P773)

    cc = []
    cc.append((-P773**2*u1**2, -P773**2*U3**2))
    cc.append((T1, R73.zero()))
    cc.append((-2*d*P773**2*u1, R73.zero()))

    for i in range(3):
        cc.append((R73.zero(), T3*V73**i))

    for i in range(2):
        cc.append(
            (R73.zero(), -2*d*P773**2*U3*V73**i)
        )

    return Matrix(
        Q73,
        8,
        8,
        lambda rr, col: (
            coeff73(cc[col][0], pivot_rows[rr])
            if pivot_rows[rr] < 5
            else coeff73(cc[col][1], pivot_rows[rr]-5)
        ),
    )


z = vector(Q73, [Q73(v) for v in seed])

for iteration in range(13):
    dcur, Z1cur, u1cur, Z3cur, U3cur, E1cur, E3cur = state73(z)

    residual = vector(
        Q73,
        [
            (
                coeff73(E1cur, row)
                if row < 5
                else coeff73(E3cur, row-5)
            )
            for row in pivot_rows
        ],
    )

    Jcur = jac73(dcur, Z1cur, u1cur, Z3cur, U3cur)
    correction = Jcur.solve_right(residual)
    z -= correction

    _, _, _, _, _, E1check, E3check = state73(z)

    vals = [
        c.valuation()
        for poly in (E1check, E3check)
        for c in poly.list()
        if c != 0
    ]
    minval = min(vals) if vals else PREC

    print(
        f"Q80TWCOMP|newton_iteration={iteration+1}|"
        f"min_residual_v73={minval}|stage=NEWTON",
        flush=True,
    )

    if minval >= 380:
        break

_, _, _, _, _, E1pad, E3pad = state73(z)
vals = [
    c.valuation()
    for poly in (E1pad, E3pad)
    for c in poly.list()
    if c != 0
]
assert not vals or min(vals) >= 300

# ---------------------------------------------------------------------------
# Recognize all 8 values in K=QQ(sqrt(-3)).
# ---------------------------------------------------------------------------


def qq_sqrt(x):
    x = QQ(x)
    if x < 0:
        return None
    num = x.numerator()
    den = x.denominator()
    if not num.is_square() or not den.is_square():
        return None
    return QQ(num.sqrt()) / QQ(den.sqrt())


def recognize_K(alpha, residue):
    dep = alpha.algebraic_dependency(2)

    if dep.degree() == 1:
        val = K0(-QQ(dep[0])/QQ(dep[1]))
        if red_K(val) == F(residue):
            return val

    if dep.degree() != 2:
        raise ArithmeticError(
            f"unexpected algebraic dependency degree {dep.degree()}: {dep}"
        )

    c0, c1, c2 = map(QQ, (dep[0], dep[1], dep[2]))
    a = -c1/(2*c2)
    norm = c0/c2
    b2 = (norm-a**2)/3
    babs = qq_sqrt(b2)

    if babs is None:
        raise ArithmeticError(
            f"algdep is not QQ(sqrt(-3))-shaped: {dep}"
        )

    for b in (babs, -babs):
        val = K0(a+b*j0)
        if red_K(val) == F(residue):
            return val

    raise ArithmeticError(
        f"could not choose Galois root of {dep} matching residue {residue}"
    )


exact = []
for index, (alpha, residue) in enumerate(zip(z, seed)):
    val = recognize_K(alpha, residue)
    exact.append(val)

    print(
        f"Q80TWCOMP|recognized={index}|value={val}|stage=ALGDEP",
        flush=True,
    )

d_exact = exact[0]
lam_exact = exact[1]
u1_exact = exact[2]
q3_exact = R0(exact[3:6])
U3_exact = R0(exact[6:8]+[K0(1)])

Z1_exact = interp1 + lam_exact*M1
Z3_exact = node7 + P7*q3_exact

W1_exact = P7**2 * Q2 * u1_exact
W3_exact = P7**2 * U3_exact

# Decisive exact identities.
assert Z1_exact**3+A0*Z1_exact+B0 == d_exact*W1_exact**2
assert Z3_exact**3+A0*Z3_exact+B0 == d_exact*W3_exact**2

A_actual = d_exact**2*A0
B_actual = d_exact**3*B0
Delta_actual = d_exact**6*Delta0

P1x = d_exact*Z1_exact
P1y = d_exact**2*W1_exact
P3x = d_exact*Z3_exact
P3y = d_exact**2*W3_exact

assert P1y**2 == P1x**3+A_actual*P1x+B_actual
assert P3y**2 == P3x**3+A_actual*P3x+B_actual

dred = red_K(d_exact)
assert (dred/dlocal).is_square()

print(
    "Q80TWCOMP|"
    f"d={d_exact}|d_mod73={int(dred)}|"
    f"P1_degrees={P1x.degree()},{P1y.degree()}|"
    f"P3_degrees={P3x.degree()},{P3y.degree()}|"
    "status=PASS_EXACT_TWIST_P1_P3",
    flush=True,
)

# ---------------------------------------------------------------------------
# Persist exact model and sections.
# ---------------------------------------------------------------------------

OUT_MODEL.write_text(
    "\n".join(
        [
            "#!/usr/bin/env sage",
            "from sage.all import PolynomialRing, QuadraticField",
            'K = QuadraticField(-3, "j")',
            "j = K.gen()",
            'R = PolynomialRing(K, "V")',
            "V = R.gen()",
            f"d = {d_exact}",
            f"A = {A_actual}",
            f"B = {B_actual}",
            "Delta = -16*(4*A^3+27*B^2)",
            "assert Delta.degree() == 24",
            'print(f"Q80ORBIT1222JACOBIAN|d={d}|field=QQ(sqrt(-3))|status=PASS_EXACT_TWIST")',
        ]
    )
    + "\n"
)

OUT_SECTIONS.write_text(
    OUT_MODEL.read_text().rstrip()
    + "\n"
    + "\n".join(
        [
            f"P1x = {P1x}",
            f"P1y = {P1y}",
            f"P3x = {P3x}",
            f"P3y = {P3y}",
            "assert P1y^2 == P1x^3+A*P1x+B",
            "assert P3y^2 == P3x^3+A*P3x+B",
            'print("Q80ORBIT1222SECTIONS|P1_height=1/7|P3_height=8/7|status=PASS_EXACT_P1_P3")',
        ]
    )
    + "\n"
)

OUT_NOTE.write_text(
    "# Q80 orbit 1222 — exact Jacobian twist and MW sections\n\n"
    "Status: **PASS_EXACT_TWIST_P1_P3**\n\n"
    "The decisive correction was to impose the full local component "
    "multiplicities at the I7 fibers, not merely node incidence. "
    "For component labels 2 or 5 of an I7 fiber, the y-coordinate has "
    "vanishing order 2 at the Weierstrass node.\n\n"
    "- P1: `W1=P7^2*Q2*u1` (I7 orders 2,2; selected I2 orders 1,1)\n"
    "- P3: `W3=P7^2*U3` (I7 orders 2,2)\n"
    "- the corrected modular coefficient system has full Jacobian rank 8\n"
    f"- exact twist representative: `{d_exact}`\n"
    "- exact P1 and P3 Weierstrass identities verified\n\n"
    "P3 is the exact height-8/7 horizontal needed for q6_7774.\n"
)

print(
    f"Q80TWCOMP|model={OUT_MODEL}|sections={OUT_SECTIONS}|note={OUT_NOTE}",
    flush=True,
)
