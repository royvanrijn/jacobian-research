#!/usr/bin/env sage
"""
Exact characteristic-zero q4_6855 compiler.

This deliberately replays the already-certified modular construction instead
of searching again:

  1. polynomialize the exact q4_1938 child at the I4 fibre reducing to S=37;
  2. reconstruct the optimal basis section P1 deterministically from the five
     reducible-fibre nodes;
  3. form H=2P1 by the exact group law;
  4. derive the connected-A3 middle-double module at the I4 fibre reducing to
     S=23 from the analytic-centre value and first jet;
  5. certify the 4D -> 2D RR kernel;
  6. compile the binary quartic/Jacobian and regress it to the pinned GF(73)
     q4_6855 child.

No p-adic lifting, algdep, or toric chart scan is used.
"""

from pathlib import Path
import json

from sage.all import (
    QQ, GF, PolynomialRing, QuadraticField, matrix, vector, sage_eval
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
ART = ROOT / "artifacts" / "generated-results" / "q80-q4-1938-char0-resolved-rr.json"
OUT = ROOT / "artifacts" / "generated-results"
QDIR = ROOT / "elkies-k3" / "data" / "fibrations" / "q80-q4-1938-char0"
OUT.mkdir(parents=True, exist_ok=True)
QDIR.mkdir(parents=True, exist_ok=True)

if not ART.exists():
    raise SystemExit(f"missing exact q4_1938 artifact: {ART}")

data = json.loads(ART.read_text())
assert data["status"] == "PASS_EXACT_Q4_1938_RESOLVED_RR"

load(str(HERE / "elliptic_neighbor_compiler.sage"))
load(str(HERE / "elliptic_neighbor_compiler_field_generic.sage"))

K = QuadraticField(-3, "j")
j = K.gen()
R = PolynomialRing(K, "S")
S = R.gen()
KF = R.fraction_field()
Sf = KF(S)

A_rat = KF(sage_eval(data["child"]["jacobian_A"], locals={"j": j, "S": Sf}))
B_rat = KF(sage_eval(data["child"]["jacobian_B"], locals={"j": j, "S": Sf}))

Fp = GF(73)
JMOD = Fp(17)
assert JMOD^2 == Fp(-3)

def red_q(c):
    c = QQ(c)
    return Fp(c.numerator()) / Fp(c.denominator())

def red_k(c):
    c = K(c)
    cc = list(c) + [QQ(0), QQ(0)]
    return red_q(cc[0]) + JMOD*red_q(cc[1])

Rp = PolynomialRing(Fp, "S")
sp = Rp.gen()

def red_poly(f):
    f = R(f)
    return Rp([red_k(c) for c in f.list()])

# ---------------------------------------------------------------------------
# 1. Exact polynomial q4_1938 parent gauge.
# ---------------------------------------------------------------------------
reducible = {}
for item in data["child"]["finite_fibres"]:
    symbol = str(item["kodaira"])
    if symbol not in ("I5", "I4", "I2") or int(item["degree"]) != 1:
        continue
    fac = R(sage_eval(item["factor"], locals={"j":j, "S":S}))
    root = K(-fac[0]/fac[1])
    reducible[int(red_k(root))] = (symbol, root)

expected_reducible = {
    8:"I5", 72:"I5", 23:"I4", 37:"I4", 27:"I2"
}
assert {r:s for r,(s,_) in reducible.items()} == expected_reducible

r37 = reducible[37][1]
den_scale = R((S-r37)^2)
A_parent = R(A_rat * KF(den_scale)^4)
B_parent = R(B_rat * KF(den_scale)^6)
Delta_parent = R(-16*(4*A_parent^3+27*B_parent^2))
assert (A_parent.degree(), B_parent.degree(), Delta_parent.degree()) == (8,12,24)

A_pin = Rp([38,22,48,52,17,58,15,32,34])
B_pin = Rp([38,43,27,37,32,42,23,53,64,3,5,6,11])

gauge_candidates = [
    g for g in Fp if g and
    red_poly(A_parent) == g^4*A_pin and
    red_poly(B_parent) == g^6*B_pin
]
if not gauge_candidates:
    raise ArithmeticError("exact polynomialized 1938 parent misses pinned GF73 gauge")

print(
    "Q806855PARENT|den_scale={}|r37={}|r37_mod73=37|gauges={}|"
    "status=PASS_EXACT_1938_POLYNOMIAL_GAUGE".format(
        den_scale, r37, tuple(int(g) for g in gauge_candidates)
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Deterministic exact reconstruction of optimal P1.
#
# P1 has nonidentity component at every reducible fibre.  Hence on the
# singular polynomial Weierstrass model it passes through all five nodes.
# Degree x(P1)<=4 therefore makes x(P1) the unique node interpolation.
# Its y-coordinate is then forced by the exact curve equation.
# ---------------------------------------------------------------------------
def node_x(root):
    ar = K(A_parent(root))
    br = K(B_parent(root))
    if ar == 0:
        raise ArithmeticError("unexpected c4=0 at old multiplicative fibre")
    x0 = K(-3*br/(2*ar))
    assert x0^3 + ar*x0 + br == 0
    assert 3*x0^2 + ar == 0
    return x0

roots = {res: pair[1] for res,pair in reducible.items()}
node_points = [(roots[res], node_x(roots[res])) for res in (8,72,23,37,27)]
P1x = R.lagrange_polynomial(node_points)
assert P1x.degree() <= 4

Fnodes = R.one()
for res in (8,72,23,37,27):
    Fnodes *= S-roots[res]

rhs = R(P1x^3 + A_parent*P1x + B_parent)
W, rem = rhs.quo_rem(Fnodes^2)
assert rem == 0
assert W.degree() == 2

lead = K(W[2])
if not lead.is_square():
    raise ArithmeticError("P1 y-quotient leading coefficient is not a square in K")
b = K(lead.sqrt())
a = K(W[1]/(2*b))
L = R(a+b*S)
assert L^2 == W

P1x_pin = Rp([54,6,47,29,9])
P1y_pin = Rp([7,0,59,26,70,67,30])

selected = None
for g in gauge_candidates:
    if red_poly(P1x) != g^2*P1x_pin:
        continue
    for P1y in (R(Fnodes*L), R(-Fnodes*L)):
        if red_poly(P1y) == g^3*P1y_pin:
            selected = (g,P1y)
            break
    if selected is not None:
        break

if selected is None:
    raise ArithmeticError("deterministic exact P1 does not match pinned modular section")
g73, P1y = selected

assert P1y^2 == P1x^3 + A_parent*P1x + B_parent

print(
    "Q806855P1|gauge={}|P1x={}|P1y={}|"
    "status=PASS_EXACT_DETERMINISTIC_P1".format(int(g73),P1x,P1y),
    flush=True,
)

# Exact H=2P1.
x1 = KF(P1x)
y1 = KF(P1y)
lam = (3*x1^2 + KF(A_parent))/(2*y1)
Hx_f = lam^2 - 2*x1
Hy_f = lam*(x1-Hx_f) - y1
Hx = R(Hx_f)
Hy = R(Hy_f)
assert Hy^2 == Hx^3 + A_parent*Hx + B_parent
assert Hx.degree() <= 4 and Hy.degree() <= 6

Hx_pin = Rp([34,12,62,23,14])
Hy_pin = Rp([51,13,25,65,24,48,26])
assert red_poly(Hx) == g73^2*Hx_pin
assert red_poly(Hy) == g73^3*Hy_pin

print(
    "Q806855H|H=2P1|Hx={}|Hy={}|status=PASS_EXACT_HORIZONTAL_2P1".format(
        Hx,Hy
    ),
    flush=True,
)

# Node-incidence regression of the historical component profile.
expected_h_hits = {8:True,23:False,27:False,37:True,72:True}
for res,expect in expected_h_hits.items():
    root = roots[res]
    hit = bool(Hx(root) == node_x(root) and Hy(root) == 0)
    assert hit == expect

# ---------------------------------------------------------------------------
# 3. Exact connected-A3 middle-double module at I4 reducing to S=23.
#
# The historical module is the value + first-jet of the pure-base chord
# along the analytic nodal centre c(S), where 3*c(S)^2+A(S)=0.
# For the chord through -H:
#
#       m_c(S) = H_y(S)/(c(S)-H_x(S)).
#
# We need only c(r) and c'(r), so no Laurent/toric expansion is necessary.
# ---------------------------------------------------------------------------
r23 = roots[23]
x023 = node_x(r23)
assert Hx(r23) != x023

c1 = K(-A_parent.derivative()(r23)/(6*x023))
den0 = K(x023-Hx(r23))
m0 = K(Hy(r23)/den0)
m1 = K(
    (
        Hy.derivative()(r23)*den0
        - Hy(r23)*(c1-Hx.derivative()(r23))
    ) / den0^2
)

assert red_k(c1) == g73^2*Fp(23)
assert red_k(m0) == g73*Fp(10)
assert red_k(m1) == g73*Fp(9)

C = matrix(K, [
    [1, r23, r23^2, m0],
    [0, 1, 2*r23, m1],
])
assert C.rank() == 2

block = {
    "name":"q4_6855 connected A3 middle-double module",
    "matrix":C,
    "quotient_basis":("value","first_jet"),
    "provenance":(
        "exact analytic nodal-centre value and first jet at the supported I4; "
        "vertical coefficients (-2,-1,-1)"
    ),
}
compiled = compile_resolved_conditions(
    ("1","S","S^2","m"), (block,), complete=True, coefficient_field=K
)
assert compiled["ambient_dimension"] == 4
assert compiled["rank"] == 2
assert compiled["kernel_dimension"] == 2
assert compiled["h0_certified"]

tail = C.matrix_from_columns([2,3])
assert tail.det() != 0
z1 = tail.solve_right(-C.column(0))
z2 = tail.solve_right(-C.column(1))
k1 = vector(K,[1,0,z1[0],z1[1]])
k2 = vector(K,[0,1,z2[0],z2[1]])
kernel = matrix(K,[k1,k2])
assert C*kernel.transpose() == matrix(K,2,2)

expected_k1 = (Fp(1),Fp(0),Fp(38),Fp(41)/g73)
expected_k2 = (Fp(0),Fp(1),Fp(45),Fp(70)/g73)
assert tuple(red_k(x) for x in k1) == expected_k1
assert tuple(red_k(x) for x in k2) == expected_k2

print(
    "Q806855RRA3|root={}|root_mod73=23|c0={}|c1={}|m0={}|m1={}|"
    "status=PASS_EXACT_CONNECTED_A3_MIDDLE_DOUBLE".format(
        r23,x023,c1,m0,m1
    ),
    flush=True,
)
print(
    "Q806855RR|ambient=4|A3_rank=2|rank=2|nullity=2|h0=2|"
    "status=PASS_EXACT_Q4_6855_RESOLVED_RR",
    flush=True,
)
print(
    "Q806855RRKERNEL|k1={}|k2={}|status=PASS_PINNED_KERNEL_LIFT".format(
        tuple(k1),tuple(k2)
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 4. Compile exact q4_6855 child from this exact RR kernel.
# ---------------------------------------------------------------------------
UR = PolynomialRing(K,"U")
U0 = UR.gen()
KU = UR.fraction_field()
RS = PolynomialRing(KU,"S")
ss = RS.gen()
U = KU(U0)

def lift_poly(f):
    return RS([KU(c) for c in R(f).list()])

def chord_coeff(row):
    aa = RS(KU(row[0]) + KU(row[1])*ss + KU(row[2])*ss^2)
    bb = KU(row[3])
    return aa,bb

a0,b0 = chord_coeff(k1)
a1,b1 = chord_coeff(k2)

hop = compile_degree_two_chord_hop(
    RS,U,a0,b0,a1,b1,
    lift_poly(Hx),-lift_poly(Hy),
    lift_poly(A_parent),lift_poly(B_parent),
)
quartic = hop["binary_quartic"]
assert quartic.degree() == 4

A_child = KU(hop["jacobian_a"])
B_child = KU(hop["jacobian_b"])

print(
    "Q806855CHILD|quartic_degree=4|status=PASS_EXACT_BINARY_QUARTIC",
    flush=True,
)

# ---------------------------------------------------------------------------
# 5. Cheap pinned GF(73) regression FIRST.
# ---------------------------------------------------------------------------
FUring = PolynomialRing(Fp,"U")
up = FUring.gen()
FU = FUring.fraction_field()

def red_u_poly(poly):
    poly = UR(poly)
    return FUring([red_k(c) for c in poly.list()])

def red_u(frac):
    frac = KU(frac)
    num = red_u_poly(frac.numerator())
    den = red_u_poly(frac.denominator())
    if den == 0:
        raise ZeroDivisionError("q4_6855 child denominator vanished mod 73")
    return FU(num)/FU(den)

A_child_pin = FUring([
    35,65,22,17,8,8,34,31,33,10,71,52,29,7,35,17
])
B_child_pin = FUring([
    3,39,8,24,18,29,17,35,8,38,58,18,59,60,20,25,19,17,53,69,65,46,5,21,49
])

A73 = red_u(A_child)
B73 = red_u(B_child)
child_scale = None
for h in Fp:
    if not h:
        continue
    if A73 == FU(h^4*A_child_pin) and B73 == FU(h^6*B_child_pin):
        child_scale = int(h)
        break
if child_scale is None:
    raise ArithmeticError("exact q4_6855 child misses pinned GF73 equation")

print(
    "Q806855REDUCTION|scale_marker={}|status=PASS_GF73_Q4_6855".format(
        child_scale
    ),
    flush=True,
)

# ---------------------------------------------------------------------------
# 6. Exact fibre certificate last; this is the potentially expensive part.
# ---------------------------------------------------------------------------
print("Q806855FIBRES|phase=exact_classification|status=RUNNING", flush=True)
classification = classify_finite_short_weierstrass_fibres(UR,A_child,B_child)

finite_totals = {}
finite_dump = []
for item in classification["finite_fibres"]:
    symbol = str(item["kodaira"])
    finite_totals[symbol] = finite_totals.get(symbol,0)+int(item["degree"])
    finite_dump.append({
        "factor":str(item["factor"]),
        "degree":int(item["degree"]),
        "kodaira":symbol,
        "orders":[int(x) for x in item["minimal_orders"]],
    })

expected_finite = {"I0*":2,"I4":2,"I2":1,"I1":2}
assert finite_totals == expected_finite

infty = classification["infinity_boundary"]
infty_orders = tuple(int(x) for x in infty["normalized_orders"])
assert infty_orders[2] == 0

print(
    "Q806855FIBRES|finite={}|infinity=smooth|"
    "status=PASS_EXACT_Q4_6855_FIBRES".format(finite_totals),
    flush=True,
)

# Persist exact parent + sections for the next hop.
parent_file = QDIR / "q80_char0_q4_1938_parent_P1_H6855.sage"
parent_file.write_text(
    "\n".join([
        "#!/usr/bin/env sage",
        "from sage.all import QuadraticField, PolynomialRing",
        'K = QuadraticField(-3, "j")',
        "j = K.gen()",
        'R = PolynomialRing(K, "S")',
        "S = R.gen()",
        f"A = {A_parent}",
        f"B = {B_parent}",
        f"P1x = {P1x}",
        f"P1y = {P1y}",
        f"Hx = {Hx}",
        f"Hy = {Hy}",
        "assert P1y^2 == P1x^3 + A*P1x + B",
        "assert Hy^2 == Hx^3 + A*Hx + B",
        'print("Q806855PARENT|status=PASS_EXACT_PARENT_P1_H")',
    ]) + "\n"
)

payload = {
    "status":"PASS_EXACT_Q4_6855_RESOLVED_RR",
    "field":"QQ(sqrt(-3))",
    "parent":{
        "polynomial_scale":str(den_scale),
        "A":str(A_parent),
        "B":str(B_parent),
        "gf73_gauge":int(g73),
    },
    "horizontal":{
        "identity":"H=2P1",
        "P1x":str(P1x),
        "P1y":str(P1y),
        "Hx":str(Hx),
        "Hy":str(Hy),
    },
    "connected_A3":{
        "root":str(r23),
        "root_mod73":23,
        "analytic_center":str(x023),
        "center_derivative":str(c1),
        "m0":str(m0),
        "m1":str(m1),
        "matrix":[[str(v) for v in row] for row in C.rows()],
    },
    "rr":{
        "ambient_dimension":4,
        "condition_rank":2,
        "kernel_dimension":2,
        "h0":2,
        "kernel":[[str(x) for x in k1],[str(x) for x in k2]],
    },
    "child":{
        "quartic":str(quartic),
        "jacobian_A":str(A_child),
        "jacobian_B":str(B_child),
        "finite_fibres":finite_dump,
        "infinity_orders":[int(x) for x in infty_orders],
        "global_pattern":"2 I0* + 2 I4 + I2 + 2 I1",
        "gf73_scale_marker":int(child_scale),
    },
    "next":"q4_candidate1 characteristic-zero propagation",
}

json_path = OUT / "q80-q4-6855-char0-resolved-rr.json"
json_path.write_text(json.dumps(payload,indent=2,default=int)+"\n")

note_path = QDIR / "Q80_CHAR0_Q4_6855_RESOLVED_RR_CERTIFICATE.md"
note_path.write_text(
    "# Q80 q4_6855 — exact characteristic-zero resolved RR certificate\n\n"
    "Status: **PASS_EXACT_Q4_6855_RESOLVED_RR**\n\n"
    "- exact q4_1938 child polynomialized at the I4 fibre reducing to `S=37`;\n"
    "- optimal `P1` reconstructed deterministically from all five reducible-fibre nodes;\n"
    "- horizontal `H=2P1` verified by the exact group law;\n"
    "- connected A3 middle-double value/jet module at the exact I4 reducing to `S=23`;\n"
    "- RR ambient dimension 4, condition rank 2, exact `h0=2`;\n"
    "- binary quartic/Jacobian compiled from that same exact kernel;\n"
    "- fibres `2 I0* + 2 I4 + I2 + 2 I1`, infinity smooth;\n"
    f"- pinned GF(73) regression, constant scale marker `{child_scale}`.\n\n"
    "Next: q4 candidate1 characteristic-zero propagation.\n"
)

print(
    "Q806855RRFNAL|json={}|parent={}|certificate={}|"
    "status=PASS_EXACT_Q4_6855_RESOLVED_RR".format(
        json_path,parent_file,note_path
    ),
    flush=True,
)
