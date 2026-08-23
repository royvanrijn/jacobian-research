#!/usr/bin/env sage
"""
Finite-field global-j probe for the Q80 orbit-1222 CM24 pencil.

No characteristic-zero normalization/maximal-order/RR is used.
For one split prime and one embedding of s=sqrt(-6), j=sqrt(-3):
  * reduce the exact residual pencil;
  * compute enough bounded finite-field fibers using Singular brnoeth;
  * interpolate the bidegree-(24,24) j-function;
  * validate on withheld fibers;
  * at p=73, (++), require exact equality with the pinned certificate.
"""

import argparse
import contextlib
import io
import json
import time
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, QQ
from sage.env import SAGE_SHARE
from sage.interfaces.singular import singular

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
MODULAR_DATA = (
    REPO_ROOT
    / "artifacts"
    / "generated-results"
    / "q80-orbit1222-char0"
    / "modular"
)
SOURCE = REPO_ROOT / "elkies-k3/scripts/derive_q80_third_q12_cm24_pencil.sage"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--sign-s", type=int, choices=(-1, 1), default=1)
parser.add_argument("--sign-j", type=int, choices=(-1, 1), default=1)
parser.add_argument("--out", type=Path, required=True)
parser.add_argument("--needed", type=int, default=56)
args = parser.parse_args()

if not SOURCE.exists():
    raise SystemExit(f"missing source: {SOURCE}")

_worker_file = __file__
__file__ = str(SOURCE)
try:
    with contextlib.redirect_stdout(io.StringIO()):
        load(str(SOURCE))
finally:
    __file__ = _worker_file

p = args.prime
finite = GF(p)
if not finite(-6).is_square() or not finite(-3).is_square():
    raise SystemExit(f"p={p} does not split both -6 and -3")

def canonical_root(value):
    r = finite(value).sqrt()
    ri = int(r)
    other = (-ri) % p
    return finite(min(ri, other))

finite_s = args.sign_s * canonical_root(-6)
finite_j = args.sign_j * canonical_root(-3)
assert finite_s**2 == -6 and finite_j**2 == -3

def pad(values, n, zero):
    values = list(values)
    return values + [zero] * (n-len(values))

def reduce_Q(value):
    value = QQ(value)
    return finite(value.numerator()) / finite(value.denominator())

def reduce_K(value):
    c = pad(K(value).list(), 2, QQ(0))
    return reduce_Q(c[0]) + finite_s*reduce_Q(c[1])

def reduce_L(value):
    c = pad(L(value).list(), 2, K(0))
    return reduce_K(c[0]) + finite_j*reduce_K(c[1])

# Pre-reduce the entire exact residual pencil to GF(p)[V] once.
RV = PolynomialRing(finite, "V")
V = RV.gen()
reduced = []
for x_coefficient in residual_cubic.list():
    w_poly = new_old_base(x_coefficient)
    row = []
    for parameter_coefficient in w_poly.list():
        parameter_coefficient = new_parameter_ring(parameter_coefficient)
        row.append(RV([reduce_L(c) for c in parameter_coefficient.list()]))
    reduced.append(tuple(row))

finite_plane = PolynomialRing(finite, names=("w", "x"), order="lex")
finite_w, finite_x = finite_plane.gens()

# Load/patch brnoeth once for the whole worker.
singular.set_ring(finite_plane._singular_())
singular.lib("brnoeth.lib")
hnoether_source = Path(SAGE_SHARE) / "singular/LIB/hnoether.lib"
hnoether_text = hnoether_source.read_text()
procedure_start = hnoether_text.index("proc extdevelop (list l, int Exaktheit)")
procedure_end = hnoether_text.index("\nexample\n", procedure_start)
patched_extdevelop = hnoether_text[procedure_start:procedure_end]
renames = {
    "ideal a=hole(lastrow);": "ideal q80row=hole(lastrow);",
    "else { ideal a=lastrow; }": "else { ideal q80row=lastrow; }",
    "a[Q]=delt;": "q80row[Q]=delt;",
    "a[Q+1]=x;": "q80row[Q+1]=x;",
    "lastrow=zurueck(a);": "lastrow=zurueck(q80row);",
    "else { lastrow=a; }": "else { lastrow=q80row; }",
}
for old, new in renames.items():
    assert patched_extdevelop.count(old) == 1
    patched_extdevelop = patched_extdevelop.replace(old, new)
singular.eval("kill extdevelop;")
singular.eval(patched_extdevelop)

def equation_at(v):
    vv = finite(v)
    F = finite_plane(0)
    for xd, w_coefficients in enumerate(reduced):
        coeff = finite_plane(0)
        for wd, parameter_poly in enumerate(w_coefficients):
            c = parameter_poly(vv)
            if c:
                coeff += c*finite_w**wd
        F += coeff*finite_x**xd
    return F

def fiber_j(v, serial):
    F = equation_at(v)
    if F.degree(finite_w) != 9 or F.degree(finite_x) != 3:
        raise ArithmeticError("degree drop")

    singular.set_ring(finite_plane._singular_())
    sf = F._singular_()
    tag = f"Q80M{serial}"
    curve = tag+"C"
    proj = tag+"P"
    local = tag+"L"
    place_name = tag+"PLACE"
    ws = tag+"WS"

    singular.eval(
        "printlevel=-1; "
        f"list {curve}=Adj_div({sf.name()}); "
        f"def {proj}={curve}[1][2]; setring {proj};"
    )
    genus = int(singular.eval(f"{curve}[2][2]"))
    if genus != 1:
        raise ArithmeticError(f"genus={genus}")

    # Same canonical infinity-branch selection as the pinned p=73 certificate:
    # penultimate local branch at the projective point at infinity = xi=-6.
    singular.eval(
        f"def {local}={curve}[5][1][1]; setring {local}; "
        f"int {tag}LI=size(POINTS)-1; setring {proj}; "
        f"int {place_name}=0; int {tag}I; "
        f"for ({tag}I=1;{tag}I<=size({curve}[3]);{tag}I={tag}I+1) "
        "{ "
        f"if ({curve}[3][{tag}I][1]==1 && "
        f"{curve}[3][{tag}I][2]=={tag}LI) "
        f"{{ {place_name}={tag}I; }} "
        "}"
    )
    place = int(singular.eval(place_name))
    if place <= 0:
        raise ArithmeticError("xi=-6 place not found")

    singular.eval(
        f"list {ws}=Weierstrass({place},3,{curve}); "
        f"int {tag}S0={ws}[1][1]; int {tag}S1={ws}[1][2]; int {tag}S2={ws}[1][3]; "
        f"poly {tag}XN={ws}[2][2][1]; poly {tag}XD={ws}[2][2][2]; "
        f"poly {tag}YN={ws}[2][3][1]; poly {tag}YD={ws}[2][3][2]; "
        f"poly {tag}F=CHI;"
    )
    semigroup = tuple(
        int(singular.eval(tag+x)) for x in ("S0", "S1", "S2")
    )
    if semigroup != (0,2,3):
        raise ArithmeticError(f"semigroup={semigroup}")

    xn = singular(tag+"XN").sage()
    xd = singular(tag+"XD").sage()
    yn = singular(tag+"YN").sage()
    yd = singular(tag+"YD").sage()
    projective_equation = singular(tag+"F").sage()

    relation_terms = (
        yn**2 * xd**3,
        xn * yn * xd**2 * yd,
        yn * xd**3 * yd,
        -xn**3 * yd**2,
        -xn**2 * xd * yd**2,
        -xn * xd**2 * yd**2,
        -xd**3 * yd**2,
    )
    remainders = tuple(
        term.reduce([projective_equation]) for term in relation_terms
    )
    monomials = sorted(
        set().union(*(remainder.dict().keys() for remainder in remainders))
    )
    M = Matrix(
        finite,
        [
            [r.dict().get(monomial, finite(0)) for r in remainders]
            for monomial in monomials
        ],
    )
    ker = M.right_kernel()
    if ker.dimension() != 1:
        raise ArithmeticError(f"relation kernel={ker.dimension()}")
    relation = ker.basis()[0]
    if relation[0] == 0 or relation[3] == 0:
        raise ArithmeticError("bad Weierstrass relation")
    relation /= relation[0]

    c = relation[3]
    a1 = relation[1]
    a2 = relation[4]
    a3 = relation[2]*c
    a4 = relation[5]*c
    a6 = relation[6]*c**2

    b2 = a1**2 + 4*a2
    b4 = a1*a3 + 2*a4
    b6 = a3**2 + 4*a6
    c4 = b2**2 - 24*b4
    c6 = -b2**3 + 36*b2*b4 - 216*b6
    delta = (c4**3-c6**2)/finite(1728)
    if delta == 0:
        raise ArithmeticError("singular Weierstrass result")
    return c4**3/delta

PINNED_VALUES_73 = (
    3,4,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,
    28,29,30,31,32,33,34,35,37,38,39,40,41,42,44,45,46,47,48,49,
    50,51,52,54,55,57,58,59,60,61,62,63,64,65
)
EXPECTED_KERNEL_73 = (
    1,53,22,17,58,36,37,57,32,17,6,49,58,33,23,30,50,51,58,33,2,25,71,1,64,
    1,29,27,37,20,35,39,40,4,9,44,62,1,3,62,50,50,40,10,4,62,68,7,3,10,
)

values = list(PINNED_VALUES_73) if p == 73 else list(range(1, p))
samples = []
failures = []
started = time.monotonic()

for serial, v in enumerate(values, start=1):
    if len(samples) >= args.needed:
        break
    t0 = time.monotonic()
    try:
        jv = fiber_j(v, serial)
        samples.append((finite(v), jv))
        elapsed = time.monotonic()-t0
        print(
            f"Q80MODJ|p={p}|s={int(finite_s)}|j={int(finite_j)}|"
            f"V={v}|jV={int(jv)}|fiber_seconds={elapsed:.3f}|"
            f"good={len(samples)}",
            flush=True,
        )
    except Exception as exc:
        failures.append((v, type(exc).__name__, str(exc)))
        print(
            f"Q80MODJ|p={p}|V={v}|status=SKIP|"
            f"type={type(exc).__name__}|message={str(exc).replace(chr(10),' ')}",
            flush=True,
        )

if len(samples) < 54:
    raise SystemExit(
        f"only {len(samples)} usable fibers at p={p}; need at least 54"
    )

training = samples[:49]
withheld = samples[49:]

def interpolation_row(value, j_value):
    return [value**i for i in range(25)] + [
        -j_value*value**i for i in range(25)
    ]

M = Matrix(finite, [interpolation_row(v,jv) for v,jv in training])
if M.rank() != 49:
    raise ArithmeticError(f"interpolation rank={M.rank()}, expected 49")
kernel_space = M.right_kernel()
if kernel_space.dimension() != 1:
    raise ArithmeticError(f"interpolation nullity={kernel_space.dimension()}")
kernel = kernel_space.basis()[0]
if kernel[0] == 0:
    raise ArithmeticError("cannot normalize j numerator constant coefficient")
kernel /= kernel[0]

N = RV(list(kernel[:25]))
D = RV(list(kernel[25:]))
for v,jv in withheld:
    if D(v) == 0 or N(v) != jv*D(v):
        raise ArithmeticError(f"withheld validation failed at V={v}")

pinned = (
    p == 73 and int(finite_s) == 33 and int(finite_j) == 17
)
if pinned and tuple(int(c) for c in kernel) != EXPECTED_KERNEL_73:
    raise ArithmeticError("p=73 (++): interpolated kernel != pinned certificate")

record = {
    "version": 1,
    "prime": p,
    "s_root": int(finite_s),
    "j_root": int(finite_j),
    "sign_s": args.sign_s,
    "sign_j": args.sign_j,
    "kernel": [int(c) for c in kernel],
    "training": [(int(v),int(jv)) for v,jv in training],
    "withheld": [(int(v),int(jv)) for v,jv in withheld],
    "failures": failures,
    "seconds": time.monotonic()-started,
    "pinned73_match": pinned,
}
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(record, indent=2, sort_keys=True, default=int)+"\n")

print(
    "Q80MODJGLOBAL|"
    f"p={p}|s={int(finite_s)}|j={int(finite_j)}|"
    f"training=49|withheld={len(withheld)}|rank=49|"
    f"seconds={record['seconds']:.3f}|"
    f"pinned73_match={pinned}|status=PASS_GLOBAL_J",
    flush=True,
)
print(f"Q80MODJGLOBAL|out={args.out}", flush=True)
