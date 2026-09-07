#!/usr/bin/env sage
"""
Fast local quadratic-character probe for the true Q80 orbit-1222 twist.

For one prime p splitting both -6 and -3:
  * reduce the exact residual genus-one pencil;
  * recover one smooth Jacobian fiber with the certified Singular/brnoeth path;
  * compare its Frobenius trace with the exact reconstructed j/fiber model;
  * repeat for both sqrt(-3) embeddings and both sqrt(-6) embeddings.

If traces are equal, the true twist relative to the base model is a square
locally (chi=+1); if traces are opposite, it is a nonsquare (chi=-1).

The two sqrt(-6) embeddings must agree.  This gives two reliable quadratic
characters per rational prime, without polynomial-section searches.
"""

import argparse
import contextlib
import io
import json
import time
from pathlib import Path

from sage.all import EllipticCurve, GF, Matrix, PolynomialRing, QQ
from sage.env import SAGE_SHARE
from sage.interfaces.singular import singular

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
ORBIT_DATA = (
    REPO_ROOT
    / "elkies-k3"
    / "data"
    / "fibrations"
    / "q80-orbit1222-char0"
)
MODULAR_DATA = (
    REPO_ROOT
    / "artifacts"
    / "generated-results"
    / "q80-orbit1222-char0"
    / "modular"
)
SOURCE = REPO_ROOT / "elkies-k3/scripts/derive_q80_third_q12_cm24_pencil.sage"
BASE_MODEL = ORBIT_DATA / "q80_char0_orbit1222_weierstrass.sage"

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--prime", type=int, required=True)
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()

if not SOURCE.exists():
    raise SystemExit(f"missing source pencil: {SOURCE}")
if not BASE_MODEL.exists():
    raise SystemExit(f"missing exact base model: {BASE_MODEL}")

# Load the exact residual pencil, then preserve every object needed from it.
_worker_file = __file__
__file__ = str(SOURCE)
try:
    with contextlib.redirect_stdout(io.StringIO()):
        load(str(SOURCE))
finally:
    __file__ = _worker_file

SourceK = K
SourceL = L
SourceResidual = residual_cubic
SourceOldBase = new_old_base
SourceParameterRing = new_parameter_ring

# Now load the reconstructed characteristic-zero j/fiber model.
_worker_file = __file__
__file__ = str(BASE_MODEL)
try:
    with contextlib.redirect_stdout(io.StringIO()):
        load(str(BASE_MODEL))
finally:
    __file__ = _worker_file

BaseK = K
BaseA = A
BaseB = B

p = args.prime
F = GF(p)
if not F(-6).is_square() or not F(-3).is_square():
    raise SystemExit(f"p={p} does not split both -6 and -3")


def canonical_root(value):
    roots = F(value).sqrt(all=True)
    return min(roots, key=lambda z: int(z))


rs0 = canonical_root(-6)
rj0 = canonical_root(-3)


def pad(values, n, zero):
    values = list(values)
    return values + [zero]*(n-len(values))


def reduce_Q(value):
    value = QQ(value)
    return F(value.numerator()) / F(value.denominator())


def reduce_source_K(value, sroot):
    coeffs = pad(SourceK(value).list(), 2, QQ(0))
    return reduce_Q(coeffs[0]) + sroot*reduce_Q(coeffs[1])


def reduce_source_L(value, sroot, jroot):
    coeffs = pad(SourceL(value).list(), 2, SourceK(0))
    return (
        reduce_source_K(coeffs[0], sroot)
        + jroot*reduce_source_K(coeffs[1], sroot)
    )


def reduce_base_K(value, jroot):
    coeffs = pad(BaseK(value).list(), 2, QQ(0))
    return reduce_Q(coeffs[0]) + jroot*reduce_Q(coeffs[1])


finite_plane = PolynomialRing(F, names=("w", "x"), order="lex")
fw, fx = finite_plane.gens()

singular.set_ring(finite_plane._singular_())
singular.lib("brnoeth.lib")

# Certified Singular 4.4 compatibility patch.
hnoether_source = Path(SAGE_SHARE) / "singular/LIB/hnoether.lib"
text = hnoether_source.read_text()
a = text.index("proc extdevelop (list l, int Exaktheit)")
b = text.index("\nexample\n", a)
patched = text[a:b]
renames = {
    "ideal a=hole(lastrow);": "ideal q80row=hole(lastrow);",
    "else { ideal a=lastrow; }": "else { ideal q80row=lastrow; }",
    "a[Q]=delt;": "q80row[Q]=delt;",
    "a[Q+1]=x;": "q80row[Q+1]=x;",
    "lastrow=zurueck(a);": "lastrow=zurueck(q80row);",
    "else { lastrow=a; }": "else { lastrow=q80row; }",
}
for old, new in renames.items():
    assert patched.count(old) == 1
    patched = patched.replace(old, new)
singular.eval("kill extdevelop;")
singular.eval(patched)


def reduced_pencil(sroot, jroot):
    reduced = []
    for xcoef in SourceResidual.list():
        wp = SourceOldBase(xcoef)
        row = []
        for parameter_coefficient in wp.list():
            parameter_coefficient = SourceParameterRing(parameter_coefficient)
            row.append(
                tuple(
                    reduce_source_L(c, sroot, jroot)
                    for c in parameter_coefficient.list()
                )
            )
        reduced.append(tuple(row))
    return reduced


def equation_at(reduced, v):
    vv = F(v)
    eq = finite_plane(0)
    for xd, wcoeffs in enumerate(reduced):
        coeff = finite_plane(0)
        for wd, vcoeffs in enumerate(wcoeffs):
            cv = F(0)
            for c in reversed(vcoeffs):
                cv = cv*vv + c
            if cv:
                coeff += cv*fw**wd
        eq += coeff*fx**xd
    return eq


def actual_fiber_curve(reduced, v, serial):
    eq = equation_at(reduced, v)
    if eq.degree(fw) != 9 or eq.degree(fx) != 3:
        raise ArithmeticError("degree drop")

    singular.set_ring(finite_plane._singular_())
    sf = eq._singular_()
    tag = f"Q80TWCHAR{serial}"
    curve = tag+"C"
    proj = tag+"P"
    local = tag+"L"
    place = tag+"PLACE"
    ws = tag+"WS"

    singular.eval(
        "printlevel=-1; "
        f"list {curve}=Adj_div({sf.name()}); "
        f"def {proj}={curve}[1][2]; setring {proj};"
    )
    genus = int(singular.eval(f"{curve}[2][2]"))
    if genus != 1:
        raise ArithmeticError(f"genus={genus}")

    singular.eval(
        f"def {local}={curve}[5][1][1]; setring {local}; "
        f"int {tag}LI=size(POINTS)-1; setring {proj}; "
        f"int {place}=0; int {tag}I; "
        f"for ({tag}I=1;{tag}I<=size({curve}[3]);{tag}I={tag}I+1) "
        "{ "
        f"if ({curve}[3][{tag}I][1]==1 && {curve}[3][{tag}I][2]=={tag}LI) "
        f"{{ {place}={tag}I; }} "
        "}"
    )
    pn = int(singular.eval(place))
    if pn <= 0:
        raise ArithmeticError("xi=-6 place not found")

    singular.eval(
        f"list {ws}=Weierstrass({pn},3,{curve}); "
        f"int {tag}S0={ws}[1][1]; int {tag}S1={ws}[1][2]; int {tag}S2={ws}[1][3]; "
        f"poly {tag}XN={ws}[2][2][1]; poly {tag}XD={ws}[2][2][2]; "
        f"poly {tag}YN={ws}[2][3][1]; poly {tag}YD={ws}[2][3][2]; "
        f"poly {tag}F=CHI;"
    )
    semigroup = tuple(int(singular.eval(tag+x)) for x in ("S0","S1","S2"))
    if semigroup != (0,2,3):
        raise ArithmeticError(f"semigroup={semigroup}")

    xn = singular(tag+"XN").sage()
    xd = singular(tag+"XD").sage()
    yn = singular(tag+"YN").sage()
    yd = singular(tag+"YD").sage()
    peq = singular(tag+"F").sage()

    terms = (
        yn**2*xd**3,
        xn*yn*xd**2*yd,
        yn*xd**3*yd,
        -xn**3*yd**2,
        -xn**2*xd*yd**2,
        -xn*xd**2*yd**2,
        -xd**3*yd**2,
    )
    rem = tuple(t.reduce([peq]) for t in terms)
    mons = sorted(set().union(*(r.dict().keys() for r in rem)))
    M = Matrix(F, [[r.dict().get(m, F(0)) for r in rem] for m in mons])
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

    E = EllipticCurve(F, [a1,a2,a3,a4,a6])
    return E


def base_fiber_curve(v, jroot):
    vv = F(v)
    Av = F(0)
    Bv = F(0)
    for c in reversed(BaseA.list()):
        Av = Av*vv + reduce_base_K(c, jroot)
    for c in reversed(BaseB.list()):
        Bv = Bv*vv + reduce_base_K(c, jroot)
    return EllipticCurve(F, [0,0,0,Av,Bv])


def probe_embedding(sroot, jroot, serial_base):
    reduced = reduced_pencil(sroot, jroot)

    # Odd values first: historically avoids a few pathological V=2 fibers.
    values = list(range(1,p,2)) + list(range(2,p,2))
    failures = []

    for offset, v in enumerate(values):
        try:
            Eactual = actual_fiber_curve(reduced, v, serial_base+offset)
            Ebase = base_fiber_curve(v, jroot)

            ja = Eactual.j_invariant()
            jb = Ebase.j_invariant()
            if ja != jb:
                raise ArithmeticError(f"j mismatch actual={ja}, base={jb}")
            if ja in (F(0), F(1728)):
                continue

            ta = int(Eactual.trace_of_frobenius())
            tb = int(Ebase.trace_of_frobenius())
            if ta == 0 or tb == 0:
                continue

            if ta == tb:
                chi = +1
            elif ta == -tb:
                chi = -1
            else:
                raise ArithmeticError(
                    f"same-j traces are not +/-: actual={ta}, base={tb}"
                )

            return {
                "V": int(v),
                "j_value": int(ja),
                "actual_trace": ta,
                "base_trace": tb,
                "chi": chi,
                "failures_before": failures,
            }
        except Exception as exc:
            failures.append(
                [int(v), type(exc).__name__, str(exc).replace("\n"," ")]
            )

    raise ArithmeticError(f"no usable nonzero-trace fiber; failures={failures[:8]}")


started = time.monotonic()
records = {}
serial = 1000

for ssign in (+1,-1):
    for jsign in (+1,-1):
        sroot = ssign*rs0
        jroot = jsign*rj0
        key = f"s{ssign:+d}_j{jsign:+d}"
        rec = probe_embedding(sroot, jroot, serial)
        serial += p + 100
        records[key] = rec
        print(
            f"Q80TWCHAR|p={p}|s={int(sroot)}|j={int(jroot)}|"
            f"V={rec['V']}|jV={rec['j_value']}|"
            f"actual_trace={rec['actual_trace']}|base_trace={rec['base_trace']}|"
            f"chi={rec['chi']}|status=PASS_LOCAL_CHARACTER",
            flush=True,
        )

# Descent consistency: the character must not depend on sqrt(-6).
chi_plus = records["s+1_j+1"]["chi"]
chi_plus_sconj = records["s-1_j+1"]["chi"]
chi_minus = records["s+1_j-1"]["chi"]
chi_minus_sconj = records["s-1_j-1"]["chi"]

if chi_plus != chi_plus_sconj or chi_minus != chi_minus_sconj:
    raise ArithmeticError(
        "twist character depends on sqrt(-6): "
        f"+j {chi_plus}/{chi_plus_sconj}, -j {chi_minus}/{chi_minus_sconj}"
    )

payload = {
    "version": 1,
    "prime": int(p),
    "sqrt_minus_6_canonical": int(rs0),
    "sqrt_minus_3_canonical": int(rj0),
    "chi_plus_j": int(chi_plus),
    "chi_minus_j": int(chi_minus),
    "records": records,
    "seconds": time.monotonic()-started,
}

args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(payload, indent=2, sort_keys=True, default=int)+"\n")

print(
    f"Q80TWCHAR|p={p}|chi_plus={chi_plus}|chi_minus={chi_minus}|"
    f"seconds={payload['seconds']:.3f}|out={args.out}|"
    "status=PASS_TWIST_CHARACTER_PRIME",
    flush=True,
)
