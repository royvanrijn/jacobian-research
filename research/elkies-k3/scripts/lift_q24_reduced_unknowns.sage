#!/usr/bin/env sage -python
"""
Reduced Hensel exactification of H92 q24.

Input: a high-precision q24 Hensel checkpoint (normally p2048) in which many
coefficients already rational-reconstruct.

Freeze every reconstructed rational coefficient EXACTLY and Hensel-lift only
the unresolved coefficients.  At p2048 this is typically:

    Z: 0 unresolved
    X: 15 unresolved
    Y: 31 unresolved
    total: 46 unknowns

Thus Newton uses a 157 x 46 Jacobian and a 46 x 46 nonsingular minor instead
of the original 157 x 156 / 156 x 156 solve.

The frozen rationals are not trusted blindly: success requires the complete
exact identity over QQ[U] and literal reduction to the certified modular q24
section.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, Zp, matrix, vector


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    h=Path.home()
    candidates += [
        h/"Documents"/"jacobian-research",
        h/"jacobian-research",
        h/"src"/"jacobian-research",
        h/"git"/"jacobian-research",
    ]
    seen=set()
    for c in candidates:
        try:
            c=c.resolve()
        except Exception:
            continue
        if c in seen:
            continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate jacobian-research")


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--input",type=Path)
parser.add_argument("--precision",type=int,default=8192)
parser.add_argument("--output",type=Path)
parser.add_argument("--exact-output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

INPUT=args.input.resolve() if args.input else LOCAL/"q24-direct-hensel-p2048.json"
OUTPUT=args.output.resolve() if args.output else LOCAL/f"q24-reduced-hensel-p{args.precision}.json"
EXACT_OUTPUT=(
    args.exact_output.resolve()
    if args.exact_output
    else LOCAL/"q8-q24-horizontal-section-qq.json"
)
MOD=LOCAL/"q24-degree46-direct-global-mod-100003.json"

q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
CHILD=next((
    p for p in q8_candidates
    if p.exists()
    and "child" in json.loads(p.read_text())
    and "minimal_A_coefficients_low_to_high"
        in json.loads(p.read_text()).get("child",{})
),None)

for path in (INPUT,MOD):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")
if CHILD is None:
    raise SystemExit("No complete corrected q8 D13 child")

seed=json.loads(INPUT.read_text())
mod=json.loads(MOD.read_text())
child=json.loads(CHILD.read_text())

assert seed["status"]=="PASS_Q24_DIRECT_HENSEL"
assert int(seed["jacobian_rank"])==156
assert int(seed["final_residual_valuation"])>=int(seed["precision"])
assert len(seed["residues"])==156
assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert child["status"]=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"

p=ZZ(seed["prime"])
assert p==ZZ(mod["prime"])
seed_prec=int(seed["precision"])
target=int(args.precision)
if target<=seed_prec:
    raise ValueError("target precision must exceed seed precision")

# Layout: Z lower 24, X 53, Y 79.
def parse_partial(values):
    return [None if v is None else QQ(v) for v in values]

Z0=parse_partial(seed["Z"])
X0=parse_partial(seed["X"])
Y0=parse_partial(seed["Y"])
assert len(Z0)==25 and Z0[24]==1
assert len(X0)==53 and len(Y0)==79

# The leading monic Z coefficient is fixed, not one of the 156 Hensel vars.
unknown=[]
for i in range(24):
    if Z0[i] is None:
        unknown.append(("Z",i))
for i in range(53):
    if X0[i] is None:
        unknown.append(("X",i))
for i in range(79):
    if Y0[i] is None:
        unknown.append(("Y",i))

nunk=len(unknown)
if not nunk:
    raise SystemExit("Checkpoint already has every coefficient reconstructed")

residues=[ZZ(v) for v in seed["residues"]]
seed_res={}
for kind,i in unknown:
    off={"Z":0,"X":24,"Y":77}[kind]
    seed_res[(kind,i)]=residues[off+i]

print(
    "Q24REDUCED_SETUP|"
    f"prime={p}|seed_precision={seed_prec}|target_precision={target}|"
    f"unknowns={nunk}|"
    f"Z_unknown={sum(k=='Z' for k,i in unknown)}|"
    f"X_unknown={sum(k=='X' for k,i in unknown)}|"
    f"Y_unknown={sum(k=='Y' for k,i in unknown)}|status=PASS",
    flush=True,
)
print(
    "Q24REDUCED_UNKNOWN|"
    + ",".join(f"{k}{i}" for k,i in unknown)
    + "|status=PASS",
    flush=True,
)

A_QQ=[QQ(v) for v in child["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ=[QQ(v) for v in child["child"]["minimal_B_coefficients_low_to_high"]]


def coeffs(poly,degree):
    d=poly.degree()
    return [
        poly[i] if d>=0 and i<=d else poly.base_ring()(0)
        for i in range(degree+1)
    ]


def build_polys(ring, unknown_values):
    U=ring.gen()
    values={key:unknown_values[j] for j,key in enumerate(unknown)}

    z=[]
    for i in range(24):
        z.append(ring.base_ring()(Z0[i]) if Z0[i] is not None else values[("Z",i)])
    Z=sum(z[i]*U**i for i in range(24))+U**24

    x=[]
    for i in range(53):
        x.append(ring.base_ring()(X0[i]) if X0[i] is not None else values[("X",i)])
    X=sum(x[i]*U**i for i in range(53))

    y=[]
    for i in range(79):
        y.append(ring.base_ring()(Y0[i]) if Y0[i] is not None else values[("Y",i)])
    Y=sum(y[i]*U**i for i in range(79))

    A=sum(ring.base_ring()(v)*U**i for i,v in enumerate(A_QQ))
    B=sum(ring.base_ring()(v)*U**i for i,v in enumerate(B_QQ))
    return U,Z,X,Y,A,B


def residual_and_jacobian(ring, unknown_values):
    U,Z,X,Y,A,B=build_polys(ring,unknown_values)
    F=Y**2-X**3-A*X*Z**4-B*Z**6

    deriv=[]
    for kind,i in unknown:
        if kind=="Z":
            deriv.append((-4*A*X*Z**3-6*B*Z**5)*U**i)
        elif kind=="X":
            deriv.append((-3*X**2-A*Z**4)*U**i)
        elif kind=="Y":
            deriv.append((2*Y)*U**i)
        else:
            raise AssertionError(kind)

    fc=coeffs(F,156)
    dc=[coeffs(d,156) for d in deriv]
    return (
        vector(ring.base_ring(),fc),
        matrix(
            ring.base_ring(),
            157,nunk,
            lambda row,col: dc[col][row],
        ),
    )


# ---------------------------------------------------------------------------
# Find a nonsingular reduced Jacobian minor modulo p.
# ---------------------------------------------------------------------------
finite=GF(p)
RF=PolynomialRing(finite,"U")
seed_mod=[finite(seed_res[key]%p) for key in unknown]
f0,J0=residual_and_jacobian(RF,seed_mod)

# Frozen exact rationals must agree with the original modular branch, hence
# the residual vanishes mod p.
if f0:
    raise ArithmeticError(
        "freezing reconstructed coefficients changed the modular solution; "
        "one of the supposedly exact coefficients is not on the q24 branch"
    )

rank=J0.rank()
print(
    f"Q24REDUCED_JAC|dimensions={J0.dimensions()}|rank={rank}|expected={nunk}|"
    f"status={'PASS' if rank==nunk else 'SINGULAR'}",
    flush=True,
)
if rank!=nunk:
    raise SystemExit(2)

rows=tuple(int(r) for r in J0.transpose().pivots())
assert len(rows)==nunk
minor0=J0.matrix_from_rows(rows)
assert minor0.is_invertible()

# ---------------------------------------------------------------------------
# Reduced p-adic Newton.
# ---------------------------------------------------------------------------
K=Zp(p,prec=target)
RK=PolynomialRing(K,"U")

# Integer representatives from p^seed_prec are valid initial approximations.
uvals=[K(seed_res[key]) for key in unknown]

def valuation_floor(v):
    vals=[x.valuation() for x in v if x]
    return target if not vals else min(vals)

for iteration in range(1,32):
    F,J=residual_and_jacobian(RK,uvals)
    val=valuation_floor(F)
    print(
        f"Q24REDUCED_NEWTON|iteration={iteration}|residual_valuation={val}|"
        f"unknowns={nunk}|status=PASS",
        flush=True,
    )
    if val>=target:
        break

    JJ=J.matrix_from_rows(rows)
    rhs=-vector(K,[F[r] for r in rows])
    delta=JJ.solve_right(rhs)
    uvals=[a+b for a,b in zip(uvals,delta)]
else:
    raise ArithmeticError("reduced Newton did not reach target precision")

Ffinal,Jfinal=residual_and_jacobian(RK,uvals)
final_val=valuation_floor(Ffinal)
assert final_val>=target

# ---------------------------------------------------------------------------
# Rational recognition of only the remaining variables.
# ---------------------------------------------------------------------------
M=p**target

def rr(value):
    try:
        return ZZ(value.lift()).rational_reconstruction(M)
    except (ArithmeticError,ValueError):
        return None

rec=[rr(v) for v in uvals]
complete=all(v is not None for v in rec)

print(
    "Q24REDUCED_RECON|"
    f"resolved={sum(v is not None for v in rec)}/{nunk}|"
    f"modulus_bits={M.nbits()}|complete={int(complete)}|status=PASS",
    flush=True,
)

# Always checkpoint the reduced lift.
payload={
    "schema":"elkies-k3.h92-q24-reduced-hensel.v1",
    "status":"PASS_Q24_REDUCED_HENSEL",
    "prime":int(p),
    "seed_precision":seed_prec,
    "precision":target,
    "unknown_layout":[[k,int(i)] for k,i in unknown],
    "selected_rows":list(rows),
    "jacobian_rank":int(rank),
    "final_residual_valuation":int(final_val),
    "reconstructed":[None if q is None else str(q) for q in rec],
    "residues":[str(ZZ(v.lift())) for v in uvals],
    "source_checkpoint":str(INPUT.relative_to(ROOT)),
}
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUTPUT}",flush=True)

if not complete:
    print(
        "Q24REDUCED_RESULT|"
        f"precision={target}|resolved={sum(v is not None for v in rec)}/{nunk}|"
        "exact=0|status=NEEDS_MORE_REDUCED_PRECISION",
        flush=True,
    )
    raise SystemExit(3)

# Fill exact coefficients.
Zv=list(Z0)
Xv=list(X0)
Yv=list(Y0)
for key,q in zip(unknown,rec):
    kind,i=key
    if kind=="Z":
        Zv[i]=q
    elif kind=="X":
        Xv[i]=q
    else:
        Yv[i]=q

assert all(v is not None for v in Zv+Xv+Yv)

RQ=PolynomialRing(QQ,"U")
Z=RQ(Zv)
X=RQ(Xv)
Y=RQ(Yv)
A=RQ(A_QQ)
B=RQ(B_QQ)

identity=Y**2-X**3-A*X*Z**4-B*Z**6
if identity:
    print(
        f"Q24REDUCED_EXACT|identity=0|degree={identity.degree()}|status=REJECT",
        flush=True,
    )
    raise SystemExit(4)

assert Z.degree()==24 and Z.leading_coefficient()==1
assert X.degree()==52 and Y.degree()==78

# Independent modular replay.
Fmod=GF(p)
Rmod=PolynomialRing(Fmod,"U")

def red(poly):
    out=[]
    for q in poly.list():
        q=QQ(q)
        d=ZZ(q.denominator())
        if d%p==0:
            raise ArithmeticError("exact coefficient denominator hits check prime")
        out.append(Fmod(ZZ(q.numerator()))/Fmod(d))
    return Rmod(out)

sec=mod["section_mod_p"]
Zm=Rmod([Fmod(int(v)) for v in sec["Z_coefficients_low_to_high"]])
Xm=Rmod([Fmod(int(v)) for v in sec["X_coefficients_low_to_high"]])
Ym=Rmod([Fmod(int(v)) for v in sec["Y_coefficients_low_to_high"]])
if red(Z)!=Zm or red(X)!=Xm or red(Y)!=Ym:
    raise ArithmeticError("exact reduced-Hensel solution fails certified mod-p replay")

max_num_bits=max(abs(ZZ(q.numerator())).nbits() for q in list(Z)+list(X)+list(Y))
max_den_bits=max(abs(ZZ(q.denominator())).nbits() for q in list(Z)+list(X)+list(Y))

exact={
    "schema":"elkies-k3.h92-q8-q24-horizontal-section-qq.reduced-hensel.v1",
    "status":"PASS_EXACT_Q24_HORIZONTAL_SECTION",
    "method":"reduced Hensel with already-recognized coefficients frozen exactly",
    "source_checkpoint":str(INPUT.relative_to(ROOT)),
    "profile":{
        "P_dot_O":24,
        "height":"52",
        "D13_local_correction":"0",
        "Z_degree":24,
        "X_degree":52,
        "Y_degree":78,
        "x_degrees":[52,48],
        "y_degrees":[78,72],
    },
    "section":{
        "Z_coefficients_low_to_high":[str(v) for v in Z.list()],
        "X_coefficients_low_to_high":[str(v) for v in X.list()],
        "Y_coefficients_low_to_high":[str(v) for v in Y.list()],
        "x_numerator_coefficients_low_to_high":[str(v) for v in X.list()],
        "x_denominator_coefficients_low_to_high":[str(v) for v in (Z**2).list()],
        "y_numerator_coefficients_low_to_high":[str(v) for v in Y.list()],
        "y_denominator_coefficients_low_to_high":[str(v) for v in (Z**3).list()],
        "x":"X/Z^2",
        "y":"Y/Z^3",
    },
    "verification":{
        "exact_weierstrass_identity":True,
        "reduction_matches_degree46_modular_section":True,
    },
    "reconstruction":{
        "seed_precision":seed_prec,
        "final_precision":target,
        "reduced_unknown_count":nunk,
        "max_numerator_bits":int(max_num_bits),
        "max_denominator_bits":int(max_den_bits),
    },
}
EXACT_OUTPUT.parent.mkdir(parents=True,exist_ok=True)
EXACT_OUTPUT.write_text(json.dumps(exact,indent=2,sort_keys=True)+"\n")
print(f"EXACT_OUTPUT|{EXACT_OUTPUT}",flush=True)
print(
    "Q24REDUCED_RESULT|"
    f"precision={target}|unknowns={nunk}|identity=PASS|modp=PASS|"
    f"max_num_bits={max_num_bits}|max_den_bits={max_den_bits}|"
    "status=PASS_EXACT_Q24",
    flush=True,
)
