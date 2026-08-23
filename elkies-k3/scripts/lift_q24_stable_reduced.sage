#!/usr/bin/env sage -python
"""
Stable reduced Hensel lift for H92 q24.

Freeze a coefficient over QQ only when its rational reconstruction is IDENTICAL
at two different p-adic precisions (default p512 and p2048).  A one-shot
reconstruction at p2048 is not trusted.

Everything else remains an unknown and is Newton/Hensel lifted from the p2048
checkpoint.  This keeps the system much smaller than 156 variables while
avoiding the false-freeze plateau.

Success still requires:
  * all 157 coefficient equations at target precision;
  * complete rational reconstruction;
  * exact QQ[U] Weierstrass identity;
  * literal reduction to the certified GF(100003) q24 section.
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
parser.add_argument("--witness",type=Path)
parser.add_argument("--seed",type=Path)
parser.add_argument("--precision",type=int,default=8192)
parser.add_argument("--output",type=Path)
parser.add_argument("--exact-output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

WITNESS=args.witness.resolve() if args.witness else LOCAL/"q24-direct-hensel-p512.json"
SEED=args.seed.resolve() if args.seed else LOCAL/"q24-direct-hensel-p2048.json"
OUTPUT=args.output.resolve() if args.output else LOCAL/f"q24-stable-reduced-hensel-p{args.precision}.json"
EXACT_OUTPUT=args.exact_output.resolve() if args.exact_output else LOCAL/"q8-q24-horizontal-section-qq.json"
MOD=LOCAL/"q24-degree46-direct-global-mod-100003.json"

q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
CHILD=next((
    p for p in q8_candidates
    if p.exists()
    and "child" in json.loads(p.read_text())
    and "minimal_A_coefficients_low_to_high" in json.loads(p.read_text()).get("child",{})
),None)

for pth in (WITNESS,SEED,MOD):
    if not pth.exists():
        raise SystemExit(f"Missing prerequisite: {pth}")
if CHILD is None:
    raise SystemExit("Missing exact D13 child")

w=json.loads(WITNESS.read_text())
s=json.loads(SEED.read_text())
mod=json.loads(MOD.read_text())
child=json.loads(CHILD.read_text())

for rec in (w,s):
    assert rec["status"]=="PASS_Q24_DIRECT_HENSEL"
    assert rec["jacobian_rank"]==156
    assert rec["final_residual_valuation"]>=rec["precision"]
    assert len(rec["residues"])==156

p=ZZ(s["prime"])
assert p==ZZ(w["prime"])==ZZ(mod["prime"])
wp=int(w["precision"])
sp=int(s["precision"])
target=int(args.precision)
if not wp < sp < target:
    raise ValueError(f"need witness_precision < seed_precision < target; got {wp} < {sp} < {target}")

def partial(rec):
    Z=[None if x is None else QQ(x) for x in rec["Z"]]
    X=[None if x is None else QQ(x) for x in rec["X"]]
    Y=[None if x is None else QQ(x) for x in rec["Y"]]
    assert len(Z)==25 and len(X)==53 and len(Y)==79 and Z[24]==1
    return Z,X,Y

Zw,Xw,Yw=partial(w)
Zs,Xs,Ys=partial(s)

# Only freeze if the SAME rational appeared independently at both precisions.
frozen={}
unknown=[]
for kind,lo,hi,a,b in (
    ("Z",0,24,Zw,Zs),
    ("X",0,53,Xw,Xs),
    ("Y",0,79,Yw,Ys),
):
    for i in range(lo,hi):
        if a[i] is not None and b[i] is not None and a[i]==b[i]:
            frozen[(kind,i)]=a[i]
        else:
            unknown.append((kind,i))

nunk=len(unknown)
print(
    "Q24STABLE_SETUP|"
    f"prime={p}|witness_precision={wp}|seed_precision={sp}|target_precision={target}|"
    f"frozen={len(frozen)}|unknowns={nunk}|"
    f"Z_unknown={sum(k=='Z' for k,i in unknown)}|"
    f"X_unknown={sum(k=='X' for k,i in unknown)}|"
    f"Y_unknown={sum(k=='Y' for k,i in unknown)}|status=PASS",
    flush=True,
)
print(
    "Q24STABLE_UNKNOWN|"+",".join(f"{k}{i}" for k,i in unknown)+"|status=PASS",
    flush=True,
)

seed_res=[ZZ(x) for x in s["residues"]]
offset={"Z":0,"X":24,"Y":77}
ures={(k,i):seed_res[offset[k]+i] for k,i in unknown}

A_QQ=[QQ(v) for v in child["child"]["minimal_A_coefficients_low_to_high"]]
B_QQ=[QQ(v) for v in child["child"]["minimal_B_coefficients_low_to_high"]]

def coeffs(poly,n):
    d=poly.degree()
    return [poly[i] if d>=0 and i<=d else poly.base_ring()(0) for i in range(n+1)]

def fixed_or_unknown(base_ring,key,unknown_values,idxmap):
    if key in frozen:
        return base_ring(frozen[key])
    return unknown_values[idxmap[key]]

idxmap={key:j for j,key in enumerate(unknown)}

def build(ring,uvals):
    br=ring.base_ring()
    U=ring.gen()

    z=[fixed_or_unknown(br,("Z",i),uvals,idxmap) for i in range(24)]
    Z=sum(z[i]*U**i for i in range(24))+U**24

    x=[fixed_or_unknown(br,("X",i),uvals,idxmap) for i in range(53)]
    X=sum(x[i]*U**i for i in range(53))

    y=[fixed_or_unknown(br,("Y",i),uvals,idxmap) for i in range(79)]
    Y=sum(y[i]*U**i for i in range(79))

    A=sum(br(v)*U**i for i,v in enumerate(A_QQ))
    B=sum(br(v)*U**i for i,v in enumerate(B_QQ))
    return U,Z,X,Y,A,B

def residual_jac(ring,uvals):
    U,Z,X,Y,A,B=build(ring,uvals)
    F=Y**2-X**3-A*X*Z**4-B*Z**6
    deriv=[]
    for k,i in unknown:
        if k=="Z":
            deriv.append((-4*A*X*Z**3-6*B*Z**5)*U**i)
        elif k=="X":
            deriv.append((-3*X**2-A*Z**4)*U**i)
        else:
            deriv.append(2*Y*U**i)
    fc=coeffs(F,156)
    dc=[coeffs(d,156) for d in deriv]
    return vector(ring.base_ring(),fc), matrix(
        ring.base_ring(),157,nunk,
        lambda row,col: dc[col][row]
    )

# Check the frozen values are compatible modulo p.
Fp=GF(p)
RF=PolynomialRing(Fp,"U")
u0=[Fp(ures[key]%p) for key in unknown]
f0,J0=residual_jac(RF,u0)
if f0:
    bad=[i for i,v in enumerate(f0) if v]
    raise ArithmeticError(
        f"stable frozen values already break the modular branch; bad rows={bad[:12]}"
    )

rank=J0.rank()
print(
    f"Q24STABLE_JAC|dimensions={J0.dimensions()}|rank={rank}|expected={nunk}|"
    f"status={'PASS' if rank==nunk else 'SINGULAR'}",
    flush=True,
)
if rank!=nunk:
    raise SystemExit(2)

rows=tuple(int(r) for r in J0.transpose().pivots())
assert len(rows)==nunk
assert J0.matrix_from_rows(rows).is_invertible()

K=Zp(p,prec=target)
RK=PolynomialRing(K,"U")
uvals=[K(ures[key]) for key in unknown]

def valuation_floor(vals):
    vv=[x.valuation() for x in vals if x]
    return target if not vv else min(vv)

for iteration in range(1,24):
    F,J=residual_jac(RK,uvals)
    full_val=valuation_floor(F)
    selected_F=vector(K,[F[r] for r in rows])
    sel_val=valuation_floor(selected_F)

    print(
        "Q24STABLE_NEWTON|"
        f"iteration={iteration}|selected_residual={sel_val}|full_residual={full_val}|"
        f"unknowns={nunk}|status=PASS",
        flush=True,
    )
    if full_val>=target:
        break

    JJ=J.matrix_from_rows(rows)
    delta=JJ.solve_right(-selected_F)
    dval=valuation_floor(delta)
    print(
        f"Q24STABLE_DELTA|iteration={iteration}|valuation={dval}|status=PASS",
        flush=True,
    )
    uvals=[a+b for a,b in zip(uvals,delta)]
else:
    raise ArithmeticError("stable reduced Newton did not reach target")

Ffinal,Jfinal=residual_jac(RK,uvals)
final_val=valuation_floor(Ffinal)
assert final_val>=target

M=p**target
def rr(v):
    try:
        return ZZ(v.lift()).rational_reconstruction(M)
    except (ArithmeticError,ValueError):
        return None

rec=[rr(v) for v in uvals]
resolved=sum(q is not None for q in rec)
print(
    f"Q24STABLE_RECON|resolved={resolved}/{nunk}|modulus_bits={M.nbits()}|"
    f"complete={int(resolved==nunk)}|status=PASS",
    flush=True,
)

payload={
    "schema":"elkies-k3.h92-q24-stable-reduced-hensel.v1",
    "status":"PASS_Q24_STABLE_REDUCED_HENSEL",
    "prime":int(p),
    "witness_precision":wp,
    "seed_precision":sp,
    "precision":target,
    "frozen_count":len(frozen),
    "unknown_layout":[[k,int(i)] for k,i in unknown],
    "selected_rows":list(rows),
    "final_residual_valuation":int(final_val),
    "reconstructed":[None if q is None else str(q) for q in rec],
    "residues":[str(ZZ(v.lift())) for v in uvals],
}
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUTPUT}",flush=True)

if resolved!=nunk:
    print(
        f"Q24STABLE_RESULT|precision={target}|resolved={resolved}/{nunk}|"
        "status=NEEDS_MORE_STABLE_PRECISION",
        flush=True,
    )
    raise SystemExit(3)

# Assemble final exact section.
vals=dict(frozen)
for key,q in zip(unknown,rec):
    vals[key]=q

RQ=PolynomialRing(QQ,"U")
Z=RQ([vals[("Z",i)] for i in range(24)]+[QQ.one()])
X=RQ([vals[("X",i)] for i in range(53)])
Y=RQ([vals[("Y",i)] for i in range(79)])
A=RQ(A_QQ)
B=RQ(B_QQ)

identity=Y**2-X**3-A*X*Z**4-B*Z**6
if identity:
    print(
        f"Q24STABLE_EXACT|identity=0|degree={identity.degree()}|status=REJECT",
        flush=True,
    )
    raise SystemExit(4)

# Independent modular replay.
sec=mod["section_mod_p"]
Rmod=PolynomialRing(Fp,"U")
def red(poly):
    return Rmod([
        Fp(ZZ(QQ(q).numerator()))/Fp(ZZ(QQ(q).denominator()))
        for q in poly.list()
    ])
Zm=Rmod([Fp(int(v)) for v in sec["Z_coefficients_low_to_high"]])
Xm=Rmod([Fp(int(v)) for v in sec["X_coefficients_low_to_high"]])
Ym=Rmod([Fp(int(v)) for v in sec["Y_coefficients_low_to_high"]])
assert red(Z)==Zm and red(X)==Xm and red(Y)==Ym

max_num_bits=max(abs(ZZ(q.numerator())).nbits() for q in list(Z)+list(X)+list(Y))
max_den_bits=max(abs(ZZ(q.denominator())).nbits() for q in list(Z)+list(X)+list(Y))

exact={
    "schema":"elkies-k3.h92-q8-q24-horizontal-section-qq.stable-reduced.v1",
    "status":"PASS_EXACT_Q24_HORIZONTAL_SECTION",
    "method":"stable-cross-precision reduced Hensel",
    "profile":{
        "P_dot_O":24,"height":"52","D13_local_correction":"0",
        "Z_degree":24,"X_degree":52,"Y_degree":78,
        "x_degrees":[52,48],"y_degrees":[78,72],
    },
    "section":{
        "Z_coefficients_low_to_high":[str(v) for v in Z.list()],
        "X_coefficients_low_to_high":[str(v) for v in X.list()],
        "Y_coefficients_low_to_high":[str(v) for v in Y.list()],
        "x_numerator_coefficients_low_to_high":[str(v) for v in X.list()],
        "x_denominator_coefficients_low_to_high":[str(v) for v in (Z**2).list()],
        "y_numerator_coefficients_low_to_high":[str(v) for v in Y.list()],
        "y_denominator_coefficients_low_to_high":[str(v) for v in (Z**3).list()],
    },
    "verification":{
        "exact_weierstrass_identity":True,
        "reduction_matches_degree46_modular_section":True,
    },
    "reconstruction":{
        "witness_precision":wp,
        "seed_precision":sp,
        "final_precision":target,
        "frozen_stable_coefficients":len(frozen),
        "lifted_coefficients":nunk,
        "max_numerator_bits":int(max_num_bits),
        "max_denominator_bits":int(max_den_bits),
    },
}
EXACT_OUTPUT.parent.mkdir(parents=True,exist_ok=True)
EXACT_OUTPUT.write_text(json.dumps(exact,indent=2,sort_keys=True)+"\n")
print(f"EXACT_OUTPUT|{EXACT_OUTPUT}",flush=True)
print(
    "Q24STABLE_RESULT|identity=PASS|modp=PASS|"
    f"frozen={len(frozen)}|lifted={nunk}|"
    f"max_num_bits={max_num_bits}|max_den_bits={max_den_bits}|"
    "status=PASS_EXACT_Q24",
    flush=True,
)
