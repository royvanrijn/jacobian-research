#!/usr/bin/env sage -python
"""
Structured Newton/Hensel lift of the H92 q24 section.

Instead of forming and solving the dense 156x156 coefficient Jacobian, exploit
the polynomial structure of

    F = Y^2 - X^3 - A X Z^4 - B Z^6.

A Newton correction satisfies

    2Y*dY - Cx*dX - Cz*dZ = -F

where
    Cx = 3X^2 + A Z^4,
    Cz = 4 A X Z^3 + 6 B Z^5.

Modulo Y:
    Cx*dX + Cz*dZ = F  (mod Y).

If Cx is invertible modulo Y then
    dX = h + c*dZ (mod Y),
    h = Cx^-1 F,
    c = -Cx^-1 Cz.

Since deg(Y)=78 and deg(dX)<=52, coefficients U^53..U^77 of that residue
must vanish.  That is 25 conditions on only 24 coefficients of dZ.
A fixed rank-24 minor gives a 24x24 p-adic solve.  Then dX is known and dY
comes from exact polynomial division by 2Y.

Thus ALL 156 section coefficients remain live, but Newton linear algebra drops
from 156x156 to 24x24.

Input defaults to q24-direct-hensel-p2048.json.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, PolynomialRing, Qp, matrix, vector


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
parser.add_argument("--seed",type=Path)
parser.add_argument("--precision",type=int,default=8192)
parser.add_argument("--output",type=Path)
parser.add_argument("--exact-output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"

SEED=args.seed.resolve() if args.seed else LOCAL/"q24-direct-hensel-p2048.json"
OUTPUT=args.output.resolve() if args.output else LOCAL/f"q24-structured-hensel-p{args.precision}.json"
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
    and "minimal_A_coefficients_low_to_high"
        in json.loads(p.read_text()).get("child",{})
),None)

for pth in (SEED,MOD):
    if not pth.exists():
        raise SystemExit(f"Missing prerequisite: {pth}")
if CHILD is None:
    raise SystemExit("Missing exact D13 child")

seed=json.loads(SEED.read_text())
mod=json.loads(MOD.read_text())
child=json.loads(CHILD.read_text())

assert seed["status"]=="PASS_Q24_DIRECT_HENSEL"
assert seed["jacobian_rank"]==156
assert seed["final_residual_valuation"]>=seed["precision"]
assert len(seed["residues"])==156
assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"

p=ZZ(seed["prime"])
assert p==ZZ(mod["prime"])
seed_prec=int(seed["precision"])
target=int(args.precision)
if target<=seed_prec:
    raise ValueError("target precision must exceed seed precision")

Aq=[QQ(v) for v in child["child"]["minimal_A_coefficients_low_to_high"]]
Bq=[QQ(v) for v in child["child"]["minimal_B_coefficients_low_to_high"]]
res=[ZZ(v) for v in seed["residues"]]

print(
    "Q24STRUCT_SETUP|"
    f"prime={p}|seed_precision={seed_prec}|target_precision={target}|"
    "newton_system=24x24|status=PASS",
    flush=True,
)


def coeff(poly,i):
    return poly[i] if i<=poly.degree() else poly.base_ring()(0)


def valuation_floor_poly(poly,cap):
    vals=[c.valuation() for c in poly.list() if c]
    return cap if not vals else min(vals)


# ---------------------------------------------------------------------------
# Determine the fixed 24 of 25 high-degree conditions over GF(p).
# ---------------------------------------------------------------------------
Fp=GF(p)
Rp=PolynomialRing(Fp,"U")
Up=Rp.gen()

sec=mod["section_mod_p"]
Zp0=Rp([Fp(int(v)) for v in sec["Z_coefficients_low_to_high"]])
Xp0=Rp([Fp(int(v)) for v in sec["X_coefficients_low_to_high"]])
Yp0=Rp([Fp(int(v)) for v in sec["Y_coefficients_low_to_high"]])
Ap=Rp([Fp(ZZ(q.numerator()))/Fp(ZZ(q.denominator())) for q in Aq])
Bp=Rp([Fp(ZZ(q.numerator()))/Fp(ZZ(q.denominator())) for q in Bq])

assert Zp0.degree()==24 and Xp0.degree()==52 and Yp0.degree()==78
assert Yp0**2-Xp0**3-Ap*Xp0*Zp0**4-Bp*Zp0**6==0

Cxp=3*Xp0**2+Ap*Zp0**4
Czp=4*Ap*Xp0*Zp0**3+6*Bp*Zp0**5

if Cxp.gcd(Yp0).degree()!=0:
    raise ArithmeticError("Cx is not invertible modulo Y at the modular solution")

invCxp=Cxp.inverse_mod(Yp0)
cp=(-invCxp*Czp)%Yp0

# Matrix mapping dZ coefficients to forbidden dX coefficients U^53..U^77.
high_degrees=list(range(53,78))
cols=[]
for j in range(24):
    col=(cp*Up**j)%Yp0
    cols.append([coeff(col,d) for d in high_degrees])

M0=matrix(Fp,25,24,lambda r,c:cols[c][r])
rank=M0.rank()
if rank!=24:
    raise ArithmeticError(f"structured modular matrix rank {rank}, expected 24")

# Pivot columns of transpose = independent rows.
rows=tuple(int(r) for r in M0.transpose().pivots())
assert len(rows)==24
assert M0.matrix_from_rows(rows).is_invertible()

print(
    "Q24STRUCT_MODULAR|"
    f"conditions=25|rank={rank}|selected_rows={','.join(map(str,rows))}|"
    "Cx_unit_mod_Y=1|status=PASS",
    flush=True,
)


# ---------------------------------------------------------------------------
# p-adic structured Newton.
# ---------------------------------------------------------------------------
K=Qp(p,prec=target)
R=PolynomialRing(K,"U")
U=R.gen()

def qpk(q):
    q=QQ(q)
    return K(ZZ(q.numerator()))/K(ZZ(q.denominator()))

A=R([qpk(q) for q in Aq])
B=R([qpk(q) for q in Bq])

def valuation_floor_coeffs(poly,cap):
    vals=[c.valuation() for c in poly.list() if c]
    return cap if not vals else min(vals)

def rem_mod(poly,modulus):
    return poly.quo_rem(modulus)[1]

# Lift the certified GF(p) inverse Cx^{-1} mod Y. It is a valid mod-p
# starting point for every p-adic iterate on this branch.
invCx_hint=R([K(int(c)) for c in invCxp.list()])

def newton_inverse_mod(cx,yy,hint,needed,outer):
    q=rem_mod(hint,yy)
    step=0
    while True:
        err=rem_mod(R.one()-cx*q,yy)
        ev=valuation_floor_coeffs(err,target)
        print(
            "Q24STRUCT_INV|"
            f"outer={outer}|step={step}|valuation={ev}|needed={needed}|status=PASS",
            flush=True,
        )
        if ev>=needed:
            return q
        if ev<=0:
            raise ArithmeticError(
                f"Newton inverse lost modular seed at outer={outer}, valuation={ev}"
            )
        q=rem_mod(q*(R(2)-cx*q),yy)
        step+=1
        if step>20:
            raise ArithmeticError("quotient inverse Newton failed to converge")

# Seed layout: Z lower 24, X 53, Y 79. Z leading coefficient = 1.
Z=R([K(res[i]) for i in range(24)]+[K.one()])
X=R([K(res[24+i]) for i in range(53)])
Y=R([K(res[77+i]) for i in range(79)])

assert Z.degree()==24 and X.degree()==52 and Y.degree()==78

def full_residual(Z,X,Y):
    return Y**2-X**3-A*X*Z**4-B*Z**6

for iteration in range(1,16):
    F=full_residual(Z,X,Y)
    fval=valuation_floor_poly(F,target)
    print(
        f"Q24STRUCT_NEWTON|iteration={iteration}|residual_valuation={fval}|status=PASS",
        flush=True,
    )
    if fval>=target:
        break

    Cx=3*X**2+A*Z**4
    Cz=4*A*X*Z**3+6*B*Z**5

    # Avoid polynomial XGCD over Qp. Refine the known mod-p inverse by
    # Newton in K[U]/(Y). Accuracy fval is sufficient because F=O(p^fval).
    needed=min(target,max(1,fval))
    # A Newton correction from residual O(p^v) only needs its linearized
    # compatibility conditions correct modulo p^(min(target,2v)).
    # Requiring final target precision here is unnecessarily strong.
    expected_next=min(target,2*fval)

    invCx=newton_inverse_mod(Cx,Y,invCx_hint,needed,iteration)
    invCx_hint=invCx

    h=rem_mod(invCx*rem_mod(F,Y),Y)
    c=rem_mod(-invCx*rem_mod(Cz,Y),Y)

    # h + c*dZ must have coefficients 53..77 equal zero.
    cols=[]
    for j in range(24):
        col=rem_mod(c*U**j,Y)
        cols.append([coeff(col,d) for d in high_degrees])

    MM=matrix(K,25,24,lambda r,cc:cols[cc][r])
    rhs=vector(K,[-coeff(h,d) for d in high_degrees])

    Ms=MM.matrix_from_rows(rows)
    bs=vector(K,[rhs[r] for r in rows])
    dzvec=Ms.solve_right(bs)

    # Verify all 25 conditions, including the redundant one.
    defect=MM*dzvec-rhs
    dval=min([v.valuation() for v in defect if v],default=target)
    print(
        "Q24STRUCT_COMPAT|"
        f"iteration={iteration}|valuation={dval}|required={expected_next}|status=PASS",
        flush=True,
    )
    if dval < expected_next:
        raise ArithmeticError(
            f"structured high-degree compatibility only to valuation {dval}; "
            f"need {expected_next}"
        )

    dZ=R(list(dzvec))
    dX=rem_mod(h+c*dZ,Y)

    # The high coefficients should now vanish to working precision.  Remove
    # numerically-zero p-adic coefficients above degree 52 by truncating to
    # the permitted correction space.
    high=[coeff(dX,d) for d in high_degrees]
    highval=min([v.valuation() for v in high if v],default=target)
    if highval < expected_next:
        raise ArithmeticError(
            f"dX high-degree tail valuation {highval}; need {expected_next}"
        )
    dX=R([coeff(dX,i) for i in range(53)])

    # Recover dY from polynomial division:
    #   2Y*dY = -F + Cx*dX + Cz*dZ.
    numer=-F+Cx*dX+Cz*dZ
    q,rem=numer.quo_rem(Y)
    remval=valuation_floor_poly(rem,target) if rem else target
    if remval < expected_next:
        raise ArithmeticError(
            f"dY division remainder valuation {remval}; need {expected_next}"
        )
    dY=q/K(2)
    if dY.degree()>78:
        tail=[coeff(dY,i) for i in range(79,dY.degree()+1)]
        tailval=min([v.valuation() for v in tail if v],default=target)
        if tailval<expected_next:
            raise ArithmeticError(
                f"dY degree overflow valuation {tailval}; need {expected_next}"
            )
        dY=R([coeff(dY,i) for i in range(79)])

    dzval=valuation_floor_poly(dZ,target)
    dxval=valuation_floor_poly(dX,target)
    dyval=valuation_floor_poly(dY,target)
    print(
        "Q24STRUCT_DELTA|"
        f"iteration={iteration}|dZ={dzval}|dX={dxval}|dY={dyval}|status=PASS",
        flush=True,
    )

    Z+=dZ
    X+=dX
    Y+=dY
else:
    raise ArithmeticError("structured Newton did not reach target precision")

Ffinal=full_residual(Z,X,Y)
final_val=valuation_floor_poly(Ffinal,target)
if final_val<target:
    raise ArithmeticError(f"final residual valuation {final_val} < {target}")


# ---------------------------------------------------------------------------
# Save checkpoint and try exact rational reconstruction.
# ---------------------------------------------------------------------------
M=p**target

def rr(v):
    try:
        return ZZ(v.lift()).rational_reconstruction(M)
    except (ArithmeticError,ValueError):
        return None

Zv=[rr(coeff(Z,i)) for i in range(24)]+[QQ.one()]
Xv=[rr(coeff(X,i)) for i in range(53)]
Yv=[rr(coeff(Y,i)) for i in range(79)]

resolved=sum(v is not None for v in Zv[:-1]+Xv+Yv)
complete=resolved==156

payload={
    "schema":"elkies-k3.h92-q24-structured-hensel.v1",
    "status":"PASS_Q24_STRUCTURED_HENSEL",
    "prime":int(p),
    "seed_precision":seed_prec,
    "precision":target,
    "final_residual_valuation":int(final_val),
    "selected_high_rows":list(rows),
    "complete":bool(complete),
    "Z":[None if v is None else str(v) for v in Zv],
    "X":[None if v is None else str(v) for v in Xv],
    "Y":[None if v is None else str(v) for v in Yv],
    "residues":[
        *[str(ZZ(coeff(Z,i).lift())) for i in range(24)],
        *[str(ZZ(coeff(X,i).lift())) for i in range(53)],
        *[str(ZZ(coeff(Y,i).lift())) for i in range(79)],
    ],
    "source_checkpoint":str(SEED.relative_to(ROOT)),
}
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUTPUT}",flush=True)

print(
    f"Q24STRUCT_RECON|resolved={resolved}/156|complete={int(complete)}|"
    f"modulus_bits={M.nbits()}|status=PASS",
    flush=True,
)

if not complete:
    print(
        "Q24STRUCT_RESULT|"
        f"precision={target}|resolved={resolved}/156|"
        "status=NEEDS_MORE_STRUCTURED_PRECISION",
        flush=True,
    )
    raise SystemExit(3)

RQ=PolynomialRing(QQ,"U")
ZQ=RQ(Zv)
XQ=RQ(Xv)
YQ=RQ(Yv)
AQ=RQ(Aq)
BQ=RQ(Bq)
identity=YQ**2-XQ**3-AQ*XQ*ZQ**4-BQ*ZQ**6
if identity:
    raise ArithmeticError(f"reconstructed exact identity nonzero degree {identity.degree()}")

# Literal modular replay.
def red(poly):
    return Rp([
        Fp(ZZ(QQ(q).numerator()))/Fp(ZZ(QQ(q).denominator()))
        for q in poly.list()
    ])
assert red(ZQ)==Zp0 and red(XQ)==Xp0 and red(YQ)==Yp0

max_num_bits=max(abs(ZZ(q.numerator())).nbits() for q in list(ZQ)+list(XQ)+list(YQ))
max_den_bits=max(abs(ZZ(q.denominator())).nbits() for q in list(ZQ)+list(XQ)+list(YQ))

exact={
    "schema":"elkies-k3.h92-q8-q24-horizontal-section-qq.structured-hensel.v1",
    "status":"PASS_EXACT_Q24_HORIZONTAL_SECTION",
    "method":"structured polynomial Newton/Hensel (24x24 corrections)",
    "profile":{
        "P_dot_O":24,"height":"52","D13_local_correction":"0",
        "Z_degree":24,"X_degree":52,"Y_degree":78,
        "x_degrees":[52,48],"y_degrees":[78,72],
    },
    "section":{
        "Z_coefficients_low_to_high":[str(v) for v in ZQ.list()],
        "X_coefficients_low_to_high":[str(v) for v in XQ.list()],
        "Y_coefficients_low_to_high":[str(v) for v in YQ.list()],
        "x_numerator_coefficients_low_to_high":[str(v) for v in XQ.list()],
        "x_denominator_coefficients_low_to_high":[str(v) for v in (ZQ**2).list()],
        "y_numerator_coefficients_low_to_high":[str(v) for v in YQ.list()],
        "y_denominator_coefficients_low_to_high":[str(v) for v in (ZQ**3).list()],
    },
    "verification":{
        "exact_weierstrass_identity":True,
        "reduction_matches_degree46_modular_section":True,
    },
    "reconstruction":{
        "seed_precision":seed_prec,
        "final_precision":target,
        "newton_linear_system":"24x24",
        "max_numerator_bits":int(max_num_bits),
        "max_denominator_bits":int(max_den_bits),
    },
}
EXACT_OUTPUT.parent.mkdir(parents=True,exist_ok=True)
EXACT_OUTPUT.write_text(json.dumps(exact,indent=2,sort_keys=True)+"\n")
print(f"EXACT_OUTPUT|{EXACT_OUTPUT}",flush=True)
print(
    "Q24STRUCT_RESULT|identity=PASS|modp=PASS|"
    f"max_num_bits={max_num_bits}|max_den_bits={max_den_bits}|"
    "status=PASS_EXACT_Q24",
    flush=True,
)
