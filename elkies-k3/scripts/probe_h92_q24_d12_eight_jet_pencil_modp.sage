#!/usr/bin/env sage -python
"""
Provisional end-to-end H3 q24 -> D12 modular compiler.

This is deliberately a DIAGNOSTIC bridge between the passing RR preflight
and the still-required resolved I9* proof.

Inputs already certified:
  * canonical D13 parent over QQ(U), reduced mod p;
  * q24 point over GF(p)(U), with X/Z^2, Y/Z^3 and deg Z=24;
  * D24 = O + P24 - 7F + V, old-fibre degree 2;
  * smooth-collision ambient 58 -> rank 48 -> dimension 10;
  * ordinary cusp jets on that 10-space have ranks 1,2,...,10.

Here we impose the FIRST EIGHT cusp jets as a provisional local quotient.
If this gives a 2D pencil whose binary-quartic Jacobian has D12 root data,
then the next task is sharply defined: prove, using the actual resolved I9*
charts, that their connected quotient cuts out this same 2D subspace.

This script does NOT promote ordinary jets to a resolved-fibre certificate.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


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
        h/"projects"/"jacobian-research",
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
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"
CORE=ROOT/"elkies-k3/scripts/elliptic_neighbor_compiler.sage"
exec(compile(CORE.read_text(),str(CORE),"exec"))

MOD=LOCAL/f"q24-degree46-direct-global-mod-{args.prime}.json"
TRANS=LOCAL/"q8-q24-physical-to-equation-translation.json"
q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((
    path for path in q8_candidates
    if path.exists()
    and json.loads(path.read_text()).get("status")
       =="PASS_EXACT_CORRECTED_Q8_D13_CHILD"
),None)
if Q8 is None:
    raise SystemExit("No passing exact D13 q8 child artifact")
for path in (MOD,TRANS,Q8,CORE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

mod=json.loads(MOD.read_text())
trans=json.loads(TRANS.read_text())
q8=json.loads(Q8.read_text())
assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert trans["status"]=="PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION"

qe=trans["q24_equation"]
assert qe["child_root_data"]==[12,264,4]
assert qe["height"]=="52"
assert qe["D13_local_correction"]=="0"
assert qe["P_dot_O"]==24
assert qe["vertical_fibre_coefficient"]==-7

p=ZZ(args.prime)
F=GF(p)
R=PolynomialRing(F,"U")
U=R.gen()
K=R.fraction_field()

def red_q(q):
    q=QQ(q)
    d=ZZ(q.denominator())
    if d%p==0:
        raise ZeroDivisionError(f"denominator divisible by {p}")
    return F(ZZ(q.numerator()))/F(d)

def red_poly(values):
    return R([red_q(QQ(v)) for v in values])

child=q8["child"]
A=red_poly(child["minimal_A_coefficients_low_to_high"])
B=red_poly(child["minimal_B_coefficients_low_to_high"])
Delta=-16*(4*A**3+27*B**2)

sec=mod["section_mod_p"]
Z=R([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
X=R([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
Y=R([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
assert (Z.degree(),X.degree(),Y.degree())==(24,52,78)
assert Z.is_monic()
assert Y**2==X**3+A*X*Z**4+B*Z**6

xP=K(X)/K(Z**2)
yP=K(Y)/K(Z**3)
assert yP**2==xP**3+K(A)*xP+K(B)

# ---------------------------------------------------------------------------
# 58-dimensional global coefficient ambient.
# ---------------------------------------------------------------------------
Amax=41
Bmax=15
ambient=[("A",i) for i in range(Amax+1)]
ambient += [("B",i) for i in range(Bmax+1)]
assert len(ambient)==58

# Smooth P.O collision:
#   a=Acoef/Z^2, b=Bcoef/Z,
#   m=(y+yP)/(x-xP)
# gives Acoef*X-Bcoef*Y = 0 mod Z^2.
modulus=Z**2
collision_cols=[]
for kind,i in ambient:
    collision_cols.append(
        (U**i*X)%modulus if kind=="A"
        else (-U**i*Y)%modulus
    )
C=matrix(F,48,58,lambda row,col:collision_cols[col][row])
assert C.rank()==48

# ---------------------------------------------------------------------------
# I9* local cusp frame.
# ---------------------------------------------------------------------------
i9=next(item for item in child["finite_fibres"] if item["kodaira"]=="I9*")
RQ=PolynomialRing(QQ,"U")
fQ=RQ(str(i9["factor"]))
f=R([red_q(c) for c in fQ.list()])
assert f.degree()==1
alpha=-f[0]/f[1]

S=PolynomialRing(F,"u")
u=S.gen()
KS=S.fraction_field()

def shift_poly(poly):
    return S(poly(alpha+u))

def shift_rf(value):
    value=K(value)
    num=shift_poly(R(value.numerator()))
    den=shift_poly(R(value.denominator()))
    if not den[0]:
        raise ZeroDivisionError("local denominator not a unit at I9*")
    return KS(num)/KS(den)

Al=shift_poly(A)
Bl=shift_poly(B)
Dl=shift_poly(Delta)
assert (Al.valuation(),Bl.valuation(),Dl.valuation())==(2,3,15)

xloc=shift_rf(xP)
yloc=shift_rf(yP)
assert xloc.valuation()>=0 and yloc.valuation()>=0
x0=F(xloc(0)); y0=F(yloc(0))
assert x0 or y0
assert y0**2==x0**3

# Chord evaluated at the singular cusp x=y=0.
mcusp=-yloc/xloc
assert mcusp.valuation()>=0

def jet_of_rf(value,n):
    value=KS(value)
    num=S(value.numerator())
    den=S(value.denominator())
    if not den[0]:
        raise ZeroDivisionError("jet denominator not a unit")
    rem=(num*den.inverse_mod(u**n))%(u**n)
    return [rem[i] for i in range(n)]

# Build the first 8 cusp-jet rows directly on all 58 ambient columns.
Jcols=[]
for kind,i in ambient:
    if kind=="A":
        g=shift_rf(K(U**i)/K(Z**2))
    else:
        g=shift_rf(K(U**i)/K(Z))*mcusp
    Jcols.append(jet_of_rf(g,8))

J=matrix(F,8,58,lambda row,col:Jcols[col][row])
assert J.rank()==8

M=C.stack(J)
combined_rank=M.rank()
kernel=M.right_kernel().basis_matrix()
print(
    "Q24D12JET_RR|"
    f"ambient=58|collision_rank={C.rank()}|jet_rank={J.rank()}|"
    f"combined_rank={combined_rank}|kernel={kernel.nrows()}|"
    f"status={'PASS_PROVISIONAL_H0_TWO' if combined_rank==56 and kernel.nrows()==2 else 'UNEXPECTED'}",
    flush=True,
)
if combined_rank!=56 or kernel.nrows()!=2:
    raise ArithmeticError("provisional 48+8 quotient does not leave dimension 2")

# ---------------------------------------------------------------------------
# Convert the two kernel rows into global coefficient pairs a+b*m.
# ---------------------------------------------------------------------------
def pair_from_row(row):
    AA=R.zero()
    BB=R.zero()
    for j,(kind,i) in enumerate(ambient):
        if kind=="A":
            AA += row[j]*U**i
        else:
            BB += row[j]*U**i
    a=K(AA)/K(Z**2)
    b=K(BB)/K(Z)
    assert (AA*X-BB*Y)%modulus==0
    # Verify provisional I9* u^8 condition.
    local=shift_rf(a)+shift_rf(b)*mcusp
    assert all(v==0 for v in jet_of_rf(local,8))
    return AA,BB,a,b

pairs=[pair_from_row(row) for row in kernel.rows()]
for idx,(AA,BB,a,b) in enumerate(pairs):
    print(
        "Q24D12JET_PENCILGEN|"
        f"i={idx}|Adeg={AA.degree()}|Bdeg={BB.degree()}|"
        f"Aterms={len([c for c in AA.list() if c])}|"
        f"Bterms={len([c for c in BB.list() if c])}|status=PASS",
        flush=True,
    )

# ---------------------------------------------------------------------------
# Degree-two chord elimination over F_p(V).
# ---------------------------------------------------------------------------
VR=PolynomialRing(F,"V")
V=VR.gen()
VF=VR.fraction_field()
UR=PolynomialRing(VF,"U")
UU=UR.gen()
UK=UR.fraction_field()

def lift_poly(poly):
    poly=R(poly)
    return UR([VF(c) for c in poly.list()])

def lift_rf(value):
    value=K(value)
    return UK(lift_poly(R(value.numerator())))/UK(
        lift_poly(R(value.denominator()))
    )

a0,b0=lift_rf(pairs[0][2]),lift_rf(pairs[0][3])
a1,b1=lift_rf(pairs[1][2]),lift_rf(pairs[1][3])
xPV,yPV=lift_rf(xP),lift_rf(yP)
AV,BV=lift_poly(A),lift_poly(B)

den=b1-VF(V)*b0
if not den:
    raise ArithmeticError("pencil basis has degenerate chord coefficient")
mval=-(a1-VF(V)*a0)/den

XR=PolynomialRing(UK,"x")
xx=XR.gen()
yline=XR(mval)*(xx-XR(xPV))-XR(yPV)
relation=yline**2-xx**3-XR(AV)*xx-XR(BV)
quadratic,remainder=relation.quo_rem(xx-XR(xPV))
assert not remainder and quadratic.degree()==2
disc=UK(quadratic[1]**2-4*quadratic[2]*quadratic[0])
assert disc

quartic,square_factor=squarefree_binary_quartic(disc,UR)
print(
    "Q24D12JET_QUARTIC|"
    f"degree={quartic.degree()}|"
    f"status={'PASS_GENUS_ONE' if quartic.degree() in (3,4) else 'FAIL'}",
    flush=True,
)
if quartic.degree() not in (3,4):
    raise ArithmeticError("provisional pencil did not produce genus one")

I,J=binary_quartic_invariants(quartic)
jacA=VF(-27)*VF(I)
jacB=VF(-27)*VF(J)

# The generic compiler's finite minimizer is ring-generic away from 2,3.
classification=classify_finite_short_weierstrass_fibres(VR,jacA,jacB)
finite=[
    {
        "factor":str(item["factor"]),
        "degree":int(item["degree"]),
        "minimal_orders":list(map(int,item["minimal_orders"])),
        "kodaira":item["kodaira"],
    }
    for item in classification["finite_fibres"]
]
root_rank=int(classification["finite_root_rank"])
euler=int(classification["finite_euler_number"])
root_det=int(classification["finite_root_determinant"])
infinity=classification["infinity_boundary"]
inf_orders=tuple(map(int,infinity["normalized_orders"]))
inf_kind="smooth"
if inf_orders[2]>0:
    ir,ie,idt,inf_kind=kodaira_data_from_short_orders(*inf_orders)
    root_rank += int(ir)
    euler += int(ie)
    root_det *= int(idt)

print(
    "Q24D12JET_CHILD|"
    f"finite={[(r['degree'],r['minimal_orders'],r['kodaira']) for r in finite]}|"
    f"infinity={inf_orders},{inf_kind}|"
    f"root_rank={root_rank}|euler={euler}|root_det={root_det}|"
    f"status={'PASS_PROVISIONAL_D12' if (root_rank,root_det)==(12,4) else 'NOT_D12'}",
    flush=True,
)

is_d12=(root_rank,root_det)==(12,4)
payload={
    "schema":"elkies-k3.h3-q24-d12-jet-pencil-modp.v1",
    "status":(
        "PASS_PROVISIONAL_D12_FROM_EIGHT_CUSP_JETS"
        if is_d12 else
        "PROVISIONAL_EIGHT_CUSP_JETS_NOT_D12"
    ),
    "proof_boundary":(
        "The first eight ordinary cusp jets are a diagnostic candidate for "
        "the missing I9* local quotient.  This artifact is NOT a resolved-"
        "component certificate.  Promotion requires proving that the actual "
        "resolved I9* chart cover cuts out the same 2D global subspace."
    ),
    "prime":int(p),
    "rr":{
        "ambient":58,
        "smooth_collision_rank":48,
        "provisional_i9star_jet_rank":8,
        "combined_rank":int(combined_rank),
        "kernel_dimension":int(kernel.nrows()),
    },
    "quartic_degree":int(quartic.degree()),
    "child":{
        "finite_fibres":finite,
        "infinity_orders":list(inf_orders),
        "infinity_kind":inf_kind,
        "root_rank":int(root_rank),
        "root_determinant":int(root_det),
        "euler":int(euler),
    },
    "next":(
        "Resolve the actual I9* fibre and compare its connected marked-chord "
        "quotient row space with the eight-jet row space."
        if is_d12 else
        "Discard the ordinary eight-jet local quotient and derive the actual "
        "resolved I9* module before further elimination."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q24-d12-jet-pencil-mod-{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24D12JET_RESULT|"
    f"kernel=2|quartic={quartic.degree()}|"
    f"root_rank={root_rank}|root_det={root_det}|"
    f"status={payload['status']}",
    flush=True,
)
