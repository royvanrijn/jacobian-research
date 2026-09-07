#!/usr/bin/env sage -python
"""
Equation-level modular compiler test for the geometrically clean q32 neighbour

    D13/MW4 --q32(current equation frame)--> D12/MW5.

Certified lattice/geometry input:
  * same horizontal P24, height 52, P.O=24, D13 correction 0;
  * old-fibre degree 2;
  * q32 current-frame factor presentation D=[16,2,w], D.O=14;
  * local vertical correction at the ACTUAL resolved I9* fibre, after assigning
    -8F to infinity:
        E01:-2 E02:-4 E03:-3 E04:-6 E05:-5 E06:-8
        E07:-7 E08:-10 E09:-9 E10a:-6 E10b:-6 E11:-11 E12:-1
    hence required component orders
        E01:2 E02:4 E03:3 E04:6 E05:5 E06:8
        E07:7 E08:10 E09:9 E10a/b:6 E11:11 E12:1.

Global degree-two chord:
    s = A/Z^2 + (B/Z) m,
    m=(y+yP)/(x-xP).

With -8F at infinity:
    deg A <= 40, deg B <= 14,
so ambient dimension is 56.
The same 24 smooth P.O collisions impose rank 48, leaving dimension 8.

Resolved local conditions:
For each exceptional component created by blow-up center E_i, pull the section
numerator
    G = A*(Z^2*x-X) + B*(Z^3*y+Y)
through the exact stored map from that infinitely-near chart to the original
I9* germ.  Vanishing of all jets below the component's required order gives
linear conditions on the post-collision 8-space.

At E10 the exceptional divisor splits into E10a/E10b, but both q32 required
orders are 6, so the common center valuation is sufficient for this first
resolved-center test.

If the resolved-center quotient has rank 6 / kernel dimension 2, eliminate
the chord to a binary quartic and classify the Jacobian.
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    GF, PolynomialRing, QQ, ZZ, matrix, sage_eval
)


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
RES=LOCAL/f"q24-i9star-resolution-mod-{args.prime}.json"
GEOCLEAN=ROOT/"elkies-k3/scripts/test_h92_q24_geometrically_clean_d12_neighbor.sage"
q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((
    path for path in q8_candidates
    if path.exists()
    and json.loads(path.read_text()).get("status")=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"
),None)
if Q8 is None:
    raise SystemExit("No passing exact D13 q8 child artifact")

for path in (MOD,RES,Q8,CORE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

mod=json.loads(MOD.read_text())
res=json.loads(RES.read_text())
q8=json.loads(Q8.read_text())
assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert res["status"]=="PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS"
assert res["actual_exceptional_irreducible_components"]==13
assert q8["status"]=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"

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
Acurve=red_poly(child["minimal_A_coefficients_low_to_high"])
Bcurve=red_poly(child["minimal_B_coefficients_low_to_high"])
Delta=-16*(4*Acurve**3+27*Bcurve**2)

sec=mod["section_mod_p"]
Z=R([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
X=R([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
Y=R([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
assert (Z.degree(),X.degree(),Y.degree())==(24,52,78)
assert Z.is_monic()
assert Y**2==X**3+Acurve*X*Z**4+Bcurve*Z**6

xP=K(X)/K(Z**2)
yP=K(Y)/K(Z**3)

# I9* base.
i9=next(item for item in child["finite_fibres"] if item["kodaira"]=="I9*")
RQ=PolynomialRing(QQ,"U")
fQ=RQ(str(i9["factor"]))
f=R([red_q(c) for c in fQ.list()])
assert f.degree()==1
alpha=-f[0]/f[1]

# ---------------------------------------------------------------------------
# 1. q32 global ambient: 56 -> collision rank48 -> dimension8.
# ---------------------------------------------------------------------------

required_inf=8
m_inf=-2
Amax=2*Z.degree()-required_inf
Bmax=Z.degree()+m_inf-required_inf
assert (Amax,Bmax)==(40,14)
ambient=[("A",i) for i in range(Amax+1)]
ambient += [("B",i) for i in range(Bmax+1)]
assert len(ambient)==56

modulus=Z**2
collision_cols=[]
for kind,i in ambient:
    collision_cols.append(
        (U**i*X)%modulus if kind=="A"
        else (-U**i*Y)%modulus
    )
C=matrix(F,48,56,lambda row,col:collision_cols[col][row])
assert C.rank()==48
K8=C.right_kernel().basis_matrix()
assert K8.dimensions()==(8,56)

print(
    "Q32D12RR_GLOBAL|"
    f"required_inf={required_inf}|Amax={Amax}|Bmax={Bmax}|ambient=56|"
    f"collision_rank={C.rank()}|post_collision={K8.nrows()}|status=PASS",
    flush=True,
)

def pair_from_ambient_row(row):
    AA=R.zero()
    BB=R.zero()
    for j,(kind,i) in enumerate(ambient):
        if kind=="A":
            AA += row[j]*U**i
        else:
            BB += row[j]*U**i
    return AA,BB

post_pairs=[pair_from_ambient_row(row) for row in K8.rows()]

# ---------------------------------------------------------------------------
# 2. Reconstruct the original-chart map at each infinitely-near center.
# ---------------------------------------------------------------------------

S=PolynomialRing(F,names=("u","x","y"),order="degrevlex")
u,x,y=S.gens()
localsS={"u":u,"x":x,"y":y}

by_label={rec["label"]:rec for rec in res["centers"]}

def parse_poly(text):
    return S(sage_eval(str(text),locals=localsS))

def center_origin_map(rec):
    if rec["path"]=="root":
        return (u,x,y)
    last=rec["path"].rsplit("/",1)[1]
    parent_label,chart,unused_direction=last.split(":",2)
    parent=by_label[parent_label]
    diag=next(d for d in parent["charts"] if d["chart"]==chart)
    return tuple(parse_poly(v) for v in diag["origin_map"])

contexts={}
for rec in res["centers"]:
    om=center_origin_map(rec)
    pt=tuple(F(v) for v in rec["point"])
    # Every infinitely-near center maps back to the original singular point.
    assert all(expr(*pt)==0 for expr in om)
    contexts[rec["label"]]=(om,pt)

# q32 local correction after placing -8F at infinity.
required_orders={
    "E01":2,
    "E02":4,
    "E03":3,
    "E04":6,
    "E05":5,
    "E06":8,
    "E07":7,
    "E08":10,
    "E09":9,
    "E10":6,   # both E10a/E10b
    "E11":11,
    "E12":1,
}
assert set(required_orders)==set(contexts)

# ---------------------------------------------------------------------------
# 3. Truncated pullback jets of the 8 post-collision basis functions.
#
# Use the polynomial numerator
#   G=A*(Z^2*x-X)+B*(Z^3*y+Y),
# which differs from s by a common unit at the I9* place.
# ---------------------------------------------------------------------------

T=PolynomialRing(F,"t")
t=T.gen()

def shift_univar(poly):
    return T(poly(alpha+t))

def trunc(poly,N):
    poly=S(poly)
    return S({
        exp:coef for exp,coef in poly.dict().items()
        if coef and sum(exp)<N
    })

def mul_trunc(a,b,N):
    return trunc(S(a)*S(b),N)

def pow_trunc(a,n,N):
    out=S.one()
    base=trunc(a,N)
    k=int(n)
    while k:
        if k&1:
            out=mul_trunc(out,base,N)
        k//=2
        if k:
            base=mul_trunc(base,base,N)
    return out

def eval_shifted_univar(poly, local_u, N):
    shifted=shift_univar(R(poly))
    result=S.zero()
    power=S.one()
    for i in range(min(N,shifted.degree()+1 if shifted else 0)):
        if shifted[i]:
            result=trunc(result+S(shifted[i])*power,N)
        power=mul_trunc(power,local_u,N)
    return trunc(result,N)

def translated_origin_map(origin_map, point, N):
    a,b,c=point
    subs=(a+u,b+x,c+y)
    return tuple(trunc(expr(*subs),N) for expr in origin_map)

def local_G(AA,BB,origin_map,point,N):
    ou,ox,oy=translated_origin_map(origin_map,point,N)
    AAl=eval_shifted_univar(AA,ou,N)
    BBl=eval_shifted_univar(BB,ou,N)
    Zl=eval_shifted_univar(Z,ou,N)
    Xl=eval_shifted_univar(X,ou,N)
    Yl=eval_shifted_univar(Y,ou,N)

    z2=mul_trunc(Zl,Zl,N)
    z3=mul_trunc(z2,Zl,N)
    first=mul_trunc(AAl,trunc(mul_trunc(z2,ox,N)-Xl,N),N)
    second=mul_trunc(BBl,trunc(mul_trunc(z3,oy,N)+Yl,N),N)
    return trunc(first+second,N)

# For each center, collect every ambient monomial coefficient below required order
# as a row on the post-collision 8-space.
local_rows=[]
center_ranks={}
for label in sorted(required_orders,key=lambda z:int(z[1:])):
    req=required_orders[label]
    origin_map,point=contexts[label]
    funcs=[
        local_G(AA,BB,origin_map,point,req)
        for AA,BB in post_pairs
    ]

    monomials=sorted({
        exp
        for ffun in funcs
        for exp,coef in ffun.dict().items()
        if coef and sum(exp)<req
    },key=lambda e:(sum(e),e))

    rows=[]
    for exp in monomials:
        row=[ffun.dict().get(exp,F.zero()) for ffun in funcs]
        if any(row):
            rows.append(row)
            local_rows.append(row)

    Mr=matrix(F,rows) if rows else matrix(F,0,8)
    center_ranks[label]=int(Mr.rank())
    print(
        "Q32D12RR_CENTER|"
        f"component={label}|required={req}|jet_rows={len(rows)}|"
        f"local_rank={Mr.rank()}|status=PASS",
        flush=True,
    )

L=matrix(F,local_rows) if local_rows else matrix(F,0,8)
resolved_rank=L.rank()
K2=L.right_kernel().basis_matrix()

print(
    "Q32D12RR_RESOLVED|"
    f"post_collision=8|raw_local_rows={L.nrows()}|"
    f"resolved_rank={resolved_rank}|kernel={K2.nrows()}|"
    f"status={'PASS_H0_TWO' if resolved_rank==6 and K2.nrows()==2 else 'DIAGNOSTIC'}",
    flush=True,
)

# If the geometry-derived quotient is not 8->2, write diagnostics and stop cleanly.
if not (resolved_rank==6 and K2.nrows()==2):
    payload={
        "schema":"elkies-k3.h3-q32-d12-resolved-center-rr-modp.v1",
        "status":"DIAGNOSTIC_Q32_RESOLVED_CENTER_RR",
        "prime":int(p),
        "global":{"ambient":56,"collision_rank":48,"post_collision":8},
        "required_orders":required_orders,
        "center_ranks":center_ranks,
        "resolved":{"raw_rows":int(L.nrows()),"rank":int(resolved_rank),"kernel":int(K2.nrows())},
        "next":(
            "If rank exceeds 6, replace ambient maximal-ideal jets by the exact "
            "surface-local quotient at the offending centers; if below 6, add "
            "the split E10 tangent-line residue conditions."
        ),
    }
    OUT=args.output.resolve() if args.output else LOCAL/f"q32-d12-resolved-center-rr-mod-{p}.json"
    OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(f"OUTPUT|{OUT}",flush=True)
    print(
        "Q32D12RR_RESULT|"
        f"resolved_rank={resolved_rank}|kernel={K2.nrows()}|status={payload['status']}",
        flush=True,
    )
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# 4. Lift the two resolved combinations back to the 56 ambient coefficients.
# ---------------------------------------------------------------------------

ambient_kernel=K2*K8
assert ambient_kernel.dimensions()==(2,56)

pairs=[]
for row in ambient_kernel.rows():
    AA,BB=pair_from_ambient_row(row)
    assert (AA*X-BB*Y)%modulus==0
    pairs.append((AA,BB,K(AA)/K(Z**2),K(BB)/K(Z)))
    print(
        "Q32D12RR_PENCILGEN|"
        f"Adeg={AA.degree()}|Bdeg={BB.degree()}|"
        f"Aterms={sum(bool(c) for c in AA.list())}|"
        f"Bterms={sum(bool(c) for c in BB.list())}|status=PASS",
        flush=True,
    )

# Recheck required center orders on the two final generators.
for gi,(AA,BB,unused_a,unused_b) in enumerate(pairs):
    for label,req in required_orders.items():
        om,pt=contexts[label]
        ffun=local_G(AA,BB,om,pt,req)
        assert not ffun  # truncated below req must vanish identically
    print(
        f"Q32D12RR_GENERATOR_CHECK|i={gi}|all_center_orders=PASS|status=PASS",
        flush=True,
    )

# ---------------------------------------------------------------------------
# 5. Degree-two chord elimination -> binary quartic -> Jacobian.
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
    return UK(lift_poly(R(value.numerator())))/UK(lift_poly(R(value.denominator())))

a0,b0=lift_rf(pairs[0][2]),lift_rf(pairs[0][3])
a1,b1=lift_rf(pairs[1][2]),lift_rf(pairs[1][3])
xPV,yPV=lift_rf(xP),lift_rf(yP)
AV,BV=lift_poly(Acurve),lift_poly(Bcurve)

den=b1-VF(V)*b0
if not den:
    raise ArithmeticError("resolved pencil basis has degenerate chord coefficient")
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
    "Q32D12RR_QUARTIC|"
    f"degree={quartic.degree()}|"
    f"status={'PASS_GENUS_ONE' if quartic.degree() in (3,4) else 'FAIL'}",
    flush=True,
)
if quartic.degree() not in (3,4):
    raise ArithmeticError("resolved q32 pencil did not produce genus one")

I,J=binary_quartic_invariants(quartic)
jacA=VF(-27)*VF(I)
jacB=VF(-27)*VF(J)

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

is_d12=(root_rank,root_det)==(12,4)

print(
    "Q32D12RR_CHILD|"
    f"finite={[(r['degree'],r['minimal_orders'],r['kodaira']) for r in finite]}|"
    f"infinity={inf_orders},{inf_kind}|"
    f"root_rank={root_rank}|euler={euler}|root_det={root_det}|"
    f"status={'PASS_D12' if is_d12 else 'NOT_D12'}",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q32-d12-resolved-center-rr-modp.v1",
    "status":(
        "PASS_Q32_RESOLVED_CENTER_RR_D12"
        if is_d12 else
        "Q32_RESOLVED_CENTER_RR_CHILD_MISMATCH"
    ),
    "prime":int(p),
    "global":{
        "infinity_required":8,
        "A_max_degree":40,
        "B_max_degree":14,
        "ambient":56,
        "collision_rank":48,
        "post_collision":8,
    },
    "resolved":{
        "required_orders":required_orders,
        "center_ranks":center_ranks,
        "raw_rows":int(L.nrows()),
        "rank":int(resolved_rank),
        "kernel_dimension":int(K2.nrows()),
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
    "proof_boundary":(
        "Resolved conditions are derived from the actual infinitely-near blow-up "
        "centers and the q32 geometric vertical correction. At the split E10 "
        "center both spinor required orders are equal, so no tangent-line "
        "asymmetry is needed for membership. A final promotion should still "
        "cross-check the two resulting generators against exact divisorial "
        "valuations on terminal smooth charts."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q32-d12-resolved-center-rr-mod-{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q32D12RR_RESULT|"
    f"resolved_rank={resolved_rank}|kernel={K2.nrows()}|quartic={quartic.degree()}|"
    f"root_rank={root_rank}|root_det={root_det}|status={payload['status']}",
    flush=True,
)
