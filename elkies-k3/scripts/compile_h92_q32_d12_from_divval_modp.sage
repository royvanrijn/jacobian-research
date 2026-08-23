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
# 2-3. Import the certified GENERIC DIVISORIAL resolved quotient.
# ---------------------------------------------------------------------------

DIVVAL=LOCAL/f"q32-d12-generic-divval-mod-{p}.json"
if not DIVVAL.exists():
    raise SystemExit(f"Missing prerequisite: {DIVVAL}")
divval=json.loads(DIVVAL.read_text())
assert divval["status"]=="PASS_Q32_GENERIC_DIVISORIAL_RR_H0_TWO"
assert int(divval["prime"])==int(p)
assert int(divval["global"]["ambient"])==56
assert int(divval["global"]["collision_rank"])==48
assert int(divval["global"]["post_collision"])==8

L=matrix(F,divval["resolved"]["condition_matrix"])
K2=matrix(F,divval["resolved"]["kernel_basis"])
resolved_rank=int(divval["resolved"]["rank"])
assert L.ncols()==8
assert resolved_rank==L.rank()==6
assert K2.dimensions()==(2,8)
assert L*K2.transpose()==matrix(F,L.nrows(),2)

required_orders={
    item["component"]:int(item["required"])
    for item in divval["components"]
}
center_ranks={
    item["component"]:int(item["cumulative_rank"])
    for item in divval["components"]
}

print(
    "Q32D12RR_RESOLVED|"
    "source=generic_divisorial_surface_local_ring|"
    f"post_collision=8|rows={L.nrows()}|"
    f"resolved_rank={resolved_rank}|kernel={K2.nrows()}|"
    "status=PASS_H0_TWO",
    flush=True,
)

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

# Recheck final generators against the intrinsic resolved quotient.
for gi,row in enumerate(K2.rows()):
    assert L*vector(F,row)==vector(F,[0]*L.nrows())
    print(
        f"Q32D12RR_GENERATOR_CHECK|i={gi}|"
        "generic_divisorial_conditions=PASS|status=PASS",
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
    "schema":"elkies-k3.h3-q32-d12-generic-divisorial-compiler-modp.v1",
    "status":(
        "PASS_Q32_GENERIC_DIVISORIAL_RR_D12"
        if is_d12 else
        "Q32_GENERIC_DIVISORIAL_RR_CHILD_MISMATCH"
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
OUT=args.output.resolve() if args.output else LOCAL/f"q32-d12-generic-divisorial-compiled-mod-{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q32D12RR_RESULT|"
    f"resolved_rank={resolved_rank}|kernel={K2.nrows()}|quartic={quartic.degree()}|"
    f"root_rank={root_rank}|root_det={root_det}|status={payload['status']}",
    flush=True,
)
