#!/usr/bin/env sage -python
"""
Scan monomial ordinary-jet codimension-eight quotients for the H3 q24 step.

Status: ACTIVE_SEARCH / diagnostic only.

The exact q24 divisor and the modular RR preflight give a 10-dimensional
post-collision space.  The failed first-eight-jet experiment retained jet
orders 8 and 9 and produced a squarefree degree-18 hyperelliptic radicand,
so ordinary u^8 vanishing is not the resolved I9* quotient.

This script uses the fact that the first ten ordinary cusp jets identify the
10-dimensional post-collision space.  For every pair of allowed jet orders
0 <= i < j <= 9, it kills the other eight orders, obtains a two-dimensional
pencil, and measures the squarefree chord-radicand degree.  Degree 3 or 4 is
a genus-one candidate and is then classified by binary-quartic invariants.

A successful pair is only a pattern hint for the true resolved module.  A
proof still requires the actual split I9* component charts, the D13 graph
matching, and the connected resolved quotient.
"""

import argparse
import json
from itertools import combinations
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
parser.add_argument(
    "--only-pair",
    nargs=2,
    type=int,
    metavar=("I","J"),
    help="scan only one allowed jet pair",
)
parser.add_argument(
    "--stop-on-genus-one",
    action="store_true",
    help="stop after the first degree-3/4 candidate",
)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"
CORE=ROOT/"elkies-k3/scripts/elliptic_neighbor_compiler.sage"
exec(compile(CORE.read_text(),str(CORE),"exec"))

MOD=LOCAL/f"q24-degree46-direct-global-mod-{args.prime}.json"
TRANS=LOCAL/"q8-q24-physical-to-equation-translation.json"
PREFLIGHT=LOCAL/f"q24-d12-rr-preflight-mod-{args.prime}.json"
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
for path in (MOD,TRANS,PREFLIGHT,Q8,CORE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

mod=json.loads(MOD.read_text())
trans=json.loads(TRANS.read_text())
preflight=json.loads(PREFLIGHT.read_text())
q8=json.loads(Q8.read_text())
assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert trans["status"]=="PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION"
assert preflight["status"]=="PASS_H3_Q24_D12_MODP_RR_PREFLIGHT"

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

# 58-dimensional global ambient and rank-48 smooth-collision condition.
ambient=[("A",i) for i in range(42)]
ambient += [("B",i) for i in range(16)]
assert len(ambient)==58
modulus=Z**2
collision_cols=[
    (U**i*X)%modulus if kind=="A" else (-U**i*Y)%modulus
    for kind,i in ambient
]
C=matrix(F,48,58,lambda row,col:collision_cols[col][row])
assert C.rank()==48
K10=C.right_kernel().basis_matrix()
assert K10.dimensions()==(10,58)

# I9* cusp and the full ten-jet evaluation.
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


def jet_of_rf(value,n):
    value=KS(value)
    num=S(value.numerator())
    den=S(value.denominator())
    if not den[0]:
        raise ZeroDivisionError("jet denominator not a unit")
    rem=(num*den.inverse_mod(u**n))%(u**n)
    return [rem[i] for i in range(n)]


Al=shift_poly(A)
Bl=shift_poly(B)
Dl=shift_poly(Delta)
assert (Al.valuation(),Bl.valuation(),Dl.valuation())==(2,3,15)

xloc=shift_rf(xP)
yloc=shift_rf(yP)
x0=F(xloc(0)); y0=F(yloc(0))
assert x0 or y0
assert y0**2==x0**3
mcusp=-yloc/xloc
assert mcusp.valuation()>=0

Jcols=[]
for kind,i in ambient:
    if kind=="A":
        local=shift_rf(K(U**i)/K(Z**2))
    else:
        local=shift_rf(K(U**i)/K(Z))*mcusp
    Jcols.append(jet_of_rf(local,10))
J10=matrix(F,10,58,lambda row,col:Jcols[col][row])
restricted=J10*K10.transpose()
assert restricted.rank()==10
assert C.stack(J10).rank()==58

# Shared new-base rings for elimination.
VR=PolynomialRing(F,"V")
V=VR.gen()
VF=VR.fraction_field()
UR=PolynomialRing(VF,"U")
UK=UR.fraction_field()


def lift_poly(poly):
    poly=R(poly)
    return UR([VF(c) for c in poly.list()])


def lift_rf(value):
    value=K(value)
    return UK(lift_poly(R(value.numerator())))/UK(
        lift_poly(R(value.denominator()))
    )


xPV,yPV=lift_rf(xP),lift_rf(yP)
AV,BV=lift_poly(A),lift_poly(B)


def pair_from_row(row):
    AA=R.zero()
    BB=R.zero()
    for column,(kind,power) in enumerate(ambient):
        if kind=="A":
            AA += row[column]*U**power
        else:
            BB += row[column]*U**power
    assert (AA*X-BB*Y)%modulus==0
    return (
        AA,
        BB,
        K(AA)/K(Z**2),
        K(BB)/K(Z),
    )


def classify_pair(allowed):
    disallowed=[index for index in range(10) if index not in allowed]
    local_rows=J10.matrix_from_rows(disallowed)
    M=C.stack(local_rows)
    combined_rank=int(M.rank())
    kernel=M.right_kernel().basis_matrix()
    result={
        "allowed_jet_orders":list(map(int,allowed)),
        "combined_rank":combined_rank,
        "kernel_dimension":int(kernel.nrows()),
    }
    if combined_rank!=56 or kernel.nrows()!=2:
        result["status"]="UNEXPECTED_RANK"
        return result

    pairs=[pair_from_row(row) for row in kernel.rows()]
    (unused_A0,unused_B0,a0raw,b0raw),(unused_A1,unused_B1,a1raw,b1raw)=pairs
    a0,b0=lift_rf(a0raw),lift_rf(b0raw)
    a1,b1=lift_rf(a1raw),lift_rf(b1raw)

    den=b1-VF(V)*b0
    if not den:
        result["status"]="DEGENERATE_CHORD_COEFFICIENT"
        return result
    mval=-(a1-VF(V)*a0)/den

    XR=PolynomialRing(UK,"x")
    xx=XR.gen()
    yline=XR(mval)*(xx-XR(xPV))-XR(yPV)
    relation=yline**2-xx**3-XR(AV)*xx-XR(BV)
    quadratic,remainder=relation.quo_rem(xx-XR(xPV))
    if remainder or quadratic.degree()!=2:
        result["status"]="ELIMINATION_FAILED"
        return result

    disc=UK(quadratic[1]**2-4*quadratic[2]*quadratic[0])
    try:
        squarefree,square_factor=squarefree_binary_quartic(disc,UR)
    except Exception as exc:
        result["status"]="SQUAREFREE_EXTRACTION_FAILED"
        result["error"]=type(exc).__name__+":"+str(exc)
        return result

    degree=int(squarefree.degree())
    genus=max(0,(degree-1)//2)
    factor_degrees=[
        {"degree":int(factor.degree()),"multiplicity":int(exponent)}
        for factor,exponent in squarefree.factor()
    ]
    result.update({
        "squarefree_degree":degree,
        "hyperelliptic_genus_if_smooth":int(genus),
        "squarefree_factor_degrees":factor_degrees,
    })
    if degree not in (3,4):
        result["status"]="HIGHER_GENUS"
        return result

    I,J=binary_quartic_invariants(squarefree)
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
    result["child"]={
        "finite_fibres":finite,
        "infinity_orders":list(inf_orders),
        "infinity_kind":inf_kind,
        "root_rank":root_rank,
        "root_determinant":root_det,
        "euler":euler,
    }
    result["status"]=(
        "GENUS_ONE_D12_CANDIDATE"
        if (root_rank,root_det,euler)==(12,4,24)
        else "GENUS_ONE_NON_D12"
    )
    return result


if args.only_pair is not None:
    i,j=sorted(args.only_pair)
    if i<0 or j>9 or i==j:
        raise ValueError("--only-pair must give two distinct orders in 0..9")
    allowed_pairs=[(i,j)]
else:
    allowed_pairs=list(combinations(range(10),2))

results=[]
for index,allowed in enumerate(allowed_pairs,1):
    result=classify_pair(allowed)
    results.append(result)
    print(
        "Q24JETPAIR|"
        f"pair={allowed[0]},{allowed[1]}|"
        f"rank={result.get('combined_rank')}|"
        f"kernel={result.get('kernel_dimension')}|"
        f"degree={result.get('squarefree_degree','NA')}|"
        f"genus={result.get('hyperelliptic_genus_if_smooth','NA')}|"
        f"status={result['status']}",
        flush=True,
    )
    if args.stop_on_genus_one and result["status"].startswith("GENUS_ONE"):
        break

candidates=[
    row for row in results
    if row["status"].startswith("GENUS_ONE")
]
d12_candidates=[
    row for row in results
    if row["status"]=="GENUS_ONE_D12_CANDIDATE"
]
status=(
    "PASS_DIAGNOSTIC_D12_JET_PAIR_CANDIDATES"
    if d12_candidates else
    "PASS_DIAGNOSTIC_GENUS_ONE_JET_PAIR_CANDIDATES"
    if candidates else
    "PASS_NO_MONOMIAL_JET_PAIR_GENUS_ONE"
)
payload={
    "schema":"elkies-k3.h3-q24-d12-ordinary-jet-pair-scan-modp.v1",
    "status":status,
    "prime":int(p),
    "post_collision_dimension":10,
    "ordinary_jet_map_rank":10,
    "pairs_scanned":len(results),
    "results":results,
    "proof_boundary":(
        "This scans monomial coordinate subspaces in the ordinary cusp-jet "
        "quotient.  Even a D12 candidate is not a resolved I9* certificate. "
        "The actual exceptional-component valuations and connected chart "
        "quotient must independently recover the same two-dimensional space."
    ),
    "next":(
        "Use the tangent-cone-aware I9* resolution to split geometric "
        "exceptional branches, match the D13 graph to the lattice vertical "
        "cycle, and derive the actual codimension-eight resolved quotient."
    ),
}
OUT=(
    args.output.resolve()
    if args.output else
    LOCAL/f"q24-d12-ordinary-jet-pair-scan-mod-{p}.json"
)
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24JETPAIR_RESULT|"
    f"scanned={len(results)}|genus_one={len(candidates)}|"
    f"d12={len(d12_candidates)}|status={status}",
    flush=True,
)
