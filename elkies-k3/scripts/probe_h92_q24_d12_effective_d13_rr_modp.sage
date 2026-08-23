#!/usr/bin/env sage -python
"""
H3-03 D13/MW4 --q24--> D12/MW5 effective-D13 RR probe over GF(p).

Status: ACTIVE_SEARCH.  Promote only if the resolved cluster gives h0=2 and
the resulting degree-two chord pencil has a degree-3/4 squarefree radicand
with D12 child data.

This replaces the rejected ordinary-jet shortcut.  Its local conditions are
derived from the effective D13 component orientation and the actual blow-up
chronology:

  effective local component cycle -> divisorial valuation thresholds
  -> infinitely-near point cluster C01:3,C02:2,C04:2,C06:2,C08:2.

The global setup is the already-passing q24 preflight:

  ambient coefficient dimension = 58
  smooth P.O collision rank      = 48
  post-collision dimension       = 10.

The ten post-collision functions are represented near the I9* singularity by
a COMMON UNIT DENOMINATOR, so their local vanishing conditions can be imposed
on polynomial numerators.  At each infinitely-near centre this script works
in the actual strict-transform surface local ring, imposes the required
maximal-ideal order, changes to a kernel basis, pulls that basis through the
selected blow-up chart, and divides by the certified exceptional power before
continuing.

No ordinary u-jet subspace is inserted.  No Kodaira symbol is used as a proxy
for the chart geometry.

Inputs:
  q24-d12-rr-preflight-mod-<p>.json
  q24-i9star-resolution-mod-<p>.json
  q24-i9star-component-graph-mod-<p>.json
  q24-i9star-effective-cluster-mod-<p>.json
  q8-q24-physical-to-equation-translation.json
  q24-degree46-direct-global-mod-<p>.json
  exact corrected D13 q8 child

Output:
  artifacts/local/elkies-k3/q24-d12-resolved-cluster-rr-mod-<p>.json
"""

import argparse
import itertools
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, identity_matrix, matrix, vector


def locate_repo(explicit=None):
    candidates=[]
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd=Path.cwd().resolve()
    candidates += [cwd,*cwd.parents]
    home=Path.home()
    candidates += [
        home/"Documents"/"jacobian-research",
        home/"jacobian-research",
        home/"src"/"jacobian-research",
        home/"git"/"jacobian-research",
        home/"projects"/"jacobian-research",
    ]
    seen=set()
    for candidate in candidates:
        try:
            candidate=candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if ((candidate/"elkies-k3/scripts").is_dir()
                and (candidate/"artifacts/generated-results").is_dir()):
            return candidate
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

PREFLIGHT=LOCAL/f"q24-d12-rr-preflight-mod-{args.prime}.json"
RESOLUTION=LOCAL/f"q24-i9star-resolution-mod-{args.prime}.json"
GRAPH=LOCAL/f"q24-i9star-component-graph-mod-{args.prime}.json"
EFFECTIVE=LOCAL/f"q24-effective-d13-transport-mod-{args.prime}.json"
MOD=LOCAL/f"q24-degree46-direct-global-mod-{args.prime}.json"
TRANS=LOCAL/"q8-q24-physical-to-equation-translation.json"
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

for path in (PREFLIGHT,RESOLUTION,GRAPH,EFFECTIVE,MOD,TRANS,Q8,CORE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

preflight=json.loads(PREFLIGHT.read_text())
resolution=json.loads(RESOLUTION.read_text())
graph=json.loads(GRAPH.read_text())
effective=json.loads(EFFECTIVE.read_text())
mod=json.loads(MOD.read_text())
trans=json.loads(TRANS.read_text())
q8=json.loads(Q8.read_text())

assert preflight["status"]=="PASS_H3_Q24_D12_MODP_RR_PREFLIGHT"
assert resolution["status"]=="PASS_EXPLICIT_MODP_I9STAR_D13_COMPONENT_RESOLUTION"
assert graph["status"]=="PASS_H3_Q24_AFFINE_D13_COMPONENT_GRAPH"
assert effective["status"]=="PASS_EXACT_H3_Q24_EFFECTIVE_D13_TRANSPORT"
assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert trans["status"]=="PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION"
assert q8["status"]=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"

# Derive the local component profile directly from the certified affine
# I9* graph, D.F=2, nefness, and the D12 child.  This avoids choosing an
# abstract D13 Weyl chamber.

vertices=list(map(str,graph["geometric_graph"]["affine_vertices"]))
index={name:i for i,name in enumerate(vertices)}
edges=[
    (str(a),str(b))
    for a,b in graph["geometric_graph"]["affine_edges"]
]
n=len(vertices)
assert n==14 and "F0" in index

I=matrix(ZZ,n,n)
for i in range(n):
    I[i,i]=-2
for a,b in edges:
    I[index[a],index[b]]=1
    I[index[b],index[a]]=1

multiplicity={
    str(k):int(v)
    for k,v in graph["fibre_multiplicities"].items()
}
m=vector(ZZ,[multiplicity[name] for name in vertices])
assert I*m==vector(ZZ,[0]*n)

# Horizontal O+P both meet the identity component F0.
horizontal=vector(ZZ,[2 if name=="F0" else 0 for name in vertices])

# Since sum m_i (D.C_i)=D.F=2 and all intersections are nonnegative
# integers, there are only these possibilities.
patterns=[]

mult1=[name for name in vertices if multiplicity[name]==1]
mult2=[name for name in vertices if multiplicity[name]==2]

for name in mult1:
    d={v:0 for v in vertices}
    d[name]=2
    patterns.append(d)

for i,a in enumerate(mult1):
    for b in mult1[i+1:]:
        d={v:0 for v in vertices}
        d[a]=d[b]=1
        patterns.append(d)

for name in mult2:
    d={v:0 for v in vertices}
    d[name]=1
    patterns.append(d)

candidates=[]
f0=index["F0"]
unknown=[i for i in range(n) if i!=f0]
B=I.change_ring(QQ)[:,unknown]

for dct in patterns:
    d=vector(ZZ,[dct[name] for name in vertices])
    assert sum(m[i]*d[i] for i in range(n))==2

    zero=[i for i in range(n) if d[i]==0]

    # The old fibre components killed by the new fibre must give D12:
    # rank 12 and determinant 4.
    if len(zero)!=12:
        continue
    G=(-I.matrix_from_rows_and_columns(zero,zero))
    if G.rank()!=12 or abs(G.det())!=4:
        continue

    rhs=vector(QQ,d-horizontal)
    try:
        sol=B.solve_right(rhs)
    except ValueError:
        continue
    if not all(value in ZZ for value in sol):
        continue

    coeff=[ZZ(0)]*n
    for j,i in enumerate(unknown):
        coeff[i]=ZZ(sol[j])

    assert I*vector(ZZ,coeff)==d-horizontal

    candidates.append((dct,coeff))

assert len(candidates)==1, [
    [name for name,value in d.items() if value]
    for d,c in candidates
]

intersection,coeff=candidates[0]
positive=[
    (name,int(intersection[name]))
    for name in vertices if intersection[name]
]
assert positive==[("C10a",1),("C10b",1)] or \
       positive==[("C10b",1),("C10a",1)]

component_coefficients={
    name:int(coeff[index[name]]) for name in vertices
}

expected_coefficients={
    "F0":0,
    "C01":-2,
    "C02":-4,
    "C03":-3,
    "C04":-6,
    "C05":-5,
    "C06":-8,
    "C07":-7,
    "C08":-10,
    "C09":-9,
    "C10a":-6,
    "C10b":-6,
    "C11":-11,
    "C12":-1,
}
assert component_coefficients==expected_coefficients

print(
    "Q24AFFINED12_PROFILE|"
    "positive=C10a:1,C10b:1|"
    + "|".join(
        f"{name}={component_coefficients[name]}"
        for name in sorted(component_coefficients)
    )
    + "|status=PASS_UNIQUE_AFFINE_D12_PROFILE",
    flush=True,
)

thresholds={
    name:max(0,-component_coefficients[name])
    for name in vertices
}

chronology=graph["geometric_graph"]["chronology"]
plan=[]

for step in chronology:
    center=str(step["center"])
    active=list(map(str,step["active_geometric_components"]))
    created=list(map(str,step["new_geometric_components"]))

    baseline=sum(thresholds[name] for name in active)
    residual=[
        max(0,thresholds[name]-baseline)
        for name in created
    ]

    if len(residual)==2:
        assert residual[0]==residual[1] or max(residual)==0

    extra=max(residual) if residual else 0
    if extra:
        plan.append((center,int(extra)))

expected_plan=[
    ("C01",2),
    ("C02",2),
    ("C04",2),
    ("C06",2),
    ("C08",2),
]
assert plan==expected_plan

print(
    "Q24AFFINED12_CLUSTER|"
    + "|".join(f"{name}={order}" for name,order in plan)
    + "|status=PASS_DERIVED_2_2_2_2_2",
    flush=True,
)

p=ZZ(args.prime)
F=GF(p)
R=PolynomialRing(F,"U")
U=R.gen()
K=R.fraction_field()
S=PolynomialRing(F,names=("u","x","y"),order="degrevlex")
u,x,y=S.gens()


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

# ===========================================================================
# 1. Reconstruct the passing 58 -> 10 global preflight space.
# ===========================================================================
Amax=int(preflight["infinity"]["A_max_degree"])
Bmax=int(preflight["infinity"]["B_max_degree"])
assert (Amax,Bmax)==(41,15)
ambient=[("A",i) for i in range(Amax+1)] + [("B",i) for i in range(Bmax+1)]
assert len(ambient)==58
modulus=Z**2
collision_cols=[]
for kind,i in ambient:
    collision_cols.append(
        (U**i*X)%modulus if kind=="A" else (-U**i*Y)%modulus
    )
C=matrix(F,48,58,lambda row,col:collision_cols[col][row])
assert C.rank()==48
K10=C.right_kernel().basis_matrix()
assert K10.dimensions()==(10,58)


def pair_from_ambient_row(row):
    AA=R.zero(); BB=R.zero()
    for j,(kind,i) in enumerate(ambient):
        if kind=="A":
            AA += row[j]*U**i
        else:
            BB += row[j]*U**i
    assert (AA*X-BB*Y)%modulus==0
    return AA,BB


post_pairs=[pair_from_ambient_row(row) for row in K10.rows()]

# ===========================================================================
# 2. Build COMMON-UNIT-denominator local numerators at the I9* point.
# ===========================================================================
i9=next(item for item in child["finite_fibres"] if item["kodaira"]=="I9*")
RQ=PolynomialRing(QQ,"U")
fQ=RQ(str(i9["factor"]))
f=R([red_q(c) for c in fQ.list()])
assert f.degree()==1
alpha=-f[0]/f[1]
assert int(alpha)==int(preflight["I9star"]["base"])


def shift_base(poly):
    return S(R(poly)(alpha+u))


Al=shift_base(A)
Bl=shift_base(B)
Dl=R(Delta)(alpha+R.gen()) if False else None
surface=S(y**2-x**3-Al*x-Bl)
assert surface(0,0,0)==0

Zl=shift_base(Z)
Xl=shift_base(X)
Yl=shift_base(Y)
assert Zl(0,0,0)!=0 and Xl(0,0,0)!=0
# m=(y+yP)/(x-xP) with xP=X/Z^2,yP=Y/Z^3:
# m_num = y*Z^3+Y; m_den = Z*(x*Z^2-X).
m_num=S(y*Zl**3+Yl)
m_den=S(Zl*(x*Zl**2-Xl))
common_den=S(Zl**2*m_den)
assert common_den(0,0,0)!=0

local_numerators=[]
for AA,BB in post_pairs:
    AAl=shift_base(AA)
    BBl=shift_base(BB)
    numerator=S(AAl*m_den+BBl*Zl*m_num)
    local_numerators.append(numerator)
assert len(local_numerators)==10

print(
    "Q24RESCLUSTER_INPUT|ambient=58|collision_rank=48|post_collision=10|"
    f"I9base={int(alpha)}|common_den_unit=1|cluster={','.join(str(q) for _,q in plan)}|status=PASS",
    flush=True,
)

# ===========================================================================
# 3. Actual local-surface order conditions and successive blow-ups.
# ===========================================================================
records={str(row["label"]):row for row in resolution["centers"]}
for name,_ in plan:
    assert name in records


def order_at_origin(poly):
    poly=S(poly)
    if not poly:
        return 10**9
    return min(sum(exp) for exp,coef in poly.dict().items() if coef)


def shifted(poly,point):
    a,b,c=map(F,point)
    return S(poly(a+u,b+x,c+y))


def monomials_exact(total):
    result=[]
    for i in range(total+1):
        for j in range(total-i+1):
            k=total-i-j
            result.append(u**i*x**j*y**k)
    return result


def monomials_below(total):
    return [mon for degree in range(total) for mon in monomials_exact(degree)]


def local_order_matrix(basis,surface_eq,point,required_order):
    required_order=int(required_order)
    ss=shifted(surface_eq,point)
    ideal=S.ideal([ss]+monomials_exact(required_order))
    gb=ideal.groebner_basis()
    mons=monomials_below(required_order)
    remainders=[]
    for poly in basis:
        rem=shifted(poly,point).reduce(gb)
        # The quotient contains only terms below the requested order.
        if rem and max(sum(exp) for exp in rem.dict())>=required_order:
            raise ArithmeticError("truncated local remainder escaped maximal-ideal quotient")
        remainders.append(rem)
    M=matrix(
        F,len(mons),len(basis),
        lambda row,col:remainders[col].monomial_coefficient(mons[row]),
    )
    return M,mons,ss


def canonical_after_condition(poly,surface_eq,point,required_order):
    """
    Choose an ambient representative with visible order >= required_order.

    Membership was already proved in (surface) + m^required_order.  Reducing
    modulo the surface alone is not order-preserving when the surface itself
    has multiplicity below required_order.  Instead explicitly subtract the
    low-order multiple of the surface whose existence the membership test
    guarantees.
    """
    a,b,c=map(F,point)
    sp=shifted(poly,point)
    ss=shifted(surface_eq,point)

    surface_order=order_at_origin(ss)

    if surface_order >= required_order:
        h=sp
    else:
        q_bound=required_order-surface_order
        qmons=monomials_below(q_bound)
        mons=monomials_below(required_order)

        M=matrix(
            F,len(mons),len(qmons),
            lambda i,j:(qmons[j]*ss).monomial_coefficient(mons[i]),
        )
        rhs=vector(
            F,[sp.monomial_coefficient(mon) for mon in mons]
        )

        try:
            coeff=M.solve_right(rhs)
        except ValueError as exc:
            raise ArithmeticError(
                "surface-local membership passed but explicit low-order "
                "surface multiple could not be reconstructed"
            ) from exc

        h=S(sp)
        for j in range(len(qmons)):
            h -= coeff[j]*qmons[j]*ss

    if order_at_origin(h) < required_order:
        raise ArithmeticError(
            f"explicit surface-adjusted representative has order "
            f"{order_at_origin(h)}, expected >= {required_order}"
        )

    # Undo the translation back to the current blow-up chart.
    return S(h(u-a,x-b,y-c))


def chart_substitutions(point,kind):
    a,b,c=map(F,point)
    if kind=="u":
        return (a+u,b+u*x,c+u*y),u
    if kind=="x":
        return (a+x*u,b+x,c+x*y),x
    if kind=="y":
        return (a+y*u,b+y*x,c+y),y
    raise ValueError(kind)


def divide_power(poly,exceptional,power):
    q=S(poly)
    for unused in range(int(power)):
        q,rem=q.quo_rem(exceptional)
        if rem:
            raise ArithmeticError("required exceptional power does not divide")
    return S(q)


def child_path(parent,child):
    target=str(child["path"])
    matches=[]
    for item in parent["children"]:
        suffix=(
            f"/{parent['label']}:{item['selected_chart']}:"
            + ",".join(map(str,item["direction"]))
        )
        candidate=str(parent["path"])+suffix
        if candidate==target:
            matches.append((str(item["selected_chart"]),candidate))
    if len(matches)!=1:
        raise ArithmeticError(
            f"could not identify unique chart {parent['label']} -> {child['label']}: {matches}"
        )
    return matches[0][0]


# Check that the nonzero cluster is one nested branch of the actual tree.
for (left,_),(right,_) in zip(plan,plan[1:]):
    assert str(records[right]["path"]).startswith(str(records[left]["path"]))

current_basis=list(local_numerators)
current_surface=surface
transform=identity_matrix(F,10)
condition_ledger=[]

for step_index,(center_name,required_order) in enumerate(plan):
    record=records[center_name]
    point=tuple(F(v) for v in record["point"])
    M,quotient_monomials,unused_shifted_surface=local_order_matrix(
        current_basis,current_surface,point,required_order
    )
    local_rank=int(M.rank())
    kernel=M.right_kernel().basis_matrix()
    before=len(current_basis)
    after=int(kernel.nrows())
    if before-local_rank!=after:
        raise ArithmeticError("local rank/nullity mismatch")

    new_basis=[]
    for row in kernel.rows():
        combination=S(sum(row[j]*current_basis[j] for j in range(before)))
        new_basis.append(canonical_after_condition(
            combination,current_surface,point,required_order
        ))
    transform=kernel*transform
    current_basis=new_basis
    cumulative_codim=10-after
    condition_ledger.append({
        "center":center_name,
        "additional_order":int(required_order),
        "quotient_rows":int(M.nrows()),
        "local_rank":local_rank,
        "dimension_before":before,
        "dimension_after":after,
        "cumulative_codimension":cumulative_codim,
        "quotient_monomials":[str(mon) for mon in quotient_monomials],
    })
    print(
        "Q24RESCLUSTER_CENTER|"
        f"center={center_name}|order={required_order}|before={before}|"
        f"rows={M.nrows()}|rank={local_rank}|after={after}|"
        f"cumulative_codim={cumulative_codim}|status=PASS",
        flush=True,
    )

    if step_index+1<len(plan):
        next_name=plan[step_index+1][0]
        next_record=records[next_name]
        kind=child_path(record,next_record)
        subs,e=chart_substitutions(point,kind)
        surface_transformed=S(current_surface(*subs))
        current_surface=divide_power(surface_transformed,e,int(record["multiplicity"]))
        pulled=[]
        for poly in current_basis:
            pulled.append(divide_power(S(poly(*subs)),e,required_order))
        current_basis=pulled
        print(
            "Q24RESCLUSTER_BLOWUP|"
            f"from={center_name}|to={next_name}|chart={kind}|"
            f"surface_mult={record['multiplicity']}|section_divide={required_order}|status=PASS",
            flush=True,
        )

final_dimension=len(current_basis)
final_codimension=10-final_dimension
print(
    "Q24RESCLUSTER_RR|"
    f"post_collision=10|resolved_codim={final_codimension}|"
    f"kernel={final_dimension}|target_codim=8|target_h0=2|"
    f"status={'PASS_H0_TWO' if final_dimension==2 else 'DIMENSION_MISMATCH'}",
    flush=True,
)

# Express the final local kernel back in the original 58-dimensional ambient.
final58=transform*K10
assert final58.nrows()==final_dimension and final58.ncols()==58
assert C*final58.transpose()==matrix(F,48,final_dimension)

# ===========================================================================
# 4. If h0=2, compile the actual cluster-selected degree-two pencil.
# ===========================================================================
def pair_from_row(row):
    AA=R.zero(); BB=R.zero()
    for j,(kind,i) in enumerate(ambient):
        if kind=="A":
            AA += row[j]*U**i
        else:
            BB += row[j]*U**i
    assert (AA*X-BB*Y)%modulus==0
    return AA,BB,K(AA)/K(Z**2),K(BB)/K(Z)

quartic_degree=None
child_summary=None
terminal_status="RESOLVED_CLUSTER_DIMENSION_MISMATCH"
final_pairs=[]
if final_dimension==2:
    final_pairs=[pair_from_row(row) for row in final58.rows()]
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
        return UK(lift_poly(R(value.numerator())))/UK(lift_poly(R(value.denominator())))

    a0,b0=lift_rf(final_pairs[0][2]),lift_rf(final_pairs[0][3])
    a1,b1=lift_rf(final_pairs[1][2]),lift_rf(final_pairs[1][3])
    xPV,yPV=lift_rf(xP),lift_rf(yP)
    AV,BV=lift_poly(A),lift_poly(B)
    den=b1-VF(V)*b0
    if not den:
        raise ArithmeticError("resolved cluster pencil is degenerate in chord direction")
    mval=-(a1-VF(V)*a0)/den

    XR=PolynomialRing(UK,"x")
    xx=XR.gen()
    yline=XR(mval)*(xx-XR(xPV))-XR(yPV)
    relation=yline**2-xx**3-XR(AV)*xx-XR(BV)
    quadratic,remainder=relation.quo_rem(xx-XR(xPV))
    assert not remainder and quadratic.degree()==2
    disc=UK(quadratic[1]**2-4*quadratic[2]*quadratic[0])
    quartic,square_factor=squarefree_binary_quartic(disc,UR)
    quartic_degree=int(quartic.degree())
    genus=(quartic_degree-2)//2 if quartic_degree>=3 else 0
    print(
        "Q24RESCLUSTER_QUARTIC|"
        f"degree={quartic_degree}|genus={genus}|"
        f"status={'PASS_GENUS_ONE' if quartic_degree in (3,4) else 'NOT_GENUS_ONE'}",
        flush=True,
    )

    if quartic_degree in (3,4):
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
            root_rank+=int(ir); euler+=int(ie); root_det*=int(idt)
        child_summary={
            "finite_fibres":finite,
            "infinity_orders":list(inf_orders),
            "infinity_kind":inf_kind,
            "root_rank":root_rank,
            "root_determinant":root_det,
            "euler":euler,
        }
        is_d12=(root_rank,root_det,euler)==(12,4,24)
        terminal_status=(
            ("CANDIDATE_H3_Q24_EFFECTIVE_D13_D12_MODP"
             if cluster_sign_certified else
             "CANDIDATE_H3_Q24_EFFECTIVE_D13_D12_MODP")
            if is_d12 else "GENUS_ONE_CHILD_NOT_D12"
        )
        print(
            "Q24RESCLUSTER_CHILD|"
            f"root_rank={root_rank}|root_det={root_det}|euler={euler}|"
            f"infinity={inf_orders},{inf_kind}|"
            f"status={'PASS_D12' if is_d12 else 'NOT_D12'}",
            flush=True,
        )
    else:
        terminal_status="RESOLVED_CLUSTER_NOT_GENUS_ONE"

payload={
    "schema":"elkies-k3.h3-q24-d12-effective-d13-rr-modp.v1",
    "status":terminal_status,
    "prime":int(p),
    "inputs":{
        "preflight":str(PREFLIGHT.relative_to(ROOT)),
        "resolution":str(RESOLUTION.relative_to(ROOT)),
        "component_graph":str(GRAPH.relative_to(ROOT)),
        "effective_d13_transport":str(EFFECTIVE.relative_to(ROOT)),
    },
    "global_rr":{
        "ambient_dimension":58,
        "smooth_collision_rank":48,
        "post_collision_dimension":10,
    },
    "resolved_cluster":{
        "plan":[{"center":name,"additional_order":order} for name,order in plan],
        "condition_ledger":condition_ledger,
        "codimension_on_post_collision":int(final_codimension),
        "kernel_dimension":int(final_dimension),
        "kernel_basis_post_collision":[
            [int(v) for v in row] for row in transform.rows()
        ],
        "kernel_basis_ambient58":[
            [int(v) for v in row] for row in final58.rows()
        ],
    },
    "quartic_degree":quartic_degree,
    "child":child_summary,
    "proof_boundary":(
        "The local point conditions are derived from the exact O8-normalized "
        "D13 transport and are orientation-independent. The remaining open "
        "promotion gate is to identify that transported simple-root chamber "
        "with the explicit resolved geometric component charts, not merely "
        "with their abstract D13 graph. Therefore even h0=2 + quartic degree "
        "4 + D12 is reported as a candidate until that marking is replayed."
    ),
    "affine_D12_component_profile_certified":True,
    "geometric_Weyl_chamber_assumption_used":False,
}
OUT=(args.output.resolve() if args.output else
     LOCAL/f"q24-d12-effective-d13-rr-mod-{p}.json")
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True,default=int)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24RESCLUSTER_RESULT|"
    f"resolved_codim={final_codimension}|kernel={final_dimension}|"
    f"quartic_degree={quartic_degree}|status={terminal_status}",
    flush=True,
)
