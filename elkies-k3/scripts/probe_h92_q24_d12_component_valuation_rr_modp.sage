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
# 1. Global degree-two ambient in the GEOMETRIC fibre trivialization.
# ===========================================================================
geometric_vertical_square=ZZ(-12)
horizontal_square=ZZ(44)
geometric_fibre_twist=ZZ(-8)
assert horizontal_square + geometric_vertical_square + 4*geometric_fibre_twist == 0

required=int(-geometric_fibre_twist)
m_inf=int(preflight["infinity"]["marked_chord_order"])
assert required==8 and m_inf==-2

Amax=int(2*Z.degree()-required)
Bmax=int(Z.degree()+m_inf-required)
assert (Amax,Bmax)==(40,14)

ambient=[("A",i) for i in range(Amax+1)] + [("B",i) for i in range(Bmax+1)]
assert len(ambient)==56

print(
    "Q24GEOM_INFINITY|"
    f"deterministic_fibre=-7|geometric_fibre={geometric_fibre_twist}|"
    f"vertical_square={geometric_vertical_square}|D_square=0|"
    f"required={required}|Amax={Amax}|Bmax={Bmax}|"
    f"ambient={len(ambient)}|status=PASS_ISOTROPIC_TRIVIALIZATION",
    flush=True,
)

modulus=Z**2
collision_cols=[]
for kind,i in ambient:
    collision_cols.append(
        (U**i*X)%modulus if kind=="A" else (-U**i*Y)%modulus
    )

C=matrix(
    F,
    modulus.degree(),
    len(collision_cols),
    lambda row,col:collision_cols[col][row],
)
collision_rank=int(C.rank())
post_dim=int(C.ncols()-collision_rank)

print(
    "Q24GEOM_COLLISION|"
    f"rows={C.nrows()}|cols={C.ncols()}|rank={collision_rank}|"
    f"post_collision={post_dim}|"
    f"status={'PASS' if collision_rank==48 and post_dim==8 else 'UNEXPECTED'}",
    flush=True,
)
assert collision_rank==48
assert post_dim==8

K10=C.right_kernel().basis_matrix()
assert K10.dimensions()==(post_dim,len(ambient))

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
assert len(local_numerators)==post_dim

print(
    "Q24RESCLUSTER_INPUT|" f"ambient={len(ambient)}|collision_rank={collision_rank}|post_collision={post_dim}|"
    f"I9base={int(alpha)}|common_den_unit=1|cluster={','.join(str(q) for _,q in plan)}|status=PASS",
    flush=True,
)

# ===========================================================================
# 3. Direct divisorial valuation conditions on the resolved exceptional cover.
# ===========================================================================
#
# The actual condition for a component E with local equation e=0 is
#
#       f in (strict_surface, e^t)
#
# in charts covering E.  No maximal-ideal jet proxy is used.

records={str(row["label"]):row for row in resolution["centers"]}
assert set(records)=={f"C{i:02d}" for i in range(1,13)}

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
            raise ArithmeticError(
                "strict-transform exceptional power does not divide"
            )
    return S(q)

path_to_label={str(rec["path"]):label for label,rec in records.items()}
parent_link={}
for parent_label,parent in records.items():
    for child in parent["children"]:
        kind=str(child["selected_chart"])
        direction=",".join(map(str,child["direction"]))
        child_path=(
            str(parent["path"])
            + f"/{parent_label}:{kind}:{direction}"
        )
        child_label=path_to_label.get(child_path)
        if child_label is None:
            raise ArithmeticError(
                f"resolution child path not found: {child_path}"
            )
        if child_label in parent_link:
            raise ArithmeticError(
                f"resolution child has two parents: {child_label}"
            )
        parent_link[child_label]=(parent_label,kind)

roots=[label for label in records if label not in parent_link]
assert roots==["C01"]
assert len(parent_link)==11

state_cache={"C01":(surface,list(local_numerators))}

def state_before(label):
    if label in state_cache:
        return state_cache[label]

    parent_label,kind=parent_link[label]
    parent_surface,parent_basis=state_before(parent_label)
    parent=records[parent_label]
    point=tuple(F(v) for v in parent["point"])

    subs,e=chart_substitutions(point,kind)
    transformed_surface=S(parent_surface(*subs))
    child_surface=divide_power(
        transformed_surface,e,int(parent["multiplicity"])
    )

    child_basis=[S(poly(*subs)) for poly in parent_basis]
    state_cache[label]=(child_surface,child_basis)
    return state_cache[label]

def quotient_matrix_for_chart(
    surface_before,basis_before,record,kind,threshold
):
    point=tuple(F(v) for v in record["point"])
    subs,e=chart_substitutions(point,kind)

    strict=divide_power(
        S(surface_before(*subs)),
        e,
        int(record["multiplicity"]),
    )
    pulled=[S(poly(*subs)) for poly in basis_before]

    exceptional_restriction=S(strict.subs({e:0}))

    # This affine chart can miss the exceptional curve completely.
    if exceptional_restriction in F and exceptional_restriction != 0:
        return None

    # Work branch-by-branch on the REDUCED exceptional divisor.
    factors=[
        (S(factor),int(multiplicity))
        for factor,multiplicity in exceptional_restriction.factor()
    ]
    if not factors:
        raise ArithmeticError(
            f"empty exceptional factorization at {record['label']}:{kind}"
        )

    branch_matrices=[]
    branch_records=[]

    for h,scheme_multiplicity in factors:
        # Prime ideal of the reduced exceptional branch in this chart.
        P=S.ideal([e,h])

        # First form the ordinary power condition in the surface ring.
        J=S.ideal([strict]) + P**int(threshold)

        # Divisorial valuation lives at the GENERIC point of this component.
        # Remove embedded/closed-point conditions by retaining only the
        # primary component supported on P.
        pieces=J.primary_decomposition()

        primary=[]
        for Q in pieces:
            radical=Q.radical()
            # Compare ideals by mutual containment; this is more robust than
            # relying on object equality/canonical Groebner presentations.
            if radical <= P and P <= radical:
                primary.append(Q)

        if len(primary)!=1:
            raise ArithmeticError(
                f"expected unique P-primary component at "
                f"{record['label']}:{kind} branch={h}; "
                f"found {len(primary)} among {len(pieces)}"
            )

        Jsym=primary[0]
        gb=Jsym.groebner_basis()

        remainders=[poly.reduce(gb) for poly in pulled]

        exponents=sorted({
            exp
            for rem in remainders
            for exp,coef in rem.dict().items()
            if coef
        })

        M=matrix(
            F,
            len(exponents),
            len(basis_before),
            lambda i,j:remainders[j].dict().get(
                exponents[i],F(0)
            ),
        )

        branch_matrices.append(M)
        branch_records.append({
            "factor":str(h),
            "scheme_multiplicity":int(scheme_multiplicity),
            "primary_components":int(len(pieces)),
            "rows":int(M.nrows()),
            "rank":int(M.rank()),
        })

    Mtotal=matrix(F,0,len(basis_before))
    for M in branch_matrices:
        if M.nrows():
            Mtotal=Mtotal.stack(M)

    return {
        "chart":kind,
        "matrix":Mtotal,
        "exceptional_restriction":str(exceptional_restriction),
        "branches":branch_records,
    }


centre_thresholds={}
for label in records:
    if label=="C10":
        left=int(thresholds["C10a"])
        right=int(thresholds["C10b"])
        assert left==right
        centre_thresholds[label]=left
    else:
        centre_thresholds[label]=int(thresholds[label])

assert centre_thresholds=={
    "C01":2,
    "C02":4,
    "C03":3,
    "C04":6,
    "C05":5,
    "C06":8,
    "C07":7,
    "C08":10,
    "C09":9,
    "C10":6,
    "C11":11,
    "C12":1,
}

all_rows=matrix(F,0,post_dim)
component_ledger=[]

for label in sorted(records):
    threshold=centre_thresholds[label]
    surface_before,basis_before=state_before(label)
    assert len(basis_before)==post_dim

    chart_results=[]
    centre_rows=matrix(F,0,post_dim)

    for kind in ("u","x","y"):
        result=quotient_matrix_for_chart(
            surface_before,basis_before,records[label],kind,threshold
        )
        if result is None:
            continue
        chart_results.append(result)
        M=result["matrix"]
        if M.nrows():
            centre_rows=centre_rows.stack(M)

    if not chart_results:
        raise ArithmeticError(
            f"no blow-up chart meets exceptional divisor {label}"
        )

    centre_rank=int(centre_rows.rank())
    before_rank=int(all_rows.rank())
    if centre_rows.nrows():
        all_rows=all_rows.stack(centre_rows)
    cumulative_rank=int(all_rows.rank())
    new_rank=cumulative_rank-before_rank

    component_ledger.append({
        "component":label,
        "threshold":int(threshold),
        "charts":[
            {
                "chart":entry["chart"],
                "rows":int(entry["matrix"].nrows()),
                "rank":int(entry["matrix"].rank()),
                "exceptional_restriction":entry[
                    "exceptional_restriction"
                ],
            }
            for entry in chart_results
        ],
        "component_cover_rank":centre_rank,
        "new_global_rank":new_rank,
        "cumulative_rank":cumulative_rank,
    })

    print(
        "Q24DIVVAL_COMPONENT|"
        f"component={label}|threshold={threshold}|"
        f"charts={len(chart_results)}|"
        f"cover_rank={centre_rank}|new_rank={new_rank}|"
        f"cumulative_rank={cumulative_rank}|status=PASS",
        flush=True,
    )

resolved_rank=int(all_rows.rank())
kernel10=all_rows.right_kernel().basis_matrix()
final_dimension=int(kernel10.nrows())
final_codimension=resolved_rank

assert final_codimension==post_dim-final_dimension

print(
    "Q24DIVVAL_RR|"
    f"post_collision={post_dim}|resolved_rank={resolved_rank}|"
    f"kernel={final_dimension}|target_rank={post_dim-2}|target_h0=2|"
    f"status={'PASS_H0_TWO' if final_dimension==2 else 'DIMENSION_MISMATCH'}",
    flush=True,
)

final58=kernel10*K10
assert final58.nrows()==final_dimension
assert final58.ncols()==len(ambient)
assert C*final58.transpose()==matrix(F,48,final_dimension)

transform=kernel10
condition_ledger=component_ledger
final_codimension=resolved_rank
final_pairs=[]


# ===========================================================================
# 4. If h0=2, compile the component-valuation-selected pencil.
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

component_valuation_certified=True
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
        "Q24DIVVAL_QUARTIC|"
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
             if component_valuation_certified else
             "CANDIDATE_H3_Q24_EFFECTIVE_D13_D12_MODP")
            if is_d12 else "GENUS_ONE_CHILD_NOT_D12"
        )
        print(
            "Q24DIVVAL_CHILD|"
            f"root_rank={root_rank}|root_det={root_det}|euler={euler}|"
            f"infinity={inf_orders},{inf_kind}|"
            f"status={'PASS_D12' if is_d12 else 'NOT_D12'}",
            flush=True,
        )
    else:
        terminal_status="RESOLVED_CLUSTER_NOT_GENUS_ONE"

payload={
    "schema":"elkies-k3.h3-q24-d12-component-valuation-rr-modp.v1",
    "status":terminal_status,
    "prime":int(p),
    "inputs":{
        "preflight":str(PREFLIGHT.relative_to(ROOT)),
        "resolution":str(RESOLUTION.relative_to(ROOT)),
        "component_graph":str(GRAPH.relative_to(ROOT)),
        "effective_d13_transport":str(EFFECTIVE.relative_to(ROOT)),
    },
    "geometric_trivialization":{"deterministic_fibre_twist":-7,"geometric_fibre_twist":int(geometric_fibre_twist),"vertical_square":int(geometric_vertical_square),"isotropic_square":0},
    "global_rr":{
        "ambient_dimension":int(len(ambient)),
        "smooth_collision_rank":int(collision_rank),
        "post_collision_dimension":int(post_dim),
    },
    "resolved_cluster":{
        "plan":[{"center":name,"additional_order":order} for name,order in plan],
        "condition_ledger":condition_ledger,
        "method":"direct_membership_in_(surface,e^threshold)_on_full_chart_cover",
        "codimension_on_post_collision":int(final_codimension),
        "kernel_dimension":int(final_dimension),
        "kernel_basis_post_collision":[
            [int(v) for v in row] for row in transform.rows()
        ],
        "kernel_basis_ambient":[
            [int(v) for v in row] for row in final58.rows()
        ],
    },
    "quartic_degree":quartic_degree,
    "child":child_summary,
    "proof_boundary":(
        "The geometric affine-I9* profile has vertical square -12. Together "
        "with P.O=24 and D^2=0 this forces fibre twist -8, so the correct "
        "infinity ambient is 56-dimensional (A<=40,B<=14), not the earlier "
        "58-dimensional deterministic-frame ambient. The I9* conditions use "
        "symbolic/generic-point component valuations. A final h0=2 plus "
        "genus-one quartic and D12 child is a modular resolved-chart "
        "certificate; characteristic-zero replay remains separate."
    ),
    "affine_D12_component_profile_certified":True,
    "geometric_Weyl_chamber_assumption_used":False,
}
OUT=(args.output.resolve() if args.output else
     LOCAL/f"q24-d12-component-valuation-rr-mod-{p}.json")
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True,default=int)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24DIVVAL_RESULT|"
    f"resolved_codim={final_codimension}|kernel={final_dimension}|"
    f"quartic_degree={quartic_degree}|status={terminal_status}",
    flush=True,
)
