#!/usr/bin/env sage -python
"""
Search the entire nearby q6 MW lattice for a cheaper q8 multisection whose
anchored D13 Abel-Jacobi class differs from q24 only by the already-explicit
G1 and G3.

Anchored D13 data:
    G1  = (1,0,0,0)   explicit
    G3  = (0,0,1,0)   explicit
    q24 = (2,-1,-1,1)

Therefore any divisor with D13 MW coordinates
    z = (*,-1,*,1)
can be converted to q24 using only G1 and G3:
    q24 = z + a*G1 + b*G3.

We enumerate actual effective sections of the preceding q6 fibration, not
formal MW vectors.  For each q6 MW word n=(n1,n2,n3), the script:

  * forms its exact Shioda vector from the already-certified q6 MW basis;
  * tries the three possible IV* component specializations
      affine, outer+, outer-
    (II* contributes no component choice);
  * reconstructs the unique integral (-2) section class;
  * verifies its q6 section properties and effective component pairings;
  * computes its q8 degree and anchored D13 MW coordinates.

This avoids a CVP call per lattice point and makes a fairly large exact search
cheap.  It reports:
  * cheapest section carrying any fourth D13 coordinate;
  * cheapest q24-compatible section z=(*,-1,*,1);
  * cheapest exact q24 section, if one appears.

Run:
  sage -python ~/Downloads/search_h92_q6_mw_for_q24_bridge.sage
Optional:
  sage -python ... --radius 100 --third-radius 3
"""

import argparse
import json
from math import isqrt
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, identity_matrix, matrix, pari, vector, xgcd
)

Q6_REFLECTIONS = (
    1, 2, 4, 3, 5, 4, 2, 6, 5, 4, 3, 1,
    7, 6, 5, 4, 2, 3, 4, 5, 6, 7,
)

H3_LIFTS = matrix(ZZ, [
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -4, 1, 0, -4, 2, -2],
    [-10, -8, -6, 0, 0, 0, 0, 0, 0, 0, 0, -8, 4, 1, -8, 5, -4],
    [-5, -4, -3, 0, 0, 0, 0, 0, 0, 0, 0, -3, 2, 0, -4, 2, -2],
])

OLD_ZERO_ROOT_SHIFTS = matrix(ZZ, [
    [5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, -1, 4],
    [12, 10, 8, 0, 0, 0, 0, 0, 0, 0, 0, 6, -1, 9],
    [5, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 4],
])

EXPECTED_Q6_HEIGHT = matrix(QQ, [
    [QQ(8)/3, QQ(1)/3, -1],
    [QQ(1)/3, QQ(8)/3, 1],
    [-1, 1, 46],
])

E6_QF_INDICES = (0, 1, 2, 3, 12, 13)
E6_SIMPLE_IN_QF = (
    (-1, -1, -1, -1, 0, 0),
    (0, 0, 0, 1, 0, 0),
    (0, 0, 1, 0, 0, 0),
    (0, 1, 0, 0, 0, 0),
    (1, 0, 0, 0, 0, 1),
    (2, 1, 0, 0, -1, 1),
)
E6_CARTAN = matrix(ZZ, [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, -1],
    [0, 0, -1, 2, -1, 0],
    [0, 0, 0, -1, 2, 0],
    [0, 0, -1, 0, 0, 2],
])
E8_CARTAN = matrix(ZZ, [
    [2, 0, -1, 0, 0, 0, 0, 0],
    [0, 2, 0, -1, 0, 0, 0, 0],
    [-1, 0, 2, -1, 0, 0, 0, 0],
    [0, -1, -1, 2, -1, 0, 0, 0],
    [0, 0, 0, -1, 2, -1, 0, 0],
    [0, 0, 0, 0, -1, 2, -1, 0],
    [0, 0, 0, 0, 0, -1, 2, -1],
    [0, 0, 0, 0, 0, 0, -1, 2],
])


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    home = Path.home()
    candidates += [
        home / "Documents" / "jacobian-research",
        home / "jacobian-research",
        home / "src" / "jacobian-research",
        home / "git" / "jacobian-research",
        home / "projects" / "jacobian-research",
    ]
    seen = set()
    for candidate in candidates:
        try:
            candidate = candidate.resolve()
        except Exception:
            continue
        if candidate in seen:
            continue
        seen.add(candidate)
        if (candidate / "elkies-k3/scripts").is_dir():
            return candidate
    raise SystemExit("Could not locate jacobian-research; pass --repo PATH")


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


def reflect(row, gram, root):
    row = vector(ZZ, list(row))
    root = vector(ZZ, list(root))
    assert root * gram * root == -2
    return row + (row * gram * root) * root


def isotropic_mate(ns, fibre):
    current = ZZ(0)
    data = [ZZ(0)] * ns.nrows()
    for index, value in enumerate(ns * fibre):
        if not value:
            continue
        divisor, left, right = xgcd(current, ZZ(value))
        data = [left * entry for entry in data]
        data[index] += right
        current = divisor
    assert abs(current) == 1
    if current == -1:
        data = [-entry for entry in data]
    mate = vector(ZZ, data)
    mate -= (mate * ns * mate // 2) * fibre
    assert mate * ns * mate == 0 and mate * ns * fibre == 1
    return mate


def child_frame_with_zero(ns, fibre, zero):
    mate = zero + fibre
    orth = matrix(
        ZZ, [list(fibre * ns), list(mate * ns)]
    ).right_kernel_matrix()
    basis = matrix(
        ZZ, [list(fibre), list(mate)] + [list(row) for row in orth.rows()]
    )
    assert abs(basis.det()) == 1
    child = -(orth * ns * orth.transpose())
    U = matrix(ZZ, ((0,1),(1,0)))
    assert basis * ns * basis.transpose() == block_diagonal_matrix(U, -child)
    return child, basis


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    half = [vector(ZZ,c) for c in matrix(ZZ,result[2]).columns()]
    roots = tuple(half + [-r for r in half])
    rb = matrix(ZZ,[list(r) for r in roots]).row_module().basis_matrix()
    rg = rb*gram*rb.transpose()
    return roots, rb, (rb.rank(),count,abs(ZZ(rg.det())))


def deterministic_simple_roots(gram):
    roots, unused, data = roots_and_data(gram)
    rank = data[0]
    regular=None
    for shift in range(1,1000):
        candidate=vector(ZZ,[
            (i+1)**2+shift*(i+1)+1 for i in range(gram.nrows())
        ])
        if all(candidate*root != 0 for root in roots):
            regular=candidate
            break
    assert regular is not None
    positive=[r for r in roots if regular*r > 0]
    pset={tuple(r) for r in positive}
    simple=[
        r for r in positive
        if not any(tuple(r-left) in pset for left in positive)
    ]
    simple=matrix(ZZ,[list(r) for r in simple])
    assert simple.nrows()==simple.rank()==rank
    return simple,simple*gram*simple.transpose()


def d13_root_adaptation(child):
    unused,root_basis,invariants=roots_and_data(child)
    assert invariants==(13,312,4),invariants
    simple,cartan=deterministic_simple_roots(child)
    smith,left,right=root_basis.smith_form()
    assert smith==left*root_basis*right
    completion=right.inverse()
    initial=simple.stack(completion[13:])
    adapted=initial*child*initial.transpose()
    root=adapted[:13,:13]
    coupling=adapted[:13,13:]
    tail=adapted[13:,13:]
    H=tail-coupling.transpose()*root.inverse()*coupling
    scale=ZZ(1)
    from sage.all import lcm
    for v in H.list():
        scale=lcm(scale,ZZ(QQ(v).denominator()))
    lll=matrix(ZZ,pari((scale*H).change_ring(ZZ)).qflllgram())
    change=block_diagonal_matrix(identity_matrix(ZZ,13),lll.transpose())
    basis=change*initial
    adapted=basis*child*basis.transpose()
    root=adapted[:13,:13]
    coupling=adapted[:13,13:]
    tail=adapted[13:,13:]
    H=tail-coupling.transpose()*root.inverse()*coupling
    assert H.det()==237
    return basis,adapted,H


def highest_root(cartan):
    half=matrix(ZZ,pari(cartan).qfminim(2)[2]).transpose().rows()
    roots=[vector(ZZ,r) for r in half]
    roots += [-vector(ZZ,r) for r in half]
    positive=[r for r in roots if all(v>=0 for v in r)]
    return max(positive,key=lambda r:sum(r))


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--threshold",type=int,default=206,
                    help="exhaustively enumerate every q6 section with q8 degree below this")
parser.add_argument("--top",type=int,default=20)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
FRAME=ROOT/"elkies-k3/data/fibrations/kumar_e7e8_mw2_frame_3.txt"
TARGET=LOCAL/"q8-target-component-nef.json"
BRANCH=LOCAL/"q8-d13-branch-anchor.json"
G3FILE=LOCAL/"q8-d13-g3-from-e77-bisection.json"
OUTPUT=(
    args.output.resolve()
    if args.output
    else LOCAL/"q6-mw-q24-bridge-search.json"
)

for path in (FRAME,TARGET,BRANCH,G3FILE):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

frame=load_gram(FRAME)
target=json.loads(TARGET.read_text())
branch=json.loads(BRANCH.read_text())
g3data=json.loads(G3FILE.read_text())
assert branch["status"]=="PASS_EXACT_D13_BRANCH_ANCHOR"
assert g3data["status"]=="PASS_EXACT_D13_G3_FROM_E77_BISECTION"

U2=matrix(ZZ,((0,1),(1,0)))
ns=block_diagonal_matrix(U2,-frame)
source_F=vector(ZZ,[1,0]+[0]*17)
source_O=vector(ZZ,[-1,1]+[0]*17)
source_simple=tuple(
    vector(ZZ,[0,0]+[ZZ(i==node) for i in range(17)])
    for node in range(15)
)

# ---------------------------------------------------------------------------
# 1. Reconstruct actual q6 fibration and the three certified basis sections.
# ---------------------------------------------------------------------------

raw_F=vector(ZZ,[3,2]+[
    0,0,-1,-1,-1,-1,-1,
    0,0,0,0,0,0,0,0,1,0
])
reflection_roots=tuple(source_simple[node-1] for node in Q6_REFLECTIONS)

q6_F=vector(ZZ,list(raw_F))
for root in reflection_roots:
    q6_F=reflect(q6_F,ns,root)
assert source_O*ns*q6_F==1

raw_zero=vector(ZZ,list(source_O))
for root in reversed(reflection_roots):
    raw_zero=reflect(raw_zero,ns,root)

raw_mate=isotropic_mate(ns,raw_F)
raw_orth=matrix(
    ZZ,[list(raw_F*ns),list(raw_mate*ns)]
).right_kernel_matrix()
raw_transport=matrix(
    ZZ,[list(raw_F),list(raw_mate)]+[list(row) for row in raw_orth.rows()]
)
raw_child=-(raw_orth*ns*raw_orth.transpose())
raw_roots=matrix(
    ZZ,pari(raw_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
raw_zero_coords=vector(ZZ,raw_zero*raw_transport.inverse())
raw_zero_lift=vector(ZZ,raw_zero_coords[2:])

basis_sections=[]
for h3_lift,shift in zip(H3_LIFTS.rows(),OLD_ZERO_ROOT_SHIFTS.rows()):
    lift=(
        raw_zero_lift
        +vector(ZZ,h3_lift)
        +vector(ZZ,shift)*raw_roots
    )
    norm=ZZ(lift*raw_child*lift)
    pole=norm//2-2
    candidate_raw=(pole+1)*raw_F+raw_mate+lift*raw_orth
    candidate=vector(ZZ,list(candidate_raw))
    for root in reflection_roots:
        candidate=reflect(candidate,ns,root)
    assert candidate*ns*candidate==-2
    assert candidate*ns*q6_F==1
    basis_sections.append(candidate)

# ---------------------------------------------------------------------------
# 2. Physical q6 E6+E8 simple components and Shioda basis.
# ---------------------------------------------------------------------------

q6_orth=matrix(
    ZZ,
    [
        list(q6_F*ns),
        list((source_O+q6_F)*ns),
    ],
).right_kernel_matrix()
physical_child=-(q6_orth*ns*q6_orth.transpose())
qf_basis=matrix(
    ZZ,pari(physical_child).qfminim(2)[2]
).transpose().row_module().basis_matrix()
assert qf_basis.rank()==14

e6_qf=matrix(ZZ,[list(qf_basis[i]) for i in E6_QF_INDICES])
e6_roots=matrix(
    ZZ,[list(vector(ZZ,row)*e6_qf) for row in E6_SIMPLE_IN_QF]
)*q6_orth
physical_roots=qf_basis*q6_orth
e8_roots=physical_roots[4:12,:]

assert -e6_roots*ns*e6_roots.transpose()==E6_CARTAN
assert -e8_roots*ns*e8_roots.transpose()==E8_CARTAN
assert e6_roots*ns*e8_roots.transpose()==matrix(ZZ,6,8)

Rroots=e6_roots.stack(e8_roots)
Groot=Rroots*ns*Rroots.transpose()
assert Groot.det()==3  # (-1)^14 does not change det

projection=identity_matrix(QQ,19)-ns*Rroots.transpose()*Groot.inverse()*Rroots

def shioda(P):
    horizontal=P-source_O-(P*ns*source_O+2)*q6_F
    assert horizontal*ns*q6_F==0
    assert horizontal*ns*source_O==0
    return vector(QQ,horizontal)*projection

phis=[shioda(P) for P in basis_sections]
Hq6=matrix(QQ,[[-a*ns*b for b in phis] for a in phis])
assert Hq6==EXPECTED_Q6_HEIGHT

# Section component patterns in IV*: affine plus the simple multiplicity-one
# components of the E6 highest root.  II* has no non-affine multiplicity-one
# simple component.
e6_high=highest_root(E6_CARTAN)
outer_indices=[i for i,v in enumerate(e6_high) if v==1]
assert len(outer_indices)==2, (e6_high,outer_indices)

pairing_patterns=[vector(ZZ,[0]*14)]
for i in outer_indices:
    p=[0]*14
    p[i]=1
    pairing_patterns.append(vector(ZZ,p))

# Verify each known basis section is reconstructed by exactly one pattern.
def candidate_for_pattern(n, p):
    n=vector(ZZ,n)
    p=vector(ZZ,p)
    phi=sum((QQ(n[i])*phis[i] for i in range(3)),vector(QQ,[0]*19))
    height=QQ(n*Hq6*n)
    correction=-QQ(p*Groot.inverse()*p)
    pole=(height+correction-4)/2
    if pole not in ZZ or pole<0:
        return None
    root_projection=vector(QQ,p)*Groot.inverse()*Rroots
    P=(
        vector(QQ,source_O)
        +(QQ(pole)+2)*vector(QQ,q6_F)
        +phi+root_projection
    )
    if not all(v in ZZ for v in P):
        return None
    P=vector(ZZ,P)
    if P*ns*P!=-2 or P*ns*q6_F!=1:
        return None
    if vector(ZZ,P*ns*Rroots.transpose())!=p:
        return None
    return P,ZZ(pole),correction,p

def section_from_word(n):
    solutions=[]
    for p in pairing_patterns:
        candidate=candidate_for_pattern(n,p)
        if candidate is not None:
            solutions.append(candidate)
    assert len(solutions)==1,(tuple(n),len(solutions))
    return solutions[0]

for i in range(3):
    n=[0,0,0]
    n[i]=1
    P,pole,corr,p=section_from_word(n)
    assert P==basis_sections[i],(i,P-basis_sections[i])

print(
    "Q6Q24SEARCH_SETUP|basis_reconstruction=PASS|"
    f"IV_outer_indices={','.join(str(i+1) for i in outer_indices)}|"
    "q6_height=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# 3. Anchored D13 coordinate reader.
# ---------------------------------------------------------------------------

F8=vector(ZZ,target["selected_q8"]["source_h3_ns_vector"])
E6q8=matrix(
    ZZ,target["selected_q8"]["E6"]["simple_root_vectors_in_source_h3_ns"]
)
E8q8=matrix(
    ZZ,target["selected_q8"]["E8"]["simple_root_vectors_in_source_h3_ns"]
)
O8=vector(ZZ,E8q8.row(0))

child8,Bzero=child_frame_with_zero(ns,F8,O8)
A13,adapted,Hd13=d13_root_adaptation(child8)
Badapt=block_diagonal_matrix(identity_matrix(ZZ,2),A13)*Bzero
Gadapt=block_diagonal_matrix(U2,-adapted)
assert Badapt*ns*Badapt.transpose()==Gadapt

def d13_coords(P):
    c=vector(QQ,P)*Badapt.inverse()
    assert all(v in ZZ for v in c)
    return vector(ZZ,c)

G1=vector(ZZ,(1,0,0,0))
G3=vector(ZZ,(0,0,1,0))
q24=vector(ZZ,(2,-1,-1,1))
assert G1*Hd13*G1==QQ(3)/4
assert G3*Hd13*G3==QQ(11)/4
assert q24*Hd13*q24==52

# ---------------------------------------------------------------------------
# 4. Derive a certified exhaustive box for every q6 section with
#    q8 degree < threshold.
# ---------------------------------------------------------------------------

# For a fixed IV* component pattern p,
#
# P(n,p)=O + (height(n)+corr(p))/2 * F6 + phi(n) + root_projection(p),
#
# because pole+2=(height+corr)/2.  Since F6.F8=2,
#
# d_q8(n,p)=n^T H n + linear*n + constant_p.
#
# Thus completing the square gives a rigorous finite ellipsoid.
assert q6_F*ns*F8==2
linear=vector(QQ,[phi*ns*F8 for phi in phis])
Hinv=Hq6.inverse()
center=-QQ(1)/2 * linear * Hinv

def ceil_sqrt_rational(value):
    value=QQ(value)
    assert value>=0
    num=ZZ(value.numerator())
    den=ZZ(value.denominator())
    q=(num+den-1)//den
    k=ZZ(isqrt(int(q)))
    if k*k<q:
        k+=1
    assert QQ(k*k)>=value
    return k

pattern_quadratics=[]
global_lower=[None]*3
global_upper=[None]*3
threshold=QQ(args.threshold)

for pattern_index,p in enumerate(pairing_patterns):
    p=vector(ZZ,p)
    correction=-QQ(p*Groot.inverse()*p)
    root_projection=vector(QQ,p)*Groot.inverse()*Rroots
    constant=QQ(source_O*ns*F8)+correction+QQ(root_projection*ns*F8)
    qminimum=constant-QQ(1)/4*(linear*Hinv*linear)
    budget=threshold-qminimum
    if budget<=0:
        continue
    bounds=[]
    for i in range(3):
        radius=ceil_sqrt_rational(budget*Hinv[i,i])
        lo=ZZ(center[i].floor())-radius-1
        hi=ZZ(center[i].ceil())+radius+1
        bounds.append((lo,hi))
        global_lower[i]=lo if global_lower[i] is None else min(global_lower[i],lo)
        global_upper[i]=hi if global_upper[i] is None else max(global_upper[i],hi)
    pattern_quadratics.append({
        "pattern_index":pattern_index,
        "pattern":p,
        "correction":correction,
        "constant":constant,
        "minimum":qminimum,
        "bounds":bounds,
    })

assert pattern_quadratics
assert all(v is not None for v in global_lower+global_upper)

print(
    "Q6Q24SEARCH_BOUND|threshold={}|center={}|box={}|patterns={}".format(
        args.threshold,
        ",".join(map(str,center)),
        ";".join(f"{global_lower[i]}..{global_upper[i]}" for i in range(3)),
        len(pattern_quadratics),
    ),
    flush=True,
)

# Verify the closed formula against the three known basis sections.
for i in range(3):
    n=vector(ZZ,[ZZ(j==i) for j in range(3)])
    P,pole,corr,p=section_from_word(n)
    formula=QQ(n*Hq6*n)+linear*n+(
        QQ(source_O*ns*F8)+corr+
        QQ((vector(QQ,p)*Groot.inverse()*Rroots)*ns*F8)
    )
    assert formula==P*ns*F8

best_fourth=[]
best_q24=[]
exact=[]
tested=0
below_threshold=0

# Enumerate the union box.  candidate_for_pattern / section_from_word certifies
# which of the three component patterns is actually integral for each MW word.
for n1 in range(int(global_lower[0]),int(global_upper[0])+1):
    for n2 in range(int(global_lower[1]),int(global_upper[1])+1):
        for n3 in range(int(global_lower[2]),int(global_upper[2])+1):
            if n1==n2==n3==0:
                continue
            n=vector(ZZ,(n1,n2,n3))
            P,pole,corr,p=section_from_word(n)
            q8deg=int(P*ns*F8)
            tested+=1
            if q8deg>=args.threshold:
                continue
            below_threshold+=1

            # Independent check against the appropriate fixed-pattern quadratic.
            root_projection=vector(QQ,p)*Groot.inverse()*Rroots
            formula=QQ(n*Hq6*n)+linear*n+(
                QQ(source_O*ns*F8)+corr+QQ(root_projection*ns*F8)
            )
            assert formula==q8deg

            c=d13_coords(P)
            z=vector(ZZ,c[-4:])
            rec={
                "q6_word":[n1,n2,n3],
                "q6_height":str(n*Hq6*n),
                "q6_P_dot_O":int(P*ns*source_O),
                "q6_component_correction":str(corr),
                "q8_degree":q8deg,
                "q8_zero_intersection":int(P*ns*O8),
                "d13_mw":list(map(int,z)),
                "source_h3_ns":list(map(int,P)),
            }

            if z[3]!=0:
                best_fourth.append(rec)

            if z[1]==-1 and z[3]==1:
                residual=q24-z
                assert residual[1]==0 and residual[3]==0
                rec["q24_correction_G1"]=int(residual[0])
                rec["q24_correction_G3"]=int(residual[2])
                best_q24.append(rec)

            if z==q24:
                exact.append(rec)

best_fourth.sort(key=lambda r:(r["q8_degree"],sum(abs(v) for v in r["q6_word"])))
best_q24.sort(key=lambda r:(r["q8_degree"],sum(abs(v) for v in r["q6_word"])))
exact.sort(key=lambda r:r["q8_degree"])

print(
    f"Q6Q24SEARCH|tested_box={tested}|below_threshold={below_threshold}|"
    f"threshold={args.threshold}|"
    f"fourth_candidates={len(best_fourth)}|"
    f"q24_compatible={len(best_q24)}|exact_q24={len(exact)}",
    flush=True,
)

def emit(prefix,records):
    for rank,rec in enumerate(records[:args.top],1):
        print(
            f"{prefix}|rank={rank}|"
            f"word={','.join(map(str,rec['q6_word']))}|"
            f"q8_degree={rec['q8_degree']}|"
            f"d13_mw={','.join(map(str,rec['d13_mw']))}|"
            f"q6_height={rec['q6_height']}|"
            f"q6_PdotO={rec['q6_P_dot_O']}|"
            +(
                f"q24_add={rec['q24_correction_G1']}*G1+"
                f"{rec['q24_correction_G3']}*G3|"
                if "q24_correction_G1" in rec else ""
            )
            +"status=CANDIDATE",
            flush=True,
        )

emit("Q6Q24SEARCH_FOURTH",best_fourth)
emit("Q6Q24SEARCH_COMPAT",best_q24)
emit("Q6Q24SEARCH_EXACT",exact)

status=(
    "PASS_FOUND_GLOBALLY_CHEAPER_Q24_BRIDGE"
    if best_q24
    else "PASS_NO_Q6_Q24_BRIDGE_BELOW_THRESHOLD"
)

print(
    "Q6Q24SEARCH_RESULT|status={}|best_compatible_degree={}|"
    "threshold={}|reference_C_degree=206|exhaustive_below_threshold=1".format(
        status,
        best_q24[0]["q8_degree"] if best_q24 else "none",
        args.threshold,
    ),
    flush=True,
)

payload={
    "schema":"elkies-k3.h92-q6-mw-q24-bridge-search.v1",
    "status":status,
    "search":{
        "threshold":args.threshold,
        "tested_box":tested,
        "below_threshold":below_threshold,
        "completed_square_center":[str(v) for v in center],
        "global_coordinate_bounds":[
            [int(global_lower[i]),int(global_upper[i])] for i in range(3)
        ],
        "pattern_quadratics":[
            {
                "pattern_index":item["pattern_index"],
                "pairing_pattern":list(map(int,item["pattern"])),
                "component_correction":str(item["correction"]),
                "constant":str(item["constant"]),
                "real_minimum":str(item["minimum"]),
                "coordinate_bounds":[[int(a),int(b)] for a,b in item["bounds"]],
            }
            for item in pattern_quadratics
        ],
        "exhaustive_for_q8_degree_below_threshold":True,
    },
    "target":{
        "q24_anchored_D13_mw":[2,-1,-1,1],
        "explicit_correction_basis":{
            "G1":[1,0,0,0],
            "G3":[0,0,1,0],
        },
        "compatibility_condition":"z[1]=-1 and z[3]=1",
        "reference_C":{
            "mw":[22,-1,45,1],
            "q8_degree":206,
            "relation":"q24=C-20*G1-46*G3",
        },
    },
    "best_fourth_direction":best_fourth[:args.top],
    "best_q24_compatible":best_q24[:args.top],
    "exact_q24":exact[:args.top],
    "boundary":(
        "This is an exact exhaustive search over every actual effective q6 section "
        "with q8 degree below the stated threshold. It identifies cheaper q8 multisection bridges but does not "
        "perform their Abel-Jacobi reduction on the q8 quartic."
    ),
}
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUTPUT}")
