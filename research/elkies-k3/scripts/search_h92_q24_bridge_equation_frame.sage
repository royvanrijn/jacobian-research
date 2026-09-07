#!/usr/bin/env sage -python
"""
Exhaustive equation-level H92 q24 bridge search.

This replaces the historical physical/raw-q6 degree search.

Parent chain is now authoritative at equation level:
    E7+E8/MW2 -> E8+E6/MW3 -> D13/MW4.

We enumerate actual effective q6 sections in the current (Weyl-transported)
q6 chamber, measure degree against the exact equation q8 fibre F8eq, and read
their Abel-Jacobi classes in the deterministic equation D13 frame.

The q24 structural-isometry audit proved that in this same deterministic
equation D13 basis the q24 horizontal class is unique up to inversion:
    +/-(-2,1,-1,1).

Known equation D13 corrections G1,G3 are reconstructed geometrically:
    G1 = IV*_E6_1 q8 section,
    AJ(old_E7_7) = G1 + G3.

A candidate is q24-compatible iff, for either q24 orientation,
    q24 - AJ(candidate) in Z*G1 + Z*G3.

The q6 degree search is rigorously exhaustive below --threshold by completing
the positive-definite height quadratic independently for each IV* component
pattern, exactly as in the corrected raw-q6 search.
"""

import argparse
import json
import sys
from math import isqrt
from pathlib import Path

from sage.all import (
    QQ, ZZ, block_diagonal_matrix, identity_matrix, lcm,
    matrix, pari, vector
)


E6_CARTAN = matrix(ZZ, [
    [2,-1,0,0,0,0],
    [-1,2,-1,0,0,0],
    [0,-1,2,-1,0,-1],
    [0,0,-1,2,-1,0],
    [0,0,0,-1,2,0],
    [0,0,-1,0,0,2],
])
E8_CARTAN = matrix(ZZ, [
    [2,0,-1,0,0,0,0,0],
    [0,2,0,-1,0,0,0,0],
    [-1,0,2,-1,0,0,0,0],
    [0,-1,-1,2,-1,0,0,0],
    [0,0,0,-1,2,-1,0,0],
    [0,0,0,0,-1,2,-1,0],
    [0,0,0,0,0,-1,2,-1],
    [0,0,0,0,0,0,-1,2],
])


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
        try: c=c.resolve()
        except Exception: continue
        if c in seen: continue
        seen.add(c)
        if (c/"elkies-k3/scripts").is_dir() and (c/"artifacts/generated-results").is_dir():
            return c
    raise SystemExit("Could not locate repo")


def run_scope(path):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    try:
        sys.argv=[str(path)]
        exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope


def roots_and_data(gram):
    qf=pari(gram).qfminim(2)
    count=ZZ(qf[0])
    if not count:
        return (),matrix(ZZ,0,gram.nrows()),(0,0,1)
    half=[vector(ZZ,c) for c in matrix(ZZ,qf[2]).columns()]
    roots=tuple(half+[-r for r in half])
    basis=matrix(ZZ,[list(r) for r in roots]).row_module().basis_matrix()
    rg=basis*gram*basis.transpose()
    return roots,basis,(basis.rank(),count,abs(ZZ(rg.det())))


def deterministic_simple_roots(gram):
    roots,unused,data=roots_and_data(gram)
    rank=data[0]
    regular=None
    for shift in range(1,1000):
        c=vector(ZZ,[(i+1)**2+shift*(i+1)+1 for i in range(gram.nrows())])
        if all(c*r != 0 for r in roots):
            regular=c
            break
    assert regular is not None
    positive=[r for r in roots if regular*r>0]
    pset={tuple(r) for r in positive}
    simple=[
        r for r in positive
        if not any(tuple(r-left) in pset for left in positive)
    ]
    simple=matrix(ZZ,[list(r) for r in simple])
    assert simple.nrows()==simple.rank()==rank
    return simple,simple*gram*simple.transpose()


def child_frame_with_zero(ns,fibre,zero):
    mate=zero+fibre
    assert fibre*ns*fibre==0
    assert zero*ns*zero==-2 and zero*ns*fibre==1
    assert mate*ns*mate==0 and mate*ns*fibre==1
    orth=matrix(ZZ,[list(fibre*ns),list(mate*ns)]).right_kernel_matrix()
    B=matrix(ZZ,[list(fibre),list(mate)]+[list(r) for r in orth.rows()])
    assert abs(B.det())==1
    child=-(orth*ns*orth.transpose())
    return child,B


def d13_root_adaptation(child):
    unused,rb,inv=roots_and_data(child)
    assert inv==(13,312,4)
    simple,cartan=deterministic_simple_roots(child)
    assert abs(cartan.det())==4

    smith,left,right=rb.smith_form()
    assert smith==left*rb*right
    assert tuple(abs(smith[i,i]) for i in range(13))==(1,)*13
    completion=right.inverse()
    initial=simple.stack(completion[13:])
    assert abs(initial.det())==1

    G=initial*child*initial.transpose()
    C=G[:13,13:]
    T=G[13:,13:]
    H=T-C.transpose()*cartan.inverse()*C

    scale=ZZ(1)
    for x in H.list():
        scale=lcm(scale,ZZ(QQ(x).denominator()))
    lll=matrix(ZZ,pari((scale*H).change_ring(ZZ)).qflllgram())
    assert abs(lll.det())==1

    change=block_diagonal_matrix(identity_matrix(ZZ,13),lll.transpose())
    A=change*initial
    G=A*child*A.transpose()

    R=G[:13,:13]
    C=G[:13,13:]
    T=G[13:,13:]
    H=T-C.transpose()*R.inverse()*C
    assert H.det()==237
    return A,G,H


def highest_root(cartan):
    half=matrix(ZZ,pari(cartan).qfminim(2)[2]).transpose().rows()
    roots=[vector(ZZ,r) for r in half]
    roots += [-vector(ZZ,r) for r in half]
    positive=[r for r in roots if all(v>=0 for v in r)]
    assert positive
    return max(positive,key=lambda r:sum(r))


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


def solve_two_span(a,b,r):
    """Return integral (x,y) with x*a+y*b=r, or None."""
    A=matrix(QQ,[list(a),list(b)]).transpose()
    try:
        sol=A.solve_right(vector(QQ,r))
    except ValueError:
        return None
    if all(v in ZZ for v in sol):
        sol=vector(ZZ,sol)
        if sol[0]*a+sol[1]*b==r:
            return (int(sol[0]),int(sol[1]))
    return None


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--threshold",type=int,default=52)
parser.add_argument("--top",type=int,default=20)
parser.add_argument("--output",type=Path)
parser.add_argument(
    "--setup-only",
    action="store_true",
    help=(
        "construct the exact q6 word machinery and named equation-D13 "
        "frame, then stop before the bounded bridge enumeration"
    ),
)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
CERT=ROOT/"elkies-k3/scripts/certify_h92_q8_equation_ns_divisor.sage"
OUT=args.output.resolve() if args.output else LOCAL/"q24-equation-bridge-search.json"

eq=run_scope(CERT)
ns=eq["ns"]
F6=vector(ZZ,eq["F6"])
O6=vector(ZZ,eq["Oold"])
Ostd=vector(ZZ,eq["Ostd"])
F8=vector(ZZ,eq["F8eq"])
H6=matrix(QQ,eq["EXPECTED_HEIGHT"])
phis=[vector(QQ,p) for p in eq["phis"]]
simple=eq["simple"]
weyl=eq["weyl_transport"]
selected=eq["selected"]

assert F6*ns*F6==0
assert O6*ns*O6==-2 and O6*ns*F6==1
assert F8*ns*F8==0 and F8*ns*F6==2

# ---------------------------------------------------------------------------
# 1. Current q6 E6+E8 roots: transport the exact physical simple roots through
#    the same 22 q6 Weyl reflections that define the current q6 chamber.
# ---------------------------------------------------------------------------

e6_raw=matrix(ZZ,selected["E6"]["simple_root_vectors_in_source_h3_ns"])
e8_raw=matrix(ZZ,selected["E8"]["simple_root_vectors_in_source_h3_ns"])
e6=matrix(ZZ,[list(weyl(vector(ZZ,r))) for r in e6_raw.rows()])
e8=matrix(ZZ,[list(weyl(vector(ZZ,r))) for r in e8_raw.rows()])

assert -e6*ns*e6.transpose()==E6_CARTAN
assert -e8*ns*e8.transpose()==E8_CARTAN
assert e6*ns*e8.transpose()==matrix(ZZ,6,8)
assert all(r*ns*F6==0 for r in e6.rows())
assert all(r*ns*F6==0 for r in e8.rows())
assert all(r*ns*O6==0 for r in e6.rows())
assert all(r*ns*O6==0 for r in e8.rows())

Rroots=e6.stack(e8)
Groot=Rroots*ns*Rroots.transpose()
assert Groot.det()==3

# The certifier's MW Shioda basis must be orthogonal to the same root lattice.
assert all(phi*ns*F6==0 and phi*ns*O6==0 for phi in phis)
assert all(
    phi*ns*r==0
    for phi in phis
    for r in Rroots.rows()
)
assert matrix(QQ,[[-a*ns*b for b in phis] for a in phis])==H6

# IV* multiplicity-one component patterns.
e6_high=highest_root(E6_CARTAN)
outer_indices=[i for i,v in enumerate(e6_high) if v==1]
assert len(outer_indices)==2
patterns=[vector(ZZ,[0]*14)]
for i in outer_indices:
    p=[0]*14
    p[i]=1
    patterns.append(vector(ZZ,p))


def candidate_for_pattern(n,p):
    n=vector(ZZ,n); p=vector(ZZ,p)
    phi=sum((ZZ(n[i])*phis[i] for i in range(3)),vector(QQ,[0]*19))
    height=QQ(n*H6*n)
    correction=-QQ(p*Groot.inverse()*p)
    pole=(height+correction-4)/2
    if pole not in ZZ or pole<0:
        return None
    root_part=vector(QQ,p)*Groot.inverse()*Rroots
    P=(
        vector(QQ,O6)
        +(QQ(pole)+2)*vector(QQ,F6)
        +phi+root_part
    )
    if not all(x in ZZ for x in P):
        return None
    P=vector(ZZ,P)
    if P*ns*P!=-2 or P*ns*F6!=1:
        return None
    if vector(ZZ,P*ns*Rroots.transpose())!=p:
        return None
    return P,ZZ(pole),correction,p


def section_from_word(n):
    sol=[]
    for p in patterns:
        c=candidate_for_pattern(n,p)
        if c is not None:
            sol.append(c)
    assert len(sol)==1,(tuple(n),len(sol))
    return sol[0]


# Regressions on actual current q6 curves.
E77=vector(ZZ,simple[6])
assert E77*ns*E77==-2 and E77*ns*F6==1

def shioda(P):
    horizontal=P-O6-(P*ns*O6+2)*F6
    assert horizontal*ns*F6==0 and horizontal*ns*O6==0
    projection=identity_matrix(QQ,19)-ns*Rroots.transpose()*Groot.inverse()*Rroots
    return vector(QQ,horizontal)*projection

def word_of(P):
    phi=shioda(P)
    pair=vector(QQ,[-phi*ns*q for q in phis])
    w=pair*H6.inverse()
    assert all(x in ZZ for x in w)
    return vector(ZZ,w)

w77=word_of(E77)
assert section_from_word(w77)[0]==E77
S3=vector(ZZ,eq["S3std"])
wS3=word_of(S3)
assert section_from_word(wS3)[0]==S3

print(
    "Q24EQSEARCH_Q6|"
    f"E77_word={','.join(map(str,w77))}|E77_degree={E77*ns*F8}|"
    f"S3_word={','.join(map(str,wS3))}|S3_degree={S3*ns*F8}|"
    "status=PASS_EFFECTIVE_Q6_RECONSTRUCTION",
    flush=True,
)

# ---------------------------------------------------------------------------
# 2. Geometrically transported equation D13 frame and q24 target.
#
# Do NOT import q24 coordinates through the old structural D13 isometry.
# Instead start from the certified PHYSICAL branch-zero D13 frame, recover the
# actual q24 (-2)-section there, and transport every geometric anchor through
#
#       tau_after_weyl = tau o Weyl.
# ---------------------------------------------------------------------------

F8phys = vector(ZZ, selected["source_h3_ns_vector"])
O8phys = vector(ZZ, e8_raw.row(0))
G1phys_curve = vector(ZZ, e6_raw.row(0))

assert F8phys*ns*F8phys == 0
assert O8phys*ns*F8phys == 1
assert G1phys_curve*ns*F8phys == 1

# Exact q6 Eichler translation old zero -> standard zero.
ov = ZZ(Ostd*ns*O6)
vtrans = Ostd - O6 - (ov + 2)*F6
assert vtrans*ns*F6 == 0
assert vtrans*ns*vtrans == -12

def tau(x):
    x = vector(ZZ, x)
    return (
        x
        + (x*ns*F6)*vtrans
        - (x*ns*vtrans)*F6
        - (vtrans*ns*vtrans//2)*(x*ns*F6)*F6
    )

assert tau(F6) == F6
assert tau(O6) == Ostd
assert tau(weyl(F8phys)) == F8

def to_equation(x):
    return tau(weyl(vector(ZZ, x)))

# Physical D13 frame, with the original branch-zero orientation.
childP, BzeroP = child_frame_with_zero(ns, F8phys, O8phys)
AP, GP, HP = d13_root_adaptation(childP)
BadaptP = block_diagonal_matrix(identity_matrix(ZZ,2), AP) * BzeroP

def d13_coords_phys(P):
    c = vector(QQ, P) * BadaptP.inverse()
    assert all(x in ZZ for x in c)
    return vector(ZZ, c)

zG1phys = vector(ZZ, d13_coords_phys(G1phys_curve)[-4:])
assert zG1phys == vector(ZZ, (1,0,0,0))

# Physical old_E7_7 bisection is the inverse-Weyl image of source E7_7.
E77phys = eq["inverse_weyl_transport"](simple[6])
assert E77phys*ns*F8phys == 2
zAJ77phys = vector(ZZ, d13_coords_phys(E77phys)[-4:])
assert zAJ77phys == vector(ZZ, (1,0,1,0))
zG3phys = zAJ77phys - zG1phys
assert zG3phys == vector(ZZ, (0,0,1,0))
assert zG3phys*HP*zG3phys == QQ(11)/4

zQ24phys = vector(ZZ, (2,-1,-1,1))
assert zQ24phys*HP*zQ24phys == 52

def class_order(dual):
    out = ZZ(1)
    for x in dual:
        out = lcm(out, ZZ(QQ(x).denominator()))
    return out

def effective_d13_section(G, z, cap=200000):
    from sage.all import IntegralLattice

    R13 = G[:13,:13]
    C13 = G[:13,13:]
    T13 = G[13:,13:]
    H13local = T13 - C13.transpose()*R13.inverse()*C13

    z = vector(ZZ, z)
    h = QQ(z*H13local*z)
    base = vector(ZZ, [0]*13 + list(z))
    pairing = vector(QQ, base*G[:,:13])
    dual = pairing*R13.inverse()
    order = class_order(dual)
    correction = {
        ZZ(1): QQ(0),
        ZZ(2): QQ(1),
        ZZ(4): QQ(13)/4,
    }[order]

    target = h + correction
    assert target in ZZ and target >= 4 and target % 2 == 0
    target = ZZ(target)

    it = IntegralLattice(R13).enumerate_close_vectors(-dual)
    chosen = None
    for unused in range(cap):
        shift = vector(ZZ, next(it))
        cand = base + vector(ZZ, list(shift) + [0]*4)
        norm = ZZ(cand*G*cand)
        if norm == target:
            chosen = cand
            break
        if norm > target:
            break
    assert chosen is not None, (tuple(z), h, correction, target)

    pole = ZZ((target - 4)//2)
    a = ZZ((target - 2)//2)
    section = vector(ZZ, [a,1] + list(chosen))
    return section, h, correction, pole

P24coordsP, h24P, corr24P, pole24P = effective_d13_section(GP, zQ24phys)
P24phys = P24coordsP * BadaptP
assert P24phys*ns*P24phys == -2
assert P24phys*ns*F8phys == 1
assert vector(ZZ, d13_coords_phys(P24phys)[-4:]) == zQ24phys

# Transport every geometric object through the same exact map.
O8 = to_equation(O8phys)
G1curve = to_equation(G1phys_curve)
E77eq = to_equation(E77phys)
P24eq = to_equation(P24phys)

assert O8*ns*F8 == 1
assert G1curve*ns*F8 == 1
assert E77eq*ns*F8 == 2
assert P24eq*ns*F8 == 1
assert P24eq*ns*P24eq == -2

# Natural equation D13 frame around the transported branch zero.
child8, Bzero = child_frame_with_zero(ns, F8, O8)
A13, G13, H13 = d13_root_adaptation(child8)
Badapt = block_diagonal_matrix(identity_matrix(ZZ,2), A13) * Bzero

def d13_coords(P):
    c = vector(QQ, P) * Badapt.inverse()
    assert all(x in ZZ for x in c)
    return vector(ZZ, c)

zG1 = vector(ZZ, d13_coords(G1curve)[-4:])
zAJ77 = vector(ZZ, d13_coords(E77eq)[-4:])
zG3 = zAJ77 - zG1
zQ24 = vector(ZZ, d13_coords(P24eq)[-4:])
zS3 = vector(ZZ, d13_coords(S3)[-4:])

assert zG1*H13*zG1 == QQ(3)/4
assert zG3*H13*zG3 == QQ(11)/4
assert zQ24*H13*zQ24 == 52
assert zAJ77 == zG1 + zG3

q24a = zQ24
q24b = -zQ24

print(
    "Q24EQSEARCH_D13|"
    f"zero_transport=weyl+tau|"
    f"AJ_S3={','.join(map(str,zS3))}|"
    f"G1={','.join(map(str,zG1))}|"
    f"AJ_E77={','.join(map(str,zAJ77))}|"
    f"G3={','.join(map(str,zG3))}|"
    f"q24={','.join(map(str,zQ24))}|"
    f"q24_phys_profile=h{h24P},corr{corr24P},O{pole24P}|"
    "status=PASS_GEOMETRIC_Q24_TRANSPORT",
    flush=True,
)

if args.setup_only:
    print(
        "Q24EQSEARCH_SETUP_RESULT|"
        "q6_words=PASS|equation_D13_marking=PASS|"
        "search=SKIPPED|status=PASS_EXACT_SETUP_ONLY",
        flush=True,
    )
    raise SystemExit(0)

# ---------------------------------------------------------------------------
# 3. Rigorous q8-degree bound in each IV* component class.
# ---------------------------------------------------------------------------

linear=vector(QQ,[phi*ns*F8 for phi in phis])
Hinv=H6.inverse()
center=-QQ(1)/2*linear*Hinv
threshold=QQ(args.threshold)

pattern_quadratics=[]
lower=[None]*3
upper=[None]*3

for pi,p in enumerate(patterns):
    correction=-QQ(p*Groot.inverse()*p)
    root_part=vector(QQ,p)*Groot.inverse()*Rroots
    constant=QQ(O6*ns*F8)+correction+QQ(root_part*ns*F8)
    minimum=constant-QQ(1)/4*(linear*Hinv*linear)
    budget=threshold-minimum
    if budget<=0:
        continue
    bounds=[]
    for i in range(3):
        radius=ceil_sqrt_rational(budget*Hinv[i,i])
        lo=ZZ(center[i].floor())-radius-1
        hi=ZZ(center[i].ceil())+radius+1
        bounds.append((lo,hi))
        lower[i]=lo if lower[i] is None else min(lower[i],lo)
        upper[i]=hi if upper[i] is None else max(upper[i],hi)
    pattern_quadratics.append({
        "index":pi,"pattern":p,"correction":correction,
        "constant":constant,"minimum":minimum,"bounds":bounds,
    })

assert pattern_quadratics
assert all(x is not None for x in lower+upper)

# Closed formula regressions.
for P,w in ((E77,w77),(S3,wS3)):
    PP,pole,corr,p=section_from_word(w)
    assert PP==P
    root_part=vector(QQ,p)*Groot.inverse()*Rroots
    formula=(
        QQ(w*H6*w)+linear*w+QQ(O6*ns*F8)
        +corr+QQ(root_part*ns*F8)
    )
    assert formula==P*ns*F8

print(
    "Q24EQSEARCH_BOUND|"
    f"threshold={args.threshold}|center={','.join(map(str,center))}|"
    f"box={';'.join(f'{lower[i]}..{upper[i]}' for i in range(3))}|"
    f"patterns={len(pattern_quadratics)}|status=PASS_EXHAUSTIVE_BOUND",
    flush=True,
)

# ---------------------------------------------------------------------------
# 4. Exhaustive enumeration.
# ---------------------------------------------------------------------------

records=[]
compatible=[]
tested=0
below=0

for a in range(int(lower[0]),int(upper[0])+1):
    for b in range(int(lower[1]),int(upper[1])+1):
        for c in range(int(lower[2]),int(upper[2])+1):
            if a==b==c==0:
                continue
            n=vector(ZZ,(a,b,c))
            P,pole,corr,p=section_from_word(n)
            tested+=1
            degree=int(P*ns*F8)
            if degree>=args.threshold:
                continue
            below+=1

            root_part=vector(QQ,p)*Groot.inverse()*Rroots
            formula=(
                QQ(n*H6*n)+linear*n+QQ(O6*ns*F8)
                +corr+QQ(root_part*ns*F8)
            )
            assert formula==degree

            z=vector(ZZ,d13_coords(P)[-4:])
            matches=[]
            for sign,target in ((+1,q24a),(-1,q24b)):
                residual=target-z
                ab=solve_two_span(zG1,zG3,residual)
                if ab is not None:
                    matches.append({
                        "orientation":sign,
                        "G1":ab[0],"G3":ab[1],
                    })

            rec={
                "q6_word":list(map(int,n)),
                "q6_height":str(n*H6*n),
                "q6_P_dot_Oold":int(P*ns*O6),
                "q6_P_dot_Ostd":int(P*ns*Ostd),
                "IVstar_correction":str(corr),
                "q8_degree":degree,
                "AJ_equation_D13":list(map(int,z)),
                "q24_matches":matches,
                "source_h3_ns":list(map(int,P)),
            }
            records.append(rec)
            if matches:
                compatible.append(rec)

records.sort(key=lambda r:(r["q8_degree"],sum(abs(x) for x in r["q6_word"])))
compatible.sort(key=lambda r:(r["q8_degree"],sum(abs(x) for x in r["q6_word"])))

print(
    "Q24EQSEARCH|"
    f"tested_box={tested}|below_threshold={below}|"
    f"compatible={len(compatible)}|threshold={args.threshold}|status=PASS",
    flush=True,
)

for rank,rec in enumerate(records[:args.top],1):
    print(
        "Q24EQSEARCH_LOW|"
        f"rank={rank}|word={','.join(map(str,rec['q6_word']))}|"
        f"degree={rec['q8_degree']}|"
        f"AJ={','.join(map(str,rec['AJ_equation_D13']))}|"
        f"corr={rec['IVstar_correction']}|"
        f"compatible={int(bool(rec['q24_matches']))}|status=CANDIDATE",
        flush=True,
    )

for rank,rec in enumerate(compatible[:args.top],1):
    corrections=";".join(
        f"orient={m['orientation']:+d},add={m['G1']}*G1+{m['G3']}*G3"
        for m in rec["q24_matches"]
    )
    print(
        "Q24EQSEARCH_COMPAT|"
        f"rank={rank}|word={','.join(map(str,rec['q6_word']))}|"
        f"degree={rec['q8_degree']}|"
        f"AJ={','.join(map(str,rec['AJ_equation_D13']))}|"
        f"q24={corrections}|status=PASS_COMPATIBLE",
        flush=True,
    )

status=(
    "PASS_FOUND_EQUATION_Q24_BRIDGE"
    if compatible
    else "PASS_NO_EQUATION_Q24_BRIDGE_BELOW_THRESHOLD"
)

payload={
    "schema":"elkies-k3.h92-q24-equation-bridge-search.v1",
    "status":status,
    "threshold":args.threshold,
    "exhaustive_below_threshold":True,
    "equation_D13":{
        "AJ_S3":list(map(int,zS3)),
        "G1":list(map(int,zG1)),
        "G3":list(map(int,zG3)),
        "q24_up_to_inversion":[list(map(int,q24a)),list(map(int,q24b))],
    },
    "search":{
        "tested_box":tested,
        "below_threshold":below,
        "box":[[int(lower[i]),int(upper[i])] for i in range(3)],
        "center":[str(x) for x in center],
    },
    "lowest":records[:args.top],
    "compatible":compatible[:args.top],
    "next":(
        "Construct the cheapest compatible q6 section as an exact rational "
        "point in the standard q6 Weierstrass model, evaluate the repaired q8 "
        "parameter, and perform the direct II*_E8_1 branch-zero AJ trace with "
        "L((d+1)O). Then add the certified G1/G3 correction."
        if compatible else
        "No equation-compatible bridge beats the threshold; use the direct "
        "S3 degree-52 route or raise the search threshold."
    ),
}

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24EQSEARCH_RESULT|"
    f"status={status}|best_degree={compatible[0]['q8_degree'] if compatible else 'none'}|"
    f"threshold={args.threshold}|exhaustive_below_threshold=1",
    flush=True,
)
