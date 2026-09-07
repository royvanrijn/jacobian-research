#!/usr/bin/env sage -python
"""
Replay how the H3 D13 --q24--> D12 neighbour was actually FOUND, then use
that discovery to target the correct equation object: the affine D12
bisection Theta.

Part A (exact ZZ lattice):
  * load the original root-adapted D13 frame;
  * rebuild the q24 fibre D=[12,2;w];
  * show D meets exactly one D13 simple component, a leaf R3, and the old
    affine component, both in degree one;
  * the other 12 old components form D12;
  * compute their D12 highest-root multiplicities;
  * define
        Theta = D - sum m_i C_i;
  * prove
        Theta^2=-2, Theta.F_old=2, Theta.O_old=10, Theta.D=0;
  * prove Theta meets only old component R12 once and misses the old affine
    component.
Thus D = Theta + sum m_i C_i is an explicit D12 fibre member.

Part B (GF(100003), diagnostic):
  * rebuild the already-passing 58 -> 10 smooth-collision space;
  * select the unique member with maximal ordinary cusp vanishing
    (collision + first 9 cusp jets, expected kernel 1);
  * interpret its generic zero divisor as a bisection via the line
        m_P = -a/b
    through -P24;
  * factor the residual quadratic discriminant.
If its squarefree degree is <=2, this is the expected rational Theta
bisection.  If not, the lattice replay still stands and tells the next task:
recover Theta using the actual old-component restrictions rather than a
generic I9* module.

No characteristic-zero q24 point lifting is used.
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    GF, PolynomialRing, QQ, ZZ, block_diagonal_matrix, matrix, pari, vector
)


Q24_WITNESS = vector(ZZ, (
    0,5,0,1,2,1,2,2,2,2,4,8,2,0,-1,1,1,
))


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


def load_gram(path):
    return matrix(ZZ,[
        [ZZ(v) for v in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ])


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"
FRAME=ROOT/"elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"

# ===========================================================================
# A. Exact replay of the D12 discovery.
# ===========================================================================

frame=load_gram(FRAME)
assert frame.dimensions()==(17,17)
U2=matrix(ZZ,((0,1),(1,0)))
ns=block_diagonal_matrix(U2,-frame)
F_old=vector(ZZ,[1,0]+[0]*17)
O_old=vector(ZZ,[-1,1]+[0]*17)
D=vector(ZZ,[12,2]+list(Q24_WITNESS))
assert D*ns*D==0 and D*ns*F_old==2

# Effective D13 simple components in this chamber are -e_i.
components=[]
for i in range(13):
    C=vector(ZZ,[0,0]+[-ZZ(j==i) for j in range(17)])
    assert C*ns*C==-2 and C*ns*F_old==0
    components.append(C)

pairings=[ZZ(D*ns*C) for C in components]
assert pairings==[0,0,1,0,0,0,0,0,0,0,0,0,0]
removed=pairings.index(1)
assert removed==2  # R3

root=frame[:13,:13]
keep=[i for i in range(13) if i!=removed]
d12=root.matrix_from_rows_and_columns(keep,keep)
roots_data=pari(d12).qfminim(2)
root_count=ZZ(roots_data[0])
assert (d12.rank(),root_count,abs(ZZ(d12.det())))==(12,264,4)

# Highest root in the effective D12 simple basis.
half=matrix(ZZ,roots_data[2]).transpose().rows()
roots=[vector(ZZ,r) for r in half]+[-vector(ZZ,r) for r in half]
positive_coordinate=[r for r in roots if all(v>=0 for v in r)]
assert positive_coordinate
highest12=max(positive_coordinate,key=lambda r:sum(r))
m13=[ZZ(0)]*13
for index,value in zip(keep,highest12):
    m13[index]=ZZ(value)

# Old D13 affine component.
d13q=pari(root).qfminim(2)
d13half=matrix(ZZ,d13q[2]).transpose().rows()
d13roots=[vector(ZZ,r) for r in d13half]+[-vector(ZZ,r) for r in d13half]
highest13=max((r for r in d13roots if all(v>=0 for v in r)),key=lambda r:sum(r))
old_affine=F_old+vector(ZZ,[0,0]+list(highest13)+[0]*4)
assert old_affine*ns*old_affine==-2

vertical=sum(
    (m13[i]*components[i] for i in keep),
    vector(ZZ,[0]*19),
)
Theta=D-vertical

assert Theta*ns*Theta==-2
assert Theta*ns*F_old==2
assert Theta*ns*O_old==10
assert Theta*ns*D==0
assert D==Theta+vertical
assert Theta*ns*old_affine==0

theta_simple_hits=[
    (i+1,int(Theta*ns*components[i]))
    for i in range(13)
    if Theta*ns*components[i]
]
assert theta_simple_hits==[(12,1)]

# Check the affine D12 incidence: Theta attaches to the node predicted by the
# D12 highest root.
print(
    "Q24D12DISCOVERY|"
    f"D_simple_pairings={','.join(map(str,pairings))}|"
    f"removed_leaf=R{removed+1}|"
    f"D12_root_data=12,{root_count},{abs(ZZ(d12.det()))}|"
    f"highest={','.join(map(str,m13))}|status=PASS",
    flush=True,
)
print(
    "Q24D12THETA|"
    f"class={','.join(map(str,Theta))}|"
    "square=-2|old_degree=2|O=10|D=0|"
    f"old_component_hits={theta_simple_hits}|old_affine=0|"
    "status=PASS_EXACT_D12_AFFINE_BISECTION",
    flush=True,
)
print(
    "Q24D12FIBRE|"
    "formula=D24=Theta+sum(m_i*C_i)|"
    f"multiplicities={','.join(map(str,m13))}|"
    "status=PASS_EXPLICIT_D12_FIBRE_MEMBER",
    flush=True,
)

# ===========================================================================
# B. Modular candidate for Theta from the existing q24 marked-chord space.
# ===========================================================================

p=ZZ(args.prime)
F=GF(p)
R=PolynomialRing(F,"U")
U=R.gen()
K=R.fraction_field()

MOD=LOCAL/f"q24-degree46-direct-global-mod-{p}.json"
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
for path in (MOD,TRANS,Q8):
    if path is None or not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

mod=json.loads(MOD.read_text())
trans=json.loads(TRANS.read_text())
q8=json.loads(Q8.read_text())
assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert trans["status"]=="PASS_EXACT_Q24_PHYSICAL_TO_EQUATION_TRANSLATION"

def red_q(q):
    q=QQ(q)
    d=ZZ(q.denominator())
    if d%p==0:
        raise ZeroDivisionError("bad reduction denominator")
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
assert Y**2==X**3+A*X*Z**4+B*Z**6
xP=K(X)/K(Z**2)
yP=K(Y)/K(Z**3)

# 58-dimensional collision-saturated global envelope.
ambient=[("A",i) for i in range(42)]+[("B",i) for i in range(16)]
modulus=Z**2
cols=[
    (U**i*X)%modulus if kind=="A" else (-U**i*Y)%modulus
    for kind,i in ambient
]
Cmat=matrix(F,48,58,lambda row,col:cols[col][row])
assert Cmat.rank()==48

# I9* place.
i9=next(item for item in child["finite_fibres"] if item["kodaira"]=="I9*")
RQ=PolynomialRing(QQ,"U")
fQ=RQ(str(i9["factor"]))
f=R([red_q(c) for c in fQ.list()])
alpha=-f[0]/f[1]

S=PolynomialRing(F,"u")
u=S.gen()
KS=S.fraction_field()

def shift_poly(poly):
    return S(poly(alpha+u))

def shift_rf(value):
    value=K(value)
    n=shift_poly(R(value.numerator()))
    d=shift_poly(R(value.denominator()))
    assert d[0]
    return KS(n)/KS(d)

xloc=shift_rf(xP)
yloc=shift_rf(yP)
mcusp=-yloc/xloc
assert mcusp.valuation()>=0

def jet(value,n):
    value=KS(value)
    num=S(value.numerator())
    den=S(value.denominator())
    assert den[0]
    rem=(num*den.inverse_mod(u**n))%(u**n)
    return [rem[i] for i in range(n)]

# First 9 cusp jets on the ORIGINAL 58 columns.  The prior preflight proved
# their image on the 10D collision kernel has rank 9.
Jcols=[]
for kind,i in ambient:
    if kind=="A":
        g=shift_rf(K(U**i)/K(Z**2))
    else:
        g=shift_rf(K(U**i)/K(Z))*mcusp
    Jcols.append(jet(g,9))
J9=matrix(F,9,58,lambda row,col:Jcols[col][row])
M=Cmat.stack(J9)
assert M.rank()==57
theta_kernel=M.right_kernel().basis_matrix()
assert theta_kernel.dimensions()==(1,58)
row=theta_kernel[0]

AA=R.zero()
BB=R.zero()
for coefficient,(kind,i) in zip(row,ambient):
    if kind=="A":
        AA += coefficient*U**i
    else:
        BB += coefficient*U**i
assert BB
assert (AA*X-BB*Y)%modulus==0

a=K(AA)/K(Z**2)
b=K(BB)/K(Z)
slope=-a/b

# m_P=slope means y+yP=slope*(x-xP), a line through -P.
XR=PolynomialRing(K,"x")
xx=XR.gen()
yline=XR(slope)*(xx-XR(xP))-XR(yP)
relation=yline**2-xx**3-XR(K(A))*xx-XR(K(B))
quadratic,remainder=relation.quo_rem(xx-XR(xP))
assert not remainder and quadratic.degree()==2
disc=K(quadratic[1]**2-4*quadratic[2]*quadratic[0])

# Squarefree part of the bisection discriminant over F[U].
num=R(disc.numerator())
den=R(disc.denominator())
sf=R(num.factor().unit()/den.factor().unit())
odd_num=[]
odd_den=[]
for factor,e in num.factor():
    if e%2:
        sf*=factor
        odd_num.append((str(factor),int(e),int(factor.degree())))
for factor,e in den.factor():
    if e%2:
        sf*=factor
        odd_den.append((str(factor),int(e),int(factor.degree())))
sf=sf.monic()
quotient=disc/K(sf)
assert quotient.is_square()
sf_degree=int(sf.degree())

print(
    "Q24D12THETA_MODP|"
    "source=maximal_9_cusp_jet_member|"
    f"Adeg={AA.degree()}|Bdeg={BB.degree()}|"
    f"slope_numdeg={R(slope.numerator()).degree()}|"
    f"slope_dendeg={R(slope.denominator()).degree()}|"
    f"bisection_squarefree_degree={sf_degree}|"
    f"status={'PASS_RATIONAL_BISECTION_CANDIDATE' if sf_degree<=2 else 'NOT_THETA'}",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q24-d12-discovery-replay.v1",
    "status":"PASS_EXACT_D12_DISCOVERY_REPLAY",
    "pinned_D13":{
        "q24_fibre":list(map(int,D)),
        "simple_pairings":list(map(int,pairings)),
        "removed_leaf":removed+1,
        "D12_root_data":[12,int(root_count),int(abs(ZZ(d12.det())))],
        "D12_highest_root_in_D13_numbering":list(map(int,m13)),
        "theta_class":list(map(int,Theta)),
        "theta_square":-2,
        "theta_old_fibre_degree":2,
        "theta_old_zero_intersection":10,
        "theta_new_fibre_intersection":0,
        "theta_old_component_hits":[list(v) for v in theta_simple_hits],
        "theta_old_affine_intersection":0,
        "fibre_decomposition":"D24 = Theta + sum m_i C_i",
    },
    "modular_theta_diagnostic":{
        "prime":int(p),
        "construction":"unique collision+first-nine-cusp-jets member",
        "A_degree":int(AA.degree()),
        "B_degree":int(BB.degree()),
        "slope_numerator_degree":int(R(slope.numerator()).degree()),
        "slope_denominator_degree":int(R(slope.denominator()).degree()),
        "bisection_squarefree_degree":sf_degree,
        "odd_numerator_factors":odd_num,
        "odd_denominator_factors":odd_den,
        "candidate_is_rational":bool(sf_degree<=2),
    },
    "next":(
        "If the modular maximal-vanishing member is rational, transport and "
        "identify it with the exact Theta class, then use its defining equation "
        "as the distinguished D12 fibre generator. Otherwise derive Theta "
        "directly from restrictions to the twelve known old D13 components; "
        "do not return to generic q24-section lifting."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q24-d12-discovery-replay-mod-{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24D12DISCOVERY_RESULT|"
    f"theta_old_degree=2|theta_O=10|"
    f"modp_bisection_sf_degree={sf_degree}|"
    "status=PASS_EXACT_D12_DISCOVERY_REPLAY",
    flush=True,
)
