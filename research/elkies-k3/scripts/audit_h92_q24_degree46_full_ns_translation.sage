#!/usr/bin/env sage -python
"""
Audit and construct the exact equation-frame degree-46 q24 bridge by applying
the already-certified q6 Eichler translation to the FULL physical NS curve.

No MW-word identification is assumed.

Inputs:
  q24-cheapest-bridge-current.json
      selected.ambient_ns_vector = exact physical raw-q6 effective section P_phys
  q8-target-component-nef.json
      physical q8 fibre F8_phys
  certify_h92_q8_equation_ns_divisor.sage
      exact equation fibre F8_eq, q6 fibre F6, old zero, standard zero section,
      and exact q6 Shioda basis/functions.

We certify:
  tau(F8_phys) = F8_eq
  P_phys.F8_phys = 46
  P_eq=tau(P_phys)
  P_eq.F8_eq = 46

Then recover P_eq's old-zero q6 MW coordinates directly from the Shioda map,
convert to standard-Weierstrass coordinates, and build the same generic q6
point from exact known points (old_zero, E7_7, affine_E7, S3) whenever the
coordinates lie in their integral span.

Finally evaluate the repaired q8 pencil and require degree 46.
"""

import argparse
import json
import sys
from pathlib import Path

from sage.all import (
    EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector
)


def locate_repo(explicit=None):
    candidates = []
    if explicit:
        candidates.append(Path(explicit).expanduser())
    cwd = Path.cwd().resolve()
    candidates += [cwd, *cwd.parents]
    h = Path.home()
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


def run_scope(path):
    saved=list(sys.argv)
    scope={"__name__":"__embedded__"}
    try:
        sys.argv=[str(path)]
        exec(compile(path.read_text(),str(path),"exec"),scope)
    finally:
        sys.argv=saved
    return scope


def polynomial(R, vals):
    return R([QQ(v) for v in vals])


def rational(K,R,data,nk,dk):
    return K(polynomial(R,data[nk]))/K(polynomial(R,data[dk]))


def monic_power_root(value, exponent):
    out=value.parent().one()
    for f,m in value.factor():
        assert m % exponent == 0
        out *= f.monic()**(m//exponent)
    return out.monic()


def rf_record(v,R):
    return {
        "numerator_coefficients_low_to_high":[str(x) for x in R(v.numerator()).list()],
        "denominator_coefficients_low_to_high":[str(x) for x in R(v.denominator()).list()],
    }


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--repo",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
GEN=ROOT/"artifacts/generated-results"
LOCAL=ROOT/"artifacts/local/elkies-k3"
CERT=ROOT/"elkies-k3/scripts/certify_h92_q8_equation_ns_divisor.sage"
CHEAP=LOCAL/"q24-cheapest-bridge-current.json"
TARGET=LOCAL/"q8-target-component-nef.json"
Q6=GEN/"elkies-k3-h92-q6-child-jacobian.json"
ZERO=GEN/"elkies-k3-h92-q6-child-zero-section.json"
COMP=GEN/"elkies-k3-h92-q6-child-e7-infinity-sections.json"
S3BR=LOCAL/"q6-third-to-q8-bridge.json"

q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((
    p for p in q8_candidates
    if p.exists()
    and "rr" in json.loads(p.read_text())
    and "kernel_polynomials" in json.loads(p.read_text()).get("rr",{})
),None)
if Q8 is None:
    raise SystemExit("No complete q8 child artifact")

for p in (CERT,CHEAP,TARGET,Q6,ZERO,COMP,S3BR,Q8):
    if not p.exists():
        raise SystemExit(f"missing {p}")

OUT=args.output.resolve() if args.output else LOCAL/"q24-degree46-equation-bridge-exact.json"

eq=run_scope(CERT)
ns=eq["ns"]
F6=vector(ZZ,eq["F6"])
Oold=vector(ZZ,eq["Oold"])
F8eq=vector(ZZ,eq["F8eq"])
Ostd=vector(ZZ,eq["Ostd"])
phis=[vector(QQ,p) for p in eq["phis"]]
H=matrix(QQ,eq["EXPECTED_HEIGHT"])

cheap=json.loads(CHEAP.read_text())
target=json.loads(TARGET.read_text())
q6=json.loads(Q6.read_text())
zero=json.loads(ZERO.read_text())
components=json.loads(COMP.read_text())
s3br=json.loads(S3BR.read_text())
q8=json.loads(Q8.read_text())

sel=cheap["selected"]
assert sel["kind"]=="raw_q6"
assert sel["q8_degree"]==46
assert sel["ambient_ns_vector"] is not None

Pphys=vector(ZZ,sel["ambient_ns_vector"])
F8phys=vector(ZZ,target["selected_q8"]["source_h3_ns_vector"])

assert Pphys*ns*Pphys==-2
assert Pphys*ns*F6==1
assert Pphys*ns*F8phys==46

# ---------------------------------------------------------------------------
# Exact q6 Eichler transvection old-zero -> standard-zero.
# ---------------------------------------------------------------------------

ov=ZZ(Ostd*ns*Oold)
v=Ostd-Oold-(ov+2)*F6
assert v*ns*F6==0
assert v*ns*Oold==0
assert v*ns*v==-12

def tau(x):
    x=vector(ZZ,x)
    return (
        x
        + (x*ns*F6)*v
        - (x*ns*v)*F6
        - (v*ns*v//2)*(x*ns*F6)*F6
    )

# Build matrix and certify integral isometry globally.
M=matrix(ZZ,[list(tau(vector(ZZ,[ZZ(i==j) for i in range(19)]))) for j in range(19)])
assert abs(M.det())==1
assert M*ns*M.transpose()==ns
assert tau(F6)==F6
assert tau(Oold)==Ostd

mapped_fibre=tau(F8phys)
assert mapped_fibre==F8eq

Peq=tau(Pphys)
assert Peq*ns*Peq==-2
assert Peq*ns*F6==1
assert Peq*ns*F8eq==46

print(
    "Q24WNS_TRANS|"
    f"physical_degree={Pphys*ns*F8phys}|"
    f"tau_fibre={int(mapped_fibre==F8eq)}|"
    f"equation_degree={Peq*ns*F8eq}|"
    f"det={M.det()}|status=PASS_EXACT_NS_TRANSLATION",
    flush=True,
)

# ---------------------------------------------------------------------------
# Read MW coordinates DIRECTLY from Shioda classes.
# ---------------------------------------------------------------------------

roots_source=eq["roots_source"]
root_gram=eq["root_gram"]
projection=eq["projection"]

def shioda(P,zero):
    horizontal=P-zero-(P*ns*zero+2)*F6
    assert horizontal*ns*F6==0 and horizontal*ns*zero==0
    return vector(QQ,horizontal)*projection

# old-group coordinates of any q6 section
phi_peq=shioda(Peq,Oold)
pair=vector(QQ,[-phi_peq*ns*q for q in phis])
mw_old_q=pair*H.inverse()
assert all(x in ZZ for x in mw_old_q)
mw_old=vector(ZZ,mw_old_q)

# Old-group coordinate of Ostd is known exactly from certifier.
mw_Ostd_old=vector(ZZ,(-2,1,0))
# standard coordinate subtracts the standard-zero old-group coordinate
mw_std=mw_old-mw_Ostd_old

# Check by rebuilding section class from old-group MW through certifier when
# identity-component reconstruction supports this vector.
rebuilt_same=None
try:
    Pcheck,unused_pole,unused_h=eq["section_from_old_mw"](mw_old)
    rebuilt_same=bool(vector(ZZ,Pcheck)==Peq)
except Exception:
    rebuilt_same=None

print(
    "Q24WNS_MW|"
    f"old_mw={','.join(map(str,mw_old))}|"
    f"standard_mw={','.join(map(str,mw_std))}|"
    f"Peq_Oold={Peq*ns*Oold}|Peq_Ostd={Peq*ns*Ostd}|"
    f"rebuild_same={'NA' if rebuilt_same is None else int(rebuilt_same)}|"
    "status=PASS_DIRECT_SHIODA_COORDINATES",
    flush=True,
)

# ---------------------------------------------------------------------------
# Exact rational q6 point from geometric points with certified std coordinates.
# We solve against old_zero=(2,-1,0), E7_7=(0,-1,0),
# affine=(0,1,0), S3=(0,0,1).
# This spans exactly the sublattice with even first coordinate.
# ---------------------------------------------------------------------------

R=PolynomialRing(QQ,"T"); T=R.gen(); K=R.fraction_field()
model=q6["minimal_short_weierstrass"]
A=polynomial(R,model["A_coefficients_low_to_high"])
B=polynomial(R,model["B_coefficients_low_to_high"])
E=EllipticCurve(K,[0,0,0,K(A),K(B)])

zd=zero["section"]
Pold=E(
    rational(K,R,zd,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high"),
    rational(K,R,zd,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high"),
)
points={}
for item in components["sections"]:
    points[item["sign"]]=E(
        rational(K,R,item,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high"),
        rational(K,R,item,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high"),
    )
Paff=points[components["source"]["affine_E7_sign"]]
Pe77=points[components["source"]["E7_7_sign"]]

s3d=s3br["third_section_canonical_q6"]
PS3=E(
    rational(K,R,s3d["x"],"numerator_coefficients_low_to_high","denominator_coefficients_low_to_high"),
    rational(K,R,s3d["y"],"numerator_coefficients_low_to_high","denominator_coefficients_low_to_high"),
)

a,b,c=map(ZZ,mw_std)
if a%2:
    print(
        "Q24WNS_RATIONAL|"
        f"standard_mw={','.join(map(str,mw_std))}|"
        "status=NEEDS_PRIMITIVE_FIRST_Q6_GENERATOR",
        flush=True,
    )
    rational_built=False
    Pgen=None
else:
    # Solve a,b using Pold=(2,-1), e77=(0,-1):
    # alpha*Pold + beta*e77 = (a,b)
    alpha=a//2
    beta=-b-alpha
    Pgen=alpha*Pold + beta*Pe77 + c*PS3
    rational_built=True
    assert Pgen in E and not Pgen.is_zero()

    # Independent coordinate regression: affine=-e77 in the second primitive direction.
    assert Pe77 + Paff == E(0)

    gx,gy=Pgen.xy()
    print(
        "Q24WNS_RATIONAL|"
        f"standard_mw={','.join(map(str,mw_std))}|"
        f"formula={alpha}*old_zero+{beta}*E7_7+{c}*S3|"
        f"x={R(gx.numerator()).degree()}/{R(gx.denominator()).degree()}|"
        f"y={R(gy.numerator()).degree()}/{R(gy.denominator()).degree()}|"
        "status=PASS_EXACT_Q6_RATIONAL_POINT",
        flush=True,
    )

# ---------------------------------------------------------------------------
# Repaired q8 parameter; if rational point built, degree MUST equal NS 46.
# ---------------------------------------------------------------------------

q8_degree=None
if rational_built:
    md=q8["marking"]["section"]
    sx=rational(K,R,md,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high")
    sy=rational(K,R,md,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high")
    Smark=E(sx,sy)

    nx,dx=R(sx.numerator()),R(sx.denominator())
    ny,dy=R(sy.numerator()),R(sy.denominator())
    h=monic_power_root(dx,2)
    assert h==monic_power_root(dy,3)
    ii=R(next(x for x in q6["finite_fibres"] if x["kodaira"]=="II*")["factor"]).monic()
    iv=R(next(x for x in q6["finite_fibres"] if x["kodaira"]=="IV*")["factor"]).monic()
    MM=(ii**2*iv**2).monic()
    normalizer=(ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
    pfun=-sy/sx
    rho=(normalizer*nx.inverse_mod(MM)).mod(MM)

    pairs=[]
    for entry in q8["rr"]["kernel_polynomials"]:
        sp=R(entry["s"]); tp=R(entry["t"])
        Bc=K(sp)/K(h)
        Ac=(
            -K(sp)*pfun/K(h)
            -K(sp)*K(normalizer)/K(nx)
            +K(sp*rho)+K(tp*MM)
        )
        pairs.append((Ac,Bc))
    (A0,B0),(A1,B1)=pairs

    gx,gy=Pgen.xy()
    mg=(gy+sy)/(gx-sx)
    Ug=K((A1+B1*mg)/(A0+B0*mg))
    un,ud=R(Ug.numerator()),R(Ug.denominator())
    assert un.gcd(ud) in QQ
    q8_degree=max(un.degree(),ud.degree())

    print(
        "Q24WNS_Q8|"
        f"ns_degree={Peq*ns*F8eq}|function_degree={q8_degree}|"
        f"numdeg={un.degree()}|dendeg={ud.degree()}|"
        f"match={int(q8_degree==46)}|status=PASS_DIAGNOSTIC",
        flush=True,
    )
    assert q8_degree==46

payload={
    "schema":"elkies-k3.h92-q24-degree46-equation-bridge-exact.v1",
    "status":(
        "PASS_EXACT_EQUATION_DEGREE46_BRIDGE"
        if rational_built and q8_degree==46
        else "PASS_NS_DEGREE46_NEEDS_PRIMITIVE_Q6_POINT"
    ),
    "physical":{
        "q6_word":sel["word"],
        "ns_vector":list(map(int,Pphys)),
        "q8_degree":46,
        "d13_mw":sel["d13_mw"],
        "q24_correction":{"G1":3,"G3":-1},
    },
    "translation":{
        "eichler_vector":list(map(int,v)),
        "matrix_rows":[list(map(int,row)) for row in M.rows()],
        "physical_fibre_to_equation":True,
        "equation_ns_vector":list(map(int,Peq)),
        "equation_q8_degree":int(Peq*ns*F8eq),
    },
    "q6_mw":{
        "old_zero_coordinates":list(map(int,mw_old)),
        "standard_coordinates":list(map(int,mw_std)),
        "rebuild_same":rebuilt_same,
    },
}
if rational_built:
    gx,gy=Pgen.xy()
    payload["rational_point"]={
        "x":rf_record(gx,R),"y":rf_record(gy,R),
        "q8_degree":int(q8_degree),
    }

OUT.parent.mkdir(parents=True,exist_ok=True)
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24WNS_RESULT|"
    f"old_mw={','.join(map(str,mw_old))}|"
    f"standard_mw={','.join(map(str,mw_std))}|"
    f"equation_degree={Peq*ns*F8eq}|"
    f"rational={int(rational_built)}|"
    f"status={payload['status']}",
    flush=True,
)
