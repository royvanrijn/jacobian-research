#!/usr/bin/env sage -python
"""
Audit actual repaired-q8 map degrees of exact q6 rational sections relevant to
the four equation-search q24 candidates.

This is ground truth: degrees are computed from the actual QQ(T) q8 parameter,
not inferred from an NS/component presentation.

In particular:
  rank 2 old word (-2,-1,1) -> std (0,-2,1) = Qmap + S3
  rank 4 old word (-4, 1,-1) -> std (-2,0,-1) = Pmap - S3

If rank 4 has actual degree 46, it is preferable because its q24 correction
uses only G1 and all q6 ingredients are already exact.
"""

import json
from pathlib import Path
from sage.all import EllipticCurve, PolynomialRing, QQ


def locate_repo():
    h=Path.home()
    for c in [
        Path.cwd().resolve(),
        h/"Documents"/"jacobian-research",
        h/"jacobian-research",
    ]:
        if (c/"elkies-k3/scripts").is_dir():
            return c
    raise SystemExit("repo not found")


def poly(R,vals):
    return R([QQ(v) for v in vals])


def rat(K,R,data,nk,dk):
    return K(poly(R,data[nk]))/K(poly(R,data[dk]))


def monic_power_root(value,e):
    out=value.parent().one()
    for f,m in value.factor():
        assert m%e==0
        out*=f.monic()**(m//e)
    return out.monic()


ROOT=locate_repo()
GEN=ROOT/"artifacts/generated-results"
LOCAL=ROOT/"artifacts/local/elkies-k3"

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
    raise SystemExit("complete q8 child missing")

q6=json.loads(Q6.read_text())
zero=json.loads(ZERO.read_text())
comp=json.loads(COMP.read_text())
s3br=json.loads(S3BR.read_text())
q8=json.loads(Q8.read_text())

R=PolynomialRing(QQ,"T"); T=R.gen(); K=R.fraction_field()
m=q6["minimal_short_weierstrass"]
A=poly(R,m["A_coefficients_low_to_high"])
B=poly(R,m["B_coefficients_low_to_high"])
E=EllipticCurve(K,[0,0,0,K(A),K(B)])

zd=zero["section"]
Oold=E(
    rat(K,R,zd,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high"),
    rat(K,R,zd,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high"),
)

entries={e["sign"]:e for e in comp["sections"]}
pts={
    sign:E(
        rat(K,R,e,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high"),
        rat(K,R,e,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high"),
    )
    for sign,e in entries.items()
}
affine=pts[comp["source"]["affine_E7_sign"]]
e77=pts[comp["source"]["E7_7_sign"]]
Pmap=e77-Oold
Qmap=e77-affine

sd=s3br["third_section_canonical_q6"]
S3=E(
    rat(K,R,sd["x"],"numerator_coefficients_low_to_high","denominator_coefficients_low_to_high"),
    rat(K,R,sd["y"],"numerator_coefficients_low_to_high","denominator_coefficients_low_to_high"),
)

md=q8["marking"]["section"]
sx=rat(K,R,md,"x_numerator_coefficients_low_to_high","x_denominator_coefficients_low_to_high")
sy=rat(K,R,md,"y_numerator_coefficients_low_to_high","y_denominator_coefficients_low_to_high")
Smark=E(sx,sy)
assert Smark==Pmap+Qmap

nx,dx=R(sx.numerator()),R(sx.denominator())
ny,dy=R(sy.numerator()),R(sy.denominator())
h=monic_power_root(dx,2)
assert h==monic_power_root(dy,3)

ii=R(next(x for x in q6["finite_fibres"] if x["kodaira"]=="II*")["factor"]).monic()
iv=R(next(x for x in q6["finite_fibres"] if x["kodaira"]=="IV*")["factor"]).monic()
M=(ii**2*iv**2).monic()

normalizer=(ny*dx*(h*dy).inverse_mod(nx)).mod(nx)
pfun=-sy/sx
rho=(normalizer*nx.inverse_mod(M)).mod(M)

pairs=[]
for entry in q8["rr"]["kernel_polynomials"]:
    sp=R(entry["s"]); tp=R(entry["t"])
    Bc=K(sp)/K(h)
    Ac=(
        -K(sp)*pfun/K(h)
        -K(sp)*K(normalizer)/K(nx)
        +K(sp*rho)+K(tp*M)
    )
    pairs.append((Ac,Bc))
(A0,B0),(A1,B1)=pairs


def map_degree(P):
    if P.is_zero():
        return None
    x,y=P.xy()
    mm=(y+sy)/(x-sx)
    U=K((A1+B1*mm)/(A0+B0*mm))
    n,d=R(U.numerator()),R(U.denominator())
    assert n.gcd(d) in QQ
    return max(n.degree(),d.degree()),n.degree(),d.degree()


tests=[
    ("old_E7_7",e77,None,None),
    ("old_E7_affine",affine,None,None),
    ("S3",S3,(-2,1,1),None),
    ("rank2_Qmap_plus_S3",Qmap+S3,(-2,-1,1),(0,-2,1)),
    ("rank4_Pmap_minus_S3",Pmap-S3,(-4,1,-1),(-2,0,-1)),
    ("Pmap_plus_S3",Pmap+S3,(-4,1,1),(-2,0,1)),
    ("Qmap_minus_S3",Qmap-S3,(-2,-1,-1),(0,-2,-1)),
]

records=[]
for name,P,oldmw,stdmw in tests:
    degree,numdeg,dendeg=map_degree(P)
    x,y=P.xy()
    rec={
        "name":name,
        "old_mw":oldmw,
        "standard_mw":stdmw,
        "degree":int(degree),
        "numdeg":int(numdeg),
        "dendeg":int(dendeg),
        "x_degrees":[int(R(x.numerator()).degree()),int(R(x.denominator()).degree())],
        "y_degrees":[int(R(y.numerator()).degree()),int(R(y.denominator()).degree())],
    }
    records.append(rec)
    print(
        "Q24ACTUALDEG|"
        f"name={name}|old_mw={oldmw or 'NA'}|std_mw={stdmw or 'NA'}|"
        f"q8={numdeg}/{dendeg}|degree={degree}|"
        f"x={rec['x_degrees'][0]}/{rec['x_degrees'][1]}|"
        f"y={rec['y_degrees'][0]}/{rec['y_degrees'][1]}|status=PASS",
        flush=True,
    )

rank2=next(r for r in records if r["name"]=="rank2_Qmap_plus_S3")
rank4=next(r for r in records if r["name"]=="rank4_Pmap_minus_S3")

print(
    "Q24ACTUALDEG_RESULT|"
    f"rank2={rank2['degree']}|rank4={rank4['degree']}|"
    f"best_accessible={'rank4' if rank4['degree']<rank2['degree'] else 'rank2'}|"
    "status=PASS_ACTUAL_Q8_DEGREES",
    flush=True,
)

out=LOCAL/"q24-equation-candidate-actual-degrees.json"
out.write_text(json.dumps({
    "schema":"elkies-k3.h92-q24-equation-candidate-actual-degrees.v1",
    "status":"PASS_ACTUAL_Q8_DEGREES",
    "records":records,
},indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{out}",flush=True)
