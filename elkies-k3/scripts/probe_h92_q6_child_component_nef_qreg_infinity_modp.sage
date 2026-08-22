#!/usr/bin/env sage -python
"""Probe the corrected component-nef q_regular finite module at infinity.

The exact finite module is

    g1=q_regular-rho,   g2=M=f_II^2*f_IV^2.

In standard coordinates
    q_regular = alpha + beta*m_inf,
with alpha=-p/h-R/Nx and beta=T^2/h.

The component-nef pencil is base-isomorphic to its pullback by tau_-P0, so
branch degrees can be computed in this standard frame.  We test the
determinant-compatible smooth-infinity scalar order (20,20) on the complete
degree bounds deg(s)<=24, deg(t)<=23 and report the exact modular kernel.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
CHILD = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-jacobian.json"
MARKING = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json"
NORMALIZER = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json"
DEFAULT_MODULE = ROOT / "artifacts/local/elkies-k3/q8-component-nef-qreg-finite-module.json"
DEFAULT_OUTPUT = ROOT / "artifacts/local/elkies-k3/q8-component-nef-qreg-infinity-modp.json"


def coefficient(finite, value):
    value = QQ(value)
    den = finite(ZZ(value.denominator()))
    if not den:
        raise ValueError("prime divides an input denominator")
    return finite(ZZ(value.numerator()))/den


def polynomial(ring, finite, coefficients):
    return ring([coefficient(finite, value) for value in coefficients])


def monic_power_root(value, exponent):
    root=value.parent().one()
    for factor,multiplicity in value.factor():
        assert multiplicity % exponent==0
        root *= factor.monic()**(multiplicity//exponent)
    return root.monic()


def degree_or_minus_one(value):
    return -1 if not value else value.degree()


parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime",type=int,default=43)
parser.add_argument("--module",type=Path,default=DEFAULT_MODULE)
parser.add_argument("--screen-levels",action="store_true")
parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
args=parser.parse_args()
if not ZZ(args.prime).is_prime() or args.prime in (2,3):
    raise ValueError("prime must be odd and !=3")
args.module=args.module.resolve()
args.output=args.output.resolve()

child=json.loads(CHILD.read_text())
marking=json.loads(MARKING.read_text())
normalizer=json.loads(NORMALIZER.read_text())
module=json.loads(args.module.read_text())
assert module["status"]=="PASS_EXACT_COMPONENT_NEF_QREG_FINITE_MODULE"

finite=GF(args.prime)
ring=PolynomialRing(finite,"T")
T=ring.gen()
field=ring.fraction_field()

sdata=marking["selected_q8"]["relative_child_section_standard_jacobian_coordinates"]
nx=polynomial(ring,finite,sdata["x_numerator_coefficients_low_to_high"])
dx=polynomial(ring,finite,sdata["x_denominator_coefficients_low_to_high"])
ny=polynomial(ring,finite,sdata["y_numerator_coefficients_low_to_high"])
dy=polynomial(ring,finite,sdata["y_denominator_coefficients_low_to_high"])
sx,sy=field(nx)/field(dx),field(ny)/field(dy)
h=monic_power_root(dx,2)
R=polynomial(ring,finite,normalizer["normalizer"]["R_coefficients_low_to_high"])
assert (R*h*dy-ny)%nx==0
p=-sy/sx

rho=polynomial(ring,finite,module["module"]["rho_coefficients_low_to_high"])
ii=polynomial(ring,finite,PolynomialRing(QQ,"T")(
    next(item for item in child["finite_fibres"] if item["kodaira"]=="II*")["factor"]
).list()).monic()
iv=polynomial(ring,finite,PolynomialRing(QQ,"T")(
    next(item for item in child["finite_fibres"] if item["kodaira"]=="IV*")["factor"]
).list()).monic()
M=ii**2*iv**2
assert M.degree()==4 and rho.degree()<4

alpha=-p/field(h)-field(R)/field(nx)
beta=field(T**2)/field(h)
alpha_num,alpha_den=ring(alpha.numerator()),ring(alpha.denominator())
beta_num,beta_den=ring(beta.numerator()),ring(beta.denominator())
assert alpha_den(0) and beta_den(0)

max_s,max_t=24,23
labels=[("s",i) for i in range(max_s+1)]+[("t",i) for i in range(max_t+1)]
pairs=[]
for kind,e in labels:
    if kind=="s":
        pairs.append((T**e,-rho*T**e))
    else:
        pairs.append((ring.zero(),M*T**e))

a_values=[C*alpha_den+B*alpha_num for B,C in pairs]
b_values=[B*beta_num for B,C in pairs]

def rows_for(values,den_degree,required_order):
    cutoff=den_degree-required_order
    top=max((degree_or_minus_one(v) for v in values),default=-1)
    return [
        [v[d] if d<=v.degree() else finite.zero() for v in values]
        for d in range(cutoff+1,top+1)
    ]

a_rows=rows_for(a_values,alpha_den.degree(),20)
b_rows=rows_for(b_values,beta_den.degree(),20)
condition=matrix(finite,a_rows+b_rows,ncols=len(labels))
kernel=condition.right_kernel_matrix()

def pair_from_row(row):
    s_poly=sum(
        row[i]*T**e for i,(kind,e) in enumerate(labels) if kind=="s"
    )
    t_poly=sum(
        row[i]*T**e for i,(kind,e) in enumerate(labels) if kind=="t"
    )
    B=s_poly
    C=-rho*s_poly+M*t_poly
    a=field(C)-field(B)*p/field(h)-field(B*R)/field(nx)
    b=field(B)/field(h)
    return a,b

branch=None
if args.screen_levels and kernel.nrows()==2:
    Acurve=polynomial(ring,finite,child["minimal_short_weierstrass"]["A_coefficients_low_to_high"])
    Bcurve=polynomial(ring,finite,child["minimal_short_weierstrass"]["B_coefficients_low_to_high"])
    a0,b0=pair_from_row(kernel.row(0))
    a1,b1=pair_from_row(kernel.row(1))
    xring=PolynomialRing(field,"x")
    x=xring.gen()
    hist={}
    good=[]
    singular=[]
    for level in finite:
        den=b1-field(level)*b0
        if not den:
            singular.append(int(level)); continue
        mvalue=-(a1-field(level)*a0)/den
        y=xring(mvalue)*(x-xring(sx))-xring(sy)
        relation=y**2-x**3-xring(Acurve)*x-xring(Bcurve)
        quad,rem=relation.quo_rem(x-xring(sx))
        assert not rem and quad.degree()==2
        disc=xring.base_ring()(quad[1]**2-4*quad[2]*quad[0])
        num,denp=ring(disc.numerator()),ring(disc.denominator())
        odd=sum(
            factor.degree()
            for value in (num,denp)
            for factor,mult in value.squarefree_decomposition()
            if mult%2
        )
        degree=int(odd+(denp.degree()-num.degree())%2)
        hist[degree]=hist.get(degree,0)+1
        if degree==4: good.append(int(level))
    branch={
        "histogram":{str(k):v for k,v in sorted(hist.items())},
        "genus_one_levels":good,
        "singular_levels":singular,
    }

payload={
    "schema":"elkies-k3.h92-q6-child-q8-component-nef-qreg-infinity-modp.v1",
    "status":"PASS_COMPONENT_NEF_QREG_INFINITY_PROBE",
    "prime":args.prime,
    "finite_rho_degree":int(rho.degree()),
    "ambient_dimension":len(labels),
    "a_rows":len(a_rows),
    "b_rows":len(b_rows),
    "rank":int(condition.rank()),
    "kernel_dimension":int(kernel.nrows()),
    "kernel_basis":[[int(v) for v in row] for row in kernel.rows()],
    "branch_screen":branch,
    "boundary":"Modular infinity intersection; characteristic-zero certification is separate.",
}
args.output.parent.mkdir(parents=True,exist_ok=True)
args.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(
    "Q8QREGINF|prime={}|rho_deg={}|ambient={}|rank={}|kernel={}|hist={}|good={}|"
    "status=PASS_COMPONENT_NEF_QREG_INFINITY_PROBE".format(
        args.prime,rho.degree(),len(labels),condition.rank(),kernel.nrows(),
        "none" if branch is None else ",".join("{}:{}".format(k,v) for k,v in branch["histogram"].items()),
        "none" if branch is None or not branch["genus_one_levels"] else ",".join(map(str,branch["genus_one_levels"])),
    ),
    flush=True,
)
