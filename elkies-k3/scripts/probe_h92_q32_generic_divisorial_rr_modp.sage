#!/usr/bin/env sage -python
"""
Intrinsic generic-component valuation test for the clean q32 D12 neighbour.

This replaces the failed maximal-ideal jet test.  For each resolved I9*
exceptional component C:

  1. choose a blow-up chart containing a generic point of C;
  2. parameterize C by a function-field parameter t;
  3. use the STRICT TRANSFORM surface equation to solve formally in a
     transverse parameter s;
  4. pull the q32 chord numerator
         G=A*(Z^2*x-X)+B*(Z^3*y+Y)
     to that formal generic point;
  5. impose v_C(G) >= required_order.

Thus we impose divisorial valuations in the surface local ring, not total
maximal-ideal order at an infinitely-near center.

Expected:
    ambient 56 -> collision rank 48 -> post-collision 8
    resolved divisorial rank 6 -> H0 dimension 2.
"""

import argparse
import json
from pathlib import Path

from sage.all import (
    FunctionField, GF, PolynomialRing, PowerSeriesRing,
    QQ, ZZ, matrix, sage_eval
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

MOD=LOCAL/f"q24-degree46-direct-global-mod-{args.prime}.json"
RES=LOCAL/f"q24-i9star-resolution-mod-{args.prime}.json"
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

for path in (MOD,RES,Q8):
    if not path.exists():
        raise SystemExit(f"Missing prerequisite: {path}")

mod=json.loads(MOD.read_text())
res=json.loads(RES.read_text())
q8=json.loads(Q8.read_text())
assert mod["status"]=="PASS_MODULAR_Q24_FROM_DIRECT_DEGREE46_BRIDGE"
assert res["status"]=="PASS_EXPLICIT_MODP_I9STAR_D13_EXCEPTIONAL_COMPONENTS"
assert q8["status"]=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"

p=ZZ(args.prime)
F=GF(p)
RU=PolynomialRing(F,"U")
U=RU.gen()
KU=RU.fraction_field()

def red_q(q):
    q=QQ(q)
    d=ZZ(q.denominator())
    if d%p==0:
        raise ZeroDivisionError(f"denominator divisible by {p}")
    return F(ZZ(q.numerator()))/F(d)

def red_poly(values):
    return RU([red_q(QQ(v)) for v in values])

child=q8["child"]
Acurve=red_poly(child["minimal_A_coefficients_low_to_high"])
Bcurve=red_poly(child["minimal_B_coefficients_low_to_high"])

sec=mod["section_mod_p"]
Z=RU([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
X=RU([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
Y=RU([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
assert (Z.degree(),X.degree(),Y.degree())==(24,52,78)
assert Y**2==X**3+Acurve*X*Z**4+Bcurve*Z**6

# I9* base.
i9=next(item for item in child["finite_fibres"] if item["kodaira"]=="I9*")
RQ=PolynomialRing(QQ,"U")
fQ=RQ(str(i9["factor"]))
f=RU([red_q(c) for c in fQ.list()])
assert f.degree()==1
alpha=-f[0]/f[1]
assert Z(alpha)!=0

# ---------------------------------------------------------------------------
# Global q32 ambient and smooth-collision quotient.
# ---------------------------------------------------------------------------

Amax=40
Bmax=14
ambient=[("A",i) for i in range(Amax+1)] + [("B",i) for i in range(Bmax+1)]
assert len(ambient)==56

modulus=Z**2
cols=[]
for kind,i in ambient:
    cols.append(
        (U**i*X)%modulus if kind=="A"
        else (-U**i*Y)%modulus
    )
C=matrix(F,48,56,lambda row,col:cols[col][row])
assert C.rank()==48
K8=C.right_kernel().basis_matrix()
assert K8.dimensions()==(8,56)

def pair_from_row(row):
    AA=RU.zero()
    BB=RU.zero()
    for j,(kind,i) in enumerate(ambient):
        if kind=="A":
            AA += row[j]*U**i
        else:
            BB += row[j]*U**i
    return AA,BB

post_pairs=[pair_from_row(row) for row in K8.rows()]

print(
    "Q32DIVVAL_GLOBAL|ambient=56|collision_rank=48|post_collision=8|"
    "status=PASS",
    flush=True,
)

# ---------------------------------------------------------------------------
# Formal generic-point machinery on resolved charts.
# ---------------------------------------------------------------------------

S=PolynomialRing(F,names=("u","x","y"),order="degrevlex")
u,x,y=S.gens()
sloc={"u":u,"x":x,"y":y}

Kr=FunctionField(F,"t")
t=Kr.gen()
PREC=16
PS=PowerSeriesRing(Kr,"s",default_prec=PREC)
ss=PS.gen()

def parse_s(text):
    return S(sage_eval(str(text),locals=sloc))

def chart_data(rec,kind):
    return next(item for item in rec["charts"] if item["chart"]==kind)

def factor_chart(factor_text,kind):
    q=parse_s(factor_text)
    if kind=="u":
        return S(q(1,x,y))
    if kind=="x":
        return S(q(u,1,y))
    if kind=="y":
        return S(q(u,x,1))
    raise ValueError(kind)

def univar_at_series(poly, ub):
    poly=RU(poly)
    out=PS.zero()
    pw=PS.one()
    base=PS(F(alpha))+ub
    # Horner is simpler and avoids high temporary powers.
    out=PS(F(0))
    for c in reversed(poly.list()):
        out=out*base+PS(F(c))
    return out

def newton_solve(strict, values_template, solve_index, root0):
    root=PS(Kr(root0))
    deriv=strict.derivative((u,x,y)[solve_index])
    for _ in range(7):
        vals=list(values_template)
        vals[solve_index]=root
        hv=PS(strict(*vals))
        dv=PS(deriv(*vals))
        if not dv or dv.valuation()!=0:
            raise ArithmeticError(
                ("nonunit implicit derivative",solve_index,root0,dv)
            )
        root=root-hv/dv
    vals=list(values_template)
    vals[solve_index]=root
    residual=PS(strict(*vals))
    if residual and residual.valuation()<PREC-2:
        raise ArithmeticError(("formal solve precision",residual.valuation()))
    return root

def component_branch(rec, factor_index=0):
    """
    Return (component_name, required_order, chart, strict, origin_map,
            chart_values(u,x,y) in PS at generic point).
    """
    label=rec["label"]
    required={
        "E01":2,"E02":4,"E03":3,"E04":6,"E05":5,"E06":8,
        "E07":7,"E08":10,"E09":9,"E10":6,"E11":11,"E12":1,
    }[label]

    factors=rec["tangent_factors"]
    if label=="E10":
        kind="x"
        fac=factors[factor_index]["factor"]
        cname=label+("a" if factor_index==0 else "b")
    else:
        kind="u"
        fac=factors[0]["factor"]
        cname=label

    data=chart_data(rec,kind)
    strict=parse_s(data["strict_transform"])
    om=tuple(parse_s(v) for v in data["origin_map"])
    q=factor_chart(fac,kind)

    if label=="E10":
        # x-chart: exceptional x=0.  q(u,y)=0 is linear in u.
        q0=S(q(u,0,y))
        cu=S(q0.coefficient({u:1}))
        # robust polynomial coefficient extraction
        qu=S(q0)
        cu=qu.monomial_coefficient(u)
        rest=S(qu-cu*u)
        assert cu and rest.degree(u)<=0
        p_u=Kr(-rest(0,F(0),t)/cu)
        template=[PS(Kr(p_u)), ss, PS(t)]
        root=newton_solve(strict,template,0,p_u)
        vals=[root,ss,PS(t)]
    else:
        qu=S(q)
        # Double-line support: q=y in the u-chart.
        if qu.degree()==1 and qu.monomial_coefficient(y):
            # component: u=0,y=0, generic x=t.
            template=[PS(Kr(0)),PS(t),ss]
            root=newton_solve(strict,template,0,Kr(0))
            vals=[root,PS(t),ss]
        else:
            # All irreducible conic cases become linear in x in u-chart:
            # q(1,x,y)=c0 + x + c2*y^2.
            cx=qu.monomial_coefficient(x)
            assert cx
            rest=S(qu-cx*x)
            p_x=Kr(-rest(F(0),F(0),t)/cx)
            template=[ss,PS(Kr(p_x)),PS(t)]
            root=newton_solve(strict,template,1,p_x)
            vals=[ss,root,PS(t)]

    # Generic component replay: strict surface vanishes.
    hv=PS(strict(*vals))
    assert (not hv) or hv.valuation()>=PREC-2

    return cname,required,kind,strict,om,vals

def eval_origin_map(om,vals):
    return tuple(PS(expr(*vals)) for expr in om)

def G_series(AA,BB,om,vals):
    ub,xb,yb=eval_origin_map(om,vals)
    AAv=univar_at_series(AA,ub)
    BBv=univar_at_series(BB,ub)
    Zv=univar_at_series(Z,ub)
    Xv=univar_at_series(X,ub)
    Yv=univar_at_series(Y,ub)
    return AAv*(Zv**2*xb-Xv) + BBv*(Zv**3*yb+Yv)

def functional_rows(values):
    vals=[Kr(v) for v in values]
    nonzero=[v for v in vals if v]
    if not nonzero:
        return []
    base_ring=nonzero[0].denominator().parent()
    common=base_ring.one()
    for v in nonzero:
        common=common.lcm(v.denominator())
    nums=[(v*common).numerator() for v in vals]
    maxdeg=max([q.degree() for q in nums if q]+[-1])
    rows=[]
    for degree in range(maxdeg+1):
        row=[
            F(q[degree]) if q and degree<=q.degree() else F(0)
            for q in nums
        ]
        if any(row):
            rows.append(row)
    return rows

def canonical_row(row):
    row=list(row)
    pivot=next(v for v in row if v)
    return tuple(v/pivot for v in row)

all_rows=[]
seen=set()
component_diag=[]

for rec in sorted(res["centers"],key=lambda r:int(r["label"][1:])):
    branches=(0,1) if rec["label"]=="E10" else (0,)
    for fi in branches:
        cname,required,kind,strict,om,vals=component_branch(rec,fi)
        series=[G_series(AA,BB,om,vals) for AA,BB in post_pairs]

        # Verify we really chose a generic parameter: not every basis function
        # may vanish identically on the branch.
        local_rows=[]
        valuations=[]
        for fser in series:
            valuations.append(PREC if not fser else int(fser.valuation()))

        for k in range(required):
            coeffs=[
                Kr(0) if (not fser or k>=fser.prec()) else Kr(fser[k])
                for fser in series
            ]
            for row in functional_rows(coeffs):
                crow=canonical_row(row)
                if crow not in seen:
                    seen.add(crow)
                    all_rows.append(list(crow))
                    local_rows.append(list(crow))

        Llocal=matrix(F,local_rows) if local_rows else matrix(F,0,8)
        Lall=matrix(F,all_rows) if all_rows else matrix(F,0,8)
        component_diag.append({
            "component":cname,
            "required":required,
            "chart":kind,
            "local_new_rows":len(local_rows),
            "local_new_rank":int(Llocal.rank()),
            "cumulative_rank":int(Lall.rank()),
            "basis_valuations":valuations,
        })
        print(
            "Q32DIVVAL_COMPONENT|"
            f"component={cname}|required={required}|chart={kind}|"
            f"new_rows={len(local_rows)}|local_rank={Llocal.rank()}|"
            f"cumulative_rank={Lall.rank()}|"
            f"basis_v={','.join(map(str,valuations))}|status=PASS",
            flush=True,
        )

L=matrix(F,all_rows) if all_rows else matrix(F,0,8)
rank=int(L.rank())
kernel=L.right_kernel().basis_matrix()

print(
    "Q32DIVVAL_RESOLVED|"
    f"post_collision=8|rows={L.nrows()}|rank={rank}|kernel={kernel.nrows()}|"
    f"status={'PASS_H0_TWO' if rank==6 and kernel.nrows()==2 else 'DIAGNOSTIC'}",
    flush=True,
)

payload={
    "schema":"elkies-k3.h3-q32-d12-generic-divisorial-valuations-modp.v1",
    "status":(
        "PASS_Q32_GENERIC_DIVISORIAL_RR_H0_TWO"
        if rank==6 and kernel.nrows()==2 else
        "DIAGNOSTIC_Q32_GENERIC_DIVISORIAL_RR"
    ),
    "prime":int(p),
    "global":{
        "ambient":56,
        "collision_rank":48,
        "post_collision":8,
    },
    "components":component_diag,
    "resolved":{
        "row_count":int(L.nrows()),
        "rank":rank,
        "kernel_dimension":int(kernel.nrows()),
        "condition_matrix":[[int(v) for v in row] for row in L.rows()],
        "kernel_basis":[[int(v) for v in row] for row in kernel.rows()],
    },
    "interpretation":(
        "Conditions are generic divisorial valuations in the resolved surface "
        "local ring. They are not maximal-ideal jets at blow-up centers."
    ),
    "next":(
        "If rank=6, lift the two kernel vectors through the post-collision "
        "basis and run the degree-two chord/binary-quartic Jacobian compiler. "
        "If rank>6, inspect only components that increase cumulative rank past "
        "6 and replace their independent traces by the connected-divisor "
        "quotient, following the Q80 A6 certificate."
    ),
}
OUT=args.output.resolve() if args.output else LOCAL/f"q32-d12-generic-divval-mod-{p}.json"
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q32DIVVAL_RESULT|"
    f"rank={rank}|kernel={kernel.nrows()}|status={payload['status']}",
    flush=True,
)
