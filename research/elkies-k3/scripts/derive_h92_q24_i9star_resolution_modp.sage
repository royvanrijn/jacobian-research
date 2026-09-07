#!/usr/bin/env sage -python
"""
Resolve the actual H3 D13 parent's I9* surface germ over GF(p) by explicit
ordinary blow-up charts.

This is the D13 analogue of derive_h92_q6_e8_resolution.sage.  It does NOT
infer geometry from the Kodaira symbol.  Starting from the actual canonical
D13 Weierstrass equation at its rational I9* place, it recursively:

  * locates singular points on each exceptional divisor;
  * blows them up in all three affine charts;
  * deduplicates overlap charts by projective tangent direction;
  * tracks the map back to the original local (u,x,y) coordinates;
  * verifies that terminal charts are smooth above every center;
  * distinguishes blow-up centres from irreducible exceptional components by
    recording the tangent cone at every centre.

A point blow-up need not create only one geometric exceptional component:
a multiplicity-two tangent cone of rank two is a pair of lines over the
algebraic closure.  Therefore the number of blow-up centres is not itself the
D13 component count.  The expected minimal resolution has 13 geometric
exceptional curves.  This artifact is still only the geometric chart tree;
marked-chord quotient conditions are deliberately deferred.
"""

import argparse
import json
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
parser.add_argument("--max-centers",type=int,default=20)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"
q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
Q8=next((
    path for path in q8_candidates
    if path.exists() and json.loads(path.read_text()).get("status")=="PASS_EXACT_CORRECTED_Q8_D13_CHILD"
),None)
if Q8 is None:
    raise SystemExit("No passing exact D13 q8 child artifact")
q8=json.loads(Q8.read_text())
child=q8["child"]

p=ZZ(args.prime)
F=GF(p)
if F.characteristic() == 2:
    raise ValueError("tangent-cone component accounting requires odd characteristic")
Q=PolynomialRing(QQ,"U")
UQ=Q.gen()
B=PolynomialRing(F,"U")
UB=B.gen()


def red_q(q):
    q=QQ(q)
    d=ZZ(q.denominator())
    if d%p==0:
        raise ZeroDivisionError(f"denominator divisible by {p}")
    return F(ZZ(q.numerator()))/F(d)


def red_poly(vals):
    return B([red_q(QQ(v)) for v in vals])


A=red_poly(child["minimal_A_coefficients_low_to_high"])
Bc=red_poly(child["minimal_B_coefficients_low_to_high"])
Delta=-16*(4*A**3+27*Bc**2)
i9=next(item for item in child["finite_fibres"] if item["kodaira"]=="I9*")
fQ=Q(str(i9["factor"]))
f=B([red_q(c) for c in fQ.list()])
assert f.degree()==1
alpha=-f[0]/f[1]

S=PolynomialRing(F,names=("u","x","y"),order="degrevlex")
u,x,y=S.gens()

Al=S(A(alpha+u))
Bl=S(Bc(alpha+u))

# Check local minimal orders in one-variable ring separately.
T=PolynomialRing(F,"t")
t=T.gen()
Al1=T(A(alpha+t)); Bl1=T(Bc(alpha+t)); Del1=T(Delta(alpha+t))
assert (Al1.valuation(),Bl1.valuation(),Del1.valuation())==(2,3,15)

surface=y**2-x**3-Al*x-Bl
assert surface(0,0,0)==0


def hessian_matrix(poly, point):
    return matrix(
        F,
        [[poly.derivative(a,b)(*point) for b in (u,x,y)] for a in (u,x,y)],
    )


def hessian_det(poly, point):
    return hessian_matrix(poly, point).det()


def shifted_polynomial(poly, point):
    a,b,c=point
    return S(poly(a+u,b+x,c+y))


def order_at_point(poly, point):
    shifted=shifted_polynomial(poly,point)
    if not shifted:
        return 10**9
    return min(sum(exp) for exp,coef in shifted.dict().items() if coef)


def tangent_cone_record(poly, point, multiplicity):
    """Record the homogeneous tangent cone and geometric component count.

    Every centre encountered here has multiplicity two.  In odd
    characteristic a quadratic tangent cone has:

      rank 1: one doubled line;
      rank 2: two geometric lines (possibly conjugate over F);
      rank 3: one smooth conic.

    Hence only rank two contributes two irreducible geometric components.
    We retain the factorization over F as a separate rationality diagnostic.
    """
    shifted=shifted_polynomial(poly,point)
    cone=S(sum(
        coefficient * u**exponent[0] * x**exponent[1] * y**exponent[2]
        for exponent,coefficient in shifted.dict().items()
        if sum(exponent)==multiplicity
    ))
    if not cone:
        raise ArithmeticError("empty tangent cone")
    factorization=cone.factor()
    factors=[
        {"factor":str(factor),"multiplicity":int(exponent)}
        for factor,exponent in factorization
    ]

    quadratic_rank=None
    geometric_components=None
    if multiplicity==2:
        quadratic_rank=int(hessian_matrix(poly,point).rank())
        if quadratic_rank not in (1,2,3):
            raise ArithmeticError(
                f"unexpected quadratic tangent-cone rank {quadratic_rank}"
            )
        geometric_components=2 if quadratic_rank==2 else 1

    return {
        "polynomial":str(cone),
        "factorization_over_base_field":factors,
        "quadratic_rank":quadratic_rank,
        "geometric_exceptional_components":geometric_components,
    }


def divide_power(poly, exceptional, power):
    q=S(poly)
    for _ in range(power):
        q,rem=q.quo_rem(exceptional)
        if rem:
            raise ArithmeticError("strict-transform exceptional power does not divide")
    return S(q)


def canonical_projective(direction):
    vals=[F(v) for v in direction]
    idx=next((i for i,v in enumerate(vals) if v),None)
    if idx is None:
        raise ArithmeticError("zero projective direction")
    inv=vals[idx]**-1
    return tuple(int(v*inv) for v in vals)


def singular_points_on_exceptional(poly, exceptional):
    ideal=S.ideal([
        exceptional,
        poly,
        poly.derivative(u),
        poly.derivative(x),
        poly.derivative(y),
    ])
    if ideal.is_one():
        return []
    dim=ideal.dimension()
    if dim != 0:
        raise ArithmeticError(
            f"exceptional singular locus not zero-dimensional: dim={dim}; "
            f"gb={ideal.groebner_basis()}"
        )
    try:
        variety=ideal.variety(ring=F)
    except TypeError:
        variety=ideal.variety()
    points=[]
    for sol in variety:
        point=tuple(F(sol.get(v,F(0))) for v in (u,x,y))
        if point not in points:
            points.append(point)
    # Strong replay.
    for pt in points:
        assert exceptional(*pt)==0
        assert poly(*pt)==0
        assert all(poly.derivative(v)(*pt)==0 for v in (u,x,y))
    return points


def chart_substitutions(point, kind):
    a,b,c=map(F,point)
    if kind=="u":
        return (a+u,b+u*x,c+u*y),u
    if kind=="x":
        return (a+x*u,b+x,c+x*y),x
    if kind=="y":
        return (a+y*u,b+y*x,c+y),y
    raise ValueError(kind)


def direction_from_chart(kind, point):
    pu,px,py=map(F,point)
    if kind=="u":
        assert pu==0
        return canonical_projective((1,px,py))
    if kind=="x":
        assert px==0
        return canonical_projective((pu,1,py))
    assert kind=="y" and py==0
    return canonical_projective((pu,px,1))


center_records=[]
leaf_records=[]
next_center=0


def blow_center(poly, origin_map, component_eqs, point, path, depth):
    global next_center
    if len(center_records)>=args.max_centers:
        raise RuntimeError("resolution exceeded max centers")

    multiplicity=order_at_point(poly,point)
    if multiplicity < 2:
        raise ArithmeticError(
            ("attempted blow-up of smooth point",path,point,multiplicity)
        )
    tangent=tangent_cone_record(poly,point,multiplicity)
    ordinary=bool(hessian_det(poly,point))
    active=[name for name,g in component_eqs.items() if g(*point)==0]

    label=f"C{next_center+1:02d}"
    next_center += 1
    record={
        "label":label,
        "path":path,
        "depth":int(depth),
        "point":[int(v) for v in point],
        "multiplicity":int(multiplicity),
        "ordinary_double_point":ordinary,
        "active_components":active,
        "tangent_cone":tangent,
        "children":[],
    }
    center_records.append(record)

    # Blow up in all affine charts and gather singular points on the NEW
    # exceptional divisor.  The chart equation `exceptional=0` can represent
    # more than one geometric component when the tangent cone has rank two;
    # that split is recorded above and must be separated by the later
    # resolved-component quotient compiler.
    candidates={}
    chart_diagnostics=[]
    for kind in ("u","x","y"):
        subs,e=chart_substitutions(point,kind)
        transformed=S(poly(*subs))
        strict=divide_power(transformed,e,multiplicity)
        new_map=tuple(S(expr(*subs)) for expr in origin_map)

        new_components={}
        for name,g in component_eqs.items():
            ordg=order_at_point(g,point)
            gt=S(g(*subs))
            if ordg>0 and ordg<10**9:
                gt=divide_power(gt,e,ordg)
            new_components[name]=gt
        new_components[label]=e

        singular=singular_points_on_exceptional(strict,e)
        chart_diagnostics.append({
            "chart":kind,
            "exceptional_coordinate":str(e),
            "strict_transform":str(strict),
            "exceptional_restriction":str(strict.subs({e:0})),
            "origin_map":[str(v) for v in new_map],
            "singular_points":[[int(v) for v in q] for q in singular],
        })
        for q in singular:
            d=direction_from_chart(kind,q)
            candidates.setdefault(d,[]).append(
                (kind,q,strict,new_map,new_components)
            )

    record["charts"]=chart_diagnostics

    # Deduplicate overlaps by projective tangent direction and recurse using a
    # deterministic preferred chart u,x,y.
    preference={"u":0,"x":1,"y":2}
    for direction,reps in sorted(candidates.items()):
        reps.sort(key=lambda item:preference[item[0]])
        kind,q,strict,new_map,new_components=reps[0]
        child_path=f"{path}/{label}:{kind}:{','.join(map(str,direction))}"
        child={
            "direction":list(direction),
            "representatives":[
                {"chart":r[0],"point":[int(v) for v in r[1]]}
                for r in reps
            ],
            "selected_chart":kind,
        }
        record["children"].append(child)
        blow_center(strict,new_map,new_components,q,child_path,depth+1)

    if not candidates:
        leaf_records.append({
            "after_center":label,
            "path":path,
            "ordinary_center":ordinary,
            "status":"SMOOTH_ABOVE_EXCEPTIONAL",
        })


# Start at the unique singular point of the Weierstrass germ.  Track the
# proper transform of the original fibre u=0 as component F0.
initial_components={"F0":u}
blow_center(surface,(u,x,y),initial_components,(F(0),F(0),F(0)),"root",0)

all_resolved=all(
    rec.get("status")=="SMOOTH_ABOVE_EXCEPTIONAL" for rec in leaf_records
)
center_count=len(center_records)
ordinary_count=sum(
    int(r["ordinary_double_point"]) for r in center_records
)
nonordinary_count=center_count-ordinary_count
if any(
    r["tangent_cone"]["geometric_exceptional_components"] is None
    for r in center_records
):
    geometric_component_count=None
else:
    geometric_component_count=sum(
        int(r["tangent_cone"]["geometric_exceptional_components"])
        for r in center_records
    )
split_centers=[
    r["label"] for r in center_records
    if r["tangent_cone"]["geometric_exceptional_components"]==2
]
passes_d13=bool(
    all_resolved
    and geometric_component_count==13
)

print(
    "Q24I9RESOLVE|"
    f"prime={p}|base={int(alpha)}|orders=2,3,15|"
    f"centers={center_count}|components={geometric_component_count}|"
    f"split_centers={','.join(split_centers)}|"
    f"ordinary={ordinary_count}|nonordinary={nonordinary_count}|"
    f"leaves={len(leaf_records)}|expected_components=13|"
    f"status={'PASS_D13_COMPONENT_COUNT' if passes_d13 else 'DIAGNOSTIC'}",
    flush=True,
)
for rec in center_records:
    tangent=rec["tangent_cone"]
    print(
        "Q24I9CENTER|"
        f"label={rec['label']}|depth={rec['depth']}|point={rec['point']}|"
        f"mult={rec['multiplicity']}|ordinary={int(rec['ordinary_double_point'])}|"
        f"tangent_rank={tangent['quadratic_rank']}|"
        f"new_components={tangent['geometric_exceptional_components']}|"
        f"active={','.join(rec['active_components'])}|"
        f"children={len(rec['children'])}",
        flush=True,
    )

payload={
    "schema":"elkies-k3.h3-q24-i9star-blowup-resolution-modp.v2",
    "status":(
        "PASS_EXPLICIT_MODP_I9STAR_D13_COMPONENT_RESOLUTION"
        if passes_d13 else
        "DIAGNOSTIC_I9STAR_BLOWUP_TREE"
    ),
    "prime":int(p),
    "base":int(alpha),
    "local_orders":[2,3,15],
    "expected_root_lattice":"D13",
    "expected_exceptional_curves":13,
    "actual_blowup_centers":int(center_count),
    "actual_geometric_exceptional_components":(
        int(geometric_component_count)
        if geometric_component_count is not None else None
    ),
    "split_tangent_cone_centers":split_centers,
    "ordinary_centers":int(ordinary_count),
    "nonordinary_centers":int(nonordinary_count),
    "centers":center_records,
    "leaves":leaf_records,
    "boundary":(
        "This artifact resolves the actual I9* surface germ by ordinary "
        "blow-up charts and distinguishes blow-up centres from irreducible "
        "geometric exceptional components using the tangent cone.  A rank-two "
        "quadratic tangent cone contributes two geometric lines.  The current "
        "chart records do not yet assign separate component names/trivializations "
        "to those two branches, identify the q24 vertical class with the full "
        "geometric D13 basis, or impose marked-chord RR quotient conditions."
    ),
    "next":(
        "Split every rank-two tangent-cone exceptional into its two geometric "
        "branches, build the full D13 intersection graph, match it to the "
        "deterministic lattice root basis, then pull the modular q24 chord and "
        "10-dimensional post-collision RR space through the resulting charts "
        "to derive the connected codimension-eight quotient."
    ),
}
OUT=(
    args.output.resolve()
    if args.output else
    LOCAL/f"q24-i9star-resolution-mod-{p}.json"
)
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24I9RESOLVE_RESULT|"
    f"centers={center_count}|components={geometric_component_count}|"
    f"expected_components=13|status={payload['status']}",
    flush=True,
)
