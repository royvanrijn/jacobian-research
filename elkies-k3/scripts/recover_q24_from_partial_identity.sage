#!/usr/bin/env sage -python
"""
Recover q24 exactly from the PARTIALLY rational-reconstructed p2048 Hensel
checkpoint by exploiting the exact Weierstrass identity as an overdetermined
algebraic system.

At p2048 we already have:
  * all 24 non-leading coefficients of monic Z,
  * most coefficients of X,
  * most coefficients of Y.

Instead of increasing p-adic precision, substitute every known rational
coefficient into

    Y^2 - X^3 - A X Z^4 - B Z^6 = 0

and peel the remaining unknown coefficients using exact coefficient equations.

As soon as X is completely recovered, Y is obtained as the exact polynomial
square root of X^3 + A X Z^4 + B Z^6; there is no need to reconstruct the
remaining Y coefficients independently.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, matrix, vector


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
parser.add_argument("--input",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

ROOT=locate_repo(args.repo)
LOCAL=ROOT/"artifacts/local/elkies-k3"
GEN=ROOT/"artifacts/generated-results"
INPUT=args.input.resolve() if args.input else LOCAL/"q24-direct-hensel-p2048.json"
OUTPUT=args.output.resolve() if args.output else LOCAL/"q8-q24-horizontal-section-qq-partial-completion.json"
MOD=LOCAL/"q24-degree46-direct-global-mod-100003.json"

q8_candidates=[
    LOCAL/"q8-corrected2cover-qq-child.json",
    GEN/"elkies-k3-h92-q6-child-q8-corrected2cover-qq-child.json",
]
CHILD=next((
    p for p in q8_candidates
    if p.exists()
    and "child" in json.loads(p.read_text())
    and "minimal_A_coefficients_low_to_high"
       in json.loads(p.read_text()).get("child",{})
),None)

if not INPUT.exists():
    raise SystemExit(f"Missing checkpoint: {INPUT}")
if CHILD is None:
    raise SystemExit("Missing exact D13 child")

hen=json.loads(INPUT.read_text())
child=json.loads(CHILD.read_text())

assert hen["status"]=="PASS_Q24_DIRECT_HENSEL"
assert int(hen["precision"])>=2048
assert int(hen["final_residual_valuation"])>=int(hen["precision"])
assert int(hen["jacobian_rank"])==156

p=ZZ(hen["prime"])
M=p**int(hen["precision"])
res=[ZZ(v)%M for v in hen["residues"]]

def parse_partial(values):
    return [None if v is None else QQ(v) for v in values]

# Prefer the checkpoint's own rational reconstructions.
Zpart=parse_partial(hen["Z"])
Xpart=parse_partial(hen["X"])
Ypart=parse_partial(hen["Y"])
assert len(Zpart)==25 and len(Xpart)==53 and len(Ypart)==79
assert Zpart[24]==1

if any(v is None for v in Zpart):
    raise SystemExit("Z is not fully reconstructed; this analytical completion expects exact Z")

ux=[i for i,v in enumerate(Xpart) if v is None]
uy=[i for i,v in enumerate(Ypart) if v is None]

print(
    "Q24PARTIAL_SETUP|"
    f"Z_exact=25/25|X_exact={53-len(ux)}/53|Y_exact={79-len(uy)}/79|"
    f"X_unknown={','.join(map(str,ux))}|Y_unknown={','.join(map(str,uy))}|status=PASS",
    flush=True,
)

# Unknown variables only for coefficients not already recognized.
names=[f"x{i}" for i in ux]+[f"y{i}" for i in uy]
if not names:
    raise SystemExit("Nothing to solve")
P=PolynomialRing(QQ,names=names,order="degrevlex")
gens=list(P.gens())
vmap={name:g for name,g in zip(names,gens)}

RU=PolynomialRing(P,"U")
U=RU.gen()

Z=RU([P(v) for v in Zpart])
X=RU([
    vmap[f"x{i}"] if Xpart[i] is None else P(Xpart[i])
    for i in range(53)
])
Y=RU([
    vmap[f"y{i}"] if Ypart[i] is None else P(Ypart[i])
    for i in range(79)
])
A=RU([P(QQ(v)) for v in child["child"]["minimal_A_coefficients_low_to_high"]])
B=RU([P(QQ(v)) for v in child["child"]["minimal_B_coefficients_low_to_high"]])

identity=Y**2-X**3-A*X*Z**4-B*Z**6
eqs=[P(identity[i]) for i in range(157)]

# Map symbolic variable back to its modular residue. Useful when a univariate
# equation has more than one rational root.
residue_by_var={}
for i in ux:
    residue_by_var[vmap[f"x{i}"]]=res[24+i]
for i in uy:
    residue_by_var[vmap[f"y{i}"]]=res[77+i]

solved={}

def modp(q):
    q=QQ(q)
    d=ZZ(q.denominator())
    if d%p==0:
        return None
    F=GF(p)
    return F(ZZ(q.numerator()))/F(d)

def residue_mod_p(var):
    return GF(p)(residue_by_var[var] % p)

def current_eq(e):
    if not solved:
        return P(e)
    return P(e.subs(solved))

def remaining_vars(e):
    return [v for v in gens if v not in solved and e.degree(v)>0]

def accept(var,val,reason):
    val=QQ(val)
    mp=modp(val)
    if mp is None or mp != residue_mod_p(var):
        return False
    solved[var]=val
    print(
        "Q24PARTIAL_SOLVE|"
        f"var={var}|reason={reason}|num_bits={abs(ZZ(val.numerator())).nbits()}|"
        f"den_bits={abs(ZZ(val.denominator())).nbits()}|"
        f"solved={len(solved)}/{len(gens)}|status=PASS",
        flush=True,
    )
    return True

def try_single_variable():
    for raw in eqs:
        e=current_eq(raw)
        if not e:
            continue
        vs=remaining_vars(e)
        if len(vs)!=1:
            continue
        v=vs[0]
        # Turn into a univariate QQ polynomial.
        QV=PolynomialRing(QQ,str(v))
        deg=e.degree(v)
        # This equation has exactly one active variable, so its coefficients
        # on v^k are rational scalars.
        ev=QV([QQ(e.monomial_coefficient(v**k)) for k in range(deg+1)])
        if ev.degree()==1:
            val=-ev[0]/ev[1]
            if accept(v,val,"univariate_linear"):
                return True
        elif 1 < ev.degree() <= 4:
            roots=ev.roots(QQ,multiplicities=False)
            matching=[r for r in roots if modp(r)==residue_mod_p(v)]
            if len(matching)==1 and accept(v,matching[0],f"univariate_deg{ev.degree()}"):
                return True
    return False

def try_linear_system():
    rem=[v for v in gens if v not in solved]
    if not rem:
        return False

    linear=[]
    for raw in eqs:
        e=current_eq(raw)
        if not e:
            continue
        if e.total_degree()<=1:
            linear.append(e)

    if not linear:
        return False

    zero_sub={v:0 for v in rem}
    rows=[]
    rhs=[]
    for e in linear:
        const=QQ(e.subs(zero_sub))
        row=[QQ(e.derivative(v)) for v in rem]
        if any(row):
            rows.append(row)
            rhs.append(-const)

    if not rows:
        return False

    AA=matrix(QQ,rows)
    bb=vector(QQ,rhs)
    aug=AA.augment(matrix(QQ,len(rhs),1,list(bb)))
    if aug.rank()!=AA.rank():
        raise ArithmeticError("exact linear subsystem became inconsistent")

    rank=AA.rank()
    print(
        f"Q24PARTIAL_LINEAR|equations={len(rows)}|variables={len(rem)}|rank={rank}|status=PASS",
        flush=True,
    )

    if rank==len(rem):
        vals=AA.solve_right(bb)
        pending=[]
        for v,val in zip(rem,vals):
            if modp(val)!=residue_mod_p(v):
                raise ArithmeticError(f"unique exact linear solution for {v} disagrees mod p")
            pending.append((v,QQ(val)))
        for v,val in pending:
            accept(v,val,"full_linear_system")
        return bool(pending)

    # RREF may still fix some pivot variables independently of free variables.
    R=aug.echelon_form()
    piv=AA.pivots()
    free=[j for j in range(len(rem)) if j not in piv]
    progress=False
    for row,pivot in enumerate(piv):
        if all(R[row,j]==0 for j in free):
            val=QQ(R[row,len(rem)])
            v=rem[pivot]
            if accept(v,val,"linear_pivot_fixed"):
                progress=True
    return progress

def build_exact_if_X_complete():
    xvals=[]
    for i in range(53):
        if Xpart[i] is not None:
            xvals.append(QQ(Xpart[i]))
        else:
            v=vmap[f"x{i}"]
            if v not in solved:
                return None
            xvals.append(QQ(solved[v]))

    RQ=PolynomialRing(QQ,"U")
    Zq=RQ(Zpart)
    Xq=RQ(xvals)
    Aq=RQ([QQ(v) for v in child["child"]["minimal_A_coefficients_low_to_high"]])
    Bq=RQ([QQ(v) for v in child["child"]["minimal_B_coefficients_low_to_high"]])
    rhs=Xq**3+Aq*Xq*Zq**4+Bq*Zq**6

    if not rhs.is_square():
        return ("X_COMPLETE_NOT_SQUARE",Zq,Xq,None)

    Yq=rhs.sqrt()
    # Choose sign from any known Y coefficient / mod-p seed.
    oriented=False
    for i,v in enumerate(Ypart):
        if v is None:
            continue
        if i<=Yq.degree() and Yq[i]==v:
            oriented=True
            break
        if i<=Yq.degree() and -Yq[i]==v:
            Yq=-Yq
            oriented=True
            break
    if not oriented:
        # fallback to mod-p leading coefficient
        F=GF(p)
        lead_seed=F(res[77+78]%p)
        ylead=modp(Yq[78])
        if ylead==lead_seed:
            oriented=True
        elif -ylead==lead_seed:
            Yq=-Yq
            oriented=True
    if not oriented:
        raise ArithmeticError("could not orient exact Y")

    assert Yq**2==rhs
    return ("PASS",Zq,Xq,Yq)


# Repeated sparse peeling.
round_no=0
while True:
    round_no+=1
    before=len(solved)

    # Single-equation peeling is cheap and tends to unlock more equations.
    while try_single_variable():
        pass

    try_linear_system()

    exact=build_exact_if_X_complete()
    if exact is not None:
        status,Zq,Xq,Yq=exact
        if status=="PASS":
            break
        print(
            "Q24PARTIAL_X|all_X_recovered=1|rhs_square=0|status=REJECT",
            flush=True,
        )

    after=len(solved)
    print(
        "Q24PARTIAL_ROUND|"
        f"round={round_no}|solved={after}/{len(gens)}|progress={after-before}|status=PASS",
        flush=True,
    )
    if after==before:
        exact=None
        break

if exact is None or exact[0]!="PASS":
    rem=[str(v) for v in gens if v not in solved]
    # Diagnostic: count equation supports after all substitutions.
    hist={}
    for raw in eqs:
        e=current_eq(raw)
        if not e:
            continue
        n=len(remaining_vars(e))
        hist[n]=hist.get(n,0)+1
    print(
        "Q24PARTIAL_RESULT|status=NEEDS_SECOND_STAGE|"
        f"solved={len(solved)}/{len(gens)}|remaining={','.join(rem)}|"
        f"support_hist={hist}",
        flush=True,
    )
    raise SystemExit(2)

unused,Zq,Xq,Yq=exact
assert Zq.degree()==24 and Zq.leading_coefficient()==1
assert Xq.degree()==52 and Yq.degree()==78

# Independent modular replay.
mod_match=None
if MOD.exists():
    mod=json.loads(MOD.read_text())
    if ZZ(mod["prime"])==p:
        F=GF(p)
        RF=PolynomialRing(F,"U")
        def red(poly):
            vals=[]
            for q in poly.list():
                q=QQ(q)
                vals.append(F(ZZ(q.numerator()))/F(ZZ(q.denominator())))
            return RF(vals)
        sec=mod["section_mod_p"]
        Zm=RF([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
        Xm=RF([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
        Ym=RF([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
        mod_match=(red(Zq)==Zm and red(Xq)==Xm and red(Yq)==Ym)
        if not mod_match:
            raise ArithmeticError("exact completion fails independent modular replay")

payload={
    "schema":"elkies-k3.h92-q8-q24-horizontal-section-qq.partial-completion.v1",
    "status":"PASS_EXACT_Q24_HORIZONTAL_SECTION",
    "method":"exact completion from partial p-adic rational reconstruction + Weierstrass identity",
    "source_checkpoint":str(INPUT.relative_to(ROOT)),
    "profile":{
        "Z_degree":24,"X_degree":52,"Y_degree":78,
        "x_degrees":[52,48],"y_degrees":[78,72],
    },
    "section":{
        "Z_coefficients_low_to_high":[str(v) for v in Zq.list()],
        "X_coefficients_low_to_high":[str(v) for v in Xq.list()],
        "Y_coefficients_low_to_high":[str(v) for v in Yq.list()],
        "x_numerator_coefficients_low_to_high":[str(v) for v in Xq.list()],
        "x_denominator_coefficients_low_to_high":[str(v) for v in (Zq**2).list()],
        "y_numerator_coefficients_low_to_high":[str(v) for v in Yq.list()],
        "y_denominator_coefficients_low_to_high":[str(v) for v in (Zq**3).list()],
    },
    "verification":{
        "exact_weierstrass_identity":True,
        "reduction_matches_degree46_modular_section":mod_match,
    },
}
OUTPUT.parent.mkdir(parents=True,exist_ok=True)
OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUTPUT}",flush=True)
print(
    "Q24PARTIAL_RESULT|X=52|Z=24|Y=78|identity=PASS|"
    f"modular_replay={int(bool(mod_match))}|status=PASS_EXACT_Q24",
    flush=True,
)
