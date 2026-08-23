#!/usr/bin/env sage -python
# Direct chord-slope compiler for the current-equation orbit42 section.

import argparse,json
from pathlib import Path
from sage.all import GF, PolynomialRing, ZZ

ROOT=Path(__file__).resolve().parents[2]
LOCAL=ROOT/"artifacts/local/elkies-k3"

parser=argparse.ArgumentParser()
parser.add_argument("--prime",type=int,default=100003)
parser.add_argument("--input",type=Path)
parser.add_argument("--output",type=Path)
args=parser.parse_args()

p=ZZ(args.prime)
F=GF(p)
IN=(args.input.resolve() if args.input else
    LOCAL/f"q24-orbit42-current-equation-section-mod-{p}.json")
if not IN.exists():
    raise SystemExit(f"missing explicit P42 artifact: {IN}")

data=json.loads(IN.read_text())
assert data["status"]=="PASS_Q24_ORBIT42_CURRENT_EQUATION_SECTION_MODP"
model=data["actual_twist_model"]
targets=data["orbit42_section_candidates"]
if not targets:
    raise SystemExit("explicit P42 artifact contains no target sections")

Ru=PolynomialRing(F,"u")
u=Ru.gen()
A=Ru(model["A"])
B=Ru(model["B"])

Rm=PolynomialRing(F,"m")
m0=Rm.gen()
Km=Rm.fraction_field()
R=PolynomialRing(Km,"u")
m=Km(m0)

def lift_poly(P):
    P=Ru(P)
    return R([Km(c) for c in P.list()])

def binary_invariants(Q):
    Q=R(Q)
    cs=[Km(0)]*5
    for i,c in enumerate(Q.list()):
        if i<5:
            cs[i]=Km(c)
    e,d,c,b,a=cs
    I=12*a*e-3*b*d+c*c
    J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c*c*c
    return I,J

def as_poly_m(v):
    v=Km(v)
    den=Rm(v.denominator())
    if den.degree()>0:
        return None
    return Rm(v.numerator())/F(den[0])

results=[]

for ti,t in enumerate(targets):
    X=Ru(t["X"])
    Y=Ru(t["Y"])
    Z=Ru(t["Z"])
    assert Z.degree()==3
    assert Y**2==X**3+A*X*Z**4+B*Z**6

    XR=lift_poly(X)
    YR=lift_poly(Y)
    ZR=lift_poly(Z)
    AR=lift_poly(A)

    raw=R(
        m**4*ZR**4
        - 6*XR*m**2*ZR**2
        - 8*YR*m*ZR
        - 3*XR**2
        - 4*AR*ZR**4
    )
    fac=raw.factor()

    reduced=R(fac.unit())
    square=R.one()
    factor_profile=[]
    for f,e in fac:
        e=int(e)
        factor_profile.append((int(f.degree()),e))
        square*=f**(e//2)
        if e%2:
            reduced*=f

    reduced=R(reduced)
    if reduced.degree() not in (3,4):
        print(
            "Q24O42CHORD_REDUCE|"
            f"target={ti}|raw_degree={raw.degree()}|"
            f"reduced_degree={reduced.degree()}|"
            f"factors={factor_profile}|status=NOT_GENUS_ONE_QUARTIC",
            flush=True,
        )
        continue

    I,J=binary_invariants(reduced)
    jacA=-27*I
    jacB=-27*J
    Delta=-16*(4*jacA**3+27*jacB**2)

    Ap=as_poly_m(jacA)
    Bp=as_poly_m(jacB)
    Dp=as_poly_m(Delta)

    fibre_profile=[]
    a11=False
    if Ap is not None and Bp is not None and Dp is not None:
        Ap=Rm(Ap); Bp=Rm(Bp); Dp=Rm(Dp)

        for f,e in Dp.factor():
            e=int(e)
            ta=Ap; va=0
            while ta and ta%f==0:
                va+=1; ta//=f
            tb=Bp; vb=0
            while tb and tb%f==0:
                vb+=1; tb//=f
            kind=(f"I{e}" if va==0 and vb==0 else f"vA{va}_vB{vb}_vD{e}")
            fibre_profile.append({
                "degree":int(f.degree()),
                "multiplicity":e,
                "kind":kind,
                "factor":[int(x) for x in f.list()],
            })
            if e==12 and va==0 and vb==0:
                a11=True

        vAinf=8-Ap.degree()
        vBinf=12-Bp.degree()
        vDinf=24-Dp.degree()
        if vDinf>0:
            kind=(f"I{vDinf}" if vAinf==0 and vBinf==0 else
                  f"vA{vAinf}_vB{vBinf}_vD{vDinf}")
            fibre_profile.append({
                "degree":"infinity",
                "multiplicity":int(vDinf),
                "kind":kind,
            })
            if vDinf==12 and vAinf==0 and vBinf==0:
                a11=True

    print(
        "Q24O42CHORD|"
        f"target={ti}|raw_degree={raw.degree()}|reduced_degree={reduced.degree()}|"
        f"square_degree={square.degree()}|"
        f"A11_I12={int(a11)}|"
        f"status={'PASS_A11_MODP' if a11 else 'QUARTIC_COMPILED'}",
        flush=True,
    )

    results.append({
        "target_index":ti,
        "X":t["X"],"Y":t["Y"],"Z":t["Z"],
        "raw_degree":int(raw.degree()),
        "raw_factor_profile":[list(x) for x in factor_profile],
        "square_factor_degree":int(square.degree()),
        "reduced_quartic_degree":int(reduced.degree()),
        "reduced_quartic_coefficients_low_to_high":[str(x) for x in reduced.list()],
        "jacobian_A":str(jacA),
        "jacobian_B":str(jacB),
        "discriminant":str(Delta),
        "fibre_profile":fibre_profile,
        "has_I12_A11":bool(a11),
    })

status=(
    "PASS_Q24_ORBIT42_A11_CHORD_MODP"
    if any(x["has_I12_A11"] for x in results)
    else "Q24_ORBIT42_CHORD_COMPILED_NEEDS_MINIMALITY_OR_MARKING_CHECK"
    if results
    else "Q24_ORBIT42_CHORD_NO_QUARTIC"
)

payload={
    "schema":"elkies-k3.h3-q24-orbit42-a11-chord-modp.v1",
    "status":status,
    "prime":int(p),
    "source":str(IN.relative_to(ROOT)),
    "results":results,
}
OUT=(args.output.resolve() if args.output else
     LOCAL/f"q24-orbit42-a11-chord-mod-{p}.json")
OUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
print(f"OUTPUT|{OUT}",flush=True)
print(
    "Q24O42CHORD_RESULT|"
    f"targets={len(results)}|A11={sum(x['has_I12_A11'] for x in results)}|"
    f"status={status}",
    flush=True,
)
