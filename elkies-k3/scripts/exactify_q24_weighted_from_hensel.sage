#!/usr/bin/env sage -python
"""
Try to exactify q24 from an existing high-precision Hensel checkpoint by
undoing the potentially disastrous monic-Z weighted normalization.

If
    x = X/Z^2, y = Y/Z^3,
then for any nonzero c,
    Z* = c Z, X* = c^2 X, Y* = c^3 Y
represents the same section and satisfies the same cleared Weierstrass identity.

We infer c from coefficients that already rational-reconstruct, rescale ALL
p-adic residues modulo p^N, center them, and test the exact QQ[U] identity.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ, gcd, lcm


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

INPUT=(
    args.input.resolve() if args.input
    else LOCAL/"q24-direct-hensel-p2048.json"
)
OUTPUT=(
    args.output.resolve() if args.output
    else LOCAL/"q8-q24-horizontal-section-qq-weighted.json"
)

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
assert hen["jacobian_rank"]==156
assert hen["final_residual_valuation"]>=hen["precision"]
assert len(hen["residues"])==156

p=ZZ(hen["prime"])
N=int(hen["precision"])
M=p**N
res=[ZZ(v)%M for v in hen["residues"]]

print(
    f"Q24WEIGHT_SETUP|prime={p}|precision={N}|modulus_bits={M.nbits()}|status=PASS",
    flush=True,
)

def rr(r):
    try:
        return ZZ(r).rational_reconstruction(M)
    except (ArithmeticError,ValueError):
        return None

recon=[rr(r) for r in res]
rz=recon[:24]
rx=recon[24:77]
ry=recon[77:]

print(
    "Q24WEIGHT_RR|"
    f"Z={sum(q is not None for q in rz)}/24|"
    f"X={sum(q is not None for q in rx)}/53|"
    f"Y={sum(q is not None for q in ry)}/79|status=PASS",
    flush=True,
)

def lcm_den(values):
    c=ZZ.one()
    for q in values:
        if q is not None:
            c=lcm(c,ZZ(q.denominator()))
    return c

# Strongest natural guess: denominator of monic Z coordinates is the
# leading coefficient of a primitive weighted-integral representative.
cZ=lcm_den(rz)

# Augment without factoring giant denominators. For weight w require d|c^w;
# multiplying c by d/gcd(d,c^w) guarantees the condition.
def augment(c,values,w):
    c=ZZ(c)
    for q in values:
        if q is None:
            continue
        d=ZZ(q.denominator())
        rem=d//gcd(d,c**w)
        if rem!=1:
            c*=rem
    return c

cCompat=augment(augment(cZ,rx,2),ry,3)

candidates=[]
for label,c in (("Z_LCM",cZ),("COMPAT",cCompat)):
    c=abs(ZZ(c))
    if c and c not in [x[1] for x in candidates]:
        candidates.append((label,c))

RQ=PolynomialRing(QQ,"U")
A=RQ([QQ(v) for v in child["child"]["minimal_A_coefficients_low_to_high"]])
B=RQ([QQ(v) for v in child["child"]["minimal_B_coefficients_low_to_high"]])

def centered(v):
    v=ZZ(v)%M
    return v-M if v>M//2 else v

def try_scale(label,c):
    cp1=c%M
    cp2=(cp1*cp1)%M
    cp3=(cp2*cp1)%M

    Zs=[centered(res[i]*cp1) for i in range(24)] + [ZZ(c)]
    Xs=[centered(res[24+i]*cp2) for i in range(53)]
    Ys=[centered(res[77+i]*cp3) for i in range(79)]

    Z=RQ(Zs)
    X=RQ(Xs)
    Y=RQ(Ys)

    maxbits=max(abs(v).nbits() for v in Zs+Xs+Ys if v)
    identity=Y**2-X**3-A*X*Z**4-B*Z**6

    print(
        "Q24WEIGHT_TRY|"
        f"label={label}|c_bits={c.nbits()}|max_scaled_bits={maxbits}|"
        f"identity={int(not identity)}|"
        f"identity_degree={-1 if not identity else identity.degree()}|status={'PASS' if not identity else 'MISS'}",
        flush=True,
    )

    if identity:
        return None

    # Normalize back to monic Z.
    cq=QQ(c)
    Zn=Z/cq
    Xn=X/(cq**2)
    Yn=Y/(cq**3)

    assert Zn.degree()==24 and Zn.leading_coefficient()==1
    assert Xn.degree()==52 and Yn.degree()==78
    assert Yn**2-Xn**3-A*Xn*Zn**4-B*Zn**6 == 0

    # Literal reduction to the certified modular q24 seed, if present.
    MOD=LOCAL/"q24-degree46-direct-global-mod-100003.json"
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
                    d=ZZ(q.denominator())
                    if d%p==0:
                        raise ZeroDivisionError
                    vals.append(F(ZZ(q.numerator()))/F(d))
                return RF(vals)
            sec=mod["section_mod_p"]
            Zm=RF([F(int(v)) for v in sec["Z_coefficients_low_to_high"]])
            Xm=RF([F(int(v)) for v in sec["X_coefficients_low_to_high"]])
            Ym=RF([F(int(v)) for v in sec["Y_coefficients_low_to_high"]])
            mod_match=(red(Zn)==Zm and red(Xn)==Xm and red(Yn)==Ym)
            if not mod_match:
                raise ArithmeticError("exact weighted solution does not match certified modular q24")

    payload={
        "schema":"elkies-k3.h92-q8-q24-horizontal-section-qq.weighted.v1",
        "status":"PASS_EXACT_Q24_HORIZONTAL_SECTION",
        "method":"weighted integral recovery from p-adic checkpoint",
        "weighted_scale":str(c),
        "source_checkpoint":str(INPUT.relative_to(ROOT)),
        "profile":{
            "Z_degree":24,"X_degree":52,"Y_degree":78,
            "x_degrees":[52,48],"y_degrees":[78,72],
        },
        "section":{
            "Z_coefficients_low_to_high":[str(v) for v in Zn.list()],
            "X_coefficients_low_to_high":[str(v) for v in Xn.list()],
            "Y_coefficients_low_to_high":[str(v) for v in Yn.list()],
            "x_numerator_coefficients_low_to_high":[str(v) for v in Xn.list()],
            "x_denominator_coefficients_low_to_high":[str(v) for v in (Zn**2).list()],
            "y_numerator_coefficients_low_to_high":[str(v) for v in Yn.list()],
            "y_denominator_coefficients_low_to_high":[str(v) for v in (Zn**3).list()],
        },
        "verification":{
            "exact_weierstrass_identity":True,
            "reduction_matches_modular_q24":mod_match,
        },
    }
    OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    OUTPUT.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    print(f"OUTPUT|{OUTPUT}",flush=True)
    print(
        f"Q24WEIGHT_RESULT|scale_bits={c.nbits()}|max_scaled_bits={maxbits}|"
        "status=PASS_EXACT_Q24",
        flush=True,
    )
    return payload

for label,c in candidates:
    result=try_scale(label,c)
    if result is not None:
        raise SystemExit(0)

print(
    "Q24WEIGHT_RESULT|status=NO_EXACT_WEIGHTED_SCALE_FROM_RECONSTRUCTED_DENOMINATORS|"
    "next=continue_hensel_or_use_AJ_S3",
    flush=True,
)
raise SystemExit(2)
