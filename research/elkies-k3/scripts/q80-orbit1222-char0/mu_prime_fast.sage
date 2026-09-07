#!/usr/bin/env sage
"""
Fast modular mu probe for Q80 orbit 1222.

Requires q80_char0_orbit1222_cs.json (exact monic C,S over QQ(sqrt(-3))).
For one split prime p:
  - use only s=+sqrt(-6);
  - compute a few classic brnoeth fibers for j=+sqrt(-3);
  - derive mu from j(V)=1728*C^3/(C^3-mu*S^2);
  - repeat for j=-sqrt(-3);
  - recover mu=a+b*sqrt(-3) mod p.

No 49-point interpolation is performed.
"""

import argparse
import contextlib
import io
import json
import time
from pathlib import Path

from sage.all import GF, Matrix, PolynomialRing, QQ
from sage.env import SAGE_SHARE
from sage.interfaces.singular import singular

HERE=Path(__file__).resolve().parent
REPO_ROOT=HERE.parents[2]
ORBIT_DATA=(
    REPO_ROOT/"elkies-k3"/"data"/"fibrations"/"q80-orbit1222-char0"
)
MODULAR_DATA=(
    REPO_ROOT/"artifacts"/"generated-results"/"q80-orbit1222-char0"/"modular"
)
SOURCE=REPO_ROOT/"elkies-k3/scripts/derive_q80_third_q12_cm24_pencil.sage"
CSFILE=ORBIT_DATA/"q80_char0_orbit1222_cs.json"

ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument("--prime",type=int,required=True)
ap.add_argument("--checks",type=int,default=2)
ap.add_argument("--out",type=Path,required=True)
args=ap.parse_args()

if args.checks < 1:
    raise SystemExit("--checks must be >=1")
if not SOURCE.exists() or not CSFILE.exists():
    raise SystemExit("missing exact source pencil or q80_char0_orbit1222_cs.json")

_worker_file=__file__
__file__=str(SOURCE)
try:
    with contextlib.redirect_stdout(io.StringIO()):
        load(str(SOURCE))
finally:
    __file__=_worker_file

p=args.prime
F=GF(p)
if not F(-6).is_square() or not F(-3).is_square():
    raise SystemExit(f"p={p} does not split both -6 and -3")

def canonical_root(x):
    r=F(x).sqrt()
    a=int(r); b=(-a)%p
    return F(min(a,b))

rs=canonical_root(-6)
rj=canonical_root(-3)

cs=json.loads(CSFILE.read_text())

def qmod(text):
    q=QQ(text)
    return F(q.numerator())/F(q.denominator())

def eval_coeff(c,jroot):
    # C,S are already certified to have zero s and s*j components.
    return qmod(c[0])+jroot*qmod(c[2])

def poly_from_exact(coeffs,jroot,name):
    R=PolynomialRing(F,name)
    return R([eval_coeff(c,jroot) for c in coeffs])

def pad(values,n,zero):
    values=list(values)
    return values+[zero]*(n-len(values))

def reduce_Q(value):
    value=QQ(value)
    return F(value.numerator())/F(value.denominator())

def reduce_K(value,sroot):
    c=pad(K(value).list(),2,QQ(0))
    return reduce_Q(c[0])+sroot*reduce_Q(c[1])

def reduce_L(value,sroot,jroot):
    c=pad(L(value).list(),2,K(0))
    return reduce_K(c[0],sroot)+jroot*reduce_K(c[1],sroot)

finite_plane=PolynomialRing(F,names=("w","x"),order="lex")
fw,fx=finite_plane.gens()

singular.set_ring(finite_plane._singular_())
singular.lib("brnoeth.lib")

# Same required Singular 4.4.x compatibility patch as the certified worker.
hnoether_source=Path(SAGE_SHARE)/"singular/LIB/hnoether.lib"
text=hnoether_source.read_text()
a=text.index("proc extdevelop (list l, int Exaktheit)")
b=text.index("\nexample\n",a)
patched=text[a:b]
renames={
    "ideal a=hole(lastrow);":"ideal q80row=hole(lastrow);",
    "else { ideal a=lastrow; }":"else { ideal q80row=lastrow; }",
    "a[Q]=delt;":"q80row[Q]=delt;",
    "a[Q+1]=x;":"q80row[Q+1]=x;",
    "lastrow=zurueck(a);":"lastrow=zurueck(q80row);",
    "else { lastrow=a; }":"else { lastrow=q80row; }",
}
for old,new in renames.items():
    assert patched.count(old)==1
    patched=patched.replace(old,new)
singular.eval("kill extdevelop;")
singular.eval(patched)


def reduced_pencil(jroot):
    reduced=[]
    for xcoef in residual_cubic.list():
        wp=new_old_base(xcoef)
        row=[]
        for parameter_coefficient in wp.list():
            parameter_coefficient=new_parameter_ring(parameter_coefficient)
            row.append(tuple(reduce_L(c,rs,jroot) for c in parameter_coefficient.list()))
        reduced.append(tuple(row))
    return reduced


def equation_at(reduced,v):
    vv=F(v)
    eq=finite_plane(0)
    for xd,wcoeffs in enumerate(reduced):
        coeff=finite_plane(0)
        for wd,vcoeffs in enumerate(wcoeffs):
            # Horner in the new-base parameter.
            cv=F(0)
            for c in reversed(vcoeffs):
                cv=cv*vv+c
            if cv:
                coeff += cv*fw**wd
        eq += coeff*fx**xd
    return eq


def fiber_j(reduced,v,serial):
    eq=equation_at(reduced,v)
    if eq.degree(fw)!=9 or eq.degree(fx)!=3:
        raise ArithmeticError("degree drop")

    singular.set_ring(finite_plane._singular_())
    sf=eq._singular_()
    tag=f"Q80FASTMU{serial}"
    curve=tag+"C"; proj=tag+"P"; local=tag+"L"; place=tag+"PLACE"; ws=tag+"WS"

    t0=time.monotonic()
    singular.eval(
        "printlevel=-1; "
        f"list {curve}=Adj_div({sf.name()}); "
        f"def {proj}={curve}[1][2]; setring {proj};"
    )
    tadj=time.monotonic()-t0

    genus=int(singular.eval(f"{curve}[2][2]"))
    if genus!=1:
        raise ArithmeticError(f"genus={genus}")

    singular.eval(
        f"def {local}={curve}[5][1][1]; setring {local}; "
        f"int {tag}LI=size(POINTS)-1; setring {proj}; "
        f"int {place}=0; int {tag}I; "
        f"for ({tag}I=1;{tag}I<=size({curve}[3]);{tag}I={tag}I+1) "
        "{ "
        f"if ({curve}[3][{tag}I][1]==1 && {curve}[3][{tag}I][2]=={tag}LI) "
        f"{{ {place}={tag}I; }} "
        "}"
    )
    pn=int(singular.eval(place))
    if pn<=0:
        raise ArithmeticError("xi=-6 place not found")

    t0=time.monotonic()
    singular.eval(
        f"list {ws}=Weierstrass({pn},3,{curve}); "
        f"int {tag}S0={ws}[1][1]; int {tag}S1={ws}[1][2]; int {tag}S2={ws}[1][3]; "
        f"poly {tag}XN={ws}[2][2][1]; poly {tag}XD={ws}[2][2][2]; "
        f"poly {tag}YN={ws}[2][3][1]; poly {tag}YD={ws}[2][3][2]; "
        f"poly {tag}F=CHI;"
    )
    tw=time.monotonic()-t0

    if tuple(int(singular.eval(tag+x)) for x in ("S0","S1","S2"))!=(0,2,3):
        raise ArithmeticError("unexpected Weierstrass semigroup")

    t0=time.monotonic()
    xn=singular(tag+"XN").sage(); xd=singular(tag+"XD").sage()
    yn=singular(tag+"YN").sage(); yd=singular(tag+"YD").sage()
    peq=singular(tag+"F").sage()

    terms=(
        yn**2*xd**3,
        xn*yn*xd**2*yd,
        yn*xd**3*yd,
        -xn**3*yd**2,
        -xn**2*xd*yd**2,
        -xn*xd**2*yd**2,
        -xd**3*yd**2,
    )
    rem=tuple(t.reduce([peq]) for t in terms)
    mons=sorted(set().union(*(r.dict().keys() for r in rem)))
    M=Matrix(F,[[r.dict().get(m,F(0)) for r in rem] for m in mons])
    ker=M.right_kernel()
    if ker.dimension()!=1:
        raise ArithmeticError("relation kernel not 1-dimensional")
    r=ker.basis()[0]
    if r[0]==0 or r[3]==0:
        raise ArithmeticError("bad Weierstrass relation")
    r/=r[0]

    c=r[3]
    a1=r[1]; a2=r[4]; a3=r[2]*c; a4=r[5]*c; a6=r[6]*c**2
    b2=a1**2+4*a2
    b4=a1*a3+2*a4
    b6=a3**2+4*a6
    c4=b2**2-24*b4
    c6=-b2**3+36*b2*b4-216*b6
    delta=(c4**3-c6**2)/F(1728)
    if delta==0:
        raise ArithmeticError("singular recovered curve")
    tj=time.monotonic()-t0
    return c4**3/delta,tadj,tw,tj


def mu_for_embedding(jroot,sign_label):
    C=poly_from_exact(cs["C"],jroot,"V")
    S=poly_from_exact(cs["S"],jroot,"V")
    reduced=reduced_pencil(jroot)

    mus=[]
    serial=0
    for v in list(range(1,p,2))+list(range(2,p,2)):
        if len(mus)>=args.checks:
            break
        serial+=1
        try:
            jv,tadj,tw,tr=fiber_j(reduced,v,serial)
            cv=C(F(v)); sv=S(F(v))
            if jv==0 or sv==0 or cv==0:
                continue
            mu=cv**3*(jv-F(1728))/(jv*sv**2)
            mus.append(mu)
            print(
                f"Q80MUFAST|p={p}|j_sign={sign_label}|V={v}|jV={int(jv)}|"
                f"mu={int(mu)}|adj={tadj:.3f}|weierstrass={tw:.3f}|"
                f"relation={tr:.3f}|good={len(mus)}",
                flush=True,
            )
        except Exception as exc:
            print(
                f"Q80MUFAST|p={p}|j_sign={sign_label}|V={v}|status=SKIP|"
                f"type={type(exc).__name__}|message={str(exc).replace(chr(10),' ')}",
                flush=True,
            )

    if len(mus)<args.checks:
        raise ArithmeticError(f"only {len(mus)} usable mu fibers for j_sign={sign_label}")
    if any(x!=mus[0] for x in mus[1:]):
        raise ArithmeticError(f"mu disagreement for j_sign={sign_label}: {mus}")
    return mus[0]


started=time.monotonic()
up=mu_for_embedding(rj,+1)
um=mu_for_embedding(-rj,-1)
a=(up+um)/2
b=(up-um)/(2*rj)

record={
    "version":1,
    "prime":p,
    "s_root_used":int(rs),
    "j_positive_root":int(rj),
    "checks_per_embedding":args.checks,
    "mu_plus_j":int(up),
    "mu_minus_j":int(um),
    "mu_a_mod_p":int(a),
    "mu_j_mod_p":int(b),
    "seconds":time.monotonic()-started,
}
args.out.parent.mkdir(parents=True,exist_ok=True)
args.out.write_text(json.dumps(record,indent=2,sort_keys=True,default=int)+"\n")
print(
    f"Q80MUFAST|p={p}|a={int(a)}|b={int(b)}|"
    f"fibers={2*args.checks}|seconds={record['seconds']:.3f}|"
    f"status=PASS_FAST_MU|out={args.out}",
    flush=True,
)
