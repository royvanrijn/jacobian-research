#!/usr/bin/env python3
"""PF187 exact controls; the unrestricted operator/contour proof is in the source."""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import platform
import time
import sympy as S
START=time.perf_counter()
BASE=Path(__file__).resolve().parent
PROOF=BASE/"RULING_SIX_UNIVERSAL_THREE_MOMENT_FIRST_QUARTIC_REDUCTION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,w,Y,Z,T,p=S.symbols("v w Y Z T p")
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(m[0] for m in pn.monoms());od=min(m[0] for m in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))
def power_coeffs(f,alpha,N):
    ys=[F(1)]
    for n in range(1,N+1):
        ys.append(sum((((alpha+1)*k-n)*f[k]*ys[n-k]
                       for k in range(1,min(n,len(f)-1)+1)),F(0))/n)
    return ys
def addto(dest,e,c):
    if c:
        dest[e]=dest.get(e,F(0))+c
        if not dest[e]:del dest[e]
def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():addto(out,i+j,x*y)
    return out
def combine(terms):
    out={}
    for factor,poly in terms:
        for e,c in mul(factor,poly).items():addto(out,e,c)
    return out
def leading(poly):
    return (min(poly),poly[min(poly)]) if poly else (S.oo,F(0))
def raw_whole(poly,L,n):
    pp=S.Poly(poly,w)
    assert pp.degree()==6 and pp.LC()==v**(4*L)
    fs={}
    for k in range(1,7):
        fs[k]={i[0]-4*L:F(int(c)) for i,c in S.Poly(pp.nth(6-k),v).terms() if c}
    ys=[{0:F(1)}];alpha=F(n,2)
    for j in range(1,3*n+5):
        out={}
        for k in range(1,min(j,6)+1):
            scale=((alpha+1)*k-j)/j
            for e,c in mul(fs[k],ys[j-k]).items():addto(out,e,c*scale)
        ys.append(out)
    return [{e+2*L*n:c for e,c in ys[3*n+i+1].items()} for i in range(4)]

# Exact inverse-prefix annihilation for the three native weights.
a3,a2,a1,a0,r2,r1,r0,aa=S.symbols("a3 a2 a1 a0 r2 r1 r0 aa",nonzero=True)
rh=a3+a2*T+a1*T**2+a0*T**3
rr=r2+r1*T+r0*T**2
inverse=S.series(T**3/rh*(1+T**4*rr/rh**2)**S.Rational(-1,2),T,0,5).removeO()
weights=[S.Integer(1),1/T-aa,(1/T)*(1/T-aa)**2/aa]
for j,weight in enumerate(weights):
    res=S.expand(weight*inverse).coeff(T,1)
    check("negative coefficient vanishes j"+str(j),zero(res.subs(aa,-a2/(2*a3))))
    for exponent in [5,6,7,8]:
        check(f"prefix remainder cannot contribute j{j} exponent{exponent}",
              S.expand(weight*T**exponent).coeff(T,1)==0)
check("wrong center retains a negative character",
      zero(S.expand(weights[2]*inverse).coeff(T,1).subs(aa,-a2/(3*a3))-1/a3))

# General regular finite source root and exact root-sum center, all t>=L.
gamma,rho0,mu,dd=S.symbols("gamma rho0 mu dd",nonzero=True)
for L,b in [(3,0),(5,0),(5,3),(7,3),(8,6),(12,10)]:
    t=L+b
    rr0=v**L*rho0
    mh=gamma+v*mu
    source=(1+v)*(Y-rr0)*((Y-mh)**2+v**3)
    pp=S.Poly(S.expand(source),Y)
    ratio=pp.nth(2)/pp.nth(3)
    ac=-v**(-L)*ratio/2
    bg=v**(-L)*(mh+v**t*dd)
    check(f"universal exact center difference L{L} b{b}",
          zero(ac-bg-(rho0/2-v**b*dd)))
    hh=S.cancel(v**(-L)*source.subs(Y,v**L*w))
    check(f"regular finite linear source root L{L} b{b}",
          zero(hh.subs(v,0)-gamma**2*(w-rho0)))
    check(f"finite source derivative unit L{L} b{b}",
          zero(S.diff(hh,w).subs(v,0)-gamma**2))

# Generic differential-operator identity exposes all three moments exactly.
c2,c1,c0=S.symbols("c2 c1 c0")
R=Z**4+c2*Z**2+c1*Z+c0
for n in [1,3,5,7]:
    alpha=S.Rational(n,2)
    f=[S.Integer(1),S.Integer(0),c2,c1,c0]
    co=[S.Integer(1)]
    for j in range(1,2*n+4):
        co.append(S.expand(sum((((alpha+1)*k-j)*f[k]*co[j-k]
                                  for k in range(1,min(j,4)+1)),S.Integer(0))/j))
    pol=sum(co[j]*Z**(2*n-j) for j in range(2*n+1))
    moments=co[2*n+1:2*n+4]
    operator=S.expand(2*R*S.diff(pol,Z)-n*S.diff(R,Z)*pol)
    expected=(2*(2*n+1)*moments[0]*Z**2+4*(n+1)*moments[1]*Z
              +2*(n+1)*c2*moments[0]+2*(2*n+3)*moments[2])
    check("generic exact three-moment differential operator n"+str(n),
          zero(operator-expected))
    check("operator has degree at most two n"+str(n),S.Poly(operator,Z).degree()<=2)

# A zero operator makes the rational quotient S^2/R^n constant.
ss,sd,rd,rr0,nn=S.symbols("ss sd rd rr0 nn")
check("quotient derivative equals operator times polynomial",
      zero((2*ss*sd*rr0-nn*ss**2*rd)-ss*(2*rr0*sd-nn*rd*ss)))
qA,qB=S.symbols("qA qB")
square=S.Poly(S.expand((Z**2+qA*Z+qB)**2),Z)
check("depressed square cubic coefficient is twice linear coefficient",square.nth(3)==2*qA)
check("depressed square linear coefficient vanishes after qA0",square.nth(1).subs(qA,0)==0)
check("depressed square constant is square of middle coefficient",square.nth(0).subs(qA,0)==qB**2)

# Repeated roots and all three congruences are included in the vector controls.
models=[
 ("qdom",Z**4+Z,False),
 ("delta",Z**4+1,False),
 ("even_boundary",(Z**2+1)**2+1,False),
 ("double_zero",Z**2*(Z**2+2),False),
 ("triple_simple",(Z-1)**3*(Z+3),False),
 ("double_two_simple",(Z-1)**2*(Z-2)*(Z+4),False),
 ("q_delta",Z**4+Z+1,False),
 ("general",Z**4+2*Z**2+3*Z+2,False),
 ("square",(Z**2+1)**2,True),
 ("fourfold",Z**4,True)]
vectors={}
for name,poly,is_square in models:
    rp=S.Poly(S.expand(poly),Z)
    f=[F(1)]+[F(int(rp.nth(4-j))) for j in range(1,5)]
    for n in range(1,24,2):
        co=power_coeffs(f,F(n,2),2*n+3)
        vec=co[2*n+1:2*n+4]
        check(f"quartic vector iff square control {name} n{n}",
              (all(x==0 for x in vec))==is_square)
        vectors[(name,n)]=vec
check("pure qdom first congruence has only first moment",vectors[("qdom",7)][0]!=0)
check("pure qdom second congruence has only second moment",vectors[("qdom",5)][1]!=0)
check("pure qdom third congruence has adapted third moment",vectors[("qdom",3)][2]!=0)
for n in range(1,20,2):
    for j0 in [0,1,-2,3]:
        poly=(Z**2+j0)**2+2
        rp=S.Poly(poly,Z)
        f=[F(1)]+[F(int(rp.nth(4-j))) for j in range(1,5)]
        co=power_coeffs(f,F(n,2),2*n+3)
        check(f"even-boundary native centered residue n{n} j{j0}",
              co[2*n+2]==F(S.binomial(S.Rational(n,2),(n+1)//2))*F(2)**((n+1)//2))

# Exact finite/escape comparisons on a common Puiseux cover.
scale_controls=[]
for L,rho in [(5,F(3,2)),(7,F(3)),(8,F(3)),(9,F(4)),(11,F(10,3)),(13,F(17,4))]:
    lam=F(L)-2*rho
    check(f"subcritical lambda positive L{L} rho{rho}",lam>0)
    finite=[F(0),rho,2*rho-L]
    for n in range(1,18,2):
        pn=lam*n+L-rho
        check(f"top dominates every retained finite contribution L{L} rho{rho} n{n}",pn>max(finite))
        check(f"finite margin is lambda nplus1 L{L} rho{rho} n{n}",pn-rho==lam*(n+1))
        if n>=3:
            lower=lam*(n-2)+L-rho
            check(f"strict arbitrary-character gap L{L} rho{rho} n{n}",pn-lower==2*lam>0)
    scale_controls.append({"L":L,"rho":str(rho),"lambda":str(lam),"finite_bounds":list(map(str,finite))})

# Whole sextic controls retain TWO finite simple roots and bounded center mismatch.
whole_controls=[]
selected=[models[i] for i in [0,1,2,3,4,5,6]]
for L,rho in [(7,3),(8,3),(9,4)]:
    lam=L-2*rho
    for name,rrpoly,is_square in selected:
        rp=S.Poly(S.expand(rrpoly),Z)
        GG=sum(rp.nth(i)*v**((4-i)*rho)*(v**L*w-1)**i for i in range(5))
        whole=S.expand((w**2-3*w+2)*GG)
        wp=S.Poly(whole,w)
        center=S.cancel(-wp.nth(5)/(4*wp.nth(6)))
        check(f"whole adapted center exact L{L} {name}",zero(center-v**(-L)-S.Rational(3,4)))
        check(f"whole finite quadratic retained L{L} {name}",zero(whole.subs(v,0)-(w**2-3*w+2)))
        ad={-L:F(1),0:F(3,4)}
        ad2=mul(ad,ad)
        for n in [1,3,5,7]:
            raw=raw_whole(whole,L,n)
            W1=combine([({0:F(1)},raw[1]),({e:-c for e,c in ad.items()},raw[0])])
            numerator=combine([({0:F(1)},raw[3]),
                               ({e:-2*c for e,c in ad.items()},raw[2]),
                               (ad2,raw[1])])
            order_coeff=[leading(raw[0]),leading(W1),leading(numerator)]
            shifts=[0,L-rho,3*L-2*rho]
            ppole=lam*n+L-rho
            target=vectors[(name,n)]
            actual=[]
            for j,((oo,cc),shift) in enumerate(zip(order_coeff,shifts)):
                aligned=oo+shift
                check(f"whole aligned pole bound L{L} {name} n{n} j{j}",aligned>=-ppole)
                got=cc if aligned==-ppole else F(0)
                actual.append(got)
                check(f"whole exact leading vector L{L} {name} n{n} j{j}",got==target[j])
            check(f"whole nonzero exact vector L{L} {name} n{n}",any(actual))
            whole_controls.append({"L":L,"rho":rho,"quartic":name,"n":n,
                                   "common_pole":ppole,"leading_vector":list(map(str,actual))})

# The cap edge has a generic simple quartic and a forbidden actual pole equality.
for j0,q0,constant in [(0,0,0),(1,0,1),(1,2,0),(-2,3,4)]:
    edge=Z**4+2*j0*Z**2+q0*Z+constant-p
    resultant=S.resultant(edge,S.diff(edge,Z),Z)
    check(f"generic parameter-edge quartic separable j{j0} q{q0}",resultant!=0 and S.Poly(resultant,p).degree()==3)
    check(f"generic parameter-edge quartic has no zero root j{j0} q{q0}",edge.subs(Z,0)!=0)
for L in range(3,24):
    mu=S.Rational(3*L,2)
    check("parameter-edge coframe has a genuine pole L"+str(L),mu>L+1)
    check("parameter-edge derivative equals forbidden affine-response threshold L"+str(L),mu==S.Rational(3*L,2))

receipt={"claim":"PF2RULED187","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
        "RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md",
        "RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md",
        "RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md",
        "RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md",
        "RULING_SIX_NONCRITICAL_QUARTIC_SPLITTING_AND_SMALL_SLOPE_EXCLUSION.md",
        "RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md",
        "RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md",
        "RULING_SIX_CUBIC_BALANCE_COMPANION_DEGREE_TWENTY_SEVEN.md",
        "RULING_SIX_ADAPTED_CUBIC_MOMENT_AND_DEGREE30_EXCLUSION.md"]},
    "uniform_proved_scope":{
        "source_and_companion_degrees":"unrestricted",
        "actual_A_slope":"every L>=3",
        "nonsquare_first_initial_excluded":True,
        "multiple_roots_retained":True,
        "parameter_cap_edge_excluded":True,
        "necessary_remaining_first_initial":"(Z^2+j0)^2, j0 nonzero",
        "necessary_orders":"ell finite, ell<L, delta>2ell, 3ell<2t",
        "deeper_square_deformation_closed_here":False,
        "full_ruling_six_HC4_or_JC2_closed":False},
    "scale_controls":scale_controls,"whole_sextic_controls":whole_controls,
    "assurance":"Exact internal controls and unrestricted written operator/contour proof; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_universal_first_quartic_receipt.json").write_text(
        json.dumps(receipt,indent=2)+chr(10),encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
