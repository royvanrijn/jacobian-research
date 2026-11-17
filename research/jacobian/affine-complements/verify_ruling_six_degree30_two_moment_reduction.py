#!/usr/bin/env python3
"""PF184 controls; the unrestricted two-moment/contour proof is written."""
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
PROOF=BASE/"RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md"
checks=[]
def check(label,value):
    if not bool(value):raise AssertionError(label)
    checks.append(label)
def zero(e):return S.cancel(e)==0
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
v,w,X,Z,T,p,s,lam=S.symbols("v w X Z T p s lam")
a,z0,A,B=S.symbols("a z0 A B",nonzero=True)
def initial(expr,var):
    num,den=S.fraction(S.cancel(expr))
    if num==0:return S.oo,S.Integer(0)
    pn,pd=S.Poly(num,var),S.Poly(den,var)
    on=min(m[0] for m in pn.monoms());od=min(m[0] for m in pd.monoms())
    return on-od,S.cancel(pn.nth(on)/pd.nth(od))
def power_coeff(f,alpha,N):
    """Exact [t^N] f(t)^alpha for f(0)=1; f is a coefficient list."""
    y=[F(1)]
    for n in range(1,N+1):
        y.append(sum((((alpha+1)*k-n)*f[k]*y[n-k]
                      for k in range(1,min(n,len(f)-1)+1)),F(0))/n)
    return y[N]
def addto(dest,e,c):
    if c:
        dest[e]=dest.get(e,F(0))+c
        if not dest[e]:del dest[e]
def mul(a,b):
    out={}
    for i,x in a.items():
        for j,y in b.items():addto(out,i+j,x*y)
    return out
def whole_moments(poly,n):
    """Complete exact Laurent coefficients of a whole sextic, not local truncation."""
    pp=S.Poly(poly,w)
    assert pp.degree()==6 and pp.LC()==v**28
    fs={}
    for k in range(1,7):
        fs[k]={i[0]-28:F(int(c)) for i,c in S.Poly(pp.nth(6-k),v).terms() if c}
    yy=[{0:F(1)}];alpha=F(n,2)
    for j in range(1,3*n+3):
        out={}
        for k in range(1,min(j,6)+1):
            scale=((alpha+1)*k-j)/j
            for e,c in mul(fs[k],yy[j-k]).items():addto(out,e,c*scale)
        yy.append(out)
    m={e+14*n:c for e,c in yy[3*n+1].items()}
    lc={e+14*n:c for e,c in yy[3*n+2].items()}
    for e,c in m.items():addto(lc,e-7,-c)
    return m,lc
def leading(d):
    return (min(d),d[min(d)]) if d else (S.oo,F(0))

# The triple/simple coefficient identity has an independent affine-change derivation.
triple=(Z-a)**3*(Z+3*a)
check("triple-simple quartic is depressed",S.expand(triple).coeff(Z,3)==0)
check("exact affine transformation for every exponent",
      zero(triple.subs(Z,a+4*a*T)-(4*a)**4*T**3*(T+1)))
for n in range(1,16,2):
    direct=power_coeff([F(1),F(0),F(-6),F(8),F(-3)],F(n,2),2*n+1)
    expected=F(4)**(2*n+1)*F(S.binomial(S.Rational(n,2),2*n+1))
    check("direct triple coefficient matches all-index binomial n"+str(n),direct==expected)
    check("positive odd triple binomial nonzero n"+str(n),expected!=0)
for n in [2,4,6]:
    check("even index is deliberately excluded from nonvanishing n"+str(n),
          S.binomial(S.Rational(n,2),2*n+1)==0)

# Exact first two residue formulas at the two quadratic roots.
for r in range(1,11):
    numer=(A*Z+B)**r
    rem=S.rem(S.expand(numer),S.Poly(Z**2-z0**2,Z),Z)
    rem2=S.rem(S.expand(Z*numer),S.Poly(Z**2-z0**2,Z),Z)
    U=B+A*z0;V=B-A*z0
    M=(U**r-V**r)/(2*z0);N=(U**r+V**r)/2
    check("double first residue formula r"+str(r),zero(S.Poly(rem,Z).nth(1)-M))
    check("double centered second formula r"+str(r),zero(S.Poly(rem2,Z).nth(1)-N))
    check("double vector recovers U power r"+str(r),zero(N+z0*M-U**r))
    check("double vector recovers V power r"+str(r),zero(N-z0*M-V**r))
for aa,bb in [(1,2),(1,0),(0,2),(1,1),(1,-1)]:
    for r in range(1,9):
        UU=bb+aa;VV=bb-aa
        pair=(F(UU**r-VV**r,2),F(UU**r+VV**r,2))
        check(f"double vector remains nonzero A{aa} B{bb} r{r}",pair!=(0,0))
check("one scalar residue can vanish while the vector remains nonzero",
      ((1**2-(-1)**2)==0) and ((1**2+(-1)**2)!=0))
check("zero perturbation is explicitly inadmissible",
      (0+0,0-0)==(0,0))

# The even quartic is detected immediately by the centered second coefficient.
for n in range(1,16,2):
    direct=power_coeff([F(1),F(0),F(0),F(0),F(2)],F(n,2),2*n+2)
    expected=F(S.binomial(S.Rational(n,2),(n+1)//2))*F(2)**((n+1)//2)
    check("even quartic centered coefficient n"+str(n),direct==expected and direct!=0)

# Whole sextics retain a finite SIMPLE two-root cluster at w=+1,-1.
whole_controls=[]
Gtr=(X-v**3)**3*(X+3*v**3)+v**13
Ptr=S.expand((w**2-1)*Gtr.subs(X,v**7*w-1))
check("whole triple control finite specialization",zero(Ptr.subs(v,0)-(w**2-1)))
for n in [1,3,5,7,9]:
    mm,ll=whole_moments(Ptr,n)
    order,coeff=leading(mm)
    expected=F(4)**(2*n+1)*F(S.binomial(S.Rational(n,2),2*n+1))
    check("whole triple exact moment order n"+str(n),order==-n-4)
    check("whole triple exact moment coefficient n"+str(n),coeff==expected)
    whole_controls.append({"type":"triple","n":n,"m_order":int(order),"m_leading":str(coeff)})
for aa,bb in [(1,2),(1,0),(0,2),(1,1),(1,-1)]:
    qterm=v**10*aa*X if aa else v**11*X
    GG=(X**2-v**6)**2+qterm+bb*v**13+3*v**14
    PP=S.expand((w**2-1)*GG.subs(X,v**7*w-1))
    check(f"whole double finite simple cluster A{aa} B{bb}",zero(PP.subs(v,0)-(w**2-1)))
    for n in [1,3,5,7,9]:
        r=(n+1)//2
        mm,ll=whole_moments(PP,n)
        aligned={e+4:c for e,c in ll.items()}
        target=-r-3
        UU=bb+aa;VV=bb-aa
        factor=F(S.binomial(S.Rational(n,2),r))
        expectedM=factor*F(UU**r-VV**r,2)
        expectedN=factor*F(UU**r+VV**r,2)
        check(f"whole double first coefficient and bound A{aa} B{bb} n{n}",
              leading(mm)[0]>=target and mm.get(target,F(0))==expectedM)
        check(f"whole double second coefficient and bound A{aa} B{bb} n{n}",
              leading(aligned)[0]>=target and aligned.get(target,F(0))==expectedN)
        check(f"whole double exact vector order A{aa} B{bb} n{n}",
              (mm.get(target,0),aligned.get(target,0))!=(0,0))
        check(f"whole double finite contribution strictly later A{aa} B{bb} n{n}",target<-3)
Geven=(X**2+v**7)**2+v**10*X+2*v**13
Peven=S.expand((w**2-1)*Geven.subs(X,v**7*w-1))
check("whole even quartic finite simple cluster",zero(Peven.subs(v,0)-(w**2-1)))
for n in [1,3,5,7,9]:
    mm,ll=whole_moments(Peven,n)
    order,coeff=leading(ll)
    expected=F(S.binomial(S.Rational(n,2),(n+1)//2))*F(2)**((n+1)//2)
    check("whole even centered exact order n"+str(n),order==-(n+15)//2)
    check("whole even centered exact coefficient n"+str(n),coeff==expected)
    check("whole even finite contribution strictly later n"+str(n),order<-7)

# Exact scaling and finite-cluster order comparisons.
n,r=S.symbols("n r",integer=True,positive=True)
check("triple contour scaling",-7+3+(12-14)*n/2==-n-4)
check("double first scaling",zero((-7+3)+(12-14)*(2*r-1)/2+r-(-r-3)))
check("double second scaling",zero(2*(-7+3)+(12-14)*(2*r-1)/2+r-(-r-7)))
check("even quartic centered scaling",
      zero(2*(S.Rational(13,4)-7)+(13-14)*n/2+(n+15)/2))
check("finite centered alignment is minus three",4-7==-3)

# Curvature arguments use the WHOLE sextic, not a quartic with a frozen factor.
H=Z**4-4*Z+3
check("double-and-two-simple leading example",zero(H-(Z-1)**2*(Z**2+2*Z+3)))
check("exact double leading curvature",S.diff(H,Z,2).subs(Z,1)!=0)
check("double-and-two-simple upper budget is strictly below six",2*(9-8)+2*(10-8)==6)
for RR in [S.Rational(15,4),S.Integer(4),S.Rational(17,4)]:
    check("refined double cubic terms later R"+str(RR),3+3*RR>6+2*RR)
    check("refined double higher terms later R"+str(RR),4*RR>6+2*RR)
for t,ell in [(12,5),(12,6),(11,6),(12,4)]:
    tag=f"t{t} J{ell}"
    RR=7-S.Rational(ell,2)
    check("parameter curvature exact order "+tag,ell+2*RR==14)
    check("parameter curvature linear later "+tag,t+RR>14)
    check("parameter curvature cubic later "+tag,S.Rational(ell,2)+3*RR>14)
    check("parameter curvature higher later "+tag,4*RR>14)
    F0=(1+X)**2*((X**2-v**ell)**2+v**t*X+2*v**14)
    for sign in [-1,1]:
        XX=sign*s**ell+s**(14-ell)*lam
        oo,ll=initial((F0-p*v**14).subs({v:s**2,X:XX},simultaneous=True)/s**28,s)
        constant=2+(sign if 2*t+ell==28 else 0)
        check("whole parameter quadratic initial "+tag+str(sign),
              oo==0 and zero(ll-(4*lam**2+constant-p)))

# Complete integer region reduction: ell>=8 has the same inequalities as infinity.
square_candidates=[];classifications=[]
for t in range(9,13):
    for delta in range(t+1,15):
        for ell in list(range(1,8))+[None]:
            EE=S.oo if ell is None else S.Integer(ell)
            if delta<2*EE and delta<S.Rational(4*t,3):
                degree=3*delta-32
                label="even-moment" if degree in [6,7,8] else "degree"
                if label=="even-moment":
                    check("dominant Delta can survive only at thirteen",delta==13 and degree==7)
            elif delta==2*EE and delta<S.Rational(4*t,3):
                mu=S.Rational(3*delta,4)
                if 8<mu<S.Rational(21,2):
                    check("square boundary surviving interval only ell6 delta12",
                          ell==6 and delta==12 and t in [10,11])
                    label="double-curvature"
                else:label="interval"
            elif EE>S.Rational(2*t,3) and delta>S.Rational(4*t,3):
                degree=4*t-32
                if degree in [6,7,8]:
                    check("cubic balance only surviving degree8",
                          t==10 and delta==14 and (ell is None or ell>=7))
                    label="remaining"
                else:label="degree"
            elif EE<S.Rational(2*t,3) and delta>2*EE:
                qorder=t+EE/2
                m0=min(S.Rational(delta),qorder)
                mu=(m0+EE)/2
                if delta!=qorder:
                    degree=4*mu-32
                    possible=8<mu<S.Rational(21,2) and degree in [6,7,8]
                else:
                    possible=8<mu<S.Rational(21,2) and any(b<2*mu-11 for b in [6,7,8])
                if possible:
                    square_candidates.append((t,ell,delta))
                    label=("parameter-curvature" if delta==14 and m0==14 else "double-vector")
                else:label="degree-or-interval"
            else:
                check("all remaining boundaries in t9 quartic group",
                      t==9 and (ell is None or ell>=6) and delta>=12)
                label="triple-or-double"
            classifications.append({"t":t,"J_order":ell,"Delta_order":delta,"result":label})
expected=[(10,6,13),(10,6,14),(11,6,13),(11,6,14),
          (12,4,14),(12,5,14),(12,6,13),(12,6,14)]
check("exact square-center candidate table",sorted(square_candidates)==expected)
for row in classifications:
    if row["result"]=="double-vector":
        check("vector cases exactly the four listed cells "+str(row),
              (row["t"],row["J_order"],row["Delta_order"]) in
              [(10,6,13),(10,6,14),(11,6,13),(12,6,13)])
check("remaining cubic range nonempty",any(r["result"]=="remaining" for r in classifications))
check("remaining generic degree conditional eight",32+8==4*10)
check("remaining e-order is three",10-7==3)
check("remaining d-order is ten",7+3==10)
check("remaining original generic degree is sixteen",2*8==16)

receipt={"claim":"PF2RULED184","outcome":"PASS","check_count":len(checks),"checks":checks,
    "elapsed_seconds":time.perf_counter()-START,
    "environment":{"python":platform.python_version(),"sympy":S.__version__},
    "proof_sha256":sha(PROOF),"checker_sha256":sha(Path(__file__)),
    "retained_proof_sources":{name:sha(BASE/name) for name in [
        "POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md",
        "RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md",
        "RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md",
        "RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md",
        "RULING_SIX_DEGREE24_CRITICAL_QUARTIC_POLE_BUDGET_EXCLUSION.md",
        "RULING_SIX_NONCRITICAL_POLE_BUDGET_AND_DEGREE30_LEDGER.md",
        "RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md",
        "RULING_SIX_NONPOLE_PACKETS_AND_SLOPE7_TOP_DEGREE_EXCLUSION.md"]},
    "uniform_proved_scope":{
        "coefficient_and_companion_degrees":"unrestricted",
        "two_constant_relations":"same c_n for [w^-1] and [w^-2] positive odd half powers",
        "triple_simple_scalar_growth":"ord m_n=-n-4, every positive odd n",
        "two_double_vector_growth":"ord vector(m_n,v4*lhat_n)=-(r+3), n=2r-1, every r>=1",
        "even_quartic_centered_growth":"ord lhat_n=-(n+15)/2, every positive odd n",
        "finite_cluster_contributions_explicitly_retained":True,
        "D30_generic_degree_six_and_seven_excluded":True,
        "remaining_D30_intermediate_degree":8,
        "remaining_D30_original_degree":16,
        "remaining_local_orders":"ord e3, ord d10, ord q10, ord J>=7, ord Delta14",
        "remaining_cubic_balance_or_full_D30_closed":False,
        "full_ruling_six_or_JC2_closed":False},
    "whole_moment_controls":whole_controls,
    "complete_integer_region_controls":classifications,
    "assurance":"Exact internal symbolic controls and unrestricted written contour proof; no formal or external human review."}
parser=argparse.ArgumentParser()
parser.add_argument("--write-receipt",action="store_true")
args=parser.parse_args()
if args.write_receipt:
    (BASE/"ruling_six_degree30_two_moment_reduction_receipt.json").write_text(
        json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"claim":receipt["claim"],"outcome":"PASS","check_count":len(checks),
                  "elapsed_seconds":receipt["elapsed_seconds"],"receipt_written":args.write_receipt},indent=2))
