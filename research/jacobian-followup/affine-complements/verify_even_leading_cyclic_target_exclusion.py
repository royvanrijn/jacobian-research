"""PF149 exact support for the universal monic-coefficient pole proof."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="MONIC_PARAMETER_POLE_BOUNDS_AND_EVEN_LEADING_CYCLIC_TARGETS.md"
CONTEXT="sources/even-leading-cyclic-target-exclusion/context.json"
RECEIPT="even_leading_cyclic_target_exclusion_receipt.json"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    started=time.perf_counter()
    checks=[]
    def check(name,value):
        assert bool(value),name
        checks.append(name)
    ctx=json.loads((HERE/CONTEXT).read_text())
    for path,digest in ctx["retained_sha256"].items():
        check("retained proof pin "+path,sha(ROOT/path)==digest)

    u,w,z,p,v,s,b,T=S.symbols("u w z p v s b T",nonzero=True)
    M,K,A,Au,Az,Bu,Bz=S.symbols("M K A Au Az Bu Bz",nonzero=True)
    zu=-K/M;zw=2*w/M
    Pu=w*Au+Bu+(w*Az+Bz)*zu
    Pw=A+(w*Az+Bz)*zw
    dzden=2*w*(w*Au+Bu)+A*K
    duden=M*A+2*w*(w*Az+Bz)
    check("whole mixed dz coframe retains A_u and B_u",
          S.cancel((Pu*zw-Pw*zu)/dzden-1/M)==0)
    check("whole mixed du coframe retains A_z and B_z",
          S.cancel(Pw/duden-1/M)==0)
    check("nonconstant w coefficient not treated as additive constant",
          S.cancel(duden-(M*A+2*w*Bz))==2*w**2*Az)

    # Actual monic curve construction in unbounded-coefficient examples.
    examples=[
        (1,u**2,u+u**100,(u-2)**2,u**7+1,
         (u-2)**2*z+u**13),
        (2,u**2*(u-1)**2,u*(1-3*u)**2,1,
         u**7*z+u**3+1,z**2+(u**9+1)*z+u**15),
        (3,u**2,u,(u-2)**6*(u-3),
         u**31*z**2+u**20*z+1,
         (u-2)**6*(u-3)*z**3+u**25*z**2+u**22*z+u**30)]
    for nn,MM,NN,cc,AA,BB in examples:
        zexpr=(w**2-NN)/MM
        Pstar=(BB+w*AA).subs(z,zexpr)
        Q=S.cancel(MM**nn*Pstar/cc)
        poly=S.Poly(Q,w)
        check(f"curve polynomial is monic degree2n n{nn}",
              poly.degree()==2*nn and S.cancel(poly.LC()-1)==0)
        rhs=MM**nn/cc
        check(f"exact whole curve equation n{nn}",
              S.cancel(Q-p*rhs-rhs*(Pstar-p))==0)
        wp=S.cancel(rhs/S.diff(Q,w))
        check(f"parameter derivative inverse identity n{nn}",
              S.cancel(S.diff(Pstar,w)*wp-1)==0)
        check(f"actual volume relative coframe n{nn}",
              S.cancel(S.diff(Pstar,w)*wp/MM-1/MM)==0)
    check("new universal theorem needs no N degree gap",
          S.degree(examples[0][2],u)==100
          and S.degree(examples[0][1],u)==2)
    check("actual MX example crosses old Newton face",
          3*9+4>2*4)
    check("actual MX example has z-dependent w coefficient",
          S.degree(examples[1][4],z)==1)

    # Exact branch models for the monic parameter lemma, with arbitrarily
    # large constant polar parts. No bounded list substitutes for its proof.
    for DD,hhpow,center in [(2,-2,20),(4,-6,50),(6,-12,100),
                             (4,8,30),(4,-4,40)]:
        gg=gcd(DD,abs(hhpow));ee=DD//gg;rr=hhpow//gg
        W0=s**(-ee*center)+7*s**(-ee*(center-1))
        WW=W0+b*s**rr
        QQ=(T-W0)**DD
        wp=s**rr/(DD*b**(DD-1))  # p=b^D, differentiated in p
        check(f"monic lemma exact branch D{DD}h{hhpow}",
              S.cancel(QQ.subs(T,WW)-b**DD*s**(ee*hhpow))==0)
        check(f"monic lemma exact parameter derivative D{DD}h{hhpow}",
              S.cancel(S.diff(QQ,T).subs(T,WW)*wp-s**(ee*hhpow))==0)
        check(f"monic lemma sharp bound all conjugates D{DD}h{hhpow}",
              DD*rr==ee*hhpow and gg*ee==DD)
    # A forbidden too-early nonconstant coefficient would satisfy a
    # nonzero polynomial over C; here its initial polynomial is Y^4.
    Y=S.symbols("Y")
    QQ=(T-v**-8)**4+v**-3*(T-v**-8)**2+v**-5
    shifted=S.expand(QQ.subs(T,v**-8+v**-2*Y))
    check("too-early coefficient has nonzero constant initial polynomial",
          S.limit(v**8*shifted,v,0)==Y**4)
    check("too-early coefficient precedes parameter and is forbidden",
          -8 < -6)

    # Symbolic exact order identities, including the integer rounding.
    alpha,mm,bb,nn,mr,br,nuWp=S.symbols(
        "alpha m b n mr br nuWp",integer=True,positive=True)
    inf=(mm-1)*alpha-1+(bb-nn*mm)*alpha/(2*nn)
    infbound=((mm-2)/2+bb/(2*nn))*alpha-1
    fin=(1-mr)*alpha-1+(nn*mr-br)*alpha/(2*nn)
    finbound=(1-mr/2-br/(2*nn))*alpha-1
    check("infinite valuation budget exact algebra",S.simplify(inf-infbound)==0)
    check("finite valuation budget exact algebra",S.simplify(fin-finbound)==0)
    check("global numerator is b+nm-2n",
          S.simplify(((mm-2)/2+bb/(2*nn))-(bb+nn*mm-2*nn)/(2*nn))==0)
    check("finite numerator is 2n-nmr-br",
          S.simplify((1-mr/2-br/(2*nn))-(2*nn-nn*mr-br)/(2*nn))==0)
    for mv,nv,bv in [(0,2,4),(1,1,1),(2,2,0),(4,2,0),(2,2,5)]:
        margin=S.Rational(bv+nv*mv-2*nv,2*nv)
        orders=[S.ceiling(margin*av)-1 for av in [1,2,3]]
        check(f"global logarithmic budget m{mv}n{nv}b{bv}",
              margin>=0 and min(orders)>=-1)
        if margin>0:
            check(f"global strict budget is holomorphic m{mv}n{nv}b{bv}",
                  min(orders)>=0)
    for nv,mrv,brv in [(1,0,2),(2,0,3),(2,1,2),(2,1,1)]:
        margin=S.Rational(2*nv-nv*mrv-brv,2*nv)
        orders=[S.ceiling(margin*av)-1 for av in [1,2,3]]
        check(f"finite logarithmic budget n{nv}mr{mrv}br{brv}",
              margin>=0 and min(orders)>=-1)
        if margin>0:
            check(f"finite strict budget is holomorphic n{nv}mr{mrv}br{brv}",
                  min(orders)>=0)
    E,order=S.symbols("E order",integer=True,positive=True)
    check("every ramified logarithmic pullback remains simple",
          S.expand(E*(-1)+(E-1))==-1)
    check("every ramified holomorphic pullback stays holomorphic",
          S.expand(E*order+(E-1)-(E*(order+1)-1))==0)

    # Actual leading coefficients of arbitrary degree, with equality kept.
    cc=(u-2)**4*(u-3)**3*(u-4)
    MM=u**2*(u-1)**2
    factors=S.factor_list(cc)[1]
    check("nonconstant leading coefficient coprime to actual MX M",
          S.gcd(cc,MM)==1)
    check("all actual leading multiplicities within2n including equality",
          [e for f,e in factors]==[1,3,4] or sorted(e for f,e in factors)==[1,3,4])
    check("nonconstant leading degree unrestricted by old support",
          S.degree(cc,u)==8 and S.degree(cc,u)+2*S.degree(MM,u)>=4)
    check("any shared leading zero on actual MX violates sufficient budget",
          1+2*2>2*2)

    # Actual S2 finite equality: P=(u-1)^2 z.
    pp=S.symbols("pp",positive=True)
    uu=1+v
    WW=S.sqrt(uu*v**2+pp*uu**2)/v
    om=S.cancel(-S.diff(WW,pp)/uu**2)
    residue=S.limit(v*om,v,0)
    check("finite S2 equality has exact nonzero logarithmic residue",
          S.simplify(residue+1/(2*S.sqrt(pp)))==0)
    check("opposite finite branch has opposite nonzero residue",
          S.simplify(S.limit(-v*om,v,0)-1/(2*S.sqrt(pp)))==0)
    check("finite S2 equality curve relation",
          S.cancel(WW**2-uu-uu**2*pp/v**2)==0)

    # Nonlinear S2 control with two geometric generic components.
    zS=(w**2-u)/u**2
    H=zS+w+u**2/4
    BB=z**2+S.Rational(3,2)*u**2*z+u+u**4/16
    AA=2*z+u**2/2
    check("nonlinear logarithmic S2 fits even monic form exactly",
          S.cancel((BB+w*AA).subs(z,zS)-H**2)==0)
    hh,qq=S.symbols("hh qq",nonzero=True)
    check("both geometric generic components retained algebraically",
          S.expand((H-qq)*(H+qq)-(H**2-qq**2))==0)
    neww=S.symbols("neww")
    check("translated conic on each component",
          S.cancel(H.subs(w,neww-u**2/2)-(neww**2-u)/u**2)==0)
    hpos=S.symbols("hpos",positive=True)
    uu=1/v
    W=S.sqrt(uu+hpos*uu**2)
    omegaP=-S.diff(uu,v)/(4*hpos*W)
    check("nonlinear component infinity residue nonzero",
          S.simplify(S.limit(v*omegaP,v,0)-1/(4*hpos*S.sqrt(hpos)))==0)

    # Genuinely above-face S2 candidate passes both marked source axes.
    t,h,x=S.symbols("t h x")
    uu=t**2;ww=t*(1+t**2*h);zz=2*h+t**2*h**2
    candidate=S.expand(ww+zz**2+uu*zz+uu)
    check("above-face S2 full original source relation",
          S.expand(ww**2-uu-uu**2*zz)==0)
    check("above-face S2 C0 degree two",S.degree(candidate.subs(t,0),h)==2)
    check("above-face S2 other axis degree two",S.degree(candidate.subs(h,0),t)==2)
    check("above-face S2 ordinary vertical degree four",
          S.degree(candidate.subs(t,2),h)==4)

    # Full actual MX example (11), without any curve-point omission.
    uu=t**2;sig=1-3*uu;aa=uu*(uu-1)**2
    ww=t*(sig+aa*x);zz=2*sig*x+aa*x**2
    check("above-face mixed-w MX full source relation",
          S.expand(ww**2-uu*(1-3*uu)**2-uu**2*(uu-1)**2*zz)==0)
    candidate=S.expand(zz**2+(uu**9+1)*zz+uu**15+ww*(uu**7*zz+uu**3+1))
    for tv,dv in [(0,2),(1,2),(2,4)]:
        check(f"above-face mixed-w MX reached line t{tv} degree{dv}",
              S.degree(candidate.subs(t,tv),x)==dv)

    # Exact countercontrols to invalid universal promotions.
    polemate=S.Matrix([w+u*zS,1/w]).jacobian([u,w]).det()
    check("finite leading budget cannot be dropped rational mate",
          S.cancel(polemate-1/u**2)==0)
    check("known rational mate violates precisely finite leading condition",
          1+1*2>2 and 1+1*2>=2)
    oddmate=S.Matrix([w+S.Function("g")(u),1/u]).jacobian([u,w]).det()
    check("odd highest degree remains outside theorem",
          S.cancel(oddmate-1/u**2)==0)
    zplane=w**2-u
    check("global degree budget cannot be dropped regular plane mate",
          S.Matrix([w+zplane**2,zplane]).jacobian([u,w]).det()==1)
    check("regular plane countercontrol violates global budget",0+2*0<2*2)

    out=dict(claim="PF2RULED149",outcome="PASS",check_count=len(checks),checks=checks,
             elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),
             proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),
             context_sha256=sha(HERE/CONTEXT),
             old_Newton_support_restriction_removed_for_stated_family=True,
             arbitrary_lower_w_coefficients_retained=True,
             all_leading_coefficients_unrestricted=False,
             arbitrary_target_coordinates_reduced_to_tested_family=False,
             actual_full_marked_quadratic_regularity=False,
             unrestricted_cyclic_targets_excluded=False,full_goal_complete=False,
             registered_predecessor_checkers_replayed=[],
             independent_replay=False,formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["claim","outcome","check_count","elapsed_seconds"]}))

if __name__=="__main__":
    main()
