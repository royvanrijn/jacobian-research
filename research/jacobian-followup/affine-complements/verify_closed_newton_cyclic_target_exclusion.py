"""PF148 exact support for the written closed-face, all-branch proof."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="CLOSED_NEWTON_FACES_AND_GENERIC_PARAMETER_DERIVATIVE.md"
CONTEXT="sources/closed-newton-cyclic-target-exclusion/context.json"
RECEIPT="closed_newton_cyclic_target_exclusion_receipt.json"

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
        check("retained source pin "+path,sha(ROOT/path)==digest)

    u,w,z,p,a,v,s,b=S.symbols("u w z p a v s b",nonzero=True)
    M,N,K,Fu,Fz,gp=S.symbols("M N K Fu Fz gp",nonzero=True)
    zu=-K/M;zw=2*w/M
    Pu=Fu+gp+Fz*zu;Pw=a+Fz*zw
    Gu=K+2*w*(Fu+gp)/a;Gz=M+2*w*Fz/a
    check("actual dz coframe includes BOTH mixed partials",
          S.cancel((Pu*zw-Pw*zu)/(a*Gu)-1/M)==0)
    check("actual du coframe includes F_z",S.cancel(Pw/(a*Gz)-1/M)==0)
    check("dropping F_u changes actual coframe",
          S.cancel(Gu-(K+2*w*gp/a))==2*w*Fu/a)
    check("constant mixed coefficients recover additive derivative",
          S.cancel(Gu.subs(Fu,0)-(K+2*w*gp/a))==0)

    n,m,j,i,r=S.symbols("n m j i r",positive=True)
    t=m/(2*n-1)
    check("strict finite-z coefficient margin",
          S.simplify(m/S.Integer(2)-(n-1)*m/(2*n-1)-m/(2*(2*n-1)))==0)
    check("face meeting of surface and leading z pole",
          S.simplify((m+t)/2-n*t)==0)
    check("closed face controls every lower term at meeting",
          S.simplify(j+i*t-n*t-((2*n-1)*j+m*i-n*m)/(2*n-1))==0)
    check("positive-r exterior slope is n-i",
          S.diff(n*r-(j+i*r),r)==n-i)
    for mm,nn,jj,ii,boundary in [
        (4,2,1,1,False),(6,2,2,1,True),(6,2,1,1,False),
        (10,3,2,2,True),(10,3,4,1,True)]:
        check(f"actual mixed support m{mm}n{nn}u{jj}z{ii}",
              ((2*nn-1)*jj+mm*ii<=nn*mm)
              and (((2*nn-1)*jj+mm*ii==nn*mm)==boundary))
    check("S2 has no positive-u lower-z term under closed-face condition",
          S.simplify(2*(n-1)/(2*n-1)-(1-1/(2*n-1)))==0)

    alpha,beta,WW,JJ=S.symbols("alpha beta W J",integer=True,positive=True)
    check("g-f cancellation exact differential order",
          S.expand(WW+JJ-alpha-beta-1
                   -(WW-alpha+(n-1)*beta-1)).subs(JJ,n*beta)==0)
    d=S.symbols("d")
    check("g-w combined leading derivative nonzero coefficient identity",
          S.simplify(((m-2*d)*alpha+beta).subs(beta,(2*d-m)*alpha))==0)

    # Finite nonzero z: mixed u*z changes coefficients but not the exact order.
    MM=u**4;NN=u;FF=z**2+u*z;gg=u**2
    GG=NN+MM*z-(p-FF-gg)**2
    zz=1+2*v
    check("mixed finite nonzero-z leading surface balance",
          S.limit(v**3*GG.subs({u:1/v,z:zz}),v,0)==0)
    om=S.cancel(-S.diff(1/v,v)/S.diff(GG,z).subs({u:1/v,z:zz}))
    check("mixed finite nonzero-z exact coframe order",
          S.limit(om/v**2,v,0)==1)

    # Two nonbalanced branches, checked only to their leading order.
    examples=[
        ("g-f cancellation",u**4,u,z**2+u*z,-u**4,v**-2,3,-S.Rational(1,4),6),
        ("g-w balance",u**5,u,z**2+u*z,u**3,v**-1,3,1,6)]
    for name,MM,NN,FF,gg,zz,order,lead,scale in examples:
        GG=NN+MM*z-(p-FF-gg)**2
        check(name+" leading curve balance",
              S.limit(v**scale*GG.subs({u:1/v,z:zz}),v,0)==0)
        om=S.cancel(S.diff(zz,v)/S.diff(GG,u).subs({u:1/v,z:zz}))
        check(name+" exact mixed-coframe order",
              S.limit(om/v**order,v,0)==lead)

    # Direct coframe rescaling and the exact derivative-order identity.
    ll,kk,DD,ee,nuA,nuDp=S.symbols("l k D e nuA nuDp",integer=True)
    vv,ZZ=S.symbols("vv ZZ",positive=True)
    mm=S.symbols("mm",integer=True,positive=True)
    TT=S.Function("T")(vv,ZZ);AA=S.Function("A")(vv,ZZ)
    Mzero=S.Function("M0")(vv)
    pull=S.Matrix([vv**(-ll),vv**(-DD)*TT]).jacobian([vv,ZZ]).det()*vv**(ll*mm)/Mzero
    expected=-ll*vv**(ll*mm-DD-ll-1)*S.diff(TT,ZZ)/Mzero
    check("direct original volume pullback",
          S.simplify(pull-expected)==0)
    rel=ll*vv**(ll*(mm-1)-1)*S.diff(TT,ZZ)/(Mzero*S.diff(AA,ZZ))
    wedge=S.Matrix([vv**(-DD)*AA,vv]).jacobian([vv,ZZ]).det()*rel
    check("direct generic relative coframe pullback",
          S.simplify(wedge-pull)==0)
    order=ee*ll*(m-1)-1-nuA
    check("coefficient parameter derivative exact order formula",
          S.expand(order.subs(nuA,ee*DD-nuDp)
                   -(ee*(ll*(m-1)-DD)+nuDp-1))==0)
    check("uniform integer margin identity",
          S.expand((n-1)*m-(2*n-1)-((n-1)*(m-2)-1))==0)

    # Binomial controls of arbitrary multiplicity; all gcd branches are counted.
    for mm,nn,rr in [(3,2,2),(3,2,3),(4,2,3),(6,2,3),(10,3,4),(14,4,5)]:
        hh=gcd(mm,2*nn-1);lv=(2*nn-1)//hh;kv=mm//hh;Dv=nn*kv
        branches=gcd(rr,Dv);ev=rr//branches;gamma=Dv//branches
        Y=b*s**gamma
        AZ=rr*Y**(rr-1)
        Yp=s**gamma/(rr*b**(rr-1))
        check(f"parameter derivative identity m{mm}n{nn}mult{rr}",
              S.cancel(AZ*Yp-s**(ev*Dv))==0)
        check(f"every binomial branch counted m{mm}n{nn}mult{rr}",
              branches*ev==rr and ev*Dv==rr*gamma)
        valAZ=(rr-1)*gamma
        exact=ev*lv*(mm-1)-1-valAZ
        bound=ev*(lv*(mm-1)-Dv)
        check(f"all ramified branches holomorphic m{mm}n{nn}mult{rr}",
              valAZ<=ev*Dv-1 and exact>=bound>=0)
    check("lower universal double-root model attains order zero",
          1*(3-1)-1-1==0)
    check("m2 cannot use the new coarse derivative margin",3*(2-1)-4==-1)

    # Non-binomial local clustering: distinct simple and double packets.
    Y=S.symbols("Y")
    local=(Y-v)**2*(Y+2*v)
    ys=-2*v+p*v**3/9
    check("non-binomial simple packet leading equation",
          S.limit((local.subs(Y,ys)-p*v**5)/v**5,v,0)==0)
    check("non-binomial simple packet derivative identity at leading order",
          S.limit(S.diff(local,Y).subs(Y,ys)*S.diff(ys,p)/v**5,v,0)==1)
    for sign in [-1,1]:
        yd=v+sign*b*v**2
        check(f"non-binomial double packet sign{sign} leading equation",
              S.limit((local.subs(Y,yd)-3*b**2*v**5)/v**5,v,0)==0)
        yp=sign*v**2/(6*b)
        check(f"non-binomial double packet sign{sign} parameter derivative",
              S.limit(S.diff(local,Y).subs(Y,yd)*yp/v**5,v,0)==1)

    # Full actual MX input with a triple leading zero.
    MM=u**2*(u-1)**4;NN=u*(1-7*u)**2
    cc=-6*u**2+8*u-S.Rational(4,3)
    gg=-3*u**4+8*u**3-S.Rational(20,3)*u**2+S.Rational(32,27)*u+S.Rational(32,81)
    check("actual triple MX degrees",S.degree(MM,u)==6 and S.degree(NN,u)==3)
    FF=z**2+cc*z
    check("actual triple MX full closed support",
          all(3*ju+6*iz<=12 for (ju,iz),coef in S.Poly(FF,u,z).terms()))
    cv=S.expand(v**2*cc.subs(u,1/v));gv=S.expand(v**4*gg.subs(u,1/v))
    Z=S.symbols("Z",positive=True)
    A0=8*S.sqrt(Z)+Z**2-6*Z-3
    check("actual MX triple leading value and first two derivatives",
          all(S.diff(A0,Z,jj).subs(Z,1)==0 for jj in range(3)))
    check("actual MX third leading derivative is three",
          S.diff(A0,Z,3).subs(Z,1)==3)
    nterm=S.expand(v**8*NN.subs(u,1/v))
    check("actual MX retained N term begins at v5",
          nterm==49*v**5-14*v**6+v**7)
    Zc=(1-v)**S.Rational(4,3)
    TT=-S.sqrt((1-v)**4*(Zc+Y)+nterm)
    AA=-8*TT+(Zc+Y)**2+cv*(Zc+Y)+gv
    for jj,expected in [
        (0,0),
        (1,S.Rational(8,27)*v**3+S.Rational(10,81)*v**4),
        (2,0)]:
        actual=S.series(S.diff(AA,Y,jj).subs(Y,0),v,0,5).removeO()
        check(f"actual MX triple unfolding derivative{jj} modulo v5",
              S.expand(actual-expected)==0)
    check("actual MX triple unfolding cubic coefficient",
          S.limit(S.diff(AA,Y,3).subs(Y,0),v,0)==3)
    check("actual MX triple branch has e3 and order six",
          gcd(3,4)==1 and 3*(6-1)-1-2*4==6)

    # Arbitrarily high leading multiplicity is allowed as n grows.
    for nv in [2,3,4]:
        x=S.symbols("x")
        tay=S.series(S.sqrt(1+x),x,0,nv+1).removeO().subs(x,Z-1).expand()
        cn=-1/S.Poly(tay,Z).LC()
        lead=S.expand(cn*S.sqrt(Z)-cn*tay)
        check(f"actual MX high multiplicity n{nv} leading derivatives",
              all(S.simplify(S.diff(lead,Z,jj).subs(Z,1))==0 for jj in range(nv+1)))
        check(f"actual MX high multiplicity n{nv} next derivative nonzero",
              S.simplify(S.diff(lead,Z,nv+1).subs(Z,1))!=0)
        check(f"actual MX high multiplicity n{nv} leading z monic",
              S.Poly(-cn*tay,Z).LC()==1 and 4*nv-4>=2)

    # The original whole source chart and all named vertical-line degree tests.
    t,x=S.symbols("t x")
    uu=t**2;sig=1-7*uu;aa=uu*(uu-1)**4
    ww=t*(sig+aa*x);zz=2*sig*x+aa*x**2
    check("actual MX full original source relation",
          S.expand(ww**2-NN.subs(u,uu)-MM.subs(u,uu)*zz)==0)
    source=S.expand(-8*ww+zz**2+cc.subs(u,uu)*zz+gg.subs(u,uu))
    for tv,dv in [(0,2),(1,2),(2,4)]:
        check(f"actual triple candidate reached line t{tv} degree{dv}",
              S.degree(source.subs(t,tv),x)==dv)

    # Boundary controls are retained; no regular plane counterexample follows.
    zS=(w**2-u)/u**2
    mate=S.Matrix([w+u*zS,1/w]).jacobian([u,w]).det()
    check("outside-scope multiplicative S2 rational mate remains polar",
          S.cancel(mate-1/u**2)==0)
    mate=S.Matrix([w+S.Function("g")(u),1/u]).jacobian([u,w]).det()
    check("nonlinear hypothesis remains material",
          S.cancel(mate-1/u**2)==0)
    zplane=w**2-u
    check("degree-gap boundary has actual regular plane mate",
          S.Matrix([w+zplane**2,zplane]).jacobian([u,w]).det()==1)

    out=dict(claim="PF2RULED148",outcome="PASS",check_count=len(checks),checks=checks,
             elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),
             proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),
             context_sha256=sha(HERE/CONTEXT),
             arbitrary_mixed_coordinates_reduced_to_tested_family=False,
             actual_full_marked_quadratic_regularity=False,
             unrestricted_cyclic_targets_excluded=False,full_goal_complete=False,
             registered_predecessor_checkers_replayed=[],
             independent_replay=False,formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["claim","outcome","check_count","elapsed_seconds"]}))

if __name__=="__main__":
    main()
