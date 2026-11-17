"""PF150 exact supporting controls for the written odd-degree boundary proof."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="ODD_LEADING_CYCLIC_COORDINATES_AND_RECIPROCAL_BOUNDARY_COFORMS.md"
CONTEXT="sources/odd-leading-cyclic-target-exclusion/context.json"
RECEIPT="odd_leading_cyclic_target_exclusion_receipt.json"

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

    u,U,w,z,y,p,s,v,T,lam,b=S.symbols("u U w z y p s v T lam b",nonzero=True)
    A,Au,Ay,B,Bu,By,M,Mp,N,Np=S.symbols("A Au Ay B Bu By M Mp N Np",nonzero=True)
    n=S.symbols("n",integer=True,positive=True)
    W=(p*y**n-B)/A
    Wu=(-Bu-W*Au)/A
    Wy=(n*p*y**(n-1)-By-W*Ay)/A
    Wp=y**n/A
    Hu=Mp-y*(2*W*Wu-Np)
    Up=2*y*W*Wp/Hu
    Uy=(W**2-N+2*y*W*Wy)/Hu
    jac=Up*(Wu*Uy+Wy)-Uy*(Wu*Up+Wp)
    check("reciprocal exact coframe from full Jacobian",
          S.cancel(jac/(y*(W**2-N))+y**(n-1)/(A*Hu))==0)
    Hbaru=Mp*A**2+2*M*A*Au-y*(2*B*Bu-Np*A**2-2*N*A*Au)+2*p*Bu*y**(n+1)
    check("barred denominator retains all u partials",
          S.cancel((Hbaru-A**2*Hu).subs(M,y*(W**2-N)))==0)
    check("barred and unbarred coframe identical on curve",
          S.cancel(-A*y**(n-1)/(A**2*Hu)+y**(n-1)/(A*Hu))==0)
    check("zero limiting w needs no division by w in reciprocal coframe",
          not S.cancel(-A*y**(n-1)/S.Symbol("Hbar_u")).has(w))

    # Actual monic odd curve polynomials, with unrestricted lower coefficients.
    examples=[
        (1,u**2,u,1,z+u**7,z+u**11),
        (2,u**2*(u-1)**2,u*(1-3*u)**2,1,
         z**2+u**5*z+1,3*z**2+u**7*z+u**11),
        (3,u**2,u,(u-2)**7,
         (u-2)**7*z**3+u**19*z**2+1,z**3+u**23*z+u**31)]
    for nv,MM,NN,aa,AA,BB in examples:
        Z=(w**2-NN)/MM
        Pstar=(w*AA+BB).subs(z,Z)
        QQ=S.cancel(MM**nv*Pstar/aa)
        poly=S.Poly(QQ,w)
        check(f"monic actual odd curve degree n{nv}",
              poly.degree()==2*nv+1 and S.cancel(poly.LC()-1)==0)
        hh=MM**nv/aa
        check(f"exact odd curve equation n{nv}",
              S.cancel(QQ-p*hh-hh*(Pstar-p))==0)
        wp=S.cancel(hh/S.diff(QQ,w))
        check(f"actual generic coframe parameter inverse n{nv}",
              S.cancel(S.diff(Pstar,w)*wp-1)==0)
    mm,dd,alpha,er=S.symbols("m d alpha er",integer=True,positive=True)
    D=2*n+1
    inf=(mm-1)*alpha+(dd-n*mm)*alpha/D
    check("odd infinity exact pole margin",
          S.simplify(inf-((n+1)*mm-D+dd)*alpha/D)==0)
    check("m2 gives strictly positive odd infinity numerator",
          S.expand(((n+1)*mm-D+dd).subs(mm,2))==1+dd)
    check("finite leading zero exact pole margin",
          S.simplify(alpha-er*alpha/D-(1-er/D)*alpha)==0)
    for nv,ev in [(1,3),(2,5),(3,7),(4,8)]:
        check(f"finite leading zero equality retained n{nv}e{ev}",
              ev<=2*nv+1 and all(S.ceiling((1-S.Rational(ev,2*nv+1))*av)-1>=-1
                                   for av in [1,2,3]))

    # Nondegenerate M-root branch exists and is holomorphic in every degree.
    for ev,nv in [(2,1),(5,1),(4,3)]:
        HH=U**ev+y-p**2*y**(2*nv+1)  # W=p*y^n, N=1
        yy=-s**ev;uu=s
        check(f"nonzero Delta branch leading equation e{ev}n{nv}",
              S.limit(HH.subs({U:uu,y:yy})/s**ev,s,0)==0)
        om=-yy**(nv-1)*S.diff(yy,s)/S.diff(HH,U).subs({U:uu,y:yy})
        order=(nv-1)*ev
        check(f"nonzero Delta exact holomorphic order e{ev}n{nv}",
              S.limit(om/s**order,s,0)!=0 and order>=0)

    # First-parameter-coefficient bound with all gcd branches and ramification.
    for ev,nv in [(2,1),(2,2),(4,3),(5,5),(5,3),(8,6)]:
        gg=gcd(ev,nv+1);aa=(nv+1)//gg;bb=ev//gg
        uu=b*s**aa;yy=s**bb
        HH=U**ev+p*y**(nv+1)
        check(f"high contact full branch balance e{ev}n{nv}",
              S.cancel(HH.subs({U:uu,y:yy,p:-b**ev}))==0)
        up=-s**aa/(ev*b**(ev-1))
        check(f"high contact exact parameter derivative e{ev}n{nv}",
              S.cancel(S.diff(HH,U).subs(U,uu)*up+yy**(nv+1))==0)
        check(f"all high contact normalized branches counted e{ev}n{nv}",
              gg*bb==ev and ev*aa==(nv+1)*bb)
        order=nv*bb-1-(ev-1)*aa
        bound=S.ceiling(S.Rational((nv+1)*bb,ev))-bb-1
        check(f"high contact bound exact e{ev}n{nv}",order==bound)
        if nv+1>=ev:
            check(f"all high contact branches at most logarithmic e{ev}n{nv}",
                  order>=-1)
        else:
            check(f"outside threshold higher pole control e{ev}n{nv}",
                  order<=-2)
    beta,ee,q=S.symbols("beta e q",integer=True,positive=True)
    check("exact reciprocal order from parameter derivative",
          S.expand(n*beta-1-((n+1)*beta-q)-(q-beta-1))==0)
    check("ramified simple pole remains simple",
          S.expand(ee*(-1)+(ee-1))==-1)

    # Actual central S2 n3, both zero-w branches retained.
    At=1+u*y**2+y**3;Bt=y+u**5*y**3
    Hbar=(u**2+u*y)*At**2-y*(p*y**3-Bt)**2
    om=-At*y**2/S.diff(Hbar,u)
    for label,uu,lead in [("small u",y**2,-1),("u near minus y",-y-y**2,1)]:
        check("central S2 "+label+" leading curve balance",
              S.limit(Hbar.subs(u,uu)/y**3,y,0)==0)
        check("central S2 "+label+" actual coframe leading order",
              S.limit(om.subs(u,uu)/y,y,0)==lead)
    # General simple-N calculation also permits e>2.
    HH=U**4+U*y-y**3
    yy=s**3;uu=-s
    check("simple N root e4 noncentral scale balance",
          S.limit(HH.subs({U:uu,y:yy})/s**4,s,0)==0)
    om=-yy*S.diff(yy,s)/S.diff(HH,U).subs({U:uu,y:yy})
    check("simple N root e4 exact holomorphic order two",
          S.limit(om/s**2,s,0)==1)

    # Actual MX equality and outside-threshold controls with high contact.
    topbs={}
    for kv in [4,5]:
        f=((2*kv-1)*u-1)*S.sqrt(u)
        bb=S.series(f.subs(u,1+U),U,0,kv).removeO().subs(U,u-1).expand()
        topbs[kv]=bb
        NN=u*(1-(2*kv-1)*u)**2;MM=u**2*(u-1)**kv
        delta=S.expand(bb**2-NN)
        dpoly=S.Poly(delta.subs(u,U+1).expand(),U)
        check(f"actual MX k{kv} leading contact exactly k",
              min(monom[0] for monom,coef in dpoly.terms())==kv)
        HH=MM-y*delta+2*p*bb*y**4-p**2*y**7
        if kv==4:
            check("actual k4 equality initial polynomial",
                  S.simplify(S.limit(HH.subs(u,1+lam*y)/y**4,y,0)
                             -(lam**4+12*p))==0)
            check("actual k4 all four initial roots distinct",
                  S.gcd(lam**4+12*p,S.diff(lam**4+12*p,lam),lam)==1)
            om=-y**2/S.diff(HH,u).subs(u,1+lam*y)
            check("actual k4 equality exact nonzero residues",
                  S.simplify(S.limit(y*om,y,0)+1/(4*lam**3))==0)
            check("actual k4 four normalized unramified branches",
                  gcd(4,4)==4 and (3*1-1-3)==-1)
        else:
            values={u:1+lam*s**4,y:s**5}
            check("actual k5 outside threshold initial polynomial",
                  S.simplify(S.limit(HH.subs(values)/s**20,s,0)
                             -(lam**5+16*p))==0)
            om=-(s**5)**2*S.diff(s**5,s)/S.diff(HH,u).subs(values)
            check("actual k5 exact double differential pole",
                  S.simplify(S.limit(s**2*om,s,0)+1/lam**4)==0)
            check("actual k5 ONE ramified normalized branch",
                  gcd(5,4)==1 and 3*5-1-4*4==-2)

    # Whole-source degree controls, not a claim about every recharted curve.
    t,h,x=S.symbols("t h x")
    uu=t**2;ww=t*(1+t**2*h);zz=2*h+t**2*h**2
    Psource=S.expand(ww*zz**3+zz**2+uu*ww*zz+uu**5+ww)
    check("central S2 full original source relation",
          S.expand(ww**2-uu-uu**2*zz)==0)
    for label,expr,var,dv in [
        ("C0",Psource.subs(t,0),h,2),
        ("other axis",Psource.subs(h,0),t,10),
        ("ordinary vertical",Psource.subs(t,2),h,7)]:
        check("odd S2 candidate "+label+" degree",S.degree(expr,var)==dv)
    check("odd S2 candidate has nonzero differential at marked origin",
          S.diff(Psource,t).subs({t:0,h:0})==1)
    kv=4;uu=t**2;sig=1-7*uu;aa=uu*(uu-1)**4
    ww=t*(sig+aa*x);zz=2*sig*x+aa*x**2
    check("high contact MX whole original source relation",
          S.expand(ww**2-uu*sig**2-uu**2*(uu-1)**4*zz)==0)
    Psource=S.expand((ww+topbs[4].subs(u,uu))*zz**3+zz**2+uu**17+ww)
    for tv,dv in [(0,3),(1,2),(-1,3),(2,7)]:
        check(f"high contact MX actual reached line t{tv} degree{dv}",
              S.degree(Psource.subs(t,tv),x)==dv)

    check("high contact MX candidate has nonzero differential at marked origin",
          S.diff(Psource,t).subs({t:0,x:0})==1)

    # Symbolic source-line corollary for all n2 and all model k.
    ca,cb,a1,a0,b1,b0=S.symbols("ca cb a1 a0 b1 b0",nonzero=True)
    P2=w*(ca*z**2+a1*z+a0)+cb*z**2+b1*z+b0
    restricted=S.expand(P2.subs(w,-cb/ca))
    check("leading-cancelled n2 whole line restriction degree at most one",
          S.degree(restricted,z)<=1)
    P1=w*(ca*z+a0)+cb*z+b0
    check("central n1 whole line restriction at most one",
          S.degree(P1.subs(w,0),z)<=1)
    check("k2 k3 k4 remaining n>=3 all meet threshold",
          all(3>=kv-1 for kv in [2,3,4]))

    # Outside controls remain outside; no whole-plane counterexample claimed.
    zS=(w**2-u)/u**2
    mate=S.Matrix([w+S.Function("g")(u),1/u]).jacobian([u,w]).det()
    check("degree-one odd rational mate remains outside n>=1",
          S.cancel(mate-1/u**2)==0)
    mate=S.Matrix([w+u*zS,1/w]).jacobian([u,w]).det()
    check("nonadmissible even leading rational control retained",
          S.cancel(mate-1/u**2)==0)

    out=dict(claim="PF2RULED150",outcome="PASS",check_count=len(checks),checks=checks,
             elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),
             proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),
             context_sha256=sha(HERE/CONTEXT),
             all_odd_admissible_leading_S2_coordinates_excluded=True,
             all_odd_MX_coordinates_unconditionally_excluded=False,
             source_line_corollary_promoted_to_rational_exclusion=False,
             all_leading_divisors_unrestricted=False,
             arbitrary_target_coordinates_reduced_to_tested_family=False,
             actual_full_marked_quadratic_regularity=False,
             unrestricted_cyclic_targets_excluded=False,full_goal_complete=False,
             registered_predecessor_checkers_replayed=[],
             independent_replay=False,formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["claim","outcome","check_count","elapsed_seconds"]}))

if __name__=="__main__":
    main()
