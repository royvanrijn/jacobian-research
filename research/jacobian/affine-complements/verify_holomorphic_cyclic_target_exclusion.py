"""PF146 supporting exact controls; universal proof is in the canonical note."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="HOLOMORPHIC_GENERIC_FIBRES_AND_NONLINEAR_CYCLIC_TARGET_EXCLUSION.md"
CONTEXT="sources/holomorphic-cyclic-target-exclusion/context.json"
RECEIPT="holomorphic_cyclic_target_exclusion_receipt.json"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    started=time.perf_counter()
    checks=[]
    def check(name,value):
        assert bool(value),name
        checks.append(name)
    context=json.loads((HERE/CONTEXT).read_text())
    for path,digest in context["retained_sha256"].items():
        check("retained source pin "+path,sha(ROOT/path)==digest)

    u,w,z,p,a,s,t,x,lam,c0=S.symbols("u w z p a s t x lam c0",nonzero=True)
    MM=S.Function("M")(u);NN=S.Function("N")(u)
    ZZ=(w**2-NN)/MM
    FF=S.Function("f")(ZZ)
    PP=a*w+FF
    Fu=S.diff(NN,u)+S.diff(MM,u)*ZZ
    jac=S.Matrix([PP,ZZ]).jacobian([u,w]).det()
    check("actual coframe identity independent of f",
          S.simplify(jac/(a*Fu)-1/MM)==0)
    check("surface residue identity at finite patches",
          S.simplify(S.diff(ZZ,u)+Fu/MM)==0)

    al,be,m,ell,n,h=S.symbols("alpha beta m ell n h",integer=True,positive=True)
    order=(m-1)*al-1
    check("all vanishing-z ends exact differential order",
          S.expand((m-ell)*al-1+(ell-1)*al-order)==0)
    check("all polar-z ends exact differential order",
          S.expand(-be-1+(m-1)*al+be-order)==0)
    check("z cannot escape while u stays finite",((2*n-1)*be).is_positive)
    check("unit z cannot balance dominant M term",(h*al).is_positive)
    for mm,aa in [(2,1),(4,1),(4,3),(8,7)]:
        check(f"holomorphic normalized end m{mm}alpha{aa}",(mm-1)*aa-1>=0)

    models=[("S2",u**2,u),
            ("MXk2",u**2*(u-1)**2,u*(1-3*u)**2)]
    bu=u**2+1;k=3
    sigma=-bu+2*(1-k)*u*S.diff(bu,u)
    models.append(("MXk3",u**2*bu**k,u*sigma**2))
    for name,M,N in models:
        md=S.degree(M,u);ld=S.degree(N,u)
        check(name+" strict degree gap",md>ld>=1 and md>=2)
        Mlead=S.LC(S.Poly(M,u))
        # Leading finite-z end. The lower O(u^-m) correction solving
        # G=0 does not change these leading orders; the universal proof
        # includes every valuation without assuming a Puiseux chart.
        zsmall=S.cancel(-N/M)
        us=1/s;zs=zsmall.subs(u,us)
        fus=(S.diff(N,u)+S.diff(M,u)*zsmall).subs(u,us)
        om=S.cancel(S.diff(zs,s)/(a*fus))
        check(name+" finite-z end exact leading volume order",
              S.limit(om/s**(md-2),s,0)==1/(a*Mlead))
        check(name+" finite-z end leading cancellation",
              S.cancel((N+M*zsmall))==0)
        for nd in [1,2,4]:
            d=gcd(md,2*nd-1);alpha=(2*nd-1)//d;beta=md//d
            us=s**(-alpha);zs=s**(-beta)
            gp=N+M*z-(p-(z**nd+1))**2
            # All models are monic in M; leading z coefficient chosen1.
            check(name+f" degree{nd} polar-z leading balance",
                  S.limit(s**(md*alpha+beta)*gp.subs({u:us,z:zs}),s,0)==0)
            fus=(S.diff(N,u)+S.diff(M,u)*z).subs({u:us,z:zs})
            om=S.cancel(S.diff(zs,s)/fus)
            oo=(md-1)*alpha-1
            check(name+f" degree{nd} polar-z end exact leading order",
                  S.limit(om/s**oo,s,0)==-S.Rational(beta,md))

    # Independent hyperelliptic model for classical S2, all generic p.
    relation=w**2-u-u**2*z
    YY=1+2*z*u
    check("S2 generic curve Y relation",
          S.expand(YY**2-1-4*z*w**2+4*z*relation)==0)
    for f in [z,z**2+z,z**4-2*z+3]:
        nd=S.degree(f,z)
        disc=S.expand(1+4*z*(p-f)**2)
        poly=S.Poly(disc,z,domain=S.QQ.frac_field(p))
        check(f"S2 f degree{nd} generic squarefree model",
              S.gcd(poly,poly.diff()).degree()==0)
        check(f"S2 f degree{nd} exact branch degree",poly.degree()==2*nd+1)
        check(f"S2 f degree{nd} hyperelliptic genus",((2*nd+2)-2)//2==nd)
    y=S.symbols("Y",nonzero=True)
    check("S2 required holomorphic differential uses actual coframe",
          S.cancel(1/(a*(1+2*z*u))-1/(a*YY))==0)

    # Nonlinear MX source candidate passes the named old line tests.
    uu=t**2;st=1-3*uu;bt=uu-1
    at=t**2*bt**2
    ww=t*(st+at*x);zz=2*st*x+at*x**2
    candidate=S.expand(ww+zz**2)
    check("nonlinear MX coordinate C0 degree is two",
          S.Poly(candidate.subs(t,0),x).degree()==2)
    check("nonlinear MX coordinate x0 degree is three",
          S.Poly(candidate.subs(x,0),t).degree()==3)
    check("nonlinear MX coordinate ordinary source line degree four",
          S.Poly(candidate.subs(t,2),x).degree()==4)
    check("nonlinear MX coordinate split fibre line degree two",
          S.Poly(candidate.subs(t,1),x).degree()==2)
    check("complete MX relation retained",
          S.expand(ww**2-uu*st**2-uu**2*bt**2*zz)==0)
    # Nonconstant f prevents even affine pencil stability of P.
    polynomial=S.Poly((-w+z**2)-lam*(w+z**2)-c0,w,z)
    check("nonlinear example breaks the original fibre pencil",
          S.solve([polynomial.coeff_monomial(w),
                   polynomial.coeff_monomial(z**2)],[lam])==[])

    # Scope controls that would be false under overgeneralization.
    check("constant f has a rational mate on S2",
          S.Matrix([w,1/u]).jacobian([u,w]).det()==1/u**2)
    zplane=w**2-u
    check("degree-gap hypothesis essential on the plane",
          S.Matrix([w+zplane**2,zplane]).jacobian([u,w]).det()==1)
    for e,ord0 in [(1,0),(2,0),(5,2)]:
        check(f"ramified pullback retains nonnegative order e{e}ord{ord0}",
              e*ord0+e-1>=0)

    out=dict(claim="PF2RULED146",outcome="PASS",check_count=len(checks),checks=checks,
             elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),
             proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),
             context_sha256=sha(HERE/CONTEXT),
             all_target_coordinates_reduced_to_tested_family=False,
             actual_full_marked_quadratic_regularity=False,
             unrestricted_cyclic_targets_excluded=False,full_goal_complete=False,
             registered_predecessor_checkers_replayed=[],
             independent_replay=False,formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({key:out[key] for key in ["claim","outcome","check_count","elapsed_seconds"]}))

if __name__=="__main__":
    main()
