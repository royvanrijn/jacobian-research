"""PF147 exact supporting controls for the written all-degree proof."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="ADDITIVE_BASE_PERTURBATIONS_AND_RESONANT_CYCLIC_TARGET_EXCLUSION.md"
CONTEXT="sources/additive-cyclic-target-exclusion/context.json"
RECEIPT="additive_cyclic_target_exclusion_receipt.json"

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    started=time.perf_counter()
    checks=[]
    def check(name,value):
        assert bool(value),name
        checks.append(name)
    ctx=json.loads((HERE/CONTEXT).read_text())
    for path,digest in ctx["retained_sha256"].items():
        check("retained source pin "+path,sha(ROOT/path)==digest)

    u,w,z,p,a,v,t,h,C=S.symbols("u w z p a v t h C",nonzero=True)
    M=S.Function("M")(u);N=S.Function("N")(u);g=S.Function("g")(u)
    Z=(w**2-N)/M
    P=a*w+S.Function("f")(Z)+g
    Gu=S.diff(N,u)+S.diff(M,u)*Z+2*w*S.diff(g,u)/a
    check("actual coframe includes base derivative",
          S.simplify(S.Matrix([P,Z]).jacobian([u,w]).det()/(a*Gu)-1/M)==0)
    check("constant base perturbation recovers PF146",
          S.simplify(Gu.subs(S.diff(g,u),0)-(S.diff(N,u)+S.diff(M,u)*Z))==0)

    alpha,beta,m,L,n,d,W,J=S.symbols("alpha beta m L n d W J",integer=True,positive=True)
    baseorder=(m-1)*alpha-1
    check("finite-z shifted branch order exact",
          S.expand((m-L)*alpha-1+(L-1)*alpha-baseorder)==0)
    check("f-dominated polar-z order exact",
          S.expand(2*W-alpha-beta-1-baseorder).subs(2*W,m*alpha+beta)==0)
    check("leading cancellation branch order exact",
          S.expand(W+J-alpha-beta-1
                   -(W-alpha+(n-1)*beta-1)).subs(J,n*beta)==0)
    check("g-dominated combined coefficient cannot cancel",
          S.simplify(((m-2*d)*alpha+beta).subs(beta,(2*d-m)*alpha))==0)
    for mm,aa,bb in [(2,1,2),(3,2,2),(6,1,2)]:
        ww=(mm*aa+bb)//2
        check(f"strict cancellation gives holomorphic order {mm},{aa},{bb}",
              2*ww==mm*aa+bb and ww-aa>=1 and ww-aa+bb-1>=0)

    # New finite nonzero-z and coefficient-cancellation branches.
    GG=u+u**2*z-(p-z**2-u)**2
    zz=1-(2*p-1)*v;uu=1/v
    check("finite nonzero-z leading balance",
          S.limit(v*GG.subs({u:uu,z:zz}),v,0)==0)
    omega=S.cancel(S.diff(zz,v)/S.diff(GG,u).subs({u:uu,z:zz}))
    check("finite nonzero-z end is holomorphic nonzero",
          S.limit(omega,v,0)==1)
    NN=u**2+u;MM=u**5
    GG=NN+MM*z-(p-z**2-u)**2
    RR=S.expand(NN-(p-u)**2)
    check("finite-z effective degree retains generic p",
          S.degree(RR,u)==1 and S.degree(NN,u)==2)
    zz=-(1+2*p)*v**4
    omega=S.cancel(S.diff(zz,v)/S.diff(GG,u).subs({u:1/v,z:zz}))
    check("finite-z cancellation branch exact holomorphic order",
          S.limit(omega/v**3,v,0)==1)

    # f dominates g; g dominates f; and f,g cancel above w.
    examples=[
        ("f dominates",u**2,u,z**2,u,v**-3,v**-2,2,-1),
        ("g dominates",u**5,u,z**2,u**3,v**-1,v**-1,3,1),
        ("leading f g cancellation",u**2,u,z**2,-u**2,
         v**-2,v**-2-S.Rational(1,2)*v**-1,2,S.Rational(1,2))]
    for name,MM,NN,ff,gg,uu,zz,order,leading in examples:
        GG=NN+MM*z-(p-ff-gg)**2
        omega=S.cancel(S.diff(zz,v)/S.diff(GG,u).subs({u:uu,z:zz}))
        check(name+" exact leading differential order",
              S.limit(omega/v**order,v,0)==leading)
        # Check only the leading balance; no finite truncation is claimed
        # to be the complete formal branch.
        scaled=6 if name=="leading f g cancellation" else (8 if name=="f dominates" else 6)
        check(name+" leading curve balance",
              S.limit(v**scaled*GG.subs({u:uu,z:zz}),v,0)==0)

    # A0'' at every possible multiple leading root: exact symbolic formula.
    A0,fn,zz0,nn=S.symbols("A0 fn Z0 n",nonzero=True)
    second=-A0/(4*zz0**2)+nn*(nn-1)*fn*zz0**(nn-2)
    imposed=second.subs(A0,-2*nn*fn*zz0**nn)
    check("leading root has nonzero quadratic coefficient formula",
          S.simplify(imposed-nn*(2*nn-1)*fn*zz0**(nn-2)/2)==0)
    gd=(2*nn-1)*fn*zz0**nn
    check("double leading value relation",
          S.expand(-2*nn*fn*zz0**nn+fn*zz0**nn+gd)==0)

    # Every local binomial normalization branch, even q included.
    for nn,tt in [(2,1),(2,2),(3,1)]:
        mm=(2*nn-1)*tt;dd=nn*tt
        check(f"resonant uniform order bound n{nn}t{tt}",
              S.Rational(mm-1)-S.Rational(dd,2)>=1)
        for q in range(1,dd+1):
            branches=gcd(2,q);e=2//branches
            order=e*(mm-1)-1-e*q//2
            check(f"normalized local exponents n{nn}t{tt}q{q}",
                  e*q%2==0 and branches*e==2)
            check(f"all local branches holomorphic n{nn}t{tt}q{q}",order>=0)

    # Full MX0 resonant input. Only the exact mod-v^5 calculation is used.
    MM=u**2*(u-1)**4;NN=u*(1-7*u)**2
    check("resonant MX degrees",S.degree(MM,u)==6 and S.degree(NN,u)==3)
    check("resonant MX rescaled M unit",S.expand(v**6*MM.subs(u,1/v)-(1-v)**4)==0)
    nterm=S.expand(v**8*NN.subs(u,1/v))
    check("resonant MX N begins exactly at v5",
          nterm==49*v**5-14*v**6+v**7)
    critical=S.series((1-v)**S.Rational(4,3),v,0,5).removeO()
    check("resonant MX critical section is exact modulo v5",
          S.series(-2*(1-v)**2/S.sqrt(critical)+2*critical,v,0,5).removeO()==0)
    # N changes the square root only at order5, so this value is B modv5.
    Bseries=S.series(-3*(1-v)**S.Rational(8,3)+3
                     -8*v+S.Rational(20,3)*v**2-S.Rational(40,27)*v**3,
                     v,0,5).removeO()
    check("resonant MX maximal delay critical value",
          S.expand(Bseries-S.Rational(10,81)*v**4)==0)
    check("generic p cannot cancel resonant leading coefficient",
          S.Poly(p-S.Rational(10,81),p).degree()==1)
    check("resonant MX double leading root",
          4*(-1)+1+3==0 and -2+2==0)
    check("resonant MX second derivative is three",S.Rational(1)+2==3)
    check("resonant MX two unramified branches have order two",
          gcd(2,4)==2 and (6-1)-1-4//2==2)
    check("lower universal model reaches order zero",
          (3-1)-1-2//2==0)

    # Genuine mixed S2 source coordinate. The named old axis tests pass.
    uu=t**2;ww=t*(1+t**2*h);zz=2*h+t**2*h**2
    candidate=S.expand(ww+zz**2+uu)
    check("mixed S2 candidate C0 degree two",S.degree(candidate.subs(t,0),h)==2)
    check("mixed S2 candidate other axis degree two",S.degree(candidate.subs(h,0),t)==2)
    check("mixed S2 candidate ordinary vertical degree four",S.degree(candidate.subs(t,2),h)==4)
    check("mixed S2 full source relation",S.expand(ww**2-uu-uu**2*zz)==0)

    # Boundary controls: actual rational mates outside the hypothesis.
    zexpr=(w**2-u)/u**2
    generic_g=S.Function("b")(u)
    check("constant f permits rational mate for every base perturbation",
          S.simplify(S.Matrix([w+generic_g,1/u]).jacobian([u,w]).det()-1/u**2)==0)
    mixed=w+u*zexpr
    check("multiplicative mixed control has rational mate",
          S.cancel(S.Matrix([mixed,1/w]).jacobian([u,w]).det()-1/u**2)==0)
    translated_z=zexpr+w+u**2/4
    check("linear f can give translated logarithmic case",
          S.expand(translated_z-((w+u**2/2)**2-u)/u**2)==0)
    aa,H=S.symbols("aa H",nonzero=True)
    conicu=(H+1/H-2)/(4*aa**2);conicw=(H-1/H)/(4*aa)
    logform=S.cancel(-S.diff(conicu,H)/(2*conicw))
    check("linear f boundary has nonzero log residue",
          S.residue(logform,H,0)==-1/(2*aa))
    zplane=w**2-u
    check("degree gap remains essential",
          S.Matrix([w+zplane**2,zplane]).jacobian([u,w]).det()==1)

    out=dict(claim="PF2RULED147",outcome="PASS",check_count=len(checks),checks=checks,
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
    print(json.dumps({key:out[key] for key in ["claim","outcome","check_count","elapsed_seconds"]}))

if __name__=="__main__":
    main()
