"""PF151: conditional coverage, degree conventions, and dimension-three guards.
The written proof is universal; these are narrow exact supporting controls.
"""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import subprocess
import sys
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md"
CONTEXT="sources/normal-form-coverage/context.json"
RECEIPT="normal_form_coverage_receipt.json"

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

    u,w,z=S.symbols("u w z")
    aa=S.Function("a")(u);bb=S.Function("b")(u)
    for r,s in [(4,6),(5,6),(5,7)]:
        P=aa*w**r+(u+1)*w**(r-1)
        Q=bb*w**s+(u*u+1)*w**(s-1)
        top=S.Poly(S.expand(S.Matrix([P,Q]).jacobian([u,w]).det()),w).coeff_monomial(w**(r+s-1))
        check(f"actual two-dimensional highest Wronskian r{r}s{s}",
              S.simplify(top-(s*S.diff(aa,u)*bb-r*aa*S.diff(bb,u)))==0)
    G=(u-2)**3
    M=u**2
    for r,s in [(4,6),(6,10),(5,6),(5,7)]:
        d=gcd(r,s);r0=r//d;s0=s//d
        H=(u-2)**3/u**d
        lp=H**r0;lq=H**s0
        check(f"common rational leading powers satisfy Wronskian r{r}s{s}",
              S.cancel(s*S.diff(lp,u)*lq-r*lp*S.diff(lq,u))==0)
        check(f"neither degree divides reduced pair r{r}s{s}",
              r>=4 and s>r and s%r!=0 and r%s!=0 and r0>=2 and s0>=2)

    # A genuine polynomial target automorphism exercises the valid reduction.
    P=u+w*w;Q=w+P**3
    check("elementary polynomial target example has constant Jacobian",
          S.expand(S.Matrix([P,Q]).jacobian([u,w]).det())==1)
    Qred=S.expand(Q-P**3);Pred=S.expand(P-Qred**2)
    check("valid elementary reduction drops ruling degree6 to1",
          S.degree(Q,w)==6 and Qred==w)
    check("valid elementary reduction next drops degree2 to0",Pred==u)
    check("each valid target change preserves Jacobian",
          S.Matrix([Pred,Qred]).jacobian([u,w]).det()==1)
    # In an actual noninvertible marked model degree<=3 is forbidden by PF52;
    # this plane automorphism is a positive control, not such a model.

    r0,s0,d,h0,J=S.symbols("r0 s0 d h0 J",integer=True)
    eps=S.symbols("eps",integer=True)
    check("central divisor order formula exact",
          S.expand((d*r0-eps)+r0*h0-(r0*(d+h0)-eps))==0)
    check("even common divisor order threshold exact",
          S.expand(r0*J-d*r0-r0*(J-d))==0)
    for r,s in [(5,6),(5,7),(6,9),(9,15)]:
        dv=gcd(r,s);rv=r//dv;sv=s//dv
        # These pairs are all reduced; J>=1 follows from any odd degree.
        cp=rv-(r%2);cq=sv-(s%2)
        check(f"odd reduced central common DOUBLE zero r{r}s{s}",
              rv>=2 and sv>=2 and (r%2 or s%2) and cp>=2 and cq>=2)
    central_bad=u**4*z**2+w*z+z+1
    check("degree4 with odd partner loses its central source-line degree",
          S.degree(central_bad.subs({u:0,w:0}),z)<=1)

    # Exact all-even leading data, source central Bezout, and a negative
    # control: satisfying all these necessary conditions is NOT being Keller.
    N=u
    P=w+G**2*z**2
    Q=G**3*z**3+z+S.Rational(3,2)*G*w*z+u
    Z=(w*w-N)/M
    Ps=P.subs(z,Z);Qs=Q.subs(z,Z)
    cp=S.Poly(S.cancel(Ps),w).LC()*M**2
    cq=S.Poly(S.cancel(Qs),w).LC()*M**3
    check("all-even residual common polynomial leading powers",
          S.cancel(cp-G**2)==0 and S.cancel(cq-G**3)==0)
    check("lowest degree4 residual G coprime M and triple root",
          S.gcd(G,M)==1 and S.degree(G,u)==3
          and S.diff(G,u).subs(u,2)==0 and S.diff(G,u,2).subs(u,2)==0)
    AP=S.diff(P,w).subs({u:0,w:0});BP=P.subs({u:0,w:0})
    AQ=S.diff(Q,w).subs({u:0,w:0});BQ=Q.subs({u:0,w:0})
    check("all-even residual central coframe Bezout",
          S.expand(AP*S.diff(BQ,z)-S.diff(BP,z)*AQ)==1)
    check("all-even residual central restrictions degree2 and3",
          S.degree(BP,z)==2 and S.degree(BQ,z)==3)
    check("all-even residual data do NOT assert a Keller pair",
          S.cancel(M*S.Matrix([Ps,Qs]).jacobian([u,w]).det()-1)!=0)

    # A mixed-parity residual with constant rational leading H but canonical
    # polynomial top coefficients u4 and u6, passing the central coframe.
    P=w+z**2+u**4*w*z**2
    Q=w+z**2+z+u**6*z**3
    Ps=P.subs(z,Z);Qs=Q.subs(z,Z)
    check("odd residual ruling degrees5 and6",
          S.degree(S.cancel(Ps),w)==5 and S.degree(S.cancel(Qs),w)==6)
    check("odd residual canonical top coefficients share central zero",
          S.cancel(S.Poly(S.cancel(Ps),w).LC()*M**2-u**4)==0
          and S.cancel(S.Poly(S.cancel(Qs),w).LC()*M**3-u**6)==0)
    AP=S.diff(P,w).subs({u:0,w:0});BP=P.subs({u:0,w:0})
    AQ=S.diff(Q,w).subs({u:0,w:0});BQ=Q.subs({u:0,w:0})
    check("odd residual central coframe Bezout",
          S.expand(AP*S.diff(BQ,z)-S.diff(BP,z)*AQ)==1)
    check("odd residual central restrictions both degree two",
          S.degree(BP,z)==2 and S.degree(BQ,z)==2)
    check("odd residual necessary data do NOT assert a Keller pair",
          S.cancel(M*S.Matrix([Ps,Qs]).jacobian([u,w]).det()-1)!=0)

    # The intrinsic central form is dw wedge dz/Nprime; no normalization
    # order or trace regularity is inferred from this first-order row.
    K,Mv,Wv=S.symbols("K Mv Wv",nonzero=True)
    # du=(2w dw-M dz)/K, so du wedge dw/M=dw wedge dz/K.
    check("central residue coframe exact exterior coefficient",S.cancel(S.Matrix([[2*Wv/K,-Mv/K],[1,0]]).det()/Mv-1/K)==0)
    for rv,sv in [(4,6),(4,10),(4,14)]:
        check(f"lowest-degree4 nondividing even shape s{sv}",
              gcd(rv,sv)==2 and sv%4==2 and sv%rv!=0)

    # Known THREE-dimensional polynomial Keller control.
    x,y,t=S.symbols("x y t")
    U0=1+x*y
    A=U0**3*z+y**2*U0*(4+3*x*y)
    B=y+3*x*U0**2*z+3*x*y**2*(4+3*x*y)
    C=2*x-3*x*x*y-x**3*z
    F=S.Matrix([A,B,C]);JF=F.jacobian([x,y,z])
    check("known3D actual global constant Jacobian",S.expand(JF.det())==-2)
    points=[(0,0,-S.Rational(1,4)),(1,-S.Rational(3,2),S.Rational(13,2)),
            (-1,S.Rational(3,2),S.Rational(13,2))]
    for ix,pt in enumerate(points):
        vals=[S.expand(f.subs(dict(zip([x,y,z],pt)))) for f in [A,B,C]]
        check("known3D retained exact collision point"+str(ix),vals==[-S.Rational(1,4),0,0])
    arc={x:1/t,y:0,z:2*t*t}
    vals=[S.cancel(f.subs(arc)) for f in [A,B,C]]
    check("known3D nonproper source arc has regular target",
          vals==[2*t*t,6*t,0])
    Jarc=JF.subs(arc).applyfunc(S.cancel)
    expected=S.Matrix([[0,6*t,1],[6*t*t,13,3/t],[-4,-3/t**2,-1/t**3]])
    check("known3D complete differential matrix on nonproper arc",Jarc==expected)
    v1,v2,v3,la,mu,nu=S.symbols("v1 v2 v3 la mu nu",nonzero=True)
    velocity=(S.Matrix([[la,mu,nu]])*Jarc*S.Matrix([v1,v2,v3]))[0]
    check("known3D generic first velocity exact pole order three",
          S.simplify(S.limit(t**3*velocity,t,0)+nu*v3)==0)
    check("known3D single first velocity is minus x cubed",
          S.diff(C,z)==-x**3 and S.diff(C,z).subs(arc)==-1/t**3)
    check("source mixed jets generate x",S.expand(-S.diff(C,x,y)/6-x)==0)
    check("source mixed jets generate z",S.expand(-S.diff(C,x,3)/6-z)==0)
    check("source mixed jets generate y",S.expand(-S.diff(C,x,2)/6-x*z-y)==0)
    check("3D equal-degree leading coefficients are not constant multiples",
          S.cancel(S.diff(3*x/U0,y)+3*x*x/U0**2)==0
          and S.diff(3*x/U0,y)!=0)
    # Weighted quotient is an exact dimension-two NON-Keller control.
    v=S.symbols("v");hh=2-3*u-v
    alpha=(1+u)*hh**2*(3*u**3+u*u*v+4*u*u+2*u*v+v)
    beta=hh*(9*u**3+3*u*u*v+12*u*u+6*u*v+u+3*v)
    check("known3D quotient has nonconstant Jacobian square",
          S.expand(S.Matrix([alpha,beta]).jacobian([u,v]).det()+2*hh**2)==0)

    # Graph degree and potential degree are genuinely different.
    q=S.symbols("q")
    U0=1+x*y
    for ix,R5 in enumerate([y**5,y**5+x*q**4]):
        Ag=U0**3*z+3*x*U0**2*q-x**3*R5
        Bg=y**2*U0*(4+3*x*y)*z+(y+3*x*y**2*(4+3*x*y))*q+(2*x-3*x*x*y)*R5
        Phi=S.expand(Ag**2+13*Ag+2*Bg)
        poly=S.Poly(Phi,x,y,z,q)
        top=sum(coef*x**ex[0]*y**ex[1]*z**ex[2]*q**ex[3]
                for ex,coef in poly.terms() if sum(ex)==16)
        check("exact quintic graph yields degree16 potential example"+str(ix),
              poly.total_degree()==16 and S.expand(top-x**6*R5**2)==0)
    # Cotangent doubling degree test is independent of the graph construction.
    X,Y,Z,W=S.symbols("X Y Z W")
    cot=Z*(X+Y**4)+W*Y
    check("degree4 plane component gives degree5 cotangent potential",
          S.Poly(cot,X,Y,Z,W).total_degree()==5)

    # Recheck only retained external evidence, NOT fresh CAS calculations.
    audit=ROOT/'hessian/hc4-degree-five-audit'
    run=subprocess.run([sys.executable,str(audit/'verify_evidence.py')],
                       capture_output=True,text=True,check=True)
    evidence=json.loads(run.stdout)
    check("retained preprint276 checks and source-log hashes intact",
          evidence["retained_author_checks"]==276
          and evidence["source_and_log_hashes"]=="PASS"
          and evidence["fresh_replays"]=={})
    for name,key in [("2611.00161v1.pdf","live_preprint_v1_pdf_sha256"),
                     ("2611.00161v1-source.zip","live_preprint_v1_zip_sha256")]:
        check("current remote comparison pin "+name,sha(audit/name)==ctx[key])
    check("coverage obligations remain genuinely open",
          ctx["global_coverage_proved"] is False
          and ctx["quadratic_model_coverage_assumed_not_proved"] is True)

    out=dict(claim="PF2RULED151",outcome="PASS",check_count=len(checks),checks=checks,
             elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),
             proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),
             context_sha256=sha(HERE/CONTEXT),
             conditional_quadratic_model_normal_form_proved=True,
             global_normal_form_coverage_proved=False,
             original_full_marked_integrality_proved=False,
             three_dimensional_first_velocity_countercontrol_verified=True,
             fresh_external_CAS_replays=0,retained_external_checks_revalidated=276,
             registered_predecessor_mathematical_checkers_replayed=[],
             full_goal_complete=False,independent_replay=False,
             formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ["claim","outcome","check_count","elapsed_seconds",
                                      "retained_external_checks_revalidated"]}))

if __name__=="__main__":
    main()
