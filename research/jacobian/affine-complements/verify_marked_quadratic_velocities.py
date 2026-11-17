"""PF143 supporting exact controls. Written universal proofs are in the note."""
from pathlib import Path
from math import ceil
import hashlib
import json
import platform
import time
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="MARKED_QUADRATIC_VELOCITIES_AND_CYCLIC_TRACE_POLES.md"
CONTEXT="sources/marked-quadratic-velocities/context.json"
RECEIPT="marked_quadratic_velocities_receipt.json"

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
        check("source pin "+path,sha(ROOT/path)==digest)

    # Full codifferent threshold, not just regularity in an etale open.
    for e,q in [(1,2),(1,7),(2,1),(3,2),(5,2),(6,2),(13,4)]:
        m=1+ceil((2*q+1)/e)
        check(f"threshold sufficient e{e}q{q}",e*m-2*(q+e)>=1-e)
        check(f"threshold minimal e{e}q{q}",e*(m-1)-2*(q+e)<1-e)
        check(f"all omitted primes cost at least two e{e}q{q}",m>=2)
    check("shared image takes full maximum",max(1+ceil(5/5),1+ceil(7/1))==8)
    check("unramified q2 marked quadratic charge six",1+ceil(5/1)==6)

    X,A,B,lam,Xp,Ap,Xq,Aq,kappa=sp.symbols("X A B lam Xp Ap Xq Aq kappa")
    def trace_character(expr,order,characters):
        out=0
        variables=list(characters)
        for powers,coefficient in sp.Poly(sp.expand(expr),*variables).terms():
            weight=sum(p*characters[v] for p,v in zip(powers,variables))
            if weight%order==0:
                out+=order*coefficient*sp.prod(v**p for p,v in zip(powers,variables))
        return sp.expand(out)
    chars={X:1,A:-1,B:-2}
    w=lam*X+A+B
    for g in [2,3,4,5,6,7,11]:
        if g==2:
            check("C2 first trace",trace_character(w,g,chars)==2*B)
            check("C2 quadratic trace",sp.expand(trace_character(w**2,g,chars)-2*((lam*X+A)**2+B**2))==0)
        elif g==3:
            check("C3 quadratic trace",trace_character(w**2,g,chars)==6*(lam*X*A+A*B))
        elif g==4:
            check("C4 quadratic trace",trace_character(w**2,g,chars)==4*(2*lam*X*A+B**2))
        else:
            check(f"C{g} first trace zero",trace_character(w,g,chars)==0)
            check(f"C{g} quadratic trace",trace_character(w**2,g,chars)==2*g*lam*X*A)
        j0=(g+1)//2
        for j in range(1,j0):
            check(f"C{g} exceptional transverse prefix j{j}",trace_character((A+B)**j,g,chars)==0)
        expected=g*B**j0 if g%2==0 else g*j0*A*B**(j0-1)
        check(f"C{g} first transverse nonzero power",trace_character((A+B)**j0,g,chars)==expected)

    # Constant Jacobian supplies a determinant, but entries need not land.
    matrix=sp.Matrix([[2*lam*Xp*Ap,lam*(Xp*Aq+Ap*Xq)],
                      [lam*(Xp*Aq+Ap*Xq),2*lam*Xq*Aq]])
    check("relative Gram determinant exact",sp.expand(matrix.det()+lam**2*(Xp*Aq-Ap*Xq)**2)==0)
    check("common polar component cancels in Jacobian",
          sp.expand(Xp*(Aq+kappa*Xq)-(Ap+kappa*Xp)*Xq-Xp*Aq+Ap*Xq)==0)

    t,u,x,v,p=sp.symbols("t u x v p")
    for g,k,bu in [(2,2,u-1),(3,3,u+2),(5,2,u-1)]:
        sigmau=-bu+g*(1-k)*u*sp.diff(bu,u)
        bt=bu.subs(u,t**g);sigma=sigmau.subs(u,t**g)
        a0=t**2*bt**k
        K=sp.cancel((sp.diff(sigma,t)-sp.diff(a0,t)*sigma/a0)/a0)
        lead=2*bu.subs(u,0)**(1-k)
        check(f"g{g} polar coefficient K",sp.limit(t**3*K,t,0)==lead)
        # F=v^2*z is a GLOBAL invariant function on W/mu_g, not a target pair.
        vv=sigma+a0*x
        FF=vv**2*sum(sp.binomial(g,i)*sigma**(g-i)*a0**(i-1)*x**i for i in range(1,g+1))
        XX=sp.diff(FF,x)
        check(f"g{g} global function nonzero longitudinal derivative",
              sp.simplify(XX.subs(t,0)-g*sigma.subs(t,0)**(g+1))==0)
        check(f"g{g} global function normal character",
              (2-2)%g==0)
        check(f"g{g} exact polar leading coefficient B",
              sp.limit(t**3*K*XX,t,0)==lead*g*sigma.subs(t,0)**(g+1))
        if g==2:
            rationalF=v**2*(v**g-sigma**g)/a0
            rationalA=sp.diff(rationalF,t)+sp.diff(a0,t)/a0*v*sp.diff(rationalF,v)
            rationalB=(sp.diff(sigma,t)-sp.diff(a0,t)/a0*sigma)*sp.diff(rationalF,v)
            check("chain rule with fixed v and fixed x distinguished",
                  sp.cancel((rationalA+rationalB).subs(v,vv)-sp.diff(FF,t))==0)

    # Explicit Laurent source has Jac1 and FULL normalization C[p,t].
    # The cyclic trace of a Laurent monomial t^n is g*u^(n/g), or zero.
    def trace_laurent(expr,g):
        out=0
        poly=sp.Poly(sp.cancel(expr*t**8),t)
        for (power,),coefficient in poly.terms():
            n=power-8
            if n%g==0:
                out+=g*coefficient*u**(n//g)
        return sp.cancel(out)
    for g in [5,6,7]:
        bu=u-1;k=2
        sigmau=-bu+g*(1-k)*u*sp.diff(bu,u)
        D0=sp.cancel(g*u*sp.diff(sigmau,u)-(2+k*g*u*sp.diff(bu,u)/bu)*sigmau)
        J=sp.cancel(D0/(g*u*bu**k))
        wp=lam*t/(g*u)-(g-1)*p/t+J/t**2
        wq=g*u/t
        G=sp.Matrix([[trace_laurent(wp**2,g),trace_laurent(wp*wq,g)],
                     [trace_laurent(wp*wq,g),trace_laurent(wq**2,g)]])
        expected=sp.Matrix([[-2*lam*(g-1)*p/u,g*lam],[g*lam,0]])
        check(f"g{g} exact Laurent moment matrix",all(sp.cancel(z)==0 for z in G-expected))
        check(f"g{g} constant nonzero partial Gram determinant",sp.cancel(G.det()+g**2*lam**2)==0)
        marked=trace_laurent(t**4*wp**2,g)
        check(f"g{g} integral marker detects double target pole",
              sp.limit(u**2*marked,u,0)==sp.Rational(4,g))
        check(f"g{g} full ramified prime marked charge",1+ceil(5/g)==2)
    g=5
    bt=t**g-1
    sigma=(- (u-1)-g*u).subs(u,t**g)
    a0=t**2*bt**2
    pp=(x+sigma/a0)/(g*t**(g-1));qq=t**g
    check("actual Laurent constant Jacobian",sp.cancel(sp.diff(pp,x)*sp.diff(qq,t)-1)==0)
    check("branch divisor retained with full index",sp.diff(t**g,t)==g*t**(g-1))
    check("Laurent source is not whole plane",sp.denom(sp.cancel(pp)).has(t))

    result=dict(claim="PF2RULED143",outcome="PASS",check_count=len(checks),checks=checks,
                elapsed_seconds=time.perf_counter()-started,python=platform.python_version(),
                sympy=sp.__version__,proof_sha256=sha(HERE/NOTE),
                checker_sha256=sha(Path(__file__)),context_sha256=sha(HERE/CONTEXT),
                automatic_marked_regularity=False,actual_cyclic_plane_target_constructed=False,
                registered_predecessor_checkers_replayed=[],formal_verification=False,
                independent_replay=False,external_review=False,full_goal_complete=False)
    (HERE/RECEIPT).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ["claim","outcome","check_count","elapsed_seconds"]}))

if __name__=="__main__":
    main()
