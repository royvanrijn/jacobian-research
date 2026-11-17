"""PF144 exact controls for the written fibre-action exclusion."""
from pathlib import Path
from math import gcd
import hashlib
import json
import platform
import time
import sympy as S

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
NOTE="CYCLIC_FIBRE_SCALING_AND_TARGET_FIELD_ESCAPE.md"
CONTEXT="sources/cyclic-fibre-field-escape/context.json"
RECEIPT="cyclic_fibre_field_escape_receipt.json"

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

    t,u,x,v,z,eps,gamma,sig,aa=S.symbols("t u x v z eps gamma sig aa",nonzero=True)
    # The two actions commute, but their effects on the field E differ.
    original=S.Matrix([t,v,z])
    deck=S.Matrix([eps*t,eps*v,z/eps**2])
    fibre=S.Matrix([t,gamma*v,z])
    left=deck.subs({t:fibre[0],v:fibre[1],z:fibre[2]},simultaneous=True)
    right=fibre.subs({t:deck[0],v:deck[1],z:deck[2]},simultaneous=True)
    check("regular fibre action commutes with free deck action",left==right)
    check("local normal coordinate is deck invariant",S.cancel(deck[0]/deck[1]-t/v)==0)
    check("local tangent coordinate is deck invariant",S.cancel(deck[1]**2*deck[2]-v**2*z)==0)
    check("local normal fibre character minus one",S.cancel(fibre[0]/fibre[1]-gamma**-1*t/v)==0)
    check("local tangent fibre character two",S.cancel(fibre[1]**2*fibre[2]-gamma**2*v**2*z)==0)
    check("fibre action scales volume",S.Matrix([gamma*v,t]).jacobian([v,t]).det()==gamma)
    check("fixed equality first coordinate",S.expand((fibre-deck)[0]-(1-eps)*t)==0)
    check("fixed equality second coordinate",S.expand((fibre-deck)[1]-(gamma-eps)*v)==0)
    check("fixed equality third coordinate",S.cancel((fibre-deck)[2]-(1-eps**-2)*z)==0)
    check("boundary tangent chart Jacobian unit",
          S.det(S.Matrix([[1/v,0],[0,v**2]]))==v)
    check("Gamma chart Jacobian unit off t0",
          S.Matrix([t**5,v/t]).jacobian([t,v]).det()==5*t**3)

    # Character arithmetic checks prime and composite orders. Smith is used
    # only after choosing a PRIME order subgroup in the written proof.
    for order in [2,3,4,5,6,8,9,15]:
        fixed_normal=int((-1)%order==0)
        fixed_tangent=int(2%order==0)
        check(f"boundary normal never fixed order{order}",fixed_normal==0)
        check(f"boundary fixed dimension order{order}",fixed_tangent==int(order==2))
        check(f"Gamma fixed dimension one order{order}",
              int(0%order==0)+int(1%order==0)==1)
        check(f"volume character consistent order{order}",(-1+2-1)%order==0)
    for p in [3,5,7]:
        check(f"odd prime has isolated boundary fixed point p{p}",(-1)%p!=0 and 2%p!=0)
    check("involution fixes whole boundary line",2%2==0 and (-1)%2!=0)

    for g,k,bu in [(2,2,u-1),(3,2,u**2+u+1),(4,3,u+2),(5,2,u-1)]:
        sigmau=-bu+g*(1-k)*u*S.diff(bu,u)
        sigma=sigmau.subs(u,t**g)
        bt=bu.subs(u,t**g)
        a=t**2*bt**k
        check(f"g{g} all bad fibres avoid sigma zeros",S.gcd(sigmau,u*bu)==1)
        check(f"g{g} t0 has nonzero v labels",sigma.subs(t,0)!=0)
        check(f"g{g} deck invariance of a",all(m[0]%g==2%g for m in S.Poly(a,t).monoms()))
        check(f"g{g} deck invariance of sigma",all(m[0]%g==0 for m in S.Poly(sigma,t).monoms()))
        vv=sigma+a*x
        zz=sum(S.binomial(g,j)*sigma**(g-j)*a**(j-1)*x**j for j in range(1,g+1))
        check(f"g{g} selected line maps isomorphically to C0",
              S.expand((vv**2*zz).subs(t,0)-g*sigma.subs(t,0)**(g+1)*x)==0)
        check(f"g{g} source line and isolated fixed point match",zz.subs({t:0,x:0})==0)
        # rho changes the label; it is not a regular source action at t0.
        H=S.cancel(sigma/a)
        check(f"g{g} rational source action has genuine double pole",
              S.limit(t**2*H,t,0)==-bu.subs(u,0)**(1-k))
        check(f"g{g} fibre scaling slice action",
              S.cancel((gamma*v)/(g*t**(g+1)*bt**k)
                       -gamma*(v/(g*t**(g+1)*bt**k)))==0)
        if g==2:
            quotientJ=S.Matrix([t**2,zz]).jacobian([x,t]).det()
            check("quadratic coarse quotient has actual critical factors",
                  S.expand(quotientJ+4*t*vv)==0)

    # The canonical-trivial classical pseudo-plane family, all-degree
    # formulas supported by three fixed transcriptions.
    h=S.symbols("h")
    for g in [2,3,5]:
        uu=t**g;ww=t*(1+t**g*h)
        vv=sum(S.binomial(g,j)*t**(g*(j-1))*h**j for j in range(1,g+1))
        check(f"classical g{g} whole plane chart lands",S.expand(uu+uu**2*vv-ww**g)==0)
        check(f"classical g{g} unit volume on entire plane",
              S.expand(S.Matrix([uu,ww]).jacobian([t,h]).det()-g*uu**2)==0)
        check(f"classical g{g} selected line isomorphism",vv.subs(t,0)==g*h)
    check("classical fixed scheme keeps the whole A1 component",S.expand(u+u**2*v-u*(1+u*v))==0)

    # General symbolic identities, not evidence of the global Smith/line inputs.
    for degree in [2,3,7]:
        f=t**degree
        check(f"polynomial line degree{degree} is not etale",
              S.diff(f,t).subs(t,0)==0)
    check("line affine coordinate is etale",S.diff(3*t+2,t)==3)

    result=dict(claim="PF2RULED144",outcome="PASS",check_count=len(checks),
                checks=checks,elapsed_seconds=time.perf_counter()-started,
                environment=dict(python=platform.python_version(),sympy=S.__version__),
                proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),
                context_sha256=sha(HERE/CONTEXT),
                target_field_invariance_automatic=False,
                all_cyclic_targets_excluded=False,full_goal_complete=False,
                registered_predecessor_checkers_replayed=[],
                independent_replay=False,formal_verification=False,external_review=False)
    (HERE/RECEIPT).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({key:result[key] for key in ["claim","outcome","check_count","elapsed_seconds"]}))

if __name__=="__main__":
    main()
