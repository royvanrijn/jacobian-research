#!/usr/bin/env python3
"""Exact controls for QEM1: fixed Moyal defects and an elliptic cohomology identity."""
from pathlib import Path
from math import comb, factorial
import hashlib
import json
import sympy as s

HERE=Path(__file__).resolve().parent
x,y,b,h,eps,tau=s.symbols("x y b h eps tau")
u=4*y**2/x**2-3*x
v=8*y**3/x**3-9*y
H=y**2-x**3


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def zero(expr,label):
    if s.cancel(s.expand(expr)) != 0:
        raise AssertionError(label+": "+str(expr))


def lam(f,g,n):
    return s.expand(sum((-1)**i*comb(n,i)*s.diff(f,y,n-i,x,i)
                        *s.diff(g,x,n-i,y,i) for i in range(n+1)))


def pb(f,g):
    return lam(f,g,1)


def mul_jet(F,G,N):
    out={}
    for i,f in F.items():
        for j,g in G.items():
            for n in range(N-i-j+1):
                k=i+j+n
                out[k]=out.get(k,0)+lam(f,g,n)/(2**n*factorial(n))
    return {k:s.cancel(a) for k,a in out.items()}


def add_jet(F,G,scale=1):
    return {k:s.cancel(F.get(k,0)+scale*G.get(k,0)) for k in set(F)|set(G)}


def curve_zero(expr,label):
    num=s.together(expr).as_numer_denom()[0]
    remainder=s.Poly(s.expand(num),y).rem(s.Poly(y**2-x**3-b,y)).as_expr()
    zero(remainder,label)


def main():
    zero(pb(v,u)-27,"constant symplectic multiplier")
    zero(v**2-u**3+27*H,"classical isogeny energy")
    A=s.factor(lam(v,u,3)/24)
    Q=s.factor((lam(v,v,2)-u*lam(u,u,2)-lam(u**2,u,2))/8)
    zero(A-96*y**2/x**8,"commutator anomaly")
    zero(Q+6*(3*x**6+8*x**3*y**2-8*y**4)/x**8,"energy anomaly")

    # Check the two anomalies by a separate coefficient assembly of star products.
    UV=mul_jet({0:u},{0:v},3)
    VU=mul_jet({0:v},{0:u},3)
    comm=add_jet(VU,UV,-1)
    zero(comm[1]-27,"assembled first commutator")
    zero(comm[3]-A,"assembled third commutator")
    cube=mul_jet(mul_jet({0:u},{0:u},2),{0:u},2)
    energy=add_jet(mul_jet({0:v},{0:v},2),cube,-1)
    zero(energy[0]+27*H,"assembled classical energy")
    zero(energy[2]-Q,"assembled second energy")

    xi=s.Function("xi")(x,y)
    eta=s.Function("eta")(x,y)
    U2=xi*s.diff(u,x)+eta*s.diff(u,y)
    V2=xi*s.diff(v,x)+eta*s.diff(v,y)
    zero(pb(V2,u)+pb(v,U2)-27*(s.diff(xi,x)+s.diff(eta,y)),
         "universal divergence linearization")
    zero(2*v*V2-3*u**2*U2+27*(xi*s.diff(H,x)+eta*s.diff(H,y)),
         "universal energy linearization")
    fac=Q/(162*H)
    div=s.diff(2*x*fac,x)+s.diff(3*y*fac,y)
    F=s.factor(-A/27-div)
    Ffiber=-(56*b**2+56*b*x**3+15*x**6)/(27*b*x**8)
    curve_zero(F-Ffiber,"fiber obstruction coefficient")
    primitive=4*y/(27*x**7)+y/(18*b*x**4)+5*y/(36*b**2*x)
    curve_derivative=s.diff(primitive,x)+3*x**2*s.diff(primitive,y)/(2*y)
    curve_zero(Ffiber/(2*y)-curve_derivative+5*x/(72*b**2*y),
               "exact elliptic obstruction differential")
    qfun=s.Function("g")(b)
    factor=qfun/(162*b)
    # Euler(b)=6b and divergence(Euler)=5.
    zero(6*b*s.diff(factor,b)+5*factor-(6*s.diff(qfun,b)-qfun/b)/162,
         "arbitrary univariate energy correction")

    # One explicit unrestricted escape jet, not an order or support census.
    f=-32*y**3/(27*x**8)
    uc=s.factor(f*s.diff(u,y))
    vc=s.factor(f*s.diff(v,y))
    lifted_comm=add_jet(mul_jet({0:v,2:vc},{0:u,2:uc},3),
                        mul_jet({0:u,2:uc},{0:v,2:vc},3),-1)
    for order in range(4):
        zero(lifted_comm.get(order,0)-(27 if order==1 else 0),
             "unrestricted jet relation order "+str(order))
    lifted_energy=add_jet(mul_jet({0:v,2:vc},{0:v,2:vc},2),
        mul_jet(mul_jet({0:u,2:uc},{0:u,2:uc},2),{0:u,2:uc},2),-1)
    G2=-18/x**2-48*y**2/x**5+112*y**4/x**8
    zero(lifted_energy[2]-G2,"explicit transverse energy")
    energy_comm=s.factor(pb(G2,H))
    if energy_comm==0:
        raise AssertionError("escape must move energy")

    # Polynomial Ore relations: ad(H)/eps is regular on the PBW generators.
    zero(pb(H,x)-2*y,"Hamiltonian action on x")
    zero(pb(H,y)-3*x**2,"Hamiltonian action on y")
    delta_x=-y
    delta_y=-s.Rational(3,2)*x**2
    zero(delta_x*s.diff(H,x)+delta_y*s.diff(H,y),"delta preserves cubic")
    # The only Weyl-relation overlap for the derivation vanishes.
    zero(pb(delta_y,x)+pb(y,delta_x),"Ore derivation relation compatibility")
    zero(2*eps*tau*delta_x+tau*eps*pb(H,x),"central C on x")
    zero(2*eps*tau*delta_y+tau*eps*pb(H,y),"central C on y")
    zero(-s.Rational(1,2)*pb(H,x)-delta_x,"diagonal reduction x relation")
    zero(-s.Rational(1,2)*pb(H,y)-delta_y,"diagonal reduction y relation")

    # Bind the normalized rational map to the retained actual Long chart.
    r,pos,p=s.symbols("r pos p")
    aa=p/pos
    SS=2*pos+4*aa-3*r*aa**2
    TT=(r*aa**3-aa**2-pos*aa)/2
    xx=8*r*pos
    yy=8*r*(3*r*p-2*pos)
    subst={x:xx,y:yy}
    zero(u.subs(subst,simultaneous=True)-(16-12*r*SS),"actual Long u chart")
    zero(v.subs(subst,simultaneous=True)-(432*r**2*TT+72*r*SS-64),
         "actual Long v chart")
    zero(s.diff(yy,p)*s.diff(xx,pos)-s.diff(yy,pos)*s.diff(xx,p)-192*r**3,
         "actual Long Planck multiplier")

    for name in ["isogeny_module_controls.json","quantum_elliptic_energy_controls.json"]:
        receipt=json.loads((HERE/name).read_text())
        if receipt["status"]!="PASS":
            raise AssertionError("retained receipt status "+name)
        for label,expected in receipt["hashes"].items():
            path=HERE/label
            if not path.exists():
                path=HERE.parent/"dixmier-audit"/label
            if sha(path)!=expected:
                raise AssertionError("retained input drift "+label)

    locked=["ELLIPTIC_QUANTUM_ENERGY_MOTION.md","check_energy_motion.py",
            "ELLIPTIC_ISOGENY_MATRIX_MODULE.md","isogeny_module_controls.json",
            "ELLIPTIC_ISOGENY_AND_QUANTUM_ENERGY.md",
            "quantum_elliptic_energy_controls.json","PROOF.md"]
    out={"claim":"QEM1","status":"PASS",
         "normalized_isogeny_and_actual_Long_chart":"PASS",
         "Moyal_defects":{"third_commutator":str(A),"second_energy":str(Q)},
         "universal_vector_field_linearizations":"PASS",
         "relative_obstruction":{"nonexact_coefficient":"-5/(72*b^2)",
             "class":"[x dx/y]","exact_primitive":str(primitive),
             "certificate":"direct rational differential identity; nonexactness is a written genus-one proof"},
         "arbitrary_univariate_correction_linearization":"PASS",
         "explicit_unrestricted_escape":{"U2":str(uc),"V2":str(vc),"energy_G2":str(G2),
             "third_energy_commutator_coefficient":str(energy_comm),
             "commutator_through_h3":"PASS"},
         "Ore_bridge_central_and_diagonal_generator_relations":"PASS",
         "retained_QIM1_QEI1_hash_chains":"PASS",
         "hashes":{name:sha(HERE/name) for name in locked},
         "coverage_boundary":"Exact finite Moyal, relative-differential and Ore-generator identities only. Nonexactness, rational logarithmic normalization, the universal h-valuation bound2..3, all-degree canonical transport and PBW flatness use written proofs. No spatial bound for a full Long lift or all-representative nonlifting is certified.",
         "independent_replay":False,"formal_verification":False,"external_review":False}
    (HERE/"energy_motion_controls.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({"claim":"QEM1","status":"PASS",
        "obstruction":"-5/(72*b^2) [x dx/y]",
        "escape_energy_commutator":str(energy_comm),
        "receipt":"energy_motion_controls.json"},indent=2))


if __name__=="__main__":
    main()
