"""Bounded exact original-residue and moving-coordinate bridge controls."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

t,z,x,b,U,S,W,Y=sp.symbols("t z x b U S W Y")


def zero(value):
    assert sp.cancel(value)==0,value


def truncate(value,order):
    return sp.expand(sum(coefficient*t**power for (power,),coefficient in sp.Poly(sp.expand(value),t).terms()
                         if power<=order))


def original_unit(phase,order):
    answer=0
    for n in range(order+1):
        value=0
        for (a,c,j),coefficient in sp.Poly(sp.expand(phase**n),z,W,Y).terms():
            if a==c and j%2==0:
                value+=coefficient*sp.factorial(a)*sp.factorial2(j-1)
        answer+=value*t**n/sp.factorial(n)
    return sp.expand(answer)


def controls():
    # Nonrational beta: e=z, g=C=1+z^2, c=1+z^2, h=-C/2.
    ell=z*(1+z*z)
    c=1+z*z
    h=-c/2
    p=z-t*h
    A=t*ell+t*t*c*c/2
    beta=(sp.sqrt(1-t*t)-1)/t
    relation=t*b*b+2*b+t
    F=1+t*b+t*x/2
    A1=t*(1+(b+x)**2)*F
    zero(sp.rem(p.subs(z,b+x)-x*F,relation,b))
    zero(sp.rem(A.subs(z,b+x)-x*A1,relation,b))
    assert sp.simplify(relation.subs(b,beta))==0
    assert sp.simplify(sp.diff(beta,t)+(1+beta*beta)/(2+2*t*beta))==0
    zero(c*(-z*c-2*z*h))  # Rmov=0, u=a=0.
    order=5
    beta_jet=sp.series(beta,t,0,order+2).removeO()
    residue=0
    for j in range(order+1):
        inv=sum((-1)**k*sp.binomial(j+k,k)*(t*beta_jet+t*x/2)**k
                for k in range(order-j+1))
        coefficient=sp.expand((1+(beta_jet+x)**2)**j*inv).coeff(x,j)
        residue+=sp.factorial2(2*j-1)*(2*t)**j*coefficient
    residue=truncate(residue,order)
    phase=ell*W*W+h*W+c*W*Y
    zero(residue-original_unit(phase,order))
    for power in range(-6,7):
        translated=sp.series((x+beta_jet)**power,t,0,order+1).removeO().expand()
        zero(translated.coeff(x,-1)-(1 if power==-1 else 0))

    # Higher raw origin order: gamma=2 and U starts as t^2*S/2.
    beta2=-t/2
    H2=2*t*(beta2+x)**2
    a10=t*beta2**2
    assert sp.Poly(a10,t).terms()==[((3,),sp.Rational(1,4))]
    xlead=t**3*S*S/2
    ulead=truncate(2*(beta2+xlead)**2*S,2)
    zero(ulead-t*t*S/2)
    high_phase=z**3*W*W+z*W*Y-W/2
    high_expected=truncate(sum(sp.factorial2(2*j-1)*sp.binomial(2*j,j)*(-t*t)**j
                               for j in range(order+1)),order)
    zero(high_expected-original_unit(high_phase,order))

    # Exact moving-coordinate identity in the degree-one pure/safe model.
    # It is deliberately NOT within the degree>=2 volume exclusion.
    q=1+2*z/t
    phi=-z/(4*t)-sp.Rational(1,8)
    zs=t*(2*S*S-sp.Rational(1,2))
    us=2*S
    zS=sp.diff(zs,S)
    zt=sp.diff(zs,t)
    def evaluate(value):
        return sp.cancel(sp.sympify(value).subs({z:zs,U:us},simultaneous=True))
    def dt(value):
        return sp.diff(value,t)-z*sp.diff(value,U)/(t*t*U)
    def dz(value):
        return sp.diff(value,z)+sp.diff(value,U)/(t*U)
    zero(us**2-q.subs(z,zs))
    zero(phi.subs(z,zs)+S*S/2)
    for value in [1,z,U,1/U,1/(z+t/2),z/U**3]:
        pull=evaluate(value)
        zero(sp.diff(pull,t)-evaluate(dt(value))-evaluate(dz(value))*zt)
        left=sp.diff(pull*zS,t)-evaluate(dt(value)+sp.diff(phi,t)*value)*zS
        zero(left-(sp.diff(pull*zt,S)-S*pull*zt))
    left=sp.diff(zS,t)-evaluate(sp.diff(phi,t))*zS
    omitted=2*S*S
    assert sp.cancel(left-(sp.diff(omitted,S)-S*omitted))!=0
    zero(zS/(t*us)-2)
    pure_phase=z*W*W+W*Y-W/2
    zero(original_unit(pure_phase,order)-1)
    return {"status":"passed","original_Wick_orders":[0,order],
            "nonrational_beta":{"phase":"z(1+z^2)W^2-(1+z^2)W/2+(1+z^2)WY",
                                "beta":"(sqrt(1-t^2)-1)/t","actual_unit_jet":str(residue)},
            "translation_residue_monomial_powers":[-6,6],
            "higher_valuation_case":{"ell":"z^3","c":"z","beta":"-t/2","A1_origin":"t^3/4",
                                     "U_leading":"t^2*S/2","actual_unit_jet":str(high_expected)},
            "exact_moving_covariant_forms":6,"omitted_beta_derivative_rejected":True,
            "degree_one_scope_control":"bridge matches pure unit1; degree>=2 volume premise is absent",
            "scope":"Bounded exact controls only; actual residue invariance, all-iterates horizontality and finite-image exclusion are written arguments."}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=controls()
    here=Path(__file__).resolve()
    files=[here,here.with_name("GAUSSIAN_MOVING_COMMON_ROOT_SCALAR_EXCLUSION.md"),
           here.with_name("GAUSSIAN_AFFINE_ISOTROPIC_CHART.md"),
           here.with_name("GAUSSIAN_SIMPLE_ORIGIN_ALTERNATE_SCALAR_EXCLUSION.md")]
    receipt={"status":"passed","checked_on":"2026-11-17",
             "environment":{"python":platform.python_version(),"sympy":sp.__version__},
             "sha256":{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in files},
             "checks":result,"independent_replay":False,"formal_verification":False,"external_review":False}
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
