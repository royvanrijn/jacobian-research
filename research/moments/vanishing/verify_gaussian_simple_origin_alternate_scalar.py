"""Exact controls for the scalar finite-seed bridge, not all-marker agreement."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp

if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

t,z,U,S=sp.symbols("t z U S")


def zero(value):
    assert sp.cancel(value)==0,value


def laurent_gaussian(value):
    answer=0
    for term in sp.Add.make_args(sp.expand(value)):
        power=term.as_powers_dict().get(S,sp.Integer(0))
        assert power.is_Integer
        coefficient=sp.cancel(term/S**power)
        assert not coefficient.has(S)
        power=int(power)
        if power%2:
            continue
        moment=(sp.factorial2(power-1) if power>=0
                else (-1)**(-power//2)/sp.factorial2(-power-1))
        answer+=coefficient*moment
    return sp.cancel(answer)


def gaussian_jet(value,order):
    series=sp.series(value,t,0,order+1).removeO().expand()
    answer=0
    for term in sp.Add.make_args(series):
        power=term.as_powers_dict().get(t,sp.Integer(0))
        coefficient=sp.cancel(term/t**power)
        assert not coefficient.has(t)
        answer+=laurent_gaussian(coefficient)*t**power
    return sp.expand(answer)


def controls():
    # Fully rational exact model ell=z(1+z), h=c=a=u=k=0.
    q=2*z*(1+z)/t
    phi=-z/(4*t*(1+z))
    zs=2*t*S*S/(1-2*t*S*S)
    us=2*S/(1-2*t*S*S)
    zS=sp.diff(zs,S)
    zt=sp.diff(zs,t)
    def evaluate(value):
        return sp.cancel(sp.sympify(value).subs({z:zs,U:us},simultaneous=True))
    def dt(value):
        return sp.diff(value,t)-U*sp.diff(value,U)/(2*t)
    def dz(value):
        return sp.diff(value,z)+sp.diff(q,z)*sp.diff(value,U)/(2*U)
    zero(us**2-q.subs(z,zs))
    zero(phi.subs(z,zs)+S*S/2)
    for exponent in range(-10,11):
        value=S**exponent
        zero(laurent_gaussian(sp.diff(value,S)-S*value))
    forms=[1,z,U,1/U,1/z,z/U**3,t/((1+z)*U)]
    gauge=t+2*t*t
    for value in forms:
        pull=evaluate(value)
        zero(sp.diff(pull,t)-evaluate(dt(value))-evaluate(dz(value))*zt)
        zero(sp.diff(pull,S)-evaluate(dz(value))*zS)
        for psi in [sp.Integer(0),gauge]:
            phase=phi-psi
            Psi=-psi
            left=(sp.diff(pull*zS,t)+sp.diff(Psi,t)*pull*zS
                  -evaluate(dt(value)+sp.diff(phase,t)*value)*zS)
            right=sp.diff(pull*zt,S)-S*pull*zt
            zero(left-right)
    # Every tested covariant iterate is evaluated from the SOURCE connection,
    # rather than obtained by differentiating a presumed scalar output.
    coefficient=1/(t*U)
    for derivative in range(4):
        actual=gaussian_jet(evaluate(coefficient)*zS/2,4)
        expected=sum(sp.factorial(2*(m+derivative))*t**m/sp.factorial(m) for m in range(5))
        zero(actual-expected)
        coefficient=sp.cancel(dt(coefficient)+sp.diff(phi,t)*coefficient)
    # The localized Jker class is zero, while its ORIGINAL marked period is1.
    stationary_w=1/(2*t*(1+z))
    zero(1-2*t*(1+z)*stationary_w)
    zero(gaussian_jet(evaluate(stationary_w/(t*U))*zS/2,4)-1/(2*t))
    original_f=sum(sp.factorial(2*m)*t**m/sp.factorial(m) for m in range(6))
    original_RW=sum(sp.factorial(2*m+1)*t**m/sp.factorial(m) for m in range(5))
    zero(original_f-2*t*original_RW-1)
    return {"status":"passed","model":"ell=z(1+z), h=c=a=u=k=0",
            "fixed_Laurent_Stein_exponents":[-10,10],
            "exact_covariant_chain_rule_forms":len(forms),"scalar_gauges":["0","t+2t^2"],
            "source_covariant_iterates":[0,3],"scalar_Wick_orders":[0,4],
            "Gysin_scope_control":{"original_Jker":"1","alternate_localized_Jker":"0",
                                   "original_W":"0","alternate_localized_W":"1/(2t)"},
            "scope":"Exact bounded controls for the written all-iterates scalar bridge. No original all-marker or source-localization agreement is asserted."}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=controls()
    here=Path(__file__).resolve()
    proof=here.with_name("GAUSSIAN_SIMPLE_ORIGIN_ALTERNATE_SCALAR_EXCLUSION.md")
    source=here.with_name("GAUSSIAN_AFFINE_ISOTROPIC_CHART.md")
    receipt={"status":"passed","checked_on":"2026-11-17",
             "environment":{"python":platform.python_version(),"sympy":sp.__version__},
             "checker_sha256":hashlib.sha256(here.read_bytes()).hexdigest(),
             "proof_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
             "h21_source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
             "checks":result,"independent_replay":False,"formal_verification":False,"external_review":False}
    if args.output:
        args.output.write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":
    main()
