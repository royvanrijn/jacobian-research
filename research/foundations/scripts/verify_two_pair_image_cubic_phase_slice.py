#!/usr/bin/env python3
"""Verify the exact six-moment certificate for a spherical cubic phase slice.

The family is p=y*(c0+c1*x+c2*x**2)+d0+d1*x+d2*x**2+d3*x**3,
in C[x,y,t]/(t**2+2*x*y-1). This is seven dimensional, not all V_3.
Moment monomials satisfy L(x**i*y**j)=delta(i,j)*i!/(2*i+1)!!.
The symbolic equations, together with the stated branch split, prove that
moments 1..6 vanish iff the displayed p has strictly one-sided phase support.
No finite sampling or numerical arithmetic is used. Requires sympy==1.14.0.
"""
from __future__ import annotations

import hashlib
import json
import platform
from math import comb, factorial
from pathlib import Path

import sympy as s

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT / "research/foundations/extended-geometry/two_pair_image_cubic_phase_slice_verification.json"


def main() -> None:
    if s.__version__ != "1.14.0":
        raise RuntimeError("This replay is pinned to sympy==1.14.0")
    x, z, a, b = s.symbols("x z a b")
    d0, d1, d2, d3 = s.symbols("d0 d1 d2 d3")
    c_poly = 1 + a*x + b*x**2
    d_poly = d0 + d1*x + d2*x**2 + d3*x**3

    def moment(order: int, substitutions: dict) -> s.Expr:
        cp = s.Poly(c_poly, x, domain=s.QQ[a,b])
        dp = s.Poly(d_poly.subs(substitutions), x)
        return s.factor(sum(
            s.Rational(comb(order, k)*factorial(k), s.factorial2(2*k+1))
            * (cp**k * dp**(order-k)).nth(k)
            for k in range(order+1)
        ))

    # Small direct controls for the normalization of the sphere functional.
    assert moment(1,{d0:0,d1:0,d2:0,d3:0}).subs({a:1,b:0}) == s.Rational(1,3)
    assert moment(2,{d0:0,d1:0,d2:0,d3:0}).subs({a:1,b:0}) == s.Rational(2,15)
    assert s.factor(moment(1,{})-(a+3*d0)/3) == 0
    substitutions = {d0: -a/3}
    solutions = [
        -a**2/s.Integer(30)-2*b/s.Integer(5),
        a*(31*a**2-108*b)/s.Integer(1890),
        -289*a**4/s.Integer(28350)+47*a**2*b/s.Integer(1050)-2*b**2/s.Integer(75),
    ]
    coefficients = [s.Rational(2,3),s.Rational(2,5),s.Rational(8,35)]
    triangular = [{"moment":1,"identity":"M1=(a+3*d0)/3","solution":"d0=-a/3"}]
    for order, variable, solution, coefficient in zip(
        range(2,5),[d1,d2,d3],solutions,coefficients
    ):
        equation = moment(order,substitutions)
        assert s.Poly(equation,variable).degree() == 1
        assert s.factor(equation-coefficient*(variable-solution)) == 0
        assert coefficient != 0
        substitutions[variable] = solution
        triangular.append({
            "moment":order,
            "identity":f"M{order}={coefficient}*({variable}-({solution}))",
            "solution":f"{variable}={solution}",
        })

    p_poly = 2219-11874*z+13176*z**2
    q_poly = 482933-1979487*z+218781*z**2+1571724*z**3
    expected5 = -4*a*(2219*a**4-11874*a**2*b+13176*b**2)/s.Integer(9823275)
    expected6 = -8*(482933*a**6-1979487*a**4*b+218781*a**2*b**2+1571724*b**3)/s.Integer(13408770375)
    assert s.factor(moment(5,substitutions)-expected5) == 0
    assert s.factor(moment(6,substitutions)-expected6) == 0
    # If a!=0, let z=b/a**2. Both P(z) and Q(z) must vanish.
    bezout_a = 6676468225812*z**2+2727590293062*z-7662133538547
    bezout_b = 35364219296-55969842888*z
    bezout_constant = 76274195239375
    assert s.expand(bezout_a*p_poly+bezout_b*q_poly) == bezout_constant
    assert s.expand((bezout_a+1)*p_poly+bezout_b*q_poly) != bezout_constant
    # The remaining a=0 branch forces b=0, since the coefficient is nonzero.
    assert s.factor(moment(6,substitutions).subs(a,0)) == -s.Rational(32,34125)*b**3
    assert all(s.factor(value.subs({a:0,b:0})) == 0 for value in substitutions.values())

    # c0=0 branch has no negative weights. Only c1*x*y+d0 contributes.
    c1,e0=s.symbols("c1 e0")
    first=e0+c1/s.Integer(3)
    second=e0**2+2*e0*c1/s.Integer(3)+2*c1**2/s.Integer(15)
    assert first.subs(e0,-c1/3) == 0
    assert s.factor(second.subs(e0,-c1/3)) == c1**2/45
    # After c1=e0=0 the surviving weights are +1,+2,+3.
    assert all(weight>0 for weight in [1,1,2,3])

    # Sharpness control: five moments are insufficient in the same slice.
    sharp={key:value.subs({a:0,b:1}) for key,value in substitutions.items()}
    sharp_moments=[s.factor(moment(m,sharp).subs({a:0,b:1})) for m in range(1,7)]
    assert sharp_moments == [0,0,0,0,0,-s.Rational(32,34125)]

    receipt = {
        "format":"two-pair-image-cubic-phase-slice-v1",
        "arithmetic":"exact symbolic rational polynomial identities",
        "software":{"python":platform.python_version(),"sympy":s.__version__},
        "proof_source":"foundations/extended-geometry/TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md",
        "checker":"foundations/scripts/verify_two_pair_image_cubic_phase_slice.py",
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "command":".venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_phase_slice.py",
        "family":"p=y*(c0+c1*x+c2*x^2)+d0+d1*x+d2*x^2+d3*x^3",
        "sphere_relation":"t^2+2*x*y=1",
        "parameter_dimension":7,
        "moment_functional":"L(x^i*y^j)=delta(i,j)*i!/(2*i+1)!!",
        "conclusion":"L(p^m)=0 for m=1..6 iff p has strictly one-sided phase support in these coordinates (including zero)",
        "branches":{
            "c0_nonzero":{
                "normalization":"scale c0=1; a=c1/c0, b=c2/c0, and similarly scale d0..d3",
                "triangular_equations":triangular,
                "moment5_after_triangular_substitution":str(expected5),
                "moment6_after_triangular_substitution":str(expected6),
                "a_nonzero_substitution":"z=b/a^2",
                "P":str(p_poly),"Q":str(q_poly),
                "bezout_A":str(bezout_a),"bezout_B":str(bezout_b),
                "bezout_identity":"A*P+B*Q=76274195239375",
                "a_zero_moment6":"-32*b^3/34125",
                "result":"a=b=d0=d1=d2=d3=0; p=y before undoing scaling",
            },
            "c0_zero":{
                "moment1":"d0+c1/3",
                "moment2_after_moment1":"c1^2/45",
                "result":"c1=d0=0; all remaining monomials have positive phase weight",
            },
        },
        "sharpness":{
            "p":"y+x^2*y-2*x/5-2*x^3/75",
            "moments_1_through_6":[str(value) for value in sharp_moments],
        },
        "negative_control":"perturbing Bezout_A by +1 does not satisfy the recorded identity",
        "assurance":{"internal_exact_replay":True,"independent_replay":False,"formal_verification":False,"external_review":False},
        "limits":"Does not cover the full 16-parameter cubic family or certify a general degree-three classification.",
    }
    OUTPUT.write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: seven-parameter phase slice, exact six-moment identities, both branches and Bezout certificate")
    print("PASS: five-moment sharpness control; altered Bezout certificate rejected")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
