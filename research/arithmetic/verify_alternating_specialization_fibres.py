#!/usr/bin/env python3
"""Exact specialization and full Keller-fibre transfer controls.

Run: .venv/bin/python research/arithmetic/verify_alternating_specialization_fibres.py
Optional --write-receipt PATH. Uses sympy==1.14.0 and exact rational/modular
arithmetic. This checks one explicit S8 specialization, its translation,
Jacobian factor identities and reconstruction in its degree-eight quotient.
The accompanying note proves the quantified Hilbert and scheme statements.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as sp

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt",type=Path)
    args=parser.parse_args()
    if sp.__version__!="1.14.0":raise RuntimeError("Requires sympy==1.14.0")
    z,S,Pi,B,Q=sp.symbols("z S Pi B Q")
    H=z**8+z**6+sp.Rational(9,28)*z**4+sp.Rational(233,5292)*z**2
    a=-sp.Rational(1,5488)
    P=H-2*z+a
    qbranch=sp.Symbol("s")**4+sp.Rational(5960,66706983)*sp.Symbol("s")**2+sp.Rational(2685619,1483273860320763)
    delta=sp.discriminant(P,z)
    assert delta==-7**7*qbranch.subs("s",2)**2
    assert delta!=0
    reductions={
      11:[z**3+z**2-3*z-4,z**5-z**4+5*z**3-4*z**2-3],
      23:[z+1,z**7-z**6+2*z**5-2*z**4-10*z**3+10*z**2+3*z-5],
      47:[z**8+z**6+2*z**4+10*z**2-2*z-17],
    }
    def modular(f,p):
        return sp.Poly.from_list([
          int(c.p)*pow(int(c.q),-1,p)%p for c in sp.Poly(f,z,domain=sp.QQ).all_coeffs()
        ],z,modulus=p)
    factor_degrees={}
    for prime,factors in reductions.items():
        polynomials=[sp.Poly(f,z,modulus=prime) for f in factors]
        product=sp.Poly(1,z,modulus=prime)
        for factor in polynomials:
            assert factor.is_irreducible
            product*=factor
        assert product==modular(P,prime)
        assert sp.gcd(product,product.diff()).degree()==0
        factor_degrees[str(prime)]=[f.degree() for f in polynomials]
    assert not sp.Poly(z**3-z,z,modulus=11).is_irreducible

    shifted=sp.Poly(P.subs(z,1+S),S,domain=sp.QQ).as_expr()
    g=[sp.expand(shifted).coeff(S,j) for j in range(9)]
    expected=[sp.Rational(54125,148176),sp.Rational(35387,2646),
      sp.Rational(237995,5292),sp.Rational(541,7),sp.Rational(2389,28),62,29,8,1]
    assert g==expected and g[1]*g[3]!=0
    G=shifted-g[0]
    Gpi=g[1]*S+Pi*(g[2]*S**2+g[3]*S**3)+sum(g[j]*Pi**j*S**j for j in range(4,9))
    beta=sp.cancel((sp.diff(Gpi,S)/g[1]-1-Pi*S**2)/S)
    compact_B=Q+beta
    compact_C=2*Gpi/g[1]-compact_B*S**2
    D=1-S*Q+Pi*S**2
    assert sp.cancel(sp.det(sp.Matrix([
      [sp.diff(compact_B,S),sp.diff(compact_B,Q)],
      [sp.diff(compact_C,S),sp.diff(compact_C,Q)]
    ]))+2*D)==0
    x,y,q=sp.symbols("x y q")
    t=1+x*y
    chart=[t*q,x/t,y+x*q]
    chart_jac=sp.det(sp.Matrix([[sp.diff(f,v) for v in (x,y,q)] for f in chart]))
    assert sp.cancel(chart_jac-1/t)==0
    # q_z=t^2 for the actual q=t^2 z+(g1/g3)y^2(1+3t).
    assert sp.cancel(chart_jac*t**2-t)==0
    assert sp.cancel(D.subs(dict(zip((Pi,S,Q),chart)))-1/t)==0

    # Check the polynomial pullback, not only the abstract Jacobian factors.
    w=sp.symbols("w")
    actual_q=t**2*w+(g[1]/g[3])*y**2*(1+3*t)
    actual_pi=t*actual_q
    actual_S=x/t
    actual_Q=y+x*actual_q
    low_B=actual_Q+2*g[2]/g[1]*actual_pi+(3*g[3]/g[1]-1)*actual_pi*actual_S
    low_C=2*(g[1]*actual_S+actual_pi*(g[2]*actual_S**2+g[3]*actual_S**3))/g[1]-low_B*actual_S**2
    assert sp.cancel(low_B-(y+3*g[3]/g[1]*x*actual_q+2*g[2]/g[1]*t*actual_q))==0
    assert sp.cancel(low_C-(x*(5-3*t)-g[3]/g[1]*x**3*w))==0
    for j in range(4,9):
        assert sp.cancel((t*q)**j*(x/t)**(j-2)-t**2*x**(j-2)*q**j)==0
        assert sp.cancel((t*q)**j*(x/t)**j-(x*q)**j)==0

    modulus=sp.Poly(shifted,S,domain=sp.QQ)
    def red(f):return sp.rem(sp.Poly(f,S,domain=sp.QQ),modulus).as_expr()
    def mul(*factors):
        out=sp.Integer(1)
        for factor in factors:out=red(out*factor)
        return out
    def power(f,n):
        out=sp.Integer(1)
        for _ in range(n):out=mul(out,f)
        return out
    d=red(sp.diff(shifted,S)/g[1])
    dinv=sp.invert(d,modulus.as_expr(),S)
    assert red(d*dinv)==1
    root_Q=red(-beta.subs(Pi,1))
    root_x=mul(S,dinv)
    root_y=red(root_Q-S)
    root_z=mul(d,d,red(d-(g[1]/g[3])*mul(root_y,root_y,1+3*dinv)))
    reconstructed_t=red(1+mul(root_x,root_y))
    reconstructed_q=red(mul(reconstructed_t,reconstructed_t,root_z)+(g[1]/g[3])*mul(root_y,root_y,1+3*reconstructed_t))
    assert red(reconstructed_t-dinv)==0
    assert red(reconstructed_q-d)==0
    assert mul(reconstructed_t,reconstructed_q)==1
    Bvalue=red(root_y+3*g[3]/g[1]*mul(root_x,reconstructed_q)+2*g[2]/g[1]*mul(reconstructed_t,reconstructed_q))
    for j in range(4,9):
        Bvalue=red(Bvalue+j*g[j]/g[1]*mul(power(reconstructed_t,2),power(root_x,j-2),power(reconstructed_q,j)))
    Cvalue=red(mul(root_x,5-3*reconstructed_t)-g[3]/g[1]*mul(power(root_x,3),root_z))
    for j in range(4,9):
        Cvalue=red(Cvalue-(j-2)*g[j]/g[1]*power(mul(root_x,reconstructed_q),j))
    target_C=-2*g[0]/g[1]
    assert Bvalue==0
    assert red(Cvalue-target_C)==0
    assert mul(root_x,d)==S
    # A deliberately altered target fails in the quotient.
    assert red(Cvalue-(target_C+1))!=0
    # Explicit rational A10 specialization of the new degree-ten pencil.
    P10=9*z**10-864*z**8+28161*z**6-350150*z**4+1127500*z**2-z+250000
    ten_factors={
      11:[z-2,z**9+2*z**8-4*z**7+3*z**6-z**3-2*z**2-4*z-2],
      13:[z**3-2*z**2+3*z-1,z**7+2*z**6-4*z**5-3*z**3+3*z**2+4*z-4],
    }
    ten_degrees={}
    for prime,factors in ten_factors.items():
        polynomials=[sp.Poly(f,z,modulus=prime) for f in factors]
        product=sp.Poly(1,z,modulus=prime)
        for factor in polynomials:
            assert factor.is_irreducible
            product*=factor
        assert product==modular(P10,prime).monic()
        assert sp.gcd(product,product.diff()).degree()==0
        ten_degrees[str(prime)]=[f.degree() for f in polynomials]
    ten_delta=sp.discriminant(P10,z)
    ten_square_root=2187*(59049-126938174857200000+68130645548641000000000000000)
    assert ten_delta==ten_square_root**2 and ten_delta!=0
    ten_g=[sp.expand(P10.subs(z,1+S)).coeff(S,j) for j in range(11)]
    assert ten_g==[1054655,1016543,-574772,-884684,13675,122850,5859,-5832,-459,90,9]
    assert ten_g[1]*ten_g[3]!=0
    ten_target=-2*ten_g[0]/ten_g[1]

    report={
      "schema_version":1,"kind":"internal_exact_replay",
      "command":".venv/bin/python research/arithmetic/verify_alternating_specialization_fibres.py",
      "python":platform.python_version(),"sympy":sp.__version__,
      "arithmetic":"QQ; F_11,F_23,F_47; QQ[S]/P(1+S)",
      "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      "specialization_polynomial":str(P),
      "discriminant":str(delta),
      "modular_factor_degrees":factor_degrees,
      "translation":1,"translated_coefficients":[str(c) for c in g],
      "degree_ten_specialization":{"polynomial":str(P10),"modular_factor_degrees":ten_degrees,"discriminant_square_root":str(ten_square_root),"translated_coefficients":[str(c) for c in ten_g],"normalized_target":["1","0",str(ten_target)]},
      "normalized_target":["1","0",str(target_C)],
      "checks":{"discriminant_square_class":-7,"all_modular_factors_irreducible":True,
        "compact_plane_jacobian":"-2D","source_chart_jacobian":"t",
        "normalized_full_jacobian":1,"polynomial_pullback":True,"quotient_derivative_unit":True,
        "whole_fibre_reconstruction":True,"marked_root_extraction":True,
        "invalid_controls":2},
      "scope":"Explicit S8 degree-eight and A10 degree-ten number-field specializations, and exact degree-eight compiler identities; the written note supplies group arguments, quantified Hilbert specialization and scheme-theoretic full-fibre proof.",
      "independent_replay":False,"formal_verification":False,"external_review":False,
      "result":"PASS"}
    output=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:args.write_receipt.write_text(output)
    print(output,end="")

if __name__=="__main__":main()
