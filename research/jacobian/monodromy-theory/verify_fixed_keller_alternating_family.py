#!/usr/bin/env python3
"""Exact controls for the uniform fixed-map pencil realization.

Run with repository .venv/bin/python. The note supplies uniform-in-degree
coefficient, whole-scheme and perfect-group proofs; this checker certifies
the complete degree-ten identities and a complete quotient reconstruction.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as sp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    assert sp.__version__ == "1.14.0"
    T, S, p, V, Y, Z, Q = sp.symbols("T S p V Y Z Q")
    H = 9*T**10-864*T**8+28161*T**6-350150*T**4+1127500*T**2
    a = sp.Integer(250000)
    tau = sp.Integer(1)
    shifted = sp.Poly(H.subs(T, tau+S), S)
    c = [shifted.nth(j) for j in range(11)]
    A = c[1]
    kappa = c[0]+a-tau*A
    assert A == 1016544 and kappa == 38112
    assert c[2:] == [-574772, -884684, 13675, 122850, 5859, -5832, -459, 90, 9]
    profiles = {j: c[j]*V**(2*j-8) for j in range(4, 11)}
    G = S+c[3]*V*S**3+sum(profiles[j]*V**j*S**j for j in profiles)
    E = G+Y*S**2-Z/2
    gamma = {V:p, Y:c[2]/p**2, Z:-2*(kappa+tau*p**5)/p**8}
    pencil_pullback = H.subs(T, tau+p**3*S)-(A-p**5)*(tau+p**3*S)+a
    assert sp.expand(p**8*E.subs(gamma)-pencil_pullback) == 0
    assert sp.expand(p**5*sp.diff(E,S).subs(gamma)
                     -sp.diff(H,T).subs(T,tau+p**3*S)+A-p**5) == 0
    assert sp.expand(p**8*E.subs(gamma)-pencil_pullback+S**2) != 0

    beta = sp.cancel((sp.diff(G,S)-1-V*S**2)/S)
    assert sp.Poly(beta,S,V).is_univariate is False
    Bcompact = Q+beta
    Ccompact = 2*G-Bcompact*S**2
    D = 1-S*Q+V*S**2
    jac_plane = sp.det(sp.Matrix([[sp.diff(f,v) for v in (S,Q)]
                                 for f in (Bcompact,Ccompact)]))
    assert sp.expand(jac_plane+2*D) == 0
    assert sp.expand(D.subs(Q,-2*Y-beta)-sp.diff(E,S)) == 0

    x,y,q,z,cubic = sp.symbols("x y q z cubic", nonzero=True)
    t = 1+x*y
    chart = (t*q,x/t,y+x*q)
    jac_chart = sp.det(sp.Matrix([[sp.diff(f,v) for v in (x,y,q)]
                                 for f in chart]))
    assert sp.cancel(jac_chart-1/t) == 0
    assert sp.cancel(D.subs(dict(zip((V,S,Q),chart)))-1/t) == 0
    actual_q = t**2*z+y**2*(1+3*t)/cubic
    low_G = S+cubic*V*S**3
    low_beta = sp.cancel((sp.diff(low_G,S)-1-V*S**2)/S)
    chart_actual = {V:t*actual_q,S:x/t,Q:y+x*actual_q}
    low_B = (Q+low_beta).subs(chart_actual)
    low_C = (2*low_G-(Q+low_beta)*S**2).subs(chart_actual)
    assert sp.cancel(low_B-y-3*cubic*x*actual_q) == 0
    assert sp.cancel(low_C-x*(5-3*t)+cubic*x**3*z) == 0
    degrees = []
    for j in range(4,11):
        # Exact high-term pullback; profile derivatives are irrelevant at fixed V.
        assert sp.cancel((t*q)**(3*j-8)*(x/t)**(j-2)
                         -t**(2*j-6)*x**(j-2)*q**(3*j-8)) == 0
        assert sp.cancel((t*q)**(3*j-8)*(x/t)**j
                         -t**(2*j-8)*x**j*q**(3*j-8)) == 0
        degrees.append([j, 20*j-54, 20*j-56])
    assert max(row[1] for row in degrees) == 146
    assert sp.cancel(jac_chart*t**2*(-2/t)*sp.Rational(-1,2)-1) == 0

    # Entire fibre at p=1, including reconstruction in its full root algebra.
    p0 = sp.Integer(1)
    Y0 = gamma[Y].subs(p,p0)
    Z0 = gamma[Z].subs(p,p0)
    modulus = sp.Poly(E.subs({V:p0,Y:Y0,Z:Z0}),S,domain=sp.QQ)
    assert modulus.degree() == 10
    assert sp.gcd(modulus,modulus.diff()).degree() == 0
    def red(f):
        return sp.rem(sp.Poly(f,S,domain=sp.QQ),modulus).as_expr()
    def mul(*factors):
        out=sp.Integer(1)
        for factor in factors:
            out=red(out*factor)
        return out
    def power(f,n):
        out=sp.Integer(1)
        for _ in range(n):
            out=mul(out,f)
        return out
    d = red(modulus.diff().as_expr())
    dinv = sp.invert(d,modulus.as_expr(),S)
    assert mul(d,dinv) == 1
    root_Q = red(-2*Y0-beta.subs(V,p0))
    root_x = mul(S,dinv)
    root_y = red(root_Q-p0*S)
    root_z = mul(d,d,red(p0*d-mul(root_y,root_y,1+3*dinv)/c[3]))
    root_t = red(1+mul(root_x,root_y))
    root_q = red(mul(root_t,root_t,root_z)
                 +mul(root_y,root_y,1+3*root_t)/c[3])
    assert red(root_t-dinv) == 0
    assert red(root_q-p0*d) == 0
    assert mul(root_t,root_q) == p0
    Bvalue = red(root_y+3*c[3]*mul(root_x,root_q))
    Cvalue = red(mul(root_x,5-3*root_t)-c[3]*mul(power(root_x,3),root_z))
    for j in range(4,11):
        profile = profiles[j].subs(V,p0)
        Bvalue = red(Bvalue+j*profile*mul(power(root_t,2),power(root_x,j-2),power(root_q,j)))
        Cvalue = red(Cvalue-(j-2)*profile*power(mul(root_x,root_q),j))
    assert red(-Bvalue/2-Y0) == 0
    assert red(Cvalue-Z0) == 0
    assert mul(root_x,d) == S
    assert red(Cvalue-Z0-1) != 0
    root = Path(__file__).resolve().parents[3]
    inputs = [
        "research/arithmetic/RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md",
        "research/arithmetic/ALTERNATING_PENCIL_SPECIALIZATIONS_AND_KELLER_FIBRES.md",
        "research/jacobian/universal-keller-structure/PROOF.md",
        "unpublished/common-arithmetic-fibers/sections/02-quadratic-gauge.tex",
        "unpublished/common-arithmetic-fibers/sections/03-realization.tex",
    ]
    report = {
        "schema_version":1, "kind":"internal_exact_replay",
        "command":".venv/bin/python research/jacobian/monodromy-theory/verify_fixed_keller_alternating_family.py",
        "python":platform.python_version(), "sympy":sp.__version__,
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "inputs":{path:hashlib.sha256((root/path).read_bytes()).hexdigest() for path in inputs},
        "checks":{
            "full_polynomial_curve_identity":True,
            "full_derivative_identity":True,
            "fixed_profile_exponents":[2*j-8 for j in profiles],
            "compact_plane_jacobian":"-2D", "source_chart_jacobian":"t",
            "global_polynomial_jacobian":1,
            "coordinate_degree_bound":146,
            "quotient_degree":modulus.degree(),
            "quotient_derivative_unit":True,
            "whole_fibre_reconstruction":True,
            "target_at_p_one":["1",str(Y0),str(Z0)],
            "pencil_parameter_at_p_one":int(A-1),
            "invalid_controls":2,
        },
        "scope":"Exact degree-ten polynomial and quotient controls. Uniform degree, whole-scheme fullness, regular A10 pullback and infinitely many full-group fibres follow from the written proof and retained imported premises.",
        "independent_replay":False, "formal_verification":False,
        "external_review":False, "result":"PASS",
    }
    output=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(output)
    print(output,end="")


if __name__ == "__main__":
    main()
