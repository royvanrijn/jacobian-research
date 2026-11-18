#!/usr/bin/env python3
"""Exact controls for the uniform real affine critical-trace obstruction."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
u, r, z, x, s, t = sp.symbols("u r z x s t")


def operator(poly):
    return sp.expand((1-z)**2*poly.subs(z, z+1) -
                     z**2*poly.subs(z, z-1) +
                     2*z*(1-z)*sp.diff(poly, z) -
                     (1-2*z)*poly)


def degree_eight_control():
    N = u**4+u**3+sp.Rational(9, 28)*u**2+sp.Rational(233, 5292)*u-sp.Rational(1, 5488)
    critical = sp.Poly((2*u*sp.diff(N, u)-N)/7, u)
    assert critical.LC() == 1 and sp.discriminant(critical.as_expr(), u) != 0
    mod = lambda expression: sp.rem(expression, critical.as_expr(), u)
    M = mod(N*N*sp.invert(u, critical.as_expr(), u))
    rho = mod(-sp.Rational(268, 2401)-sp.Rational(9261, 2)*M)
    a = mod(N*sp.invert(u-rho, critical.as_expr(), u))
    assert sp.degree(a, u) > 0
    theta = mod(2*rho+M*sp.invert(a*a, critical.as_expr(), u))
    assert mod(theta+sp.Rational(1, 6)+2*rho) == 0
    assert mod(u*u-theta*u+rho*rho) == 0
    assert mod(N-a*(u-rho)) == 0
    characteristic = sp.Poly(sp.resultant(critical.as_expr(), r-rho, u), r).monic()
    q = characteristic.exquo(sp.gcd(characteristic, characteristic.diff())).monic()
    assert q.degree() == 2 and characteristic == q*q and q.nth(0) != 0

    intercept = -sp.Rational(1, 6)
    normalizedN = sp.expand(N.subs(u, intercept*u)/intercept**4)
    normalizedS = sp.expand(critical.as_expr().subs(u, intercept*u)/intercept**4)
    normalizedq = sp.expand(q.as_expr().subs(r, intercept*r)/intercept**2)
    Q = sp.expand(normalizedq.subs(r, z*(1-z)))
    assert sp.degree(Q, z) == 4
    assert sp.expand(normalizedS.subs(u, z*z)-Q*Q.subs(z, -z)) == 0
    reflection = ((1-z)*normalizedN.subs(u, z*z) +
                  z*normalizedN.subs(u, (1-z)**2))
    assert sp.expand(reflection+7*Q*Q) == 0
    assert operator(Q) == 0
    return {
        "N_coefficients_ascending": [str(sp.Poly(N, u).nth(j)) for j in range(5)],
        "canonical_rho_as_polynomial_in_M": "-268/2401-9261*M/2",
        "trace_intercept": "-1/6", "trace_slope": "-2",
        "secant_is_nonconstant": True,
        "scale_is_primitive": True,
        "normalized_q_coefficients_ascending":
            [str(sp.Poly(normalizedq, r).nth(j)) for j in range(3)],
        "critical_factorization_and_reflection_identity": True,
        "operator_vanishes_in_degree4": True}


def run():
    theta = s+t*r
    delta = theta**2-4*r*r
    assert sp.expand((theta+2*r)*sp.diff(delta, r)-2*t*delta -
                     4*r*(theta+2*r)*(t-2)) == 0

    gamma = [sp.Integer(0), sp.Integer(1)]
    first_kind = [sp.Integer(2), x]
    for j in range(2, 33):
        gamma.append(sp.expand(t*gamma[-1]-gamma[-2]))
        first_kind.append(sp.expand(x*first_kind[-1]-first_kind[-2]))
    controls = []
    for h in range(1, 17):
        m = 2*h
        L = 2*m-1
        numerator = sp.Poly(1+sum(first_kind[1:m])-L, x, domain=sp.ZZ)
        B, rem = numerator.div(sp.Poly(x-2, x, domain=sp.ZZ))
        assert rem.is_zero
        assert sp.expand(gamma[m]-gamma[m-1]+L -
                         (t+2)*B.as_expr().subs(x, -t)) == 0
        assert B.eval(2) == sum(j*j for j in range(1, m))
        value = operator(z**m)
        polynomial = sp.Poly(value, z)
        assert polynomial.is_zero or polynomial.degree() <= m-1
        coefficient = sp.Rational(m*(m-2)*(m-4), 3)
        assert polynomial.nth(m-1) == coefficient
        controls.append({"h": h, "pencil_degree": 4*h,
                         "slope_witness_identity": True,
                         "operator_leading_coefficient": str(coefficient)})

    # Direct leading/constant trace-basis checks independent of the
    # first-kind witness construction.
    C = [sp.Integer(0), sp.Integer(1)]
    T = [sp.Integer(2), theta]
    for j in range(2, 11):
        C.append(sp.expand(theta*C[-1]-r*r*C[-2]))
        T.append(sp.expand(theta*T[-1]-r*r*T[-2]))
    for j in range(1, 11):
        D = sp.Poly(T[j]-(theta-2*r)*C[j], r)
        assert D.degree() <= j and D.nth(0) == 0
        assert sp.expand(D.nth(j)-2*(gamma[j]-gamma[j-1])) == 0

    q0, q1, q2 = sp.symbols("q0 q1 q2")
    R = z*(1-z)
    Q12 = R**3+q2*R**2+q1*R+q0
    assert sp.expand(operator(Q12)-8*(1-2*z)*R*R) == 0
    n = sp.symbols("n")
    assert sp.expand(n-n*(n-1)+n*(n-1)*(n-2)/3 -
                     n*(n-2)*(n-4)/3) == 0

    degree8 = degree_eight_control()
    paths = [HERE/"UNIFORM_REAL_AFFINE_TRACE_OBSTRUCTION.md",
             Path(__file__).resolve(),
             HERE/"RATIONAL_ALTERNATING_PENCIL_REALIZATION.md"]
    return {
        "schema_version": 1, "result": "PASS",
        "python": platform.python_version(), "sympy": sp.__version__,
        "uniform_claim": "No good real even-numerator paired pencil ofdegree4h>=12 has trace affine in canonical scale.",
        "uniform_corollary": "No such real pencil has scalar secant and affine scale.",
        "uniform_written_proof_required": True,
        "finite_control_scope": "h1,...,16; trace-basis degrees1,...,10; exact degree8 boundary control.",
        "critical_jet_identity": True,
        "trace_constant_and_leading_basis_identities": True,
        "controls": controls,
        "degree12_operator": "8*(1-2*z)*z^2*(1-z)^2",
        "degree8_boundary_control": degree8,
        "scope_limits": ["Even numerator degree only.", "Trace affine in scale, not every singular spectral cubic.",
                         "No restriction that critical points be real.",
                         "Nonaffine trace/nonconstant secant general models remain open."],
        "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in paths},
        "assurance": {"formal_verification": False,
                      "independent_external_replay": False,
                      "external_human_review": False}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    result = run()
    payload = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()

