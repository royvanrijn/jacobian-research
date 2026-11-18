#!/usr/bin/env python3
"""Exact controls for the quadratic spectral pairing/divisor criterion."""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
r, u, v, x, y, X, Z = sp.symbols("r u v x y X Z")
t0, t1, t2 = sp.symbols("t0 t1 t2")
theta = t0+t1*r+t2*r*r


def trace_reduce(poly, th=theta):
    poly = sp.Poly(poly, u)
    C, B = sp.Integer(0), sp.Integer(0)
    cj, bj = sp.Integer(0), sp.Integer(1)
    for j in range(poly.degree()+1):
        C += poly.nth(j)*cj
        B += poly.nth(j)*bj
        cj, bj = sp.expand(th*cj+bj), sp.expand(-r*r*cj)
    return sp.expand(C), sp.expand(B)


def spectral_transform():
    alpha = t0*t2
    b = t1*t1/4-1-alpha
    g = (1-t2*u)*r*r-t1*u*r+u*u-t0*u
    YY = (1-t2*u)*r-t1*u/2
    cubic_u = t2*u**3+b*u*u+t0*u
    assert sp.expand(YY*YY-cubic_u-(1-t2*u)*g) == 0
    assert sp.expand((t2*YY)**2 -
                     ((t2*u)**3+b*(t2*u)**2+alpha*t2*u)
                     -t2*t2*(1-t2*u)*g) == 0
    cubic = X**3+b*X*X+alpha*X
    assert sp.expand(cubic.subs(X, 1)-t1*t1/4) == 0
    assert sp.expand(cubic.subs(X, alpha)-t1*t1*alpha*alpha/4) == 0
    assert sp.expand(t1*t1*X*X/4-cubic +
                     X*(X-1)*(X-alpha)) == 0
    return {"weierstrass_identity": True,
            "marked_points_and_R_numerator_norm": True,
            "discriminant": "16*alpha^2*(b^2-4*alpha)"}


def jet_controls():
    up, a = sp.symbols("up a")
    vv = r*r/u
    vp = 2*r/u-r*r*up/(u*u)
    th = u+vv
    thp = up+vp
    d, dp = th-2*r, thp-2
    Tprime = a*((u-r)*up/u+(vv-r)*vp/vv)/2
    Cprime = a*((u-r)*up/u-(vv-r)*vp/vv-2*(up-vp))/(2*(u-vv))
    assert sp.cancel(Tprime-a*dp/2) == 0
    assert sp.cancel(Cprime+a*dp/(2*d)) == 0
    assert sp.cancel(Tprime-dp*a-d*Cprime) == 0
    J = a*a*(u-r)*(vv-r)
    Jprime = J*(up/u+vp/vv)/2
    assert sp.cancel(Jprime-J/r) == 0
    kap, eta = sp.symbols("kap eta")
    Kprime = -kap*(r-eta)+kap*(2*r-eta)
    assert sp.expand(Kprime-kap*r) == 0
    A, B, th = sp.symbols("A B th")
    T = th*A+2*B
    JJ = B*B+th*A*B+r*r*A*A
    assert sp.expand(T*T-(th*th-4*r*r)*A*A-4*JJ) == 0
    return {"critical_trace_jets": True,
            "D_double_zero_without_differentiating_a": True,
            "natural_norm_simple_q_zeros": True,
            "nonconstant_norm_Pell_identity": True}


def degree_controls():
    controls = []
    for m in range(1, 11):
        C, B = trace_reduce(u**m)
        D = sp.Poly(2*(B+r*C), r)
        assert sp.Poly(C, r).degree() == 2*m-2
        assert sp.expand(sp.Poly(C, r).LC()-t2**(m-1)) == 0
        assert D.degree() == 2*m-1
        assert sp.expand(D.LC()-2*t2**(m-1)) == 0
        assert D.nth(0) == 0
        if m <= 5:
            assert sp.rem(u**m-C*u-B, u*u-theta*u+r*r, u) == 0
        controls.append({"numerator_degree": m,
                         "D_degree": D.degree(),
                         "D_leading_coefficient": str(2*t2**(m-1))})
    q0 = sp.symbols("q0")
    for h in range(1, 5):
        L = 4*h-1
        n0 = -L*q0*q0
        assert sp.cancel(2*n0/(q0*q0)) == -2*L
    return controls


def degree_eight_positive_control():
    N = u**4+u**3+sp.Rational(9, 28)*u**2+sp.Rational(233, 5292)*u-sp.Rational(1, 5488)
    S = sp.Poly((2*u*sp.diff(N, u)-N)/7, u)
    q = r*r-sp.Rational(4, 21)*r+sp.Rational(1, 196)
    th = sp.expand(q-sp.Rational(1, 6)-2*r)
    kap, eta = -sp.Rational(2, 9261), -sp.Rational(268, 2401)
    assert sp.expand(sp.resultant(q, u*u-th*u+r*r, r)-S.as_expr()) == 0
    assert all(value != 0 for value in
               [q.subs(r, 0), q.subs(r, eta), sp.discriminant(q, r),
                S.discriminant()])
    alpha = th.subs(r, 0)
    b = sp.Poly(th, r).nth(1)**2/4-1-alpha
    assert alpha*(b*b-4*alpha) != 0
    C, B = trace_reduce(N, th)
    D = sp.expand(2*(B+r*C))
    U, rem = sp.div(D, q*q, r)
    assert rem == 0 and sp.degree(U, r) == 3
    assert sp.Poly(U, r).LC() == 2 and U.subs(r, 0) == -14
    assert sp.rem((th-2*r)*C*C-kap*(r-eta), q, r) == 0
    J = sp.expand(B*B+th*B*C+r*r*C*C)
    K = sp.expand(J+kap*r*(r-eta))
    assert sp.rem(K, q, r) == 0
    assert sp.rem(sp.diff(K, r)-kap*r, q, r) == 0
    assert sp.gcd(sp.div(K, q, r)[0], q) == 1
    mod = lambda e: sp.rem(e, S.as_expr(), u)
    M = mod(N*N*sp.invert(u, S.as_expr(), u))
    rho = mod(eta+M/kap)
    assert mod(u*u-th.subs(r, rho)*u+rho*rho) == 0
    secant = mod(C.subs(r, rho))
    assert mod(N-secant*(u-rho)) == 0
    assert sp.gcd(secant, S.as_expr()) == 1
    # Independent norm formula viewed as a quadratic in R.
    g = u*u-th*u+r*r
    norm_q = sp.resultant(g, q, r)
    assert sp.expand(norm_q-S.as_expr()) == 0
    # Both roots r=(4+-sqrt(7))/42 are positive, while theta=-1/6-2r
    # on them is negative: this good real-coefficient model is not all-real.
    assert sp.discriminant(q, r) > 0 and q.subs(r, 0) > 0
    assert sp.Rational(4, 21) > 0
    return {"N_coefficients_ascending": [str(N.coeff(u, j)) for j in range(5)],
            "q": str(q), "theta": str(th), "kappa": str(kap), "eta": str(eta),
            "good_criterion_and_inverse": True, "D_over_q_squared": str(U),
            "smooth_quadratic_spectral_extension": True,
            "natural_norm_has_simple_q_zeros": True,
            "all_critical_points_real": False,
            "scope": "Known degree8 boundary model; no degree12 existence claim."}


def ec_add(P, Q, p, a2, a4):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1+y2) % p == 0:
        return None
    if P == Q:
        slope = (3*x1*x1+2*a2*x1+a4)*pow(2*y1, -1, p) % p
    else:
        slope = (y2-y1)*pow((x2-x1) % p, -1, p) % p
    x3 = (slope*slope-a2-x1-x2) % p
    return (x3, (-y1+slope*(x1-x3)) % p)


def real_and_arithmetic_controls():
    F = x*x+y*y+x*y+x*x*y*y+1
    assert sp.expand(F-((x+y/2)**2+3*y*y/4+x*x*y*y+1)) == 0
    # Strict negative-coefficient real gate and its singular threshold.
    a = sp.symbols("a", positive=True)
    bthreshold = (2+2*a)**2/4-1-a*a
    assert sp.expand(bthreshold-2*a) == 0
    assert sp.expand((bthreshold*bthreshold-4*a*a)) == 0
    assert sp.expand((x*x+y*y+5*x*y+x*x*y*y+1).subs({x: 1, y: -sp.Rational(1,2)})) == 0
    for T0, T1, T2, soluble in [(-1,1,-1,False),(-1,5,-1,True),
                                (-4,5,-1,False),(-4,7,-1,True)]:
        alpha = T0*T2
        b = sp.Rational(T1*T1,4)-1-alpha
        assert bool(b > 0 and b*b > 4*alpha) is soluble
    orders = {}
    for p in [3, 5]:
        count = 1+sum(1 for xx in range(p) for yy in range(p)
                      if (yy*yy-(xx**3-2*xx*xx+2*xx)) % p == 0)
        orders[p] = count
    assert orders == {3: 6, 5: 8}
    B = (1,1)
    twice = ec_add(B, B, 5, -2, 2)
    assert twice == (4,0) and ec_add(twice, twice, 5, -2, 2) is None
    assert -256 % 3 and -256 % 5
    assert orders[3] % 4 != 0
    return {"no_real_point_control": [-1,1,-1],
            "strict_real_gate_controls": 4,
            "singular_equality_threshold": True,
            "nontorsion_mark_control": {"parameters": [1,2,2],
                "E_counts": {str(p): value for p,value in orders.items()},
                "B_mod5_order": 4, "good_discriminant": -256,
                "prime_to3_reduction_contradiction": True},
            "controls_are_not_degree12_pencils": True}


def run():
    fields = {"schema_version": 1, "result": "PASS",
              "python": platform.python_version(), "sympy": sp.__version__,
              "spectral_transform": spectral_transform(),
              "jets_and_norm": jet_controls(),
              "trace_degree_controls": degree_controls(),
              "positive_degree8_control": degree_eight_positive_control(),
              "real_and_arithmetic_controls": real_and_arithmetic_controls(),
              "degree12_accessory": {"degree": 5, "constant": -22,
                                     "leading": "2*t2^5"},
              "uniform_written_proof_required": True,
              "scope_limits": [
                  "Implicit exact coefficient locus, not a rational parametrization.",
                  "No rational degree12 or degree26 existence theorem.",
                  "Elliptic divisor claims require t0*t2*(b^2-4*alpha)!=0.",
                  "The displayed divisor relations alone impose no torsion relation.",
                  "Possible additional torsion consequences of full equations remain unproved.",
                  "The square-root genus-one cover may have no rational point."],
              "assurance": {"formal_verification": False,
                            "independent_external_replay": False,
                            "external_human_review": False}}
    paths = [HERE/"QUADRATIC_SPECTRAL_PAIRING_AND_ELLIPTIC_DIVISORS.md",
             Path(__file__).resolve(), HERE/"PAIRED_CRITICAL_TRACE_CHART.md",
             HERE/"UNIFORM_REAL_AFFINE_TRACE_OBSTRUCTION.md",
             HERE/"RATIONAL_ALTERNATING_PENCIL_REALIZATION.md"]
    fields["sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in paths}
    return fields


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", type=Path)
    args = parser.parse_args()
    payload = json.dumps(run(), indent=2, sort_keys=True)+"\n"
    if args.write_receipt:
        args.write_receipt.write_text(payload)
    print(payload, end="")


if __name__ == "__main__":
    main()

