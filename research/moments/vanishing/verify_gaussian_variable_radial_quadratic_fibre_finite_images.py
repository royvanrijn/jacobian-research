"""Bounded exact controls for variable radial quadratic Gaussian fibres.

The written proof establishes the all-degree relative-source and full
period-comparison theorem. These finite controls do not extrapolate a prefix.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp

s, y, z, w, t, r, V = sp.symbols("s y z w t r V")


def zero(expr):
    assert sp.simplify(sp.powsimp(sp.expand(expr), force=True)) == 0, expr


def gamma_expect(poly):
    return sp.expand(sum(c*factorial(i) for (i,), c in sp.Poly(sp.expand(poly), s).terms()))


def radial_wick(poly):
    answer = 0
    for (i, j), coefficient in sp.Poly(sp.expand(poly), s, y).terms():
        if j % 2 == 0:
            answer += coefficient*factorial(i)*sp.Rational(
                factorial(j), 2**(j//2)*factorial(j//2))
    return sp.expand(answer)


def wick(poly):
    answer = 0
    for (i, j, k), coefficient in sp.Poly(sp.expand(poly), z, w, y).terms():
        if i == j and k % 2 == 0:
            answer += coefficient*factorial(i)*sp.Rational(
                factorial(k), 2**(k//2)*factorial(k//2))
    return sp.expand(answer)


def exp_coefficients(u_coefficients, order):
    result = [sp.Integer(1)]
    for n in range(1, order+1):
        result.append(sp.expand(sum(j*u_coefficients[j]*result[n-j]
                                    for j in range(1, n+1))/n))
    return result


def completed_coefficients(a, B, g, order):
    u = [sp.Integer(0)]*(order+1)
    if order:
        u[1] = g
    for n in range(2, order+1):
        u[n] = sp.expand(B**2*(2*a)**(n-2)/2)
    exponent = exp_coefficients(u, order)
    determinant = [sp.binomial(2*j, j)*a**j/2**j for j in range(order+1)]
    return [sp.expand(sum(determinant[j]*exponent[n-j] for j in range(n+1)))
            for n in range(order+1)]


def y_moment_coefficients(k, a, B, order):
    result = [sp.Integer(0)]*(order+1)
    for ell in range(k//2+1):
        start = k-2*ell
        prefactor = sp.Rational(factorial(k),
                               2**ell*factorial(ell)*factorial(k-2*ell))*B**start
        for extra in range(order-start+1):
            result[start+extra] += prefactor*sp.rf(k-ell, extra)*(2*a)**extra/factorial(extra)
    return [sp.expand(value) for value in result]


def source_controls():
    order = 4
    cases = [
        ("quadratic-root", s**2, 1+s, s**2),
        ("triple-root", (s-1)**3, s, 1+s**2),
        ("affine-cancelling-top", 2*s+1, s**2+s+1, s**3/8+3*s**2/16+s),
    ]
    markers = [(0,0), (1,0), (2,0), (0,1), (0,2), (1,1), (0,3)]
    records = []
    for name, a, B, g in cases:
        phase = g+B*y+a*y**2
        q = completed_coefficients(a, B, g, order)
        powers = [sp.expand(phase**n) for n in range(order+1)]
        for i, k in markers:
            yy = y_moment_coefficients(k, a, B, order)
            for n in range(order+1):
                original = radial_wick(s**i*y**k*powers[n])/factorial(n)
                completed = gamma_expect(s**i*sum(q[n-j]*yy[j] for j in range(n+1)))
                zero(original-completed)

        # Relative boundary condition is checked with rational d^-1 markers.
        for relative, power in [(True,0), (True,1), (False,2)]:
            f = [sp.expand((s if relative else 1)*sp.rf(power,n)*(2*a)**n/factorial(n))
                 for n in range(order+1)]
            for n in range(order+1):
                full = sp.expand(sum(f[j]*q[n-j] for j in range(n+1)))
                lhs = gamma_expect(sp.diff(full,s)-full)
                zero(lhs+full.subs(s,0))
                if relative:
                    zero(lhs)
                elif n == 0:
                    assert lhs == -1
        records.append(dict(name=name, a=str(a), B=str(B), g=str(g),
                            original_balanced_markers=[f"S^{i}Y^{k}" for i,k in markers],
                            Wick_orders=[0,order], relative_endpoint_retained=True))
    return records


def full_source_operator_controls():
    aa = sp.Function("a")(s)
    bb = sp.Function("B")(s)
    gg = sp.Function("g")(s)
    d = 1-2*t*aa
    U = t*gg+t**2*bb**2/(2*d)
    vlogs = sp.diff(d,s)/(2*d)
    vlogt = -1/(2*t*d)
    clogt = -1/(2*t)
    zero(sp.diff(vlogs,t)-sp.diff(vlogt,s))
    zero(clogt-vlogt-aa/d)
    K = -1+sp.diff(U,s)-vlogs
    J = sp.diff(U,t)-vlogt
    zero(sp.diff(K,t)-sp.diff(J,s))
    eta_J = sp.diff(U,t)+aa/d
    zero(sp.diff(K,t)-sp.diff(eta_J,s))

    # Exact scalar gauge product rule, without asserting a false generic
    # scalar identity for a nonfinite phase.
    p, psi = 1+t, t+t**2
    zero(sp.diff(sp.exp(-psi)/p,t)/(sp.exp(-psi)/p)
         +sp.diff(psi,t)+sp.diff(p,t)/p)
    return dict(fixed_S_V_log_derivative="-1/(2td)",
                eta_log_derivative="a/d",
                spatial_time_flatness="exact with arbitrary symbolic a(s),B(s),g(s)",
                gauge_product_rule=True)


def endpoint_controls():
    a1 = sp.symbols("a1", nonzero=True)
    a0 = sp.symbols("a0")
    beta = (1/(2*t)-a0)/a1
    Shat = beta+V**2/a1
    lift = lambda value: sp.diff(value,t)+sp.diff(value,V)/(4*t**2*V)
    zero(lift(Shat))
    zero(lift(Shat/V)+Shat/(4*t**2*V**3))
    naive = sp.diff(Shat/V,t)
    zero(naive+1/(2*a1*t**2*V))
    assert naive != 0
    zero((-1)/(-V)-1/V)  # Logarithmic residue is deck invariant.
    zero(Shat.subs(V,-V)-Shat)
    zero((Shat/V).subs(V,-V)+Shat/V)
    return dict(moving_endpoint="S=0, V^2=a0-r/2",
                exact_lift="partial_t|V+r^2/(4V)partial_V",
                rejected_naive_fixed_V_derivative=str(naive),
                anti_relative_endpoint_line="retained",
                regular_branch_Gysin="dV/V is deck invariant")


def volume_and_affine_controls():
    m = sp.symbols("m", positive=True, integer=True)
    theta, eps, N = sp.symbols("theta eps N", positive=True)
    ell, k = sp.symbols("ell k", nonnegative=True, integer=True)
    subcritical = sp.Rational(1,2)+(1-m/2)*theta/m
    zero(subcritical-((1-theta)/2+theta/m))
    zero(sp.Rational(1,2)+(sp.Rational(1,2)+eps)-(m-1)/m-(1/m+eps))
    zero(sp.Rational(1,2)+((2+ell)/2-1)*N-(sp.Rational(1,2)+ell*N/2))

    a1, bm = sp.symbols("a1 bm", nonzero=True)
    a0, bprev, gnext = sp.symbols("a0 bprev gnext")
    beta = (r/2-a0)/a1
    topg = bm**2/(4*a1)
    next_coefficient = (((2*m-1)*topg*beta+gnext)/(r*a1**(2*m-2))
                        -(2*m*bm**2*beta+2*bm*bprev)/(4*r*a1**(2*m-1)))
    C = gnext/a1**(2*m-2)-bm*bprev/(2*a1**(2*m-1))+bm**2*a0/(4*a1**(2*m))
    zero(next_coefficient-(C/r-bm**2/(8*a1**(2*m))))
    zero(sp.diff(r*next_coefficient,r)+bm**2/(8*a1**(2*m)))

    # Representative top cancellation has actual degree4 with valuation0.
    B = s**2+s+1
    g = s**3/8+3*s**2/16+s
    beta0 = (r/2-1)/2
    ss = beta0+V**2/2
    phase = sp.expand((-s+g/r-B**2/(4*r*V**2)).subs(s,ss))
    zero(phase.coeff(V,6))
    zero(phase.coeff(V,4)+sp.Rational(1,128))
    zero(phase.coeff(V,-2)+B.subs(s,beta0)**2/(4*r))

    # Generic M>1/2 gap and volume identities (all D=2+2k, no integer prefix).
    M = sp.Rational(1,2)+theta
    DD = 2+2*k
    zero(2*M-1-2*theta)
    for nu in [-1,0]:
        action = DD*M-nu
        zero((sp.Rational(1,2)-M+action)
             -(sp.Rational(1,2)+(DD-1)*M-nu))

    # Descending T=(-C_D)V^D; its valuation is independent of the unit sign.
    T, Cdesc = sp.symbols("T Cdesc", positive=True)
    D = sp.symbols("D", positive=True, integer=True)
    inverse = (T/Cdesc)**(1/D)
    zero(sp.diff(inverse,T)-Cdesc**(-1/D)*T**(1/D-1)/D)
    return dict(all_multiplicities="convex bound min(1/2,1/m); cancellation >1/m",
                higher_degree_infinity="1/2+(L-2)N/2 >=1/2 for ALL L>=2",
                generic_affine_top_cancellation="C/r-b_m^2/(8a1^(2m))",
                affine_actual_top_control="degree4 coefficient -1/128, after degree6 cancellation",
                affine_all_M_cases=["M<1/2", "M=1/2 finite-multiplicity epsilon>0",
                                    "M>1/2 top gap>=2M-1"],
                descending_infinity_tail="T=-C_D V^D, scalar order1/2-v(C_D)/D>=1/2")


def finite_marker_controls():
    plus = z**3*w**2+z*y**2+z**2*w*y
    phase = plus+2*y+3
    markers = [sp.Integer(1), z, w, y, w**2, z*w, w*y]
    order = 7
    epsi = exp_coefficients([0,3,2]+[0]*(order-2),order)
    normalized = []
    for marker in markers:
        degree = sp.Poly(marker,z,w,y).total_degree()
        translated = sp.expand(sum(t**n*wick(
            (marker*plus**n).subs(y,y+2*t))/factorial(n) for n in range(degree+1)))
        for n in range(order+1):
            original = wick(marker*sp.expand(phase**n))/factorial(n)
            predicted = sum(epsi[n-j]*translated.coeff(t,j) for j in range(n+1))
            zero(original-predicted)
            if n > degree:
                zero(wick(marker*sp.expand(plus**n)))
        normalized.append(str(translated))
    return dict(original_markers=[str(q) for q in markers],
                original_Wick_orders=[0,order], normalized_periods=normalized,
                safe_C1_controls="moments n>deg Q vanish within the bounded window")


def quartic_support_controls():
    data = {}
    for degree in [4,5]:
        monomials = {(i,j) for i in range(degree//2+1) for j in range(degree-2*i+1)}
        remaining = {(i,j) for i,j in monomials if i > 0 or j <= 1}
        expected = ({(0,0),(0,1),(1,0),(1,1),(1,2),(2,0)} if degree == 4 else
                    {(0,0),(0,1),(1,0),(1,1),(1,2),(1,3),(2,0),(2,1)})
        assert remaining == expected
        data[str(degree)] = [f"S^{i}Y^{j}" for i,j in sorted(remaining)]
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = dict(
        status="passed",
        actual_radial_source_controls=source_controls(),
        full_source_covariant_controls=full_source_operator_controls(),
        relative_endpoint_controls=endpoint_controls(),
        symbolic_volume_and_top_controls=volume_and_affine_controls(),
        finite_original_marker_controls=finite_marker_controls(),
        boundary_structural_support_controls=quartic_support_controls(),
        scope=("Bounded exact controls only. All multiplicities and infinity scales "
               "are covered by the written proof, not finite prefixes. The quartic "
               "corollary uses full image finiteness; no pure-quartic necessity is inferred."),
        independent_replay=False, formal_verification=False, external_review=False)
    directory = Path(__file__).resolve().parent
    names = [
        "GAUSSIAN_VARIABLE_RADIAL_QUADRATIC_FIBRE_FINITE_IMAGES.md",
        Path(__file__).name,
        "GAUSSIAN_RADIAL_ONE_SIDED_FINITE_IMAGES.md",
        "GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md",
        "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
        "GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md",
        "GAUSSIAN_AFFINE_ISOTROPIC_CHART.md",
        "GAUSSIAN_VARIABLE_W_SQUARE_FINITE_IMAGES.md",
        "GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md",
        "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md",
        "GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md",
        "../../../unpublished/gaussian-boundary/analytic_core.tex",
        "../../../unpublished/gaussian-boundary/seed_period_comparison.tex",
        "../../../unpublished/gaussian-boundary/sources/laurent_primary_manifest.json",
    ]
    report.update(python=platform.python_version(), sympy=sp.__version__,
                  sha256={name:hashlib.sha256((directory/name).read_bytes()).hexdigest()
                          for name in names})
    text=json.dumps(report,indent=2)+"\n"
    if args.output is not None:
        args.output.write_text(text,encoding="utf-8",newline="\n")
    print(text,end="")


if __name__ == "__main__":
    main()
