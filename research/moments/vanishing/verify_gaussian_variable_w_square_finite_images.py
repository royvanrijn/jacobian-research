"""Bounded exact controls for the unrestricted variable-W-square theorem.

The written original-source, all-face and complete period-comparison proof
establishes the unrestricted result. These are finite algebraic controls.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp

z, w, y, t, r, tau = sp.symbols("z w y t r tau", nonzero=True)


def zero(expression):
    assert sp.cancel(sp.expand(expression)) == 0, expression


def wick(expression):
    answer = 0
    for (iz, iw, iy), coefficient in sp.Poly(sp.expand(expression), z, w, y).terms():
        if iz == iw and iy % 2 == 0:
            answer += coefficient*factorial(iz)*sp.Rational(
                factorial(iy), 2**(iy//2)*factorial(iy//2))
    return sp.expand(answer)


def periods(phase, markers, order):
    result = [[] for _ in markers]
    power = sp.Integer(1)
    for n in range(order+1):
        if n:
            power = sp.expand(power*phase)
        for row, marker in zip(result, markers):
            row.append(sp.expand(wick(marker*power)/factorial(n)))
    return result


def exponential_coefficients(exponent, order):
    result = [sp.Integer(1)]
    for n in range(1, order+1):
        result.append(sp.expand(sum(j*exponent.coeff(t, j)*result[n-j]
                                    for j in range(1, n+1))/n))
    return result


def checks():
    alpha, kappa, aa, h0, h1 = sp.symbols("alpha kappa aa h0 h1", nonzero=True)
    face_cases = []
    leading_controls = []
    tail_b, tail_T = sp.symbols("b T", positive=True)
    for m in [3, 5, 7, 9]:
        sigma_power = (m-3)//2
        U = sp.sqrt(2*kappa*(r-2*aa))*w**m
        H = 2*r/sp.sqrt(2*kappa*(r-2*aa))
        zero(r*2*w*w**(2*sigma_power)/U-H/w**2)
        hh = h0+h1*(z-alpha)
        phi = (-((r*z-hh)**2)/(4*kappa*r*(z-alpha)**m)
               +z**2/r+(1+z)**2/(2*r*(r-2*aa)))
        pulled = sp.expand(phi.subs(z, alpha+w**2))
        leading_pole = pulled.coeff(w, -2*m)
        zero(leading_pole+(r*alpha-h0)**2/(4*kappa*r))

        face = -(alpha-h1*tau**2)**2/(4*kappa*tau**(2*m))
        derivative = ((alpha-h1*tau**2)*(m*alpha-(m-2)*h1*tau**2)
                      /(2*kappa*tau**(2*m+1)))
        zero(sp.diff(face, tau)-derivative)
        # The two possible critical levels are distinct; every critical
        # point is simple as a derivative zero, hence shifted roots <=2.
        x_other = m*alpha/((m-2)*h1)
        zero((alpha-h1*x_other)+2*alpha/(m-2))
        assert sp.cancel(alpha-h1*x_other) != 0

        inverse_w = (tail_b/tail_T)**sp.Rational(1, 2*m)
        expected_tail = (-tail_b**sp.Rational(-1, 2*m)
                         *tail_T**(sp.Rational(1, 2*m)-1)/(2*m))
        assert sp.simplify(sp.diff(inverse_w, tail_T)/inverse_w**2-expected_tail) == 0
        assert sp.Rational(1, 2)-sp.Rational(1, 2*m) > 0
        assert sp.Rational(m-1, 2) > 0
        face_cases.append(dict(m=m, sigma_power=sigma_power,
                               critical_displacement_bound=str(sp.Rational(m-1, 2)),
                               worst_deep_tail_order=str(sp.Rational(1, 2)-sp.Rational(1, 2*m))))

        outside = z**m*w**2+(2+z)*w+y**2/3+(1+z**2)*y+z**3
        leading = periods(outside, [z**sigma_power], sigma_power)[0]
        assert leading[:-1] == [0]*sigma_power
        assert leading[-1] == 2**sigma_power
        leading_controls.append(dict(m=m, original_marker=str(z**sigma_power),
                                     leading_coefficient=int(leading[-1])))

    face3 = -(1-tau**2)**2/(4*tau**6)
    zero(face3+sp.Rational(1, 27)
         -(tau**2-3)**2*(4*tau**2-3)/(108*tau**6))
    zero((face3+sp.Rational(1, 27)).subs(tau, sp.sqrt(3)))
    zero(sp.diff(face3, tau).subs(tau, sp.sqrt(3)))
    assert sp.simplify(sp.diff(face3, tau, 2).subs(tau, sp.sqrt(3))) != 0

    reductions = []
    for ell in [z**3, (z-1)**9, (z-1)**8,
                z**3*(z-1)**5, z**4*(z-1)**3,
                z**2*(z-1)**3*(z+2)**5]:
        factors = sp.sqf_list(ell, z)[1]
        rho = sp.prod(factor**((int(m)-1)//2) for factor, m in factors)
        reduced = sp.cancel(ell/rho**2)
        residual_degree = int(sp.degree(reduced, z))
        roots_weight = sum(int(sp.degree(factor, z))*(1 if int(m)%2 else 2)
                           for factor, m in factors)
        assert residual_degree == roots_weight
        single_odd = (len(factors) == 1 and sp.degree(factors[0][0], z) == 1
                      and int(factors[0][1])%2 == 1)
        assert (residual_degree == 1) == single_odd
        reductions.append(dict(ell=str(sp.factor(ell)), rho=str(sp.factor(rho)),
                               reduced_degree=residual_degree, odd_single_exception=bool(single_odd)))

    plus = z**3*w**2+z**2*w+z*y+z
    phase = plus+2*y+3
    markers = [sp.Integer(1), z, w, y, w**2, z*w, w*y]
    order = 6
    rows = periods(phase, markers, order)
    exp_psi = exponential_coefficients(3*t+2*t**2, order)
    normalized = []
    for marker, row in zip(markers, rows):
        degree = sp.Poly(marker, z, w, y).total_degree()
        polynomial = sp.expand(sum(t**n*wick(
            sp.expand((marker*plus**n).subs(y, y+2*t)))/factorial(n)
            for n in range(degree+1)))
        expected = [sp.expand(sum(exp_psi[n-j]*polynomial.coeff(t, j)
                                  for j in range(n+1))) for n in range(order+1)]
        assert row == expected
        assert polynomial == 0 or sp.degree(polynomial, t) <= 2*degree
        if marker == w:
            zero(polynomial-(9*t+2*t**2))
        if marker == w**2:
            zero(polynomial-(524*t**2+64*t**3+4*t**4))
        normalized.append(str(polynomial))
    safe_rows = periods(plus, markers, order)
    for marker, row in zip(markers, safe_rows):
        degree = sp.Poly(marker, z, w, y).total_degree()
        assert all(value == 0 for value in row[degree+1:])

    degree_one = periods(z*w**2, [sp.Integer(1)], 6)[0]
    assert degree_one == [1, 0, 0, 0, 0, 0, 0]
    double_origin = periods(z**2*w**2, [sp.Integer(1)], 4)[0]
    assert double_origin == [sp.Rational(factorial(2*n), factorial(n)) for n in range(5)]

    return dict(status="passed", exact_face_and_tail_cases=face_cases,
                critical_cubic_shift="F+1/27=(tau^2-3)^2(4tau^2-3)/(108tau^6); multiplicity2",
                original_zero_mass_marker_controls=leading_controls,
                fixed_cover_reductions=reductions,
                finite_branch_original_markers=[str(v) for v in markers],
                finite_branch_normalized_periods=normalized,
                original_Wick_orders=[0, order],
                safe_marker_cutoff="m>deg Q through moment6 for seven markers",
                scope_controls=["degree-one ell=z is pure/safe and outside degell>=2",
                                "common double origin has factorial original unit; do not transport it through a clean-source argument"],
                scope="Finite exact controls only. All scalar shifts, all scales and full marked-curve gluing are proved in the written argument, including the permanently irregular double pole.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    names = ["GAUSSIAN_VARIABLE_W_SQUARE_FINITE_IMAGES.md", Path(__file__).name,
             "GAUSSIAN_QUADRATIC_W_MIXED_FINITE_IMAGES.md",
             "GAUSSIAN_AFFINE_ISOTROPIC_CHART.md",
             "GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md",
             "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
             "../../../results/gaussian-boundary/analytic_core.tex",
             "../../../results/gaussian-boundary/seed_period_comparison.tex"]
    result.update(python=platform.python_version(), sympy=sp.__version__,
                  sha256={name: hashlib.sha256((directory/name).read_bytes()).hexdigest()
                          for name in names},
                  independent_replay=False, formal_verification=False, external_review=False)
    report = json.dumps(result, indent=2)+"\n"
    if args.output is not None:
        args.output.write_text(report, encoding="utf-8", newline="\n")
    print(report, end="")


if __name__ == "__main__":
    main()
