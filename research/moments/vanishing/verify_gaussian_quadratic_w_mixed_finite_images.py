"""Exact bounded controls for the quadratic-W mixed finite-image theorem.

The written proof, not this bounded replay, supplies every root multiplicity,
every original marker and the complete period-comparison argument.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp

z, w, y, t, r = sp.symbols("z w y t r")
ss = sp.symbols("s")


def zero(value):
    assert sp.cancel(sp.expand(value)) == 0, value


def wick(value):
    answer = 0
    for (iz, iw, iy), scalar in sp.Poly(sp.expand(value), z, w, y).terms():
        if iz == iw and iy % 2 == 0:
            answer += scalar*factorial(iz)*sp.Rational(
                factorial(iy), 2**(iy//2)*factorial(iy//2))
    return sp.expand(answer)


def periods(phase, markers, order):
    answer = [[] for _ in markers]
    power = sp.Integer(1)
    for n in range(order+1):
        if n:
            power = sp.expand(power*phase)
        for row, marker in zip(answer, markers):
            row.append(sp.expand(wick(marker*power)/factorial(n)))
    return answer


def exponential_coefficients(exponent, order):
    result = [sp.Integer(1)]
    for n in range(1, order+1):
        result.append(sp.expand(sum(j*exponent.coeff(t, j)*result[n-j]
                                    for j in range(1, n+1))/n))
    return result


def radical_marker(ell, c):
    g = sp.gcd(sp.Poly(ell, z, domain=sp.QQ),
               sp.Poly(c**2, z, domain=sp.QQ)).monic().as_expr()
    rho = sp.Integer(1)
    multiplicities = []
    for factor, multiplicity in sp.sqf_list(g, z)[1]:
        rho *= factor**((multiplicity-1)//2)
        multiplicities.append(multiplicity)
    rho = sp.expand(rho)
    ell_bar = sp.cancel(ell/rho**2)
    c_bar = sp.cancel(c/rho)
    g_bar = sp.gcd(sp.Poly(ell_bar, z), sp.Poly(c_bar**2, z)).monic().as_expr()
    zero(g_bar-g/rho**2)
    assert all(m <= 2 for _, m in sp.sqf_list(g_bar, z)[1])
    assert sp.denom(ell_bar) == sp.denom(c_bar) == 1
    derivative_gcd = sp.gcd(sp.Poly(c, z), sp.Poly(sp.diff(c, z), z)).monic().as_expr()
    assert sp.rem(derivative_gcd, rho, z) == 0
    return g, rho, ell_bar, c_bar, multiplicities


def charge_one_period(marker):
    """Exact original-marker formula for Pdag, using Gamma times Gaussian."""
    R = 1+ss+ss**2+ss*y
    answer = 0
    for (iz, iw, iy), scalar in sp.Poly(sp.expand(marker), z, w, y).terms():
        if iw < iz:
            continue
        n = iw-iz
        value = 0
        for (is_, jy), coefficient in sp.Poly(sp.expand(ss**iw*y**iy*R**n), ss, y).terms():
            if jy % 2 == 0:
                value += coefficient*factorial(is_)*sp.Rational(
                    factorial(jy), 2**(jy//2)*factorial(jy//2))
        answer += scalar*value*t**n/factorial(n)
    return sp.expand(answer)


def checks():
    a = sp.symbols("a")
    cases = [
        ("high_odd_root_linear_reduced_ell", (z-1)**9, (z-1)**5),
        ("high_even_root", (z-1)**8, (z-1)**5),
        ("two_high_roots", z**5*(z-1)**7, z**3*(z-1)**4),
        ("finite_positive_origin", z**3*(z-1)**6, z**2*(z-1)**4),
        ("constant_nonzero_ell", sp.Integer(2), (z-1)**6*(z+2)**3),
        ("linear_ell_nonconstant_c", z+1, (z+1)**2),
        ("constant_nonzero_c", z**7*(z-1)**4, sp.Integer(1)),
        ("W_affine_boundary", sp.Integer(0), z**3*(z-1)**5),
    ]
    reductions = []
    for label, ell, c in cases:
        g, rho, ell_bar, c_bar, multiplicities = radical_marker(ell, c)
        original_q = c**2+(2*r-4*a)*ell
        reduced_q = c_bar**2+(2*r-4*a)*ell_bar
        zero(original_q-rho**2*reduced_q)
        effective_degree = max(sp.degree(ell_bar, z), 2*sp.degree(c_bar, z))
        assert effective_degree >= 2
        if sp.degree(c, z) >= 1:
            assert sp.degree(rho, z) <= sp.degree(c, z)-1
            assert sp.degree(c_bar, z) >= 1
        else:
            assert rho == 1
        reductions.append(dict(case=label, g_multiplicities=multiplicities,
                               rho=str(sp.factor(rho)),
                               reduced_ell_degree=str(sp.degree(ell_bar, z)),
                               reduced_c_degree=int(sp.degree(c_bar, z)),
                               reduced_raw_degree=int(effective_degree)))

    # Exact algebra in the formerly exceptional constant-e moving source.
    ee, aa, CC = sp.symbols("ee aa CC", nonzero=True)
    J0 = -2*ee+t*(4*aa*ee-CC)
    inverse_jet = -sum((t*(4*aa*ee-CC)/(2*ee))**j
                      for j in range(5))/(2*ee)
    zero(J0*inverse_jet-1+(t*(4*aa*ee-CC)/(2*ee))**5)

    # A fixed original polynomial-time representative of the escaping kernel.
    ell = z*(1+z)
    c, h, u = ell, sp.Integer(2), 1+z**2*(1+z)/4
    C, e = ell, sp.Integer(1)
    moving_R = C*(-z*C-2*e*h)+4*e**2*c*u
    zero(moving_R)
    phase = ell*w**2+w*(h+c*y)+u*y
    kappa_1 = -sp.diff(C, z)*z/2
    kappa_2 = sp.diff(C, z)*(ell*w+u)
    rows = periods(phase, [kappa_1, kappa_2], 4)
    for n in range(6):
        coefficient = (rows[0][n-1] if 1 <= n <= 5 else 0)
        coefficient += rows[1][n-2] if 2 <= n <= 5 else 0
        zero(coefficient)

    # Leading original rho period can be nonzero despite zero initial mass.
    ell = z**5*(z-1)**7
    c = z**3*(z-1)**4
    _, rho, _, _, _ = radical_marker(ell, c)
    outside = ell*w**2+w*(2+z+c*y)+y**2/3+(1+z)*y+z**2
    leading = periods(outside, [rho], 2)[0]
    assert leading == [0, 0, -4]

    # A genuinely weighted MOVING source: raw gcd z^4 is outside the old
    # unweighted volume range; rho=z reduces it to z^2.
    moving_phase = z**5*w**2-w/2+z**2*w*y
    moving_marker_coefficients = periods(moving_phase, [z], 5)[0]
    assert moving_marker_coefficients == [0, -sp.Rational(1, 2), 0, 0, 0,
                                         sp.Rational(5, 8)]
    S = sp.symbols("S")
    # x=2t(beta+x)^4 S^2, beta=-t/2: this is its exact five-jet.
    selected_z = -t/2+t**5*S**2/8
    selected_J = 1-8*t*selected_z**3*S**2
    selected_jet = sp.series(selected_z/selected_J, t, 0, 6).removeO()
    source_marker_period = wick(selected_jet.subs(S, y))
    zero(source_marker_period-sum(value*t**n for n, value
                                  in enumerate(moving_marker_coefficients)))

    # Finite branch with BOTH a triple origin and a multiplicity-six other root.
    ell = z**3*(z-1)**6
    c = z**2*(z-1)**4
    h = z**2*(z-1)**3
    _, rho, _, _, _ = radical_marker(ell, c)
    plus = ell*w**2+w*(h+c*y)+z*y+z**2
    assert all(iz-iw >= 1 for (iz, iw, iy), scalar
               in sp.Poly(sp.expand(plus), z, w, y).terms())
    order = 6
    u0, k0 = sp.Integer(2), sp.Integer(3)
    phase = plus+u0*y+k0
    markers = [sp.Integer(1), w, y, z*w, w**2, w*y, rho]
    rows = periods(phase, markers, order)
    exp_psi = exponential_coefficients(k0*t+u0**2*t**2/2, order)
    normalized = []
    for marker, row in zip(markers, rows):
        degree = sp.Poly(marker, z, w, y).total_degree()
        polynomial = sp.expand(sum(t**n*wick(
            sp.expand((marker*plus**n).subs(y, y+u0*t)))/factorial(n)
            for n in range(degree+1)))
        assert polynomial == 0 or sp.degree(polynomial, t) <= 2*degree
        expected = [sp.expand(sum(exp_psi[n-j]*polynomial.coeff(t, j)
                                  for j in range(n+1))) for n in range(order+1)]
        for actual, target in zip(row, expected):
            zero(actual-target)
        normalized.append(str(polynomial))
    assert normalized[0] == "1" and normalized[-1] == "0"
    assert rows[-1] == [0]*(order+1)

    safe_rows = periods(plus, markers[:-1], order)
    for marker, row in zip(markers[:-1], safe_rows):
        degree = sp.Poly(marker, z, w, y).total_degree()
        assert all(value == 0 for value in row[degree+1:])

    # Exact charge-one benchmark; the all-degree proof is positivity plus
    # angular balance, while these diagonal entries are direct raw controls.
    benchmark = z+z**2*w+z**3*w**2+z**2*w*y
    benchmark_rows = periods(benchmark, [w**n for n in range(5)], 6)
    diagonal = []
    for n, row in enumerate(benchmark_rows):
        assert all(value == 0 for m, value in enumerate(row) if m != n)
        assert row[n].is_Integer and row[n] >= 1
        diagonal.append(int(row[n]))
    assert diagonal == [1, 9, 535, 89593, 29527341]
    benchmark_markers = [z**iz*w**iw*y**iy
                         for iz in range(4) for iw in range(4-iz)
                         for iy in range(4-iz-iw)]
    benchmark_controls = periods(benchmark, benchmark_markers, 4)
    for marker, row in zip(benchmark_markers, benchmark_controls):
        polynomial = charge_one_period(marker)
        assert row == [polynomial.coeff(t, n) for n in range(5)]
    kernel_element = 2*w-9*w*y
    zero(charge_one_period(kernel_element))
    zero(charge_one_period(w*kernel_element)+280*t**2)
    assert charge_one_period(w*kernel_element) != 0

    return dict(status="passed", reduced_cover_cases=reductions,
                constant_e_inverse="generic five-term t-geometric identity",
                original_escaping_kernel_coefficients_0_to_5="all zero",
                high_multiplicity_nonzero_marker_coefficients=[int(v) for v in leading],
                weighted_moving_marker_coefficients=[str(v) for v in moving_marker_coefficients],
                weighted_moving_original_vs_selected_source="exact through order5, including zero initial mass",
                finite_branch_original_markers=[str(v) for v in markers],
                finite_branch_normalized_periods=normalized,
                original_Wick_moment_orders=[0, order],
                safe_original_cutoff="six markers vanish for m>deg Q through moment6",
                charge_one_benchmark_diagonal_periods=diagonal,
                charge_one_original_monomial_controls=len(benchmark_markers),
                charge_one_kernel_control="2W-9WY maps to0; W(2W-9WY) maps to -280t^2",
                scope_controls=["rho may vanish identically in the finite positive-angular branch",
                                "linear reduced ell occurs and retains effective cover degree2"],
                scope="Bounded exact algebraic and original-Wick controls; universal classification uses the complete written source and period-comparison proof.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    names = ["GAUSSIAN_QUADRATIC_W_MIXED_FINITE_IMAGES.md", Path(__file__).name,
             "GAUSSIAN_AFFINE_ISOTROPIC_CHART.md",
             "GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md",
             "GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md",
             "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
             "GAUSSIAN_UNIT_ELL_FORMAL_SOURCE_FINITE_IMAGE_EXCLUSION.md",
             "GAUSSIAN_SIMPLE_ORIGIN_ALTERNATE_SCALAR_EXCLUSION.md",
             "GAUSSIAN_MOVING_COMMON_ROOT_SCALAR_EXCLUSION.md",
             "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md"]
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
