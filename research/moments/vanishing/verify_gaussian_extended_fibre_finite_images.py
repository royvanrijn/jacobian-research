"""Bounded exact source and Wick controls for extended Gaussian fibres.

The universal classifications use the written comparison proofs. These
controls verify original-marker identities and representative boundary cases.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp

z, w, y, t = sp.symbols("z w y t")


def zero(expression):
    assert sp.cancel(sp.expand(expression)) == 0, expression


def wick(poly):
    answer = 0
    for (iz, iw, iy), scalar in sp.Poly(sp.expand(poly), z, w, y).terms():
        if iz == iw and iy % 2 == 0:
            answer += scalar * factorial(iz) * sp.Rational(
                factorial(iy), 2 ** (iy // 2) * factorial(iy // 2))
    return sp.expand(answer)


def periods(phase, markers, order):
    answer = [[] for _ in markers]
    power = sp.Integer(1)
    for n in range(order + 1):
        if n:
            power = sp.expand(power * phase)
        for row, marker in zip(answer, markers):
            row.append(sp.expand(wick(marker * power) / factorial(n)))
    return answer


def coefficients(expression, order):
    return [sp.expand(expression).coeff(t, n) for n in range(order + 1)]


def exponential_coefficients(exponent, order):
    # Formal recurrence n a_n=sum_j j psi_j a_(n-j).
    answer = [sp.Integer(1)]
    for n in range(1, order + 1):
        answer.append(sp.expand(sum(j * exponent.coeff(t, j) * answer[n-j]
                                    for j in range(1, n+1)) / n))
    return answer


def multiply_truncated(left, right, order):
    return [sp.expand(sum(left[j] * right[n-j] for j in range(n+1)))
            for n in range(order+1)]


def checks():
    order = 6
    eta, gamma = sp.Integer(2), sp.Integer(-1)
    h, c = eta*z + z**2, gamma*z + z**3
    g0 = y**3 + y**2 + y
    phase = w*(h+c*y) + g0 + z*y**3 + z**2*y
    source_markers = [sp.Integer(1), y, y**2, y**3]
    actual = periods(phase, source_markers, order)
    ordinary = periods(g0, source_markers[:3], order)
    for j in range(3):
        for n in range(order+1):
            lhs = actual[j][n]
            if n:
                lhs -= eta*actual[j][n-1] + gamma*actual[j+1][n-1]
            zero(lhs-ordinary[j][n])
    assert actual[0][1] != ordinary[0][1]  # Omitting J0 must fail.

    # Nonquadratic, variable-Y-degree fibre: positive angular expansion.
    plus = w*(z**2+z**3*y) + z*y**4 + z**2*y**3
    u, k = sp.Integer(2), sp.Integer(3)
    phase = plus+u*y+k
    markers = [sp.Integer(1), w, y, z*w, w**2, w*y]
    actual = periods(phase, markers, order)
    psi = k*t+u**2*t**2/2
    exp_psi = exponential_coefficients(psi, order)
    normalized = []
    for marker, row in zip(markers, actual):
        degree = sp.Poly(marker, z, w, y).total_degree()
        amplitude = 0
        for n in range(degree+1):
            shifted = sp.expand((marker*plus**n).subs(y, y+u*t))
            amplitude += t**n*wick(shifted)/factorial(n)
        amplitude = sp.expand(amplitude)
        expected = multiply_truncated(exp_psi, coefficients(amplitude, order), order)
        for a, b in zip(row, expected):
            zero(a-b)
        assert sp.degree(amplitude, t) <= 5*degree if degree else amplitude == 1
        normalized.append(str(amplitude))
    zero(actual[0][0]-1)
    safe_rows = periods(plus, markers, order)
    for marker, row in zip(markers, safe_rows):
        degree = sp.Poly(marker, z, w, y).total_degree()
        assert all(value == 0 for value in row[degree+1:])
    assert safe_rows[1][1] != 0

    # Constant nonzero W^2: direct t-adic determinant localization.
    lam, aa, cc = sp.symbols("lam aa cc", nonzero=True)
    delta = -2*t*lam+t**2*(4*aa*lam-cc**2)
    inverse = -sum((-t*(cc**2-4*aa*lam)/(2*lam))**j
                   for j in range(6))/(2*t*lam)
    remainder = sp.expand(delta*inverse-1)
    assert all(remainder.coeff(t, n) == 0 for n in range(6))
    hh, uu, kk = sp.symbols("hh uu kk")
    precision = sp.Matrix([[-2*lam*t, -cc*t], [-cc*t, 1-2*aa*t]])
    source = sp.Matrix([t*hh-z, t*uu])
    numerator = (z**2+t*(-2*z*hh-2*aa*z**2)
                 +t**2*(hh**2+4*aa*z*hh-2*cc*z*uu-4*lam*kk)
                 +t**3*(-2*aa*hh**2+2*cc*hh*uu-2*lam*uu**2
                        +2*(4*lam*aa-cc**2)*kk))
    zero(precision.det()-delta)
    zero(2*delta*(t*kk+(source.T*precision.inv()*source)[0]/2)-numerator)

    # An isotropic square with lambda=1 and nonconstant h,u,k.
    ell = w+y-z/2
    quadratic = ell**2+y+3
    covariance = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    xs = sp.Matrix([z, w, y])
    matrix = sp.hessian(quadratic, (z, w, y))
    assert matrix[1, 1] == 2
    assert covariance*matrix != sp.zeros(3)
    assert (covariance*matrix)**2 == sp.zeros(3)
    vector = sp.Matrix([sp.diff(quadratic, x).subs({z: 0, w: 0, y: 0})
                        for x in xs])
    mean = t*(sp.eye(3)+t*covariance*matrix)*covariance*vector
    tilted_covariance = (sp.eye(3)+t*covariance*matrix)*covariance
    exponent = 3*t+t**2*(vector.T*(sp.eye(3)+t*covariance*matrix)
                         *covariance*vector)[0]/2
    zero(exponent-(3*t+t**2/2+t**3))
    qmarkers = [sp.Integer(1), z, w, y, w**2, z*y]
    qrows = periods(quadratic, qmarkers, order)
    qamps = [1, mean[0], mean[1], mean[2],
             tilted_covariance[1, 1]+mean[1]**2,
             tilted_covariance[0, 2]+mean[0]*mean[2]]
    qexp = exponential_coefficients(sp.expand(exponent), order)
    for row, amplitude in zip(qrows, qamps):
        expected = multiply_truncated(qexp, coefficients(sp.sympify(amplitude), order), order)
        for a, b in zip(row, expected):
            zero(a-b)

    quartic = w**2+2*w*y+4*z*y-4*z**2+z
    quartic_b = covariance*sp.hessian(quartic, (z, w, y))
    assert quartic_b**3 == sp.zeros(3)
    assert quartic_b**2 != sp.zeros(3)
    quartic_vector = sp.Matrix([1, 0, 0])
    quartic_exponent = sp.expand(t**2*(quartic_vector.T*
        (sp.eye(3)+t*quartic_b+t**2*quartic_b**2)*covariance*quartic_vector)[0]/2)
    zero(quartic_exponent-(t**3+2*t**4))
    quartic_row = periods(quartic, [sp.Integer(1)], order)[0]
    assert quartic_row == exponential_coefficients(quartic_exponent, order)

    # Boundary control: the unit is 1/(1-t), not an exponential polynomial.
    bad_quadratic = w**2+z*w
    bad = periods(bad_quadratic, [sp.Integer(1)], order)[0]
    assert bad == [1]*(order+1)
    assert (covariance*sp.hessian(bad_quadratic, (z, w, y)))**3 != sp.zeros(3)

    return dict(status="passed", order=order,
                original_source_cancellation="3 markers, including cubic G(0,Y), through order6",
                omitted_J0_negative_control="rejected",
                arbitrary_positive_fibre_normalized_periods=normalized,
                positive_fibre_original_markers=len(markers),
                safe_fibre_cutoff="m>deg Q vanishes for six markers through order6",
                determinant_localization="symbolic inverse through relative order5",
                schur_numerator="symbolic precision inverse matches all coefficients",
                primitive_quadratic_original_markers=len(qmarkers),
                primitive_quadratic_exponent=str(exponent),
                maximal_degree_primitive_exponent=str(quartic_exponent),
                nonnilpotent_quadratic_unit_coefficients=[int(value) for value in bad],
                scope="Bounded exact controls only; universal conclusions use retained written comparison proofs.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    names = ["GAUSSIAN_W_AFFINE_ARBITRARY_FIBRE_FINITE_IMAGES.md",
             "GAUSSIAN_CONSTANT_W_QUADRATIC_FINITE_IMAGES.md",
             "GAUSSIAN_AFFINE_ISOTROPIC_CHART.md",
             "GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md",
             "GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md",
             "MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md",
             "../../../unpublished/gaussian-boundary/polynomial_stationary_balls.tex",
             "../../../unpublished/gaussian-boundary/repeated_line_exclusion.tex",
             Path(__file__).name]
    result.update(python=platform.python_version(), sympy=sp.__version__,
                  sha256={name: hashlib.sha256((directory/name).read_bytes()).hexdigest()
                          for name in names},
                  independent_replay=False, formal_verification=False, external_review=False)
    report = json.dumps(result, indent=2)+"\n"
    if args.output is not None:
        args.output.write_text(report, encoding="utf-8", newline="\n")
    print(report, end="")
