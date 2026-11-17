"""Finite exact replay of the proved original-kernel/filtration formulas.

These checks compare independent original Wick counts and heat formulas.
They do not prove unrestricted safety, image rank or the universal formulas.
Run with --output PATH to retain a finite internal verification receipt.
"""

from fractions import Fraction
from math import factorial
import argparse
import hashlib
import json
from pathlib import Path
import platform

import sympy as sp

z, w, y, t = sp.symbols("z w y t")
variables = (z, w, y)
inverse = {z: z, w: w-z*y-z**3/2, y: y+z**2}
forward = {z: z, w: w+z*y-z**3/2, y: y-z**2}


def monomials(bound):
    for r in range(bound+1):
        for s in range(bound+1-r):
            for u in range(bound+1-r-s):
                yield (r, s, u), z**r*w**s*y**u


def heat(poly, sign=1):
    total = sp.expand(poly)
    derivative = sp.expand(poly)
    j = 0
    while derivative != 0:
        j += 1
        derivative = sp.expand(2*sp.diff(derivative, z, w)+sp.diff(derivative, y, 2))
        if derivative != 0:
            total += sign**j*derivative/(2**j*factorial(j))
    return sp.expand(total)


def standard_heat(poly, variable, sign=1):
    total = sp.expand(poly)
    derivative = sp.expand(poly)
    j = 0
    while derivative != 0:
        j += 1
        derivative = sp.expand(sp.diff(derivative, variable, 2))
        if derivative != 0:
            total += sign**j*derivative/(2**j*factorial(j))
    return sp.expand(total)


def dfac_odd(exponent):
    # exponent is the even Gaussian power; return (exponent-1)!!.
    if exponent % 2:
        return 0
    return factorial(exponent)//(2**(exponent//2)*factorial(exponent//2))


def direct_coefficient(r, s, u, m):
    # Independent multinomial Wick expansion of (w+zy-z^3/2)^m/m!.
    answer = Fraction(0)
    for k in range(m+1):
        numerator = m-r+s-4*k
        if numerator % 2:
            continue
        j = numerator//2
        i = m-j-k
        if min(i, j) < 0 or r+j+3*k != s+i:
            continue
        moment_y = dfac_odd(u+j)
        if not moment_y:
            continue
        answer += Fraction((-1)**k*factorial(r+j+3*k)*moment_y,
                           2**k*factorial(i)*factorial(j)*factorial(k))
    return sp.Rational(answer.numerator, answer.denominator)


def direct_period(poly, order):
    coefficients = [sp.Rational(0) for _ in range(order+1)]
    for (r, s, u), scalar in sp.Poly(sp.expand(poly), *variables).terms():
        for m in range(order+1):
            coefficients[m] += scalar*direct_coefficient(r, s, u, m)
    return sp.Poly(sum(sp.expand(c)*t**m for m, c in enumerate(coefficients)), t)


def predicted_period(poly):
    pulled = sp.expand(poly.subs(inverse, simultaneous=True))
    return sp.Poly(heat(pulled).subs({z: t, w: 0, y: 0}), t)


def checks():
    originals = list(monomials(4))
    for _, marker in originals:
        predicted = predicted_period(marker)
        actual = direct_period(marker, 13)
        assert actual == predicted, (marker, actual, predicted)

    filtration = []
    for d in range(5):
        rows = [[predicted_period(marker).nth(j) for j in range(3*d+1)]
                for _, marker in monomials(d)]
        rank = sp.Matrix(rows).rank()
        assert rank == 3*d+1
        top = predicted_period(w**d).nth(3*d)
        assert top == (-sp.Rational(1, 2))**d
        filtration.append({"marker_degree": d, "image_dimension": rank,
                           "max_time_degree": 3*d})

    hermite_count = 0
    max_transformed_degree = 0
    for (r, s, u), marker in monomials(3):
        H = heat(marker, -1)
        transformed = sp.expand(H.subs(forward, simultaneous=True))
        max_transformed_degree = max(max_transformed_degree,
                                     sp.Poly(transformed, *variables).total_degree())
        expected = sp.Poly(t**r if s == u == 0 else 0, t)
        assert direct_period(transformed, 28) == expected
        hermite_count += 1

    inverse_count = 0
    for _, marker in monomials(5):
        assert heat(heat(marker, -1)) == marker
        assert heat(heat(marker), -1) == marker
        inverse_count += 1

    af, bf = sp.symbols("af bf")
    f = af*z+bf
    image_w, image_y = w+f*y-z*f**2/2, y-z*f
    assert sp.expand(2*z*image_w+image_y**2-2*z*w-y**2) == 0
    assert sp.expand((w-f*y-z*f**2/2).subs({w: image_w, y: image_y}, simultaneous=True)-w) == 0

    # Original ideal images for nonisotropic P=z and independent standard y.
    # The ideal (z+y) has content one: (z+y)*y*z^m lifts a Hermite polynomial.
    for m in range(6):
        h_m = standard_heat(z**m, z).subs(z, t)
        actual = 0
        for k in range(11):
            actual += sp.Rational(dfac_odd(m+k), factorial(k))*t**k
        predicted = sp.series(sp.exp(t**2/2)*h_m, t, 0, 11).removeO()
        assert sp.expand(actual-predicted) == 0
        # E_y[(z+y)*y*z^m]=z^m, including the mass-zero weighted seed.
        conditional = 0
        for (iz, iy), scalar in sp.Poly((z+y)*y*z**m, z, y).terms():
            conditional += scalar*z**iz*dfac_odd(iy)
        assert sp.expand(conditional-z**m) == 0

    # The ideal (z^2) has a two-dimensional defect under the heat map.
    ideal_polys = [standard_heat(z**m, z) for m in range(2, 9)]
    assert sp.Matrix([[poly.coeff(z, j) for j in range(9)]
                      for poly in ideal_polys]).rank() == 7
    for m, poly in enumerate(ideal_polys, 2):
        assert standard_heat(poly, z, -1) == z**m
    assert sp.rem(standard_heat(sp.Integer(1), z, -1), z**2, z) == 1
    assert sp.rem(standard_heat(z, z, -1), z**2, z) == z

    return {
        "status": "passed",
        "original_monomial_markers": len(originals),
        "marker_degrees_tested": [0, 4],
        "direct_wick_orders_tested": [0, 13],
        "exact_filtration_ranks": filtration,
        "transformed_hermite_markers": hermite_count,
        "hermite_input_degree": [0, 3],
        "max_transformed_marker_degree": max_transformed_degree,
        "transformed_hermite_wick_orders": [0, 28],
        "heat_inverse_monomials": inverse_count,
        "heat_inverse_degree": [0, 5],
        "generic_affine_shear": "quadratic preservation and inverse checked symbolically",
        "original_ideal_image_controls": {
            "transverse_content_one": "P=z,J=(z+y); fixed markers (z+y)*y*z^m, m=0..5, original Wick orders 0..10",
            "phase_content_degree_two": "P=z,J=(z^2); heat images of z^2..z^8 have rank 7 in degrees <=8; inverse-heat divisibility checked",
        },
        "scope": "Finite exact replay; universal theorem follows from the written identities.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    result.update(
        python=platform.python_version(), sympy=sp.__version__,
        proof_sha256=digest(directory / 'GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md'),
        checker_sha256=digest(Path(__file__)),
        independent_replay=False, formal_verification=False, external_review=False,
    )
    report = json.dumps(result, indent=2) + '\n'
    if args.output is not None:
        args.output.write_text(report)
    print(report, end='')
