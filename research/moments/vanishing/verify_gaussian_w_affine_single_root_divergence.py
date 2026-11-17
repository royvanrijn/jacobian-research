"""Exact algebraic controls for the single-nonzero-root divergence proof.

The selected-sheet/contour argument is a written proof, not certified by
these finite controls. No moment prefix is used as an exclusion criterion.
"""

from fractions import Fraction
from math import comb, factorial
import argparse
import hashlib
import json
from pathlib import Path
import platform

import sympy as sp

z, w, y, t, xi, a = sp.symbols("z w y t xi a")


def zero(expression):
    assert sp.simplify(sp.expand(expression)) == 0, expression


def dfac_even(exponent):
    if exponent % 2:
        return 0
    return factorial(exponent)//(2**(exponent//2)*factorial(exponent//2))


def wick(poly):
    answer = 0
    for (iz, iw, iy), scalar in sp.Poly(sp.expand(poly), z, w, y).terms():
        if iz == iw:
            answer += scalar*factorial(iz)*dfac_even(iy)
    return sp.expand(answer)


def prototype_lagrange(order):
    values = []
    for N in range(order+1):
        value = Fraction(0)
        for n in range(N+1):
            for k in range(n+1):
                j = n+k-N
                if j < 0 or j > k or j % 2:
                    continue
                value += (Fraction((-1)**n*comb(n, k)*comb(k, j)*comb(2*n+2*k, n))
                          *Fraction(-1, 4)**(k-j)*dfac_even(j))
        values.append(sp.Rational(value.numerator, value.denominator))
    return values


def checks():
    c, h, u, k = sp.symbols("c h u k")
    d = 1-2*a*t
    square = t*h*w+t*k+(t*c*w+t*u+xi)**2/(2*d)
    A = t**2*c**2/(2*d)
    B = t*h+t*c*(t*u+xi)/d
    sigma = t*k+(t*u+xi)**2/(2*d)
    zero(square-(A*w**2+B*w+sigma))
    zero((t*c/sp.sqrt(d))**2-2*A)

    alpha, beta = sp.symbols("alpha beta", nonzero=True)
    for m in range(2, 8):
        inverse = z/(alpha+beta*z)**m
        zs = alpha/((m-1)*beta)
        rho = sp.Rational((m-1)**(m-1), m**m)
        gamma = beta*alpha**(m-1)
        zero(sp.diff(inverse, z).subs(z, zs))
        zero(inverse.subs(z, zs)-rho/gamma)
        zero(sp.diff(inverse, z, 2).subs(z, zs)
             +(m-1)*beta*(alpha*m/(m-1))**(-m-1))
        for n in range(7):
            zero(sp.expand((alpha+beta*z)**(m*n)).coeff(z, n)
                 -sp.binomial(m*n, n)*gamma**n)

    p = (z-1)**2*w+(z-1)**4*w*y-y/4
    order = 6
    power = sp.Integer(1)
    raw = 0
    for n in range(order+1):
        if n:
            power = sp.expand(power*p)
        raw += wick(power)*t**n/factorial(n)
    normalized = sp.series(sp.exp(-t**2/32)*raw, t, 0, order+1).removeO().expand()
    lagrange = prototype_lagrange(order)
    assert [normalized.coeff(t, n) for n in range(order+1)] == lagrange
    assert lagrange[:4] == [1, -2, 35, -sp.Rational(775, 2)]

    wc = sp.symbols("wc", nonzero=True)
    prototype_inverse = z/(z-1)**4-t/(z-1)**2+t**2/4
    critical_t = (3*wc+4)/(2*wc**2)
    zero(sp.diff(prototype_inverse, z).subs({z: wc+1, t: critical_t}))
    zero(prototype_inverse.subs({z: wc+1, t: critical_t})
         -(9*wc+16)/(16*wc**3))

    n = sp.symbols("n", positive=True, integer=True)
    bn = sp.factorial(4*n)/(sp.factorial(n)*sp.factorial(3*n))
    ratio = sp.combsimp(bn/bn.subs(n, n-1))
    rho = sp.Rational(27, 256)
    zero(sp.limit(ratio, n, sp.oo)-1/rho)
    zero(sp.limit(n*(rho*ratio-1), n, sp.oo)+sp.Rational(1, 2))

    return {
        "status": "passed",
        "generic_source": "Y square completion and exact linear-Gaussian kernel coefficient",
        "critical_geometry_controls_m": [2, 7],
        "base_lagrange_controls_n": [0, 6],
        "prototype_independent_wick_lagrange_orders": [0, 6],
        "prototype_first_coefficients": [str(value) for value in lagrange[:4]],
        "prototype_critical_ratio": "rho=27/256; first correction -1/(2n)",
        "prototype_critical_curve": "2t*w^2-3w-4=0; Xstar=(9w+16)/(16w^3)",
        "scope": "Finite algebraic controls. Unrestricted divergence is the written selected-sheet contour/convolution proof.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    result.update(
        family="W-affine single-nonzero-root divergence and one-root finite classification",
        python=platform.python_version(), sympy=sp.__version__,
        proof_sha256=digest(directory / 'GAUSSIAN_W_AFFINE_SINGLE_ROOT_DIVERGENCE.md'),
        checker_sha256=digest(Path(__file__)),
        independent_replay=False, formal_verification=False, external_review=False,
    )
    report = json.dumps(result, indent=2) + '\n'
    if args.output is not None:
        args.output.write_text(report, encoding='utf-8', newline='\n')
    print(report, end='')
