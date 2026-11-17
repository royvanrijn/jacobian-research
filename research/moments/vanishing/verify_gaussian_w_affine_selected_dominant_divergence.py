"""Exact bounded controls for selected-dominant W-affine divergence.

The unrestricted theorem is the written selected-sheet, contour, and Gaussian
convolution proof. These algebraic controls do not certify its analytic
hypotheses for arbitrary inputs or infer divergence from moment prefixes.
"""

import argparse
import hashlib
import json
from math import factorial, gcd
from pathlib import Path
import platform

import sympy as sp

z, w, y, t = sp.symbols("z w y t")


def zero(expression):
    assert sp.simplify(sp.expand(expression)) == 0, expression


def ordinary_moment(exponent):
    if exponent % 2:
        return 0
    return factorial(exponent) // (2 ** (exponent // 2) * factorial(exponent // 2))


def wick(poly):
    result = 0
    for (iz, iw, iy), scalar in sp.Poly(sp.expand(poly), z, w, y).terms():
        if iz == iw:
            result += scalar * factorial(iz) * ordinary_moment(iy)
    return sp.expand(result)


def circular_wick(poly):
    result = 0
    for (iz, iw), scalar in sp.Poly(sp.expand(poly), z, w).terms():
        if iz == iw:
            result += scalar * factorial(iz)
    return sp.expand(result)


def positive_aperiodic_representative(poly):
    """Exact positive-rational representative test, not a general classifier."""
    polynomial = sp.Poly(poly, z, domain=sp.QQ)
    if polynomial.degree() < 2 or polynomial.nth(0) <= 0:
        return False
    if any(value < 0 for value in polynomial.all_coeffs()):
        return False
    support = [j for j in range(1, polynomial.degree() + 1)
               if polynomial.nth(j) != 0]
    period = 0
    for j in support:
        period = gcd(period, j)
    return period == 1


def root_order(poly, center):
    polynomial = sp.Poly(poly, z)
    for j in range(polynomial.degree() + 1):
        if sp.simplify(sp.diff(poly, z, j).subs(z, center)) != 0:
            return j
    raise AssertionError("zero polynomial has no finite root order")


def original_wick_coefficients(phase, order):
    power = sp.Integer(1)
    result = []
    for n in range(order + 1):
        if n:
            power = sp.expand(power * phase)
        result.append(sp.expand(wick(power) / factorial(n)))
    return result


def completed_coefficients(c, h, fibre_a, u, k, order):
    # After exact Y completion the exponent is
    # t(k+hW) + t^2(cW+u)^2/[2(1-2at)].
    q = [sp.Integer(0), k + h * w]
    for j in range(2, order + 1):
        q.append(sp.expand((c * w + u) ** 2 * (2 * fibre_a) ** (j - 2) / 2))
    exponential = [sp.Integer(1)]
    for n in range(1, order + 1):
        exponential.append(sp.expand(sum(j * q[j] * exponential[n - j]
                                         for j in range(1, n + 1)) / n))
    # (1-2at)^(-1/2) has this exact binomial coefficient.
    prefactor = [sp.binomial(2 * j, j) * (fibre_a / 2) ** j
                 for j in range(order + 1)]
    return [sp.expand(sum(prefactor[j] * circular_wick(exponential[n - j])
                          for j in range(n + 1)))
            for n in range(order + 1)]


def checks():
    c = (1 + z) * (1 + 2 * z)
    inverse = z / c
    z_plus, z_minus = sp.sqrt(2) / 2, -sp.sqrt(2) / 2
    x_plus, x_minus = 3 - 2 * sp.sqrt(2), 3 + 2 * sp.sqrt(2)
    zero(c - z * sp.diff(c, z) - (1 - 2 * z ** 2))
    for center, value in ((z_plus, x_plus), (z_minus, x_minus)):
        zero(sp.diff(inverse, z).subs(z, center))
        zero(inverse.subs(z, center) - value)
        assert sp.simplify(sp.diff(inverse, z, 2).subs(z, center)) != 0
    assert bool(x_plus > 0) and bool(x_minus > x_plus)
    assert positive_aperiodic_representative(c)

    equality_c = (1 + z) ** 4 * (1 + 2 * z)
    equality_h = (1 + z) ** 2
    tau = (sp.sqrt(41) - 3) / 16
    chi = equality_c - z * sp.diff(equality_c, z)
    zero(chi - (1 + z) ** 3 * (1 - 3 * z - 8 * z ** 2))
    zero((1 - 3 * z - 8 * z ** 2).subs(z, tau))
    assert bool(tau > 0)
    assert root_order(equality_c, -1) == 4
    assert root_order(equality_h, -1) == 2
    assert root_order(equality_c, -sp.Rational(1, 2)) == 1
    assert positive_aperiodic_representative(equality_c)
    # Generic derivative identity gives G''=-z*c''/c^2 at a genuine critical.
    g = z / equality_c
    zero(equality_c ** 3 * sp.diff(g, z, 2)
         + z * equality_c * sp.diff(equality_c, z, 2)
         + 2 * sp.diff(equality_c, z) * chi)
    assert sp.simplify(sp.diff(chi, z).subs(z, tau)) != 0

    rejected = {
        "c0_zero": z + z ** 2,
        "degree_below_two": 1 + z,
        "positive_support_gcd_two": 1 + z ** 2,
    }
    for poly in rejected.values():
        assert not positive_aperiodic_representative(poly)
    # The period-two example also has tied critical moduli.
    tied = z / (1 + z ** 2)
    zero(sp.diff(tied, z).subs(z, 1))
    zero(sp.diff(tied, z).subs(z, -1))
    zero(tied.subs(z, 1) + tied.subs(z, -1))

    fibre_a = sp.Rational(2, 3)
    h, u, k = 1 - z + z ** 3, 2 + z, z ** 2 - z + sp.Rational(1, 5)
    phase = w * (h + c * y) + fibre_a * y ** 2 + u * y + k
    order = 4
    direct = original_wick_coefficients(phase, order)
    completed = completed_coefficients(c, h, fibre_a, u, k, order)
    assert direct == completed

    # Generic Y completion and the exact linear-Gaussian auxiliary coefficient.
    cs, hs, us, ks, aa, xi = sp.symbols("cs hs us ks aa xi")
    d = 1 - 2 * aa * t
    square = t * hs * w + t * ks + (t * cs * w + t * us + xi) ** 2 / (2 * d)
    A = t ** 2 * cs ** 2 / (2 * d)
    B = t * hs + t * cs * (t * us + xi) / d
    amplitude = t * ks + (t * us + xi) ** 2 / (2 * d)
    zero(square - (A * w ** 2 + B * w + amplitude))
    zero((t * cs / sp.sqrt(d)) ** 2 - 2 * A)

    n = sp.symbols("n", positive=True, integer=True)
    bn = sp.binomial(2 * n, n) / 4 ** n
    bn_ratio = sp.combsimp(bn / bn.subs(n, n - 1))
    zero(bn_ratio - (2 * n - 1) / (2 * n))
    zero(sp.limit(n * (bn_ratio - 1), n, sp.oo) + sp.Rational(1, 2))
    radius = sp.symbols("radius", positive=True)
    ln = (sp.factorial(2 * n) / (2 ** n * sp.factorial(n))
          * radius ** (-2 * n) * sp.binomial(4 * n, 2 * n) / 4 ** (2 * n))
    ln_ratio = sp.combsimp(ln / ln.subs(n, n - 1))
    zero(sp.limit(ln_ratio / n, n, sp.oo) - 2 / radius ** 2)

    eta, aux = sp.symbols("eta aux", nonzero=True)
    zero(eta * (t / eta) * aux - t * aux)

    return {
        "status": "passed",
        "written_theorem_scope": (
            "c(0) nonzero, degree at least two; selected unique dominant simple "
            "inverse singularity. Sufficient: simple genuine finite criticals "
            "with distinct critical-value moduli, or phase-scaled nonnegative "
            "aperiodic coefficients. Arbitrary constant fibre a and polynomial h,u,k."
        ),
        "bounded_controls": {
            "distinct_root_example": "c=(1+z)(1+2z); critical values 3+-2sqrt(2)",
            "equality_example": "c=(1+z)^4(1+2z), h=(1+z)^2; orders (4,2) at -1",
            "equality_positive_critical": "(sqrt(41)-3)/16",
            "original_wick_vs_completion_orders": [0, order],
            "arbitrary_perturbation_example": {
                "a": str(fibre_a), "h": str(h), "u": str(u), "k": str(k),
                "unit_coefficients": [str(value) for value in direct],
            },
            "generic_algebra": [
                "Y completion and auxiliary Gaussian coefficient",
                "central-binomial exact ratio and first asymptotic correction",
                "Gaussian-majorant ratio limit",
                "critical-value reparametrization argument identity",
            ],
        },
        "negative_controls": list(rejected),
        "verification_scope": (
            "Finite exact algebraic controls only. The unrestricted selected-sheet "
            "continuation, uniform contour bounds, generic density, and divergence "
            "are the written proof; no moment-prefix inference."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    result.update(
        family="W-affine selected-dominant multiroot divergence",
        python=platform.python_version(),
        sympy=sp.__version__,
        proof_sha256=digest(directory / "GAUSSIAN_W_AFFINE_SELECTED_DOMINANT_DIVERGENCE.md"),
        checker_sha256=digest(Path(__file__)),
        source_proof_sha256=digest(directory / "GAUSSIAN_W_AFFINE_SINGLE_ROOT_DIVERGENCE.md"),
        independent_replay=False,
        formal_verification=False,
        external_review=False,
    )
    report = json.dumps(result, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(report, encoding="utf-8", newline="\n")
    print(report, end="")


if __name__ == "__main__":
    main()
