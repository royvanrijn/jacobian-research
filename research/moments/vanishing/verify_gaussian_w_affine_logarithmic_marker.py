"""Finite exact controls for the original logarithmic-marker argument.

The written fixed-puncture comparison proof, not these bounded controls,
establishes the unrestricted finite-image classification. Run --output PATH
to retain an internal receipt pinned to the current proof and checker.
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
    assert sp.simplify(sp.expand(expression)) == 0, expression


def even_gaussian_power(power):
    if power % 2:
        return 0
    return factorial(power)//(2**(power//2)*factorial(power//2))


def wick(poly):
    answer = 0
    for (iz, iw, iy), scalar in sp.Poly(sp.expand(poly), z, w, y).terms():
        if iz == iw:
            answer += scalar*factorial(iz)*even_gaussian_power(iy)
    return sp.expand(answer)


def period_coefficients(marker, phase, order):
    power = sp.Integer(1)
    answer = []
    for n in range(order+1):
        if n:
            power = sp.expand(power*phase)
        answer.append(sp.expand(wick(marker*power)/factorial(n)))
    return answer


def radical_data(c):
    poly = sp.Poly(c, z, domain=sp.QQ)
    rho = sp.gcd(poly, poly.diff()).monic().as_expr()
    red = sp.cancel(c/rho)
    zero(rho*red-c)
    zero(rho/c-1/red)
    assert sp.Poly(red, z).degree() >= 1
    assert sp.gcd(sp.Poly(red, z), sp.Poly(sp.diff(red, z), z)).degree() == 0
    assert sp.Poly(red, z).degree()-2 >= -1
    return rho, red


def checks():
    cc, dd, rrho, JJ, hh, aa, uu, kk = sp.symbols(
        "cc dd rrho JJ hh aa uu kk", nonzero=True)
    zeta_y = t*cc/JJ
    zero(zeta_y/(t*cc)-1/JJ)
    zero((zeta_y/(t*dd)).subs(cc, rrho*dd)-rrho/JJ)
    ss = (z/t-hh)/cc
    phase = -ss**2/2+t*(aa*ss**2+uu*ss+kk)
    zero((phase+y**2/2-t*(aa*y**2+uu*y+kk)).subs(z, t*(hh+cc*y)))

    # Generic leading-root identity for every integer s>=0.
    s = sp.symbols("s", nonnegative=True, integer=True)
    h0, rho_s, q1, q2 = sp.symbols("h0 rho_s q1 q2", nonzero=True)
    leading_quotient = (h0+t*q1)**s*(rho_s+t*q2)
    zero(leading_quotient.subs(t, 0)-rho_s*h0**s)

    cases = []
    for label, c in [
        ("mixed_nonzero_root_multiplicities_4_2", (z-1)**4*(z+2)**2),
        ("origin_order4_and_nonzero_order2", z**4*(z-1)**2),
        ("squarefree", (z-1)*(z+2)),
        ("one_multiple_root", (z-1)**4),
    ]:
        rho, red = radical_data(c)
        cases.append({"case": label, "c": str(sp.expand(c)),
                      "rho": str(sp.factor(rho)), "c_red": str(sp.factor(red)),
                      "infinity_order_of_dz_over_c_red": int(sp.degree(red, z)-2)})

    c = (z-1)**4*(z+2)**2
    h = (z-1)**2*(z+2)
    rho, _ = radical_data(c)
    mixed_initial = wick(rho)
    zero(mixed_initial+2)
    # The nonzero root 1 has precisely the formerly exceptional orders (4,2).
    assert sp.Poly(c, z).diff((z, 3)).eval(1) == 0
    assert sp.Poly(c, z).diff((z, 4)).eval(1) != 0
    assert sp.diff(h, z).subs(z, 1) == 0
    assert sp.diff(h, z, 2).subs(z, 1) != 0

    # Original Wick control where rho(0)=0 but the source period is nonzero.
    H0 = sp.symbols("H0", nonzero=True)
    c = z**4*(z-1)**2
    rho, _ = radical_data(c)
    phase = w*(H0+z+c*y)
    leading = period_coefficients(rho, phase, 3)
    assert leading[:3] == [0, 0, 0]
    zero(leading[3]+H0**3)

    # Shared zero origin is compatible with a nonzero (4,2) root.
    c = z**2*(z-1)**4
    h = z**2*(z-1)**2
    rho, _ = radical_data(c)
    phase = w*(h+c*y)
    assert all(iz-iw >= 1 for (iz, iw, iy), scalar
               in sp.Poly(sp.expand(phase), z, w, y).terms())
    assert all(iz-iw >= 1 for (iz, iw, iy), scalar
               in sp.Poly(sp.expand(rho), z, w, y).terms())
    assert period_coefficients(1, phase, 4) == [1, 0, 0, 0, 0]
    assert period_coefficients(rho, phase, 4) == [0, 0, 0, 0, 0]

    normalization_cases = []
    for p in [sp.Integer(2), 1+t, t**3, t**3*(1+2*t)]:
        HH = sp.cancel(1/(t*p))
        numerator, denominator = sp.fraction(HH)
        order = sp.degree(denominator, t)-sp.degree(numerator, t)
        assert order == 1+sp.degree(p, t)
        normalization_cases.append({"p": str(p), "infinity_order_r_over_p": int(order)})

    return {
        "status": "passed",
        "generic_identities": [
            "eta pulls back to dY/J",
            "rho*eta pulls back to rho(zeta)*dY/J",
            "phase pulls back to Gaussian phase plus original fibre exponent",
            "all-s leading coefficient rho_s*h0^s",
        ],
        "radical_and_logarithmic_cases": cases,
        "mixed_4_2_initial_marker_period": str(mixed_initial),
        "origin_order4_h0_nonzero_original_coefficients_0_to_3": [str(v) for v in leading],
        "shared_origin_compatible_4_2_example": "positive angular support checked; original unit/marker Wick controls through order4",
        "normalizations_including_zero_initial_mass": normalization_cases,
        "scope": "Finite exact algebraic and original Wick controls. The written fixed-log comparison proof establishes unrestricted root orders.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    result.update(
        python=platform.python_version(), sympy=sp.__version__,
        proof_sha256=digest(directory / 'GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md'),
        checker_sha256=digest(Path(__file__)),
        fixed_log_import_sha256=digest(directory / 'GAUSSIAN_AFFINE_ISOTROPIC_CHART.md'),
        independent_replay=False, formal_verification=False, external_review=False,
    )
    report = json.dumps(result, indent=2) + '\n'
    if args.output is not None:
        args.output.write_text(report, encoding='utf-8', newline='\n')
    print(report, end='')
