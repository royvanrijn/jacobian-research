#!/usr/bin/env python3
"""Exact finite checks supporting SIMPLEX_INTERFACE_SYMMETRY.md.

No all-k induction theorem, analytic kernel inequality or global optimum
is certified by this script. The optional quadrature is exploration only.
"""

import argparse
from fractions import Fraction as Q
import hashlib
import itertools
import json
import math
from pathlib import Path
from statistics import NormalDist


def parity(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def cycle_type(p):
    seen, lengths = set(), []
    for i in range(len(p)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = p[j]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def edge_trace(p, complement_sign=False):
    """Induced action on an oriented reference edge and its complement."""
    n, result = len(p), 0
    for i, j in itertools.combinations(range(n), 2):
        if {p[i], p[j]} != {i, j}:
            continue
        sign = 1 if p[i] < p[j] else -1
        if complement_sign:
            remaining = [a for a in range(n) if a not in (i, j)]
            sign *= parity([remaining.index(p[a]) for a in remaining])
        result += sign
    return result


def character_checks():
    # Actual edge traces are compared with independently specified irreducibles.
    s4 = {
        (1, 1, 1, 1): {"31": 3, "22": 2, "211": 3, "1111": 1},
        (2, 1, 1): {"31": 1, "22": 0, "211": -1, "1111": -1},
        (2, 2): {"31": -1, "22": 2, "211": -1, "1111": 1},
        (3, 1): {"31": 0, "22": -1, "211": 0, "1111": 1},
        (4,): {"31": -1, "22": 0, "211": 1, "1111": -1},
    }
    checked = 0
    for p in itertools.permutations(range(4)):
        chi = s4[cycle_type(p)]
        assert edge_trace(p) == chi["31"] + chi["211"]
        assert edge_trace(p, True) == chi["22"] + chi["211"] + chi["1111"]
        checked += 1
    for p in itertools.permutations(range(3)):
        standard = sum(p[i] == i for i in range(3)) - 1
        assert edge_trace(p) == standard + parity(p)
        checked += 1
    # Negative control: forgetting the orientation is already wrong at (12).
    p = (1, 0, 2, 3)
    unsigned = sum({p[i], p[j]} == {i, j} for i, j in itertools.combinations(range(4), 2))
    assert unsigned != edge_trace(p)
    return {"permutations_checked": checked, "orientation_negative_control": "passed"}


def explicit_gap_check():
    t = Q(15, 16)
    assert t - t**3 / 6 > Q(4, 5)
    assert 6 * 20**2 < 49**2
    assert 21 * 11**2 > 50**2
    upper = Q(64, 9) * (Q(13, 24) - Q(1240, 2695) + Q(451, 2750))
    assert 1 < upper < Q(7, 4)
    assert not upper < Q(3, 2)  # Detect an unsupported stronger bound.
    return {
        "H_half_rational_upper": str(upper),
        "H_half_upper_less_than_7_over_4": True,
        "scope": "Rational arithmetic in S21; uses analytic Gaussian integral and classical 3<pi<22/7 from the proof.",
        "gap": "ell*(2/sqrt(3)-1), sign sector at rho=1/2 only",
    }


def elementary_half_point_and_reflections():
    """The finite polynomial/rational arithmetic used in S32 and S36."""
    x = Q(93, 100)
    assert x - x**3/6 + x**5/120 - x**7/5040 > Q(4, 5)
    x = Q(101, 500)
    assert x - x**3/6 > Q(1, 5)
    x = Q(411, 1000)
    assert x - x**3/6 + x**5/120 < Q(2, 5)
    x = Q(253, 1000)
    assert x - x**3/6 > Q(1, 4)
    assert 16 * (Q(1, 5) - Q(1, 5)**3/3) - Q(4, 239) > Q(25, 8)
    assert 3 * Q(577, 1000)**2 < 1
    assert 15 * Q(258, 1000)**2 < 1
    assert 6 * Q(409, 1000)**2 > 1
    assert 21 * Q(219, 1000)**2 > 1
    R = Q(4, 3)*(1+2*Q(93, 100)/Q(25, 8))
    R += 2*Q(409, 1000)*(1+2*Q(101, 500)/Q(25, 8))
    R += Q(219, 1000)*(1-2*Q(411, 1000)/Q(22, 7))
    L = Q(4, 3)*Q(577, 1000) + Q(258, 1000)*(1-2*Q(253, 1000)/Q(25, 8))
    upper = Q(16, 9)*R - Q(16, 3)*L*L
    assert R == Q(2650234883, 825000000)
    assert L == Q(4619803, 4687500)
    assert upper == Q(32058163928093, 60424804687500)
    assert 0 < upper < Q(9, 16)
    # With u=w-1>=0, compare the exact two difference polynomials.
    plus = [Q(math.comb(4, j)) + (2 if j <= 1 else 0) - (3 if j == 0 else 0) for j in range(5)]
    minus = [Q(math.comb(4, j)) - (4 if j <= 1 else 0) + (3 if j == 0 else 0) for j in range(5)]
    assert plus == [0, 6, 6, 4, 1]
    assert minus == [0, 0, 6, 4, 1]
    assert all(x >= 0 for x in plus + minus)
    # Replacing 3/2 by 2 in K- >= c*D is false for every w>1.
    assert 4*(1-Q(2)) < 0
    return {
        "R_half_upper": str(R), "L_half_lower": str(L),
        "G_half_upper": str(upper), "G_half_less_than_9_over_16": True,
        "reflection_difference_coefficients_in_u": {"plus": [str(x) for x in plus], "minus": [str(x) for x in minus]},
        "scope": "Exact finite algebra in S32/S36. The all-correlation theorem additionally uses the written OU parity spectra, ground-state identities, interaction monotonicity and extra-dimensional Hermite decomposition; this script does not verify those analytic arguments.",
    }


def iv(x):
    return (Q(x[0]), Q(x[1])) if isinstance(x, tuple) else (Q(x), Q(x))


def add(x, y):
    x, y = iv(x), iv(y)
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    x = iv(x)
    return -x[1], -x[0]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    values = [a * b for a in iv(x) for b in iv(y)]
    return min(values), max(values)


def div(x, y):
    y = iv(y)
    if y[0] * y[1] <= 0:
        raise ValueError("interval denominator must exclude zero")
    return mul(x, (1 / y[1], 1 / y[0]))


def round_out(x, scale=10**10):
    lo, hi = iv(x)
    return Q((lo * scale).__floor__(), scale), Q((hi * scale).__ceil__(), scale)


def sqrt_point(x, scale=10**12):
    x = Q(x)
    lo = Q(math.isqrt(x.numerator * scale**2 // x.denominator), scale)
    hi = lo + Q(1, scale)
    assert lo * lo <= x < hi * hi
    return lo, hi


def sqrt_iv(x):
    lo, hi = iv(x)
    return sqrt_point(lo)[0], sqrt_point(hi)[1]


def atan_point(x, count=20):
    x = Q(x)
    assert 0 < x < 1
    value = sum(((-1)**n * x**(2*n+1) / (2*n+1) for n in range(count)), Q(0))
    next_term = (-1)**count * x**(2*count+1) / (2*count+1)
    return min(value, value + next_term), max(value, value + next_term)


def pi_interval():
    # tan(4 atan(1/5)-atan(1/239))=1, with the angle in (0,pi/2).
    # Thus Machin's identity gives pi=16 atan(1/5)-4 atan(1/239).
    return round_out(sub(mul(16, atan_point(Q(1, 5))), mul(4, atan_point(Q(1, 239)))))


def asin_point(x, count=64):
    x = Q(x)
    if x < 0:
        return neg(asin_point(-x, count))
    assert 0 <= x < 1
    def coefficient(n):
        return Q(math.comb(2*n, n), 4**n * (2*n+1))
    lower = sum((coefficient(n) * x**(2*n+1) for n in range(count)), Q(0))
    # Positive decreasing coefficients bound the omitted tail geometrically.
    tail = coefficient(count) * x**(2*count+1) / (1 - x*x)
    return round_out((lower, lower + tail))


def asin_iv(x):
    lo, hi = iv(x)
    return asin_point(lo)[0], asin_point(hi)[1]


def operator_bounds(lower, upper):
    """Rational interval extensions of the analytic functions H and G."""
    rho, pi = (Q(lower), Q(upper)), pi_interval()
    rho2 = mul(rho, rho)
    sigma2 = sub(1, rho2)
    def quadrant(x):
        return add(Q(1, 4), div(asin_iv(x), mul(2, pi)))
    def J(u):
        base = add(1, rho2)
        delta = sub(mul(base, base), mul(rho2, u*u))
        return div(quadrant(div(mul(rho, u), base)), sqrt_iv(delta))
    j2, jhalf, jneg = J(Q(2)), J(Q(1, 2)), J(Q(-1))
    H = mul(div(Q(16, 9), rho2), add(sub(j2, mul(2, jhalf)), jneg))
    R = add(add(mul(4, j2), mul(4, jhalf)), jneg)
    L = add(div(mul(2, quadrant(rho)), sqrt_iv(sigma2)),
            div(quadrant(mul(Q(-1, 2), rho)), sqrt_iv(sub(1, mul(Q(1, 4), rho2)))))
    G = sub(mul(div(Q(4, 9), rho2), R),
            mul(div(mul(Q(16, 9), sigma2), rho2), mul(L, L)))
    return H, G


def interval_gap_check():
    # Check interval primitives on exact cases, including their rejection path.
    assert mul((-2, -1), (3, 4)) == (-8, -3)
    assert div((2, 4), (-2, -1)) == (-4, -1)
    try:
        div(1, (-1, 1))
    except ValueError:
        pass
    else:
        raise AssertionError("invalid denominator accepted")
    half_H, half_G = operator_bounds(Q(1, 2), Q(1, 2))
    assert half_H[1] < Q(7, 4)
    assert half_G[1] < Q(9, 16)
    maxima = [Q(0), Q(0)]
    for j in range(250, 600):
        H, G = operator_bounds(Q(j, 1000), Q(j+1, 1000))
        assert H[1] < Q(19, 10), (j, "H")
        assert G[1] < Q(4, 5), (j, "G")
        maxima = [max(maxima[0], H[1]), max(maxima[1], G[1])]
    rounded = [round_out(x)[1] for x in maxima]
    return {
        "interval": ["1/4", "3/5"],
        "closed_cover_mesh": "1/1000",
        "closed_cover_intervals": 350,
        "H_upper": str(rounded[0]),
        "G_upper": str(rounded[1]),
        "targets": {"H_less_than": "19/10", "G_less_than": "4/5"},
        "half_correlation_G_upper": str(round_out(half_G)[1]),
        "arithmetic": "Python Fraction, integer-square-root enclosures; positive arcsin series with geometric tail, alternating atan series and Machin pi identity; outward rational rounding.",
        "scope": "Uniform H/G bounds on every point of the stated interval, via the analytic operator reduction S24-S26. No sampling argument.",
        "full_three_cell_gap": "ell*(1-3/sqrt(10)) on [1/4,3/5]",
        "full_three_cell_gap_at_half": "ell*(2/sqrt(3)-1)",
    }


def quadrature():
    """Test S23 at four chosen correlations, without error certification."""
    import numpy as np
    phi0 = 1 / math.sqrt(2 * math.pi)
    cdf = np.vectorize(lambda x: (1 + math.erf(x / math.sqrt(2))) / 2)
    output = []
    for count in (180, 360):
        grid, weights = np.polynomial.legendre.leggauss(count)
        r, dr = 5 * (grid + 1), 5 * weights  # Finite [0,10] cutoff.
        w = dr * np.exp(-r * r / 2) / (2 * math.pi)
        rr, ss = r[:, None], r[None, :]
        for rho in (0.1, 0.5, 0.9, 0.99):
            sigma = math.sqrt(1 - rho * rho)
            u = rho * r / sigma
            a = rho / sigma * (2 * phi0 * cdf(u) + phi0 * np.exp(-3 * u * u / 8) * cdf(-u / 2))
            # Include quadrature weights before exponentiation to avoid overflow.
            base = (np.log(w)[:, None] + np.log(w)[None, :]) / 2
            base -= rho * rho * (rr * rr + ss * ss) / (2 * sigma * sigma)
            same = np.exp(base + rho * rr * ss / (sigma * sigma)) / (sigma * sigma)
            other = np.exp(base - rho * rr * ss / (2 * sigma * sigma)) / (sigma * sigma)
            scaled = (2 * same + other) / np.sqrt(a[:, None] * a[None, :])
            constraint = np.sqrt(w / a)
            # Euclidean coordinates are sqrt(a*w)*h, so the constraint is nu-mean.
            basis, _ = np.linalg.qr(constraint[:, None], mode="complete")
            reduced = basis[:, 1:].T @ scaled @ basis[:, 1:]
            maximum = float(np.linalg.eigvalsh(reduced)[-1])
            standard_basis, _ = np.linalg.qr(np.sqrt(w)[:, None], mode="complete")
            sign_basis, _ = np.linalg.qr((np.sqrt(w)*r)[:, None], mode="complete")
            def restricted_min(matrix, basis):
                return float(np.linalg.eigvalsh(basis[:, 1:].T @ matrix @ basis[:, 1:]) [0])
            standard_min = restricted_min(np.diag(a/rho)-(2*same+other)/rho, standard_basis)
            sign_min = restricted_min(np.diag(a/rho)-2*(same-other)/rho, sign_basis)
            mixed_rotation_min = restricted_min(np.diag(a/rho)-(2*same+other), standard_basis)
            ell = 1.5*phi0
            gap = ell*(1/sigma-math.sqrt(3)/2) if rho <= .5 else ell*(1-rho)*(3*rho-1)/(3*rho*sigma)
            assert standard_min >= gap-1e-8
            assert sign_min >= ell*sigma-1e-8
            assert mixed_rotation_min >= ell*sigma-1e-8
            output.append({"nodes": count, "cutoff": 10, "rho": rho,
                           "largest_constrained_generalized_eigenvalue": maximum,
                           "standard_L2_minimum": standard_min, "sign_L2_minimum_off_rotations": sign_min,
                           "extra_degree_one_standard_minimum_off_rotations": mixed_rotation_min,
                           "analytic_full_gap_lower": gap, "analytic_sign_gap_lower": ell*sigma})
    return {"numpy_version": np.__version__, "scope": "Uncertified Nystrom exploration on [0,10]; no tail, discretization or rounding error enclosure; not an all-rho or fixed-rho proof.", "cases": output}


def perimeter_controls():
    import numpy as np
    q = NormalDist().inv_cdf(2/3)
    beta = 3/(4*math.pi**2)-1/(4*math.pi*math.sqrt(3))
    c = (2*math.exp(-q*q/2)-1.5)/math.pi
    x, w = np.polynomial.legendre.leggauss(40)
    cases = []
    for theta in (.1, .05, .025):
        u = theta*(x+1)/2
        Cq = float(theta/2*np.sum(w*np.exp(-q*q/(1+np.cos(u))))/math.pi)
        NT = 2/3-3/(4*math.pi**2)*((math.pi-theta)**2-math.acos(math.cos(theta)/2)**2)
        J_upper = math.exp(-q*q/(1-math.cos(theta)))
        normalized_without_jump = (2*Cq-NT)/theta
        cases.append({"theta": theta, "simplex_beta_quotient": (3*theta/(2*math.pi)-NT)/theta**2,
                      "slab_normalized_deficit_without_jump": normalized_without_jump,
                      "slab_normalized_limit": c,
                      "jump_correction_upper": 2*J_upper/theta})
        assert normalized_without_jump-2*J_upper/theta > c
    for rho in (.1, .5, .9):
        def I(t):
            return (.25+math.asin(rho*t)/(2*math.pi))/(2*math.pi*math.sqrt(1-rho*rho*t*t))
        derivative = 3*math.acos(-rho)/(2*math.pi**2*math.sqrt(1-rho*rho))
        derivative += 3*math.acos(rho/2)/(4*math.pi**2*math.sqrt(1-rho*rho/4))
        assert abs(derivative-(6*I(1)+3*I(-.5))) < 1e-12
    return {"scope": "Uncertified floating normalization/expansion controls, not a proof of D0-D9 or monotonicity failure.",
            "q": q, "beta3": beta, "cases": cases}


def main():
    if not __debug__:
        raise RuntimeError("Run without -O: the exact certificate uses assertions.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quadrature", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    receipt = {
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_source_sha256": hashlib.sha256((here / "SIMPLEX_INTERFACE_SYMMETRY.md").read_bytes()).hexdigest(),
        "perimeter_proof_source_sha256": hashlib.sha256((here / "THREE_CELL_DEFICIT_AND_PERIMETER.md").read_bytes()).hexdigest(),
        "exact_checks": {"characters": character_checks(), "half_correlation_bound": explicit_gap_check(),
                         "elementary_half_point_and_reflections": elementary_half_point_and_reflections(), "interval_gap": interval_gap_check()},
        "analytic_claims_certified_by_this_script": False,
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    if args.quadrature:
        receipt["standard_sector_exploration"] = quadrature()
        receipt["perimeter_controls"] = perimeter_controls()
    rendered = json.dumps(receipt, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
