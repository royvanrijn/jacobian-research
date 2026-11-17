"""Exact finite certificates for binary pure vanishing over QQ.

Completeness uses the Hall/shifted-ray theorem in BINARY_GVC_FINITE_CERTIFICATE.md.
Positive weight certificates themselves use only elementary support separation.
No moment prefix is used to decide an input. Coefficients and coordinate changes
are rational; the unique Hall direction descends to QQ when degree(P) >= r.
"""
from __future__ import annotations

from fractions import Fraction
from math import comb, gcd
from typing import Mapping

import sympy as sp

Exponent = tuple[int, int]
Polynomial = dict[Exponent, Fraction]


def exact_polynomial(terms: Mapping[Exponent, object]) -> Polynomial:
    result = {}
    for exponent, coefficient in terms.items():
        if len(exponent) != 2 or any(type(i) is not int or i < 0 for i in exponent):
            raise ValueError("exponents must be pairs of nonnegative integers")
        if isinstance(coefficient, float):
            raise ValueError("floating-point coefficients are not exact inputs")
        value = Fraction(coefficient)
        if value:
            result[tuple(exponent)] = value
    return result


def encode(terms: Polynomial) -> list[list[object]]:
    return [[i, j, str(c)] for (i, j), c in sorted(terms.items())]


def decode(rows: list[list[object]]) -> Polynomial:
    result = {}
    for i, j, coefficient in rows:
        exponent = (i, j)
        if exponent in result:
            raise ValueError("duplicate input monomial")
        result[exponent] = coefficient
    return exact_polynomial(result)


def degree(terms: Polynomial) -> int:
    if not terms:
        raise ValueError("the zero polynomial has no degree in this interface")
    return max(sum(exponent) for exponent in terms)


def transform(terms: Polynomial, chart: dict, *, symbol: bool) -> Polynomial:
    """Apply P(M z) or lambda(M^(-T) X), with the same matrix M.

    Finite root a: M = [[1,0],[-a,1]], so the substitutions are
    P(x,y-a*x) and lambda(X+a*Y,Y). The infinity chart swaps coordinates.
    """
    if chart["kind"] == "swap":
        return {(j, i): c for (i, j), c in terms.items()}
    a = Fraction(chart["root"])
    result: Polynomial = {}
    for (i, j), coefficient in terms.items():
        length = i if symbol else j
        for k in range(length + 1):
            if symbol:
                exponent = (k, i + j - k)
                value = coefficient * comb(i, k) * a ** (i - k)
            else:
                exponent = (i + j - k, k)
                value = coefficient * comb(j, k) * (-a) ** (j - k)
            result[exponent] = result.get(exponent, Fraction(0)) + value
    return {exponent: c for exponent, c in result.items() if c}


def rational_directions(symbol: Polynomial, r: int) -> tuple[list[dict], dict]:
    """Factor the lowest symbol and retain its rational derivative directions.

    Nonlinear irreducible factors are recorded, not silently treated as linear
    directions. A Hall direction would be unique, so cannot belong to them.
    """
    t = sp.Symbol("t")
    top = {exponent: c for exponent, c in symbol.items() if sum(exponent) == r}
    univariate = sp.Poly(sum(sp.Rational(c.numerator, c.denominator) * t**i
                             for (i, _j), c in top.items()), t, domain=sp.QQ)
    content, factors = sp.factor_list(univariate)
    reconstructed = sp.Poly(content, t, domain=sp.QQ)
    directions = []
    factor_rows = []
    for factor, multiplicity in factors:
        reconstructed *= factor**multiplicity
        factor_rows.append({
            "coefficients_descending": [str(c) for c in factor.all_coeffs()],
            "multiplicity": multiplicity,
            "degree": factor.degree(),
        })
        if factor.degree() == 1:
            root = -Fraction(factor.nth(0)) / Fraction(factor.nth(1))
            directions.append({"kind": "finite", "root": str(root),
                               "multiplicity": multiplicity})
    assert reconstructed == univariate
    directions.sort(key=lambda row: Fraction(row["root"]))
    infinity_multiplicity = r - univariate.degree()
    if infinity_multiplicity:
        directions.append({"kind": "swap", "multiplicity": infinity_multiplicity})
    return directions, {"content": str(content), "factors": factor_rows,
                        "infinity_multiplicity": infinity_multiplicity}


def strict_slope(symbol: Polynomial, polynomial: Polynomial) -> dict:
    """Solve every strict inequality a*s+b>0, together with s>1.

    Return exact infeasibility data or a small integer weight. The mediant
    bound is proved separately, rather than inferred from bounded searches.
    """
    lower = Fraction(1)
    upper = None
    lower_pair = None
    upper_pair = None
    for alpha in sorted(symbol):
        for beta in sorted(polynomial):
            a, b = alpha[0] - beta[0], alpha[1] - beta[1]
            pair = {"operator": list(alpha), "polynomial": list(beta)}
            if a == 0:
                if b <= 0:
                    return {"feasible": False, "reason": "constant inequality",
                            "pair": pair, "constant": b}
                continue
            endpoint = Fraction(-b, a)
            if a > 0 and endpoint > lower:
                lower, lower_pair = endpoint, pair
            if a < 0 and (upper is None or endpoint < upper):
                upper, upper_pair = endpoint, pair
    interval = {"lower": str(lower), "upper": None if upper is None else str(upper),
                "lower_pair": lower_pair, "upper_pair": upper_pair}
    if upper is not None and lower >= upper:
        return {"feasible": False, "reason": "empty open interval", **interval}
    if upper is None:
        u, v = lower.numerator + lower.denominator, lower.denominator
    else:
        u = lower.numerator + upper.numerator
        v = lower.denominator + upper.denominator
    divisor = gcd(u, v)
    u, v = u // divisor, v // divisor
    assert u > v > 0
    operator_min = min(u*i + v*j for i, j in symbol)
    polynomial_max = max(u*i + v*j for i, j in polynomial)
    gap = operator_min - polynomial_max
    assert gap > 0
    assert max(u, v) <= degree(symbol) + degree(polynomial)
    return {"feasible": True, **interval, "weight": [u, v],
            "operator_min": operator_min, "polynomial_max": polynomial_max,
            "gap": gap}


def classify(symbol: Mapping[Exponent, object], polynomial: Mapping[Exponent, object],
             multiplier: Mapping[Exponent, object] | None = None) -> dict:
    """Decide all-positive-power pure vanishing for a rational binary pair.

    A negative result rejects the *premise*, not GVC. It uses the stated
    mathematical completeness theorem. Every positive nontrivial result
    includes a directly verifiable weight certificate.
    """
    symbol = exact_polynomial(symbol)
    polynomial = exact_polynomial(polynomial)
    multiplier = None if multiplier is None else exact_polynomial(multiplier)
    result = {"format": "binary-gvc-finite-certificate-v1", "field": "QQ",
              "input_symbol": encode(symbol), "input_polynomial": encode(polynomial),
              "input_multiplier": None if multiplier is None else encode(multiplier),
              "completeness_source": "GVC2SC"}
    if not symbol or not polynomial:
        result.update(pure_status="all_powers_zero", method="zero input")
        if multiplier is not None:
            result["mixed_cutoff"] = 1
        return result
    if (0, 0) in symbol:
        return {**result, "pure_status": "pure_premise_false",
                "method": "nonzero constant symbol", "first_nonzero_pure_power": 1}

    r = min(sum(exponent) for exponent in symbol)
    d, R = degree(polynomial), degree(symbol)
    result.update(lowest_operator_order=r, highest_operator_order=R, polynomial_degree=d)
    if d < r:
        transformed_symbol, transformed_polynomial = symbol, polynomial
        chart = {"kind": "identity"}
        separation = {"feasible": True, "weight": [1, 1], "operator_min": r,
                      "polynomial_max": d, "gap": r - d}
        result["method"] = "ordinary degree gap"
    else:
        directions, factorization = rational_directions(symbol, r)
        result["lowest_symbol_factorization"] = factorization
        trials, candidates = [], []
        for chart in directions:
            transformed_symbol = transform(symbol, chart, symbol=True)
            transformed_polynomial = transform(polynomial, chart, symbol=False)
            e = min(i for i, j in transformed_symbol if i + j == r)
            t = max(i for i, j in transformed_polynomial if i + j == d)
            assert e == chart["multiplicity"]
            trials.append({"chart": chart, "operator_x_multiplicity": e,
                           "top_polynomial_max_x": t, "hall_candidate": t < e})
            if t < e:
                candidates.append((chart, transformed_symbol, transformed_polynomial))
        assert len(candidates) <= 1, "distinct Hall candidates contradict multiplicity counting"
        result["direction_trials"] = trials
        if not candidates:
            return {**result, "pure_status": "pure_premise_false",
                    "method": "no Hall direction; uniqueness rules out nonrational directions"}
        chart, transformed_symbol, transformed_polynomial = candidates[0]
        separation = strict_slope(transformed_symbol, transformed_polynomial)
        result.update(chart=chart, transformed_symbol=encode(transformed_symbol),
                      transformed_polynomial=encode(transformed_polynomial),
                      separation=separation)
        if not separation["feasible"]:
            return {**result, "pure_status": "pure_premise_false",
                    "method": "unique Hall chart has no strict positive separator"}
        result["method"] = "unique rational Hall chart and strict support separation"

    result.update(pure_status="all_powers_zero", chart=chart,
                  transformed_symbol=encode(transformed_symbol),
                  transformed_polynomial=encode(transformed_polynomial), separation=separation)
    if multiplier is not None:
        transformed_q = multiplier if chart["kind"] == "identity" else transform(
            multiplier, chart, symbol=False)
        u, v = separation["weight"]
        q_weight = max((u*i + v*j for i, j in transformed_q), default=0)
        cutoff = q_weight // separation["gap"] + 1
        q_degree = degree(multiplier) if multiplier else 0
        bound = (R + d) * q_degree + 1
        assert cutoff <= bound
        result.update(transformed_multiplier=encode(transformed_q), multiplier_weight=q_weight,
                      mixed_cutoff=cutoff, degree_only_cutoff=bound)
    return result
