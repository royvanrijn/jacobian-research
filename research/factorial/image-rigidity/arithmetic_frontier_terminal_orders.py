"""Exact certificates for the no-carry, single-digit terminal-order theorem.

Written proof: ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex.
A finite replay checks implementation, not the unbounded theorem.
"""
from fractions import Fraction
from math import comb, factorial, gcd, prod
from pathlib import Path
import hashlib
import json
import platform
import sympy as s

from prime_power_windows import (expected_edges, factorial_valuation, newton_residuals,
                                 residue, valuation)

X = s.Symbol("x")


def exact_integer(value, minimum):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError("exact integer >= " + str(minimum) + " required")
    return value


def parameters(p, r, exponent, m, offset):
    p = exact_integer(p, 2)
    if not s.isprime(p):
        raise ValueError("prime p required")
    r = exact_integer(r, 0)
    exponent = exact_integer(exponent, 1)
    m = exact_integer(m, 1)
    offset = exact_integer(offset, 0)
    if m >= p:
        raise ValueError("1 <= m < p required")
    place = p ** r
    gap = p ** exponent
    digit = offset // p ** (r + exponent) % p
    if digit + m >= p:
        raise ValueError("unsupported carry at place r+s")
    return p, place, gap, digit


def class_certificate(p, r, exponent, m, offset):
    p, place, gap, digit = parameters(p, r, exponent, m, offset)
    height_step = (gap * place - 1) // (p - 1)
    epsilon = (-1) ** height_step % p
    residual = [comb(m, j) * prod(range(digit + 1, digit + j + 1))
                * pow(epsilon, j, p) % p for j in range(m + 1)]
    return {
        "prime": p, "r": r, "s": exponent, "m": m, "C": offset,
        "terminal_N": m * place, "gap": gap, "offset_digit_t": digit,
        "edge": {"left": [0, 0], "right": [m * place, m * height_step],
                 "e": place, "h": height_step,
                 "residual_coefficients": residual, "residual_squarefree": True},
        "squarefree": True, "irreducible": m == 1,
        "proof": "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
        "all_terminal_orders": False, "unrestricted_binomial_SFC": False,
        "finite_extrapolation": False, "external_review": False,
        "independent_replay": False, "formal_verification": False,
    }


def coefficients(terminal, offset, gap):
    return [comb(terminal, k) * factorial(offset + gap * k) // factorial(offset)
            for k in range(terminal + 1)]


def polynomial(terminal, offset, gap):
    return s.Poly.from_list(list(reversed(coefficients(terminal, offset, gap))), X)


def literal_certificate(p, r, exponent, m, offset):
    result = class_certificate(p, r, exponent, m, offset)
    actual = newton_residuals(coefficients(result["terminal_N"], offset, result["gap"]), p)
    assert actual == [result["edge"]], (result, actual)
    result["actual_literal_hull_and_residual_match"] = True
    residual = s.Poly.from_list(list(reversed(result["edge"]["residual_coefficients"])),
                               X, modulus=p)
    epsilon = (-1) ** result["edge"]["h"] % p
    digit = result["offset_digit_t"]
    expression = (epsilon * X ** 2 * residual.diff().diff().as_expr()
                  + (1 - epsilon * (m - digit - 2) * X) * residual.diff().as_expr()
                  - epsilon * m * (digit + 1) * residual.as_expr())
    assert s.Poly(expression, X, modulus=p).is_zero
    assert s.gcd(residual, residual.diff()).degree() == 0
    assert residual.degree() == m
    return result


def lower_one_certificate(p, r, exponent, m, beta=1):
    beta = exact_integer(beta, 1)
    terminal = m * p ** r
    if terminal < 2:
        raise ValueError("terminal N >= 2 required for two positive orders")
    result = class_certificate(p, r, exponent, m, terminal + beta - 2)
    result["consecutive_pair"] = {
        "orders": [terminal - 1, terminal], "lower_exponent": 1,
        "integer_Gamma_shape": beta, "all_complex_coefficients": True,
        "mechanism": "squarefree tilted polynomial and exact lower-one bridge"}
    return result


def general_lower_certificate(p, r, exponent, lower):
    lower = exact_integer(lower, 0)
    if r < 1:
        raise ValueError("r >= 1 required for terminal N >= 2")
    terminal = p ** r
    result = class_certificate(p, r, exponent, 1, lower * terminal)
    result["consecutive_pair"] = {
        "orders": [terminal - 1, terminal], "lower_exponent": lower,
        "integer_Gamma_shape": 1, "all_complex_coefficients": True,
        "mechanism": "irreducibility of degree N and lower degree N-1"}
    return result



def factorial_unit_residue(number, prime):
    """Unit part of number!, computed from digits without constructing number!."""
    valuation = 0
    remaining = number
    digit_product = 1
    while remaining:
        remaining, digit = divmod(remaining, prime)
        digit_product = digit_product * factorial(digit) % prime
        valuation += remaining
    return pow(-1, valuation, prime) * digit_product % prime


def coefficient_unit(terminal, offset, gap, index, prime):
    """Unit part of binom(N,k)*(C+dk)!/C!, without large factorials."""
    numerator = (factorial_unit_residue(terminal, prime)
                 * factorial_unit_residue(offset + gap * index, prime))
    denominator = (factorial_unit_residue(index, prime)
                   * factorial_unit_residue(terminal - index, prime)
                   * factorial_unit_residue(offset, prime))
    return numerator * pow(denominator, -1, prime) % prime


def minus_one_class_certificate(p, r, exponent, m=1):
    p = exact_integer(p, 2)
    if not s.isprime(p):
        raise ValueError("prime p required")
    r = exact_integer(r, 1)
    exponent = exact_integer(exponent, 1)
    m = exact_integer(m, 1)
    if m >= p:
        raise ValueError("1 <= m < p required")
    terminal, gap = m * p ** r - 1, p ** exponent
    edges = []
    if r < exponent:
        base_slope = (gap - 1) // (p - 1)
        for row in expected_edges(terminal, p, (-1) ** base_slope % p):
            result = dict(row)
            result["left"] = [row["left"][0], row["left"][1]
                              + base_slope * row["left"][0]]
            result["right"] = [row["right"][0], row["right"][1]
                               + base_slope * row["right"][0]]
            result["h"] = row["h"] + base_slope * row["e"]
            edges.append(result)
        mechanism = "existing fixed-shape prime-power theorem, Gamma shape N<gap"
    else:
        q = p ** (r - exponent)

        def falling_edge(left, right, step, residual_degree, left_height):
            height_step = (gap * step - 1) // (p - 1)
            epsilon = (-1) ** height_step % p
            falling = 1
            residual = [1]
            for j in range(1, residual_degree + 1):
                falling *= residual_degree - j + 1
                residual.append(falling * pow(epsilon, j, p) % p)
            return {"left": [left, left_height],
                    "right": [right, left_height + residual_degree * height_step],
                    "e": step, "h": height_step,
                    "residual_coefficients": residual, "residual_squarefree": True}

        first_right = (p - m) * q
        edges.append(falling_edge(0, first_right, q, p - m, 0))
        for i in range(2, exponent + 1):
            step = p ** (i - 1) * q
            left = (p ** (i - 1) - m) * q
            right = (p ** i - m) * q
            left_height = (gap * left + m - 1) // (p - 1) - (i - 1)
            edges.append(falling_edge(left, right, step, p - 1, left_height))
        last_middle = (gap - m) * q
        final_left = m * (gap - 1) * q
        if m > 1:
            left_height = (gap * last_middle + m - 1) // (p - 1) - exponent
            edges.append(falling_edge(last_middle, final_left, gap * q, m - 1,
                                      left_height))
        width = m * q - 1
        if width:
            left_height = gap * final_left // (p - 1) - exponent
            height = (gap * width - (m - 1)) // (p - 1) + exponent
            divisor = gcd(width, height)
            ratio = (coefficient_unit(terminal, terminal - 1, gap, terminal, p)
                     * pow(coefficient_unit(terminal, terminal - 1, gap, final_left, p),
                           -1, p)) % p
            residual = [1] + [0] * (divisor - 1) + [ratio]
            edges.append({
                "left": [final_left, left_height],
                "right": [terminal, (gap * terminal - (m - 1)) // (p - 1)],
                "e": width // divisor, "h": height // divisor,
                "residual_coefficients": residual, "residual_squarefree": True})
        mechanism = "digit-sum polygon, reversed truncated exponentials, final binomial"
    return {
        "prime": p, "r": r, "s": exponent, "m": m, "terminal_N": terminal,
        "gap": gap, "Gamma_shape_of_tilted_polynomial": terminal,
        "edges": edges, "squarefree": True, "mechanism": mechanism,
        "positive_consecutive_orders": [terminal - 1, terminal] if terminal >= 2 else None,
        "proof": "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
        "all_terminal_orders": False, "unrestricted_binomial_SFC": False,
        "finite_extrapolation": False}


def minus_one_literal_certificate(p, r, exponent, m=1):
    result = minus_one_class_certificate(p, r, exponent, m)
    terminal, gap = result["terminal_N"], result["gap"]
    actual = newton_residuals(coefficients(terminal, terminal - 1, gap), p)
    assert actual == result["edges"], (result, actual)
    return result


def gap_three_two_minus_one_class_certificate(r):
    return minus_one_class_certificate(3, r, 1, 2)


def gap_three_two_minus_one_literal_certificate(r):
    return minus_one_literal_certificate(3, r, 1, 2)





def native_prime_certificate(lower, gap, terminal, prime):
    lower = exact_integer(lower, 0)
    gap = exact_integer(gap, 1)
    terminal = exact_integer(terminal, 2)
    prime = exact_integer(prime, 2)
    if not s.isprime(prime) or prime <= gap * terminal:
        raise ValueError("prime p>dN required")
    index = (-lower * terminal) % prime
    if not 1 <= index <= gap:
        raise ValueError("p must divide a member of the initial length-d block")
    if valuation(lower * terminal + index, prime) != 1:
        raise ValueError("selected initial block member must have valuation exactly one")
    return {
        "lower_exponent": lower, "gap": gap, "terminal_N": terminal,
        "prime": prime, "initial_block_index": index,
        "normalized_constant": 1, "all_nonconstant_coefficient_valuations": 1,
        "reciprocal_is_Eisenstein": True, "irreducible_over_Qp": True,
        "consecutive_orders": [terminal - 1, terminal], "coprime": True,
        "arbitrary_integer_gap_allowed": True,
        "proof": "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
        "all_terminal_orders": False, "unrestricted_binomial_SFC": False,
        "finite_extrapolation": False}


def native_prime_literal_certificate(lower, gap, terminal, prime):
    result = native_prime_certificate(lower, gap, terminal, prime)
    normalized = coefficients(terminal, lower * terminal, gap)
    assert normalized[0] == 1
    assert all(valuation(c, prime) == 1 for c in normalized[1:])
    reciprocal = s.Poly.from_list(normalized, X)
    assert reciprocal.LC() == 1
    assert valuation(int(reciprocal.TC()), prime) == 1
    assert reciprocal.is_irreducible
    preceding = polynomial(terminal - 1, lower * (terminal - 1), gap)
    current = polynomial(terminal, lower * terminal, gap)
    assert s.gcd(current, preceding).degree() == 0
    return result



def forward_native_prime_certificate(lower, gap, central_order, remainder):
    lower = exact_integer(lower, 2)
    gap = exact_integer(gap, 1)
    central_order = exact_integer(central_order, 3)
    remainder = exact_integer(remainder, 1)
    if lower < gap or remainder > gap:
        raise ValueError("a>=d and 1<=j<=d required")
    prime = lower * central_order + remainder
    if not s.isprime(prime):
        raise ValueError("aN+j must itself be prime")
    native = native_prime_certificate(lower, gap, central_order, prime)
    next_offset = lower * (central_order + 1)
    next_end = (lower + gap) * (central_order + 1)
    assert next_offset >= prime and next_end < 3 * prime
    assert prime > central_order + 1 and prime > 3
    if 2 * prime <= next_end:
        first_index = (2 * prime - next_offset + gap - 1) // gap
        assert central_order <= first_index <= central_order + 1
        nonzero_width = central_order + 2 - first_index
        root_valuation_profile = ["0", "-1/" + str(nonzero_width)]
    else:
        first_index = None
        nonzero_width = None
        root_valuation_profile = ["0"]
    native.update({
        "forward_pair": [central_order, central_order + 1],
        "backward_and_forward_pairs_coprime": True,
        "native_root_valuation": "-1/" + str(central_order),
        "next_moment_possible_root_valuations": root_valuation_profile,
        "next_first_divisible_coefficient_index": first_index,
        "next_nonzero_edge_width": nonzero_width})
    return native


def forward_native_prime_literal_certificate(lower, gap, central_order, remainder):
    result = forward_native_prime_certificate(lower, gap, central_order, remainder)
    prime = result["prime"]
    current_coefficients = coefficients(central_order, lower * central_order, gap)
    next_coefficients = coefficients(central_order + 1, lower * (central_order + 1), gap)
    current_edges = newton_residuals(current_coefficients, prime)
    assert len(current_edges) == 1
    assert current_edges[0]["left"] == [0, 0]
    assert current_edges[0]["right"] == [central_order, 1]
    next_edges = newton_residuals(next_coefficients, prime)
    first_index = result["next_first_divisible_coefficient_index"]
    if first_index is None:
        assert all(valuation(c, prime) == 0 for c in next_coefficients)
        assert len(next_edges) == 1
        assert next_edges[0]["left"] == [0, 0]
        assert next_edges[0]["right"] == [central_order + 1, 0]
    else:
        assert [valuation(c, prime) for c in next_coefficients] == (
            [0] * first_index + [1] * (central_order + 2 - first_index))
        assert len(next_edges) == 2
        assert next_edges[0]["right"] == [first_index - 1, 0]
        assert next_edges[1]["right"] == [central_order + 1, 1]
        assert next_edges[1]["e"] == result["next_nonzero_edge_width"]
        assert next_edges[1]["h"] == 1
    current = s.Poly.from_list(list(reversed(current_coefficients)), X)
    following = s.Poly.from_list(list(reversed(next_coefficients)), X)
    assert s.gcd(current, following).degree() == 0
    return result


def split_literal_certificate(p, exponent, terminal, offset):
    p = exact_integer(p, 2)
    if not s.isprime(p):
        raise ValueError("prime p required")
    exponent = exact_integer(exponent, 1)
    terminal = exact_integer(terminal, 1)
    offset = exact_integer(offset, 0)
    gap = p ** exponent
    base_slope = (gap - 1) // (p - 1)
    gamma = offset // gap + 1
    actual = [Fraction(c, p ** (base_slope * k))
              for k, c in enumerate(coefficients(terminal, offset, gap))]
    base = [comb(terminal, k) * prod(range(gamma, gamma + k))
            * (-1) ** (base_slope * k) for k in range(terminal + 1)]
    for actual_coefficient, base_coefficient in zip(actual, base):
        ratio = actual_coefficient / base_coefficient
        assert valuation(ratio, p) == 0 and residue(ratio, p) == 1
    actual_edges = newton_residuals(actual, p)
    base_edges = newton_residuals(base, p)
    assert actual_edges == base_edges
    return {"prime": p, "s": exponent, "terminal_N": terminal, "C": offset,
            "first_newton_hull_and_residuals_match": True,
            "all_residuals_squarefree": all(row["residual_squarefree"]
                                          for row in actual_edges)}


def common_quotient_certificate(p, exponent, lower, terminal):
    p = exact_integer(p, 2)
    if not s.isprime(p):
        raise ValueError("prime p required")
    exponent = exact_integer(exponent, 1)
    lower = exact_integer(lower, 0)
    terminal = exact_integer(terminal, 2)
    gap = p ** exponent
    upper_q = lower * terminal // gap
    lower_q = lower * (terminal - 1) // gap
    if upper_q != lower_q or upper_q + terminal >= p:
        raise ValueError("same quotient q and q+N<p required")
    base_slope = (gap - 1) // (p - 1)
    resultant_valuation = (
        (terminal - 1) * factorial_valuation(lower * terminal, p)
        + terminal * factorial_valuation(lower * (terminal - 1), p)
        + base_slope * terminal * (terminal - 1))
    return {
        "prime": p, "s": exponent, "lower_exponent": lower, "gap": gap,
        "terminal_N": terminal, "quotient_q": upper_q,
        "consecutive_orders": [terminal - 1, terminal],
        "coprime": True, "resultant_valuation": resultant_valuation,
        "all_complex_phase_coefficients": True,
        "proof": "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
        "all_terminal_orders": False, "unrestricted_binomial_SFC": False,
        "finite_extrapolation": False}


def common_quotient_literal_certificate(p, exponent, lower, terminal):
    result = common_quotient_certificate(p, exponent, lower, terminal)
    gap = result["gap"]
    base_slope = (gap - 1) // (p - 1)
    epsilon = (-1) ** base_slope
    reduced = []
    for degree in (terminal, terminal - 1):
        normalized = coefficients(degree, lower * degree, gap)
        scaled = [Fraction(c * epsilon ** k, p ** (base_slope * k))
                  for k, c in enumerate(normalized)]
        expected = [comb(degree, k) * prod(range(result["quotient_q"] + 1,
                                                result["quotient_q"] + k + 1)) % p
                    for k in range(degree + 1)]
        actual = [residue(c, p) for c in scaled]
        assert actual == expected
        assert actual[-1] != 0
        reduced.append(s.Poly.from_list(list(reversed(actual)), X, modulus=p))
    assert int(s.resultant(*reduced)) % p
    high = polynomial(terminal, lower * terminal, gap)
    low = polynomial(terminal - 1, lower * (terminal - 1), gap)
    normalized_resultant = int(s.resultant(high, low))
    assert valuation(normalized_resultant, p) == base_slope * terminal * (terminal - 1)
    literal_raw_valuation = (valuation(normalized_resultant, p)
        + (terminal - 1) * factorial_valuation(lower * terminal, p)
        + terminal * factorial_valuation(lower * (terminal - 1), p))
    assert literal_raw_valuation == result["resultant_valuation"]
    return result


def gap_one_discriminant_formula(terminal, gamma):
    terminal = exact_integer(terminal, 1)
    gamma = exact_integer(gamma, 1)
    return ((-1) ** (terminal * (terminal - 1) // 2)
            * prod(j ** j for j in range(1, terminal + 1))
            * prod(prod(range(gamma, gamma + i)) for i in range(1, terminal)))


def global_precision_obstruction_certificate(h):
    h = exact_integer(h, 0)
    if (h + 1) % 3 == 0:
        raise ValueError("this exact precision-three family requires 3 not dividing h+1")
    terminal, gamma = 9 * h + 4, 3 * h + 2
    disc_valuation = sum(
        j * valuation(j, 3) + (j - 1) * valuation(gamma + terminal - j, 3)
        for j in range(1, terminal + 1))
    assert disc_valuation >= terminal + 1
    return {
        "h": h, "terminal_N": terminal, "gap": 3, "base_Gamma_shape": gamma,
        "best_unit_variable_scaling_coefficient_precision": 3,
        "base_discriminant_valuation": disc_valuation,
        "global_precision_test_can_certify": False,
        "actual_polynomial_has_repeated_root": None,
        "interpretation": "failure of a sufficient method, not an SFC counterexample",
        "proof": "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
        "finite_extrapolation": False}


def global_precision_obstruction_literal_certificate(h):
    result = global_precision_obstruction_certificate(h)
    terminal, gamma = result["terminal_N"], result["base_Gamma_shape"]
    q = 3 * h + 1
    preferred_unit = (3 * q + 1) * (3 * q + 2)
    actual = coefficients(terminal, terminal - 1, 3)
    base = [comb(terminal, k) * prod(range(gamma, gamma + k))
            for k in range(terminal + 1)]
    precisions = []
    for unit in [preferred_unit] + [preferred_unit * (1 + 3 ** j) for j in range(5)]:
        difference = [Fraction(c, (3 * unit) ** k) - base[k]
                      for k, c in enumerate(actual)]
        precision = min(valuation(c, 3) for c in difference if c)
        precisions.append(precision)
        assert precision <= 3
    assert precisions[0] == 3
    if terminal <= 13:
        reversed_base = s.Poly.from_list(base, X)
        actual_discriminant = int(s.discriminant(reversed_base))
        assert actual_discriminant == gap_one_discriminant_formula(terminal, gamma)
        assert valuation(actual_discriminant, 3) == result["base_discriminant_valuation"]
    return result


def replay():
    literal_count = bridge_count = pair_count = rejected_count = 0
    for p in (2, 3, 5, 7):
        for r in range(3):
            for exponent in (1, 2):
                for m in range(1, p):
                    terminal = m * p ** r
                    if terminal > 75:
                        continue
                    max_digit = p - m - 1
                    possible_offsets = {0, terminal - 1,
                                        max_digit * p ** (r + exponent)
                                        + p ** r - 1}
                    for offset in sorted(possible_offsets):
                        if offset > 500:
                            continue
                        literal_certificate(p, r, exponent, m, offset)
                        literal_count += 1
                    if terminal < 2:
                        continue
                    result = lower_one_certificate(p, r, exponent, m)
                    gap = result["gap"]
                    if terminal <= 18 and gap <= 9:
                        tilted = polynomial(terminal, terminal - 1, gap)
                        upper = polynomial(terminal, terminal, gap)
                        lower = polynomial(terminal - 1, terminal - 1, gap)
                        theta = s.Poly(X * tilted.diff().as_expr(), X)
                        assert (upper - tilted - s.Rational(gap, terminal) * theta).is_zero
                        assert (lower - tilted + s.Rational(1, terminal) * theta).is_zero
                        assert s.gcd(upper, lower).degree() == 0
                        bridge_count += 1
    gamma_bridge_count = 0
    for p, r, exponent, m, beta in ((3, 2, 1, 1, 2), (3, 1, 1, 2, 4),
                                    (5, 1, 1, 2, 20), (5, 0, 1, 3, 4),
                                    (2, 2, 1, 1, 3)):
        result = lower_one_certificate(p, r, exponent, m, beta)
        terminal, gap = result["terminal_N"], result["gap"]
        offset = terminal + beta - 2
        gamma = offset + 1
        literal_certificate(p, r, exponent, m, offset)
        tilted = polynomial(terminal, offset, gap)
        upper = polynomial(terminal, offset + 1, gap)
        lower = polynomial(terminal - 1, offset, gap)
        theta = s.Poly(X * tilted.diff().as_expr(), X)
        assert (upper - tilted - s.Rational(gap, gamma) * theta).is_zero
        assert (lower - tilted + s.Rational(1, terminal) * theta).is_zero
        assert s.gcd(upper, lower).degree() == 0
        gamma_bridge_count += 1
    for p in (2, 3, 5):
        for r in (1, 2):
            terminal = p ** r
            for exponent in (1, 2):
                for lower in (1, 2, 4, 7, 11):
                    digit = lower // p ** exponent % p
                    if digit == p - 1:
                        try:
                            general_lower_certificate(p, r, exponent, lower)
                        except ValueError:
                            rejected_count += 1
                        else:
                            raise AssertionError("carry case was accepted")
                        continue
                    result = general_lower_certificate(p, r, exponent, lower)
                    literal_certificate(p, r, exponent, 1, lower * terminal)
                    pair_count += 1
                    if terminal <= 9 and result["gap"] <= 9:
                        high = polynomial(terminal, lower * terminal, result["gap"])
                        low = polynomial(terminal - 1, lower * (terminal - 1), result["gap"])
                        assert s.gcd(high, low).degree() == 0
    for args in ((4, 1, 1, 1, 0), (3, 1, 0, 1, 0),
                 (3, 1, 1, 3, 0), (3, 1, 1, 1, 18),
                 (3, 1, 1, 2, 9)):
        try:
            class_certificate(*args)
        except ValueError:
            rejected_count += 1
        else:
            raise AssertionError("unsupported input accepted")
    minus_one_count = two_minus_one_count = 0
    for p in (2, 3, 5, 7):
        for r in range(1, 4):
            for exponent in (1, 2, 3):
                if p ** exponent > 50:
                    continue
                for m in range(1, p):
                    if m * p ** r > 150:
                        continue
                    minus_one_literal_certificate(p, r, exponent, m)
                    minus_one_count += 1
    for r in range(1, 5):
        gap_three_two_minus_one_literal_certificate(r)
        two_minus_one_count += 1
    # Unsupported exact parameters are not accepted by the new terminal families.
    for invalid in ((4, 1, 1), (3, 0, 1), (3, 1, 0)):
        try:
            minus_one_class_certificate(*invalid)
        except ValueError:
            rejected_count += 1
        else:
            raise AssertionError("unsupported minus-one input accepted")
    native_prime_count = 0
    for args in ((3, 2, 4, 13), (5, 3, 2, 11), (4, 3, 3, 13),
                 (5, 3, 4, 23), (7, 5, 4, 29), (8, 6, 3, 29),
                 (12, 10, 4, 53), (100, 3, 2, 67)):
        native_prime_literal_certificate(*args)
        native_prime_count += 1
    for args in ((12, 2, 2, 5), (3, 2, 4, 17), (3, 2, 4, 7),
                 (3, 2, 4, 15), (1, 3, 2, 5)):
        try:
            native_prime_certificate(*args)
        except ValueError:
            rejected_count += 1
        else:
            raise AssertionError("unsupported native-prime input accepted")
    forward_native_count = 0
    forward_profiles = set()
    for args in ((2, 2, 3, 1), (3, 3, 3, 2), (3, 3, 4, 1),
                 (5, 3, 4, 3), (3, 2, 4, 1), (8, 6, 3, 5),
                 (6, 6, 3, 1), (6, 6, 3, 5)):
        forward_result = forward_native_prime_literal_certificate(*args)
        forward_profiles.add(forward_result["next_nonzero_edge_width"])
        forward_native_count += 1
    assert forward_profiles == {None, 1, 2}
    for args in ((1, 1, 3, 1), (2, 2, 2, 1), (2, 3, 3, 1),
                 (3, 3, 5, 1)):
        try:
            forward_native_prime_certificate(*args)
        except ValueError:
            rejected_count += 1
        else:
            raise AssertionError("unsupported forward native-prime input accepted")
    split_count = 0
    split_failures_detected = 0
    for args in ((3, 1, 22, 21), (3, 1, 63, 62), (5, 2, 9, 17),
                 (2, 3, 7, 27), (7, 1, 11, 83)):
        split_result = split_literal_certificate(*args)
        split_count += 1
        split_failures_detected += not split_result["all_residuals_squarefree"]
    assert split_failures_detected >= 2
    common_quotient_count = 0
    for args in ((5, 1, 1, 3), (5, 2, 13, 3), (7, 1, 3, 4),
                 (7, 1, 5, 4), (7, 2, 17, 5), (11, 1, 5, 6),
                 (11, 1, 1, 7)):
        common_quotient_literal_certificate(*args)
        common_quotient_count += 1
    for args in ((3, 1, 1, 3), (5, 1, 2, 3), (7, 2, 25, 5),
                 (2, 1, 0, 2), (4, 1, 1, 2)):
        try:
            common_quotient_certificate(*args)
        except ValueError:
            rejected_count += 1
        else:
            raise AssertionError("unsupported common-quotient input accepted")
    precision_count = 0
    for h in (0, 1, 3, 4, 6, 7):
        global_precision_obstruction_literal_certificate(h)
        precision_count += 1
    for h in (2, 5):
        try:
            global_precision_obstruction_certificate(h)
        except ValueError:
            rejected_count += 1
        else:
            raise AssertionError("unsupported precision-obstruction input accepted")
    # Symbolic identities check the coefficient algebra used by the unbounded proof.
    symbolic_gamma = s.Symbol("gamma")
    previous = s.Poly(1, X)
    symbolic_disc_count = 0
    for degree in range(1, 7):
        current = s.Poly(sum(s.binomial(degree, k) * s.rf(symbolic_gamma, k)
                             * X ** (degree - k) for k in range(degree + 1)), X)
        assert (current.diff() - degree * previous).is_zero
        exact_formula = ((-1) ** (degree * (degree - 1) // 2)
            * prod(j ** j for j in range(1, degree + 1))
            * prod(s.rf(symbolic_gamma, i) for i in range(1, degree)))
        assert s.expand(s.discriminant(current) - exact_formula) == 0
        previous = current
        symbolic_disc_count += 1
    symbolic_q = s.Symbol("q")
    def unit_block(j):
        return (3 * j - 2) * (3 * j - 1)
    preferred = unit_block(symbolic_q + 1)
    difference = unit_block(symbolic_q + 2) * unit_block(symbolic_q + 3) - preferred ** 2
    identity = 9 * ((6 * symbolic_q + 8) * preferred
                    + 9 * (2 * symbolic_q + 2) * (4 * symbolic_q + 6))
    assert s.expand(difference - identity) == 0
    # Deliberate repeated-factor countercontrol for the imported residual checker.
    repeated = newton_residuals([1, 2, 1], 3)
    assert not all(row["residual_squarefree"] for row in repeated)
    # A carry can change the endpoint; rejection is mathematically meaningful.
    carry_hull = newton_residuals(coefficients(3, 18, 3), 3)
    assert carry_hull[-1]["right"][1] != 4
    return {
        "status": "pass", "arithmetic": "exact rational and finite-field",
        "literal_hull_residual_cases": literal_count,
        "lower_one_bridge_and_pair_cases": bridge_count,
        "additional_integer_Gamma_bridge_and_pair_cases": gamma_bridge_count,
        "general_lower_exponent_certificates": pair_count,
        "prime_power_minus_one_literal_cases": minus_one_count,
        "gap_three_two_minus_one_literal_cases": two_minus_one_count,
        "common_quotient_resultant_cases": common_quotient_count,
        "arbitrary_gap_native_prime_cases": native_prime_count,
        "forward_native_prime_slope_cases": forward_native_count,
        "forward_native_prime_width_profiles": ["no_nonzero_edge", 1, 2],
        "literal_factorial_split_cases": split_count,
        "literal_split_inseparable_residual_controls": split_failures_detected,
        "global_precision_obstruction_cases": precision_count,
        "symbolic_gap_one_discriminant_cases": symbolic_disc_count,
        "unsupported_inputs_rejected": rejected_count,
        "repeated_factor_countercontrol_detected": True,
        "carry_endpoint_countercontrol_detected": True,
        "lower_one_gap_three_terminal_family": "N=m*3^r or m*3^r-1, m=1,2, N>=2",
        "lower_one_prime_power_gap_extra_family": "gap=p^s, terminal N=m*p^r-1>=2, 1<=m<p",
        "ordinary_arbitrary_lower_family":
            "gap=p^s, terminal N=p^r>=2, digit_s(a)<=p-2",
        "ordinary_arbitrary_lower_extra_family":
            "gap=p^s, q=floor(aN/gap)=floor(a(N-1)/gap), q+N<p",
        "global_precision_method_obstruction_family": "N=4 or13 mod27",
        "arbitrary_gap_native_prime_family":
            "prime p>dN dividing some aN+j, 1<=j<=d, with valuation exactly1",
        "forward_native_prime_family":
            "a>=d, a>=2, N>=3, prime aN+j with 1<=j<=d certifies both adjoining pairs",
        "unrestricted_binomial_SFC": False, "all_terminal_orders": False,
        "finite_run_proves_unbounded_quantifiers": False,
        "proof": "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
        "python_version": platform.python_version(),
        "sympy_version": s.__version__,
        "independent_replay": False, "external_review": False,
        "formal_verification": False, "priority_claim": False}


if __name__ == "__main__":
    report = replay()
    here = Path(__file__).resolve().parent
    report["input_sha256"] = {
        name: hashlib.sha256((here / name).read_bytes()).hexdigest()
        for name in ("ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex",
                     "arithmetic_frontier_terminal_orders.py", "prime_power_windows.py",
                     "PRIME_POWER_TWO_TERM_SFC.tex", "HIGHER_GAP_ADJOINT_RANK.tex")}
    destination = here / "arithmetic_frontier_terminal_orders_verification.json"
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
