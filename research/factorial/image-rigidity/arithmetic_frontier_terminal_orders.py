"""Exact arithmetic criteria and method-boundary certificates for binomial moments.

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




def cubic_factorial_unit_mod(quotient, tail, index, modulus):
    """Literal nonmultiples of three in the factorial interval."""
    start = 3 * quotient + tail
    value = 1
    for number in range(start + 1, start + 3 * index + 1):
        if number % 3:
            value = value * number % modulus
    return value


def cubic_unit_mod81_formula(quotient, tail, index):
    alpha, beta = ((quotient + 1, 1), (7 * quotient + 6, 7),
                   (quotient - 7, 1))[tail]
    return (1 + 9 * (alpha * comb(index, 2) + beta * comb(index, 3))) % 81


def fraction_mod_unit(value, modulus):
    value = Fraction(value)
    assert value.denominator % 3
    return value.numerator * pow(value.denominator, -1, modulus) % modulus


def cubic_unit_log_truncation(quotient, tail, index, precision):
    """Exact rational representative of the retained diagonal eigenvalue."""
    answer = Fraction(0)
    base = comb(quotient + 1, 2)
    # Beyond these finite ranges the proved valuation bounds are >= precision.
    for order in range(1, precision):
        if 2 * order - valuation(order, 3) < precision:
            summed = sum(comb(quotient + j, 2) ** order
                         for j in range(1, index + 1))
            answer += Fraction((-1) ** (order + 1) * 9 ** order
                               * (summed - index * base ** order), order)
    for order in range(2, 2 * precision):
        if order - valuation(order, 3) < precision:
            endpoint = sum((Fraction(1, (3 * quotient + a) ** order)
                            for a in range(1, tail + 1)), Fraction(0))
            answer += (Fraction((-1) ** (order + 1) * 3 ** order
                                * (index ** order - index), order) * endpoint)
    return answer


def cubic_unit_exp_truncation(logarithm, precision):
    return sum((logarithm ** order / factorial(order)
                for order in range(precision)
                if 2 * order - factorial_valuation(order, 3) < precision),
               Fraction(0))


def cubic_unit_operator_replay():
    mod81_count = 0
    for quotient in list(range(41)) + [80, 81, 82, 242, 243, 244, 728, 729, 730]:
        for tail in range(3):
            first = cubic_factorial_unit_mod(quotient, tail, 1, 81)
            for index in list(range(81)) + [242, 243, 244, 728, 729, 730]:
                actual = (cubic_factorial_unit_mod(quotient, tail, index, 81)
                          * pow(first, -index, 81)) % 81
                assert actual == cubic_unit_mod81_formula(quotient, tail, index)
                mod81_count += 1
    logarithm_count = 0
    exponential_count = 0
    for precision in range(1, 10):
        modulus = 3 ** precision
        for quotient in [0, 1, 2, 7, 22, 42, 67]:
            for tail in range(3):
                for index in [0, 1, 2, 3, 5, 9, 17]:
                    first = cubic_factorial_unit_mod(quotient, tail, 1, modulus)
                    actual = (cubic_factorial_unit_mod(quotient, tail, index, modulus)
                              * pow(first, -index, modulus)) % modulus
                    direct_log = sum(
                        (Fraction((-1) ** (order + 1) * (actual - 1) ** order, order)
                         for order in range(1, precision)
                         if 2 * order - valuation(order, 3) < precision), Fraction(0))
                    formula_log = cubic_unit_log_truncation(
                        quotient, tail, index, precision)
                    assert fraction_mod_unit(direct_log - formula_log, modulus) == 0
                    logarithm_count += 1
                    reconstructed = cubic_unit_exp_truncation(formula_log, precision)
                    assert fraction_mod_unit(reconstructed, modulus) == actual
                    exponential_count += 1
    polynomial_count = 0
    for quotient in [0, 1, 7, 22]:
        for tail in range(3):
            for terminal in [1, 2, 5, 9]:
                literal = coefficients(terminal, 3 * quotient + tail, 3)
                for precision in [1, 2, 3, 5, 8]:
                    modulus = 3 ** precision
                    first = cubic_factorial_unit_mod(quotient, tail, 1, modulus)
                    for index, coefficient in enumerate(literal):
                        divided, remainder = divmod(coefficient, 3 ** index)
                        assert remainder == 0
                        actual = divided * pow(first, -index, modulus) % modulus
                        log_value = cubic_unit_log_truncation(
                            quotient, tail, index, precision)
                        unit = cubic_unit_exp_truncation(log_value, precision)
                        base = comb(terminal, index) * prod(
                            range(quotient + 1, quotient + index + 1))
                        assert fraction_mod_unit(base * unit, modulus) == actual
                    polynomial_count += 1
    isometry_count = 0
    block = lambda j: 2 + 9 * j * (j - 1)
    for quotient in range(20):
        for depth in range(21):
            other = quotient + 3 ** depth
            difference = (Fraction(block(other + 2), block(other + 1))
                          - Fraction(block(quotient + 2), block(quotient + 1)))
            assert valuation(difference, 3) == depth + 2
            isometry_count += 1
    # Countercontrols reject a frozen offset and a corrupted cubic correction.
    first = cubic_factorial_unit_mod(42, 0, 1, 81)
    actual = cubic_factorial_unit_mod(42, 0, 2, 81) * pow(first, -2, 81) % 81
    assert actual != cubic_unit_mod81_formula(22, 0, 2)
    first = cubic_factorial_unit_mod(0, 1, 1, 81)
    actual = cubic_factorial_unit_mod(0, 1, 3, 81) * pow(first, -3, 81) % 81
    corrupted = (1 + 9 * (6 * comb(3, 2) + comb(3, 3))) % 81
    assert actual != corrupted
    return {
        "scope": "all q,k>=0,c=0,1,2,N>=1 by the written operator proof",
        "mod81_literal_unit_cases": mod81_count,
        "tested_precision_exponents": list(range(1, 10)),
        "logarithm_truncation_cases": logarithm_count,
        "exponential_truncation_cases": exponential_count,
        "literal_factorial_operator_polynomial_cases": polynomial_count,
        "offset_isometry_cases": isometry_count,
        "frozen_offset_countercontrol_rejected": True,
        "corrupted_cubic_countercontrol_rejected": True,
        "squarefreeness_preservation_claimed": False,
        "finite_tests_prove_unbounded_scope": False,
    }



def prime_gap_cyclic_determinant(terminal, prime, gamma, z, u):
    """Literal integral recurrence, unsigned cofactors, and descending minor."""
    rows = []
    auxiliary = s.Poly(1, z)
    for degree in range(terminal + prime - 2):
        auxiliary = s.Poly(s.expand(
            (z + prime * degree + gamma + prime - 1) * auxiliary.as_expr()
            - z * auxiliary.diff().as_expr()), z)
        if degree + 1 >= terminal:
            rows.append(auxiliary)
    cyclic = s.Matrix(prime - 1, prime, lambda i, j: sum(
        rows[i].nth(j + prime * k) * u ** k
        for k in range((terminal + prime - 2) // prime + 1)))
    minors = [cyclic[:, [k for k in range(prime) if k != j]].det(method="domain-ge")
              for j in range(prime)]
    W = s.Poly(s.expand(sum(z ** (prime - 1 - j)
        * minors[j].subs(u, (-1) ** prime * z ** prime)
        for j in range(prime))), z)
    assert W.degree() == (prime - 1) * terminal and W.LC() == 1
    descending = W.all_coeffs()
    at = lambda j: descending[j] if 0 <= j < len(descending) else s.Integer(0)
    H = s.Matrix(terminal - 1, terminal - 1,
                 lambda i, j: at(prime * (j + 1) - (i + 1)))
    return s.Poly(H.det(method="domain-ge"), gamma)


def prime_gap_unit_replay():
    gamma, u, lam, z = s.symbols("gamma u lambda z")
    matrix_count = 0
    periodicity_count = 0
    determinant_count = 0
    constant_count = 0
    ordinary_pair_count = 0
    def reduce_prime(expr, prime, *variables):
        return s.Poly(s.expand(expr), *variables, modulus=prime).as_expr()
    for prime in [2, 3, 5, 7]:
        operator = s.diag(*(gamma - 1 - j for j in range(prime)))
        for j in range(prime - 1):
            operator[j + 1, j] = 1
        operator[0, prime - 1] = u
        Q = u + gamma ** prime - gamma
        target = lam ** prime - lam - Q
        assert reduce_prime(operator.charpoly(lam).as_expr() - target,
                            prime, lam, gamma, u) == 0
        e0 = s.eye(prime)[:, 0]
        columns = [e0]
        for j in range(1, prime):
            columns.append(operator * columns[-1])
        basis = s.Matrix.hstack(*columns)
        assert s.expand(basis.det()) == 1
        matrix_count += 1
        inverse = basis.inv()
        reduced = operator.subs(u, -gamma ** prime + gamma).applyfunc(
            lambda value: reduce_prime(value, prime, gamma))
        powers = [e0]
        for j in range(3 * prime + 2):
            powers.append((reduced * powers[-1]).applyfunc(
                lambda value: reduce_prime(value, prime, gamma)))
        for terminal in [1, 2, prime - 1, prime, prime + 1, 2 * prime + 3]:
            for j in range(prime - 1):
                coordinate = (inverse * powers[terminal + j]).applyfunc(
                    lambda value: reduce_prime(value, prime, gamma))
                expected = 1 + ((terminal + j - 1) % (prime - 1))
                assert coordinate == s.eye(prime)[:, expected]
                periodicity_count += 1
    for prime, bound in [(2, 7), (3, 7), (5, 4), (7, 3)]:
        for terminal in range(1, bound + 1):
            determinant = prime_gap_cyclic_determinant(
                terminal, prime, gamma, z, u)
            sign = (-1) ** (terminal * (terminal - 1) // 2)
            assert (determinant - s.Poly(sign, gamma)).trunc(prime).is_zero
            determinant_count += 1
    for prime in [2, 3, 5, 7, 11]:
        for terminal in range(1, 13):
            N = terminal * (terminal - 1) // 2
            constant = Fraction(factorial(prime * terminal - 1) ** (terminal - 1),
                prime ** N * prod(k ** k * factorial(prime * k - 1)
                                  for k in range(1, terminal)))
            units = [prod(number for number in range(1, prime * k)
                          if number % prime) for k in range(1, terminal + 1)]
            telescoped = Fraction(units[-1] ** (terminal - 1), prod(units[:-1]))
            assert constant == telescoped
            assert valuation(constant, prime) == 0
            assert residue(constant, prime) == (-1) ** N % prime
            constant_count += 1
    for prime in [2, 3, 5, 7]:
        for terminal in range(2, 10):
            upper = polynomial(terminal, terminal, prime)
            lower = polynomial(terminal - 1, terminal - 1, prime)
            assert s.gcd(upper, lower).degree() == 0
            ordinary_pair_count += 1
    # Wrong characteristic constant and composite-gap generalization must fail.
    prime = 3
    operator = s.diag(*(gamma - 1 - j for j in range(prime)))
    for j in range(prime - 1):
        operator[j + 1, j] = 1
    operator[0, prime - 1] = u
    wrong = lam ** prime - lam - u - gamma ** prime
    assert reduce_prime(operator.charpoly(lam).as_expr() - wrong,
                        prime, lam, gamma, u) != 0
    composite = prime_gap_cyclic_determinant(2, 4, gamma, z, u).trunc(2)
    assert composite.is_zero or composite.degree() > 0
    return {
        "scope": "all prime gaps p, all n>=1, rational p-integral positive gamma by written proof",
        "characteristic_polynomial_and_unimodular_basis_cases": matrix_count,
        "positive_power_periodicity_column_cases": periodicity_count,
        "symbolic_cyclic_determinant_congruence_cases": determinant_count,
        "constant_normalization_and_unit_telescope_cases": constant_count,
        "literal_ordinary_lower_one_pair_cases": ordinary_pair_count,
        "wrong_characteristic_polynomial_rejected": True,
        "prime_unit_congruence_composite_extension_rejected": True,
        "all_terminal_prime_gap_lower_one": True,
        "all_positive_real_Gamma_shapes_claimed": False,
        "finite_tests_prove_unbounded_scope": False,
    }


def prime_power_frobenius_replay():
    x,X,g=s.symbols('x X g')
    checks=0;profiles=[]

    def vp(v,p):
     v=Fraction(v);n=v.numerator;d=v.denominator
     if not n:return 10**9
     z=0
     while n%p==0:n//=p;z+=1
     while d%p==0:d//=p;z-=1
     return z

    def rank_mod(rows,p):
     if not rows:return 0
     A=[[int(a)%p for a in row] for row in rows];r=0
     for col in range(len(A[0])):
      pivot=next((k for k in range(r,len(A)) if A[k][col]),None)
      if pivot is None:continue
      A[r],A[pivot]=A[pivot],A[r]
      inv=pow(A[r][col],-1,p);A[r]=[(a*inv)%p for a in A[r]]
      for k in range(len(A)):
       if k!=r and A[k][col]:
        fac=A[k][col];A[k]=[(a-fac*b)%p for a,b in zip(A[k],A[r])]
      r+=1
      if r==len(A):break
     return r

    def rowmul(row,A,p):
     return [sum(row[k]*A[k][j] for k in range(len(row)))%p for j in range(len(row))]

    def krylov(p,c,t,N,kappa=1):
     A=[[0]*p for _ in range(p)]
     for i in range(p):
      A[i][i]=(c-i)%p
      if i<p-1:A[i][i+1]=1
      else:A[i][0]=t%p
     out=[];row=[kappa%p]+[0]*(p-1)
     for _ in range(N+1):out.append(row);row=rowmul(row,A,p)
     return out

    for p in [2,3,5]:
     start=checks;q=p*p
     # AH coefficients via zH'=-sum(z^(p^r))H; no floating arithmetic.
     N=2*p*p+2;hh=[Fraction(1)]
     powers=[];k=1
     while k<=N:powers.append(k);k*=p
     for n in range(1,N+1):hh.append(-sum(hh[n-k] for k in powers if k<=n)/n)
     for a in hh:assert vp(a,p)>=0;checks+=1
     for j in range(p):
      assert hh[j]==Fraction((-1)**j,factorial(j));checks+=1
      for k in range(1,(N-j)//p+1):
       assert vp(hh[p*k+j]*p**k,p)>=k;checks+=1
     # Closed p-typical potential class; p2 included.
     coeff=[p*(k+2) for k in range(5)];b=coeff[0]-1
     assert b%p;checks+=1
     for r in range(1,4):
      a=Fraction((coeff[r]-1)*p**(p**r-1),b**(p**r))
      assert vp(a,p)>=p**r-1>=1;checks+=1
     # Explicit local Hermite polynomials and q-1 degree-restricted basis.
     H={(j,r):s.Poly((x-j)**r*(1-(x-j)**(p-1))**p,x) for j in range(p) for r in range(p)}
     for (j,r),poly in H.items():
      for a in range(p):
       value=s.Poly(poly.as_expr().subs(x,a+p*X).expand(),X)
       if a==j:
        scaled=s.Poly(value.as_expr()/p**r-X**r,X)
        assert all(vp(c,p)>=1 for c in scaled.all_coeffs());checks+=1
       else:
        assert all(vp(c,p)>=p for c in value.all_coeffs());checks+=1
     rows=[H[j,r] for j in range(p) for r in range(p-1)]
     rows += [s.Poly(H[j,p-1].as_expr()-H[0,p-1].as_expr(),x) for j in range(1,p)]
     assert len(rows)==q-1 and all(h.degree()<=q-2 for h in rows);checks+=1
     assert rank_mod([[h.nth(k) for k in range(q-1)] for h in rows],p)==q-1;checks+=1
     # Normalized one-step block matrix: all p residue values for c_j,t;
     # large n>=p included to reject any accidental n<p argument.
     ns=sorted(set([1,2,p,p+1,2*p+3]))
     for n in ns:
      for c in range(p):
       for t in range(p):
        ks=[krylov(p,(c+j)%p,t,n+p,((-1)**j*pow(factorial(j),-1,p))%p) for j in range(p)]
        mat=[]
        for j in range(p):
         for r in range(p-1):
          row=[0]*q
          vec=ks[j][n+r] if j==0 else [(pow(j,n,p)*a)%p for a in ks[j][r]]
          row[j*p:(j+1)*p]=vec;mat.append(row)
        for j in range(1,p):
         row=[0]*q
         row[j*p:(j+1)*p]=[(pow(j,n,p)*a)%p for a in ks[j][p-1]]
         row[0:p]=[(-a)%p for a in ks[0][n+p-1]]
         mat.append(row)
        assert rank_mod(mat,p)==q-1;checks+=1
        assert rank_mod(ks[0][n:n+p-1],p)==p-1;checks+=1
     # Primitive polynomial reduction-degree contraction controls.
     for height in [p,p*p]:
      P=s.Poly((x**(p-1)-1)**height*x**(height-2)+p*(1+x**(p*height+2)),x)
      unit_indices=[k for k in range(P.degree()+1) if vp(P.nth(k),p)==0]
      r0=min(unit_indices);assert r0<=height-2;checks+=1
      M=min(vp(P.nth(k),p)+k for k in range(P.degree()+1))
      Q=s.Poly(P.as_expr().subs(x,p*x)/p**M,x)
      assert min(vp(c,p) for c in Q.all_coeffs())==0;checks+=1
      reduced_degree=max(k for k in range(Q.degree()+1) if vp(Q.nth(k),p)==0)
      assert reduced_degree<=height-2;checks+=1
     profiles.append({'p':p,'prime_square_normalized_rows':q-1,'n_values':ns,'checks':checks-start})

    # Exact diagnostic of the sharper valuation conjecture.
    valuation_cases=0
    for p,ss in [(2,2),(2,3),(3,2)]:
     q=p**ss;e=(q-1)//(p-1)-ss
     for n in range(2,8):
      N=n*(n-1)//2
      for gamma in range(1,10):
       T=sum(comb(n,k)*s.rf(gamma,q*k)*x**k for k in range(n+1))
       disc=s.discriminant(T,x)
       fac=q**N*s.prod(k**k for k in range(1,n+1))*s.prod(s.rf(gamma,q*k) for k in range(1,n))
       assert vp(disc,p)-vp(fac,p)==e*N;checks+=1;valuation_cases+=1
    prime_power_pair_cases=0
    for p,stage,beta in [(2,2,s.Rational(1,3)),(2,3,s.Integer(1)),
                          (3,2,s.Rational(1,2)),(3,2,s.Integer(1))]:
        q=p**stage
        for n in [5,7,10,13]:
            upper=s.Poly(sum(comb(n,k)*s.rf(beta+n,q*k)*x**k
                             for k in range(n+1)),x)
            lower=s.Poly(sum(comb(n-1,k)*s.rf(beta+n-1,q*k)*x**k
                             for k in range(n)),x)
            assert s.gcd(upper,lower).degree()==0
            checks+=1;prime_power_pair_cases+=1
    parameter_closure_cases=0
    eta=s.Symbol("eta")
    for p in [2,3,5]:
        for stage in [2,3,4]:
            for center in range(p+2):
                for j in range(p):
                    residue=(center-j)%p
                    descendant=(center+p**stage*eta-j-residue)/p
                    expected=(center-j-residue)//p+p**(stage-1)*eta
                    assert s.expand(descendant-expected)==0
                    checks+=1;parameter_closure_cases+=1
    exact_center_cases=0
    u_center=s.Symbol("u_center")
    for p,stage in [(2,2),(2,3),(3,2)]:
        q=p**stage;C=(q-1)//(p-1)
        for n in [2,3]:
            delta=prime_gap_cyclic_determinant(n,q,g,x,u_center)
            N=n*(n-1)//2
            for m in range(q):
                formula=Fraction(factorial(q*n-m-1)**(n-1),
                    q**N*prod(k**k*factorial(q*k-m-1) for k in range(1,n)))
                assert Fraction(int(delta.eval(-m)))==formula
                assert vp(formula,p)==(C-stage)*N
                checks+=1;exact_center_cases+=1
    return {"all_recorded_checks":checks,"profiles":profiles,
        "prime_power_valuation_cases":valuation_cases,
        "literal_prime_power_lower_one_pairs":prime_power_pair_cases,
        "stage_parameter_closure_cases":parameter_closure_cases,
        "negative_integer_center_formula_cases":exact_center_cases,
        "valuation_status":"proved by root-free integer disks and exact center values",
        "proof_scope":"all s>=1: stage-dependent Artin-Hasse closure, Gauss-degree induction, exact valuations on integer disks",
        "finite_tests_prove_unbounded_scope":False,
        "all_real_positive_Gamma_shapes_claimed":False}



def quintic_common_pair_radius_replay():
    """Complete residue calculations supporting the uniform graded proof."""
    x, y = s.symbols("x y")
    e_table = {1: [0, 1, 1, 1, 3], -1: [0, 0, 1, 1, 1],
               -3: [0, 0, 0, 2, 4]}
    table_cases = 0
    vanished_sections = 0
    for delta, row in e_table.items():
        for state, e in enumerate(row):
            ceiling = (2 * state + delta + 4) // 5
            numerator = s.Poly((1 - x) ** (5 * ceiling - 2 * state - delta)
                * (1 - x + 2 * x ** 2) ** (4 - state), x, modulus=5)
            for residue in range(e):
                target = (state - residue) % 5
                assert all(numerator.nth(j) == 0
                           for j in range(target, numerator.degree() + 1, 5))
                vanished_sections += 1
            assert numerator.nth(state - e) % 5
            table_cases += 1
    def section(poly, residue):
        return s.Poly(sum(poly.nth(j) * y ** ((j - residue) // 5)
                         for j in range(residue, poly.degree() + 1, 5)),
                      y, modulus=5)
    further_flag_cases = 0
    for state, expression in [(2, 1 - 2 * x + x ** 3 - x ** 4),
        (3, 1 + x - 2 * x ** 2 + 2 * x ** 4 - 2 * x ** 5)]:
        numerator = s.Poly(expression, x, modulus=5)
        assert section(numerator, state).is_zero
        assert (2 * section(numerator, state - 1)
                - section(numerator, state - 2)).is_zero
        further_flag_cases += 1
    corrupted = s.Poly(1 - x + x ** 3 - x ** 4, x, modulus=5)
    assert not (2 * section(corrupted, 1) - section(corrupted, 0)).is_zero
    def raise_coefficients(row, constant):
        return [(row[j - 1] if j else 0)
                + (constant - j) * (row[j] if j < len(row) else 0)
                for j in range(len(row) + 1)]
    def auxiliary_coefficients(degree, gamma, gap):
        row = [1]
        for k in range(degree):
            row = raise_coefficients(row, gap * k + gamma + gap - 1)
        return row
    physical_bessel_cases = 0
    bessel = [1]
    for n in range(26):
        actual = auxiliary_coefficients(n, 2 * n + 1, 3)
        assert all((a - b) % 5 == 0 for a, b in zip(actual, bessel))
        bessel = raise_coefficients(bessel, 2 * n)
        physical_bessel_cases += 1
    flag_coefficient_cases = 0
    for m in range(2, 49):
        if m % 5 not in [2, 3]:
            continue
        ell, state = divmod(m, 5)
        L = 5 * ell + 1
        B = auxiliary_coefficients(m, 2 * m + 1, 3)
        C = raise_coefficients(B, 5 * m + 3)
        F = [C[j] - (5 * m + 3 - L) * (B[j] if j < len(B) else 0)
             for j in range(len(C))]
        assert F[5 * ell + 3] % 5
        for j, coefficient in enumerate(F):
            if coefficient:
                k = ell - j // 5
                assert valuation(coefficient, 5) >= k + (j % 5 < 3)
        flag_coefficient_cases += 1
    z, u = s.symbols("z u")
    basic = [s.Integer(1), z, z ** 2 + z,
             z ** 3 + 3 * z ** 2 + 3 * z, z ** 4 + z ** 3]
    expected = [s.Integer(1), s.Integer(1), u, u, u ** 2]
    unit_origin_cases = 0
    for index, B in enumerate(basic):
        C = s.Poly(s.expand((z + 3) * B - z * s.diff(B, z)), z)
        B = s.Poly(B, z)
        def cyclic(poly, residue):
            return sum(poly.nth(j) * u ** ((j - residue) // 3)
                       for j in range(residue, poly.degree() + 1, 3))
        cr = [cyclic(C, j) for j in range(3)]
        br = [cyclic(B, j) for j in range(3)]
        minors = [s.Poly(cr[i] * br[j] - cr[j] * br[i], u, modulus=5)
                  for i, j in [(0, 1), (0, 2), (1, 2)]]
        common = minors[0]
        for minor in minors[1:]:
            common = s.gcd(common, minor)
        assert common.monic() == s.Poly(expected[index], u, modulus=5)
        unit_origin_cases += 1
    return {
        "complete_annulus_numerator_cases": table_cases,
        "vanishing_residue_sections": vanished_sections,
        "complete_further_flag_numerator_cases": further_flag_cases,
        "corrupted_further_flag_identity_rejected": True,
        "physical_bessel_congruence_cases": physical_bessel_cases,
        "literal_further_flag_coefficient_cases": flag_coefficient_cases,
        "complete_unit_argument_cyclic_ideal_cases": unit_origin_cases,
        "uniform_claim": "outer argument excludes v5(u)<3/5",
        "finite_tests_prove_unbounded_scope": False,
        "boundary_exclusion_in_separate_profile_certificate": True,
        "all_terminal_quintic_SFC": False,
    }



def quintic_inner_endpoint_replay():
    """Exact coefficient controls for the uniform valuation-one restriction."""
    n, w = s.symbols("n w")
    modulus = 125

    def reduce_coefficients(expression):
        return s.Poly(expression, n, domain=s.ZZ).trunc(modulus).as_expr().expand()

    def coefficient_formula(index):
        if index == 0:
            return s.Integer(1)
        ceiling = (3 * index + 4) // 5
        factorial_part = factorial(index)
        order = index - ceiling - factorial_valuation(index, 5)
        if order >= 3:
            return s.Integer(0)
        value = 5 ** order * pow(
            factorial_part // 5 ** factorial_valuation(index, 5), -1, modulus)
        for j in range(index - ceiling):
            value = reduce_coefficients(value * (n - ceiling - j))
        for j in range(3 * index):
            if j % 5:
                inverse = pow(j, -1, modulus)
                reciprocal = -inverse * sum(
                    (5 * inverse * n) ** a for a in range(3 - order))
                value = reduce_coefficients(value * reciprocal)
        return value

    expected = [
        1, -25 * n ** 2 + 35 * n - 62, -13,
        25 * n ** 2 + 5 * n + 15, 25 * n ** 2 + 5 * n + 10,
        5 * n ** 2 - 35 * n + 60,
        50 * n ** 3 + 15 * n ** 2 - 60 * n + 50,
        -15 * n ** 2 + 40 * n + 50,
        25 * n ** 3 + 50 * n ** 2 + 50 * n,
        25 * n ** 3 - 25 * n ** 2 + 25 * n - 25,
        50 - 50 * n ** 4,
        -25 * n ** 4 - 25 * n ** 3 - 25 * n ** 2 - 25 * n,
        25 * n ** 4 + 50 * n ** 3 - 25 * n ** 2 - 50 * n,
    ]
    actual = [coefficient_formula(j) for j in range(13)]
    assert all(reduce_coefficients(a - b) == 0 for a, b in zip(actual, expected))
    truncated = sum(a * w ** j for j, a in enumerate(actual))
    assert s.Poly(truncated, n, w).trunc(5).as_expr() == 1 - 2 * w + 2 * w ** 2
    evaluated_two = reduce_coefficients(truncated.subs(w, 2))
    evaluated_nine = reduce_coefficients(truncated.subs(w, 9))
    assert evaluated_two == 25 * n ** 2 - 50 * n + 25
    assert evaluated_nine == 50 * n ** 2 - 25 * n - 50
    assert reduce_coefficients(evaluated_two.subs(n, n + 1) - evaluated_two) == (
        50 * n - 25)
    assert reduce_coefficients(evaluated_nine.subs(n, n + 1) - evaluated_nine) == (
        25 - 25 * n)
    for root in (2, 4):
        assert (1 + 3 * root + 2 * root ** 2) % 5 == 0
        assert (3 + 4 * root) % 5
    tail_exceptions = 0
    for j in range(13, 18):
        assert j - (3 * j + 4) // 5 - factorial_valuation(j, 5) == 3
        tail_exceptions += 1
    # e_j >= (3j-11)/20 > 2 for j>=18 proves the remaining infinite tail.
    assert Fraction(3 * 18 - 11, 20) > 2
    assert max((3 * j + 4) // 5 for j in range(13)) == 8

    def literal_coefficient(index, j):
        if j > index:
            return Fraction(0)
        return Fraction(comb(index, j) * 5 ** j * factorial(5 * index - 3 * j),
                        factorial(5 * index))

    def modular_fraction(value):
        assert value.denominator % 5
        return value.numerator * pow(value.denominator, -1, modulus) % modulus

    literal_cases = 0
    for index in range(25, 50):
        for j in range(13):
            assert modular_fraction(literal_coefficient(index, j)) == (
                int(actual[j].subs(n, index)) % modulus)
            literal_cases += 1
    threshold_cases = 0
    for index in range(8, 25):
        finite = sum(modular_fraction(literal_coefficient(index, j)) * w ** j
                     for j in range(index + 1))
        table = truncated.subs(n, index)
        assert s.Poly(finite - table, w, domain=s.ZZ).trunc(modulus).is_zero
        threshold_cases += 1

    complete_residue_cases = 0
    surviving_disks = set()
    for index in range(25, 50):
        rows = [
            [modular_fraction(literal_coefficient(degree, j))
             for j in range(degree + 1)]
            for degree in (index, index + 1)
        ]
        for unit in range(1, modulus):
            if unit % 5 == 0:
                continue
            values = [sum(a * pow(unit, j, modulus) for j, a in enumerate(row))
                      % modulus for row in rows]
            if values == [0, 0]:
                condition = (unit % 25 == 2 and index % 5 == 3) or (
                    unit % 25 == 9 and index % 5 == 1)
                assert condition
                surviving_disks.add((index % 5, unit % 25))
            complete_residue_cases += 1
    assert surviving_disks == {(3, 2), (1, 9)}

    # The canceled infinite family is not the finite moment at small indices.
    assert literal_coefficient(1, 2) == 0
    assert int(actual[2].subs(n, 1)) % modulus != 0
    corrupted = truncated + 25 * w ** 3
    assert reduce_coefficients(corrupted.subs(w, 2) - evaluated_two) != 0
    return {
        "scope": "all n>=8, actual consecutive quintic moments, v5(u)=1",
        "coefficient_table_identities": len(actual),
        "root_disk_evaluation_and_difference_identities": 4,
        "simple_residual_roots": 2,
        "tail_bound_exception_indices": tail_exceptions,
        "tail_bound_all_j_ge18": "(3j-11)/20>2",
        "literal_factorial_coefficient_cases": literal_cases,
        "finite_degree_threshold_cases": threshold_cases,
        "complete_precision_three_unit_cases": complete_residue_cases,
        "necessary_disks": [
            {"u_over_5_mod25": 2, "n_mod5": 3},
            {"u_over_5_mod25": 9, "n_mod5": 1},
        ],
        "excluded_index_classes_mod5": [0, 2, 4],
        "false_infinite_interpolation_control_rejected": True,
        "corrupted_coefficient_table_rejected": True,
        "other_inner_valuations_resolved": False,
        "all_terminal_quintic_SFC": False,
        "finite_tests_prove_unbounded_scope": False,
    }


def quintic_boundary_profiles_replay():
    """Literal complete boundary certificates; no gcd calculation is trusted."""
    x, y, t, v = s.symbols("x y t v")
    def falling(a, count):
        value = 1
        for j in range(count):
            value = value * (a - j) % 5
        return value
    def gcoeff(a, state, delta):
        ceiling = (2 * state + delta + 4) // 5
        numerator = s.Poly((1 - x) ** (5 * ceiling - 2 * state - delta)
            * (1 - x + 2 * x ** 2) ** (4 - state), x, modulus=5)
        H = s.Poly(s.series((1 - y) ** (-(2 * a + ceiling))
            * (1 - y + 2 * y ** 2) ** (-(a + 1)),
            y, 0, 5).removeO().expand(), y, modulus=5)
        return {j: sum(numerator.nth(j - 5 * k) * H.nth(k)
                       for k in range(5) if j >= 5 * k) % 5
                for j in range(24)}
    def leading(a, state, delta, e):
        G = gcoeff(a, state, delta)
        scale = factorial(state) * pow(factorial(e), -1, 5) % 5
        return s.Poly(sum(((-1) ** k * scale * falling(a, k)
            * G.get(5 * k + state - e, 0)) % 5 * t ** k
            for k in range(5)), t, modulus=5)
    def further(a, state, e, ee):
        G = gcoeff(a, state, 1)
        scale = factorial(state) * pow(factorial(ee), -1, 5) % 5
        return s.Poly(sum(((-1) ** k * scale * falling(a, k)
            * (ee * G.get(5 * k + state - ee + 1, 0)
               - (ee - e) * G.get(5 * k + state - ee, 0))) % 5 * t ** k
            for k in range(5)), t, modulus=5)
    def coefficients_mod5(poly):
        if poly.is_zero:
            return [0]
        return [int(poly.nth(j)) % 5 for j in range(poly.degree() + 1)]
    def sections(poly):
        return [s.Poly(sum(poly.nth(j) * v ** ((j - residue) // 3)
            for j in range(residue, poly.degree() + 1, 3)), v, modulus=5)
            for residue in range(3)]
    def from_coefficients(row):
        return s.Poly(sum(int(c) * v ** j for j, c in enumerate(row)),
                      v, modulus=5)
    def raise_coefficients(row, constant):
        return [(row[j - 1] if j else 0)
            + (constant - j) * (row[j] if j < len(row) else 0)
            for j in range(len(row) + 1)]
    def actual_auxiliary(degree):
        row = [1]
        for k in range(degree):
            row = raise_coefficients(row, 3 * k + 2 * degree + 3)
        return row
    def literal_profile(row, block, residue):
        values = []
        for k in range(block + 1):
            j = 5 * (block - k) + residue
            coefficient = row[j] if j < len(row) else 0
            assert coefficient % 5 ** k == 0
            values.append((coefficient // 5 ** k) % 5)
        assert all(value == 0 for value in values[5:])
        values = values[:5]
        while len(values) > 1 and values[-1] == 0:
            values.pop()
        return values
    certificate_path = Path(__file__).resolve().parent / (
        "arithmetic_frontier_quintic_boundary_profiles.json")
    data = json.loads(certificate_path.read_text())
    assert len(data["profiles"]) == 25
    flags = [0, 1, 1, 1, 3]
    seen = set()
    bezout_count = 0
    literal_count = 0
    profiles = []
    for record in data["profiles"]:
        state, a = record["m_mod5"], record["ell_mod5"]
        assert state in range(5) and a in range(5)
        assert (state, a) not in seen
        seen.add((state, a))
        e = flags[state]
        b = leading(a, state, 1, e)
        if state < 4:
            ee = e + (1 if state in [0, 1] else 2)
            c = further(a, state, e, ee)
            kind = "F=C-(5m+3-L)B"
        else:
            ee = 5
            c = leading((a + 1) % 5, 0, -1, 0)
            kind = "C"
        assert (e, ee, kind) == (record["B_exponent"],
            record["other_exponent"], record["other_row"])
        assert coefficients_mod5(b) == record["B_leading_t_coefficients"]
        assert coefficients_mod5(c) == record["other_leading_t_coefficients"]
        bp = s.Poly(s.expand(x ** (20 + e)
            * b.as_expr().subs(t, x ** -5)), x, modulus=5)
        cp = s.Poly(s.expand(x ** (20 + ee)
            * c.as_expr().subs(t, x ** -5)), x, modulus=5)
        br, cr = sections(bp), sections(cp)
        minors = [br[i] * cr[j] - br[j] * cr[i]
                  for i, j in [(0, 1), (0, 2), (1, 2)]]
        assert [coefficients_mod5(f) for f in minors] == (
            record["cyclic_minors_v_coefficients"])
        multipliers = [from_coefficients(row)
                       for row in record["Bezout_coefficients"]]
        target = s.Poly(v ** record["Bezout_target_origin_power"],
                        v, modulus=5)
        assert sum((a * b for a, b in zip(multipliers, minors)),
                   s.Poly(0, v, modulus=5)) == target
        bezout_count += 1
        profiles.append((multipliers, minors, target))
        # Independent reconstruction from literal integer auxiliary coefficients.
        ell = a + 5
        m = 5 * ell + state
        B = actual_auxiliary(m)
        C = raise_coefficients(B, 5 * m + 3)
        if state < 4:
            L = 5 * ell + e
            other = [C[j] - (5 * m + 3 - L)
                * (B[j] if j < len(B) else 0) for j in range(len(C))]
            block, residue = ell, ee
        else:
            other, block, residue = C, ell + 1, 0
        assert literal_profile(B, ell, e) == record["B_leading_t_coefficients"]
        assert literal_profile(other, block, residue) == (
            record["other_leading_t_coefficients"])
        literal_count += 1
    assert seen == {(state, a) for state in range(5) for a in range(5)}
    factorial_count = 0
    for ell in range(10):
        for state in range(5):
            n = 5 * ell + state
            for k in range(min(ell, 6) + 1):
                for i in range(5):
                    j = 5 * (ell - k) + i
                    if j > n + 1:
                        continue
                    ratio = Fraction(factorial(n), factorial(j) * 5 ** k)
                    assert ratio.denominator % 5
                    actual = ratio.numerator * pow(ratio.denominator, -1, 5) % 5
                    expected = ((-1) ** k * factorial(state)
                        * pow(factorial(i), -1, 5) * falling(ell, k)) % 5
                    assert actual == expected
                    factorial_count += 1
    # Changing a nonzero Bezout coefficient must destroy the retained identity.
    aa, mm, target = profiles[0]
    altered = aa[:]
    index = next(j for j, poly in enumerate(altered) if not poly.is_zero)
    altered[index] = altered[index] + s.Poly(1, v, modulus=5)
    assert sum((a * b for a, b in zip(altered, mm)),
               s.Poly(0, v, modulus=5)) != target
    return {
        "complete_profiles": len(seen),
        "literal_Bezout_identities": bezout_count,
        "independent_actual_auxiliary_profiles": literal_count,
        "factorial_block_reduction_cases": factorial_count,
        "corrupted_Bezout_control_rejected": True,
        "boundary_excluded": True,
        "strict_uniform_common_pair_radius": "v5(u)>3/5",
        "finite_reduction_exhausts_all_m": True,
        "finite_prefix_extrapolation": False,
        "all_terminal_quintic_SFC": False,
        "certificate": certificate_path.name,
    }


def integer_bareiss_determinant(matrix):
    """Fraction-free determinant with exact-division checks and row pivoting."""
    entries = [list(row) for row in matrix]
    size = len(entries)
    previous, sign = 1, 1
    for k in range(size - 1):
        if entries[k][k] == 0:
            pivot_row = next((j for j in range(k + 1, size) if entries[j][k]), None)
            if pivot_row is None:
                return 0
            entries[k], entries[pivot_row] = entries[pivot_row], entries[k]
            sign = -sign
        pivot = entries[k][k]
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                numerator = entries[i][j] * pivot - entries[i][k] * entries[k][j]
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                entries[i][j] = quotient
            entries[i][k] = 0
        previous = pivot
    return sign * entries[-1][-1]


def evaluated_moment_coefficients(terminal, gap, lower):
    """Descending coefficients of sum binom(j,k)*(j*a+1)_(d*k)*t^k."""
    running, ascending = 1, [1]
    for index in range(1, gap * terminal + 1):
        running *= terminal * lower + index
        if index % gap == 0:
            ascending.append(comb(terminal, index // gap) * running)
    return list(reversed(ascending))


def integer_sylvester_resultant(terminal, gap, lower):
    high = evaluated_moment_coefficients(terminal, gap, lower)
    low = evaluated_moment_coefficients(terminal - 1, gap, lower)
    size = 2 * terminal - 1
    rows = []
    for shift in range(terminal - 1):
        rows.append([0] * shift + high + [0] * (size - shift - len(high)))
    for shift in range(terminal):
        rows.append([0] * shift + low + [0] * (size - shift - len(low)))
    return integer_bareiss_determinant(rows)


def integer_horner(descending, value):
    result = 0
    for coefficient in descending:
        result = result * value + coefficient
    return result


def scaled_gaussian_taylor(descending, real_numerator, imag_numerator, denominator):
    """Return B_k=D^(degree-k)*[w^k]H((X+iY)/D+w), as Gaussian integers."""
    result = [(descending[0], 0)]
    power = 1
    for coefficient in descending[1:]:
        power *= denominator
        updated = []
        for k in range(len(result) + 1):
            real, imag = result[k] if k < len(result) else (0, 0)
            previous_real, previous_imag = result[k - 1] if k else (0, 0)
            updated.append((
                real_numerator * real - imag_numerator * imag + previous_real,
                imag_numerator * real + real_numerator * imag + previous_imag))
        updated[0] = (updated[0][0] + coefficient * power, updated[0][1])
        result = updated
    return result


def integer_rouche_bounds(descending, real_numerator, imag_numerator,
                          denominator, radius_denominator):
    assert denominator > 0 and radius_denominator > 0
    scale, remainder = divmod(denominator, radius_denominator)
    assert remainder == 0 and scale > 0
    taylor = scaled_gaussian_taylor(descending, real_numerator, imag_numerator,
                                   denominator)
    linear_lower = scale * max(abs(taylor[1][0]), abs(taylor[1][1]))
    norms = [abs(real) + abs(imag) for real, imag in taylor]
    tail = 0
    for k in range(len(norms) - 1, 1, -1):
        tail = tail * scale + norms[k]
    remainder_upper = norms[0] + tail * scale * scale
    return linear_lower, remainder_upper


def resultant_hurwitz_counterexample_certificate():
    """Verify the literal degree-280 resultant and one right-half-plane zero.

    No floating point, Arb, Routh recursion, or symbolic resultant is used here.
    Degree-bound interpolation connects 281 exact integer determinants to the
    literal coefficient input. An integer Rouché inequality then certifies a
    simple root in the displayed rational disk.
    """
    here = Path(__file__).resolve().parent
    path = here / "arithmetic_frontier_terminal_orders_resultant_hurwitz_counterexample.json"
    certificate = json.loads(path.read_text(encoding="utf-8"))
    assert certificate["schema_version"] == 1
    assert certificate["source_resultant"]["terminal_N"] == 8
    assert certificate["source_resultant"]["gap"] == 5
    assert certificate["coefficient_order"] == "descending"
    assert certificate["degree_upper_bound"] == 280
    assert certificate["identity_verification_nodes"] == {"first": 0, "last": 280}
    descending = list(map(int, certificate["primitive_coefficients"]))
    content = int(certificate["resultant_content"])
    assert content > 0 and len(descending) == 281
    assert all(coefficient > 0 for coefficient in descending)
    assert gcd(*descending) == 1
    for lower in range(281):
        assert integer_sylvester_resultant(8, 5, lower) == (
            content * integer_horner(descending, lower))
    center = certificate["center"]
    real_numerator = int(center["real_numerator"])
    imag_numerator = int(center["imag_numerator"])
    denominator = int(center["denominator"])
    radius_denominator = int(certificate["radius_denominator"])
    margin = int(certificate["rouche_strict_margin_multiplier"])
    assert real_numerator * radius_denominator > denominator
    assert imag_numerator * radius_denominator > denominator
    lhs, rhs = integer_rouche_bounds(descending, real_numerator, imag_numerator,
                                     denominator, radius_denominator)
    assert rhs > 0 and margin >= 1 and lhs > margin * rhs
    # Corrupt coefficients must fail the independent determinant input check.
    corrupted = list(descending)
    corrupted[-1] += 1
    assert integer_sylvester_resultant(8, 5, 0) != content * integer_horner(corrupted, 0)
    # A disk at the positive real point 1 is not certified as a root disk.
    wrong_lhs, wrong_rhs = integer_rouche_bounds(descending, denominator, 0,
                                                denominator, radius_denominator)
    assert wrong_lhs <= wrong_rhs
    # Independent small controls for the integer determinant and Taylor helpers.
    assert integer_bareiss_determinant([[1, 2], [3, 4]]) == -2
    assert integer_bareiss_determinant([[0, 1], [1, 0]]) == -1
    assert scaled_gaussian_taylor([1, 0, 1], 1, 1, 1) == [(1, 2), (2, 2), (1, 0)]
    return {
        "terminal_N": 8, "gap": 5, "resultant_degree": 280,
        "coefficient_identity_nodes_checked": 281,
        "literal_primitive_coefficients_all_positive": True,
        "exact_integer_rouche_zero_count_in_disk": 1,
        "rouche_zero_is_simple": True, "disk_strictly_in_open_right_half_plane": True,
        "rouche_linear_to_remainder_floor": str(lhs // rhs),
        "required_Arb_or_Routh_dependency": False,
        "corrupted_coefficient_and_wrong_disk_controls_rejected": True,
        "unrestricted_unshifted_resultant_Hurwitz_stability": False,
        "refutes_adjacent_gap_resultant_stability": False,
        "refutes_coefficient_positivity": False, "SFC_counterexample": False,
        "certificate": path.name}


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
    cubic_unit_operator = cubic_unit_operator_replay()
    prime_gap_units = prime_gap_unit_replay()
    prime_power_frobenius = prime_power_frobenius_replay()
    quintic_pair_radius = quintic_common_pair_radius_replay()
    quintic_boundary = quintic_boundary_profiles_replay()
    quintic_inner_endpoint = quintic_inner_endpoint_replay()
    hurwitz_counterexample = resultant_hurwitz_counterexample_certificate()
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
        "cubic_factorial_unit_operator": cubic_unit_operator,
        "prime_gap_unit_determinant": prime_gap_units,
        "prime_power_frobenius_noncollision": prime_power_frobenius,
        "quintic_five_adic_common_pair_radius": quintic_pair_radius,
        "quintic_boundary_profiles": quintic_boundary,
        "quintic_five_adic_inner_endpoint": quintic_inner_endpoint,
        "unshifted_resultant_hurwitz_counterexample": hurwitz_counterexample,
        "unsupported_inputs_rejected": rejected_count,
        "repeated_factor_countercontrol_detected": True,
        "carry_endpoint_countercontrol_detected": True,
        "lower_one_gap_three_terminal_family": "all N>=2 for ordinary factorials",
        "prime_gap_lower_one_all_terminal_orders": True,
        "prime_power_lower_one_all_terminal_orders": True,
        "ordinary_binomial_degree_at_most_four": True,
        "all_positive_real_Gamma_shapes_at_higher_prime_gaps": False,
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
                     "PRIME_POWER_TWO_TERM_SFC.tex", "HIGHER_GAP_ADJOINT_RANK.tex",
                     "arithmetic_frontier_terminal_orders_resultant_hurwitz_counterexample.json",
                     "sources/dlmf-1-10-rouche.html",
                     "gap_three_hurwitz_structure.tex", "ALL_UNIVARIATE_TWO_TERM_SFC.tex",
                     "UNIFORM_PADE_PARABOLA_SFC.tex", "CUBIC_ADJOINT_STABILITY_SFC.tex",
                     "DYADIC_BINOMIAL_SFC.tex", "lower_exponent_initial_and_boundary.tex",
                     "arithmetic_frontier_quintic_boundary_profiles.json")}
    destination = here / "arithmetic_frontier_terminal_orders_verification.json"
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
