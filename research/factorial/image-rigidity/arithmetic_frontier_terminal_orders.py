"""Exact certificates for the no-carry, single-digit terminal-order theorem.

Written proof: ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex.
A finite replay checks implementation, not the unbounded theorem.
"""
from fractions import Fraction
from math import comb, factorial, prod
from pathlib import Path
import hashlib
import json
import platform
import sympy as s

from prime_power_windows import newton_residuals

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
        "unsupported_inputs_rejected": rejected_count,
        "repeated_factor_countercontrol_detected": True,
        "carry_endpoint_countercontrol_detected": True,
        "lower_one_gap_three_terminal_family": "N=3^r or 2*3^r, N>=2",
        "ordinary_arbitrary_lower_family":
            "gap=p^s, terminal N=p^r>=2, digit_s(a)<=p-2",
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
