#!/usr/bin/env python3
"""Exact identities, complete F11 Bezout certificates, and actual-jet guards.

The all-order argument is in STEPANOV_LOWER_FOUR_ELEVENTH_GAP_CLOSURE.md.
The finite jet examples are transcription controls, not its proof.
"""

import argparse
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path
import platform


PRIME = 11


def trim(poly):
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result or [0]


def add(left, right):
    return trim([(left[i] if i < len(left) else 0)
                 + (right[i] if i < len(right) else 0)
                 for i in range(max(len(left), len(right)))])


def scale(poly, scalar):
    return trim([scalar * coefficient for coefficient in poly])


def sub(left, right):
    return add(left, scale(right, -1))


def mul(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return trim(result)


def derivative(poly):
    return trim([i * poly[i] for i in range(1, len(poly))])


def evaluate(poly, x):
    result = 0
    for coefficient in reversed(poly):
        result = result * x + coefficient
    return result


def reduction(poly, prime=PRIME):
    return trim([coefficient % prime for coefficient in poly])


def rising_polynomial(length):
    result = [1]
    for j in range(length):
        result = mul(result, [j, 1])
    return result


def divmod_field(poly, divisor, prime=PRIME):
    dividend, divisor = reduction(poly, prime), reduction(divisor, prime)
    assert divisor != [0]
    quotient = [0] * max(1, len(dividend) - len(divisor) + 1)
    inverse = pow(divisor[-1], -1, prime)
    while dividend != [0] and len(dividend) >= len(divisor):
        shift = len(dividend) - len(divisor)
        coefficient = dividend[-1] * inverse % prime
        quotient[shift] = coefficient
        dividend = reduction(sub(dividend, [0] * shift + scale(divisor, coefficient)), prime)
    return trim(quotient), dividend


def extended_gcd(left, right, prime=PRIME):
    old_r, r = reduction(left, prime), reduction(right, prime)
    old_s, s, old_t, t = [1], [0], [0], [1]
    while r != [0]:
        q, remainder = divmod_field(old_r, r, prime)
        old_r, r = r, remainder
        old_s, s = s, reduction(sub(old_s, mul(q, s)), prime)
        old_t, t = t, reduction(sub(old_t, mul(q, t)), prime)
    inverse = pow(old_r[-1], -1, prime)
    return tuple(reduction(scale(poly, inverse), prime) for poly in (old_r, old_s, old_t))


def verify_bezout(left, right, a, b, prime=PRIME):
    # This recomposes the polynomial identity rather than trusting a gcd flag.
    assert reduction(add(mul(a, left), mul(b, right)), prime) == [1]


def symbolic_top_coefficients():
    # Coefficients are integer polynomials in C. Construct the X coefficients
    # of -(X+C-1)(-1)^j(X)_j directly for j=9,8,7, corresponding to w1,w2,w3.
    top_rows = {}
    for w_index, j in enumerate([9, 8, 7], start=1):
        rising = rising_polynomial(j)
        sign = -((-1) ** j)
        columns = {}
        for degree in [10, 9, 8]:
            previous = rising[degree - 1] if degree - 1 < len(rising) else 0
            current = rising[degree] if degree < len(rising) else 0
            columns[degree] = scale([previous - current, current], sign)
        top_rows[w_index] = columns

    alpha = [546, 69, 3]
    beta = [30, 2]
    expected_scalar = {1: mul([2, 1], [3, 1]), 2: [3, 1], 3: [1]}
    for w_index, columns in top_rows.items():
        actual = add(sub(mul(alpha, columns[10]), mul(beta, columns[9])), columns[8])
        assert actual == expected_scalar[w_index]

    h = rising_polynomial(10)
    assert h[10:7:-1] == [1, 45, 870]
    numerator = add(sub(scale(alpha, h[10]), scale(beta, h[9])), [h[8]])
    assert numerator == [66, -21, 3]
    denominator_mod_121 = evaluate(h, 1) % 121
    constant_mod_121 = 66 * pow(denominator_mod_121, -1, 121) % 121
    assert denominator_mod_121 % 11 == 10
    assert constant_mod_121 == 55

    # C=15N-3 and C+3=15N: exact denominator cancellation.
    assert add([-3, 15], [3]) == [0, 15]
    native_residue = 3 * pow(15, -1, 11 ** 3) % (11 ** 3)
    assert native_residue == 1065
    return {
        "coefficient_convention": "ascending integer coefficients in C",
        "top_coefficients_by_w": {
            str(w_index): {"f" + str(degree): value for degree, value in columns.items()}
            for w_index, columns in top_rows.items()},
        "scalar_multipliers_for_f10_f9_f8": [alpha, scale(beta, -1), [1]],
        "H_top_coefficients": [1, 45, 870],
        "constant_numerator": numerator,
        "H_at_one_mod121": denominator_mod_121,
        "S0_mod121_when_C_divisible121": constant_mod_121,
        "remaining_terminal_mod1331": native_residue
    }


# Dual numbers a+b*u over F11, with u^2=0.
def dadd(left, right):
    return ((left[0] + right[0]) % PRIME, (left[1] + right[1]) % PRIME)


def dscale(value, scalar):
    return (value[0] * scalar % PRIME, value[1] * scalar % PRIME)


def dmul(left, right):
    return (left[0] * right[0] % PRIME,
            (left[0] * right[1] + left[1] * right[0]) % PRIME)


def dpow(value, exponent):
    result = (1, 0)
    for _ in range(exponent):
        result = dmul(result, value)
    return result


def dual_number_identity():
    rho = (1, 10)
    u = (0, 1)
    assert dadd(dadd(dpow(rho, 11), dscale(rho, -1)), dscale(u, -1)) == (0, 0)
    assert dpow(dadd(rho, (-1, 0)), 2) == (0, 0)
    h = [dpow(rho, 10 - degree) for degree in range(11)]
    h[0] = dadd(h[0], (-1, 0))
    product = [(0, 0)] * 12
    for degree, coefficient in enumerate(h):
        product[degree] = dadd(product[degree], dmul(dscale(rho, -1), coefficient))
        product[degree + 1] = dadd(product[degree + 1], coefficient)
    expected = [(0, 0)] * 12
    expected[0], expected[1], expected[11] = (0, 10), (10, 0), (1, 0)
    assert product == expected
    h_at_one = (0, 0)
    for coefficient in h:
        h_at_one = dadd(h_at_one, coefficient)
    assert h_at_one == (10, 0)
    f10, f9, f8 = [dscale(h[degree], -1) for degree in [10, 9, 8]]
    scalar = dadd(dadd(dscale(f10, 546), dscale(f9, -30)), f8)
    assert scalar == (0, 5)
    return {"ring": "F11[u]/(u^2)", "rho": rho,
            "quotient_product_verified": True, "H_at_one": h_at_one,
            "F_top_coefficients": [f10, f9, f8], "S": scalar}


def finite_field_certificates():
    ell = [0] * 9
    alternate = [0] * 9
    for k in range(9):
        ell[8 - k] = (-1) ** k * comb(8, k) * factorial(k + 1) % PRIME
        alternate[8 - k] = comb(8, k) * prod(range(10 - k, 10)) % PRIME
    assert ell == alternate == [1, 4, 1, 6, 7, 9, 3, 6, 1]
    assert evaluate(ell, -1) % PRIME == 10
    squarefree_gcd, squarefree_a, squarefree_b = extended_gcd(ell, derivative(ell))
    assert squarefree_gcd == [1]
    verify_bezout(ell, derivative(ell), squarefree_a, squarefree_b)

    witnesses = []
    for residue in range(1, PRIME):
        quadratic = [(-6 * residue) % PRIME, 1, 1]
        gcd, a, b = extended_gcd(ell, quadratic)
        assert gcd == [1]
        verify_bezout(ell, quadratic, a, b)
        witnesses.append({"d_mod11": residue, "quadratic": quadratic,
                          "A": a, "B": b,
                          "identity": "A*ell+B*quadratic=1 in F11[x]"})
    assert {record["d_mod11"] for record in witnesses} == set(range(1, PRIME))

    # A damaged witness must fail recomposition.
    damaged = add(witnesses[0]["A"], [1])
    rejected = False
    try:
        verify_bezout(ell, witnesses[0]["quadratic"], damaged, witnesses[0]["B"])
    except AssertionError:
        rejected = True
    assert rejected
    return {"prime": PRIME, "coefficient_convention": "ascending, residues 0 through 10",
            "ell": ell, "ell_at_minus_one": 10,
            "squarefree_witness": {"A": squarefree_a, "B": squarefree_b,
                                   "identity": "A*ell+B*ell'=1 in F11[x]"},
            "complete_nonzero_residue_witnesses": witnesses,
            "deliberately_corrupted_witness_rejected": rejected}


def actual_lower_reduction(terminal, ell):
    n = terminal - 1
    gamma = 4 * n + 1
    c = gamma + 11 * n + 11
    assert c % 121 == 0 and gamma % 11 == 0 and n % 11 == 8
    normalized_rising = 1
    coefficients = []
    for k in range(11):
        if k:
            block = k - 1
            # The actual product of eleven consecutive factors, divided by
            # exactly one factor 11. No Euler/cyclic quotient is involved.
            first = gamma + 11 * block
            assert first % 11 == 0
            normalized_block = (first // 11) * prod(first + j for j in range(1, 11))
            normalized_rising = normalized_rising * normalized_block % 11
        coefficient = comb(n, k) * normalized_rising % 11
        coefficients.append(coefficient)
    assert coefficients[:9] == list(reversed(ell))
    assert coefficients[9:] == [0, 0]
    assert normalized_rising == 0
    # For all k>=10, the integer normalized block recurrence multiplies this
    # zero by further integers, proving the entire remaining tail is zero.
    return {"N": terminal, "n": n, "gamma": gamma, "C": c,
            "coefficient_k_0_through_10_mod11": coefficients,
            "tail_from_k10_vanishes_by_integer_block_recurrence": True,
            "reduction": "x^(n-8)*ell(x)"}


def inverse_matrix(matrix, modulus):
    size = len(matrix)
    rows = [[value % modulus for value in row]
            + [int(i == j) for j in range(size)] for i, row in enumerate(matrix)]
    determinant_mod11 = 1
    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row][column] % PRIME)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            determinant_mod11 = -determinant_mod11
        value = rows[column][column]
        determinant_mod11 = determinant_mod11 * value % PRIME
        inverse = pow(value, -1, modulus)
        rows[column] = [entry * inverse % modulus for entry in rows[column]]
        for row in range(size):
            if row != column:
                factor = rows[row][column]
                rows[row] = [(a - factor * b) % modulus
                             for a, b in zip(rows[row], rows[column])]
    assert all(rows[i][j] == int(i == j) for i in range(size) for j in range(size))
    inverse = [row[size:] for row in rows]
    assert all(sum(matrix[i][k] * inverse[k][j] for k in range(size)) % modulus == int(i == j)
               for i in range(size) for j in range(size))
    assert determinant_mod11 != 0
    return inverse, determinant_mod11


def matrix_vector(matrix, vector, modulus):
    return [sum(a * b for a, b in zip(row, vector)) % modulus for row in matrix]


def valuation11(number):
    assert number != 0
    result = 0
    while number % PRIME == 0:
        result += 1
        number //= PRIME
    return result


def actual_jet(terminal, precision=6, u_degree=3):
    n = terminal - 1
    gamma = 4 * n + 1
    c = gamma + 11 * n + 11
    assert n >= 2 and gamma % 11 == 0 and c % 121 == 0
    modulus = PRIME ** precision
    size = PRIME * (u_degree + 1)
    b = [1] + [0] * (size - 1)
    lower = None
    for j in range(n + 1):
        if j == n:
            lower = b[:]
        # Original polynomial coefficients, before projection. The Euler
        # derivative multiplies by full k, never k modulo eleven.
        b = [((b[k - 1] if k else 0) + (11 * j + gamma + 10 - k) * b[k]) % modulus
             for k in range(size)]
    ls = [b]
    for j in range(9):
        previous = ls[-1]
        ls.append([((k - j) * previous[k] - (previous[k - 1] if k else 0)) % modulus
                   for k in range(size)])
    matrices, rhs = [], []
    for degree in range(u_degree + 1):
        matrices.append([[-ls[10 - m][11 * degree + row] % modulus
                          for m in range(1, 11)] for row in range(1, 11)])
        rhs.append([-lower[11 * degree + row] % modulus for row in range(1, 11)])
    inverse, determinant = inverse_matrix(matrices[0], modulus)
    kernel = []
    for degree in range(u_degree + 1):
        vector = rhs[degree][:]
        for i in range(1, degree + 1):
            contribution = matrix_vector(matrices[i], kernel[degree - i], modulus)
            vector = [(a - b) % modulus for a, b in zip(vector, contribution)]
        kernel.append(matrix_vector(inverse, vector, modulus))
    for degree in range(u_degree + 1):
        lhs = [0] * 10
        for i in range(degree + 1):
            contribution = matrix_vector(matrices[i], kernel[degree - i], modulus)
            lhs = [(a + b) % modulus for a, b in zip(lhs, contribution)]
        assert lhs == rhs[degree]
    scalar = [((c + 2) * (c + 3) * vector[0] + (c + 3) * vector[1] + vector[2]) % modulus
              for vector in kernel]
    assert scalar[0] % 121 == 55
    assert scalar[1] % 11 == 5
    constant_modulus = PRIME ** min(n, precision)
    denominator = prod(1 - c + j for j in range(10))
    assert denominator % 11 != 0
    predicted_constant = (3 * (c * c - 7 * c + 22)
                          * pow(denominator, -1, constant_modulus)) % constant_modulus
    assert scalar[0] % constant_modulus == predicted_constant
    return {"N": terminal, "n": n, "actual_lower_shape_gamma": gamma,
            "actual_upper_moment_shape": gamma + 4, "C": c, "v11_C": valuation11(c),
            "coefficient_modulus": modulus, "u_degree": u_degree,
            "minor_determinant_mod11": determinant,
            "normalized_last_component": 1, "kernel_first_ten_coefficients": kernel,
            "all_retained_rows_recomposed": True, "S_coefficients": scalar,
            "S0_over_11_mod11": scalar[0] // 11 % 11, "S1_mod11": scalar[1] % 11,
            "constant_formula_modulus": constant_modulus,
            "constant_formula_value": predicted_constant}


def run(repo_root):
    symbolic = symbolic_top_coefficients()
    dual = dual_number_identity()
    finite = finite_field_certificates()
    terminals = [97, 1065, 2396, 11713]
    lower = [actual_lower_reduction(terminal, finite["ell"]) for terminal in terminals]
    jets = [actual_jet(terminal) for terminal in terminals]
    assert [record["v11_C"] for record in jets] == [2, 3, 3, 4]
    source_dir = Path("research/factorial/image-rigidity")
    names = ["STEPANOV_LOWER_FOUR_ELEVENTH_GAP_CLOSURE.md",
             "verify_stepanov_lower_four_eleventh_gap.py"]
    pins = [{"path": str(source_dir / name),
             "sha256": sha256((repo_root / source_dir / name).read_bytes()).hexdigest()}
            for name in names]
    return {
        "status": "PASS", "base_commit": "19eb9912fb60142f8b620dac24896931d4ff3e46",
        "python_version": platform.python_version(), "dependencies": "Python standard library only",
        "scope": "universal polynomial identities, complete ten-case F11 certificates, and finite actual-jet guards",
        "all_order_closure_established_by_finite_probes": False,
        "all_order_proof": "STEPANOV_LOWER_FOUR_ELEVENTH_GAP_CLOSURE.md, Sections 2-6",
        "source_pins": pins, "symbolic_top_coefficient_checks": symbolic,
        "dual_number_check": dual, "complete_finite_field_certificate": finite,
        "actual_lower_reduction_guards": lower, "actual_normalized_jet_guards": jets,
        "assurance": "internal exact replay; no external review or formal verification claimed"
    }


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O or -OO: exact assertions must execute.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run(args.repo_root.resolve())
    output = args.output or Path(__file__).with_name("stepanov_lower_four_eleventh_gap_verification.json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "complete_F11_bezout_cases": 10,
                      "actual_jet_probes": len(report["actual_normalized_jet_guards"]),
                      "source_pins": len(report["source_pins"]), "receipt": str(output)}))
