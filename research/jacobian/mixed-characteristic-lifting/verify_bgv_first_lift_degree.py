#!/usr/bin/env python3
"""Replay the sharp bounds and the Frobenius-shifted branch obstruction modulo 81."""

from verify_bgv_odd_prime_lifting import (
    ONE,
    add,
    degree,
    family,
    jac,
    mod,
    quotient,
    radial_substitute,
    scale,
    sub,
)


CERTIFICATE = {
    (38, 9): 1,
    (40, 10): 2,
    (42, 11): 1,
    (67, 16): 1,
    (69, 17): 1,
    (73, 19): 1,
    (75, 20): 1,
    (85, 25): 2,
}

# Weights on the linear correction problem over Z/9 starting from
# the specific modulo-9 radial lift, with total degree capped at 81.
RADIAL_BRANCH_CERTIFICATE = {
    (41, 13): 6,
    (43, 14): 6,
    (47, 16): 6,
    (49, 17): 6,
    (51, 13): 3,
    (53, 19): 3,
    (55, 20): 3,
    (60, 15): 6,
    (62, 16): 6,
    (66, 18): 3,
    (69, 22): 3,
    (72, 21): 6,
    (76, 18): 3,
    (78, 19): 6,
    (80, 20): 2,
    (82, 21): 3,
    (84, 22): 6,
    (89, 22): 2,
    (91, 23): 2,
    (114, 27): 1,
}

# Second digits added to the coefficient representatives 0,...,8
# of the constructed mod-9 pair. Keys are exponents of x and y.
SECOND_DIGIT_P = {
    (3, 1): 2,
    (8, 1): 1,
    (10, 2): 2,
    (11, 0): 2,
    (13, 1): 2,
    (14, 4): 1,
    (15, 2): 1,
    (16, 0): 2,
    (18, 1): 2,
    (19, 4): 2,
    (20, 2): 1,
    (21, 5): 1,
    (23, 6): 1,
    (25, 2): 1,
    (26, 0): 1,
    (26, 5): 2,
    (31, 5): 2,
    (32, 3): 1,
    (32, 8): 1,
    (35, 7): 1,
    (37, 8): 1,
    (41, 10): 1,
    (44, 9): 1,
    (50, 12): 1,
}

SECOND_DIGIT_Q = {
    (2, 2): 2,
    (7, 2): 1,
    (11, 4): 1,
    (17, 2): 2,
    (20, 1): 2,
    (29, 8): 1,
    (32, 2): 2,
    (32, 7): 2,
    (35, 1): 2,
    (38, 10): 1,
    (41, 4): 1,
    (44, 8): 2,
    (47, 7): 2,
    (53, 10): 2,
    (56, 14): 2,
    (59, 13): 1,
    (65, 16): 2,
}


# Exact certificates at modulus 81. The cap-108 functional tests all
# later digits and all capped Frobenius shifts of the first lift.
FROBENIUS_CAP108_DUAL = {
    (76, 18): 3, (78, 19): 3, (87, 21): 3, (91, 23): 6,
    (93, 24): 6, (94, 27): 3, (96, 28): 3, (97, 26): 3,
    (99, 27): 3, (103, 24): 3, (103, 29): 6, (109, 27): 3,
    (114, 27): 1, (118, 29): 2, (145, 35): 2,
}

# These digits are multiplied by 9 and added to the radial first lift.
DEGREE109_WITNESS_P = {
    (1, 0): 7, (5, 2): 6, (7, 3): 3, (8, 1): 1,
    (10, 2): 2, (11, 0): 7, (13, 1): 6, (14, 4): 1,
    (16, 0): 3, (16, 5): 3, (17, 3): 8, (20, 2): 7,
    (22, 3): 3, (23, 1): 3, (23, 6): 1, (25, 2): 1,
    (25, 7): 3, (26, 0): 1, (26, 5): 5, (28, 1): 6,
    (31, 5): 5, (32, 3): 1, (32, 8): 4, (34, 4): 3,
    (35, 2): 3, (35, 7): 4, (37, 8): 1, (38, 1): 3,
    (38, 6): 6, (39, 9): 2, (41, 10): 1, (43, 6): 6,
    (43, 11): 3, (44, 4): 6, (44, 9): 4, (46, 0): 6,
    (46, 10): 6, (47, 8): 6, (49, 9): 6, (50, 7): 3,
    (50, 12): 7, (52, 3): 6, (52, 13): 3, (53, 11): 6,
    (55, 12): 3, (56, 10): 6, (59, 14): 6, (61, 15): 3,
    (62, 13): 3, (64, 9): 3, (68, 16): 6, (70, 12): 3,
    (82, 18): 6, (88, 21): 6,
}

DEGREE109_WITNESS_Q = {
    (0, 1): 2, (2, 2): 2, (7, 2): 1, (8, 0): 6,
    (10, 1): 8, (11, 4): 1, (12, 2): 8, (13, 0): 4,
    (16, 4): 5, (17, 2): 5, (18, 5): 2, (19, 3): 3,
    (20, 1): 2, (23, 5): 3, (26, 4): 6, (28, 5): 6,
    (29, 8): 7, (32, 2): 2, (32, 7): 8, (35, 1): 8,
    (38, 5): 3, (38, 10): 1, (41, 4): 4, (44, 8): 5,
    (47, 7): 8, (50, 1): 6, (53, 5): 3, (53, 10): 2,
    (56, 14): 5, (59, 8): 6, (59, 13): 7, (62, 7): 6,
    (65, 11): 3, (65, 16): 8, (71, 14): 3, (77, 17): 6,
    (80, 16): 6, (83, 20): 3, (86, 19): 3,
}


def evaluate(polynomial, weights, modulus):
    return sum(
        weight * polynomial.get(exponent, 0)
        for exponent, weight in weights.items()
    ) % modulus


def check_frobenius_obstruction(lift_p, lift_q, error, cap, weights):
    assert evaluate(error, weights, 9) == 3
    for i in range(cap + 1):
        for j in range(cap + 1 - i):
            monomial = {(i, j): 1}
            left = jac(monomial, lift_q)
            right = jac(lift_p, monomial)
            assert evaluate(left, weights, 9) == 0
            assert evaluate(right, weights, 9) == 0
            if i % 3 == j % 3 == 0:
                assert evaluate(quotient(left, 3), weights, 9) == 0
                assert evaluate(quotient(right, 3), weights, 9) == 0


def main():
    _, seed_p, seed_q = family(3)
    error = mod(quotient(sub(jac(seed_p, seed_q), ONE), 3), 3)
    assert evaluate(error, CERTIFICATE, 3) == 1

    residue_p = mod(seed_p, 3)
    residue_q = mod(seed_q, 3)
    for i in range(81):
        for j in range(81 - i):
            monomial = {(i, j): 1}
            assert evaluate(jac(monomial, residue_q), CERTIFICATE, 3) == 0
            assert evaluate(jac(residue_p, monomial), CERTIFICATE, 3) == 0

    digit = mod({
        exponent: -coefficient * pow(sum(exponent) + 2, -1, 3)
        for exponent, coefficient in error.items()
    }, 3)
    multiplier = add(ONE, scale(digit, 3))
    lift_p = radial_substitute(seed_p, multiplier, 9)
    lift_q = radial_substitute(seed_q, multiplier, 9)
    assert mod(lift_p, 3) == residue_p and mod(lift_q, 3) == residue_q
    assert mod(jac(lift_p, lift_q), 9) == ONE
    assert (degree(lift_p), degree(lift_q)) == (51, 81)
    assert lift_q[(65, 16)] == 6

    lift_p_27 = add(lift_p, scale(SECOND_DIGIT_P, 9))
    lift_q_27 = add(lift_q, scale(SECOND_DIGIT_Q, 9))
    assert mod(lift_p_27, 9) == lift_p and mod(lift_q_27, 9) == lift_q
    assert mod(jac(lift_p_27, lift_q_27), 27) == ONE
    assert (degree(lift_p_27), degree(lift_q_27)) == (62, 81)
    assert not any(i == 0 for i, j in lift_p_27)
    assert {(i, j): c for (i, j), c in lift_q_27.items() if i == 0} == {(0, 1): 1}

    error_81 = mod(quotient(sub(jac(lift_p, lift_q), ONE), 9), 9)
    check_frobenius_obstruction(lift_p, lift_q, error_81, 81, RADIAL_BRANCH_CERTIFICATE)
    check_frobenius_obstruction(lift_p, lift_q, error_81, 108, FROBENIUS_CAP108_DUAL)

    lift_p_81 = add(lift_p, scale(DEGREE109_WITNESS_P, 9))
    lift_q_81 = add(lift_q, scale(DEGREE109_WITNESS_Q, 9))
    assert mod(lift_p_81, 9) == lift_p and mod(lift_q_81, 9) == lift_q
    assert mod(jac(lift_p_81, lift_q_81), 81) == ONE
    assert (degree(lift_p_81), degree(lift_q_81)) == (109, 105)

    print("Mod 9: cap 80 obstructed on all 6642 correction monomials.")
    print("Mod 9: witness degrees (51, 81), reduction and determinant verified.")
    print("Mod 27: witness degrees (62, 81), reduction and determinant verified.")
    print("Mod 81: radial and all 3^812 Frobenius-shifted mod-9 choices obstructed at cap 81.")
    print("Mod 81: cap 108 fails for all Frobenius shifts (11990 regular and 1406 divided columns).")
    print("Mod 81: witness degrees (109, 105), radial first lift and determinant verified.")


if __name__ == "__main__":
    main()
