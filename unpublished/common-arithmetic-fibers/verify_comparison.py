#!/usr/bin/env python3
"""Check the Ulam coefficient specialization in exact rational arithmetic.

This verifies the displayed coefficient formula for arbitrary supplied
rational coefficients and selected controls, not either full paper proof.
"""
from fractions import Fraction as Q


def ulam_coefficients(c):
    assert len(c) >= 4 and c[1] and c[-1]
    c = list(map(Q, c))
    eta = {1: 3*(c[3]/c[1]-1)}
    eta.update({k-2: k*c[k]/c[1] for k in range(4, len(c))})
    target_q, target_r = -2*c[2]/c[1], -2*c[0]/c[1]
    omega = [-target_r, Q(2), -target_q, Q(2)] + [Q(0)]*(len(c)-4)
    for j, value in eta.items():
        omega[j+2] += Q(2, j+2)*value
    return omega


def main():
    checked = 0
    for n in range(3, 25):
        for sign in (-1, 1):
            c = [Q(sign*(k*k-3*k+7), k+1) for k in range(n+1)]
            expected = [2*x/c[1] for x in c]
            assert ulam_coefficients(c) == expected
            bad = expected.copy()
            bad[3] += 1
            assert ulam_coefficients(c) != bad
            checked += 1
    # The formula remains valid when the translated cubic coefficient is zero.
    c = [Q(7), Q(3), Q(-2), Q(0), Q(5)]
    assert ulam_coefficients(c) == [2*x/c[1] for x in c]
    print(f'PASS: {checked+1} coefficient controls, including c3=0; altered coefficients rejected.')


if __name__ == '__main__':
    main()
