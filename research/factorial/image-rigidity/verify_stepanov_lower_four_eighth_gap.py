#!/usr/bin/env python3
"""Exact algebra and bounded coefficient controls for the a=4,d=8 proof."""
from fractions import Fraction
from math import comb, factorial
from pathlib import Path
import json


def mul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def valuation(x, p=2):
    if not x:
        return None
    x = Fraction(x)
    a, b = abs(x.numerator), x.denominator
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    while b % p == 0:
        b //= p
        v -= 1
    return v


def upper_row(n, c5, c6, c7):
    c, C = 4*n+8, 12*n+9
    ell = []
    for j in range(7):
        f = [8*n, 1]
        for r in range(j):
            f = mul(f, [c-r, -1])
        ell.append(f)
    b6 = c7
    b5 = Fraction(c6-b6*ell[6][6], ell[5][6])
    b4 = Fraction(c5-b6*ell[6][5]-b5*ell[5][5], ell[4][5])
    assert b5 == -c6-(16*n+33)*c7
    assert b4 == c5+(12*n+30)*c6+(144*n*n+480*n+545)*c7
    result = b4+(C+3)*b5+(C+2)*(C+3)*b6
    expected = c5+18*c6+(96*n*n+168*n+281)*c7
    assert result == expected
    return result


def main():
    identity_values = []
    for n in [0, 1, 2]:
        row = [upper_row(n, *basis) for basis in [(1, 0, 0), (0, 1, 0), (0, 0, 1)]]
        assert all(x.denominator == 1 for x in row)
        identity_values.append({"n": n, "row": [int(x) for x in row]})
    # The leading-three coefficient expansion has degree at most two
    # in n, so the preceding values are a complete identity certificate.
    # A corrupted constant is rejected by that same certificate.
    assert any(item["row"][2] != 96*item["n"]**2+168*item["n"]+280
               for item in identity_values)
    controls, boundary = 0, []
    for n in range(1, 129):
        C = 12*n+9
        original = Fraction(8*(n+1), C*(C+1)*(C+2)*(C+3))
        simplified = Fraction(2, 3*(12*n+9)*(12*n+10)*(12*n+11))
        assert original == simplified and valuation(original) == 0
        shape = 4*n+1
        rising = 1
        for k in range(1, n+1):
            for j in range(8*(k-1), 8*k):
                rising *= shape+j
            v = valuation(comb(n, k)*rising)
            assert v >= 7*k
            if k >= 2:
                assert v > 7*k
            if k == 1 and n % 4 == 1:
                assert v == 7
                boundary.append(n)
            controls += 1
    boundary_coefficient = factorial(12)//factorial(4)
    assert valuation(boundary_coefficient) == 7
    report = {
        "status": "PASS",
        "scope": "complete bounded-degree identity certificate plus finite coefficient controls",
        "symbolic_identity_degree_bound_in_n": 2,
        "symbolic_identity_values": identity_values,
        "corrupted_upper_row_constant_rejected": True,
        "kappa_normalization_controls": 128,
        "actual_factorial_coefficient_valuation_controls": controls,
        "radius_seven_boundary_indices_checked": boundary,
        "n1_exact_reciprocal_root_u_equals_minus_one_over_t": boundary_coefficient,
        "n1_exact_reciprocal_root_valuation": 7,
        "universal_proof": "STEPANOV_LOWER_FOUR_EIGHTH_GAP_CLOSURE.md, Sections 1-2",
        "imported_generic_stage_lemma": "STEPANOV_LOWER_THREE_EIGHTH_GAP_CLOSURE.md, Section 2",
        "all_order_theorem_not_inferred_from_a_finite_prefix": True,
    }
    Path(__file__).with_name("stepanov_lower_four_eighth_gap_verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O or -OO: exact assertions must execute.")
    main()
