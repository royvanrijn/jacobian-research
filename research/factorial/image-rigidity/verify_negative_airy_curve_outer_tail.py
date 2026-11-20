#!/usr/bin/env python3
"""Replay exact elementary certificates for the negative Airy outer tail.

This checks algebra and rational guards, not a numerical Airy quadrature.
The joint Airy lower bounds are mathematical imports identified in the note.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

if not __debug__:
    raise SystemExit("optimized execution is refused")

HERE = Path(__file__).resolve().parent


def trim(p):
    p = list(map(Q, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    out = [Q(0)] * max(len(a), len(b))
    for i, v in enumerate(a): out[i] += v
    for i, v in enumerate(b): out[i] += v
    return trim(out)


def scale(a, q):
    return trim([Q(q) * v for v in a])


def mul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return trim(out)


def power(a, n):
    out = [Q(1)]
    for _ in range(n): out = mul(out, a)
    return out


def shift(a, q):
    out = [Q(0)]
    for j, v in enumerate(a):
        out = add(out, scale(power([Q(q), Q(1)], j), v))
    return out


def guarded_transform(a, bound, degree):
    # (1+z)^degree * a(bound*z/(1+z)).
    out = [Q(0)]
    for j, v in enumerate(a):
        term = [Q(0)] * j + scale(power([Q(1), Q(1)], degree-j), v*Q(bound)**j)
        out = add(out, term)
    return out


def exact_certificates():
    y = [Q(0), Q(1)]
    ym1 = [Q(-1), Q(1)]
    B = [Q(2), Q(-3), Q(0), Q(4)]
    C = [Q(1), Q(2), Q(3), Q(8), Q(4)]
    P = [Q(20), Q(-1), Q(-16), Q(52), Q(64), Q(16)]
    # Clearing 81*y^4 from T1.
    first = scale(mul(mul(power(ym1, 4), power([1, 2], 4)), mul([2, 1], B)), 4)
    left = add(first, scale(power(add(B, scale(y, -3)), 4), -1))
    right = scale(mul(mul(y, power(ym1, 4)), P), 12)
    assert left == right, "cubic phase-loss identity"
    F = add(scale(mul([2, 1], power([-1, 0, 0, 4], 2)), 4), scale(mul(C, B), -1))
    assert F == list(map(Q, [6, 3, 0, -75, -24, 0, 96, 48]))
    right_coefficients = shift(F, 1)
    assert right_coefficients == list(map(Q, [54, 594, 2079, 3429, 3096, 1584, 432, 48]))
    left_coefficients = guarded_transform(F, Q(3, 8), 7)
    expected = [Q(6), Q(345,8), Q(531,4), Q(114135,512), Q(110697,512), Q(60273,512), Q(260283,8192), Q(394833,131072)]
    assert left_coefficients == expected
    assert all(c > 0 for c in right_coefficients + left_coefficients)
    # x'^2+1 and the exact radial derivative simplify to C*B/denominator.
    assert add(power([1, 1, 1], 2), scale(mul(power(y, 3), [2, 1]), 3)) == C
    # Exact imaginary phase: 3*y*x^2 - y^3+3*y = 2, after clearing 3*y.
    x2_numerator = mul(power(ym1, 2), [2, 1])
    assert add(x2_numerator, add(scale(power(y, 3), -1), scale(y, 3))) == [Q(2)]
    return {"phase_loss_identity": True, "arc_ratio_identity": True,
            "positive_right_coefficients": [str(x) for x in right_coefficients],
            "positive_left_coefficients": [str(x) for x in left_coefficients]}


def rational_guards():
    rho = Q(3, 8)
    guards = {
        "left_interval_value": Q(4,3)*Q(3,8)**2 - 1 + Q(2,3)/Q(3,8) == Q(139,144),
        "left_interval_below_one": Q(139,144) < 1,
        "P_left_positive": 20-16-1 >= 3,
        "P_right_positive": 52-16-1 >= 35,
        "ratio_at_nine_eighths": Q(1,3)*(1-Q(8,9))**2 == Q(1,243),
        "rational_cubic_loss": 365**2 > 2*243**2,
        "positive_absolute_loss": rho**3/Q(365) > Q(1,7000),
        "compact_tail_coefficient": 4*450*36000*365/(3*rho**2) == 56064000000,
        "compact_coefficient_below_limit": 56064000000 < 10**11,
        "large_tail_coefficient": Q(64*365, 9*rho**2) < 19000,
        "zero_zeta_coefficient": 72000/rho**2 == 512000,
        "zero_zeta_decay": rho**3/3 > Q(1,14000),
    }
    # An actual upper bound, not just a truncated-series lower estimate:
    # e^(16/3), terms through 40; later term ratios <=(16/3)/42<1.
    x = Q(16,3)
    terms = [x**j/Q(__import__('math').factorial(j)) for j in range(42)]
    guards["exp_sixteen_thirds_upper"] = sum(terms[:41]) + terms[41]/(1-x/Q(42)) < 225
    assert all(guards.values()), [k for k, v in guards.items() if not v]
    return guards


def negative_controls():
    # The P coefficients and the arc bound are materially used.
    p_bad = [Q(20),Q(-1),Q(-16),Q(52),Q(64),Q(15)]
    y = [Q(0),Q(1)]
    ym1 = [-1,1]
    B = [2,-3,0,4]
    first = scale(mul(mul(power(ym1,4),power([1,2],4)),mul([2,1],B)),4)
    left = add(first,scale(power(add(B,scale(y,-3)),4),-1))
    assert left != scale(mul(mul(y,power(ym1,4)),p_bad),12)
    # At y=1 the squared ratio is 2, so a claimed bound 1/2 is false.
    C1, B1, denominator1 = Q(18), Q(3), Q(3)*Q(3)**2
    assert C1*B1 > Q(1,4)*denominator1
    assert Q(27,186880) < Q(1,6000), "a stronger asserted absolute loss is rejected"
    return 3


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "negative_airy_curve_outer_tail_verification.json")
    args = parser.parse_args()
    result = {"schema_version": 1,
              "scope": "exact elementary negative-sector canonical Airy outer-tail certificates; not a full-period remainder",
              "proof_sha256": hashlib.sha256((HERE / "NEGATIVE_AIRY_CURVE_OUTER_TAIL.md").read_bytes()).hexdigest(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "certificates": exact_certificates(),
              "rational_guards": rational_guards(),
              "negative_controls_rejected": negative_controls(),
              "phase_loss_constant": "1/243", "physical_loss_constant": "1/7000",
              "outer_tail_bound": "10^11*m^(-1/3)*exp(-m/14000)*joint_Airy_norm",
              "formal_verification": False, "external_review": False, "independent_replay": False}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"passed": True, "rational_guards": len(result["rational_guards"]),
                      "negative_controls_rejected": 3, "output": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
