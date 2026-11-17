#!/usr/bin/env python3
"""Complete exact depth-one replay through 100000; Python standard library."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import platform

HERE = Path(__file__).resolve().parent
SOURCE = "ELLIPTIC_EXCEPTION_BOUNDED_FATAL_PRIME_CERTIFICATE.md"
CERTIFICATE = "elliptic_exception_bounded_fatal_prime_certificate.json"
A, B = -10758382, 6600759744
LONG_A, LONG_B = -10656, 27091730


def primes_to(bound):
    sieve = bytearray(b"\x01")*(bound+1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(bound)+1):
        if sieve[p]:
            sieve[p*p::p] = b"\x00"*(((bound-p*p)//p)+1)
    return [p for p in range(2, bound+1) if sieve[p]]


def legendre(a, p):
    a %= p
    return 0 if not a else 1 if pow(a, (p-1)//2, p) == 1 else -1


def eligible(p):
    return (p >= 11 and p % 3 == 2 and p not in (1889, 2709173)
            and legendre(51730, p) == 1
            and (legendre(2, p) == 1 or legendre(13545865, p) == 1))


def cardinality(p):
    """Complete affine x scan, plus the point at infinity."""
    chi = [-1]*p
    chi[0] = 0
    for x in range(1, (p+1)//2):
        chi[x*x % p] = 1
    return p+1+sum(chi[x*(x*(x+LONG_A)+LONG_B) % p] for x in range(p))


def modq(q, modulus):
    return q.numerator*pow(q.denominator, -1, modulus) % modulus


def reduce_point(point, p):
    if point is None:
        return ("f", 0)
    x, y = point
    modulus = p*p
    if x.denominator % p:
        return ("a", modq(x, modulus), modq(y, modulus))
    return ("f", modq(-x/y, modulus))


def add(point1, point2, p, branches=None):
    """Group law on E(Q_p)/E_2, in affine and formal charts."""
    modulus = p*p
    if point1[0] == "f" and point2[0] == "f":
        if branches is not None:
            branches["formal plus formal"] += 1
        return ("f", (point1[1]+point2[1]) % modulus)
    if point1[0] == "f":
        point1, point2 = point2, point1
    if point2[0] == "f":
        if branches is not None:
            branches["affine plus formal"] += 1
        x, y = point1[1:]
        t = point2[1]
        return ("a", (x+2*y*t) % modulus, (y+(3*x*x+A)*t) % modulus)
    x, y = point1[1:]
    u, v = point2[1:]
    if (u-x) % p:
        if branches is not None:
            branches["unit x difference"] += 1
        slope = (v-y)*pow((u-x) % modulus, -1, modulus) % modulus
    elif (y+v) % p:
        if branches is not None:
            branches["same non-two-torsion reduction"] += 1
        slope = (x*x+x*u+u*u+A)*pow((y+v) % modulus, -1, modulus) % modulus
    elif (v-y) % p:
        if branches is not None:
            branches["opposite non-two-torsion reductions"] += 1
        return ("f", (u-x)*pow((v-y) % modulus, -1, modulus) % modulus)
    else:
        if branches is not None:
            branches["two-torsion reduction"] += 1
        return ("f", (y+v)*pow((3*x*x+A) % modulus, -1, modulus) % modulus)
    xx = (slope*slope-x-u) % modulus
    yy = (-y+slope*(x-xx)) % modulus
    return ("a", xx, yy)


def multiply(point, n, p):
    result = ("f", 0)
    while n:
        if n & 1:
            result = add(result, point, p)
        point = add(point, point, p)
        n //= 2
    return result


def rational_add(point1, point2):
    if point1 is None:
        return point2
    if point2 is None:
        return point1
    x, y = point1
    u, v = point2
    if x == u and y == -v:
        return None
    slope = (3*x*x+A)/(2*y) if point1 == point2 else (v-y)/(u-x)
    xx = slope*slope-x-u
    return xx, -y+slope*(x-xx)


def exact_rational_controls(points):
    samples = [None]
    for generator in points.values():
        multiples = [None]
        for n in range(1, 9):
            multiples.append(rational_add(multiples[-1], generator))
        samples += multiples[1:]+[(x, -y) for x, y in multiples[1:]]
    branches = Counter()
    cases = 0
    for p in [11, 17, 23, 29, 41, 53, 101]:
        for point1 in samples:
            for point2 in samples:
                expected = reduce_point(rational_add(point1, point2), p)
                actual = add(reduce_point(point1, p), reduce_point(point2, p), p, branches)
                if actual != expected:
                    raise AssertionError(("rational group-law crosscheck", p, point1, point2))
                cases += 1
    return cases, dict(branches)


def replay():
    cert = json.loads((HERE/CERTIFICATE).read_text())
    checks = []

    def check(label, value):
        if not value:
            raise AssertionError(label)
        checks.append(label)

    check("inclusive bound", cert["bound_inclusive"] == 100000)
    check("literal curve", cert["curve_short"] == {"A": A, "B": B})
    check("integral model translation",
          A == LONG_B-LONG_A**2//3 and B == 2*LONG_A**3//27-LONG_A*LONG_B//3)
    points = {name: tuple(F(q) for q in point) for name, point in cert["points_short"].items()}
    check("literal points", points == {"P": (F(199200, 49), F(59459616, 343)),
                                       "Q": (F(1563680, 529), F(304670016, 12167))})
    for name, (x, y) in points.items():
        check("exact rational curve point "+name, y*y == x*x*x+A*x+B)
    cases, branches = exact_rational_controls(points)
    check("complete independent rational addition domain", cases == 7623)
    check("all six modular chart branches exercised", len(branches) == 6 and min(branches.values()) > 0)
    primes = primes_to(cert["bound_inclusive"])
    check("complete prime sieve count", len(primes) == 9592)
    active = [p for p in primes if p >= 11 and p % 3 == 2 and p not in (1889, 2709173)]
    split = [p for p in active if legendre(51730, p) == 1]
    selected = [p for p in primes if eligible(p)]
    check("active-prime domain count", len(active) == 4804)
    check("branch-split domain count", len(split) == 2336)
    check("eligible-union domain count", len(selected) == cert["eligible_prime_count"] == 1722)
    rows = cert["rows"]
    check("literal row schema", cert["columns"] == ["prime", "curve_order",
          "P_formal_parameter_div_p_mod_p", "Q_formal_parameter_div_p_mod_p"])
    check("complete selected-prime coverage", [row[0] for row in rows] == selected)
    check("common witness label", cert["witness_for_every_row"] == "Q")
    p_nonwitness = []
    total_counted_x = 0
    for row in rows:
        p, n, up, uq = row
        actual_n = cardinality(p)
        total_counted_x += p
        if actual_n != n:
            raise AssertionError(("complete point count", p, n, actual_n))
        if n % 2 or n % p == 0 or (n-p-1)**2 > 4*p:
            raise AssertionError(("point-count consistency", p, n))
        expected_units = {"P": up, "Q": uq}
        for name, point in points.items():
            result = multiply(reduce_point(point, p), n, p)
            if result[0] != "f" or result[1] % p:
                raise AssertionError(("multiple is formal", p, name, result))
            if result[1]//p != expected_units[name]:
                raise AssertionError(("formal coefficient", p, name, result))
        if not up:
            p_nonwitness.append(p)
        if not 1 <= uq < p:
            raise AssertionError(("Q exact depth one", p, uq))
    check("complete point-count and formal-row replay", len(rows) == 1722)
    check("P higher-depth row list", p_nonwitness == [401])
    check("401 P/Q contrast", next(row for row in rows if row[0] == 401) == [401, 408, 0, 288])
    check("all rows witnessed without a generator basis", all(0 < row[3] < row[0] for row in rows))
    return {"passed": True, "named_controls": len(checks), "checks": checks,
            "bound_inclusive": cert["bound_inclusive"], "primes_in_domain": len(primes),
            "active_good_noncore_primes": len(active), "branch_split_primes": len(split),
            "eligible_union_primes": len(selected), "complete_point_counts": len(rows),
            "affine_x_values_counted": total_counted_x, "modular_scalar_multiplications": 2*len(rows),
            "Q_depth_one_witnesses": len(rows), "P_higher_depth_rows": p_nonwitness,
            "exact_rational_addition_controls": cases, "addition_chart_case_counts": branches,
            "no_simultaneous_higher_depth_rows_in_domain": True,
            "full_Mordell_Weil_basis_assumed": False, "claim_beyond_bound": False,
            "new_all_places_Hasse_seed_claimed": False, "finite_checks_prove_unbounded_depth_claim": False,
            "formal_verification": False, "external_review": False, "independent_replay": False,
            "python": platform.python_version(), "arithmetic": "standard-library exact integers and fractions"}


if __name__ == "__main__":
    report = replay()
    names = [SOURCE, CERTIFICATE, Path(__file__).name]
    report["input_sha256"] = {name: sha256((HERE/name).read_bytes()).hexdigest() for name in names}
    (HERE/"elliptic_exception_bounded_fatal_prime_checks.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "checks"}, indent=2))
