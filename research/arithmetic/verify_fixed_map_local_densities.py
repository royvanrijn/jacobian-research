#!/usr/bin/env python3
"""Finite exact controls for FIXED_MAP_LOCAL_DENSITIES.md; not an asymptotic proof."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform

ROOT = Path(__file__).resolve().parents[2]
EXPECTED = {
    10: (11, 2, 2, 0),
    32: (104, 4, 4, 0),
    100: (1038, 16, 10, 6),
    300: (9121, 89, 22, 67),
}


def icbrt(n: int) -> int:
    assert n >= 0
    lo, hi = 0, 1
    while hi**3 <= n:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid
    return lo


def valuations(n: int) -> dict[int, int]:
    assert n > 0
    answer = {}
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            answer[p] = e
        p += 1
    if n > 1:
        answer[n] = 1
    return answer


def local_integer_info(n: int) -> tuple[bool, int, int]:
    facts = valuations(n)
    inert_ok = all(e % 3 == 0 for p, e in facts.items() if p % 3 == 2)
    e3 = facts.get(3, 0)
    return inert_ok, e3, n // 3**e3


def local_condition(m: int, n: int, info: dict | None = None) -> bool:
    if m == 0:
        return False
    gm, em, um = info[abs(m)] if info else local_integer_info(abs(m))
    gn, en, un = info[n] if info else local_integer_info(n)
    if m < 0:
        um = -um
    return gm and gn and (em - en) % 3 == 0 and (um * pow(un, -1, 9)) % 9 in (1, 8)


def projective_height(m: int, n: int) -> int:
    assert n > 0 and math.gcd(m, n) == 1
    coordinates = (9*n, -9*n, 32*m, 24*m + 3*n)
    content = math.gcd(*coordinates)
    assert content == math.gcd(n, 8) * math.gcd(m, 3)
    return max(map(abs, coordinates)) // content


def cube_count(bound: int) -> dict:
    """Enumerate primitive roots r/s, separately from the (m,n) direct loop."""
    counts = [0, 0]
    max_r, max_s = icbrt(3*bound // 4), icbrt(8*bound // 3)
    for s in range(1, max_s + 1):
        for r0 in range(1, max_r + 1):
            if math.gcd(r0, s) != 1:
                continue
            for r in (-r0, r0):
                if projective_height(r**3, s**3) <= bound:
                    counts[r > 0] += 1
    return {"height_bound": bound, "negative": counts[0], "positive": counts[1],
            "total": sum(counts), "root_box": [max_r, max_s]}


def direct_count(bound: int) -> dict:
    max_m, max_n = 3*bound // 4, 8*bound // 3
    info = {v: local_integer_info(v) for v in range(1, max(max_m, max_n) + 1)}
    cubes = {r**3 for r in range(icbrt(max(max_m, max_n)) + 1)}
    totals = [0, 0, 0, 0]
    local_sign, failure_sign = [0, 0], [0, 0]
    for n in range(1, max_n + 1):
        for m in range(-max_m, max_m + 1):
            if math.gcd(m, n) != 1 or projective_height(m, n) > bound:
                continue
            globally_soluble = m != 0 and abs(m) in cubes and n in cubes
            locally_soluble = local_condition(m, n, info)
            failure = locally_soluble and not globally_soluble
            assert not globally_soluble or locally_soluble
            totals[0] += 1
            totals[1] += locally_soluble
            totals[2] += globally_soluble
            totals[3] += failure
            if locally_soluble:
                local_sign[m > 0] += 1
            if failure:
                failure_sign[m > 0] += 1
    assert tuple(totals) == EXPECTED[bound], (bound, totals)
    assert totals[1] == totals[2] + totals[3]
    root_count = cube_count(bound)
    assert root_count["total"] == totals[2]
    return {"height_bound": bound, "all": totals[0], "local": totals[1],
            "global": totals[2], "hasse_failure": totals[3],
            "local_negative_positive": local_sign,
            "failure_negative_positive": failure_sign,
            "parameter_box": [max_m, max_n]}


def poly_mul(a: list, b: list) -> list:
    result = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def algebra_controls() -> dict:
    quadratic = [Fraction(1), Fraction(1), Fraction(1)]
    for c in (Fraction(-3, 2), Fraction(-1), Fraction(1, 3), Fraction(1), Fraction(2)):
        cubic = [-c**3, 0, 0, 1]
        assert poly_mul([-c, 1], [c*c, c, 1]) == cubic
        assert poly_mul(cubic, quadratic) == poly_mul([-c, 1], poly_mul([c*c, c, 1], quadratic))
    assert poly_mul([-1, 0, 0, 1], quadratic) == poly_mul([-1, 1], poly_mul(quadratic, quadratic))
    assert poly_mul([0, 0, 0, 1], quadratic) == [0, 0, 0, 1, 1, 1]
    for a, expected in [(0, False), (1, True), (-1, True), (2, False), (3, False),
                        (7, False), (19, True), (91, True), (152, True),
                        (Fraction(19, 37), True), (64, True)]:
        a = Fraction(a)
        assert local_condition(a.numerator, a.denominator) == expected
    return {"cubic_factorizations": 5, "quintic_factorizations": 5,
            "exceptional_factorizations": 2, "local_parameter_controls": 11}


def local_weights() -> dict:
    d2 = Fraction(1, 2) + Fraction(4, 8) + Fraction(16, 16) + Fraction(64, 16)
    d3 = Fraction(2, 3) + 9*Fraction(2, 9)
    assert d2 == 6 and d3 == Fraction(8, 3)
    residue_counts = {(e, f): 0 for e in (0, 1) for f in (0, 1)}
    for r in range(6):
        for s in range(6):
            if math.gcd(math.gcd(r, s), 6) == 1:
                residue_counts[int(s % 2 == 0), int(r % 3 == 0)] += 1
    assert residue_counts == {(0, 0): 12, (0, 1): 4, (1, 0): 6, (1, 1): 2}
    # Since 8^(2/3)=4, retain the cube weight exactly as A+B*3^(2/3).
    weights = [Fraction(0), Fraction(0)]
    for (e, f), count in residue_counts.items():
        weights[f] += Fraction(count * 4**e, 36)
    assert weights == [Fraction(1), Fraction(1, 3)]
    assert Fraction(11, 1728) / Fraction(143, 864) == Fraction(1, 26)
    assert Fraction(71, 20736) + Fraction(72, 20736) == Fraction(143, 20736)
    return {"ambient_D2": str(d2), "ambient_D3": str(d3),
            "primitive_cube_residue_counts_00_01_10_11": [12, 4, 6, 2],
            "cube_weight_1_and_3_to_2_over_3": list(map(str, weights))}


def numerical_constants() -> dict:
    # A numerical control only, not a certified quadrature or proof of (5).
    count = 10000
    lo, hi = 36**(-1/3), 32**(-1/3)
    h = (hi - lo) / count
    def f(u):
        return ((1 - 24*u**3) / 3)**(1/3)
    integral = h/3*(f(lo) + f(hi)
                   + 4*sum(f(lo+i*h) for i in range(1, count, 2))
                   + 2*sum(f(lo+i*h) for i in range(2, count, 2)))
    negative = 288**(-1/3)
    positive = 324**(-1/3) + integral
    area = negative + positive
    zeta2 = math.pi**2/6
    return {"area_negative": negative, "area_positive": positive,
            "area_total": area, "global_count_constant": area*(3+3**(2/3))/(2*zeta2),
            "all_target_constant": 143/(864*zeta2),
            "scope": "uncertified Simpson control; not used for exact tests"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    algebra_count = algebra_controls()
    weights = local_weights()
    direct = [direct_count(bound) for bound in EXPECTED]
    cubes = [cube_count(bound) for bound in (1000, 1000000)]
    assert [row["total"] for row in cubes] == [44, 4693]
    source_paths = ["research/arithmetic/FIXED_MAP_LOCAL_DENSITIES.md",
                    "research/arithmetic/verify_fixed_map_local_densities.py",
                    "unpublished/fixed-map-hasse-failures/sections/02-fixed-map.tex",
                    "unpublished/fixed-map-hasse-failures/sections/03-local-global.tex",
                    "unpublished/fixed-map-hasse-failures/sections/04-rational-counting.tex"]
    receipt = {
        "status": "passed", "python": platform.python_version(),
        "command": "python3 research/arithmetic/verify_fixed_map_local_densities.py --output research/arithmetic/fixed_map_local_densities.verification.json",
        "scope": "Internal finite exact controls only; no asymptotic, formal, independent-replay, or external-review certificate.",
        "algebra_and_local_controls": algebra_count, "local_weights": weights,
        "direct_target_enumeration": direct, "primitive_cube_enumeration": cubes,
        "numerical_controls": numerical_constants(),
        "source_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths},
    }
    if args.output:
        args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
