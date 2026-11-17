#!/usr/bin/env python3
"""Exact arithmetic controls for the finite dependency core; standard library."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import prod
from pathlib import Path
import json
import platform
import runpy

HERE = Path(__file__).resolve().parent
SOURCE = "ELLIPTIC_EXCEPTION_DEPENDENCY_CLOSED_CORE.md"
BASELINE = "verify_elliptic_exception_finite_obstruction_baseline.py"
BASELINE_INPUT = "elliptic_exception_finite_obstruction_baseline.json"
ARITHMETIC = "verify_elliptic_exception_one_prime_positivity.py"


def replay():
    baseline = runpy.run_path(str(HERE / BASELINE), run_name="core_baseline")
    prior = baseline["replay"]()
    dep = runpy.run_path(str(HERE / ARITHMETIC), run_name="core_arithmetic")
    cert = json.loads((HERE / BASELINE_INPUT).read_text())
    a, b = cert["curve"]["a2"], cert["curve"]["a4"]
    checks = []

    def check(label, value):
        if not value:
            raise AssertionError(label)
        checks.append(label)

    check("Q/U baseline exact replay", prior["passed"])
    check("curve coefficients", (a, b) == (-10656, 27091730))
    q0 = 2709173
    core = {2, 3, 5, 1889}
    dependencies = core | {7}
    for p in sorted(dependencies | {q0}):
        check("dependency or ramified-class prime " + str(p), dep["prime"](p))
    disc = 16 * b * b * (a * a - 4 * b)
    check("complete curve discriminant factorization",
          disc == 2**9 * 5**2 * 7**3 * 1889 * q0**2)
    check("bad-prime list complete", set([2, 5, 7, 1889, q0]) ==
          {p for p in [2, 3, 5, 7, 1889, q0] if disc % p == 0})
    check("q0 automatic class-two residue", q0 % 8 == 5)
    check("q0 automatic class-five-q valuation", dep["valuation"](5*q0, q0) == 1)
    check("seven automatic residue", 7 % 3 == 1)
    check("five node b valuation", dep["valuation"](b, 5) == 1)
    check("five node a unit", a % 5 != 0)
    check("five split tangent", a % 5 == 4)
    check("five discriminant valuation below twelve", dep["valuation"](disc, 5) == 2)
    r, delta, p = 5328, 1295854, 1889
    check("exact nodal cubic coefficient identity",
          [0, b, a, 1] == [0, r*r-delta, -2*r, 1])
    check("1889 node discriminant factor", delta == p*686)
    check("1889 node discriminant valuation", dep["valuation"](delta, p) == 1)
    check("1889 node x unit", r % p != 0)
    check("1889 nonsplit tangent", pow(r, (p-1)//2, p) == p-1)
    check("1889 discriminant valuation below twelve", dep["valuation"](disc, p) == 1)
    counts = {}
    for p in [5, 1889]:
        # Complete reduction scan: one point at infinity, one y for zero,
        # two y for a nonzero quadratic residue, and zero otherwise.
        total = 1
        nodes = []
        for x in range(p):
            rhs = (x*x*x + a*x*x + b*x) % p
            total += 1 if rhs == 0 else 2 if pow(rhs, (p-1)//2, p) == 1 else 0
            if rhs == 0 and (3*x*x+2*a*x+b) % p == 0:
                nodes.append(x)
        expected_node = 0 if p == 5 else r % p
        check("unique reduced node " + str(p), nodes == [expected_node])
        expected_ns = 4 if p == 5 else 1890
        check("complete nonsingular reduction count " + str(p), total-1 == expected_ns)
        counts[str(p)] = {"all_projective_points": total,
                          "nonsingular_points": total-1, "node_x": expected_node}
    quotient5 = 2*counts["5"]["nonsingular_points"]
    quotient1889 = counts["1889"]["nonsingular_points"]
    check("five quotient-order arithmetic", quotient5 == 8)
    check("1889 quotient-order complete factorization", quotient1889 == 2*3**3*5*7)
    finite_support = {2, 3, 5, 7}
    check("finite core quotient primary support",
          prod([2**3]) == quotient5 and prod([2, 3**3, 5, 7]) == quotient1889)
    check("formal and finite primary dependency closure",
          finite_support | core == dependencies)
    points = {name: tuple(F(q) for q in point)
              for name, point in cert["points"].items()}
    unit_data = {}
    for name in ["Q", "U"]:
        point = points[name]
        translated = dep["add"](point, points["T"], a)
        check("literal T-translation identity " + name,
              translated == (F(b)/point[0], -b*point[1]/point[0]**2))
        check("T-translation remains on curve " + name, dep["on_curve"](translated, a, b))
        f = 128*point[0]-25865
        check("1889 F-unit core witness " + name, dep["valuation"](f, 1889) == 0)
        unit_data[name] = {"F_valuation_1889": 0,
                           "F_numerator_mod1889": f.numerator % 1889,
                           "F_denominator_mod1889": f.denominator % 1889,
                           "class": cert["classes"][name]}
    signs = list(product([-1, 1], repeat=4))
    # Character order is (-3, 2, 5*q0, d), as in EPE3.
    eligible2 = {s for s in signs if s[0] == -1 and s[1] == 1 and s[3] == 1}
    eligible5q = {s for s in signs if s[0] == -1 and s[2] == 1 and s[3] == 1}
    check("four-character sign space", len(signs) == 16)
    check("class-two numerator eligibility sign count", len(eligible2) == 2)
    check("class-five-q numerator eligibility sign count", len(eligible5q) == 2)
    check("common eligible sign count", len(eligible2 & eligible5q) == 1)
    check("eligible union sign count", len(eligible2 | eligible5q) == 3)
    check("class-specific safe complement density", F(16-len(eligible2), 16) == F(7, 8))
    check("common safe complement density", F(16-len(eligible2 | eligible5q), 16) == F(13, 16))
    check("Kummer sign pattern in both eligible sets", (-1, 1, 1, 1) in eligible2 & eligible5q)
    return {"passed": True, "exact_controls": len(checks), "checks": checks,
            "baseline_exact_controls_replayed": prior["exact_controls"],
            "test_core": sorted(core), "primary_dependency_set": sorted(dependencies),
            "complete_reduction_scan_primes": [5, 1889],
            "reduction_counts": counts, "core_witnesses": unit_data,
            "finite_checks_prove_induction": False,
            "unrestricted_all_finite_sets_claimed": False,
            "new_everywhere_local_seed_claimed": False,
            "independent_replay": False, "external_review": False,
            "formal_verification": False, "python": platform.python_version(),
            "arithmetic": "standard-library exact integers and fractions"}


if __name__ == "__main__":
    report = replay()
    names = [SOURCE, Path(__file__).name, BASELINE, BASELINE_INPUT, ARITHMETIC]
    report["input_sha256"] = {name: sha256((HERE/name).read_bytes()).hexdigest()
                              for name in names}
    (HERE / "elliptic_exception_dependency_closed_core_checks.json").write_text(
        json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ["passed", "exact_controls",
          "baseline_exact_controls_replayed", "test_core", "primary_dependency_set",
          "core_witnesses", "input_sha256"]}, indent=2))
