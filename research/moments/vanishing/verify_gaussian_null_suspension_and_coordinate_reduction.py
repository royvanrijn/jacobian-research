"""Bounded exact QQ controls for null suspension and its original map.

The unrestricted statements are the companion written proofs. No finite
moment packet is used to infer all-marker safety or polynomiality.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from math import factorial
from pathlib import Path
import platform

import sympy as sp

if not __debug__:
    raise RuntimeError("Run without -O so exact assertions execute")

x, v, y, Z, W, t = sp.symbols("x v y Z W t")
variables = (x, v, y, Z, W)


def zero(value):
    assert sp.expand(value) == 0, value


def wick(polynomial):
    answer = sp.Integer(0)
    for (a, b, c, d, e), coefficient in sp.Poly(sp.expand(polynomial), *variables).terms():
        if a == b and d == e and c % 2 == 0:
            answer += coefficient*factorial(a)*factorial(d)*sp.Rational(
                factorial(c), 2**(c//2)*factorial(c//2))
    return sp.expand(answer)


def period(phase, marker, order):
    return sp.expand(sum(wick(marker*phase**n)*t**n/factorial(n)
                         for n in range(order+1)))


def truncate(value, order):
    return sp.series(value, t, 0, order+1).removeO().expand()


def timed_period(phase, marker, order):
    result = 0
    for (power,), coefficient in sp.Poly(sp.expand(marker), t).terms():
        if power <= order:
            result += t**power*period(phase, coefficient, order-power)
    return sp.expand(result)


def degree_control(degree):
    # Formal old monomial/time pairs: no truncated period dependence enters.
    actual, expected = set(), set()
    for a, b, c in itertools.product(range(degree+1), repeat=3):
        old_degree = a+b+c
        if old_degree > degree:
            continue
        for az in range(degree-old_degree+1):
            for bw in range(degree-old_degree-az+1):
                if bw >= az:
                    assert sp.Rational(factorial(bw), factorial(bw-az)) != 0
                    actual.add(((a, b, c), bw-az))
        for power in range(degree-old_degree+1):
            expected.add(((a, b, c), power))
    assert actual == expected
    return {"marker_degree": degree, "formal_old_marker_time_pairs": len(actual)}


def checks():
    old = x**2*v**2+x+y**2/3+x*y
    suspended = old+Z
    zero(sp.diff(suspended, Z)-1)
    zero((Z+old).subs(Z, Z-old)-Z)
    zero((Z-old).subs(Z, Z+old)-Z)
    markers = [1, x, v, y, y**2, x*v]
    for marker in markers:
        for a, b in itertools.product(range(4), repeat=2):
            actual = period(suspended, marker*Z**a*W**b, 4)
            expected = 0 if b < a else truncate(
                sp.Rational(factorial(b), factorial(b-a))*t**(b-a)
                *period(old, marker, 4), 4)
            zero(actual-expected)
    for marker in markers+[Z, W, Z*W, Z**2*W, Z*W**3]:
        lifted = W*marker-sp.diff(marker, Z)
        zero(period(suspended, lifted, 4)-truncate(t*period(suspended, marker, 4), 4))
    # At P=0 the exact auxiliary formula gives all 1,t,...,t^D.
    for degree in range(9):
        for power in range(degree+1):
            zero(period(Z, W**power, degree)-t**power)
    filtered = [degree_control(degree) for degree in range(9)]
    # Nonlinear null-variable elimination: direct Wick versus finite jets.
    simple_old = x*v+y
    polynomials = [simple_old+Z*y+Z**2*x, (simple_old+Z)**2+(simple_old+Z),
                   (simple_old+Z)**3+2, sp.Integer(3)]
    nonlinear_controls = 0
    for phase in polynomials:
        phase0 = phase.subs(Z, 0)
        hs = [sp.Integer(1)]
        for _ in range(4):
            hs.append(sp.expand(sp.diff(hs[-1], Z)+t*sp.diff(phase, Z)*hs[-1]))
        for a, b in itertools.product(range(3), repeat=2):
            for marker in (sp.Integer(1), x, v, y**2):
                expected = 0 if b < a else timed_period(
                    phase0, marker*sp.Rational(factorial(b), factorial(b-a))
                    *hs[b-a].subs(Z, 0), 3)
                zero(period(phase, marker*Z**a*W**b, 3)-expected)
                nonlinear_controls += 1
    # The scalar-only boundary is itself a submersion: f=1-t, M=1, N=6.
    U, A, z, w, eta, V, T = sp.symbols("U A z w eta V T")
    cn = sum(sp.binomial(6, j+1)*eta**j for j in range(6))
    polynomial_seed_phase = A+w*z+V*(1+eta)-(1+eta)*cn*T**2/2
    zero(sp.diff(polynomial_seed_phase, A)-1)
    seed = 1-t
    zero(sp.expand((1-t)*sum(t**j for j in range(6))).series(t, 0, 6).removeO()-1)
    assert seed != 1
    # Its full all-marker formulas are imported, not inferred by this test.
    return {"status": "passed", "coefficient_field": "QQ",
            "suspension_marker_count": len(markers)*16,
            "suspension_time_orders": [0, 4],
            "Stein_lift_marker_count": len(markers)+5,
            "nonlinear_null_elimination_controls": nonlinear_controls,
            "nonlinear_time_orders": [0, 3],
            "filtered_formal_map": filtered,
            "scalar_seed_boundary": "partial_A P=1, f=1-t; the full imported image is nonfinite",
            "scope": "Bounded exact original contractions and formal filtered-map controls; unrestricted conclusions use the written proof.",
            "independent_replay": False, "formal_verification": False,
            "external_review": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    directory = Path(__file__).resolve().parent
    names = {"proof": "GAUSSIAN_NULL_SUSPENSION_AND_COORDINATE_REDUCTION.md",
             "angular_core_import": "GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md",
             "submersion_import": "GAUSSIAN_ANGULAR_CORE_LITERAL_REGULAR_SECTOR.md",
             "seed_import": "POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md",
             "safe_recomposition_import": "GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md",
             "vacuum_construction_import": "SAFE_SEED_DIVISIBILITY_AND_PHASE_MARKER_FILTRATION.md"}
    pins = {key+"_sha256": hashlib.sha256((directory/name).read_bytes()).hexdigest()
            for key, name in names.items()}
    pins["checker_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result.update(python=platform.python_version(), sympy=sp.__version__, **pins)
    output = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
