#!/usr/bin/env python3
"""Exact finite character checks for the prime coset 3 modulo 10.

This certifies the quartic Fourier algebra, character centering coefficients,
fixed residue factors, and rational spectral-gap arithmetic only. It does not
prove analytic estimates, the comparison theorem, positivity at large scales,
or an onset of representability. Standard library only; no floating arithmetic.

python3 research/arithmetic/verify_prime_coset_characters.py \
  --output research/arithmetic/prime_coset_character_checks.json
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


# A Gaussian rational is an exact pair (real, imaginary).
ZERO = (Fraction(0), Fraction(0))
ONE = (Fraction(1), Fraction(0))
I = (Fraction(0), Fraction(1))
UNITS = (1, 3, 7, 9)
GENERATOR_POWERS = (1, 3, 9, 7)


def require(condition, detail):
    if not condition:
        raise AssertionError(detail)


def add(z, w):
    return (z[0] + w[0], z[1] + w[1])


def mul(z, w):
    return (z[0] * w[0] - z[1] * w[1],
            z[0] * w[1] + z[1] * w[0])


def scale(z, q):
    return (z[0] * q, z[1] * q)


def conj(z):
    return (z[0], -z[1])


def power(z, n):
    require(n >= 0, "only nonnegative powers are used")
    out = ONE
    for _ in range(n):
        out = mul(out, z)
    return out


def total(values):
    out = ZERO
    for value in values:
        out = add(out, value)
    return out


def encoded(z):
    return {"real": str(z[0]), "imaginary": str(z[1])}


def check():
    chi = {a: power(I, j) for j, a in enumerate(GENERATOR_POWERS)}
    fourier = [scale(conj(power(I, k)), Fraction(1, 4))
               for k in range(4)]
    singular = [scale(power(I, t), Fraction(1, 4))
                for t in range(4)]

    first_euler_cases = []
    for t in range(4):
        for a in UNITS:
            direct = power(chi[a], t) if a == 3 else ZERO
            factored = total(mul(fourier[k], power(chi[a], t+k))
                             for k in range(4))
            require(direct == factored, f"twisted Euler coefficient t={t}, a={a}")
            first_euler_cases.append({"twist_power": t, "prime_residue": a,
                                      "coefficient": encoded(direct)})
        # The zeta constituent occurs at k=-t modulo four.
        require(fourier[(-t) % 4] == singular[t],
                f"singular exponent conjugation at twist {t}")

    # A model sum bar(psi(v))*M_psi(v) has exactly the indicated twisted
    # mean under the uniform unit measure. This is finite Fourier algebra.
    centering = []
    for t in range(4):
        row = []
        for k in range(4):
            entry = scale(total(mul(power(chi[a], t),
                                    conj(power(chi[a], k)))
                                for a in UNITS), Fraction(1, 4))
            require(entry == (ONE if t == k else ZERO),
                    f"centering coefficient t={t}, k={k}")
            row.append(encoded(entry))
        centering.append(row)

    # Negative controls deliberately corrupt the conjugation or omit both
    # complex branches; the resulting exact Fourier identities must fail.
    wrong = [scale(power(I, k), Fraction(1, 4)) for k in range(4)]
    wrong_value = total(mul(wrong[k], power(chi[3], k)) for k in range(4))
    require(wrong_value == ZERO, "wrong-conjugation coefficient should vanish")
    require(wrong_value != ONE, "negative control accepted wrong conjugation")
    truncated = total(mul(fourier[k], power(chi[3], k)) for k in (0, 2))
    require(truncated == (Fraction(1, 2), Fraction(0)),
            "expected precise missing-complex-branch defect")
    require(truncated != ONE, "negative control accepted dropped complex branches")
    require(singular[1] != ZERO and singular[3] != ZERO,
            "complex twisted singular means must survive")

    residue_cases = []
    for target in range(10):
        roots = [[u, v] for u in UNITS for v in UNITS
                 if (u + v - target) % 10 == 0]
        gamma = Fraction(10 * len(roots), len(UNITS)**2)
        expected = (Fraction(0) if target % 2 else
                    Fraction(5, 2) if target % 5 == 0 else Fraction(15, 8))
        require(gamma == expected, f"fixed residue factor at target {target}")
        residue_cases.append({"target_mod_10": target, "root_count": len(roots),
                              "Gamma": str(gamma), "ordered_roots": roots})

    rho = Fraction(1, 4)
    gaps = [rho - z[0] for z in singular[1:]]
    require(gaps == [Fraction(1, 4), Fraction(1, 2), Fraction(1, 4)],
            "branch real-part gaps")
    relative = [Fraction(3, 5) * gap for gap in gaps]
    require(relative == [Fraction(3, 20), Fraction(3, 10), Fraction(3, 20)],
            "rough-model relative exponents")
    delta = min(Fraction(1, 8192 * 4), Fraction(3, 10) * min(gaps))
    require(delta == Fraction(1, 32768), "final stated logarithmic saving")

    return {
        "status": "passed",
        "scope": "Exact finite character and residue algebra only; no analytic proof or onset bound.",
        "arithmetic": "Gaussian rational pairs using fractions.Fraction; no floating-point operations.",
        "modulus": 10,
        "prime_coset": [3],
        "generated_residue_group": list(UNITS),
        "character_generator": {"residue": 3, "value": encoded(I)},
        "fourier_L_exponents": [encoded(z) for z in fourier],
        "twisted_singular_exponents": [encoded(z) for z in singular],
        "twisted_first_Euler_coefficient_cases": first_euler_cases,
        "twisted_first_Euler_coefficients_checked": len(first_euler_cases),
        "centering_matrix": centering,
        "centering_coefficients_checked": 16,
        "fixed_residue_cases": residue_cases,
        "fixed_residue_cases_checked": len(residue_cases),
        "real_part_gaps": [str(gap) for gap in gaps],
        "relative_rough_model_exponents": [str(value) for value in relative],
        "logarithmic_saving": str(delta),
        "negative_controls_rejected": [
            "Fourier L exponents without conjugation (coefficient at prime residue 3 becomes zero)",
            "Omitting both complex branches (coefficient at prime residue 3 becomes 1/2 instead of 1)",
        ],
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = check()
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
        print(f"Passed 16 twisted Euler identities, 16 centering coefficients, "
              f"10 residue factors, exact gaps, and 2 negative controls; wrote {args.output}")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
