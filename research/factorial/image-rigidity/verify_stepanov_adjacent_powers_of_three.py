#!/usr/bin/env python3
"""Exact algebra guards; the unrestricted proof is in the companion note.

Only Python's standard library is used. Finite n/r checks are explicitly
reported as guards, not a replacement for the written universal arguments.
"""

import argparse
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, factorial, prod
from pathlib import Path
import platform


def trim(poly):
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result or [0]


def add(left, right):
    size = max(len(left), len(right))
    return trim([(left[i] if i < len(left) else 0)
                 + (right[i] if i < len(right) else 0) for i in range(size)])


def scale(poly, scalar):
    return trim([scalar * value for value in poly])


def sub(left, right):
    return add(left, scale(right, -1))


def mul(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return trim(result)


def power(poly, exponent):
    result = [1]
    while exponent:
        if exponent & 1:
            result = mul(result, poly)
        poly = mul(poly, poly)
        exponent //= 2
    return result


def derivative(poly):
    return trim([i * poly[i] for i in range(1, len(poly))])


def evaluate(poly, value):
    result = 0
    for coefficient in reversed(poly):
        result = result * value + coefficient
    return result


def mod(poly, prime):
    return trim([coefficient % prime for coefficient in poly])


def monic_remainder(poly, divisor):
    assert divisor[-1] == 1
    result = trim(poly)
    while result != [0] and len(result) >= len(divisor):
        shift = len(result) - len(divisor)
        coefficient = result[-1]
        result = sub(result, [0] * shift + scale(divisor, coefficient))
    return result


def remainder_mod(poly, divisor, prime):
    divisor = mod(divisor, prime)
    assert divisor != [0]
    result = mod(poly, prime)
    inverse = pow(divisor[-1], -1, prime)
    while result != [0] and len(result) >= len(divisor):
        shift = len(result) - len(divisor)
        coefficient = result[-1] * inverse % prime
        result = mod(sub(result, [0] * shift + scale(divisor, coefficient)), prime)
    return result


def gcd_mod(left, right, prime):
    left, right = mod(left, prime), mod(right, prime)
    while right != [0]:
        left, right = right, remainder_mod(left, right, prime)
    return mod(scale(left, pow(left[-1], -1, prime)), prime)


def model(n):
    result = [0] * (n + 1)
    for k in range(n + 1):
        result[n - k] = (-1) ** k * comb(n + 2 * k, 3 * k)
    return result


def compose_cube(poly):
    result = [0] * (3 * (len(poly) - 1) + 1)
    for i, coefficient in enumerate(poly):
        result[3 * i] = coefficient
    return result


def valuation3(number):
    assert number != 0
    result = 0
    while number % 3 == 0:
        result += 1
        number //= 3
    return result


@lru_cache(maxsize=None)
def unit_factorial(m):
    return prod((3 * j - 2) * (3 * j - 1) for j in range(1, m + 1))


# Sparse multivariate polynomials for universal low-degree identities.
def mconst(dimension, value):
    return {(0,) * dimension: value} if value else {}


def mvar(dimension, index):
    exponent = [0] * dimension
    exponent[index] = 1
    return {tuple(exponent): 1}


def madd(*polynomials):
    result = {}
    for poly in polynomials:
        for exponent, coefficient in poly.items():
            result[exponent] = result.get(exponent, 0) + coefficient
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def mscale(poly, scalar):
    return {exponent: scalar * coefficient for exponent, coefficient in poly.items()
            if scalar * coefficient}


def mmul(left, right):
    result = {}
    for a, av in left.items():
        for b, bv in right.items():
            exponent = tuple(x + y for x, y in zip(a, b))
            result[exponent] = result.get(exponent, 0) + av * bv
    return {exponent: coefficient for exponent, coefficient in result.items()
            if coefficient}


def mpow(poly, exponent):
    dimension = len(next(iter(poly)))
    result = mconst(dimension, 1)
    for _ in range(exponent):
        result = mmul(result, poly)
    return result


def mderivative(poly, index):
    result = {}
    for exponent, coefficient in poly.items():
        if exponent[index]:
            new = list(exponent)
            new[index] -= 1
            result[tuple(new)] = coefficient * exponent[index]
    return result


def structural_identities():
    one = mconst(2, 1)
    x, t = mvar(2, 0), mvar(2, 1)
    xt = mmul(x, t)
    b = madd(one, mscale(xt, -1))
    denominator = madd(mpow(b, 3), t)
    numerator = mpow(b, 2)
    identity = madd(mmul(madd(one, mscale(xt, 2)), numerator),
                    mscale(denominator, -1),
                    mmul(t, mderivative(denominator, 1)))
    assert identity == {}

    y1, y2, y3 = [mvar(3, i) for i in range(3)]
    total = madd(y1, y2, y3)
    pairs = madd(mmul(y1, y2), mmul(y1, y3), mmul(y2, y3))
    product = mmul(mmul(y1, y2), y3)
    identity = madd(mpow(total, 3), mscale(mmul(total, pairs), -3),
                    mscale(product, 3), mscale(mpow(y1, 3), -1),
                    mscale(mpow(y2, 3), -1), mscale(mpow(y3, 3), -1))
    assert identity == {}

    x, b, c, a = [mvar(4, i) for i in range(4)]
    one = mconst(4, 1)
    r = madd(b, mmul(x, c), mmul(mpow(x, 2), a))
    dr = madd(c, mscale(mmul(x, a), -1))
    rn = madd(mmul(madd(x, mscale(one, -1)), b),
              mmul(mpow(x, 2), c), mmul(mpow(x, 3), a))
    identity = madd(mmul(madd(x, mscale(one, -1)), r), mmul(x, dr),
                    mscale(rn, -1), mscale(mmul(mpow(x, 2), a), -1))
    assert all(coefficient % 3 == 0 for coefficient in identity.values())

    # Exact identity for each block of the factorial unit part.
    assert sub(mul([-2, 3], [-1, 3]), add([2], scale(mul([0, 1], [-1, 1]), 9))) == [0]
    return ["formal logarithmic-derivative numerator identity over Z[x,t]",
            "third-power Newton identity over Z[y1,y2,y3]",
            "common-root derivative identity over F3[x,B,C,a0]",
            "factorial unit-block identity over Z[j]"]


def normalization_guards():
    cases = []
    count = 0
    for r in range(1, 5):
        a = 3 ** r
        for n in range(1, 9):
            coefficients = []
            for k in range(n + 1):
                numerator = comb(n, k) * factorial(a * (n + 2 * k)) // factorial(a * n)
                b, remainder = divmod(numerator, 3 ** (a * k))
                assert remainder == 0
                unit_numerator = unit_factorial(k) * prod(
                    unit_factorial(3 ** j * (n + 2 * k)) for j in range(r))
                unit_denominator = prod(unit_factorial(3 ** j * n) for j in range(r))
                assert b * unit_denominator == comb(n + 2 * k, 3 * k) * unit_numerator
                assert valuation3(b) == valuation3(comb(n + 2 * k, 3 * k))
                assert b % 9 == ((-1) ** k * comb(n + 2 * k, 3 * k)) % 9
                coefficients.append(b)
                count += 1
            assert coefficients[0] == 1
            assert coefficients[-1] == unit_factorial(a * n)
            assert coefficients[-1] % 3 != 0
            # Hash binary integers, avoiding Python's optional decimal digit cap.
            payload = b"|".join(value.to_bytes((value.bit_length() + 7) // 8 or 1, "big")
                                for value in coefficients)
            cases.append({"r": r, "a": a, "n": n,
                          "coefficients_in_k_order_sha256": sha256(payload).hexdigest(),
                          "constant_mod9": coefficients[-1] % 9})
    # r=0 is deliberately outside this common-model congruence.
    outside_coefficient = factorial(3) // factorial(1) // 3
    assert outside_coefficient == 2 and outside_coefficient % 9 != (-1) % 9
    return count, cases


def sequence_guards(maximum=64):
    models = [model(n) for n in range(maximum + 2)]
    assert models[:3] == [[1], [-1, 1], [1, -4, 1]]
    for n in range(3, len(models)):
        recurrence = add(add(mul([-1, 3], models[n - 1]),
                             mul([0, 0, -3], models[n - 2])),
                         mul([0, 0, 0, 1], models[n - 3]))
        assert models[n] == recurrence

    aa = [[1]]
    for n in range(1, maximum + 2):
        aa.append(add(scale(aa[n - 1], -1), mul([0, 1], aa[n - 3]) if n >= 3 else [0]))
    for n in range(1, maximum + 1):
        am2 = compose_cube(aa[n - 2]) if n >= 2 else [0]
        am1 = compose_cube(aa[n - 1])
        an = compose_cube(aa[n])
        residual = add(add(an, mul([0, 1], am1)), mul([0, 0, 1], am2))
        assert mod(residual, 3) == mod(models[n], 3)
        residual_derivative = sub(am1, mul([0, 1], am2))
        assert mod(residual_derivative, 3) == mod(derivative(models[n]), 3)
        identity = sub(add(mul([-1, 1], models[n]),
                           mul([0, 1], derivative(models[n]))), models[n + 1])
        assert mod(identity, 3) == mod(mul([0, 0, 1], am2), 3)

    ss = [[3], [-1], [1]]
    for m in range(3, 3 * (maximum + 1) + 1):
        ss.append(add(scale(ss[m - 1], -1), mul([0, 1], ss[m - 3])))
    vv = [[3], [0], [0, 2]]
    for m in range(3, maximum + 2):
        vv.append(add(mul([0, 1], vv[m - 2]), mul([0, 0, 1], vv[m - 3])))
    for n in range(1, maximum + 2):
        assert ss[3 * n] == add(models[n], mul([0, 2], models[n - 1]))
        newton = add(sub(power(ss[n], 3), scale(mul(ss[n], vv[n]), 3)), [0] * n + [3])
        assert ss[3 * n] == newton

    residual_records = []
    for n in range(1, maximum + 1):
        lower, upper = models[n:n + 2]
        common = gcd_mod(lower, upper, 3)
        assert gcd_mod(common, derivative(lower), 3) == [1]
        if common == [1]:
            continue
        assert gcd_mod(common, derivative(common), 3) == [1]
        lower_rem = monic_remainder(lower, common)
        upper_rem = monic_remainder(upper, common)
        assert all(value % 3 == 0 for value in lower_rem + upper_rem)
        obstruction = sub(mul([value // 3 for value in upper_rem], derivative(lower)),
                          mul([value // 3 for value in lower_rem], derivative(upper)))
        obstruction = remainder_mod(obstruction, common, 3)
        assert gcd_mod(common, obstruction, 3) == [1]
        residual_records.append({"lower_n": n, "common_residual_degree": len(common) - 1,
                                 "common_coefficients_ascending_mod3": common,
                                 "first_hensel_obstruction_mod_g": obstruction})

    assert gcd_mod(models[3], models[4], 3) == [2, 1]
    assert evaluate(models[3], 7) % 9 == 0
    assert evaluate(models[4], 7) % 9 == 3
    return residual_records


def run(repo_root):
    structural = structural_identities()
    coefficient_count, normalization = normalization_guards()
    residual_records = sequence_guards()
    relative_dir = Path("research/factorial/image-rigidity")
    source_names = ["STEPANOV_ADJACENT_POWERS_OF_THREE_SFC.md",
                    "verify_stepanov_adjacent_powers_of_three.py"]
    sources = [{"path": str(relative_dir / name),
                "sha256": sha256((repo_root / relative_dir / name).read_bytes()).hexdigest()}
               for name in source_names]
    return {
        "status": "PASS",
        "base_commit": "19eb9912fb60142f8b620dac24896931d4ff3e46",
        "python_version": platform.python_version(),
        "dependencies": "Python standard library only",
        "scope": "exact algebra identities and finite guards for the companion unrestricted proof",
        "universal_theorem_proved_by_finite_sampling": False,
        "universal_theorem_source": "STEPANOV_ADJACENT_POWERS_OF_THREE_SFC.md, Sections 2-5",
        "source_pins": sources,
        "structural_identities": structural,
        "actual_normalization_domain": {"r": [1, 2, 3, 4], "lower_n_min": 1, "lower_n_max": 8,
                                        "all_k": "0 <= k <= n", "coefficient_count": coefficient_count},
        "actual_normalization_cases": normalization,
        "sequence_guard_domain": {"recurrence_lower_n_max": 65,
                                  "trace_terminal_N_min": 1, "trace_terminal_N_max": 65,
                                  "residual_adjacent_lower_n_min": 1,
                                  "residual_adjacent_lower_n_max": 64},
        "nontrivial_residual_gcd_and_first_lift_controls": residual_records,
        "negative_controls": {
            "plain_mod3_coprimality_fails": {"lower_n": 3, "gcd_ascending_mod3": [2, 1]},
            "that_disk_is_excluded_mod9": {"lower_root_mod9": 7, "upper_value_mod9": 3},
            "r_zero_does_not_have_the_stated_model": {"actual_n1_k1_mod9": 2,
                                                       "model_n1_k1_mod9": 8},
            "ramified_integrality_is_insufficient": "v3(3^(1/3))=1/3; the written proof first proves unramifiedness"
        },
        "assurance": "internal exact replay; no formal verification or external review claimed"
    }


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Run without -O or -OO: exact assertions must execute.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = run(args.repo_root.resolve())
    serialized = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(serialized)
    print(serialized, end="")
