"""Replay the first wild binary row exclusions in F_9 with a local order.

Run with Python 3.12 and Singular 4.3.2. The arithmetic witnesses below
compare the two native coefficient constructions at every retained order.
"""

from dataclasses import dataclass
from functools import lru_cache
from math import comb, factorial
import subprocess


@dataclass(frozen=True)
class Profile:
    kind: str
    x: int = 0
    y: int = 0
    twice: int = 0


def digits_can_double(n):
    return n % 3 <= 1 and n // 3 <= 1


def profiles():
    for x in range(9):
        for y in range(9):
            yield Profile("ordinary", x, y)
            if digits_can_double(x) and digits_can_double(y):
                yield Profile("ordinary", x, y, -1)
                yield Profile("ordinary", x, y, 1)
    for kind in ("x-low", "x-high", "y-low", "y-high", "separate"):
        yield Profile(kind)


def residue_coefficient(profile, x, y, m):
    if profile.kind == "ordinary":
        if m > 2 or m == 2 and profile.twice == 0:
            return 0
        factor = profile.twice if m == 2 else 1
        return factor * comb(x + m * profile.x, x) * comb(y + m * profile.y, y) % 3

    if profile.kind == "separate":
        return 1
    # One slot apart. Only the low count's higher digits use the high slot.
    axis_count = x if profile.kind.startswith("x") else y
    low, high = (m, axis_count) if profile.kind.endswith("low") else (axis_count, m)
    return comb(high + low // 3, high) % 3


def arithmetic_example(profile):
    if profile.kind == "ordinary":
        if profile.twice == 0:
            return 3, 3, 3 * profile.x + 2, 3 * profile.y + 1
        i, j = 3 * profile.x + 1, 3 * profile.y + 1
        weight = 2 + sum(n % 3 + n // 3 for n in (profile.x, profile.y))
        if (-1) ** (weight - 1) != profile.twice:
            i += 27
        return 3, 3, i, j
    return {
        "x-low": (3, 3, 1, 0),
        "x-high": (1, 3, 3, 0),
        "y-low": (3, 3, 0, 1),
        "y-high": (3, 1, 0, 3),
        "separate": (9, 3, 1, 0),
    }[profile.kind]


factorial = lru_cache(None)(factorial)


def check_arithmetic(profile):
    qx, qy, i, j = arithmetic_example(profile)
    checks = 0
    for n in range(9):
        for m in range(n + 1):
            for x in range(n - m + 1):
                y = n - m - x
                numerator = factorial(qx * x + i * m) * factorial(qy * y + j * m)
                denominator = (
                    factorial(x) * factorial(y) * factorial(m)
                    * factorial(qx)**x * factorial(qy)**y
                    * (factorial(i) * factorial(j))**m
                )
                coefficient, remainder = divmod(numerator, denominator)
                if remainder or coefficient % 3 != residue_coefficient(profile, x, y, m):
                    raise ArithmeticError((profile, n, x, y, m, coefficient, remainder))
                checks += 1
    return checks


def polynomial_case(profile):
    if profile.kind != "ordinary":
        return False
    insertions = (1, 2) if profile.twice else (1,)
    return all(
        m + (2 - m * profile.x) % 9 + (2 - m * profile.y) % 9 <= 4
        for m in insertions
    )


CONSTANT = ("1", "k", "-k^2", "U", "k*U", "-k^2*U", "-U^2", "-k*U^2", "k^2*U^2")
SURVIVOR = Profile("ordinary", 0, 0, -1)


def native_polynomial(profile, n):
    terms = []
    for m in range(n + 1):
        for x in range(n - m + 1):
            y = n - m - x
            c = residue_coefficient(profile, x, y, m)
            if c:
                terms.append(f"{c}*u^{m}*b^{y}")
    return " + ".join(terms)


def singular_case(index, profile):
    define = [
        f"  poly R{n}={native_polynomial(profile, n)};"
        for n in range(9)
    ]
    equations = [
        " + ".join(f"R{n-j}*({CONSTANT[j]})" for j in range(n + 1))
        for n in range(5, 9)
    ]
    generators = ",\n    ".join(["k^3", *equations])
    expected = "u^2" if profile == SURVIVOR else "0"
    verify = f"""  ideal I={generators};
  ideal G=std(I);
  if (reduce(k^2,G)=={expected}) {{ print("PASS"); }}
  else {{ print("FAIL: {profile.kind} {profile.x} {profile.y} {profile.twice}"); }}
"""
    return f"proc check{index}()\n{{\n" + "\n".join(define) + "\n" + verify + f"}}\ncheck{index}();\n"


def main():
    cases = list(profiles())
    checks = sum(map(check_arithmetic, cases))
    algebra_cases = [profile for profile in cases if not polynomial_case(profile)]
    code = (
        "ring R=(3,i),(B,V,k,u),ds;\n"
        "minpoly=i^2+1;\n"
        "poly b=i+B;\n"
        "poly U=i-1+V;\n"
        + "".join(singular_case(index, profile) for index, profile in enumerate(algebra_cases))
        + "quit;\n"
    )
    run = subprocess.run(
        ["Singular", "-q"], input=code, text=True, capture_output=True, timeout=60
    )
    if run.returncode or run.stderr or run.stdout.splitlines() != ["PASS"] * len(algebra_cases):
        raise RuntimeError(run.stdout + run.stderr)
    print(f"{checks} exact native coefficients across {len(cases)} profiles")
    print(f"{len(cases) - len(algebra_cases)} polynomial cases")
    print(f"{len(algebra_cases) - 1} local k^2 containments and one survivor")


if __name__ == "__main__":
    main()
