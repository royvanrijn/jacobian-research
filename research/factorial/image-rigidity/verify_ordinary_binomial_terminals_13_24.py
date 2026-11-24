#!/usr/bin/env python3
"""Exact annulus complements and Bézout replay at ordinary terminals 13–24."""
import hashlib
import json
import math
import platform
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import flint


FIRST, LAST = 13, 24
PRIMES = (1000003, 1000033)
GAP_CUTOFFS = (73, 81, 89, 97, 105, 114, 122, 131, 139, 148, 157, 166)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def complement(n, gap):
    radii = []
    for m in (n - 1, n):
        c = next(c for c in range(1, 2 * m) if (c + 1) ** m < 2 * c ** m)
        radii.append(c)
    b0 = math.factorial((gap + 1) * (n - 1)) // math.factorial(n - 1)
    b1 = math.factorial((gap + 1) * n) // math.factorial(n)
    margin = b1 ** (n - 1) - (radii[0] * radii[1]) ** (n * (n - 1)) * b0 ** n
    require(margin > 0, f"Strict annulus inequality at N={n}")
    cutoff_zero = (4 * n * n - 1).bit_length()
    require(2 ** cutoff_zero >= 4 * n * n, f"Constant-family cutoff at N={n}")
    return radii, cutoff_zero


def factorial_tables(p, last):
    factorial = [1] * (last + 1)
    inverse = [1] * (last + 1)
    for j in range(1, last + 1):
        factorial[j] = factorial[j - 1] * j % p
    inverse[last] = pow(factorial[last], -1, p)
    for j in range(last, 0, -1):
        inverse[j - 1] = inverse[j] * j % p
    return factorial, inverse


def moment(m, a, d, p, tables, binomial):
    factorial, inverse = tables
    start = m * a
    normalization = inverse[start]
    coefficients = [
        normalization * factorial[start + d * k] * binomial[k] % p
        for k in range(m + 1)
    ]
    result = flint.nmod_poly(coefficients, p)
    require(result.degree() == m and result[0] == 1, f"Exact degree for {(m, a, d)}")
    return result


def domain(n, gap, cutoff_zero):
    for d in range(1, gap):
        for a in range(1, d * (n - 2)):
            yield a, d
    for d in range(1, cutoff_zero):
        yield 0, d


def verify_region(job):
    n, gap = job
    radii, cutoff_zero = complement(n, gap)
    last_argument = max(n * ((gap - 1) * (n - 1) - 1), n * (cutoff_zero - 1))
    require(n >= FIRST and n <= LAST and 1 < cutoff_zero < gap, "Terminal domain")
    for p in PRIMES:
        require(p > last_argument and all(p % j for j in range(2, math.isqrt(p) + 1)),
                f"Verified degree-preserving prime {p} at N={n}")
    tables = {p: factorial_tables(p, last_argument) for p in PRIMES}
    binomial = [tuple(math.comb(m, k) for k in range(m + 1)) for m in (n - 1, n)]

    def pair(a, d, p):
        return (moment(n - 1, a, d, p, tables[p], binomial[0]),
                moment(n, a, d, p, tables[p], binomial[1]))

    for p in PRIMES:
        for a, d in ((0, 2), (1, 1), (4, 10)):
            for m, polynomial in zip((n - 1, n), pair(a, d, p)):
                expect = [
                    (math.comb(m, k) * math.prod(range(a*m + 1, a*m + d*k + 1))) % p
                    for k in range(m + 1)
                ]
                require(list(polynomial) == expect, "Independent rising-factorial coefficient check")

    p0, p1 = PRIMES
    stream = hashlib.sha256()
    got = 0
    fallbacks = []
    for a, d in domain(n, gap, cutoff_zero):
        p = p0
        lower, upper = pair(a, d, p)
        g, u, v = lower.xgcd(upper)
        if g != 1:
            fallbacks.append({"a": a, "d": d, "first_gcd_degree": g.degree()})
            p = p1
            lower, upper = pair(a, d, p)
            g, u, v = lower.xgcd(upper)
        require(g == 1 and u * lower + v * upper == 1, f"Full Bezout at {(n, a, d)}")
        stream.update(f"{a},{d},{p}\n".encode())
        got += 1
    want = (n - 2) * gap * (gap - 1) // 2 - (gap - 1) + cutoff_zero - 1
    require(got == want, f"Complete domain at N={n}")
    return {
        "terminal": n,
        "annulus_radii": radii,
        "gap_cutoff": gap,
        "constant_gap_cutoff": cutoff_zero,
        "maximum_factorial_argument": last_argument,
        "cases": got,
        "fallbacks": fallbacks,
        "case_stream_sha256": stream.hexdigest(),
    }


def main():
    jobs = list(zip(range(FIRST, LAST + 1), GAP_CUTOFFS, strict=True))
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(verify_region, jobs))
    result = {
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "python_flint": flint.__version__,
        "primes": PRIMES,
        "total_cases": sum(row["cases"] for row in rows),
        "terminals": rows,
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
