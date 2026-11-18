#!/usr/bin/env python3
"""Exact finite controls for the asymmetric sign-transfer analysis.

Standard library only; run from any working directory. The finite checks
cover convolution/sign/Rosser identities and least-squares projection.
They do not verify asymptotic means, dispersion or binary coverage.
No file is written unless --output is explicitly supplied.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from math import prod
from pathlib import Path
import sys

if not __debug__:
    raise RuntimeError("Run these finite controls without Python -O.")

ROOT = Path(__file__).resolve().parents[3]
SOURCE_PATHS = (
    "research/arithmetic/ASYMMETRIC_SIGN_TRANSFER.md",
    "research/arithmetic/ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md",
    "research/arithmetic/ASYMMETRIC_SEMILINEAR_SIEVE.md",
    "research/arithmetic/ASYMMETRIC_AUXILIARY_AVERAGE.md",
    "research/arithmetic/sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf",
    "research/arithmetic/sources/iwaniec-half-dimensional-sieve-1976.pdf",
    "research/arithmetic/sources/loughran-matthiesen-1904.12845v2.pdf",
)
EXPECTED_PRIMARY = {
    SOURCE_PATHS[4]: "6cb425e37285b1a09aeeade57210948f885d40f1d4d7b337671cfbd0a4a4711f",
    SOURCE_PATHS[5]: "ca6e84aafb32ee62dcd5ce76679c65652bf677ccd4e36d200576024fac1eca67",
    SOURCE_PATHS[6]: "ffbb0c2e6a6a39e090ebd6102f2a3b8936208df33c11bdf7c892744756c10203",
}
COUNTS = Counter()


def require(ok, label):
    if not ok:
        raise RuntimeError(label)
    COUNTS[label] += 1


@lru_cache(maxsize=None)
def factors(n):
    if n < 1:
        raise ValueError("Only positive integers are factored.")
    out = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            e += 1
            n //= p
        if e:
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return tuple(out)


@lru_cache(maxsize=None)
def divisors(n):
    out = [1]
    for p, e in factors(n):
        out = [d * p ** j for d in out for j in range(e + 1)]
    return tuple(out)


def tau(n):
    return prod(e + 1 for _, e in factors(n))


def mu(n):
    return 0 if any(e > 1 for _, e in factors(n)) else (-1) ** len(factors(n))


def p7(n):
    return all(p % 12 == 7 for p, _ in factors(n))


def chi4(n):
    return 0 if n % 2 == 0 else 1 if n % 4 == 1 else -1


def chi3(n):
    return 0 if n % 3 == 0 else 1 if n % 3 == 1 else -1


def u(n):
    return tau(n) ** 2 if all(p % 12 == 1 for p, _ in factors(n)) else 0


def v(n):
    return tau(n) if all(p % 3 == 1 for p, _ in factors(n)) else 0


def m7(n):
    return int(mu(n) != 0 and p7(n))


def b(n):
    fs = factors(n)
    return 2 ** sum(e == 1 for _, e in fs) if p7(n) and all(e <= 2 for _, e in fs) else 0


def b_sf(n):
    return 2 ** len(factors(n)) * m7(n)


def lower_rosser(d, D):
    if mu(d) == 0:
        return False
    ps = sorted((p for p, _ in factors(d)), reverse=True)
    return all(prod(ps[:j - 1]) * ps[j - 1] ** 2 < D
               for j in range(2, len(ps) + 1, 2))


def check_convolutions():
    for n in range(1, 3001):
        require(sum(m7(d) * m7(n // d) for d in divisors(n)) == b(n),
                "squarefree_convolution")
        target = tau(n) if all(p % 12 == 1 for p, _ in factors(n)) else 0
        require(sum(b(d) * chi4(n // d) * v(n // d) for d in divisors(n)) == target,
                "signed_convolution")


def check_binary_signs():
    for N in range(2, 500, 12):
        for a in range(1, N):
            if not u(a):
                continue
            for c in divisors(N - a):
                if not p7(c):
                    continue
                t = (N - a) // c
                if v(t):
                    require(chi4(c) == chi4(t), "binary_sign_transfer")
                if t % 6 == 1:
                    for q in divisors(t):
                        r = t // q
                        require(chi3(q) * chi4(q) * chi4(r) == chi4(c) * chi3(q),
                                "reverse_hyperbola_phase")
        for c in (1, 7, 19, 31, 49, 133):
            actual = sum(u(N - c * t) * v(t) for t in range(1, (N - 1) // c + 1))
            twisted = sum(u(N - c * t) * chi4(t) * v(t)
                          for t in range(1, (N - 1) // c + 1))
            require(twisted == chi4(c) * actual, "coefficient_resonance")


def check_rosser_and_compression():
    for D in (10, 100, 1000, 10000):
        z = D - 1
        for d in range(1, min(3001, D)):
            if mu(d) and p7(d) and d ** 3 < D ** 2:
                require(lower_rosser(d, D), "automatic_rosser_support")
        for c in range(1, 3001):
            if not b(c):
                continue
            fs = factors(c)
            k = prod(p for p, e in fs if e == 2)
            ell = prod(p for p, e in fs if e == 1)
            d = k * ell
            require(c == k * k * ell and d == prod(p for p, _ in fs),
                    "cubefree_unique_decomposition")
            for K in (1, 7, 100):
                if k > K:
                    continue
                selected = lower_rosser(d, D) and all(p < z for p, _ in factors(d))
                lam = mu(d) * int(selected)
                psi = b(c) * int(selected)
                require(mu(k) * 2 ** len(factors(ell)) * lam == chi4(c) * psi,
                        "compressed_sign_coefficient")


def check_least_squares():
    for fun in (b, b_sf):
        for X in (12, 37, 80, 151):
            vals = {n: F(fun(n)) for n in range(X + 1, 2 * X + 1)}
            for m in (1, 2, 3, 4, 6, 12, 17):
                classes = defaultdict(list)
                for n in vals:
                    classes[n % m].append(n)
                means = {r: sum(vals[n] for n in ns) / len(ns) for r, ns in classes.items()}
                total = sum(x * x for x in vals.values())
                projection = sum(len(classes[r]) * a * a for r, a in means.items())
                best = sum((vals[n] - means[n % m]) ** 2 for n in vals)
                require(best == total - projection, "least_squares_projection_identity")
                # Complex amplitudes represented by exact rational real/imaginary parts.
                real = {r: F(3 * r - m + X, 2) for r in classes}
                imag = {r: F(r - X, 3) for r in classes}
                error = sum((vals[n] - real[n % m]) ** 2 + imag[n % m] ** 2 for n in vals)
                excess = sum(len(ns) * ((real[r] - means[r]) ** 2 + imag[r] ** 2)
                             for r, ns in classes.items())
                require(error == best + excess and error >= best,
                        "least_squares_complex_pythagoras")
                S = {r: sum(vals[n] for n in ns) for r, ns in classes.items()}
                require(projection <= max(F(1, len(ns)) for ns in classes.values()) * sum(S.values()) ** 2,
                        "positive_projection_bound")


def check_negative_controls():
    require(133 % 12 == 217 % 12 and lower_rosser(133, 1000)
            and not lower_rosser(217, 1000), "control_selector_not_periodic")
    require(b_sf(133) == b_sf(217) == 4, "control_equal_unselected_coefficients")
    # N=62, a=13, c=7, v=1, q=1, r=7: the missing complementary twist is detected.
    require(u(13) > 0 and 13 + 7 * 1 * 1 * 7 == 62 and 62 % 12 == 2,
            "control_valid_resonance_pair")
    require(chi3(1) * chi4(1) != chi4(7) * chi3(1),
            "control_omitted_complementary_twist_fails")
    require(sum(b(d) * v(7 // d) for d in divisors(7)) != 0
            and sum(b(d) * chi4(7 // d) * v(7 // d) for d in divisors(7)) == 0,
            "control_omitted_signed_convolution_twist_fails")
    vals = [b(n) for n in range(13, 25)]
    require(sum(x * x for x in vals) - F(sum(vals) ** 2, len(vals)) > 0,
            "control_constant_replacement_not_exact")


def provenance():
    hashes = {s: hashlib.sha256((ROOT / s).read_bytes()).hexdigest() for s in SOURCE_PATHS}
    for s, expected in EXPECTED_PRIMARY.items():
        require(hashes[s] == expected, "immutable_primary_pin")
    return {
        "checker": Path(__file__).resolve().relative_to(ROOT).as_posix(),
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": hashes,
        "python_version": sys.version.split()[0],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional JSON receipt; otherwise print only.")
    args = parser.parse_args()
    if args.output is not None:
        protected = {(ROOT / rel).resolve() for rel in SOURCE_PATHS}
        protected.add(Path(__file__).resolve())
        if args.output.resolve() in protected:
            raise ValueError("The output path would overwrite a source or the checker.")
    pins = provenance()
    check_convolutions()
    check_binary_signs()
    check_rosser_and_compression()
    check_least_squares()
    check_negative_controls()
    receipt = {
        "result": "PASS",
        "scope": "Finite exact identities and least-squares controls only; no asymptotic mean, dispersion, critical remainder, or binary coverage is verified.",
        "provenance": pins,
        "checks": dict(sorted(COUNTS.items())),
        "negative_controls": [
            "same residue133/217 has different Rosser selection",
            "omitting complementary chi4(r) changes the reverse phase",
            "omitting the signed convolution twist fails at7",
            "a constant periodic replacement has nonzero finite residual",
        ],
    }
    serialized = json.dumps(receipt, indent=2) + "\n"
    if args.output is not None:
        args.output.write_text(serialized, encoding="utf-8")
    print(serialized, end="")


if __name__ == "__main__":
    main()
