#!/usr/bin/env python3
"""Finite-field audit of the standalone D1/F2 coefficient-space map.

The normalized seeds are

    H(W) = W^2 (1-W) Q(W),  deg(Q) = N-3,  Q(1) = 1.

For each requested (N, p), this script exhausts the F_p-rational seed
space, retains the open where H/W^2 and H'' are squarefree of their expected
degrees, and groups seeds by

* the marked Hessian divisor [H''] (the selected root remains W=1), and
* the unmarked divisor modulo W -> aW for a in F_p^*.

The marked grouping detects failures of affine-mark faithfulness.  These
occur in small characteristic because Frobenius terms are killed by the
second derivative.  The calculation is deliberately independent of the
weighted-map reconstruction formulas used elsewhere in the repository.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import product


def trim(f: tuple[int, ...], p: int) -> tuple[int, ...]:
    values = list(f)
    while values and values[-1] % p == 0:
        values.pop()
    return tuple(x % p for x in values)


def degree(f: tuple[int, ...], p: int) -> int:
    return len(trim(f, p)) - 1


def derivative(f: tuple[int, ...], p: int) -> tuple[int, ...]:
    return trim(tuple((i * f[i]) % p for i in range(1, len(f))), p)


def divmod_poly(
    numerator: tuple[int, ...], denominator: tuple[int, ...], p: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    num = list(trim(numerator, p))
    den = trim(denominator, p)
    if not den:
        raise ZeroDivisionError("polynomial division by zero")
    quotient = [0] * max(0, len(num) - len(den) + 1)
    inv_lead = pow(den[-1], -1, p)
    while len(num) >= len(den):
        shift = len(num) - len(den)
        coefficient = num[-1] * inv_lead % p
        quotient[shift] = coefficient
        for i, value in enumerate(den):
            num[i + shift] = (num[i + shift] - coefficient * value) % p
        while num and num[-1] == 0:
            num.pop()
    return trim(tuple(quotient), p), trim(tuple(num), p)


def gcd_poly(f: tuple[int, ...], g: tuple[int, ...], p: int) -> tuple[int, ...]:
    left, right = trim(f, p), trim(g, p)
    while right:
        _, remainder = divmod_poly(left, right, p)
        left, right = right, remainder
    if not left:
        return ()
    scale = pow(left[-1], -1, p)
    return tuple(scale * x % p for x in left)


def squarefree(f: tuple[int, ...], p: int) -> bool:
    f = trim(f, p)
    return bool(f) and degree(gcd_poly(f, derivative(f, p), p), p) == 0


def hessian(h: tuple[int, ...], p: int) -> tuple[int, ...]:
    return derivative(derivative(h, p), p)


def projective_key(f: tuple[int, ...], p: int) -> tuple[int, ...]:
    f = trim(f, p)
    first = next(value for value in f if value)
    inverse = pow(first, -1, p)
    return tuple(inverse * value % p for value in f)


def scale_variable(f: tuple[int, ...], a: int, p: int) -> tuple[int, ...]:
    return tuple(value * pow(a, i, p) % p for i, value in enumerate(f))


def unmarked_key(f: tuple[int, ...], p: int) -> tuple[int, ...]:
    return min(projective_key(scale_variable(f, a, p), p) for a in range(1, p))


def scaling_stabilizer(f: tuple[int, ...], p: int) -> tuple[int, ...]:
    key = projective_key(f, p)
    return tuple(
        a
        for a in range(1, p)
        if projective_key(scale_variable(f, a, p), p) == key
    )


def normalized_seeds(n: int, p: int):
    """Yield coefficient tuples for exact-degree normalized seeds."""
    q_degree = n - 3
    # Q has q_degree+1 coefficients and Q(1)=1.  Eliminate its top one.
    for initial in product(range(p), repeat=q_degree):
        top = (1 - sum(initial)) % p
        q = initial + (top,)
        if q[0] == 0 or top == 0:
            continue
        h = [0] * (n + 1)
        for i, value in enumerate(q):
            h[i + 2] = (h[i + 2] + value) % p
            h[i + 3] = (h[i + 3] - value) % p
        yield tuple(h)


def clean_seed(h: tuple[int, ...], n: int, p: int) -> bool:
    primitive = trim(h[2:], p)
    k = hessian(h, p)
    return (
        degree(h, p) == n
        and degree(primitive, p) == n - 2
        and squarefree(primitive, p)
        and degree(k, p) == n - 2
        and k[0] != 0
        and squarefree(k, p)
    )


def format_poly(f: tuple[int, ...], variable: str = "W") -> str:
    terms = []
    for exponent, coefficient in enumerate(f):
        if not coefficient:
            continue
        monomial = "" if exponent == 0 else variable if exponent == 1 else f"{variable}^{exponent}"
        terms.append(str(coefficient) if not monomial else f"{coefficient}*{monomial}")
    return " + ".join(terms) if terms else "0"


def audit(n: int, p: int) -> dict[str, object]:
    clean = []
    marked = defaultdict(list)
    unmarked = defaultdict(list)
    symmetric = []
    for h in normalized_seeds(n, p):
        if not clean_seed(h, n, p):
            continue
        clean.append(h)
        k = hessian(h, p)
        marked[projective_key(k, p)].append(h)
        unmarked[unmarked_key(k, p)].append(h)
        stabilizer = scaling_stabilizer(k, p)
        if len(stabilizer) > 1:
            symmetric.append((h, stabilizer))

    marked_collisions = [fiber for fiber in marked.values() if len(fiber) > 1]
    return {
        "degree": n,
        "prime": p,
        "clean": len(clean),
        "marked_fibers": len(marked),
        "marked_collisions": len(marked_collisions),
        "largest_marked_fiber": max(map(len, marked.values()), default=0),
        "largest_unmarked_rational_fiber": max(map(len, unmarked.values()), default=0),
        "symmetric": len(symmetric),
        "marked_witness": marked_collisions[0][:2] if marked_collisions else None,
        "symmetry_witness": symmetric[0] if symmetric else None,
    }


DEFAULT_CASES = ((4, 3), (4, 5), (5, 3), (5, 7), (6, 3), (6, 5), (6, 7), (8, 3), (8, 5), (8, 7), (8, 11))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "cases",
        nargs="*",
        metavar="N:P",
        help="degree/prime pairs (default: the built-in degree 4,5,6,8 audit)",
    )
    args = parser.parse_args()
    cases = tuple(tuple(map(int, item.split(":"))) for item in args.cases) or DEFAULT_CASES

    for n, p in cases:
        result = audit(n, p)
        print(
            f"N={n} p={p}: clean={result['clean']} "
            f"marked_keys={result['marked_fibers']} "
            f"marked_collision_keys={result['marked_collisions']} "
            f"max_marked_fiber={result['largest_marked_fiber']} "
            f"max_unmarked_Fp_fiber={result['largest_unmarked_rational_fiber']} "
            f"Fp_symmetric={result['symmetric']}"
        )
        witness = result["marked_witness"]
        if witness:
            left, right = witness
            print(f"  marked collision H1={format_poly(left)}")
            print(f"                   H2={format_poly(right)}")
            print(f"                   [H1'']=[H2'']={projective_key(hessian(left, p), p)}")
        symmetry = result["symmetry_witness"]
        if symmetry:
            h, stabilizer = symmetry
            print(f"  symmetry H={format_poly(h)} stabilizer={stabilizer}")


if __name__ == "__main__":
    main()
