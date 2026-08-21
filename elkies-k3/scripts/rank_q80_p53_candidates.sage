#!/usr/bin/env sage
"""Reproduce the retracted standard-P1 p=53 ranking experiment.

The source residues are now known to have the wrong residual-I2 component for
P2.  This is a historical heuristic diagnostic, not target data.  For every
possible residue of the remaining ambient
parameter, it finds short projective representatives

    p = (a + b*sqrt(-6))/c

compatible with the fully marked p=7,11,29 data and one of the two embeddings
at 29 and 53.  The exact P1 identity then forces
q=18-2p and e=(p-42)^2/36, so ranking p ranks the whole ambient surface.
"""

from sage.all import *
from itertools import product as itertools_product


primes = (7, 11, 29, 53)
modulus = prod(primes)


def short_values(residues, roots, coefficient_bound=7):
    root_crt = ZZ(CRT_list(list(roots), list(primes)))
    residue_crt = ZZ(CRT_list(list(residues), list(primes)))
    lattice = matrix(ZZ, [
        [modulus, 0, 0],
        [-root_crt, 1, 0],
        [residue_crt, 0, 1],
    ]).LLL()
    values = {}
    for coefficients in itertools_product(
        range(-coefficient_bound, coefficient_bound + 1), repeat=3
    ):
        if coefficients == (0, 0, 0):
            continue
        a, b, c = map(ZZ, vector(ZZ, coefficients) * lattice)
        if c == 0 or any(c % prime == 0 for prime in primes):
            continue
        common = gcd((a, b, c))
        a, b, c = a // common, b // common, c // common
        if c < 0:
            a, b, c = -a, -b, -c
        score = a*a + 6*b*b + c*c
        values[(a, b, c)] = min(score, values.get((a, b, c), infinity))
    return sorted(values, key=lambda triple: values[triple])


ranked = []
for root29, root53 in itertools_product((9, 20), (10, 43)):
    roots = (1, 4, root29, root53)
    for residue53 in range(53):
        values = short_values((4, 2, 5, residue53), roots)
        if not values:
            continue
        a, b, c = values[0]
        score = a*a + 6*b*b + c*c
        ranked.append((score, residue53, roots, a, b, c))

for rank, item in enumerate(sorted(ranked)[:32], 1):
    score, residue53, roots, a, b, c = item
    print(
        f"Q80P53RANK|rank={rank}|score={score}|p53={residue53}|"
        f"roots={roots}|p=({a}+({b})*sqrt(-6))/{c}",
        flush=True,
    )

print(
    "Q80P53RANK|SUMMARY|status=HEURISTIC_PRIORITY_ONLY|"
    f"candidates={len(ranked)}",
    flush=True,
)
