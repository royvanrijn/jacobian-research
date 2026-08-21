#!/usr/bin/env python3
"""Verify determinant obstructions to rank-zero MW fibrations.

For this K3, |disc(NS)| = 948 = 2^2 * 3 * 79.  If an elliptic
fibration had Mordell--Weil rank zero, Shioda's determinant formula would
give det(R) = 948 * |MW_tors|^2 for its rank-17 ADE root lattice R.  But
every irreducible ADE factor of rank at most 17 has determinant supported
on primes at most 17, so det(R) cannot be divisible by 79.

The reconstructed rational specialization instead has Picard rank 20 and
|disc(NS)|=43.  Rank-zero would require a rank-18 ADE determinant divisible
by 43, which is impossible for the same reason.  The enumeration below
independently checks both lattices.
"""

from functools import lru_cache


factors: list[tuple[str, int, int]] = []
for rank in range(1, 19):
    factors.append((f"A{rank}", rank, rank + 1))
for rank in range(4, 19):
    factors.append((f"D{rank}", rank, 4))
factors.extend([("E6", 6, 3), ("E7", 7, 2), ("E8", 8, 1)])
factors.sort(key=lambda item: (item[1], item[0]))


@lru_cache(maxsize=None)
def determinant_products(rank: int, start: int) -> frozenset[int]:
    if rank == 0:
        return frozenset({1})
    products: set[int] = set()
    for i in range(start, len(factors)):
        _, factor_rank, determinant = factors[i]
        if factor_rank > rank:
            break
        for tail in determinant_products(rank - factor_rank, i):
            products.add(determinant * tail)
    return frozenset(products)


def main() -> None:
    cases = (("generic_rank19", 17, 948, 79), ("picard20", 18, 43, 43))
    for name, root_rank, discriminant, prime in cases:
        assert discriminant % prime == 0
        assert all(det % prime for _, rank, det in factors if rank <= root_rank)
        products = determinant_products(root_rank, 0)
        assert products
        assert all(det % prime for det in products)
        assert all(det % discriminant for det in products)
        print(
            "MW0_OBSTRUCTION|status=PASS"
            f"|case={name}|root_rank={root_rank}"
            f"|ADE_determinants={len(products)}"
            f"|disc_NS={discriminant}|obstructing_prime={prime}"
        )


if __name__ == "__main__":
    main()
