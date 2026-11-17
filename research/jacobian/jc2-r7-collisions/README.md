# All affine singularity strata of the polynomial gap-21 family

14 November 2026. PF2R7COLL1 classifies every parameter in the polynomial
coefficient family L: twelve ordinary nodes generically, with exactly
34 exceptional parameters split into triple-point, tacnode and cusp
strata. Three small modular gcd certificates make the exceptional counts
uniform over their algebraic parameter sets.

The old sharp d=0 braid sample has three ordinary triple points and
three nodes. Its historical `generic` label must not be used to extend
that sample's braid result to the generic twelve-node family.

PF2R7COLL2 also forces every unknown additional omitted image to avoid
these triple values in a residual degree-six Keller candidate. The
source would omit those values entirely; this is a necessary restriction,
not a global impossibility proof.

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-r7-collisions/verify.py
```

[Proof](PROOF.md) · [strata](strata.json) · [certificates](certificates.json) ·
[replay](receipt.json) · [validation](validation.json) ·
[two complete coefficient families](../jc2-r7-image-family/README.md).

The [Laurent-family result](../jc2-r7-laurent-collisions/README.md) now
completes the other component as well, giving268 nonnodal coefficient
points across the complete normalized union. Actual second-image equations,
infinity comparison and global Keller realization remain open. Authored written geometry and exact
algebra; no formalization or independent external review.
