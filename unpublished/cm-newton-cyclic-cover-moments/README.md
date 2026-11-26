# Nonzero sphere moments from nonordinary cyclic covers

The draft proves that one coherent Frobenius line can certify a
nonzero sphere moment even when the primitive Jacobian factor of the
associated cyclic cover is nonordinary. The first construction gives
infinitely many congruences on each fixed cover at its certifying prime,
with unbounded conductor and least sphere degree across the family.
The second fixes the cover and phase. It gives an exact residue-class
criterion for the distinguished coefficient and proves that every good
nonsplit prime is nonordinary. At conductor \(35\) the certifying
nonordinary classes have density \(7/24\). At the fixed prime five, the
explicit digit formula gives \(5\)-rank \(42\) in dimension \(60\) for
conductor \(31\), and the rank fraction tends to one half as the conductor
grows.
Both constructions give primitive factors that are not potentially CM.

The conductor-18 reconstruction theorem covers every distinct
nonzero branch pair, quartic odd and radial even lower terms,
and arbitrary characteristic-zero coefficient fields. An
effective finite window determines the full reconstruction
kernel on fixed support. The conductor-162 example measures a
rank-two Hasse–Witt block and gives an isolated two-marker
zero whenever both moment equations vanish.

[Manuscript PDF](cm-newton-cyclic-cover-moments.pdf) ·
[LaTeX source and replay ZIP](cm-newton-cyclic-cover-moments-source.zip) ·
[Metadata](metadata.json) · [Verification record](verification.json) ·
[Checksums](SHA256SUMS.txt)

The proof sources are the [fixed-prime construction](../../research/foundations/extended-geometry/CM_NEWTON_NONSPLIT_NONCM_PRIME_POWER_MOMENTS.md)
and the [twisted quotient construction](../../research/foundations/extended-geometry/CM_NEWTON_UNBOUNDED_TWISTED_MOMENT_QUOTIENTS.md).
This manuscript accompanies the [Newton polygon paper](../cm-newton-polygon-constraints/),
which develops density theorems and global constraints on Frobenius norms.
The reconstruction arguments are supported by the
[moving-branch theorem](../../research/foundations/extended-geometry/CM_NEWTON_MOVING_BRANCH_RATIO_CERTIFICATES.md),
[effective cutoff](../../research/foundations/extended-geometry/CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#effective-bounds-for-the-reconstruction-window),
and [block-rank and rigidity calculations](../../research/foundations/extended-geometry/CM_NEWTON_BLOCK_RANK_MOMENTS.md).
The [sphere period note](../../research/foundations/extended-geometry/CM_NEWTON_SPHERE_PERIOD_CYCLES.md)
computes the limiting periods of two sphere variations in a
separable conductor-18 phase.

## Build and replay

From the extracted ZIP root, run

```sh
SOURCE_DATE_EPOCH=1795651200 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=cm-newton-cyclic-cover-moments main.tex
python3 -I verification/verify_cm_newton_nonsplit_noncm_prime_power_moments.py
python3 -I verification/verify_cm_newton_unbounded_twisted_moment_quotients.py
sage -python verification/verify_cm_newton_moving_branch_ratio_certificates.py
python3 -I verification/verify_cm_newton_block_rank_moments.py
```

The first replay evaluates individual coherent Frobenius maps and
full rational sphere moments on four examples, and checks the exact
rank at conductor \(31\) using the coherent orbit matrices. The second
computes twisted quotient dimensions, enumerates the residue classes,
and checks selected coefficient residues at conductors \(35\) and \(133\).
The moving-branch checker selects pivots modulo a prime and
verifies the three parameter identities by exact rational
multiplication. The last checker replays the rank counts,
rank-two determinant and rigidity congruence. External
specialist review and priority assessment remain open.
