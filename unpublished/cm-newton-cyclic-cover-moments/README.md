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
zero whenever both moment equations vanish. An arithmetic
attachment lemma uses unmarked denominator witnesses to locate
finite vanishing cycles in the sphere-variation span. Two
conductor-18 pencils show how the resulting periods detect both
coherent inputs when escaping periods contribute. For the
degree-\((5,1)\) support, rows \(0,\ldots,7\) recover all eight
marker inputs wherever their determinant is nonzero. On every
line varying just \(B(0)\), this determinant has degree 12
with the same nonzero leading coefficient.

For the same conductor-18 top and lower support, an explicit nonempty
boundary stratum gives all-order coherent detection with six correction
columns. With three correction columns, shifts \(0,\ldots,13\) detect
both coherent inputs over every characteristic-zero coefficient field.
The window is optimal for the formal response, where an algebraic phase
misses a coherent direction through shift 12. On the specified
four-parameter slice, the ninth row is sharp. Realization of the formal
sharpness witnesses by actual Frobenius inputs remains open.

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
separable conductor-18 phase. The
[escaping-variation note](../../research/foundations/extended-geometry/CM_NEWTON_SPHERE_ESCAPING_VARIATIONS.md)
gives the attachment lemma, the two pencil calculations, and the
eight-row determinant. The
[tail-period note](../../research/foundations/extended-geometry/CM_NEWTON_CYCLIC_LIMIT_TAIL_PERIODS.md)
records the local escaping scalings. The
[fixed-boundary stratum](../../research/foundations/extended-geometry/CM_NEWTON_SPHERE_ATTACHMENT_STRATUM.md)
gives the geometric detection theorem and the sharp ninth-row slice.
The [fourteen-row window](../../research/foundations/extended-geometry/CM_NEWTON_FOURTEEN_ROW_LOCUS.md)
combines the affine integral-place calculation and the boundary charts
with two exact slope exclusions. It uses the
[seventeen-row calculation](../../research/foundations/extended-geometry/CM_NEWTON_SEVENTEEN_ROW_LOCUS.md)
for the underlying finite algebra and the formal sharpness witness.
The source archive includes these proof notes and their exact replay inputs.

## Build and replay

From the extracted ZIP root, run

```sh
SOURCE_DATE_EPOCH=1795651200 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=cm-newton-cyclic-cover-moments main.tex
python3 -I verification/verify_cm_newton_nonsplit_noncm_prime_power_moments.py
python3 -I verification/verify_cm_newton_unbounded_twisted_moment_quotients.py
sage -python verification/verify_cm_newton_moving_branch_ratio_certificates.py
python3 -I verification/verify_cm_newton_block_rank_moments.py
sage -python verification/verify_cm_newton_sphere_escaping_variations.py
sage -python verification/verify_cm_newton_cyclic_limit_tail_periods.py
sage -python verification/verify_cm_newton_sphere_attachment_stratum.py
sage -python verification/verify_cm_newton_seventeen_row_locus.py
sage -python verification/verify_cm_newton_fourteen_row_locus.py
```

The first replay evaluates individual coherent Frobenius maps and
full rational sphere moments on four examples, and checks the exact
rank at conductor \(31\) using the coherent orbit matrices. The second
computes twisted quotient dimensions, enumerates the residue classes,
and checks selected coefficient residues at conductors \(35\) and \(133\).
The moving-branch checker selects pivots modulo a prime and
verifies the three parameter identities by exact rational
multiplication. The block-rank checker replays the rank counts,
rank-two determinant and rigidity congruence. The two escaping
period checkers replay the exact critical-coordinate substitutions,
irreducibility tests, finite moment ranks, and local scaling algebra.
The stratum checker verifies the critical-value basis, the ninth-row
unit certificate and the first-eight blind example. The two window
checkers are for different scopes. The seventeen-row program replays
rank and finite-algebra checks. The fourteen-row program checks the
affine and weighted boundary reductions and the exact slope exclusions
using the shorter prefix. Its exact rational matrix calculations can
take substantially longer than the finite-field checks. See the
[verification record](verification.json) for completed replays.
External specialist review and priority assessment remain open.
