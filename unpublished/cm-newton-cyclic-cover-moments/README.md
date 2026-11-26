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

[Manuscript PDF](cm-newton-cyclic-cover-moments.pdf) ·
[LaTeX source and replay ZIP](cm-newton-cyclic-cover-moments-source.zip) ·
[Metadata](metadata.json) · [Verification record](verification.json) ·
[Checksums](SHA256SUMS.txt)

The proof sources are the [fixed-prime construction](../../research/foundations/extended-geometry/CM_NEWTON_NONSPLIT_NONCM_PRIME_POWER_MOMENTS.md)
and the [twisted quotient construction](../../research/foundations/extended-geometry/CM_NEWTON_UNBOUNDED_TWISTED_MOMENT_QUOTIENTS.md).
This manuscript accompanies the [Newton polygon paper](../cm-newton-polygon-constraints/),
which develops density theorems and global constraints on Frobenius norms.

## Build and replay

From the extracted ZIP root, run

```sh
SOURCE_DATE_EPOCH=1795651200 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=cm-newton-cyclic-cover-moments main.tex
python3 -I verification/verify_cm_newton_nonsplit_noncm_prime_power_moments.py
python3 -I verification/verify_cm_newton_unbounded_twisted_moment_quotients.py
```

The first replay evaluates individual coherent Frobenius maps and
full rational sphere moments on four examples, and checks the exact
rank at conductor \(31\) using the coherent orbit matrices. The second
computes twisted quotient dimensions, enumerates the residue classes,
and checks selected coefficient residues at conductors \(35\) and \(133\).
The manuscript
proves the parameterized statements by symbolic arguments. External
specialist review and priority assessment remain open.
