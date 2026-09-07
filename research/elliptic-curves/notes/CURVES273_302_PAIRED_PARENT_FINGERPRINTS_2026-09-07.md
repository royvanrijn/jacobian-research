# Curves 273 and 302: paired determinant-1092 and strict-local audit

This is a paired, fail-closed comparison of the public record curves 273 and
302. It does not identify the discoverers' construction.

The cheapest exact test is negative.  For the explicit arithmetic-MW17
determinant-1092 parent recovered through curve 302, the numerator of

\[
j_{1092}(t)-j_{273}
\]

has degree 24, is irreducible over `QQ`, and has no rational root.  Hence
curve 273 is not a rational fibre of that particular parent, even after a
quadratic twist.  The test says nothing about a different parent, a higher
degree base change, or the original construction.

## Frozen rank-17 core on 273

The same frozen numerical-height selector and exact integral-point shell
recovery used for the curve-302 determinant-1092 lead was applied to the
already stored primitive rank-17 candidate subspace of the 30 public points
of 273.  It selects 1,036 antipodal rays and retains 450 integral rays.  The
450 constant-norm equations have a one-dimensional rational solution space;
normalizing the common norm to four gives an even positive rank-17 form
with

```text
determinant       1020 = 2^2 * 3 * 5 * 17
Smith factors     1^16, 1020
minimum           4
signed norm-4 shell size  2526
```

All conclusions after the frozen numerical ray selection are exact.  In
particular this is not the determinant-1092 form.

It also does not give a substitute Shimura/K3 lead: the eight even ternary
genera of signature `(2,1)` and determinant `-1020` contain none with the
discriminant form required to complement `U + (-G)`.  Thus this particular
recovered form cannot be the primitive Mordell--Weil lattice of a rootless
rank-17 elliptic K3.  This does **not** exclude a different rank-17
subspace of curve 273 or a non-K3 source.

## Strict-local cubic comparison

For each curve, take the integral cubic 2-division model in `z=4x`, use every
bad rational prime of the public minimal model and every real embedding, and
compute the kernel of localization on the certified public point span.  This
is the exact *known-public strict-local Kummer kernel*; it is not a complete
strict Selmer group or an S-class-group calculation.

| curve | public mod-2 span | real embeddings of cubic field | local image rank | strict-local kernel |
| --- | ---: | ---: | ---: | ---: |
| 273 | 30 | 1 | 16 | 14 |
| 302 | 31 | 3 | 21 | 10 |

Both public spans therefore have a nonzero strict-local kernel, but the
supports, archimedean types and exact dimensions differ.  That is a shared
qualitative feature, not a matching strict-local fingerprint.  In particular,
the calculation does not justify treating either kernel as the full strict
class space, nor does it produce a cover or a new rational point.

A proof-certified full class-group computation for the reduced cubic field of
273 was also attempted as a separate check, but exceeded PARI's configured
1 GiB stack.  It is not used here, and no class-group comparison is inferred
from that resource failure.

## Replay

From `research/`:

```sh
sage -python elliptic-curves/cas/compare_curve273_curve302_parent_fingerprints.sage
```

The compact replay certificate is
[`curve273_curve302_parent_fingerprints_v1.json`](../../artifacts/generated-results/elliptic-curves/curve273_curve302_parent_fingerprints_v1.json).
