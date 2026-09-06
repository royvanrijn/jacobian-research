# Complete factor-of-13 scaling classification for the six compact R17 families

For every primitive integer parameter pair `(n,d)`, the displayed homogeneous
coefficients admit division by `13^4,13^6` exactly in one projective residue
cell per family. In that cell precisely one factor of 13 is removable from
the Weierstrass coordinates. A second is impossible. For nonsingular fibres,
the resulting model is minimal at 13 and still has bad reduction there.

The [aggregate proof](../../artifacts/generated-results/elliptic-curves/r17_13_scaling_geometry_v2.json)
combines exact polynomial identities, complete finite residue exclusions and
independent symbolic and exhaustive checks. This extends the earlier bounded
[local-feature audit](LOCAL_FEATURE_RANK25_2026-09-06.md) to all primitive
parameters of these six displayed families. It produces no new point or curve.

| Family | Unique nonminimal cell modulo 13 | Parameter matrix M |
| --- | --- | --- |
| 103b2 | d = 0 | (1, 0; 0, 13) |
| 11952 | n = 6d | (1, 5; -2, 3) |
| 074d9 | n = 10d | (3, 1; -1, 4) |
| 07ca9 | d = 0 | (1, 0; 0, 13) |
| 08234 | n = 2d | (2, -3; 1, 5) |
| 08f72 | n = 5d | (3, 2; -2, 3) |

Here `(n,d)^T = M (u,v)^T`. Every matrix has determinant 13 and its columns
span the stated index-13 residue lattice. A primitive old pair has primitive
new coordinates and avoids the one-dimensional kernel of M modulo 13.
Conversely, primitive `(u,v)` outside that kernel gives a primitive old pair:
any common divisor of the old coordinates divides `det M = 13`.

## Polynomial identities and generic group

Write the original degree-eight and degree-twelve binary forms as `A_h,B_h`.
For the displayed matrices the forms

```
A_new(u,v) = A_h(M(u,v)) / 13^4,
B_new(u,v) = B_h(M(u,v)) / 13^6
```

have integral coefficients. All 84 cells—14 projective cells in each family—
were checked, including every nonintegral outcome. Independent Sage symbolic
substitution verifies every binary-form identity.

On each corresponding homogeneous fibre, `x = 13^2 X`, `y = 13^3 Y`
gives the explicit Weierstrass isomorphism. Over the rational function field,
the invertible rational base change, with the usual degree-four and degree-six
homogenization of point coordinates, preserves the existing generic
rank-17 subgroup. These are new coordinate choices for existing families,
not new generic families or additional independent sections.

## Why the classification is complete

For an affine residue cell, any rational parameter with denominator prime to
13 gives a compatible sequence of residues modulo powers of 13. At depth k,
the audit retains every lift satisfying

```
A = 0 mod 13^min(k,4),    B = 0 mod 13^min(k,6).
```

The infinity cell uses the reversed polynomials and `d/n`. All cells other
than the integral one terminate with an empty residue set. The same test on
the divided forms rules out another factor of 13, after excluding only the
kernel corresponding to nonprimitive old pairs.

Although the declared limit permits depth six and 4,096 live residues,
every actual exclusion closes by depth three; no live set exceeds 13.
A second checker independently enumerates the complete relevant residue sets
modulo at most `13^3`, rather than trusting the lifting implementation.
There are no UNKNOWN branches in this result.

Since 13 divides neither 48 nor 864, the short-model invariants
`c4 = -48 A`, `c6 = -864 B` make these coefficient conditions necessary for
any integral model change with a removable factor of 13. Direct coordinate
scaling makes them sufficient. Exhaustive discriminant evaluation on every
eligible divided-model residue shows that the minimal model remains bad at
13. Thus there is no hidden good-reduction trace at 13 to restore in these
families, and no deeper power-of-13 model-scale population to search.

The version-two aggregate also checks the original discriminant in all 84
projective residue cells directly. Each is zero modulo 13. Outside the
unique nonminimal cell, the original model is already minimal at 13; inside
it, the once-divided model is minimal and bad. This supplies the original-model
part of the universal bad-reduction claim explicitly. The version-one
certificate and its exact scaling proof remain retained.

## Search implications remain limited

The integral charts suggest an exact way to organize parameter lattices, but
Euclidean reduction of such a lattice does not automatically reduce elliptic
coefficients. Among these six Gauss-basis square boxes, only the 074d9 chart
improves the conservative weighted coefficient bound, by about 4.2%; the other
five bounds worsen. The bounded signed-permutation checks are not a complete
base-automorphism classification.

No new parameter scan, prime-score extension or rational-point search was
launched from this audit. Any transformed population still needs its own
coefficient geometry, finite limits, exact curve deduplication and exposure
protocol. Global minimality at other primes, conductor, whole-curve rank and
near-record discovery remain separate questions.

Replay with `report_r17_13_scaling_geometry_v2.py --check`,
`verify_r17_integral13_charts.sage`, and
`verify_r17_13_scaling_classification.py`. Their source files, result JSONs and
supervision records are retained; an isolated portable replay for this new
supplement has not yet been claimed.
