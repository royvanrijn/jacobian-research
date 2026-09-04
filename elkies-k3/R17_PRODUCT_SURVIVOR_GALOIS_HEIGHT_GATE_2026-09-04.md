# Product-survivor Galois height gate and constructor priority

**Update (2026-09-05):** the former first target `19bad:083ad` has
[arithmetic rank zero](R17_PRODUCT_19BAD_083AD_ARITHMETIC_RANK_ZERO_2026-09-05.md).
The regulator squareclasses from its existing reductions at 131 and 137
are incompatible. Its point-Kummer and Tate quotients are zero, and both
requested height boxes, as well as all higher heights, are empty.
The original five-target bounds and constructor metrics below are retained
as source evidence; the active arithmetic worklist now has four targets.

## Result

The complete two-prime product-twist classification leaves exactly five
targets.  At each of their ten irreducible quadratic branch places, the
original smooth fibre's two-division cubic remains irreducible over the
quadratic residue field.  The three nonidentity component classes of the
resulting `I0*` fibre are indexed by the three nonzero two-torsion points.
Therefore none is fixed by residue-field Galois, and a `QQ(u)`-rational
section must meet the identity component at all four `I0*` fibres.

For every rational section on any of these five product twists, the Shioda
height formula consequently reduces to

```text
height(P) = 8 + 2(P.O).                             (1)
```

In particular, a rational height-eight product-character section is
necessarily disjoint from the zero section and lies in the direct polynomial
box

```text
deg X <= 8,       deg Y <= 12.                      (2)
```

The `P.O=1,2` height-eight cases allowed by the unspecialized geometric height
inequality do not occur over `QQ(u)` on these targets.  This is a Galois
component obstruction, not another bisection-carrier scan.

The prior norm-eight and deep-trace inversions exclude the zero Tate class in
(2).  Hence any height-eight section that exists on these five targets has a
nonzero Tate class.  The first genuinely higher-pole search is instead

```text
P.O=1,  height=10,  deg numerator(X)<=10,  deg numerator(Y)<=15,
```

with a linear pole denominator.  A `P.O=2` section has height twelve and the
corresponding bounds `(12,18)` with a quadratic denominator.  Neither layer
has been searched here.

There is also an arithmetic rank gate not visible from the geometric Tate
degree alone.  At `p=137`, every survivor's complete normalized elliptic
quotient has Tate factor exactly

```text
(Z-1)(Z+1).
```

A `QQ(u)`-rational section specializes to a Frobenius-fixed divisor class, so
the normalized `Z=1` multiplicity proves

```text
0 <= rank E^(d)(QQ(u)) <= 1
```

for each survivor, while the geometric bound remains at most two.  Thus one
constructed nonzero product-character section would certify arithmetic rank
exactly one, the precise missing contribution in the `17+1+1+1` objective.

The exact certificate is
[`elkies-k3-r17-product-survivor-galois-height-gate-v1.json`](../artifacts/generated-results/elkies-k3-r17-product-survivor-galois-height-gate-v1.json),
SHA-256
`d6f5c2bb74531506b2eecdcec6462c44ec79addc5c2fa8e9c12bb7a3cad7fe85`.

## Constructor-aware target order

All five genus-one bases have exact Jacobian rank one and degree-four maps to
the alternate parameter `u`.  Ranking only by the Jacobian generator height
misses a large finite-height effect in the map.  The exact starting metrics
are:

| product pair | generator height (approx.) | bits in `u(P)` | bits in `u(2P)` |
|---|---:|---:|---:|
| `19bad:083ad` | 136.090 | 124 | 913 |
| `11ee2:0c36e` | 90.018 | 635 | 1,667 |
| `0c10b:17a1a` | 93.417 | 654 | 1,735 |
| `0f82c:025be` | 111.935 | 779 | 2,067 |
| `11ae6:0f82c` | 114.531 | 792 | 2,110 |

The bit counts are exact; the canonical heights are the stored Sage
approximations. Thus `19bad:083ad` was the primary explicit-construction
target because its first integer parameter is dramatically smaller.
Explicitly, the first image is

```text
u(P) = 15723338214416440692688694148490520697
       / 4764204189425575887883321058627395288.
```

`11ee2:0c36e` is the secondary target because its generator has the smallest
canonical height and therefore the best asymptotic growth among the five
degree-four maps.  The certificate also records exact `u(nP)` bit sizes for
`1 <= n <= 4` and the rational values themselves; no rational-base claim is
made.  These remain genus-one
integer-parameter constructors

```text
n -> nP -> u(nP) -> E_(u(nP)),
```

not `QQ(v)` families.

## Third-prime test of the primary target

The complete degree-28 product-twist Frobenius quotient for
`19bad:083ad` was independently recomputed at `p=151`.  A separate fibrewise
`n=1,2` character-sum audit agrees with the toric output.  The normalized Tate
factor is again

```text
(Z-1)(Z+1),
```

exactly as at `p=131` and `p=137`.  The unconditional conclusions are

```text
0 <= rank E^(d)(QQbar(u)) <= 2,
0 <= rank E^(d)(QQ(u)) <= 1.                       (3)
```

The persistent unique `Z-1` direction is consistent with a rational section but does
not produce one; the `Z+1` direction is likewise only a reduction pattern.
The certificates are
[`alternate-orbit-19bad--alternate-orbit-083ad-p151-v1.json`](../artifacts/generated-results/elkies-k3-r17-product-extra-prime-audits/alternate-orbit-19bad--alternate-orbit-083ad-p151-v1.json)
and
[`elkies-k3-r17-product-alternate-orbit-19bad--alternate-orbit-083ad-p151-toric-frobenius-v1.json`](../artifacts/generated-results/elkies-k3-r17-product-alternate-orbit-19bad--alternate-orbit-083ad-p151-toric-frobenius-v1.json),
with SHA-256 values
`f5a0c13f2162efb121caf7cf4b47bba55bbafd4ff2cd7f5de9069262322d7df9`
and
`86fd05853673d2518e0c3d869e131b95ff11053e6d10ae2dc95c7deb9b9152a0`.

## Historical next gate (superseded for `19bad:083ad`)

The original handoff was to compute the target-specific two-Selmer/Kummer quotient for `19bad:083ad` and
retain only the nonzero classes modulo the restricted invariant `R17` span.
Then solve the class-sliced form of (2).  If that finite collection of
projectively complete systems is empty, move to the height-ten linear-pole
box; do not enlarge the exhausted integral coboundary atlas.

The later rank-zero theorem closes this target without those solver runs.
For any of the other four targets, a nonzero rational solution would close
the arithmetic product rank at exactly one by (3).

The multisection geometry in Garbagnati--Salgado motivates this use of
suitably ramified covers.  Its hypotheses and Picard-rank-three theorem are
not being applied to these Picard-rank-nineteen surfaces, and it supplies no
rank jump here.

## Replay

```bash
sage -python \
  elkies-k3/scripts/certify_r17_product_survivor_galois_height_gate.sage
sage -python \
  elkies-k3/scripts/certify_r17_product_survivor_galois_height_gate.sage \
  --check

sage -python \
  elkies-k3/scripts/audit_r17_product_twist_extra_prime.sage \
  --pair-key 'alternate-orbit-19bad:alternate-orbit-083ad' \
  --prime 151
sage -python \
  elkies-k3/scripts/audit_r17_product_twist_extra_prime.sage \
  --pair-key 'alternate-orbit-19bad:alternate-orbit-083ad' \
  --prime 151 --check

elkies-k3/scripts/run_r17_product_toric_frobenius_extra_prime.sh \
  'alternate-orbit-19bad:alternate-orbit-083ad' 151
```

The extra-prime runner uses the same pinned toric pipeline as the original
runner, first creates the displayed independent moment audit, and passes that
audit explicitly to the final verifier.

## Claim boundary

The ten residue-field factorizations, rational component conclusion, formula
(1), degree boxes, arithmetic rank-at-most-one bounds, constructor metrics,
and the full `p=151` Frobenius quotient are exact.  No product-character
section, full two-Selmer group, nonzero Kummer representative, or
characteristic-zero Picard-rank-twenty class has been constructed.
Higher-height sections remain `UNKNOWN` for the other four targets;
`19bad:083ad` is excluded at every rational height by the linked follow-up.

<!-- status-consumer: EC-K3-R17-PRODUCT-19BAD-083AD-ARITHMETIC-RANK-ZERO fe572bd5979b5d2c -->
