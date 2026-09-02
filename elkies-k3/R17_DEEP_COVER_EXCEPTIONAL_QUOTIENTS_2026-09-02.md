# R17 low-genus covers at the rank-25--28 fibres

<!-- status-consumer: EC-K3-ELKIES-2026-LOW-GENUS-COVER-QUOTIENTS 31a6363906ad0ac0 -->

## Result

The tested degree-three and degree-four equation universes add no exceptional
direction beyond the complete degree-two bisection atlas.  In the ordered
public complements defining `L_t/M_t`, the cumulative captured ranks are

| parameter `t` | known exceptional rank | `R_t(2)` | `R_t(3)` | `R_t(4)` |
|---|---:|---:|---:|---:|
| `-2/377` | 8 | 5 | 5 | 5 |
| `-308/251` | 9 | 3 | 3 | 3 |
| `2456/135` | 10 | 2 | 2 | 2 |
| `-9529/5471` | 11 | 1 | 1 | 1 |

The rank-28 fibre was evaluated first.  Its sole captured direction is still

```text
(0,1,0,-1,1,0,0,-1,0,1,0)
 = Q2-Q4+Q5-Q8+Q10.
```

It is exposed by the degree-two orbit `orbit-15a68`, whose minimum trace norm
is 10.  Its deterministic equation-cost tuple begins with 11 group additions,
support 11, dependency closure 11, 3,751 serialized input bits, maximum
coefficient 2, and coefficient `l1` norm 12.  No tested degree-three or
degree-four cover splits there, so the other ten directions have no finite
minimum norm or equation cost in this tested universe.

The degree-four layer gains zero rank after the degree-three layer.  This is
the requested stopping condition, so degrees five and six were not run.  The
negative conclusion is bounded exactly as described below; it is not a
non-existence theorem for all higher-degree covers or for all of `M/4M`.

## Tested cover universe

The equation calculations use the compact published R17 model and its exact
17-section basis.

1. Degree two is the existing complete atlas of 39,120 rational bisection
   translation classes, all with minimum norm 10.
2. The degree-three frontier is a deterministic inversion-closed sample of
   1,025 cosets.  It contains 138 norm-20 rational vertices; one representative
   of each of the 69 inversion pairs was constructed.  Every residual cubic is
   irreducible over `QQ(t)`.  None has a rational component at any of the four
   controls.
3. The degree-three deep layer is complete: the full `3^17` coset census has
   exactly 320 norm-26 rational translation cosets, or 160 inversion pairs.
   All 160 residual cubics are irreducible over `QQ(t)`, with squarefree
   discriminant degree four.  None has a rational component at a control.
4. The degree-four layer uses the pinned deterministic inversion-closed sample
   of 1,025 cosets.  Its 106 norm-34 rational vertices give 53 inversion pairs.
   The residual factor patterns are 45 irreducible quartics and eight `2+2`
   decompositions.  The squarefree discriminant degrees are respectively six
   and four.  None has a rational component at a control.

Thus the new equation universe performs

```text
4 * (69 + 160 + 53) = 1,128
```

exact control factorizations and finds no new split cover.  The complete deep
degree-three statement is exact for all 320 deep translation cosets.  The
norm-20 frontier and degree-four statements are exact only for their declared
deterministic samples.  The complete lattice census contains 18,024,296
norm-20 degree-three cosets, and `M/4M` contains `4^17` cosets; those larger
equation universes were not constructed.

## Exact construction

For a norm-26 trace section

```text
tau = (Nx/h^2, Ny/h^3),  deg(h)=11,
```

the unique relation

```text
f0 + f1*x + f2*y + f3*x^2 = 0
```

has coefficient-degree bounds `(20,16,14,12)`.  The exact interpolation matrix
has size `65 x 66` and rank 65.  Eliminating `y` gives a quartic in `x`; division
by the known factor `x-x(tau)` leaves the certified residual cubic.  The
norm-20 construction is identical with bounds `(16,12,10,8)` and a rank-49
`49 x 50` matrix.

For a norm-34 degree-four trace, the relation

```text
f0 + f1*x + f2*y + f3*x^2 + f4*x*y = 0
```

has bounds `(25,21,19,17,15)`.  Its `101 x 102` interpolation matrix has rank
101.  Eliminating `y` and removing the trace factor leaves the residual
quartic.  Every stored relation, division, generic factorization, and control
factorization is exact over `QQ` or `QQ(t)`.

## Exact quotient coordinates and circuits

For every split degree-two cover, the specialized point was transported to
the pinned minimal and short models.  A full relation block was discovered
from canonical heights and then verified by exact rational group addition.
Solving that block gives integral coordinates in the ordered known basis, so
the quotient vectors in the certificate are exact elements of `L_t/M_t`, not
only finite-quotient fingerprints.

There are no parallel duplicates at any of the four controls.  The rank-25
fibre has one fundamental circuit in complexity order:

```text
orbit-1cb25 - orbit-0cff7 + orbit-1ea09 - orbit-0d4ca = 0
```

The rank-26, rank-27, and rank-28 displayed class sets are independent and
have no circuit.  The certificate also records every exact quotient vector,
the complexity-ordered rank increments, and the first tested cover whose span
contains each ordered public basis direction.  `null` in the latter table
means that direction is not individually exposed by the tested span.

## Reproduction

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/complete_lattice_foundry_degree3_spectrum.py \
  --frame-id NS0001-F001 --workers 8 --chunk-size 1000000 \
  --retain-norm 26 \
  --output artifacts/generated-results/elkies-k3-r17-degree3-deep-cosets-v1.json \
  --checkpoint artifacts/generated-results/elkies-k3-r17-degree3-deep-cosets-v1.json.partial

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/construct_elkies_2026_sampled_frontier_trisections.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/construct_elkies_2026_deep_trisections.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/construct_elkies_2026_sampled_quadrisections.sage

python3 \
  elliptic-curves/scripts/analyze_elkies_2026_deep_cover_quotients.py
```

The first command takes several minutes and visits all `3^17` degree-three
cosets using inversion.  The three equation builders and quotient analyser
support `--check` against their pinned outputs.

## Claim boundary

This is a computation about displayed cover classes and the known subgroup
`L_t`.  It does not prove that `L_t` is the full specialized Mordell--Weil
group, does not give a rank upper bound, and does not exclude a split cover in
the untested norm-20 or degree-four cosets.  The degree-five/six non-run is a
stopping decision forced by zero rank gain at two successive requested
layers, not evidence that all degree-five/six covers fail.
