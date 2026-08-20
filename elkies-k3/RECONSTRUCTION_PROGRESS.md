# Elkies rank-17 K3 reconstruction progress

Status checkpoint: 2026-08-20, consolidated after the E6/MW3 finite-field
breakthrough later the same day.

## Goal

Recover an explicit elliptic K3 fibration realizing the already-recovered
rank-17 Mordell--Weil lattice, together with explicit generic sections.  Once
that model exists, feed real specializations into the ignition/cascade search in
`RANK_GROWTH_SEARCH.md`.

This document separates three evidence levels:

- **exact**: certified by integer/lattice/algebraic computation in the repo;
- **numerical**: high-accuracy solutions of the reconstruction equations;
- **pending**: not yet exactified or proved.

## Exact starting point

The recovered rank-17 lattice is stored in `data/lattice/rank17_gram.txt` and
`data/lattice/short_vector_basis_gram.txt`.  The current identification work
matches the Elkies `X(6,79)` example: quaternion discriminant `D=6`, level
`M=79`; the corrected local criterion at the odd ramified prime `p=3` is
satisfied.  See `results/IDENTIFICATION-NOTES.md`.

The norm-4 shell is unusually rich:

- 1311 unsigned `+/-` pairs;
- 2622 signed minimal vectors;
- 184242 exact additive identities `a+b=c` among signed minimal vectors;
- a 17-vector seed generates the complete signed shell under these additive
  relations.

These files are the principal reconstruction data:

```text
data/lattice/rank17_gram.txt
data/lattice/short_vector_basis_gram.txt
data/relations/all_2622_signed_short_basis.npy
data/relations/minimal_additive_triples.npy
```

Public Elkies lecture material gives the working model expected for the
reconstruction: after a suitable change of variables one may work with

```text
y^2 = x^3 + A(T) x + B(T)
```

with `deg A <= 8`, `deg B <= 12`; minimal height-4 sections have quartic
`x(T)` and sextic `y(T)`.  The canonical final representation is still to be
stored under `data/k3-model/` once exact coefficients are recovered.

## Consolidated primary route: 17-by-17 frame to E6/MW3

The direct Coxeter continuation described below remains a useful diagnostic,
but it is no longer the primary reconstruction route. Exact neighbor searches
from the 17-by-17 positive frame have produced lower-Mordell--Weil-rank
fibrations with larger reducible root systems. The most promising node has

```text
ADE = E6 + A3^2 + A1^2
MW rank = 3
reduced MW determinant = 79/16
fibers = IV* + I4 + I4 + I2 + I2 + 4 I1.
```

Putting `IV*` at infinity reduces the short Weierstrass degrees to
`deg A <= 5`, `deg B <= 8`. Exact triangular elimination and a rational
parametrization reduce the one-section locus to eight fiber equations in
eleven variables, of expected dimension three.

The first complete point on this locus has been reconstructed exactly over
`GF(31)`. Its full section identity and local fiber conditions vanish, and the
reduced `8 x 8` Jacobian is nonsingular. Thus it is a smooth modular seed with
a convergent local Hensel chart. It is not yet the target rank-3 seed: the
first polynomial companion satisfies `P1=-2*P3`, and the canonical
one-denominator `P2` is absent on that first surface.

The exhaustive reduced-core scan is nevertheless finite and concrete. It
tested all 893,730 points on the declared branch and left seven split-root
cores. A degree-parity argument forces the remaining `P2` numerator to have
the form `X2=C+q0*F`; hence each fixed surface needs only 837 `(r,q0)` tests.
See [`E6_P2_REDUCTION_2026-08-20.md`](E6_P2_REDUCTION_2026-08-20.md).

## Why direct reconstruction was reduced

The naive symbolic system (`build_rank17_system.sage`) introduces `A`, `B`, and
17 quartic/sextic sections simultaneously and is too large for direct Groebner
work.

The recovered shell contains a 9-vector clique with Gram

```text
4 on the diagonal, 2 off the diagonal.
```

Write its oriented sections as `V_0,...,V_8`.  Since
`<V_i,V_j>=2`, every difference

```text
D_ij = V_i - V_j
```

is again minimal.  Thus the clique immediately gives 45 known minimal
x-coordinate classes: 9 vertices plus 36 differences.

For every minimal additive triple `P+Q=R`, the chord slope is quadratic and

```text
x_P + x_Q + x_R = m_PQR^2.
```

For the Coxeter subsystem this gives exactly 120 square relations on 45
quartic x-polynomials.  `build_coxeter9_x_reconstruction.py` constructs the
120x45 incidence matrix and verifies every relation against the global 184242
triple catalog.

Observed exact/numerical rank data:

```text
incidence_shape = 120 x 45
numeric_rank = 45
rank mod 5,7,11,... = 45
exact full-column rank over Q = certified
left nullity = 75
```

Hence all 225 x-coordinate coefficients can be eliminated linearly once the
quadratic slopes are known.

## Coherent-slope reduction

The 84 triangle slopes are not independent.  If `m_ij` is the slope through
`V_i` and `-V_j`, then

```text
slope(D_ij,D_jk) = m_ik - m_ij - m_jk.
```

Therefore the nonlinear reconstruction can be expressed using only the 36 pair
slopes, i.e. 108 scalar coefficients before gauge fixing.  The numerical gauge
fixes four coefficients, leaving 104 free Coxeter variables.

Given the slopes, the implementation reconstructs in sequence:

```text
36 pair slopes
    -> 9 quartic x_i
    -> 9 sextic y_i
    -> one common degree-8 A(T)
    -> one common degree-12 B(T).
```

`solve_coxeter9_slopes_numeric.py` repeatedly finds roots with algebraic
residuals around `1e-13` to `1e-16`.  Many roots lie on or close to the
cuspidal/isotrivial component, so residual size alone is not a useful ranking
criterion.

## Coxeter-root diagnostics

`analyze_coxeter9_numeric_roots.py` ranks roots using:

- scale-aware `A/B` and discriminant strength;
- effective degrees of `A`, `B`, and the discriminant;
- the polynomial non-isotriviality invariant

  ```text
  3 A' B - 2 A B';
  ```

- numerical Jacobian nullity.

Two roots became especially important:

### `root-000001`

```text
raw ~ 6.97e-15
non-isotrivial
(deg A, deg B, deg Delta) = (8,12,24)
nullity at 1e-8 ~ 9
```

This is an excellent generic Coxeter-scaffold point.  However, attempts to
impose the first rank-10 extension on its initial lattice labeling drove the
surface strongly toward the isotrivial boundary.

### `root-000029`

```text
raw ~ 8.14e-13
non-isotrivial
(deg A, deg B, deg Delta) = (8,12,24)
discriminant diagnostic ~ 0.18
```

This root became the best testbed for the first extension, but the first
apparently machine-precision rank-10 hit on it was later found to be a
near-collision and is not yet an independent-section certificate.

## Exact rank-17 extension chain

Starting with the raw 9-vector Coxeter clique, the combinatorial extension
search selected the following signed minimal vectors:

```text
2313 2525 307 1303 1859 1441 683 2351 2143
961 2402 1642 1300 1023 2216 2392 2610
```

The final matrix has:

```text
rank = 17
coordinate determinant = 2
77 pairings with absolute value 2
additive closure = 1396/2622 before explicit saturation bridging
```

The determinant `2` is **optimal**, not a failure.  The raw Coxeter clique has
Gram determinant 5120 while its saturated rank-9 lattice has determinant 1280,
so the raw clique itself has saturation index

```text
sqrt(5120/1280) = 2.
```

Any 17-vector coordinate matrix retaining those nine raw generators must
therefore have determinant divisible by 2.  The selected chain attains this
minimum, and its saturation is the complete recovered rank-17 lattice.

`audit_rank17_extension_chain.py` exposes the index-2 parity class, searches a
minimal saturation bridge, and exports exact `|pairing|=2` continuation anchors.

## Rank-10 continuation: first attempts

The first continuation used section `961` with only its `|pairing|=2` links to
the nine raw Coxeter generators.  It found excellent algebraic fits, often
`1e-9` or better, but these collapsed toward the isotrivial boundary.  Adding
all nine available anchors against the full 45-section Coxeter shell improved
the algebraic fit further but showed a clear Pareto split:

- algebraic residual near `1e-10`/`1e-11` -> discriminant and j-variation
  collapse by many orders of magnitude;
- healthy non-isotrivial branch -> algebraic residual only around `1e-5`.

The line/group-law formulas were checked; the remaining issue was not a sign
bug but a symmetry/labeling ambiguity and, as found later, a collision escape.

## Coxeter S9 ambiguity

The raw clique Gram is invariant under permutation of its nine generators.
Thus a numerical Coxeter solution has no intrinsic labeling telling us which
numerical `V_i` corresponds to which full rank-17 lattice generator.

For lattice section `961`, the original generator-pairing signature is

```text
-1 0 -1 0 -1 -2 -2 -2 -1
```

Its multiset has 1260 distinct `S9` permutations; including the opposite
orientation of the new section gives 2520 distinct fingerprints.  These were
scanned exhaustively on fixed Coxeter roots by
`scan_rank10_coxeter_fingerprints.py`.

### Root-000001 scan

Best fixed-surface result:

```text
score ~ 4.412e-7
```

No compelling rank-10 hit appeared.

### Root-000029 scan: apparent hit, later downgraded

The exhaustive scan found:

```text
fingerprint index = 2452
V pairings = 0 0 2 2 2 1 1 1 1
anchors = 9
curve residual = 5.220e-14
line residual = 5.850e-13
score = 5.850e-13
```

The winning fingerprint is a permutation of the **negative** of the original
section-961 signature.  This left 288 compatible full-lattice Coxeter mappings.

A joint refinement then reached

```text
base residual  ~ 3.36e-13
curve residual ~ 7.73e-16
line residual  ~ 4.11e-14
healthy discriminant / non-isotrivial diagnostics
```

However, `validate_rank10_independence.py` subsequently found

```text
nearest_x = V3
x_distance ~ 2.52e-7
```

and the intended quadratic-line identities were not numerically separated from
many accidental non-target line fits.  Six intended edges had absolute errors
around `3e-14` to `4e-14`, but several non-target sections had comparable
absolute errors.  The original validator's strict relative threshold therefore
marked only 3 of 9 target edges; simply relaxing that threshold would create
many false extra edges.

The key conclusion is that the original `5.85e-13` fixed-surface hit and its
joint refinement are **near-collision candidates, not yet evidence of an
independent tenth section**.  Their small residuals can be explained by the
ill-conditioned gauge together with `Q` approaching an existing Coxeter
section.

The rank-10 tangent/Jacobian experiments are also not being used as evidence at
this checkpoint.  Several versions showed that strongest-SVD-gap nullities are
parameterization- and scaling-sensitive, and an incremental-codimension test
was initially invalid because numerical tangent leakage was amplified by
column normalization.  These diagnostics remain useful for conditioning work,
but they are secondary to direct section/incidence validation.

## Current experiment: nondegenerate rank-10 refinement

`refine_rank10_nondegenerate.py` now repeats the winning-fingerprint refinement
with explicit structural guards:

1. `Q` must stay a scale-free positive distance from all 45 known Coxeter
   sections, using the median pairwise Coxeter x-distance as normalization;
2. the discriminant and non-isotriviality branch must remain healthy;
3. all nine required target line identities are imposed;
4. after solving, required edges must be numerically separated from all
   non-target and wrong-sign quadratic-line fits.

The principal threshold-free validation statistic is

```text
edge_gap = min(non-target line error) / max(required-edge line error).
```

The previous near-collision solution has edge gap of order one.  The new solver
requires, by default, both target/non-target and correct/wrong-sign gaps of at
least `1e2`, in addition to the collision and branch guards.

A successful result should therefore simultaneously have:

```text
small algebraic residual
healthy discriminant / j-variation
distinct_ratio >= 0.02
edge_gap >= 1e2
wrong_sign_gap >= 1e2
```

If no such solution exists for the current fingerprint/root, the next step is
to rescan the 2520 Coxeter fingerprints with the same non-collision and
incidence-separation criteria rather than interpreting the old machine-small
residual as rank growth.

## Next steps

Primary E6/MW3 route:

1. Reconstruct every admissible `P1` chart over the seven complete split-root
   `GF(31)` cores and both `lambda/mu` labelings.
2. Exhaust the reduced 837-case `P2` search on each reconstructed surface.
3. Reject all group-law dependencies against `P1` and `P3`; retain only a
   genuinely independent second section.
4. Lift a successful modular seed at several good primes and perform rational
   or p-adic coefficient recognition.
5. Verify the characteristic-zero Weierstrass identities and Kodaira fibers
   exactly, then use neighbor transformations toward the rootless rank-17
   fibration.

Fallback direct-Coxeter route:

6. If the seven-core E6 tranche is empty, rescan the 2520 Coxeter fingerprints
   with collision guards and edge-gap scoring rather than reviving the old
   near-collision.
7. Use the index-2 saturation bridge as a gluing constraint for any genuine
   continuation through sections
   `961,2402,1642,1300,1023,2216,2392,2610`.

Final certification remains unchanged: exactify all coefficients, reconstruct
the complete 1311-pair shell, store the canonical model under `data/k3-model/`,
and only then resume specialization searches.

## Current claim boundary

What is exact today:

- the rank-17 lattice and its 1311-pair minimal shell;
- the 184242 additive relations;
- the Coxeter-9 decomposition and incidence-rank reduction;
- the coherent-slope algebra;
- the determinant-2 optimal rank-17 extension chain and its saturation
  interpretation.

What is numerical today:

- reconstruction of many Coxeter-9 elliptic K3 points;
- the non-isotrivial `root-000029` scaffold;
- a very small-residual **near-collision** candidate for the first extension,
  now explicitly downgraded pending nondegenerate validation.

What is exact over a finite field today:

- a complete E6/P1 surface and section over `GF(31)`;
- nonsingularity of its reduced `8 x 8` Jacobian;
- dependency `P1=-2*P3` on that residue branch;
- completeness of the seven split-root reduced cores and the degree reduction
  `X2=C+q0*F` for the canonical `P2` search.

What remains unproved:

- an independent generic tenth section realizing the intended fingerprint;
- an exact Weierstrass model for the intended rank-17 fibration;
- exact generic coordinates for sections 10 through 17;
- any new rank-21 specialization arising from this family.
