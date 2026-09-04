# Arithmetic-marking classifier for rank-19 NS lattices

Date: 2026-09-04.

Status: **ACTIVE, fail-closed infrastructure**.

<!-- status-consumer: EC-K3-DET1236-GENUS2-RATIONAL-POINTS 5a3c84eb9f7f0604 -->
<!-- status-consumer: EC-K3-DET1236-MARKED-SHIMURA-CURVE e482668e1208f764 -->
<!-- status-consumer: EC-K3-DET1236-RATIONAL-CM-LOCUS bd6ab0e86ca70ab2 -->

<!-- status-consumer: EC-K3-GOLAY-DET720-QQ-MARKING-OBSTRUCTION 972f591d2885f9ba -->
<!-- status-consumer: EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER eec5710ee1b498ab -->
<!-- status-consumer: EC-K3-ARITHMETIC-FIRST-MARKED-T-FOUNDRY 9e9e0a1a8ac7c088 -->
<!-- status-consumer: EC-K3-DET500-DET750-QQ-MARKING-OBSTRUCTIONS 14498ad134ffa60e -->

## Outcome

The rank-seven foundry catalogue now has a separate arithmetic pre-screen.
Among the 827 catalogue surfaces, 66 carry an exact even rootless rank-17
frame and enter the classifier.  The current exact decisions are

```text
ARITHMETICALLY_POSSIBLE    1
ARITHMETICALLY_EXCLUDED    5
UNKNOWN                   60
```

The positive row is the already-realized determinant-948 `NS0001` control.
The excluded rows are determinants 500 and 750, determinant-720 Golay,
determinant-950 `NS0024`, and determinant-1184 `NS0031`.
Every other row remains `UNKNOWN`; no bounded search, coarse modular curve, or
formal local branch is promoted. One of those rows, determinant 1236, now has
the more precise Phase-2 certificate `UNRESOLVED_FOR_EXPLICIT_REASON`: its
exact marked curve is known, but its rational non-CM locus is not.

The generated equation-agent handoff contains only new different-NS rows with
an exact `ARITHMETICALLY_POSSIBLE` certificate.  It is currently empty.  In
particular no determinant-720 row is sent to Riemann--Roch equation work.

## Pipeline

For each exact rootless-MW17 candidate the classifier records:

1. `NS = U orthogonal_sum W(-1)` from every qualifying frame and checks rank,
   signature, determinant, evenness, minimum, and Smith invariants;
2. the exact rank-three primitive K3 complement `T=NS^perp` imported from the
   catalogue's opposite-discriminant-form construction;
3. the replayed rational even-Clifford algebra and integral order from the
   T-arithmetic ledger;
4. the abstract full-marking group
   `O^+(T)^* = kernel(O^+(T) -> O(A_T))`, separately from any coarse
   norm-one modular or Shimura curve;
5. certified easy quotient maps and rational-point or Frobenius tests when a
   decision record exists;
6. one of the three exact public classifications and the next arithmetic gate
   for every unknown.

The separation in step 4 is essential.  The existing T-arithmetic ledger may
replace a literal ternary form by a similar primitive integral quadratic
form.  This preserves the rational orthogonal group and is enough to compute
a coarse Clifford curve, but it can discard marking level.  Even when no
rescaling occurs, the stable discriminant-kernel subgroup still has to be
computed.  Thus a genus-zero `X_0(N)` row is not automatically an arithmetic
source.

The machine-readable outputs are

- [`../artifacts/generated-results/elkies-k3-rank19-arithmetic-marking-classifier-v1.json`](../artifacts/generated-results/elkies-k3-rank19-arithmetic-marking-classifier-v1.json),
- [`../artifacts/generated-results/elkies-k3-rank19-arithmetic-marking-equation-survivors-v1.json`](../artifacts/generated-results/elkies-k3-rank19-arithmetic-marking-equation-survivors-v1.json).

The small auditable decision registry is
[`data/arithmetic/rank19-arithmetic-marking-decisions-v1.json`](data/arithmetic/rank19-arithmetic-marking-decisions-v1.json).
New registered decisions, including exact-curve unresolved records, must name
exact certificate assertions and any external theorem inputs.

## The six controls

### Determinant 948: possible

The exact ternary lattice has division even-Clifford algebra and the clean
coprime squarefree Eichler pair

```text
(D,N) = (6,79),       reduced discriminant D*N = 474.
```

The quotient `X_0^6(79)/<w_474>` has genus two and model

```text
u^2 = 16*t^6 - 19*t^4 + 88*t^2 - 48.
```

It has the explicit non-CM rational point
`(t,u)=(14/13,64*251/13^3)`.  More importantly, the exact rootless equation,
seventeen rational sections, two-prime Picard-rank proof, saturation audit,
and torsion exclusion prove that this point really lifts to a full rational
rank-19 marking.  This is an existence certificate, not merely absence of an
obstruction.

### Determinant 950: excluded

Here `T = U + <950>` is split and the primitive `2E8/MW1` frame forces a
degree-475 cyclic Inose isogeny.  A full marking would give a noncuspidal,
non-CM rational point on `X_0^+(475)`.  Momose's theorem excludes it, using
the prime `19`.  The classifier imports the exact numerical certificate and
keeps the theorem input explicit.

### Determinant 1184: excluded

The split order gives the exact norm-one curve

```text
X_ns(4) x_{X(1)} X_0(37),       genus 23.
```

Forgetting level four maps to `X_0(37)`.  Vélu's rational-point
classification leaves two noncuspidal points, and both have Frobenius pair
`(trace,determinant)=(2,3) mod 4` at `19`, which is absent from the unramified
non-split Cartan.  Neither point lifts, including after quadratic twist.

### Determinant 720: excluded

For the literal transcendental lattice with
`A_T=Z/2+Z/6+Z/60`, the primitive similarity order has norm-one group
conjugate to `Gamma_0(15)`. Its spin action on `A_T` has image `S3`; imposing
the stable discriminant kernel adds the identity condition modulo `2` and
gives a group rationally conjugate to `Gamma_0(60)`. Thus the exact marked
curve is `X_0(60)`, of genus seven. The Mazur--Kenku cyclic-isogeny
classification excludes rational noncuspidal points. The rational `X_0(15)`
points do not lift, and the known rational `3A5` model instead saturates from
determinant `720` to determinant `20` with index six.

### Determinants 500 and 750: excluded

The two rootless rows have

```text
T_N = U(5) + <10N> = 5(U+<2N>),      N=2,3.
```

Their coarse curves `X_0(2)` and `X_0(3)` omit the literal projective mod-five
marking. The exact spin image is `PSL_2(F_5)=A5`; the Fricke and full
orthogonal images are `S5` and `S5 x C2`, and no non-spin coset is stable.
The full marked curves are the genus-four `X_H(50)` and genus-nine `X_H(75)`,
both degree two over `X_0(50)` or `X_0(75)`. Mazur--Kenku excludes rational
noncuspidal points on those quotients, while exact cusp Galois calculations
leave four rational cusps on each marked curve. Thus neither row has a K3
period point. See
[`DET500_DET750_QQ_MARKING_OBSTRUCTIONS_2026-09-04.md`](DET500_DET750_QQ_MARKING_OBSTRUCTIONS_2026-09-04.md).

### Determinant 1236: exact curve, rational lift unresolved

For

```text
T = [[-2,0,1],[0,4,0],[1,0,154]],       A_T=Z/1236Z,
```

the literal content-one Clifford order is the Eichler order with
`(D,N)=(6,103)`. Its Atkin--Lehner action exhausts the eight-element group
`O(A_T)`, and only the full Fricke class can be made stable after central
inversion. Thus the exact marked curve is

```text
X_0^6(103)/<w_618>,       genus 6, cusps 0, e2=12, e3=2.
```

It has exactly ten rational CM points—two of discriminant `-3` and four each
of discriminants `-43` and `-67`—all rank-20 specializations. It maps with
degree two to the explicit genus-two curve

```text
y^2 = 1944*x^6 + 441*x^4 - 90*x^2 + 9,
```

which has exactly fourteen rational points and a rank-one `618f1` quotient.
Bielliptic quadratic Chabauty and an exact Mordell--Weil sieve prove
the complete rational-point set, with elliptic images
`+/-G`, `+/-3G`, `+/-4G`, and `+/-10G`; it separates the two fixed fibers into
the rational discriminant-`-3` CM pair and a quadratic discriminant-`-24` CM
pair.  An exact characteristic-zero cover candidate now evaluates all twelve
non-fixed fibres: four fibres are rational for the displayed twist, while its
`-3` twist has none. The complete CM residue-field calculation proves that
all eight conditional lifts are CM and rules out the no-lift twist once the
candidate branch divisor is identified. Its Jacobian splits as the six rank-one
factors `618a1` through `618f1`; the genus-two part is `618e1 x 618f1` and
the cover Prym is `618a1 x ... x 618d1`. Classical Chabauty therefore misses
its strict rank bound, while the quadratic-Chabauty dimension screen passes
but still awaits a characteristic-zero identification of the candidate's
cubic CM branch orbit. That identification would give arithmetic exclusion,
not a positive handoff. The row therefore stays out of the equation handoff. See
[`DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md`](DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md).

## Why 948 looks hospitable

The current evidence points to a structural explanation, but not yet a
classification theorem.  Determinant `948` combines a clean Eichler order
with a very low-genus Atkin--Lehner quotient and an actual non-CM rational
lift.  The split-Clifford failures show that splitness itself is not
favorable: their markings force rigid modular quotients where global
rational-point theorems or a single Frobenius class rule out every lift.

The useful predictor is therefore not determinant size or split versus
division.  It is the tuple

```text
(integral order type, stable discriminant kernel,
 low-genus quotient, rational non-CM lift).
```

The last entry is decisive.  The first three only identify where to look.

## Next arithmetic targets

The construction search now starts from the full rank-three `T` ledger, not
from the remaining rootless-frame list. The generated arithmetic-first queue
contains all 827 transcendental rows and does not use rootless-frame data in
its priority. It propagates the classifier's five exclusions, the separate
split determinant-378 exclusion, one realized positive control, and an
820-row arithmetic research queue. Twenty-one rows have coarse genus at
most two, but each still needs its literal stable discriminant kernel and
rational-point decision.

The split determinant-378 row is now excluded because its literal stable curve
is `X_0(63)`, whose rational points are cusps. The next exact-coarse rows are
currently

```text
det 256: X_0(2),
det 512: X_0(4).
```

These labels remain coarse diagnostics. Compute each literal stable kernel,
then seek a rational noncuspidal non-CM point. Only an exact positive marking
may trigger `NS=T^perp`, rootlessness, or equation work. The full queue is in
[`../artifacts/generated-results/elkies-k3-arithmetic-first-marked-t-foundry-v1.json`](../artifacts/generated-results/elkies-k3-arithmetic-first-marked-t-foundry-v1.json).

Inside the already-rootless subqueue, determinant 1236 is now the strongest
near-exclusion target because its full marking curve, low-genus quotient
tower, complete fourteen-point rational locus, and complete ten-point
rational CM locus are exact. Its only next task is the modular-CM
identification of the exact cover candidate's branch orbit, not K3 equation
compilation.

## Replay

```bash
sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage
sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage --check
sage -python elkies-k3/scripts/certify_golay_det720_qq_marking_obstruction.sage --check
sage -python elkies-k3/scripts/certify_det500_det750_qq_marking_obstructions.sage --check
sage -- elkies-k3/scripts/certify_det1236_genus2_rational_points.sage --check
sage -python elkies-k3/scripts/certify_det1236_marked_shimura_curve.sage --check
python3 elkies-k3/scripts/build_rank7_determinant_aware_ranking.py --check
python3 elkies-k3/scripts/build_arithmetic_first_marked_t_foundry.py --check
```

The checkers do not reprove Momose, Vélu, Mazur--Kenku, the split-Eichler
normalizer theorem, the Inose correspondence, or the rank-three period/spin correspondence. Those inputs
remain named in the decision registry and canonical obstruction notes.
