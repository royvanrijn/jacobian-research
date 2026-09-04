# Arithmetic-marking classifier for rank-19 NS lattices

Date: 2026-09-04.

Status: **ACTIVE, fail-closed infrastructure**.

<!-- status-consumer: EC-K3-GOLAY-DET720-QQ-MARKING-OBSTRUCTION 972f591d2885f9ba -->
<!-- status-consumer: EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER 93e6c5626d369572 -->
<!-- status-consumer: EC-K3-ARITHMETIC-FIRST-MARKED-T-FOUNDRY d94b3dbddf5cb529 -->

## Outcome

The rank-seven foundry catalogue now has a separate arithmetic pre-screen.
Among the 827 catalogue surfaces, 66 carry an exact even rootless rank-17
frame and enter the classifier.  The current exact decisions are

```text
ARITHMETICALLY_POSSIBLE    1
ARITHMETICALLY_EXCLUDED    3
UNKNOWN                   62
```

The positive row is the already-realized determinant-948 `NS0001` control.
The excluded rows are determinant-720 Golay, determinant-950 `NS0024`, and
determinant-1184 `NS0031`.
Every other row remains `UNKNOWN`; no bounded search, coarse modular curve, or
formal local branch is promoted.

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
New non-`UNKNOWN` decisions must name exact certificate assertions and any
external theorem inputs.

## The four controls

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

## Why 948 looks hospitable

The current evidence points to a structural explanation, but not yet a
classification theorem.  Determinant `948` combines a clean Eichler order
with a very low-genus Atkin--Lehner quotient and an actual non-CM rational
lift.  The two nearby failures are split-Clifford cases, so splitness itself
is not favorable: their markings force rigid modular quotients where global
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
its priority. It propagates the classifier's three exclusions, the separate
split determinant-378 exclusion, one realized positive control, and an
822-row arithmetic research queue. Twenty-three rows have coarse genus at
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

## Replay

```bash
sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage
sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage --check
sage -python elkies-k3/scripts/certify_golay_det720_qq_marking_obstruction.sage --check
python3 elkies-k3/scripts/build_rank7_determinant_aware_ranking.py --check
python3 elkies-k3/scripts/build_arithmetic_first_marked_t_foundry.py --check
```

The checkers do not reprove Momose, Vélu, Mazur--Kenku, the Inose
correspondence, or the rank-three period/spin correspondence. Those inputs
remain named in the decision registry and canonical obstruction notes.
