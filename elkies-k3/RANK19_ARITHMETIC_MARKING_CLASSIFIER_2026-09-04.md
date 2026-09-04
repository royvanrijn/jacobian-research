# Arithmetic-marking classifier for rank-19 NS lattices

Date: 2026-09-04.

Status: **ACTIVE, fail-closed infrastructure**.

<!-- status-consumer: EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER 40745008c2fe2a80 -->

## Outcome

The rank-seven foundry catalogue now has a separate arithmetic pre-screen.
Among the 827 catalogue surfaces, 66 carry an exact even rootless rank-17
frame and enter the classifier.  The current exact decisions are

```text
ARITHMETICALLY_POSSIBLE    1
ARITHMETICALLY_EXCLUDED    2
UNKNOWN                   63
```

The positive row is the already-realized determinant-948 `NS0001` control.
The excluded rows are determinant-950 `NS0024` and determinant-1184 `NS0031`.
Every other row remains `UNKNOWN`; no bounded search, coarse modular curve, or
formal local branch is promoted.

The generated equation-agent handoff contains only new different-NS rows with
an exact `ARITHMETICALLY_POSSIBLE` certificate.  It is currently empty.  In
particular determinant `720` is sent to arithmetic curve identification, not
to Riemann--Roch equation work.

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

## The three controls

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

## Next arithmetic target

The determinant-720 row has

```text
A_NS = Z/2 x Z/6 x Z/60,
T = [ 8  -2   2 ]
    [-2  -4  10 ]
    [ 2  10  -4 ],
```

and split even-Clifford order of reduced discriminant `45`.  Its isotropic
vector has divisibility `3`; the order is not certified as `Gamma_0(45)`.
The next calculation is therefore the NS0031-style one: embed the order in
`M_2(QQ)`, derive its exact congruence conditions and signature, compute the
stable discriminant kernel for `Z/2 x Z/6 x Z/60`, and only then look for a
cheap quotient and rational points.

## Replay

```bash
sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage
sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage --check
python3 elkies-k3/scripts/build_rank7_determinant_aware_ranking.py --check
```

The checker does not reprove Momose, Vélu, the Inose correspondence, or the
rank-three period/spin correspondence.  Those inputs remain named in the
decision registry and canonical obstruction notes.
