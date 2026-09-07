# Historic +12 controls from retained equation factors

Both historic +12 rows now have complete masked generic-block and local
boundary computations. Curve 356 requires at least **11 additional strict
rational directions**, whereas the present constraints on curve 385 require
at least **4**. The arithmetic uses only their equations, marked generic
sections and independently checked factor hints. Existing rank labels are
joined afterward.

The [panel accounting note](FRESH_STRICT_BLOCK_NECESSITIES.md) incorporates
these rows, taking coverage from nine to eleven of the sixteen frozen
fibres. This supplement documents the new input certificates and the
distinct implications of the two historic cases.

## What was recovered and verified

The [protocol](FRESH_RETAINED_FACTOR_PROTOCOL.json) projects only:

* curve 356's discriminant-prime list from
  `record_prime_factor_proofs_20260904.json`;
* curve 385's `factor_hint_primes` from the historical input embedded in
  `exceptional_selmer_feasibility_v1.json`;
* retained local primes and remaining cofactors for the unfinished fresh
  rows, plus 2 and 3 for the change to the frozen integral short model.

The [hint payload](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_retained_factor_hints_v1.json)
contains only tokens, hint integers and provenance. It does not export
exceptional coordinates, Kummer masks, CT matrices or rank labels. The
curve 385 source records a **failed** class-group attempt; its failure
supplies no class-group information. Its input integers remain usable
as hints because the new worker verifies them independently.

For each frozen monic cubic f, the worker proves primality of usable hints,
extracts their exact valuations in 16 disc(f), and checks the complete
product. Both historic products close: 13 primes for curve 356 and 16
for curve 385. This verifies coverage in the current integral model,
including any scale primes, without assuming that old displayed-model
exponents can be copied unchanged. No general factorization routine runs.

The five remaining rows retain verified composite cofactors:

| Frozen case | Role | Remaining cofactor bits |
|---|---|---:|
| new-40, 074d9 at 2818/1535 | +10 | 364 |
| new-72, 11952 at 2012/211 | +10 | 387 |
| new-186, 11952 at 4286/1881 | +10 | 352 |
| new-90, MW16 at −1867/270 | +11 | 387 |
| MW16 at −3187/3697 | observed +1 control | 429 |

These sizes are coverage metadata, not rank features. The experiment did
not factor those residuals or replace the cases.

## The historic block comparison

Use the same definitions as the panel: m=dim G, k=dim(G∩U), c_S=dim U,
ell the local point-product dimension, and a=(ell−1)−(m−k).
On each completed row a nonzero equation-defined derivative functional
gives

\[
 \dim\operatorname{Sel}_2(E)/G=(c_S-k)+e,
 \qquad 0\le e\le a.
\]

| Quantity | Curve 356 | Curve 385 |
|---|---:|---:|
| Retained rank lower bound / marked generic rank | 29 / 17 | 29 / 17 |
| Recorded gain | +12 | +12 |
| Inherited strict dimension k | 1 | 0 |
| Inherited Artin matrix | [1] | empty |
| Inherited −1-twist CT-switch rank | 0 | 0 |
| Local point-product dimension ell | 18 | 26 |
| Additional boundary capacity a | 1 | 8 |
| Additional strict rational dimension forced by rank ≥29 | ≥11 | ≥4 |
| Total strict rational dimension, hence c_S, forced by rank ≥29 | ≥12 | ≥4 |

For curve 356, the cubic field has one real root. The derivative class
beta=−disc(f)f'(theta) lies outside the full generic local point image at
**2**, certified in local squareclass coordinates. Local self-duality
therefore gives a nonzero reciprocity functional. For curve 385, the
three-real-root sign calculation supplies the witness at infinity.
Both derivatives have norm disc(f)^4 and verified even valuations at
all omitted good polynomial-discriminant primes.

For curve 356 the exact interval is

\[
 \dim\operatorname{Sel}_2(E)/G=(c_S-1)+e,\quad0\le e\le1,
 \qquad \operatorname{rank}E\le c_S+17.
\]

Thus at least eleven of its twelve recorded quotient directions are
strict modulo G. At least twelve independent strict rational classes
exist in total. Their quadratic character compositum over K has degree
at least 4096, but this is an encoding of the forced classes, **not a
geometric carrier explaining their rational solubility**. No independent
twelve-class basis was constructed.

For curve 385,

\[
 \dim\operatorname{Sel}_2(E)/G=c_S+e,\quad0\le e\le8,
 \qquad \operatorname{rank}E\le c_S+25.
\]

This forces only four strict directions. It does not prove that the
actual number is four or that the boundary contribution reaches eight.
The two controls therefore have different *permitted decompositions*,
not a measured difference of seven in their actual strict rational ranks.

The inherited Artin value [1] on curve 356 is compatible with its zero
alternating switch: A+A^T=[0] over F2. One scalar Artin entry is neither
an extra rational point nor a nonzero alternating CT obstruction. Inherited
CT again provides no explanation of the exceptional quotient.

## Consequences and remaining gates

The historic +12 comparisons sharpen the same lesson as the fresh panel:
one cannot identify a successful jump block by jump size alone, and counting
the inherited structure misses the relevant excess. The meaningful
**incidence** target is c_S−k together with the actual boundary contribution
e. Local constraints bound the latter; they do not determine it.

Independent upper certificates c_S≤12 for curve 356 or c_S≤4 for curve 385
would match the respective retained rank ≥29 bounds, prove exact rank 29,
and force Sha[2]=0. These are concrete conditional targets, not obtained
upper bounds. A class computation below those thresholds would contradict
the retained rank certificates and must be investigated as an error.

The next **solubility** gate remains an independently constructed additional
class basis and CT information on that basis, followed by a sufficient
global rational-solubility criterion. The inherited block, class-dimension
necessities, and vanished inherited switches do not provide it. The fresh
MW16 +11 comparison is still incomplete and must not be represented as a
negative result.

## Reproducibility

The [supplement](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_retained_factor_supplement_v1.json)
retains factors, generic Kummer classes, half ideals, Artin evaluations and
the boundary witnesses. The [label-joined comparison](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_retained_factor_comparison_v1.json)
is separate. Existing workers are imported with an isolated checkpoint
directory; none of their files or earlier artifacts is edited.

The [verification](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_retained_factor_verification_v1.json)
passes: all prime products and composite residuals, both 17-dimensional
generic independence checks, full generic local matrices, 30 separate local
power tests for curve 356's strict representative, its exact half-ideal
square and Artin reduction/Jacobi calculation, derivative norm/support,
and the rank-label deductions. Coordinate replay uses the same local
character backend, while local-power and integer Jacobi checks use distinct
interfaces or arithmetic.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_fresh_retained_factors.py check
```

No exceptional point was supplied to a worker. No new point search,
parameter sweep, class-group campaign, or active-search modification ran.
