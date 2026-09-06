# MW16 at 3/17: class-2-rank 16 and curve rank 22 under GRH

Mathematical authority: `MATH_STATUS.json`, entries `EC-SMALL-CONDUCTOR-CLASS-TARGET-20260906`, `EC-SMALL-CONDUCTOR-CLASS-LOWER16-20260906` and `EC-SMALL-CONDUCTOR-CLASS-CHARACTERS-20260906`.

**Current class-2-rank bounds: 16 <= g <= 16.** The lower bound is unconditional;
the upper bound depends on GRH for the stated ideal-class characters.
The unconditional curve-rank lower bound is **22**. Matching bounds prove **exact rank 22 under GRH**.

The [preceding descent study](SMALL_CONDUCTOR_DESCENT_SHORTCUT_2026-09-06.md) establishes
the cubic field, local correction `rank <= dim Sel_2 <= g+7`, proved even
Selmer parity and the interval-certified generating cutoff 37,638 under GRH.
Thus `g<=16` suffices for rank 22 under the same assumption.

The [completion proof](SMALL_CONDUCTOR_CLASS_COMPLETION_PROOF_2026-09-06.md)
closes the bound using a quadratic-character explicit formula. The formal
relation quotient remains 18; exact memberships of small prime classes in
the sixteen-anchor span suffice to prove that span generates under GRH.
At cutoff 50,000 the worst-case corrected margin exceeds 17.16. A second
interval calculation with a weaker archimedean bound is also positive.
The [certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_class_completion_v1.json) lists every membership
and unresolved prime. No twenty-third point or missing norm relation is claimed.

## Unconditional lower bound and independent ideal anchors

The [lower-bound certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_class_lower16_v1.json) uses the known points
to form `beta_i=4*x(P_i)-theta`, with norm `(8*y(P_i)+4*x(P_i))^2`.
Residue characters at a fixed set of 128 rational primes prove independence
of all 22 field square classes. Their valuation-parity matrix at every bad
prime has rank 4. Its kernel supplies 18 independent products with even
valuations everywhere: away from `2*Delta(E)`, separability and the square
norm prove the assertion, including when x has a pole. Adjoining -1 adds
an independent class because its cubic norm is negative. The field Selmer
group has dimension at least 19; the unit square-class dimension is 3.
The [field Selmer exact sequence](https://arxiv.org/html/1606.07178#S4.SS3)
therefore proves **g >= 19-3 = 16** unconditionally.

The [character certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_class_characters_v1.json) imposes dyadic unit
Hilbert symbols and positive signs at the three real embeddings on these
19 classes. The constraint rank is 3, leaving 16 independent everywhere-
unramified quadratic extensions. Each dyadic residue field is F2; the
units `1+pi^k`, `1<=k<=2e`, generate units modulo squares because successive
unit-filtration quotients are F2 and `U^(2e+1)` consists of squares by Hensel.
Class field theory gives 16 ordinary ideal-class characters. Their values
on the listed 16 prime ideals form an invertible matrix, proving those
ideal classes independent. Protected waves keep these anchors free during
elimination and never target them. A purported relation supported only
on the anchors causes an explicit failure.

## Authorized relation continuation

The user set the goal of reducing the quotient to 16. Each wave freezes a
finite target list, candidate region, smoothness bound, one-worker time
and memory limits, source hashes and checkpoint policy before execution.
Targets are free prime-ideal columns after supported row reduction.
Square boxes use exact Hessian-reduced index-prime lattices. Near-root
strips use the three real norm roots transformed into those lattices:
exact algebraic arithmetic floors each slope times `2^96`, and each
positive denominator v tests the integer center and its two neighbors.
Primitive-coordinate filters and previous search regions are replayed.
Protected waves also skip targets that have become elimination pivots;
the checker reconstructs every such adaptive decision. Increasing the
smoothness cutoff permits reexamining a previously searched region.

Every retained witness has its norm and principal ideal checked exactly,
including the nonmonic norm identity's fixed square factor. All prime
coordinates above 37,638 remain until exact elimination cancels them.
Supported rank is independently checked as `rank(all)-rank(outside)`.
Rejected candidates are unnecessary for a rank upper bound; counts and
digests describe the search, without an exhaustive-miss assertion.

Residual-representative waves also allow a pivot prime ideal, or an outside-
base ideal whose outside normal form vanishes, when its reduced
normal form contains unresolved nonanchor coordinates. A greedy spanning
set is selected in ascending prime order, then filled with inexpensive
eligible ideals. A target is skipped only when its normal form lies
entirely in the certified anchor span. These normal forms guide selection;
they are not treated as newly proved class-group relations.

The capped implementation bounds every product-tree node by `M+1`, where M
is the primorial. Capping commutes with positive multiplication, and a
remainder at most M is unchanged modulo any product above M. Thus every
leaf remainder stays exact. The [fixed benchmark](../../artifacts/generated-results/elliptic-curves/small_conductor_capped_remainders_benchmark_v1.json)
checks agreement against both the full tree and scalar division on 58,819
actual values plus 200 edge-case lists. Timing is specific to that target.

## Audited waves

| Wave | Targets checkpointed | Region | Retention cutoff | Candidate occurrences | Relation occurrences | Independent supported gain | Norm-matrix dimension | Worker seconds |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| [box 1](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_wave_001_v1.json) | 512 | box 128 | 37638 | 9891662 | 824 | 588 | 2291 | 130.5 |
| [box 2](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_wave_002_v1.json) | 512 | box 128 | 37638 | 9901917 | 807 | 434 | 1857 | 128.32 |
| [box 3](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_wave_003_v1.json) | 512 | box 128 | 50000 | 9918007 | 1248 | 346 | 1511 | 143.13 |
| [box 4](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_wave_004_v1.json) | 512 | box 128 | 50000 | 9847135 | 1181 | 345 | 1166 | 140.91 |
| [strip 1](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_strip_wave_001_v1.json) | 64 | v <= 1024 | 50000 | 358109 | 55 | 24 | 1142 | 4.09 |
| [strip 2](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_strip_wave_002_v1.json) | 512 | v <= 1024 | 50000 | 2866338 | 431 | 201 | 941 | 34.0 |
| [strip 3](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_strip_wave_003_v1.json) | 485 | v <= 1024 | 50000 | 2715170 | 407 | 98 | 843 | 31.65 |
| [strip 4](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_strip_wave_004_v1.json) | 512 | v <= 4096 | 50000 | 8594390 | 472 | 228 | 615 | 126.82 |
| [strip 5](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_strip_wave_005_v1.json) | 247 | v <= 4096 | 50000 | 4147007 | 211 | 76 | 539 | 60.15 |
| [protected 1](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_protected_wave_001_v1.json) | 64 | v <= 4096 | 100000 | 1430177 | 334 | 4 | 535 | 31.23 |
| [protected 2](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_protected_wave_002_v1.json) | 450 | v <= 4096 | 100000 | 9375401 | 1982 | 48 | 487 | 201.8 |
| [protected 3](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_protected_wave_003_v1.json) | 465 | v <= 16384 | 50000 | 19513829 | 530 | 240 | 247 | 437.05 |
| [capped 1](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_capped_wave_001_v1.json) | 225 | v <= 32768 | 50000 | 12218826 | 227 | 102 | 142 | 111.48 |
| [capped 2](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_capped_wave_002_v1.json) | 123 | v <= 65536 | 50000 | 12798226 | 162 | 65 | 77 | 114.79 |
| [residual 1](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_residual_wave_001_v1.json) | 64 | v <= 1024 | 50000 | 358297 | 63 | 10 | 67 | 3.07 |
| [residual 2](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_residual_wave_002_v1.json) | 512 | v <= 1024 | 50000 | 1916673 | 284 | 49 | 18 | 16.99 |

The capped phase seeds the known-point principal parity relations directly.
They supply **3 independent rows** beyond its norm-only seed; this gain
is already included in the displayed capped-phase dimensions.

From the starting bound 2,879, the formal relation quotient has fallen to **18**.
The explicit-formula completion proves that the sixteen independent anchor classes generate, reaching **16 under GRH**.
The curve-rank upper bound is **22 under GRH**.
These are ideal relations, not new rational points or an algebraic-rank parity claim.

## Replay

```bash
sage -python elliptic-curves/cas/certify_small_conductor_class_completion.sage check
```

The checker replays the inherited curve and field proofs, principal-ideal
witnesses, target selection, matrix transitions and applicable point relations.
It also verifies all prime memberships in the anchor span with a second
elimination order and checks both positive explicit-formula margins.
Raw protocols and checkpoints are under the matching `small-conductor-class-target*`
directories in `artifacts/local/elliptic-curves/`. Each generated certificate
pins its exact sources and witness chunks.
The [portable replay](../../artifacts/generated-results/elliptic-curves/small_conductor_class_target_portable_replay_v1.json) passes from a fresh extracted directory.
