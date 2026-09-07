# Restricted-prime point exposure and an invariant coefficient floor

The new restricted-prime mapping passes the full blind27→28 regression.
Exact history comparison then avoids98 already completed boxes on two retained27
curves. A separate98-box refinement on the two determinant1092 fibres completes
with no returned finite points and no added direction.

The further obstruction is partly intrinsic: every integral completed-square
binary quartic with the two fibres' j-invariants needs a coefficient of at least
255 or438 bits. Their best refined charts have258 or444 bits. Coordinate
minimisation cannot put these fibres into the83–92-bit quartic range of the two
retained27 curves. This is a coefficient bound, not a rank or point-height bound.

## Mapping policy and blind gate

The original [factor-free control and point portfolio](BLIND_FACTOR_FREE_CONTROL_AND_PROSPECTIVE_EXPOSURE_2026-09-07.md)
remain unchanged. The new mapping first constructs their Gauss/hyperellred
quartic, then invokes `hyperellminimalmodel` with the explicit list of25 primes
from2 through97, followed by `hyperellred`. PARI documents that the optional
prime list restricts the minimality guarantee to that list; the installed
PARI2.15.4 supports the option. No unrestricted minimisation or complete
factorization of a large discriminant is requested. [PARI reference](https://pari.math.u-bordeaux.fr/dochtml/ref/Hyperelliptic_curves.html)

Every returned transformation is checked by an exact polynomial identity and
rational square scale. Point correctness does not depend on asserting the
minimality of the output. This is a new finite coordinate policy, not a claim
that the old and new search boxes are equivalent.

`blind_bounded_prime_28_control_v2.py` regenerates all2048 masks and49 centres
from only the original27-point seed of the known public28 curve. The artifact
read whitelist again excludes the public exceptional point, previous maps and
winning-chart data. All49 maps precede points, and all49 height125000,
ten-second boxes complete without rank stopping. The full49-point cloud
certifies28 modulo2,3,5; standalone Sage and exact geometry/history checks pass.
All49 refined boxes have an explicit new-coordinate witness relative to the
factor-free boxes. The known-curve recovery is not a new rank28 discovery.

The preserved first version stopped during preparation because serialized
matrix entries had not been converted back to rationals. It searched zero
boxes. Its2.166648820 seconds are included in the control's67.711791708
supervised seconds. The [control certificate](../../artifacts/generated-results/elliptic-curves/blind_bounded_prime_28_control_v2.json)
retains both versions and their input-access evidence.

## Prospective history gate and results

The frozen four-curve protocol retains every prior seed, parity mask and
centre. It selects ID186 from the R17 branch, ID90 from MW16, and the already
fixed determinant1092 fibres T=1 and1009/101. No public exceptional point,
new parameter population, score or validation statistic enters selection.
All196 refined maps freeze before any point search.

For each map, compare the new box with the completed factor-free box and any
explicitly declared same-centre global-minimal box. ID186's latter box completed
at height1000000; ID90's completed at125000. A curve with no proved new-history
coordinate witness receives zero point boxes, with no refill. Otherwise its
full49-chart roster runs under the unchanged125000/ten-second limits.

| Curve | Refined maps | Proved covered by prior boxes | New boxes searched | Initial → final bound |
|---|---:|---:|---:|---:|
| ID186, R17 | 49 | 49 | 0 | 27 → 27 |
| ID90, MW16 | 49 | 49 | 0 | 27 → 27 |
| determinant1092, T=1 | 49 | 0 | 49 | 17 → 17 |
| determinant1092, T=1009/101 | 49 | 0 | 49 | 17 → 17 |

The independent replay strengthens the skip decision: every retained27
transition to the old completed global-minimal chart is diagonal with entries
1 and ±1. Thus all98 skipped boxes are proved contained in their declared
completed history, rather than merely lacking a discovered exposure witness.
The comparison does not claim to cover every different-centre historical chart
or every possible translated representative.

Both new-parent rosters have an explicit coordinate witness outside their old
box on every chart. Such witnesses need not lift to rational curve points.
All98 searches complete but return zero finite points. Both full point clouds
therefore consist of the17 input sections, verified independently modulo2,3,5
and by standalone Sage. No rank upper bound follows.

The [completed portfolio report](../../artifacts/generated-results/elliptic-curves/bounded_prime_point_portfolio_v1.json)
accounts for136.746749954 supervised seconds, including the preserved failed
independent-proof attempt and its correction. Together with the blind control,
the cost is204.458541662 seconds. Nested driver time is not counted twice.

## Universal coefficient bound

Let j=N/D be reduced, with D>0. Write bN and bD for the bit lengths of |N| and D.

For an integral short model y²=x³+Ax+B whose largest coefficient has k bits,
|A|,|B|<2^k. Its invariant is

```
j = 6912 A^3 / (4 A^3 + 27 B^2).
```

Before cancellation, numerator and denominator have at most3k+13 and3k+5
bits respectively:6912<2^13 and4+27<2^5. Cancellation can only decrease them.
Hence every such model has

```
k >= max(1, ceil((bN-13)/3), ceil((bD-5)/3)).
```

For an integral binary quartic with coefficients a,b,c,d,e of absolute value
less than2^k, use the classical invariants

```
I = 12ae - 3bd + c^2
J = 72ace + 9bcd - 27ad^2 - 27b^2e - 2c^3
j = 6912 I^3 / (4I^3 - J^2).
```

These conventions agree with [Cremona–Fisher, *On the equivalence of binary quartics*](https://www.dpmms.cam.ac.uk/~taf1000/papers/quarticequiv.pdf).
The following size estimate is the elementary deduction used here:
|I|<16·2^(2k), |J|<137·2^(3k),6912·16³<2^25 and4·16³+137²<2^16.
Thus, after any cancellation,

```
k >= max(1, ceil((bN-25)/6), ceil((bD-16)/6)).
```

This applies to integral completed-square models Y²=F(X,Z). It is not a bound
on the separate coefficients of a generalized model y²+Q(x)y=P(x), whose
completed square has F=4P+Q². It also places no lower bound on the height of an
individual rational search coordinate.

| Fibre | j numerator/denominator bits | Integral short coefficient floor | Integral quartic coefficient floor | Observed refined quartic bits |
|---|---:|---:|---:|---:|
| determinant1092, T=1 | 1555 / 1533 | 514 | 255 | 258–350 |
| determinant1092, T=1009/101 | 2648 / 2629 | 879 | 438 | 444–615 |

[Independent Sage verification](../cas/verify_bounded_prime_history_and_height_v2.sage)
recomputes all196 rational centres, polynomial map identities, quartic
invariants, j-values, explicit new-coordinate witnesses and prior-box
inclusions. The [certificate](../../artifacts/generated-results/elliptic-curves/bounded_prime_history_and_height_v1.json)
also gives floors79/80 and observed ranges83–92/85–90 for IDs186/90.
The first proof attempt used matrix iteration where a flat list was required;
that failed API conversion and its0.692040725 seconds remain recorded.

## Exact visibility audit

All6,664 tests of the17 original signed section representatives in the two
old and two refined rosters lie outside height125000. There is no
within-completed-box omission among these witnesses. The smallest heights are:

| Fibre | Before refinement | After refinement |
|---|---:|---:|
| T=1 | 49542743827981 | 49542697072701 |
| T=1009/101 | 6524843559948928288229712911 | 13049687119897856576459425822 |

The [exact visibility certificate](../../artifacts/generated-results/elliptic-curves/det1092_refined_visibility_v1.json)
retains each section's minimum and the ordered observation digests. This is a
retrospective check of these representatives, not a masked calibration or an
exclusion of other representatives in their cosets. The earlier known28
control already showed why representative-only visibility cannot establish
exceptional-direction sensitivity.

The universal coefficient floor explains why further normalization cannot
make these two fibres as small as the retained27 comparison curves. It does
not exclude their having further independent points. The completed results
support an exposure/calibration or arithmetic-height change, not an automatic
increase in parameter-population size.

## Recheck

```sh
python3 elliptic-curves/cas/report_blind_bounded_prime_28.py --check
python3 elliptic-curves/cas/report_bounded_prime_exposure.py --check
python3 elliptic-curves/cas/audit_det1092_refined_visibility.py check
```

The fresh proof producers preserve existing outputs. To regenerate them, use
an explicitly new run directory/version; do not overwrite the frozen records.
