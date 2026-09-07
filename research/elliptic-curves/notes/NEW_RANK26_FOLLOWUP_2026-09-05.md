# Follow-up and remaining gaps for the new rank-at-least-26 curve

Update: the [completed gap and exact-geometry audit](EXACT_PARITY_AND_COORDINATE_AUDIT_2026-09-06.md) supersedes the unlaunched-tail and omitted-class status below. All 301 tails and 48 extra chart attempts now replay; the complete retained cloud has 2,338 points and still certifies rank at least 26. The original experiment and its frozen evidence remain unchanged.

The [new curve](NEW_COMPACT_RANK26_CURVE_2026-09-05.md) has 26 independently
certified points. Its first adaptive follow-up completed all **301 planned
chart attempts** in 1,742.69 seconds under its 1,800-second / 1.5-GiB limit.
The worker found no further independently certified direction. Independent
replay passed for all 301 chart/admission histories and the complete
1,947-point cloud through prime 997, still with lower bound 26.

These charts use only the curve's 26 discovered points: its original
17 generic directions and nine further independent directions. The
predeclared class pool combines 301 generic parity representatives with
cyclic nonzero nine-bit quotient words, then orders the resulting centres
by numerical specialized norms. PARI canonical heights use explicit
384-bit precision. The numerical geometry schedules searches; it is not
an independence certificate or a proof of optimal CVP representatives.

Each chart uses height 100,000 and four seconds of search time, with
admission primes through 997. All 301 individual height boxes timed out;
mean denominator-range coverage is **0.8080451827242525**. Thus completing
the declared attempt is not exhausting these height boxes, and no upper
rank bound follows. Full raw point clouds are retained, including points
that did not enlarge the independently certified basis.

The [replayer](../cas/replay_compact_r17_new26_followup.py) verifies the
generic transport, integral parity representatives, rounded-metric norms,
chart ordering, every retained chart and point admission, and the final
state. The [complete-cloud certificate](../../artifacts/generated-results/elliptic-curves/compact_r17_new26_followup_recorded_mod2_v1.json)
is a separate exact check through prime 997. Verification logs and limits
are under `artifacts/local/elliptic-curves/compact-r17-new26-followup-v1/verification/`.

## Exact family-incidence extension

The original incidence proof covered the first 32 curves only. The seven
subsequent additions, stable IDs 33–39, now have their own
[exact certificate](../../artifacts/generated-results/elliptic-curves/latest7_cross_family_j_incidence_v1.json)
and [Sage-free replay](../../artifacts/generated-results/elliptic-curves/latest7_cross_family_j_incidence_replay_v1.json).
All **84 pairs** with the same twelve presentations are resolved:

- 77 have no rational projective `j`-preimage, by good-reduction image exclusions;
- seven have complete rational-preimage factorizations, with just the original
  parameter in the original family and no preimage at infinity.

Consequently none of these seven curves obtains another generic subgroup
by transferring sections from those other presentations. Together with the
original proof this covers all **39 × 12 = 468** pairs, without extending
the conclusion to future fibres or untested families. Equal `j` alone is
not a rational-isomorphism or point-independence certificate; the unique
remaining own-family fibres already have exact transports in their original
point certificates.

## Exact orbit compression for future workers

For a pointed quartic with pole divisor `O+C`, the two points over one
parameter satisfy `P+Q=C`. When `C` is an explicitly verified integral word
in the current subgroup, adjoining either point also generates the other.
The optional [compression helper](../cas/research_runtime/pointed_orbit_compression.py)
keeps the first point of each involution orbit `P -> C-P` and records the
exact partner relation for every omitted point. Raw point clouds remain
intact. Fixed points and unpaired points are retained.

On the already replayed 43 initial charts of the new curve, the
[exact algebraic audit](../../artifacts/generated-results/elliptic-curves/pointed_orbit_compression_audit_v1.json)
reduces **384 admission events to 192**, with 192 checked relations and all
43 centres recomputed from integral subgroup words. Targeted tests cover a
new independent pair, a fixed point, invalid centre words and off-curve
inputs. This proves preservation of the generated integral subgroup; it
does not claim identical observation histories or a measured wall-time
speedup. No frozen or running worker was changed to use this helper.

## Conductor and reproduction limits

A 60-second bounded local audit trial-divides the minimal discriminant by
primes through 10,000. All bad primes except 3 have `c4` a unit and hence
multiplicative reduction; exact Tate local data at 3 give its conductor
exponent. The remaining cofactor has 100 digits and is not factored or
asserted prime. The exact conductor divides

```text
664058614956913086062441294281223548751518341016391732166558605492361906820548369081057425270473176425309691030.
```

This upper bound does not beat the listed rank-at-least-26 conductor
minimum. It neither establishes a conductor record nor excludes a smaller
actual conductor. No full factorization or descent campaign was launched.
The separate rank-at-least-22 small-conductor result remains the current
proved low-conductor near-record comparison.

The initial experiment's 4,897-file evidence archive passed member-hash
checks, and its isolated rank26, 24-result, 43-chart and six-census replays
passed. The first isolated 39-curve inventory replay reached its
300-second cap after printing 24 verified rows. A separately logged
600-second retry completed the entire exact inventory check in 110.62
seconds. Both attempts are retained; local full-inventory replay also passed.

## Concrete remaining search gaps

All 301 adaptive boxes have explicitly retained unsearched denominator tails.
A bounded continuation could cover those intervals using the same curve and
chart maps. Their endpoints and original state keys are recorded in
`compact-r17-new26-followup-v1/unsearched-tails.json`; no tail continuation
has been launched. New limits and exact chart bindings are required.

Families `11952` and `08f72` each have 49 representatives of **computed**
generic norm 12 in the fresh census. The balanced initial experiment's
fixed 43-chart cap omits six of those classes per family. Its eight retained
candidates therefore admit a concrete possible 48-chart extension. This
gap is recorded, but that extension has **not been launched**. It needs its
own frozen limits and replay plan; computed norms do not prove exact coset
minima, new rational points or rank jumps.

The [follow-up evidence manifest](../../artifacts/generated-results/elliptic-curves/new_rank26_followup_evidence_v1.json)
and adjacent ZIP retain the finite experiment and these diagnostics.
Exact rank and new rank-at-least-28/32 curves remain open.
