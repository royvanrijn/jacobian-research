# Million-height control and the saved MW16 candidate population

`../MATH_STATUS.json` remains authoritative. This note records bounded visibility,
backend and selection experiments following the
[wider-retention discoveries](WIDER_RETENTION_DISCOVERIES_2026-09-06.md).
Independent rational points certify lower bounds; completed boxes never certify
an upper bound or the absence of more points.

## The native known29 control now recovers its last direction

The prior exact translated-point audit placed a representative of the remaining
known direction at height918522 in original generic chart12. A separately frozen
benchmark searched that chart at height1000000 for60 seconds, with two unchanged
engines and the same independently certified28-point seed. PARI completed in
38.012909 seconds and recovered29; GMP was censored after60.092315 seconds,
with117838 completed denominators, retaining28. Both point/admission histories
and exact independent-point proofs replay.

This is **retrospective calibration**, not a blind validation or a new curve.
The chart and budget were chosen using the earlier public-point visibility
result. Neither worker consumed the public point or its translation word.
It supplies a finite cost/visibility gate for a separately declared pilot;
it does not guarantee that arbitrary new curves have a missing point in this box.

A second exact diagnostic tested98 earlier CVP proposals plus/minus every row
of the existing28-row LLL basis:5488 translations, without another numerical
metric calculation or adaptive enlargement. All lie outside the completed
125000 boxes. Their best height is21823475806, worse than the earlier918522.
All5488 full rational group words and coordinates replay. This rules out an
improvement among these proposals only, not a better translate in general.

Certificates:
[`native29_million_chart_benchmark_v1.json`](../../artifacts/generated-results/elliptic-curves/native29_million_chart_benchmark_v1.json),
[`native29_cvp_neighbours_replay_v1.json`](../../artifacts/generated-results/elliptic-curves/native29_cvp_neighbours_replay_v1.json).

## A correct sieve prototype is not a general speed improvement

The separate `pointed_quartic_sieve_v2.cpp` adds a necessary primitive-residue
filter, skips empty denominator rows and replaces repeated remainder operations
in the block loop. Its exact primitive square hits agree with Python brute force
on48 fixed quartics and two intervals each (96 checks), and with the unchanged
worker on three full real100000-height boxes.

Measured old/new times are5.253/6.308,5.223/6.090 and5.057/3.622 seconds.
Two of three charts regress. The prototype is **not promoted**; the production
worker and its frozen certificates remain unchanged. These three timing pairs
are not a general performance theorem.

Certificate:
[`pointed_sieve_v2_benchmark_v1.json`](../../artifacts/generated-results/elliptic-curves/pointed_sieve_v2_benchmark_v1.json).

## Complete the older MW16 exposure before extending its selector

Both previously new MW16 curves with25 certified points still had incomplete
old metric-coordinate boxes. A fixed follow-up kept their exact43 generic16
centres, original order, own25-point seed and height100000, but used PARI
horizontal maps/backend and a10-second cap. All86 boxes and both admission
histories completed. The union of old and new retained clouds contains1241 and
650 points respectively; exact mod2,3,5 checks still give25 on each. This changes
coordinate/backend/time cap together and makes no single-factor causal claim.

The five-family MW16 selector had also retained only562-prime scores. A fixed
five-address benchmark completed5978 additional prime traces per address, with
40 direct Python character sums agreeing. All1280 saved signed H4096 finalists
were then extended, reusing the existing102000390-address population without a
new scan. All562 original short scores per candidate and every extension roster
replay. The extension ran in199.060 seconds; its replay took6.912 seconds.

Four candidates per family were selected by the combined score through32749,
then selection-band good count, denominator and signed numerator. Validation
primes32771..65521, catalogue records, public points and earlier ranks do not
enter selection, including tie-breaking. All20 receive the same generic16-only
point-search policy, with no replacement of previously tested candidates.

The43 generic parity labels per MW16 family are the previously recorded computed
choices, **not an exact certificate of every generic or specialized maximum**.
This distinction matters: the separate exact MW17 maximum-class theorem does
not automatically apply to these lattices.

All860 point boxes completed; all20 admission histories, exact maps and complete
point-cloud proofs passed independent replay. Post-batch comparison identifies
known curves542 (26 recovered points) and548 (24), which are rediscoveries.
Twelve candidates have no prior-equation or catalogue match: three26s, five24s,
two20s and two18s. The eight fresh curves at or above22 enter the inventory,
alongside the previously measured MW16 curve2407/532 whose bound improves from20
to23. The other previously measured examples receive no inventory improvement.

| Inventory ID | MW16 frame | Parameter | Certified lower bound |
|---|---|---|---:|
| new-20260906-90 | a1-fibration-01 | -1867/270 | 26 |
| new-20260906-91 | a1-fibration-01 | -557/3572 | 26 |
| new-20260906-92 | a1-fibration-02 | 3161/432 | 26 |

All three26s have proved globally minimal integral equations, exact point
transports and standalone Sage point files. Their invariant gcds are169,169,9;
small exact local checks establish minimality without a discriminant-factor
campaign. The [minimal proof](../../artifacts/generated-results/elliptic-curves/extended20_mw16_high_rank_minimal_v1.json)
is canonical for the equations and26 independent points.

A separate post-batch audit transports and combines retained point witnesses
from every matching earlier equation, using exact rational isomorphisms. It
uses no catalogue points. No union increases the initial cohort lower bounds.
This closes a possible loss from comparing only the separate subgroup sizes.

The [V8 inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v8.json)
and its CSV retain all89 earlier IDs and contain98 distinct j-invariants:
five27s, eleven26s, twenty-one25s, twenty-four24s, twenty23s and seventeen22s.
All98 independent-point certificates and catalogue exclusions replay. Absence
from the pinned593 equations is not universal novelty.

The nine new inventory entries have108 exactly replayed incidences across the
same twelve presentations:99 exclusions and only their nine original parameter
preimages. Combined with the prior cohorts, all1176 pairs for98 curves are
accounted for; the same21 R17 duplicate presentations remain the same proved
generic subgroup. No extra high-rank presentation was found in this atlas.

## Exact MW16 parity geometry closes the maximum-class gap

The independent generic-lattice audit checks all327680 parity upper witnesses
and all42 proposed maximum classes. Scaling each Gram by2 gives an integral
lattice; a checked unimodular basis change and exhaustive rational LDL
ellipsoid enumeration exclude every vector of scaled norm<=22 in each queried
class. Their displayed norm23 witnesses therefore prove exact minimum23/2.
All other parities have directly checked witnesses of norm<=23/2.

The exact discrete maxima are23/2 in all five frames, with12,4,8,10,8 maximum
classes respectively. None was omitted by the old43-chart cap. This resolves
that specific generic maximum-class concern; it does not certify specialized
maxima or rule out useful lower-norm classes. It is not the continuous covering
radius or a saturation theorem.

The first producer failed before enumeration because Sage ZZ would not directly
coerce Python Fraction. Version2 checks integral denominators and explicitly
converts them to integers. Both sources, protocols and the failure are retained;
all five version2 queries and independent Python replays pass. The first new
incidence replayer also retained an old cohort cutoff; its version2 corrects89
versus70 and passes against unchanged exact incidence output. An outer reporting
driver used the wrong elapsed-time key after successful arithmetic calls; the
separately launched missing stages all passed. No failed attempt is silently
rewritten or counted as evidence.

## The deeper new27 pilot completed without a28th direction

Exactly one prospective curve was selected by the smallest largest global-minimal
coefficient among the five new27s:11952 at2012/211, inventory72. The same49
original generic maps, original order, own27-point subgroup, height1000000 and
60 seconds per chart were frozen before execution. No known-control chart
choice or oracle entered this worker. All49 boxes completed and its history
replayed, within1820.181 seconds total for the supervised pilot.

The union of initial,301 adaptive125000 and million-height retained clouds
contains2590 points. Exact mod2,3,5 finite-quotient checks still give27. The
[result](../../artifacts/generated-results/elliptic-curves/new27_million_retention_11952_coverage_v1.json)
is a bounded visibility outcome, not rank27 exactly or an exclusion of rank28.
The broader goal of a new near-record/record curve remains open.
