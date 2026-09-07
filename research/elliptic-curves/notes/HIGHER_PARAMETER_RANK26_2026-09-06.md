# Higher-parameter pilot: a new rank-at-least26 curve

The fixed higher-parameter R17 pilot supplies two new catalogue-unmatched
inventory entries. The stronger one, `new-20260906-99`, is the `11952` fibre at
`7460/32309`. Its initial49 charts certify25; a separate301-centre follow-up
raises the exact lower bound to26. Its globally minimal equation is

```
y^2 + x*y = x^3
 - 67032315925132515518574925613377641806054416249435904230*x
 + 211585624000252210792909737627250671997115932131400571389543578375708225639461940452.
```

The [minimal-model proof](../../artifacts/generated-results/elliptic-curves/higher26_minimal_proof_v1.json)
contains26 exactly independent rational points. The invariant gcd is1, which
closes global minimality. The [Sage point file](../../artifacts/generated-results/elliptic-curves/new_higher_rank26_curve_11952.sage)
exports the equation and witnesses. This is a lower bound, not an exact rank.

The other addition, `new-20260906-100`, is `103b2` at `4271/6508`, with22
certified independent points. Both are absent up to rational isomorphism from
the pinned593 catalogue equations and402 earlier measured equations. All24
new attempts are mutually nonisomorphic and unmatched in those comparisons;
only these two meet the inventory threshold22. No universal novelty is proved.

The [V11 inventory](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v11.json)
now contains100 distinct curves: six lower bounds27, eleven26, twenty-one25,
twenty-four24, twenty23 and eighteen22. Existing IDs are preserved.
All100 proofs and theCSV replay. The two additions have only their own rational
parameters in the twelve recorded presentations. The aggregate now accounts
for1200 incidence pairs:1079 exclusions and121 rational preimages, including
the same21 duplicate presentations of the existing R17 subgroup.

## Population and point exposure

The [previous parameter-height audit](NEW_MW16_RANK27_2026-09-06.md) placed all
eight mapped public reported-rank-at-least28 examples outside parameter
height16384. The separately verified periodic Nagao accumulator made a higher
population affordable without changing the original integer score policy.

For each of six families and two signs, the frozen pilot chooses one denominator
residue modulo64 by SHA256 of a literal salt and the frame/sign. Preassigned
parities give six odd and six even slices. Numerators and denominators are
bounded by32768. This is122368792 primitive addresses, close to the earlier
122400468 full4096 addresses, but it is neither uniform sampling nor full32768
coverage. Public parameters, equations, ranks and points do not choose slices.

Each slice retains512 candidates using the562 cached primes5..4093, giving6144
addresses. All receive the same5978 additional traces4099..65521. The selector
uses the combined sum through32749, then good count, denominator and signed
numerator; validation32771..65521 never enters ties. All short and extended
scores and the fixed top-four-per-family roster replay. The extended phase
includes36728832 new prime traces; no missing trace is treated as a good score.

The24 point attempts start from only17 generic sections and use every43/49
exact generic maximum parity class. All1080 boxes complete at height125000
and ten seconds per chart, with at most two workers. The initial bounds are
13x17,3x18,2x19,2x20,2x21,1x22,1x25. All histories and full retained clouds replay.
The parameter scan took21.774 seconds and its replay2.776; the prime extension
563.697 and replay31.088; the point batch543.273 and its verification67.768.
These are recorded local timings, not portable performance guarantees.

The unique fresh25-point result receives a separate301-centre follow-up.
Its eight discovered directions beyond generic17 had not entered the initial
centres. The protocol pairs301 largest computed generic census norms with
nonzero8-bit quotient words, ordered by Hamming weight then integer, cycling
after255 words. These are301 pairings, not an exhaustive quotient search.
Numerical heights choose representatives; rational centres, parity, maps,
raw points and independence are checked exactly. Chart81 adds the26th
point, with generic mask1044 and quotient word213. All301 boxes complete at
unchanged125000/10-second limits. The combined677-point cloud certifies26
modulo2,3,5. Further centre experiments have separate protocols and evidence.

## What the larger population exposed

The [exact descriptive audit](../../artifacts/generated-results/elliptic-curves/higher24_visibility_cost_v1.json)
compares the initial24 outcomes with the earlier24 wider-retention attempts,
using identical per-chart budgets. The median largest minimal coefficient
increases from227 to288 bits; the median retained cloud decreases from573 to
49.5 points. All48 normalized models pass the existing cheap global-minimality
proof. The eight mapped public reported>=28 examples have a median277.5-bit
largest minimal coefficient. Their ranks are catalogue metadata in this audit.

The new25 curve itself has277-bit largest minimal coefficient. Its adaptive
gain shows that the initial generic-centre exposure was incomplete for its
rational point subgroup. The comparison does not establish why the other
curves yielded fewer points: coefficient size, score quality, true rank and
point visibility are confounded. No bounded result excludes a higher rank.

There is also a selector question. Elkies–Klagsbrun use
`sum_good log(#E(F_p)/p)` and discuss the deterioration of score quality as the
search region grows; they do not give a universally calibrated height/prime
cutoff. Their experimental split-multiplicative bonuses are also heuristic.
[Primary paper, sections2 and7](https://arxiv.org/pdf/2003.00077).
Watkins et al. compare weighted and unweighted scores and warn that extending
a prime range alone need not improve discrimination.
[Primary paper, section5](https://www.dpmms.cam.ac.uk/~taf1000/papers/rankcongr.pdf).

Our earlier `compare_bounded_prime_selectors.py` already compared the two sums
on a small ordinary-fibre diagnostic panel. It did not promote a policy or
calibrate the present higher-height population. A retained-trace comparison
here would be a new bounded experiment, still unable to recover addresses
removed by the short-prime retention step. Neither alternative score nor a
coefficient-size cut is an established remedy for the current misses.

New rank-at-least28/32 curves, exact ranks and universal novelty remain open.
