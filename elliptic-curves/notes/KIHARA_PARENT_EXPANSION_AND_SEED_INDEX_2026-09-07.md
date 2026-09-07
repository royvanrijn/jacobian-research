# Four further Kihara parents, corrected seeds, and completed exposure

**Four additional K3 parents, pairwise nonisomorphic over Q and distinct
from the previous seven, now have tested fibres.** This gives at least
**eleven Q-distinct tested parents**. The four fixed Kihara parameters
3/2, 5/2, 11/3 and 1009/101 each received 49 completed point-search boxes.
All four certified lower bounds remain **14**: 196 completed boxes,
**zero discovered directions**. No near-record curve was obtained.

Authority: `EC-KIHARA-FOUR-PARENT-EXPANSION-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json).
The [report](../../artifacts/generated-results/elliptic-curves/kihara_parent_expansion_report_v1.json)
binds the retained transcripts, complete-cloud proofs, model transports,
independent replay and catalogue comparisons.

The [subsequent first-parent proof](KIHARA_FIRST_PARENT_PICARD_AND_RANK_2026-09-07.md)
now resolves the t=3/2 parent: rational/geometric Picard17/18, generic
MW12/13, full rational NS determinant756, and an all-Q-fibration MW15
upper bound. Exact ranks and NS types for the other three remain unknown.

## Which objects are parents?

The [earlier portfolio audit](PARENT_PORTFOLIO_AND_SECTION_LABEL_AUDIT_2026-09-07.md)
found that all 201 high-rank inventory curves, despite twelve family
labels, came from the production K3 X948. Six fixed Mestre parents then
expanded the tested portfolio to seven. Their
[full NS calculation](MESTRE_FULL_NS_GRAMS_2026-09-07.md) identifies one
additional geometric lattice type, of discriminant 468.

Kihara's rank-at-least-14 **one-parameter path is not a K3 surface**.
Replaying the retained symbolic quartic, removing every polynomial
fourth/sixth-power common factor from its short Jacobian, gives

```
deg A = 144, deg B = 216, deg Delta = 426, chi(O) = 36.
finite fibres: 398 I1 + 8 I2 + 3 I4; infinity: I6.
Euler number = 432; geometric reducible-fibre root rank = 22.
```

The finite discriminant is coprime to c4; the displayed factor
multiplicities and the infinity valuation therefore give the complete
semistable configuration. The Hodge bound on geometric MW rank is 336,
which is only an upper bound and offers no discovery claim. The
[global certificate](../../artifacts/generated-results/elliptic-curves/kihara_global_parent_geometry_v1.json)
retains the exact invariant polynomials and removed scale.

The path varies both a six-root K3 parent and its fibre. For fixed p,q,
put a1=0 and

```
a2 = (2p^2+pq+2q^2)^2
a3 = 2(p+q)^2(2p^2+pq+q^2)
a4 = q^2(4p^2-pq+4q^2)
a5 = p(2p-q)(2p^2+4pq+5q^2)
a6 = 4p^4+8p^3q+9p^2q^2-2pq^3+2q^4.
```

Normalize the six roots by a2, and use T=U/a2. The monic degree-six
square approximation G to
`F=product_i (x-a_i/a2-T)(x-a_i/a2+T)` yields the quartic
`y^2=(G^2-F)/T^2`. Its twelve root points and one further explicit
Kihara point supply twelve Jacobian sections when the first is the origin.
For each of the five fixed path parameters, exact specialization sends
these sections to `P2-P1,...,P13-P1` in a certified rank-14 group.
These twelve differences are independent, hence the generic Q-rank is
at least twelve. This argument does **not** assert generic rank fourteen.

All five fixed parents have degrees (8,12,20), squarefree finite
Delta, and split I4 at infinity: **20 I1 + I4**, hence K3 surfaces.
Their exact generic MW ranks, geometric Picard ranks, full NS lattices,
and whether these represent a new geometric lattice type are **UNKNOWN**.
At the published t=2 control, reduced p=5,q=-18 gives T=453/608.
The other four configurations use
`p=t^2(8+3t^2)`, `q=-6(2+t^2)(4+t^2)` at the fixed parameters above.

## Parent separation with a fixed prime pool

The [five-parent certificate](../../artifacts/generated-results/elliptic-curves/kihara_five_parent_distinctness_v1.json)
uses only primes 131,239,251, without extending the pool after a miss.
Each accepted prime passes the exact semistable reduction gate. Smooth
fibre counts are independently recomputed by Sage; singular fibres are
counted directly, with the split I4 resolution contributing 3p.

| Path t | #X(F131) | #X(F239) | #X(F251) |
|---|---:|---:|---:|
| 2, published control | — | — | 67178 |
| 3/2 | — | 61424 | 66948 |
| 5/2 | 19176 | 60803 | — |
| 11/3 | 18888 | — | 66611 |
| 1009/101 | 19164 | 60983 | 67260 |

A dash means that this fixed displayed-model gate did not pass; it is not
a proof that no good model exists at that prime. Different counts at a
common good prime prove Q-nonisomorphism of the smooth proper surfaces.
Using the previous seven parents' replayed certificates gives **44 of 45**
new comparisons separated. The sole unresolved pair is control t=2 versus
t=5/2, which has no common accepted prime in the fixed pool. All four
fresh pilot parents are mutually separated and separated from all seven
previous parents. Eleven is thus a proved lower bound; twelve is not yet
proved by these witnesses. Q-nonisomorphism is not literature novelty.

## A seed subgroup gap and excessive denominator clearing

The original fourteen displayed Kihara points use the fifteenth quartic
point as origin. If K is its image under the hyperelliptic involution,
the twelve product-root points satisfy the exact relation

```
6 K = P1 + ... + P12.
```

The divisor of G-y is the sum of the twelve product-root points minus
six times each infinity point. Also the two infinity points are linearly
equivalent to the origin plus its involution. This proves the relation;
the replay checks it directly in each specialized group as well.
Replacing P1 by K enlarges the supplied subgroup by **exact index six**:
`P1=6K-P2-...-P12`, and the basis change has determinant six.
The new fourteen points are independently certified, so this determinant
is an actual subgroup index. It is **not a rank gain**.

The initial intake's finite mod-two rank was thirteen on all four rows.
That result remains preserved as `UNRESOLVED_SEED_QUOTIENT`; it never
supplied a rational rank upper bound. The corrected seeds have finite
quotient rank fourteen modulo 2,3,5. This does not claim full saturation
at every prime or an exact rational rank.

The original intake also overcleared denominators. Using the natural
weights of the degree-(144,216) model and denominator(t)^36, with exact
coordinate transports, gives:

| t | First intake coefficient bits | Corrected bits |
|---|---:|---:|
| 3/2 | 1587 | 585 |
| 5/2 | 1673 | 671 |
| 11/3 | 2054 | 885 |
| 1009/101 | 9455 | 2264 |

The frozen V1 sources, seeds and failed quotient diagnostics remain intact.
V2 is a separate intake. No point search ran on the excessive V1 models.
The corrected adapter is
[`intake_kihara_fresh_fibres_v2.sage`](../cas/intake_kihara_fresh_fibres_v2.sage).
The published Kihara source and baseline certificate are preserved.

## Completed four-fibre experiment

The protocol froze these four retained fibres, their corrected own14
seeds, 2048 distinct nonzero SHA256 parity masks each, and 49 charts
selected by largest computed CVP norm. Canonical heights used 384-bit
precision, rounded at 10^6; LLL was unimodular and parity/norm checks exact.
Factor-free Gauss/hyperellred maps avoided quartic-content factorization.
All 196 maps completed before the first point search.

Each chart received H=125000 and a ten-second PARI cap. Each row had
180 seconds for preparation, 600 for the point worker and 300 for history
replay, at most 2 GiB and one worker. There was no scoring, validation-prime
selection, public-target selection, rank stopping, refill or next wave.

All 196 searches report `bounded_search_complete`. Each full retained
cloud contains exactly the fourteen seed points. History replays, exact
geometry, and complete-cloud modulo 2,3,5 audits pass. An independent Sage
implementation enumerates complete finite elliptic groups and their
quotients, checks basis independence, and excludes rational two-torsion.
The [rank replay](../../artifacts/generated-results/elliptic-curves/kihara_fresh_point_pilot_sage_replay_v1.json)
is separate from the producer. Completed exposure trusts the pinned PARI
transcripts; these bounded misses do not exclude higher rational ranks.

The total recorded supervised time is **222.899432532 seconds**, including
both seed intakes, preparation, workers, history, finite audits, parent
geometry and standalone replay. The enclosing cloud-audit wrapper is
reported separately and not added on top of its constituent jobs.

Post-search exact j/isomorphism comparisons find no matches for any of
the four curves against the pinned 620-equation catalogue or 201-equation
inventory. These are four mutually different j-invariants. This limited
comparison does not justify “never seen before”; the rank-14 curves are
not added to the high-rank inventory. No larger sweep follows the pilot.

## Reproducibility and remaining gap

The [portable bundle](../../artifacts/generated-results/elliptic-curves/kihara_parent_replay_bundle_v1.json)
contains the retained global quartic, five parent equations, thirteen
quartic points and twelve section coordinates per parent, old/new models,
rank witnesses, and prior-parent count certificates. Copy it with
[`verify_kihara_parent_bundle.sage`](../cas/verify_kihara_parent_bundle.sage)
to an empty directory and run:

```sh
sage -python verify_kihara_parent_bundle.sage --input kihara_parent_replay_bundle_v1.json
```

The fresh-directory replay passed in **1.406585056 seconds**, under a
120-second/2-GiB cap. It imports no repository code. It reconstructs the
fixed root quartics and point maps, checks global minimal geometry from
the retained symbolic quartic, independently proves control independence,
checks both seed-model transports and index-six relations, and recomputes
ten accepted new surface counts. Earlier seven-parent counts remain
explicit dependencies of the separation conclusion.

```sh
python3 elliptic-curves/cas/report_kihara_parent_expansion.py --check
sage -python elliptic-curves/cas/verify_kihara_fresh_point_rank.sage --check
```

The construction is grounded in [Kihara (2001), pp.50–51](https://doi.org/10.3792/pjaa.77.50)
and the retained exact source replay; the surface geometry uses the usual
minimal Weierstrass and Shioda–Tate framework
([Schütt–Shioda](https://arxiv.org/abs/0907.0298)).

The search now reaches more parents. It still lacks evidence that these
new parents deliver near-record specialization gains, and lacks a useful
high-generic-rank alternative fibration on these Kihara parents. Resolving
those geometric or point-visibility gaps is a distinct experiment, not
permission to repeat the completed boxes or launch a larger population scan.
