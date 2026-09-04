# Blind multi-stratum half-lattice jump ladder

<!-- status-consumer: EC-K3-R17-REFRESH-BLIND-JUMP-LADDER -->

Status: **passing fixed-panel detector experiment; exact blind rank gains;
exact post-freeze displayed-jump analysis; no rank-32 promotion without the
separate residual-Selmer gate**.

## Outcome

The 573-curve refresh supplied sixteen quotient-eligible new norm-twelve
fibres spanning displayed jumps `+3,+4,+5,+6,+8,+10,+11,+12`.  Curve 499 was
excluded before searching: its exact specialization audit proves that the
displayed subgroup does not contain specialized MW17, so the requested
displayed-subgroup quotient is not defined.

The blinded input contains, for each eligible curve, only its short equation,
the seventeen exactly transported generic points, their exact generic height
Gram, curve id, and atlas class.  It contains no displayed rank, jump, public
point beyond the reconstructed MW17, or public-complement coordinate.

The frozen response was only

> exact quotient rank recovered before opening the public complement.

The sealed responses are:

| curve | exact blind rank gain | displayed jump opened later | initial gain | adaptive gain |
| ---: | ---: | ---: | ---: | ---: |
| 478 | 6 | 4 | 6 | 0 |
| 498 | 6 | 6 | 6 | 0 |
| 531 | 11 | 11 | 9 | 2 |
| 532 | 3 | 3 | 3 | 0 |
| 534 | 11 | 11 | 4 | 7 |
| 535 | 10 | 11 | 3 | 7 |
| 536 | 11 | 11 | 4 | 7 |
| 537 | 10 | 10 | 7 | 3 |
| 538 | 5 | 5 | 5 | 0 |
| 539 | 6 | 6 | 6 | 0 |
| 540 | 8 | 8 | 8 | 0 |
| 541 | 8 | 8 | 8 | 0 |
| 543 | 12 | 12 | 3 | 9 |
| 544 | 0 | 11 | 0 | 0 |
| 545 | 11 | 11 | 6 | 5 |
| 546 | 8 | 8 | 8 | 0 |

Thus the blind search recovers `126/136` displayed-jump dimensions in
aggregate.  This fraction is descriptive: curve 478 recovers six independent
directions beyond MW17 even though its displayed subgroup supplied only four.

## Frozen policy and cross-class correction

The policy gives every fibre the same maximum opportunity:

1. completely enumerate `MW17/2MW17` using the exact generic height form and
   select the first 43 parity classes by decreasing exact generic norm, then
   mask;
2. choose shortest representatives and order them by the specialized
   canonical-height form at scale `10^6`, audited independently at `10^5`;
3. classify every returned point by exact group law and finite-reduction
   independence;
4. when the initial quotient rank is nonzero, recompute the complete current
   lattice state and search a fixed 301-chart adaptive pool; and
5. use height bound 100,000, a 15-second wall timeout, a 1 GB PARI stack, and
   no retry for every chart.

Every initial and adaptive order is bound to a hash of the current basis,
height Gram, generic coordinates, quotient complement, chart universe, and
ordered chart identities.  Fifteen fibres consume all `43+301=344` charts.
Curve 544 has initial gain zero, so the predeclared structural rule has no
quotient coordinate in which to construct an adaptive lift and stops after
43 charts.  There are no timeouts or PARI failures.

Protocol v1 mistakenly required every determinant-948 native MW17 lattice to
have exactly 43 *maximum*-norm parity classes.  Curve 478 completed with score
6, but curve 498 contradicted that lattice-transfer assertion before any of
its charts was searched.  Protocol v2 therefore changes only the class-set
wording to the first 43 in exact generic-depth order.  This selects the same
43 classes in the same order on curve 478, preserves all budgets and
confirmatory endpoints, discloses the already known curve-478 score, and was
frozen before the other fifteen recovery outcomes.  Curve 478 was rerun from
the redacted input and reproduced score 6.

## Confirmatory tests

Before the v2 recovery outcomes, the retained tests were:

- positive Kendall tau-b between exact blind gain and displayed jump, passing
  only if `tau_b >= 0.35` and the exact one-sided tied-margin permutation
  `p <= 0.05`; and
- enrichment of the displayed `+10/+11/+12` tail by detector score at least
  10, passing only if the tail-minus-body risk difference is at least `0.25`
  and one-sided Fisher exact `p <= 0.05`.

Both pass.  Kendall `tau_b=0.7503122325921043`; the exact permutation tail is
`60,852 / 2,421,619,200`, or `p=2.5128641200069772e-5`.  The detector-positive
table is `7/8` in the true upper tail and `0/8` outside it, for risk difference
`0.875`, infinite sample odds ratio, and Fisher `p=1/1430`.

The joint predeclared decision is therefore
`PASS_USABLE_EXTREME_JUMP_DETECTOR`.  Half-lattice recovery has evidence here
as an extreme-jump detector rather than only a point finder.  This does not
authorize a serious rank-32 point search by itself: the existing completed
residual 2-Selmer quotient gate remains mandatory on the same minimal curve.

## Public-complement opening and new rank information

The blind v2 artifact was sealed with SHA-256
`0699b53c2bc7d77673231bc0d377dc725880efd26a18d1eb2af613d28578c165`
before the public complement was opened.  Fourteen final blind bases are
integrally contained in their displayed public subgroups.  Curve 539 has the
same blind and displayed rank 23 but three blind basis points are not recovered
as integral combinations of the displayed basis; no finite-index equality or
inequality is inferred from that bounded relation failure.

For curve 478 the blind final basis has exact finite-reduction rank 23.  Hence
the experiment improves its unconditional lower bound from the displayed
`rank >= 21` to

`rank E_478(Q) >= 23`.

This is a lower bound, not an exact Mordell--Weil rank or saturation theorem.

## Replay

The expensive blind run is preserved rather than included in the normal test
suite.  Preparation and protocol freezes are:

```sh
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/prepare_r17_refresh_jump_ladder_inputs.sage

python3 elliptic-curves/cas/build_r17_refresh_jump_ladder_protocol.py
python3 elliptic-curves/cas/build_r17_refresh_jump_ladder_protocol_v2.py
```

The v2 blind run and post-freeze checks are:

```sh
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/run_r17_refresh_jump_ladder_blind_v2.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/verify_r17_refresh_jump_ladder.sage \
  --blind artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_blind_v2.json \
  --protocol artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_protocol_v2.json \
  --output artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_verification_v2.json

python3 elliptic-curves/cas/analyze_r17_refresh_jump_ladder_v2.py --check
.venv/bin/python -m unittest \
  elliptic-curves/tests/test_r17_refresh_jump_ladder.py
```

## Claim boundary

The exact scores are ranks of discovered groups modulo specialized MW17.
Displayed jumps are exact ranks of the certified displayed-subgroup quotients,
not assertions that those subgroups are full.  Bounded misses prove no point
absence, Selmer structure, saturation, or rank upper bound.  The exact
permutation calculation conditions on this fixed atlas-refresh panel; it is
not a population-sampling theorem or evidence for unrelated fibrations.
