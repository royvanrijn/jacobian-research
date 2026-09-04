# Blind multi-stratum half-lattice jump ladder

<!-- status-consumer: EC-K3-R17-REFRESH-BLIND-JUMP-LADDER a2d7034fb8977c18 -->

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
prove that any selected candidate has rank 32.  In production, however, a
completed residual 2-Selmer computation is an exclusion/closure tool rather
than a universal search prerequisite: only a certified upper bound below 32
vetoes a fibre, while incomplete descent affects scheduling only.

## Fibration and `j`-class sensitivity

Write `S` for the exact blind rank of the discovered subgroup modulo the
specialized MW17, and `q` for the exact free rank of the later-opened displayed
public subgroup modulo MW17.  The pooled confirmatory endpoint above is the
predeclared `S >= 10` versus `q >= 10` test.  The `q >= 11` cut and all splits
in this section were computed after unblinding, so they are sensitivity
analyses rather than additional confirmatory endpoints.

Pooled over all sixteen fibres, `S >= 10` occurs on `6/7` rows with `q >= 11`
and `1/9` rows below eleven.  The risk difference is `0.7460`, the sample odds
ratio is `48`, and the one-sided Fisher probability is `4/715`.

The fibration split is:

| fibration frame | rows | `tau_b(S,q)` | exact ordinal `p` | `S>=10` in `q>=10` vs below | `S>=10` in `q>=11` vs below |
| --- | ---: | ---: | ---: | --- | --- |
| published R17 | 13 | `0.7033` | `2656/4324320` | `7/8` vs `0/5`, `p=2/429` | `6/7` vs `1/6`, `p=43/1716` |
| alternate Q80 | 3 | `1` | `1/6` | no `q>=10` row | no `q>=11` row |

Conditioning the ordinal randomization on those two frame labels gives the
block-restricted `tau_b=0.7163` and exact one-sided
`p=7082/25945920`.  Thus the pooled ordinal signal is not created solely by
the frame imbalance.  It does not, however, validate the extreme-tail rule on
alternate Q80: its three rows are exactly `(S,q)=(3,3),(6,6),(8,8)`.

At the finer rational-`PGL2` `j`-map level, only class `08234` contains both
tail and non-tail rows.  Its seven pairs are

```text
(6,4), (8,8), (10,10), (10,11), (11,11), (11,11), (11,11).
```

They give `tau_b=0.8767`, exact ordinal `p=2/210`, and the `q>=10` table
`5/5` versus `0/2` with Fisher `p=1/21`.  At the stricter `q>=11` cut the
table is `4/4` versus `1/3`, with `p=1/7`; that threshold is not independently
resolved by this small within-class panel.  Class `07ca9` has three tail-only
rows `(12,12),(0,11),(11,11)`, including the explicit false negative, and the
other four represented classes have no `q>=10` row.

An exact randomization that permutes `q` only within the six `j`-classes gives
block-restricted `tau_b=0.8804` and ordinal `p=2/2520`.  The corresponding
conditional high-`S` tail probabilities are `3/63=1/21` for `q>=10` and
`9/63=1/7` for `q>=11`.  These calculations exclude every cross-class pair;
their limitation is support, not an observed reversal.

The operational decision is therefore one-sided and family-scoped.  High `S`
may participate in scheduling among candidates in the calibrated
norm-twelve R17 setting, including within `08234`; it is not merely a
post-selection point finder there.  Low `S` must not veto a candidate, because
the searches are bounded and curve 544 has `q=11,S=0`.  Nor does this panel
authorize transporting the rule to alternate Q80, an unrepresented `j`-class,
or a changed lattice/basis.  Those settings need their own high-`q` blinded
controls.  Every follow-up still needs a declared finite search budget and
exact point-independence certification; a completed residual-Selmer upper bound
is required for exact-rank closure, not for proving rank at least 32 from
points.

The completed v2 artifact is not retroactively changed to repair its bootstrap
asymmetry.  Production MW17-jump-v2 instead has a separately frozen,
outcome-blind one-in-eight rescue arm.  Assigned clean-zero fibres receive the
next 301 generic half-classes in seven batches, and after the first certified
escape any unused slots switch to the original quotient-adaptive policy.  This
keeps the 344-chart total cap while ensuring that a treated zero-gain fibre
does not need a pre-existing quotient direction to receive full exposure.

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
not a population-sampling theorem or evidence for unrelated fibrations.  The
fibration/`j`-class and `q>=11` calculations are explicitly post-freeze
sensitivity analyses.
