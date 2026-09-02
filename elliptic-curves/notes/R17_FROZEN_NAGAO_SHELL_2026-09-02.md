# Prospective frozen-Nagao R17 height shell

Date: 2026-09-02

## Outcome

The current compact-R17 weakest-block rule was frozen completely and applied
without a population presieve to the disjoint shell

```text
10,000 < H(t)=max(|a|,b) <= 30,000,   t=a/b,   b>0,   gcd(a,b)=1.
```

The shell contains exactly `972,697,152` primitive parameters. Every parameter
was scored; only the leading prefixes were materialized. The complete stream
took `529.446` seconds on the recorded host.

The frozen target lane, pooled ordinary-Nagao control, and deterministic random
control each contain 128 parameter-disjoint fibres. Complete evaluation of the
preexisting 39,120-bisection atlas gives:

| lane | rows | rows with certified quotient gain | gain sum beyond generic 17 |
| --- | ---: | ---: | ---: |
| frozen weakest block | 128 | 3 | 3 |
| pooled ordinary Nagao | 128 | 3 | 3 |
| deterministic random | 128 | 1 | 1 |

Thus this prospective shell certifies seven new rank-at-least-18 fibres, one
new quotient direction on each. It does **not** show an advantage over pooled
ordinary Nagao in this 128-row comparison. The descriptive success fractions
are `3/128`, `3/128`, and `1/128`; they are not population rank-jump
probabilities.

## Frozen rule and controls

The target comparator is copied verbatim from
[`elkies_2026_compact_t_nagao_positive_control_h10000_v1.json`](../../artifacts/generated-results/elliptic-curves/elkies_2026_compact_t_nagao_positive_control_h10000_v1.json):

1. use the same 102 primes from 19 through 599 in the same three pairwise
   disjoint round-robin blocks;
2. center and population-standardize every good-fibre contribution on
   `P^1(F_p)`;
3. mean-impute a singular local fibre by standardized contribution zero;
4. divide each block sum by the square root of its prime count;
5. rank by descending weakest block, then descending block mean, good-prime
   count, bad-prime count, height, denominator, and numerator.

No prime, window, normalization, imputation rule, tie-break, or cutoff changed.
The old development box is `H<=10,000`; the new shell starts at `H=10,001`, so
the two parameter sets are disjoint by construction.

The ordinary control ranks by the mean of the same three standardized blocks,
then the frozen comparator. The random control consists of the smallest
predeclared deterministic `splitmix64` keys. The lanes are deduplicated in
target, ordinary, random order before any bisection outcome is evaluated.

## Exact gains

| lane | lane rank | `t` | bisection orbit | quotient gain | rank lower bound | strongest reduction bound |
| --- | ---: | ---: | --- | ---: | ---: | ---: |
| frozen | 9 | `28781/3342` | `orbit-11aa9` | 1 | 18 | 997 |
| frozen | 38 | `-26851/25276` | `orbit-1b011` | 1 | 18 | 499 |
| frozen | 121 | `29111/2186` | `orbit-18753` | 1 | 18 | 199 |
| ordinary | 8 | `-11677/17012` | `orbit-0acda` | 1 | 18 | 199 |
| ordinary | 103 | `-19909/14703` | `orbit-12f39` | 1 | 18 | 199 |
| ordinary | 110 | `-25961/9577` | `orbit-01362` | 1 | 18 | 199 |
| random | 119 | `-2337/28069` | `orbit-1037b` | 1 | 18 | 199 |

For every row the two bisection branches lie exactly on the specialized
source fibre, their sum equals the stored trace, and that trace equals its
published R17 basis word. Finite-reduction replay proves that the displayed
point escapes the specialized generic 17-dimensional subgroup.

The frozen-rank-91 fibre `t=-12229/7635` has two nonzero split bisections but is
censored: the base prime ensemble sees generic finite-quotient rank 16, not
17. No gain is promoted for it.

These seven directions live on seven different fibres. They do not combine to
give one rank-24 fibre. A zero atlas outcome says only that the already-known
bisection mechanism supplied no certified direction; it is not a rank-17
claim or an upper bound.

## Progressive depth and the point-search boundary

Depth was allocated only by frozen target rank:

- all 384 rows received the complete bisection-atlas pass and finite-reduction
  replay through prime 199;
- frozen ranks 17--64 were replayed through prime 499;
- frozen ranks 1--16 were replayed through prime 997.

The stronger tiers confirm the rank-9 and rank-38 gains. They do not change
the ranking and were selected before point outcomes were opened.

No unrestricted `ratpoints`, eclib/mwrank, normalized-slope, or two-cover point
search was run. The repository's fail-closed R17 policy requires a completed
same-fibre residual quotient

```text
Sel_2(E_t)/<P_1,...,P_17>
```

of dimension at least 15 before any such search. None of the shell candidates
has that computation, and Magma is absent on this host. This is an
authorization/backend boundary, not a negative search result or rank upper
bound. The newly added relative-2-Selmer suite likewise records backend
unavailability and does not cover these shell candidates.

## Artifacts and replay

The complete ranking and compact exact outcome are:

- [`r17_frozen_nagao_shell_h10001_30000_v1.json`](../../artifacts/generated-results/elliptic-curves/r17_frozen_nagao_shell_h10001_30000_v1.json), SHA-256
  `b8b4b533f69b08680beffb4384302d3b95b04856670e1f781f121682fe50f112`;
- [`r17_frozen_nagao_shell_search_v1.json`](../../artifacts/generated-results/elliptic-curves/r17_frozen_nagao_shell_search_v1.json), SHA-256
  `7f9824e5f3552d929bff9e5380a082cfb11538c46b96a5b4dc66ad9a97bb7de6`.

Generate the ranking and matched cohort with:

```bash
.venv/bin/python elliptic-curves/scripts/run_r17_frozen_nagao_shell.py
```

Run the base exact atlas evaluation with:

```bash
.venv/bin/python elliptic-curves/scripts/label_r17_training_bisections.py \
  --input artifacts/local/elliptic-curves/r17-frozen-shell-h10001-30000-cohort.jsonl \
  --output artifacts/local/elliptic-curves/r17-frozen-shell-h10001-30000-bisection-labels.jsonl \
  --summary artifacts/local/elliptic-curves/r17-frozen-shell-h10001-30000-bisection-labels-summary.json \
  --workers 4 --prime-bound 199
```

The exact tier inputs are the first 16 and rows 17--64 of the already frozen
target lane. Their commands and hashes are recorded in the compact result.
After those replays, rebuild the compact result with:

```bash
.venv/bin/python elliptic-curves/scripts/summarize_r17_frozen_nagao_shell.py
```

This is a new exact instance of the existing R17 bisection-specialization
theorem, not a new general theorem, so `MATH_STATUS.json` is unchanged.
