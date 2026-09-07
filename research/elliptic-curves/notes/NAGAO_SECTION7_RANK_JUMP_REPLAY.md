# Nagao section-7 quotient fingerprint and score replay

## Exact quotient response

Nagao's printed parameter is `t=5081/94`; the repository constructor uses
`T=2t=5081/47`. The arithmetic generic rank is exactly 12 and the fibre has
a certified rank-20 subgroup. The generic basis used here is the eleven
visible sections with indices 0 through 10 followed by the `plus-7/27`
section. Exact group-law replay embeds those twelve points in the saturated
rank-20 basis.

The resulting displayed quotient has:

- free rank 8;
- Smith factors `1,2,2,2,2,2,2,2,2,2,2,2` and saturation index `2048`;
- tensor dimensions `19,8,8` over `F_2,F_3,F_5`;
- eight complete numerical projected-height successive minima.

The extra mod-2 dimension is Smith torsion in the quotient of the displayed
lattices; it must not be identified with nineteen free Mordell--Weil
directions.

The archived bounded degree-two searches returned 224 distinct exact
relations in the rank-20 basis. Modulo the generic subgroup they span all
eight free quotient directions over `Q` and all nineteen tensor-quotient
dimensions over `F_2`. They occupy 170 quotient classes; 40 returned points
are in the zero mod-2 quotient class. This is complete for the returned
bounded searches, not a complete atlas of every degree-two cover. No
degree-three atlas is available. The archived degree-four/`ell2cover` attempt
ended with no output at its strict timeout, so its value is missing, not zero.

The machine-readable fingerprint is
[`../../artifacts/generated-results/elliptic-curves/nagao_section7_rank_jump_fingerprint_v1.json`](../../artifacts/generated-results/elliptic-curves/nagao_section7_rank_jump_fingerprint_v1.json).

## Complete frozen-box replay

The historical global scan exhausted all 18,244,819 positive primitive
parameters `T=a/b` with `1<=a<=30000` and `1<=b<=1000`. Its ordinary Nagao
score used primes 5 through 199 for training and the disjoint band 211 through
397 for validation. The score formula and integer scale were frozen in the
archived scanner and do not use point-search or Mordell--Weil labels.

Replaying the known rank-20 fibre against every parameter gives:

| score band | exact position | population fraction |
| --- | ---: | ---: |
| training 5--199 | 9,041,935 | 0.4955 |
| validation 211--397 | 755,065 | 0.04138 |

Thus the two bands disagree substantially. The validation placement is a
useful retrospective development signal, but even it misses a one-percent
candidate budget. The training placement is essentially median. This replay
is not a prospective holdout: the fibre was known historically, although the
archived global frontier explicitly excluded it before selection. Every
other fibre remains censored; the replay supplies no negative rank labels.

The compact replay artifact is
[`../../artifacts/generated-results/elliptic-curves/nagao_section7_rank_jump_replay_v1.json`](../../artifacts/generated-results/elliptic-curves/nagao_section7_rank_jump_replay_v1.json).

## Reproduction

```sh
PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_nagao_section7_rank_jump_fingerprint.py --check

PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_nagao_section7_rank_jump_replay.py --check
```

The rank statements, embedding, Smith structure, and returned cover-point
relations are exact. Canonical heights are numerical at 80-digit PARI
precision. The quotient is relative to the certified rank-20 subgroup; no
claim is made that it is the full `E_T(Q)`.
