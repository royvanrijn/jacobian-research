# Fermigier quotient fingerprints and frozen-score replay

## Outcome

The two certified Fermigier/Mestre controls now have quotient-first
fingerprints relative to the proved generic rank-12 subgroup.

| fibre | free quotient lower bound | dimensions mod 2/3/5 | generic-subgroup index in displayed group | first / last quotient minimum | bounded degree-two visible span |
| --- | ---: | ---: | ---: | ---: | ---: |
| E22, `u=19754/39` | 10 | `22/11/10` | 24,576 | `9.267916 / 12.049441` | 10 |
| rank-20 near miss, `u=28917/20` | 8 | `8/8/8` | 1 | `22.171961 / 34.049852` | 8 |

The E22 Smith factors are eleven copies of `2` and one `12`. This is the
reason its mod-2 quotient dimension is 22: twelve dimensions are finite Smith
torsion, not twelve additional free rank directions. The rank-20 embedding is
primitive. In both fibres the declared signed weight-one/two direction ball
spans the whole displayed free quotient. Modulo sign its intersection matroid
is the cycle matroid of `K_11` for E22 and `K_9` for the rank-20 fibre; the
artifact records the exact circuit census. Degree-three and degree-four cover
visibility remain missing, not zero.

The complete height Grams, successive minima, short vectors, embeddings,
Smith data, and cover-intersection records are in
[`fermigier_rank_jump_fingerprints_v1.json`](../../artifacts/generated-results/elliptic-curves/fermigier_rank_jump_fingerprints_v1.json).

## Complete frozen-score replay

The historical global scan supplies a deterministic common population:

```text
every primitive T=a/b with 0 <= a <= 100000 and 1 <= b <= 1000,
with T and -T identified: 60,815,684 parameters.
```

Four already-declared local score orderings were replayed exhaustively:
discovery-band rank, discovery-band rank plus the fixed repeated-discriminant
term, and the corresponding two held-band orderings. No point-search result or
rank label enters a score. The exact positive positions are:

| score ordering | E22 position | rank-20 position |
| --- | ---: | ---: |
| discovery rank | 2,755,127 | 3,070,200 |
| discovery composite | 2,755,324 | 3,070,421 |
| held rank | 3,556,159 | 12,535,323 |
| held composite | 3,556,250 | 12,535,549 |

Thus every ordering has zero recall at budgets 10, 100, 1,000, 8,000, 16,133,
and 100,000. This is a useful negative result: these frozen local scores do
not retrieve escape from the generic lattice in the Fermigier family and
must not be used to launch a new candidate search.

This replay is labelled development/retrospective, not held out. E22 was a
known calibration when the historical global score was designed. The held
prime band is disjoint from the discovery band, but that fact alone does not
make the known fibre a prospective holdout. All 60,815,682 other fibres remain
censored; none is assigned quotient rank zero, so no accuracy, specificity,
ROC, or precision claim is made.

The compact replay is
[`fermigier_rank_jump_replay_v1.json`](../../artifacts/generated-results/elliptic-curves/fermigier_rank_jump_replay_v1.json).

## Replay

```bash
PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_fermigier_rank_jump_fingerprints.py

PYTHONPATH=elliptic-curves:elliptic-curves/cas python3 \
  elliptic-curves/cas/build_fermigier_rank_jump_replay.py

.venv/bin/python elliptic-curves/scripts/run_rank_jump_laboratory.py
```

Each builder also accepts `--check`. The score replay recompiles its small
C++17 enumerator and repeats the complete box, so its check is intentionally
slower than the fingerprint check.
