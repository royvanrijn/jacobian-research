# Broad matched MW17 search

Status: **implemented prospective CPU campaign; no search result is claimed here.**

The controller `elliptic-curves/cas/run_broad_mw17_search.py` compares certified generic-rank-17 parents on the same frozen rational-address window. It currently ships with two immediately runnable X1092 presets:

- `x1092-curve302`: the recovered curve302 MW17 parent, automatically adapted through its certified short reduced parameter chart while retaining the proved height Gram;
- `x1092-class1`: the independently replayed new class-1 MW17 realization.

Additional compatible short-Weierstrass MW17 parent packets can be supplied as `--parent NAME=/absolute/path.json` without changing the controller.

## What is frozen

For each parent the existing prospective search implementation freezes the parent, source tree, Sage/GP executables, score tables, and exact address window. The broad controller then freezes three disjoint expensive-search arms:

1. **ranked**: highest current Nagao/Mestre-style aggregate score;
2. **control**: score-independent hash-selected controls inherited from the existing prospective protocol;
3. **diversity**: deterministic extrema in standardized cheap feature space (`model_bits`, `discriminant_bits`, `smooth_prime_count`), excluding ranked and control rows.

The diversity arm is an exploratory arm, not a theorem or learned predictor. It exists to gather data away from one scalar score.

The same address interval is scored on every parent. Historical exceptional parameters and points are not construction or selection inputs.

## Default campaign

The defaults intentionally start beyond the two previous 65,536-address class-1 windows:

- address offset: `131072`
- addresses per parent: `262144`
- score primes: through `997`
- ranked expensive fibres per parent: `256`
- independent controls per parent: `64`
- diversity fibres per parent: `64`
- global expensive-search workers: `4`

The existing worker starts at 24 point-search calls and adaptively raises its cumulative allowance at certified rank lower bounds 20, 23, 25 and 27. Every completed batch receives full-cloud reconciliation and two finite-rank certificate implementations. Censored work remains UNKNOWN.

## Start it

From `research/`:

```sh
python elliptic-curves/cas/run_broad_mw17_search.py prepare \
  --folder artifacts/local/elliptic-curves/broad-mw17-search-v1 \
  --workers 4

python elliptic-curves/cas/run_broad_mw17_search.py run \
  --folder artifacts/local/elliptic-curves/broad-mw17-search-v1
```

Preparation performs the cheap scoring and freezes all expensive-search selections **before** `run` dispatches a fibre.

Inspect progress with:

```sh
python elliptic-curves/cas/run_broad_mw17_search.py status \
  --folder artifacts/local/elliptic-curves/broad-mw17-search-v1
```

Graceful stop:

```sh
touch artifacts/local/elliptic-curves/broad-mw17-search-v1/STOP
```

No new fibres are dispatched after `STOP`; active bounded fibres drain. Remove the file and invoke `run` again to resume uncompleted selected fibres.

For a smaller smoke run, use a new folder:

```sh
python elliptic-curves/cas/run_broad_mw17_search.py prepare \
  --folder /tmp/broad-mw17-smoke \
  --window 256 --offset 131072 \
  --ranked 4 --controls 2 --diversity 2 --workers 2
python elliptic-curves/cas/run_broad_mw17_search.py run --folder /tmp/broad-mw17-smoke
```

Run the pure controller tests with:

```sh
python -m pytest elliptic-curves/tests/test_broad_mw17_search.py -q
```

## Interpretation

The primary outputs are per-parent/per-arm certified lower-bound distributions and CPU cost. A null bounded panel is operational evidence only, never a rank upper bound. The scientifically useful comparisons are:

- score-selected vs independent controls at equal expensive-search exposure;
- diversity arm vs the scalar score;
- curve302 frame class 6 vs class 1 on the same X1092 and the same address protocol;
- CPU per certified >=20, >=23, >=25 and >=27 event.

Full cheap score rows remain in each frozen child runtime, so later selector work can use the complete feature population rather than only successful fibres. Any V2 learned selector must be frozen and tested on a disjoint future address window; this campaign must not be retrospectively retuned.
