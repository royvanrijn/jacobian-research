# Relative 2-Selmer pipeline for the compact R17 family

## Current outcome

The basis-level inputs cover the held-out rank-21 mechanism control, the
rank-25--28 controls, and the first ten candidates in the frozen height-10000
weakest-block Nagao ranking. The original Magma supervisor records all fifteen
jobs as `backend_unavailable`. An open-source Sage/PARI replacement is now
implemented and validated end to end on a small curve; the large R17 controls
remain resource-bounded computations. No incomplete run supplies a Selmer
dimension, unrealized Selmer class, exact-rank statement, or point-search
success.

The first pinned open run on the rank-21 control used a 300-second wall limit,
a 4 GB RSS limit, a 2 GB PARI stack, and all twelve proved discriminant-prime
hints. It timed out inside `ellrankinit` at 440,283,136 bytes peak observed
RSS, before BNF certification or `ell2cover`. This is a measured backend
bottleneck, not a Selmer result.

The same frozen method was applied without public-point hints to the top
high-Nagao candidate `t=-5643/6760`. Its 120-second diagnostic run also timed
out inside `ellrankinit`, at 230,608,896 bytes peak observed RSS. No candidate
cover search was reached.

The exact input manifest is
[`elkies_2026_relative_2selmer_suite_inputs_v1.json`](../../artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json).
The host/backend audit is
[`elkies_2026_relative_2selmer_suite_run_v1.json`](../../artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_run_v1.json).

| case | `t` | certified known rank | held-out quotient directions |
| --- | ---: | ---: | ---: |
| rank-21 mechanism control | `3/8` | 21 | 4 |
| rank-25 control | `-2/377` | 25 | 8 |
| rank-26 control | `-308/251` | 26 | 9 |
| rank-27 control | `2456/135` | 27 | 10 |
| rank-28 control | `-9529/5471` | 28 | 11 |

These are lower-bound controls.  The displayed ranks are not asserted to be
exact until a matching upper bound exists.

## Frozen computation

### Open-source path

`run_elkies_2026_relative_2selmer_open.py` uses PARI's `ell2cover`, whose
output is a basis of everywhere locally soluble binary-quartic 2-covers with
maps to the elliptic curve. For each selected case it:

1. starts an isolated Sage worker with the minimal model but no generic or
   exceptional points;
2. calls `ellrankinit` and requires `bnfcertify` to return one on every
   cubic-field BNF in the rank context;
3. calls `ell2cover`, records every basis quartic and its map, and runs a
   bounded `hyperellratpoints` search without public point hints;
4. reloads the exact generic and held-out points only after the worker exits;
5. uses exact good-prime quotients to express recovered cover images in the
   known Mordell--Weil basis; and
6. emits point-to-Selmer rows when the blindly recovered basis classes span
   the known image.

A full `ell2cover` return plus successful BNF certification is an actual
2-Selmer basis. A bounded quartic-search miss does not prove that its class is
unrealized. PARI supplies the basis quartics but does not expose addition of
arbitrary cover classes, so explicit covers for non-basis quotient
combinations still use the repository's cubic-etale intersection-of-quadrics
layer.

### Optional Magma cross-check

For each specialization the generated Magma job performs the following steps.

1. Reconstruct the global minimal fibre and the seventeen specialized generic
   sections from the pinned compact R17 model.
2. Construct multiplication by two and call
   `SelmerGroup([2] : Bound := -1, Raw := true)`, the shared-map form of
   `TwoSelmerGroup`. `Bound=-1` requests unconditional class-group data; no
   GRH class-group bound is installed.
3. Use `DescentMaps([2])` and the returned `AtoS` map to compute the actual
   Selmer coordinate `AtoS(mu(P_i))` of every generic section.  The job aborts
   unless these rows have rank 17.
4. Extend those rows to a basis of the full 2-Selmer group and thereby obtain
   an explicit basis of

   ```text
   Sel_2(E/Q) / image(<P_1,...,P_17>).
   ```

5. Before declaring any public exceptional points, materialize quotient
   representatives with `TwoCover`, record their quartic equations, and run a
   bounded `Points` search.  The generic x-coordinates are the only Selmer
   hints.  Thus the measured recovered quotient span is genuinely blind to the
   public extra-point coordinates.
6. Declare the held-out control points only after the blind phase, map them
   into the same quotient basis, compute their span, and label every stored
   cover class as inside or outside that known rational span.

If a quotient has at most 255 nonzero classes, every nonzero class is built
and searched.  Above that threshold the job builds and searches a canonical
quotient basis.  This preserves a spanning set of unexplained directions
without pretending that an exponential all-class search is feasible.  The
parser explicitly records whether enumeration was exhaustive.

The parser rejects missing stages, source-hash changes, inconsistent
dimensions, incomplete generic/quotient bases, and class-count mismatches.  A
bounded cover-search miss remains a negative experiment, not evidence for a
Tate--Shafarevich class.

## High-Nagao application

The frozen method has been instantiated on these first ten prospective
candidates; the open backend has been executed on the first candidate only:

```text
-5643/6760, 1452/7817, 4298/8873, -7634/2859, -841/8544,
461/4420, 6695/1353, 1217/151, 9783/7559, 9446/3605.
```

Their exact minimal models, specialized generic sections, Nagao records, and
program hashes are in the input manifest.  They remain heuristic candidates
until the complete descent returns.  In particular, no candidate is promoted
to point search merely because it was selected by Nagao score.

## Backend calibration

Input reconstruction and generation of all fifteen source-pinned jobs took
about 13.1 seconds in the first local run.  This is an input-generation
benchmark, not a descent benchmark.

Earlier open-source entry points did not provide the requested basis-level
result on the controls:

- eclib/mwrank failed on the rank-21 minimal model in about one second with
  `lower bound on c too large`, before returning a Selmer rank;
- PARI `ellrank`, with all 21 certified points supplied, did not return within
  a strict 60-second calibration;
- PARI `ellrankinit` with all twelve proved rank-21 factor hints did not return
  within the pinned 300-second open-suite run and peaked at 440,283,136 bytes
  observed RSS;
- the already pinned rank-28 eclib and PARI attempts both timed out after 300
  seconds, and the factor-supplied PARI attempts timed out after 600 and 1800
  seconds.

None of these outcomes is a Selmer upper bound. PARI 2.17.3 additionally
exposes `ell2cover`, which supplies the previously missing full locally
soluble cover basis. The open runner separately certifies the field BNF and
reconstructs the point embedding from blindly recovered cover images. Magma's
raw interface remains an optional independent cross-check.

The mathematical reduction through a cubic etale algebra and an `S`-class/
`S`-unit computation follows Schaefer--Stoll,
[*How to do a p-descent on an elliptic curve*](https://mathe2.uni-bayreuth.de/stoll/papers/p-descent-long.pdf).
The explicit-cover minimization and reduction layer is consistent with
Cremona--Fisher--Stoll,
[*Minimisation and reduction of 2-, 3- and 4-coverings of elliptic curves*](https://arxiv.org/abs/0908.1741).

## Replay

Generate the five controls and first ten candidates:

```bash
python3 elliptic-curves/cas/build_elkies_2026_relative_2selmer_suite.py \
  --output-dir artifacts/local/elliptic-curves/elkies-2026-relative-2selmer-suite-v1 \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --candidate-count 10 --search-bound 1000 \
  --enumerate-class-limit 255 --overwrite
```

Run the rank-21 control with the open-source backend and explicit wall/RSS
limits:

```bash
python3 elliptic-curves/cas/run_elkies_2026_relative_2selmer_open.py \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --output-dir artifacts/local/elliptic-curves/elkies-2026-relative-2selmer-open-v1 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_open_rank21_300s_v1.json \
  --case control-r21-t3_8 --timeout-per-case 300 \
  --rss-limit-bytes 4000000000 --pari-stack-bytes 2000000000 \
  --search-bound 1000 --certificate-prime-bound 1000 --overwrite
```

Repeat `--case` for controls or candidates, or use `--controls-only`. The
output is self-classifying: only `COMPLETE_CERTIFIED_PARI_TWO_SELMER_BASIS`
contains a full basis and quotient calculation.

For an optional licensed Magma cross-check, supervise every case with explicit
wall/RSS limits:

```bash
python3 elliptic-curves/cas/run_elkies_2026_relative_2selmer_suite.py \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --log-dir artifacts/local/elliptic-curves/elkies-2026-relative-2selmer-suite-v1/logs \
  --output artifacts/local/elliptic-curves/elkies_2026_relative_2selmer_suite_run.json \
  --timeout-per-case 86400 --rss-limit-bytes 16000000000 --overwrite
```

Only after every job completes, parse the source-matched transcripts:

```bash
python3 elliptic-curves/cas/parse_elkies_2026_relative_2selmer_suite.py \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --log-dir artifacts/local/elliptic-curves/elkies-2026-relative-2selmer-suite-v1/logs \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_results_v1.json
```

The parser cannot turn a timeout, backend failure, or partial log into a
result.  `MATH_STATUS.json` therefore remains unchanged.
