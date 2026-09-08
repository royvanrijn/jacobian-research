# Four MW17 bases through the curve302 adaptive cascade

**Experiment in progress; no terminal conclusion yet.**

Starting commit: `5257d634a390b3ce6d65790d93a7817cfea74b3b`.
The checked-out head subsequently changed through concurrent repository work;
this experiment continues with its original frozen source and fixture hashes.
Mathematical status remains authoritative through `research/MATH_STATUS.json`.

This is the follow-up to the [blinded generic reconstruction](DET1092_BLIND_MW16_RECONSTRUCTION_2026-09-08.md).
The question is whether alternative bases of the same recovered integral MW17
subgroup materially change the **specialized 17→31 search on curve302**.
The experiment does not select new fibres or construct a different elliptic
curve. All four starting groups are equal, with exact unimodular certificates.
Their common certified conductor has163 digits, approximately9.058040×10^162;
see the [user-requested conductor diagnostic](../../artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1/conductor_diagnostic.json)
and authority `ECR31`. A basis change cannot improve this invariant.

## Frozen comparison

The [roster](../../artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1/roster.json)
and [protocol](../../artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1/protocol.json)
were frozen before starting any arm:

| Arm | Generic basis change |
|---|---|
| Original | Identity |
| Recovered P14 | Replace P14 by P14−P16 |
| Recovered P11 | Replace P11 by −P2−P3+P4+P7−P9+P11−P12+P17 |
| Random unimodular | Twelve elementary row additions, deterministic seed10921731 |

The two recovered sections are checked exactly over `Q(t)` against their
previously frozen reconstruction certificates before specialization. The
random operations and all four exact matrices/inverses are in the fixtures.
These choices use completed generic reconstruction results; this is not a
prospectively validated rule for choosing a productive unknown family.

Every worker starts with only generic points specialized at `t=0`. Each uses
its own transformed generic orbit table, own discoveries and own checkpoint
history. A Python artifact-access guard rejects other arms and historical
exceptional-point artifacts. This is an audited process input boundary, not a
claim of designer ignorance or an operating-system sandbox. The V3 selector
was developed by retrospective calibration; that boundary is preserved.

Generic orbit words are transported by `U^{-1}`. In each epoch the height metric
is rounded once, at scale10^6 and384-bit arithmetic, in the original generic17
basis followed by that arm's discoveries; the integer form is then transported
exactly by `diag(U,I)`. Thus a rebase does not accidentally change the orbit
population or the underlying rounded metric when the extension points agree.
Full extension enumeration, exact CVP refinement, chart selection and actual
centre deduplication otherwise follow V3. Babai and LLL choices can still depend
on the displayed basis. Exact CVP refers to the rounded integer form, not an
exact canonical-height optimum.

Common limits are4096 charts,20 epochs, target lower bound31, height125000,
10 seconds per PARI chart,2,000,000 CVP nodes per refinement,7200 seconds of
worker CPU and14400 seconds elapsed per arm. Each process has a12GiB address
space limit. The four arms run concurrently with single-threaded BLAS. Only
finite group tables are cached within an arm; there is no cross-arm search
cache. Each gain is checked using exact finite group certificates with torsion
exclusion, then stale charts are discarded and the landscape rebuilt. Complete
no-gain epochs stop; timeouts and failures remain censored evidence.

## Reproduction and evidence

The preparation script is
[`det1092_basis_cascade_prepare.sage`](../cas/det1092_basis_cascade_prepare.sage);
the frozen search engine is
[`det1092_basis_cascade_v1.sage`](../cas/det1092_basis_cascade_v1.sage).
The [supervisor](../cas/det1092_basis_cascade_supervise.py) refuses a second
launch of an already launched arm. Preparation also refuses to overwrite an
existing package. Fresh reproduction therefore requires an isolated checkout
and fresh output directories, while preserving the current experiment.

Raw checkpoints, full extension arrays, exact CVP records, all selected and
searched charts, transcripts, failed/censored statuses and rank audits are
under `research/artifacts/local/elliptic-curves/det1092-basis-cascade-v1/`.
The checked package is under
`research/artifacts/generated-results/elliptic-curves/det1092_basis_cascade_v1/`.
Post-run packaging compares new equations against historical charts only after
all four arms are terminal. It cannot influence any worker.

The independent [raw replay](../cas/det1092_basis_cascade_replay.sage)
checks seed/orbit/metric transport, every retained point/map witness, the
literal chart prefix, and every promoted finite rank certificate. The
[portable checker](../cas/det1092_basis_cascade_compact_check.sage) checks
compressed evidence without requiring the raw run directory and proves the
universal inverse chord identities used to classify the quartic presentations.
Neither replay reruns the finite point search or asserts exhaustive CVP
coverage beyond the frozen worker's retained evidence.
