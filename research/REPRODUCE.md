# Reproducing the results

Choose a claim in [MATH_STATUS.json](MATH_STATUS.json) or the
[research catalogue](index/README.md), read its canonical source, and identify
the smallest check that establishes the required endpoint.

## Navigation and small regression checks

These commands use the Python standard library and can run from the repository root:

```sh
make render-navigation
make check-navigation
python3 research/scripts/research.py show EC-CURVE302-RECOVERED-MW17-PARENT
```

They validate navigation, claim metadata, checker-source hashes, lesson references
and preserved snapshots. Small regression controls include the literal NS0031
period-group counter-witness. They do not rebuild curve inventories, replay
large certificates or run CAS searches.

## Exact replay

The [full command catalogue](replay/CATALOGUE.md) preserves the former long guide.
Its EC/K3 commands run from `research/`. Other programmes moved to
[their archive](archive/non-elliptic/README.md); their original relative paths
now resolve from `research/archive/non-elliptic/`. Commands are preserved as
historical evidence, not a current execution queue.
The [elliptic-curve replay guide](elliptic-curves/REPRODUCE.md) and canonical
proof notes specify the relevant inputs, versions and limits.

The level-474 H3 rational base has a current Sage 10.9
quadratic-Chabauty and finite Mordell--Weil-sieve certificate:

```sh
sage elkies-k3/scripts/certify_h3_level474_rational_points_qc.sage --check
```

It fetches only a SHA-pinned upstream implementation when its local cache is
empty. The historical Magma output remains unavailable; the source-family
replay still has its separately recorded input boundary.
Do not replace those missing inputs by downloading a newer external copy.

For NS0031, the [current proof](elkies-k3/NS0031_FULL_STABLE_MARKING_OBSTRUCTION_2026-09-14.md)
excludes full rational rank-19 markings using the full stable group, its canonical
genus-10 curve over `Q`, and a complete rational-point calculation on a genus-2
quotient. Its two replay commands run from the repository root:

```sh
sage -python research/elkies-k3/scripts/certify_ns0031_stable_marking.sage --check
magma research/elkies-k3/scripts/certify_ns0031_genus2_rational_points.m
```

The Sage command recomputes the finite group, modular-symbol and Sturm checks
and audits the retained Magma proof. The Magma command reruns the unconditional
global descent and elliptic Chabauty; its recorded run used Magma 2.29-10.
The [packet](artifacts/generated-results/elkies-k3-ns0031-stable-marking-v1/README.md)
preserves the raw output and input/output hashes. The written normalizer proof,
canonical period-map theorem and CM identification remain explicit theorem inputs.

The [earlier note](elkies-k3/NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
preserves the invalidated norm-one argument and its exact counter-witness.
`make verify-ns0031-arithmetic` checks only the retained finite modular arithmetic;
`make verify-ns0031-period-group` checks the missing reflection. Neither of those
older checks alone proves the corrected rational-marking obstruction.

For the three retained six-fibre cohorts, `make verify-fresh6-retained-ranks`
checks 37 subgroup endpoints using saved points and sufficient reduction
primes. The [replay boundary](elliptic-curves/notes/FRESH6_RETAINED_SEED_COHORT_2026-09-09.md#retained-rank-replay)
distinguishes this small check from full search/map replay, a cached receipt
check and conductor factorization.

1. Read the full claim scope and its checker/software-lock fields.
2. Check that the retained certificate and input paths exist.
3. Read the replay command and its implementation before running it. A
   `verify` or `check` name alone does not guarantee a cheap witness replay.
4. Run only the selected check. Missing local artifacts remain missing until
   reconstruction is explicitly part of the task.

The full catalogue remains the canonical source for the two historical H3 q24
certificates `EC-K3-H3-Q24-EXACT-SECTION` and `EC-K3-H3-Q24-QQ-D12`.
Its status-consumer markers and relative links moved with the material.

## Environments and broad suites

For a version-matched Python replay, first select the interpreter recorded in
[.python-version](.python-version), then create the environment from `research/`:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install --requirement requirements.txt
```

The default `make check` and navigation checks use system Python. Explicit
mathematical replay targets use `research/.venv/bin/python` when present and
otherwise the existing root `.venv/bin/python`. Direct commands in the historical
catalogue assume the former path.  Check the interpreter version before
claiming a matching replay; a compatible environment or a passing metadata
check is not evidence that the recorded runtime was reproduced.

Sage, PARI and native workers are needed only by their named replays.
The active [Makefile](Makefile) exposes navigation checks and explicit EC replay.
The full old [Makefile](archive/non-elliptic/Makefile) is preserved in the archive.
Broad verification and research campaigns require their own scope and budget.

<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 71290f86a475074b -->
