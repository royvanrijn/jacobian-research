# Reproducing the results

Choose a claim in [MATH_STATUS.json](MATH_STATUS.json) or the
[research catalogue](index/README.md), read its canonical source, and identify
the smallest check that establishes the required endpoint.

## Navigation checks: no research calculations

These commands use the Python standard library and can run from the repository root:

```sh
make render-navigation
make check-navigation
python3 research/scripts/research.py show EC-CURVE302-RECOVERED-MW17-PARENT
```

They validate navigation, claim metadata, checker-source hashes, lesson references
and preserved snapshots. They do not replay mathematical certificates, discover
factors, rebuild curve inventories or run CAS searches.

## Exact replay

The [full command catalogue](replay/CATALOGUE.md) preserves the former long guide.
Its EC/K3 commands run from `research/`. Other programmes moved to
[their archive](archive/non-elliptic/README.md); their original relative paths
now resolve from `research/archive/non-elliptic/`. Commands are preserved as
historical evidence, not a current execution queue.
The [elliptic-curve replay guide](elliptic-curves/REPRODUCE.md) and canonical
proof notes specify the relevant inputs, versions and limits.

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
