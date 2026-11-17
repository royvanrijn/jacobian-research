# Verification and replay

Run maintenance from the repository root:

```sh
make check-navigation
make check
```

These check the ledger, source references, navigation and small tooling controls. They do not prove every statement or repeat the research calculations.

Each [published paper](../results/README.md) includes its deposited PDF, source and metadata. Its README states the result and verification evidence. [Unpublished manuscripts](../unpublished/README.md) state their remaining obligations.

Research Python dependencies are pinned in [requirements.txt](requirements.txt). A standard Python 3 installation is sufficient for navigation. Particular proof checks may require SymPy, Sage, Singular, Lean or Julia; read their prerequisites first.

The [foundational library](foundations/README.md) retains explicit proof targets in its Makefile. Use `make -C research/foundations TARGET` for a deliberately selected target. The long [replay catalogue](replay/CATALOGUE.md) records historical commands; inspect current paths, inputs and claim scopes first.

Do not rerun successful sealed checkers, large searches or absent local caches merely to refresh navigation. [Source-review coverage](knowledge/PARTIAL_REVIEW.md) exposes missing and stale reviews.
