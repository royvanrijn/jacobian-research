# Verification

From the repository root:

```sh
make check
make render-navigation
```

`make check` validates current result sources, checker hashes, declared current input pins, the acyclic implication graph and its proof-source freshness, cached primary-source digests, deposited-file checksums, paper metadata, local Markdown links and the small tooling controls. It does not repeat the mathematical research calculations.

The research index, topic README pages, scope audit, dependency diagrams and status overview are generated from `MATH_STATUS.json`. Run `make render-navigation` after changing that catalogue; `make check` rejects stale summaries.

Use `python3 research/scripts/research.py list`, `show RESULT-ID` or `search QUERY` to inspect the compact current catalogue. Topic pages link the main proofs and retained verification inputs.

The [published bundles](../results/README.md) retain the exact deposited PDF/source bytes. [Draft metadata](../unpublished/README.md) records the build command and remaining review work. Formal projects have their own pinned toolchains and scope pages. Python proof-checker dependencies are in [requirements.txt](requirements.txt).

Read a checker and its requirements before selecting a mathematical replay. Retained receipts of internal checks are separate from independent replay and external review.

The [dependency guide](DEPENDENCY_GUIDE.md) explains how to record a proof so its
consequences propagate. Use `explain ID` for blockers and proof provenance,
`impact ID` for descendants and `simulate --prove ID` for a read-only hypothetical
improvement. Graph checks do not upgrade formal, replay, review or publication flags.

## Two-pair cubic Image strata

The [current cubic note](foundations/extended-geometry/TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md)
records the full sixteen-parameter coordinate model, parity coverage, exact mixed
slices, and the two missing full-family implications. With the pinned environment:

```sh
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_structure.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_phase_slice.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_affine_slice.py
```

These are exact internal replays of structural identities and the complete seven-
and ten-parameter slice certificates. The odd-parity proof imports the existing
internal Gaussian classification; the checkers do not prove that theorem or cover
the entire mixed cubic family. The sphere-shear argument supplies all-order mixed
vanishing for its classified cases; finite power regressions have only their
recorded bounds.
