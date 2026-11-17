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

## Complete two-pair cubic Image theorem

The [canonical proof](foundations/extended-geometry/TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md)
proves the full sixteen-parameter cubic stratum, the sharp mixed bound
`m > 3 deg_z Q`, and the least failing bihomogeneous balanced degree four.
It combines explicit boundary volume extension, a sourced proper V-filtration
lemma, and the exact eleven-parameter polynomial certificate.

```sh
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_structure.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_affine_slice.py
.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_phase_slice.py
```

The first two replays support the complete written argument; they do not
computer-verify its imported V-filtration theorems. The phase-slice check and the
optional repeated-profile checkers have their stated smaller arithmetic scopes.
The [primary-source audit](foundations/extended-geometry/two_pair_image_cubic_primary_sources.json)
and [completion audit](foundations/extended-geometry/two_pair_image_cubic_completion_audit.json)
record dependencies, actual replays and the distinction between internal proof,
external review and formal verification.

The current graph also indexes source declarations, cited resources and exact-domain
bounds. Use `python3 research/scripts/research.py coverage`, `references`,
`bounds --opportunities` and `opportunities` to locate omissions and affected work.
Equivalent statements retain separately proved directions in acyclic components;
proofs never acquire assurance flags through propagation.
