# Verification

From the repository root:

```sh
make check
make render-navigation
```

`make check` validates current result sources, checker hashes, deposited-file checksums, paper metadata, local Markdown links and the small tooling controls. It does not repeat the mathematical research calculations.

The research index, topic README pages and status overview are generated from `MATH_STATUS.json`. Run `make render-navigation` after changing that catalogue; `make check` rejects stale summaries.

Use `python3 research/scripts/research.py list`, `show RESULT-ID` or `search QUERY` to inspect the compact current catalogue. Topic pages link the main proofs and retained verification inputs.

The [published bundles](../results/README.md) retain the exact deposited PDF/source bytes. [Draft metadata](../unpublished/README.md) records the build command and remaining review work. Formal projects have their own pinned toolchains and scope pages. Python proof-checker dependencies are in [requirements.txt](requirements.txt).

Read a checker and its requirements before selecting a mathematical replay. Retained receipts of internal checks are separate from independent replay and external review.
