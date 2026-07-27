# AGENTS.md

This is a research repository of theorem notes, exact certificates, exploratory
calculations, formalizations, and papers about polynomial Keller maps. Keep
changes narrow, reproducible, and explicit about what has actually been proved.

## Repository map

- `README.md` is the main mathematical overview; `REPRODUCE.md` is the command
  catalogue.
- `verified/`, `cancellation/`, `extended-geometry/`, and `plane-jc/` contain
  active mathematical notes. `formal/` contains Lean developments.
- `scripts/` contains checkers and experiments. Prefer extending a nearby
  script over introducing a second implementation of the same calculation.
- `papers/` contains manuscripts; `archive/` is historical and should not be
  treated as a current source.
- Generated results belong in `artifacts/generated-results/`, with their
  reproducing command documented.

## Mathematical status

- `MATH_STATUS.json` is the sole status authority. `STATUS.md` is generated;
  update it with `python3 scripts/render_status.py`, never by hand.
- Preserve the distinction between theorem, conditional result, computation,
  experiment, and open problem. A successful bounded search is not a proof.
- Put proofs in their canonical note or paper. Other documents should link to
  that source instead of copying the argument or maintaining a second status.
- Do not silently refresh pinned certificates or generated artifacts. Record
  the command, parameters, software assumptions, and any changed hash.

## Editing and checks

- Preserve unrelated work in the tree and avoid broad mechanical rewrites.
- Keep Markdown links relative and run the narrowest relevant verifier first.
- Python verification normally uses `.venv/bin/python`; dependency-free
  independent checks deliberately use `python3`.
- For Python or documentation changes, run `make check`. For the foundational
  certificate, also run `make verify-minimal`. Use the targeted commands in
  `REPRODUCE.md` for specialized or expensive calculations rather than running
  the full suite by default.
- Report commands that were not run, especially when they require Singular,
  Macaulay2, Lean, Julia, LaTeX, or a long symbolic computation.
