# Repository workflow

Keep one clear account of each result and enough evidence to check it. The directory
structure separates mathematical work from manuscript editions and deposited files.
Read [the root instructions](../AGENTS.md) before starting; use
[the mathematics guide](MATHEMATICS.md) when changing a claim or its evidence.

## Where work belongs

| Location | Purpose | What belongs here |
|---|---|---|
| `research/` | Current mathematics by topic | Canonical statements and proofs, open obligations, necessary checkers and exact inputs |
| `research/MATH_STATUS.json` | Structured current mathematical status | Stated scope, canonical sources, retained evidence and separate assurance flags |
| `research/STATUS.md` | Generated navigation | Output of `make render-navigation`; edit its structured source |
| `unpublished/` | Manuscript bundles awaiting revision, review or deposit | Paper source, current PDF, build instructions, metadata, verification scope and next review tasks |
| `results/` | Deposited editions | Exact deposited files, checksums, identifiers, submission metadata and verification scope |
| `working-guide/` | Reusable working standards | Policies and reusable procedures, without per-result progress records |

A research note can contain incomplete arguments and failed approaches when they are
needed to explain an open obligation. Mark their status where they occur. A paper
presents a coherent mathematical contribution that can be read without the worktree.
A deposited edition is immutable: prepare corrections or later editions under
`unpublished/`, retaining the original deposit and its provenance.

Use the topic catalogue to find existing work before adding a new directory or proof.
Do not create competing descriptions named `final`, `final-v2`, `latest` or dated
handoffs. Improve the canonical source and let Git retain superseded text. Preserve
a failed approach when it establishes a useful obstruction or counterexample; state
that conclusion directly instead of retaining an unbounded discovery diary.

## Start each task

1. Read `AGENTS.md`, the relevant topic README, the canonical proof and the applicable
   manuscript metadata. Inspect the actual checker and inputs before interpreting
   an old receipt.
2. Inspect `git status --short` and relevant diffs. Existing changes may belong to
   another task or person; do not reset, reformat, stage or delete them.
3. Identify the intended change, its authoritative source, its consumers and the
   checks it requires. During organization work, preserve mathematical scope and
   avoid starting new research calculations.
4. Reuse the current build and dependency instructions. Do not install a second
   environment or update toolchains just to reorganize files.

Run these commands from the repository root in its Linux environment:

```sh
python3 research/scripts/research.py list
python3 research/scripts/research.py show RESULT-ID
python3 research/scripts/research.py search QUERY
```

Replace `RESULT-ID` and `QUERY` with the relevant identifier or search words.
The commands read the compact catalogue; follow its source references for the proof.

## Keep sources and evidence together

Each principal result needs a canonical statement with its hypotheses and scope,
a complete proof or explicit open obligations, precise primary references, and the
necessary checkers and inputs. Update the structured catalogue when that mathematical
picture changes. Keep theorem status separate from publication status and from
internal checking, independent replay, formal verification and external review.

For a computation, preserve what another person needs to rerun and assess it:

- The checker source and its exact input or a deterministic input-generation method.
- The command, working directory, software versions and relevant dependency pins.
- Any seed, bounds, arithmetic domain, tolerance, timeout or other parameter that
  limits the result.
- A compact result record giving the actual outcome and what it establishes.
  Include source or input hashes when required by the existing schema.

Use the existing metadata and verification record for the result. Extend it only
when necessary; do not add a second status ledger. A hash identifies bytes and
detects changes; it does not certify an argument. After changing a checker, inspect
the change and rerun the relevant checks before updating its recorded hash. Merely
replacing the hash would conceal stale evidence.

Keep long stdout, checkpoints, caches and exploratory outputs in ignored local
storage such as `.scratch/` or `research/artifacts/local/`. Bound retained logs to
the information needed to diagnose or reproduce an outcome. Promote an artifact
into retained evidence only when a proof, checker or review actually needs it, with
a documented role and provenance. Never discard the only exact input, certificate
or source needed to verify a retained claim.

## Keep changes small and inspectable

Use descriptive filenames and relative links. Follow the repository's UTF-8 and
LF conventions. Avoid unrelated formatting, duplicate utilities, copied dependency
trees, editor files and absolute machine-specific paths.

Before moving or deleting a file, find references in prose, LaTeX includes, imports,
scripts, metadata, manifests and tests. Check untracked files as well as tracked
ones. Update the affected references in the same change. Confirm that a proposed
duplicate is reproducible or superseded before removing it; similar names do not
establish identical content. Preserve exact upstream mirrors identified by source
ledgers, including their checksums and provenance.

Use `git diff --stat`, focused diffs and `git diff --check` to review the change.
Avoid broad cleanup commands such as `git clean` or recursive deletion in a shared
dirty worktree. Stage only the intended paths when staging is part of the task.

LaTeX auxiliaries, compiled caches and environments are ignored. Some paper PDFs
and deposited bibliographies are intentionally retained source artifacts; consult
[the ignore rules](../.gitignore) and manuscript metadata before removing them.
Build drafts with their recorded command, direct disposable output to ignored
storage where supported, and retain the designated review PDF.

## Run checks that cover the change

The current top-level commands are:

```sh
make render-navigation
make check-navigation
make check
```

Run `make render-navigation` after changing mathematical status, then
`make check-navigation`. The renderer updates `research/STATUS.md`, `research/README.md`, and the topic
README pages from `research/MATH_STATUS.json`. Edit that structured source for
result labels and statuses so the summaries stay aligned.
For documentation edits that leave status unchanged, run `make check-navigation`
without regenerating unrelated files.

`make check-navigation` validates catalogue references and checker hashes, manuscript
metadata, exact deposited-file checksums, local Markdown links and tooling controls.
`make check` runs that target and compiles the Python tooling sources. These commands
do **not** replay mathematical computations, build manuscripts, validate every LaTeX
cross-reference, check external links or establish theorem correctness. Their
implementation lives in [the research Makefile](../research/Makefile) and
[the catalogue tool](../research/scripts/research.py).

Select mathematical replay from the result's actual proof and checker instructions.
Read the code before running an expensive or stateful command. Record interrupted
or timed-out runs as such; an incomplete run supplies no completed verification.
For a manuscript change, build the manuscript and inspect the resulting PDF as
described in [the paper guide](PAPERS.md). Use a formal project's recorded toolchain
and build command when formal sources change.

If an existing check fails, determine whether the failure predates or results from
the task. Fix failures caused by the task. Report unrelated failures precisely and
keep them separate from successful checks; do not rewrite unrelated mathematics or
weaken a validator to obtain a pass.

## Prepare a paper bundle

Use the existing manuscript conventions in [unpublished papers](../unpublished/README.md).
A new bundle needs its own README and `metadata.json`, source and designated PDF,
an exact build command, verification scope, and concrete next review tasks.
Record the actual mathematical result, its limits and which checks occurred.

Keep internal planning, repository paths, run labels and remaining project tasks in
the bundle's README or metadata. Apply [the paper guide](PAPERS.md) to the manuscript
itself. Scientific limitations needed to understand the theorem remain in the paper;
they must not be hidden as internal housekeeping.

A submission bundle must reproduce the reviewed PDF from its included sources,
bibliography and figures. Freeze the reviewed bytes and record the edition before
submission. A directory name or passing checker does not mean the paper has been
reviewed, accepted or deposited. Create a `results/` edition only after an actual
public deposit is verified, preserving its exact files and checksums. Submit or
publish only when the user has authorized that action.

## Finish with an evidence based handoff

Review the full diff against the task and report the changed canonical paths, the
resulting behavior or claim, the checks actually run and their outcomes, and any
remaining obligation. Distinguish a new proof from a consistency check and a local
run from independent review. If work is incomplete, identify the next concrete
action in the existing result or manuscript record.

Keep this handoff in the conversation unless durable task information is missing
from the canonical record. Future agents should be able to continue from the
repository's current sources without reconstructing a chain of process diaries.
