# Working in this repository

This is AI-assisted mathematical research. The front page is a portfolio of
released and submitted papers. Keep the current mathematical picture clear for human readers.

Read [the research and paper working guide](working-guide/README.md) before
developing results, changing their evidence, or preparing manuscripts. It covers
repository hygiene, mathematical validity, paper writing and review requirements.

Keep a single source of truth for this project. Do all new research in the
canonical Ubuntu-24.04 WSL checkout at `/home/royvanrijn/src/jacobian-research`
or a linked Git worktree that shares its Git directory. Record the results and
their supporting evidence in this repository. Use `royvanrijn/jacobian-research`
as the project's sole GitHub repository. Follow the
[canonical workspace instructions](working-guide/REPOSITORY.md#canonical-workspace)
to verify the checkout before starting.

Before starting a new mathematical approach, whenever possible check the primary
literature for an exact solution or an argument that can be adapted. Inspect the
actual theorem statements and hypotheses, identify the closest applicable result
and remaining gap, and use that review to guide new proofs and computations.
Do not infer novelty or an unsolved status merely from an unsuccessful search.

- results/: current editions of released or submitted papers. Record pending submissions in paper metadata.
- unpublished/: manuscripts in preparation.
- research/: current principal results, exact proof support and open questions by topic.
- research/MATH_STATUS.json: current mathematical status and proof dependency DAG; publication is separate.

Start a DAG assignment with `python3 research/scripts/research.py gaps --unresolved`, then use
`explain ID` to inspect the exact scope, proof routes and blockers. The `gaps --all` inventory
also includes resolved and falsified endpoints and endpoints outside topic summaries.

Preserve the exact deposited files and their checksums in published/ inside each paper directory.
Keep the best current source, PDF and review record at the top level of its paper directory.
Merge a revision from unpublished/ into its released paper rather than maintaining
two current homes.
Do not infer publication, independent replay, formal verification or external
review from an internal checker passing.

Keep the stated theorem scope and its necessary proof sources, primary references,
checker and inputs together. A finite experiment does not prove an unrestricted
statement. Record open obligations plainly; never upgrade assurance during cleanup.

Prefer current result documents over handoffs, process diaries, preupdate snapshots
or duplicated status narratives. Git already preserves the development history.
Do not introduce new research calculations during organization work.

Run make render-navigation after status changes, then make check-navigation.
These commands check current source references, checker hashes, acyclic implications,
graph evidence freshness, exact publication
files, metadata and local links. They do not replay mathematical computations.
