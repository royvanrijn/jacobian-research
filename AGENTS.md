# Working in this repository

This is AI-assisted mathematical research. The front page is a portfolio of
published results. Keep the current mathematical picture clear for human readers.

Read [the research and paper working guide](working-guide/README.md) before
developing results, changing their evidence, or preparing manuscripts. It covers
repository hygiene, mathematical validity, paper writing and review requirements.

- results/: current editions of released papers; exact deposits and submission metadata in each paper's published/ directory.
- unpublished/: papers without a public deposit, awaiting revision, review or publication.
- research/: current principal results, exact proof support and open questions by topic.
- research/MATH_STATUS.json: current mathematical status and proof dependency DAG; publication is separate.

Preserve the exact deposited files and their checksums under results/<paper>/published/.
Keep the best current source, PDF and review record at results/<paper>/ top level.
Merge a revision from unpublished/ into its released paper rather than maintaining
two current homes; unpublished/ is for papers without a public deposit.
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
