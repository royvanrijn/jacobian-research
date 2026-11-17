# Working in this repository

This is AI-assisted mathematical research. The front page is a portfolio of
published results. Keep the current mathematical picture clear for human readers.

- results/: deposited papers, exact published files and submission metadata.
- unpublished/: complete manuscript bundles still needing revision, review or publication.
- research/: current principal results, exact proof support and open questions by topic.
- research/MATH_STATUS.json: compact current mathematical status; publication is separate.

Preserve the exact deposited files and their checksums. Keep later revisions in
unpublished/. Do not infer publication, independent replay, formal verification
or external review from an internal checker passing.

Keep the stated theorem scope and its necessary proof sources, primary references,
checker and inputs together. A finite experiment does not prove an unrestricted
statement. Record open obligations plainly; never upgrade assurance during cleanup.

Prefer current result documents over handoffs, process diaries, preupdate snapshots
or duplicated status narratives. Git already preserves the development history.
Do not introduce new research calculations during organization work.

Run make render-navigation after status changes, then make check-navigation.
These commands check current source references, checker hashes, exact publication
files, metadata and local links. They do not replay mathematical computations.
