# Research memory

Use the [algorithmic lessons](knowledge/ALGORITHMS.md) and
[scoped failed routes](knowledge/FAILED_ROUTES.md) before planning a calculation.
Both link back to evidence. Mathematical status comes only from
[MATH_STATUS.json](MATH_STATUS.json).

## Find the prior result

Run from either checkout directory using the appropriate script path:

```sh
# From the repository root:
python3 research/scripts/research.py search "pointed quartic"
python3 research/scripts/research.py search "q8" --history
python3 research/scripts/research.py search "annihilator"
python3 research/scripts/research.py show EC-CURVE302-RECOVERED-MW17-PARENT
python3 research/scripts/research.py routes --area elkies-k3
python3 research/scripts/research.py work --area elliptic-curves
python3 research/scripts/research.py resources --kind mechanism
```

Search reads current scopes, method records, proposed work, maintained JSON
records and source text, including canonical TeX/code and historical EC/K3 Markdown.
`--history` adds other archived programmes, cleanup snapshots and inherited work items. Results
label their authority. `show` returns the complete selected claim, work item or
resource, with its evidence and current references.
Use `--kind historical` or `--kind lesson` to focus a search. Current claims and
lessons are ranked ahead of archived prose.

The [catalogue](index/README.md) covers all registered claims. The
[document inventory](index/documents.tsv) also includes source notes with no
registered claim; these must not be silently promoted.

## Unknowns and suggested work

The [work ledger](knowledge/WORK_LEDGER.md) covers every active open problem,
partial result and parked problem. Proposed actions state a next step,
completion evidence and prerequisites. Their claim fingerprints require review
when a supporting scope changes; rendering cannot approve that review.
All actions remain unscheduled until work is explicitly taken up.

The [source-level partial review](knowledge/PARTIAL_REVIEW.md) distinguishes
actually reviewed obligations from merely indexed claims. Its source and code
hashes make later drift visible; `research.py check --require-partial-review`
fails until every current partial result has a fresh review. The navigation
check now requires that coverage. This reviews the surviving obligations;
it does not independently reprove every theorem in the registry.
Maintained proof sources remain mandatory. Hash receipts for ignored local
calculation files are checked when those files are present; `show ID` reports
missing local evidence separately. Such receipts preserve an inspection, not
a promise that a fresh checkout contains every replay input.

The [historical checklist reconciliation](knowledge/LEGACY_WORK_REVIEW.md)
retains every unchecked retrospective item, including completed subsets and old
blocked labels. Use `research.py show LEGACY-20260904-556` to see why completed
parent recovery and unknown original provenance now have different destinations.

[Structured retrieval](index/resources.md) reads the K3 mechanism/process
ledger and curve database directly; archived support gates require `--history`. For example,
`research.py resources --kind curve --unknown-only` lists missing exact
conductor certificates without starting factorization.

## What to retain from every experiment

Add or update one record in [knowledge/lessons.json](knowledge/lessons.json) when
work changes how a later calculation should be done. Record:

- the question and conditions where the method applies;
- the reusable implementation and exact source;
- the failed approach and what it actually excludes;
- the certificate boundary and what would justify revisiting it.

Supporting claims have reviewed state/title/scope fingerprints. If one changes,
search marks the lesson for review and `check-navigation` fails until the lesson
is reviewed against the new scope. `research.py show LESSON-ID` prints both the
current fingerprints and the full supporting claims; do not refresh fingerprints
without reviewing applicability. Rendering never approves a changed lesson.

Keep empirical timings tied to the measured component, input and environment.
A faster map worker does not imply the same speedup for an entire search.
Keep bounded misses, timeouts and missing prerequisites distinct.

## Durable preflight

Reuse exact model/field contexts, certified point signatures, sealed landscapes
and completed chart receipts before reconstructing them. Deduplicate the actual
fibration, subgroup and bounded box at the equivalence relation being used.
Check theorem gates before expensive equation work; heuristic scores only
schedule explicitly bounded work.

After a result changes, update its canonical note and typed claim, repair its
consumers, add the method lesson, then run:

```sh
make render-navigation
make check-navigation
```

The lightweight navigation workflow runs these checks on relevant pull requests;
it requires no CAS installation or research calculations.

The [earlier knowledge-base snapshot](archive/repository-cleanup-2026-09-12/research__KNOWLEDGE_BASE.md.txt)
is preserved as history. It contains dated counts and pre-recovery frontiers.

<!-- status-consumer: EC-K3-R17-KUMMER-CLASSGROUP-PRESSURE-COMPARISON 74b1dae24470b531 -->
