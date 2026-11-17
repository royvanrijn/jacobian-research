# Research and paper working guide

Use this guide whenever you develop mathematics, change its supporting files, or prepare a paper. The goal is a repository whose current claims can be checked and papers that a human reviewer can understand without knowing how the work was produced.

Start with the repository [agent instructions](../AGENTS.md). This guide supplies the detailed workflow; it does not change the status of any existing theorem, paper, or verification record.

## Essential rules

1. State exactly what is proved, under which hypotheses, and where the proof is. Keep conjectures, conditional results, experiments, and proved statements distinguishable.
2. Check the difficult mathematical steps and every dependency. Fluent prose, agreement between agents, a successful build, and a passing checker each provide limited evidence.
3. Verify citation identity and support. Attribute authoritative secondary sources accurately; inspect the actual statement and justification of essential imported theorems.
4. Write around the contribution and argument. Keep research updates, discovery chronology, repository paths, task names, and agent instructions out of manuscript prose. Explain failed natural approaches when they clarify an obstruction; state mathematical limitations directly.
5. Explain the central idea before technical machinery and give enough detail to check the hard step. Every definition, lemma, and example must earn its place; remove repetition and routine detail that conceal the argument.
6. Maintain one current account of a result and its evidence. Preserve necessary proof sources and exact replay inputs; use Git for superseded development history.
7. Preserve deposited bytes and checksums in each released paper's `published/` archive. Keep its current revision at the paper directory's top level, and describe publication, mathematical review, computation, and formalization with their actual scopes.
8. End with evidence: what changed, which checks ran, their outcomes, and any remaining obligations. Never claim a check, independent review, or publication that did not happen.

## Find the relevant instructions

| Guide | Use it for |
|---|---|
| [Repository workflow](REPOSITORY.md) | Choosing locations, preserving evidence, avoiding clutter, changing status, and running the existing tools. |
| [Mathematical validity](MATHEMATICS.md) | Statement design, proof dependencies, counterexamples, computations, and formal verification. |
| [Writing papers](PAPERS.md) | A coherent standalone article, reliable citations, readable proofs, and avoiding documented failures in AI-assisted papers. |
| [Review and release checklist](CHECKLIST.md) | Deciding whether a draft can move to review or an authorized submission. |
| [Reusable records and outline](TEMPLATES.md) | Compact claim, reference, and verification records and a paper outline. Adapt these inside existing records. |
| [Sources and rationale](SOURCES.md) | Primary research and expert guidance behind the editorial rules, with their limits. |

## Work in stages

**Orient.** Read the current result, metadata, proof sources, and applicable instructions. Inspect the worktree and identify existing changes before editing. Decide the actual deliverable and its audience.

**Establish the mathematics.** Write the proposed statement and its dependencies before writing an abstract that announces it. Test vulnerable cases. Resolve gaps or narrow the statement. Record an unresolved implication as an obligation, not as a completed theorem.

**Make the evidence reproducible.** Keep exact inputs and a small checker where computation matters. Explain what the computation certifies and how it connects to the proof. Pin the relevant environment and distinguish internal execution from independent replay.

**Build the argument for a reader.** Choose the contribution and compare exact statements with the closest theorem, testing specializations and known reductions before claiming novelty. Explain the central idea and a representative example before machinery. Write and revise using [the paper guide](PAPERS.md).

**Review in separate passes.** Audit mathematical dependencies, source attribution, exposition, computation where applicable, and the rendered document. Apply [the PDF-only blind-reader procedure](PAPERS.md#review-the-actual-article), recording supplied context, findings, and revisions. An agent with workspace context does not qualify; isolated AI review remains internal. Record external human review only when it occurs.

**Freeze the appropriate edition.** Use [the checklist](CHECKLIST.md), build the exact sources, inspect the PDF, and align the abstract and metadata with the final theorem. A draft may be ready for mathematical review while still awaiting that review. Submission and public deposit follow the user's authorization and the target's current requirements.

**Maintain the result.** When a proof changes, revisit dependent statements, summaries, abstracts, metadata, checkers, and formal claims. Preserve the identity of earlier deposited editions and describe corrections explicitly.

## How to apply the rules

The repository and mathematical requirements are working rules. Editorial choices such as section lengths, notation tables, and appendices depend on the result and audience. Use them when they reduce the reader's effort; do not pad a short paper to fit a template.

The research in [Sources and rationale](SOURCES.md) documents failure modes and established exposition advice. It is not a survey of what most researchers believe. Judge drafts by their mathematics, evidence, attribution, and readability. Style alone does not establish AI authorship or mathematical validity.
