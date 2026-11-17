# Writing readable mathematical papers

A paper should teach a mathematician something they can understand, verify, and build upon. Make the question, exact contribution, relationship to prior work, and proof mechanism clear. Follow [Mathematical validity](MATHEMATICS.md); editing cannot repair an unsupported theorem.

## The human-reader test

Before considering a manuscript finished, verify that:

- A mathematician can understand the contribution from the introduction alone.
- The central insight is explained without requiring every technical detail.
- The hardest argument is visible and verifiable.
- Every mathematical object, lemma, and example serves a clear purpose.
- A reader unfamiliar with the research workspace can follow the paper without additional explanations.

This tests communication alongside correctness: does reading the paper convey mathematics the reader can use?

## Plan a contribution that can be assessed

Before drafting, record the precise problem and compare with the closest known theorem in the existing literature or claim record:

| Question | Required answer |
|---|---|
| What did the previous result establish? | Its exact statement and hypotheses. |
| What does ours establish? | The exact mathematical difference in hypotheses or conclusion. |
| Could ours be an immediate corollary? | Check specialization, changes of notation, equivalences, and known reductions. |
| What genuinely new argument is needed? | Identify the precise proof step and its location. |

Challenge the claimed advance using the strongest relevant version of prior work. If it reduces to a known theorem, describe its value accurately as a consequence, exposition, or a new proof. Do not build a novelty claim around a change of formulation.

Choose the audience and its shared background. State enough for readers to check applications of standard mathematics without reproving it. Outline the theorem, its inputs, and its proof in the order that aids understanding. [Exposition guidance](SOURCES.md#exposition-guidance) supports early examples and a direct route to the idea.

## Keep the manuscript independent of the workspace

Describe established results accurately, including partial results, conditional theorems, and rigorous counterexamples. The article must make sense with only the submitted PDF, bibliography, and any declared supplement.

Exclude research updates, agent plans, task identifiers, internal nicknames, branch names, local paths, command transcripts, audit receipts, and changing repository status. Replace internal labels with ordinary mathematical definitions or descriptive lemma names. A mathematical explanation of why a natural approach fails can clarify an obstruction; present that explanation without recounting the discovery chronology or agent activity.

Include every new essential argument in the article or an identified submitted appendix. Cite public literature precisely. A stable, versioned code or data citation may belong in an availability statement. Explain the algorithm, input class, and proof role in the article; keep operational paths and commands in the supplement's README. An artifact link cannot replace the explanation of why a method proves the theorem.

Keep mathematical qualifications explicit. Put unproved assumptions in the theorem and name any restricted family. State unresolved boundaries that explain scope. An open question can be appropriate; a promise of a forthcoming proof is not evidence. Required acknowledgments, AI-use disclosures, and verification statements must accurately describe dependencies and material assistance.

## Title and abstract

Use a descriptive title identifying the object and result. Avoid promotional claims and coined labels without mathematical value. A bounded or conditional result must not appear to resolve the unrestricted problem.

Give the problem, exact main result with decisive restrictions, mechanism if useful, and mathematical consequence within the target's limit. Minimize undefined notation, abbreviations, citations, and background that delays the result.

Check each sentence against the final theorem. “All,” “sharp,” “optimal,” “first,” “complete,” “effective,” and “classifies” each create obligations. Optimality needs a stated sense and a matching witness or bound; effectiveness needs a justified procedure. Finding no earlier paper does not establish priority.

Finalize the abstract after the argument stabilizes. Metadata and submission abstracts must preserve the same mathematical meaning, hypotheses, and qualifications.

## Introduction and relationship to prior work

Start with the question and why it matters. State the principal result early, with enough definitions to understand its content. A technical version may appear later, but the early version must preserve scope. Explain the closest relevant work and the exact advance; attribute borrowed theorems, ideas, and constructions where they enter.

Before substantial technical machinery, explain the central idea in a few paragraphs of ordinary mathematical language. Why is the result plausible? What is the obstacle? What insight overcomes it? Identify the key step and final deduction so readers can recognize where the idea might be reused. Where possible, work through the simplest nontrivial example. A section roadmap or list of lemmas does not supply this explanation.

Compare hypotheses and conclusions explicitly, using the internal novelty comparison. Avoid catalogues of weakly related citations or comparisons with obsolete or artificially restricted versions. Technical comparisons may wait for their notation.

## Definitions, notation, and proof flow

Introduce notation before use and keep one meaning per symbol in context. Specify relevant fields, coefficient rings, characteristic, topology, dimension, and conventions. Distinguish objects from their coordinates and representatives. Avoid indices that merely encode development history.

Explain mathematical relationships in prose. Display equations readers need to inspect or reference, integrate them into sentences, and say what they establish. Define constants and their dependence, especially for uniform or asymptotic claims.

Begin a long proof with a short plan and give the hard idea room. State useful intermediate claims with meaningful hypotheses and explain how they assemble. Explain why cases are exhaustive and cover degeneracies. End by connecting the established facts to the conclusion.

Distinguish intuition, proof sketches, and proofs explicitly. An informal explanation can convey the mechanism but cannot silently replace an essential argument. Readers must be able to tell what is established and what is illustrative; label sketches and locate the complete proof of every essential claim.

Check cited hypotheses in the present setting. Replace “standard arguments,” “it is clear,” or “a routine calculation” when they conceal an essential inference. Compress routine algebra once its mechanism is visible. [Tao's discussion](SOURCES.md#mathematical-practice) describes how excessive routine detail can obscure the difficult step.

Keep main ideas in the main text. Move lengthy routine expansions, exhaustive tables, and implementation detail to a referenced appendix or supplement, explaining their role. Appendices remain part of mathematical review; moving a gap does not resolve it.

## Check mathematical density

Every definition must be used and earn its complexity. Every lemma must carry a useful step; combine trivial lemmas when that improves the narrative. Prefer familiar terminology to names for one-off constructions. Remove duplicated motivation and theorem summaries, retaining reminders only where readers need them.

Never inflate a modest result to make it appear substantial. Do not compress a difficult inference into an unexplained sentence either. Aim for maximum mathematical insight per unit of reader effort. Audit each object, lemma, and example for its purpose, and each hard step for sufficient explanation. Length follows the argument's needs.

## Examples, figures, and computations

Use examples to introduce an object, demonstrate the mechanism, expose a boundary, or prove sharpness. Work a representative case through enough detail to explain the construction. An illustration does not prove a universal claim.

Each figure or table must answer a mathematical question in the text. Define symbols, conditions, and units; use readable labels and captions; explain what readers should notice. Reference it from the argument. Remove decoration, code screenshots, unexplained matrices, and tables that repeat prose.

For computations state the searched class, bounds, arithmetic, completeness argument, and verification method. Explain any essential computation's logical role and provide a certificate or reproducible procedure. Numerical observations remain evidence unless rigorous error bounds certify the conclusion. Follow [the computation rules](MATHEMATICS.md).

## Verify citations and attribution

For substantive citations, prefer an inspected primary source. Record the supporting theorem or passage and assumptions in the existing literature record. Confirm authors, title, year, version, and stable identifier against the work: a real article or resolving DOI can still fail to support the claim. [An empirical audit](SOURCES.md#citation-reliability) motivates checking identity and support separately.

For historical claims or inaccessible originals, an authoritative secondary source can be legitimate when accurately attributed. State which source was actually inspected; do not imply the original was read. An essential imported theorem still requires inspection of its actual mathematical statement and justification in a source that supplies them. A summary or abstract cannot settle that dependency. Seek an accessible version, supply the proof, or leave the dependency unresolved.

Use search results and generated bibliographies to locate candidates, then verify them. Identify preprints accurately without implying peer review. Check corrections, superseded versions, or retractions and explain material differences in the version used. Attribute borrowed ideas and text, including translations and AI-mediated paraphrases. Never invent authors, quotations, acknowledgments, review outcomes, or priority.

## Recognize recurring failures in AI-assisted papers

These are review targets, not authorship tests. The [source notes](SOURCES.md) distinguish empirical findings, expert advice, and position papers.

| Failure a reviewer encounters | Required revision |
|---|---|
| Fluent transitions whose conclusions do not follow | Trace claims to proofs and check hypotheses and objects at each transition. |
| Elaborate setup without a visible central idea | Explain plausibility, obstacle, and insight before machinery; remove setup that serves no purpose. |
| Plausible citations without relevant support | Inspect the supporting passage; repair or remove unsupported claims. |
| An introduction promising more than the proof establishes | Narrow the title, abstract, introduction, theorem, and metadata together. |
| Routine calculations expanded while the central step is asserted | Compress routine work and supply the essential argument. |
| Workspace narration, template residue, TODOs, or inconsistent notation | Revise from the reader's perspective and inspect the complete compiled article. |

Do not optimize for an AI detector or remove supposedly suspicious words mechanically. Use concrete claims and precise verbs; remove inflated praise and empty transitions so the mathematics is easier to assess.

## Review the actual article

Review logic, citations, exposition, computation where applicable, and layout in separate passes. Read the title, abstract, introduction, theorem, proof overview, and conclusion together for consistency. Follow the proof linearly and note backward searches for definitions or hidden assumptions.

Build the declared source bundle with its recorded toolchain. Resolve missing citations and cross-references, inspect warnings, and visually check the PDF: equations, fonts, page breaks, figures, tables, links, and appendices. Repeat affected checks after substantial edits. Compilation alone establishes neither readability nor correctness.

Run a blind-reader test on the compiled PDF before calling the manuscript finished. Give a reviewer only that PDF, without research history, repository access, this conversation, or accompanying agent explanations. Ask them to identify:

- What exactly is new?
- Why is it true at a conceptual level?
- Which step is hardest, and where is it proved?
- What are the main limitations?
- What prevents them from verifying the result?

Record the PDF edition, reviewer, supplied context, answers, and obstacles in the existing verification record. Revise the manuscript to address those obstacles and recheck affected passages. This tests what the PDF communicates; review of declared supplements and proof dependencies remains necessary separately. If the test has not occurred, record it as unperformed.

An AI with workspace or research context does not qualify as a blind reader. A fresh AI review counts only when context is actually isolated to the PDF and must be recorded as internal review, separately from external human review. Use [the checklist](CHECKLIST.md) to record outcomes. Human authors remain responsible for understanding the proof and describing the assistance and review that occurred.
