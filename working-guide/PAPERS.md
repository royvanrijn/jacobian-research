# Writing readable mathematical papers

A reviewer should be able to identify the question, the exact contribution, its relationship to prior work, and the mechanism of the proof. Every section should help the reader understand or verify one of those things. Follow [Mathematical validity](MATHEMATICS.md) while composing; editing cannot repair an unsupported theorem.

## Plan a contribution that can be assessed

Before drafting, write a short internal answer to these questions: What is the precise problem? What was already known? What does this paper add? What idea makes the new step possible? What hypotheses or limitations matter? Locate the proof or other evidence for each answer. If the answer is only a reformulation of a known result, identify it as exposition or a new proof and explain its value.

Choose the audience and assume a reasonable shared background. Define specialized objects needed beyond that background. A paper need not reprove standard mathematics, but it must state enough for a reader to check applications of it. An internal note is not a substitute for a missing proof.

Build a logical outline with the central theorem, its inputs, and its proof. Order the exposition for understanding. The order in which experiments, failures, or conversations happened rarely supplies that order. Established [exposition guidance](SOURCES.md#exposition-guidance) supports early examples, a clear contribution, and a direct route to the idea.

## Keep the manuscript independent of the workspace

The article describes completed mathematics. It must make sense to a reviewer who has only the submitted PDF, its bibliography, and any declared supplement.

Exclude ongoing research updates, plans for future agents, task identifiers, internal result nicknames, branch names, local paths, command transcripts, audit receipts, changing repository status, and accounts of what agents tried. Do not refer to “our current pipeline,” “the next research sprint,” or a lemma available only in a research folder. Replace an internal label with an ordinary mathematical definition or a descriptive lemma name.

Bring every new essential argument into the article or a clearly identified submitted appendix. Cite public mathematical literature precisely. A stable, versioned citation to code or data may belong in an artifact or availability statement when needed for reproducibility. Describe the algorithm, input class, and proof role in the article; keep operational paths and commands in the supplement's README. An artifact link must not conceal the only explanation of why the method proves the theorem.

Keep actual mathematical qualifications. If an argument assumes an unproved statement, put that hypothesis in the theorem. If the result covers a restricted family, name that family. State an unresolved mathematical boundary when it explains the scope, without narrating the authors' ongoing work. A brief open question can be appropriate; a progress report or promise of a forthcoming proof is not evidence.

Required acknowledgments, AI-use disclosures, and statements about verification belong in their appropriate places and must remain accurate. Removing workspace narration is not permission to hide a proof dependency or material assistance.

## Title and abstract

Use a descriptive title that identifies the object and result. Avoid promotional claims, coined labels without mathematical value, and announcements that outrun the theorem. If only a bounded or conditional result is proved, the title must not imply a complete resolution.

The abstract should give the problem, the exact main result with its decisive restrictions, the main mechanism if useful, and the mathematical consequence. Use as much space as needed to express these clearly within the target's limit. Minimize undefined notation, abbreviations, bibliography references, and background that delays the result.

Check each sentence against the final theorem. “All,” “sharp,” “optimal,” “first,” “complete,” “effective,” and “classifies” each create obligations. An optimal bound needs a stated sense of optimality and a matching example or lower bound. An existence proof does not supply an effective algorithm without a justified procedure. A search that finds no earlier paper does not establish a priority claim.

Finalize the abstract after the argument stabilizes. Keep metadata and submission abstracts identical in mathematical meaning; any difference in notation or length must preserve hypotheses and qualifications.

## Introduction and relationship to prior work

Start with the mathematical question and why the distinction being studied matters. Give a small example when it makes the question concrete. Explain the closest relevant work before claiming an advance over it. Attribute an existing theorem, proof idea, or construction at the point where it enters the story.

State the principal result early, with enough definitions to understand its content. A more technical version may appear later, but the early version must preserve scope. Explain the proof strategy in ordinary mathematical language: the reduction, the obstacle it leaves, the lemma that removes it, and the final deduction. A list of section numbers does not provide this explanation.

Compare hypotheses and conclusions explicitly. Explain what the new result adds and what existing work already provides. Avoid a long catalogue of weakly related citations and comparisons with obsolete or artificially restricted versions of the closest result. A technical comparison can wait until its notation is introduced.

## Definitions, notation, and proof flow

Introduce notation before use and keep one meaning per symbol in the same context. Specify the field, coefficient ring, characteristic, topology, dimension, and conventions when relevant. Distinguish objects from their coordinates and representatives. Keep notation short enough to read; avoid several indices that exist only to encode the development history.

Use paragraphs to explain mathematical relationships. Display an equation when readers need to inspect or refer to it, and integrate it into a grammatical sentence. State what an equation establishes before moving to the next step. Define constants and explain their dependence, especially in asymptotic statements and uniform estimates.

Begin a long proof with a short plan. Give the hard idea room. State useful intermediate claims as lemmas with meaningful hypotheses, and explain how they assemble into the main theorem. When splitting into cases, say why the cases exhaust the possibilities and where degenerate cases are treated. End the proof by connecting the established facts to the stated conclusion.

Use a citation only after checking its hypotheses in the present setting. Replace “standard arguments,” “it is clear,” and “a routine calculation” with the needed reasoning whenever they conceal an essential inference. Conversely, compress routine algebra once the mechanism is visible. Excessive explanation of elementary steps can make an unproved difficult step harder to notice; [Tao's discussion](SOURCES.md#mathematical-practice) describes this concern in AI-assisted mathematics.

Keep main proof ideas in the main text. Put lengthy routine expansions, exhaustive tables, and implementation detail in an appendix or supplement with an explicit reference and a description of their role. Appendices remain part of the mathematical review; moving a gap there does not resolve it.

## Examples, figures, and computations

Use examples for a purpose: introduce an object, demonstrate the mechanism, show a boundary case, or prove sharpness. Work at least one representative case through enough detail that the reader can follow the general construction. Distinguish an illustration from evidence for a universal claim.

Every figure or table must answer a mathematical question raised in the text. Define its symbols, explain the conditions and units, and state what the reader should notice. Use readable labels and captions that explain the content. Refer to each figure from the argument. Avoid decorative diagrams, screenshots of code, giant unexplained matrices, and tables that repeat the prose.

For computational results state the searched class, bounds, arithmetic, completeness argument, and verification method. If a computation is essential to a proof, explain its logical role and provide the certificate or reproducible procedure. Report a numerical observation as numerical evidence unless a rigorous error bound certifies the required conclusion. See [the computation rules](MATHEMATICS.md).

## Verify citations and attribution

For every substantive citation, open the primary source, identify the relevant theorem or passage, and record its assumptions and exact support in the existing literature record. Confirm title, authors, year, version, and DOI or other stable identifier against the actual work. A resolving DOI can still point to the wrong article. An existing article can still fail to support the cited claim. [An empirical citation audit](SOURCES.md#citation-reliability) motivates checking both.

Use search results and generated bibliographies to locate candidates, then verify them. Do not cite a source as checked when only an abstract or another author's summary was available. Seek an accessible primary version; if the exact needed result remains inaccessible, record that internally and resolve the dependency before presenting it as established. A public preprint is a legitimate source when identified accurately; do not call it peer reviewed without evidence.

Check corrected, superseded, or retracted versions when relying on a result. Cite the version actually used and explain material differences. Give credit for borrowed ideas and text, including translations and AI-mediated paraphrases. Avoid copying a related-work paragraph and changing its wording without inspecting its sources. Never invent authors, acknowledgments, quotations, review outcomes, or priority.

## Recognize recurring failures in AI-assisted papers

These are review targets, not tests of who wrote a paper. The [source notes](SOURCES.md) distinguish empirical findings from expert advice and position papers.

| Failure a reviewer encounters | Required revision |
|---|---|
| Fluent paragraphs whose conclusions do not follow from earlier sections | Trace each main claim to its proof and check that the same hypotheses and objects survive each transition. |
| Repeated motivation, duplicated theorem summaries, or long generic conclusions | Keep the needed explanation once, with a short reminder only where it helps a reader resume the argument. |
| Impressive terminology or a named framework without a mathematical purpose | Define the object, explain its use, and remove unnecessary terminology. |
| Plausible citations with no verifiable or relevant source | Verify the bibliographic identity and the exact supporting passage; repair or remove the unsupported claim. |
| A polished introduction promising more than the proofs establish | Narrow the title, abstract, introduction, theorem, and metadata together. |
| Easy calculations explained at length while the central step is asserted | Compress routine work and supply the missing central argument. |
| Topic lists presented as a literature synthesis | Explain relationships, distinctions, and the specific unresolved question; cite selectively. |
| Agent commentary, template residue, TODOs, or inconsistent notation | Rewrite from the intended reader's perspective and inspect the complete compiled article. |

Do not edit merely to defeat an AI detector or remove a list of supposedly suspicious words. Prefer concrete claims and precise verbs; remove inflated praise such as “groundbreaking” and empty transitions. A restrained voice is useful because it makes the actual claim easier to assess.

## Examples of revisions

These examples are illustrative, not claims about results in this repository. Substitute a real theorem only after checking its evidence.

| Draft wording | Revision and reason |
|---|---|
| “Our ongoing research pipeline has now closed the degree gate.” | “Theorem 2 excludes maps of degree at most four satisfying hypotheses (i)–(iii).” This names the mathematical scope; the proof must establish it. |
| “The proof is in the repository under the latest closure note.” | Include the proof as Lemma 3 or cite a stable public result with an exact locator. A reviewer needs an inspectable argument. |
| “Extensive tests prove the identity for all parameters.” | “Exact evaluation verifies the listed parameter values.” Claim the general identity only with a separate symbolic or other complete proof. |
| “The remaining cases will be resolved in our next paper.” | “The argument requires the nondegeneracy hypothesis.” State the boundary that affects the present result. |
| “We introduce a powerful unified framework.” | State the new construction and the theorem it enables. Remove the praise. |

## Review the actual article

Do separate readings for logic, citations, exposition, and layout. During the exposition pass, read the title, abstract, introduction, main statement, proof overview, and conclusion together. They should tell the same mathematical story. Then follow the proof linearly and note every backward search for a definition or hidden assumption.

Build from the declared source bundle with the declared toolchain. Resolve missing citations and cross-references, inspect warnings, and visually inspect the resulting PDF. Check equations, fonts, page breaks, figures, table widths, links, and appendices. After substantial edits, repeat the affected checks. A successful compiler exit does not establish readable pages or a valid argument.

Use [the review checklist](CHECKLIST.md) to record the outcome. Internal agents may help locate errors and confusing passages. Human authors remain responsible for understanding the proof and accurately describing the review and assistance that actually occurred.
