# Review and release checklist

Apply this to the exact edition being reviewed. Keep the outcome in the paper's existing verification record; use [the record template](TEMPLATES.md#verification-and-review-record) as needed. Each completed item needs a location, observation, or command result. Mark an inapplicable item with a reason. An unknown or unperformed check remains open.

## Mathematics

- [ ] The main statement specifies the objects, hypotheses, quantifiers, parameter ranges, and conclusion. Its scope matches the title, abstract, introduction, and metadata.
- [ ] Every essential lemma is proved or precisely sourced; the dependency chain is complete and has no circular step or hidden unproved input.
- [ ] The hard step has been checked directly. Exceptional, degenerate, boundary, and empty cases are covered where relevant.
- [ ] Reductions establish the needed implications. Claims of equivalence prove both directions; claims of completeness justify coverage.
- [ ] Finite, generic, local, numerical, and conditional evidence has not been promoted to a stronger statement without a proof.
- [ ] Uniformity, termination, effectiveness, and sharpness claims have the required bounds, arguments, or witnesses.
- [ ] Any discovered error has been propagated through dependent results and summaries.

Use [Mathematical validity](MATHEMATICS.md) to decide what evidence each item needs. An unresolved essential step prevents presenting the unconditional theorem as established. Narrow it, make the dependency explicit, or keep it as a conjecture or research note.

## References and contribution

- [ ] The closest prior work has been read sufficiently to assess the stated advance. Known results, adapted ideas, and new contributions are distinguished.
- [ ] Each substantive citation has verified bibliographic identity and an inspected passage that supports its use under the required assumptions.
- [ ] Preprint, published, corrected, and retracted versions are described accurately. Quotations and borrowed material are attributed.
- [ ] Novelty and priority claims use language justified by the literature check. There are no fabricated references or review outcomes.

## Human readability

- [ ] A reader can find the problem, main result, contribution, and proof idea early without knowing the repository or discovery history.
- [ ] Specialized notation is defined before use and stays consistent. Examples explain the construction or boundary where useful.
- [ ] The proof overview identifies the obstacle and the key step; the proof itself supplies that step.
- [ ] Necessary reasoning has enough detail; repeated motivation, routine expansions, promotional language, and unused terminology have been removed.
- [ ] The manuscript contains no ongoing research reports, repository paths, task IDs, agent instructions, draft comments, or promised future proofs. Required disclosures and substantive mathematical limits remain.
- [ ] Appendices, figures, tables, and supplements are referenced and have a clear purpose.

An internal review can prepare a draft for human review. Record whether an unfamiliar reader actually assessed it and what that reader checked; do not substitute a generated approval for that event.

## Computation and formalization

- [ ] Essential computation has its exact inputs, declared arithmetic and bounds, environment, command, output, and logical role available.
- [ ] Checkers reject malformed or incorrect controls and verify certificates or conclusions rather than merely trusting producer output.
- [ ] Replays identify the person or process, source edition, scope, and outcome. Internal replay, independent implementation, and external review are distinguished.
- [ ] Formal claims name the actual checked theorem, assumptions and axioms, imported inputs, toolchain, and correspondence to the paper. Placeholders and unformalized bridges are disclosed.

## Source bundle and rendered paper

- [ ] The declared bundle builds from the recorded environment without undeclared workspace files. Source and PDF describe the same edition.
- [ ] References and citations resolve. Warnings have been inspected and material layout problems fixed.
- [ ] The PDF has been visually inspected, including equations, page breaks, fonts, tables, figures, and appendices.
- [ ] Abstract, title, authors, identifiers, source manifests, hashes, and verification metadata match the reviewed edition. No placeholder or template field remains.
- [ ] Local links and repository checks have run; failures are reported accurately, including unrelated existing failures. Relevant mathematical checks are reported separately.

## Submission or deposit

- [ ] Required mathematical review and substantive corrections are complete for the intended action. Any remaining obligation is compatible with the exact stated status.
- [ ] The target's current submission format, authorship, attribution, AI disclosure, and artifact requirements have been checked. The disclosure describes the assistance actually used, including proof or code assistance where applicable.
- [ ] The author has authorized the external action. The verified bundle and metadata are ready for that action.
- [ ] Any existing deposited edition remains byte-for-byte intact; the new edition has its own identity and accurate relationship to previous versions.

## Record a specific outcome

Choose the narrowest supported status: research note, draft with open obligations, ready for mathematical review, revised after documented review, or ready for the authorized submission. Record public deposit only after verifying the actual public record and files. Keep publication, external review, formal coverage, and computational checks as separate facts.

For a failure, state the affected claim, evidence, and required correction. For an unperformed check, state it as unperformed. Never turn unchecked boxes into a general assertion that the paper is “fully verified.”
