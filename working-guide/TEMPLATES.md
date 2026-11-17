# Reusable records and paper outline

Use these as prompts for reasoning, not boilerplate to paste into papers. Fill the relevant fields in the existing proof, literature, metadata, or verification record. Do not create a second status system or one file per review pass. Preserve existing schemas and adapt the level of detail to the result.

## Claim and dependency record

```text
Claim: descriptive name and manuscript theorem number
Exact statement: objects, hypotheses, quantifiers, conclusion
Status: proved / conditional / conjectural / finite observation / refuted
Closest theorem: exact prior statement and hypotheses
Contribution: exact difference in hypotheses or conclusion
Corollary check: specialization, notation changes, equivalences, known reductions
New argument: precise proof step and location, or accurate alternative contribution
Proof: current source and precise section or lemma
Dependencies: each essential lemma or cited theorem and its hypotheses
Critical step: inference most in need of scrutiny and its justification
Coverage: why all required cases or parameters are included
Computation: role and limits, or none
Formalization: exact theorem and assumptions, or none
Open obligations: specific gaps, or none found in the stated scope
```

Keep the statement precise enough that a later agent can detect a weakened hypothesis or strengthened conclusion. A status label never replaces a proof. If a dependency changes, trace the implications before reusing the conclusion.

## Reference evidence record

```text
Source identity: verified authors, title, date, identifier, version
Accessed source: inspected URL or retained version with provenance
Source type: primary / authoritative secondary; reason and attribution if secondary
Inspected location: theorem, page, section, or passage; justification for essential theorems
Used for: exact claim in our text
Source statement: mathematical content and hypotheses in our own words
Applicability: why those hypotheses hold in our setting
Attribution: theorem, construction, proof idea, quotation, or background
Version issues: corrections or differences that affect this use
Unresolved access or support: precise issue, if any
```

Verification checks identity and support. Attribute secondary accounts accurately without implying the original was read. An essential imported theorem requires inspection of its actual statement and justification; a summary alone leaves that dependency unresolved.

## Verification and review record

```text
Edition: source revision plus exact relevant file hashes
Claim scope reviewed: statements and exclusions
Method: proof reading / exact checker / numerical exploration /
        formal build / editorial review / PDF inspection / blind PDF review
Reviewer or executor: actual person or process and date
Supplied context: exact materials and prior workspace or research exposure
Independence: internal agent, author replay, separate implementation,
              or actual external reviewer; explain shared dependencies
Environment: tool versions and required inputs
Commands: working directory, exact invocation, bounds, random seed if used
Observed result: exit status and concise substantive output
Controls: invalid cases rejected and known cases checked, if relevant
Supports: exact conclusion justified by this check
Does not cover: relevant proof or review obligations outside this method
Blind-reader answers if applicable: novelty, conceptual reason, hardest step and
    proof location, main limitations, verification obstacles
Findings and resolution: issue, affected claims, correction, recheck
Remaining obligations: specific unperformed or failed checks
Outcome: precise draft or review status supported by the evidence
```

For visual inspection, record the actual PDF edition and findings. For blind review, document PDF-only context isolation; shared workspace context disqualifies the test and isolated AI review remains internal. For a proof review, state which dependencies were read. A proposed command has no observed result. Record failures and unperformed checks as well as successes.

## Paper outline

Adapt this order when the mathematics calls for another structure. A short note may combine sections. A computational proof needs enough method and certification detail to be checked.

1. **Title and abstract.** Identify the problem, exact result, decisive restrictions, and mathematical consequence. Keep the title consistent with the proved scope.
2. **Introduction.** Explain the question, prior work, exact contribution, and main theorem. Before machinery, give plausibility, obstacle, insight, and a simple nontrivial example where possible. Explain failed natural approaches when useful; avoid discovery chronology.
3. **Definitions and setup.** Give the objects, conventions, hypotheses, and external results needed for the proof. Verify the applicability of those results.
4. **Key argument.** Present the lemmas in a useful logical order. Explain the difficult step and prove the required coverage.
5. **Main theorem and consequences.** Assemble the argument, identify the scope, and show consequences and sharpness examples that are actually established.
6. **Further discussion if useful.** Explain a mathematical limitation or a precise open question. Omit progress reports and promises about ongoing work.
7. **Appendices or supplement if needed.** Supply lengthy supporting arguments, certificates, or implementation detail. Keep the main mechanism visible in the article.
8. **References, acknowledgments, and required disclosures.** Use verified entries and truthful statements of assistance, support, and review. Follow the target's required placement.

## Standard AI attribution

Use the [standard AI attribution](PAPERS.md#ai-assistance-and-review) verbatim in each manuscript. Keep mathematical provenance, acknowledgments, and the actual scope of verification separate from that paragraph. Do not add draft or review-status notes.

Follow the target's required disclosure placement. Keep all statements factual; never invent a reviewer, an affiliation, an author contribution, or an endorsement, or claim that a human audited every step unless that happened.

## Completion note for the user

```text
Changed: current files and purpose
Mathematical effect: changed scope or none
Verified: checks actually run and the conclusions they support
Open: relevant failures, unperformed checks, and remaining obligations
Next action: concrete next step if the requested work is incomplete
```

Report the current result concisely. Keep internal paths here when needed for action; keep them out of manuscript prose.
