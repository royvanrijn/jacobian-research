# Reusable records and paper outline

Use these as prompts for reasoning, not boilerplate to paste into papers. Fill the relevant fields in the existing proof, literature, metadata, or verification record. Do not create a second status system or one file per review pass. Preserve existing schemas and adapt the level of detail to the result.

## Claim and dependency record

```text
Claim: descriptive name and manuscript theorem number
Exact statement: objects, hypotheses, quantifiers, conclusion
Status: proved / conditional / conjectural / finite observation / refuted
Contribution: what is new; closest prior result and exact difference
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
Accessed source: primary URL or retained version with provenance
Inspected location: theorem, page, section, or passage
Used for: exact claim in our text
Source statement: mathematical content and hypotheses in our own words
Applicability: why those hypotheses hold in our setting
Attribution: theorem, construction, proof idea, quotation, or background
Version issues: corrections or differences that affect this use
Unresolved access or support: precise issue, if any
```

Verification has two parts: locating the real work and checking what it establishes. If only a summary was inspected, record that limitation and resolve the substantive dependency before treating it as checked.

## Verification and review record

```text
Edition: source revision plus exact relevant file hashes
Claim scope reviewed: statements and exclusions
Method: proof reading / exact checker / numerical exploration /
        formal build / editorial review / PDF inspection
Reviewer or executor: actual person or process and date
Independence: internal agent, author replay, separate implementation,
              or actual external reviewer; explain shared dependencies
Environment: tool versions and required inputs
Commands: working directory, exact invocation, bounds, random seed if used
Observed result: exit status and concise substantive output
Controls: invalid cases rejected and known cases checked, if relevant
Supports: exact conclusion justified by this check
Does not cover: relevant proof or review obligations outside this method
Findings and resolution: issue, affected claims, correction, recheck
Remaining obligations: specific unperformed or failed checks
Outcome: precise draft or review status supported by the evidence
```

For visual inspection, record the actual PDF edition and findings. For a proof review, state which dependencies were read. A command that was proposed but not executed has no observed result. Record reproducible failures as well as successful controls.

## Paper outline

Adapt this order when the mathematics calls for another structure. A short note may combine sections. A computational proof needs enough method and certification detail to be checked.

1. **Title and abstract.** Identify the problem, exact result, decisive restrictions, and mathematical consequence. Keep the title consistent with the proved scope.
2. **Introduction.** Explain the question, a motivating example, relevant prior work, the contribution, the main theorem, and the proof idea. Avoid development history.
3. **Definitions and setup.** Give the objects, conventions, hypotheses, and external results needed for the proof. Verify the applicability of those results.
4. **Key argument.** Present the lemmas in a useful logical order. Explain the difficult step and prove the required coverage.
5. **Main theorem and consequences.** Assemble the argument, identify the scope, and show consequences and sharpness examples that are actually established.
6. **Further discussion if useful.** Explain a mathematical limitation or a precise open question. Omit progress reports and promises about ongoing work.
7. **Appendices or supplement if needed.** Supply lengthy supporting arguments, certificates, or implementation detail. Keep the main mechanism visible in the article.
8. **References, acknowledgments, and required disclosures.** Use verified entries and truthful statements of assistance, support, and review. Follow the target's required placement.

## Disclosure drafting prompts

Describe the actual tools and their roles: language editing, literature discovery, conjecture generation, proof development, code, data analysis, figures, or formalization. State the checking actually performed by the authors. Do not describe mathematical assistance as mere proofreading when it affected the argument. Do not claim a human audited every step unless that happened.

Consult the target's current policy for wording and placement. Keep this disclosure factual and separate from the mathematical narrative. Never invent a reviewer, an affiliation, an author contribution, or an endorsement.

## Completion note for the user

```text
Changed: current files and purpose
Mathematical effect: changed scope or none
Verified: checks actually run and the conclusions they support
Open: relevant failures, unperformed checks, and remaining obligations
Next action: concrete next step if the requested work is incomplete
```

Report the current result concisely. Keep internal paths here when needed for action; keep them out of manuscript prose.
