# Constructive class endpoint and the new 302 parent

Historical construction checkpoint. The later
[positive construction report](TWO_CONSTRUCTED_STRICT_CLASSES_AND_302.md)
supersedes the unmet calibration endpoint below: two additional strict classes
have now been independently certified on the +6 reference. Their rational
solubility and the transfer to302 remain unresolved. This earlier report is
retained for its frozen protocols and evidence archive.

The primary research endpoint is an explicit strict cubic Kummer class
outside the marked generic subgroup, constructed from the equation and
generic sections without exceptional points. **That endpoint is not yet
reached.** The completed MW16-05 calibration attempts below are construction
attempts, not another capacity theorem or a rank predictor.

The newly certified determinant1092 MW17 parent makes curve302 the primary
mechanism target. The +6 fibre remains a smaller calibration control; no
further automatic enlargement of its completed experiments is planned here.

## What the new parent adds

The [source theorem](../notes/CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md)
supplies a full saturated arithmetic generic rank17 basis and literal
specialization at t=0 to302. The
[intake certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_new_parent_intake_v1.json)
checks the source proof's file bindings, evaluates all seventeen sections,
checks smoothness at zero and the squarefree degree24 discriminant, and
computes an exact unimodular change of basis of determinant **-1** from the
previously audited primitive17-core to the new specialization image.

Thus the old core-relative retrospective arithmetic applies to an actual
generic subgroup. Its quotient in the displayed independent group is
**D/sp(M17)=Z^14**. The whole rational group need not equal D, and exact
rank31 of302 remains UNKNOWN. The fourteen quotient rank directions are
not automatically fourteen *strict* classes; that requires the complete
local-boundary accounting.

The new information is the parameter direction preserving all seventeen
sections. The old halving-field, pair-character and component computations
need not be repeated. They do not identify the class-creation event.

The [constructor input](../../artifacts/generated-results/elliptic-curves/rank_jump_curve302_new_parent_inputs_v1.json)
contains only the family equation, seventeen generic sections, zero
parameter and their specialized coordinates. It omits the public31 basis,
embedding matrices and fourteen residual points. The separate intake audit
uses the embedding matrices only to identify the old core retrospectively.
The parent itself was reconstructed retrospectively from curve302 data;
this projection does not turn its provenance into a prospective experiment.

## Completed +6 construction calibration

All arithmetic workers use only the MW16-05 equation at3/17, its sixteen
generic classes, and the eligible early principal-relation pool. Protected,
capped and residual waves with exceptional-point character anchors are
excluded. The reference is known retrospectively to have six extra rational
strict directions; that knowledge is not a class input or stopping target.

| Construction | Bounded outcome |
|---|---|
| Default compact BNF,300 seconds | Timed out; no BNF returned |
| Seeded compact BNF,256MiB | Accepted2378 equation-only seeds; first exact matrix reduction overflowed |
| Seeded relation-only BNF,256MiB | Same matrix-stage memory failure |
| Seeded relation-only BNF,1GiB,300 seconds | Matrix reduction succeeded; six further principal elements; timeout |
| One free-column norm wave,512 targets,8192 strip limit | Completed in175.71seconds;11,452,396 candidate occurrences;383 principal-relation occurrences |

The separate memory protocol records why1GiB replaced the initial cap.
Time, worker count and inputs were unchanged; failed runs were preserved.
The first seeded log parser missed markers prefixed by PARI's progress
counter, so its captured accepted-element count is not reliable. The later
variant starts markers on fresh lines and preserves their exact columns.
This affects logging, not the reported memory failure or class endpoint.

The targeted wave selects from the matrix after precisely four early box
and five early root-strip waves. It never reads the protected continuations.
It uses no target16 stopping rule. Every completed target is checkpointed.
The [witness replay and extraction](../../artifacts/generated-results/elliptic-curves/rank_jump_reference_targeted_class_extraction_v1.json)
rechecks all retained norm identities and6943 individual ideal valuations.
Deduplication leaves207 additional primitive elements beyond the original
4134, including the six from the seeded run.

The extractor retains full outside-S prime blocks, applies
pi(alpha)=Norm(alpha)alpha, and computes its parity kernel. It then combines
those products with G, imposes local squares at every S place and positivity
at the three real embeddings, and checks fixed good-prime characters.
The combined dictionary has4341 elements and parity rank4340. Its seven
strict coefficient dependencies have detected squareclass rank6, all
already supplied by G. There is **no additional independence certificate**.
No conclusion about absence of other fibre classes follows.

This distinguishes three different outputs that must not be conflated:
principal relations were obtained; extra strict classes were not obtained;
global rational solubility of an extra constructed class was therefore not
tested. No elliptic point search or new parameter was run.

## Next mathematical endpoint on302

1. Use the frozen equation/generic projection to compute the complete local
   boundary and generic strict subgroup. Reuse retained equation-only
   factorization evidence where available, but do not infer complete S from
   the list of nontrivial component groups.
2. Construct one beta outside G with square norm, an exact square principal
   ideal and all strict local conditions. A partial class/unit computation
   can supply an individually verifiable beta; full class-group completion
   is not required for that certificate.
3. Compare the resulting ideal/cover identity with the *actual new family*.
   Identify which ramification, principalization or local condition changes
   at zero. Until this is done, neither a specialization condition nor a
   new soluble block is established.
4. Classify rationality separately. A nonzero Cassels--Tate pairing against
   a certified Selmer class obstructs rationality. Zero pairing alone is
   insufficient; a rational point or another complete global argument is
   needed. Public exceptional points may evaluate a frozen construction
   retrospectively, never select its inputs.

No new condition on t, additional strict class, simultaneous-solubility
carrier or Sha classification is claimed in this update. The new parent
closes the missing *family and generic-subgroup identification* for302;
the missing implication is still from its specialized arithmetic to an
explicit new strict class.

## Replay and ownership

The new scripts and `rank_jump_*` artifacts are separate from live search
work. Active search protocols, worker settings, candidates, outputs and
MATH_STATUS entries are unchanged. Construction commands refuse to overwrite
their artifacts; they are not instructions to rerun completed searches.

The principal construction scripts are
[reference_class_targeted_relations.py](reference_class_targeted_relations.py)
and [reference_targeted_class_extraction.py](reference_targeted_class_extraction.py).
The302 intake is [curve302_rank_jump_parent_intake.py](curve302_rank_jump_parent_intake.py).
The protocols alongside them retain their exact input and resource gates.

The [portable checkpoint archive](../../artifacts/generated-results/elliptic-curves/rank_jump_reference_constructive_evidence_v1.zip)
and [manifest](../../artifacts/generated-results/elliptic-curves/rank_jump_reference_constructive_evidence_v1.json)
retain547 input, log and checkpoint files, including the failures.
[Independent local-power checks](../../artifacts/generated-results/elliptic-curves/rank_jump_reference_constructive_checkpoint_verification_v1.json)
verify every strict candidate product through a separate PARI interface.
After restoring the checkpoint files to an empty replay checkout, run:

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_reference_constructive_checkpoint.py check
```
