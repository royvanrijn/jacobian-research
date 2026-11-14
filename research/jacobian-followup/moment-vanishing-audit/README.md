# Moment and vanishing review

Reviewed through 14 November 2026, following the user's explicit request.
The shared [claim registry](../../MATH_STATUS.json) remains authoritative.
The scope here is this follow-up; archived work items remain unscheduled.

The useful additions are [quadratic finite determination, sharper per-input
cutoffs and a complete ternary quadratic normal form](QUADRATIC_GAUSSIAN.md),
the [compact connected group Mathieu classification](COMPACT_GROUPS.md),
and its [complete disconnected compact Lie group extension](DISCONNECTED_COMPACT_GROUPS.md).
The [complete affine-W classification](CUBIC_AFFINE_W_CLASSIFICATION.md) now
closes an entire twelve-coefficient cubic family, classifies the complete
sixteen-coefficient quartic family and gives a sufficient 5D-3 pure-moment
gate in every degree D>=2 within that architecture. Its checks are in [affine_w_verification.json](affine_w_verification.json).
The earlier exact bounded checks are in [verification.json](verification.json);
scoped metadata checks and repository-wide blockers are in [audit.json](audit.json).

## Literature reconciliation

| Primary source | Current use and boundary |
|---|---|
| [Derksen--van den Essen--Zhao, 1506.05192](https://arxiv.org/abs/1506.05192) | Universal GMC implies JC; a Gaussian counterexample does not reverse this implication. |
| [Long, 2607.18186v1](https://arxiv.org/abs/2607.18186v1) | Explicit failure for all n>=3; one must also read subsequent positive results for n=2. |
| [Wilson, 2607.23887v1](https://arxiv.org/abs/2607.23887v1) | Independent two-variable positive theorem; strict one-sided support gives m>deg Q. Local G2T already closes this question and credits Long's overlapping multitorus theorem. Wilson is an additional public source, not an external review of G2T. |
| [Long's multitorus source](https://github.com/octonion/mathematics/blob/main/gmc/gmc2_stronger_arbitrary_torus.tex) | Already imported as FWT1; one factorial radial variable with arbitrary angular rank. It does not authorize a many-radial extension. |
| [Long, 2607.19012v1](https://arxiv.org/abs/2607.19012v1) | The xz and SU(2) failures were already reproduced. The [public Lean release](https://github.com/long-mathematics/xz-mathieu-su2-counterexamples) reports core completion with supporting omissions; no local Lean build was performed here. |
| [Long, 2609.16178v1](https://arxiv.org/abs/2609.16178v1) | Adds the complete compact connected group classification; see the scoped source audit and derived disconnected extension. |
| [Liu, coherent-state circuits v2.2](https://zenodo.org/records/21689760) | Public record inspected only. Its virtually toral reduction is closed by our existing finite-trace theorem; the new proof supplies the finite-type component argument directly. |
| [Dvorsky, 2608.07338v1](https://arxiv.org/abs/2608.07338v1) | Five-variable homogeneous cubic operator witness already recorded as DVG1; it does not improve the local ternary Delta^6 endpoint. |
| [Liu, 21688251 v3](https://zenodo.org/records/21688251) | Full source read for the quadratic gate and mixed cutoff. The older [21674142 v1](https://zenodo.org/records/21674142) redirects via its newer-version link to this v3; retain version identity rather than infer recency from record numbers. |
| [Liu, 21688601 v1](https://zenodo.org/records/21688601) | Separate earlier formulation with rank-specific pure-threshold sharpness. Do not substitute its older limitations for the v3 statement. |
| [Our GVC paper](https://doi.org/10.5281/zenodo.21803325) | Local GVC2ENV and GVC3HC establish the binary/ternary boundary with their recorded proof assurances; the Delta^6 witness does not settle Delta itself. |

Primary abstract/HTML pages and the quadratic TeX source were inspected.
Searches included Gaussian Moments Conjecture, Generalized Vanishing
Conjecture, factorial torus moments and compact Mathieu results, with recent
date filters and author follow-ups. No later directly relevant theorem was
identified in this pass; this is a search-coverage statement, not exhaustive
bibliographic certification. Some displayed publication and identifier dates
are inconsistent; the links and version labels above are the identities used.

## What is closed and what survives

| Question | Evidence and next gate |
|---|---|
| GMC(2) or binary GVC existence | Closed locally by G2T/GVC2ENV; do not recreate their historical chart searches. |
| Quadratic GMC in any covariance rank | Finite matrix classification and mixed cutoff now audited; no quadratic counterexample search is needed. |
| Compact connected nonabelian Mathieu failure | Closed by Long's all-group construction, including central forms. The componentwise finite-trace deduction now covers disconnected compact Lie groups; minimal witness degrees remain open. |
| Minimum GMC degree in three variables | Quadratics are now completely described. The full affine-W/univariate-coefficient cubic architecture is now excluded with a sufficient twelve-moment gate. A failure needs genuinely omitted coupling terms or an unproved frame reduction; the global question stays open. |
| Homogeneous ternary Delta^k, k<=5 | OP-GVC3-MIN remains open. The degree-eight/ten modular units still require characteristic-zero promotion; GVC3D10F5G closes D=0 only. |
| Ordinary-Laplacian compression | OP-GVC2-RP remains open. DVG4 excludes the canonical one-pair completion; neither the compact transfer nor the quadratic theorem supplies a replacement. |
| Finite pure thresholds beyond quadratics | Source-only Liu v1 support/width claims need their own exact certificates and a version comparison before adoption. No all-degree numerical cutoff is inferred from a bounded prefix. |
| Exact moment-scheme audit | WORK-MOMENT-SCHEMES still needs characteristic, saturation, scheme thickness and multiplier quantifiers on the remaining charts. |

## Reproduction and limits

From the repository root:

```sh
python3 research/jacobian-followup/moment-vanishing-audit/verify_moment_transfers.py
python3 research/jacobian-followup/moment-vanishing-audit/verify_affine_w_classification.py
```

The checker uses only rational arithmetic and is independent of downloaded
producer code. It verifies the Hopf/sphere transfer, ternary quadratic mixed
examples, all-rank chain examples at stated finite bounds and a rank-one
finite-prefix trap. The adjacent written proofs establish the quantified
claims. The imported all-rank sharp *pure* threshold and general Lie theory
are not computationally certified by this receipt.

The latest ordinary rendering and navigation-check attempts fail on KDSQ6s
missing retained artifact consumer; ordinary research search still fails on
missing EC evidence. Earlier checks also reported unrelated source links and
checker drift. See the dated audit receipt for coverage. This follow-up does
not reconstruct unrelated outputs or treat those failures as mathematical
obstructions.
