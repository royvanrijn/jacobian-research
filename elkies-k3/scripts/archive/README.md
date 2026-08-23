# Archived Elkies K3 scripts

This directory preserves exploratory and superseded scripts that explain how the current
H3 and Q80 constructions were reached. It is historical evidence, not a proof entry-point
directory.

The archive was imported on 2026-08-23 from material predating the current q8 and Q80
closeouts. The audit at commit `4eac04a442a132b696a85384f08b81569870a940`
enumerated the complete batch and grouped it by research question. Unique failed scripts
remain because their assumptions, local calculations, or negative results are useful.

When an identically named file exists in `../`, the root copy is authoritative and this
copy is an older snapshot.

## A. H3 q8 marking, chamber, and q24-bridge audits

These scripts established or challenged the exact divisor marking before the corrected
q8 compiler existed:

```text
analyze_h92_q8_component_nef_s3_target.sage
analyze_h92_q8_visible_zero_profiles.sage
analyze_h92_q8_visible_zero_profiles_v2.sage
audit_h3_d13_mw_basis_profiles.sage
audit_h3_d13_q24_pullback_invariants.sage
audit_h3_d13_q24_section_pullback*.sage
audit_h92_q6_child_q8_marking_height.sage
audit_h92_q6_child_q8_point_heights.sage
audit_h92_q6_third_literal_support.sage
audit_h92_q8_d13_zero_anchor*.sage
audit_h92_q8_explicit_old_curves.sage
audit_h92_q8_old_q6_mw_images.sage
audit_h92_q8_original_h3_sections.sage
audit_h92_q8_representative_selection.sage
certify_h92_q6_third_to_q8_bridge*.sage
certify_h92_q8_equation_ns_divisor.sage
derive_h92_q8_d13_branch_anchor.sage
derive_h92_q8_d13_g3_from_e77_bisection.sage
derive_h92_q8_global_orientation*.sage
probe_h92_q6_child_q8_corrected2cover_global_modp*.sage
probe_h92_q6_third_p2_trace_modp*.sage
probe_h92_q8_complete_chamber_reduction.sage
search_h92_q6_mw_for_q24_bridge*.sage
search_h92_source_mw_for_q24_bridge*.sage
```

Outcome: the divisor and chamber questions were real, but the final q8 obstruction was
not geometric. It was two implementation errors: treating a binary-quartic covariant map
as a primitive MW map rather than a 2-cover map, and omitting the denominator factor
`Dx` in q-frame CRT normalization. The authoritative replacement is
`../derive_h92_q6_child_q8_corrected2cover_qq.sage`; see
[`../../H3_Q8_REAUDIT_2026-08-22.md`](../../H3_Q8_REAUDIT_2026-08-22.md).

The archived `derive_h92_q6_child_q8_corrected2cover_qq.sage` is an older snapshot and
must not be used in place of the root copy.

## B. Superseded `corrected1278` q8 route

```text
probe_h92_q8_corrected1278_all_nodes_direct.sage
probe_h92_q8_corrected1278_e7_identity.sage
probe_h92_q8_corrected1278_e7_identity_full_jets.sage
probe_h92_q8_corrected1278_global.sage
probe_h92_q8_corrected1278_global_v2.sage
probe_h92_q8_corrected1278_marked_caps.sage
probe_h92_q8_corrected1278_two_translated_divisors.sage
probe_h92_q8_corrected1278_y_root_full_jets.sage
probe_h92_q8_corrected1278_zu_root_full_jets.sage
```

This hand-built degree-16 / `q6^8` lane attempted to impose the resolved E7 conditions
through explicit local jets and translated divisors. It was useful for discovering which
local modules were missing, but it solved the wrong enlarged marking problem. The exact
corrected q8 RR system is only `13 -> 2` over `QQ`.

## C. Superseded `true1600` and old q8 normalization routes

```text
probe_h92_q8_true1600_all_nodes_direct.sage
probe_h92_q8_true1600_e75_full_normal_jets.sage
probe_h92_q8_true1600_first_finite_node*.sage
probe_h92_q8_true1600_full_chart_power_modp.sage
probe_h92_q8_true1600_p59_global.sage
probe_h92_q8_true1600_restricted_node_direct.sage
probe_h92_q8_true1600_restricted_node_fast.sage
probe_h92_q8_true1600_translated_divisor.sage
probe_h92_q8_true1600_two_translated_divisors*.sage
probe_h92_q8_true1600_y_root_full_jets.sage
probe_h92_q8_true1600_zu_root_full_jets.sage
probe_h92_q8_horizontal_fixed_reduction.sage
probe_h92_q8_lower_twist_multiplication.sage
probe_h92_q8_lower_twist_tail.sage
probe_h92_q8_s3_trace_modp*.sage
probe_h92_q8_source_nef_multisections.sage
probe_h92_q8_true_twist_ladder.sage
q8_global_lattice_intersection_diagnostic.sage
q8_global_lattice_intersection_modp_v2.sage
q8_global_true_e7_all_nodes_gate.sage
q8_global_true_e7_generic_gate*.sage
q8_global_true_e7_node43_gate.sage
q8_qframe_e7_node_probe*.sage
q8_qframe_finite_transition_probe.sage
repair_h92_q6_child_q8_marking_2cover.sage
```

The `true1600` lane used a source-side degree-18/large-cap presentation; other scripts
used the old degree-46 marking or a q-frame normalizer that silently left a vertical pole.
They remain diagnostic evidence for the corrected denominator and 2-cover rules, not
alternative current compilers.

## D. q6 third-section and singular Hensel experiments

```text
check_h92_q6_third_rr_stability.sage
inspect_h92_q6_third_hensel.sage
lift_h92_q6_third_hensel.sage
rebase_h92_q6_third_iistar_ivstar*.sage
recover_h92_q6_standard_zero_translation_v2.sage
recover_h92_q6_third_common_denominator.sage
recover_h92_q6_third_rational_adaptive.sage
recover_h92_q6_third_trace_modp.sage
recover_h92_q6_third_trace_samples_modp*.sage
```

These tried to recover the difficult third q6 section directly in singular residue
coordinates. They exposed non-transversality and useful trace/translation machinery, but
the direct Hensel lane was not the robust construction. The current route uses exact
section descent and marked transport; the surviving trace ideas are reused by the active
q24 interpolation script.

## E. Guessed E6/MW3 chart and component searches

```text
build_mw3_semistable_p1_checkpoint.sage
build_mw3_semistable_p1_system.sage
check_e6_rank2_components_gf31.sage
e6_i4_blowup_classifier*.sage
extend_mw3_p1_p110_norm_state.sage
mw3_p1_stage2_chain.sage
mw3_semistable_component_profiles.sage
probe_mw3_i11_formal_root_state.sage
probe_mw3_p1_residuals.sage
probe_mw3_p1_stage3_state.sage
rank_mw3_p1_i2_five_from_state.sage
scan_e6_p3_first_order_31adic.sage
scan_e6_p3_second_order_31adic.sage
search_e6_third_polynomial_gf31.sage
```

The E6 frame itself is a valid lattice backtrack. The split Weierstrass chart attacked by
these scripts was inferred from Kodaira data and was never obtained by applying the
certified neighbour transport. Later exact glue recovery added the missing condition
`(P1+P2).O=1`, equivalently `<P1,P2>=-5/6`; all old finite-field survivors fail it.
These scripts therefore preserve a valuable negative result, not a source reconstruction.

## F. q8 root-lattice/module certification experiments

```text
certify_q8_e7_round10_pure_m.sage
certify_q8_e8_resolved_reduced_lattice.sage
certify_q8_e8_saturated_lattice.sage
```

These are historical local-module checks that helped distinguish a root-lattice basis
from the resolved quotient conditions actually needed by the compiler. The durable rules
have moved into the field-generic compiler and its verifier in `../`.

## Exact duplicates removed by the audit

The following files had the same Git blob SHA as the corresponding unnumbered file and
were deleted:

```text
probe_h92_q8_corrected1278_e7_identity (1).sage
probe_h92_q8_corrected1278_global_v2 (1).sage
probe_h92_q8_true1600_translated_divisor (1).sage
probe_h92_q8_true_twist_ladder (1).sage
probe_mw3_p1_residuals (1).sage
```

No unique historical script was deleted or moved.

## How to use this archive

Use it to answer questions such as:

- Which local condition first exposed a bad marking?
- Was a failed lift transverse or singular?
- Which finite-field or jet calculation ruled out an apparent shortcut?
- Where did a later compiler rule originate?

Do not cite an archived PASS line as a current theorem. Reproduce the claim through the
root-level verifier named by a current note, and treat archive output as provenance only.