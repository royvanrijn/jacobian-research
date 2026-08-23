# Q80 final q6 — exact characteristic-zero CM24 certificate

Status: **PASS_EXACT_Q80_FINAL_Q6_CHILD**

- source: selected exact q4 candidate1 parent;
- horizontal: exact `P2-P3`, recovered as a difference of easier exact MW sections;
- resolved RR: ambient `4`, condition rank `2`, `h0(D)=2`;
- exact quartic degree: `4`;
- finite fibres: `4 I3 + I4 + I6 + 2 I1`;
- infinity: smooth;
- root lattice: `4A2+A3+A5`;
- root rank: `16`;
- root determinant: `1944`;
- CM24 MW rank: `2`;
- Euler number: `24`.

This is the characteristic-zero CM24 specialization of the separately certified generic final neighbour `A1/MW16 -> rootless/MW17`.

Reproduce the terminal chain from the repository root with

```bash
sage elkies-k3/scripts/trace_q80_candidate1_marked_transport.sage
sage elkies-k3/scripts/recover_q80_final_q6_via_basis_sections.sage
sage elkies-k3/scripts/certify_q80_final_q6_char0_rr_from_basis.sage
sage elkies-k3/scripts/compile_q80_final_q6_char0_child.sage
```

The final child model is `q80_char0_final_q6_child.sage` in this directory. For the reconstruction history, modular-regression boundary, and superseded direct-Hensel approaches, see [`../../../Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md`](../../../Q80_FINAL_Q6_CLOSEOUT_2026-08-23.md).
