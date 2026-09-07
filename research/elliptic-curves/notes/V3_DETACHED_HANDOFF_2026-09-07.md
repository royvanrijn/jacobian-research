# Frozen V3 detached handoff

V2 remains unchanged (4585-file original seal passed). V3 protocol:
`eed77c4ac724c0132bad7cd39348da75e6b7ebace2a9986f6715b51ff18d5123`.

**Verified result:** fixed V2 M30 → 31 on chart 101, with independent full
schedule replay and mod-2/3/5 certificates. The successful chart entered through
coset coverage at Babai rank 2286/8192, witness `-42239/7324`.

**Continuing:** fresh M17 replay reached every rank through M28; approximately
590 of 1503 M28 charts had completed without another gain at detachment.
The independent checker has verified the completed chain through M28 and its
entire new landscape (65536 cosets, 1652 exact refinements). The current epoch
is not yet terminal. No claim of autonomous M17→31 is made yet.

Supervisor PID at launch: **2742303**, detached with `nohup setsid`, no terminal.
The previous foreground worker/checker were stopped; the supervisor resumes
their atomic checkpoints without changing policy. It then runs full witness,
rank, metric, preservation and transcript-binding checks and packages the result.
Any failure stops the supervisor; it does not retune or restart a failed policy.
Budgets remain the frozen V3 budgets, not an unbounded new campaign.

From the repository root:

```bash
tail -f research/artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3/replay-M17-run.log
tail -f research/artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3/detached-supervisor.log
```

Final result (created only after all checks pass):
`research/artifacts/generated-results/elliptic-curves/adaptive_visibility_cascade_v3.json`.
Read `runs`, `fixed_M30_target_achieved`, and `M17_target_achieved`; a completed
bounded stall is not rank 31. If the supervisor fails, inspect the corresponding
step log before resuming. Do not launch a second worker while this one is alive.

Policy and calibration details: [V3 note](ADAPTIVE_HALF_LATTICE_V3_2026-09-07.md).
Supervisor: [script](../cas/run_visibility_v3_detached.sh).
After completion, update that note and the mathematical-status registry from
the verified result and regenerate STATUS.md with the repository renderer.
