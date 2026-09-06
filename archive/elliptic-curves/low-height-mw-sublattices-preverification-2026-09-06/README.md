# Historical preverification outputs

Use the [active result](../../../elliptic-curves/notes/LOW_HEIGHT_MW_SUBLATTICES_2026-09-06.md).
These snapshots are retained as execution history, not active certificates.

The initial GP scorer needed an empty-kernel case at full ambient rank and a
row/column correction for minimum witnesses. Subsequent verification found
that bounded `qfminim`'s second output was being reported as the minimum;
it can instead be the largest returned norm. The active script explicitly
minimizes returned norms and selects a corresponding witness. Proposal and
candidate checkpoints, rankings and failed calibration outcomes were unchanged.

The first successful arithmetic replay and its summary are also retained
here because the final verifier adds independent 100-digit ball enumeration
and sparse-pool reconstruction. Their mathematical checks were not retracted.
Historical hashes refer to the files as they existed during those runs.
