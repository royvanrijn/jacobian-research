# Determinant-1092 intrinsic low-shell cascade pilot

The terminal certificate is
[`det1092_low_shell_cascade_v1.json`](../../artifacts/generated-results/elliptic-curves/det1092_low_shell_cascade_v1.json).

Three score-stratified determinant-1092 fibres were selected before point
execution, one from each pre-existing score stratum.  On each certified
subgroup the protocol sampled 2,048 SHA-addressed nonzero parity classes,
selected the 49 smallest exact rounded canonical-height representatives, froze
all maps, and searched/replayed/certified the complete wave.  All 147 boxes
completed with rank 17 modulo 2, 3, and 5, so the policy correctly stopped
without an adaptive enlargement.

The completed 48-fibre determinant-1092 deep-centre exposure (2,352 boxes) and
the completed 104-fibre target-free A1/MW16 exposure are retained matched
negative controls.  Neither was read by the new selection or execution.

The requested basis-randomization check is negative: a unimodular shear
preserves the coset norms but only 21 of the 49 sampled selected cosets.  The
SHA parity schedule is coordinate-indexed, so this version must **not** be
described as basis invariant.  A follow-up needs a canonical, basis-free parity
addressing scheme before a wider cascade campaign.
