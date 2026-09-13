# Relative 2-Selmer pipeline — historical boundary

The complete resource-bounded R17 protocol, command catalogue and timeout
receipts are preserved byte-for-byte in [the archive](../../archive/elliptic-curves/notes/ELKIES_R17_RELATIVE_2SELMER_PIPELINE.md.txt)
(`sha256: b3ef0d3d1cdd87c80cd4688f97869c4af621007a4fc1e4c68e8260d943bbce34`).
It does not authorize a fresh BNF, Selmer, cover, or point-search campaign.
Mathematical status remains in [MATH_STATUS.json](../../MATH_STATUS.json).

## Retained evidence and boundary

The class-group-free known-subgroup audit certifies mod-two dimensions
`21,25,26,27,28` on the five public controls and dimension 17 for MW17 on all
ten frozen high-Nagao candidates. Those are exact lower-bound labels, not full
Selmer quotients. The [input manifest](../../artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json)
and [run receipt](../../artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_run_v1.json)
record that every original Magma job was unavailable and the first PARI/Sage
controls timed out before a certified BNF or `ell2cover` basis. No run returns
a Selmer dimension, rank upper bound, unrealized class, or candidate result.

For record fibres 356 and 385, retain a certified global squareclass upper
envelope, all 29 exact Kummer rows, and certified local maps. Extend the MW29
rows to the global envelope before forming local conditions. A zero kernel
then proves the residual quotient is zero; a nonzero kernel with missing
places or global data is only a monotone upper bound. The odd MW29 dimension
and certified even total Selmer parity let a certified raw upper bound of one
sharpen to zero. No current certified BNF or equivalent global provider exists,
so both residual groups remain `UNKNOWN`.

The current reusable interfaces are
[`build_mw29_relative_2selmer_matrix.py`](../cas/build_mw29_relative_2selmer_matrix.py),
[`audit_mw29_relative_selmer_witness_bound.py`](../cas/audit_mw29_relative_selmer_witness_bound.py),
and the [higher-2-power programme](R17_RECORD_PAIR_HIGHER_2POWER_SELMER_PROGRAM.md).
Read [ALGORITHMS.md](../../knowledge/ALGORITHMS.md) before any separately
scoped arithmetic provider; a bounded relation plateau or timeout is not a
global envelope.

The complete bisection control corpus shows that the known mechanism recovers
all four rank-21 quotient directions but spans only `5,3,2,1` directions on
the rank-25--28 controls. It is therefore not a complete exceptional-direction
detector. Do not restart translated-bisection shells or generic backend retries
without either a certified product-twist point or a complete residual
class-field computation.
