# Historical half-lattice chart calibration — 4 September 2026

The September 4 blind rank-28 replay, its fixed-basis class-order ablation,
and its CRT continuation are historical. The full source, tables, commands,
input hashes and bounded-search receipts are retained in
[`archive/elliptic-curves/notes/HALF_LATTICE_FAKE_DESCENT_REPLAY_2026-09-04.md.txt`](../../archive/elliptic-curves/notes/HALF_LATTICE_FAKE_DESCENT_REPLAY_2026-09-04.md.txt)
(SHA-256 `3bcd8e39223c615c439d3aa76dbc1e9e75175d037b8f5024c1a9219804e1c1ab`).

The current mathematical boundary is
[`EC-K3-R17-074D9-HALF-LATTICE-PROMOTION-GATE`](../../elkies-k3/R17_PROSPECTIVE_CRT_RANK_JUMP_EXPERIMENT_2026-09-04.md):
the frozen protocol cannot promote a rank-32 candidate. Its current
[promotion certificate](../../artifacts/generated-results/elkies-k3-r17-prospective-crt-half-lattice-promotion-gate-v1.json)
and [equal-budget ablation](../../artifacts/generated-results/elliptic-curves/half_lattice_search_ablation_summary_v1.json)
are the relevant compact evidence. Exact maximum-class geometry and the
former 49-versus-43 cap are in
[EXACT_PARITY_AND_COORDINATE_AUDIT_2026-09-06.md](EXACT_PARITY_AND_COORDINATE_AUDIT_2026-09-06.md).

## Preserved boundaries

A pointed quartic here is a birational search chart of the same elliptic
curve, not a nontrivial two-covering torsor or Selmer-class representative.
Its local solubility and a chart miss therefore imply neither a Selmer fact,
rational-point absence, saturation, nor a rank bound. Keep the birational
point map distinct from the actual covering map as specified in
[ALGORITHMS.md](../../knowledge/ALGORITHMS.md).

The historical `depth`, `old-deep-43`, and quotient-weight fields order charts
only in their pinned lattice presentation. Any basis, finite-index, coordinate,
or independent-point change invalidates that order and requires a separately
blinded calibration. Bounded recovery performance is not a prospective
rank-jump predictor. No command in the archived record authorizes a restart.

<!-- status-consumer: EC-K3-R17-074D9-HALF-LATTICE-PROMOTION-GATE 9a1f080523c9ecae -->
