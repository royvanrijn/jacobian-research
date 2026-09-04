# Historical audit input

The gzip file preserves the exact zero-gain rescue protocol before the
2026-09-04 external audit rejected infinite and NaN search budgets in
`production_search_gates.py`. The arm was frozen but unrun. Its active
protocol was regenerated with only that implementation hash and the derived
protocol hash changed; candidates, assignments, ordering and limits are identical.

Original uncompressed SHA256: `4ad41a89a43205f657e72e0cedd222777e7906c02b57b435aa60794214685188`.

The historical [budget validator](production_search_gates.before.py.txt)
matches the old protocol's source hash. The historical
[finite-reduction helper](rank_certification.before.py.txt) matches the
Fermigier generic-rank artifact's recorded producer dependency. Both copies
were checked against the audit's initial SHA256 inventory. They are replay
inputs, not active implementations.

See the [audit](../../../elliptic-curves/notes/EXTERNAL_AUDIT_2026-09-04.md).
