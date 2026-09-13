# R17 small-field class-quotient laboratory: frozen Phase 0

Phase 0 is a frozen, rank-blind cohort of 100 ordinary R17/MW17 fibres.  Its
[cohort artifact](../../artifacts/generated-results/elkies-k3-r17-small-field-class-quotient-cohort-v1.json)
contains no class-quotient feature and no point-search outcome.  Phase 1 was
never run; every feature value and detector result remains `UNKNOWN`.  This is
not a claim in [`MATH_STATUS.json`](../../MATH_STATUS.json) or authorization
for a BNF, Selmer, or point-search campaign.

For a future separately scoped experiment, the retained feature is

\[
 Q_t=\bigl(\operatorname{Cl}(K_t)/\langle S_t\rangle\bigr)[2]
 /\langle c_{S_t}(G_1),\ldots,c_{S_t}(G_{17})\rangle.
\]

It is deliberately different from quotienting the class group by doubles:
the latter loses the embedded 2-torsion information when 4-torsion is present.
The corresponding exact Smith-form, localization, unconditional-BNF, blinding
and censoring requirements are the current
[localized-2-torsion lesson](../../knowledge/ALGORITHMS.md#method-ec-localized-2torsion-quotient).

The full selection rule, fixtures, abandoned shard protocol, expansion boundary
and historical commands are preserved byte-for-byte in the
[archive](../../archive/elliptic-curves/notes/R17_SMALL_FIELD_CLASS_QUOTIENT_LAB.md.txt).
Do not launch Phase 1 merely because the cohort or code exists.  It requires
complete unconditional feature records for the whole frozen cohort, a
feature-file-bound detector protocol, and a separate evidence plan.
