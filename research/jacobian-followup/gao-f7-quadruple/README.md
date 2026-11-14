# F7: complete quadruple-root strata

**GAOF7FOUR2:** all normalized F7 targets with a quadruple-or-higher
eliminant root have type4+1^8,4+2+1^6, or5+1^7. Exactly **912** have
type4+2+1^6 and six-point fibres; exactly **118** have type5+1^7 and
seven-point fibres. All others have eight-point fibres. The complete
nonordinary locus therefore contains **1030** normalized targets.

The three-point type **(4,3,2) is excluded**. Only **(3,2,2,2)**
remains for global three-point F7 existence. Full F7 strata, F5 and JC2
remain open. These counts are complex geometric targets after weighted
normalization; original target sets include their C* scaling curves.

[Complete proof](COMPLETE_PROOF.md) · [combined replay](verify_complete.py) ·
[receipt](complete-receipt.json) · [exact finite presentation](finite-degenerations.json).

The proof certifies an integer norm using a determinant height bound,
345 primes and8052300 determinant comparisons, then uses its squarefree
degree912 factor and a local norm-multiplicity lemma. All240
presentation-pole targets are covered, including216 newly replayed
points in the remaining charts. See the
[discriminant receipt](discriminant-receipt.json),
[norm receipt](norm-grid-receipt.json) and
[chart receipt](remaining-charts-receipt.json).

The [GAOF7FOUR1 predecessor lemma](PROOF.md) proves generic degree eight,
finite fibres and the finite degeneration locus; its
[original pole lemma](EXCEPTIONAL_PROOF.md) covers24 targets.
The [campaign record](campaign.json) preserves inconclusive attempts.
[Dependencies](requirements.txt) include a C compiler and matched FLINT
headers/library for the norm replay. Independent end-to-end, formal,
external and novelty review remain outstanding.
