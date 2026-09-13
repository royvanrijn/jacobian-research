# Exceptional ancestry, equation-only arithmetic, and strict descent classes

The full frozen September 12 study, including every carrier, arithmetic row,
strict word, producer protocol, and replay history, is preserved byte-for-byte
in [the archive](../../archive/elliptic-curves/notes/RANK_TRIANGLE_ANCESTRY_AND_DESCENT_2026-09-12.md.txt)
(`sha256: e51982bb641e8a7e7dc0131100025d11ca48dd85005a867dcf4c54f0458e2927`).
It does not authorize a new search, BNF run, descent, or point computation.

The [principal 11952 rank-28 comparison](PRINCIPAL28_EXCEPTIONAL_ANCESTRY_2026-09-12.md)
is the current principal control. This source remains canonical for the bounded
[rank-triangle claim](../../MATH_STATUS.json): the earlier 302-versus-11952
carriers, frozen equation-only panel, and finite strict-class dictionary.

## Certified finite results

- Every one of the selected 14 exceptional 302 points and 8 selected 11952
  points has an explicit degree-two genus-one carrier from its matched
  alternate fibration. Certified generic rank 17 and independence exclude
  degree one, so the minimum multisection degree is two for these chosen
  representatives. Minimum genus and a complete ancestry classification remain
  `UNKNOWN`.
- The exact records cover 504 plus 288 carriers. All 792 carrier identities
  replay, and all 792 altered-carrier controls reject. Within either selected
  fibre the tested quadratic covers are pairwise distinct; this does not rule
  out another shared cover or identify covers after arbitrary base changes.
- The frozen arithmetic panel contains 17 of 2,080 cases plus 302. Ten search
  cases have equation-only Brumer--Kramer offsets 5--9; 302 has offset 12.
  Seven selected search cases remain `UNKNOWN`, and all ten provisional BNF
  probes were unresolved. There is no independent class-group, Selmer, or rank
  upper bound, no enrichment conclusion, and no change to the rank-at-least-31
  or rank-at-least-25 statements.
- Exact Kummer/local replay gives `dim W=31`, `dim G=17`,
  `rank loc(W)=21`, `rank loc(G)=17`, hence a known strict kernel of dimension
  10 and a four-dimensional localized quotient. The retained ten strict words
  and half ideals use the fixed basis; only E3 and E7 individually have a
  generic strict correction in that complement. This is not the full strict
  Selmer kernel or the full Selmer group.

The primary replay evidence is the
[targets](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/targets.json),
[geometry](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/geometry.json),
[combined verification](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/verified.json),
[strict dictionary](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/strict.json),
and [manifest](../../artifacts/generated-results/elliptic-curves/rank_triangle_v1/manifest.json).
The companion [Curve302 descent anatomy](CURVE302_DESCENT_ANATOMY_2026-09-10.md)
contains the ideal-square and independence evidence.

## Reuse boundary

Reuse the already certified carriers and strict dictionary; do not recompute
them or promote this finite panel to a rank mechanism. A new result would need
an explicitly scoped common-cover or genus-zero construction, a stated
equivalence/generation argument, or an independent exact upper-bound
certificate. The [ancestry method record](../../knowledge/ALGORITHMS.md#method-ec-ancestry-reuse-the-bounded-ancestry-panel-without-promoting-it-to-a-rank-mechanism)
holds the operational rule.

With the recorded Sage environment, the narrow verifier is
`sage -python research/elliptic-curves/cas/verify_rank_triangle.sage`; the
archive retains the full command sequence and its historical failure boundary.
