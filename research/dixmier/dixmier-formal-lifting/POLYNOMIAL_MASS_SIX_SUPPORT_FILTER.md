# A polynomial Weyl mass filter for candidate pairs

[Shaik, Hexagon 2610.00273v1, Theorem 1.1](sources/shaik-2610.00273v1.pdf)
states: for any characteristic-zero field `K` and actual
`P,Q∈A₁(K)` with `[Q,P]=1`, if `P` has mass at most six, then
`K⟨P,Q⟩=A₁(K)`. The same conclusion holds with `Q` in place of
`P`. Here mass is the number of occupied components in the
`Z`-grading `deg X=1`, `deg Y=−1`: in normal PBW support it is
`#{i−j : coefficient of X^iY^j is nonzero}`. It counts grades, not
PBW monomials, differential order, or total degree.

Consequently a genuine nongenerating *polynomial* Weyl pair must
have at least seven occupied grades in **each** member. If a Weyl
automorphism sends such a pair to another polynomial pair, the
image remains nongenerating and both image masses must again be at
least seven. This is a necessary support filter only. The paper
also states a prime-outer-exponent face exclusion (Theorem 1.3),
whose separate hypotheses need checking before any application.

The paper reports a Lean formalization of Theorem 1.1. The primary
statement and argument were inspected; the Lean project and its
transitive axioms were not replayed here. Thus this is a scoped
imported theorem, with no independent formal-verification claim.

Apply the filter to
[actual simultaneous polynomial subrectangular candidates](CENTERED_46_RECTANGLE_NEWTON_FACE_GATE.md)
*after* constructing both full operators. A face calculation alone
does not determine mass of the operator. The displayed order-(4,6)
rectangle allows grades from `−4` through `8` for `P` and `−6`
through `12` for `Q`, before coefficient constraints. Its surviving
square/cube top faces can have only three and four grades, but the
unfixed lower coefficients may occupy more than six. The existing
rectangle therefore is not removed by the theorem. Its total-degree
gcd six also lies below the separate published degree bound noted
there.

For the active
[subprincipal constraint](PBW_TRANSPORT_AND_FINAL_COEFFICIENT_COMPARISON.md),
the actionable check is: once a simultaneous pair has actual finite
PBW support in the original polynomial `A₁(K)`, compute both exact
grade sets, including any admissible polynomial automorphism used in
normalization. A mass at most six then excludes that candidate as a
nongenerating pair. Rational/Laurent pairs, matrix carriers, formal
series and pairs over localized coefficient rings require a separate
polynomial admission argument before this filter applies. Nothing in
the mass theorem evaluates the centered Weyl corner `A_W`.
