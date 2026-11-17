# Absolute Sunada pairs of Keller maps: research audit

This note records a literature and construction audit completed on
2026-07-27.  It starts from the verified
[global Davenport--Sunada construction](GLOBAL_SUNADA_KELLER_COVERS.md)
and asks for two polynomial maps

\[
F_P,F_L:\mathbb A^n_K\longrightarrow\mathbb A^n_K
\]

with constant nonzero Jacobian, generic degree seven, a common
\(\Gamma=\operatorname{GL}_3(\mathbb F_2)\) Galois closure, and generic
inverse covers corresponding to the point and line stabilizers.  The desired
maps should remain stably polynomially left--right inequivalent.

## 1. Bottom line

The group theory and arithmetic comparison are complete.  The missing
theorem is a **prescribed-monodromy polynomial algebraization theorem**:
realize the two fixed degree-seven function-field extensions as the
regular-reconstruction opens of polynomial Keller maps with affine-space
source and target.

The July 2026 counterexample to the Jacobian conjecture changes the context
but does not close this gate.  Absolute Keller maps now exist in dimension
three in every generic degree at least three.  Their known marked-root
families have generic symmetric monodromy \(S_N\).  For \(N=7\), that does
not produce the nonconjugate index-seven Gassmann pair: the symmetric and
alternating degree-seven actions have no such partner.  A degree-seven
Keller map is therefore no longer the existential difficulty; forcing the
specific nonsymmetric \(\operatorname{GL}_3(\mathbb F_2)\) inverse cover is.

The current Davenport attack remains honest and sharply localized.  Its
leading unresolved gate is still the translated determinant incidence from
the [post-coordinate audit](DAVENPORT_POST_COORDINATE_ATTACKS.md), or a
genuinely new monodromy core whose boundary already matches one of the known
absolute marked-root algebraizations.

No paper located in this audit realizes a nontrivial Gassmann pair as the
generic inverse covers of two polynomial Keller self-maps of affine space.

## 2. What an absolute theorem must say

Let

\[
k=K(y_1,\ldots,y_n),\qquad
L_P=K(x_1,\ldots,x_n),\qquad
L_L=K(x'_1,\ldots,x'_n)
\]

be the function fields induced by \(F_P,F_L\).  The exact generic theorem
requires:

1. \([L_P:k]=[L_L:k]=7\);
2. the two normal closures coincide over \(k\);
3. their common Galois group is
   \(\Gamma=\operatorname{GL}_3(\mathbb F_2)\);
4. the two subfields are fixed by the nonconjugate point and line
   stabilizers;
5. \(\det DF_P,\det DF_L\in K^*\); and
6. the distinguished affine opens in the two finite normalizations are
   actually affine spaces.

Items 1--4 give the Sunada/Gassmann conclusion.  Item 5 says that all finite
ramification of the normalization lies outside the affine source.  Item 6
is the algebraization gate.

The known absolute degree-seven maps do not meet items 2--4.  The known
Davenport Cox maps meet items 1--5 but fail item 6.

## 3. Arithmetic scope: good fibers versus every fiber

The common permutation character gives more than equal point counts.  Fixed
counts for every power of a Frobenius element recover its complete cycle
partition, hence the residue-degree multiset and zero-dimensional zeta
function.  Therefore an absolute Gassmann pair automatically has identical
zeta functions at every closed point of a common good finite-etale open.

This does **not** by itself compare fibers on the nonproperness boundary.
A polynomial Keller map is etale on affine source but is not finite; sheets
can escape at infinity above a target hypersurface.  Two Gassmann-equivalent
normalizations can have different distinguished affine opens, different
nonproperness sets, or different surviving sheets on those sets.

Accordingly there are two possible final claims:

- **good-fiber indistinguishability:** equality on one common dense open;
- **all-target indistinguishability:** equality for every good-reduction
  target point, including points on the nonproperness boundary.

The first follows from the Gassmann pair.  The second needs a
boundary-decorated Gassmann theorem: the omitted primes, their residue
degrees, and their incidence over each boundary stratum must also match.

The phrase “arithmetic point-count data cannot distinguish the maps” should
state explicitly which version is proved.

## 4. Stable inequivalence is a separate obligation

The current relative weighted maps are separated by their normalized
Hessian divisors.  A new absolute affine modification need not preserve the
same Hessian polynomial.  Retaining only the generic point/line function
fields is not enough to transfer the Hessian obstruction automatically.

An absolute construction must therefore do one of two things:

1. remain inside a suspension formalism for which the existing Hessian
   functoriality theorem applies; or
2. replace the Hessian argument by an intrinsic invariant of the finite
   normalization together with its distinguished affine open.

The decorated normalization boundary is the more robust candidate.  It can
record ramification index, residue degree, omitted versus affine sheets,
conductor data, and boundary intersections, all of which are stable under
adjoining identity variables.

## 5. Why the new absolute maps do not solve the problem formally

The foundational threefold map and its all-degree weighted extensions have
inverse equations of the form

\[
H(W)-BCW+cAC^2=0.
\]

For an admissible degree-\(N\) seed, this gives a polynomial Keller map of
\(\mathbb A^3\) with generic degree \(N\).  The generic monodromy is
\(S_N\).

Three tempting deductions fail:

1. **Choosing \(N=7\).**  This gives \(S_7\), not the Fano
   \(\operatorname{GL}_3(\mathbb F_2)\) action.
2. **Restricting to the Davenport pencil.**  The restriction recovers the
   desired covers only as pullbacks along a positive-codimension target
   slice.  It does not change the full inverse cover of the absolute map.
3. **Finite base change reducing \(S_7\) to \(\Gamma\).**  The resulting
   resolvent base is a new finite cover of affine space.  Proving that its
   source, target, and pulled-back reconstruction open are affine spaces is
   another form of the original algebraization problem.

The July 2026 constructions are nevertheless important design evidence:
ramification can be moved entirely to the normalization boundary by a
polynomial reciprocal or weighted chart.  The missing ingredient is a chart
adapted to the Fano cover rather than to a generic one-variable seed.

## 6. Certified Davenport boundary

The present point cover has

\[
\Delta(T,g_T(Y))=E_{3,g}E_{6,g}J_g^2,\qquad J_g=g_T'(Y),
\]

and the line cover has the conjugate factorization.  Thus the height-one
ledger is

\[
\begin{array}{c|ccc}
 &E_3&E_6&J\\ \hline
\pi^*\Delta&1&1&2\\
\operatorname{Jac}(\pi)&0&0&1.
\end{array}
\]

The integral lattice is saturated and admits a one-row unimodular
completion.  The obstruction is geometric:

- stable straightening preserves the nonzero unit ranks;
- every suspension retaining \((T,g_T(Y))\) has determinant divisible by
  \(J_g\);
- the natural \(E_3\)-coordinate introduces a fourth coprime divisor;
- the two tangent charts cover the marking line but their standard weighted
  lift has a forced \(1/C\) overlap pole;
- the unique affine-plane node charts glue through a complete exceptional
  \(\mathbb P^1\);
- the direct determinant splice has Jacobian \(D^3/2\);
- the unshifted determinant output has a critical axis; and
- the translated incidence reduces to a plane map unless both the old
  \(T\)-coordinate and linear \(U\)-dependence are genuinely absorbed.

The local checkers for the global cover, Cox boundary, proportional tangent
atlas, and post-coordinate attacks all pass.

## 7. New audit: the LaMacchia sign family

The literature contains another two-parameter
\(\operatorname{GL}_3(\mathbb F_2)\) family.  Write

\[
\begin{aligned}
f_{s,t}(X)={}&X^7+(-6t+2)X^6+(8t^2+4t-3)X^5\\
&+(-s-14t^2+6t-2)X^4\\
&+(s+6t^2-8t^3-4t+2)X^3\\
&+(8t^3+16t^2)X^2+(8t^3-12t^2)X-8t^3.
\end{aligned}
\]

LaMacchia proved that its generic Galois group is
\(\operatorname{GL}_3(\mathbb F_2)\).  Bosma--de Smit proved that
\(f_{s,t}\) and \(f_{-s,t}\) give the two degree-seven Gassmann actions.

Because \(s\) occurs linearly,

\[
f_{s,t}(X)=f_{0,t}(X)+sX^3(1-X).
\]

On \(X(1-X)\ne0\), its root incidence is the rational plane map

\[
(t,X)\longmapsto(t,\phi(t,X)),\qquad
\phi=-\frac{f_{0,t}(X)}{X^3(1-X)}.
\]

Exact factorization gives

\[
\phi=
\frac{(X+1)A(t,X)B(t,X)}{X^3(X-1)},
\]

where

\[
\begin{aligned}
A={}&2t^2-2tX^2+2tX+X^3-X^2,\\
B={}&-4tX^2+8tX-4t+X^3+2X^2-2X.
\end{aligned}
\]

This initially looks friendlier than the Davenport chart: it has two
explicit pole divisors \(X=0,1\).  The exact audit shows that it does not
remove the essential ledger:

\[
\frac{\partial\phi}{\partial X}
=\frac{J_8(t,X)}{X^4(X-1)^2},
\]

and the degree-seven discriminant is

\[
\operatorname{Disc}_X(f_{s,t})
=64t^6\,\mathcal B(s,t)^2.
\]

Pulling back the non-coordinate reduced branch factor gives

\[
\mathcal B(\phi(t,X),t)
=\frac{E_4(t,X)E_8(t,X)J_8(t,X)^2}
       {X^{12}(X-1)^4}.
\]

Thus the same structural pattern \((1,1,2)\) reappears inside that branch,
with degrees \((4,8,8)\) instead of \((3,6,6)\).  The full reduced branch
also contains \(t=0\), where three roots move onto the omitted \(X=0\)
boundary, so the total boundary is not simpler than the Davenport one.

There is also a decisive stable-equivalence warning:

\[
(t,X)\longmapsto(t,-\phi(t,X))
\]

is obtained from the first chart by the polynomial target involution
\((t,s)\mapsto(t,-s)\).  The raw \(f_{s,t}/f_{-s,t}\) presentation is
therefore unsuitable by itself for proving stable left--right inequivalence,
even though the two covers are nonisomorphic over the fixed target.

The LaMacchia family is valuable as a simpler testbed for a universal
three-column algebraizer.  It is not a shortcut to the requested
inequivalent pair.

Run

```bash
.venv/bin/python scripts/audit_lamacchia_sunada_core.py
```

for the exact factorization and target-sign checks.

## 8. Ranked research routes

### A. Generalized marked-line algebraization of the Davenport core

This is the best route that directly uses the 2026 absolute-map mechanism.
The standard weighted model forces the incidence slope to be a product
\(BC\), which creates the \(1/C\) transition pole.  The marked-line identity

\[
(W,\sigma)\longmapsto
\bigl(\sigma,Y(W)-\sigma X(W)\bigr),
\qquad
J=-(Y'-\sigma X'),
\]

shows that a product presentation is not intrinsic.

The exact target is to choose polynomial horizontal and vertical
coordinates \(X(W,s),Y(W,s)\), plus a reciprocal source chart, such that:

1. the two Davenport marking charts differ by polynomial source and target
   automorphisms on their overlap;
2. the additive slope/intercept cocycle is absorbed without inverting a
   boundary coordinate;
3. the transversality factor is the required derivative divisor;
4. the full determinant ledger is constant; and
5. the inverse equation remains the Davenport degree-seven cover, not a
   generic \(S_7\) pencil.

This strictly contains the standard weighted-product ansatz already ruled
out.

### B. Translated determinant incidence

Continue the leading attack in
[DAVENPORT_POST_COORDINATE_ATTACKS.md](DAVENPORT_POST_COORDINATE_ATTACKS.md).
The first live affine-in-\(U\) coefficient pencil must have an alternating
triangular coordinate of polydegree at least two.  Otherwise the invariant
unit gate factors the Jacobian through a nonunit.  A nonlinear-\(U\)
coupling is the other live class.

Every candidate must be rejected immediately if it:

- retains \(T\) as a target coordinate;
- makes \(S\) separately recoverable;
- exposes \(T+Y^2\);
- is affine-linear in mask variables;
- is one-triangular; or
- reduces after translation to a plane Keller pair.

### C. Fano marked-incidence slice

The foundational cubic map can be understood as “mark a root and forget the
mark” on one exceptional normalized factorization slice whose source is
\(\mathbb A^3\).  The analogous group-theoretic search is:

1. build a \(\Gamma\)-cover from a Fano-plane configuration space;
2. form its point-marked and line-marked incidence spaces;
3. choose one common affine slice of the unmarked quotient;
4. test whether both regular-reconstruction opens have class
   \(\mathbb L^n\), trivial units, and polynomial coordinate rings; and
5. compute the residue Jacobians before attempting explicit coordinates.

This bypasses the Davenport derivative curve completely.  It is
speculative, but it aligns the monodromy and algebraization mechanisms at
the start instead of imposing \(\Gamma\) on a generic \(S_7\) map later.

### D. A different \(\Gamma\)-cover with an \(\mathbb A^1\) or
\(\mathbb G_m\) critical normalization

The degree-seven Gassmann pair is essentially unique at the group-action
level, but the \(\Gamma\)-torsor realizing it is not.  Search LaMacchia,
rigid, and other parametric \(\operatorname{PSL}_2(7)\) families for one
whose selected ramification normalization already has logarithmic type
\(\mathbb A^1\) or \(\mathbb G_m\), matching the two known absolute
algebraizations.

The fast prefilter is the reduced branch pullback.  A candidate with a
three-puncture normalization or another \((1,1,2)\) ledger has not improved
the essential boundary geometry.  The LaMacchia calculation above fails
this prefilter and also carries a target-sign equivalence.

## 9. Concrete next computations

1. **Independent marked-line reciprocal ring.**  The slope/intercept
   cocycle is now solved by a determinant-one triangular automorphism.  The
   elementary reciprocal modification has unavoidable \(D^{-3}\) and
   \(D^{-4}\) poles; search non-elementary affine modifications of
   \(K[s,W,D,D^{-1}]\) that contain the Davenport target algebra.
2. **Beyond length-two Jung coordinates.**  A three-coefficient Newton
   gate now excludes every quadratic-first length-two coordinate in all
   second-shear degrees; degree dominance excludes the other orientation
   and all higher first shears.  The first constant-direction
   quadratic-in-\(U\) suspension also fails: its top equation reduces to
   the centralizer of \(g-\lambda T\), and its middle equation has an
   unavoidable degree-five contradiction.  The next affine-in-\(U\)
   screen is Jung length three; the next nonlinear screen needs
   nonconstant quadratic coefficients or two auxiliaries.  See
   [the independent marked-line opening](DAVENPORT_INDEPENDENT_MARKED_LINE_OPENING.md)
   and [the quadratic auxiliary audit](DAVENPORT_QUADRATIC_U_INCIDENCE.md).
   For the full quadratic vector, the \(U^5\) equation is projective
   developability.  Every affine-linear rank-two direction is impossible;
   the non-origin rank-one direction is also impossible after an integrated
   conic/cubic and component analysis.  Only the through-origin affine
   rank-one PDE and nonlinear developable directions survive.  The
   through-origin PDE is now a Laurent zero-Jacobian equation; all primitive
   monomial common-parameter branches fail.  For nonlinear directions the
   next equation is \(J(q,\phi^2/\lambda)=0\).  The pure Davenport conic
   \(\mathbf b=(1,g,g^2)\) integrates through \(J_2\); its factor-through
   component and its complete \(R=0\) branch fail, while any genuine
   \(R\ne0\) survivor must generate the full degree-seven Davenport field.
   For \(R\in K^*\) and the natural generator \(v=Y\), the quartic-square
   equation has no reduced Davenport-basis root of degree at most four;
   degree four ends in an irreducible quintic over \(K(T)\).  Only basis
   degrees five and six reach the remaining \(J_1\) differential gate for
   general \(c\).  On the \(R=1\) slice, a necessary \((T,z)=(1,0)\)
   \(J_1\)-fiber has unit Gröbner basis, closing all basis degrees.
   See
   [the quadratic-vector developability audit](DAVENPORT_QUADRATIC_VECTOR_DEVELOPABILITY.md)
   [the quadratic-survivor calculation](DAVENPORT_QUADRATIC_SURVIVORS.md),
   and [the constant-normal basis screen](DAVENPORT_CONSTANT_NORMAL_BASIS_SCREEN.md).
3. **Fano incidence class calculation.**  Compute equivariant finite-field
   point counts for the point- and line-marked normalized slices.  Reject a
   slice as soon as either class differs from \(\mathbb L^n\) or has
   nonconstant units.
4. **Boundary-decorated arithmetic check.**  Decide whether the intended
   theorem is good-fiber or all-target indistinguishability.  For the
   stronger form, add omitted-prime cycle data to the construction from the
   beginning.
5. **Stable invariant transfer.**  State in advance whether the absolute
   map will retain the weighted Hessian divisor or use the intrinsic
   decorated normalization boundary.

## 10. Literature

- L. Alpöge, announced three-dimensional Keller counterexample (2026);
  exact map and provenance are summarized in the repository's foundational
  verification.
- A. Gallagher, [*The Jacobian counterexample, explained*](https://doi.org/10.5281/zenodo.21479195)
  (2026), including absolute weighted lifts in every generic degree
  \(N\ge3\).
- T. Shaska,
  [*Graded Keller maps and the Jacobian Conjecture*](https://arxiv.org/abs/2607.20210)
  (2026), for the mixed-sign quotient geometry and boundary escape.
- Z. Jelonek and M. Lasoń,
  [*Quantitative properties of the non-properness set of a polynomial map*](https://arxiv.org/abs/1411.5011),
  for hypersurface and uniruledness constraints on nonproperness sets.
- S. Kaliman and M. Zaidenberg,
  [*Affine modifications and affine hypersurfaces with a very transitive automorphism group*](https://arxiv.org/abs/math/9801076),
  for the affine-modification framework \(A[I/f]\).
- W. Bosma and B. de Smit,
  [*On arithmetically equivalent number fields of small degree*](https://www.math.ru.nl/~bosma/pubs/ANTS2002.pdf),
  for the degree-seven classification and the LaMacchia sign family.
- T. Sunada,
  [*Riemannian coverings and isospectral manifolds*](https://annals.math.princeton.edu/1985/121-1/p04),
  for the almost-conjugate subgroup mechanism.

The literature establishes the group-theoretic and affine-modification
ingredients separately.  The prescribed-monodromy polynomial
algebraization required here remains an open synthesis.
