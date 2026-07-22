# Audit of the “Six Questions” take-aways

This is a provenance and research-scope audit of a supplied summary attributed
to a recent paper of Bartosz Naskręcki.  As of 22 July 2026, the repository
has not located a complete public source, bibliographic identifier, version,
or stable URL for that paper.  The summary's claims retain the status
**external claims under source audit**.  Some are now also independent
repository theorems because the necessary algebra has been reconstructed and
checked below; independent reproduction does not settle bibliographic
provenance or constitute external review.

The author's
[public papers list](https://bnaskrecki.faculty.wmi.amu.edu.pl/doku.php/papers)
did not list this item when checked on that date.  This is a search-status
statement, not evidence that the manuscript does not exist.

The summary is valuable because it separates several finite questions from
the already-understood foundational construction.  This note records which
ones are closed here, which have now been proved, and which require the
paper's exact input data.

## 1. Foundational factorization: closed

The compact
[normalized factorization model](../verified/NORMALIZED_FACTORIZATION_MODEL.md)
now gives the complete chain

\[
 (L,Q)\longmapsto(LQ,\operatorname{Res}(L,Q)),\qquad
 X\simeq\mathbb A^3,\qquad
 F_{\rm orig}=B\circ G\circ A.
\]

It includes the global polynomial inverse across `a=0`, the tangent-space
proof of étaleness, `det DTheta=-Res^2`, `det DG=-1`, and the exact linear
comparison.  Further coordinate checks should be added only when they test a
new reusable mechanism.

## 2. Vertical derivation plus slice: proved for the cubic

On

\[
 X=\{a^2e-abd+b^2c=1,\ ad+bc=1\},
\]

the saturated primitive vertical derivation is

\[
 \boxed{D=a^2\partial_c-ab\partial_d-2b^2\partial_e.}
\]

It is locally nilpotent, preserves both equations, and has

\[
 y=2bd-ae\in\ker D,\qquad
 z=2d^2+ce+6bd^2+3bce-\frac92e,\qquad D(z)=1
\]

in `k[X]`.  Consequently

\[
 k[X]=\ker(D)[z]=k[a,y,z].
\]

This is the requested “vertical derivation plus polynomial slice” proof of
the cylinder.  The formulas and ideal reductions are part of the exact
factorization verifier.

With `D_0=D/2`, the action has the intrinsic factor form

\[
 Q\longmapsto Q+\frac{t}{2}L(aT-2bS).
\]

It preserves the resultant because the residual factor is unchanged modulo
`L`, and it preserves the selected cubic coefficient by direct
multiplication.  The normalization equation gives `b-ay=1`, and integrating

\[
 D_0(-e)=b^2,\qquad
 D_0(4yd)=-2aby,\qquad
 D_0(2cy^2)=a^2y^2
\]

produces the alternative slice

\[
 \widetilde z=4yd+2cy^2-e,\qquad D_0(\widetilde z)=1.
\]

It satisfies `tilde z=2z` in `k[X]`.  Thus the “integrated Bézout identity”
and the existing inverse coordinate are the same construction with different
normalizations.

The primitive saturated vertical generator is unique up to a nonzero scalar
after fixing the projection.  It is not the only locally nilpotent
derivation: `f(a,y)D` is locally nilpotent for every invariant `f`.  Only
invariant units retain a global everywhere-free slice action; nonconstant
factors add fixed fibers.

For higher weighted marked-root maps, the affine source is already presented
as `A^3_(x,y,z)`.  What remains open is not the existence of an arbitrary
locally nilpotent derivation on that affine space, but whether normalization
and reconstruction select a **canonical vertical derivation and slice** that
can be expressed intrinsically on the marked-root incidence.  Merely
transporting `partial_z` through the known source coordinates would be
tautological and is not counted as progress.

## 3. Weighted invariant-coordinate reduction: proved

For source weights `(1,-1,-k)`, put

\[
 u=xy,\qquad v=x^kz.
\]

If `F_i=x^(w_i)A_i(u,v)` on `x!=0`, then the exact chain-rule identity is

\[
 \det JF=x^{w_1+w_2+w_3+k}
 \det
 \begin{pmatrix}
 w_1A_1&(A_1)_u&(A_1)_v\\
 w_2A_2&(A_2)_u&(A_2)_v\\
 w_3A_3&(A_3)_u&(A_3)_v
 \end{pmatrix}.
\]

In the balanced case `sum(w_i)=-k`, the three-variable Keller equation is
exactly the displayed two-variable determinant equation.  The proof,
polynomial-extension scope, and exact checks are in the
[weighted invariant-coordinate lemma](WEIGHTED_INVARIANT_JACOBIAN_REDUCTION.md).
This closes the reusable differential reduction without relying on the later
sixteen-monomial Gröbner computation.

## 4. Weighted-support rigidity: independently reproduced

The fuller summary supplies enough structural conditions to reconstruct the
support without copying an external coefficient list.  Source weights
`(1,-1,-2)`, output weights `(-2,-1,1)`, degree bounds `(7,6,4)`, and
linearity in `z` give exactly `7+6+3=16` monomials.

After fixing the linear part, the general determinant lemma produces sixteen
coefficient equations in thirteen unknowns.  On `pq!=0`, the effective
two-dimensional diagonal action uniquely normalizes

\[
 p=[x^2y]F_3=-3,\qquad q=[x^3z]F_3=-1.
\]

An exact Gröbner computation then gives the complete triangular basis,
including `(B_01-3)^2`, and proves

\[
 \mathcal O_{\rm ansatz/gauge}
 \simeq\mathbb Q[\varepsilon]/(\varepsilon^2).
\]

The theorem, all eleven triangular relations, and the precise gauge scope are
in the
[foundational weighted coefficient scheme](FOUNDATIONAL_WEIGHTED_COEFFICIENT_SCHEME.md).
This is classification inside the displayed weighted support, not
classification of all degree-profile `(7,6,4)` Keller maps.

## 5. Dual-number deformation: explicit and polynomial-orbit trivial

The Gröbner basis now gives an explicit

\[
 F_\varepsilon=F+\varepsilon H,\qquad\varepsilon^2=0.
\]

The exact checker proves:

- the linearized Jacobian determinant vanishes;
- `det J(F+tH)=-2+t^2 Omega` for an explicit nonzero `Omega`;
- the full normalized ansatz has no reduced second-order lift;
- `H` is not tangent to the affine source--target left--right orbit;
- nevertheless, the
  [polynomial-orbit calculation](CUBIC_DUAL_NUMBER_ORBIT_TANGENCY.md) gives
  `H=DF\,V` for the polynomial divergence-free field
  `V=(DF)^(-1)H`, and an explicit decomposition into 133 locally nilpotent
  shears realizes it as the tangent of a reduced polynomial source-orbit
  curve.

Thus the claimed local algebra

\[
 \mathbb Q[\varepsilon]/(\varepsilon^2)
\]

is now a proved repository result.  Its nilpotent direction is therefore an
artifact of the normalized finite-support slice, not an infinitesimal class
modulo the full polynomial orbit.  The
[formal-orbit theorem](FORMAL_ORBIT_TRIVIALITY.md) strengthens this: every
Artin deformation of a Keller map is uniquely source-trivial, and in
dimension three every finite determinant-one jet has a reduced polynomial
special-automorphism representative.  Thus polynomial-orbit tangency and
finite-order/formal triviality are closed already before stabilization;
stabilization adds nothing at those levels.

What remains is the formal-to-algebraic gap.  Formal triviality does not
produce one regular polynomial automorphism family over a reduced base, and
therefore does not decide exact reduced polynomial left--right equivalence,
with or without stabilization.  The literal line `F+tH` is not a Keller
family beyond first order, and the normalized ansatz has no reduced
second-order lift, so that global question cannot be inferred from this
dual-number point alone.

## 6. The `pq=0` boundary: primary decomposition complete

The
[coefficient-scheme computation](FOUNDATIONAL_WEIGHTED_COEFFICIENT_SCHEME.md)
generates the full unsaturated thirteen-variable Keller ideal explicitly by
the sixteen coefficients of the invariant determinant.  Both punctured
coordinate divisors `p=0,q!=0` and `q=0,p!=0` are empty.  On `p=q=0`, exact
primary decomposition gives two nonreduced primary components with no
embedded associated prime.  Their radicals are the affine three-spaces of
triangular automorphisms

\[
 (A,B,C)=(v+\alpha _2u^2+\alpha _3u^3+\alpha _4u^4,u,2)
\]

and

\[
 (A,B,C)=(v+\alpha _2u^2+\alpha _3u^3,
          u+\lambda(v+\alpha _2u^2+\alpha _3u^3),2).
\]

They meet in the affine plane obtained by setting `lambda=alpha_4=0`, and
every reduced point on either component is a polynomial automorphism.  The
previously known line `(A,B)=(v,u+lambda v)` is only a slice of the second
component.

The reduced global scheme is also determined.  It is the union of these two
boundary components and the degree-ten toric closure `T` of the reduced
`pq!=0` orbit.  The toric component attaches along exactly the two explicit
lines `A=v+lambda u^2, B=u` and `A=v, B=u+lambda v`.  The only unresolved
boundary information is the nilpotent gluing of the generically doubled
toric component to the two nonreduced boundary components along those known
lines; no reduced component or reduced attachment remains to be found.

## 7. Hyperplane classification: already proved, with stronger obstructions

For cubic forms, hyperplanes have the three contact types

\[
 (1,1,1),\qquad(2,1),\qquad(3).
\]

The active
[marked-point dimension barrier](MARKED_POINT_DIMENSION_BARRIER.md)
and its exact factorization calculation prove:

\[
\begin{array}{c|c}
\text{contact type}&[X_h]\\ \hline
(1,1,1)&L^3-L\\
(2,1)&L^3\\
(3)&L^3-L^2.
\end{array}
\]

Thus the two unsuccessful types are excluded by their Grothendieck classes
and good-reduction point counts, which are stronger than Euler
characteristic.  The `(2,1)` orbit is `A^3` by the global polynomial model.
The `PGL_2`-equivariance and two normalization scalars give explicit
isomorphisms between every two hyperplanes in that orbit; the proof is now
spelled out in the active classification note.

Their finer isomorphism types, Makar--Limanov invariants, and fundamental
groups remain optional secondary questions because they do not affect the
affine-space classification.

## 8. Five different equivalence actions

The actions in the summary do not form one undifferentiated group.

\[
\begin{array}{ccl}
\mathbb G_m^{\rm factor}
&:&(L,Q)\mapsto(\lambda L,\lambda^{-1}Q),\\
\mathbb G_m^{\rm source}
&:&(x,y,z)\mapsto(\lambda x,\lambda^{-1}y,\lambda^{-2}z),\\
G_{\rm diag}^{\rm ansatz}
&\subset&
\operatorname{Aut}_{\rm poly}(\mathbb A^3)_{\rm src}
\times\operatorname{Aut}_{\rm poly}(\mathbb A^3)_{\rm tgt}\\
&\subset&\text{stable polynomial left--right groupoid}.
\end{array}
\]

The factor gauge is a redundancy in the factor-pair presentation and is
removed when passing to the normalized source.  The source torus is the
one-dimensional kernel of the three-parameter diagonal coordinate changes on
the coefficient support.  The effective quotient producing the dual-number
point is only `G_diag^ansatz`.  The affine-orbit calculation enlarges this to
all affine source and target transformations and still does not contain the
nilpotent tangent.  Those two quotient calculations alone do not control the
full polynomial or stable orbit.  The separate polynomial-orbit calculation
does: the tangent is already source-orbit trivial, and the formal-orbit
theorem trivializes every Artin deformation before stabilization.  Neither
statement crosses the remaining formal-to-algebraic gap for one exact family
over a reduced base.

## 9. Priority

The degreewise decorated-normalization theorem proves an `(N-3)`-dimensional
family of stable classes for every `N>=4`.  The formerly open deformation
question

\[
 \boxed{\text{Is the decorated-normalization map generically unramified?}}
\]

is now resolved affirmatively by the Fitting divisor with the marks zero and
infinity.  If a normalized tangent deformation is trivial modulo scaling of
the normalization line, then
`dot H''=alpha H''+beta rH'''`.  Integration and the normalized endpoint
conditions force first `beta=0` and then `alpha=0`.  The same argument,
together with exact rerooting, shows that the map is generically etale of
degree `N-2` onto the normalization of its reduced image.  Node pairing and
the conductor are unnecessary for this differential statement.

The remaining external-source task is bibliographic: locate the complete
manuscript and compare its conventions and certificate with this independent
reconstruction.  The separate internal problem is cross-stratum generator
rigidity for the full intrinsic affine cover; that refinement could distinguish
the finite rerootings but is not needed for generic unramifiedness.
