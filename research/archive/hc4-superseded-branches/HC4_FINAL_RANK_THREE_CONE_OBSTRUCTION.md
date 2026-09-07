# Cone obstruction in the final rank-three `[4]` HC4 stratum

> **Archived route.** The identifier `HC4RSD75` is assigned in the active
> chain to the final Frobenius closure.  This cone-only argument is retained
> for provenance.

## Status

This note continues `HC4RSD71--74` and removes every cone-type gradient-image
geometry from the only remaining relative-nilpotent HC4 branch.

> **Theorem HC4RSD75 — cone images force a linear invariant.**
> Let `F=grad A` have generic Hessian rank three in four variables, and suppose
> the projective closure of
> \[
> Y=\overline{\operatorname{im}F}\subset\mathbb A^4
> \]
> is a cone.  Then the associated Hessian-kernel quasi-translation has a
> constant linear invariant.  Consequently such a packet is already closed by
> `HC4RSD65` and the corrected `HC4RSD70`.
>
> In particular the Piontkowski Gauss-rank-two focal type `(3,1)` cannot occur
> in the genuinely linearly-independent final `[4]` branch.

## 1. Vertex at infinity

If the cone vertex contains a point at infinity which supplies the affine Gauss
ruling direction, then the affine gradient image is cylindrical in a constant
direction.  Equivalently its defining equation is independent of one constant
linear coordinate.  Its normal/kernel vector therefore has a constant linear
relation, and this is exactly the branch closed by `HC4RSD65/70`.

The only apparently new case is therefore a finite cone vertex.

## 2. Translate a finite vertex to the origin

Let `v` be a finite cone vertex and put

\[
\widetilde A=A-v\cdot x,
\qquad
\nabla\widetilde A=F-v.
\]

Translation does not change the Hessian:

\[
\operatorname{Hess}\widetilde A=T.
\]

Because `Y-v` is a cone, its irreducible defining equation may be chosen
homogeneous:

\[
R(F-v)=R(\nabla\widetilde A)=0,                     \tag{2.1}
\]

with `R` homogeneous and irreducible.

Since `rank T=3` generically, the algebraic relations among the four components
of `grad \widetilde A` form a height-one prime ideal.  Hence the irreducible
hypersurface equation `R` generates that ideal and is, in particular, a
nonzero relation of minimum degree.

## 3. Apply the four-variable singular-Hessian relation theorem

De Bondt's Theorem 4.6 in *Quasi-translations and singular Hessians* states the
following in precisely this nonhomogeneous four-variable setting:

> If `h` is a polynomial in four variables and `R` is a homogeneous relation
> of `grad h` of minimum degree, then `R` can be expressed as a polynomial in
> three constant linear forms in its target variables.

Apply this to `h=\widetilde A` and (2.1).  There is therefore a nonzero constant
vector `p` such that

\[
D_pR=p\cdot\nabla R=0.                              \tag{3.1}
\]

The associated primitive Hessian-kernel vector is, up to a scalar invariant,

\[
k=(\nabla R)(F-v).
\]

Equation (3.1) gives immediately

\[
\boxed{p\cdot k=0.}                                 \tag{3.2}
\]

Thus the components of `k` are linearly dependent over the base field.

## 4. Consequence

The final unresolved `[4]` stratum was defined by the opposite condition: the
four components of its associated quasi-translation are linearly independent.
Therefore a cone image is impossible.

Combining this with the previous global reductions leaves only **non-conical,
affinely singular** Gauss-rank-two developable gradient images.  Piontkowski's
rank-two classification reduces the remaining geometry to the non-cone focal
types `(1,1)`, `(1,2)`, `(2,1)`, and `(2,2)` (with any trivial cone factors
already removed).

## External input

The algebraic input is Michiel de Bondt, *Quasi-translations and singular
Hessians*, especially Theorem 4.6: a homogeneous minimum-degree relation of a
four-variable polynomial gradient depends on only three linear forms.  The
projective focal nomenclature is from Jens Piontkowski, *Developable Varieties
of Gauss rank 2*.
