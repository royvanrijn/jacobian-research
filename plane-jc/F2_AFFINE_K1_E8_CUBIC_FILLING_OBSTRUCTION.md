# Homological filling obstruction for the cubic E8 equality row

> **Status.**  Exact topological exclusion.  Assume the F2 affine
> nonproperness set has one component, the E8 cusp `P^5-Q^3=0`, and the
> geometric degree is six.  The logarithmic budget and cusp-group atlas
> leave one higher-inertia possibility: the natural degree-six `A_6` action
> with meridian type `3+1+1+1`.  The corresponding six-sheet cover of the
> cusp complement has first homology `Z^4`.  Filling the three peripheral
> components fixed by the meridian—exactly the three affine pullback
> components—leaves first homology `Z`.  But filling all affine pullback
> components recovers the affine source `A^2`, whose first homology is zero.
> Therefore the cubic equality row is impossible.  Combined with the
> simple-inertia exclusion and the degree-six permutation audit, this closes
> every one-component E8 completion at geometric degree six.  It does not
> address larger geometric degree, additional nonproperness components, or
> non-E8 target curves.

The lifted cellular boundary matrices and their integral Smith forms are
replayed by
[`verify_f2_affine_k1_e8_cubic_filling_obstruction.py`](../scripts/verify_f2_affine_k1_e8_cubic_filling_obstruction.py).

## 1. Complement and cover

Let

\[
 C:\quad P^5-Q^3=0,
 \qquad U=\mathbb A^2\setminus C.                \tag{1.1}
\]

Weighted radial retraction identifies `U` with the complement of the
`(3,5)` torus knot up to the harmless radial factor.  Thus

\[
 \pi_1(U)=G=\langle a,b\mid a^3=b^5\rangle,      \tag{1.2}
\]

with geometric meridian

\[
 m=a^{-1}b^2.                                    \tag{1.3}
\]

If `F:A^2->A^2` is Keller and `C` is its entire nonproperness set, then

\[
 V=F^{-1}(U)\longrightarrow U                   \tag{1.4}
\]

is a connected finite étale cover of degree `d`: away from the
nonproperness set the quasi-finite map is proper, hence finite, and the
Jacobian condition makes it étale.

At `d=6`, the unique cubic-inertia action is

\[
\begin{aligned}
 A&=(1\ 2\ 3)(4\ 5\ 6),\\
 B&=(2\ 4\ 3\ 5\ 6),\\
 M&=A^{-1}B^2=(1\ 3\ 5).
\end{aligned}                                    \tag{1.5}
\]

Its image is the natural `A_6`, so the cover is connected.  The meridian
has one three-cycle and three fixed points.

## 2. Which peripheral components are filled

Peripheral orbits classify the components of the inverse image of a small
boundary torus around `C`.  Since the preferred longitude is trivial in
(1.5), the orbit data are

\[
 (e,f)=(3,1),\qquad(1,1),\qquad(1,1),\qquad(1,1). \tag{2.1}
\]

The three fixed orbits are the affine components of `F^{-1}(C)`.  The
three-cycle is the nonproper boundary escape and is absent from the affine
source.

Put

\[
 A_F=F^{-1}(C)\subset\mathbb A^2.                \tag{2.2}
\]

Then `V=A^2-A_F`.  Adding one irreducible affine curve component to its
complement kills the corresponding meridian normally.  This follows
directly from van Kampen applied to a tubular neighborhood of its smooth
locus; singular points add cells of dimension at least two and do not
restore first homology.  Since the Keller map is locally biholomorphic,
distinct fixed sheets give distinct local pullback components.  Therefore

\[
 H_1(\mathbb A^2;\mathbb Z)
 =H_1(V;\mathbb Z)/\langle
 \text{the three fixed-sheet meridians}\rangle.  \tag{2.3}
\]

The left side is zero.

## 3. Lifted cellular complex

Use the standard two-complex of (1.2): one vertex, two oriented one-cells,
and one two-cell with attaching word

\[
 a^3b^{-5}.                                      \tag{3.1}
\]

The degree-six cover has six vertices and twelve oriented one-cells: one
`a`-edge and one `b`-edge from every sheet.  Its connected one-skeleton has
cycle rank

\[
 12-6+1=7.                                       \tag{3.2}
\]

Lift (3.1) from all six vertices.  With respect to the seven fundamental
cycles determined by one spanning tree, the six lifted relators have Smith
diagonal

\[
 (1,1,1,0,0,0).                                 \tag{3.3}
\]

Hence

\[
 H_1(V;\mathbb Z)\simeq\mathbb Z^4.             \tag{3.4}
\]

Now append the three closed lifted meridian paths based at the three fixed
sheets.  The resulting `7 x 9` relation matrix has Smith diagonal

\[
 \boxed{(1,1,1,1,1,1,0).}                       \tag{3.5}
\]

There is no torsion, but one free class remains:

\[
 \boxed{
 H_1(V;\mathbb Z)/\langle
 \text{fixed meridians}\rangle\simeq\mathbb Z.} \tag{3.6}
\]

The checker performs the calculation in both left- and right-action path
conventions; they give the same Smith form.

## 4. Contradiction and scope

Equations (2.3) and (3.6) give

\[
 0=H_1(\mathbb A^2;\mathbb Z)\simeq\mathbb Z,
\]

a contradiction.  Therefore

\[
 \boxed{
 \text{the degree-six cubic-inertia E8 equality row does not occur.}} \tag{4.1}
\]

This obstruction is independent of the contracted-chain self-intersections
and explains why the locally consistent Smith packet cannot be globally
completed.  The local normal form remains a valid local model; it simply
cannot be the boundary of an affine-plane cover with the prescribed
monodromy and only one nonproperness component.

The argument uses the one-component hypothesis in identifying the base of
the finite cover with `A^2-C`.  If another nonproperness component is
present, its meridians and their lifted fillings add further relations and
must be included in the cellular matrix.  Likewise, actions of larger
degree require their own permutation representation.  No claim about those
rows is made here.

## Literature inputs

- S. Yu. Orevkov,
  [*On three-sheeted polynomial mappings of C2*](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  for the finite-cover compactification, Euler-multiplicity method, and
  linear affine-dicritical chains.
- Nguyen Van Chau,
  [*Non-proper value set and the Jacobian condition*](https://arxiv.org/abs/math/0305088),
  for the one-point-at-infinity structure of the nonproperness curve.

The group presentation, lifted CW complex, Smith reduction, and van Kampen
filling argument are written out above rather than imported as a black box.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_e8_cubic_filling_obstruction.py
```
