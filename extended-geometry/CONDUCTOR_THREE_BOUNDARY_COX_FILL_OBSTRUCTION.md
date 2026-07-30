# The rational conductor three-boundary Cox-fill obstruction

> **Status.** The smallest symmetric three-boundary escape from the
> conductor one-chart theorem is closed exactly.  For both the nodal and
> cuspidal rational conductor algebras, the Cox fill \(xyz=p\) makes the
> incidence coordinates polynomial and removes the localization-unit
> obstruction.  It does not produce a new Keller source.  The descended
> dualizing form pulls back with one conductor pole.  Independently, the
> normalized nodal fill is smooth but has a non-affine Hodge--Deligne
> polynomial, while the normalized cuspidal fill is normal but singular.
> Thus neither modified source is affine three-space, even after ordinary
> polynomial stabilization.

This is a scoped conductor obstruction theorem, not an obstruction to every
three-boundary construction.  It tests the canonical symmetric Cox fill
suggested by the failure of the separated localization.

Work over a characteristic-zero field \(k\).

## 1. The node and cusp incidence algebras

Write the nodal conductor algebra as

\[
 A_{\rm n}
 =k[p,q]/(q^2-pq-p^3)
 \subset k[t],
 \qquad
 p=t(t-1),\quad q=t^2(t-1).                         \tag{1.1}
\]

Its normalization parameter satisfies

\[
 T^2-T-p=0,\qquad pT-q=0.                            \tag{1.2}
\]

Away from \(p=0\), reconstruction is \(T=q/p\).  At \(p=q=0\),
equation (1.2) has the two points \(T=0,1\), which are the two branches
identified by the node.

For the cusp use

\[
 A_{\rm c}
 =k[p,q]/(q^2-p^3)
 =k[t^2,t^3]\subset k[t],
 \qquad p=t^2,\quad q=t^3.                            \tag{1.3}
\]

Here the normalization equations are

\[
 T^2-p=0,\qquad pT-q=0,                               \tag{1.4}
\]

and the conductor fiber is \(T^2=0\).  It records the missing first
normalization jet rather than two distinct branches.

In both cases the conductor in \(A\) is \((p,q)\), and its pullback to
\(k[t]\) is \(p\,k[t]\).  Explicitly, the normalization quotient is cyclic
on the class \(\bar T\), with

\[
 (k[t]/A)\simeq A/(p,q)\,\bar T.                      \tag{1.5}
\]

Indeed \(p\bar T=0\); also \(qT=q+p^2\) for the node and \(qT=p^2\) for
the cusp.  Reduction to the conductor fiber shows that there is no larger
annihilator.

## 2. The minimal three-boundary fill

Replace the separated open \(k[t,1/p]\) by

\[
 \widetilde X_c
 =
 \{xyz=c(t)\}
 \subset\mathbb A^4_{x,y,z,t},                       \tag{2.1}
\]

where

\[
 c_{\rm n}(t)=t(t-1),\qquad c_{\rm c}(t)=t^2.         \tag{2.2}
\]

The three factors \(x,y,z\) are the three boundary directions.  Descending
the normalization gives

\[
 X_{\rm n}
 =
 \{q^2-xyz\,q-(xyz)^3=0\}
 \subset\mathbb A^4_{x,y,z,q},                       \tag{2.3}
\]

and

\[
 X_{\rm c}
 =
 \{q^2-(xyz)^3=0\}
 \subset\mathbb A^4_{x,y,z,q}.                       \tag{2.4}
\]

The finite maps

\[
 \pi_c:\widetilde X_c\longrightarrow X_c             \tag{2.5}
\]

send \(p\) to \(xyz\) and \(q\) to \(t\,c(t)\).  Every displayed coordinate
is polynomial.  On \(xyz\ne0\), the inverse reconstruction is

\[
 \boxed{t=\frac q{xyz}.}                              \tag{2.6}
\]

Thus this ambient coupling genuinely escapes the elementary divisibility
contradiction in the one-chart theorem: the pole belongs to an integral
normalization coordinate, while all incidence outputs descend
polynomially.

The conductor of (2.5) is the extension of \((p,q)\); on the normalization
its divisor is

\[
 p=xyz.                                               \tag{2.7}
\]

For the node it has order one on each of the six components
\(D_{x,0},D_{y,0},D_{z,0},D_{x,1},D_{y,1},D_{z,1}\).
For the cusp it has order two on the three generic boundary components
\(D_x,D_y,D_z\).

## 3. The determinant valuation ledger

Let \(F_c(x,y,z,q)\) denote the equation (2.3) or (2.4).  On the descended
hypersurface use the dualizing residue

\[
 \Omega_{\rm desc}
 =
 \frac{dx\wedge dy\wedge dz}{\partial F_c/\partial q}.
                                                               \tag{3.1}
\]

On (2.1), a residue generator is

\[
 \Omega_{\rm norm}
 =
 \frac{dx\wedge dy\wedge dt}{xy}.                    \tag{3.2}
\]

Holding \(x,y\) fixed gives

\[
 dz=\frac{c'(t)}{xy}\,dt+\text{terms in }dx,dy.       \tag{3.3}
\]

In both conductor cases,

\[
 \left.\frac{\partial F_c}{\partial q}\right|_{\widetilde X_c}
 =c(t)c'(t).                                         \tag{3.4}
\]

Consequently

\[
 \boxed{
 \pi_c^*\Omega_{\rm desc}
 =\frac{\Omega_{\rm norm}}{c(t)}
 =\frac{\Omega_{\rm norm}}{xyz}.
 }                                                    \tag{3.5}
\]

Polynomiality has therefore succeeded, but the determinant ledger retains
exactly one inverse conductor.  In valuation form the right side has pole
order one along every nodal boundary component and pole order two along
every cuspidal boundary component.

This is the new obstruction exposed by the three-boundary coupling.  The
one-chart unit has disappeared; its character has moved into the
normalization/dualizing comparison.

## 4. Affine-space recognition

The nodal normalized fill is smooth.  Indeed, a singular point of

\[
 xyz-t(t-1)=0
\]

would have \(yz=xz=xy=0\) and \(2t-1=0\).  But \(t=1/2\) makes
\(xyz=-1/4\), so none of \(x,y,z\) can vanish.

Its Grothendieck class is obtained by separating the two zero fibers of
\(c_{\rm n}\):

\[
\begin{aligned}
[\widetilde X_{\rm n}]
&=(\mathbb L-2)(\mathbb L-1)^2
 +2\bigl(\mathbb L^3-(\mathbb L-1)^3\bigr)\\
&=\boxed{\mathbb L^3+2\mathbb L^2-\mathbb L}.         \tag{4.1}
\end{aligned}
\]

After base change to \(\mathbb C\), its compactly supported
Hodge--Deligne polynomial is

\[
 (uv)^3+2(uv)^2-uv\ne(uv)^3.                          \tag{4.2}
\]

Hence

\[
 \widetilde X_{\rm n}\not\simeq\mathbb A^3.           \tag{4.3}
\]

Multiplication by \((uv)^r\) preserves the nonzero defect, so ordinary
polynomial stabilization cannot repair it.

For the cusp,

\[
[\widetilde X_{\rm c}]
 =(\mathbb L-1)^3+
 \bigl(\mathbb L^3-(\mathbb L-1)^3\bigr)
 =\mathbb L^3.                                       \tag{4.4}
\]

The class alone does not recognize affine space.  The exact singular locus
is

\[
 V(t,y,z)\ \cup\ V(t,x,z)\ \cup\ V(t,x,y).            \tag{4.5}
\]

It has codimension two, so the irreducible hypersurface
\(xyz=t^2\) is normal by Serre's criterion, but it is not smooth and hence
is not affine three-space.  The singular locus persists after adjoining
polynomial variables.

## 5. Obstruction theorem

> **Theorem.** Let \(A\subset k[t]\) be the rational node
> \(k+t(t-1)k[t]\) or the rational cusp \(k[t^2,t^3]\), with conductor
> generator \(p=t(t-1)\) or \(p=t^2\).  Form the symmetric three-boundary
> Cox fill by imposing \(xyz=p\), and descend it through the conductor
> incidence algebra.  Then:
>
> 1. the marked-root equations and reconstruction are exact;
> 2. all descended incidence coordinates are polynomial;
> 3. the normalization/dualizing ledger retains the factor \(p^{-1}\);
> 4. in the nodal case the normalized fill is smooth but not stably affine
>    space by its Hodge--Deligne polynomial; and
> 5. in the cuspidal case the normalized fill is normal but singular, with
>    the three-axis locus (4.5).
>
> Therefore this symmetric three-boundary conductor construction cannot be
> the modified affine-space source of a polynomial Keller family.

The scope is important.  The theorem does not exclude:

1. an asymmetric sequence of affine modifications;
2. a three-boundary relation not of Cox-product form \(xyz=p\); or
3. a distributed source/target ledger in which an additional target
   conductor factor cancels (3.5).

The next viable search should alter one of those three features rather than
increase the degree of the same symmetric product.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_conductor_three_boundary_cox_fill.py
```

The driver invokes Singular 4, verifies the two elimination ideals,
smoothness of the nodal fill, and the exact cuspidal singular locus.  It
then checks the marked-root reconstruction, polynomial descent,
dualizing-form pullback, boundary valuations, and motivic affine-space
tests, and writes
`artifacts/generated-results/conductor_three_boundary_cox_fill.json`.
