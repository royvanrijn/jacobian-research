# Generic exclusion on the smooth-quartic polar conic

## Status

This note proves `HC4NHM22`.  It continues `HC4NHM20` by imposing the full
reciprocal-Hessian equations at the generic point of the polar conic
\(\Delta=0\).  The set-theoretic support is the determinant-zero boundary,
so the generic point of the first exceptional divisor left by `HC4NHM16` is
empty.  The theorem does not classify every lower denominator stratum or the
non-generic two-line fibers.

Replay the exact calculation with

```bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_polar_conic_generic_gate.py
```

The checker constructs the same 81 equations as `HC4NHM16`.  Singular 4.4.1
performs the staged calculation over a rational function field.

## 1. A rational parametrization of the polar

Retain

\[
\begin{aligned}
\Delta={}&(3p^2-qr)\tau^5+(9pr-q^2)\tau^4
 +(18r^2-6pq)\tau^3\\
&+(18p^2-6qr)\tau^2+(9pr-q^2)\tau+(3r^2-pq).
\end{aligned}
\tag{1.1}
\]

The point

\[
 [p:q:r]=[1:3:1]
\tag{1.2}
\]

lies on every fiber of (1.1).  Geometrically it is the quadratic
\(H=s^2-st+t^2\), which contains the two nonrational roots of
\(s^3+t^3\).  In the affine chart \(p=1\), take the line

\[
 q=3+\theta,\qquad r=1+m\theta.
\tag{1.3}
\]

Direct substitution gives

\[
 \Delta(1,3+\theta,1+m\theta)
   =\theta(D\theta-N),
\tag{1.4}
\]

where

\[
\begin{aligned}
D={}&18m^2\tau^3+3m^2-m\tau^5-6m\tau^2-\tau^4-\tau,\\
N={}&3m\tau^5-9m\tau^4-36m\tau^3+18m\tau^2-9m\tau-6m\\
 &+\tau^5+6\tau^4+6\tau^3+6\tau^2+6\tau+1.
\end{aligned}
\tag{1.5}
\]

Thus the second intersection is \(\theta=N/D\).  Restoring the homogeneous
scale gives the generic parametrization

\[
 p=c,\qquad q=c(3+N/D),\qquad r=c(1+mN/D).
\tag{1.6}
\]

It is birational away from the base point and the usual line-parametrization
denominators.  In particular, calculation over
\(\mathbf Q(\tau,m,c)\) in (1.6) is calculation at the generic point of the
smooth polar conic, not a collection of specialized fibers.

## 2. The reduced reciprocal-Hessian ideal

Use the complete 81-equation system of `HC4NHM16`, with active deformation
coordinates

\[
 b_0,\ldots,b_{14},u,v,w
\tag{2.1}
\]

and parameters

\[
 b_{15},b_{16},b_{17},c,m,\sigma,\tau.
\tag{2.2}
\]

After (1.6), the ten linear equations retain rank ten.  The first relation
is

\[
 u+\tau v=0.
\tag{2.3}
\]

The other nine eliminate

\[
 b_{12},b_{10},b_9,b_7,b_6,b_4,b_3,b_1,b_0
\tag{2.4}
\]

over the parameter function field.  Eight selected coefficients of the
remaining reciprocal system suffice.  Their reduced standard basis has 11
elements; recombination with the linear basis has 21 elements and dimension
zero.  Most importantly, for every coordinate \(a\) in (2.1), the exact
membership calculation gives

\[
 \boxed{a^6\in I_{\mathrm{polar}}.}
\tag{2.5}
\]

Consequently

\[
 \sqrt{I_{\mathrm{polar}}}
 \supseteq
 (b_0,\ldots,b_{14},u,v,w).
\tag{2.6}
\]

This is a set-theoretic statement, which is exactly what is needed for the
exclusion.  It does not assert that the generic polar fiber is reduced.

## 3. Determinant-zero support

On (2.6), the reciprocal quadratic matrix is

\[
A=\begin{pmatrix}
0&0&-y^2\\
0&0&x^2\\
-y^2&x^2&p x^2+qxy+r y^2
       +z(b_{15}x+b_{16}y+b_{17}z)
\end{pmatrix}.
\tag{3.1}
\]

Its determinant is identically zero.  Therefore the prospective quotient
\(Q=\det(A)/(z\ell)\) vanishes and cannot be a smooth nonzero quartic.

> **Theorem `HC4NHM22` -- Generic polar-conic exclusion.**  In the
> basepoint-free squarefree-line row \(d_0=(x^2,y^2,0)\), the generic point
> of the first exceptional divisor \(\Delta=0\) admits no solution of the
> necessary reciprocal-Hessian equations with nonzero quartic quotient.
> Indeed, the polar conic has the rational parametrization (1.6), and over
> its parameter function field the sixth power of every active deformation
> coordinate belongs to the reciprocal ideal.  Its set-theoretic support is
> therefore the determinant-zero matrix (3.1).

This removes the generic smooth-conic component requested after `HC4NHM20`.
What remains on this first divisor is confined to:

1. the line-parametrization denominator and secondary-basis denominator
   strata;
2. after the Fermat-symmetry transport `HC4NHM23`, the two line components
   in the two remaining six-point slope normal forms; the size-three orbit
   \(\tau^3=-1\) is covered generically and on its first secondary strata;
3. intersections of those strata with further hidden denominators.

The complementary residual-line chart and the other reciprocal boundary
types remain separate.

See
[`HC4_SMOOTH_QUARTIC_FERMAT_SYMMETRY_ORBITS.md`](HC4_SMOOTH_QUARTIC_FERMAT_SYMMETRY_ORBITS.md)
for the exact covariance and orbit calculation.
