# The final two smooth-quartic line-fiber normal forms

## Status

This note proves `HC4NHM24`.  It continues `HC4NHM20` and `HC4NHM23` at
the generic points of the two six-point Fermat-symmetry orbits left open by
the orbit reduction.  One irreducible quartic coefficient field treats both
orbits simultaneously.  On each of the two polar-line components, an exact
eleven-row certificate forces every active reciprocal-Hessian deformation
coordinate to vanish.  The surviving boundary matrix has zero determinant.

Replay the characteristic-zero certificate with

```bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_final_line_orbits.py
```

This is a generic-component theorem over an exact algebraic function field.
It is not inferred from the auxiliary finite-field search that located the
small witness.  The lower linear-pivot and eleven-row determinant strata are
not claimed here.

For repository maintenance, `--audit-existing-only` verifies the exact hash
of the imported 81-equation builder without constructing the exact field or
coefficient matrices.  It is a provenance check, not a replacement proof.

## 1. A quartic field containing both orbit types

Recall the quadratic direction from `HC4NHM20`,

\[
K_\tau=(\tau^5+6\tau^2)s^2
 +3(\tau^4+\tau)st+(6\tau^3+1)t^2.
\tag{1.1}
\]

Fix the Fermat root \([s:t]=[-1:1]\).  Its vanishing polynomial is

\[
K_\tau(-1,1)=(\tau+1)f_4(\tau),
\qquad
f_4=\tau^4-4\tau^3+10\tau^2-4\tau+1.
\tag{1.2}
\]

The root \(\tau=-1\) belongs to the size-three orbit already closed by
`HC4NHM17` and `HC4NHM23`.  The four roots of \(f_4\) meet the other two
orbits.  Indeed,

\[
f_4=
\bigl(\tau^2-(2+2i)\tau+1\bigr)
\bigl(\tau^2-(2-2i)\tau+1\bigr).
\tag{1.3}
\]

If \(a=\tau+\tau^{-1}\), then
\(\tau^3+\tau^{-3}=a^3-3a\).  The two factors in (1.3) therefore give

\[
\tau^3+\tau^{-3}=-22+10i,
\qquad
\tau^3+\tau^{-3}=-22-10i,
\tag{1.4}
\]

which are exactly the two quotient values in `HC4NHM23`.  Since \(f_4\) is
irreducible over \(\mathbf Q\), computations in

\[
E=\mathbf Q[\tau]/(f_4)
\tag{1.5}
\]

prove the identities for all four conjugate slopes and hence for both
six-point orbits.

## 2. The two polar-line components

Modulo \(f_4\), the visible pivot factors as

\[
\Delta=(3p+q+3r)\,Q_\tau,
\tag{2.1}
\]

where

\[
\begin{aligned}
Q_\tau={}&(6\tau^3-30\tau^2+15\tau-4)p\\
&+(-4\tau^3+10\tau^2-5\tau+1)q
 +(6\tau^3+1)r.
\end{aligned}
\tag{2.2}
\]

Thus the exceptional polar conic is the union of the resultant line
\(3p+q+3r=0\) and the residual-polar line \(Q_\tau=0\).  Write

\[
A=6\tau^3-30\tau^2+15\tau-4,
\quad B=-4\tau^3+10\tau^2-5\tau+1,
\quad C=6\tau^3+1.
\tag{2.3}
\]

Generic affine coordinates on the two components are

\[
(p,q,r)=c(3,3m,-3-m)
\tag{2.4}
\]

and

\[
(p,q,r)=c(C,Cm,-A-Bm),
\tag{2.5}
\]

respectively.  Both calculations below take place over

\[
L=E(m,c,\sigma,b_{15},b_{16},b_{17}).
\tag{2.6}
\]

## 3. The eleven-row certificate

Specialize the complete 81 coefficient equations of `HC4NHM16` by either
(2.4) or (2.5).  In both cases the ten linear equations have pivots

\[
b_0,b_1,b_3,b_4,b_6,b_7,b_9,b_{10},b_{12},u.
\tag{3.1}
\]

The quotient by this linear layer is the polynomial ring over \(L\) in

\[
b_2,b_5,b_8,b_{11},b_{13},b_{14},v,w.
\tag{3.2}
\]

Number the original 81 equations in the deterministic order produced by
`build_equations("squarefree-line", False)`.  Reduce the eleven rows

\[
28,34,46,49,52,56,62,64,67,71,76
\tag{3.3}
\]

modulo (3.1).  For each line component, these rows contain exactly eleven
distinct monomials.  Their coefficient matrix over \(L\) is \(11\) by
\(11\) and has exact rank \(11\).  In particular, its row space contains
each of the eight degree-one monomials in (3.2).  Hence all eight free
coordinates lie in the specialized reciprocal-Hessian ideal.  The linear
relations (3.1) then put all eighteen active coordinates

\[
b_0,\ldots,b_{14},u,v,w
\tag{3.4}
\]

in that ideal.

The denominator-cleared witness determinants also give a finite description
of the excluded lower strata.  Factorization over \(E\) has degree/exponent
profiles

\[
c^{10}\ell_1(m)^{17}g_3(m)
\tag{3.5}
\]

on the resultant line, and

\[
c^{10}q_2(m)^{17}g_7(m)
\tag{3.6}
\]

on the residual-polar line, where the subscripts are the irreducible degrees
over \(E\); in (3.5), \(\ell_1=m+40/11\).  The factor \(c=0\) is outside the
projective line chart.  Thus the remaining determinant boundary is a finite
linear-plus-cubic locus on the first component and a finite
quadratic-plus-septic locus on the second.

Only \(b_{15},b_{16},b_{17}\) survive.  They alter the bottom-right entry
of the boundary matrix but do not change its determinant, which remains
identically zero.

## 4. Result

> **Theorem `HC4NHM24` -- final six-point line-orbit generic gate.**
> Let \(\tau\) belong to either of the two six-point Fermat-symmetry orbits
> of non-generic polar fibers not covered by \(\tau^3=-1\).  At the generic
> point of each of the resultant and residual-polar line components, the
> complete normalized 81-equation reciprocal-Hessian packet has only the
> active-coordinate origin as support.  The surviving boundary matrix has
> determinant zero.  Consequently neither of the two remaining line-fiber
> normal forms supports a smooth irreducible quartic quotient generically.

Together with `HC4NHM23`, this removes every generic line-component type
from the fifteen degenerate polar fibers.  What remains inside those fibers
is lower-dimensional: linear-pivot and witness-determinant strata, including
the lower secondary strata already separated in the first orbit.  The
generic-conic parameterization denominators, hidden denominator divisors,
the complementary residual-line chart, and the other reciprocal boundary
types also remain separate.
