# Polar geometry of the smooth-quartic pivot

## Status

This note proves `HC4NHM20`.  It identifies the first visible denominator
\(\Delta\) of `HC4NHM16` invariantly as a first polar of a binary resultant
and gives its complete generic/degenerate fiber geometry.  The subsequent
theorem `HC4NHM22` closes the generic smooth polar-conic component, but the
lower denominator strata remain.  The later symmetry theorem `HC4NHM23`
reduces the fifteen non-generic two-line fibers to three automorphism types
and transports the `tau=-1` certificates across the first type.

Replay every identity with

```bash
.venv/bin/python scripts/verify_hc4_smooth_quartic_pivot_polar_geometry.py
```

## 1. The displayed pivot

Retain the squarefree-line normalization

\[
 F(s,t)=s^3+t^3,
 \qquad
 H(s,t)=p s^2-\frac q3st+r t^2.
 \tag{1.1}
\]

The Hesse cubic in the coefficients of \(H\) is exactly the binary
resultant

\[
 \mathscr R_F(H)
 =\operatorname{Res}(s^3+t^3,H)
 =p^3+\frac{q^3}{27}+r^3-pqr.
 \tag{1.2}
\]

Define the \(\tau\)-dependent quadratic

\[
 \begin{aligned}
 K_\tau(s,t)={}&(\tau^5+6\tau^2)s^2
 +3(\tau^4+\tau)st\\
 &+(6\tau^3+1)t^2.
 \end{aligned}
 \tag{1.3}
\]

Because the middle coefficient of \(H\) is \(-q/3\), variation by
\(K_\tau\) sends

\[
 (p,q,r)\longmapsto
 (p+\epsilon(\tau^5+6\tau^2),
 q-9\epsilon(\tau^4+\tau),
 r+\epsilon(6\tau^3+1)).
 \tag{1.4}
\]

Taking the first variation of (1.2) gives

\[
 \boxed{
 \Delta=\left.\frac{d}{d\epsilon}\right|_{\epsilon=0}
 \mathscr R_F(H+\epsilon K_\tau).
 }
 \tag{1.5}
\]

Expanding (1.5) is precisely

\[
\begin{aligned}
\Delta={}&(3p^2-qr)\tau^5+(9pr-q^2)\tau^4
 +(18r^2-6pq)\tau^3\\
&+(18p^2-6qr)\tau^2+(9pr-q^2)\tau+(3r^2-pq).
\end{aligned}
\tag{1.6}
\]

Thus the three repeated coefficient combinations in (1.6) are not
accidental: they are the scaled gradient of the Hesse resultant cubic.

## 2. Root-evaluation normal form

Over the algebraic closure, let \(\rho_1,\rho_2,\rho_3\) be the three roots
of \(F\), and put

\[
 h_i=H(\rho_i),\qquad k_i=K_\tau(\rho_i).
 \tag{2.1}
\]

Up to a nonzero scalar,

\[
 \mathscr R_F(H)=h_1h_2h_3
 \tag{2.2}
\]

and its first polar is

\[
 \boxed{
 \Delta=k_1h_2h_3+k_2h_1h_3+k_3h_1h_2.
 }
 \tag{2.3}
\]

For fixed \(\tau\), equation (2.3) is a conic in the coefficient plane
\(\mathbf P(S^2U)\).  The determinant of its Hessian matrix is

\[
 2k_1k_2k_3.
 \tag{2.4}
\]

Consequently the polar conic is smooth exactly when \(K_\tau\) is coprime
to \(F\).  If, say, \(k_1=0\), then

\[
 \Delta=h_1(k_2h_3+k_3h_2),
 \tag{2.5}
\]

so the conic splits into two lines.  The first is the resultant hyperplane
\(H(\rho_1)=0\); the second is its residual polar line.

## 3. Complete degeneration parameter

The degeneration resultant is

\[
\begin{aligned}
\operatorname{Res}(F,K_\tau)
={}&(\tau+1)(\tau^2-\tau+1)\\
&\cdot(\tau^4-4\tau^3+10\tau^2-4\tau+1)\\
&\cdot(\tau^8+4\tau^7+6\tau^6+32\tau^5+83\tau^4\\
&\hspace{3.2em}+32\tau^3+6\tau^2+4\tau+1).
\end{aligned}
\tag{3.1}
\]

It is squarefree of degree fifteen.  Hence:

1. away from the fifteen finite roots of (3.1), \(\Delta=0\) is a smooth
   polar conic in the \((p,q,r)\)-plane;
2. at each root, exactly one \(k_i\) vanishes and the fiber is the two-line
   arrangement (2.5);
3. no fiber has a double common root or a more degenerate polar conic.

The factors at \(\tau=-1\) seen in `HC4NHM17` are therefore one rational
fiber of a uniform polar degeneration, rather than an isolated coefficient
accident.

## 4. What contact geometry this does and does not encode

The cubic \(\mathscr R_F=0\) in the coefficient plane is the union of the
three hyperplanes of quadratics sharing one root with the fixed binary cubic
\(F\).  Equation (1.5) is its first polar at \(K_\tau\).  This is the exact
resultant/contact meaning of \(\Delta\).

It is **not** by itself the tangent, bitangent, flex, or hyperflex
discriminant of the prospective smooth quartic \(Q=\det(A)/(z\ell)\).  Those
contacts enter the reciprocal frontend through the separate condition

\[
 \gcd(F_x,F_y,G)^2\mid(Q|_L)(\ell|_L).
 \tag{4.1}
\]

Thus a direct classification of the components of \(\Delta\) by contact
types of \(Q\) would conflate two different geometries.  The correct finite
stratification of the visible pivot is instead:

\[
 \text{smooth polar conic}\quad\text{or}\quad
 \text{one resultant line plus one residual polar line}.
 \tag{4.2}
\]

## 5. Result and next calculation

> **Theorem `HC4NHM20` -- Polar-resultant pivot geometry.**  The visible
> pivot \(\Delta\) of the first basepoint-free smooth-quartic row is the
> first polar (1.5).  Its generic coefficient-space fiber is a smooth conic.
> Its complete degeneration locus is the squarefree degree-fifteen
> polynomial (3.1), and every degenerate fiber is the two-line arrangement
> (2.5).

The subsequent theorem `HC4NHM22` performs the first calculation proposed
here.  A rational parametrization from the universal resultant point
\([p:q:r]=[1:3:1]\), followed by an exact radical-membership calculation in
the 81-equation reciprocal-Hessian ideal, excludes the generic polar conic.
See
[`HC4_SMOOTH_QUARTIC_POLAR_CONIC_GENERIC_GATE.md`](HC4_SMOOTH_QUARTIC_POLAR_CONIC_GENERIC_GATE.md).
The next invariant calculation is therefore confined to its parametrization
and secondary-basis denominator strata and, after `HC4NHM23`, to the two
remaining six-point normal-form types of universal line components.  See
[`HC4_SMOOTH_QUARTIC_FERMAT_SYMMETRY_ORBITS.md`](HC4_SMOOTH_QUARTIC_FERMAT_SYMMETRY_ORBITS.md).
