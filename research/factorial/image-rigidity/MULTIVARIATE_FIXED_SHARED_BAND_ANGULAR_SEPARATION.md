# Eventual angular separation on fixed shared curvature bands

Remove equal rows from two distinct ordinary exponent columns, and
assume all remaining rows are shared. Use \(W,\rho,w\) for the shared
curvature, real stationary curvature, and chord-centered field in the
[effective fixed-band proof](MULTIVARIATE_EFFECTIVE_SHARED_REAL_BOUNDARIES.md).

**Theorem.** Given positive rational \(W_-\le W_+\) and rational
\(0<\kappa<1\), one can compute a terminal threshold \(J\) with the
following property. For \(N\ge J\) and \(W_-\le W\le W_+\), if
\(\max\rho\le1-\kappa\) or \(\max\rho\ge1+\kappa\), the orders
\(N-1,N\) have disjoint root arguments under independent positive
rescalings. The threshold is uniform in the ordinary profiles and
is unevaluated.

The field rectangle, limiting saddle geometry, and actual contour
witnesses are those of the cited proof. The changes concern the
radial comparison after the two fields are put at the same height.
In the subcritical fold chart we use a complex interpolation of the
local Airy integral on a fixed relative order disk.

## The paired field geometry

Write \(n=N-1\). Each field is the negative logarithm of the
corresponding root, up to a real profile and order shift. Equality of
arguments therefore makes \(w_N-w_n\) real. The independent
root-annulus bound in the fixed-band proof puts both fields in
\[
             |\Re w|\le M,\qquad -\pi\le\Im w\le0
\]
after conjugation. Use the same equal-height cells in the compact
profile family as in its [regular contour cover](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md).

For each limiting pair in the supercritical band, if both fields are
on the active equal-modulus graph then they have the same field by
uniqueness at a given height. The corresponding statement in
the subcritical band follows from strict height monotonicity on the
closed active arc, including its simple fold. An equal-height
point off the active graph or arc has a dominant-term witness.
On the lower boundary the real-gap pair has one balance scale.
At the upper boundary the sole two-term chart is the supercritical
coexistence scale, with one-term exclusion at every other real
field in these bands. These statements and the contour witnesses
are the finite witness alternatives used by the fixed-band proof.
A common chart is therefore available at every limiting pair
where both fields can be roots.

## Regular interior and lower boundary charts

A common convex regular chart supplies a phase difference \(D\)
with \(|\Re D'|\) bounded below and \(|e^D-1|\) bounded below.
The analytic nonvanishing Gaussian amplitude ratio \(L\) is the
same at both orders, by the uniform expansion in the regular
cover. At the exact roots its two-term equation gives
\[
                  e^{jD(w_j)}=L(w_j)(1+O(1/j)).
\]
The relative error is uniform with computable constants. This
equation also applies at a common two-gap chart on the lower
boundary. The exact reflection factor there belongs to the
exponential phase and gives the multiplier \(-1\) at the balance.

With \(\delta=w_N-w_n\in\mathbb R\), the modulus equations give
\(|\Re D(w_j)|=O(1/n)\) and then \(|\delta|=O(1/n)\) by the real
slope. Subtract the two real equations. The difference of log
amplitudes is \(O(|\delta|)\), the extra real phase at order \(n\)
is \(O(1/n)\), and the residual is \(O(1/n)\).
Consequently
\[
                        |\delta|=O(n^{-2}).
\]
This is the real comparison of the
[weak-curvature angular proof](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md#regular-modulus-comparison)
with its \(O(1/j)\) version of the expansion.
Dividing the complex zero equations gives
\(e^{D(w_n)}=1+O(1/n)\), contradicting the phase margin.

## Supercritical coexistence at the upper boundary

The [upper boundary expansion](MULTIVARIATE_EFFECTIVE_SHARED_REAL_BOUNDARIES.md#upper-boundary)
on a common convex coexistence chart is
\[
           e^{jD(w_j)}=L(w_j)(1+B(w_j)/j+O(j^{-2})).
 \tag{1}
\]
Here \(D,L,B\) are analytic functions for the fixed profile at
both orders. The stationary difference \(D\) is real on the real
axis, its real derivative is bounded away from zero, and the
negative amplitude ratio \(L\) at the real center gives
\(\Re L\le-\eta<0\) on the chart. Choose the chart small enough
that \(|D|\le1/4\). The chart can be chosen symmetric about the real
axis for these analytic functions, using only the lower half for
the two order fields.

As in the regular case, the modulus equations give
\(\delta=w_N-w_n=O(n^{-2})\).
Now the quotient of the amplitude ratios is \(1+O(n^{-2})\).
Also \(B(w_N)/N-B(w_n)/n=O(n^{-2})\) by analyticity and
bounded derivative. Dividing (1) gives
\[
           \exp(D(w_n)+N(D(w_N)-D(w_n)))=1+O(n^{-2}).
\]
The exponent has modulus less than \(1/2\) for sufficiently
large \(n\), so it is itself \(O(n^{-2})\).
Write \(y=\Im w_n=\Im w_N\). Since \(D'\) is real on the real
axis, \(|\Im D'(x+iy)|\le C|y|\) on the joining segment.
Thus \(|\Im D(w_n)|\le C|y|/n+O(n^{-2})\).
On the vertical segment from \(\Re w_n\) to \(w_n\), the
real derivative of \(D\) has a fixed sign and positive lower
bound. It follows that \(|\Im D(w_n)|\ge c|y|\).
Together these bounds give \(y=O(n^{-2})\).
The real part of \(e^{nD(w_n)}\) is now positive, contrary to
(1) and the negative real part of \(L\).
All constants follow from the common chart margins and analytic
bounds in the upper boundary witness.

## The subcritical fold and complex order

Use a certified elementary Airy chart from the
[effective simple-fold proof](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_FOLD.md).
It has the local normal form
\[
           \phi(q;w)=C(w)+u^3/3-\lambda(w)u
\]
and \(|\Im\lambda'|\) bounded away from zero on a common convex
disk. This derivative bound uses the rotation chosen by the
actual oriented Airy contour. In the exceptional one-row case
the horizontal equal-modulus ray is inactive, as in the cited proof.

The analytic finite-order interpolation needed here can be
constructed on the two fixed truncated Airy rays. For complex
order \(\nu\) with \(|\nu-n|\le\eta n\), use the same local
Gamma amplitude \(A_\nu(q(u,w))q_u(u,w)\), after the phase and
nonzero scalar have been factored out. The chart is compact in
the upper saddle domain. At its saddles and the finite rays
the Gamma arguments divided by the order stay away from zero
and have a strict sector margin. The same holds upon small
relative complex changes of the order, by decreasing \(\eta\).
The sectorial Stirling remainder in the regular cover thus gives
uniform analytic amplitude bounds with all fixed derivatives on
smaller disks. This is the local complex-order analytic amplitude
used in the
[real-order fold interpolation](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_FOLD.md#bounds-on-each-certified-chart).
Only the finite local rays need to be interpolated. The actual
integer-order omitted contour, including the joins to these rays,
has strict exterior loss.

Choose the ray radius inside the analytic amplitude domain.
Absorb a nonzero scalar and \(e^{\nu C(w)}\nu^{-1/3}\) in the
interpolation normalization. With \(z=\nu^{2/3}\lambda(w)\)
the normalized finite integral has the form
\[
 I(\nu,w)=P(\nu,w)\operatorname{Ai}(z)
     +\nu^{-1/3}Q(\nu,w)\operatorname{Ai}'(z)+E(\nu,w),\qquad
 |E|\le \frac{C_0}{n}
       \left(|\operatorname{Ai}(z)|+
                \frac{|\operatorname{Ai}'(z)|}{\sqrt{1+|z|}}\right).
 \tag{2}
\]
Here \(P\) stays away from zero, and \(P,Q\) are analytic and
bounded on the chart, with bounded field derivatives.
To obtain (2), divide the amplitude by \(u^2-\lambda\) with
the constant and linear remainder given by Cauchy integrals.
One integration by parts of the divided term gains \(1/\nu\).
After \(v=\nu^{1/3}u\), the rays have a small common rotation.
The
[effective Airy joint norm](AIRY_EFFECTIVE_ELEMENTARY_JOINT_NORM.md)
applies to these rotated finite rays when the physical bound
on \(|\lambda|\) is small enough. Its lower bound absorbs the
fixed-radius boundary terms and the strict-loss exterior chains
after shrinking \(|\lambda|\) further. At the integer orders the
difference from the original normalized polynomial is
\(O(e^{-c_0 n})\) times the joint Airy norm. This uses only the
integer contour identity and its exterior loss. The choices
of the complex-order disk, the amplitude bounds, and these
loss bounds are computable local tests.

At a true integer-order zero, (2) and the Airy integral
equation give \(z_j=a_{k_j}+O(n^{-1/3})\), with an Airy zero
\(a_{k_j}<0\) and \(|\Im\lambda(w_j)|=O(1/n)\).
The horizontal derivative lower bound yields \(|\delta|=O(1/n)\).
Therefore \(|z_N-z_n|=O(n^{-1/3})\), including the change in the
factor \(j^{2/3}\). As in the
[weak-curvature label comparison](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md#airy-labels-at-the-same-imaginary-field),
the Airy ODE spacing lower bound \(c_1(1+|a_k|)^{-1/2}\)
forces the same label \(k\) after shrinking the physical chart.
In detail, the constant in the \(O(n^{-1/3})\) difference is
bounded on the original chart. Then
\((1+|a_k|)n^{-2/3}\) can be made sufficiently small by decreasing
the admitted \(|\lambda|\) and increasing \(n\).

For that label substitute
\(w=\lambda^{-1}(\nu^{-2/3}z)\) in (2) on a disk of radius
\(C_2 n^{-1/3}\) about \(a_k\), using a fixed smaller relative
complex order disk. The substitution stays in the complex field
chart by another fixed decrease of the admitted \(|\lambda|\).
Divide by \(P\operatorname{Ai}'(a_k)\).
The local Airy integral equation bounds this expression by
the affine model \(z-a_k\) with the controlled
\(\nu^{-1/3}(Q/P)\) correction and \(O(1/n)\) remainder.
Choose \(C_2\) before shrinking \(|\lambda|\) so the affine
term dominates on the circle. The same Rouché and derivative
bounds in the weak-curvature Airy proof give one analytic zero
\(z=a_k+b_k(\nu)\), with \(b_k=O(n^{-1/3})\) and
\(b_k'=O(n^{-4/3})\) by Cauchy in order.
They also compare the true integer roots to this branch with
error \(O(e^{-c n})\) in the physical field.
Its derivative there is
\[
  w'(\nu)=(\lambda^{-1})'
       \bigl(-\tfrac23 a_k \nu^{-5/3}+O(n^{-2})\bigr).
\]
The inverse has imaginary derivative with fixed sign and
magnitude bounded below. As \(|a_k|\) stays bounded away from
zero, the imaginary physical displacement over one order has
magnitude at least a positive multiple of \(|a_k|N^{-5/3}\)
for sufficiently large \(n\). The exponential errors cannot
remove it. This contradicts equality of physical arguments.

## Computable compact cover

The regular interior, lower boundary, and upper boundary
comparisons use the same finite witnesses as the effective
fixed-band proof. Supply at each subcritical fold a local
elementary Airy witness as above, choosing its local order disk
and analytic losses before reducing the admitted field radius
for the Airy labels and branch comparison. This gives a
computable local threshold also at every fold pair.
By the active-arc and equal-height geometry, all remaining
pair cells have a one-term exclusion at one order. Thus every
admitted limiting pair has a strict open certificate. The same
dovetailing search over profile and field nets and finite
contour witnesses terminates by compactness. Its largest local
threshold proves the stated uniform angular comparison.
