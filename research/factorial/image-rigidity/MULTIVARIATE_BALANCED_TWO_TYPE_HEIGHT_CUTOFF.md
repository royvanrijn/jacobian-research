# Effective all-order cutoff for balanced opposed factorial types

Fix any ordinary mirrored core, possibly empty. A mirrored core
is a multiset of pairs \((a,b),(b,a)\), with equal rows immaterial.
Append \(r_+\) copies of \((L,L+d)\) and \(r_-\) copies of
\((Y+e,Y)\), with fixed positive integer gaps \(d,e\).
The ordinary smaller exponents \(L,Y\) are nonnegative integers
and the counts are positive integers with \(dr_+=er_-\).

**Theorem.** A computable cutoff depending on the fixed core and
the two gaps excludes exact common roots in every order once both
heights exceed it, uniformly over all balanced counts. Moreover the
finite set of possible collision profiles in this family is
computable. The core is allowed to have a consecutive collision
on its own. Both assertions concern the profile with the nonempty
tail appended.

A physical collision forces a relative-center error proportional
to the product of tail curvature and relative gap. For two
balanced types this implies a bound on the integer difference
\(Y-L\). Exact affine-height separation then treats the finite
set of offsets with no bound on the tail counts.

## The radial estimate forces bounded height difference

Use a computable absolute terminal cutoff \(J\) and a rational
diffuse tolerance \(\eta_0\le1/2\) as in the
[mirrored signed-tail theorem](MULTIVARIATE_MIRRORED_CORE_SIGNED_DIFFUSE_RADIAL_SEPARATION.md).
In its linear radial estimate let
\[
              g=\min\!\left(1,\frac{4}{J(2+1/\eta_0)}\right),
       \qquad B=J3^J(2/g)^{J-1},\qquad K=32BJ^2 .
\]
For the tail write
\[
 \tau=r_+\frac{d^2}{L+d}+r_-\frac{e^2}{Y+e},\qquad
 \eta=\max\!\left(\frac d{L+d},\frac e{Y+e}\right),\qquad
 \mathcal L=r_+\log(1+d/L)-r_-\log(1+e/Y).
\]
Here \(L,Y>0\), which is automatic above the proposed height
cutoff. The total signed degree difference is zero, by the
mirrored core and the exact tail balance. For
\(\eta\le\min(\eta_0,g/(54B))\), apply the
[linear radial criterion](MULTIVARIATE_MIRRORED_CORE_SIGNED_DIFFUSE_RADIAL_SEPARATION.md#linear-radial-stability-and-the-all-degree-criterion)
with the positive rational parameter \(\kappa=K\eta\).
Its required inequality \(\eta\le\kappa/(32BJ^2)\) is an equality.
Hence a collision at any order requires
\[
                              |\mathcal L|<K\eta\tau .           \tag{1}
\]

Write \(G=dr_+=er_->0\), \(m=\min(L,Y)\), and \(D=\max(d,e)\).
With
\[
                f_d(L)=\frac{\log(1+d/L)}d
                           =\int_0^1\frac{dt}{L+dt},
\]
one has \(\mathcal L/G=f_d(L)-f_e(Y)\),
\(\tau/G=d/(L+d)+e/(Y+e)\), and \(\eta\le D/m\).
Also
\[
              \left|f_d(L)-\frac1L\right|
                  =\int_0^1\frac{d\,t}{L(L+d t)}\,dt
                  \le \frac d{2L^2}.
\]
Thus (1) gives
\[
          \left|\frac1L-\frac1Y\right|
                   \le \frac{C}{m^2},\qquad
          C=(d+e)(KD+1/2).
\]
If \(m\ge2C\) and \(M=\max(L,Y)\), then
\(1-m/M\le C/m\le1/2\), so \(M\le2m\).
Multiplying the reciprocal bound by \(mM\) yields
\[
                                  |Y-L|\le2C .                 \tag{2}
\]
All these constants are effective and independent of the
balanced count \(G\). The choice of the initial lower height
to make the relative-gap condition hold is effective as well.

## Finite affine tests

At a possible collision with both heights large, (2) makes
\(Y=L+k\) for one of the integers \(|k|\le\lceil2C\rceil\).
Apply the
[effective opposed affine-height consequence](MULTIVARIATE_OPPOSED_INDEPENDENT_HEIGHT_RATIONAL_BIAS.md#first-marks-exclude-the-forced-affine-family)
to each pair \((L,L+k)\); it covers arbitrary independent
counts above a computable height, allowing unequal gaps and
every fixed ordinary core. The finite maximum of these
thresholds and the radial threshold is the asserted cutoff.
In particular it does not require regularity of the core.

Any remaining collision has a bounded value for one height.
That value and its fixed positive gap bound the corresponding
count by the [absolute curvature cutoff](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md).
Exact balance bounds the other count. Absorb the bounded-height
row into the core and use the
[all-directions affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md)
for the remaining height with both counts fixed. Its effective
finite list can be checked by exact resultants through the
absolute terminal cutoff, which proves the second assertion.

For the empty core an independent, broader conclusion is the
[finite reduction for two distinct unequal row types](MULTIVARIATE_SOFT_TYPE_FINITE_REDUCTION.md).
It allows either orientation, arbitrary signed difference,
and unbounded gaps without strengthening the height
hypotheses. The mirrored-core assertion above goes beyond that
two-type case by retaining all the fixed mirrored rows.

## A finite circle-root test on the core

The same conclusion holds for a fixed ordinary core that
passes the following finite test, regardless of whether its
rows are literally mirrored. At each adjacent pair
\(2\le n<N=n+1<J\), require the core endpoint scales
\(G_j^C=([z^j]Q_j^C/Q_j^C(0))^{1/j}\) to satisfy
\(G_n^C=G_N^C\), and require both endpoint-normalized core
polynomials to have simple unit-circle roots. These are
algebraic tests on the exact core coefficients.

For a core passing the test, the Gaussian symbol applied
to each core parent preserves the circle location and a
positive circular mesh. This is the degreewise [mesh
preservation under Schur–Szegő composition](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#actual-mirrored-reference-profiles),
using the core factor's mesh and a circle-rooted Gaussian
factor. The mesh is uniformly bounded below at the
finitely many tested orders by the corresponding core
meshes, which are computable from algebraic roots.
The mirrored radial proof of (1) uses only this positive
mesh, the log-coefficient perturbation of the tail, and the
equality of the core scales. Its coefficient estimate
uses factorial log convexity and applies to this same
core since \(P_j^C\) is endpoint-normalized with positive
coefficients. In the relative center comparison the core
scales again cancel, just as for a mirrored core. Thus the
Rouché estimate gives (1) with a computable constant
obtained from the core mesh. The offset and finite decision
arguments then apply unchanged.
