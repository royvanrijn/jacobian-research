# Quadratic angular stability for nearly mirrored private unit profiles

Consider an ordinary two-monomial binomial with one private unit
row in each orientation and any shared rows. Discard equal rows and
let each remaining shared row with smaller exponent \(\ell\) and gap
\(d\) give its \(d\) Gauss slots \((\ell/d,h/d)\), \(1\le h\le d\).
Suppose there are equally many shared slots in the two orientations
and choose a pairing of opposite slots. In pair \(i\) write the two
slots as \((u_i+\delta_i,x_i+\xi_i)\) and
\((u_i-\delta_i,x_i-\xi_i)\). Thus \(u_i>0\) and \(0<x_i\le1\).
Define
\[
 \begin{gathered}
 \kappa=2\sum_i\frac{1}{u_i+1+x_i/2},\qquad
 q=(1+\kappa/4)^{-1},\quad b=1-q,\\
 g_i(s)=\delta_i+s\xi_i,\qquad
 a_i=\max\{|g_i(0)|,|g_i(1/2)|\},\qquad
 \eta=\max_i a_i/u_i,\\
 C=\sum_i\frac{a_i}{u_i^2},\qquad D=\sum_i\frac{a_i^2}{u_i^3},\\
 B=\frac14\sup_{0\le s\le1/2}\int_0^\infty
         y^2\left|\sum_i g_i(s)e^{-(u_i+s x_i)y}\right|\,dy .
 \end{gathered}
\]
A mirrored core can be included by choosing identical opposite
slots, which contribute to \(\kappa\) and have \(a_i=0\).
Assume \(\kappa>0\). Let \(J\ge4\) be an integer at least the
[effective absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

**Theorem.** The ordinary profile has disjoint consecutive root
arguments in every degree if there is a real \(R\) satisfying
\[
 \begin{gathered}
 R\ge\max\{2,3000b\sqrt J\},\qquad R\eta\le\frac12,\\
 \min\{4RC,\ RB+8R^2D\}\le\frac{b^5}{2^{16}}.                 \tag{1}
 \end{gathered}
\]
In particular SFC(2) holds for this class over every
characteristic-zero field. The angular conclusion holds under
independent positive rescalings in the two consecutive degrees.
For curvature in a fixed positive bounded interval, the first
test in (1) admits skew costs on the order of \(J^{-1/2}\).

The second test retains cancellation of signed slot deviations.
For example \(B=0\) if the sums of \(\delta_i\) and \(\xi_i\) both
vanish in each group with a common \((u_i,x_i)\). The elementary
bound \(B\le\frac12\sum_i a_i/u_i^3\) is available without
testing signs in the integral.

The joint matched tanh data determine both coordinates of each
slot. In terms of the lower and upper values \(T,T'\) at orders
\(n,N=n+1\), the coordinates are
\[
 \frac12\left(\frac n{T'}-\frac{n-1}T-1\right),\qquad
 \frac12\left(1+\frac{n^2-1}T-\frac{n^2}{T'}\right).
\]
Thus both tests in (1) depend on this joint signed data at a
single consecutive pair. They apply to pairings across different
integer row gaps as well as across rows of equal gap.
[Exact matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

## The mirrored reference slope

In degree \(m\), cancel the two private factorials against the
binomial coefficient and normalize the endpoint coefficients to
one by a positive scalar and a positive change of variable.
By the Gauss multiplication formula, an increasing slot \((u,x)\)
contributes the rising factors \((mu+x)_k\), and a decreasing
slot contributes \((mu+x)_{m-k}\), up to the factors removed by
normalization. Replace the two slots in each pair by their common
mean \((u_i,x_i)\) to obtain a mirrored reference \(P_m\) with
coefficients \(c_{m,k}\). Each reference pair represents opposite
real exponent rows of unit gap with common positive Gamma shape
\(x_i\). In degree \(m\ge2\) the continuous log coefficient has
second derivative at least \(\kappa/m\). Indeed each slot contributes
\(\psi'(mu_i+x_i+k)\), or its reflected expression. The trigamma
series and its integral comparison give \(\psi'(y)\ge1/y\) for \(y>0\).
The affine endpoint normalization leaves
second derivatives unchanged. The resulting quadratic chord
bound gives, with \(h_k=\min(k,m-k)\),
\[
 0<c_{m,k}\le \exp[-\kappa k(m-k)/(2m)]\le q^{h_k}.
                                                                  \tag{2}
\]
[Mirrored reference interlacing](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md)

There is a useful additional slope reserve when (2) attenuates
the endpoint coefficients. Use the mirrored resolvent notation
\(K_m=A_m-B_m e^{-i\theta}\) and
\(\psi_m(\theta)=m\theta/2-\arg K_m\). At a reference root,
\[
 |P_m'|=2|K_m|\psi_m',\qquad \psi_m'\ge m/4 .
\]
This includes the real negative root in odd degree by continuity.
[Mirrored resolvent slope bounds](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_UNIFORM_ANGULAR_RESERVE.md#the-resolvent-and-its-phase-slopes)

Write \(a=m/2\), let \(\rho\) be the folded midpoint law of the
absolute shared log ratio, and reweight it proportionally to
\(\cosh(as)\). Denote expectation in this reweighted probability
law by \(E\). Then
\[
 \begin{aligned}
 c_{m,1}&=E[\cosh((a-1)s)/\cosh(as)],\\
 A_m-B_m&=\frac12E\left[
  \frac{\cosh((a+1)s)/\cosh(as)-1}{\cosh s-\cos\theta}
                                    \right].
 \end{aligned}
\]
For \(x=\cosh s-1\) and \(y=\tanh(as)\sinh s\) one has
\(0\le y-x\le1\). The lower bound follows from
\(\cosh((a-1)s)\le\cosh(as)\) and the upper from
\(y\le\sinh s\). Also
\[
       \frac{x+y}{x+2}\ge \frac{y-x}{2},
       \qquad 2(x+y)-(x+2)(y-x)=x(4+x-y)\ge0.
\]
Consequently \(A_m-B_m\ge(1-c_{m,1})/4\). Since
\(|K_m|\ge A_m-B_m\), the derivative bound becomes
\[
                           |P_m'|\ge mb/8.                     \tag{3}
\]

Put \(\epsilon_m=b^2/(128m)\). In \(|z|\le1+2\epsilon_m\),
(2) gives \(|P_m''(z)|\le4m^2/b\), using
\(\sum_k q^{h_k}\le2/b\) and
\((1+2\epsilon_m)^m<2\).
On a circle of radius \(\epsilon_m\) centered at a reference zero
\(\zeta\), (3) and Taylor's formula give
\[
                         |P_m(z)|\ge mb\epsilon_m/16.           \tag{4}
\]
The linear Taylor term also dominates the remainder on the
circle of radius \(2\epsilon_m\). Hence each larger disk contains
exactly one reference zero, and the smaller disks are disjoint.

## An analytic deformation respecting reversal

In each pair replace the common mean slot by
\((u_i+t\delta_i,x_i+t\xi_i)\) in the increasing orientation and by
\((u_i-t\delta_i,x_i-t\xi_i)\) in the decreasing orientation.
Apply the same chord normalization to the complex family,
using analytic logs of the individual rising factors with
positive real part. Let \(F_m(t,z)\) be its polynomial.
It agrees with \(P_m\) at \(t=0\). At \(t=1\) the normalized
ordinary polynomial is recovered.

For a Gauss slot of base \(A\), the normalized increasing log
coefficient is
\[
 f_{m,k}(A)=\log(A)_k-(k/m)\log(A)_m .
\]
The decreasing slot uses \(f_{m,m-k}\). If \(\Re A\ge a_*>0\),
\[
 |f_{m,k}'(A)|\le \frac{k(m-k)}{2a_*^2},\qquad
 |f_{m,k}''(A)|\le \frac{k(m-k)}{a_*^3}.
\]
For \(0<k<m\), the first derivative is \(k(m-k)/m\) times the
difference of the means of \((A+j)^{-1}\) over \(0\le j<k\)
and \(k\le j<m\). Each paired difference is at most the index
distance divided by \(a_*^2\), and the mean distance is \(m/2\).
The second derivative uses the negative squared reciprocals.
Their paired differences are at most twice the index
distance divided by \(a_*^3\), by integration between the
indices. Both endpoint derivatives vanish.

Write \(A_i=mu_i+x_i\). The two deformed bases are
\(A_i\pm tm g_i(1/m)\), so their real parts are at least
\(m u_i/2\) for \(|t|\le R\) by (1).
Integrating the first derivative along the segment from
\(0\) to \(t\) bounds the total log coefficient change by
\(4RC\,k(m-k)/m\le 4RC h_k\). Taylor's formula with the
second derivative bounds the error after the linear term
by \(8R^2D\, k(m-k)/m\le 8R^2D h_k\).

To bound the signed linear term, put \(s=1/m\) and use
\((A+j)^{-1}=\int_0^\infty e^{-(A+j)v}\,dv\) at \(A=A_i\).
The first variation in \(t\) is
\[
 \int_0^\infty G_{m,k}(y/m)
               \sum_i g_i(s)e^{-(u_i+s x_i)y}\,dy ,
 \quad
 G_{m,k}(v)=\frac{m-2k}{m}\sum_{j=0}^{m-1}e^{-jv}
                              -\sum_{j=k}^{m-k-1}e^{-jv}
\]
for \(k\le m/2\). For \(k>m/2\) use
\(G_{m,k}=-G_{m,m-k}\).
For \(0<k<m/2\), \(G_{m,k}\) equals \(2k(m-2k)/m\)
times the difference of the mean of \(e^{-jv}\) on the end
blocks and its mean on the middle block. Both sets of
indices have the same mean \((m-1)/2\). The Taylor remainder
of the convex function \(e^{-jv}\) about this center is
between zero and \(m^2v^2/8\) on the entire index interval.
It follows that \(|G_{m,k}(y/m)|\le h_k y^2/4\),
including the endpoint and middle indices where \(G=0\).
The linear variation has modulus at most \(Bh_k\).
This proves both log coefficient bounds in (1).

Put \(M=\min\{4RC,RB+8R^2D\}\).
From (2) and \(|e^w-1|\le |w|e^{|w|}\), on
\(|z|\le1+\epsilon_m\) the evaluation difference is at most
\(16M/b^2\). In detail \(M\le b/2\),
\(q e^{b/2}\le1-b/2\), and
\(\sum_k h_k(1-b/2)^{h_k}\le8/b^2\).
The remaining factor \((1+\epsilon_m)^m\) is less than two.
The upper bound is at most \(b^3/4096\) by (1),
whereas the right side of (4) equals \(b^3/2048\).
Rouché therefore gives a unique simple root \(z_\zeta(t)\)
in each disk when \(|t|<R\). The roots depend holomorphically
on \(t\) by the implicit function theorem and uniqueness.

The exact reciprocity identity is
\(F_m(-t,z)=z^mF_m(t,1/z)\). The coefficients also obey
complex conjugation in \(t\). For sufficiently small real
\(t\), uniqueness in the root disk gives
\(z_\zeta(-t)=1/\overline{z_\zeta(t)}\).
Define \(L_\zeta(t)=\log(z_\zeta(t)/\zeta)\) using its analytic
branch with \(L_\zeta(0)=0\). The reflection identity extends
to the parameter disk as
\(L_\zeta(-\bar t)=-\overline{L_\zeta(t)}\).
Thus \(\Im L_\zeta(t)\) for real \(t\) is the imaginary part of
the even function \((L_\zeta(t)+L_\zeta(-t))/2\). That function
vanishes to order at least two and has modulus at most
\(2\epsilon_m\) throughout the parameter disk. The Schwarz
lemma applied twice yields the angular displacement bound
\[
                    |\Im L_\zeta(1)|\le 2\epsilon_m/R^2 .
                                                                  \tag{5}
\]

## Consecutive degrees

The mirrored reference has angular separation greater than
\(1/(150000000N^2)\) at terminal \(N=n+1\).
For \(2\le n<N<J\), the sum of both bounds in (5) is
below \(b^2/(32nR^2)\). This is smaller than the reference gap:
use \(N^2/n\le 3J/2\), \(R^2\ge9000000b^2 J\), and
\(225000000<288000000\).
Terminal two follows from the real linear root and strict positivity
of the second moment on the real axis. Terminals \(N\ge J\) use
the absolute cutoff. Since the ordinary polynomials have integer
coefficients, complex coprimality gives a rational Bezout identity
and hence coprimality over every characteristic-zero field.
[Quantitative reference gap](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_UNIFORM_ANGULAR_RESERVE.md#a-quadratic-lower-gap-in-the-order)

The criterion includes ordinary mixed profiles at intermediate
curvature. For even \(L\ge2\), take \(L/2\) copies of
\((L,L+1),(L+2,L+1)\) in addition to the private unit pair.
Their shared absolute curvature
\((L/2)[1/(L+1)+1/(L+2)]\) tends to one, and their reference
\(\kappa=L/(L+2)\) is at least \(1/2\).
Here \(C=L/(2L+1)^2\) and \(\eta=1/(2L+1)\).
For each fixed admissible \(J\), (1) holds for every
sufficiently large even \(L\).

Pairing can also unfold a shared row into more than one row.
For \(4\mid L\), take \(L/4\) blocks, each having the
increasing row \((2L,2L+2)\) and two decreasing rows
\((L+1,L)\). In each block, one slot \((L,1)\) pairs
with an identical opposite slot and the other pairs
\((L,1/2)\) with \((L,1)\). Thus \(\eta=1/(8L)\),
\(C=1/(32L)\), and both the reference and ordinary shared
curvatures tend to one. These unbalanced row counts
also satisfy (1) for sufficiently large \(L\).
