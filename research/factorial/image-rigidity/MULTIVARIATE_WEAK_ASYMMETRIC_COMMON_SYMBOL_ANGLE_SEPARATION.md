# Weak-curvature angular separation with independent factorial-symbol blends

For distinct fully shared ordinary exponent columns let
\(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\), and
\(W=\sum_i d_i^2/m_i>0\), omitting equal rows. Let \(L_j\)
be the endpoint-normalized log-factorial symbol per degree.
At consecutive orders \(n,N=n+1\), choose arbitrary
\(\lambda_n,\lambda_N\in[0,1]\). Use the polynomials
\[
 P_m(z)=\sum_{k=0}^m\binom mk
 \exp\!\left(m[(1-\lambda_m)L_n(k/m)+\lambda_m L_N(k/m)]\right)z^k,
                       \qquad m=n,N .
\]

**Theorem.** There are absolute computable constants
\(\eta>0,c>0\) such that if \(W\le\eta\), then for every
\(n\ge1\) the two sets of root arguments have circular
distance at least \(c\sqrt{W/2}\,n^{-5/3}\).
The cap and gap are uniform in the number and signs of
the rows and in both blend parameters. Either symbol can
therefore be frozen on both grids. The conclusion also
holds for real smaller exponents at least one with
integer gaps. The constants are unevaluated.

The [matched midpoint hierarchy](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md#proof-of-the-curvature-hierarchy)
compares \(C=L_n''(1/2)\) and \(W\) by
\(W/2\le C\le W\). Thus \(C\le\eta/2\) suffices. The signed hierarchy controls the bounded
degrees. For growing degrees the exact midpoint law
for the blends produces one stationary phase and a
common Airy branch without a circular-root assumption.

## The exact positive midpoint laws

Put \(\varepsilon=\sqrt{W/2}\), \(e_i=d_i/\varepsilon\),
and \(u_i=e_i/m_i\). The
[ordinary midpoint Gamma moment function](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md#exact-midpoint-normalization-and-effective-density-bounds)
is
\[
 T_j(h)=\prod_i
  \frac{\Gamma(jm_i+1+e_i h)}{\Gamma(jm_i+1)}
                         (jm_i)^{-e_i h}.
\]
For both integer orders define
\[
 D_m(h)=m\!\left[
   (1-\lambda_m)\frac{\log T_n(nh/m)}{n}
            +\lambda_m\frac{\log T_N(Nh/m)}{N}\right].
                                                               \tag{1}
\]
This is a log moment generating function of a real
probability law. In fact for \(p>0\),
\[
 \log\frac{\Gamma(p+z)}{\Gamma(p)}
   =z\psi(p)+\int_0^\infty
       (e^{-zu}-1+zu)\,\frac{e^{-pu}}{u(1-e^{-u})}\,du
                    \qquad \Re z>-p .
\]
Integrating the trigamma series twice verifies the
identity, including the constant and derivative at zero.
Its positive compensated measure defines a law by
truncation to compound Poisson sums and convergence
of the second moment near zero. Real dilation and
nonnegative linear combinations of the log moment
functions preserve this construction. This is the
classical infinite divisibility of the log-Gamma
law. [Jurek, *A note on gamma random variables and Dirichlet series*](https://doi.org/10.1016/S0167-7152(00)00071-7)

If \(X_m\) has the law in (1) and
\(g_\varepsilon(v)=2\sinh(\varepsilon v/2)/\varepsilon\),
then \(\mathbb E[g_\varepsilon(X_m+\xi)^m]\) is the
polynomial value at \(z=-e^{\varepsilon\xi}\), up to a
nonzero factor and a positive radial scaling.
To see this, for \(q=k/m\) substitute
\(h=\varepsilon(k-m/2)\) in (1). Each
\(\log T_j(j\varepsilon(q-1/2))/j\) has the same
endpoint-normalized chord as \(L_j(q)\), since its
subtracted terms are affine in \(q\).
Expanding the power of the hyperbolic sine uses
exactly these \(m+1\) exponential moments.
The affine discrepancy is real and gives only a
scalar and a positive change of \(z\).
Thus one may choose \(|\Im\xi|\le\pi/\varepsilon\),
and imaginary fields compare physical angles
independently of those scalings.

## A shared analytic interpolation in order

For complex \(m\) near \(n\), use
\(\lambda(m)=\lambda_n+(m-n)(\lambda_N-\lambda_n)\)
in (1), with the analytic Gamma logarithms.
The period \(n\) or \(N\) inside each \(T_j\) is
fixed. The weak-density functions are
\[
 \begin{split}
 K(t)&=\sum_i m_i[(1+u_i t)\log(1+u_i t)-u_i t],\\
 J_1(t)&=\tfrac12\sum_i\log(1+u_i t),\qquad
 b_m=m/n-m\lambda(m)/(nN).
 \end{split}
\]
For \(h=mt\), Stirling gives
\[
                     D_m(h)=m K(t)+b_mJ_1(t)+R_m(t).
                                                               \tag{2}
\]
For a fixed sufficiently small relative disk
\(|m-n|\le\eta_1 n\), \(b_m\) is uniformly bounded.
The remainder is \(O(W)\) on that disk and
\(O(W/n)\) at the two integer endpoints. At the
endpoints \(b_m=1+O(1/n)\). The distinction comes
from the at most \(O(n)\) size of \(\lambda(m)\)
on the complex-order disk, while the \(jK(t)\)
terms in the two Gamma logarithms combine exactly.

Here are the global estimates needed to apply the
existing density computation to (2). On
\(h=n(\alpha+iv)\) for \(\alpha\) in a fixed bounded
interval and \(\varepsilon\) small,
\[
      \Re(m+u_i h)\ge c n .
\]
Hence \(1+u_i h/m\) has modulus bounded below and
argument in a sector of width strictly less than
\(2\pi\) about the positive ray, uniformly bounded
away from the negative real ray. The logarithms
for the two \(T_j\)'s can both use this sector.
The row Gamma arguments are of order at least
\(n m_i\) in modulus, with the same sector margins,
so each \(T_j(jt)\) has a uniform \(O(W/n)\)
absolute logarithmic Stirling remainder on the
whole vertical line, by \(\sum_i1/m_i\le W\).
This proves the stated bounds for \(R_m\).
Integrating the \(u_i\)-log derivative on the
segment from \(1\) to \(1+u_i h/m\) gives
\(|J_1(t)|\le C\varepsilon(1+|v|)\), using
\(\sum_i |u_i|\le2\varepsilon\). These segments
stay in the rotated right half-plane
\(\Re(m z)\ge c n\).

The real second \(v\)-derivative of \(m K(h/m)\)
is
\[
   -n^2\sum_i e_i^2\,
          \Re\frac1{m m_i+e_i h}
              \le-\frac{c n}{1+\varepsilon^2 v^2}.
\]
This is the same estimate as in the
[effective density calculation](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md#exact-midpoint-normalization-and-effective-density-bounds).
Indeed the numerator in each reciprocal is
at least a constant times \(n m_i\) and its
denominator is bounded by
\(C n^2 m_i^2(1+\varepsilon^2 v^2)\).
Use \(\sum_i e_i^2/m_i=2\).
Along the vertical line the \(J_1,R_m\) bounds
therefore give the same integrated growth bound
\(\exp(C n(1+|v|)-c' n v^2/(1+\varepsilon |v|))\)
for the inverse Fourier integrand on bounded
complex \(s\)-disks. This proves holomorphy and
the validity of the same vertical shifts.

In particular, with \(t(s)\) solving \(K'(t)=s\)
near \(s/2\), let \(A(s)=K(t(s))-s t(s)\).
The Gaussian density calculation with the exact
analytic amplitude \(e^{b_mJ_1+R_m}\) now gives
\[
 p_m(s)=\sqrt m\,e^{mA(s)}a_m(s),\qquad
 A(s)=-s^2/4+O(\varepsilon).
                                                               \tag{3}
\]
On every fixed bounded \(s\)-disk all the fixed
derivative bounds in that calculation are uniform.
On the complex-order disk, \(a_m\) is analytic,
bounded, and close to \(1/\sqrt{4\pi}\) by
\(O(\varepsilon+1/n)\). At the integer endpoints
it satisfies
\[
       a_m(s)=\frac{e^{J_1(t(s))}}
                      {\sqrt{2\pi K''(t(s))}}+O(1/n).
                                                               \tag{4}
\]
In detail the fixed-disk amplitude of the shifted
vertical saddle is
\(e^{b_mJ_1(t(s))+R_m(t(s))}/\sqrt{2\pi K''(t(s))}\);
bounded derivatives of this amplitude give the
same relative Gaussian \(O(1/n)\) error after the
odd terms are integrated out. Complex \(m\) stays
in a small relative disk, so the strict Gaussian
quadratic loss and the principal square root remain
valid. The real tails are likewise uniform:
at each integer order, (2) at \(t=\pm c_1\)
gives \(\mathbb E[e^{m c_1|X_m|}]\le2e^{C_{c_1}m}\)
for every fixed \(c_1>0\) after shrinking the cap.

## Localization of both off-circle fields

The Gaussian large-field exclusion can be run
without an annulus theorem. For
\(|\Im\xi|\le\pi/\varepsilon\), \(\xi\ne0\), the
identity for the modulus of \(\sinh\) gives
\[
            \varepsilon|\coth(\varepsilon\xi/2)|
                        \le\varepsilon+\pi/|\xi|.
\]
Indeed \(|\sinh(\varepsilon\xi/2)|^2\ge
\varepsilon^2|\xi|^2/\pi^2\) on this strip, and
\(|\cosh v|^2\le1+|\sinh v|^2\).
For \(|\xi|\ge Y\), the ratio
\[
 \frac{g_\varepsilon(s+\xi)}{g_\varepsilon(\xi)}
  =\cosh(\varepsilon s/2)+
          \coth(\varepsilon\xi/2)\sinh(\varepsilon s/2)
\]
therefore differs from one by \(O_R(\varepsilon+1/Y)\)
on each bounded complex \(s\)-disk of radius \(R\).
On the real line its modulus is bounded by
\(\exp((\varepsilon/2+C/Y)|s|)\), as follows by
bounding \(\cosh v+b\sinh |v|\) by
\(\exp(\max(1,b)|v|)\).
The proof of the
[large-field one-saddle exclusion](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md#uniform-regular-comparison-off-a-fixed-fold-neighborhood)
now applies without restricting \(\Re\xi\).
Choose the real tail cutoff using the uniform
exponential moments, then \(Y\) large and the
curvature cap small. On the compact retained
disk the phase of the ratio plus \(A(s)\) is
close with derivatives to \(-s^2/4\).
Its single Gaussian contribution dominates the
other pieces and is nonzero.

To localize the remaining \(|\xi|\le Y\), fix any
small \(\tau>0\) and first use \(\Re\xi=x\ge\tau\).
The Gaussian model has one strictly dominant saddle
on a continued horizontal contour. Put
\[
 s_*=\frac{2}{(\xi+\sqrt{\xi^2+8})/2}.
\]
The square root has positive real part and
\(\Re s_*>0\), since \(x>0\). The saddle satisfies
\(s_*(s_*+\xi)=2\). On the horizontal through
\(s_*\), the real phase, as a function of the real
coordinate \(u\), is up to a constant
\[
            -u^2/4+\tfrac12\log((u+x)^2+B^2),\qquad
                B=\Im(s_*+\xi).
\]
For \(u>0\) this is strictly greater than its value
at \(-u\). Its derivative vanishes at precisely
one positive point: the equation
\((u+x)^2+B^2=2+2x/u\) has strictly increasing
left and decreasing right sides. Its derivative
changes strictly through zero there and the
second derivative is negative. Thus the global
maximum is the right saddle, with strict loss
on any retained exterior piece. The left half
case \(\Re\xi\le-\tau\) follows by reversing the
model coordinate. The saddle Hessian, local
transversal contours, and strict exterior losses
are uniform on these compact sets by the explicit
model formulas. Use the modulus base
\(e^{\Re A}|g_\varepsilon|\) elsewhere, including
possible base zeros. Choose a fixed large real
cutoff using the exponential tails, and move
the finite interval to the continued horizontal
saddle with strict-loss endpoint connectors.
For small \(\varepsilon\), the actual phase
and amplitude from (3) preserve these bounds.
One-saddle dominance excludes zeros.

Hence for all sufficiently large orders every
root has \(|\Im\xi|\le Y\) and \(|\Re\xi|<\tau\).
The width \(\tau\) can be any fixed positive
width chosen before narrowing the cap. No
profile symmetry is used.

## Angular separation in the regular and Airy charts

Take the effective regular and fold charts from
the [weak angular comparison](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md).
Their model losses and strict slopes near the
imaginary fields persist on a sufficiently thin
fixed strip \(|\Re\xi|<\tau\), by compactness on
the regular field intervals and by the open
fold neighborhood. Choose \(\tau\) small enough
after fixing those charts. The one-saddle
exclusions and the distinct regular intervals
overlap the fold chart as in that proof.
The exact phases are unchanged, by (2)–(3).
The regular leading amplitudes at the two
integer orders agree by (4), with the same
relative \(O(n^{-1/2})\) root equations. At the
fold, the amplitude on the analytic complex
order disk is uniformly close to the positive
Gaussian value, so the same Airy division,
exponential integer residual and derivative
bounds apply. These checks use analytic
complex order and a nonzero central amplitude
for the finite fold integral. They do not
require that the complex-order amplitude be
constant in order.

For the quantitative separation, assume two
imaginary fields differ by at most \(g n^{-5/3}\)
with \(g\) fixed small. In a regular two-saddle
chart, first move one field vertically by this
amount to a common imaginary height. The
absolute-value equations give real displacement
\(O(n^{-1})+O(g n^{-5/3})\) by the real slope,
and subtraction improves its bound to
\(O(n^{-3/2})+O(g n^{-5/3})\).
Here the vertical move changes the leading
real phase by \(O(g n^{-5/3})\), with uniform
derivative bounds. Division of the two complex
root equations then forces \(e^D=1+O(n^{-1/2})\),
against the strictly positive phase margin
for the exact saddle difference \(D\).
In the fold chart, the imaginary Airy parameter
at each root is \(O(1/n)\). Its nonzero imaginary
derivative along a horizontal segment therefore
gives the same first real displacement bound,
including \(O(g n^{-5/3})\) for the vertical move.
The scaled Airy coordinates differ by
\(O(n^{-1/3})\). The
[Airy label argument](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md#airy-labels-at-the-same-imaginary-field)
then forces the same label on a fixed smaller
chart. Its analytic order branch has imaginary
motion of size at least \(c_2 n^{-5/3}\),
with \(c_2>0\) uniform, by the same Airy inverse
derivative formula. The exponential integer
residual cannot remove this margin. Choose
\(g<c_2/2\) after increasing the terminal
threshold. Conjugation handles the negative
fold. Around zero imaginary field one uses
the regular chart across the negative physical
ray. All field ranges are covered by the
fixed chart overlaps.

Finally take \(\varepsilon Y<\pi/2\). The physical
angles all lie in the same open half-circle about
the negative ray, and the imaginary-field
distance is converted to physical angular distance
by multiplication by \(\varepsilon\). This proves
the gap for all large adjacent orders.

## Bounded orders from the signed midpoint hierarchy

For bounded \(m\) the
[matched Taylor bounds](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md#proof-of-the-curvature-hierarchy)
give, for either blended symbol \(\widehat L\),
\[
 \widehat L(1/2+x)-\widehat L(1/2)-x\widehat L'(1/2)
     =\varepsilon^2x^2+\varepsilon^4 a_2x^2
                   +\sum_{r\ge3}\varepsilon^{2r-2}a_r x^r .
\]
The signed coefficients satisfy
\(|a_r|\le 2^{r-1}/[r(r-1)]\) for \(r\ge3\), by
bounding the absolute row Taylor coefficients
using the midpoint kernel and
\(\max |d_i|/m_i\le W\).
Also \(|a_2|\le1\): the midpoint curvature is
\(W+O(W^2/j)\), with absolute error
at most \(W^2/(2j)\) from
\(0\le1-u/(e^u-1)\le u/2\).
The estimates are preserved by convex blending.
The discarded constant and linear coefficients
only change the real scalar and radial center.

Put \(x_k=k/m-1/2\). Replace \(\varepsilon\) by an
auxiliary complex \(\rho\) in the displayed
expansion and in the angular factor
\(e^{m\rho \xi x_k}\). On a fixed small
\(\rho\)-disk, the alternating binomial sum of
the resulting exponentials is analytic and
uniformly bounded for each bounded degree and
bounded \(\xi\). At weight \(w\) in \(\rho\)
the coefficient polynomial in \(x_k\) has
degree at most \(w\). Every correction past
the leading linear and quadratic factors has
strictly lower degree than its weight.
The \(m\)-th finite difference thus kills all
weights below \(m\). Its weight \(m\) coefficient
is
\[
                          (-1)^m m![v^m]e^{\xi v+v^2/m}.
\]
Cauchy's bound for the higher weights gives a
uniform \(O(\varepsilon)\) error after
division by \(\varepsilon^m\) at the physical
parameter. The limiting roots are the purely
imaginary scaled Hermite roots. Their consecutive
imaginary coordinates are strictly distinct by
[scaled Hermite interlacing](MULTIVARIATE_SCALED_HERMITE_INTERLACING.md).
For the finitely many bounded adjacent orders,
choose fixed root disks with disjoint imaginary
projections across degrees, use Rouché, and
shrink the cap once more. Each degree has its
full number of roots in those disks, which
also supplies the claimed positive gap after
possibly decreasing \(c\). The choices are effective:
first choose the model fold and regular chart margins
and the strip width, then certify the off-axis model
losses on the compact complement by interval subdivision
using the explicit saddle. Decrease the curvature cap
for these contours and choose the large terminal cutoff
from the effective density and Airy bounds. The remaining
scaled Hermite roots are algebraic and their finitely
many disks can be isolated by rational root intervals and
the bounded-order Taylor estimates.

For completeness, the same matched hierarchy
also extends every available
[paired-half Gaussian phase lift](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md#using-the-midpoint-hierarchy-in-a-paired-phase-lift)
to independent blends on the two grids.
With \(K=L_n''\), \(G=L_N''-L_n''\), and \(C=K(1/2)\),
the correction has
\(E_m''=m(K-C+\lambda_mG)\) and
\(E_N''-E_n''=K-C+(N\lambda_N-n\lambda_n)G\).
Since \(|N\lambda_N-n\lambda_n|\le N\), the
existing tolerance in that lift is unchanged.
This applies, for example, to the explicit
matched skew and increment certificate in
the strip \(16/3\le2A\le6\) there.
