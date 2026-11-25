# Stationary sheets and boundary contributions with private variables

Let \(a,b\in\mathbb Z_{\ge0}^r\) be distinct. Discard equal
rows, put \(d_i=b_i-a_i\), and define the two private degrees by
\[
 P=\sum_{a_i=0} b_i,\qquad Q=\sum_{b_i=0} a_i.
\]
The degree \(P\) controls the boundary at the pure \(a\)
monomial, and \(Q\) the boundary at the pure \(b\) monomial.
Use the ordinary factorial functional and the scaled external
parameter \(\zeta=n^\Delta t\), where \(\Delta=\sum_i d_i\).

The stationary map has exactly
\(\mathbf1_{\{P=0\}}+\mathbf1_{\{Q=0\}}\) upper-half-plane
sheets over \(-\pi<\Im X<0\). It is injective in the upper
half-plane when there is a private variable, apart from the
constant map for two private linear monomials.
On the integral side the actual boundary contributions have
explicit amplitudes. This also gives local consecutive-root
exclusion at stable real coexistence points and in every bounded
scaled window of the continuous one-sided boundary transitions.

## The degree on every open horizontal strip

Write \(x_i(q)=a_i+d_iq\), \(D_+=\sum_{d_i>0}d_i\), and
\(D_-=\sum_{d_i<0}|d_i|\). On the upper half-plane define
\[
 s(q)=\sum_i d_i\log x_i(q),\qquad
 X(q)=s(q)+\log(1-q)-\log q
\]
using principal logarithms. The rational stationary parameter
is \(\zeta(q)=e^{-X(q)}\).
Let \(s_{\rm sh}\) be the sum restricted to shared rows
\(a_i b_i>0\). The private rows give the exact identity
\[
 X(q)=C+s_{\rm sh}(q)+(P-1)\log q+(1-Q)\log(1-q)
 \tag{1}
\]
for a real constant \(C\). In particular a degree-one private
side makes the corresponding endpoint removable in \(X\).

For each integer \(j\), restrict \(X\) to the inverse image
in the upper half-plane of \(j\pi<\Im X<(j+1)\pi\).
This map is proper. Its degree, evaluated at
\(h=\Im X/\pi\notin\mathbb Z\), is
\[
 \begin{aligned}
 &\mathbf1_{\{P=0,\,-1<h<0\}}
 +\mathbf1_{\{Q\ge2,\,0<h<Q-1\}}\\
 &\qquad+\mathbf1_{\{Q-1<h<D_--1\}}
 +\mathbf1_{\{\Delta>0,\,D_--1<h<D_+-1\}} .
 \end{aligned}                                                     \tag{2}
\]
Equivalently the negative strip has the claimed degree, each
positive strip with \(0<h<\max(D_+,D_-)-1\) has degree one,
and the remaining strips have degree zero. In particular every
nonempty inverse image of a strip of degree one is a single
component mapped conformally onto that strip.

Here is the argument-principle proof. Along the real axis
the height \(\Im X/\pi\) is an integer. On \((0,1)\) it is
zero. Immediately to the left of zero it is \(P-1\),
with each shared left pole contributing its positive gap
to the height on the side farther left. Immediately to the
right of one it is \(Q-1\), and each shared right pole
adds its absolute gap as the abscissa increases.
Use a large upper half-disk with small semicircular
indentations around real logarithmic singularities. Count
upward crossings of the ray \(\Im X=\pi h,\ \Re X>\Re w\)
for a target \(w\) at that height. Straight segments have
no crossing. On an indentation the real part escapes with
the sign of its logarithmic divergence and the imaginary
part interpolates monotonically between the endpoint heights
for small radii whenever those heights straddle \(h\).
This follows from the affine angular term of the singular
logarithm and a uniformly vanishing analytic error, also
after angular differentiation.

At zero the real part escapes right only for \(P=0\),
giving the first term of (2). At one it escapes right
only for \(Q\ge2\), giving the second term. It escapes
left at every shared left pole, and right at every shared
right pole. The successive elementary height intervals
of the shared right poles together give the third term.
Finally, on \(q=R e^{i\theta}\) the asymptotic is
\[
 X(q)=\Delta\log R+C_\infty
          +i\{\Delta\theta+\pi(D_--1)\}+O(R^{-1}),
 \tag{3}
\]
with \(C_\infty\) real. The error and its angular derivative
are uniform, by factoring out \(q\) in the logarithms.
The arc contributes only for \(\Delta>0\), which gives
the last term. When \(\Delta=0\), the limiting height is
an integer and there is no crossing for a noninteger \(h\).
This proves (2). If the image stays in a compact subset
of the strip, finite boundary approaches are excluded by
their integer heights or divergent real part, and escape
to infinity is excluded by (3). This proves properness.

If \(P+Q>0\), every strip degree is at most one.
The flat profile \(P=Q=1\) without varying shared rows
has \(X\) constant. In every other private profile \(X\)
is nonconstant, by (1) and the locations of the shared
poles. It is injective on the whole upper half-plane.
Indeed two distinct preimages of a target on an
integer-height line would give overlapping open images
of disjoint neighborhoods, and hence two preimages of
a nearby noninteger-height target. A critical point is
also excluded by local ramification. Conjugation
excludes nonreal zeros of the rational function \(X'\).
For \(P=Q=0\), (2) recovers the degree in the
[shared-support argument](MULTIVARIATE_SADDLE_MAP_TWO_SHEETS.md).
It also locates its unique possible nonreal critical point
in the negative strip: local ramification cannot occur
on an integer-height boundary or in a strip of degree
less than two.

The conditional strip \(|\Im s|<\pi/2\) remains valid
with private rows. The proof by the vertical domain
\(\Re x_i(q)>0\) is exactly the same as in the
[shared-support conditioning argument](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#an-exact-conditional-integral).
At a finite side with a nonzero imaginary part, including
the new sides \(0,1\), a touching row contributes absolute
height at least \(\pi/2\). At the real side its logarithm
makes the real part escape. On the upper part of the
conditional strip, \(\Im X\in(-\pi,\pi/2)\).
Thus for upper-half-plane external \(\zeta\) the only
possible lift is \(X=-\operatorname{Log}\zeta\).
In particular if both sides have private variables there
is no upper conditional saddle for such a \(\zeta\).
This count leaves the actual contributing contours to be
identified away from the real regimes below.

## The global image of the conditional coordinate

For every profile, the function \(s(q)\) itself is
injective on the upper half-plane. Its image is the
strip
\[
 0<\Im s<\pi\max(D_+,D_-)
\]
with finitely many horizontal closed rays removed.
Every bounded gap between successive left poles
supplies a leftward ray, and every bounded gap between
successive right poles supplies a rightward ray.
The interval beyond the leftmost pole supplies another
leftward ray if \(D_+>0\) and \(\Delta<0\). The interval
beyond the rightmost pole supplies a rightward ray if
\(D_->0\) and \(\Delta>0\). The rays have the boundary heights
of the corresponding real intervals, all integer
multiples of \(\pi\). Two opposite rays at the same
height have disjoint closures.

To prove injectivity, perform the same winding count
without \(\log(1-q)-\log q\). All left poles, including
zero if \(P>0\), escape left. The right poles, including
one if \(Q>0\), escape right and contribute precisely
for \(0<\Im s/\pi<D_-\) at noninteger heights.
The outer arc has limiting heights \(D_-\) and \(D_+\)
and contributes precisely the adjacent interval
\(D_-<\Im s/\pi<D_+\) if \(\Delta>0\).
Thus the degree at noninteger heights is one throughout
the stated strip and zero elsewhere. The open mapping
argument gives injectivity also on the integer levels.

On an internal integer level, a finite real-boundary
value cannot be in the image. Otherwise continuity of
the inverse at that target would give an interior
limit of the preimages approaching the real boundary.
Conversely take a missing target on an internal level
and approximate it by targets of noninteger heights.
Their preimages have a real limit along a subsequence:
an interior limit would put the target in the image,
a pole would send the real part to infinity, and escape
to infinity would send the real part to infinity or,
when \(\Delta=0\), reach the top height. Thus the
missing targets are precisely the finite real-boundary
values at that height.

On each bounded left gap the real part of \(s\) tends
to \(-\infty\) at both ends, and on each bounded right
gap it tends to \(+\infty\). The same respective behavior
holds on the two stated far intervals. The image of
each such gap along the real boundary is therefore a
closed ray starting at its finite extremum, by the
intermediate value theorem and attainment of the
extremum in the interval. These are exactly the
real intervals at internal heights. If opposite rays
at the same height touched, their union would remove
the whole horizontal line and separate the image into
two nonempty components. This contradicts connectedness.
The saddle-centered integral with a positive transverse
quadratic form applies on \(s({\cal V})\), where
\({\cal V}=\{q:\Re x_i(q)>0\text{ for every }i\}\).

The conditional density itself has a larger continuation.
Let \(v_i=\operatorname{sign}(d_i)/\sum_j |d_j|\).
For real \(s\), set \(\log y_i=v_i s+u_i\) with
\(u\in d^\perp\). Since \(d\cdot v=1\), the log-coordinate
Jacobian is \(1/\|d\|\), and the real conditional density
is
\[
 \frac1{\|d\|}\int_{d^\perp}
  \exp\!\left(
   \sum_i(na_i+1)(v_i s+u_i)
             -n\sum_i e^{v_i s+u_i}\right)\,du .
\]
This is holomorphic throughout
\(|\Im s|<(\pi/2)\sum_i |d_i|\).
On compact subsets each \(\cos(v_i\Im s)\) is bounded
below by a positive constant. At \(u_i\to+\infty\)
the negative exponential dominates, and at
\(u_i\to-\infty\) the coefficient \(na_i+1>0\)
gives an integrable linear decay. Analytic continuation
identifies this density with the saddle-centered formula
on \(s({\cal V})\). A complex transverse saddle with
some \(\Re x_i\le0\) still requires a deformation of its
integration cycle. For integer \(n\) the scalar moment
contour can be shifted to any horizontal line in this
full strip. Indeed rotate the original positive \(y_i\)-ray
by \(v_i h\) with \(h=\Im s\). Each rotation preserves the
integral term by term, since the integrand is a polynomial
times \(e^{-n\sum y_i}\) and \(|v_i h|<\pi/2\).
In the log coordinates these rotations give precisely
the affine conditional formula on \(s\in\mathbb R+ih\).

## Actual boundary contributions

Let \(a\) denote the boundary \(y_i=a_i\) after the
original rescaling \(X_i=ny_i\), and assume \(P>0\).
Write \(r_a=\#\{i:a_i>0\}\) and
\(K_a=\prod_{a_i>0}a_i^{d_i}\).
For a small fixed neighborhood of that boundary, its
local contribution to the unscaled moment is
\[
 \Lambda(X^{na})\bigl(g_a(\zeta)+O(n^{-1})\bigr),\qquad
 g_a(\zeta)=
 \begin{cases}
  1,&P\ge2,\\
  (1-\zeta K_a)^{-1},&P=1 .
 \end{cases}                                                      \tag{4}
\]
For \(P\ge2\) the estimate is uniform on compact sets of
\(\zeta\). For \(P=1\) use compact subsets of
\(\Re(\zeta K_a)<1\). The \(b\) boundary has the reciprocal
formula with \(Q,\zeta^{-1},K_b=\prod_{b_i>0}b_i^{-d_i}\)
and pure factor \(t^n\Lambda(X^{nb})\).

For the proof, tilt by the pure \(a\) moment. The active
variables \(X_i/n\) converge to \(a_i\) under their Gamma
laws of shape \(na_i+1\). The other \(X_i\) are independent
unit exponentials and
\[
 Y_n=\prod_i(X_i/n)^{d_i}
   =n^{-P} \prod_{a_i=0}X_i^{b_i}
                \prod_{a_i>0}(X_i/n)^{d_i}.
\]
The multiplier is \((1+\zeta Y_n)^n\).
For \(P\ge2\), conditional on bounded private variables,
it tends to one. With \(S_0=\sum_{a_i=0}X_i\), on
\(S_0/n\le\eta\) and the active neighborhood,
\(n|Y_n|\le C\eta^{P-1}S_0\); take \(\eta\) small for
integrable exponential domination. For \(P=1\) there is one
private variable \(X_0\) with exponent one. The limit is
\(e^{\zeta K_a X_0}\). The bound
\(\log|1+z|\le\Re z+C|z|^2\) for small \(z\)
gives uniform domination with a strictly positive decay
margin against its exponential density. This yields (4)
at leading order.

Both stable boundary packets also have an analytic
\(1/n\) correction with relative error \(O(n^{-2})\).
To see this, work before dividing by the pure moment.
Use \(X_i/n=a_i+\epsilon u_i\) for the active variables,
\(\epsilon=n^{-1/2}\), and leave the private \(X_i\)
unscaled. The active exponent has a nondegenerate
Gaussian expansion, because
\(\sum_{a_i>0}(a_i\log y_i-y_i)\) has its strict
maximum at \(a\). In the private multiplier write
\[
 n\log(1+\zeta Y_n)
 =\epsilon^{-2}\log\!\left(
      1+\epsilon^{2P}\zeta W(a+\epsilon u,X_{\rm private})
                    \right).
\]
Here \(W\) is analytic in the active coordinates and a
fixed monomial in the private coordinates. The apparent
factor \(\epsilon^{-2}\) is removed by factoring
\(\log(1+w)=w\cdot(\log(1+w)/w)\). Take smooth cutoffs
supported where the active variables are near \(a\) and
\(S_0/n\le\eta\). There \(w\) stays in a small disk,
and derivatives through order four after scaling are
bounded by polynomials in the scaled variables times
the same Gaussian and private exponential majorant.
Choose the active cutoff even in \(y-a\).
Odd powers of \(\epsilon\) integrate to zero by
\(u\mapsto-u\); the private cutoff is even in
\(\epsilon\). Taylor expansion therefore gives the
\(1/n\) term with error \(O(n^{-2})\).
The same argument for the pure \(a\) integral, or
Stirling's expansion, gives the correction relative to its
exact pure moment. The argument for the other boundary
uses the reversed monomials.

In the scaled integral \(I_n\) the \(a\) contribution
therefore has exponential phase
\(A_a=\sum_{a_i>0}a_i\log a_i-\sum a_i\)
and power \(n^{r_a/2-r}\).
The \(b\) contribution has phase
\(\log\zeta+A_b\) and power \(n^{r_b/2-r}\).
These powers can differ from one another and from the
power \(n^{-r/2}\) of a nondegenerate interior saddle.

## Real geometry and consecutive comparison at coexistence

For real \(q\in(0,1)\) the curvature function is
\[
 \rho(q)=P(1-q)+Qq
   +q(1-q)\sum_{a_i b_i>0}\frac{d_i^2}{x_i(q)}.
 \tag{5}
\]
It is concave, strictly so if there is a varying shared
row. As before
\((\log\zeta(q))'=(1-\rho(q))/(q(1-q))\).
For real positive \(\zeta\), the derivative of the
conditional phase with respect to \(s\) has the sign of
\(\log\zeta-\log\zeta(q)\) on \(0<q<1\).
This follows by comparing the log odds of
\(\zeta e^s/(1+\zeta e^s)\) and \(q\).
Where \(q\le0\) the phase derivative is positive, and
where \(q\ge1\) it is negative.

When both \(P,Q\ge1\), the interior stationary points
are minima, except for the flat case \(P=Q=1\) without
varying shared rows. In the nonflat case the pure boundary
rates cross once at \(\log\zeta=A_a-A_b\).
When \(P=1\), the endpoint stability threshold is
\(\zeta_a^*=1/K_a\), obtained by canceling the log at zero
in (1). Analogously when \(Q=1\) it is
\(\zeta_b^*=K_b\). The crossing is strictly in both
stability ranges where those restrictions apply.
For instance at the finite scale \(\zeta_a^*\) the
function \(\zeta(q)\) is strictly decreasing from its
endpoint value, so the phase is strictly increasing
through \((0,1)\) and the \(b\) boundary is higher.
At the other finite threshold the reasoning reverses.
In the flat case the two varying monomials are
independent linear variables, whose moment polynomials
are truncated geometric series times \(n!\), so the
consecutive orders are coprime directly.

For one private side, orient it so \(P>0,Q=0\).
If \(P\ge2\), or if \(P=1\) and \(\rho'(0)>0\),
there is exactly one interior fold by (5). Beyond
the turning value there are a local minimum and a
nondegenerate interior maximum, in increasing \(q\)
order. The maximum and the \(a\) boundary exchange
dominance once: the difference of their phase values
has derivative \(q_{\rm int}>0\) with respect to
\(\log\zeta\). At the fold the incipient saddle
is strictly below the boundary, as the phase is
strictly decreasing between them. For \(P\ge2\),
at large enough \(\zeta\) the interior dominates,
which follows by testing any fixed positive point.
For \(P=1\) and \(\rho'(0)>0\), at
\(\zeta=\zeta_a^*\) the phase increases from the
boundary, so the unique crossing is strictly inside
the stability range of (4).
If \(P=1,Q=0,\rho'(0)\le0\), there is instead a
continuous transition at \(\zeta_a^*\): below this
scale the boundary is the unique maximum and above
it the interior maximum is unique by monotonicity
of \(\zeta(q)\).

At each of the coexistence points with two stable
contributions, consecutive moments have no common
zero in a sufficiently small fixed complex neighborhood
of the positive real scale, for all large orders.
Here is the comparison including the different powers of
\(n\). Denote the two leading terms by
\[
 n^{-\mu_j} A_j(\zeta)e^{n\Phi_j(\zeta)}
       (1+B_j(\zeta)/n+O(n^{-2})),\qquad j=1,2.
\]
Their amplitudes are positive on the positive real
coexistence point. For an interior saddle use the
localized conditional contour of the
[shared-support analysis](MULTIVARIATE_SHARED_SUPPORT_CONDITIONAL_SADDLE.md#a-local-exclusion-around-the-positive-real-axis);
for a boundary use (4) and its first correction.
The rest of the real integral has a strict gap.
Indeed the phase has compact upper level sets on the
closed orthant where at least one monomial is nonzero,
using logarithmic growth at infinity. Excise fixed
neighborhoods of its distinct global maxima. On the
remaining region the real integral is bounded by a finite
one-order integral times an exponentially smaller
factor. This also bounds it for nearby complex \(\zeta\)
using \(|y^a+\zeta y^b|\le
(1+|\zeta-\zeta_0|/\zeta_0)(y^a+\zeta_0y^b)\).
Shrink the disk so the connecting segments of any
interior local contour retain this gap. These estimates
justify the two-term formula with the stated errors.

Set \(\delta=\Phi_2-\Phi_1\). At a zero,
\(e^{n\delta}=-n^{\mu_2-\mu_1}(A_1/A_2)(1+O(n^{-1}))\).
At the consecutive order use
\(\zeta'=(1+1/n)^\Delta\zeta\) and divide the equations.
For \(\Delta\ne0\) the exponential limit on the left at
the coexistence point is
\(\exp(\Delta(q_2-q_1))\), where an \(a\) boundary
has \(q=0\) and a \(b\) boundary has \(q=1\).
It differs from one because \(q_2\ne q_1\).
For \(\Delta=0\), division gives
\[
 e^\delta=(1+1/n)^{\mu_2-\mu_1}(1+O(n^{-2})).
\]
Thus \(\delta=O(1/n)\) in a sufficiently small disk.
If the powers differ, taking absolute values in the
first zero equation instead requires
\(n\Re\delta=(\mu_2-\mu_1)\log n+O(1)\).
If they agree, division gives \(\delta=O(n^{-2})\);
as \(\delta\) has a simple zero at coexistence,
the first zero equation then has left side tending
to one and right side tending to a negative value.
All cases exclude a simultaneous zero.

## Order increments at complex boundary cancellations

The same comparison applies conditionally at a complex scale
where a regular two-term expansion with the actual boundary
and interior contributions is known to exhaust the integral.
For an \(a\) boundary and an interior saddle in \(s({\cal V})\),
equal leading moduli and equal consecutive multipliers when
\(\Delta\ne0\) would force \(\Re q=0\).
Indeed the boundary phase is independent of \(\zeta\), while
the increment of the interior phase gains \(\Delta q\).
But \(P>0\) implies that some \(x_i=d_i q\) with \(d_i>0\),
so \(\Re q>0\) in \({\cal V}\).
With a \(b\) boundary the analogous requirement
\(\Re q=1\) is incompatible with \(\Re q<1\).

For \(\Delta=0\), a boundary and an interior saddle
always carry unequal powers of \(n\), by (4) and the
interior Gaussian power. Adjacent cancellation would
require \(\Re\delta=O(1/n)\), even if the imaginary
part of \(\delta\) tends to a nonzero multiple of
\(2\pi\). The unequal powers require
\(n\Re\delta=(\mu_2-\mu_1)\log n+O(1)\), ruling this out.
For two stable boundaries the multiplier moduli differ
when \(\Delta\ne0\). If \(\Delta=0\), equality of the
leading multipliers forces the positive real crossing
\(\zeta=e^{A_a-A_b}\), already treated above.
These tests require a two-contribution expansion in the
first place. They do not establish which contributions
occur on a general complex contour.

## Continuous endpoint windows

Suppose \(P=1,Q=0,\rho'(0)\le0\).
Orient the sole \(b\)-private variable as coordinate zero,
so \(x_0(q)=q\) and \(d_0=1\). Write
\(\beta=\log\zeta\), \(\beta_0=-\log K_a\), and
\(\beta(q)=\log\zeta(q)\). The finite log-scale at zero
has the expansion
\[
 \beta(q)=\beta_0+c_1 q+c_2 q^2+c_3 q^3+O(q^4),\qquad
 c_j=\frac{1-(-1)^{j+1}\sum_{i>0}d_i^{j+1}/a_i^j}{j}.
 \tag{6}
\]
Indeed \(\beta(q)=-\log(1-q)-\sum_{i>0}d_i\log x_i(q)\).
Here \(c_1=-\rho'(0)\ge0\); if \(c_1=0\), then
\(c_2=-\rho''(0)/4>0\) by strict concavity.
Let \(m=2\) if \(c_1>0\) and \(m=3\) if \(c_1=0\);
put \(\ell=c_{m-1}^{-1/m}\).
For every fixed \(H<\infty\), the consecutive moments
have no common zero \(t\) with
\[
 \left|\ell n^{(m-1)/m}
           \bigl(\log(n^\Delta t)-\beta_0\bigr)\right|\le H
 \tag{7}
\]
for all sufficiently large \(n\). The logarithm is the
one near the positive scale. Reversing the monomials
gives the corresponding \(b\) boundary result.

### Uniform half-line limits

Use \(q\) as the integration variable near zero.
In the exact conditional integral, the Jacobian cancels
the vanishing product factor:
\[
 s'(q)\prod_i x_i(q)
   =\bigl(1+q s_{\rm sh}'(q)\bigr)\prod_{i>0}x_i(q).
\]
Here \(s(q)=\log q+s_{\rm sh}(q)\) since \(Q=0\).
The phase
\[
 \phi(q,\beta)
  =\log(1+e^\beta q\,e^{s_{\rm sh}(q)})
        +\sum_{i>0}a_i\log x_i(q)-\sum_i x_i(q)
\]
is analytic at zero, with \(\phi(0,\beta)=A_a\).
For \(q>0\) its derivative equals
\(s'(q)(L_q-q)\), where the log odds of \(L_q\)
exceed the log odds of \(q\) by \(\beta-\beta(q)\).
Therefore
\[
 \phi(q,\beta_0)-A_a
             =-\frac{c_{m-1}}m q^m+O(q^{m+1}),\qquad
 \phi_{q\beta}(0,\beta_0)=1.
 \tag{8}
\]
For example the prefactor in the linearized
log-odds difference is
\(s'(q)q(1-q)=(1+q s_{\rm sh}'(q))(1-q)\).
When \(c_1=0\) this has no term linear in \(q\).

The transverse quadratic form
\(\operatorname{diag}(x_i(q))|_{\{\sum_i d_i u_i=0\}}\)
is positive definite through \(q=0\): a vector supported
on coordinate zero cannot lie in that hyperplane.
The local amplitude after the transverse Gaussian
integration is analytic and nonzero to leading order at
zero, proportional to
\[
 A(q)=
 \frac{s'(q)\prod_i x_i(q)}
 {\sqrt{\det(\operatorname{diag}(x_i(q))|_{\{\sum_i d_i u_i=0\}})}} .
 \tag{9}
\]
Only the local transverse integral is needed at zero.
Its smooth cutoff can be chosen even on the hyperplane.
After scaling \(u=\tau v\), its normalized amplitude
is a smooth even function of real \(\tau=n^{-1/2}\),
including zero, uniformly for \(q\) near zero on the
real half-line. Indeed all fixed derivatives after the
scaling factor out the quadratic term and are bounded
by polynomials times a uniform Gaussian. The cutoff
localizes in \(u\), and the real part of the transverse
quadratic form is uniformly positive there. The evenness
follows from changing \(v\) to \(-v\). Thus the transverse
corrections start at order \(\tau^2=1/n\).

Now put \(\epsilon=n^{-1/m}\),
\(q=\ell\epsilon z\) and
\(\beta=\beta_0+\epsilon^{m-1}h/\ell\).
Normalize by the boundary exponential, the leading
Gaussian amplitude, and the Jacobian powers. The resulting
period \(J_\epsilon(h)\) tends to
\[
 F_m(h)=\int_0^\infty e^{-z^m/m+hz}\,dz .
 \tag{10}
\]
The normalized period has a real \(C^3\) extension in
\(\epsilon\) to zero, locally uniformly in complex \(h\)
with its holomorphic derivatives. For details, choose
fixed smooth cutoffs near \(q=0,u=0\).
After rescaling, the normalized base phase has the
factor \(z^m\) times a smooth quotient with negative real
part; the parameter difference has the factor \(hz\)
times an analytic quotient. Thus the derivatives in
\(\epsilon\) through order three of the local integrand
are dominated by polynomials times
\(e^{-c z^m+C|z|}\) after transverse integration.
For the latter, the even smooth expansion in \(\tau\)
can be taken through any fixed number of terms before
inserting \(\tau=\epsilon^{m/2}\), so this substitution
also has the asserted derivatives. The cutoff endpoints
give exponentially small terms. The complement of these
cutoffs has an exponential gap at \(\zeta_a^*\), since
the boundary at \(a\) is its unique maximum. The local
coordinates \((q,u)\) extend regularly on the real
half-space at that boundary: \(y_0=q e^{u_0}\), and
the active logarithmic coordinates and the relation
\(u_0=-\sum_{i>0}d_i u_i\) determine the other variables
with a nonsingular Jacobian. Thus the excised local
neighborhood leaves a genuine gap in the original orthant.
Differentiation of the normalized original remainder
introduces only powers of \(\epsilon^{-1}\) and powers
of the base phase, which are integrable against the
one-order exponential density. This proves the same
smoothness for the complete period.

### Consecutive exclusion in the scaled windows

For \(m=2\), integration by parts gives \(F_2'=hF_2+1\),
so every zero is simple.
At any simple model zero \(h_0\) in either model,
the actual branch has \(h(\epsilon)=h_0+O(\epsilon)\)
with bounded derivative. Its log scale
\(\beta(n)=\beta_0+n^{-(m-1)/m}h(n^{-1/m})/\ell\)
then satisfies
\[
 \beta'(n)=-\frac{m-1}{m\ell}\,h_0 n^{-(2m-1)/m}
                    +O(n^{-2}).
 \tag{11}
\]
Here \(h_0\ne0\) since \(F_m(0)>0\).
Under a common root the exact change in the log scale
would be \(\Delta\log(1+1/n)\).
The change in \(h\) is \(O(n^{-1/m})\), so the two roots
must approach the same model zero. In a bounded window
this follows from Rouché's theorem and separation of the
finitely many model zeros. For simple branches (11)
contradicts the exact change when \(\Delta\ne0\) by
size, and when \(\Delta=0\) by its nonzero leading term.

For \(m=3\) one can also handle a hypothetical multiple
zero without using any global simplicity statement about
the Scorer function. The identity \(F_3''=hF_3+1\)
shows that such a zero is exactly double. Suppose \(h_0\)
is one. Then
\[
 F_3(h_0)=F_3'(h_0)=0,\quad
 F_3''(h_0)=1,\quad F_3'''(h_0)=0,\quad
 F_3^{(4)}(h_0)=h_0,\quad F_3^{(5)}(h_0)=3.
\]
The order-\(\epsilon\) correction in \(J_\epsilon\) is
\[
 \ell\,\frac{A'(0)}{A(0)}F_3'
                       -\frac{c_3\ell^4}{4}F_3^{(4)}.
 \tag{12}
\]
Indeed the next base phase term in (8) is
\(-c_3 q^4/4\); the mixed \(h z^2\) term vanishes
by the observation after (8).
The transverse correction starts at \(\epsilon^3\).
Thus if \(c_3\ne0\), (12) has a nonzero value at
\(h_0\). Substituting \(h=h_0+\sqrt\epsilon\,z\)
and dividing \(J_\epsilon\) by \(\epsilon\) yields a
quadratic limit with two distinct nonzero roots.
The implicit function theorem in \(\sqrt\epsilon\)
gives branches \(h=h_0+O(\sqrt\epsilon)\) with
\(dh/d\epsilon=O(\epsilon^{-1/2})\).
Their separation is of order \(\sqrt\epsilon\),
larger than the adjacent-order field change
\(O(\epsilon)\). Hence a common root must belong to
the same branch. Its derivative in (11) retains the
same leading term with error \(O(n^{-11/6})\), again
contradicting the exact change for either value of
the degree difference. These branch constructions
use the real \(C^3\) extension established above,
so derivatives of an unproved error term are not
needed.

The only exception \(c_3=0\) is elementary to identify.
The cubic condition says \(\sum_{i>0}d_i^2/a_i=1\).
Since
\(d_i^4/a_i^3=(d_i^2/a_i)^2/a_i\le d_i^2/a_i\),
equality in \(c_3=0\) requires exactly one term,
\(a_i=1\) and \(|d_i|=1\).
The assumption \(Q=0\) gives \(d_i=1\).
The profile is therefore \(X\) versus \(X^2Y\), up to
discarded equal rows. Here \(\Delta=2,\ \ell=1\), and
the exact reduced phase and leading amplitude are
\[
 \phi(q,\beta)-A_a
    =\log(1+q)+\log(1+e^\beta q(1+q))-2q,\qquad
 A(q)/A(0)=\sqrt{1+2q}.
\]
Using \(q=\epsilon z,\ \beta=\epsilon^2h\), the
normalized expansion through order two is
\[
 \begin{aligned}
 J_\epsilon(h)={}&F_3(h)+\epsilon F_3'(h)\\
 &+\epsilon^2\left[
  -\tfrac12 F_3''+\tfrac12h^2 F_3'
              -hF_3'''+\tfrac25 F_3^{(5)}\right](h)
        +O(\epsilon^3).
 \end{aligned}
\]
At a double zero \(h_0\), substitute \(h=h_0+\epsilon z\)
and divide by \(\epsilon^2\). The limit is
\(z^2/2+z+7/10\), whose two roots are
\(-1\pm i\sqrt{2/5}\).
An adjacent-order field change at fixed \(t\) is
\(2\epsilon+o(\epsilon)\) in this scaled field.
Differences of the limiting \(z\) offsets are zero or
purely imaginary, so none can match this real change.
This covers the exceptional profile and completes (7).

The sheet count and the real and boundary results are
for fixed profiles. They do not identify every complex
saddle contribution and do not give an all-order
coprimality assertion for general profiles.
