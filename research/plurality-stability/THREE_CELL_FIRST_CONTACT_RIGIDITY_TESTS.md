# Shape restrictions at three-cell first contact

Let \(A\) be a balanced three-cell Gaussian maximizer at \(0<r<1\)
and \(T\) the centered regular three-sector partition. Work in
the active plane and put \(D(s)=S_s(T)-S_s(A)\) for the fixed
partition \(A\). The arguments use the regularity, complete
adjacency and full rank of the cell moments of actual maximizers.
We use \(U_i=T_r1_{A_i}\), prices \(\beta_i\), and adjusted scores
\(U_i-\beta_i\).
[Maximizer structure and interface forms](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

## Reflected ties have at most second-order contact

Suppose \(A\) has a central reflection exchanging two labels and
is nonsimplex. If \(A\) ties \(T\) at a correlation \(r\) approached
by distinct simplex-optimal correlations and \(D'(r)=0\), then
\(D''(r)>0\). The bound is uniformly positive on a compact
correlation band over nonsimplex ties with these hypotheses.

Use the [reflection reduction](THREE_CELL_REFLECTION_REDUCTION.md)
with \(A_3=\{y>h(|x|)\}\) and \(A_1=\{x>0,\ y<h(x)\}\).
Write \(F(x,y)=U_3(x,y)-U_1(x,y)\), \(\sigma=\sqrt{1-r^2}\)
and \(\beta=\beta_3-\beta_1\). On the right graph \(F=\beta\).
If \(D''(r)=0\), the mixed inequality (G10) forces
\(g(x)=\partial_r F(x,h(x))\) to be constant there.
Indeed one can test any compact speed of zero Gaussian mean
on this pair interface, with zero speed on the others. For
either direction of approach in (G10), use the signed
correlation increment \(\epsilon_j=s_j-r\) and a shape
displacement \(c\epsilon_j f\) with arbitrary real \(c\).
The quadratic inequality after division by \(\epsilon_j^2\)
is the same. The derivative in \(g\) is at fixed observation
and fixed cells.
We prove that constancy of \(g\) forces the simplex.

Let \(\phi,\Phi\) be the standard normal density and distribution,
and let \(p_\sigma\) be the density with variance \(\sigma^2\).
Conditioning in the horizontal source coordinate gives
\[
 F(x,y)=-m_0(x)+I(x,y),\qquad m_0(x)=\Phi(rx/\sigma),\qquad
 I(x,y)=\int p_\sigma(t-rx)c(t)
     \Phi((ry-h(|t|))/\sigma)\,dt,\quad c(t)=1+1_{\{t>0\}}.
\]
The vertical derivative is everywhere positive. The limits of
\(F(x,y)\) as \(y\to-\infty,+\infty\) are \(-m_0(x),1\)
respectively. Thus the level \(F(x,y)=\beta\) extends to the
analytic graph \(y=Y(x)\) on \(\{x:m_0(x)>-\beta\}\).
Analyticity in \(x\) continues the constant \(\partial_r F\)
along the extended graph.

If \(\beta<0\), the left endpoint is the finite \(x_c\) with
\(m_0(x_c)=-\beta\), and \(Y(x)\to-\infty\) as \(x\downarrow x_c\).
At fixed \(x,y\), differentiation of \(p_\sigma(t-rx)\) in \(r\)
has logarithmic derivative a quadratic in \(t-rx\) with leading
coefficient \(-r/\sigma^4\). It is bounded above for \(x\) near
\(x_c\). Write \(L=(ry-h(|t|))/\sigma\). Since \(h\ge h(0)\),
for large negative \(y\) both \(L\) and
\(\partial_r L=(y-rh(|t|))/\sigma^3\) are bounded above by
\(-c|y|\) for a constant \(c>0\). The Gaussian lower-tail bound
\(\phi(L)/\Phi(L)\ge |L|\) follows by integrating
\((-u)\phi(u)\) for \(u\le L<0\). Hence
\(I_r/I\le C-c' y^2\). Also \(I_r\to0\) along this limit.
For the horizontal derivative use Gaussian domination. For the
vertical derivative use \(y-rh(|t|)=\sigma^2 y+r\sigma L\) and
\(|y|\le C(1+|L|)\) on the negative tail.
Constancy of \(g\) now gives \(g=-\partial_r m_0(x_c)\), so
\[
 I_r(x,Y(x))=\partial_r m_0(x)-\partial_r m_0(x_c)
             =O(x-x_c)=O(I(x,Y(x))).
\]
This contradicts \(I_r/I\to-\infty\).

If \(\beta>0\), the level extends to every real \(x\).
Use the positive-tail proof of (RS5)–(RS6) in the
[reflected shape analysis](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md).
For constant \(g\) it gives \(h(t)=m t+b+o(1)\),
\(t(h'(t)-m)\to0\), and
\[
 b=\frac{\sigma\sqrt{1+m^2}\,z}{r-1},\qquad
 g=-\frac{2z\phi(z)}{\sigma^2},\qquad z=\Phi^{-1}((1+\beta)/2).
\]
Put \(u=\Phi^{-1}(\beta)\). At \(x=-t\to-\infty\), the contribution
of the positive horizontal source, including its \(r\) derivative
at \(y=m t+O(1)\), is exponentially small up to a polynomial
factor in \(t\). Writing \(S=-rt+\sigma Z\) with \(Z\) standard
normal, the remaining conditional argument is
\(L=(ry-h(-S))/\sigma\) on \(S<0\).
For \(y=m t+c\) its smoothed limit is
\(\Phi((rc-b)/(\sigma\sqrt{1+m^2}))\), locally uniformly in \(c\).
Monotonicity therefore yields
\(Y(-t)=m t+c_*+o(1)\), where \(rc_*-b=\sigma\sqrt{1+m^2}\,u\).
Differentiation at fixed \(t,y\) gives
\[
 \partial_r L=\frac{y+(-t-rZ/\sigma)h'(-S)}{\sigma}
                  +\frac{rL}{\sigma^2}.
\]
On \(|Z|\le c_0t\) choose \(c_0>0\) so that \(-S\ge rt/2\).
Then \(t(h'(-S)-m)\to0\) uniformly, and outside this event
the Gaussian tail controls the derivative using the bounded
slope. In the limit \(L\to\sqrt{1+m^2}\,u+mZ\) and the \(Z\)
terms in \(\partial_r L\) cancel. Integration gives
\[
 g=\frac{\phi(u)}{r\sigma^2}\bigl[u-(1+r)z\bigr].
\]
Comparing the two expressions for \(g\), and writing
\(H(v)=\phi(v)/(1-\Phi(v))\), yields
\[
 z-u=rz\left(\frac{H(z)}{H(u)}-1\right).
\]
Here \(u<z\), \(z>0\), and \(H(z)-H(u)<z-u\) because
\(H'(v)=1-\operatorname{Var}(Z\mid Z>v)<1\), by integrating the
first and second moments of a truncated normal. Moreover
\(H(u)>z\). In fact \(z\) is the median of \(Z\mid Z>u\).
At displacements \(0<t<z-u\) the left density \(\phi(z-t)\)
is larger than the right density \(\phi(z+t)\); beyond \(z-u\)
only the right density remains. The signed difference of the
right and left densities integrates to zero by the median
condition. Weighting it by \(t\) makes its integral strictly
positive, since one may subtract \(z-u\) from the weight.
This proves the asserted mean inequality. The right side of
the displayed equation is thus strictly less than \(r(z-u)\),
a contradiction.

For \(\beta=0\), (RS5) gives \(g=0\) and (RS7) gives the
simplex. All three cases show that a nonsimplex reflected
maximizer has nonconstant \(g\). At a tangential tie,
Taylor expansion along the approaching optimal correlations
gives \(D''(r)\ge0\), so it is strictly positive.

To see uniformity, rotate and relabel the reflected maximizers
to fix the reflection and label transposition. Strong maximizer
compactness on the band and uniform local simplex isolation
keep subsequential limits maximizing, reflected and nonsimplex.
The limits remain simplex-optimal accumulation points: if
correlations in the sequence take infinitely many distinct
values, use those correlations; if eventually constant along a
subsequence, use the corresponding hypothesis. The Hermite
multipliers give joint continuity of \(D'\) and \(D''\) under
strong convergence on a compact band inside \((0,1)\). Thus
a sequence with \(D''\to0\) would produce the excluded limit.
[Compactness and local comparison](POSITIVE_CORRELATION_CONTINUATION.md)

## Semidefinite affine flows

A nonsimplex balanced planar maximizer is invariant under no
flow \(x\mapsto p+e^{tB}(x-p)\) with \(B\) a nonzero symmetric
semidefinite matrix. Invariance is understood up to null sets
for all real \(t\). Reversing time permits \(B\ge0\).
A positive scalar \(B\) is excluded by the conical classification
about arbitrary centers in the maximizer analysis.

For unequal positive eigenvalues, use coordinates \((x,y)\)
about \(p\) in the eigenbasis and write \(q=\lambda_2/\lambda_1>1\).
On each sign half-plane of \(x\), the label depends only on
\(y/|x|^q\). This follows by changing variables to
\((\log|x|,y/|x|^q)\) and using invariance under all translations
of the first coordinate. On a generic vertical slice the
analytic adjusted pair contrasts have discrete zeros. Such a
slice exists since each contrast is nonzero as an analytic
function, and for each pair one can fix a horizontal slice
on which it is nonzero as a function of \(x\). Its discrete
zero set contains every exceptional vertical coordinate.
Consequently the profiles in \(y/|x|^q\) have locally discrete
label transitions.

Each transition \(c\) yields a two-label interface arc
\(y=c|x|^q\) tending to the origin in these coordinates.
By (G20) its pair score gradient does not vanish at the apex.
The tangent there is horizontal, so the analytic implicit
function theorem gives one graph for that pair near the apex.
Thus in each half-plane there is at most one transition for
each unordered pair. If there is a positive transition take
the largest one, and if there is a negative transition take
the smallest one. At observations \((x,c|x|^q)\) with
\(|x|\to\infty\), the conditional sample, expressed relative
to \(p\), has the same horizontal sign eventually and its
invariant coordinate tends almost surely to \(r^{1-q}c\).
The bounded drift \((r-1)p\) and Gaussian displacement do not
affect this limit. The limit lies beyond the extreme transition,
so one of the two adjacent labels has conditional probability
tending to one. Its contrast with the other tends to one.
This is incompatible with the constant stationary contrast
on the arc, which lies strictly between \(-1\) and \(1\).

With no nonzero transition, the labels are constant on the four
coordinate quadrants. Since all three labels occur, one label
is used twice. If the quadrants with that label are adjacent,
their union is a half-plane, excluded by the full rank of each
individual cell. If the quadrants are opposite, all four
incident edges carry genuine jumps. At the apex, the adjusted
score gradient for the label in quadrant \((s,t)\) has the
form \((u_s,v_t)\), where \(s\) is the horizontal sign and
\(t\) the vertical sign. This follows by matching tangential
gradient components across each interface. The repetition
of one label on opposite quadrants then forces \(u_+=u_-\)
and \(v_+=v_-\), contradicting (G20) on a genuine edge.

For rank one, after an orthogonal change the label on each side
of \(x=p_x\) depends only on \(y\). If both profiles agree,
the whole partition is cylindrical. Otherwise some open
interval is a genuine vertical \(i,j\) interface, by the same
generic-slice argument. Set \(w=\Pr[Y_x>p_x\mid X_x=p_x]\)
and denote the half-plane labels at height \(y\) by \(I_+(y)\)
and \(I_-(y)\). The adjusted pair equality extends along the
whole line by analyticity. Tangential Gaussian deconvolution
then makes \(w(e_{I_+(y)})+(1-w)(e_{I_-(y)})\) have a constant
\(i,j\) contrast in \(y\). Two different active pair types
would give two independent constant contrasts. The resulting
constant probability vector contains at most two labels,
contrary to the three positive masses. Thus all genuine
jumps across the line use the same unordered pair \(i,j\).
If \(w\ne1/2\), the constant contrast from one such jump is
\(\pm(2w-1)\); it is neither zero nor \(\pm1\). Hence the
two profiles cannot agree at any height of positive measure,
and the third label is absent. If \(w=1/2\), the third label
must occur in both profiles or in neither, so its cell is
cylindrical. Both cases are excluded. Deconvolution here uses
the nonvanishing Fourier multiplier of Gaussian convolution,
by division on compact frequency tests.

One consequence for (CT7) is that the flux kernel on symmetric
matrices detects a nonzero constrained mixed covariance derivative
whenever it contains a nonzero semidefinite matrix. If it
contains none, there is a positive definite linear combination of
the centered second-moment matrices \(\int_{A_i}(xx^{\mathsf T}-I)d\gamma\).
Indeed project the compact convex trace-one positive semidefinite
set onto the orthogonal complement of the flux kernel. Its image
excludes zero, so strict separation supplies the indicated
positive definite matrix in that complement. This dichotomy
concerns anisotropic covariance derivatives.
[Anisotropic dilation identity and flux convention](THREE_CELL_CORRELATION_CONTACT.md#anisotropic-dilation-and-mass-tangent-directions)

## Restrictions on an indefinite linear symmetry

For a centered indefinite flow, normalize the generator to
\(\operatorname{diag}(1,-q)\), \(q>0\), and use the profiles
in \(u=y|x|^q\) on the two sign half-planes of \(x\).
If such a flow preserves a balanced planar maximizer, there
are no active axis jumps, \(q\) is rational, and on each sign
half-plane at least one of the two tails of its profile fails
to be eventually constant.

An active axis jump on \(x=0\) between \(i,j\) continues its
stationary score equality along the full line. Deconvolve in
\(y\) to obtain
\(\int (1_{A_i}-1_{A_j})(x,y)p_\sigma(x)\,dx=c_{ij}\).
On each sign of \(y\), flow invariance lets one write the
contrast as a function of \(x|y|^{1/q}\), so its even part
has constant Gaussian averages at every width, by continuity
of those averages in the width. Substituting the square of
the normal coordinate and applying Laplace uniqueness shows
that this even part is the constant \(c_{ij}\) almost
everywhere. The genuine jump makes it zero on a neighborhood
of the line inside the jump patch, hence \(c_{ij}=0\).
Thus reflection swaps \(i,j\) globally and fixes the third cell.
The argument for an active jump on the other axis is identical
after reordering axes and reversing time. The reflected graph
description (RF1) would then have \(h(\lambda x)=\lambda^{-q}h(x)\)
for \(\lambda>0\), using \(q\) for the ratio in those coordinates.
Continuity and finiteness at \(0\) would give \(h=0\),
incompatible with balance.

A nonzero transition \(c\) of one of the profiles exists
because constant coordinate-sign labels cannot be balanced
with three labels. On the corresponding arc, the entire
adjusted score contrast \(F_{ij}\) satisfies
\(F_{ij}(\pm e^z,c e^{-qz})=0\) for real \(z\), hence for
all complex \(z\). The entire extension in the two arguments
follows directly from Gaussian smoothing of bounded indicators.
If \(q\) is irrational, shifting real \(z\) by \(2\pi i n\)
gives infinitely many roots on a circle in the second argument
at each fixed horizontal coordinate of that sign. The identity
theorem would make the contrast vanish identically, impossible
by Gaussian injectivity.

Finally suppose one of the profiles has eventually constant
tails. Write \(f_i(u)\) for its indicator of label \(i\),
\(t_i^\pm\) for its two tail values, and
\(d_i=\int(f_i(u)-t_i^{\operatorname{sgn}u})du\).
At an observation \((x,c|x|^{-q})\) in its half-plane, as
\(|x|\to\infty\) the score has expansion
\[
 U_i=\frac{t_i^++t_i^-}{2}+
 \frac{\phi(0)}{\sigma |x|^q}
       [rc(t_i^+-t_i^-)+r^{-q}d_i]+o(|x|^{-q}).
\]
To verify it, condition on the horizontal source coordinate.
It has the prescribed sign and magnitude at least \(r|x|/2\)
apart from an exponentially small probability. The step
at \(u=0\) gives the stated normal-CDF expansion. The compactly
supported difference contributes by the change of vertical
variable from \(y\) to \(u=y|X_{\rm horiz}|^q\). After
multiplying by \(|x|^q\), dominated convergence gives
\(r^{-q}\phi(0)d_i/\sigma\).

All labels occurring in this profile must have equal leading
adjusted score, by taking \(c\) in the interior of each
winning interval in turn. Thus the next affine terms in \(c\)
must support its winners. If the tail labels coincide and
another label occurs, the common tail label has \(d_i<0\)
and another has \(d_j>0\), while both slopes are zero.
This is impossible. If the two tail labels differ, the
slopes of the used labels are ordered and distinct, so
there are two or three intervals. With two intervals and
switch at \(c\), the intercept integrals of the left and
right labels are \(c,-c\); the tie equation gives
\((r^{-q}-r)c=0\). But a switch at \(c=0\) is an active axis
jump. With three intervals and switches \(a<b\), the
intercept integrals are \(a,b-a,-b\); the equations on
the switches give \(b=(2-r^{q+1})a\) and
\(a=(2-r^{q+1})b\), inconsistent with \(a<b\).
The sole remaining profile with eventual constant tails
would give one label a whole half-plane, contradicting its
mass of \(1/3\). This leaves the possibility of unbounded
transitions in an indefinite invariant profile.

## No centered trace-zero symmetric flow

For \(q=1\) even the profiles with unbounded transitions are
impossible. In this case the invariant variable \(u\) in one
horizontal half-plane is the product of the signed horizontal
coordinate (made positive in that half-plane) and the vertical
coordinate. Consider a genuine transition \(u=c\) between
labels \(i,j\), and let \(F(u)\) be their bounded indicator
contrast on this half-plane, now defined on the whole real
\(u\)-axis. The actual contrast at observation
\((x,c/x)\) with \(x>0\), using the positive coordinate for
the chosen half-plane, is constant by stationarity. Write
\(Y_1,Y_2\) for independent normals of variance \(\sigma^2\)
and means \(rx,rc/x\). Its difference from
\(\mathbb E F(Y_1Y_2)\) is at most \(2\Pr[Y_1<0]\).
This is \(O(e^{-V x^2}/x)\) for \(V=r^2/(2\sigma^2)\).
For the other half-plane use the reflected horizontal
coordinate in the same calculation.

Completing the square after conditioning on \(Y_1\) gives
the characteristic function of the product at frequency \(\xi\):
\[
 \mathbb E e^{i\xi Y_1Y_2}
    =g(\xi)e^{i\alpha(\xi)c}e^{-v(\xi)s},\qquad
 s=x^2+c^2/x^2,
\]
where
\[
 g(\xi)=(1+\sigma^4\xi^2)^{-1/2},\qquad
 \alpha(\xi)=\frac{r^2\xi}{1+\sigma^4\xi^2},\qquad
 v(\xi)=\frac{r^2\sigma^2\xi^2}{2(1+\sigma^4\xi^2)}.
\]
Here \(v\) is even, and increases from \(0\) to \(V\) on the
positive frequency axis. Interpret \(\widehat F\) as a
distribution with Fourier convention \(e^{-i\xi u}\).
It pairs continuously with \(H^1(\mathbb R_\xi)\): by
Plancherel and Cauchy–Schwarz the inverse Fourier transform of
an \(H^1\) function has \(L^1\) norm bounded by its \(H^1\)
norm times a fixed constant. For each fixed \(c\) define the
compactly supported distribution
\[
 \langle\mu_c,\psi\rangle=
  (2\pi)^{-1}\langle\widehat F,\,g e^{i\alpha c}\psi(v)\rangle .
\]
The \(H^1\) norm of the test factor is controlled by the
\(C^1([0,V])\) norm of \(\psi\), since \(g,g'\) lie in \(L^2\)
and \(\alpha',v'\) are bounded. Thus the distribution has
order at most one. In particular the inverse characteristic
function is integrable, so the pairing with \(\psi(v)=e^{-sv}\)
equals \(\mathbb E F(Y_1Y_2)\).

Let \(\kappa\) be the stationary pair score on the transition.
The Laplace transform of \(\mu_c-\kappa\delta_0\) is
\(O(e^{-Vs})\) as \(s\to+\infty\), since \(s-x^2=O(x^{-2})\)
and every sufficiently large \(s\) is attained. Push this
distribution forward by \(z=e^{-v}\) onto \([e^{-V},1]\).
Its moments of integer order \(n\) are \(O(e^{-Vn})\).
Consequently its Cauchy transform, initially expanded at
infinity using these moments, extends analytically to the
region of modulus larger than \(e^{-V}\). This forces the
pushforward distribution to vanish on \((e^{-V},1]\).
To see it, subtract the two Cauchy values just below and
above the real axis. Their difference is the convolution
of the distribution with the Poisson kernel times \(2\pi i\).
Test against a smooth function compactly supported in that
region of the real axis. Holomorphic continuation makes
the difference tend to zero, while the approximate identity
converges against the compact distribution in \(C^1\).
Thus \(\mu_c\) vanishes on the open interval \((0,V)\).

Since \(v\) is a smooth bijection from \((0,\infty)\) onto
\((0,V)\), the last conclusion gives the distributional equation
on \(\xi>0\)
\[
 e^{i\alpha(\xi)c}\widehat F(\xi)
    +e^{-i\alpha(\xi)c}\widehat F(-\xi)=0.
\]
Here \(g>0\) has been canceled on compact frequency intervals.
For two distinct genuine transitions \(c_1,c_2\) of the same
pair in the half-plane, the two equations force the Fourier
support on positive frequencies into the zeros of
\(\sin(\alpha(\xi)(c_1-c_2))\). There are finitely many such
positive frequencies because \(\alpha\) is positive, has one
maximum and tends to zero at both ends of \((0,\infty)\).
The same holds on negative frequencies. Hence \(\widehat F\)
has finite support, including at most an additional point at
zero. A distribution supported on a finite set is a finite
linear combination of derivatives of point masses. Its inverse
Fourier transform is a continuous exponential polynomial,
whereas \(F\) is finite-valued almost everywhere. Continuity
then makes \(F\) constant, contradicting the jumps.

It follows that each pair has at most one genuine transition
on each sign half-plane. The resulting profiles have eventually
constant tails, already excluded. Thus no centered nonzero
trace-zero symmetric flow preserves a balanced three-cell
maximizer.

In particular, a nonsimplex maximizer at which all zero-flux
symmetric covariance directions fail to detect a constrained
mixed derivative would have a one-dimensional flux kernel.
If that kernel had dimension at least two, it would contain
a nonzero trace-zero matrix, which has nonzero normal speed
by the result just proved. A generator of the one-dimensional
kernel must be indefinite of nonzero trace by the two flow
exclusions and (CT7). The indefinite profile restrictions
then apply to that generator.

## A quantitative restriction on centered straight walls

Suppose \(E\) is the tangential coordinate set of the \(i,j\)
interface on a line through the origin, with a fixed choice
of normal pointing into \(j\). Use coordinates \((s,z)\) with
this normal coordinate \(z\), and let \(k\) denote the remaining
label here. Let \((T_0,S_0)\) be a standard scalar Gaussian pair
of correlation \(r\) and \(Z\) be independent with variance
\(\sigma^2=1-r^2\). Put
\(B=\{(s,z)\in A_i:z>0\}\cup\{(s,z)\in A_j:z<0\}\).
Then
\[
 \mathbb E[|Z|(1_{A_k}+2\,1_B)(S_0,Z)1_{\{T_0\in E\}}]
   \le \sigma\sqrt{2/\pi}\Pr[T_0\in E,\ S_0\notin E].        \tag{SW1}
\]
In particular \(\gamma_1(E)\le1-\delta\) with \(\delta>0\)
uniform over a compact correlation band.

At \((t,0)\) the inward \(j\) derivative of
\(T_r(1_{A_j}-1_{A_i})\) equals
\[
 2r p_\sigma(0)-
 \frac r{\sigma^2}
   \mathbb E[|Z|(1_{A_k}+2\,1_B)(S_0,Z)\mid T_0=t].
\]
This follows by differentiating the conditional normal mean
and comparing the integrand to the step \(\operatorname{sgn} Z\).
Test (G1) with the indicator speed of this portion of the pair
wall, using the weighted-space extension of (G17).
After division by the line-density factor \(\phi(0)\), the
interaction term is \(2r p_\sigma(0)\Pr[T_0,S_0\in E]\):
the incidence norm is two and the normal Mehler factor at the
two zeros is \(1/\sigma\). Substituting the derivative and
rearranging gives (SW1). The marginal expectation
\(\mathbb E[|Z|1_{A_k}(S_0,Z)]\) has a positive lower bound
depending only on the band. Indeed choose a fixed interval
\(a\le|Z|\le b\) with \(a>0\) so that its standard-normal
complement has mass less than \(1/6\). It intersects \(A_k\)
in standard planar Gaussian mass at least \(1/6\), and the
normal density ratios there are uniformly positive. The part
of this expectation with \(T_0\notin E\) is at most
\(\mathbb E|Z|\Pr[T_0\notin E]\) by independence. Inequality
(SW1) then gives the positive lower bound on \(1-\gamma_1(E)\).

If a centered line contains a straight \(i,j\) arc, the
adjusted pair equality extends to the entire line. If the
third adjusted score is never strictly above these scores on
the line, and it is not identically tied there, then off
isolated third-score ties and isolated zeros of the normal
gradient the line is an active \(i,j\) wall. The normal
derivative is nonzero on the given arc by (G20), so its
zeros are discrete by analyticity. Any such zero would be
approached by an actual wall segment and contradict the
reciprocal-gradient condition (G20). Its sign is thus fixed,
giving \(\gamma_1(E)=1\) and contradicting (SW1).

The three scores also cannot tie identically on this line.
In that case consider, for each pair, the first normal Taylor
coefficient of its adjusted contrast that is not identically
zero along the line. Their zero sets are discrete. Away
from them the strict winners immediately on either side are
constant on each complementary interval. On each such
interval at least one pair has different first normal
derivatives, since the given regular arc has a nonzero pair
gradient. Consequently the side winners are distinct.

If a side winner changes between two labels at an exceptional
point, their first normal derivatives coincide there. The
tangential derivatives coincide by the line equality. If
all three normal derivatives coincide, any incident wall
arc accumulating there contradicts (G20). Otherwise the
third derivative is strictly below the common maximum on
the positive side, or strictly above the common minimum
on the negative side. Continuity of the first derivatives
and the exact line equality exclude the third label on a
small collar on the given side. Both changing labels occur
in arbitrarily small such collars, so by the analytic
interface stratification a jump arc for this pair accumulates
at the exceptional point. Its vanishing gradient again
contradicts (G20). Thus the side winners never change,
giving a whole consistently oriented wall in contradiction
to (SW1). Every centered line containing a straight pair
arc therefore has an open segment dominated by the third
adjusted score.

## A limitation of the signed interface spectral identities

The translation and rotation identities and the balanced
spectral gap can all hold for an abstract planar Mehler
matrix with two different normals for each pair type. A
proof of flatness requires further geometric information.
The following six-point example gives this limitation for
every \(0<r<1\).

Take two triples of points \(p_{sj}\) of radius \(R\), where
\(j=0,1,2\) has angle \(2\pi j/3\) in the first triple and
\(\theta+2\pi j/3\) in the second, with \(0<\theta<\pi/3\).
Give point \(p_{sj}\) cyclic pair type \(j\), with incidence
vector \(d_j\) of squared norm two and mutual dot products
\(-1\). Use counting measure and multiply the relative Mehler
kernel by \(d_j\cdot d_l\). The common positive Mehler factor
at these radii can be dropped from the matrix and its
gradient weights. With \(t=rR^2/(1-r^2)\) and
\(\omega=e^{2\pi i/3}\), put
\[
 A=e^{t\cos\theta},\quad B_\pm=e^{t\cos(\theta\pm2\pi/3)},\quad
 B=2A-B_+\omega-B_-\omega^2,\quad C=2A-B_+-B_-,
\]
\[
 D=2e^t+e^{-t/2},\qquad E=2e^t-2e^{-t/2}.
\]
The cyclic Fourier block at mode zero has diagonal \(E\)
and off-diagonal \(C\); at mode one the entries are \(D\)
on the diagonal and \(B,\bar B\) off the diagonal. The
remaining block is the complex conjugate. Choose \(t>0\)
with \((E+C)/(D+|B|)=r\), which exists by continuity and
the limits zero at \(t=0\) and one at infinity.
The common weight is \(a=E+C=r(D+|B|)>0\).

Identify the plane with \(\mathbb C\) and prescribe unit normals
\[
 N_{1j}=e^{i\chi}\omega^j,\qquad
 N_{2j}=e^{i(\chi-\arg B)}\omega^j,\qquad
 \chi=\pi/2+(\arg B+\theta)/2 .
\]
They are the two real coordinate copies of a mode-one
eigenvector of the kernel matrix \(K\), so \(rKN=aN\).
The two normals of a given type differ because \(B_+\ne B_-\)
and \(\operatorname{Re}B>0\). For a right-angle rotation \(J\),
the six speeds \(Jp_{sj}\cdot N_{sj}\) are the same nonzero
constant \(R\cos((\arg B+\theta)/2)\). They form the
mode-zero common eigenvector of eigenvalue \(a\).

Every other eigenvalue is strictly below \(a\), except for
the two translation eigenvalues \(a/r\): \(C>0\) and
\(D-|B|\le D-\operatorname{Re}B<E+C\), since
\(D-E=3e^{-t/2}<C+\operatorname{Re}B\).
Balance is \(\sum_{s,j}d_j f_{sj}=0\). It includes both
mode-zero vectors, and on the nonconstant modes it imposes
opposite coefficients in the two triples. Their restricted
interaction Rayleigh quotient is \(D-\operatorname{Re}B<a\).
Thus the balanced form \(aI-K\) has only the rotation
zero mode and a gap on its balanced complement. Its
unrestricted negative index is two, and \(aI-rK\) has
exactly the two translation zero modes. The example consists
of point weights and assigned normals, so the missing input
is their realization as interfaces of cells.
