# Joint rotations of Gaussian angular slices

At fixed angular templates on the circles about the Gaussian
origin, rotation on each circle preserves every cell mass.
The rotations enter the noise-stability objective through all
angular harmonics. A deterministic positive semidefinite form
can certify a joint optimum. The certificate uses the
stationary torque of the full objective and controls the
difference between each higher harmonic and its first-harmonic
bound.

## A certificate for all harmonics

Let \(\nu\) be a probability measure on the radii under
consideration. For each \(m\ge1\), let \(C_m(s,t)\) be a
measurable Hermitian kernel, with \(C_m(t,s)=\overline{C_m(s,t)}\),
and suppose
\[
                  \iint\sum_{m\ge1}m^2|C_m(s,t)|
                             \,d\nu(s)d\nu(t)<\infty .
\]
For measurable unit phases \(z\), consider the real energy
\[
 \mathcal E(z)=\sum_{m\ge1}\iint
       \overline{z(s)^m}C_m(s,t)z(t)^m\,d\nu(s)d\nu(t).       \tag{AS1}
\]
Assume \(x\) is stationary for all bounded real angular
variations \(x(s)e^{i\epsilon h(s)}\). Gauge \(x\) to the
constant phase by putting
\[
 c_m(s,t)=\overline{x(s)^m}C_m(s,t)x(t)^m=a_m(s,t)+i b_m(s,t),
 \qquad \beta(s,t)=\sum_{m\ge1}m b_m(s,t).
\]
The derivative at \(x\) in the direction \(h\) is
\(2\int h(s)\big[\int\beta(s,t)d\nu(t)\big]d\nu(s)\).
The differentiation is justified by the integrability
hypothesis. Hence the skew real kernel \(\beta\) has zero
row integrals almost everywhere.

For a competing phase \(z\), put \(y=\overline{x}z\) and
\(u=\overline{y(s)}y(t)\) on each ordered pair. Taking real
parts in (AS1) gives the exact pair-deficit identity
\[
 \mathcal E(x)-\mathcal E(z)=
  \iint\sum_{m\ge1}
    \bigl[a_m(1-\operatorname{Re}u^m)+b_m\operatorname{Im}u^m\bigr]
                                  \,d\nu(s)d\nu(t).         \tag{AS2}
\]
With \(a_m^-=\max(-a_m,0)\), define the real symmetric weight
\[
 \omega=a_1-\sum_{m\ge2}
       \bigl[m^2 a_m^-+m(m-1)|b_m|\bigr].
\]
Write \(K_\omega,K_\beta\) for the integral kernels and
\(D_\omega\) for multiplication by the row integral of
\(\omega\). Then the Hermitian form
\(\mathcal L=D_\omega-K_\omega-iK_\beta\) satisfies
\[
              \mathcal E(x)-\mathcal E(z)
                             \ \ge\ \langle y,\mathcal L y\rangle .   \tag{AS3}
\]
The form is well defined for bounded test functions by the
integrability hypothesis. If the absolute row integrals of
\(\sum_m m^2|C_m|\) are essentially bounded, the Schur bound
makes these operators bounded on \(L^2(\nu)\).

To prove (AS3), use
\[
 1-\operatorname{Re}u^m\le m^2(1-\operatorname{Re}u),\qquad
 \big|\operatorname{Im}u^m-m\operatorname{Im}u\big|
                   \le m(m-1)(1-\operatorname{Re}u).
\]
The first inequality follows by factoring \(1-u^m\) and
squaring its modulus. For the second,
\[
 u^m-1-m(u-1)=(u-1)^2
                      \sum_{\ell=0}^{m-2}(m-1-\ell)u^\ell .
\]
Since \(|u-1|^2=2(1-\operatorname{Re}u)\), the result follows.
On each pair the sum in (AS2) is therefore at least
\(\omega(1-\operatorname{Re}u)+\beta\operatorname{Im}u\).
Its integral is the right side of (AS3).

If \(\mathcal L\) is positive semidefinite, \(x\) is a
global maximizer among the allowed rotations. Stationarity
makes \(\mathcal L1=0\). A form bound by \(\lambda>0\) on
the orthogonal complement of the constants gives the deficit
bound \(\lambda(1-|\int\overline{x}z\,d\nu|^2)\) and uniqueness
up to a common rotation. With a single harmonic and finitely
many uniformly weighted radii, (AS3) is an equality and
\(\mathcal L\), up to scalar normalization and diagonal
unitary conjugation, is
\(\operatorname{Diag}(\operatorname{Re}((C_1x)\overline x))-C_1\),
the rank-one dual certificate in Bandeira, Boumal and
Singer's Lemmas 4.3 and 4.4. Their deterministic result
identifies this matrix's positivity with rank-one optimality
in the lifted semidefinite program. The higher-harmonic
form above is a sufficient certificate for the joint rotation
problem.
[Bandeira, Boumal and Singer, published paper, Lemmas 4.3–4.4](https://web.math.princeton.edu/~amits/publications/Tightness_MPSA_Published.pdf)

## The Gaussian angular kernels

For label \(j\), let \(\widehat f_{j,s}(m)\) be the
\(m\)-th complex Fourier coefficient of its indicator
at radius \(s\) before rotation, using \(e^{-im\theta}\)
and normalized circular measure. Denote the joint radial
density relative to \(\nu\otimes\nu\) by \(w_\rho(s,t)\),
restricting to pairs in the radial set being optimized.
Put \(\kappa=\rho st/(1-\rho^2)\). The angular Fourier
eigenvalues of the conditional kernel are \(I_m(\kappa)/I_0(\kappa)\).
For the energy above, take
\[
 C_m(s,t)=2w_\rho(s,t)\frac{I_m(\kappa)}{I_0(\kappa)}
                \sum_j\widehat f_{j,s}(m)
                            \overline{\widehat f_{j,t}(m)},
                         \qquad m\ge 1.                    \tag{AS4}
\]
A rotation by angle \(\theta(s)\) multiplies each Fourier
coefficient by \(e^{-im\theta(s)}\), so (AS1) follows with
\(z=e^{i\theta}\) after dropping the zero harmonic, which
is independent of rotations. The coefficients are factorially
small in \(m\) uniformly on any compact annulus because
\(I_m(\kappa)/I_0(\kappa)\le(\kappa/2)^m/m!\), as follows
from the power series for \(I_m\). Thus the bounded row
condition holds there.

For the full positive radial axis take \(d\nu(s)=s e^{-s^2/2}ds\).
The joint density \(w_\rho(s,t)\,d\nu(s)d\nu(t)\) is
\[
 \frac{st}{1-\rho^2}
 e^{-(s^2+t^2)/(2(1-\rho^2))} I_0(\kappa)\,ds\,dt .
\]
The second derivative at zero of \(e^{\kappa\cos\theta}\)
gives \(\sum_{m\ge1}m^2 I_m(\kappa)=\kappa e^\kappa/2\).
The Fourier coefficients of indicators are bounded by one,
so the assumed integrated \(m^2\) bound follows from
integrability of a polynomial times
\(e^{-(s^2+t^2)/(2(1+\rho))}\).

## Stability of a common reference orientation

Suppose \(C_m^0\) are real symmetric reference kernels on a
probability space, with \(\sum_m m^2\|C_m^0\|_\infty<\infty\).
Assume a uniform pair margin
\[
 \sum_{m\ge1} C_m^0(s,t)(1-\cos(mu))\ge c(1-\cos u)
                      \qquad(c>0) .                        \tag{AS5}
\]
If \(\sum_m m^2\|C_m-C_m^0\|_\infty\) is sufficiently small,
the joint phase problem has a unique global maximizer up
to common rotation. Its mean-zero real phase has
\(L^\infty\) norm at most a constant times the size of the
perturbation. Uniqueness is among all measurable phases.

To find the stationary phase, consider the real mean-zero
map
\[
 T_C(\theta)(s)=\int \sum_m m\operatorname{Im}
       \bigl(C_m(s,t)e^{im(\theta(t)-\theta(s))}\bigr)\,d\nu(t).
\]
It is continuously differentiable on \(L^\infty\), with the
kernel parameter in the displayed weighted summable norm.
The symmetry gives \(\int T_C=0\). The linearization at
\((C^0,0)\) is \(-L_W\), where
\[
 W=\sum_m m^2 C_m^0,\qquad
 L_W f(s)=\int W(s,t)(f(s)-f(t))\,d\nu(t).
\]
Taking \(u\to0\) in (AS5) gives \(W\ge c\). The real form
on mean-zero functions satisfies
\(\langle f,L_W f\rangle
=\tfrac12\iint W(s,t)(f(s)-f(t))^2\ge c\|f\|_2^2\).
The inverse on \(L^2\) mean-zero maps \(L^\infty\) to
\(L^\infty\) boundedly. Indeed, writing
\(d(s)=\int W(s,t)d\nu(t)\ge c\), a solution of \(L_W f=g\)
satisfies \(d f=g+K_W f\) and
\(\|K_W f\|_\infty\le\|W\|_\infty\|f\|_2\).
The implicit function theorem on mean-zero \(L^\infty\)
therefore gives a small stationary \(\theta_*\).

Now gauge by \(x=e^{i\theta_*}\), and put
\(\eta=\sum_m m^2\|c_m-C_m^0\|_\infty\).
Then \(\eta=o(1)\) with the perturbation. This follows for
each finite set of harmonics by the small phase, and for
the tail by the summable reference bound. Apply
\(u^m-1=m(u-1)+[u^m-1-m(u-1)]\) to each complex coefficient
error in (AS2). After retaining \(\beta\operatorname{Im}u\),
the real linear term and the complex remainder cost at
most \(\eta(1-\operatorname{Re}u)\) together, since
\(m+m(m-1)=m^2\). The integrated \(\beta\) term has absolute
value at most \(\|K_\beta\|\,\|y-\int y\|_2^2\) because
\(K_\beta 1=0\). The Schur bound gives \(\|K_\beta\|\le\eta\).
Using (AS5) proves
\[
 \mathcal E(x)-\mathcal E(z)
       \ge (c-2\eta)\left(1-\left|\int\overline{x}z\,d\nu\right|^2\right).
                                                             \tag{AS6}
\]
The comparison uses a strict pair margin at the reference
and is independent of signs of individual higher harmonics.

On a compact annulus, (AS5) holds for a common template in
which two incident arcs have equal lengths \(a(s)\), with
\(a(s)\in[\alpha,\pi-\alpha]\) for some \(\alpha>0\).
These lengths can vary measurably across radii. The
centered overlap of the three labels is even, nonincreasing
in circular distance and nonconstant. Decompose it into
centered superlevel arcs and convolve with the conditional
von Mises kernel, as in (RU1). The result has strict maximum
at zero and strictly negative second derivative there.
In fact each nontrivial superlevel arc contributes
\(k'_\kappa(\ell)-k'_\kappa(-\ell)<0\) at zero, where
\(0<\ell<\pi\) is the arc half-length. Compactness of
the two widths, the angular parameter and the angle
gives (AS5), using the curvature for \(u\) near zero.
The positive radial density on the compact annulus is
bounded below. Perturbing the two incident widths
independently by a small uniform amount changes the kernels
by a small amount in the weighted norm, by the Fourier
coefficient bound and factorial decay in (AS4). Thus the
adapted phase remains the unique global optimum among
rotations even when the two incident widths differ.
[The equal-adjacent-arc comparison (RU1)](THREE_CELL_RADIAL_UNTWISTING.md#alignment-on-a-pair-of-circles)

## Exact optimization over slice shapes near equal masses

On a fixed compact annulus, prescribe measurable circular mass
vectors \(p(s)\) uniformly close to \((1/3,1/3,1/3)\).
For the objective that includes pairs of radii within the
annulus, every maximizing family consists of consecutive
arcs in one common cyclic order. The angles are the unique
phases of (AS6), up to an overall rotation after choosing
the cyclic order. Reversing the common cyclic order gives
the reflected maximizer with the same value. In particular,
the near-balanced joint optimization has this exact solution
over arbitrary measurable slices.

To see the shape assertion, first allow fractional assignments
with the prescribed circular masses almost everywhere. They
form a weak-star compact set. The objective is continuous
on this set because the compact-annulus kernel is uniformly
approximable by finite sums of products. Hence there is a
maximizer. If the \(L^\infty\) mass error tends to zero,
any weak-star limiting maximizer maximizes the equal-mass
problem. Testing consecutive arcs with a common orientation
gives the lower bound on its value. The equal-mass case of
(JA2) bounds the value pair by pair, and its equality case
for fractional slices forces the same ordinary equal arcs
at almost every pair of radii. Hence the limit is a common
equal-arc assignment, and weak-star convergence to indicators
upgrades to \(L^1\). This argument is uniform on compact
correlation bands by also taking a convergent subsequence
of the correlations.

At a maximizing assignment the scores for a slice are the
angular convolutions against all radii in the annulus,
including the joint radial weights. The scores of every
sufficiently near-balanced optimizer are uniformly close,
in the angular \(C^1\) norm over the annulus, to those of
its limiting equal arcs, after a common isometry. Indeed the
kernel and its angular derivative are bounded and the
preceding compactness gives \(L^1\) closeness. The reference
scores have three simple winning crossings, one at each
sector boundary, with the third score strictly below.
The derivatives of the crossing contrasts are negative in
the direction from the ending to the starting label, and
uniformly separated from zero. These facts follow from the
strictly decreasing equal-arc convolution and the strictly
positive radial weights on the compact annulus.

Here is a linear best-response argument that also handles
fractional assignments. Write the endpoints with their given
widths as \(a_i+x_i(s)+u\), where \(a_i=2\pi i/3\),
\(\sum_i x_i(s)=0\) and
\(p_i(s)=1/3+(x_i(s)-x_{i-1}(s))/(2\pi)\).
For each radius, add the three incident score contrasts
evaluated at these tentative endpoints. This sum has a
unique zero near \(u=0\), because its derivative there is
uniformly negative and its value at zero is small. The
root is measurable, for instance by bisection. The three
contrasts at this root sum to zero, so they are consecutive
differences of three small label prices. Subtract the prices
from the scores. The priced winner changes exactly at the
three chosen endpoints, by the sign and simple-crossing
bounds. Thus the resulting arcs, which have the required
masses, are the unique linear best response. Integrating
the response supplies a measurable competitor. At the
quadratic maximizer the derivative toward this competitor is
nonpositive, whereas strict priced winner inequalities make
it positive whenever the optimizer differs on a set of
positive measure. The competitor therefore equals the
optimizer. This proves the shape assertion and allows (AS6)
to give the global optimum for the prescribed masses.
[Equal-mass equality for fractional slices and (JA2)](THREE_CELL_JOINT_ANGULAR_ENVELOPE.md#exact-pair-envelope-near-balanced-masses)

## The quartic loss from optimizing pairs separately

Let \(x(s)\) be the mean-zero endpoint increments from the
equal-arc template, \(\epsilon=\|x\|_\infty\), and use a
probability radial measure \(\nu\) on the compact annulus.
Let \(w_\rho\) be the corresponding radial pair density.
Write \(F_\kappa(x,y,u)\) for the pair interaction when the
first arc template is rotated by the angle \(u\) relative to
the second. The quadratic torque from (RC1) is
\[
 \begin{split}
 W(s,t)&= w_\rho(s,t)\frac{6\Delta_\kappa}{(2\pi)^2},\\
 G(s,t)&=\frac{w_\rho(s,t)}{(2\pi)^2}\,
                      x(s)^{\mathsf T} B_\kappa x(t),\\
 g(s)&=\int G(s,t)\,d\nu(t),\qquad
                    \vartheta=L_W^{-1}g .                  \tag{AS7}
 \end{split}
\]
Here \(B_\kappa\) is the skew matrix in (RC1) and
\(\Delta_\kappa=k_\kappa(0)-k_\kappa(2\pi/3)\). The inverse
uses the mean-zero gauge, since \(G\) is antisymmetric.
The optimizing phase in this gauge is
\(\theta_*=\vartheta+O(\epsilon^3)\) in \(L^\infty\).
Indeed the derivative of \(w_\rho F_\kappa\) at \(u=0\)
is \(G+O(\epsilon^3)\), its second derivative is
\(-W+O(\epsilon)\), and the bounded inverse argument for
\(L_W\) first gives \(\theta_*=O(\epsilon^2)\). Applying the
same inverse to the stationarity equation gives the expansion.

With ordered radial pair integration, the gain over the
common-angle family is
\[
      \mathcal E(e^{i\theta_*})-\mathcal E(1)
          =\langle g,L_W^{-1}g\rangle+O(\epsilon^5).         \tag{AS8}
\]
The separate pair envelope exceeds the exact joint
optimum by
\[
 \frac12\iint W(s,t)
  \left[\frac{G(s,t)}{W(s,t)}
               -\bigl(\vartheta(s)-\vartheta(t)\bigr)\right]^2
                   d\nu(s)d\nu(t)+O(\epsilon^5).            \tag{AS9}
\]
For a phase of size \(O(\epsilon^2)\), Taylor expansion gives
the phase-dependent term
\(2\langle g,\theta\rangle-\langle\theta,L_W\theta\rangle
+O(\epsilon^5)\). This proves (AS8). At an individually
optimal relative angle the pair gain is
\(G^2/(2W)+O(\epsilon^5)\), by (RC1) and the same Taylor
estimate. Subtract (AS8), or complete the square, to obtain
(AS9). The bracket is exactly the weighted orthogonal
residual from projecting the preferred pair shifts onto
differences of one common phase function. This gives the
leading loss for the actual pair envelope and the actual
joint optimum in the near-balanced regime.
[Preferred pair rotation (RC1)](THREE_CELL_ANGULAR_ROTATION_CYCLES.md#preferred-relative-rotation)

## Balanced widths and the determinant of their covariance

Write \(B_\kappa=b_\kappa M\) as in (RC1), and put
\[
 D(s,t)=x(s)^{\mathsf T}Mx(t),\qquad
 \alpha(s,t)=\frac{b_\kappa}{6\Delta_\kappa}>0 .
\]
Then \(G/W=\alpha D\) in (AS7). If
\(\int x\,d\nu=0\), the skew kernel \(D\) has zero row
integrals. Hence \(\iint D(s,t)(v(s)-v(t))\,d\nu(s)d\nu(t)=0\)
for every bounded \(v\). Weighted Cauchy–Schwarz gives
\[
 \iint W\big[\alpha D-(v(s)-v(t))\big]^2
 \ \ge\ \frac{(\iint\alpha D^2)^2}{\iint D^2/W}
 \ \ge\ W_{\min}\alpha_{\min}^2\,\iint D^2,                \tag{AS10}
\]
where the integrals use \(\nu\otimes\nu\) and the minima
are strictly positive bounds on the compact annulus. The
inequality is immediate when \(D=0\). In orthonormal
coordinates on the mean-zero label plane, \(M\) is a skew
rotation matrix with singular values \(\sqrt 3\). For the
centered \(x\) this gives
\(\iint D^2=6\det(\operatorname{Cov}_{\nu}x)\).
Thus the quartic coefficient in (AS9) is bounded below
by three times \(W_{\min}\alpha_{\min}^2\) times the
determinant, with \(x\) evaluated at unit perturbation
scale. For a centered profile of rank one, \(D=0\), so
the quadratic torque and the quartic loss both vanish.

A positive quartic bound also holds for each bounded
rank-two direction on the full radial axis. Let \(\mu\)
be ordinary Gaussian radial probability, let \(X(s)\) be
a bounded profile in the mean-zero label plane with
\(\int X\,d\mu=0\) and positive definite covariance on
that plane, and take \(x(s)=\epsilon X(s)\). Choose a
compact annulus \(H\) so large that, with
\(\nu=\mu|_H/\mu(H)\),
\[
 A_H=\int XX^{\mathsf T}\,d\nu,\qquad
 c_H=A_H^{-1}\int X\,d\nu,\qquad
             \|c_H\cdot X\|_{L^\infty(H)}\le\tfrac12 .
\]
The inverse acts on the label plane. Such an \(H\) exists
because \(A_H\) tends to the positive definite global
covariance and the conditional mean tends to zero.
The weight \(a(s)=1-c_H\cdot X(s)\) on \(H\) has
\(\int aX\,d\nu=0\) and \(1/2\le a\le3/2\).

Now use \(a(s)a(t)D(s,t)\) in (AS10) in place of \(D\)
as the test against the bracket. This kernel again has
zero row integrals under \(\nu\). The lower bound becomes
\[
 \iint_H W\big[\alpha D-(v(s)-v(t))\big]^2
       \ \ge\ \frac{W_{\min}\alpha_{\min}^2}{81}
                                       \iint_H D^2,
\]
where \(\iint_H\) integrates over \(H^2\). Indeed the
numerator has the lower bound
\(\alpha_{\min}^2(\iint_H D^2)^2/16\), and the
denominator is at most \(81\iint_H D^2/(16W_{\min})\).
Here \(\iint_H D^2=6\epsilon^4\det A_H>0\). Apply
(AS9) to the exact localized optimum over all measurable
slice shapes with the given masses on \(H\). For sufficiently
small \(\epsilon>0\), its pair-envelope deficit is at least
\(c_{\rho,X}\epsilon^4\). Every full-plane realization of
the mass profile has at least this deficit on \(H^2\),
with the radial restriction factors included in
\(w_\rho\), and the deficits for all other pairs are
nonnegative. Hence the full integrated pair envelope
exceeds the noise stability of every realization by
\(c_{\rho,X}\epsilon^4\), after absorbing the restriction
factors into the constant.

## Optimizing an annulus with the exterior fixed

Fix the annulus \(X=[d,D]\subset(0,\infty)\). Consider a
common-angle reference family as in the equal-adjacent
case of (AS5), with widths \(a^0(s)\in[\alpha,\pi-\alpha]\)
on \(X\), and arbitrary widths \(a^0(t)\in[0,\pi]\) outside.
Assume \(a^0\) is bounded below by a positive constant
on a compact set of positive radial measure in \(X^c\).
Choose the reference phase zero. Perturb the templates
uniformly in angular \(L^1\) on \(X\) and in integrated
angular \(L^1\) on \(X^c\). The exterior is held fixed at
this perturbed configuration. When the sum of these two
errors is small, rotations of the annular templates admit
a unique global optimum, pinned by the fixed exterior.
Its real phase is small in \(L^\infty(X)\).

Use the radial probability
\(\nu=\mu|_X+\mu(X^c)\delta_*\), with the exterior collapsed
to a point of phase one. For \(s\in X\), set
\[
 C_m(s,*)=\frac{1}{\mu(X^c)}\int_{X^c} C_m(s,t)d\mu(t)
\]
using the Gaussian coefficients from (AS4) before collapse.
Extend by the Hermitian relation. This retains both orders
of the annulus–exterior interaction exactly. The
exterior–exterior energy is constant, so for purposes of
(AS5) one may set \(C_1(*,*)=c\) and all its higher
coefficients equal to zero. This adds the same constant
for every phase. Pairs from \(X^2\) satisfy (AS5) at the
reference. Pairs with the exterior point do as well, by
retaining the positive-measure compact exterior set and
using the nonnegative loss on the remaining exterior
circles from (RU1).

The smallness of the kernel perturbation in (AS5)
follows by bounding each Fourier coefficient error by
the angular \(L^1\) error. For \(s\in X\) and \(t\in X^c\),
the radial factor in (AS4), omitting the factor two,
obeys
\[
 w_\rho(s,t)\frac{I_m(\kappa)}{I_0(\kappa)}
 \le \frac{e^{-\rho^2t^2/(2(1-\rho^2))}}{1-\rho^2}
                    I_m\!\left(\frac{\rho D t}{1-\rho^2}\right).
\]
The sum of the right side with weights \(m^2\) is uniformly
bounded in \(t\), by the Bessel sum in (AS4). Consequently
the two angular \(L^1\) errors control the weighted
summable kernel norm on \(X\cup\{*\}\). Apply (AS5)–(AS6)
and subtract the phase at the exterior point. Among
competitors that fix this phase, the variance in (AS6)
dominates \(\mu(X^c)\int_X|z(s)-x(s)|^2\,d\mu(s)\).
This pins the phase uniquely.

There is an exact shape conclusion near the balanced
reference, even when the exterior is perturbed. Prescribe
near-balanced widths on \(X\), and hold fixed any exterior
assignment whose integrated angular \(L^1\) distance from
the centered equal-arc simplex is sufficiently small.
Then the maximizer of the full stability over assignments
on \(X\) at the prescribed widths is the consecutive-arc
family in the reference cyclic order, with the pinned
phase just constructed. The shapes and phase are unique
almost everywhere. To prove this, optimize over fractional
assignments on \(X\) as above. If the width and exterior
errors tend to zero, the exterior scores converge uniformly
in angular \(C^1\) on \(X\). This follows from the Gaussian
kernel and the same estimate with the first angular
derivative. Every limiting optimizer thus maximizes the
balanced problem with exterior exactly the reference. The
equal-mass pair inequality bounds every pair meeting \(X\),
and equality on the mixed pairs pins every limiting slice
to the reference. Weak-star convergence to this ordinary
assignment gives \(L^1\) convergence inside \(X\). Now apply
the priced linear-response proof on \(X\), including the
fixed exterior in the scores. The optimizer has consecutive
arcs in the reference order with the prescribed masses,
and their joint phase is the unique solution above.

In this pinned setting an actual maximizer with the three
exact reference pair-score lines and equal prices is the
reference simplex on \(X\). At each radius the three
ordinary score contrasts are uniformly close in angular
\(C^1\) to the reference ones. In a neighborhood of
each reference endpoint, only the incident pair can win
and its contrast is strictly decreasing in the reference
orientation. Exact equality on the reference ray fixes
the unique crossing at the reference endpoint. The reference
winner persists away from these neighborhoods. This
corollary applies to compact discrepancies satisfying
the stated annular width and fixed-exterior smallness
hypotheses. It leaves open how to enforce those hypotheses
for an arbitrary maximizer.
[Full reference pair-score lines for a compact discrepancy](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#exact-germ-and-uniform-radius)

## Circular connectivity propagates the simplex germ

Let \(A\) be a balanced planar Gaussian noise-stability
maximizer with equal prices and all three exact centered
reference score lines. A circular partition into arcs means
that each label is a single connected arc up to angular null
sets. Empty arcs and a full circle are allowed. The three
labels together partition the circle, without a prescribed
cyclic order. Suppose almost every circle of radius \(s<R\)
has this property, for some \(R>0\). Then \(A\) agrees with
its centered simplex or antipodal germ throughout \(B_R\).

Use strict score winners to represent \(A\), resolving
score ties by any deterministic rule. The full reference
lines and the nonzero gradient differences at the origin
give the exact germ by the pair factorization (CG1).
For \(s>0\), call a radius exceptional if one of the pair
contrasts vanishes identically on its entire circle. These
radii form a locally finite set. Indeed, for each pair
one can choose a fixed radial direction on which its
contrast is nonzero in the germ. Its restriction to that
ray is analytic, so its zeros, including the exceptional
radii for this pair, are isolated. There are also no
exceptional radii in a sufficiently small disk by the
same argument near the origin.
[Exact centered germ, CG1](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#exact-germ-and-uniform-radius)

At any nonexceptional radius, every pair contrast is a
nonzero analytic function of the angle. Thus only finitely
many angles tie. The winner partition is continuous there
as a function of radius in angular \(L^1\), by pointwise
persistence of the strict winner angles and dominated
convergence. In particular all three arc lengths, interpreted
as \(2\pi\) times angular masses, are continuous there.

Consider a nonexceptional radius where each label is
a single arc of width strictly between zero and \(\pi\).
For a pair \(i,j\), the closed arc of spectator \(k\)
misses at least one of the antipodal intersections of
this circle with its exact reference line. At the missed
point, the two equal pair scores are maximal by continuity
from the available labels in an angular neighborhood.
The normal derivative of \(U_i-U_j\) on the reference line
is analytic and nonzero at the origin, so it can vanish
at only a discrete set of parameters. At radii avoiding
these zeros, the angular derivative of that contrast is
nonzero at the missed point. Since the spectator is absent
nearby on the circle, the two labels switch there and
their unique common endpoint lies on the reference line.

Thus at almost every such radius the three endpoints
choose one of the two antipodal rays on each reference
line. The three forward reference rays are \(120^\circ\)
apart. Flipping exactly one or exactly two rays produces
angular gaps \(\pi/3,\pi/3,4\pi/3\), so one width would
exceed \(\pi\). With zero or three flips, the labels are
exactly the reference or antipodal simplex. Indeed the
arc between two pair-endpoints has label their pairs have
in common.

Let \(I=(0,s_0)\) be the interval ending at the first
exceptional radius in \((0,R)\), or at \(R\) if none
exists. The set of radii in \(I\) where all widths are
strictly between zero and \(\pi\) is relatively open
and contains a neighborhood of zero. On that set all
widths equal \(2\pi/3\) almost everywhere by the endpoint
argument, hence everywhere by continuity. The set is
therefore also closed in \(I\) and equals \(I\).
The angular winner map is continuous there and takes
values in the two-point set of labeled reference and
antipodal partitions on a dense set of radii. It must
equal its germ throughout \(I\).

If \(s_0<R\), choose a pair whose contrast vanishes
identically on the exceptional circle and approach its
intersection with the active germ half-ray from below
in radius. Along this interval the two incident labels
are adjacent angular winners by the preceding propagation.
The spectator contrast \(h=U_i-U_k\) on the germ half-ray
is nonnegative there. It is analytic and nonzero near the
origin by the germ, so \(h>0\) sufficiently close to and
below \(s_0\). This gives an actual pair wall, and the
normal derivative of \(U_i-U_j\) cannot vanish in its
closure. Indeed away from the isolated zeros of this
normal derivative the wall is regular, and if it vanished
at the endpoint, the reciprocal of the pair-gradient
magnitude would have a nonintegrable \(1/(s_0-s)\) lower
bound on the regular wall. This violates (G20). But at
the intersection with the exceptional circle the pair
contrast vanishes on both the reference line and the
circle, so both independent tangential derivatives vanish.
The contradiction excludes an exceptional radius in
\((0,R)\) and proves the disk assertion.
[Active reference rays and reciprocal-gradient control](THREE_CELL_COMPACT_DEFECT_CENTRAL_TRIPOD.md#genuine-triple-points-on-the-reference-lines)

In particular, if a maximizing compact simplex discrepancy
has connected circular arcs almost everywhere in the plane,
take arbitrarily large \(R\). The exterior fixes the germ
orientation and the whole partition is the simplex. This
includes arbitrary initial widths and cyclic orders.

## A uniform obstruction to connected circular slices at right contact

Fix a compact correlation band \(J\subset(0,1)\).
There are constants \(R,\varepsilon>0\) depending only
on \(J\) such that every nonsimplex balanced planar maximizer
satisfying the exact-line and derivative hypotheses of (SO2)
has label indicators on \(B_R\) at least
\(\varepsilon\) away, in summed Gaussian \(L^1\), from
every partition into circular arcs on almost every
radius below \(R\). In particular this applies to a
nonsimplex maximizing compact discrepancy at a tie
approached from the right by simplex-optimal correlations.

Take \(R=L+1\), where \(L\) is the uniform excursion
window in (SO2). If the distance tended to zero, take
a sequence of maximizing partitions and approximating
circular-arc partitions on \(B_R\). Strong maximizer
compactness gives a limiting balanced maximizer with
correlation in \(J\). Rotate and relabel the sequence to
fix the reference lines. Uniform measurable local simplex
isolation keeps the limit nonsimplex. The exact lines,
equal prices and the derivative inequality pass to the
limit by strong convergence, as in (SO2).

The circular-arc class on \(B_R\) is closed in summed
Gaussian \(L^1\). To see this, pass to a subsequence
with angular \(L^1\) convergence on almost every radius
by Fubini. On one circle the indicators of arcs form
a compact set in angular \(L^1\), by parametrization
by center and length in \(S^1\times[0,2\pi]\), including
empty and full arcs. For triples of these indicators,
the condition that they are disjoint and sum to one
passes to the limit. Hence the limiting maximizer has
arc slices throughout \(B_R\) in the stated almost
everywhere sense. The disk assertion makes it the germ
simplex on \(B_R\). All its germ-arm scores then have
nonnegative spectator contrasts there, contradicting
the strict excursion within \(L\) from (SO2).

Equivalently, the integral over \((0,R)\) of the angular
\(L^1\) distance to the compact family of arc triples is
at least \(\varepsilon\). A countable dense subset gives
measurability and approximates this integral by actual
arc-slice comparisons: at each radius choose the first
member within a prescribed small error of the infimum.
The pointwise distance is at most two. Removing a small
central disk thus gives a fixed compact annulus on which
the distance is bounded below on a set of radii of
uniformly positive radial measure. The conclusion is
uniform over compact discrepancy radii.
[Uniform right-contact excursion, SO1–SO2](THREE_CELL_SPECTATOR_SOURCE_ORDER.md#a-uniform-excursion-on-one-of-the-three-germ-arms)
