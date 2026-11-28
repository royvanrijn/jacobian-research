# Gaussian asymptotics for three-cell ends isolated in wedges

Fix \(0<r<1\), put \(\sigma^2=1-r^2\), and let \(p_\sigma\)
denote the centered normal density with variance \(\sigma^2\).
Consider a stationary planar partition. In an orthonormal
coordinate system about the Gaussian origin, suppose almost
every vertical section of \(\{t>S,\ |z|<m t\}\) has one switch
between labels \(i,j\) at \(z=h(t)\), with
\(\operatorname*{ess\,sup}_{t>R}|h(t)|/t\to0\).
The sublinear graph theorem first gives an analytic
representative with bounded height. Write \(\beta\) for the pressure
of the lower label minus the upper, and let \(b\) solve
\(\beta=2\Phi((1-r)b/\sigma)-1\).

With \(R=r^2\) and \(\Lambda=R/(2\sigma^2)\), define
\[
    \kappa=\Lambda\frac{m^2}{1+m^2},
    \qquad L_\kappa=\frac{\kappa}{1+R-2\sigma^2\kappa}.
\]
Then for every \(0<\lambda<L_\kappa\),
\[
            h(t)-b=O(e^{-\lambda t^2}).                   \tag{WG1}
\]
In particular the nonlinear residual in the transformed
height equation decays at a strictly faster Gaussian rate
than \(\kappa\). In a network with one bounded affine-direction
graph per pair, this separates the leading effect of an acute
neighboring arm from the nonlinear height terms.
[Boundedness of sublinear isolated wedge graphs](THREE_CELL_ISOLATED_END_RIGIDITY.md#every-sublinear-isolated-graph-is-bounded)

## A wedge error and the rate it permits

Extend the two-label graph from \(s>S\) to all transverse
heights in that half-plane, and extend it by a bounded contrast
for \(s\le S\). The difference between this fictitious
contrast and the true contrast is supported in
\(\{s\le S\}\cup\{s>S,\ |v|\ge m s\}\), apart from null sets.
At an observation \((t,z)\) with \(z\) bounded, the
conditional mean of \(V-mS'\), where \(S',V\) are source
coordinates, is \(r(z-mt)\), and its variance is
\(\sigma^2(1+m^2)\). The analogous formula for \(-V-mS'\)
bounds the other side. Each side therefore contributes
\(O(e^{-k t^2})\) for every \(k<\kappa\), uniformly in
bounded \(z\). The contribution of \(S'\le S\) has
quadratic rate \(\Lambda>\kappa\), with a possible
linear shift. The same bounds hold after any fixed
number of score derivatives. Each derivative inserts a
polynomial in standard normal conditional coordinates,
and Hölder's inequality and Gaussian moments absorb
an arbitrarily small loss in the exponent.

The proof of (IE1)–(IE2) consequently holds with a
Gaussian-small error and its first derivatives on
bounded transverse observation sets. More explicitly,
the independent conditional Gaussian integral of the
fictitious two-label contrast gives the right-side
integral in (IE1). Its transverse derivative has uniform
strict sign on the graph since \(h\) is bounded. The
error derivative tends to zero. The implicit function
theorem gives analyticity and uniform derivative bounds
on a tail, by the bounds on the smoothed score derivatives.
Integration by parts gives (IE2), with the wedge error
still Gaussian-small. The extremal-slope and pressure
arguments of (IE2) then give \(h'\to0\) and \(h\to b\).
[The corresponding bounded half-plane computation](THREE_CELL_ISOLATED_END_RIGIDITY.md#bounded-ends-converge-to-the-pressure-offset)

Put \(\psi(u)=2\Phi(u/\sigma)-1\), \(d=(1-r)b\),
and \(v(t)=\psi(d+h(t)-b)-\psi(d)\).
The argument (TF1)–(TF6) gives one positive Gaussian
rate for \(v\). In (TF1), the new wedge error adds a
term \(O(e^{-k t^2})\) to the already present Gaussian
error. It is uniform for all late \(t\) in the section
where the moving longitudinal cutoff is used. The
bounded extension of \(v\) still bounds the residual
globally. The estimates (TF3)–(TF6) only change the
constants and the positive Gaussian error exponent,
so the same induction yields a positive rate.
[Initial Gaussian rate, TF1–TF6](THREE_CELL_TAIL_FLATNESS.md)

Use a fixed cutoff \(Q=v1_{\{t>S_0\}}\), with \(S_0\)
large enough for the estimates. Suppose that \(v\) has
rate \(\lambda>0\). For \(f=T_rQ-rQ\), the stationary
residual bound is a sum of a bounded-left term, the
longitudinal leakage \(C T_r1_{(-\infty,S_0]}\), the
nonlinear term \(C e^{-\mu(\lambda)t^2}\), and the wedge
error \(C_k e^{-k t^2}\) on the right. Here
\[
 \mu(\lambda)=\min\left\{2\lambda,\
            \lambda+\frac{2R\lambda}{1+4\sigma^2\lambda}\right\}.
\]
This is (TF7), with the extra error obtained by using
the fictitious graph in the pressure equation. In
(TF8)–(TF14), set \(d_k=\min\{\mu(\lambda),k\}\)
in place of \(\mu\). The telescoping step gives the
positive-real bound for \(P_1Q\) with rate
\(B=d_k/[R(1+2d_k)]\). The complex bound still has
rate \(a=\lambda/(1+2\lambda)\) on real parts.
Whenever \(\lambda<L_k=k/(1+R-2\sigma^2k)\), we have
\(d_k>\lambda\) and \(a<B<1/2\), so the same sector
and contour shift apply. In that range the shifted
rate for \(T_rQ\) strictly exceeds \(\lambda\), because
\[
              d_k(1+2\sigma^2\lambda)>(1+R)\lambda.
                                                               \tag{WG2}
\]
For the \(\mu\) branch this is (TF14). For the \(k\)
branch it is exactly \(\lambda<L_k\).
The residual rate \(d_k\) also exceeds \(\lambda\), so
the estimate for \(Q\) improves. The update from (TF8)
with \(d_k\) is continuous below \(L_k\) and strictly
increases the exponent there. Iteration therefore
reaches every exponent strictly below \(L_k\).
Letting \(k\uparrow\kappa\) proves (WG1).
[The shifted contour and update, TF8–TF14](THREE_CELL_TAIL_FLATNESS.md#improving-the-gaussian-rate-by-shifting-the-inversion-contour)

Finally \(\mu(L_\kappa)>\kappa\). Indeed
\(L_\kappa(1+R)=\kappa(1+2\sigma^2L_\kappa)\),
while (TF14) is strict with \(\lambda=L_\kappa\).
By continuity, some exponent \(\lambda<L_\kappa\)
already gives \(\mu(\lambda)>\kappa\).

## Leading interaction between the affine directions

Suppose outside a compact set the network consists of
exactly three ends, one for each pair, with distinct
outward directions. Suppose each is a bounded transverse
graph about its direction, in orthonormal coordinates
about the origin. Each has a sufficiently narrow wedge
as above. In the coordinates of a reference end, write
\(\theta_j\in(-\pi,\pi]\) for the relative directions of
the other two ends. Only those with \(0<|\theta_j|<\pi/2\)
enter its far forward half-plane. A perpendicular
direction has bounded longitudinal projection, and a
backward direction has projection tending to
\(-\infty\). Assume at least one acute direction enters,
and put
\[
 \theta_*=\min_{0<|\theta_j|<\pi/2}|\theta_j|,\qquad
 \kappa_*=\Lambda\sin^2\theta_*.
\]

Each entering arm \(j\) is a graph over the reference
longitudinal coordinate \(s\) on a late tail. The
longitudinal derivative in its own graph parameter tends
to \(\cos\theta_j>0\). Its limiting line in the reference
coordinates is
\[
    v\cos\theta_j-s\sin\theta_j=b_j,
\]
where \(b_j\) is its own pressure offset, with upward
normal \((-\sin\theta_j,\cos\theta_j)\).
For each sufficiently large \(s\), the difference of the
reference lower-minus-upper contrast from its binary graph
extension is a sum of single steps at the entering graphs.
An upper arm contributes with sign \(\eta_j=+1\) on its
far side, and a lower arm contributes with sign \(\eta_j=-1\).
To verify the signs, start between the reference graph
and the first upper graph, where the contrast is \(-1\).
Crossing toward larger \(v\) first reaches the spectator
contrast zero and, if another upper arm occurs, the
lower reference label with contrast \(+1\).
Starting below the reference graph and crossing toward
smaller \(v\) gives the reverse signs. Each statement uses
that each of the three pairs occurs exactly once outside
the compact set.

At reference observation \((t,b)\), the conditional
probability of the affine far side for arm \(j\) is
\[
    \Phi\!\left(-\frac{r|\sin\theta_j|t+c_j}{\sigma}\right),
    \qquad c_j=\operatorname{sgn}(\theta_j)
                   (b_j-r\cos\theta_j\,b).                \tag{WG3}
\]
Replacing an entering graph by this line produces a
strict improvement in the quadratic rate of its error.
In coordinates \((\ell,w)\) along that arm the discrepancy
is a strip about \(w=b_j\) of width \(O(e^{-\lambda_j\ell^2})\)
for some \(\lambda_j>0\), by (WG1). The two conditional
coordinates have independent variance \(\sigma^2\).
Their means at a bounded transverse observation are
\(r\cos\theta_j t+O(1)\) and \(-r\sin\theta_j t+O(1)\).
Thus the transverse density on the strip is
\(O(\exp[-\Lambda\sin^2\theta_j\,t^2+O(t)])\),
and integrating the width against the longitudinal
density gives a further strictly positive quadratic
rate. The extension of the affine far side through
\(s\le S_0\) only contributes the longitudinal Gaussian
error of rate \(\Lambda\). These estimates apply also
to fixed-order score derivatives, by the same polynomial
absorption. Evaluating at \(h(t)\) instead of \(b\) changes
(WG3) by a faster Gaussian term, by (WG1) and the mean
value theorem.

Take the reference wedge aperture arbitrarily close to
\(\theta_*\) from below. Formula (WG1) and strictness of
\(\mu(L_{\kappa_*})>\kappa_*\) give a reference exponent
\(\lambda\) with \(\mu(\lambda)>\kappa_*\).
In the stationary residual formula the affine
contributions (WG3) have the opposite sign, because the
two-label score variation must cancel the exterior
contribution. Dividing by the observation Taylor
prefactor \(1+O(h(t)-b)\) only gives a further
Gaussian improvement. Hence for some \(\epsilon>0\),
\[
 f(t)=-
       \sum_{|\theta_j|=\theta_*}
             \eta_j\Phi\!\left(-\frac{r\sin\theta_*\,t+c_j}{\sigma}\right)
       +O(e^{-(\kappa_*+\epsilon)t^2})
                    \quad(t\to+\infty).                   \tag{WG4}
\]
The remainder contains the faster-angle terms,
the nonlinear terms, and the discrepancy of the
affine replacements.

For \(P_1\) the variance-one normal convolution, the
Mehler telescoping identity (TF2) now gives
\[
 \begin{split}
 P_1Q(t)
    &=-\sum_{|\theta_j|=\theta_*}\eta_j
             \Phi\!\left(
              -\frac{\sin\theta_*\,t+c_j}
                     {\sqrt{\sigma^2+r^2\sin^2\theta_*}}
                  \right)
       +O(e^{-(B_*+\epsilon')t^2}),\\
 B_*&=\frac{\kappa_*}{R(1+2\kappa_*)},
 \qquad \epsilon'>0 .                                    \tag{WG5}
 \end{split}
\]
Indeed convolving each affine tail in (WG4) adds
\(r^2\sin^2\theta_*\) to the variance in its denominator.
The term of (TF2) with argument \(t/r\) gives the
displayed tails. Arguments \(t/r^k\) for \(k\ge2\) give
strictly faster quadratic rates. Extending the
remainder to the left contributes a bounded-left
Gaussian term under \(P_1\), and its positive-side
part has strictly faster rate by square completion.
The geometric sum in (TF2) preserves these bounds.
For a unique nearest arm, or two equally near arms
with different shifts \(c_j\), this determines the
eventual sign and the exact quadratic rate of
\(P_1Q\). The conclusion concerns the smoothed height
profile.

## A signed derivative obstruction at a wedge end

Continue to assume this three-end geometry, and use the
upper-minus-lower reference contrast \(F\) and pressure
\(\beta\) in this section. The affine offset \(b\) is
unchanged upon reversing the contrast. Assume \(\beta\ne0\)
and there is a unique nearest acute direction. Let \(E(t,z)\)
be the outside-wedge contribution from the entering arms,
defined by their step differences from the fictitious
binary graph in the far longitudinal half-plane. Let
\(R_0\) be the remaining bounded-left contribution to the
difference from the binary affine score. The sign of \(E\)
on the tail is given by the nearest affine step in (WG3),
with the sign now reversed from the preceding section.

Let \(g(t)=\left.\partial_\rho T_\rho
(1_{\rm upper}-1_{\rm lower})(t,h(t))\right|_{\rho=r}\),
holding the observation and source partition fixed.
Let \(g_\infty\) be this derivative for the affine binary
halfspace at its pressure line. Suppose the eventual sign
of \(E\) is opposite to \(\operatorname{sgn}\beta\).
Then no such end can satisfy either
\[
 \operatorname{sgn}(\beta)(g(t)-g_\infty)\le0
     \quad\text{eventually},
 \qquad\text{or}\qquad
 |g(t)-g_\infty|=O(e^{-(\kappa_*+\delta)t^2})
                  \quad\text{for some }\delta>0.          \tag{WG6}
\]

To prove this, set \(u=h-b\), \(\psi(a)=2\Phi(a/\sigma)-1\),
\(G(z)=\psi(rz-b)\), and \(U=F-G=W+E+R_0\), where \(W\)
is the right strip contribution in (CR2).
The differential equation (CR5) holds because (CR1)
uses only \(F(t,h(t))=\beta\) and \(LT_r=-r\partial_rT_r\).
Let \(a_0=r\psi'((r-1)b)>0\). We will use this equation with
its constant limiting coefficients before Gaussian smoothing.
[The signed strip coefficient and scalar equation, CR1–CR5](THREE_CELL_REFLECTED_CONSTANT_CORRELATION_RIGIDITY.md#the-signed-quadratic-term)

Take a rate \(\lambda\) from (WG1) for which
\(\mu(\lambda)>\kappa_*\). It can be chosen with
\(\lambda<\kappa_*\). Put \(A(\eta)=R\eta/(1+2\sigma^2\eta)\)
and \(\zeta=\min(\lambda,A(2\lambda),\kappa_*)\).
Then \(\lambda+\zeta>\kappa_*\).
The bound on \(u\) also holds for each fixed derivative
with an arbitrarily small loss in rate. In general, if
\(u=O(e^{-\lambda t^2})\) and all its fixed derivatives are
bounded, Taylor expansion at \(n\) consecutive points
\(t+j\epsilon\), \(0\le j<n\), gives
\[
 |u^{(m)}(t)|\le C_{m,n}
   \bigl(\epsilon^{-m} e^{-\lambda t^2}
                 +\epsilon^{n-m}\|u^{(n)}\|_\infty\bigr)
 \qquad(m<n).
\]
Use \(\epsilon=e^{-\lambda t^2/n}\) and increase \(n\) for
the desired rate. The derivatives of \(h\) are bounded
by the implicit score equation and the uniform transverse
gradient lower bound on the tail.

The strip \(W\) has linear term
\(-\psi'(rz-b)T_r\widetilde u(t)\), where
\(\widetilde u=u1_{\{t>S_0\}}\), and its remainder is bounded
by \(C T_r(\widetilde u^2)(t)\). The same estimates hold
after any fixed number of score derivatives, with a
small rate loss for horizontal derivatives due to
Gaussian polynomial factors. For the linear term,
use \(\widetilde u=v/\psi'((1-r)b)+O(v^2)\) on the
cutoff tail in the notation of (WG1).
Since \(T_rQ=rQ+f\) and (WG4) gives rate arbitrarily close
to \(\kappa_*\) for \(f\), \(T_r\widetilde u\) has rate
arbitrarily close to \(\zeta\). Its derivatives follow by
the same Taylor argument and uniform convolution
derivative bounds. Thus \(U_z,U_{tz},U_{zz}\) have
rates arbitrarily close to \(\zeta\), uniformly for
\(z\) in a fixed bounded neighborhood of \(b\).

In (CR5), \(F_z(t,h(t))-a_0=O(|u|+|U_z|)\).
The formulas for its coefficients give
\(p+t=O(|U_{tz}|+|u'|)\) and
\(q-1=O(|u|+|U_z|)\). Hence their products with
\(u,u',u''\) have rate strictly larger than \(\kappa_*\).
Also, the cancellation of the strip's linear term
in (CR4) holds for every bounded \(z\) near \(b\),
so the \(z\) derivative of \(\mathcal B_z W\) has the
same bound by \(C T_r(\widetilde u^2)\). Replacing
\(z=h(t)\) by \(b\) in \(\mathcal B_z W\) and
\(\mathcal B_zE\) therefore costs only a faster Gaussian.
The left term has rate \(\Lambda\). It follows that
for some \(\delta_1>0\),
\[
 -a_0(\partial_{tt}-t\partial_t+1)u(t)
   =\mathcal B_bW(t,b)+\mathcal B_bE(t,b)
       -r(g(t)-g_\infty)+O(e^{-(\kappa_*+\delta_1)t^2}).
                                                               \tag{WG7}
\]

Smooth this identity in \(t\) by \(P_1\), with \(u\)
smoothly extended by zero beyond a left cutoff. The
identity need only hold on a tail, since everything
to the left contributes at variance-one Gaussian
rate \(1/2\), with a possible linear shift.
The smoothed left side is
\(-a_0[-t(P_1u)'(t)+P_1u(t)]\). Indeed
\(P_1(t u')=t(P_1u)'+(P_1u)''\).
Moreover \(P_1u\) is a constant multiple of \(P_1Q\)
up to a smoothed \(O(u^2)\) error and terms from the
left cutoff. By (WG5), its rate exceeds
\(\kappa_* /(1+2\kappa_*)\), since
\(B_*=\kappa_* /[R(1+2\kappa_*)]\) and \(R<1\).
The smoothed quadratic error is also faster since
\(2\lambda>\kappa_*\).
A small rate loss from the derivative of \(P_1u\)
and the polynomial factor \(t\) preserves the strict
inequality of rates.

The signed terms on the right of (WG7) give a
contradiction to this rate. By (CR4),
\(\mathcal B_bW\) has the weak sign of \(\beta\) if
the strip cutoff is sufficiently large. The leading
nearest-arm contribution makes \(\mathcal B_bE\)
strictly opposite in sign to \(E\). In fact the
second vertical derivative of the affine unsigned
probability is a positive multiple of \(a\phi(a)\)
with \(a=(r\sin\theta_*\,t+O(1))/\sigma\).
Its coefficient contains
\(r^2\cos^2\theta_*/\sigma^2>0\), and it dominates the
first derivative and zeroth-order terms of the operator.
The other-arm and affine replacement errors have
strictly faster quadratic rates. Thus
\(\operatorname{sgn}(\beta)\mathcal B_bE(t,b)\)
is bounded below for large \(t\) by a positive
constant times \(e^{-(\kappa_*+\eta)t^2}\) for every
\(\eta>0\). After \(P_1\), retaining the late
half-line in this lower bound gives a positive term
of rate \((\kappa_*+\eta)/(1+2(\kappa_*+\eta))\), by
square completion. Take \(\eta\) sufficiently small
that the smoothed left side and remainder in (WG7)
are strictly faster than this rate. The left-cutoff
errors are faster as well. The signed \(g\) term in
(WG6) either has the same weak sign before smoothing
or contributes another faster error. This proves the
obstruction.

The signed obstruction applies to a tied maximizer
with zero deficit slope and curvature whenever
the mixed contact inequality makes \(g\) constant
along this interface. It does not require the acute
arm to be a reflected copy of the reference arm.
[Mixed contact inequality, G10](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)

## An angle condition for a rapidly converging correlation derivative

Assume a unique nearest acute direction as above, but allow any
pressure and either sign of the nearest-arm contribution.
Set \(m_0=2\) for nonzero pressure and \(m_0=3\) for zero
pressure. If
\[
                m_0r^2\cos^2\theta_*>1,                  \tag{WG8}
\]
then the end cannot satisfy
\(|g(t)-g_\infty|=O(e^{-(\kappa_*+\delta)t^2})\)
for any \(\delta>0\). In particular its correlation derivative
cannot be constant.

The derivation of (WG7) did not use a pressure sign. At
\(b=0\), the fixed-height strip integrand
\(\mathcal B_0[\psi(-u)-\psi(0)]\) in (CR4), with the
operator applied before setting \(z=0\), is odd as a
function of the source displacement \(u\).
Indeed \(\psi,\psi''\) are odd and the \(z\partial_z\)
term vanishes at \(z=0\). The linear term still cancels,
leaving an \(O(|u|^3)\) bound. At arbitrary nonzero
\(b\) the corresponding remainder is \(O(|u|^2)\).
Thus \(|\mathcal B_bW(t,b)|\le C T_r|\widetilde u|^{m_0}(t)\).
The coefficient estimates leading to (WG7) continue
to hold whenever \(\lambda\) is increased, since
\(\mu(\lambda)\) increases.

For a fixed \(\varepsilon>0\), if \(u\) currently has
rate \(\lambda\), (WG7) gives
\[
       |(\partial_{tt}-t\partial_t+1)u(t)|
         \le C_\varepsilon
                 e^{-\min(A(m_0\lambda),\,\kappa_*-\varepsilon)t^2},
       \qquad A(a)=\frac{R a}{1+2\sigma^2 a}.             \tag{WG9}
\]
This holds after dropping the faster terms and
absorbing the linear shifts in the nearest-arm Gaussian.
The scalar operator in (WG9) preserves any positive
Gaussian decay exponent of the forcing on the bounded
decaying solution used here. For \(t>0\), put \(u(t)=t v(t)\)
and \(F_0=(\partial_{tt}-t\partial_t+1)u\). The identity
\[
     (e^{-t^2/2} t^2 v'(t))'=e^{-t^2/2}t F_0(t)
\]
integrates twice from infinity, since \(u,u'\) already
have positive Gaussian rates. Explicitly,
\[
  u(t)=t\int_t^\infty \frac{e^{z^2/2}}{z^2}
                      \int_z^\infty e^{-s^2/2}s F_0(s)\,ds\,dz.
\]
If \(F_0=O(e^{-\eta t^2})\) for \(\eta>0\), this formula
bounds \(u\) by \(O(e^{-\eta t^2})\).

Under (WG8), \(A(m_0\lambda)>\lambda\) for
\(0<\lambda\le\kappa_*\), because
\(m_0R>1+2m_0\sigma^2\kappa_*\).
Iteration of (WG9) therefore reaches every rate strictly
below \(\kappa_*\) by taking \(\varepsilon\) sufficiently small.
In particular choose such a rate with \(A(m_0\lambda)>\kappa_*\).
The strip term in (WG7) now has faster quadratic rate
than the nearest-arm term. The smoothing argument following
(WG7) contradicts its signed leading term independently of
the pressure and leakage signs.

## Geometry when at least one wedge end is straight

Suppose the three-end network above belongs to an actual balanced
global maximizer. A straight reference arm cannot have any acute
or antipodal neighbor. Consequently, if the three directions
are not all in a closed semicircle, either they are all straight
with centered simplex exterior, or the partition has a global
reflection exchanging two labels. In the latter case one arm
is straight and the other two are curved and acute to each other.
A competing geometry with directions in a closed semicircle must
have all directions in an open semicircle, both successive gaps
acute, and all three arms curved.

To exclude an acute neighbor of a straight reference, use the
lower-minus-upper convention of (WG4)–(WG5). The transformed
height \(Q\) is zero with a late cutoff. Thus (WG5) forces two
acute neighbors with equal relative angles in magnitude and
equal shifts \(c_+=c_-\). With only one leading term, the
Gaussian tail cannot be canceled by the faster remainder.
With two unequal shifts, the smaller shift dominates by the
Gaussian tail ratio. The two angles must therefore be
\(\theta,-\theta\), with \(0<\theta<\pi/2\). In increasing
transverse order on a far forward slice, the label pattern
is \(k,i,j,k\), where \(i,j\) are the lower and upper reference
labels and \(k\) is the spectator.

The common shift gives \(b_++b_-=2r\cos\theta\,b\).
On the other hand the three lower-minus-upper pressures
sum to zero by their cyclic incidences on this slice.
Since \(x\mapsto 2\Phi((1-r)x/\sigma)-1\) is odd and strictly
increasing, the sign of the sum of its values at \(b_+,b_-\)
equals the sign of \(b_++b_-\).
The pressure sum then gives the opposite sign to \(b\).
Hence \(b=0\), and the reference pressure is zero.

At a straight interface the pair score is constant on its
entire supporting line. Injectivity in the longitudinal
observation gives the conditional slice equation (ER1).
On a far positive slice the neighboring heights
\(v_-(s)<0<v_+(s)\) cut the sequence \(k,i,j,k\).
The zero slice contrast is equivalent to
\(\Phi(v_-(s)/\sigma)+\Phi(v_+(s)/\sigma)=1\).
Thus \(v_-=-v_+\), with pointwise equality by continuity.
This gives reflection across the reference line for all
sufficiently large positive longitudinal coordinates.
For sufficiently negative longitudinal coordinates the
spectator fills the entire half-plane, since all three
directions point forward. In the intermediate bounded
longitudinal slab the far upper and lower parts also have
the spectator label. These facts follow from the exact
three-arm description and the cyclic label pattern.
So the partition is reflected up to a compact discrepancy.
[The slice identity and compact contrast uniqueness, ER1–ER3](THREE_CELL_ONE_STRAIGHT_END_REFLECTION.md)

If both remaining arms are straight, either of them has
the original reference direction as its unique nearest acute
neighbor. The angle to the other arm is \(2\theta>\theta\).
Then (WG5) on that straight arm is impossible. If the remaining
arms are curved, the compact contrast uniqueness argument
applies on both of them to the discrepancy between the
partition and its reflected, relabeled copy. Its hypotheses
hold by (WG1) and the exterior reflection just proved.
The two contrasts determine all three indicators, so the
partition is globally reflected. The single-reflection reduction
then places the straight \(i/j\) arm on the reflection axis and
the other two arms on the graph \(y=h(|x|)\) with \(h\) increasing,
in coordinates in which the axis arm points down.
Neither of their outward directions can project positively
onto the axis arm. This contradicts the acute angles.
[The global reflected graph](THREE_CELL_REFLECTION_REDUCTION.md)

Now suppose the straight reference direction and another
are antipodal. The third direction cannot be acute to the
straight reference. If it is perpendicular to both, all three
are isolated in full half-planes, so all three are straight.
The antipodal straight-line pressure contradiction then
applies. In the remaining case take the straight reference
eastward, its antipodal arm westward, and the third
upper-left after reflecting if needed. On far-left physical
slices, write \(H(s)\to p\) for the west arm's height and
\(V(s)\to+\infty\) for the upper-left arm's height. In increasing
transverse order the labels are \(i,k,j\), where \(i,j\) are the
lower and upper labels of the straight reference. The slice
equation reads
\[
     \Phi((H(s)-rb)/\sigma)
       -\bigl[1-\Phi((V(s)-rb)/\sigma)\bigr]=\beta.
\]
Taking \(s\to-\infty\) gives \(\beta=\Phi((p-rb)/\sigma)\).
The positive upper tail therefore forces \(H(s)>p\)
on all sufficiently late slices. In oriented west-arm
coordinates \((\ell,w)=(-s,-v)\), its displacement from
the pressure line is strictly negative. Its transformed
height \(Q_{\rm west}\) can be chosen nonpositive by taking
a late cutoff. The upper-left arm is a unique lower
acute neighbor in these west-arm coordinates. Formula
(WG5) gives \(P_1Q_{\rm west}>0\) eventually, since its
step coefficient is \(-1\). This contradicts positivity
of Gaussian convolution on nonnegative functions and
completes the antipodal exclusion.
[Full half-plane straightening](THREE_CELL_TAIL_FLATNESS.md#reaching-the-straightness-threshold)
[Antipodal straight-line pressures](THREE_CELL_ONE_STRAIGHT_END_REFLECTION.md#the-antipodal-case)

For the resulting geometry, first suppose the directions
are not contained in a closed semicircle. Their three
successive gaps are all less than \(\pi\), sum to \(2\pi\),
and at most one is acute. With no acute gap, each reference
forward half-plane eventually contains no other arm, so all
three are straight. Then (PG10) gives the centered simplex
exterior. If one gap is acute, the other two gaps are
strictly obtuse, and the opposite arm is isolated in
a full forward half-plane. It is straight by (TF7)–(TF14).
In its far reverse direction the spectator fills a cone,
so the straight offset and pressure are zero.
Both other arms project strictly backward and cross far
negative longitudinal slices at opposite transverse ends.
The cyclic labels make those slices consist of opposite
reference-label tails with spectator between.
The slice calculation (ER2) gives exact reflection on
those slices. It holds also in the straight forward
half-plane and in the far upper and lower parts of the
intermediate slab, because their labels are fixed there
by the exterior graph arrangement. Thus it holds outside
a compact set. The other two arms must both be curved by
the acute exclusion, and compact contrast uniqueness
makes the reflection global as before.
[Centered simplex exterior, PG10](THREE_CELL_PRESSURE_GRADIENT_TAIL.md#exterior-geometry-when-all-three-ends-straighten)

If the directions belong to an open semicircle and an
arm is straight, it has no acute or antipodal neighbor.
That arm must be at an extreme of the angular interval
containing the three directions, since a middle arm would
otherwise have both gaps at least \(\pi/2\), making the
interval at least \(\pi\). The reverse of an extreme
direction lies strictly in the exterior angular sector
adjacent to it through the gap greater than \(\pi\).
This sector has one of its two incident labels.
Analytic continuation of the constant score along the
straight supporting line would then make its pressure
\(+1\) or \(-1\) from the reverse limit, impossible.
Thus no arm is straight in the open-semicircle case.
If either successive gap were at least \(\pi/2\), the
extreme beyond that gap would have no acute neighbor and
would straighten by half-plane isolation. So both gaps
are acute. Finally, directions contained in a closed but
no open semicircle have an antipodal pair. At least
one member of that pair has no acute neighbor, including
when the third direction is perpendicular to both.
It straightens by half-plane isolation, contradicting
the antipodal exclusion.
## Pressure compatibility at higher-order contact

Suppose an actual balanced planar maximizer \(A\) ties the
regular simplex \(T\) at \(r\), with simplex-optimal
correlations accumulating at \(r\). Put
\(D(s)=S_s(T)-S_s(A)\). Assume \(D'(r)=D''(r)=0\)
and suppose there is a bounded isolated wedge graph end
for each of the three label pairs. This pressure conclusion
does not need a global count of the ends. Either \(A\)
is a regular simplex or exactly two labels have the same
stationary price.

By (G10), the mixed fixed-partition score derivative \(g\)
annihilates all compact balanced regular interface speeds.
The mass-flux map is onto the centered label space by
complete adjacency and regular spanning-tree bumps.
Thus \(g\) is a difference of constant label numbers on
every regular pair interface. In cyclic label orientations
its three pair values sum to zero. The stationary pressure
differences \(\beta\) also sum to zero.
[Mixed contact inequality (G10)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)

Along any of the wedge ends, the constant derivative \(g\)
equals the affine binary limit
\[
             g_\infty(\beta)=-\frac{2z\phi(z)}{1-r^2},
             \qquad z=\Phi^{-1}((1+\beta)/2).               \tag{WG10}
\]
To justify the limit, compare the fixed partition to the
binary labeling along its affine pressure line at height
\(b\). The source differences are a Gaussian-thin strip
in the wedge by (WG1), and the complement of the wedge.
Both have Gaussian-small conditional probability at the
graph on a tail, by the same strip and wedge estimates
as above. Differentiation of the conditional density
inserts a polynomial with conditional moments growing at
most polynomially in the observation. Cauchy–Schwarz gives
vanishing derivative error on the tail. On the affine
line, with \(b\) held fixed while the correlation varies,
the binary contrast at correlation \(s\) equals
\(2\Phi((1-s)b/\sqrt{1-s^2})-1\). Its derivative at
\(s=r\) is (WG10).

The odd function \(F(\beta)=2z\phi(z)\) vanishes at zero
and has
\[
      F'(\beta)=1-z^2,\qquad F''(\beta)=-z/\phi(z).
\]
It is strictly concave on the positive half of its
domain. Consequently \(F(a)+F(b)>F(a+b)\)
for \(a,b>0\) with \(a+b<1\). Three nonzero cyclic
pressures summing to zero would, up to a common sign,
have the form \(a,b,-(a+b)\). Equation (WG10) would
then make their derivatives sum to a nonzero number,
a contradiction.

If all pressures vanish, (WG10) gives \(g=0\)
on the full regular interface. The dilation inequality
(G6) forces the dilation speed to vanish, hence
centered conicality and the simplex. Otherwise the
pressure triple has one zero and two opposite nonzero
values. The zero-pressure wedge end has \(b=0\)
and is subject to the cubic derivative test (WG9).
[Dilation inequality (G6)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[Centered conical rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#convex-and-conical-maximizers)

In the exact-three-end geometry above, a nonsimplex
higher-order tie has all directions in an open semicircle
and all three ends curved. Indeed, positively spanning
directions force a simplex exterior or global reflection.
The simplex exterior has all pressures zero, already
excluded at a nonsimplex higher-order tie by the preceding
argument. Global reflection is excluded by the
reflected contact curvature theorem.
[Higher-order ties in the reflected class](THREE_CELL_FIRST_CONTACT_RIGIDITY_TESTS.md#reflected-ties-have-at-most-second-order-contact)
