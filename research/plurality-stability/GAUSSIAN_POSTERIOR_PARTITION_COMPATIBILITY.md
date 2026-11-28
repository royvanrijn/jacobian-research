# Conditioning and compatibility of Gaussian partitions

Fix \(k\ge3\) and \(0<s<\rho<1\). Reveal a Gaussian component of variance
\(s\) shared by two samples of correlation \(\rho\). The posterior
partitions of a centered simplex are translated simplices, with label
masses that vary with the revealed point. Optimizing over a separate
partition at each revealed point, even with those exact posterior masses,
gives strictly greater agreement than the centered simplex. This holds
for every reveal time and every positive correlation. A compatibility
constraint across revealed points is therefore essential in using
posterior optimization to compare fixed partitions.

## The conditional comparison

Work first in dimension \(k-1\), the active space of the simplex. Let
\(U\sim N(0,sI)\) be independent of a pair \(Z,Z'\) of standard Gaussians
with correlation
\[
                       r=\frac{\rho-s}{1-s}.
\]
Then \(X=U+\sqrt{1-s}\,Z\) and \(Y=U+\sqrt{1-s}\,Z'\) are standard
Gaussians with correlation \(\rho\). For a fixed partition \(A\), put
\(A^u=(A-u)/\sqrt{1-s}\), meaning the same operation on each cell. Direct
conditioning gives
\[
                  S_\rho(A)=\mathbb E_U S_r(A^U).                \tag{PC1}
\]
For the centered simplex \(T\), the cone property gives
\(T^u=T-u/\sqrt{1-s}\). Write \(p(u)\) for its mass vector. This is a
continuous function of \(u\), since the polyhedral interfaces are
Gaussian-null and their translations converge in measure.

Use coordinates \(H=\mathbf 1^\perp\subset\mathbb R^k\), with the cones
\(T_i=\{x\in H:x_i\ge x_j\ \forall j\}\). For a translation
\(a\in H\), the masses of \(T+a\) are the winning probabilities of
independent standard normals with scores \(G_i-a_i\), since projection
of \((G_i)\) onto \(H\) changes all scores by the same amount. If
\(a_i<a_j\), swapping \(G_i,G_j\) sends every \(j\)-winning point to
an \(i\)-winning point. The inclusion is strict in probability: take both
scores before the shifts sufficiently large, the raw coordinates \(G_i,G_j\)
sufficiently close, and all others sufficiently small. Thus the masses
are all equal exactly when \(a=0\).

Fix \(u_0\ne0\). By Heilman–Mossel–Neeman, Theorem 2.6, a shifted flat
partition into \(k\) cells in dimension \(k-1\) with unequal masses admits
a strictly more stable competitor at \(0<r<1\) and the same masses.
In particular there is an ordinary \(B\) of masses \(p(u_0)\) with
\(S_r(B)>S_r(T^{u_0})\).
[Heilman–Mossel–Neeman, Theorem 2.6](https://arxiv.org/abs/1403.0885)

To make a measurable family on a neighborhood of \(u_0\), repair the
masses of \(B\) to \(p(u)\). For target masses \(p\), send
\[
     \frac{(p_i(u_0)-p_i)_+(p_j-p_j(u_0))_+}
          {\sum_\ell (p_\ell-p_\ell(u_0))_+}
\]
units of mass from each surplus cell \(i\) to each deficit cell \(j\).
Keep \(B\) when the denominator is zero. The transfers can be made by
consecutive intervals of the Gaussian rank of the first coordinate
within each source cell, with the recipients in fixed label order.
The first coordinate has a nonatomic distribution in each cell, so this
rank assigns mass exactly and the family is jointly measurable. Its
Gaussian \(L^1\) distance from \(B\), summed over the indicators, is
\(\|p-p(u_0)\|_1\). Changing a partition on a set of mass \(b\) changes
its stability by at most \(2b\), by the union bound at the two Gaussian
endpoints. The repaired family thus beats \(T^u\) by a positive uniform
amount on a sufficiently small neighborhood of \(u_0\).

Use the repaired family on that neighborhood and \(T^u\) elsewhere.
This family \(B_u\) has the exact posterior mass vector \(p(u)\) at all
revealed points. The rules \(B_U(Z)\) and \(B_U(Z')\), interpreted as label
values, have uniform label marginals by (PC1) and Gaussian invariance.
Their agreement is strictly greater than \(S_\rho(T)\), since the
neighborhood has positive \(U\)-probability. The construction can be
extended cylindrically to every \(d\ge k-1\). The extra information
available to each rule is the revealed \(U\).

## The degree-one obstruction

There is also a direct obstruction at the level of conditional centroid
energy. For any partition \(A\), let \(m_i^A=\int_{A_i}x\,d\gamma(x)\) and
\(W_1(A)=\sum_i|m_i^A|^2\). For every nonzero \(a\in H\), the translated
simplex \(T+a\) has a strictly improving same-mass competitor for \(W_1\).

Let \(w_{ij}(a)\) be the Gaussian surface area of the shared facet of
\(T_i+a,T_j+a\). The divergence theorem for Gaussian density gives
\[
                  m_i^{T+a}
     =\frac1{\sqrt2}\sum_{j\ne i}w_{ij}(a)(e_i-e_j).
                                                               \tag{PC2}
\]
This can be applied first in a ball and then taken to the limit, since the
Gaussian-weighted area of the outer sphere tends to zero. For three
distinct labels \(i,j,\ell\), swapping \(j,\ell\) carries the unshifted
\(i,j\) facet to the unshifted \(i,\ell\) facet. On the former facet,
\[
  a\cdot x-a\cdot \sigma_{j\ell}x=(a_j-a_\ell)(x_j-x_\ell).
\]
Here \(x_j>x_\ell\) on its relative interior. After factoring out the
common Gaussian density \(e^{-|x|^2/2-|a|^2/2}\), the two remaining
weights are \(e^{-a\cdot x}\) and \(e^{-a\cdot\sigma_{j\ell}x}\).
Thus \(a_j>a_\ell\) implies \(w_{ij}(a)<w_{i\ell}(a)\).

If \(T+a\) maximized \(W_1\) at its own masses, the linear score
\((m_i^{T+a}-m_j^{T+a})\cdot x\) would be constant on every open shared
facet. Otherwise one can exchange equal small Gaussian masses in two
one-sided neighborhoods of the facet where the score values are strictly
ordered. Assign the higher score to \(i\) and the lower score to \(j\).
The resulting first-order change \(2\sum_i m_i^{T+a}\cdot\Delta m_i\)
is positive, and the quadratic change \(\sum_i|\Delta m_i|^2\) is
nonnegative. Constancy on the facet implies that \(m_i^{T+a}-m_j^{T+a}\)
is parallel to \(e_i-e_j\), so (PC2) gives \(w_{i\ell}(a)=w_{j\ell}(a)\)
for each third label \(\ell\). The strict comparison of facet weights
with the common label \(\ell\) forces \(a_i=a_j\). This holds for every
pair, and \(a\in H\) then forces \(a=0\).

In the Brownian representation, let \(B_t\) be standard Brownian motion
and \(P_t=(\Pr(B_1\in A_i\mid B_t))_i\). Differentiating its Gaussian
conditional probabilities gives the vector martingale identity
\[
       S_\rho(A)=\|\mathbb E P_0\|^2+
          \mathbb E\int_0^\rho
             \frac{W_1((A-B_t)/\sqrt{1-t})}{1-t}\,dt .
                                                               \tag{PC3}
\]
Indeed \(\nabla_x P_{t,i}\) at \(x=B_t\) is the conditional-cell centroid
divided by \(\sqrt{1-t}\), so the trace of the quadratic variation is
the integrand. The derivatives are bounded on any time interval ending
strictly before 1. This is the coordinate sum of Eldan's martingale
identity. The strict centroid comparison shows that the simplex does
not maximize this instantaneous trace at its nonuniform posterior
masses, even though its initial masses are balanced.
[Eldan, Section 2.2, equations (10)–(13)](https://arxiv.org/abs/1307.2781)

## A quantitative compatibility bound

For an arbitrary measurable family \(B_u\) of ordinary partitions set
\(L(U,X)=e_{B_U((X-U)/\sqrt{1-s})}\), and put \(F(X)=\mathbb E[L(U,X)\mid X]\).
Use the same definitions at \(Y\). Define
\[
                 \mathcal E=\mathbb E(1-|F(X)|^2).
\]
This quantity is zero exactly when the adaptive rule is a deterministic
label of \(X\) almost surely. Indeed \(\mathcal E\) is the conditional
variance of the unit-vector label integrated over \(X\). The Gaussian
density of \((U,X)\) is everywhere positive. Consequently, in terms of
the family \(B_u\), zero defect is equivalent to the equality of
\(B_u((x-u)/\sqrt{1-s})\) and \(B_{u'}((x-u')/\sqrt{1-s})\) for almost
every \(u,u',x\).

Writing \(S^{\rm ad}\) for the adaptive agreement, one has
\[
                  |S^{\rm ad}-S_\rho(F)|\le 2\mathcal E .
                                                               \tag{PC4}
\]
To prove this, expand the difference of the pairwise dot products as
\((L(U,X)-F(X))\cdot L(U,Y)+F(X)\cdot(L(U,Y)-F(Y))\).
A centered vector paired with a simplex-valued vector has absolute dot
product at most half its \(\ell^1\)-norm. Conditional on \(X\), this half
norm for \(L(U,X)-F(X)\) has expectation \(1-|F(X)|^2\). The same holds
at \(Y\), proving (PC4).

In particular, if the adaptive family has balanced expected masses,
\(F\) is a feasible fractional balanced partition. The optimum among
fractional balanced partitions equals the ordinary optimum \(V_\rho\)
in the same finite dimension, by weak compactness, compactness of the
positive OU kernel, and the strict-convexity split on fractional
coordinates. Thus \(S^{\rm ad}\le V_\rho+2\mathcal E\).
[Existence and ordinary extremizers](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md#optimization-and-volume-repair)

When the family has the specific posterior masses \(p(u)\) of \(T\),
compatibility determines the endpoint partition, with a quantitative
version. Put \(f^T=(1_{T_i})_i\). Since the conditional law of \(X\)
given \(U=\sqrt{s}\,u\) is \(N(\sqrt{s}\,u,(1-s)I)\), the mass error of
the averaged endpoint rule satisfies
\[
 \|T_{\sqrt{s}}(F-f^T)\|_{L^2(\gamma)}^2
   =\mathbb E\big|\mathbb E[F(X)-L(U,X)\mid U]\big|^2
   \le \mathcal E.                                             \tag{PC5}
\]
Thus \(\mathcal E=0\) implies \(F=f^T\), by injectivity of the
positive OU multipliers. For \(\mathcal E>0\), use Hermite
projections \(P_j\) and set \(\theta=\log(\rho)/\log(s)\in(0,1)\).
Hölder's inequality and the squared diameter \(2\) of the simplex
give
\[
 \begin{split}
 \|T_{\sqrt{\rho}}(F-f^T)\|_2^2
 &=\sum_j \rho^j\|P_j(F-f^T)\|_2^2\\
 &\le
 \left(\sum_j s^j\|P_j(F-f^T)\|_2^2\right)^\theta
       \|F-f^T\|_2^{2(1-\theta)}
 \le \mathcal E^\theta 2^{1-\theta}.
 \end{split}
\]
Both \(T_{\sqrt{\rho}}F\) and \(T_{\sqrt{\rho}} f^T\) have \(L^2\)
norm at most one. Subtracting their squared norms and using (PC4)
therefore yields
\[
       |S^{\rm ad}-S_\rho(T)|
             \le 2\mathcal E+
                  2^{(3-\theta)/2}\mathcal E^{\theta/2}.         \tag{PC6}
\]
The bound uses only the observed posterior masses and the endpoint
compatibility defect. It applies throughout the positive-correlation
range, independently of simplex optimality for fixed balanced
partitions.

## An adaptive martingale satisfying the full posterior perimeter comparison

The gain in (PC1) can be realized by a continuous posterior
martingale with the correct balanced initial state and ordinary
terminal labels. Each instantaneous covariance comes from an
ordinary conditional Gaussian partition, and every two-time sharp
perimeter-profile inequality holds. These conditions still permit
strictly larger squared posterior mass than the centered simplex at
any specified positive correlation.

Fix \(0<s<\rho<1\). Let \(X_t\) be standard Brownian motion in
the simplex active space, with its usual filtration
\(\mathcal F_t\). Use the family \(B_u\) constructed after (PC1):
it has the exact conditional simplex masses \(p(u)\), equals
\(T^u\) outside a bounded ball \(W\) of revealed points, and beats
\(T^u\) on \(W\) at correlation \((\rho-s)/(1-s)\).
Set \(\widetilde B_u=u+\sqrt{1-s}B_u\), meaning that the
transformation is applied to every cell, and let
\[
 L=e_{\widetilde B_{X_s}(X_1)},\qquad
 P_t=\mathbb E[L\mid\mathcal F_t].
\]
The label \(L\) is an ordinary label in the full Brownian filtration.
For \(t\le s\) the posterior \(P_t\) is exactly the simplex
posterior by the tower property, since the label masses conditional
on \(X_s=u\) are the simplex masses there. In particular
\(P_0=(1/k,\ldots,1/k)\). For \(s\le t<1\), conditioned on
\(X_s=u\) and \(X_t=x\), its coordinates are
\[
 P_{t,i}(u,x)=
     \int 1_{\widetilde B_{u,i}}\bigl(x+\sqrt{1-t}\,z\bigr)
                                                         \,d\gamma(z).
                                                                  \tag{PC7}
\]
At \(t=s\) the two descriptions agree by exact conditional mass
repair. Gaussian convolution in (PC7) is continuous in \(x,t\) for
each fixed \(u\) at times below one. The paths therefore match
continuously at \(s\), and the martingale converges to \(L\) at one,
since \(\mathcal F_{1-}=\mathcal F_1\) for continuous Brownian paths.

To calculate the covariance density at a time \(t<1\) away from
the restart, use the ordinary conditional partition
\((D-X_t)/\sqrt{1-t}\), where \(D=T\) before the restart and
\(D=\widetilde B_{X_s}\) after it. If \(M_t\) is its matrix of
Gaussian centroids, differentiation of the conditional Gaussian
formula gives \(dP_t=M_t\,dX_t/\sqrt{1-t}\). Indeed the conditional
probability maps satisfy the backward heat equation, so the drift
in Itô's formula vanishes. This proof is local to intervals away
from the restart and terminal time. Continuity leaves no separate
quadratic variation atom at the restart. In particular all
pointwise tests on ordinary centroid matrices at masses \(P_t\)
hold, including the weighted perimeter test (WC1).
[Weighted centroid inequality](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#a-sharp-weighted-centroid-inequality-at-arbitrary-masses)

Write \(I\) for the Gaussian multi-bubble perimeter profile.
For all deterministic \(0\le v<w<1\),
\[
 \mathbb E[I(P_w)\mid\mathcal F_v]\ \ge\
           \sqrt{\frac{1-w}{1-v}}\,I(P_v).                       \tag{PC8}
\]
If both times are before the switch, or both after it, the future
endpoint uses one fixed partition conditional on \(\mathcal F_v\).
At \(v=s\) the restarted partition is already known. Apply the
sharp posterior perimeter comparison (PF4) to the partition
recentered at \(X_v\) and scaled by \(1/\sqrt{1-v}\). Its posterior
at \(w\) uses Gaussian smoothing parameter
\(\sqrt{(w-v)/(1-v)}\), so (PF4) gives the displayed factor.
If \(v<s<w\), apply the two estimates from \(v\) to \(s\) and
from \(s\) to \(w\), then use the tower property and the
multiplicative scale factors. At \(w=1\) the corresponding
comparison has both sides zero, since \(P_1\) is a vertex and
the perimeter of that one-cell assignment is zero.
[Sharp comparison along the indicator flow](GAUSSIAN_SCORE_SUPPORT_FLOW.md#sharp-comparison-along-the-indicator-flow)

On the other hand, conditional on the reveal \(X_s=u\), two
independent continuations after time \(\rho\) yield the stability
of \(B_u\) at correlation \((\rho-s)/(1-s)\). Their agreement is
the mean squared norm of the time-\(\rho\) posterior. Thus
\[
               \mathbb E|P_\rho|^2
       =\mathbb E_U S_{(\rho-s)/(1-s)}(B_U)>S_\rho(T).             \tag{PC9}
\]
For any fixed partition the endpoint label would instead be a
function only of \(X_1\), and its posterior would depend only on
the current position \(X_t\). The adaptive martingale satisfies
all the two-time comparisons (PC8) and the ordinary centroid
condition at the same time. Even their combination with
martingale evolution therefore cannot prove the desired squared
posterior comparison without a further restriction on how
conditional partitions at different histories arise from the
same endpoint rule.

### Detecting the endpoint restriction at a single time

For \(s<t<1\), define the history defect
\[
             R_t=\mathbb E\big|P_t-\mathbb E[P_t\mid X_t]\big|^2 .
\]
For a family with unchanged branches off \(W\), \(R_t=0\) at
any one such time if and only if \(\widetilde B_u=T\) for
almost every revealed \(u\), up to null sets. Indeed the joint
Gaussian density of \((X_s,X_t)\) is everywhere positive.
Vanishing defect means that the Gaussian-smoothed indicators
in (PC7) agree as functions of \(x\) for almost every \(u\).
A branch off \(W\) identifies this common function as the
smoothed reference simplex. Injectivity of the Euclidean
Gaussian convolution on bounded functions, by Fourier
multipliers, makes the endpoint indicators themselves agree.
The reverse implication is immediate. At the restart time
this test gives zero even for a nontrivial restart, because
the conditional masses agree. The distinction is detected
at every later time.

There is a quantitative bound for every bounded restart ball
\(W\) of radius \(R\). Let \(\tau_t=\sqrt{s(t-s)/t}\) and
\(q_t=2\Phi(-R/\tau_t)>0\), where \(\Phi\) is the standard
normal CDF. The conditional law \(X_s\mid X_t=x\) has
covariance \(\tau_t^2 I\). Its probability of being outside
\(W\) is at least \(q_t\), by bounding the ball by a slab
and using that a Gaussian interval of given length has
largest mass when centered at the mean. Put
\(P_t^T=\mathbb E[e_{T(X_1)}\mid X_t]\).
At fixed \(X_t=x\) the vector \(P_t-P_t^T\) is zero off the
restart ball. Hence Cauchy–Schwarz gives
\[
 \big|\mathbb E[P_t-P_t^T\mid X_t=x]\big|^2
 \le (1-q_t)\mathbb E[|P_t-P_t^T|^2\mid X_t=x].
\]
Subtract this squared mean and integrate to obtain
\[
 R_t\ge q_t\,\mathbb E|P_t-P_t^T|^2,\qquad
 \big|\mathbb E|P_t|^2-S_t(T)\big|
                    \le 2\sqrt{R_t/q_t}.                        \tag{PC10}
\]
Thus a strict increase in agreement after the bounded restart
requires a positive single-time defect from the fixed-endpoint
description. For the adaptive example that defect is distinct
from the perimeter-profile and centroid constraints, both of
which remain valid.
