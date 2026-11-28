# A distant tie despite strict measurable simplex stability for a modified Gaussian objective

On balanced three-cell partitions of the Gaussian plane there is a
one-parameter family of compact strictly positive definite quadratic
objectives invariant under orthogonal maps and relabeling. The regular
simplex starts as the unique global maximizer and remains uniformly
strictly locally maximizing in Gaussian label distance throughout
the family. Nevertheless at an interior parameter it ties a
nonsimplex global maximizer at positive label distance, and at a
larger parameter it is strictly suboptimal.

The construction keeps a positive Ornstein–Uhlenbeck term and adds
a rank-one radial term with zero score at the simplex. The objective
is a modified quadratic functional. Its coexistence point gives
a counterexample to a continuation argument based solely on
positive definiteness, compactness and strict simplex local growth.

## The objective and a distant competitor

Fix \(r_0=1/4\), let \(q(x)=|x|^2/2\), and put
\[
 g(x)=1_{\{q(x)\ge\log 3\}}-\frac13,\qquad
 \mathcal W_g(A)=\sum_i\left(\int g1_{A_i}\,d\gamma_2\right)^2.
\]
Write \(W_1(A)=\sum_i|\int x1_{A_i}\,d\gamma_2|^2\) and define
\[
      \mathcal E_\lambda(A)=W_1(A)+\lambda\mathcal W_g(A)
                               +\frac1{100}S_{r_0}(A),
           \qquad 0\le\lambda\le\frac{11}{2}.                  \tag{RK1}
\]
The kernel operator is the sum of the degree-one Hermite projection,
\(\lambda g\otimes g\), and \(T_{r_0}/100\). It is compact on
\(L^2(\gamma_2)\) and its quadratic form is strictly positive on
every nonzero vector. The centered regular simplex \(T\) has
\(\mathcal W_g(T)=0\) because its cells have the same mass on every
centered circle.

At \(\lambda=0\), the simplex uniquely maximizes (RK1). Here is the
three-cell first-moment argument. For \(m_i=\int x1_{A_i}\,d\gamma_2\),
write \(W=\sum_i|m_i|^2\). Pointwise comparison with a maximum gives
\(W\le\mathbb E\max_i m_i\cdot X\). For three real numbers, half
the sum of their three pairwise absolute differences is the
maximum minus the minimum. Gaussian sign symmetry therefore
gives
\[
 \mathbb E\max_i m_i\cdot X
 =\frac1{2\sqrt{2\pi}}\sum_{i<j}|m_i-m_j|
 \le \frac{3\sqrt W}{2\sqrt{2\pi}},
\]
where \(\sum_i m_i=0\) was used in Cauchy–Schwarz. Thus
\(W\le9/(8\pi)\). The regular sector partition attains equality
by direct integration. Equality requires all three moment
differences to have equal lengths and each label to select its
largest linear score, which makes it the centered simplex.
Also \(S_{r_0}\) is globally maximized by \(T\).
[Exact three-cell noise stability on the lower interval](THREE_CELL_CENTRAL_SHIFT_BOOTSTRAP.md)

For the radial competitor \(B\), take the three intervals for \(q\)
cut at \(\log(3/2)\) and \(\log3\). Since \(q\) is exponential
with mean one, \(B\) is balanced. Its first moments vanish and
the \(g\)-moments of its three cells are \((-1/9,-1/9,2/9)\).
It follows that, at \(\lambda=11/2\),
\[
 \mathcal E_{11/2}(B)-\mathcal E_{11/2}(T)
     \ge \frac{11}{27}-\frac{9}{8\pi}-\frac1{100}>0.             \tag{RK2}
\]
The last strict inequality uses \(\pi>3\).

## The simplex interface gap

Orient the three interface rays cyclically, let \(e_\ell\) be their
unit directions and \(d_\ell\) their label incidence vectors, and
write \(d\nu(t)=(2\pi)^{-1}e^{-t^2/2}\,dt\) for \(t>0\).
Put \(I_0=\int d\nu=1/(2\sqrt{2\pi})\).
The linear first-moment score has normal slope \(a_1=3I_0\)
on each ray. For real speeds \(f_\ell\in L^2(\nu)\), the
local quadratic cost is
\[
                A_1(f)=3I_0\sum_\ell\int f_\ell^2\,d\nu.
\]
For \(z_\ell=\int t f_\ell\,d\nu\), the degree-one interaction
is the squared norm of \(\sum_\ell d_\ell\otimes e_\ell z_\ell\).
The Gram matrix has diagonal \(2\) and off-diagonal \(1/2\),
because the two dot products for different rays are \(-1\)
and \(-1/2\). The radial interaction is the squared norm of
\(\sum_\ell d_\ell\int g(te_\ell)f_\ell\,d\nu\), whose Gram
matrix has diagonal \(2\) and off-diagonal \(-1\).

Write \(u=(f_1+f_2+f_3)/3\) and \(v_\ell=f_\ell-u\).
Balanced normal flux means that the three integrals
\(\int f_\ell\,d\nu\) agree, so \(\int v_\ell\,d\nu=0\).
The closest-regular-moment rotation gauge is
\(\sum_\ell\int t f_\ell\,d\nu=0\). Indeed differentiating the
regular moment matrix under a rotation gives its infinitesimal
rotation, and the rotated cyclic moment difference pairs with
each unit ray \(e_\ell\) in the same nonzero scalar. Thus the
first-moment interaction of \(u\) vanishes in the gauge, and
its radial interaction vanishes by \(d_1+d_2+d_3=0\).
For the standard component the two interaction eigenvalues
are respectively \(3/2\) and \(3\).

Under \(\nu/I_0\) let \(R\) denote the half-normal radius and
\(p=\mathbb P(R\ge\sqrt{2\log3})\). Center \(t\) and \(g\) in
the integrals against \(v_\ell\), and apply Cauchy–Schwarz.
The total degree-one plus weighted radial interaction satisfies
\[
 B_1(f)+\lambda B_g(f)
       \le \left(\frac{1-2/\pi}{2}+\lambda p(1-p)\right) A_1(f).
                                                               \tag{RK3}
\]
The factor is bounded away from one for \(0\le\lambda\le 11/2\).
For an elementary bound on \(p\), set
\(M(x)=e^{x^2/2}\int_x^\infty e^{-y^2/2}\,dy\).
The function \(U(x)=(x^2+2)/(x^3+3x)\) satisfies
\(U'-xU+1=-6/(x^2(x^2+3)^2)<0\) for \(x>0\).
Since \(M'-xM+1=0\), integration with the factor
\(e^{-x^2/2}\) gives \(M<U\). The function \(U\) decreases
and \(\sqrt{2\log3}>7/5\). Therefore
\[
 p<\frac4{15}U(7/5)<\frac4{25},
\]
using \(\sqrt{2/\pi}<4/5\) and
\(U(7/5)=495/868<3/5\). In (RK3) the factor is at most
\(2/11+(11/2)(84/625)<1\) by \(\pi<22/7\).

## From the interface gap to measurable local growth

For a partition \(A\) align \(T\) by its closest regular moment
matrix. Near the simplex orbit this choice also aligns the
partitions in label distance, by moment continuity and
the full rank of the regular matrix. Write
\(h_{ij}=1_{T_i\cap A_j}\) for \(i\ne j\). For the pair \(i,j\)
use orthogonal coordinates \((z,t)\) with \(t\ge0\) toward
\(T_i\) and \(z\) in the full pair line. The first-moment
linear loss is
\[
          L_1=a_1\sum_{i\ne j}\int t h_{ij}\,d\gamma_2.
\]
The exact energy loss for \(W_1+\lambda\mathcal W_g\) is \(2L_1\)
minus its degree-one and radial quadratic transfer interactions.
Here the radial linear term is zero.

Suppose a uniform positive lower multiple of \(L_1\) for this
loss fails arbitrarily close to the simplex orbit. For such
a sequence, put \(s^2=L_1\). Moment convergence gives \(s\to0\).
For each ordered pair, the projected densities
\[
             \theta_{ij}(z)
          =\frac1{s\phi_1(0)}\int_0^\infty h_{ij}(z,t)
                                                  \phi_1(t)\,dt
\]
are bounded in \(L^2(\nu_{\rm line})\), where
\(d\nu_{\rm line}=\phi_1(0)\,d\gamma_1(z)\).
Indeed if \(a(z)=\int h_{ij}(z,t)\phi_1(t)\,dt\), then
\(\int t h_{ij}(z,t)\phi_1(t)\,dt\ge a(z)^2/(2\phi_1(0))\).
This is the one-dimensional minimization of a first moment
with density bounded by \(\phi_1(0)\). It also implies
that the total Gaussian transfer mass is \(O(s)\), so
\(L_1\) controls the square of the label distance.

Pass to subsequential weak limits of all the projected
densities. The normalized normal first moment is \(O(s)\),
so the limit transfer measure is supported on the pair line.
Its support is the actual ray because \(h_{ij}\) is supported
in \(T_i\) and that cone meets the line on the ray. The
weak \(L^2\) bound makes the vertex a null set of the limit.
Tangential first moments are tight uniformly by Cauchy–Schwarz
and the Gaussian second moment on the line. Normal first
moments vanish by the normalized bound. The first moments
therefore converge to those of the limiting ray fluxes.
Integration of the bounded radial feature \(g\) also passes
to the limit. Its only discontinuity is on one centered
circle, and each limiting projected density gives zero mass
to the intersection with that circle. Equivalently one can
cut off a neighborhood of the finitely many intersections
and use the uniform \(L^2\) bounds there, with the normalized
mass outside any fixed normal tube tending to zero.

The signed difference of the two directed limits for a pair
gives the speed \(f_\ell\) on that ray. Balance and the moment
gauge pass to these fluxes. The normal first moment inequality
and weak lower semicontinuity give
\[
              1\ge \frac{a_1}2\sum_{i\ne j}\int\theta_{ij}^2\,d\nu
                \ge \frac12 A_1(f).
\]
For the last inequality use nonnegativity of each directed
limit. Normalized quadratic interactions converge by the
moment convergence just established. After also taking a
subsequence of \(\lambda\) to its limit, (RK3) bounds that
interaction strictly below \(2\), contrary to failure of a
positive multiple of \(L_1\). This proves uniform quadratic
growth for \(W_1+\lambda\mathcal W_g\) near the simplex orbit.
The positive term \(S_{r_0}/100\) can only increase the loss
relative to \(T\), so (RK1) has the same local growth.
[Normal-flux comparison for the original Mehler objective](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)

## An interior coexistence parameter

Balanced fractional assignments form a weakly compact convex
set in \(L^2(\gamma_2;\mathbb R^3)\). The quadratic objective
in (RK1) is weakly continuous by operator compactness.
Strict positivity of the quadratic form makes the objective
strictly convex. A maximizing fractional assignment is
ordinary: if two label coordinates are bounded away from zero
on a set of positive mass, split a subset into two equal-mass
pieces and make opposite small exchanges between those
coordinates on the two pieces. The original assignment is
the midpoint of distinct balanced feasible assignments, and
cannot maximize. The same reasoning gives strong compactness
for ordinary maximizing sequences with converging parameters:
every weak limit is maximizing and ordinary, and all ordinary
assignments have the same squared \(L^2\) norm.

The value minus \(\mathcal E_\lambda(T)\) is continuous and
nondecreasing in \(\lambda\), since \(\mathcal W_g(T)=0\)
and all \(\mathcal W_g(A)\ge0\). It vanishes on an initial
positive interval. Otherwise maximizers beating \(T\) at
parameters tending to zero have a strong limit that is a
simplex by uniqueness at zero. Uniform local growth then
forbids the strict improvement near that orbit. The value
difference is positive at \(11/2\) by (RK2). Its last zero
is thus interior. Take ordinary maximizers at parameters
decreasing to that zero from above. Uniform local growth
keeps them a fixed positive label distance from the simplex
orbit. Strong compactness gives an ordinary nonsimplex limit
attaining the same value as \(T\) at the coexistence parameter.

## A strictly positive reversible Markov kernel

The construction can also be made with a pointwise strictly positive
kernel and constant row sums. For every sufficiently small fixed
\(r>0\), substitute the objective
\[
      \mathcal F_{r,\lambda}(A)=r^{-1}S_r(A)
                                     +\lambda\mathcal W_g(A)
                        \qquad (0\le\lambda\le11/2).             \tag{RK4}
\]
Adding the constant kernel \(2\) changes the value by the same
constant for all balanced partitions. Its row sum is then
\(r^{-1}+2\), since \(\int g\,d\gamma_2=0\). Dividing by this
row sum gives a reversible Markov operator. The kernel is pointwise
strictly positive: the negative part of
\(\lambda g(x)g(y)\) is at most \(11/9\) and the Mehler
kernel is positive. The addition and normalization
preserve compactness and strict spectral positivity.

Here the same uniform local argument uses the normal-flux
limit for the small-correlation Mehler kernel. If local
growth fails uniformly as \(r_j\to0\) and the orbit distance
tends to zero, align the moments and put
\(s_j^2=L_{r_j}(A_j;T)/r_j\). The simplex pair contrast
vanishes on its full pair line and its spatial gradient is
bounded by \(C r_j\) by Gaussian convolution. Thus
\(s_j^2\le C\int_{\{A_j\ne T\}}|x|\,d\gamma_2\to0\)
by Gaussian integrability. The weighted transfer bound
(C8) gives uniform tightness of the normalized transfer
measures including their first moments.
The ordered flux compactness (C9)–(C10) gives the
balanced gauged ray speeds and, by (C11), \(1\ge A_1(f)/2\).
The slope in this limit is \(a_1\) as in (RK3). By (C12),
the normalized centered Mehler interaction converges to
\(B_1(f)\). These estimates involve the reference loss and
the transfer measures, so they apply unchanged to the
sequence for (RK4). In particular adding the bounded feature
does not change the localization or normalization.
Its only discontinuity is the same centered circle as
above, which is null for the limiting directed ray
measures. Hence
\(\mathcal W_g(A_j)/s_j^2\to B_g(f)\), since the
reference \(g\)-moments are zero. The exact deficit in
(RK4), divided by \(s_j^2\), tends to
\(2-B_1(f)-\lambda B_g(f)\) along a subsequence in
\(\lambda\). By (RK3) this is bounded below by a positive
constant. This proves local growth uniformly for all small
positive \(r\) and all the stated \(\lambda\). The global
bound on normalized transfer mass gives a quadratic
lower bound in Gaussian label distance.
[The fixed-cell small-correlation normal-flux limit, C8–C12 and A8](SMALL_POSITIVE_CORRELATION_INTERVAL.md#normal-flux-contradiction-and-dimension-reduction)

At \(\lambda=0\), exact simplex uniqueness holds for all
sufficiently small \(r\). On the other hand the Hermite
expansion and balance give
\[
 r^{-1}(S_r(T)-S_r(B))\le W_1(T)+\tfrac23 r .
\]
The radial competitor \(B\) still beats \(T\) at \(11/2\)
for small \(r\), by the same strict inequality as in
(RK2). Keeping one such \(r\) fixed, the compactness
and coexistence argument above gives an interior
tilt with both a simplex and a distant ordinary
nonsimplex maximizer for this strictly positive
reversible Markov kernel.

The kernel can also be chosen real analytic. Replace the step
in \(g\) by
\([1+\exp(-\epsilon^{-1}(q-\log3))]^{-1}\), then subtract its
Gaussian mean. These features are uniformly bounded and
converge almost everywhere to \(g\) as \(\epsilon\downarrow0\).
The radial competitor's three features and the half-normal
variance in (RK3) converge by dominated convergence, so the
two strict margins persist for a fixed sufficiently small
\(\epsilon>0\). The normal-flux proof then uses a continuous
feature. A constant kernel \(6\) instead of \(2\) suffices
for strict pointwise positivity and row normalization, since
the centered feature has absolute value at most one.
