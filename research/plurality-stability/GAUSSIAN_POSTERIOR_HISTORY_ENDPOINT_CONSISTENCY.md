# Brownian posterior history determines the Gaussian endpoint

For a fixed ordinary balanced Gaussian partition, a finite
history of posterior label probabilities still leaves a
strict gap in a supermartingale bound at the regular simplex.
The same holds with finitely many heat-time derivatives at
each observation. In contrast, observing the posterior along
one Brownian path at geometrically spaced times accumulating
at an interior time determines the endpoint indicators. The
uniqueness holds simultaneously for all bounded endpoint
functions on almost every path.

## Finite posterior history leaves a strict gap

Fix \(k\ge3\) and \(0<\rho<1\). Let \(B,C\) be standard
Brownian motions in the simplex active space with cross
variation \(\rho I\,dt\), and let
\[
 p_A(s,z)=P_{1-s}1_A(z),\qquad
 V_A(t,x,y)=
 \mathbb P\bigl(A(B_1)=A(C_1)\mid B_t=x,\ C_t=y\bigr).
\]
The notation \(1_A\) denotes the vector of cell indicators.
At a deterministic time \(0<t<1\), allow the bound's state
to retain the locations and posterior vectors of both paths
at any finite list of deterministic observation times up to
\(t\), including \(t\). It may also retain finitely many
heat-time derivatives of the posterior vectors at those
locations. Assume evaluation at each fixed input partition
is continuous in the observed locations, and that along
each fixed balanced partition the resulting process is an
integrable supermartingale whose terminal value dominates
the agreement indicator. Then its initial value on the
centered simplex \(T\) is strictly larger than \(S_\rho(T)\).
These are the supermartingale hypotheses used in (BL1).
[The local state and conditional agreement formula](GAUSSIAN_BELLMAN_LOCAL_STATE_RIGIDITY.md#supermartingale-bounds-and-conditional-agreement)

At time \(t\), conditional agreement given the full path
history is still \(V_A(t,B_t,C_t)\), by independence of
the future Brownian increments. Thus the same conditional
expectation argument as in (BL1) bounds \(V_A(t,x,y)\)
from above by the evaluated state, now including the
past locations. This inequality holds pointwise in the
locations by continuity and positive joint Gaussian density
of the observations at distinct positive times. Deterministic
time-zero positions and repeated observations are counted
once.

To obtain strictness, use the compact normal variation
on one open simplex facet from (BL3), with \(v\) an
interior direction of the facet. Each added probability
constraint has first-variation weight given by a Gaussian
density of covariance \((1-s)I\), where \(s\le t\).
Along \(z=\ell v\), the weights for its heat-time derivatives
are the same density multiplied by polynomials in \(\ell\)
of distinct even degrees. The global volume weight is the
standard Gaussian density. Choose observation locations
with distinct projections on \(v\) at each positive time,
and with the two current observations satisfying (BL2).
These are open conditions attained at tuples in the
support of the joint law.

The weights are linearly independent on the facet. A
linear relation restricted to the ray extends along
its full supporting line by analyticity. At \(\ell\to
-\infty\), first divide by the Gaussian having largest
variance, and among those the smallest projected mean.
The remaining Gaussian ratios tend to zero faster than
any polynomial. The polynomial for that Gaussian must
vanish, and repetition eliminates the others. The
distinct degrees within each density group eliminate
its coefficients.

The conditional agreement derivative \(g\) from (BL3)
lies outside this span. At \(\ell\to-\infty\), the
observation ray points into the region of the spectator
labels, away from both facet labels. The two conditional
pair contrasts in \(g\) therefore have an extra Gaussian
decay factor. Each of its two marginal density factors
has variance \(1-t\), no larger than the variance of
any probe density. Thus \(g\) divided by any probe density
tends to zero faster than any polynomial on the negative
ray. The same elimination forces all coefficients of a
putative span identity to be zero. The positive-ray
limit (BL4) excludes \(g\equiv0\).
[The conditional facet derivative and the two ray limits, BL2–BL4](GAUSSIAN_BELLMAN_LOCAL_STATE_RIGIDITY.md#a-variation-with-three-gaussian-volume-constraints)

Choose a compact smooth facet speed annihilating all the
probe derivatives and the full volume derivative, with
strictly positive conditional agreement derivative.
Linear independence supplies compact correction speeds
with an invertible probe Jacobian, as in (BL5). The
implicit function theorem preserves the full volumes
and every observed posterior coordinate exactly. A fixed
small positive amplitude improves \(V_A(t,x,y)\) by a
positive uniform amount on a neighborhood of the chosen
tuple. Each corrected partition is ordinary, and the
same partition is used at both terminal observations.
The pointwise bound must therefore exceed \(V_T\) by
that amount on the neighborhood. Along \(T\), this
neighborhood has positive probability. The supermartingale
bound from time zero to \(t\) gives the strict initial
gap.

With finitely many spatial derivatives instead, the
probe first-variation weights still lie in the same
polynomial-density span. Vanishing to first order and
exact repair of the full volumes give the obstruction
to local Lipschitz sharpness, as in (BL7). This assertion
does not require preservation of spatial derivatives
to all orders in the perturbation parameter.
[The local Lipschitz obstruction and the exact first-jet constraint](GAUSSIAN_BELLMAN_LOCAL_STATE_RIGIDITY.md#finite-spatial-derivatives-and-local-lipschitz-bounds)

## Geometric observations identify the endpoint

Fix \(0<t_0<1\), \(0<\eta<t_0\), and \(0<\lambda<1\).
Set \(t_n=t_0-\eta\lambda^n\) and
\[
 G_n=\frac{B_{t_n}-B_{t_0}}{\sqrt{t_0-t_n}} .
\]
The \(G_n\) form a stationary standard Gaussian sequence
with cross-covariance \(\lambda^{|n-m|/2}I\). Almost surely,
they visit every open ball with rational center and
radius infinitely often. Indeed for such a ball \(E\),
the covariance of its two indicators is at most
\(\lambda^{|n-m|/2}\operatorname{Var}(1_E(G_0))\) by
the largest positive-degree Ornstein–Uhlenbeck multiplier.
The variance of their empirical mean is \(O(N^{-1})\).
Chebyshev's inequality and Borel–Cantelli along \(N=2^j\)
give convergence to \(\gamma(E)>0\) along that subsequence.
Take the countable intersection over balls. This is a
property of the path alone.

Fix such a path and bounded measurable endpoint
functions \(f,g\). Suppose their posterior values agree
at every sample \((t_n,B_{t_n})\). Write
\(h(t,x)=P_{1-t}(f-g)(x)\) and \(x_0=B_{t_0}\).
Gaussian smoothing makes \(h\) smooth in \(t,x\) for
\(t<1\) and real analytic in \(x\) at \(t=t_0\).
Suppose there is a nonzero mixed derivative at
\((t_0,x_0)\), and let \(m\) be the least weight
\(2j+|\beta|\) of such a derivative
\(\partial_t^j\partial_x^\beta h\).
Writing \(\delta_n=t_0-t_n\), Taylor expansion at the
accumulation point gives, on every subsequence with
bounded \(G_n\),
\[
 \frac{h(t_n,B_{t_n})}{\delta_n^{m/2}}
   =\sum_{2j+|\beta|=m}
       \frac{(-1)^j\partial_t^j\partial_x^\beta h(t_0,x_0)}
            {j!\,\beta!}\,G_n^\beta+o(1).                 \tag{PH1}
\]
An ordinary Taylor expansion through total order \(m\)
suffices, since \(|t_n-t_0|=O(\delta_n)\) and
\(|B_{t_n}-x_0|=O(\sqrt{\delta_n})\) on the subsequence.
In the displayed polynomial, each spatial monomial
has only one possible time index \(j\). By definition
of \(m\) the polynomial is nonzero. Every point of
space is a limit along a bounded subsequence of the
\(G_n\), whereas the left side of (PH1) vanishes.
The polynomial would be zero everywhere, a contradiction.

All spatial derivatives at the accumulation point
therefore vanish. Real analyticity makes \(h(t_0,\cdot)\)
identically zero, and Gaussian convolution is injective
on bounded functions, so \(f=g\) almost everywhere.
The argument used only the path event above and
properties shared by all bounded \(f,g\), which gives
the simultaneous assertion. It applies componentwise
to fractional as well as ordinary label indicators.

## Finite prefixes give qualitative endpoint control

Fix a path in the uniqueness event and an ordinary
reference partition \(T\) of prescribed Gaussian
volumes. For every \(\delta>0\) there exist a finite
prefix length \(N\) and a tolerance \(\tau>0\) such
that every fractional assignment \(f=(f_i)\) with
these volumes and
\[
 \max_{n\le N,\ i}
   \left|P_{1-t_n}(f_i-1_{T_i})(B_{t_n})\right|\le\tau
\]
satisfies \(1-\sum_i\int_{T_i}f_i\,d\gamma<\delta\).
The simplex-valued assignments of prescribed volumes
are weak-* compact. Each probe is a continuous linear
functional in this topology, because its Gaussian
density is integrable against Lebesgue measure and
has an \(L^1(\gamma)\) density relative to \(\gamma\).
The displayed disagreement quantity is also continuous.
If no finite prefix and tolerance worked, the closed
sets of assignments with disagreement at least
\(\delta\) and first \(N\) probe discrepancies at most
\(1/N\) would be nonempty and nested. Their intersection
would contradict the geometric-observation uniqueness.

For ordinary assignments the disagreement quantity
is the Gaussian mass of the changed labels. In
particular, exact compatibility with increasingly
long prefixes of a fixed simplex history makes every
eligible partition converge to \(T\) in Gaussian
label distance, uniformly over the eligible family.
For a fixed later observation, the difference of
conditional agreements then tends to zero by the
union bound under its two Gaussian marginals. The
finite-history obstruction gives a strict gap at each
finite stage and the compactness result gives no
positive bound independent of the prefix length.
