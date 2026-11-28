# Reflection tests for a centered three-cell height graph

Fix \(0<r<1\) and an actual balanced planar three-cell maximizer
with a separating axis \(x=0\) and a finite continuous graph
\(y=h(x)\). Thus cell 3 is above the graph and cells 1 and 2
are below it on the right and left, respectively. The
[affine separation argument](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md#finite-heights-and-centered-separation)
obtains this geometry from an affine separating line for
\(r\ge2/3\). The following reflection tests use the
geometry and stationarity at the fixed correlation.

Write \(h_+(s)=h(s)\) and \(h_-(s)=h(-s)\) for \(s\ge0\), and
\(h_0=h(0)\). The two profiles agree if either one is monotone.
They also agree if the inverse images of some closed
neighborhood of \(h_0\) are compact for both profiles and
the first nonzero height term at each endpoint has odd order.
Any asymmetry requires infinitely many strict sign changes
of \(h_+-h_-\) on the open half-line.

## The common height distribution

Put \(\sigma=\sqrt{1-r^2}\) and let \(p_\sigma\) be the
centered Gaussian density of standard deviation \(\sigma\).
Stationarity on the \(1/2\) ray and analyticity in the height
give \(V_1(0,y)=V_2(0,y)\) at every \(y\). At high \(y\)
both unadjusted scores tend to zero, so their prices agree.
Vertical deconvolution as in (ASR3) shows that \(h_+\) and
\(h_-\) have the same distribution under \(p_\sigma(s)\,ds\)
on \(s>0\). Each of \(\Phi(h_+)\) and \(\Phi(h_-)\) also has
the same weighted integral \(1/3\) under the standard
Gaussian measure on that half-line, by label balance.

Suppose one of the profiles is nondecreasing. Its superlevel
sets are upper intervals in \(s\). Among sets of a prescribed
\(p_\sigma\)-mass, an upper interval uniquely maximizes
integration against \(p_1\). Indeed the likelihood ratio
\[
              p_1(s)/p_\sigma(s)
            =\sigma\exp\!\left(\frac{r^2s^2}{2\sigma^2}\right)
\]
is strictly increasing. Apply this to the superlevel sets
of \(\Phi(h_\pm)\), which have matching \(p_\sigma\)-mass,
and use layer-cake integration. Equality of the two total
weighted integrals forces equality for almost every
superlevel set. The profiles agree almost everywhere
and then everywhere by continuity. If one profile is
nonincreasing, use the minimizing lower intervals in the
same argument.

## The density germ at the junction level

Each profile is real analytic on \(s>0\) and has an analytic
extension at zero. On its side the stationary score
\(V_3-V_i\) vanishes on the graph and has strictly positive
vertical derivative, so the implicit function theorem applies
up to the common junction. Suppose now that both profiles
have nonzero derivatives at zero and their preimages of
\([h_0-\epsilon,h_0+\epsilon]\) are compact for some
\(\epsilon>0\).

For either profile, the preimages of \(h_0\) are finite.
The derivative at zero precludes an identically constant
profile. Each interior preimage can be put in a local
real-analytic coordinate \(u\) increasing with \(s\), so
that the height is \(h_0\pm u^m\) for some integer \(m\ge1\).
If \(\chi(u)\) is the analytic local primitive of the
measure \(p_\sigma(s)\,ds\) in this coordinate, the
contribution of that neighborhood to the cumulative
height distribution, apart from a constant, has the
following form for a height offset \(v\) near zero.

- For even \(m\), there is a contribution on just one
  side of \(v=0\), involving
  \(\chi(|v|^{1/m})-\chi(-|v|^{1/m})\) up to an
  overall sign. The odd powers of \(|v|^{1/m}\) in this
  expression are all noninteger powers of \(|v|\).
- For odd \(m\), the contribution involves
  \(\chi(\operatorname{sgn}(v)|v|^{1/m})\), up to
  sign changes in the argument and the value according
  to the orientation. Every integer power in its
  one-sided expansions is part of the same analytic
  germ on the two sides, since for \(k=mq\) the factor
  \(\operatorname{sgn}(v)^k |v|^q\) is \(v^q\).

The rest of the cumulative distribution is constant
near \(v=0\). To see this, take disjoint neighborhoods of
all preimages, including the endpoint zero. On the compact
preimage band outside those neighborhoods the height stays
a positive distance from \(h_0\), and points outside that
band cannot contribute a change of mass near \(h_0\).

In the convergent one-sided Puiseux expansions of the
cumulative distribution, take the integer-power germ
on the right minus the integer-power germ on the left,
expressing both as powers of signed \(v\). The interior
contributions cancel by the two cases above. At the
simple endpoint \(s=0\), write \(s(v)\) for the analytic
local inverse of its height \(h_0+v\), extended as an
analytic function across zero, and set
\(\Psi(s)=\int_0^s p_\sigma(t)\,dt\).
If the endpoint derivative is positive, the cumulative
increment is \(\Psi(s(v))\) on the right only. If it is
negative, the cumulative increment is \(-\Psi(s(v))\)
on the left only. Thus in both cases the right-minus-left
integer germ contributed by the endpoint is
\(\Psi(s(v))\).

The two profiles have identical cumulative distributions,
so the uniqueness of their one-sided Puiseux expansions
gives equality of these endpoint germs. Since
\(\Psi'(0)>0\), their inverse height germs are equal.
Their profiles therefore agree near zero and hence on the
whole half-line by real analyticity. This argument permits
interior critical points at the same height as the junction.

At least one endpoint derivative is positive for every
centered affine-separated graph. At \((0,h_0)\), use the
notation
\[
 a=\partial_x U_1>0,\qquad b=-\partial_x U_2>0,\qquad
 c=\partial_y(V_3-V_1)=\partial_y(V_3-V_2)>0 .
\]
The first two signs follow from the strict one-sided
supports of the corresponding labels. Implicit differentiation,
using \(\sum_i U_i=1\), gives
\[
 h_+'(0)=(2a-b)/c,\qquad h_-'(0)=(2b-a)/c .
\]
In particular both slopes are positive if the two endpoint
derivatives are nonzero and the reflection criterion
applies.

Suppose the compact-preimage assumption holds and both endpoint
orders are odd. The sum of the endpoint derivatives is positive,
so at least one endpoint is simple. If the other has odd order
\(m\ge3\), examine the coefficient of \(|v|^{1/m}\) in the
two one-sided cumulative distributions and add those coefficients.
An interior preimage with even order \(n\) would contribute to
that exponent only at numerator \(n/m\), which is even if it
is an integer. Its Puiseux difference contains only odd
numerators, so there is no contribution. For an interior
preimage with odd \(n\), any integer numerator \(n/m\) is odd.
The coefficients in the two one-sided expansions then have
opposite signs and cancel. The simple endpoint has no term at
that exponent. The nonsimple endpoint, however, has a nonzero
leading coefficient at exponent \(1/m\), on just one side.
This contradicts equality of the cumulative distributions.
Therefore both endpoints are simple and the preceding proof
gives reflection.

In particular, if both tails eventually stay on one side
of the common pressure offset, the height has compact
preimages of every bounded interval by the linear escape
theorem. The odd endpoint condition then suffices for
the same reflection conclusion.
[Properness from one-sided tails](THREE_CELL_ONE_SIDED_END_ESCAPE.md#properness-of-a-centered-affine-separated-graph)

## The folded conditional scores

For \(x>0\) set \(u_\pm(x,y)=1_{\{y>h_\pm(x)\}}\) in folded
coordinates and \(d=u_+-u_-\). Let \(P_y\) denote smoothing
in the vertical variable by the one-dimensional Mehler
operator at \(r\), and put \(f_y(s)=(P_y d)(s,y)\).
For every \(y\), the sign of \(f_y(s)\) at each \(s>0\) is
the sign of \(d(s,\cdot)\) wherever that function is nonzero,
which is the sign of \(h_-(s)-h_+(s)\). The two folded
adjusted score contrasts are
\[
 G_+(x,y)=V_3(x,y)-V_1(x,y),\qquad
 G_-(x,y)=V_3(-x,y)-V_2(-x,y).
\]
Their difference satisfies
\[
 D(x,y)=G_+(x,y)-G_-(x,y)
     =\int_0^\infty
        [2p_\sigma(s-rx)-p_\sigma(s+rx)] f_y(s)\,ds .
                                                               \tag{AFR1}
\]
The coefficients come from the source contrasts
\(2u_\pm-1\) on the same side of the observation and
\(u_\mp\) on its other side. The equal prices cancel.
At \(x=0\), \(D(0,y)=0\). At every actual disagreement
\(d(x,y)\ne0\), the sign of \(D(x,y)\) equals that of
\(d(x,y)\), because the respective two folded scores
have opposite strict signs between their height roots.
The common height law also gives
\[
                    \int_0^\infty p_\sigma(s)d(s,y)\,ds=0
                                                               \tag{AFR2}
\]
for almost every \(y\) and ensures both signs at almost
every height with a positive-measure disagreement.

The folded energy yields another necessary condition.
Let \(M\) and \(M_R\) denote the full Mehler kernel
between two points in the positive horizontal half-plane,
with the second point unchanged and reflected, respectively.
Keep \(u_++u_-\) fixed and vary their allocations on the
set where that sum is one. The sum of the two half-plane
volumes is then fixed, and their difference stays fixed
exactly when \(\int d\,d\gamma_2\) is preserved over the
folded half-plane. The portion of the noise stability
that depends on the allocation equals
\(\frac12\langle d,(2M-M_R)d\rangle\).
In fact the same-side terms are
\(2\sum_\pm\langle u_\pm,M u_\pm\rangle\) plus
terms linear in their sum, and across the sides only
common label 3 contributes, giving
\(2\langle u_+,M_R u_-\rangle\). Thus the actual
allocation maximizes the stated positive-kernel quadratic
form subject to the prescribed signed mass.

## The number of asymmetric crossings

Remove the positive Gaussian factors in the kernel in
(AFR1). With \(\lambda=r/\sigma^2>0\), the resulting
kernel on \(x\ge0,\ s>0\) is
\[
 L(x,s)=2e^{\lambda xs}-e^{-\lambda xs}
       =\sum_{n\ge0}\frac{(2-(-1)^n)\lambda^n}{n!}x^ns^n .
\]
It is strictly totally positive for increasing distinct
nodes in each variable. For finite truncations,
Cauchy–Binet expands each ordered minor as a sum of
nonnegative products of generalized Vandermonde
determinants. Each such determinant on positive
increasing nodes is positive. One can see its
nonvanishing by dividing a generalized polynomial
by its lowest power and repeatedly applying Rolle's
theorem, and its sign by coalescing the nodes at a
positive point and taking the Wronskian. Nodes with
first coordinate zero follow by continuity, and the
ordinary Vandermonde term with consecutive powers
starting at zero gives a strictly positive contribution
even then. Passing to the infinite series proves the
claim.

If \(h_--h_+\) has a finite number \(N\) of sign changes
and is not identically zero, \(D(\cdot,y)\) has at most
\(N\) distinct zeros on \([0,\infty)\) for each \(y\).
For otherwise choose \(N+1\) zero rows of \(L\).
Use its determinant on those rows with \(N\) fixed
columns at the sign-change locations of \(f_y\) and
one variable column \(s\). Strict total positivity shows
that the determinant has exactly the prescribed
alternating signs between the fixed locations, up to
one global sign. Its integral against \(f_y(s)\) and
the positive Gaussian factor is therefore nonzero.
Expansion along the variable column instead makes the
integral zero by the chosen zero rows and (AFR1).
This contradiction gives the bound. It includes
the zero at \(x=0\).

For \(N=0\), the common height law already forces equality.
For \(N=1\), choose a height with disagreement in the
first sign interval and with (AFR2). At this height
there is disagreement of the opposite sign in the
second interval. By (AFR1), the score difference has
opposite signs at those two points. Its positive zero
between them and its zero at the axis violate the bound.

For \(N=2\), again choose a height with disagreement in
the first interval and with (AFR2). The opposite sign
occurs in the middle interval. For large \(x\), the
score difference at this height has the sign of the
last interval, which is also the sign of the first
interval. Indeed \(f_y\) is nonzero on a compact
subinterval of the last sign interval, strictly to the
right of all sources with opposite sign. In the
transform with kernel \(L\), the contribution of that
subinterval eventually dominates all opposite
contributions by exponential tilting in \(s\).
There are therefore two positive zeros as well as the
zero at the axis, again impossible. This excludes
every asymmetric graph with at most two sign changes.

## Eventual ordering forces reflection

Suppose the two half-profiles are eventually ordered. Then
they agree. The proof first uses the common height law to
reduce any possible ordered asymmetric graph to bounded
profiles, then compares their stationary equations on the tail.

First suppose one profile is unbounded above. Take \(Y\) above
the two profiles on the compact part preceding the eventual order.
At any height \(y>Y\), the difference of the two sublevel
indicators has one weak sign everywhere. Its integral with
respect to \(p_\sigma(s)\,ds\) is zero by the common height law.
Thus the sublevel indicators agree almost everywhere for every
rational \(y>Y\). Using these levels, the profiles agree on
the open set where either is above \(Y\), first almost
everywhere and then everywhere by continuity. This set is
nonempty, so analyticity gives equality throughout.
Unboundedness below is handled by negating both heights
in this argument.

It remains to treat two bounded profiles. Both half-profiles
give bounded isolated ends, and they converge
to the same pressure height \(b\). Their outer pair pressures
are equal because the two inner labels have the same price.
[Bounded isolated ends and the pressure offset](THREE_CELL_ISOLATED_END_RIGIDITY.md#bounded-ends-converge-to-the-pressure-offset)

Assume for definiteness that \(\delta(x)=h_+(x)-h_-(x)\ge0\)
on a tail \(x>S\). If \(\delta\) vanishes identically there,
analyticity proves equality everywhere. Otherwise choose a
compact interval \(I=[a,a']\) with \(a>S\) on which
\(\delta\) has a positive lower bound. Define
\[
   J(y,z)=\Phi((r y-z)/\sigma),\qquad
   k(y,s)=\int_0^1
       p_\sigma(r y-h_-(s)-\theta\delta(s))\,d\theta .
\]
The folded score identity (AFR1), evaluated at \(y=h_+(x)\),
and the two stationary root conditions give
\[
\begin{split}
 A(x)\delta(x)&=\int_0^\infty
       [\,2p_\sigma(s-rx)-p_\sigma(s+rx)\,]\,
              k(h_+(x),s)\delta(s)\,ds,\\
 A(x)&=\int_0^1
      \partial_yG_-(x,h_-(x)+\theta\delta(x))\,d\theta .
                                                               \tag{AFR3}
\end{split}
\]
Here \(J(y,h_-(s))-J(y,h_+(s))=k(y,s)\delta(s)\).
The left vertical derivative comes from the folded conditional
representation of \(G_-\). On bounded \(y\) tending to \(b\),
its same-side term tends to \(2r p_\sigma((1-r)b)\) and
its opposite-side term tends to zero. The first limit follows
by concentrating \(p_\sigma(s-rx)\) at \(s\to+\infty\) and
using \(h_-(s)\to b\). Thus \(A(x)\to2r p_\sigma((1-r)b)\).
Also \(k(h_+(x),s)\to p_\sigma((1-r)b)\) uniformly when
both \(x\) and \(s\) tend to infinity.

The possibly negative contribution in (AFR3) is confined
to \(s\le S\), where its magnitude is bounded by
\(C\exp(-(rx-S)^2/(2\sigma^2))\) for large \(x\).
The fixed interval \(I\) contributes at least
\(c\exp(-(rx-a)^2/(2\sigma^2))\). In fact all heights
are bounded, and \(2p_\sigma(s-rx)-p_\sigma(s+rx)\ge
p_\sigma(s-rx)\) for \(s,x\ge0\). Since \(a>S\), this
positive contribution absorbs the earlier one with a strict
remainder for large \(x\). In particular \(\delta(x)>0\)
there.

Choose \(0<\ell<r<u<1\). The interval \(I\) is disjoint from
\([\ell x,ux]\) for large \(x\), so it can be used to absorb
the negative error before retaining only the integral on
\([\ell x,ux]\). On the retained interval the ratio
\(p_\sigma(s+rx)/p_\sigma(s-rx)\) is exponentially small
in \(x^2\). The coefficient limits in (AFR3), and
\(r<1\), therefore give
\[
        \delta(x)\ge\gamma\int_{\ell x}^{ux}
                     \delta(s)p_\sigma(s-rx)\,ds
                             \qquad(x\ \text{large})
                                                               \tag{AFR4}
\]
for some \(\gamma>1\). The integral kernel on the retained
interval has mass tending to one. Choose \(R\) so large that
this mass times \(\gamma\) is at least \(\gamma_0>1\)
for \(x>R\), and \(\delta\) is strictly positive throughout
\([\ell R,R]\). Iterate (AFR4) using the normalized kernel
until first entry into that compact interval. Every retained
step multiplies the coordinate by a factor between \(\ell\)
and \(u\), so the stopping time is finite and is at least
\(\log(x/R)/\log(1/\ell)\). The positive minimum of \(\delta\)
on \([\ell R,R]\) then gives a lower bound growing at least
as a positive power of \(x\). This contradicts boundedness.

Exchanging the two profiles proves the result for the opposite
eventual order. In particular, an asymmetric maximizing graph
must have both signs of the discrepancy arbitrarily far out
and hence infinitely many strict sign changes.

## The common height law and balance can both hold without reflection

The common height law alone, even together with the exact
Gaussian masses and simple positive endpoint slopes, allows
asymmetry when heights near the junction level occur arbitrarily
far out. For any \(0<r<1\), put
\[
 q(s)=2\Phi(s/\sigma)-1,\qquad
 v(q)=q(1-q)-\mu,\qquad
 \mu=2\int_0^\infty q(s)(1-q(s))p_1(s)\,ds .
\]
The variable \(q\) is uniform on \((0,1)\) for the normalized
measure \(2p_\sigma(s)\,ds\). Choose a nonzero linear combination
\(g\) of \(\sin(2\pi q)\) and \(\sin(4\pi q)\) satisfying
\[
                       \int_0^\infty g(q(s))p_1(s)\,ds=0 .
\]
This is one linear constraint on a two-dimensional space.
The function \(g\) is antisymmetric under \(q\mapsto1-q\)
and vanishes at both endpoints. Choose \(a>0\) so that
\(a|g'(0)|<1/2\), and choose \(\epsilon>0\) sufficiently
small that
\[
 F_\pm(s)=\frac23+
              \epsilon\,[v(q(s))\pm a g(q(s))]
\]
take values in \((0,1)\). Set \(h_\pm=\Phi^{-1}(F_\pm)\).
Reflection of \(q\) exchanges the two expressions, so the
profiles have the same \(p_\sigma\)-height distribution.
Also \(\int_0^\infty F_\pm p_1=1/3\), by the choice of
\(\mu\) and \(g\). Thus the epigraph and its two complementary
half-plane pieces are balanced.

The profiles are analytic, have a common endpoint height,
and have positive derivatives there, since \(v'(0)=1\) and
\(q'(0)>0\). They differ because \(g\) is nonzero. Both
approach their common endpoint height as \(s\to\infty\),
so the compact-preimage hypothesis in the density-germ
argument fails. These profiles satisfy the height law and
balance. Their construction imposes no stationary score
equations along the two graph interfaces.
