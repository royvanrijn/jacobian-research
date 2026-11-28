# Sharp pressure-gradient test for bounded isolated ends

The [shifted-contour straightening
theorem](THREE_CELL_TAIL_FLATNESS.md#improving-the-gaussian-rate-by-shifting-the-inversion-contour)
makes the straight alternative unconditional for the bounded
isolated graphs considered here. The conditional gradient identity
and the exterior geometry below also identify the straight-end
decay rate and the moment constraints on a compact discrepancy.

Fix \(0<r<1\) and a balanced planar three-cell Gaussian noise-stability
maximizer. Suppose an interface is a bounded graph \(z=h(t)\) on a full
forward half-plane \(t>S\), with one label on each side and the third
label absent there. Orient \(f\) as the indicator of the upper label
minus that of the lower label, and write \(F=T_r f\).
Set \(\sigma=\sqrt{1-r^2}\) and let \(p_\sigma\) denote the normal
density with standard deviation \(\sigma\). The interface pressure
\(\beta=F(t,h(t))\) determines \(b\) by
\[
 \beta=2\Phi((r-1)b/\sigma)-1.
\]
The [bounded-end theorem](THREE_CELL_ISOLATED_END_RIGIDITY.md#bounded-ends-converge-to-the-pressure-offset)
and [Gaussian flattening](THREE_CELL_TAIL_FLATNESS.md) give
\(h(t)=b+O(e^{-ct^2})\). Define
\[
 \delta(t)=2r p_\sigma((1-r)b)-|\nabla F(t,h(t))|.
\]
Then \(\delta(t)>0\) on a sufficiently late tail. The end is straight
if and only if, for some finite \(A,C\),
\[
  \delta(t)\le C e^{-(rt-A)^2/(2\sigma^2)}                  \tag{PG1}
\]
there. In the curved case, \(\delta(t)e^{(rt-A)^2/(2\sigma^2)}
\to+\infty\) for every finite \(A\). The test includes zero pressure.
For a straight end, the coefficient of the linear correction to the
Gaussian rate gives the maximal longitudinal extent of the source
discrepancy from the binary pressure halfspace.

At nonzero pressure, the fixed-partition correlation derivative
gives a complementary signed condition. A curved end forces
\(g(t)-g_\infty\) to have the sign of \(\beta\) and exceed an
arbitrarily shifted critical Gaussian along subsequences, where
\[
 g(t)=\left.\partial_\rho T_\rho f(t,h(t))\right|_{\rho=r},
 \qquad
 g_\infty=\left.\partial_\rho
       \bigl[2\Phi((\rho b-b)/\sqrt{1-\rho^2})-1\bigr]\right|_{\rho=r}.
\]
These assertions hold without reflection or a condition on
correlation contact. With three bounded isolated ends, one for
each pair, exact straightness of all three forces agreement with
a centered regular simplex outside a compact set.

## Conditional gradient defect

Write \(x=(t,h(t))\), and let \(Y\) have law \(N(rx,\sigma^2 I_2)\).
At a large regular interface point set \(n=\nabla F(x)/|\nabla F(x)|\),
oriented toward the upper label, and put \(d=(1-r)b\). The affine binary
contrast \(f_x^*(Y)=\operatorname{sgn}(n\cdot(Y-rx)-d)\) has the same
conditional mean \(\beta\) as \(f\). Differentiating the Gaussian
density in the observation and using this equality of means gives
\[
 \begin{split}
 \delta(t)&=\frac r{\sigma^2}\,
   \mathbb E\bigl[(n\cdot(Y-rx)-d)(f_x^*(Y)-f(Y))\bigr]\\
           &=\frac r{\sigma^2}\,
   \mathbb E\bigl[|n\cdot(Y-rx)-d|\,|f_x^*(Y)-f(Y)|\bigr].
                                                               \tag{PG2}
 \end{split}
\]
Indeed the directional derivative for the comparison contrast,
viewed as a fixed function of the source, is \(2r p_\sigma(d)\).
The sign identity uses \(-1\le f\le1\). This is the elementary
halfspace comparison behind the sharp reverse Gaussian gradient
bound. The third label has positive measure and gives zero
contrast, so the deficit is strictly positive.

Write \(u=h-b\) and \(a=h'(t)\), so \(n=(-a,1)/\sqrt{1+a^2}\).
For \(Y=(s,v)\), the comparison interface is
\[
 \ell_t(s)=r h(t)+a(s-rt)+d\sqrt{1+a^2}.
\]
On \(s>S\), the contribution inside the longitudinal integral
in (PG2) equals
\[
 \frac{2}{\sqrt{1+a^2}}
       \int_{\min(h(s),\ell_t(s))}^{\max(h(s),\ell_t(s))}
        |v-\ell_t(s)|p_\sigma(v-rh(t))\,dv .              \tag{PG3}
\]
In particular, for any fixed \(0<\lambda<r<\Lambda\), this
contribution is comparable to \((h(s)-\ell_t(s))^2\)
uniformly on \([\lambda t,\Lambda t]\) for all sufficiently large
\(t\). To verify boundedness of the arguments, the bounded second
derivative of \(h\) and Gaussian bounds on \(u\) give a Gaussian
bound on \(a\). For example, choose an increment
\(\epsilon=e^{-ct^2/4}\) and bound the difference between
\(u'(t)\) and \((u(t+\epsilon)-u(t))/\epsilon\) by
\(\|u''\|_\infty\epsilon/2\). Thus \(ta\to0\).
The same comparison in (PG3) holds on any fixed compact source
interval in \(s>S\) because \(\ell_t(s)\to b\) uniformly there.

If \(u\) is nonzero arbitrarily far out, for each \(A\) choose
\([a_0,a_1]\subset(\max(A,S),\infty)\) with \(|u|\)
bounded below by a positive constant throughout the interval.
On that interval the contribution (PG3) is bounded below for
all sufficiently large \(t\). Keeping only this contribution
in (PG2) yields
\[
    \delta(t)\ge c_{a_0} e^{-(rt-a_0)^2/(2\sigma^2)}
                  \qquad(rt>a_1,\ t\text{ sufficiently large}).
                                                               \tag{PG4}
\]
Since \(a_0>A\), this gives the asserted limit for a curved end.

If \(h=b\) on a tail, then at all late observations \(n=(0,1)\)
and \(f_x^*(s,v)=\operatorname{sgn}(v-b)\).
Formula (PG2) becomes
\[
 \begin{split}
 \delta(t)&=\frac r{\sigma^2}\int H(s)p_\sigma(s-rt)\,ds,\\
 H(s)&=\int |v-b|\,|\operatorname{sgn}(v-b)-f(s,v)|
                           p_\sigma(v-rb)\,dv .           \tag{PG5}
 \end{split}
\]
The nonnegative bounded function \(H\) is zero when the
two contrasts agree almost everywhere on the source slice.
Its essential support has a finite right endpoint \(B\)
and is nonempty because of the third label.
The Gaussian upper tail bound gives (PG1) with \(A=B\).
For every \(\epsilon>0\) the integral of \(H\) on
\([B-\epsilon,B]\) is positive, giving the lower bound
\(c_\epsilon e^{-(rt-B+\epsilon)^2/(2\sigma^2)}\).
Consequently
\[
   \lim_{t\to\infty}
      \frac{\log\delta(t)+r^2t^2/(2\sigma^2)}{t}
                  =\frac{rB}{\sigma^2}.                   \tag{PG6}
\]

## Nonzero-pressure correlation-derivative test

Assume \(\beta\ne0\) and that the end is curved.
The derivative \(g\) tends to \(g_\infty\) since the source
difference from the affine binary contrast consists of a
Gaussian-thin strip on the right and bounded data on the
left. Differentiating the conditional density inserts a
polynomial, and integrating against this strip and the
left tail preserves Gaussian decay.
[Reflection-free application of the strip calculation](THREE_CELL_SECOND_ORDER_PARAMETER_DETECTION.md#higher-order-ties-with-three-isolated-ends)

Write \(G(z)=2\Phi((rz-b)/\sigma)-1\) and \(U=F-G\).
Equations (CR1)–(CR5) in the purely two-label forward half-plane
give a scalar equation for \(u\),
\[
 u''+p(t)u'+q(t)u
      =\frac{\mathcal B_z U(t,h(t))-r(g(t)-g_\infty)}
             {-F_z(t,h(t))},
 \qquad p(t)=-t+O(1),\quad p'(t)=O(1),\quad q(t)\to1 .
                                                               \tag{PG7}
\]
Here \(\mathcal B_z\) is the operator in (CR4).
In the derivation the generator identity
\(LT_r=-r\partial_rT_r\) holds with the observation fixed.
The denominator is \(C'(u)-U_z=-F_z\) in the notation of
(CR1). It is negative on the tail. The coefficient of the
source displacement squared in \(\mathcal B_z U\) at the
pressure line is
\(-b r p_\sigma((r-1)b)/(1+r)^2\), with the sign of
\(\beta\). For any large cutoff \(S\), the remaining left
error is \(O(e^{-(rt-S)^2/(2\sigma^2)})\).
[Derivative equation and signed Taylor coefficient, CR1–CR5](THREE_CELL_REFLECTED_CONSTANT_CORRELATION_RIGIDITY.md#the-signed-quadratic-term)

Take \(S\) large enough to fix the sign of the strip
integrand. Any fixed interval on which \(u^2\) is bounded
below has a contribution with this sign and dominates
the left error if its left endpoint exceeds \(S\).
In particular for each \(A\) there exist \(a_0>\max(A,S)\)
and \(c_{a_0}>0\) with
\[
  \operatorname{sgn}(\beta)\mathcal B_z U(t,h(t))
      \ge c_{a_0}e^{-(rt-a_0)^2/(2\sigma^2)}                 \tag{PG8}
\]
on a sufficiently late tail.

By the one-sided isolated-end theorem, \(u\) has positive
and negative values arbitrarily far out. Set
\(v=e^{\frac12\int^t p(s)\,ds}u\).
Its conjugated equation has zeroth-order coefficient
\(q-p'/2-p^2/4<0\) on a tail. The conjugated left side is strictly
negative at positive maxima of \(v\) on bounded positive
components and strictly positive at negative minima on bounded
negative components. The denominator in (PG7) is negative.
Consequently \(r(g-g_\infty)-\mathcal B_zU\) takes both
signs arbitrarily far out. For \(\beta>0\), take negative
minima of \(v\) tending to infinity. For \(\beta<0\),
take positive maxima. In both cases (PG8) gives
\[
 \limsup_{t\to\infty}
    e^{(rt-A)^2/(2\sigma^2)}\operatorname{sgn}(\beta)
                   \bigl(g(t)-g_\infty\bigr)=+\infty
           \qquad\text{for every }A\in\mathbb R.           \tag{PG9}
\]
Thus the upper bound on the signed deviation alone at
the critical Gaussian rate is enough for straightness.

## Exterior geometry when all three ends straighten

Assume a bounded isolated graph end exists for every pair.
Their outward directions have pairwise nonpositive dot
products, as otherwise one end carries the absent label
into another's full forward half-plane. If none are
antipodal, the reverse of each direction belongs locally
at infinity to its third label.
[Reverse third-label cones](THREE_CELL_SECOND_ORDER_PARAMETER_DETECTION.md#eliminating-the-reverse-cone-hypothesis-for-r12)

A straight end at offset \(b\) has \(F(t,b)=\beta\) for
all \(t\) by analytic continuation along the supporting
line. Its reverse third-label cone then forces \(\beta=0\)
and \(b=0\). An antipodal pair of straight ends is
impossible. In the coordinates of the antipodal case
in the linked argument, let \(p,q\) be the upper and
lower vertical line offsets, with upper pair contrast
\(U_2-U_3\) and lower pair contrast \(U_1-U_3\).
Continuation along the upper line gives
\(\psi((r-1)p)=-\Phi((q-rp)/\sigma)\), where
\(\psi(v)=2\Phi(v/\sigma)-1\). It follows that \(p>0\)
and \(q<p\), since \(2\Phi((1-r)p/\sigma)-1<
\Phi((1-r)p/\sigma)\). Continuation along the lower
line gives \(p<q\).

Fix the first of the three centered straight directions
as \((1,0)\). The other two directions have nonpositive
horizontal components and transverse components with
opposite signs. At least one of the horizontal components
is strictly negative, since otherwise these two directions
are antipodal. For nonvertical directions, write the upper and lower
straight lines in this coordinate frame as \(v=-m_+s\)
and \(v=m_-s\), respectively, for \(s<0\), with
\(m_+,m_->0\).

For a fixed sufficiently negative \(s\), the two other
forward half-planes cover the entire transverse line
\(\{(s,v):v\in\mathbb R\}\). Indeed the upper half-plane
contains \(v>c_++d_+s\) and the lower contains
\(v<c_--d_-s\), where \(d_+,d_-\ge0\) and at least
one is positive. Each half-plane's label assignment
is exact beyond its cutoff. Overlap with the first arm's
forward half-plane at large positive or negative \(v\)
identifies the outer label of the upper or lower arm,
respectively. On the reverse slice the upper and lower
labels of the first pair therefore
occupy the two opposite tails beyond the other two
supporting lines, with the spectator between. When
one of the other directions is vertical, its supporting
line is \(s=0\). On \(s<0\) the corresponding upper or
lower first-pair label is entirely absent, while the
other has positive slice mass.

On the other hand, \(F(t,0)=0\) holds for every real \(t\).
Injectivity of longitudinal Gaussian convolution gives
\[
    \int (1_{A_i}-1_{A_j})(s,v)p_\sigma(v)\,dv=0
                     \quad\text{for almost every }s.    \tag{PG10}
\]
In particular the case of a vertical other direction
is impossible. With both slopes finite, equal Gaussian
tails in (PG10) force \(m_+=m_-\). Repeating for each
end proves that all three angular gaps are \(2\pi/3\).
Their three forward half-planes cover the complement
of a sufficiently large disk. They agree there with
the centered regular simplex, up to null sets.

With connected strict-winner cells and exactly three proper
ends, the existing network argument extends this conclusion
to the whole plane. Without that interior network condition,
the conclusion here is exterior agreement.
[Network rigidity from straight ends](THREE_CELL_ISOLATED_END_RIGIDITY.md#three-pairwise-obtuse-ends)

In the non-antipodal three-end class, even one straight end
must have zero pair pressure. With a bounded isolated graph
for each pair, shifted-contour straightening gives zero
pressure for all three pairs.

## Moment constraints on the compact discrepancy

Under the three-straight-end conclusion, let \(T\) be the
exterior simplex and \(\mu=\gamma_{\sigma^2 I_2}\). Write
\(e_{ij}\) for its outward interface direction for the pair
\(i,j\). In the coordinates of (PG10), multiplying the
slice identity by \(p_\sigma(s)w(s)\) and integrating gives
\[
 \int w(e_{ij}\cdot x)(1_{A_i}-1_{A_j})(x)\,d\mu(x)=0
                                                               \tag{PG11}
\]
for every \(w\) integrable against \(N(0,\sigma^2)\).
In particular, define the symmetric tensor
\[
  M_n(A_i)=\int x^{\otimes n}1_{A_i}(x)\,d\mu(x)
  \qquad(n\ge0).
\]
For every \(n\), the pair moment difference is perpendicular
to \(e_{ij}^{\otimes n}\). The signs are those of the contrast
in (PG11), with the first label's tensor minus the second's.
The same equations hold for \(T\), and hence for the moments
of the compact discrepancies \(1_{A_i}-1_{T_i}\).

In degree zero this gives \(\mu(A_i)=1/3\) for all \(i\).
Let \(w_i\) be the unit center of \(T_i\), so \(w_1+w_2+w_3=0\)
and \(w_i-w_j\) is normal to \(e_{ij}\).
The degree-one equations give
\(M_1(A_i)-M_1(A_j)=c_{ij}(w_i-w_j)\).
Summing cyclically and using independence of two edge vectors
shows that the three coefficients \(c_{ij}\) are equal.
Since the three vector moments sum to zero, there is a real
\(c\) with \(M_1(A_i)=c w_i\) for all \(i\). The orientation
of \(c\) is left free by the projection conditions.

If \(M_1(T_i)=c_T w_i\), then \(c_T>0\) and
\[
   -3c_T=\int \min_i w_i\cdot x\,d\mu(x)
    \ \le\ 3c=\sum_i \int_{A_i} w_i\cdot x\,d\mu(x)
    \ \le\ \int \max_i w_i\cdot x\,d\mu(x)=3c_T.
\]
The upper equality forces \(A=T\) by the pointwise score
comparison. The lower equality forces \(A=-T\), impossible
with the given exterior. Thus any nontrivial compact
discrepancy has a strict degree-one loss.

These are moments in the exact stability expansion. The
joint density of the source and observation gives
\[
 S_r(A)=\sigma^2\sum_{n\ge0}\frac{(r/\sigma^2)^n}{n!}
                         \sum_i \|M_n(A_i)\|^2 .          \tag{PG12}
\]
The norm is the Euclidean tensor norm, with degree zero
interpreted as mass. To obtain (PG12), factor the joint
density as \(\sigma^2 e^{r x\cdot y/\sigma^2}\,d\mu(x)d\mu(y)\)
and expand the exponential. Absolute integrability follows
from \(r|x\cdot y|\le r(|x|^2+|y|^2)/2\) and \(r<1\).
The degree-zero terms for \(A\) and \(T\) agree. If \(A\ne T\),
the loss in degree one is \(3r(c_T^2-c^2)>0\). A tie
would have to recover this amount in degrees \(n\ge2\),
subject to the separate pair projection constraint in each degree.

The projection constraints permit nonzero compact discrepancies.
For example, in a fixed annulus \(0<R_1<|x|<R_2<\infty\)
assign \(x\) to the label of \(-x\) in \(T\) and keep the
simplex assignment elsewhere. Every annular slice has equal
angular masses, so the modified partition is balanced.
Reflection in \(e_{ij}\) interchanges labels \(i,j\), because
the antipodal map commutes with that reflection and the
annulus is radial. The contrast is therefore odd in the
normal coordinate of (PG10), giving every constraint (PG11).
This example tests the projection equations in isolation.

[Bounded components and fixed-shape stationary correlations for a compact discrepancy](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md)
