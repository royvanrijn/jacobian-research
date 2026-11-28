# A Mehler triple lift for coherent Gaussian ribbon comparisons

On the regular three-cell Gaussian tripod, the cross-ray Mehler
kernel is the two-coordinate marginal of a positive exchangeable
measure on three copies of the outward ray. This gives a nonnegative
cross-ray form whenever the three local ribbon labels can be jointly
assigned with the three hard labels available as anchors.

In particular, fix \(L<\infty\). For correlations sufficiently close
to one, no nonzero bounded ribbon profile supported in
\(s\le L\sqrt{1-\rho^2}\) has a negative second-order deficit if its
agreement kernel passes either the pointwise five-point label test
or the four-point test with shared three-point marginals. The
feature direction can vary along every ray. The comparison concerns
the ribbon perturbations below. If each ray has vector mean zero,
there is no noise-scale support restriction and the conclusion
holds for \(9/10\le\rho<1\). General agreement kernels and
unrestricted support with only total vector mean zero remain open.
[The three-ribbon construction with independent local laws](GAUSSIAN_AGREEMENT_MULTI_FACET_COUNTEREXAMPLE.md)

## A positive exchangeable lift

For \(0<\rho<1\) and \(u,v\ge0\), put
\[
 K_{-,\rho}(u,v)=e^{-(u^2+v^2+\rho uv)/2}
              =g(u)g(v)h(u+v),\qquad
 g(t)=e^{-(2-\rho)t^2/4},\quad h(t)=e^{-\rho t^2/4}.
\]
There is a continuous nonnegative function \(H\) on the half-line
such that the measure with density
\[
            g(u)g(v)g(w)H(u+v+w)\,du\,dv\,dw
                                                               \tag{RT1}
\]
has all three pair marginals \(K_{-,\rho}(u,v)\,du\,dv\).

To construct \(H\), let \(W_1,W_2,\ldots\) be independent with
probability density \(-g'\) on the half-line. Its integral is one
by the endpoint values of \(g\). Set \(S_0=0\),
\(S_j=W_1+\cdots+W_j\), and
\[
                    H(t)=\sum_{j\ge0}\mathbb E[-h'(t+S_j)].
\]
The summands are nonnegative. For any fixed \(c>0\) one has
\(-h'(t)\le C_c e^{-ct}\), and
\(\mathbb E e^{-cW_1}<1\). Hence the series converges uniformly
and is exponentially bounded. It satisfies
\(H(t)-\int_0^\infty(-g'(w))H(t+w)\,dw=-h'(t)\).
Writing \(I(t)=\int_t^\infty g(x-t)H(x)\,dx\) and differentiating
the translated \(g\) factor gives \(I'(t)=h'(t)\).
The exponential bound justifies the differentiation and gives
\(I(t)\to0\). Since \(h(t)\to0\), \(I(t)=h(t)\), which proves
(RT1). The measure is finite because its pair marginal is finite.

## What the local label laws enforce

Take three orthogonal unit label features \(e_1,e_2,e_3\).
A ribbon state on a pair \(ij\) has feature
\((e_i+e_j)/2+u_{ij}/\sqrt2\), with \(u_{ij}\) a unit vector
orthogonal to the labels. No common choice of \(u_{ij}\) is
required for different states on that ray. The hard anchor labels
are mutually distinct, and the agreements of a ribbon state with
them are \(1/2,1/2,0\) on its two labels and the absent label.

A five-point law on two hard anchors and one state from each of
the three ribbons fixes the binary label choices for all three
states. To see this, name the two anchored labels \(1,2\), and
name the remaining label \(3\). The ribbon \(12\) chooses one of
its two anchored labels with equal probabilities. The ribbon
\(13\) chooses label \(1\) with probability \(1/2\), never
chooses label \(2\), and chooses label \(3\) otherwise. The
case of \(23\) is the same. These binary signs are thus jointly
realized in the five-point law. With the pairs cyclically
oriented as \(12,23,31\), let \(R_e=+1\) mean that the
ribbon chooses the first label of \(e\). Two different
ribbons meet in one label, which requires opposite signs
in those two pairs. By their zero sign means, for distinct
ribbons \(e,f\) their agreement is
\[
 \frac{1-\mathbb E R_eR_f}{4}
               =\frac14+\frac12\langle u_e,u_f\rangle .
                                                               \tag{RT2}
\]

Compatible four-point laws already give these joint binary
signs. Here compatibility means that the laws on all
four-point sets have the same restrictions on their shared
three-point sets. Consider ribbon states \(A_{12},B_{13},C_{23}\)
with respective pair agreements
\[
 p=\Pr(A=B),\qquad q=\Pr(A=C),\qquad r=\Pr(B=C).
\]
A law with anchors \(1,2\) and states \(A,B\) pins their
two label choices as above. It gives probability \(p\) to
both choosing anchor \(1\), and also probability \(p\) to
both avoiding that anchor, by the unbiased marginals. In
the latter event they are different from each other.
By compatibility these assertions also hold in the law with
anchor \(1\) and all three ribbons. In that law \(C\) never
chooses anchor \(1\). When both \(A,B\) avoid the anchor,
\(C\) agrees with exactly one of them. Whenever they do not
both choose the anchor it can agree with at most one, and
when they both choose it \(C\) agrees with neither. This gives
\(p\le q+r\le1-p\). Repeating at the other two anchors gives
\[
                p\le q+r,\quad q\le p+r,\quad r\le p+q,
                       \qquad p+q+r\le1.                    \tag{RT3}
\]
These are exactly the unbiased binary feasibility
inequalities. For example let \(X,Y\) denote the signs of
\(A,B\) choosing label \(1\), and \(Z\) the sign of \(C\)
choosing label \(2\). The prescribed pair correlations are
\(4p-1,\ 1-4q,\ 4r-1\) on \(XY,XZ,YZ\), respectively.
The law on \(\{\pm1\}^3\) with mass
\((1+\sum_{i<j}r_{ij}\eta_i\eta_j)/8\) at \(\eta\) is
nonnegative by (RT3) and has these correlations and zero
means. Converting to the cyclic conventions gives (RT2).
This argument needs no opposite-sign companion state on
any ribbon.

## The variable-direction form

On ray \(e=ij\), let \(W_e(s)\) be a bounded profile with
values in a common real Hilbert space orthogonal to the hard
labels. For each \(s\) with nonzero profile, give its
symmetric strip about the ray width determined by
\[
 \int_{-h_\epsilon(s)}^{h_\epsilon(s)}e^{-t^2/2}\,dt
                         =2\epsilon\|W_e(s)\|,
\]
and use feature
\((e_i+e_j)/2+W_e(s)/(\sqrt2\|W_e(s)\|)\) inside the strip.
Use the hard simplex features elsewhere. For actual strips
take compact radial supports away from the origin, so the
strips are disjoint for sufficiently small \(\epsilon\).
The comparison permits any lower support cutoff. A zero
Gaussian arclength mean of the total profile ensures exact
row balance, by integrating the strip features. We use no
mean condition for the quadratic estimate.

With \(d\nu(s)=(2\pi)^{-1}e^{-s^2/2}\,ds\), write \(M_+,M_-\)
for the same-ray and distinct-ray Mehler kernels, and
\(a_\rho(s)=\rho(2M_+1(s)+M_-1(s))\) for the ordinary
simplex score coefficient. The exact second-order coefficient
of the relaxed objective deficit is
\[
 \begin{split}
 Q_\rho(W)={}&\sum_e\int a_\rho(s)\|W_e(s)\|^2\,d\nu(s)
      -2\sum_e\iint M_+(s,t)\langle W_e(s),W_e(t)\rangle\,d\nu(s)d\nu(t)\\
 &\hspace{6mm}-4\sum_{e<f}\iint
          M_-(s,t)\langle W_e(s),W_f(t)\rangle\,d\nu(s)d\nu(t).
 \end{split}                                                 \tag{RT4}
\]
Indeed the doubled soft score regret on a strip is
\(\epsilon^2\int a_\rho\|W_e\|^2\,d\nu+o(\epsilon^2)\).
The transverse integral of its feature change is
\(\sqrt2\epsilon W_e(s)\), since the soft label changes
cancel by normal symmetry. In the quadratic interaction
restricting the kernel to the rays costs \(o(\epsilon^2)\)
on compact supports. This is the same expansion as for the
constant-direction ribbon, with the inner products retained.
[Raw score coefficient and ribbon expansion (MW2)–(MW3)](GAUSSIAN_AGREEMENT_MULTI_FACET_COUNTEREXAMPLE.md#the-raw-three-ray-variation)

Let \(\sigma=\sqrt{1-\rho^2}\), write
\(\widetilde W_e(u)=W_e(\sigma u)\), and set
\(K_{+,\rho}(u,v)=e^{-(u^2+v^2)/2+\rho uv}\).
Under this scaling \(M_\pm(s,t)d\nu(s)d\nu(t)\) becomes
\((2\pi)^{-2}K_{\pm,\rho}(u,v)\,du\,dv\).
For fixed \(u_1,u_2,u_3\), put \(a_e=\|\widetilde W_e(u_e)\|\).
When all are nonzero, (RT2) gives
\[
 \sum_e a_e^2
    -4\sum_{e<f}\langle \widetilde W_e(u_e),\widetilde W_f(u_f)\rangle
                  =\mathbb E\left(\sum_e a_eR_e\right)^2\ge0.
\]
If some vanish, just use the pair law for the nonzero
coefficients, or the single remaining square. The laws
need not be selected measurably: their existence proves the
pointwise inequality for the displayed measurable form.
Integrate against (RT1). By the three pair margins this proves
\[
 D_-:=\sum_e\int K_{-,\rho}1(u)\|\widetilde W_e(u)\|^2\,du
    -4\sum_{e<f}\iint K_{-,\rho}(u,v)
             \langle\widetilde W_e(u),\widetilde W_f(v)\rangle\,du\,dv
 \ \ge\ 0.                                                  \tag{RT5}
\]
Here \(K_{\pm,\rho}1(u)=\int_0^\infty K_{\pm,\rho}(u,v)\,dv\).
Thus the cross-ray form is positive without assuming a global
joint law or compatible *weighted* pair marginals. The profile
amplitudes occur inside the square.

The remaining part of \((2\pi)^2Q_\rho(W)\) is the sum of
the same-ray terms
\[
 \sum_e\iint K_{+,\rho}(u,v)
          \|\widetilde W_e(u)-\widetilde W_e(v)\|^2\,du\,dv
\]
minus \((1-\rho)\sum_e\int (2K_{+,\rho}1+K_{-,\rho}1)
\|\widetilde W_e\|^2\). If the scaled supports lie in
\([0,L]\), restrict one argument of the difference terms to
\([L+1,L+2]\), where the profiles vanish, and the other to
\([0,L]\). Counting both orders gives a lower bound
\(2\exp[-(L^2+(L+2)^2)/2]\sum_e\|\widetilde W_e\|_2^2\).
Both kernel row integrals are at most \(\sqrt{2\pi}\),
since \(K_{+,\rho}\le e^{-(u-v)^2/2}\) and
\(K_{-,\rho}\le e^{-(u^2+v^2)/2}\). Together with (RT5),
\[
 (2\pi)^2 Q_\rho(W)\ge
  \left[2e^{-(L^2+(L+2)^2)/2}-3\sqrt{2\pi}(1-\rho)\right]
                           \sum_e\|\widetilde W_e\|_2^2 .
                                                               \tag{RT6}
\]
The coefficient is positive for \(\rho\) sufficiently close
to one depending only on \(L\). Thus this compact
noise-scale class of variable-direction perturbations loses
strictly at second order under either of the local label-law
conditions. A separate comparison below removes the noise-scale
support condition when each ray has vector mean zero. Both
statements concern the ribbon model.

The standard mean repair does not change the conclusion.
Suppose the actual profile is obtained from one with scaled
support in \([0,L]\) by subtracting on each ray its vector
arclength mean times a fixed compactly supported scalar bump
of arclength integral one, away from the origin. For small
\(\sigma\), the local label laws on the original supports
are unchanged by this repair, so (RT6) applies to the
unrepaired profiles \(W^0_e\). Their arclength means have
norm at most \(C_L\sigma\|\widetilde W_e^0\|_2\).
The score row \(a_\rho\) is uniformly \(O(\sigma^{-1})\);
thus the squared \(L^2(a_\rho\,d\nu)\) norm of the repair
is at most \(C_L\sigma\sum_e\|\widetilde W_e^0\|_2^2\).
The original profiles have squared weighted norm at most
\(C\sum_e\|\widetilde W_e^0\|_2^2\), by the scaled kernel
row bounds. Cauchy–Schwarz with the positive symmetric
kernel in each interaction of (RT4) bounds it by the product
of the weighted norms, using \(M_\pm1\le a_\rho/\rho\).
Therefore the change in \(Q_\rho\) is
\(O_L(\sqrt\sigma)\sum_e\|\widetilde W_e^0\|_2^2\), which is
absorbed by the uniform positive coefficient in (RT6).

## Per-ray mean zero removes the noise-scale support condition

Assume \(\int W_e\,d\nu=0\) on each ray separately. Under
either local label-law hypothesis, every nonzero bounded compact
profile with this property has \(Q_\rho(W)>0\) for
\(9/10\le\rho<1\), with no bound on its support in units of
\(\sqrt{1-\rho^2}\).

Use the normalized half-Gaussian measure
\(d\eta(s)=\sqrt{2/\pi}\,e^{-s^2/2}ds\) on the positive
half-line. In these coordinates all terms in \(Q_\rho\)
have the common factor
\((2\sqrt{2\pi})^{-2}\sigma^{-1}\), since
\(d\nu=(2\sqrt{2\pi})^{-1}d\eta\).
Write \(m_c(s,t)=\sigma M_{\rho,c}^{[2]}(s,t)\) for the
Mehler kernel with angle cosine \(c\) after removing that
factor. The kernel \(m_1\) is the standard one-dimensional
Mehler kernel and \(m_{-1}(s,t)=m_1(s,-t)\).
On \(L^2(\eta)\), for a single vector profile \(f\), put
\[
 \mathcal E(f)=\frac12\iint(m_1+m_{-1})
                \|f(s)-f(t)\|^2\,d\eta(s)d\eta(t),\qquad
 \mathcal R(f)=\frac12\iint(m_1-m_{-1})
                \|f(s)-f(t)\|^2\,d\eta(s)d\eta(t).
\]
The same-ray difference form is \(\mathcal E+\mathcal R\),
and \(\mathcal R\) is nonnegative. The even extension of
\(f\) identifies \(\mathcal E\) with
\(2\sum_{j\ge0}(1-\rho^{2j})\|f_{2j}\|_2^2\), where
\(f_{2j}\) denotes its degree-\(2j\) Hermite projection.

For \(x=\rho s/\sigma\), direct Gaussian integration gives
\(m_1 1(s)=2\Phi(x)\) and
\(m_{-1/2}1(s)=2e^{-3x^2/8}\Phi(-x/2)\). The latter is at
most \(4\Phi(-x)\). Indeed, write the Gaussian tail as
\(\Phi(-x)=\phi(x)\int_0^\infty e^{-xy-y^2/2}\,dy\);
in the tail with \(x/2\), the change of variable \(y=2z\)
gives the bound since \(e^{-2z^2}\le e^{-z^2/2}\).
Hence \(2m_1 1+m_{-1/2}1\le4\).
Dropping the nonnegative cross-ray term (RT5), the
normalized \(Q_\rho\) is therefore at least the sum, over
the three rays, of
\[
                 \mathcal E(f)+\mathcal R(f)
                            -4(1-\rho)\|f\|_2^2 .
\]

Put \(\delta=1-\rho\) and split the even, mean-zero profile
as \(f=g+h\) into degree two and degrees at least four.
For \(g(s)=c(s^2-1)\) with a constant vector \(c\),
\(\mathcal E(g)=2\sigma^2\|g\|_2^2\).
The reflected difference term with kernel \(m_{-1}\) is
at most \(2\sigma^5\|g\|_2^2\). Indeed \(s=\sigma u\),
\(t=\sigma v\) gives its bound as
\((2/\pi)\sigma^5\|c\|^2\) times
\(\iint_{u,v\ge0}(u^2-v^2)^2 e^{-(u^2+v^2)/2}\,du\,dv=2\pi\),
and \(\|g\|_2^2=2\|c\|^2\). Therefore
\(\mathcal R(g)\ge2\sigma^2(1-\sigma^3)\|g\|_2^2\).

Also \(\mathcal E(h)\ge2(1-\rho^4)\|h\|_2^2\) and
\(\mathcal R(h)\le\mathcal E(h)\). The seminorm triangle
inequality and Cauchy–Schwarz give
\(\mathcal R(g)\le4[\mathcal R(f)+\mathcal E(h)/3]\).
Combining these with the orthogonal splitting of
\(\mathcal E\) shows that each ray contributes at least
\[
 \left[\frac{\sigma^2(1-\sigma^3)}2-2\delta^2\right]\|g\|_2^2
             +\left[\frac43(1-\rho^4)-4\delta\right]\|h\|_2^2 .
                                                               \tag{RT7}
\]
Both brackets are strictly positive for \(9/10\le\rho<1\).
For the first use \(\sigma^3<1/10\), and for the second
use \(1+\rho+\rho^2+\rho^3>3\). This proves the claim.

## Finite coherence does not imply a global representation

For every finite level \(\ell\ge2\) there is a balanced
unit-feature kernel on a single ribbon with coherent
\(\ell\)-point label laws that has no global representation
as a mixture of three-label agreement kernels. This is a
representation obstruction, with no claimed objective gain.

Choose an odd \(N>\ell\) and divide a compact subinterval
of one ray into \(N\) states of equal Gaussian arclength
measure. Use unit residual directions at the vertices
of a centered regular \(N\)-simplex. Their off-diagonal
inner products are \(-1/(N-1)\). The widths are all equal
on these states, so the residual mean is zero, and
normal symmetry gives the usual three equal soft masses.
Outside the ribbon keep the three hard anchors.

On any tuple use the deterministic hard labels for its
anchors and assign its distinct ribbon states binary signs.
For those signs use the corresponding restriction of an
exchangeable unbiased law on \(\ell\) signs with pair
correlation \(-1/(N-1)\). Such a law is a mixture of
independent signs and the uniform sign law restricted to
strings of least absolute sum. The latter correlation is
\(-1/(\ell-1)\) for even \(\ell\) and \(-1/\ell\) for odd
\(\ell\), by expanding the squared sign sum. Both are
at most \(-1/(N-1)\). The restrictions are coherent by
exchangeability, and repeated occurrences of a state
have the same sign. Same-ribbon pair agreements are
\((1+\mathbb E R_iR_j)/2\), and the hard-anchor
agreements come from the unbiased marginals. Thus the laws
realize the feature kernel on every such tuple.

A global representation would have three distinct anchor
labels almost surely. Relative to these anchors, its
\(N\) ribbon states would have binary signs with the
specified correlations. It would follow that
\(\mathbb E(\sum_iR_i)^2=N-N(N-1)/(N-1)=0\), impossible
since \(N\) is odd. Thus finite coherence alone cannot
extend the constant-direction ribbon representation to
arbitrary feature directions.
