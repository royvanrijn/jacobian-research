# Exact stability comparison for balanced reflected affine wedges

For every \(0<\rho<1\), the regular simplex uniquely maximizes
Gaussian noise stability within the family
\[
\begin{aligned}
A_3(m)&=\{(x,y):y>m|x|+b(m)\},\\
A_1(m)&=\{(x,y):x>0,\ y<mx+b(m)\},\\
A_2(m)&=\{(x,y):x<0,\ y<-mx+b(m)\},
\end{aligned}
\qquad m\ge0 ,
\]
where \(b(m)\) gives the third cell Gaussian mass \(1/3\).
In fact the stability is strictly increasing up to
\(m_*=1/\sqrt3\) and strictly decreasing thereafter. The proof
uses the exact score on the entire moving boundary, including
the contribution of the reflected half-plane.

## Balance and the stability derivative

Write \(Q=1-\Phi\), where \(\Phi\) is the standard Gaussian
distribution function, and let \(\phi=\Phi'\). The balancing
equation is
\[
2\int_0^\infty \phi(x)Q(mx+b(m))\,dx=\frac13 .
\tag{1}
\]
Its left side is continuous and strictly decreasing in \(b\)
from one to zero. Differentiation and the implicit function
theorem give
\[
b'(m)=-\frac{\int_0^\infty x w_m(x)\,dx}
                  {\int_0^\infty w_m(x)\,dx},\qquad
w_m(x)=\phi(x)\phi(mx+b(m)).
\tag{2}
\]
The equations hold at \(m=0\) as right derivatives, or using
the smooth implicit extension to negative \(m\).
Differentiating (1) is justified by a constant multiple of
\((1+x)\phi(x)\) on each compact \(m\)-interval.
In particular \(b\) is strictly decreasing. At \(m_*\) the
centered wedge has opening \(120^\circ\), so \(b(m_*)=0\).

Put \(F_m=T_\rho(1_{A_3(m)}-1_{A_1(m)})\) and
\(f_m(x)=F_m(x,mx+b(m))\) for \(x>0\). Raising the right
boundary releases mass from label 3 to label 1, and raising
the left boundary makes the reflected transfer. If
\(S(m)=\sum_i\langle1_{A_i(m)},T_\rho 1_{A_i(m)}\rangle\),
the first variation is therefore
\[
S'(m)=-4\int_0^\infty f_m(x)(x+b'(m))w_m(x)\,dx .
\tag{3}
\]
For a direct justification, differentiate the two identical
endpoint contributions to the joint Gaussian integral.
In a compact \(m\)-interval the absolute boundary velocity
is bounded by \(C(1+x)\), and each conditional score is
bounded by one. The resulting integrand is dominated by a
constant multiple of \((1+x)\phi(x)\). Equivalently, in the
quadratic expansion in indicator differences the first
variations are the fiber endpoint derivatives, and the
quadratic term is \(o(|\Delta m|)\). This latter bound follows
from the Mehler bound
\(M_\rho(u,v)\le (1-\rho^2)^{-1}
\exp[\rho(|u|^2+|v|^2)/(2(1+\rho))]\)
and the \(O(|\Delta m|)\) Gaussian weighted measure of the
changed strips. Here the exponent coefficient is strictly
below \(1/2\), so the weighted strip estimate is again
dominated by a constant multiple of
\((1+x)e^{-c x^2}\).

After using (2), the right side of (3) is
\(-4\int_0^\infty w_m\) times the covariance of \(f_m(x)\)
and \(x\) under the probability density proportional to
\(w_m\). It remains to determine the monotonicity of \(f_m\).

## Monotonicity of the folded score

Let \(\sigma=\sqrt{1-\rho^2}\), \(z=\rho x/\sigma\),
\(d=(1-\rho)b(m)/\sigma\), and \(q=mz-d\). The unreflected
binary half-plane score of the line \(y=mx+b(m)\),
evaluated at its boundary, is
\(\lambda=2\Phi(-d/\sqrt{1+m^2})-1\).
The source on \(x<0\) contributes the exact correction
\[
f_m(x)=\lambda+C(z),\qquad
C(z)=2A(z)-B(z),
\]
\[
A(z)=\int_0^\infty\phi(z+t)Q(q+mt)\,dt,\qquad
B(z)=\int_0^\infty\phi(z+t)Q(q-mt)\,dt.
\tag{4}
\]
Indeed the binary source takes the values \(+1\) above
and \(-1\) below the extension of the right affine line.
The true source on the negative half-plane takes the
values \(+1\) above the reflected line and zero below
it, since label 1 is absent. Conditioning at horizontal
source \(-\sigma t\) gives the difference
\(2Q(q+mt)-Q(q-mt)\). These are also the folded affine
integrals used in the
[tail comparison for reflected maximizing graphs](THREE_CELL_REFLECTED_SHAPE_RIGIDITY.md#global-convexity-or-concavity-forces-the-simplex).

Changing variables from \(t\) to \(u=z+t\) in (4) gives
\[
C'(z)=-\phi(z)Q(q)+2m\int_z^\infty
                      \phi(u)\phi(2mz-d-mu)\,du .
\]
The derivative follows from the endpoint term and the
Gaussian integrable bound on the derivatives of the
vertical tail. Define the strictly decreasing positive
function \(I(a)=\int_0^\infty e^{-at-t^2/2}\,dt\), for
which \(Q(a)=\phi(a)I(a)\). Factoring out
\(\phi(z)\phi(q)\) from the last integral and rescaling
the remaining integration variable yields
\[
\frac{C'(z)}{\phi(z)Q(q)}
=\frac{2m}{\sqrt{1+m^2}}\,
 \frac{I(((1-m^2)z+md)/\sqrt{1+m^2})}{I(mz-d)}-1 .
\tag{5}
\]
If \(0\le m<m_*\), then \(d>0\), the prefactor in (5)
is less than one, and
\[
\frac{1-m^2}{\sqrt{1+m^2}}\ge m,\qquad
\frac{md}{\sqrt{1+m^2}}\ge-d .
\]
Thus the \(I\)-ratio is at most one and \(C'<0\).
If \(m>m_*\), then \(d<0\), both argument comparisons
reverse and the prefactor exceeds one, so \(C'>0\).
The slope comparison for \(m\ge1\) follows immediately
from its nonpositive left side. For \(m<1\) it follows
by squaring the nonnegative sides and comparing \(1\)
with \(3m^2\). At \(m=m_*\) all comparisons are equalities
and \(C'=0\).

Since \(z\) strictly increases with \(x\), \(f_m(x)\) is
strictly decreasing before \(m_*\) and strictly increasing
after it. The covariance in (3) has the same strict sign
as that monotonicity, by writing it as half the mean of
\((x-x')(f_m(x)-f_m(x'))\) for two independent draws from
the normalized density. This proves the claimed stability
comparison. The result compares affine folded graphs at
finite amplitude. General reflected maximizing graphs
can have nonconstant slope, so their exclusion requires
an additional argument.

## A one-sign obstruction for a curved graph

Suppose a continuous reflected balanced stationary graph
has height \(h(x)=\ell(x)+u(x)\) on the right, where
\(\ell(x)=mx+b\), \(m\ge0\), and \(u(x)\to0\) as \(x\to\infty\).
At a nonsimplex graph the following two possibilities are excluded:
\[
m\ge m_* \ \text{and}\ u(x)\ge0\ \text{for every }x\ge0,
\qquad
m\le m_* \ \text{and}\ u(x)\le0\ \text{for every }x\ge0.
\tag{6}
\]
In these comparisons the residual need not be monotone.

Let \(F_h=T_\rho(1_{A_3(h)}-1_{A_1(h)})\) for the reflected
partition with height \(h\). If \(p_\sigma(t)=\phi(t/\sigma)/\sigma\),
conditioning on the horizontal source gives
\[
F_h(x,y)=
\int_0^\infty [2p_\sigma(s-\rho x)+p_\sigma(s+\rho x)]
     Q((h(s)-\rho y)/\sigma)\,ds
     -\Phi(\rho x/\sigma).
\tag{7}
\]
The horizontal source is in the right half-plane with
probability \(\Phi(\rho x/\sigma)\), where the pair contrast
is \(2\) times the third-cell probability minus one.
On the left the contrast is just the third-cell probability.
As \(x\to\infty\), formula (7) on the graph tends to the
binary affine value \(\lambda=2\Phi(-(1-\rho)b/
(\sigma\sqrt{1+m^2}))-1\). The reflected source contribution
vanishes, and a typical right horizontal source tends
to infinity where \(u\to0\). Stationarity makes this
limiting value the score on the entire graph.

For the affine graph, write its score as
\(F_\ell(x,\ell(x))=\lambda+C(\rho x/\sigma)\) as in (4),
now without assuming balance of \(\ell\). Put \(\ell_t=\ell+tu\).
Subtracting the two scores and integrating the derivative
of (7) along \(\ell_t\) gives
\[
Z(x)\bigl[\rho u(x)-Pu(x)\bigr]=-C(\rho x/\sigma),
\quad Pu(x)=\frac{1}{Z(x)}\int_0^\infty L(x,s)u(s)\,ds,
\quad Z(x)=\int_0^\infty L(x,s)\,ds ,
\tag{8}
\]
where
\[
L(x,s)=[2p_\sigma(s-\rho x)+p_\sigma(s+\rho x)]
       \int_0^1 p_\sigma(\ell_t(s)-\rho\ell_t(x))\,dt .
\]
The kernel \(L\) is continuous and strictly positive. The
normalization \(Z(x)\) has a uniform positive lower bound:
restrict to \(s\in[\rho x,\rho x+1]\), on which
\(m(s-\rho x)+(1-\rho)b+t(u(s)-\rho u(x))\) is uniformly
bounded. Thus the transition density of \(P\) is bounded
above by a constant times \(p_\sigma(s-\rho x)\), since
\(p_\sigma(s+\rho x)\le p_\sigma(s-\rho x)\) for \(x,s\ge0\).
In particular a chain with transition \(P\) satisfies
\(\mathbb E[S_{n+1}\mid S_n=x]\le \rho x+C_0\), uniformly
in \(x\). For any fixed starting point its successive means
are bounded. Markov's inequality therefore places at least
half of each successive law in one fixed compact interval.
By strict positive continuity of the kernel, that interval
has a common positive transition probability into any fixed
positive-length subinterval. If \(u\ge0\) is continuous and
nonzero, choose that subinterval where \(u\) is bounded below
by a positive number. It follows that \(P^n u(x)\) is bounded
below by a positive number for every \(n\ge1\) at the fixed
starting point.

Now suppose \(m\ge m_*\) and \(u\ge0\). Balance of \(h\)
forces \(b\le b(m)\le0\), since replacing \(h\) by \(\ell\)
weakly increases the third-cell measure. The first inequality
is strict if \(u\ne0\). Formula (5), which remains valid with
this intercept, gives \(C'>0\) unless \(m=m_*\) and \(b=0\).
As \(C(z)\to0\), this gives \(C<0\) for every finite \(z\).
Equation (8) implies \(Pu<\rho u\), hence
\(P^n u(x)\le \rho^n u(x)\) for every \(n\). This contradicts
the preceding lower bound when \(u\ne0\). If \(u=0\),
stationarity itself forces \(C=0\), hence \(m=m_*\) and
\(b=0\) under balance. For \(m\le m_*\) and \(u\le0\) the
intercept inequalities reverse, \(C>0\) outside the simplex,
and the same argument applies to \(-u\). This proves (6).

## An isolated affine tail of slope at most one

For an actual reflected maximizer, suppose the graph has an
affine asymptote \(h(x)=mx+b+o(1)\) with \(0\le m\le1\).
Rotate coordinates so the longitudinal coordinate is
\(t=(x+my)/\sqrt{1+m^2}\) along the right asymptote.
On cell 2, where \(x<0\) and \(y<h(|x|)\), one has
\[
t\le\frac{-(1-m^2)|x|+mb+o(1)}{\sqrt{1+m^2}}
\qquad (x\to-\infty).
\]
Consequently cell 2 and the downward axis interface have
longitudinal coordinate bounded above. On the right graph
the coordinate \(t=(x+mh(x))/\sqrt{1+m^2}\) is strictly
increasing to infinity by the
[reflected graph theorem](THREE_CELL_REFLECTION_REDUCTION.md).
Thus in a full forward half-plane every transverse fiber has
exactly one switch between labels 1 and 3. Its transverse
height is \((h(x)-mx)/\sqrt{1+m^2}\), which is bounded on
the tail. This verifies the hypotheses of the
[exact straightening theorem for bounded isolated
ends](THREE_CELL_TAIL_FLATNESS.md).
The right tail is exactly straight, so analyticity makes the
entire right graph affine. The stability comparison above
then identifies the simplex.

For example, a one-sided bound on the fixed-shape correlation
derivative along the right graph forces a linear tail phase
and such an affine asymptote by
[(CR9) and the pressure correction](THREE_CELL_REFLECTED_CONSTANT_CORRELATION_RIGIDITY.md#the-two-possibilities-for-the-derivative-at-infinity).
This implication still requires the resulting slope to be
at most one. The cases of a nonlinear tail phase and a
steep affine tail remain open in this comparison.
