# Uniform small-parameter contraction of the posterior simplex lift

Fix \(k\ge2\) and a compact band \(J=[a,b]\subset(0,1)\).
Let \(T\) be the balanced regular simplex in \(\mathbb R^d\),
\(d\ge k-1\), and let \(A\) be any balanced ordinary partition
in that space. Write
\[
 L_\rho(A;T)=\sum_{i\ne j}
   \int_{T_i\cap A_j}T_\rho(1_{T_i}-1_{T_j})\,d\gamma_d
\]
for the linear score regret. For the posterior simplex lift
\(\mathcal L_s\), there are \(C<\infty\) and \(s_*>0\), depending
only on \(k,a,b\), such that
\[
 0\le L_\rho(\mathcal L_s A;\mathcal L_s T)
     \le C s^2 L_\rho(A;T)
 \qquad(0<s\le s_*,\ \rho\in J).                             \tag{SC1}
\]
In particular the estimate is uniform in \(d\). The lift of \(T\)
is itself a regular simplex. If \(D_\rho(A)=S_\rho(T)-S_\rho(A)\)
denotes the stability deficit, with reference simplices taken in
the appropriate dimension, then
\[
               D_\rho(\mathcal L_s A)\le 2Cs^2 L_\rho(A;T).
                                                               \tag{SC2}
\]
[Posterior simplex lift and model prices](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#a-sharp-posterior-entropy-and-centroid-tradeoff)

## Score regret in a posterior fiber

Identify the active space of \(T\) with \(H=\mathbf 1^\perp\subset
\mathbb R^k\), write \(x=(u,y)\) with \(u\in H\), and let \(T\)
assign the leader of \(u\). Write \(p=T_s1_A\), \(q=T_s1_T\), and
\(\sigma=\sqrt{1-s^2}\). In the fiber Gaussian \(z\in H\), the
lift of a partition of posterior masses \(v(x)\) assigns
\(\arg\max_i (z_i-\beta_i(v(x)))\), with \(\beta(v)\in H\) the
unique price vector. For \(q\) the prices are \(-su/\sigma\):
the reference lift assigns the leader of \(su+\sigma z\).
The latter map from \((x,z)\) to \(H\) is a coisometry, so
\(\mathcal L_s T\) is a regular simplex.

Let \(I(v)\) be the Gaussian cluster perimeter of the model at
masses \(v\). Its assignment maximizes the conditional expectation
of the chosen coordinate \(z_i\) at the prescribed masses, with
value \(\sqrt2 I(v)\) and gradient \(\beta(v)=\sqrt2\nabla I(v)\)
on \(H\). For a fixed \(x\), the conditional regret of using the
\(p(x)\) assignment in place of the \(q(x)\) assignment for the
linear scores \(su_i+\sigma z_i\) is exactly
\[
 \begin{split}
  \sigma\{\sqrt2[I(q)-I(p)]+\beta(q)\cdot(p-q)\}
       &=\sqrt2\,\sigma\,D_{-I}(p\Vert q),\\
 D_{-I}(p\Vert q)&=-I(p)+I(q)+\nabla I(q)\cdot(p-q).
 \end{split}                                                   \tag{SC3}
\]
The identities for the value, its gradient and the Hessian
\(-D^2 I(v)=L_v^{-1}\), where \(L_v\) is the model facet Laplacian
on \(H\), are the posterior metric formulas.
[Model assignment value and price Hessian](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#posterior-gradients-in-the-perimeter-metric)

For a regular simplex on \(H\), the smoothed contrast
\(T_\rho(1_{T_i}-1_{T_j})\) vanishes on the full pair hyperplane
by reflection. Its directional derivatives are bounded by
\(\rho\mathbb E|G|/\sqrt{1-\rho^2}\) for \(G\sim N(0,1)\), by the
Gaussian convolution derivative. Therefore the contrast in the
active coordinate \(su+\sigma z\), on the side where \(i\) wins,
is at most \(C_b[(su+\sigma z)_i-(su+\sigma z)_j]\). Its
nonnegativity is the positive simplex score margin. Integrating
(SC3) bounds the score regret of the lift by
\(C_b\mathbb E D_{-I}(p\Vert q)\).

## Uniform control of the model Hessian

Set \(m(v)=\min_i v_i>0\). If \(i,j\) achieve the largest and
smallest model prices, the event that \(i\) wins is contained
in \(\{z_i-z_j\ge\beta_i-\beta_j\}\). The difference is Gaussian of
variance two, so \(m(v)\le\exp[-(\beta_i-\beta_j)^2/4]\).
Since \(\beta\) is centered, \(|\beta(v)|^2\le4k\log(1/m(v))\).
The model is the translated regular simplex with translation
\(\beta(v)\). Each of its facets contains the translate of a fixed
unit disk of the corresponding facet of the centered model
(with the zero-dimensional disk interpreted as a point).
The standard Gaussian density on that disk is at least
\(c_k\exp(-|\beta(v)|^2)\), using
\((|\beta|+R_k)^2/2\le |\beta|^2+R_k^2\) for a bound \(R_k\) on
the original disk. Thus all model facet areas are at least
\(c_k m(v)^{4k}\). The complete-graph Laplacian formula on \(H\)
gives
\[
 \|L_v^{-1}\|\le C_k m(v)^{-M},\qquad M=4k .
\]
With \(m_*=\min\{\min_i p_i,\min_i q_i\}\), the segment joining
\(p\) and \(q\) has all coordinates at least \(m_*\). Taylor's
formula and Cauchy–Schwarz now give
\[
 \mathbb E D_{-I}(p\Vert q)
 \le C_k(\mathbb E m_*^{-2M})^{1/2}
                   (\mathbb E|p-q|^4)^{1/2}.                  \tag{SC4}
\]

For \(s^2\le[2(1+2M)]^{-1}\), reverse Gaussian
hypercontractivity with input exponent \(1/2\) and output exponent
\(-2M\) gives
\[
          \big(\mathbb E (T_s1_E)^{-2M}\big)^{-1/(2M)}
                    \ge \gamma_d(E)^2 .
\]
One may first apply it to \(1_E+\epsilon\) and let
\(\epsilon\downarrow0\) by monotone convergence of the negative
moments. Each of the \(2k\) coordinates in \(m_*\) is the posterior
of a set of mass \(1/k\). Hence
\(\mathbb E m_*^{-2M}\le 2k\,k^{4M}\).
[Chen–Dafnis–Paouris, Gaussian hypercontractivity and its reverse, equations (3.13)–(3.14)](https://arxiv.org/pdf/1306.2410v2)

## Smoothing a balanced transfer

Let \(e_i=1_{A_i}-1_{T_i}\). For each spectator value \(y\), use the
weight \(w_\alpha\) of (B2) on the active space, with an upper
correlation bound \(\max\{b,1/2\}\), and set
\(F_i(y)=\int_H w_\alpha(u)|e_i(u,y)|\,d\gamma_H(u)\).
Apply the global version of (B5) in each active slice to each of
the ordered transfers \(T_\ell\to A_j\). Its normal variable is
\(u\cdot (\mathbf e_\ell-\mathbf e_j)/\sqrt2\ge0\). The
uniform saturated score margin of (A4) controls the cost
\(\min\{u\cdot(\mathbf e_\ell-\mathbf e_j)/\sqrt2,1\}\)
with a positive constant depending only on \(k,a,b\).
Square the sum of the finitely many weighted transfer masses
using Cauchy–Schwarz and integrate in \(y\). This proves
\[
                  \sum_i\|F_i\|_2^2\le C L_\rho(A;T).
                                                               \tag{SC5}
\]
[Positive score margin and global weighted mass bound (B5)](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)

For \(0<v\le1/2\), split \(e_i=e_i^0+e_i^\perp\) into
its spectator average and the spectator-mean-zero part.
These subspaces are orthogonal and invariant under \(T_v\).
The active Gaussian mean of \(e_i^0\) is zero by balance.
The centered Mehler bound (B4) on \(H\) gives
\[
  0\le\langle e_i^0,T_v^H e_i^0\rangle
       \le C v\left(\int_H w_\alpha|e_i^0|\,d\gamma_H\right)^2
       \le C v\|F_i\|_2^2 .
\]
On the mean-zero spectator space \(0\le T_v^\perp\le v\,\mathrm{Id}\)
by the Hermite expansion. Since \(T_v^H\) is positive
semidefinite, the product operator obeys the corresponding order.
The uncentered active Mehler bound in the proof of (B4) gives
\[
 \begin{split}
 0\le \langle e_i^\perp,T_v e_i^\perp\rangle
 &\le v\langle e_i^\perp,(T_v^H\otimes\mathrm{Id})e_i^\perp\rangle\\
 &\le C v\|F_i+\mathbb E F_i\|_2^2
 \le 4 C v\|F_i\|_2^2 .
 \end{split}
\]
All the Mehler constants here are in dimension \(k-1\). Summing
and using (SC5) gives
\(\sum_i\langle e_i,T_v e_i\rangle\le C v L_\rho(A;T)\).
[Centered and uncentered active Mehler bounds near (B4)](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)

For \(3s^2\le1/2\), square the \(L^2\to L^4\) hypercontractive
bound and use \(T_s=T_{1/\sqrt3}T_{\sqrt3 s}\). The triangle
inequality in \(L^2\) for the sum of squares gives
\[
 (\mathbb E|p-q|^4)^{1/2}
 \le \sum_i\|T_s e_i\|_4^2
 \le \sum_i\|T_{\sqrt3s}e_i\|_2^2
 =\sum_i\langle e_i,T_{3s^2}e_i\rangle
 \le C s^2 L_\rho(A;T).
\]
Together with (SC3)–(SC4) this proves (SC1). Both lifts are
balanced, so (B1) and the positive semidefiniteness of \(T_\rho\)
give (SC2).

## Strict uniform comparison on a band with uniqueness and a local gap

Assume the balanced simplex is the unique unrestricted global
maximizer at every \(\rho\in J\), up to the simplex orbit and
null sets. Also assume the dimension-uniform local regret growth
(DL1) on \(J\), which holds for \(k=3\). For the moment-aligned
reference \(T_A\) there is \(c_J>0\), independent of \(d\), such
that
\[
                      D_\rho(A)\ge c_J L_\rho(A;T_A)
                            \qquad(\rho\in J).                \tag{SC6}
\]
Indeed, (DL1) provides this estimate for \(d_k(A)<\eta\). If the
infimum of \(D_\rho(A)\) on \(d_k(A)\ge\eta\) were zero, take a
sequence with deficits tending to zero and pass to a subsequence
with correlations tending to a point of \(J\). Dimension-free
near-maximizer compactness and uniqueness at the limiting
correlation give \(d_k(A)\to0\) along a further subsequence,
a contradiction. This supplies a positive deficit uniformly
on the complement. There \(L_\rho(A;T_A)\le1\), since each point
contributes at most one contrast of magnitude at most one.
[Dimension-free local growth (DL1)](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md)
[Near-maximizer compactness with varying dimensions and correlations](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)

Choose \(s_0\le s_*\) with \(2Cs_0^2<c_J\). For every nonsimplicial
balanced \(A\), (SC2) and (SC6) imply
\[
 S_\rho(\mathcal L_s A)>S_\rho(A)
 \qquad (0<s\le s_0,\ \rho\in J).                              \tag{SC7}
\]
For a simplex equality holds because the lift is a simplex.
Without the uniqueness hypothesis, (SC1)–(SC2) still hold at
every \(\rho\in J\). In a fixed dimension the same argument for
(SC7) only needs uniqueness and (B11) in that dimension, using
the compactness of \(T_{\rho_*}\) in finite dimension. In detail,
weakly convergent indicators have simplex-valued, balanced
limits, while compactness and the continuity of \(T_\rho\) in
operator norm at \(\rho_*\in J\) preserve their stability
limits. Strict positive definiteness of \(T_{\rho_*}\)
forces any maximizing fractional tuple to be ordinary, by the
usual two-label mass-preserving splitting. When the weak limit
is ordinary, the fixed squared \(L^2\) norms and the norm of the
limit coincide, giving strong convergence to the simplex orbit.
[Quadratic expansion (B1) and local regret theorem (B11)](POSITIVE_CORRELATION_CONTINUATION.md)
