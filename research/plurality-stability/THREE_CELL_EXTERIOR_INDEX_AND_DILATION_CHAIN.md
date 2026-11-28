# Exterior stability and scale-linked limits for Gaussian interfaces

For an actual planar three-cell Gaussian noise-stability maximizer,
the index computation gives exterior positivity up to a vanishing
error in the weighted space of interface speeds. Genuine exterior
coercivity follows when the critical reciprocal-gradient mass in
unit balls vanishes at infinity. That condition is also necessary
for compactness of the interface interaction. The stationary score
equations at infinity have a sequence of scale-linked local limits.
The first variation in a fixed list of these windows alone admits
periodic three-label profiles with infinitely many parallel walls.

In the local Allen–Cahn equation the implication from finite index
to exact exterior stability precedes the end argument. Wang–Wei
assume a bounded scalar solution with a double-well potential in
Theorem 2.1. Their exterior translation limits are entire stable
solutions of the same autonomous equation. The scalar planar
classification yields transition-gradient separation in Lemma 3.2.
Their curvature estimates and the nodal-domain argument for
directional derivatives in Proposition 5.1 then control the ends.
For a local second variation, separating the supports also removes
their mutual interaction. The following Gaussian statements make
the corresponding weighted-space and translation-limit issues
explicit.
[Kelei Wang and Juncheng Wei, Finite Morse Index Implies Finite Ends, Theorem 2.1, Lemma 3.2 and Proposition 5.1](https://arxiv.org/pdf/1705.06831)

## Exterior stability up to a vanishing error

Write \(r\in(0,1)\), \(\sigma^2=1-r^2\) and orient each pair
interface \(e=(i,j)\) from \(i\) to \(j\), with label vector
\(d_e=e_i-e_j\). Let \(d\nu_e=\phi_2\,d\mathcal H^1|_{\Sigma_e}\)
and \(a_e>0\) be the magnitude of the pair-score gradient on the
regular interface. For brevity, integrate over the tagged union of
the three pair interfaces with measure \(\nu\). The weighted norm
is \(\|f\|_a^2=\int a f^2\,d\nu\), and the ordinary Hessian is
\(Q_0(f)=\|f\|_a^2-\langle f,Kf\rangle\) with
\(K(e,x;e',y)=(d_e\cdot d_{e'})M_r(x,y)\).

Let \(V\) be the span of the two translation speeds and \(P_V\)
its weighted orthogonal projection. The translation equation (G2)
and index decomposition (G17) give
\[
 Q_0(f)\ge -\frac{1-r}{r}\,\|P_V f\|_a^2 .                       \tag{EI1}
\]
For an exterior speed supported on \(\{|x|>R\}\), duality in
\(V\) gives \(\|P_V f\|_a\le\|f\|_a\sup_{\|v\|_a=1,v\in V}
\|1_{|x|>R}v\|_a\). This supremum is at most
\(C\exp(-cR^2)\), uniformly on compact interior bands of
correlation for balanced maximizers. Indeed translation speeds
are normal components of constant planar vectors, the weighted
translation Gram matrix is uniformly positive by compactness
and (G24), and \(a\) is bounded above by the conditional Gaussian
derivative estimate. The polynomial interface-length bound
then gives the claimed Gaussian tail. Thus exterior speeds have
\(Q_0(f)\ge-Ce^{-cR^2}\|f\|_a^2\). No volume condition is imposed
on these speeds.
[Translation index and compactness of its weighted Gram matrix, G17 and G24](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)
[Polynomial interface-length bound](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#polynomial-growth-of-the-network)

## The critical measure for compactness

Put \(\kappa=(1-r)/(2(1+r))\) and set
\[
 d\lambda_e(x)=e^{r|x|^2/(1+r)}\frac{d\nu_e(x)}{a_e(x)} .    \tag{EI2}
\]
This measure has density \((2\pi)^{-1}e^{-\kappa|x|^2}/a_e(x)\)
against arclength. The isometry
\(g_e(x)=e^{-r|x|^2/(2(1+r))}a_e(x)f_e(x)\)
from \(L^2(a\,d\nu)\) to \(L^2(d\lambda)\) turns the interaction
into the bounded positive operator \(\mathcal K\) with kernel
\[
 \mathcal K(e,x;e',y)=\sigma^{-2}(d_e\cdot d_{e'})
             \exp\left[-\frac{r|x-y|^2}{2\sigma^2}\right].  \tag{EI3}
\]
This follows from the identity
\(M_r(x,y)=\sigma^{-2}e^{r(|x|^2+|y|^2)/(2(1+r))}
e^{-r|x-y|^2/(2\sigma^2)}\). In particular, the spatial kernel
after conjugation is translation invariant, and the relevant
weights are the reciprocal-gradient measures (EI2).

The total measure \(\lambda\) has uniformly bounded mass in unit
balls. For this, apply (G1) after conjugation to the indicator of
a unit ball on one pair. The mass is finite by (G20). On its
Cartesian square the same-pair kernel in (EI3) has a positive
lower bound independent of the ball center. Thus \(\lambda_e(B)^2\)
times a positive constant is at most \(\lambda_e(B)\).
Summing over the three pairs proves the bound.

The interaction \(\mathcal K\) is compact if and only if
\[
               \lim_{R\to\infty}\sup_{|z|>R}\lambda(B_1(z))=0 . \tag{EI4}
\]
For sufficiency, cover the plane by a fixed unit mesh. With
both coordinates restricted to \(|x|,|y|>R\), the row integrals
of the absolute kernel in (EI3) are bounded by a constant
times the largest exterior unit-ball mass, with a fixed shift
in \(R\). Indeed the Gaussian supremum on mesh cells summed
over all cells is uniformly bounded. The Schur estimate gives
operator-norm decay of the exterior compression. The interior
compression is Hilbert–Schmidt, by finite measure there and
bounded kernel. For mixed interior and exterior terms use
Cauchy–Schwarz for the positive interaction form and the
decay of the exterior compression. This makes the whole
operator a norm limit of compact interior compressions. The
same proof shows \(Q_0(f)\ge(1-o(1))\|f\|_a^2\) on exterior
supports if (EI4) holds.

For necessity, if (EI4) fails take unit balls with escaping
centers and with a common lower bound on the mass of one
pair. Their indicators normalized in \(L^2(d\lambda_e)\) form
a weakly null sequence while the same-pair kernel lower bound
makes their interaction stay positive. This contradicts
compactness. One sufficient geometric hypothesis for (EI4)
is \(a_e(x)\ge c\exp[-(\kappa-\eta)|x|^2]\) on all far interfaces
for some \(\eta>0\), by the polynomial length bound.
[First extra-coordinate bound, G1, and reciprocal-gradient local integrability, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

## The dilation chain at infinity

Let \(p_n\) escape to infinity in the fixed active plane and
translate around all the sites of its contraction orbit. For each
integer \(k\) put \(F_{n,k,i}(u)=1_{A_i}(r^k p_n+u)\). A diagonal
weak-* subsequence in \(L^\infty\) gives \(F_{k,i}\), where
\(0\le F_{k,i}\le1\) and \(\sum_iF_{k,i}=1\) almost everywhere.
With \(G_\sigma\) the Gaussian density of variance \(\sigma^2 I\)
and \(\beta_i\) the mass prices, the translated adjusted scores
are exactly
\[
 U_{n,k,i}(u)=(G_\sigma*F_{n,k+1,i})(ru)-\beta_i .          \tag{EI5}
\]
The reason is that the conditional source mean at \(r^k p_n+u\)
is \(r^{k+1}p_n+ru\). Every Gaussian density derivative is
integrable and the translated inputs are uniformly bounded.
Weak-* convergence gives pointwise convergence in (EI5) and
of every fixed derivative. Uniform derivative bounds and
equicontinuity give convergence in \(C^m_{\rm loc}\) for every
finite \(m\). The score winner rule therefore passes to the
limit as
\[
 F_{k,i}(u)>0\quad\Longrightarrow\quad
 U_{k,i}(u)=\max_j U_{k,j}(u)\quad\text{a.e.},\qquad
 U_{k,i}=(G_\sigma*F_{k+1,i})(r\,\cdot)-\beta_i .          \tag{EI6}
\]
To see this, on each compact subset of a strict losing region
the uniform score gap eliminates that label for all sufficiently
large \(n\). Testing against functions supported there gives its
weak-* limit zero. If all adjusted pair contrasts in chart \(k\)
are nonidentically zero, their analytic zero sets are null, so
the limit in that chart is ordinary and the label convergence
is strong in \(L^1_{\rm loc}\).

There is a joint limit statement when (EI4) fails. Take
\(p_n\) so the translated reciprocal-gradient measures
\(\lambda_{n,e}(du)=\lambda_e(d(p_n+u))\) have a nonzero vague
subsequential limit \(\lambda_{*,e}\). Existence uses the uniform
unit-ball mass bound, and nonzero mass follows by testing with
a function equal to one on the unit ball of each chart with mass
bounded below. Extract the score chain on the same subsequence.
The limiting measure on pair \(e=(i,j)\) is supported in
\[
 \{u: U_{0,i}(u)=U_{0,j}(u)=\max_\ell U_{0,\ell}(u),
                 \quad \nabla(U_{0,i}-U_{0,j})(u)=0\}.      \tag{EI7}
\]
A strict score gap in the limit excludes that pair on a
neighborhood for large \(n\). If its limiting gradient is
nonzero, actual \(a_e\) is bounded below on the corresponding
interface portion by (EI5). Formula (EI2), the Gaussian radial
factor and the polynomial arclength bound make the translated
mass of that neighborhood tend to zero. This proves (EI7).

The limiting measures also satisfy, for every continuous compactly
supported collection of speeds \(g_e\),
\[
 \sum_e\int g_e^2\,d\lambda_{*,e}\ \ge\
 \sum_{e,e'}\iint \mathcal K(e,u;e',v)g_e(u)g_{e'}(v)
                 \,d\lambda_{*,e}(u)d\lambda_{*,e'}(v).     \tag{EI8}
\]
Pull each test back to the translated window and use (EI1).
Its negative error vanishes because its support escapes to
infinity and the local masses of the measures are bounded.
Compactly supported continuous kernel tests pass to the vague
limit, including the product integrals, by approximation with
sums of separated functions. This proves (EI8). In particular,
if all such score-chain limits lack critical maximal pair ties in
the observed chart, then the weighted interaction is compact.
This is a sufficient nondegeneracy condition, and no such
classification of the chain limits is assumed here.

## Nondegenerate critical points and the critical measure

Let \(h_n=U_{n,0,i}-U_{n,0,j}\) and \(h=U_{0,i}-U_{0,j}\).
Center the observed chart at a critical maximal pair tie by a fixed
translation, so \(h(0)=0\) and \(Dh(0)=0\). If the Hessian of \(h\)
there is definite, the point cannot belong to
\(\operatorname{supp}\lambda_{*,(i,j)}\).

Change the sign of \(h_n\) if necessary. Smooth convergence gives
\(D^2 h_n\ge mI\) on a small fixed ball, for some \(m>0\), and a
unique nearby critical point \(q_n\to0\). On a ray from \(q_n\)
inside the ball the outward derivative at radius \(t>0\) is
at least \(mt\). There is at most one positive radial root
of \(h_n=0\). In radial coordinates its regular reciprocal-gradient
arclength is
\[
 \frac{d\mathcal H^1}{|Dh_n|}
   =\frac{t}{|\partial_t h_n|}\,d\theta
   \le \frac{1}{m}\,d\theta .                                      \tag{EI9}
\]
The same bound applies on the active subset of the level curve.
The critical measure of a sufficiently small fixed ball is at
most \(C e^{-\kappa(|p_n|-C')^2}\) by (EI2), so it tends to zero.

At a nondegenerate saddle suppose the third score is strictly
lower at the limiting tie. It stays strictly below the two pair
scores on a neighborhood for all large \(n\). Let \(q_n\to0\)
be the nearby critical point and \(\delta_n=h_n(q_n)\). Uniformly
bi-Lipschitz local coordinates \((x,y)\) give
\(h_n=\delta_n+xy\). To construct them, make a fixed linear change
taking the limiting Hessian to signed diagonal form. Taylor's
formula at \(q_n\) writes \(h_n-\delta_n=z^{\mathsf T}H_n(z)z\).
The matrix \(H_n\) is symmetric, has entries uniformly close
to those of a fixed diagonal matrix with one positive and one
negative entry on small balls, and has uniformly bounded first
derivatives. Complete the square in the first variable and take
the positive square roots of the resulting absolute coefficients.
The resulting change to \(u^2-v^2\) has derivative uniformly
close to an invertible constant on a sufficiently small ball.
Set \(x=u+v\), \(y=u-v\).

If \(\delta_n=0\), all regular points of the crossing axes in
this neighborhood belong to the actual pair interface. The
integral of \(1/|Dh_n|\) on either axis diverges at the crossing,
contrary to (G20). Thus \(\delta_n\ne0\). The coarea density
of \(xy=-\delta_n\) is \(dx/|x|\), and in a fixed coordinate
square \(|x|\) ranges from \(|\delta_n|\) divided by the side
bound to the side bound. The coordinate change multiplies the
coarea density by its bounded inverse Jacobian determinant.
Formula (EI2) gives on a sufficiently small fixed chart ball \(B\)
\[
 \lambda_{n,(i,j)}(B)
 \le C e^{-\kappa(|p_n|-C')^2}
                    \bigl(1+\log^+(1/|\delta_n|)\bigr).             \tag{EI10}
\]
If the saddle belongs to the support of the limiting pair measure,
a nonnegative continuous test supported in \(B\) gives a positive
lower bound on the left. For some \(c>0\) and all large \(n\),
\[
                 |\delta_n|\le
            \exp\!\left[-c\,e^{\kappa(|p_n|-C')^2}\right].
\]
This necessary condition at a strictly two-winner saddle leaves
degenerate critical points and saddles at three-way ties open.

## A periodic chain with three labels

The exact rule (EI6), in isolation, does not bound the number of
walls in a chart. Index the label centers on a circle of period
\(2\pi\) by \(0,2\pi/3,4\pi/3\) and use their nearest-center
intervals. In chart \(k\), label \(u\) by its nearest center at
the angular coordinate \(c r^{-k}u_1\bmod 2\pi\), with \(c>0\).
Take all the prices equal. Gaussian smoothing of any of these
interval indicators on the circle is an even function of the
distance to its center, strictly decreasing between distances
\(0\) and \(\pi\). For completeness, the wrapped heat kernel
has this same monotonicity. Approximate its point source at zero
by normalized \(\exp(C\cos\theta)\), \(C\to\infty\). The derivative
of each heat evolution is nonpositive on \((0,\pi)\) by the
Dirichlet heat maximum principle and oddness about the endpoints.
The limit derivative has this sign and is strictly negative for
positive time by the strong maximum principle and the nonvanishing
first harmonic. Convolution of the wrapped kernel with an interval
of half-width \(\pi/3\) preserves strict radial monotonicity:
its derivative is the difference of the kernel at the two interval
endpoints, and the upper endpoint's circular distance is larger
for angles strictly between \(0\) and \(\pi\).

The convolved scores of the three intervals therefore still select
the nearest center. After smoothing chart \(k+1\), evaluation
at \(ru\) gives exactly the angular coordinate of chart \(k\).
Thus all these profiles solve (EI6), with infinitely many parallel
walls of the three pair types in each complete chart. Their
winning pair gradients are nonzero on the walls for each fixed
\(k\) by the same strict monotonicity. This is a solution of the
limit equations and is not asserted to arise from a single actual
Gaussian maximizer. In a finite list of translated windows where
the actual gradients are uniformly bounded below, the local
\(\lambda\)-masses tend to zero by (EI2). The interaction among
these windows vanishes in the weighted norm by (EI3), so their
localized stability is automatically positive. Establishing finite
ends still requires geometry relating the windows along the
contraction orbit, beyond (EI6) and such finite-window tests.

## Fixed exterior windows and a central window of the same partition

Requiring all fixed windows to come from a single balanced
partition still leaves the central window free if the winner
equations are imposed only in the limit. To see this, take an
ordinary chain \(F_k\) solving (EI6), for example the stripes
above. Fix a unit vector \(v\) and set \(p_n=r^{-n^3}v\). Define a
single partition in the balls \(B_n(r^kp_n)\), \(|k|\le n\), by
\[
              1_{A_i}(r^k p_n+u)=F_{k,i}(u),\qquad |u|<n .
                                                               \tag{EI11}
\]
After discarding finitely many rows, all balls are disjoint.
The radial exponents \(n^3-k\) are distinct integers, in disjoint
blocks indexed by \(n\). Two distinct centers with larger exponent
\(j\) are separated by at least \((1-r)r^{-j}\), which exceeds the
sum of their radii for all sufficiently large rows. The total
Gaussian mass of the prescribed balls can be made arbitrarily
small, since the smallest radius of their centers in row \(n\)
is \(r^{-(n^3-n)}\). The balls also avoid any fixed bounded set.
Prescribe any labels on that set with each prescribed mass less
than \(1/3\), then use the nonatomic remaining Gaussian measure
to give all three labels total mass \(1/3\).

For each fixed \(k\), both chart \(k\) and its conditional source
chart \(k+1\) match the assigned profiles on balls of radius
tending to infinity. Thus all these translated profiles have
limit \(F_k\), and their exact smoothed scores converge to (EI5)
with every fixed derivative on compact sets. Gaussian derivative
tails bound the error from sources outside the assigned balls.
Yet the chart at the moving index \(k_n=n^3\) is centered at \(v\)
and has the prescribed local central geometry, independent of
the exterior chain. This construction imposes no exact winner
rule or maximality at finite \(n\). It shows why linking the
exterior chain to the central geometry requires information on
a growing part of the orbit or another global restriction from
maximality.
