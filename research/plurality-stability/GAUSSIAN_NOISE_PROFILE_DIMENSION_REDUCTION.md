# Dimension reduction for a Gaussian noise profile

For a Gaussian partition, consider its noise stabilities at several positive
correlations as a vector. Maximizing a smooth increasing function of this
vector is an exact finite-dimensional problem. A linear objective needs at
most \(k-1\) dimensions for \(k\) cells. For a nonlinear objective, the
additional dimensions are bounded by the negative index of its Hessian on
the relevant space of translation derivatives. In particular, every
maximum of the minimum stability across \(m\) fixed positive correlations
can be attained in at most \(k+m-2\) dimensions, with arbitrary fixed
thresholds subtracted from the individual stabilities.

The argument combines the bulk Hermite translation identity (SD3) with
simultaneous first-order constraints. For a single correlation the
dimension bound is the one in Heilman–Tarter, Theorem 1.9. Their Lemma 5.5
proves the corresponding positive translation form for a single
correlation using interfaces. Here the bulk identity applies to the
weighted score chosen by the differential of the objective.
[Heilman–Tarter](https://doi.org/10.1017/fms.2021.56)
[Bulk spectral translation identity](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md#the-bulk-translation-identity)

## Smooth objectives

Fix \(k\ge2\) prescribed masses \(p_i>0\) summing to one, \(m\ge1\)
correlations \(0<\rho_j<1\), and a finite ambient dimension \(d\ge1\).
Write \(T_j=T_{\rho_j}\) for the Gaussian Ornstein–Uhlenbeck operator.
For an ordered partition with indicator tuple \(f=(f_i)\), define
\[
 F(f)=(S_1(f),\ldots,S_m(f)),\qquad
 S_j(f)=\sum_i\langle f_i,T_jf_i\rangle_{\gamma_d}.
\]
Let \(\Phi\) be \(C^1\) on a neighborhood of \([0,1]^m\) with
nonnegative partial derivatives and nonzero gradient there. A maximizing
partition exists. For each such maximizer, let
\[
 M=\operatorname{span}\left\{\int x f_i(x)\,d\gamma_d(x):1\le i\le k\right\},
 \qquad U=M^\perp,\qquad \xi=\nabla\Phi(F(f)).
\]
For \(v\in U\), translate by \(f_i^{t,v}(x)=f_i(x-tv)\) and let
\[
 A(v)=\left.\frac{d}{dt}F(f^{t,v})\right|_{t=0}.
\]
This is a linear map \(U\to\xi^\perp\). The partition is cylindrical over
the orthogonal complement of \(\ker A\) in \(\mathbb R^d\). In particular
its active dimension is at most
\[
                  \dim M+\operatorname{rank}A\le k+m-2.       \tag{NP1}
\]

If \(\Phi\) is \(C^2\) at the maximizing profile, define the real
quadratic form on \(U\)
\[
 q(v)=D^2\Phi(F(f))[A(v),A(v)].
\]
Then \(D^2\Phi(F(f))\) is negative definite on \(\operatorname{im}A\)
whenever that image is nontrivial. The negative index \(\nu_-(q)\),
the maximal dimension of a subspace on which \(q\) is negative
definite, equals \(\operatorname{rank}A\). The active dimension is
at most
\[
                     \dim M+\nu_-(q).                         \tag{NP2}
\]
In particular \(\nu_-(q)\) is at most the negative index of the
Hessian restricted to \(\xi^\perp\). If the Hessian is nonnegative on
\(\operatorname{im} A\), the partition is cylindrical over the
centroid span \(M\). The same centroid conclusion holds when \(\Phi\)
is convex and \(C^1\), by the supporting-plane inequality below.

### Existence and simultaneous rounding at the optimum

Relax the indicators to simplex-valued tuples with means \(p_i\).
This set is closed, bounded and convex in the product Gaussian
\(L^2\) space, so it is weakly compact. Each \(T_j\) is compact
because its Hermite multiplier \(\rho_j^n\) tends to zero and its
degree-\(n\) Hermite subspaces are finite-dimensional for \(d<\infty\).
The profile is therefore weakly continuous, and the relaxed problem
attains its maximum.

If a relaxed maximizer is fractional, two coordinates \(f_a,f_b\)
are uniformly bounded below by a positive number on some
positive-measure set. Pick a bounded nonzero real \(h\) supported there
with
\[
 \int h\,d\gamma_d=0,\qquad
 \langle h,T_j(f_a-f_b)\rangle_{\gamma_d}=0
       \quad(1\le j\le m).
\]
To get such an \(h\), take a linear combination of the indicators of
\(m+2\) disjoint positive-measure subsets and solve the homogeneous
system of \(m+1\) linear constraints. Gaussian measure is nonatomic.
For sufficiently small positive \(\epsilon\), the tuple
\(f+\epsilon h(e_a-e_b)\) is feasible. At every correlation its
stability increases by \(2\epsilon^2\langle h,T_jh\rangle\), which is
strictly positive by the strictly positive Hermite multipliers.
The dot product of this increase with \(\xi\) is strictly positive,
so differentiability of \(\Phi\) contradicts maximality for small
\(\epsilon\). Thus all relaxed maximizers are ordinary, and every
ordinary maximizer also solves the relaxed problem.

### Prices and small exact mass repairs

Differentiate \(\Phi\) at the optimizing tuple along a feasible
segment to any other relaxed tuple of the same masses. It follows
that \(f\) maximizes the linear assignment with fixed scores
\(Tf_i\), where \(T=\sum_j\xi_jT_j\). There exist mass prices \(c_i\)
for which \(f\) chooses pointwise winners of \(s_i=Tf_i-c_i\)
almost everywhere. Indeed the optimal value of the linear assignment
as a function of its masses is finite and concave. A supergradient
\(c\) exists at the prescribed interior mass vector, and the
supporting inequality makes \(f\) an unconstrained optimizer of
the assignment with the adjusted scores.

For \(v\in U\), the translated masses differ from \(p_i\) by \(O(t^2)\).
Their first derivatives are \(-v\cdot\int x f_i\,d\gamma_d\), as seen
by translating the Gaussian density. The second derivatives of the
mass integrals are locally bounded by differentiation of that density.
Each \(S_j(f^{t,v})\) is also \(C^2\) in \(t\). For example,
\(S_j(f^{t,v})=\sum_i\|T_{\sqrt{\rho_j}}f_i^{t,v}\|_2^2\) and
\(T_{\sqrt{\rho_j}}f_i^{t,v}(x)=(T_{\sqrt{\rho_j}}f_i)(x-tv/\sqrt{\rho_j})\).
The smoothed indicator has bounded derivatives of every fixed order,
by differentiation of the Gaussian convolution. This also gives
linearity of \(A(v)\) in \(v\).

The essential closures of the finitely many cells cover \(\mathbb R^d\),
and their intersection graph is connected. Otherwise the unions over
two complementary groups of graph components would separate
\(\mathbb R^d\) into disjoint nonempty closed sets. At a point in two
essential closures the continuous scores of both cells are tied for
the largest score. Fix a spanning tree in the graph. For each oriented
tree edge \(i\to l\) and \(\eta>0\), choose a bounded source patch
within cell \(i\) of positive measure where \(s_l-s_i\ge-\eta\).
Such a patch exists by the tying point and continuity.

Send the translated mass errors along the tree and transfer from the
translated source patches to realize those amounts exactly. The tree
flow over an edge is the sum of mass errors on one side of its cut,
so the total mass transferred is \(O(t^2)\) with a constant
independent of \(\eta\). For each fixed \(\eta\) and small enough
\(|t|\), all translated patches still have positive measure bounded
away from zero. The transferred subsets can be chosen disjoint in
each source cell because the total mass requested tends to zero.
Denote the repaired tuple by \(f^{t,v,\eta}\). Every individual
stability changes by \(O(t^2)\) under this repair. Indeed a label
change on a set of Gaussian measure \(\delta\) changes the agreement
event only when one of its two Gaussian arguments lies in that set,
of joint probability at most \(2\delta\).

For the weighted Lagrangian
\[
 \mathcal L(u)=\sum_i\langle u_i,Tu_i\rangle
                         -2\sum_i c_i\int u_i\,d\gamma_d
\]
the lower error bound of the repair is
\[
 \mathcal L(f^{t,v,\eta})-\mathcal L(f^{t,v})
                    \ge -C\eta t^2+o_\eta(t^2),               \tag{NP3}
\]
where \(C\) does not depend on \(\eta\). To check it, expand the
quadratic Lagrangian in the repair change. Each transferred mass
element contributes twice the adjusted weighted score difference
to the linear term. These score differences on the translated
patches are bounded below by \(-\eta-o_\eta(1)\), because the scores
of the translated tuple converge to \(s_i\) locally uniformly.
This follows also from the displayed translation formula for
Gaussian smoothing. The quadratic term in the change is
nonnegative since \(T\) is a nonnegative mixture of positive
operators. The constant \(C\) in (NP3) comes only from the
tree mass bound.

### Translation directions that force cylindrical cells

For any \(v\in U\),
\[
 F(f^{t,v,\eta})-F(f)=tA(v)+O_\eta(t^2).
\]
Maximality for both signs of \(t\), at any fixed \(\eta\), gives
\(\xi\cdot A(v)=0\). This proves the asserted image condition.
If \(v\in\ker A\), the entire repaired profile change is \(O_\eta(t^2)\).
Differentiability of \(\Phi\) at the optimum then gives
\[
            \xi\cdot(F(f^{t,v,\eta})-F(f))\le o_\eta(t^2).
\]
When \(\Phi\) is \(C^2\), Taylor expansion for any \(v\in U\) instead gives
\[
  \xi\cdot(F(f^{t,v,\eta})-F(f))
                   +\tfrac12t^2q(v)\le o_\eta(t^2).
\]
The repaired masses equal the original masses, so the linear
profile change is exactly
\(\mathcal L(f^{t,v,\eta})-\mathcal L(f)\).
Using (NP3), dividing by \(t^2\), then letting \(t\to0\) and
\(\eta\downarrow0\) bounds the second derivative of
\(\mathcal L(f^{t,v})\) above by zero for \(v\in\ker A\), and
above by \(-q(v)\) in the \(C^2\) case. Its first derivative is
zero since \(\xi\cdot A(v)=0\) and \(v\in U\).

Apply the bulk identity (SD3) to the same winner scores. On
Hermite degree \(n\) the multiplier of \(T\) is
\(\lambda_n=\sum_j\xi_j\rho_j^n\). The degree-zero multiplier may
be removed from \(T\) while replacing \(c_i\) by \(c_i-\lambda_0p_i\).
This changes the Lagrangian on the translated tuple, up to a
constant, by \(\lambda_0\sum_i(\int f_i^{t,v}d\gamma_d-p_i)^2=O(t^4)\).
In the identity the \(n=0\) term vanishes since
the degree-zero Gaussian coefficient of the weak derivative
\(D_v f_i\) is \(v\cdot\int x f_i\,d\gamma_d=0\).
Thus
\[
 \frac12\left.\frac{d^2}{dt^2}\mathcal L(f^{t,v})\right|_{t=0}
   =\sum_i\sum_{n\ge1}(\lambda_n-\lambda_{n+1})
                            \|P_nD_v f_i\|_2^2 .              \tag{NP4}
\]
The bulk identity needs the pointwise winner condition and
two derivatives of the translated Lagrangian, established above.
It therefore applies even though \(\Phi\) need not be convex and
\(f\) need not maximize the quadratic energy of \(T\) globally.
The multipliers decrease strictly in positive degree since
\(\xi_j\ge0\) and some \(\xi_j>0\). Whenever the second derivative
is nonpositive, every Hermite coefficient of \(D_v f_i\) vanishes.
In coordinates with first axis \(v/|v|\) for nonzero \(v\), the
Gaussian adjoint derivative raises that coordinate's Hermite degree.
The vanishing coefficients therefore say that \(f_i\) has no positive
Hermite degree in that coordinate, so each cell is cylindrical
along \(v\).

For (NP1) take a basis of \(\ker A\), whose codimension in
\(\mathbb R^d\) is \(\dim M+\operatorname{rank}A\). The separate
cylindrical conclusions give dependence only on the projection
to the remaining orthogonal subspace, as is seen in a common
orthonormal Hermite basis. The vectors generating \(M\) sum to
zero, giving \(\dim M\le k-1\). For (NP2), any \(v\) with
\(q(v)\ge0\) is a direction of invariance, so \(A(v)=0\).
Every nonzero vector in \(\operatorname{im} A\) accordingly has
strictly negative Hessian quadratic form. On the orthogonal
complement of \(\ker A\) in \(U\), \(q\) is negative definite.
This proves the claimed equality of index and rank.

For convex \(C^1\) objectives the supporting-plane inequality gives
\(\xi\cdot(F(u)-F(f))\le0\) for every feasible tuple \(u\), by
maximality of \(\Phi(F(f))\). Thus \(f\) also maximizes the weighted
quadratic objective and the same translation argument works on
all of \(U\). This proves the centroid conclusion for this case.

## A minimum across correlations

Fix real thresholds \(b_j\). The maximum of
\(\min_j(S_j(f)-b_j)\) in finite ambient dimension \(d\) is attained
by a partition using at most \(\min(d,k+m-2)\) dimensions.
For \(\beta>0\) use the smooth objective
\[
 \Phi_\beta(s)=-\beta^{-1}\log
                   \sum_{j=1}^m\exp[-\beta(s_j-b_j)] .
\]
Its partial derivatives are positive and it converges uniformly to
the minimum, with error at most \(\log(m)/\beta\). Hence by (NP1) the
optimal maximum values in dimensions \(d\) and \(\min(d,k+m-2)\)
agree in the limit \(\beta\to\infty\). The exact minimum attains
its maximum on the relaxed compact class in that smaller
dimension. The simultaneous split above strictly increases all
criteria at a fractional optimizer, which strictly increases
their minimum. Thus a maximizing partition in the smaller
dimension is ordinary. Extension cylindrically preserves masses
and every stability. This proves an attainment bound for the
nonsmooth objective.

**Simultaneous domination.** For any partition \(f\) in any finite
dimension, there is a partition \(g\) with the same masses in at most
\(\min(d,k+m-2)\) dimensions such that
\[
                         S_j(g)\ge S_j(f)\quad(1\le j\le m).
                                                               \tag{NP5}
\]
Apply the minimum theorem with \(b_j=S_j(f)\). Its optimal value
is at least zero because \(f\) is admissible. The reduced maximizer
therefore has the claimed inequalities. In particular, every
coordinatewise nondecreasing objective of the finite noise profile
has the same supremum in all dimensions at least \(k+m-2\), with
attainment in that dimension whenever the supremum is attained
in any finite dimension.

For a compact interval \(I\subset(0,1)\) and continuous thresholds
\(b(\rho)\), a finite grid gives a dimension bound for a partition
that is within any given tolerance of the best value of
\(\inf_{\rho\in I}(S_\rho(f)-b(\rho))\) over all finite dimensions.
Indeed \(\|T_\rho-T_\sigma\|\le C_I|\rho-\sigma|\) by the Hermite
multipliers, with \(C_I=\sup_{n\ge1}n(\max I)^{n-1}<\infty\).
Therefore the grid error is controlled uniformly by \(C_I\) and
the modulus of continuity of the thresholds. Apply the exact
finite-grid result in at most \(k+m-2\) dimensions for a grid of
size \(m\). A positive value for thresholds \(b_j=S_{\rho_j}(T)\)
witnesses one partition that improves on the simplex at all
the selected correlations together.
On the same grid, (NP5) gives a partition whose entire profile on
\(I\) is within an arbitrarily small one-sided error of dominating
the profile of any prescribed input partition. At a radius between
grid points, apply the operator-norm continuity bound once to each
partition and use their ordering at the closest grid point.

## A two-correlation test at a simplex contact

Let \(r\in(0,1)\) be a correlation at which the balanced simplex \(T\)
is globally optimal. For \(h>0\) with \(r+h<1\) and
\(\varepsilon\ge0\), define
\[
 J(h,\varepsilon)=
   \max\{S_{r+h}(A): A\text{ balanced},\
                         S_r(A)\ge S_r(T)-\varepsilon\},
\]
where one may allow all ambient dimensions at least \(k-1\). By (NP5)
the maximum is unchanged by restricting to dimension \(k\). In that
dimension the relaxed constraint is weakly closed and the objective is
weakly continuous. A fractional maximizer could be split to increase
both stabilities, so there is an ordinary maximizer. Compactness gives
\(\lim_{\varepsilon\downarrow0}J(h,\varepsilon)=J(h,0)\): weak limits
remain feasible at the limiting constraint. For \(h\downarrow0\), the
Hermite multiplier Taylor expansion is uniform over partitions, since
the multipliers have uniformly bounded second derivatives on compact
subintervals of \((0,1)\). It gives
\[
 \lim_{h\downarrow0}\frac{J(h,0)-S_r(T)}h
     =\max_{A:\,S_r(A)=S_r(T)} \partial_r S_r(A).
                                                               \tag{NP6}
\]
For the equality, expand \(S_{r+h}\) on the fixed tie set and use the
uniform remainder. The derivative maximum exists by compactness of
that set in the relaxed class in dimension \(k\); its members are
ordinary by the positive single-correlation split. The original
single-correlation dimension reduction also places every global
maximizer in at most \(k-1\) dimensions. The constrained
two-correlation problem gives a fixed-dimensional test for a tie
whose outgoing slope exceeds the simplex slope.

## A continuum of correlations with finite negative index

Let \(I\subset(0,1)\) be compact, and use
\(F(f)(\rho)=S_\rho(f)\) as a profile in \(C(I)\) with the uniform norm.
Suppose \(\Phi\) is twice continuously Fréchet differentiable on a
neighborhood of the profiles with values in \([0,1]\). Assume its
differential at each such profile is a nonzero positive functional.
If the Hessian at a maximizing profile has negative index at most
\(q\) on the kernel of its differential, that maximizing partition
uses at most \(k-1+q\) dimensions. The bound can also use the index on
the image of centroid-orthogonal translations, as in (NP2). A uniform
index bound gives the same dimension bound in every finite ambient
dimension.

Existence follows from the weak compactness argument above. The profile
map is continuous into \(C(I)\) for weak convergence in the relaxed
class: a finite grid and the uniform bound \(C_I\) reduce this to weak
continuity at finitely many correlations. By the Riesz representation
theorem the differential at a maximizer is integration against a finite
nonnegative measure \(\mu\) on \(I\), with \(\mu(I)>0\). The weighted
assignment is thus given by \(T=\int_I T_\rho\,d\mu(\rho)\).

The weighted scores are real analytic. For a bounded function \(h\),
the translated Gaussian kernel formula for \(T_\rho h\) extends to
complex input \(x\). The kernel's absolute value is the real Gaussian
density centered at \(\rho\operatorname{Re} x\), times
\(\exp(\rho^2|\operatorname{Im}x|^2/(2(1-\rho^2)))\). Compactness of
\(I\) supplies locally uniform domination for integration in both the
Gaussian variable and \(\rho\). In particular the weighted scores are
analytic after subtracting mass prices.

If two coordinates of a relaxed maximizer are both uniformly positive
on a positive-measure set, their adjusted assignment scores agree
there. Analyticity makes them agree everywhere. The positive-degree
multipliers of \(T\) are strictly positive, so the two coordinates
differ by a constant. A bounded mean-zero perturbation \(h\) between
this pair, supported on the positive-measure set, then has zero linear
term in every \(S_\rho\). Its quadratic increase is strictly positive
at each \(\rho\in I\), and integrating that increase against \(\mu\)
is strictly positive. Differentiability rules out fractional
maximizers. This obtains the simultaneous cancellation across \(I\)
from the winner rule.

Finally the profile translation has a uniform second-order expansion
on \(I\), using the smoothed-translation formula and compactness of
\(I\). The same tree repair changes the profile by \(O(t^2)\) in the
uniform norm, and (NP3) applies to the mixture. Thus the Taylor
argument and the bulk identity (NP4) apply with \(A:U\to C(I)\) and
\(\lambda_n=\int_I \rho^n\,d\mu(\rho)\). This multiplier is exponentially
bounded and strictly decreasing in positive degree. The Hessian must
be negative definite on the image of \(A\), so
\(\operatorname{rank} A\le q\). The proof gives exact reduction for
smooth objectives of finite negative index, while the grid estimate
above covers the nonsmooth continuum minimum to prescribed accuracy.

These are written dimension reduction results with an internal
source audit. A comparison with the simplex on an unresolved
correlation interval still requires control of the optimizing
partition in the reduced dimension.
