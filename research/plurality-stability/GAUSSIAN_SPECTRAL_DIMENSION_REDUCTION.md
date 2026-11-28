# Dimension reduction for Gaussian spectral energies

Let \(P_j\) be the \(j\)-th Gaussian Hermite projection on
\(L^2(\gamma_d)\). For a real sequence \((\lambda_j)_{j\ge0}\)
define the self-adjoint multiplier \(T=\sum_j\lambda_jP_j\) and
\[
 \mathcal S_T(f)=\sum_{i=1}^k\langle f_i,Tf_i\rangle_{\gamma_d}.
\]
The competitors are ordinary measurable partitions with prescribed
positive masses \(p_1,\dots,p_k\). Addition of \(\lambda_0P_0\)
affects the objective by a mass-dependent constant, so we use
\(\lambda_0=0\) when convenient.

Suppose \(\lambda_j>0\) for \(j\ge1\) and
\(\lambda_j\le C b^j\) for some \(b<1\). If
\(\lambda_j\ge\lambda_{j+1}\) for \(j\ge1\), with strict
inequality except at finitely many \(j\), every maximizer
is cylindrical over the span of its cell centroids
\(\int x f_i(x)d\gamma_d(x)\). In particular its active
dimension is at most \(k-1\). More generally for a bounded
nonnegative nonincreasing sequence on \(j\ge1\), the optimal
value is attained by a partition using at most \(k-1\)
dimensions in every ambient dimension \(d\ge k-1\).

The case \(\lambda_j=\rho^j\), \(0<\rho<1\), is the
dimension reduction of
[Heilman–Tarter, Theorem 1.9 and Lemma 5.5](https://arxiv.org/pdf/2011.05583).
Their translation identity at the interfaces motivates the
bulk identity below. The bulk proof uses neither interface
regularity nor an individual boundary normal for each part
of the multiplier. A positive correlation mixture supported
in \([0,b]\), \(b<1\), with positive mass in \((0,b]\), fits
the first statement. There are also examples outside
nonnegative correlation mixtures. For \(0<a<1/2\),
\(\lambda_j=(j+1)a^j\) satisfies the hypotheses. Its positive
degree subsequence is strictly log-concave, whereas moments
of a nonnegative measure on \([0,1]\) are log-convex by
Cauchy–Schwarz.

## Optimization and volume repair

The relaxed admissible class of simplex-valued \(f\) with
the same means is weakly compact in \(L^2(\gamma_d)^k\).
The multiplier is compact since each degree is finite
dimensional and \(\lambda_j\to0\). Hence the quadratic
objective is weakly continuous and attains its maximum.
It is strictly convex on the fixed-mean class since all
positive-degree multipliers are positive. A maximizer
must therefore be extreme. Any nonordinary competitor
admits a nontrivial balanced split: on a positive-measure
set two coordinates have values at least some \(\delta>0\),
so opposite changes in these coordinates on two equal-measure
subsets give the split. Thus every relaxed maximizer is an
ordinary partition, and every ordinary maximizer is also
a relaxed maximizer.

Linearizing along feasible convex segments shows that a
maximizer \(f\) also maximizes the linear objective with scores
\(Tf_i\) among competitors of its masses. For those fixed
scores the maximum linear objective, regarded as a function
of the mass vector, is finite and concave on the probability
simplex. It has a supporting supergradient \(c\) at the
prescribed positive mass vector. Consequently \(f\) globally
maximizes the unrestricted linear objective with scores
\(s_i=Tf_i-c_i\). It chooses pointwise winners of the scores.
This argument needs only the elementary supporting-plane
theorem for a finite concave function at an interior point.

The scores are smooth. To see the needed estimates choose
\(\sqrt b<s<1\) and factor \(T=T_s B T_s\), where \(T_s\) is
Ornstein–Uhlenbeck at correlation \(s\) and \(B\) has
bounded multipliers \(\lambda_j/s^{2j}\). The Gaussian kernel
of \(T_s\) and each of its spatial derivatives are continuous
as \(L^2(\gamma_d)\)-valued functions of the input point,
bounded uniformly there on compact sets. This follows
directly from the conditional Gaussian formula and its
derivatives. In particular \(T\) maps strongly convergent
\(L^2\) functions to scores converging uniformly on each
compact set.

Let \(A_i\) be the cells and \(C_i\) their essential closures:
\(x\in C_i\) if every ball about \(x\) meets \(A_i\) in
positive measure. The finite collection of closed sets
\(C_i\) covers \(\mathbb R^d\), since the cells partition
it up to null sets. Its intersection graph is connected.
Otherwise its unions over two complementary groups of
components would disconnect \(\mathbb R^d\). If
\(x\in C_i\cap C_j\), continuity and the winner condition
give \(s_i(x)=s_j(x)=\max_l s_l(x)\). For each direction
\(i\to j\) of an edge in a fixed spanning tree and each
\(\epsilon>0\), there is a bounded positive-measure source
patch in cell \(i\) where \(s_j-s_i\ge-\epsilon\).

Fix a vector \(v\) perpendicular to all cell centroids and
translate the partition by \(t v\), writing \(f^t_i=f_i(\cdot-tv)\).
The mass derivative of cell \(i\) is
\(-\int f_i(x)(v\cdot x)d\gamma_d(x)=0\); differentiating
the translated Gaussian density gives mass errors \(O(t^2)\).
They can be repaired by routing their signed imbalances
along the tree, moving \(O(t^2)\) total mass. For each edge
transfer use a subset of its translated source patch of the
exact required mass. These patches remain uniformly positive
in measure for \(t\) small. The subsets can be chosen disjoint
because there are finitely many requests and their total mass
tends to zero. Finite-dimensional nonatomic Gaussian measure
allows every requested fraction.

Write \(\mathcal L(f)=\mathcal S_T(f)-2\sum_i c_i\int f_i d\gamma_d\).
The quadratic interaction of all corrections is nonnegative,
and the linear change is at least
\(-C_1\epsilon t^2+o_\epsilon(t^2)\). Indeed \(Tf_i^t-c_i\)
converges to \(s_i\) uniformly on the fixed compact set
of patch translates, and the total transferred mass has a
bound independent of \(\epsilon\). The repaired partition
has the original masses. Thus
\[
 \mathcal L(f^t)-\mathcal L(f)\le C_1\epsilon t^2+o_\epsilon(t^2).
 \tag{SD1}
\]
In particular \((d^2/dt^2)\mathcal L(f^t)|_{t=0}\le0\)
after taking the limit and then letting \(\epsilon\downarrow0\).

## The bulk translation identity

Let \(D_v=\partial_v\) and let \(D_v^*=v\cdot x-D_v\) be its
Gaussian adjoint. Write \(g_{ij}=P_j D_v f_i\), interpreted
against Hermite polynomials by integration by parts. The
raising and lowering relations give
\(\|g_{ij}\|_2^2\le |v|^2(j+1)\|P_{j+1}f_i\|_2^2\).
Hence its exponentially weighted square sums converge.

The function \(M(x)=\max_i s_i(x)\) is locally Lipschitz and
\(M=\sum_i f_i s_i\) almost everywhere. At differentiability
points of \(M\), the gradient of each active smooth score
equals that of \(M\) by first-order contact. Differentiating
this equality in distributions gives
\(\sum_i s_i D_v f_i=0\). In particular
\[
 \sum_i\langle D_v f_i,(v\cdot x)s_i\rangle_{\gamma_d}=0 .
 \tag{SD2}
\]
One can first test the distributional identity against
\((v\cdot x)\gamma_d(x)\) times a smooth cutoff and then remove
the cutoff. Indeed \(D_v^*((v\cdot x)s_i)\) and
\((D_v^*)^2s_i\) are in \(L^2(\gamma_d)\) by the raising and
lowering bounds and exponential smoothing. First-order
cutoff terms vanish by the same bounds and boundedness of
\(f_i\).

To justify differentiating the objective, use the factorization
\(T=T_s B T_s\) above. For bounded \(f\),
\(T_s f^t(x)=(T_s f)(x-tv/s)\). All fixed-order spatial
derivatives of \(T_s f\) are bounded in \(x\), by
differentiation of the translated Gaussian convolution.
Thus \(T_s f^t\) is twice differentiable in \(L^2(\gamma_d)\).
The mass terms are twice differentiable by translation of the
Gaussian density. Using the factorization, the translations
have first derivative \(-D_v f_i\) and second derivative
\(D_v^2 f_i\) in the resulting smoothed pairings. Hence
\[
 \frac12\frac{d^2}{dt^2}\mathcal L(f^t)\bigg|_{t=0}
   =\sum_{i,j}\lambda_j\|g_{ij}\|_2^2
     +\sum_i\langle D_v^2 f_i,s_i\rangle_{\gamma_d}.
\]
Gaussian integration by parts, (SD2) and Hermite lowering give
\[
 \sum_i\langle D_v^2 f_i,s_i\rangle_{\gamma_d}
   =-\sum_i\langle D_v f_i,D_v Tf_i\rangle_{\gamma_d}
   =-\sum_{i,j}\lambda_{j+1}\|g_{ij}\|_2^2 .
\]
We obtain the bulk identity
\[
 \frac12\frac{d^2}{dt^2}\mathcal L(f^t)\bigg|_{t=0}
       =\sum_{i,j}(\lambda_j-\lambda_{j+1})\|g_{ij}\|_2^2 .
 \tag{SD3}
\]
In the single-correlation case its right side is
\((1-\rho)\sum_i\|T_{\sqrt\rho} D_v f_i\|_2^2\).

Our choice of \(v\) gives \(g_{i0}=0\). By (SD1) all the
nonnegative terms in positive degrees in (SD3) vanish.
Under eventual strict decrease, each \(D_v f_i\) therefore has
only finitely many possible nonzero Hermite components.
It is a polynomial as a distribution: apply \(T_s\) and let
\(s\uparrow1\) in the identity \(D_v T_s f_i=sT_sD_v f_i\),
using local \(L^2\) convergence for \(f_i\).
The directional derivative of a bounded function cannot be
a nonzero polynomial. In coordinates running along \(v\),
its primitive on almost every line is a polynomial plus
a constant and is bounded, so it is constant. Hence \(D_v f_i=0\).
Taking a basis of the common centroid complement proves the
cylindrical conclusion. The centroids sum to zero, so their
span has dimension at most \(k-1\).

For the optimal-value extension to a bounded nonnegative
nonincreasing sequence, subtract its limit from every degree.
This subtracts a constant on the ordinary admissible class,
since \(\sum_i\|f_i\|_2^2=1\). Modify degree zero separately
if necessary. The resulting positive-degree sequence
\(a_j\ge0\) tends to zero. The multipliers
\(\lambda_j^{(s)}=s^j(a_j+1-s)\) for \(j\ge1\), \(s<1\),
satisfy the strict exponentially smoothing hypotheses and
converge uniformly in \(j\) to \(a_j\). The optimal values
in dimensions \(d\) and \(k-1\) coincide for every
approximant, hence also in the limit by uniform convergence
of the objectives on ordinary partitions.
The compact relaxed problem in dimension \(k-1\) has an
ordinary maximizer even if the residual form is semidefinite.
Among its relaxed maximizers, maximize the strictly convex
compact quadratic objective for \(T_{1/2}\). A nonordinary
maximizer of this secondary problem would admit the balanced
split above. Convexity forces both parts to maximize the
residual form and strict convexity prevents them both from
having no increase in the secondary form. This proves
attainment of the reduced value by an ordinary partition.

## Stability under mixtures of correlation tests

Fix any positive correlation \(r\) with a proved planar
three-cell score-relative margin. Thus every balanced
ordinary planar \(A\) admits an orientation of the tripod
\(T\) such that
\(S_r(T)-S_r(A)\ge c_r L_r(A;T)\), with \(c_r>0\).
The [sharp variance comparison](THREE_CELL_SHARP_VARIANCE_ABSORPTION.md)
supplies such anchors, including its upper endpoint.

Fix \(b\in[r,1)\). For any probability distribution \(\mu\)
on \([-r,b]\), consider sampling correlation from
\[
                  (1-\theta)\delta_r+\theta\mu .
 \tag{SD4}
\]
There is \(\theta_0(r,b)>0\) such that for every
\(0\le\theta\le\theta_0(r,b)\), the regular simplex uniquely
maximizes this averaged stability among balanced Gaussian
three-cell partitions, in every dimension \(d\ge2\), up to
orthogonal transformations, labels and null sets.
The support of \(\mu\) may include positive correlations
outside an individually proved interval.

In the plane, on fixing the orientation chosen for \(r\), set
\(h_i=1_{A_i}-1_{T_i}\). The exact expansion at any
\(t\in[-r,b]\) is \(\Delta_t=2L_t-Q_t\) with
\(Q_t=\sum_i\langle h_i,T_t h_i\rangle\); the transfer
quantity at negative \(t\) is interpreted by the same
signed score formula. The score difference for the two
labels vanishes on their reflecting plane for every \(t\).
Uniform boundedness and the Gaussian convolution gradient
bound for \(|t|\le b\) control its absolute value by
\(C_b\min(\operatorname{dist}(x,H_{ij}),1)\).
The positive \(r\)-score saturated margin bounds this same
minimum from below on each transfer, so
\(|L_t|\le C_{r,b} L_r\). Also \(|Q_t|\le Q_b\) by the
Hermite expansion and balance, and \(Q_b\le C_{r,b}L_r\)
by the [score-relative interaction estimate (B3)](POSITIVE_CORRELATION_CONTINUATION.md)
and comparison of the two positive score margins.
Thus \(|\Delta_t|\le C'_{r,b}L_r\) uniformly. Averaging gives a
strict deficit away from an oriented tripod as soon as
\((1-\theta)c_r-\theta C'_{r,b}>0\).

The eigenvalues of (SD4) decay exponentially. Their successive
differences for degrees \(j\ge1\) are at least
\(r^j[(1-\theta)(1-r)-\theta(1+r)]>0\) if
\(\theta<(1-r)/2\): this bounds the negative-correlation
contribution, and positive correlations have nonnegative
successive differences. Since the eigenvalues tend to zero,
the strict differences also give positivity. Applying (SD3)
to an actual maximizer reduces it to at most two dimensions,
where the preceding comparison gives the claimed value and
equality cases. No estimate on single-correlation optimizers
outside the certified bands is needed.
