# Finite reflection tests for balanced positive-noise maximizers

Fix \(k\ge3\), a compact band \([a,b]\subset(0,1)\), and the active
dimension \(n=k-1\). Write \(T_\rho\) for Gaussian noise and \(\Pi_j\)
for projection onto Hermite degree \(j\). For a balanced maximizer \(f\)
and \(g=R\pi_{ij}f\), the defect \(d_\gamma(f,g)=\sum_l\|f_l-g_l\|_1\)
has a finite Hermite approximation with a uniform \(O(N^{-1/2})\) error.
At sufficiently high finite degree, a small projected spanning-tree
defect forces a four- or five-cell maximizer to be a simplex.
A strict comparison on the balanced moment body suffices to verify
this premise for every maximizer on a compact band.

Heilman–Tarter's [Theorem 1.9](https://arxiv.org/pdf/2011.05583v2)
reduces ordinary maximizers in every finite ambient \(d\ge k-1\)
to \(n\) active coordinates, up to rotations and null sets. Noise
stability is unchanged under cylindrical extension. The local
step uses [reflection rigidity](REFLECTION_SYMMETRY_OF_MAXIMIZERS.md#approximate-reflection-symmetries)
and the effective-neighborhood argument of
[(B10)–(B11)](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison).

## Balanced replacement supplies a uniform signal

Among balanced fractional assignments, a maximum exists by weak
compactness and compactness of \(T_\rho\) in finite dimension. It must
be ordinary. If two labels overlap on positive measure, restrict to
a set where both are bounded below by a common positive number and
split that set into two equal-measure pieces. Transferring opposite
small signed amounts on these pieces gives two balanced feasible
fractional assignments whose midpoint is the original assignment.
Since all eigenvalues \(\rho^j\) of \(T_\rho\) are strictly positive,
the stability is strictly convex along the transfer. Thus the
fractional and ordinary maxima agree. The first variation at a
maximizer \(f\) gives \(\langle f-g,T_\rho f\rangle\ge0\) for every
balanced fractional \(g\), with summation over labels.

For a pair \(i,j\), put \(U=A_i\cup A_j\), \(p=2/k\) and
\(h=1_{A_i}-1_{A_j}\). If \(Q_\rho(v)=\langle v,T_\rho v\rangle\),
the pair contributes \((Q_\rho(1_U)+Q_\rho(h))/2\).
Repartition \(U\) by the median \(m\) of any unit-direction projection
within \(U\). This has exactly the prescribed masses and leaves all
other labels fixed. The contrast of the replacement has its first
moment in that direction equal to \(\int_U |x\cdot e-m|\,d\gamma_n\).
Writing \(\phi_0=(2\pi)^{-1/2}\), the strip of width \(2s\) about
\(m\) has mass at most \(2\phi_0s\). Integrating the lower bound
\((p-2\phi_0s)_+\) gives a first moment of at least \(p^2/(4\phi_0)\).
The other Hermite contributions to the replacement are nonnegative.
Global maximality gives, for every pair,
\[
 Q_\rho(h)\ge \rho\,\frac{p^4}{16\phi_0^2}\ge \kappa,
 \qquad \kappa=\frac{a p^4}{16\phi_0^2}.
 \tag{FR1}
\]
This uses a finite balanced replacement of the two labels.

The exact scores \(s_i=T_\rho f_i-\lambda_i\) can be chosen so that
\(f\) assigns to their pointwise winner. To see this, minimize
\(\int\max_i(T_\rho f_i-\lambda_i)\,d\gamma_n+k^{-1}\sum_i\lambda_i\)
over the mean-zero price gauge. The expression is continuous and
coercive there since the unpriced scores are nonnegative. Each score
difference is real analytic and nonconstant, by injectivity of
\(T_\rho\) and the zero mean and nonzero \(L^2\) norm of \(f_i-f_j\).
Thus for every price vector ties are null. Differentiation at a
minimizing vector gives winning masses \(1/k\). The resulting winner
assignment attains the linear maximum at these prices. The first
variation makes \(f\) a linear maximizer too, so it coincides almost
everywhere with this assignment. The unpriced scores lie in \([0,1]\)
and both cells win on positive measure, hence
\(|\lambda_i-\lambda_j|\le1\).

## Uniform Gaussian perimeter and a Hermite tail

There is a finite, effective \(C=C(a,b,k)\) such that every
actual maximizer on the band, in its active space, satisfies
\[
 \sum_i P_{\gamma_n}(A_i)\le C,\qquad
 \sum_{j>N}\sum_i\|\Pi_j1_{A_i}\|_2^2
       \le \frac{\sqrt\pi\,C}{2(1-e^{-1})\sqrt{N+1}}\quad(N\ge1).
 \tag{FR2}
\]
Here \(P_{\gamma_n}\) is the Gaussian-weighted perimeter. The
perimeter bound follows by counting zeros on rays from several nearby
points where a given pair score is separated from zero.

Fix a pair score \(F=s_i-s_j\). Balance gives
\(\int hF=Q_\rho(h)\), and \(|F|\le2\). Choose \(R_0\) with
\(2^{n/2}e^{-R_0^2/4}\le\kappa/4\). The Gaussian radial exponential
moment then gives a point \(x_0\) of norm at most \(R_0\) with
\(|F(x_0)|\ge\kappa/2\).
On every line with real base point \(x\) and unit direction \(v\),
\(F(x+zv)\) extends to an entire function of \(z\), with the bound
\[
 |F(x+zv)|\le1+\exp\!\left(\frac{b^2(\operatorname{Im}z)^2}{2(1-b^2)}\right).
\]
The real part of the complex Gaussian convolution exponent is
the real displaced kernel exponent plus the indicated
imaginary-square term. The bound \(|h|\le1\) and domination
on compact sets justify the extension and the estimate.

Differentiating in a unit direction in the Gaussian convolution
also gives a global Lipschitz bound \(L_b=b/\sqrt{1-b^2}\) for \(F\).
Use the anchors \(x_0\) and \(x_j=x_0+h e_j\) for the coordinate
unit vectors \(e_j\), with \(h=\min\{1,\kappa/(4L_b)\}\).
Every anchor \(w\) has \(|F(w)|\ge\kappa/4\).
On a unit ray from \(w\), Jensen's formula on the complex disk
of radius \(2A\), for \(A=R+R_0+1\), bounds the number of zeros
on \([0,A]\) by
\[
 Z_R=\frac{\log\!\left(1+e^{2b^2A^2/(1-b^2)}\right)
                  +\log(4/\kappa)}{\log 2}.
\]
The entire bound gives the numerator and every such root has
radius at most \(A\). If a root is on the external circle one
uses the limit through nearby radii.

For the sign indicator \(u=1_{\{F>0\}}\), the distribution
\((x-w)\cdot Du\) acts on a smooth test \(\varphi\) supported
in the ball \(B_R\) by
\[
 \big\langle (x-w)\cdot Du,\varphi\big\rangle
   =-\int_{\mathbb S^{n-1}}\int_0^\infty
          u(w+tv)\frac{d}{dt}\big(t^n\varphi(w+tv)\big)\,dt\,dv.
\]
The identity follows by polar integration of the divergence of
\((x-w)\varphi\). Integration by parts on each ray bounds its
absolute value by \(|\mathbb S^{n-1}|A^nZ_R\|\varphi\|_\infty\):
the sign has at most \(Z_R\) jumps where the test is supported.
The terms at zero and infinity vanish. Taking the difference
between the expressions based at \(x_0\) and at \(x_j\) shows
that the variation of \(D_j u\) on \(B_R\) is at most
\(2|\mathbb S^{n-1}|A^nZ_R/h\).
Every winner cell is a finite intersection of pair sign sets
up to a null set. Its perimeter on \(B_R\) is bounded by the
sum of their perimeters, by the product rule for bounded
variation, or by first mollifying each indicator in \([0,1]\).
Summing in pairs, directions and integer Gaussian shells gives
\(C(a,b,k)<\infty\) in (FR2). The bound is effective since
\(A^nZ_R=O_{a,b,k}((1+R)^{n+2})\) and the Gaussian density on
\(B_{\ell+1}\setminus B_\ell\) is at most
\((2\pi)^{-n/2}e^{-\ell^2/2}\).

For the Hermite estimate let \(X,Y\) be independent standard
Gaussians. Along \(X_s=X\cos s+Y\sin s\), the derivative
\(\dot X_s\) is an independent standard Gaussian for each \(s\).
For smooth scalar \(v\), integrating its derivative along this path
bounds \(\mathbb E|v(X_\theta)-v(X)|\) by
\(\theta\sqrt{2/\pi}\int|\nabla v|\,d\gamma_n\).
For indicators of finite Gaussian perimeter, the same bound
follows by cutting off at large radius and mollifying, with
convergence in Gaussian \(L^1\) and the weighted variation upper
bound. Indeed cutoff contributes at most the Gaussian mass of
its outer strip times its derivative bound. The compactly
supported mollification bound follows by convolving the
Gaussian weight against the variation measure and taking the
mollifier radius to zero. Summing over the indicators gives
\[
 1-S_t(A)\le \frac{\arccos(t)}{\sqrt{2\pi}}\sum_i P_{\gamma_n}(A_i).
\]
Take \(t=1-1/(N+1)\), for which \(1-t^j\ge1-e^{-1}\) for \(j>N\)
and \(\arccos(t)\le\pi/\sqrt{2(N+1)}\). The Hermite expansion
gives the second assertion of (FR2).

## A finite sufficient comparison on a positive correlation band

For \(k=4,5\), let \(\varepsilon>0\) suffice for the spanning-tree
hypothesis of [reflection rigidity](REFLECTION_SYMMETRY_OF_MAXIMIZERS.md#approximate-reflection-symmetries)
on \([a,b]\). Let \(\mathcal K_N\) be the compact convex set of
first-through-\(N\) Hermite coefficients of balanced fractional
\(k\)-label assignments in \(\mathbb R^{k-1}\).
For a tuple \(P\) in this space, write \(P_j\) for its degree-\(j\)
part and set
\[
 \begin{split}
 H_N(P,\rho)&=1/k+\sum_{j=1}^N\rho^j\|P_j\|_2^2,\\
 \tau_N(P)&=1-1/k-\sum_{j=1}^N\|P_j\|_2^2,\\
 D_N(P)&=\min_{\mathcal T}\ \max_{ij\in\mathcal T}\
       \min_{\|v\|=1}\|P-R_v\pi_{ij}P\|_2^2 .
 \end{split}
\]
The first minimum in \(D_N\) ranges over label spanning trees and
\(R_v\) is linear reflection in the plane perpendicular to \(v\).
The group action preserves each Hermite degree. For an indicator
maximizer with tuple \(P\), \(\tau_N(P)\) is its exact omitted
Hermite weight. For each pair and each reflection, the squared
norm of the full reflected difference equals
\(d_\gamma(A,R_v\pi_{ij}A)\). Its omitted squared norm is at most
\(4\tau_N(P)\), by the triangle inequality on the tail subspace.
If \(D(A)\) denotes the raw spanning-tree defect, taking the
minima and maximum yields
\[
 D_N(P)\le D(A)\le D_N(P)+4\tau_N(P),\qquad
 S_\rho(A)\le H_N(P,\rho)+\rho^{N+1}\tau_N(P).
 \tag{FR3}
\]
In particular (FR2) controls the error in the finite reflection
defect by \(O_{a,b,k}(N^{-1/2})\), uniformly in the original
ambient dimension.

An explicit sufficient comparison for simplex optimality with
equality rigidity throughout the band is
\[
 \sup_{\substack{P\in\mathcal K_N,\ \tau_N(P)\le
                 \sqrt\pi C/(2(1-e^{-1})\sqrt{N+1})\\
                   D_N(P)+4\tau_N(P)\ge\varepsilon}}
         \bigl(H_N(P,\rho)+\rho^{N+1}\tau_N(P)\bigr)<S_\rho(1_T)
 \qquad\text{for every }\rho\in[a,b].
 \tag{FR4}
\]
Every maximizer satisfies the first constraint by (FR2). The
inequality and (FR3) force \(D(A)<\varepsilon\). The reflection
bridge and the local gap in the active space finish the
comparison. Dimension reduction gives the ambient conclusion.

For specifying the finite comparison the support function of
\(\mathcal K_N\) is
\[
 \sigma_{\mathcal K_N}(c)
  =\inf_{\lambda\in\mathbb R^k}
       \left\{\frac1k\sum_i\lambda_i+
       \int\max_i\left(\sum_{|\alpha|=1}^N c_{i,\alpha}H_\alpha(x)
                                  -\lambda_i\right)d\gamma_n(x)\right\}.
 \tag{FR5}
\]
Here \(H_\alpha\) are orthonormal Hermites. Compactness of the convex
image of the feasible fractional assignments under the joint map
of masses and coefficients ensures this dual formula. More
explicitly, maximize the linear polynomial-score integral with
prescribed masses \(p\). Its optimum is an upper semicontinuous
concave function on the full mass simplex. It has a supporting
supergradient at the interior point \(p_i=1/k\). Those prices
make the prescribed-mass optimum equal to the unrestricted
maximum of the priced integral plus the price term in (FR5).
Pointwise maximization computes the unrestricted value even if
some polynomials tie identically.

The inequalities in (FR4) are finite-dimensional. A strict
uniform margin admits verification by support halfspaces,
uniform direction nets for the reflections, and a correlation
mesh. Continuity of \(D_N\) follows from the finite tree family
and compactness of the unit sphere. In particular, the
maximizing head value on a decreasing sequence of compact outer
polytopes and decreasing tolerances converges to the value on the
closed bad-defect subset of \(\mathcal K_N\). Polynomial Gaussian
tails bound truncation in the support integrals. The dependence on
correlation is uniformly Lipschitz on the band for both the head
and \(S_\rho(1_T)\), using the Hermite series. This is a
reflection-based sufficient certificate complementary to
[(B17)](POSITIVE_CORRELATION_CONTINUATION.md#a-finite-test-for-distant-competitors-on-a-correlation-band).
No sign in (FR4) is established here.

The finite comparison is exhaustive in the limit under strict
global uniqueness on the band. If it failed along \(N\to\infty\),
take fractional realizers \(f_N\) of the witnessing tuples at
\(\rho_N\to\rho_*\). Their actual stability is at least the
head and thus at least \(S_{\rho_N}(1_T)-b^{N+1}\). A weak limit
is optimal at \(\rho_*\) by compactness of the noise operators
and operator norm continuity in correlation. It is an ordinary
simplex by the uniqueness premise. Since \(\|f_N\|_2^2\le1\)
and the limiting norm squared is one, convergence is strong
along that subsequence. It follows that \(\tau_N(P)\to0\) and,
by using the fixed transposition reflections of the limiting
simplex, \(D_N(P)\to0\). This contradicts the bad-defect
constraint.

## Balanced copies across a bisecting reflection

A second exact condition on a maximizer comes from copying a
halfspace by a reflected transposition. For fixed \(i,j\) there
is a central plane bisecting every spectator cell \(A_l\) with
\(l\notin\{i,j\}\). The \(k-2\) signed halfspace imbalances form
a continuous odd map of the unit sphere in the active space
into \(\mathbb R^{k-2}\), so the
[Borsuk–Ulam theorem](https://doi.org/10.4064/fm-20-1-177-190)
gives a zero.
The combined mass of \(A_i,A_j\) in either halfspace is then
\(1/k\).

For any such reflection \(R\), write \(F(x)=f(x)\) and
\(G(x)=\pi_{ij} f(Rx)\) on its positive halfspace \(H^+\).
Using \(F\) on \(H^+\) and \(\pi_{ij}F\circ R\) on \(H^-\)
gives a balanced assignment \(B_F\), and \(B_G\) is defined
the same way using \(G\). Spectator masses agree by bisection,
and the two swapped labels each have the combined pair mass
of the positive halfspace. With \(M_\rho\) the Mehler kernel
relative to the Gaussian measure, direct expansion gives
\[
 \frac{S_\rho(B_F)+S_\rho(B_G)}2-S_\rho(A)
  =\iint_{(H^+)^2}M_\rho(x,Ry)
             \langle F(x)-G(x),\,\pi_{ij}(F(y)-G(y))\rangle
                            \,d\gamma(x)d\gamma(y)\le0 .
 \tag{FR6}
\]
The same-halfspace terms cancel. The two cross terms are
equal by symmetry of \(M_\rho(x,Ry)\), so their difference
is the displayed quadratic form. Sorting the two copies
by (RG1) is a further feasible improvement of the left side.
This is a finite variation test on a bisecting plane. A sign
bound yielding small reflection defect from it remains open.

## A limitation of first variation alone

For any \(k\ge3\), \(n\ge1\), and \(0<\rho<1\), maximize stability
over balanced fractional assignments measurable with respect to
\(|x|\). The radial subspace is weakly closed and invariant under
\(T_\rho\), and a maximum exists by compactness. The Gaussian radius
is nonatomic, so the same two-sided overlap perturbation proves
that a radial maximum \(A\) is ordinary. For any balanced
fractional \(B\), the conditional expectation \(\bar B\) given
\(|x|\) is feasible in the radial class. Radial maximality gives
\[
 2\langle B-A,T_\rho A\rangle
  =2\langle\bar B-A,T_\rho A\rangle
  \le -\langle\bar B-A,T_\rho(\bar B-A)\rangle .
 \tag{FR7}
\]
If \(\bar B=A\), the fact that \(A\) takes values at the vertices
of the label simplex implies \(B=A\) almost everywhere. Thus (FR7)
is strict for every other \(B\). Each centered reflection fixes
\(A\) cellwise, so \(d_\gamma(A,R\pi_{ij}A)=4/k\) for every
reflection and every transposition. Each reflection-orbit flip
that preserves unordered labels is trivial on this partition.

When \(n\ge k-1\), this radial partition has strictly smaller
stability than the simplex for \(0<\rho<F_k/(1-1/k)\), where
\(F_k=\sum_i\|\Pi_1 1_{T_i}\|_2^2>0\). Indeed its first weight
vanishes and its higher contribution is bounded by
\(\rho^2(1-1/k)\). The simplex has first contribution \(\rho F_k\).
Thus even strict linear optimality and all fixed-unordered-label
reflection-orbit tests allow large reflection defects. The
finite replacement bound (FR1) and the stability comparison
(FR4) use global maximality beyond those linear tests.


## Cross-interface normals require internal interactions

An actual maximizer with \(k\ge3\) has a nonbipartite interface
adjacency graph. Here two labels are adjacent when their common
reduced boundary has positive Gaussian surface measure. In
particular a four-cell maximizer has an interface triangle. A
triangle-free five-cell maximizer would have a chordless
five-cycle. This graph records contacts and can differ from
the reflection graph.

To prove the restriction, lift the active maximizer to one extra
Gaussian coordinate \(z\). The lift still maximizes by
Heilman–Tarter dimension reduction, since an improvement in the
lifted dimension would have a representative in the original
active dimension. Use their regularity and second variation
[Lemmas 2.4, 4.2 and 4.8](https://arxiv.org/pdf/2011.05583v2).
Orient the interfaces crossing any fixed label split \(C\,|\,D\)
from \(C\) to \(D\). Test normal speed \(z\) on these interfaces
and zero on the internal interfaces. For the compact tests
multiply by smooth active cutoffs exhausting their regular
two-cell faces and replace \(z\) by smooth compactly supported
odd \(z\chi(z)\). Extend the face speeds by local normal fields
away from the junctions. Exact volume corrections are allowed
since each face speed has zero first-order mass change by
oddness. First remove the \(z\) cutoff, using the bounded noise
operator in \(L^2(\gamma_1)\) and \(T_\rho z=\rho z\).
Then increase the nonnegative active cutoffs to one. The
cross-edge interaction is monotone because shared-label
incidence signs agree. The local term is bounded by the
uniform gradient bound for the smoothed scores and (FR2).

Write \(\nu\) for Gaussian surface measure in the active
space, \(N_l\) for the exterior normal on cell \(l\), and
\(a_{ij}=\|\nabla T_\rho(1_{A_i}-1_{A_j})\|\) on its interface.
At a point \(x\) on a cross-face let
\[
 U_l(x)=\int_{\partial^* A_l,\ {\rm cross}} M_\rho(x,y)\,d\nu(y),
 \quad
 V_l(x)=\int_{\partial^* A_l,\ {\rm cross}} M_\rho(x,y)N_l(y)\,d\nu(y),
 \quad
 I_l(x)=\int_{\partial^* A_l,\ {\rm internal}} M_\rho(x,y)\,d\nu(y).
\]
The second variation gives
\(\sum_{\rm cross}\int (U_i+U_j)\,d\nu\le
  \sum_{\rm cross}\int a_{ij}/\rho\,d\nu\).
Indeed the cross-edge incidence inner products sum to
\(U_i+U_j\), with factor \(\rho\) from the transverse
linear speed. Gaussian divergence in the Mehler formula gives
\(\nabla T_\rho1_{A_l}(x)
  =-\rho\int_{\partial^* A_l}M_\rho(x,y)N_l(y)\,d\nu(y)\).
For fixed \(x\) the kernel is bounded as a function of \(y\).
The finite Gaussian perimeter justifies the boundary integral.
All the double integrals are finite as well. The Mehler kernel
times the two Gaussian surface densities is bounded by
\(C_{b,n}\exp[-(|x|^2+|y|^2)/(2(1+b))]\), and the
unweighted local perimeters grow at most polynomially by the
ray estimate in (FR2). Taking norms yields the necessary cut
inequality
\[
 \sum_{\rm cross}\int
       \bigl(U_i+U_j-\|V_i-V_j\|\bigr)\,d\nu
       \ \le\ \sum_{\rm cross}\int (I_i+I_j)\,d\nu .
 \tag{FR8}
\]
If the adjacency graph were bipartite, the right side would
vanish for its bipartition. The integrands on the left are
nonnegative. Equality and strict positivity of Mehler on each
boundary force the exterior normal of each incident cell to be
a single oriented unit vector on its entire reduced boundary,
by equality in the norm of the integral of unit vectors.
Each cell has positive boundary measure since it has measure
strictly between zero and one. Each is then a halfspace,
since the distributional derivative vanishes tangentially and
has one sign in the fixed normal direction. At most two
nontrivial halfspaces can be pairwise disjoint up to null
sets. This contradicts \(k\ge3\).

## Excluding constant thresholds between two reflection components

Suppose the exact reflection graph of a maximizer has two
components \(C,D\) and the threshold in
[(RG4)](REFLECTION_SYMMETRY_OF_MAXIMIZERS.md#two-components-give-parallel-halfspace-fibers)
is constant. Positive masses make the constant \(c\) finite.
Then for each \(i\in C\) the label indicator is a product of its
internal simplex cell \(T_i^C\) with the halfspace \(z>c\), and
for \(j\in D\) the product uses \(z<c\).
On \(z=c\), with the internal coordinates in the respective
open winner chambers, both labels occur on opposite sides.
Their priced scores therefore agree there by continuity. The
unpriced scores there are positive constants times
\(T_\rho^{V_C}1_{T_i^C}(p_C)\) and
\(T_\rho^{V_D}1_{T_j^D}(p_D)\), respectively. The chamber
coordinates vary independently, so each response must be
constant on its chamber. At least one component has at least
two labels. In that component the winning response is real
analytic, equals the reciprocal of the component size at the
origin, and tends to one along any interior winner ray.
It cannot be constant on an open chamber. This excludes
constant-height product candidates in the \(3+1,2+2,4+1\)
and \(3+2\) cases.


## Relation to the near-one localization

The existing all-cell near-one compactness argument in
[the fourth-score localization proof](THREE_CELL_EXTERIOR_TUBES.md)
already makes \(D(A)\to0\) as \(\rho\uparrow1\) for each fixed
\(k\), uniformly over ambient dimensions for actual maximizers.
Indeed a simplex \(T'\) at distance \(d_\gamma(A,T')\) supplies
its exact transposition reflections as candidate reflections
for \(A\), each with error at most \(2d_\gamma(A,T')\) by the
triangle inequality. Near the endpoint the local radius varies
with \(\rho\), and the available compactness limit does not
give a rate relative to that radius. On a compact interior
band the sufficient comparison (FR4) instead uses a fixed
positive local threshold.
