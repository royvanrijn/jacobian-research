# The regular simplex loses for many equicorrelated samples

Fix \(k\ge3\) labels, pairwise correlation \(0<\rho<1\), and dimension
\(d\ge1\). For a balanced Gaussian partition \(A\), let \(M_m(A)\) be
the probability that \(m\) jointly Gaussian observations all receive
the same label. Each observation is standard and different
observations have cross covariance \(\rho I_d\). Write \(H_m\) for
the probability that all observations fall in a single halfspace of
Gaussian measure \(1/k\). Then

\[
 \sup_A M_m(A)=(2+o(1))H_m,\qquad
 H_m=m^{-\alpha+o(1)},\quad \alpha=\frac{1-\rho}{\rho}.
 \tag{MS1}
\]

The error is uniform over dimensions. The \(k\)-slab partition of a
single coordinate into Gaussian quantile intervals achieves the
asymptotic optimum. For the regular simplex \(T\) when \(d\ge k-1\),

\[
 M_m(T)=m^{-\beta\alpha+o(1)},\qquad \beta=\frac{2(k-1)}{k}>1.
 \tag{MS2}
\]

In particular the slab-to-simplex agreement ratio is
\(m^{(k-2)(1-\rho)/(k\rho)+o(1)}\). This refutes the many-sample
extension proposed in the conclusion of [Isaksson–Mossel,
Section 8](https://arxiv.org/pdf/0903.3362#page=26) even for
equicorrelated observations and exactly equal label masses. See also
[Isaksson's thesis, Section 3.1, pp. 16–17](https://publications.lib.chalmers.se/records/fulltext/120378.pdf#page=30).
The two-observation inequality has objective \(M_2\) and remains a
separate problem.

## High posterior probabilities for a simplex and a halfspace

Put \(s=\sqrt{\rho}\), \(\sigma=\sqrt{1-\rho}\), and \(a=\sigma/s\).
Write the observations as \(sX+\sigma Z_j\) for independent
standard Gaussian vectors \(X,Z_1,\ldots,Z_m\). Set
\(p_i^A(x)=\Pr(sx+\sigma Z_1\in A_i)\).
Then \(M_m(A)=\sum_i\mathbb E[p_i^A(X)^m]\).
Write \(\Phi\) for the one-dimensional Gaussian distribution function
and \(\bar\Phi=1-\Phi\).

For the halfspace \(B=\{x_1>h\}\), where \(h=\Phi^{-1}(1-1/k)\),
its posterior is \(\Phi((sX_1-h)/\sigma)\). Its upper tail is

\[
 \Pr(p^B(X)>1-\epsilon)=\bar\Phi(a t+h/s),\qquad
 t=\Phi^{-1}(1-\epsilon).
 \tag{MS3}
\]

For \(T\) use the active space \(\mathbf1^\perp\subset\mathbb R^k\),
where cell \(i\) is the region with largest \(i\)-th coordinate.
Its pairwise inward unit normals are \(n_{ij}=(e_i-e_j)/\sqrt2\).
If \(p_i^T(x)>1-\epsilon\), each containing halfspace has posterior
greater than \(1-\epsilon\), so \(x\cdot n_{ij}>at\) for all
\(j\ne i\). Their sum has norm \(\sqrt{k(k-1)/2}\). With
\(u_i=(\sum_{j\ne i}n_{ij})/|\sum_{j\ne i}n_{ij}|\), it follows that
\(x\cdot u_i>\sqrt\beta\,at\). Thus for every \(\epsilon\in(0,1/2)\)

\[
 \Pr(p_i^T(X)>1-\epsilon)\le \bar\Phi(\sqrt\beta\,at).
 \tag{MS4}
\]

Conversely \(n_{ij}\cdot u_i=1/\sqrt\beta\) for all \(j\ne i\).
Put \(t'=\Phi^{-1}(1-\epsilon/(k-1))\). On the event
\(X\cdot u_i>\sqrt\beta\,(at'+1)\) and
\(|X-(X\cdot u_i)u_i|\le1\) in the active space, every pair
margin is greater than \(at'\). The union bound gives posterior
greater than \(1-\epsilon\). The two conditions in this event
are independent and the second has fixed positive probability.
Together with \(\log\bar\Phi(t)=-t^2/2+o(t^2)\) and \(t'/t\to1\),
this proves the logarithmic tails

\[
 \Pr(p^B(X)>1-\epsilon)=\epsilon^{\alpha+o(1)},\qquad
 \Pr(p_i^T(X)>1-\epsilon)=\epsilon^{\beta\alpha+o(1)} .
 \tag{MS5}
\]

If \(P\in[0,1]\) has \(\Pr(P>1-\epsilon)=\epsilon^{\eta+o(1)}\)
for fixed \(\eta>0\), then \(\mathbb E P^m=m^{-\eta+o(1)}\).
For the lower bound, restrict to \(P>1-1/m\). For the upper bound,
split at \(1-C\log(m)/m\), using \(C>\eta\). This proves the power
laws in (MS1) and (MS2).

The halfspace term also has the equivalent

\[
 H_m\sim \Gamma(1+\alpha)\,
             \bar\Phi(a u_m+h/s),\qquad
 u_m=\Phi^{-1}(1-1/m).
 \tag{MS5a}
\]

Indeed \(H_m=\mathbb E\bar\Phi(a V_m+h/s)\), where \(V_m\) is the
maximum of \(m\) independent standard normals from the resampling
coordinates, after a sign change. The variable
\(w=m\bar\Phi(V_m)\) has density \((1-w/m)^{m-1}\) on \((0,m)\)
because \(\bar\Phi(V_m)\) is the minimum of \(m\) uniform variables.
Put \(v=\bar\Phi^{-1}(w/m)\). The logarithm of the ratio
\(\bar\Phi(a v+h/s)/\bar\Phi(a u_m+h/s)\) is
\(\int_v^{u_m}a\lambda(a z+h/s)\,dz\), where
\(\lambda=\varphi_1/\bar\Phi\) is the Gaussian hazard. Also
\(\log w=\int_v^{u_m}\lambda(z)\,dz\). The Mills asymptotic
\(\lambda(z)/z\to1\) gives the pointwise ratio limit \(w^\alpha\)
for fixed \(w>0\). For \(w\le1\) the ratio is at most one. For
\(1\le w\le m/2\) it is at most \(w^{C_0}\), since
\(a\lambda(az+h/s)/\lambda(z)\) is bounded by some \(C_0\) on
\([0,\infty)\). On those two ranges the densities are bounded
by \(e^{-w/2}\) for \(m\ge2\). The remaining integrated ratio
is bounded by \(2^{-m}/\bar\Phi(a u_m+h/s)\), which tends to zero
by the Gaussian tail asymptotic. Dominated convergence now gives
\(\int_0^\infty w^\alpha e^{-w}\,dw=\Gamma(1+\alpha)\).

There is also a consequence for the law of the whole posterior
vector. For a balanced partition with one halfspace cell, the
expectation of the symmetric convex function
\(p\mapsto\sum_i(p_i-b)_+\) strictly exceeds the simplex value for
all \(b\) sufficiently close to \(1\). Indeed (MS3) eventually
exceeds \(k\) times (MS4), by their different tail exponents, and
integrating the tail inequality from \(b\) to \(1\) gives the
claim. A symmetric convex-order comparison with the simplex at
the top therefore fails even on ordinary balanced partitions.

## A sharp dimension-independent upper bound

For an arbitrary partition and \(0<\epsilon<1/2\), set
\(D_i=\{x:p_i^A(x)>1-\epsilon\}\), and keep \(t=\Phi^{-1}(1-\epsilon)\).
If \(x\in D_i\), \(y\in D_j\), and \(i\ne j\), then

\[
                         |x-y|\ge 2at.                          \tag{MS6}
\]

In fact the two Gaussian laws with means \(sx,sy\) and covariance
\(\sigma^2I_d\) have respective probabilities greater than
\(1-\epsilon\) of the disjoint sets \(A_i,A_j\). The sum of the two
errors of the test \(A_i\) is at most \(2\epsilon\).
For the two densities, the smallest sum of errors is achieved
by assigning each point to its larger density. The comparison
is a halfspace test with threshold at the midpoint of the means,
and the sum is \(2\bar\Phi(s|x-y|/(2\sigma))\) by projection.
This gives (MS6).

For any three nonempty sets \(D_i,D_j,D_\ell\), take one point
from each and denote their average by \(z\). Writing \(x_1,x_2,x_3\)
for the three points gives

\[
 \tfrac13\sum_{\nu=1}^3|x_\nu|^2-|z|^2
 =\tfrac19\sum_{\mu<\nu}|x_\mu-x_\nu|^2\ge (2at)^2/3 .
\]

Let \(\varphi_d\) be the standard Gaussian density. This identity
implies that \(e^{-2a^2t^2/3}\varphi_d(z)\) dominates the geometric
mean of the three Gaussian densities restricted to their respective
sets. The pointwise domination also holds for triples with a point
outside its prescribed set, because the geometric mean is then zero.
Prékopa–Leindler with three weights \(1/3\) shows

\[
 \big(\gamma_d(D_i)\gamma_d(D_j)\gamma_d(D_\ell)\big)^{1/3}
                    \le e^{-2a^2t^2/3}.
 \tag{MS7}
\]

The functions are Borel and integrable. The three-function form
with arbitrary positive weights summing to one appears in
[Prékopa, Theorem 3](https://acta.bibl.u-szeged.hu/14411/1/math_034_335-343.pdf#page=3).

At most two labels have \(\gamma_d(D_i)>e^{-2a^2t^2/3}\).
For these labels the Exchangeable Gaussians Theorem bounds each
entire contribution \(\mathbb E[p_i^A(X)^m]\) by \(H_m\). Its
hypotheses hold because the \(m\) standard Gaussian observations
are equicorrelated and each tested set is the same cell of mass
\(1/k\). [Isaksson–Mossel, Theorem 1.2](https://arxiv.org/pdf/0903.3362#page=2).
For every other label, splitting at \(1-\epsilon\) bounds the
contribution by \((1-\epsilon)^m+e^{-2a^2t^2/3}\). Select two
labels including all those above the measure threshold. This
gives the finite-\(m\) upper bound

\[
 M_m(A)\le 2H_m+(k-2)\big((1-\epsilon)^m+e^{-2a^2t^2/3}\big).
 \tag{MS8}
\]

Choose \(\epsilon=C\log(m)/m\) with \(C>\alpha=a^2\) and \(m\)
large. The two terms in parentheses are bounded by \(m^{-C}\)
and \(m^{-4\alpha/3+o(1)}\), respectively, so they are \(o(H_m)\).
The lower bound \(2H_m\) comes from the two opposite extreme
halfspaces of the balanced quantile slab partition. This proves
(MS1) uniformly in \(d\), with an explicit comparison (MS8).

## A balanced low-influence voting counterexample for many copies

For each voter, first choose a common uniform ballot in \([k]\).
Independently for each of the \(m\) copies, retain that ballot
with probability \(\sqrt\rho\) and otherwise replace it by an
independent uniform ballot. Do this independently across voters.
The \(m\) indicators of whether that voter chose candidate \(1\)
have pairwise correlation \(\rho\). For each copy, center and
normalize the count of votes for candidate \(1\). The joint
central limit theorem for the \(m\) counts gives the scalar
equicorrelated model above. Apply the \(k\)-quantile slab rule to
these normalized counts. For fixed \(m\) its agreement converges
to \(M_m\) of the slab partition, and its label probabilities
tend to \(1/k\), since all the Gaussian thresholds have null
probability.

Changing a single vote moves the normalized count by at most
a constant times \(n^{-1/2}\) for \(n\) voters. It can change
the result only near one of the finitely many thresholds,
whose combined probability tends to zero by the central limit
theorem and continuity of the Gaussian law. This bound is
uniform over voters. For exact balance, relabel surplus profiles
to deficient labels until every label receives \(k^n/k\) profiles.
The fraction \(\delta_n\) changed tends to zero. Coordinate
resampling changes its probability of altering the outcome by
at most \(2\delta_n\), by the two uniform marginals. The
agreement probability over \(m\) copies changes by at most
\(m\delta_n\). Thus the repaired rules are exactly balanced and
their maximum influence tends to zero.

For plurality, the joint central limit theorem applied to the
full vote count vectors gives Gaussian input in the label
hyperplane \(\mathbf1^\perp\), after a common positive scaling.
The cross covariance between distinct copies on that hyperplane
is \(\rho I\) after normalization. The plurality winner converges
to the simplex rule since pairwise Gaussian ties are null.
Therefore for every sufficiently large fixed \(m\) the slab
rules have strictly greater limiting probability of unanimous
agreement than plurality in this model.
