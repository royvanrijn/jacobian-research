# Generic ordinary cores are stable under diffuse tails on rational slopes

Fix a finite list of rational numbers \(\mathcal S\subset[-1,1]\).
Let \(J\ge4\) be any integer at least the
[absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Fix positive integral gaps \(d_i\), orientations
\(\epsilon_i\in\{+1,-1\}\), and positive integral multiplicities \(m_i\)
for \(r\ge2(J-1)\) groups of ordinary rows. The smaller exponent
\(A_i\ge0\) varies independently in each group. All \(m_i\) copies
in that group have smaller exponent \(A_i\).

Append any number of unequal shared ordinary rows. For the tail set
\[
  \tau=\sum_{\rm tail}\frac{(b-a)^2}{\max(a,b)},\qquad
  \lambda=\sum_{\rm tail}\frac{b-a}{\max(a,b)},\qquad
  \eta=\max_{\rm tail}\frac{|b-a|}{\max(a,b)} .
\]
The empty tail has all three quantities zero.

**Theorem.** There is a proper algebraic subset of the \(r\)-dimensional
base space, defined over \(\mathbb Q\), such that each ordinary core
outside it has positive computable tolerances \(\delta_C,\eta_C\)
with the following property. Every ordinary augmented profile
satisfying
\[
                 \eta\le\eta_C,\qquad
                 |\lambda-\sigma\tau|\le\delta_C
                  \quad\text{for some }\sigma\in\mathcal S
\]
satisfies SFC(2) in every degree. There is no restriction on the
number of tail rows or on their total signed gap. The tolerances
are for the fixed core. The exceptional fraction in
\(\{0,\ldots,H\}^r\) is \(O((H+1)^{-1})\) for the fixed
gaps, orientations, multiplicities and slope list.

In particular one may use \(\mathcal S=\{0\}\). If all tail
rows have a common larger exponent \(M\), then
\(\lambda=h/M\), where \(h=\sum_{\rm tail}(b-a)\). A possible
failure has bounded \(h\) by the physical root annuli in the
proof below. Thus for each core covered by the theorem, every
sufficiently diffuse such tail is regular, with arbitrary
multiplicities and orientations of its own rows.

## A single terminal, integer twist and rational slope

Fix \(n=N-1\), where \(4\le N<J\), an integer \(h\), and a
slope \(\sigma=e/b\) with \(b>0\).
Strip the factors \((jA_i)!\) from the factorial rows at
orders \(j=n,N\), including their multiplicities. The resulting
core polynomials \(Q_j(A,z)\) have coefficients polynomial
in \(A\) over \(\mathbb Q\). Work on the open set where all their
sampled rising factors \(jA_i+\ell\) are nonzero. This includes
all ordinary nonnegative integral bases. Normalize each constant
coefficient to one.

The diffuse reference on this slope is
\[
 \begin{split}
 F_n(A,t,z)
   &=\sum_{k=0}^{n}[z^k]Q_n(A,z)\,
              t^{bN k(n-k)} z^k,\\
 F_N(A,t,z)
   &=\sum_{k=0}^{N}[z^k]Q_N(A,z)\,
              t^{bn k(N-k)+e k} (N/n)^{hk} z^k,
 \qquad t=\exp[-\tau/(2bnN)] .
 \end{split}                                                     \tag{1}
\]
Here \(Q_j\) denotes the constant-normalized stripped core.
For \(t\in\mathbb C^*\), these formulas define an algebraic family.

Consider the incidence of two distinct nonzero common roots of
the pair in (1). For each fixed \(t\), its image in the base
space has dimension at most \(r-2\).
To see this, use the
[paired coefficient rank argument](MULTIVARIATE_EXCEPTIONAL_CODIMENSION_TWO.md).
The logarithmic differential column for each base is multiplied by
its positive multiplicity \(m_i\). This multiplication leaves every
rank bound in that argument intact, including the bound on the
rank-defective locus. On the full-rank chart, fixing the lower
coefficients leaves the first two upper coefficients free.
Each of the two distinct nonzero common roots imposes one independent
linear condition on those two coefficients. In fact their
evaluation determinant is \(\zeta_1\zeta_2(\zeta_2-\zeta_1)\).
When one unit-gap orientation accounts for all but at most one of
the base groups, reverse both polynomials if that orientation is
decreasing. Reversal is an invertible transformation on the
nonzero roots, preserves the Gaussian multipliers, and inverts
the positive rescaling, replacing \(h,e\) by \(-h,-e\).
The dominant-orientation chart then fixes a ratio of the last two
upper coefficients after fixing the other bases and the lower
coefficients. Those last two indices are at least three,
so the first two upper coefficients are free in this chart as well.
Multiplication of the coefficients by the nonzero \(t\)- and
\(h\)-dependent factors in (1) preserves each of these arguments.

Consequently the incidence image in \((A,t)\) has dimension
at most \(r-1\): the parameter \(t\) has dimension one, and every
fiber has dimension at most \(r-2\). Projection onto \(A\) and
Zariski closure preserve this upper bound. This gives a proper
algebraic exceptional set for the chosen \(N,h,\sigma\).

At a real positive \(t\le1\) in (1), a physical common root gives
two distinct nonzero common roots by conjugation. The even
reference is strictly positive on the real axis. At \(\tau=0\)
this is the exponential integral of an even power of the real
binomial. At positive \(\tau\), the Gaussian coefficient
multiplier averages that polynomial over a positive lognormal
rescaling of its argument. The positive upper argument scale
in (1) preserves this property.

Terminal two separates by real positivity. At terminal three,
the even polynomial is quadratic, so a common nonreal root would
imply divisibility. For the core's binomial-stripped coefficients
\(c_{j,k}\), the endpoint invariant
\(c_{j,1}c_{j,j-1}/(c_{j,0}c_{j,j})\) decreases strictly with \(j\)
by the [endpoint nondivisibility proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients).
The Gaussian multiplier multiplies this invariant by
\(\exp[-\tau(j-1)/j]\), so the strict decrease persists.
That proof uses only this decrease and positivity and is
invariant under positive argument rescaling. Thus there is
no reference contact at terminal three either.

## The finite union and the ordinary tail

Let \(D=\sum_i m_i d_i\), and put \(\Lambda=N/n\).
The core endpoint scale is
\(G_j=([z^j]Q_j/Q_j(0))^{1/j}\). Its adjacent drift obeys
\[
                        |\log(G_N/G_n)|\le D\log\Lambda.          \tag{2}
\]
For one row write \(g(x)=\log\Gamma(x+1)\) and
\(H(x)=(g(xb)-g(xa))/x\). Then \(xH'(x)\) is the
oriented integral of \(xu\,g''(xu)\) from \(a\) to \(b\).
For \(y>0\) the trigamma series gives
\(0\le y\,g''(y)\le1\), by comparison of the decreasing sum
with its integral. The inequality at zero follows by continuity.
Integrating \(|xH'(x)|\le|b-a|\) proves (2).

For \(0\le\tau\le12\), both an actual common root with
\(\eta\le1/2\) and a diffuse reference contact must have
\[
                            |h|\le D+N(n+N+24).                 \tag{3}
\]
Indeed the
[physical endpoint annuli and center comparison](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#uniform-approximation-of-the-physical-pair)
give
\[
            |h|\log\Lambda
                  \le (n+N+24)+|\log(G_N/G_n)|.
\]
Substitute (2) and use \(\log\Lambda\ge1/N\) to get (3).
Here \(h\) is the total signed gap
of the tail. Thus the union of the exceptional sets just proved,
over \(N<J\), (3) and \(\mathcal S\), is a finite union of
proper algebraic subsets. It is defined over \(\mathbb Q\).
Its integer grid has \(O((H+1)^{r-1})\) points, by the
nonzero polynomial grid bound.

Fix an ordinary core outside that union. For each terminal,
allowed twist and slope, the physical reference resultant is
nonzero on the compact segment \(0\le\tau\le12\).
The finite family has a positive absolute minimum for the fixed
core. This minimum has a computable positive lower bound:
on each segment the resultant is a rational Laurent polynomial in
\(t\) on \([\exp(-12/(2bnN)),1]\), so one can isolate its roots
and bound its values by rational intervals. By continuity on a
compact neighborhood, a computable positive tolerance \(\delta_C\)
in \(|\lambda-\sigma\tau|\) keeps these resultants bounded away
from zero. For instance, bound their first derivatives in the
upper scale coordinate \(v=\exp[-\lambda/(2nN)]\) on a rational
box containing the compact interval of \(v\)'s.

The same physical center comparison gives coefficient error
\(O_{C,J}(\eta\tau)\) between the exact tail augmentation and the
reference at its actual \(\tau,\lambda,h\). Resultants are
Lipschitz in coefficients on the bounded comparison boxes.
For a sufficiently small computable \(\eta_C\), this keeps
the ordinary resultants nonzero at the remaining low terminals.
The endpoint \(\tau=0\) was included in the compact segments, so
there is no positive lower curvature requirement.
Curvature \(\tau\ge12\) is covered by the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md),
and the absolute cutoff handles \(N\ge J\).
Rational Bézout identities give the characteristic-zero claim.

For the common-exponent example and any failure in the
low-terminal range, (3) bounds \(|\lambda|\) by a constant times
\(1/M\). A nonempty integer tail has \(1/M\le\eta\).
Reducing the fixed core's tolerance therefore brings every
possible failure into the covered neighborhood of slope zero.

## Prescribed algebraic curves and compact parameter sets

The incidence argument also works with an arbitrary prescribed
complex algebraic curve \(C_N\subset(\mathbb C^*)^2\) for each
\(4\le N<J\), or a finite union of such curves. In lieu of (1) use
the coefficient multipliers \(q^{N k(n-k)}\) and
\(q^{n k(N-k)}\Lambda^{hk}v^k\). For fixed \((q,v)\) the same
codimension-two bound applies. Allowing one algebraic dimension in
\((q,v)\) therefore gives the same proper exceptional set after
projection to the bases.

More precisely, let \(K\) be a compact subset of
\(\{0\le\tau\le12,\ |\lambda|\le\tau\}\) such that its image
\[
     (q,v)=\bigl(e^{-\tau/(2nN)},e^{-\lambda/(2nN)}\bigr)
\]
lies on \(C_N\) at each such terminal. For ordinary cores outside a
proper algebraic subset, the same minimum and perturbation argument
gives a neighborhood of \(K\) and a positive diffuseness tolerance.
Every ordinary tail with its parameters in that neighborhood and
with diffuseness below the tolerance gives all-degree SFC(2).
Both choices can depend on the fixed core and compact set.
When the compact set admits effective approximations with error
bounds, interval tests of its finite list of resultants give
computable tolerances.

If the curves are defined over \(\mathbb Q\), the exceptional set
can be taken over \(\mathbb Q\) directly. For arbitrary complex
curves, its rational points still lie in a proper rational algebraic
subset. Indeed, decompose any nonzero vanishing polynomial in a
\(\mathbb Q\)-basis of the finite-dimensional space spanned by its
coefficients. At a rational point where it vanishes, every rational
coefficient polynomial in that decomposition vanishes, including
at least one nonzero one. The same integer grid bound follows.

For example, fix any one value \(\tau_0\in[0,12]\). The full segment
\(\{(\tau_0,\lambda):|\lambda|\le\tau_0\}\) has image on the vertical
line \(q=e^{-\tau_0/(2nN)}\) at every terminal. A generic ordinary
core therefore has a diffuse stability neighborhood of that
entire segment, with the signed bias free to vary.
