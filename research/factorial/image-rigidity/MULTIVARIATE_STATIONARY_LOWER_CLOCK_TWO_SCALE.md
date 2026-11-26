# Two scales for a same-sign lower-order clock response

Fix adjacent orders \(n,N=n+1\), \(n\ge2\), an ordinary finite
two-monomial core and a simple phase contact at positive Gaussian
strength \(W_*\) and signed bias \(\lambda_*\).
Use the [refined phase expansion](MULTIVARIATE_STATIONARY_ZERO_CLOCK_DICHOTOMY.md#the-refined-phase),
with \(g(W_*+\delta)=a\delta^{2r}+O(\delta^{2r+1})\), \(a\ne0\).
Assume the cubic response satisfies
\(A(W_*+\delta)=b\delta^r+O(\delta^{r+1})\), allowing \(b=0\),
and the clock response satisfies
\[
   B(W_*+\delta)=c_s\delta^s+O(\delta^{s+1}),\qquad
             0<s<r,\qquad s\ {\rm even},\qquad ac_s>0 .
\]
Thus the lower-order clock term has the same phase sign as the
even Gaussian contact. The opposite sign and odd clock orders
have their own unfolding tests in the cited note.

Write \(\sigma=\operatorname{sgn}(a)\), and let
\(\mathcal C=\mathcal C(W_*,\lambda_*)\) be that note's
[exact compact moment set](MULTIVARIATE_STATIONARY_ZERO_CLOCK_DICHOTOMY.md#the-exact-compact-moment-set).
It records \(P=K/\sqrt{B_4}\), \(Q=Z/\sqrt{B_4}\), and
\(\mu=V/B_4\). Set
\[
 \begin{aligned}
 L(P,Q,\mu)&=\sigma(C+E\mu+D_{11}P^2+D_{12}PQ+D_{22}Q^2),\\
 \chi&=\frac{s}{2r}-\frac13,\qquad
 L_0=\min_{\mathcal C} L(P,Q,\mu).
 \end{aligned}
\]
For the second scale take \((P,0,\mu)\in\mathcal C\) and
\(X\in\mathbb R\). If \(\chi<0\) impose \(\mu=0\).
On this domain define
\[
 L_{\rm slow}(X,P,\mu)=
     |a|X^{2r}+\sigma bP X^r+L(P,0,\mu)
    +\begin{cases}
         |c_s|X^s|\mu|^{2/3},&\chi=0,\\
         0,&\chi\ne0 ,
      \end{cases}
\]
and let \(L_1\) be its minimum. This domain is nonempty
and \(L_1\) is attained since the moment set is compact
and the polynomial is uniformly coercive in \(X\).

**Theorem.** If \(L_0,L_1>0\), the indicated angular contact
is excluded for all sufficiently diffuse ordinary shared tails
in a neighborhood of the reference. The tolerances are uniform
in tail row count, gaps and heights. This includes saturated bias.
If \(\min(L_0,L_1)<0\) and \(|\lambda_*|<W_*\), the reference
is a limit of simple transverse interior physical contacts
of finite ordinary anchored cores, preserving the adjacent
orders, the integer signed gap and every strict total curvature
cap. Both strict conditions have finite semialgebraic certificates
with enclosures of the reference constants. The zero minimum
case remains outside the theorem.

## A necessary two-scale limit

Let \(\beta=B_4\) for a nonempty ordinary tail.
At any putative angular contact with
\((W,\lambda,\eta)\to(W_*,\lambda_*,0)\), one has
\(X=\delta/\beta^{1/(2r)}=O(1)\). To see this, use
\(\sigma B(W)Z\ge0\) for small real \(\delta\),
\(K,Z=O(\sqrt\beta)\), and the uniform expansion
\[
 \begin{split}
 0={}&\sigma g(W)+\sigma A(W)K+\sigma B(W)Z\\
      &+\beta L(P,Q,\mu)+o(\beta).
 \end{split}
\]
The contact implies
\(|\delta|^{2r}=O(\sqrt\beta|\delta|^r+\beta)\), which gives
the bound. It also bounds the nonnegative clock contribution
\(\sigma B(W)Z/\beta\). Pass to a subsequence with limits for
\(X,P,Q,\mu\). If \(X\to0\), the other Gaussian and cubic
terms divided by \(\beta\) tend to zero. Then \(L_0>0\)
contradicts the nonnegative clock contribution.

If \(X\) tends to a nonzero limit, the clock contribution is
asymptotic to
\[
              |c_s|X^s\,Q\,\beta^{(s-r)/(2r)}
\]
in the sense of a relative \(1+o(1)\) error on the whole term.
Its boundedness implies \(Q\to0\). Independently the ordinary
row values \(x_i=d_i/m_i\) satisfy
\[
                         |V|=\left|\sum_i x_i^3\right|
                                \le Z^{3/2}.
\]
Thus the clock contribution is bounded below, up to its relative
\(1+o(1)\) error, by
\[
                   |c_s|X^s|\mu|^{2/3}\beta^\chi.            \tag{1}
\]
If \(\chi<0\), its boundedness forces \(\mu\to0\).
If \(\chi=0\), the right side of (1) yields the additional
nonnegative term of \(L_{\rm slow}\).
If \(\chi>0\), nonnegativity suffices. In all cases a strict
positive \(L_1\) excludes the remaining limit. The inequalities
and the finite moment set give computable tolerances. For example,
with certified positive margins, one first chooses a small fixed
threshold for \(|X|\) using the \(L_0\) margin. At \(|X|\) above
that threshold, the bounded clock contribution gives explicit
bounds tending to zero for \(Q\) and, when \(\chi<0\), \(|\mu|\).
Compact semialgebraic enclosures of these shrinking sections
certify the \(L_1\) margin.

## Ordinary anchors realizing the negative limits

The anchors below have limiting strength \(tW_*\) and bias
\(t\lambda_*\), where \(t<1\) is close enough to one to preserve
a strict negative value. In the three matrix triples describing
\(\mathcal C\), this replaces \((A_q,B_q,C_q)\) by
\((tA_q,\sqrt t B_q,C_q)\). Consequently \(P,Q\) scale by
\(\sqrt t\) and \(\mu\) stays fixed. The remaining strength
tends to \((1-t)W_*\) and the remaining bias required by radial
matching tends to \((1-t)\lambda_*\), which is strictly feasible.
The total signed gap is obtained by subtracting the exact anchor
gap from the original integer tail gap.

Suppose first \(L_0<0\). Choose a feasible negative point with
\(Q>0\), possible by mixing its matrix data with data of positive
\(Q\) in an arbitrarily small proportion. To produce the latter,
allocate some positive mass at either \(q=1\) or \(q=-1\) with
\(C_q>0\) and choose \(B_q\) nonzero of the matching sign, keeping
\(B_q^2\le A_q C_q\). This respects the two mass constraints by
starting from \(A_+=(W_*+\lambda_*)/2\),
\(A_-=(W_*-\lambda_*)/2\), \(A_0=0\).
Use the exact ordinary realization of the scaled moment point.
With these anchor coefficients held fixed, vary the remaining
Gaussian strength so \(W-W_*=\beta^{1/(2s)}X\).
The signed phase divided by \(\beta\) converges with one real
derivative on bounded \(X\)-intervals to
\[
                    L(P,Q,\mu)+|c_s|Q X^s .
\]
Here and below the moments are those after the fixed scaling
by \(t\). The Gaussian and cubic terms vanish at this scale
since \(r>s\). The limit has a simple positive zero, so the
exact anchored references cross transversely nearby.

Now suppose \(L_1<0\) and the preceding case does not apply.
Choose a negative point with \(X\ne0\) and \(Q=0\).
The sign constraints on the matrices imply \(B_+=B_-=0\).
One can replace their third coordinates by
\(C_{\operatorname{sgn}\mu}=|\mu|\) and zero for the other
sign, then take \(C_0=1-|\mu|\).
When \(\mu=0\) give both signed triples third coordinate zero.
These changes preserve \(P,Q,\mu\) and feasibility because
\(C_0\) can only increase. Decompose the triples as in the
refined phase proof. The signed triplet with \(C_q=|\mu|>0\)
now has one pure third-coordinate contribution, which can be
realized by a single row with signed unit gap. With the
realization parameter \(\varepsilon\) and
\(\beta/\varepsilon^2\to1\), this row has
\(Z/\beta^{2/3}\to|\mu|^{2/3}\).

All other contributions to \(Z\) can be made
\(O(\beta^\nu)\) for any prescribed \(2/3<\nu<1\) as
follows. For a gap-zero-label triple (label \(q=0\)) with
positive mass, the realization has one actual row with
\(x=O(\varepsilon)\), so \(Z=O(\varepsilon^2)\).
For a pure third-coordinate triple at \(q=0\), choose integer
gap \(d=\lceil\varepsilon^{-\kappa}\rceil\) with
\(0<\kappa<1\) and \((2+\kappa)/3>\nu\).
Then use \(x\sim(C_0'\varepsilon^2/d)^{1/3}\), where
\(C_0'\) is that triple's third coordinate. Midpoint rounding
preserves this rate, its curvature tends to zero, and its
contribution to \(Z\) is \(O(\varepsilon^{2(2+\kappa)/3})\).
Any signed pure mass triple uses the unit cloud with
\(|x|=O(\varepsilon^2)\) and \(Z=O(\varepsilon^2)\).
The correction to the desired integer signed gap can grow. At
sufficiently large heights its contribution is negligible, as in
the refined phase proof.

At the slow scale \(W-W_*=\beta^{1/(2r)}X\) the signed
phase divided by \(\beta\) now converges with one derivative
to \(L_{\rm slow}(X,P,\mu)\). For \(\chi>0\) the single unit
row's clock term tends to zero. For \(\chi=0\) its coefficient
tends to \(|c_s||\mu|^{2/3}\). For \(\chi<0\) the chosen
domain has \(\mu=0\), and one takes
\(\nu>1-s/(2r)\) to make the clock term tend to zero.
In all cases the fixed log-coefficient corrections are analytic
in the varying Gaussian strength after the strength \(W_{\rm A}\)
of the anchors has been removed, and the derivative convergence
follows from (2)–(4) of the refined phase expansion.

At the selected argument the limiting polynomial is strictly
negative. On the same side of zero and further from zero it
is strictly positive, by its positive leading coefficient.
Choose two finite such values of \(X\). Their opposite signs
persist for sufficiently small anchors along the realizing
sequence. The [ordinary regularization of sign changes](MULTIVARIATE_MULTIPLE_CONTACT_ANGULAR_RANK_REVERSAL.md#ordinary-perturbations-with-simple-transverse-angular-contacts)
allows arbitrarily weak additional integer rows to preserve
both endpoint signs and make the intervening crossing transverse.
In the lower and upper isolating root disks take perturbations
small enough also to preserve strict feasibility of the radial
bias. This produces the claimed finite cores even if the
limiting polynomial had a higher odd-order zero. The center
matching and the curvature cap follow from the same argument
as in the [zero-clock unfolding](MULTIVARIATE_STATIONARY_ZERO_CLOCK_DICHOTOMY.md#exclusion-and-transverse-unfolding).

The two minima use compact semialgebraic moment data and one
real clock variable, which can be cut off effectively by
coercivity. In the critical case \(|\mu|^{2/3}\) is represented
by \(v\ge0,\ v^3=\mu^2\). Strict margins therefore admit
finite certificates from rational enclosures of the reference
constants by [real quantifier elimination](https://arxiv.org/pdf/1609.02879). These conclusions classify an additional stationary
regime at every adjacent degree pair. The transverse contacts
from the negative branch still require exact arithmetic
exclusion for arbitrary ordinary completions.

