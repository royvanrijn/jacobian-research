# Free-mass minimizers in the negative voting interval

Let \(f=(f_1,\ldots,f_k)\) be a measurable simplex-valued function on
\((\mathbb R^d,\gamma_d)\), where \(p_i=\int f_i\,d\gamma_d\), \(u=1/k\),
and \(\delta=\|p-u\mathbf 1\|_2\). Use the same Ornstein–Uhlenbeck
operator and noise stability \(S_{-r}\) as in the
[negative voting interval](NEGATIVE_VOTING_INTERVAL.md). Write
\(b_i=\int x f_i(x)\,d\gamma_d(x)\) and \(W_1(f)=\sum_i\|b_i\|^2\).

Throughout \(0<r\le1/(k-1)\), all sufficiently large \(k\), and
\(d\ge k-1\), the following estimates hold with absolute constants.

* Every \(f\) satisfying \(S_{-r}(f)\le S_{-r}(T_k)\) has
  \[
   \|p-u\mathbf 1\|_2=O(r\log k/\sqrt{k}),
   \qquad W_1(f)\le F_k+O(r(\log k)^2/k).
   \tag{M1}
  \]
  In particular, each mass differs from \(1/k\) by at most
  \(O(\log(k)/k^{3/2})\) at the endpoint of the voting interval.
* A global minimizer over fractional assignments exists in each
  fixed finite dimension. Every minimizer is an ordinary partition
  up to Gaussian null sets. Its masses satisfy
  \[
   \max_i |p_i-1/k|=O(r\log k/k).
   \tag{M2}
  \]
  This is \(O(\log(k)/k^2)\) at the endpoint. The existence and purity
  assertions hold for all \(k\ge2,\ d\ge1,\ 0<r<1\), outside the
  asymptotic domain of the mass estimate.

## Rebalancing possible improvers

The [full-interval estimates (V3a), (V6)](NEGATIVE_VOTING_INTERVAL.md#v4-first-moments-of-possible-improving-assignments)
give absolute \(b,B>0\) with \(b/k\le p_i\le B/k\) and
\[
 \delta^2/r \le W_1(f)-F_k+C r/k
 \tag{M3}
\]
for every possible improver. The halfspace first-moment comparison
in the same note and the normal Mills bounds give
\(\|b_i\|\le \phi(\Phi^{-1}(p_i))\le C p_i\sqrt{\log k}\).
Indeed \(t=\Phi^{-1}(1-p_i)>0\), \(t\le\sqrt{2\log(1/p_i)}\),
and \(\phi(t)\le p_i(t+1/t)\). Here \(t\) is bounded away from zero
for large \(k\). In particular \(\max_i\|b_i\|\le C\sqrt{\log k}/k\).

For each surplus label \(i\), remove the fraction
\(\alpha_i=(p_i-u)/p_i\) of its density. Distribute the resulting
nonnegative density among the deficit labels in the proportions
of their deficits. This gives a balanced fractional assignment \(g\).
Writing \(b_i(g)\) for its moments, the triangle inequality and
the proportions summing to one give
\[
 \sum_i\|b_i(g)-b_i\|
 \le2\sum_{p_i>u}\alpha_i\|b_i\|
 \le C\sqrt{\log k}\,\|p-u\mathbf 1\|_1 .
\]
For \(p=u\mathbf1\), simply take \(g=f\). The balanced first-moment
inequality [in its fractional form](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#external-endpoint-and-application-audit)
gives \(W_1(g)\le F_k\). Expanding the squared moments and dropping
the nonnegative squared differences gives
\[
 W_1(f)-F_k\le
 2\sum_i \|b_i\|\,\|b_i-b_i(g)\|
 \le C\frac{\log k}{\sqrt{k}}\,\delta .
 \tag{M4}
\]
Combining (M3) and (M4) with the quadratic inequality
\(\delta^2\le C r(\log k)\delta/\sqrt{k}+C r^2/k\)
proves the first assertion of (M1). Substitution into (M4) gives the
second assertion.

## Existence, purity, and pointwise optimality

For fixed \(d\ge1,\ k\ge2,\ 0<r<1\), the feasible assignments are a
bounded, closed convex subset of \(L^2(\gamma_d)^k\), so they are
weakly compact. In finite \(d\) the degree-\(m\) Hermite eigenspace is
finite dimensional and the eigenvalues \((-r)^m\) tend to zero. Thus
\(T_{-r}\) is compact, and the quadratic objective is continuous on
weakly convergent sequences of feasible assignments. A minimizing
sequence has a weakly convergent subsequence attaining the infimum.

At a minimizer, transfer any small fraction of \(f_i\) into \(f_j\)
on a measurable set where \(T_{-r} f_i>T_{-r} f_j\).
The first variation is strictly negative if \(f_i\) is positive
there on a set of positive Gaussian measure. Consequently
\[
 f_i(x)>0\quad\Longrightarrow\quad
 (T_{-r}f_i)(x)=\min_j (T_{-r} f_j)(x)
 \quad\text{for almost every }x .
 \tag{M5}
\]

Suppose two labels \(i,j\) overlap on a set of positive measure. There
is \(\varepsilon>0\) such that \(E=\{f_i,f_j\ge\varepsilon\}\) has
positive Gaussian measure. Pick distinct points \(x,y\) in the support
of Gaussian measure restricted to \(E\). For sufficiently small
\(\eta>0\), let \(h_1,h_2\) be the indicators of \(E\cap B(x,\eta)\)
and \(E\cap B(y,\eta)\), each divided by the Gaussian measure of its
support. For \(x_1=x,\ x_2=y\), the matrix
\[
 G_{\eta,ab}=\langle h_a,T_{-r} h_b\rangle
\]
converges to \(G_{ab}=M_{-r}(x_a,x_b)\), by continuity of the Mehler
density relative to \(\gamma_d\otimes\gamma_d\). Its entries satisfy
\[
 \frac{M_{-r}(x,x)M_{-r}(y,y)}{M_{-r}(x,y)^2}
  =\exp\!\left(-\frac{r|x-y|^2}{1-r^2}\right)<1.
\]
Thus \(G\) and all sufficiently close \(G_\eta\) have negative
determinant. Choose real coefficients with
\(h=c_1h_1+c_2h_2\) such that \(\langle h,T_{-r}h\rangle<0\).
Both signs of a sufficiently small transfer \(f_i\mapsto f_i+t h\),
\(f_j\mapsto f_j-t h\) are feasible since \(h\) is bounded and
supported in \(E\). Its linear term vanishes by (M5), and its
quadratic term is \(2t^2\langle h,T_{-r}h\rangle<0\), a contradiction.
No two labels overlap on a positive-measure set, which proves purity.

## A per-label mass bound at a minimizer

For a global minimizer in the large-alphabet voting regime, (V3a)
applies since the simplex competes. In particular \(p_i\asymp1/k\)
uniformly. Let \(s=1/\log k\), and let \(P_m\) denote the orthogonal
projection onto Hermite degree \(m\). For \(r\le s\), the exact
Gaussian expansion and Cauchy–Schwarz give
\[
 \begin{aligned}
 \left|\langle f_i,T_{-r}f_j\rangle-p_ip_j\right|
 &\le (r/s)\sum_{m\ge1}s^m\|P_m f_i\|_2\|P_m f_j\|_2\\
 &\le (r/s)\|T_{\sqrt s}f_i\|_2\|T_{\sqrt s}f_j\|_2 .
 \end{aligned}
\]
Each squared norm equals \(\langle f_i,T_s f_i\rangle\).
The functional Borell comparison at positive correlation bounds this
by \(J_s(p_i,p_i)\), including fractional \(f_i\), as in
[(V4), (V6)](NEGATIVE_VOTING_INTERVAL.md#v2-hermite-remainder-at-the-voting-scale).
The small-correlation bivariate normal density bound used there gives
\(J_s(p_i,p_i)\le C/k^2\) when \(p_i\le B/k\).
The total error is therefore \(O(r\log(k)/k^2)\).

By (M5), \(0\ge\langle f_i,T_{-r}(f_i-f_j)\rangle\), so
\[
 p_i(p_i-p_j)\le C r\log(k)/k^2 .
\]
Divide by \(p_i\ge b/k\) and exchange \(i,j\). The full spread of the
mass vector is at most \(C r\log(k)/k\), and its average is \(1/k\).
This proves (M2).

## Matched-mass pairs

For the paired functional
\[
 P_r(f,g)=\sum_i\langle f_i,T_{-r}g_i\rangle ,
 \qquad \int f_i\,d\gamma_d=\int g_i\,d\gamma_d ,
\]
a fractional minimizer exists by the same weak compactness argument.
After reflecting the second assignment, this is Heilman's paired
problem at positive correlation. [Lemma 7.3](https://arxiv.org/html/2209.11216#S7)
supplies existence of an ordinary minimizing pair. In fact every
fractional minimizer is ordinary, even if their common masses are
held fixed.

If \(f_i,f_j\) overlap, choose \(\varepsilon>0\) such that
\(E=\{f_i,f_j\ge\varepsilon\}\) has positive measure. All bounded
zero-mean transfers \(h\) between these labels supported on \(E\)
are admissible in both signs for \(f\) alone, so
\(\langle h,T_{-r}(g_i-g_j)\rangle=0\). Thus \(T_{-r}(g_i-g_j)\)
is constant almost everywhere on \(E\). The Gaussian smoothing of
a bounded function is real analytic: viewing the shift in the
conditional Gaussian kernel as complex, every derivative converges
locally uniformly by domination against a Gaussian times an
exponential linear in its integration variable. A real analytic
function on \(\mathbb R^d\) which vanishes on a set of positive
Lebesgue measure vanishes identically, as follows locally by
induction on dimension from the one-dimensional identity theorem
and Fubini's theorem. Hence \(T_{-r}(g_i-g_j)\) is globally constant,
and injectivity on Hermite degrees gives
\(g_i-g_j=p_i-p_j\) almost everywhere. The matched masses of both
labels are positive, so \(g_i,g_j\) also overlap on a set
\(E'=\{g_i,g_j\ge\varepsilon'\}\) of positive measure.

Pick a nonzero bounded zero-mean function \(h\) supported on \(E\).
Its OU image is nonconstant on \(E'\), by analyticity and
injectivity. Choose a bounded zero-mean function \(j_0\) on \(E'\)
such that \(\langle h,T_{-r}j_0\rangle\ne0\), for example the
restriction of \(T_{-r}h\) to \(E'\) minus its average there.
The simultaneous transfers \((h,-h)\) on \(f\) and
\((j_0,-j_0)\) on \(g\) preserve both mass vectors, so both linear
terms vanish. Choose the sign of \(j_0\) to make the quadratic
term \(2t^2\langle h,T_{-r}j_0\rangle<0\), contradicting minimality.
This proves the claim for \(f\), and the same argument with the
roles interchanged proves it for \(g\).

The mass estimate (M2) also holds for a global pair minimizer when
the common masses are free. First, the functional Borell comparison
gives \(P_r(f,g)\ge\sum_i J_{-r}(p_i,p_i)\), with the same
probabilities \(p_i\) on both sides. Comparing this with
\(P_r(T_k,T_k)\) proves \(p_i\asymp1/k\) by the same (V1)–(V3a)
argument. Next let \(\lambda\in\mathbf1^\perp\) be the multiplier
for the equality of masses in the linear first variation.
For completeness, minimize the linear functional
\[
 (F,G)\longmapsto \sum_i
 \big(\langle F_i,T_{-r}g_i\rangle+\langle G_i,T_{-r}f_i\rangle\big)
\]
over all fractional pairs having a given difference of mass vectors.
Its value is a finite convex function of the difference vector in a
neighborhood of zero. Indeed small centered difference vectors are
realized by constant assignments near \(k^{-1}\mathbf1\), and the
functional is uniformly bounded. A subgradient at zero, with sign
reversed, is the required \(\lambda\). The first-order variational
inequality states that the original pair minimizes this linear
objective at mass difference zero. After adding \(\lambda\) times
the mass difference, the original pair minimizes the objective over
all fractional pairs, with no remaining coupling constraint. The
pointwise selection condition therefore holds for scores
\(T_{-r}g_i+\lambda_i\) on \(f\) and \(T_{-r}f_i-\lambda_i\) on \(g\).

Testing the respective \(i\)-versus-\(j\) conditions against \(f_i\)
and \(g_i\) and adding cancels
\(p_i(\lambda_i-\lambda_j)\), giving
\[
 0\ge \langle f_i,T_{-r}(g_i-g_j)\rangle+
       \langle g_i,T_{-r}(f_i-f_j)\rangle .
\]
The zeroth-degree term is \(2p_i(p_i-p_j)\). The higher terms
satisfy the same \(O(r\log(k)/k^2)\) bound as in the proof of
(M2), by applying the cross-term Cauchy–Schwarz estimate. This
establishes (M2) for pairs.

## Approximate selection by affine scores

For a same-partition minimizer define affine scores
\(A_i(x)=b_i\cdot x-(p_i-u)/r\). Then the exact objective's
selection condition gives the expected affine score loss bound
\[
 \int \left(\max_i A_i(x)-\sum_i f_i(x)A_i(x)\right)
 \,d\gamma_d(x)
 \le C r(\log k)^2/k .
 \tag{M6}
\]
For a matched-mass minimizing pair, the same estimate holds on
each side, with the opposite side's moment vectors and prices
\((p_i-u+\lambda_i)/r\) on \(f\) and
\((p_i-u-\lambda_i)/r\) on \(g\).

To check the bound, put
\[
 R_i(x)=r^{-1}\big(T_{-r} f_i(x)-p_i+r b_i\cdot x\big).
\]
For \(q=\log k>2\), \(t=(q-1)^{-1/2}\), and \(s=1/\log k\),
one has \(r/t\le \sqrt{s}\) on the voting interval for large \(k\).
Factor \(T_{-r}=T_t T_{-r/t}\) and remove the first two degrees.
By Nelson's Gaussian hypercontractivity
[as stated in Gross, Theorem 1.1](https://doi.org/10.1007/BF02392573),
\(T_t\) maps \(L^2(\gamma_d)\) contractively into \(L^q(\gamma_d)\).
The Hermite expansion and the previously established bound
\(\|T_{\sqrt s}f_i\|_2=O(1/k)\) give
\[
 \|R_i\|_q\le
 \frac1r\left(\frac{r}{t\sqrt s}\right)^2
 \|T_{\sqrt s}f_i\|_2
 \le C r(\log k)^2/k .
\]
The expectation of the maximum of their absolute values
is at most \(k^{1/q}\max_i\|R_i\|_q=e\max_i\|R_i\|_q\).
The actual score \(T_{-r}f_i\) equals \(u-rA_i+rR_i\).
On any label selected by \(f\), (M5) makes its actual score minimal,
so its affine loss against the largest affine score is at most
\(2\max_i|R_i|\). Integration proves (M6). The paired version
uses the same factorization with the other assignment's remainders
and the two price shifts from the paired first variation.

If \(\mathrm{gap}_A\) is the largest affine score minus the second
largest, any measurable selector for the largest affine score
differs from the minimizing assignment on measure at most
\[
 \gamma_d\{\mathrm{gap}_A\le h\}+C r(\log k)^2/(kh)
 \qquad(h>0).
\]
Indeed, outside the small-gap set, an incorrect label costs at
least \(h\) in (M6). Thus the expected score bound reduces
closeness to affine cells to an explicit affine gap estimate.
