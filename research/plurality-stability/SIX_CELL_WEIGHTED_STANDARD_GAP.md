# A directed standard-sector gap for six cells

The balanced standard block of the six-cell Gaussian simplex has
gap at least \(\ell/1000\) for \(0<\rho\le7/10\), where \(\ell\)
is the distance between two Gaussian cell moments. The other six-cell
blocks and the dimension-independent ambient lift give a full
balanced interface gap modulo rotations throughout the same range.
This is an interface theorem at the simplex.

The [directed checker](check_six_cell_weighted_gap.py) also
certifies the sufficient reversible-chain ratio above one at
\(\rho=1/5,1/2,13/20\).
The normalization, the ordinary-balance correction and the
uniform mass mismatch bound are in the
[reversible-chain note](SIX_CELL_REVERSIBLE_CHAIN.md).
Use \(\mu\) for ordinary conditional Gaussian measure on the
facet, \(p=\gamma_H(F)\) for its Gaussian probability, \(D\)
for the normalized row, \(\pi=D\mu/\int D\,d\mu\),
\(\gamma\) for the gap of the row-normalized chain on invariant
functions, and \(Z=(\int D\,d\mu)(\int D^{-1}d\mu)\) for the
mass mismatch. The chain ratio is \(\gamma/[(1-\rho)Z]\).

## A sharper mass mismatch bound at small correlation

There is a stronger bound when \(\rho\) is small:
\[
       Z-1<\frac{15}{8}\frac{\rho^2}{1-\rho^2}.
\]
Indeed \(Z-1=\frac12\mathbb E[(D(X)-D(Y))^2/(D(X)D(Y))]
\le\operatorname{Var}_\mu D\) for independent \(X,Y\sim\mu\).
Because \(D\) is permutation-invariant, the Gaussian Neumann
gap \(2\) from (SN10) gives
\(\operatorname{Var}_\mu D\le\frac12\int|\nabla D|^2d\mu\).
Use (FV3) with four gaps \(g_i=x_2-x_{i+2}\) and
\(s_i=\rho g_i/\sigma\). Its success probability is \(3pD\).
The derivative with respect to \(s_i\) is bounded by the marginal
density at \(-s_i\) of \(W=X_i+N/\sqrt2\), where \(X_i\) is
normal of variance \(3/2\) and \(N\) is independent Rayleigh.
That density on the nonpositive half-line is bounded by its
value \(3/(4\sqrt{3\pi})\) at zero, by the convolution integral
\(\int_0^\infty n e^{-n^2/2}\phi_{\sqrt{3/2}}(-s_i-n/\sqrt2)\,dn\).
All the gradient inner products of the \(g_i\) are nonnegative,
and the squared norm of their sum is 12. The gradient estimate
therefore gives
\[
 \operatorname{Var}_\mu D
       \le \frac9{8\pi(3p)^2}\frac{\rho^2}{1-\rho^2}
       <\frac{15}{8}\frac{\rho^2}{1-\rho^2}.
\]
The rational comparison uses \(\pi>3\) and \(p>747/5000\).
In particular \(Z<69/64\) at \(\rho=1/5\).

## Directed finite-rank certificate in six cells

The finite comparison extends (FV2)--(FV7) of the
[five-cell proof](FIVE_CELL_STANDARD_GAP.md).
The four gap coordinates are \(g_i=x_2-x_{i+2}\), \(1\le i\le4\),
with covariance \(C_0=I+\mathbf1\mathbf1^\top/2\) before conditioning
on the positive orthant. Their representation is
\(x_1=x_2=\sum_i g_i/6\) and \(x_{i+2}=x_2-g_i\).
For \(s_i=\rho g_i/\sigma\), the exact probability from (FV3) is
\(3pD=\Pr\{X_i+N/\sqrt2\ge-s_i,\ 1\le i\le4\}\), with
\(X\sim N(0,C_0)\) and independent Rayleigh \(N\). Conditioning on
the common Gaussian variable and \(N\) leaves independent coordinate
comparisons, all increasing functions of the same scalar. Applying
the increasing-function covariance inequality inductively shows
that the joint probability is at least the product of the marginals.
The same one-dimensional bound as (FV4) gives
\[
 D_\rho^{-1}\le E_\rho=3p\prod_{i=1}^4 (1+Z_i/3),\qquad
 Z_i=(e^{-s_i^2/3}+e^{-4s_i^2/9}+e^{-4s_i^2/3})/3 .
\]
The \(E_\rho\) here is the reciprocal multiplier envelope;
the mass mismatch \(Z\) defined above is a separate scalar. Use
\(E_0=1\), and use \(E_a\) for \(\rho\in[a,c]\) because the row
probability is increasing in \(\rho\).

In the interaction \(\sigma A_\rho=B_\rho/\rho\), ordinary balance
removes the constant ambient Hermite coefficient. All remaining
ambient coefficient forms are positive by (S11), so upper
coefficient weights on each correlation slab suffice. For the
finite expansion choose the orthogonal coordinate directions
\[
 (2,2,-1,-1,-1,-1),\quad(0,0,1,-1,0,0),\quad
 (0,0,1,1,-2,0),\quad(0,0,1,1,1,-3),\quad(1,-1,0,0,0),
\]
with squared lengths \(12,2,6,12,2\). Restricted to the facet
and its \((23)\)-transport, their linear coefficients on \(g\) are
respectively
\[
 \begin{split}
 &(1,1,1,1),\,(-1,1,0,0),\,(-1,-1,2,0),\,
                  (-1,-1,-1,3),\,(0,0,0,0),\\
 &(-2,1,1,1),\,(0,1,0,0),\,(0,-1,2,0),\,
                  (0,-1,-1,3),\,(1,0,0,0).
 \end{split}
\]
Permutational invariance reduces the adjacent contributions to
four copies of this transport. The self term has weight two. In the
basis of normalized monomial orbit averages (with constants omitted),
the exact common denominator is \(6\cdot12^m m!\) at ambient degree
\(m\). The analogous planar degree-\(m\) denominator is \(12^m m!\).
For example the first standard coefficient on the orbit of \(g_1\)
is \(48/(6\cdot12)=2/3\), while the planar coefficient is
\(16/12=4/3\). This agrees with the initial \(1/2\) ratio in (S77).
The same tail bound (FV5) from zero extension to the pair plane
has cap \(t=c^N/p\) at cutoff \(N\), with its exact finite planar
subtractions. It requires only pointwise invariance and ordinary
mean zero and holds in dimension four by the same argument.

For a slab constant \(f\in(0,1]\) and reserve \(\epsilon=1/1000\),
let \(K\) be the resulting rational finite kernel, so that
\(\sigma A_\rho\preceq K+tI\). The inverse lower envelope for
\(fD_\rho-t-\epsilon\) is \(e=E_a/(f-t-\epsilon)\), since
\(D_\rho\ge1\). Let \(v_i\) denote the finite monomial orbit features.
The Gram matrix in the weighted Hilbert space, projecting the
representing vectors onto the ordinary-balanced space, is
\[
    G_{ij}=\int ev_i v_j\,d\mu
          -\frac{\int ev_i\,d\mu\int ev_j\,d\mu}{\int e\,d\mu}.
\]
The real nonzero eigenvalues of \(QG\), for the rational
coefficient matrix \(Q\) of \(K\), are the eigenvalues of the
finite form in the weighted space. Cauchy--Schwarz or
congruence by the positive square root of \(G\) makes them
real even with the signed finite tail corrections. Therefore
\(\|K\|_{e^{-1}}\le(\operatorname{tr}((QG)^m))^{1/m}\) for
even \(m\). A strict value below one proves the claimed comparison
with at least \(\epsilon\) reserve. Multiplying out the Gram
denominator actually gives, for a certified norm upper \(n_*\),
\[
  \sigma A_\rho\preceq [\,t+n_*(f-t-\epsilon)\,]D_\rho
                 \quad\text{on ordinary-balanced functions.}
\]
This follows from \(E_a^{-1}\le D_\rho\) and \(D_\rho\ge1\).
When \(f=1\) it proves the full balanced standard-sector sign,
with gap at least \(\ell/1000\) in the facet norm as in (FV2).
The other six-cell sectors and the ambient lift are already
settled in the dependency graph.

The weighted moment calculation expands the product defining \(E_a\),
reducing it to 35 exchange-equivalence classes of diagonal Gaussian
tilts of \(C_0\). Each tilt changes the precision by the nonnegative
diagonal rates. The positive density ratio determinant and the
factor \(3p\) give the same cancellation of \(p\) as in (FV6).
For any tilt the multivariate polynomial moments use the integration
by parts recursion (FV7). Its only new base case is the four-dimensional
orthant probability. Standardize the covariance to the correlation
matrix \(R\), and integrate its derivative along \(I+t(R-I)\) from
\(t=0\) to \(t=1\). The derivative, proved by the same Fourier calculation
as for the equicorrelation mass above, is
\[
 \sum_{i<j}
   \frac{R_{ij}}{2\pi\sqrt{1-t^2R_{ij}^2}}\,
      \left(\frac14+\frac{\arcsin c_{ij}(t)}{2\pi}\right),
\]
where \(c_{ij}(t)\) is the off-diagonal correlation of the remaining
two coordinates conditioned on coordinates \(i,j\) equal to zero.
The Schur complement uses exact covariances and directed roots.
Acb integration of this analytic expression returns an enclosing
ball. The square-root branch checks use the analytic flag and the
arcsine input is verified inside the unit disk whenever Acb requests
analyticity. Covariances are positive definite along the real path.

For a pointwise or interval certificate of the stationary
comparison, choose \(f=1-(Z_0-1)(1/a-1)\) on a slab starting at
\(a>0\) for a positive target multiplier \(Z_0\). More generally
any verified interaction factor \(b_*=t+n_*(f-t-\epsilon)<1\)
gives \(\gamma\ge1-\rho b_*\). Indeed, subtract the ordinary
mean of \(h\) to get \(h_0\) and use exactly the same Dirichlet
form. Since \(\int D h_0^2d\mu\ge \tau\operatorname{Var}_\pi h\),
its \(1-\rho b_*\) margin bounds the stationary gap. Dividing
by the upper bound on \(Z\) above gives the rigorous lower
ratio in each receipt. On a slab this uses the left endpoint
for \((1-\rho b_*)/(1-\rho)\), which is increasing in \(\rho\),
and the right endpoint for the mass upper bound.

## Certificate and replay

The [receipt](six_cell_weighted_gap_verification.json) retains
the rational cover, the stationary-ratio points and their strict
Arb trace bounds. The five adjoining slabs for the balanced gap use the following
ambient degrees and strict upper trace bounds.

| Slab | Degree | Trace power | Upper trace |
|---|---:|---:|---:|
| \([0,1/3\)] | 6 | 2 | 0.343297 |
| \([1/3,1/2\)] | 8 | 4 | 0.734444 |
| \([1/2,3/5\)] | 10 | 4 | 0.837674 |
| \([3/5,13/20\)] | 12 | 4 | 0.831743 |
| \([13/20,7/10\)] | 14 | 8 | 0.814791 |

All estimates use
512-bit Arb, an exact monomial coefficient expansion and
\(\epsilon=10^{-3}\). The checker requires each even normalized
trace to be strictly below one. The finite matrix includes the
ordinary-balance correction in its Gram form.

Exact rational comparison points give the following outward-rounded
lower bounds.

| Correlation | Degree | Strict lower bound for the stationary ratio |
|---|---:|---:|
| \(1/5\) | 7 | \(1.04904\) |
| \(1/2\) | 10 | \(1.06910\) |
| \(13/20\) | 14 | \(1.02122\) |

Run from the repository root with the commands recorded in the
receipt. The tests at a single value of \(\rho\) certify the
stationary ratio at that exact rational point. The stationary-ratio bound also follows on the entire first
slab \(0<\rho\le1/3\). Its normalized interaction upper bound
is below \(0.5892\). The rational estimates \(\pi>25/8\)
from the inscribed regular 24-gon and \(p>747/5000\) sharpen the
gradient estimate to \(Z-1<(9/5)\rho^2/(1-\rho^2)\).
For \(0<\rho\le59/200\),
\[
             \frac{9\rho}{5(1+\rho)}<1-0.5892.
\]
For \(59/200\le\rho\le1/3\), use \(Z<1.17\) and
\[
          \frac{\rho(1-0.5892)}{1-\rho}>0.17 .
\]
In each case \(1-0.5892\rho>(1-\rho)Z\). The retained directed
point receipts still use the coarser bound \(15/8\) where
applicable.
