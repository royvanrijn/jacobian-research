# All exterior hooks from a truncated Gaussian determinant

**Theorem (All higher exterior hooks).** For every \(k\ge4\),
\(d\ge k-1\), and \(0<\rho<1\), comparison (S91f) holds on the ordered
facet, strictly in the interior in every degree \(2\le p\le k-2\).
In particular, for every \(3\le p\le k-2\), the full isotypic block \(V_{(k-p-1,1^{p+1})}\) has
the ambient-dimension-independent bound
\[
 \mathfrak q_\rho(f)\ge
 \frac{\ell}{\sqrt{1-\rho^2}}
 \bigl(1-\rho^{p(p-1)/2-1}\bigr)\|f\|^2 .
 \tag{EH1}
\]
Strict comparison also gives a positive gap at \(p=2\).
The [interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md) constructs
the positive harmonic fields and moment formulas in (S91) and
combines the adjacent conditional terms in (S92a). Equations (S73)
and (S92) give separate proofs for the first two of these degrees.

## Truncated Gaussian moment derivatives

Fix a finite truncation level \(u\). For a real mean \(m\), let \(P_m\),
\(p_m\) and \(F_m\) denote, respectively, \(\Phi(u-m)\), the density of
\(N(m,1)\) conditional on being at most \(u\), and its CDF. Let
\(f(m)=\phi(u-m)/P_m\) and \(M_r(m)=\int y^r p_m(y)\,dy\).
For \(y<u\) set
\[
 w(m,y)=-\partial_m F_m(y)
        =p_m(y)-f(m)F_m(y).
\]
Set \(w(m,y)=0\) for \(y\ge u\). Differentiation under the Gaussian
integral or integration by parts gives
\[
 M_r'(m)=\int_{-\infty}^u r y^{r-1} w(m,y)\,dy,\qquad
 \lim_{y\uparrow u}\frac{w(m,y)}{u-y}=f'(m).
 \tag{EH2}
\]
Indeed \(w(m,u^-)=0\) and \(\partial_y w(m,u^-)=-f(m)(u-m+f(m))
=-f'(m)\). The Gaussian tails remove the lower integration endpoint.
Since \(m-f(m)=M_1(m)<u\), the expression
\(f'(m)=f(m)(u-m+f(m))\) is positive. Differentiating the conditional
mean in the exponential family also gives
\(1-f'(m)=M_1'(m)=\operatorname{Var}_{p_m}(Y)\ge0\).

For strictly decreasing \(x_1,\ldots,x_j\) write
\(\Delta(x)=\prod_{a<b}(x_a-x_b)\). We prove the following difference
bound for every \(p\ge2\) and every \(t_1>\cdots>t_p\):
\[
 0\le s_p\det[(f'(t_i),M'_1(t_i),\ldots,M'_{p-1}(t_i))]_{i=1}^p
 \le f'(t_1)(p-1)!\Delta(t_2,\ldots,t_p),
 \qquad s_p=(-1)^{(p-1)(p-2)/2}.
 \tag{EH3}
\]

## Derivative kernel of the conditioned Brownian flow

Start standard one-dimensional Brownian motion at \(m\) at time
zero and condition on its time-one endpoint being at most \(u\).
For \(t<1\) set \(h_t(x)=\Phi((u-x)/\sqrt{1-t})\).
The conditioned process \(X^m\) has generator
\(\frac12\partial_x^2+b_t(x)\partial_x\), with
\(b_t=\partial_x\log h_t\). This follows by applying the heat equation
for \(h_t\) to its transition densities
\(\phi_{\sqrt{t-s}}(y-x)h_t(y)/h_s(x)\) for \(s<t<1\).
By the derivative bound just proved,
\(-1/(1-t)\le \partial_x b_t(x)\le0\). In particular on every
\([0,T]\), \(T<1\), the drift is smooth, globally Lipschitz in space,
and of at most linear growth.

Using the same driving Brownian path for different means,
differentiation in the initial value gives
\[
 \partial_m X^m_T=J^m_T
   =\exp\left(\int_0^T\partial_x b_t(X^m_t)\,dt\right)\in[1-T,1].
\]
For the derivative formula, subtract the two integral equations with
the same additive noise, divide by the initial separation, and pass
to the linear scalar equation using the bounded drift derivative.
For a bounded smooth test function with compactly supported
derivative, \(J^m_T\le1\) justifies differentiating its expectation.
Thus the weighted transition measure
\(\mathbb E[J^m_T\mathbf 1_{\{X^m_T\in dy\}}]\) has continuous
density
\[
 w_T(m,y)=-\partial_m\int_{-\infty}^y
              \phi_{\sqrt T}(z-m)\frac{h_T(z)}{h_0(m)}\,dz.
 \tag{EH4}
\]
This follows by integrating the differentiated density against the
test function and integrating by parts in \(y\). Dominated
convergence in (EH4) gives \(w_T(m,y)\to w(m,y)\) as \(T\uparrow1\)
for \(y<u\) and a zero limit for \(y>u\): for fixed \(m\) and
\(T\in[1/2,1)\) the integrand and its \(m\)-derivative are bounded
by an integrable Gaussian times a polynomial.

For decreasing starting points \(x_i\) and decreasing endpoints
\(y_i\), \(\det[w_T(x_i,y_j)]\) is the density of independent
conditioned processes at the ordered endpoints, restricted to no
collision up to time \(T\), with path weight \(\prod_i J_T^{x_i}\).
Indeed, expand the determinant as the signed sum over endpoint
permutations. On paths that meet, swap the futures of the first
meeting pair, breaking ties lexicographically from the histories
at that time. The common generator and the strong Markov
property preserve the law,
and the sum of the path integrals in the weight is preserved.
The swapped endpoints change permutation sign. On paths with no
collision only the identity endpoint permutation is
possible. This is the Karlin--McGregor cancellation, with the
time-dependent drift and the product weights included explicitly.
The classical Brownian case is discussed in the
[retained Grabiner paper](sources/grabiner-weyl-chamber-1999.pdf),
pp. 179--180 and Section 4.

The determinant is strictly positive for distinct ordered endpoints.
The path weights are bounded below by \((1-T)^j\). Given endpoints,
the joint path law is that of independent Brownian bridges because
the conditioning factor \(h_T(y)/h_0(x)\) depends only on the endpoints.
Small uniform tubes around the linear interpolations have positive
bridge probability, and the interpolations are separated throughout.
The positive lower bound persists in some neighborhood of the
endpoints, so continuity of (EH4) gives strict positivity there
even if the density formula is first obtained almost everywhere.
This proves strict total positivity of \(w_T\) on decreasing rows
and decreasing columns, by applying the same argument to every
selected submatrix.

The Vandermonde is harmonic because its Laplacian is alternating
of degree two less. On the strictly ordered chamber, the drift
contribution to its time-dependent generator is
\[
 \Delta(x)\sum_{a<b}\frac{b_t(x_a)-b_t(x_b)}{x_a-x_b}\le0.
\]
Up to the first collision, Itô's formula and localization therefore
give a nonnegative supermartingale after stopping at the collision,
where \(\Delta=0\). Dropping the path weights, which are at most
one, gives
\[
 \int_{y_1>\cdots>y_j}\det[w_T(x_a,y_b)]\,\Delta(y)\,dy
       \le\Delta(x).
\]
For the limiting kernel \(w\), pointwise convergence of the entries,
nonnegativity of all the ordered determinants and Fatou give the
same bound
\[
 \int_{u>y_1>\cdots>y_j}\det[w(x_a,y_b)]\,\Delta(y)\,dy
       \le\Delta(x).
 \tag{EH5}
\]
For \(j=1\) the Vandermonde is the empty product.

Apply the leading-minor bound (S63) to the strictly totally positive
matrix \([w_T(t_i,y_j)]_{i=1,\ldots,p;\,j=0,\ldots,p-1}\), where
\(u>y_0>y_1>\cdots>y_{p-1}\). First take \(T\uparrow1\), then divide
by \(u-y_0\) and take \(y_0\uparrow u\). Equations (EH2) and
(EH4) yield
\[
 0\le \det[f'(t_i),\,w(t_i,y_1),\ldots,w(t_i,y_{p-1})]_{i=1}^p
 \le f'(t_1)\det[w(t_i,y_j)]_{i=2,\ldots,p;\,j=1,\ldots,p-1}.
 \tag{EH6}
\]
By multilinearity and antisymmetrization in the integration variables,
the determinant on the left of (EH3), including its sign \(s_p\),
is the integral of the determinant in (EH6) times
\((p-1)!\Delta(y)\) on \(u>y_1>\cdots>y_{p-1}\). To see the factor,
apply the moment formula (EH2) to each column. The determinant of
the moment monomials \([r y_j^{r-1}]_{j,r=1}^{p-1}\) is
\(s_p(p-1)!\Delta(y)\). Splitting the integral over all \(y\)
into order chambers gives the ordered integral with no other
factor. The Gaussian tails justify the integrations.
Now (EH5) and (EH6) give (EH3).

## Applying the derivative bound on the facet

Use the notation of (S92a), including \(m_1=m_2=0>m_3>\cdots>m_k\),
and set \(S=2f(0)+\sum_{i=3}^k f(m_i)\). For the increasing
Helmert index list \(2\le j_1<\cdots<j_p\le k-1\) define
\[
 L_j[v]=\sum_{i=1}^j v(m_i)-jv(m_{j+1}),\qquad
 L_{j+1}[v]-L_j[v]
    =(j+1)\int_{m_{j+2}}^{m_{j+1}}v'(t)\,dt.
 \tag{EH7}
\]
The positive Helmert row factors \(\alpha_j\) can be removed from
both sides of the coordinate comparison. After multiplying (S92a)
by \(f(0)\), the interaction coordinate is \(s_p\) times the
determinant with rows \(L_{j_a}[v]\), for
\(v=(f,M_1,\ldots,M_{p-1})\). The constants subtracted in the \(U_r\)
terms do not affect a Helmert row. The comparison coordinate is
\(S W_p\), whose determinant rows, after the same normalization,
are \((1,L_{j_a}[m],\ldots,L_{j_a}[m^{p-1}])\) with sign
\(s_p\).

Subtract consecutive rows from bottom to top in the interaction
determinant. Use \(L_1[v]=0\) and expand each resulting row by
(EH7). The row variables now lie in successively lower disjoint
intervals. Inside each integral apply (EH3). The integral over the
first row variable and its integer weights is \(L_{j_1}[f]<S\):
the defining partial sum omits strictly positive terms and the
subtracted term is positive.
For the other \(p-1\) rows the remaining integrand is
\((p-1)!\Delta(t_2,\ldots,t_p)\). This is exactly the determinant
formula for the normalized coordinate of \(W_p\): subtract its
consecutive rows, expand on the first column, and use (EH7) with
the powers. Thus the conditional interaction coordinate is at
most \(L_{j_1}[f]W_p\) in the chosen coordinate, strictly less than
\(S W_p\) there.

Integrating the conditional comparison against the common density
in (S92a) proves (S91f) for every exterior degree. Use the exact
scaling and normalized positive operator argument after (S91f),
and the ambient Hermite decomposition underlying (S42). In each active
block the normalized compact, self-adjoint, positive-kernel operator
\(T\) satisfies \(T\psi<\rho^{p(p-1)/2-1}\psi\) for the interior
positive field \(\psi=\sqrt{a_\rho}W_p\). Pairing
\(|Th|\le T|h|\) with \(\psi\) gives the corresponding strict bound
on every nonzero eigenvalue in magnitude. Compactness gives the
norm bound, and \(a_\rho/\rho\ge\ell/\sqrt{1-\rho^2}\) gives the
gap. Extra Hermite degree \(m\) multiplies the normalized interaction
by \(\rho^m\). This gives (EH1) and also the strict gap when \(p=2\).

The diffusion and determinant comparisons are written proofs. They
have received internal mathematical review.

The [reflection-chamber comparison](REFLECTION_CHAMBER_DIFFUSION_COMPARISON.md)
extends the derivative-kernel estimate to general Gaussian endpoint tilts
and to every character of a finite real reflection group.
