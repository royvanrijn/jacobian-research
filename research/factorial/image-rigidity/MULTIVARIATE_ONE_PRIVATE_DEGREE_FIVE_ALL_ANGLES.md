# All-degree angular separation at private degree five

Let \(a,b\) be ordinary distinct nonnegative integer exponent columns
with private rows on exactly one side, of total private exponent
\(P=5\). The factorial moment polynomials
\[
 Q_m(z)=\sum_{k=0}^m\binom mk\prod_i((m-k)a_i+kb_i)!\,z^k
\]
have disjoint root arguments at consecutive orders. In particular
they are coprime over every characteristic-zero field.

The new range is \(N=n+1\ge15\) and
\(S_{\rm sh}=\sum_{\rm sh}(b_i-a_i)^2/\max(a_i,b_i)<41/18\).
The [stronger entropy bound](MULTIVARIATE_ONE_PRIVATE_ENTROPY_SHARPER_RESERVE.md)
handles \(S_{\rm sh}\ge41/18\). The terminal
[fourteen theorem](MULTIVARIATE_TERMINAL_FOURTEEN_ALL_PROFILES.md)
and the [earlier angular proofs](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#matching-the-long-phase-hypotheses)
give the conclusion for \(N\le14\).

The argument uses the actual first coefficient \(B_m\) to set the
scale of the half at the private endpoint. Cancellation of one
private Gamma factor gives a coefficient envelope on the full half
and bounds the interacting part of its logarithm. The adjacent
logarithms control the phase after the radii have been matched.
Edo and van den Essen's Theorem 3.2 gives the coprimality reduction
in two variables, and their Proposition 3.5 covers two specified
support templates. Filaseta and Rocks give a density bound for
irreducibility on the single-power template \(1+\lambda X^m\).
The shared rows in the present theorem can have arbitrary exponents.
[Edo–van den Essen, §3](https://arxiv.org/pdf/1304.3956)
[Filaseta–Rocks, Theorem 1.1](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf)

## Actual coefficient and first-factor bounds

Remove equal rows and orient from the side with private exponents
\((e_i,0)\). Set \(d_i=b_i-a_i\) on shared rows and put
\[
 W=4+\sum_{\rm sh}d_i^2/a_i,\qquad
 L=\max(1,\max_{\rm sh}|d_i|/a_i),\qquad S=S_{\rm sh}.
\]
The maximum over an empty set is omitted. For a shared row with
smaller exponent \(\ell\) and gap \(d=|d_i|\), \(d\le3\ell\).
Otherwise \(d\ge4\) and \(d^2/(\ell+d)>3d/4>41/18\).
Consequently \(\sum_{\rm sh}d_i^2/\ell_i\le4S\), and
\(4\le W<14\). If \(L>1\) the maximizing row increases and
contributes at least \(L^2\) to \(W-4\), so \(L/W\le1/4\).
This inequality also holds for \(L=1\).

Normalize the two endpoint coefficients to one by a positive
rescaling of the variable and multiplication by a positive scalar.
For real \(m\in[n,N]\) and \(0\le x\le m\), the normalized
interpolation is
\[
 c_{m,x}=\frac{\Gamma(1+m)}{\Gamma(1+x)\Gamma(1+m-x)}
 \prod_i\frac{\Gamma(1+ma_i+d_i x)}
 {\Gamma(1+ma_i)^{1-x/m}\Gamma(1+mb_i)^{x/m}},
\]
where \(d_i=b_i-a_i\) on all rows. At integer \(m\) split the
normalized polynomial as
\[
 H_m^{\rm e}(z)+z^m H_m^{\rm s}(1/z).
\]
The first half starts from the private exponents, the second from
the opposite endpoint. Each half takes half of the middle
coefficient if \(m\) is even. Bounds on the inclusive coefficients
apply to either half. Write \(B_m=c_{m,1}\) and
\(f_m(x)=\log(c_{m,x}\Gamma(1+x)/B_m^x)\).
The [cancelled Gamma expression](MULTIVARIATE_ONE_PRIVATE_CANCELLED_GAMMA_BUDGET.md#cancellation-and-the-derivative-bounds)
shows \(f_m''\ge0\). It gives, for integer \(2\le k<n\)
in the first half, with \(\eta_k=n-k-1/2\),
\[
 0\le f_m(k)\le F_k:=\frac{Wk(k-1)}{2\eta_k},
 \qquad -F_k/\eta_k\le\partial_m f_m(k)\le0.                 \tag{1}
\]
The first upper bound follows by twice summing forward
differences bounded by \(W/\eta_k\), using \(f_m(0)=f_m(1)=0\).

Write \(C_m=(1+\frac12\log m)/m\) and
\(\bar q(m)=\exp(-5+5 C_m)\).
In the oriented envelope retain all shared rows and take
its private chord parameter zero. Write \(Q_{\rm sh}(w)\)
for that envelope and \(q_{\rm sh}=Q_{\rm sh}(0)\).
The interaction budget with \(h=1\) gives
\(W q_{\rm sh}\le(4+S)e^{-S/2}\). The shared-row increment
bound gives \(Q_{\rm sh}(1/m)\le q_{\rm sh}e^{3S/m}\) since
\(\sum_{\rm sh}d_i^2/\ell_i\le4S\).
[Interaction budget (4)](MULTIVARIATE_ONE_PRIVATE_CANCELLED_GAMMA_BUDGET.md#tight-endpoint-budgets-for-the-interaction-terms)
[Oriented Gamma chord and increment](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#a-gamma-chord-bound-in-both-orientations)

In either orientation the private product is bounded by
\(f_{m,j}^5\), where \(f_{m,j}=\Gamma(1+j)/\Gamma(1+m)^{j/m}\)
and \(j=m-k\) or \(k\), respectively. This is the convexity
bound for \(\log\Gamma(1+ex)-e\log\Gamma(1+x)\).
At \(k=1\) in the first orientation the Gamma shift gives
\(f_{m,m-1}=\Gamma(1+m)^{1/m}/m\le e^{-1+C_m}\).
Thus \(B_m/m\le\bar q(m)Q_{\rm sh}(1/m)\). For \(m\ge14\),
\((4+S)\exp[-(1/2-3/m)S]\) decreases with \(S\ge0\).
Together with the interaction estimate this proves
\[
 B_m/m\le\bar q(m),\qquad WB_m/m\le4\bar q(m),\qquad
 LB_m/m\le\bar q(m).                                        \tag{2}
\]
These bounds hold for the actual row product at real orders.
[Private factorial convexity and Stirling bounds](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#the-private-factorial-estimates)
[Signed Stirling remainder](https://dlmf.nist.gov/5.11.ii)

For the first-factor comparison, the binomial contributes \(1/m\)
to \((\log B_m)'\). The loss from the shared rows is nonnegative
and at most \(3\sum_{\rm sh}d_i^2/\ell_i/(4m^2)\).
The loss from each private row is nonnegative by the same
positive-start derivative formula. Directly differentiating
that row's exact factor
\(\Gamma(1+em)^{1/m}/\prod_{j=0}^{e-1}(em-j)\)
also bounds its loss by
\[
 \frac{1+\frac12\log(em)+(e-1)/(2(1-1/m))}{m^2}.
\]
Indeed \(\psi(1+y)\ge\log y\) for \(y>0\), by comparison
of \(\sum_{j=1}^M(y+j)^{-1}\) with
\(\int_0^M(y+t)^{-1}\,dt\) in the digamma series limit.
Apply this and the upper Stirling bound to the numerator
\(em\psi(1+em)-\log\Gamma(1+em)-em\). The other derivative
terms have total loss \(\sum_{j=0}^{e-1}j/[m(em-j)]\).
Use \(j/(em-j)\le(j/e)/(m-1)\) and
\(\sum_{j=0}^{e-1}j/e=(e-1)/2\) to obtain the displayed bound.
[Positive-start derivative formula](MULTIVARIATE_ABSOLUTE_CURVATURE_15.md#exact-derivative-at-a-positive-start)
[Digamma series](https://dlmf.nist.gov/5.7.E6)

Summing the losses and using \(\log e_i\le e_i-1\) gives
\[
 0\le \frac1m-(\log B_m)'\le
 \frac{5(1+\frac12\log m)+2/(1-1/m)+41/6}{m^2}
 <\frac2m .                                                   \tag{3}
\]
Each term in the numerator divided by \(m\) decreases for
\(m\ge14\), and the scalar test at 14 gives the last inequality.
In particular \(|(\log B_m)'|\le1/m\).
Let \(D=\max_{n\le m\le N}B_m\) and \(U=(n+1)\bar q(n)\).
Then \(D\le U,\ WD\le4U,\ LD\le U\) by (2), since \(C_m\)
decreases on \(m\ge14\).

### Tail bounds in the two orientations

For the first half define the exact real-index envelope
\(U_m(x)=(c_{m,x}\Gamma(1+x))^{1/x}\), \(0<x\le m/2\).
For integral \(k\) one has \(c_{m,k}^{1/k}\le e U_m(k)/k\).
The maximum of \(e U_m(x)/x\) on any closed subinterval of
\((0,m/2]\) occurs at an endpoint. To prove this, it suffices
to show that \(x(\log U_m)'\) is increasing. After cancellation,
the nonaffine terms in \(\log(c_{m,x}\Gamma(1+x))\) are
\(g(x)=\log\Gamma(s+dx)\) with either \(d<0\) or \(d>0,\ s\ge1\).
The derivative of \(x[(g(x)-g(0))/x]'\) is
\(\int_0^1 t[g''(xt)+xt\,g'''(xt)]\,dt\).
For decreasing arguments both terms in the integrand are
nonnegative. For an increasing argument set \(y=s+dx\).
The sign follows from \(\psi_1(y)+(y-s)\psi_2(y)\ge0\).
In fact \(v\psi_1(1+v)=\int_0^\infty e^{-u}h(u/v)\,du\)
increases for \(v>0\), where \(h(t)=t/(e^t-1)\) decreases.
Apply this with \(v=y-1\ge y-s\) and \(\psi_2(y)<0\).
Thus \(x(\log(eU_m(x)/x))'\) is increasing, giving the
endpoint maximum. The integral for \(\psi_1\) follows by
differentiating the digamma integral at positive argument.
[DLMF 5.9.16](https://dlmf.nist.gov/5.9.E16)

At \(x=m/2=t\), the private product is bounded by \(f_{m,t}^5\)
and the shared product by one. Substituting the Gamma binomial
and using Stirling bounds gives
\[
 \frac{e U_m(t)}{t}
 \le \frac e t
       \left(\frac{\Gamma(1+t)^4}{\Gamma(1+2t)^{3/2}}\right)^{1/t}
 \le \frac{e^{4 C_t}}8 .
\]
For any integer \(2\le K<n/2\), the other endpoint bound in
the same envelope follows from (1). For the second half reserve
one \(f_{m,k}\) to cancel \(\Gamma(1+k)\) in the binomial
as in the [entropy proof](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md#matching-the-long-phase-hypotheses).
The remaining four factors give a separate bound.
\[
\begin{split}
 c_{m,k}^{1/k}
 &\le \max\!\left(
       \frac{eD}{K}e^{W(K-1)/(2\eta_K)},\frac{e^{4 C_{m/2}}}{8}
            \right)
     &&\text{at the private start, } K\le k\le m/2,\\
 c_{m,m-k}^{1/k}
 &\le e f_{m,k}^{4/k}
 \le \min\!\left(e^5(k!)^{4/k}/m^4,\frac e{16}e^{4 C_{m/2}}\right)
     &&\text{at the opposite start, } 1\le k\le m/2 .
\end{split}                                                    \tag{4}
\]
The second midpoint estimate uses the increase of
\(\Gamma(1+k)^{1/k}\) with \(k\), by convexity of \(\log\Gamma\).
For small \(m\) the replay also uses the inclusive bound
\(\binom mk f_{m,j}^5\) directly in each orientation.

## From half-log bounds to phase separation

For holomorphic power series let \(\|\cdot\|_R\) be the absolute
coefficient norm at \(R\). Let \(\mathcal E=z\partial_z\).
Suppose both halves have logarithms on \(|z|\le R\), real
on the positive ray. Let \(M,F,C\) bound the sum of the two
individual log norms and their first and second Euler norms,
respectively, in each order. Let \(A,E\) bound the sum of the
two adjacent log difference norms and their Euler norms,
evaluated at a common argument of the series. If \(M<n\log R\),
Rouché on the two circles in the original and reciprocal variables
confines every root to \(|\log|z||<\log R\). At \(|z|=R\),
the modulus ratio \(|H_m^{\rm e}(z)/H_m^{\rm s}(1/z)|\)
is bounded by \(e^M<R^m\), and \(z^m H_m^{\rm s}(1/z)\)
has all \(m\) zeros inside this circle by nonvanishing of the
half on its disk. The reciprocal argument is identical.

Write \(u_m=\log|z_m|\) for roots at a common angle \(\theta\)
and \(J_m(u,\theta)=\log H_m^{\rm e}(e^{u+i\theta})
-\log H_m^{\rm s}(e^{-u-i\theta})\).
Their moduli satisfy \(m u_m=\Re J_m(u_m,\theta)\). Hence, if
\(F<n\), subtraction gives
\[
 n|u_N-u_n|\le v=\frac{\log R+A}{1-F/n}.
\]
At fixed radius integrate the adjacent Euler difference from
the positive ray. For the change of radius, integrate the
mixed derivative in angle and log radius, using that the
imaginary part vanishes at \(\theta=0\) for each radius. This
yields
\[
 |\Im(J_N(u_N,\theta)-J_n(u_n,\theta))|
                 \le (E+Cv/n)\theta.                           \tag{5}
\]
If \(E+Cv/n<1/2\), the phase equations would give
\(\theta-(\Im J_N-\Im J_n)\in2\pi\mathbb Z\) strictly between
\(\theta/2\) and \(3\theta/2\) for \(0<\theta<\pi\).
Both ranges below meet this same phase bound.

## The range \(D\le3/5\)

Use \(R=1+1/n\). Set \(K=4\) for \(14\le n<100\) and \(K=8\)
for \(n\ge100\). Extract \(e^{B_mz}\) from the first half by
putting \(J_m^0=H_m^{\rm e}-e^{B_mz}\) and
\(T_m=e^{-B_mz}J_m^0\).
The coefficients of \(J_m^0\) vanish at indices zero and one.
At \(2\le k\le K\), (1) and (3) bound the coefficient and its
real-order derivative in absolute value by
\[
 \frac{B_m^k}{k!}(e^{F_k}-1),\qquad
 \frac{B_m^k}{k!}
       \left(\frac{k}{n}(e^{F_k}-1)+e^{F_k}F_k/\eta_k\right).
                                                               \tag{6}
\]
For \(n<100\), partition \(4\le W<14\) at 6 and 9. On an
interval \([w_-,w_+]\) apply (6) with
\(B_m\le\min(3/5,U,4U/w_-)\) and \(W\le w_+\).
Take the maximum for each coefficient and derivative.
For the coefficient tails use the sum of the two integer
order bounds \(\binom mk f_{m,j}^5\). To calculate them, the
replay bounds the factorial mean \((m!)^{1/m}\) below
by a rational with denominator \(10000\), by integer
comparison of its \(m\)-th power with \(m!\).

For \(n\ge100\) bound \(B_m\le3/5\), \(W\le14\), use (4)
with \(K=8\), and use \(R c_{m,k}^{1/k}<9/25\) from \(k=9\).
The opposite half uses the second bound in (4) for \(k\le4\)
with \((k!)^{1/k}\le9/4\), and \(R c_{m,m-k}^{1/k}<11/50\)
from \(k=5\). The scalar estimates for this range can all use
the upper constants at \(n=100\). The coefficient and derivative
bounds decrease as \(n\) increases, and the annulus test uses
the increasing lower bound \(n\log(1+1/n)\ge1-1/(2n)\).

Let \(j_0,j_1,j_2\) bound the individual Euler norms of \(J_m^0\)
and \(a_0,a_1\) the adjacent difference norms. Sum (6), with
weights \(k^jR^k\), and add bounds for the polynomial and
exponential tails in both calculations. For each polynomial
or exponential tail, the difference of the two order coefficients
is bounded by their common upper bound because they are
nonnegative. Put \(x=R\min(3/5,U)\) (use \(x=3R/5\) for
the uniform range) and \(b=x/n\). Then \(R|B_N-B_n|\le b\).
Product differentiation gives the following bounds for the
individual Euler norms of \(T_m\) and its adjacent differences:
\[
\begin{array}{ll}
 t_0=e^x j_0,& d_0=e^x(a_0+bj_0),\\
 t_1=e^x(j_1+x j_0),&
 d_1=e^x\{a_1+x a_0+b[j_1+(1+x)j_0]\},\\
 t_2=e^x\{j_2+2xj_1+(x+x^2)j_0\}.&
\end{array}
\]
For the differences use \(\|e^{-B_Nz}-e^{-B_nz}\|_R\le e^x b\)
and its Euler bound \(e^x b(1+x)\), obtained on the segment
between \(B_N\) and \(B_n\).
Let \(l_0,l_1,l_2\) be the coefficient and Euler bounds for
\(H_m^{\rm s}-1\), also bounding the adjacent difference for
orders zero and one by the nonnegativity of coefficients.
The replay gives \(t_0,l_0<1\). The log series for the halves
and integration along segments between their residual series
give, with \(h=1-t_0,\ s=1-l_0\),
\[
\begin{split}
 M&=x+t_0/h+l_0/s,& F&=x+t_1/h+l_1/s,\\
 C&=x+t_2/h+(t_1/h)^2+l_2/s+(l_1/s)^2,\\
 A&=b+d_0/h+l_0/s,&
 E&=b+d_1/h+t_1d_0/h^2+l_1/s+l_1l_0/s^2 .
\end{split}
\]
For the adjacent Euler norms use the quotient rule on the
two logarithms. The replay verifies \(M<1-1/(2n)\),
\(F<n\), and (5) with \(E+Cv/n<1/2\), using
\(\log R\le1/n\) in the bound on \(v\).

## The range \(D>3/5\)

Here \(n\ge70\) by (2) and the scalar tests at \(14\le n<70\).
Put \(q_*=43/5000\), \(K=\lceil10D\rceil+8\),
\(\eta=n-K-1/2\), and \(\eta_*=777/1000\).
The bounds (2) give \(D/n,L D/n<q_*\), and \(WD/n<4q_*\).
Also \(K\le25D\), \(\eta>\eta_* n\), and \(K+1<n/2\).
The first inequality of (4) gives \(c_{m,k}^{1/k}<19/40\)
throughout the tail from \(K\).

For \(2\le r\le K\) define \(\nu_r=\Delta_x^r f_m(0)\)
and take \(h_2=W/(n-5/2)\), \(\alpha=L e^{1/\eta}/\eta\).
Apply the cancelled Gamma bounds with the index cutoff \(r\).
Since
\(\log((n-5/2)/(n-r-1/2))\le(r-2)/\eta\), the bounds become
\[
 |\nu_r|\le (r-2)!\alpha^{r-2}h_2=:H_r,\qquad
 |\partial_m\nu_r|\le (r-1)H_r/\eta.
\]
With \(u=\alpha K<283/1000\) one has \(H_r\le h_2\).
On the first \(K\) coefficients the exponential formula for
connected hypergraphs gives \(\ell_k(m)=B_m^k\Sigma_k/k!\)
for the formal log coefficient, where \(\Sigma_k\) is the
sum of connected products with edge weight \(e^{\nu_r}-1\)
for an edge of size \(r\). This follows by expanding
\(e^{f_m(k)}=\prod_{r=2}^k(e^{\nu_r})^{\binom kr}\)
over labelled edges and decomposing each hypergraph into
connected components.

Bound the weights by \(e^{h_2}H_r\). In the derivative of
the product the sizes contribute \(X/\eta\), with
\(X=\sum_{\rm edges}(r-1)\). Since
\(X\le(8/e)k\,e^{X/(8k)}\), the derivative can use the same
inflated edge bounds \((r-2)!\alpha^{r-2}\beta\) as the
undifferentiated sum, where \(\beta=e^{h_2+1/8}h_2\),
with prefactor \(8k/(e\eta)\).
By the discovery-tree sum, if \(\alpha/\beta\le1/3\),
the connected bound is \((k-1)!\Xi^{k-1}\), with
\[
 \Xi=\frac{(256/81)\beta}{1-u}
          \exp\!\left(\beta K\left(\frac12+\frac{u}{6(1-u)}\right)\right).
\]
Here \(256/81\) bounds the discovery-tree base
\((1+t)^{(1+t)/t}\) at \(t=\alpha/\beta\le1/3\), and the
exponent bounds additional edges by
\(\sum_{r=2}^k\binom kr(r-2)!\alpha^{r-2}\beta
\le\beta k(k-1)(1/2+u/[6(1-u)])\).
[Discovery-tree sum](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum)
[Differentiated bound](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md#a-smaller-differentiated-connected-logarithm-bound)

Indeed \(h_2D\le4q_*(28/27)\),
\(\alpha/\beta\le e^{1/\eta}/(4\eta_*)<1/3\), and
\(\beta D<43/1000\). The replay gives \(\Xi B_m<p=7/20\)
and \(1+8/(e\eta_*)<5\). Combining with (3), the log
coefficients through \(K\) and their derivatives satisfy
\[
 |\ell_k(m)|\le B_m p^{k-1}/k,\qquad
 |\ell_k'(m)|\le 5B_m p^{k-1}/n .
\]
At \(k=2\) use the single edge of size two to improve the
first bound to \(|\ell_2(m)|\le(19/500)B_m/2\), since
\(e^{h_2}h_2D<19/500\).

Write \(A_m(z)=\sum_{k=1}^K\ell_k(m)z^k\) and take
\(r_i=6/5,\ r_o=21/10\).
Summing with \(1/k\le1/3\) for \(k\ge3\) gives
\(\|A_m\|_{r_i}<(7/5) B_m\), \(\|A_m\|_{r_o}<(37/10) B_m\).
Also \(\|A_N-A_n\|_{r_i}\le5(D/n)r_i/(1-pr_i)\).
The series \(e^{A_m}\) and \(H_m^{\rm e}\) agree through \(K\).
On \(r_i\) the norm of
\((H_m^{\rm e}-e^{A_m})e^{-A_m}\) is bounded by
\[
 (4/7)^{K+1}e^{(51/10)D}
       +\frac{(57/100)^{K+1}}{1-57/100}e^{(7/5)D}.
\]
The first term uses the outer radius for the exponential
tail, the second the actual coefficient tail from (4).
The first is below \(e^{-519/100}\) by \(K+1\ge10D+9\)
and \(\log(7/4)>11/20\). For the second use
\(D\le(K-8)/10\), \(K\ge14\), and \((57/100)e^{7/50}<1\).
The log correction \(\|\log H_m^{\rm e}-A_m\|_{r_i}\)
is below \(1/125\) by the log series. The opposite half
uses (4) with \(n\ge70\), head indices up to four and
ratio \(29/100\) on \(r_i\) from index five.
It gives \(\|\log H_m^{\rm s}\|_{r_i}<1/250\).

On \(r_i\) the sum of individual log norms divided by \(n\)
is below \(1/80\), and the sum of adjacent log norms is
below \(23/200\). At \(R=21/20\) the first and second Euler
multipliers from \(r_i\) are bounded by \(3\) and \(33\),
using \(\log(r_i/R)>13/100\) and maximizing
\(x e^{-x\log(r_i/R)}\) and \(x^2 e^{-x\log(r_i/R)}\).
Thus in (5) one can use \(M/n<1/80<\log R\),
\(F/n<3/80\), \(C/n<33/80\), \(A<23/200\) and
\(E<3(23/200)\). With \(\log R<1/20\) this gives
\[
 E+Cv/n<
 3\frac{23}{200}
 +\frac{33}{80}\frac{1/20+23/200}{1-3/80}
 <\frac12.
\]

Conjugation handles the lower half-plane. The positive real
ray has no roots, and the even-order factorial moment is
strictly positive on the real line by its product-exponential
integral, since the two monomials are distinct. The conclusion
is invariant under the independent positive normalizations.
The integer polynomials have a rational Bézout identity,
which transfers coprimality to every characteristic-zero field.

## Exact scalar replay

Run the following from the repository root.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_one_private_degree_five.py
~~~

The [checker](verify_multivariate_one_private_degree_five.py)
uses rational bounds on elementary series and integer
comparisons for factorial means. Its only order loop is
\(14\le n<100\) in the direct range. The three intervals for
\(W\) enclose the whole interaction budget and the remaining
estimates are uniform over the parameter ranges. The checker
passed with CPython 3.12.3 in the internal replay.
