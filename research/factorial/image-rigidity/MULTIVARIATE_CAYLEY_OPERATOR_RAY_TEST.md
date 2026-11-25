# A Cayley operator test for angular separation

A single operator carries the grouped logarithmic coefficient corrections
through the Cayley transform in any degree. This gives a sufficient ray
test for a pair of orders once their references have disjoint angular
sector counts. The test applies to ordinary profiles with shared or
private rows.

## Exact change of operator

Take \(f_j(t)=\sum_{k=0}^j a_{j,k}t^k\) and a reference
\(g_j(t)=\sum_{k=0}^j b_{j,k}t^k\), both with positive real coefficients.
For \(\lambda>0\) write
\[
 \Phi(z)=\frac{\lambda z-1}{1+\lambda z},\qquad
 \mathcal C_j f(z)=\lambda^{-j}(1+\lambda z)^j f(\Phi(z)).
\]
Let \(D_j(v)\) be the real polynomial of degree at most \(j\) with
\(D_j(k-j/2)=\log(a_{j,k}/b_{j,k})\) for \(0\le k\le j\).
On polynomials of degree at most \(j\) define
\[
 \mathcal V_j=\frac{\lambda^2z^2-1}{2\lambda}\frac{d}{dz}
                   -\frac{j\lambda z}{2}.
\]
Then, with \(G_j=\mathcal C_j g_j\), the exact identity is
\[
                 \mathcal C_j f_j
                    =\exp\!\big(D_j(\mathcal V_j)\big)G_j.                 \tag{1}
\]
Indeed, \(\Phi'(z)=2\lambda/(1+\lambda z)^2\), and the derivative
of the prefactor in \(\mathcal C_j\) gives
\[
 \mathcal C_j(t\partial_t-j/2)\mathcal C_j^{-1}
   =\frac{\Phi}{\Phi'}\partial_z
      -\frac{j\lambda\Phi}{(1+\lambda z)\Phi'}-\frac j2
   =\mathcal V_j.
\]
Thus \(\mathcal C_j(t^k)\) is an eigenvector with eigenvalue \(k-j/2\),
and applying the exponential to that basis proves (1). The values
on the centered integer grid determine the operator, so the formula
also applies to a convergent series of polynomial corrections evaluated
on the same spectrum.

An affine term in \(D_j\) can be removed before an angular comparison.
For a real term \(c+d v\), its exponential multiplies \(f_j\) by a
positive constant and changes its variable by a positive dilation.

## A sufficient ray inequality

For \(L>0\) use the weighted coefficient norm
\(\|\sum_m c_m z^m\|_L=\sum_m|c_m|L^m\) and its induced operator
norm. The action
\[
 \mathcal V_j z^m=-\frac{m}{2\lambda}z^{m-1}
                  +\frac{(m-j)\lambda}{2}z^{m+1}
\]
gives the exact column-sum norm
\[
                  \|\mathcal V_j\|_L
             =\frac j2\max\{(\lambda L)^{-1},\,\lambda L\}.
\]
Put \(A_j=D_j(\mathcal V_j)\) and \(d_j=\|A_j\|_L\).
The exponential series implies, for \(0\le s\le1\),
\[
 \big\|\exp(sA_j)G_j-G_j\big\|_L
       \le \frac{e^{d_j}-1}{d_j}\,\|A_jG_j\|_L
       \le (e^{d_j}-1)\|G_j\|_L,                                         \tag{2}
\]
where the first factor has its continuous value \(1\) at \(d_j=0\).
The first bound permits cancellation in the action on the reference.
Either bound also controls the pointwise error on \(|z|\le L\).

For orders \(n,N=n+1\), choose finitely many directions \(\Theta\)
whose open angular sectors each contain roots of at most one of
\(g_n,g_N\). Require that neither reference has a root on the rays.
At each order choose a scale \(\lambda=\lambda_j\) and radius \(L_j\)
such that all \(j\) zeros of \(G_j\) are in \(|z|<L_j\). One possible
certificate for this last condition is
\(\lvert G_{j,j}\rvert L_j^j>\sum_{m<j}\lvert G_{j,m}\rvert L_j^m\)
for the coefficients \(G_{j,m}\). Let \(K_j\) consist of the disk
boundary and the closure in the closed disk of the portions of the
inverse ray images \(\Phi^{-1}\{r e^{i\theta}:r>0\}\),
\(\theta\in\Theta\). The following inequality at both orders suffices:
\[
           \frac{e^{d_j}-1}{d_j}\|A_j G_j\|_{L_j}
                     < \min_{z\in K_j}|G_j(z)|.                           \tag{3}
\]
All \(j\) roots stay in the disk along the exponential homotopy by
Rouché's theorem on its boundary. The homotopy is the Cayley image
of the polynomial with coefficients \(b_{j,k}^{1-s}a_{j,k}^s>0\),
so its roots avoid \(\Phi^{-1}(0)\) and the pole of \(\Phi\).
Inequality (3) keeps the roots off every inverse separator ray.
Continuity of the multiset of roots then preserves the count in
each original angular sector. Consequently the \(f_j\) have
disjoint root arguments.

## Ordinary factorial moments

For two ordinary exponent columns \((a_i),(b_i)\), let
\[
 p_{j,k}=\binom jk\prod_i((j-k)a_i+kb_i)!,
 \qquad a_{j,k}=\frac{p_{j,k}}{p_{j,0}^{1-k/j}p_{j,j}^{k/j}}.
\]
This endpoint normalization preserves root arguments. On an unequal
row use \(\ell=\min(a_i,b_i)\), \(d=|a_i-b_i|\), and
\(x=d/(j\ell+h)\) for \(1\le h\le d\). Each such factor on an
increasing row contributes
\[
 H(k)=\sum_{l<k}\log(1+lx)-\frac{k}{j}\sum_{l<j}\log(1+lx)
\]
to \(\log(a_{j,k}/\binom jk)\), with indices in each sum starting
at zero. Grouping the \(dk\) factorial factors by their residues
modulo \(d\) gives this formula also at \(\ell=0\). A decreasing row
uses \(H(j-k)\).

For \(j\ge2\) put \(T=(j-1)x/(2+(j-1)x)\), and let \(U_j\)
be the sum of \(\operatorname{artanh}T\) over the factors. Use the
Gaussian reference \(b_{j,k}=\binom jk\exp[-U_j k(j-k)/(j-1)]\).
For each integer \(r\ge1\), the polynomial
\[
 S_{j,r}(k)=\sum_{l<k}\left(1-\frac{2l}{j-1}\right)^r
                -\frac{k}{j}\sum_{l<j}\left(1-\frac{2l}{j-1}\right)^r
\]
has degree at most \(r+1\), by polynomial summation. Reflection
gives \(S_{j,r}(j-k)=(-1)^{r+1}S_{j,r}(k)\) on the grid.
Let \(V_r\) be the sum of the \(T^r\)'s and \(M_r\) their sum
with sign positive for increasing rows and negative for decreasing
rows. Expanding \(\log(1-(1-2l/(j-1))T)\) and using
\(S_{j,1}(k)=k(j-k)/(j-1)\) gives the exact correction
\[
 D_j(k-j/2)=
 \sum_{\substack{r\ge3\\r\ {\rm odd}}}
       \frac{V_r}{r}\big(S_{j,1}(k)-S_{j,r}(k)\big)
 -\sum_{\substack{r\ge2\\r\ {\rm even}}}
       \frac{M_r}{r}S_{j,r}(k).                                            \tag{4}
\]
Every \(T\) is less than one, and the sums converge absolutely
on the finite grid. Conjugation to the eigenbasis of \(\mathcal V_j\)
then gives norm convergence of the corresponding matrix series.
Equation (4) therefore supplies the grouped directions in (1)
directly from the ordinary exponents. Separation for the
normalized pair implies coprimality of the original rational
polynomials over every characteristic-zero field by Bézout.

The [matched transport inequality](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
bounds how the scalar moments in (4) change across orders.
A uniform application of (3) still requires an upper bound on the
transformed error relative to a lower bound on the reference ray
margin, throughout the parameter range. The
[scaled Hermite comparison](MULTIVARIATE_SCALED_HERMITE_INTERLACING.md)
provides strict interlacing for the limiting references across the
physical deficit ratios. Extending the ray bounds beyond that limit
requires estimates at the angular gap scale. The
[near-flat coefficient relaxation](MULTIVARIATE_MATCHED_TANH_MOMENTS.md#a-failure-of-independent-errors-in-the-near-flat-angular-comparison)
exhibits why replacing grouped directions with independent
coefficient tails of fixed weight can lose that scale.
