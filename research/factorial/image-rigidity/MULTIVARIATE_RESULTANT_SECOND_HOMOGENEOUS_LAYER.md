# Positive second homogeneous layer for unbalanced factorial resultants

Fix a nonempty list of unequal exponent rows with gaps \(d_i\ge1\)
and orientations \(\epsilon_i\in\{1,-1\}\). Let \(A_i\) denote the
independent smaller exponents, and define the normalized ordinary
moment polynomials by
\[
 F_j(A,z)=\sum_{k=0}^j\binom jk z^k
      \prod_i\prod_{h=1}^{d_i\ell_{j,i}(k)}(jA_i+h),
 \qquad
 \ell_{j,i}(k)=
 \begin{cases}
 k&\epsilon_i=1,\\
 j-k&\epsilon_i=-1.
 \end{cases}
\]
Take \(N=n+1\ge2\), \(R(A)=\operatorname{Res}_z(F_n,F_N)\),
\(D_\pm=\sum_{\epsilon_i=\pm1}d_i\), \(D=D_++D_-\), and
\(\Delta=D_+-D_-\). Assume \(\Delta\ne0\).

The [leading homogeneous part](MULTIVARIATE_ALL_ORDER_EXCEPTIONAL_LOCUS.md#the-highest-homogeneous-part-of-one-resultant)
of \(R\) is the monomial
\[
 R_{\rm top}(A)=
  (nN)^{nND_-}(n^\Delta-N^\Delta)^{nN}
           \prod_i A_i^{nNd_i}.
\]
The homogeneous component of degree \(nND-1\) is exactly
\[
 R_{\rm top}(A)\sum_i\frac{K_i}{A_i},                         \tag{1}
\]
where, with \(p=-(n/N)^\Delta/(1-(n/N)^\Delta)\) and
\(p_i=p\) for \(\epsilon_i=1\), \(p_i=1-p\) for
\(\epsilon_i=-1\),
\[
 K_i=\frac{d_i^2}{2}\bigl(nN+(n+N-2nN)p(1-p)\bigr)
                         +\frac{d_i}{2}(n+p_i)>0.            \tag{2}
\]
Thus all coefficients in this layer's support are positive.
The formula is for all adjacent orders and all unbalanced signed
gap patterns. It determines no coefficients in lower layers.

## Proof

Fix \(s_i>0\), put \(A_i=s_i/u\), and write
\(B_s=\prod_{\epsilon_i=-1}s_i^{d_i}\),
\(U_s=\prod_i s_i^{\epsilon_i d_i}\). Normalize \(F_j\) to
\[
 P_j(u,z)=
 \frac{u^{jD_-}}{j^{jD_-}B_s^j}
                   F_j(s/u,u^\Delta z/U_s).
\]
Its coefficient at \(z^k\) is
\[
 \binom jk j^{\Delta k}\prod_i
       \prod_{h=1}^{d_i\ell_{j,i}(k)}(1+u h/(j s_i)).
\]
At \(u=0\) the polynomial is \((1+j^\Delta z)^j\).
The logarithmic derivative at \(u=0\) of the factor multiplying
\(\binom jk j^{\Delta k}\) in that coefficient is
\[
 \sum_i\frac{d_i^2\ell_{j,i}(k)^2+d_i\ell_{j,i}(k)}{2j s_i}.
                                                               \tag{3}
\]
The resultant is nonzero at \(u=0\), so its first logarithmic
derivative is defined there.

Differentiate the resultant in each polynomial separately.
Holding \(P_N\) fixed at \(u=0\), its root \(-N^{-\Delta}\)
has multiplicity \(N\). The product formula for the resultant
shows that the contribution from \(P_n\) is
\(N(\partial_u P_n)/P_n\) at that root. The contribution from
\(P_N\) when \(P_n\) is fixed is
\(n(\partial_u P_N)/P_N\) at \(-n^{-\Delta}\).
No derivative of individual roots at a binomial power is needed.

Dividing a binomial sum at \(w=j^\Delta z\) by \((1+w)^j\)
gives the formal binomial generating polynomial
\((1-q+q x)^j\), where \(q=w/(1+w)\). Thus the first and
second moments of the index are \(jq\) and
\(j^2q^2+jq(1-q)\), valid for every \(w\ne-1\).
For the lower resultant contribution \(q=p\), and for the
upper contribution \(q=1-p\). A decreasing row substitutes
\(j-k\) and therefore swaps these two parameters.
Apply these moment identities to (3). The coefficient of
\(1/s_i\) in the total logarithmic derivative is
\[
 \frac{d_i^2}{2}
       \bigl(nN[p_i^2+(1-p_i)^2]+(N+n)p_i(1-p_i)\bigr)
     +\frac{d_i}{2}\bigl(Np_i+n(1-p_i)\bigr)=K_i .
\]
The variable and scalar normalizations multiply the resultant by
a constant in \(s\) times \(u^{nND}\). Hence this derivative
is the ratio of its degree-\((nND-1)\) layer to its top
layer evaluated at \(s\), which proves (1) on the positive
orthant and thus as a polynomial identity.

Finally, \(p(1-p)<0\) and \(n+N-2nN=1-2n^2<0\).
Since \(\Delta\) is a nonzero integer, the smaller of \(p\)
and \(1-p\) is at least \(-n\): it has the form
\(-r/(1-r)\) for \(r=(n/N)^{|\Delta|}\le n/N\).
This proves the positivity in (2).
