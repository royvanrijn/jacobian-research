# Shifted local correlations and additive energy in two semigroup models

The binary prescribed-prime singular series and the local intensity of
a random prime semigroup share exactly the same factor at an excluded
prime when \(k=2\). This gives an explicit mixed shifted correlation
between different prescribed-prime semigroups. In the random model,
the same calculation identifies the exact energy constant needed for
a second-moment route to the conditional intensity law.

## Local correlation of the random-prime intensity

Use the notation of the [critical random-prime model](RANDOM_PRIME_SEMIGROUP_LOCAL_INTENSITIES.md).
Fix \(k\ge2\), \(\gamma=1-1/k\) and \(\lambda=1\), with \(2\) admitted.
In a fixed realization of the prime selection, put \(c_p=b_p^k\), where
\(b_p=(p t_p-1)/(p-1)\). Write \(C_p(h)=\mathbb E_U[g_p(U_p)g_p(U_p+h)]\),
where \(U_p\) is Haar uniform on \(\mathbb Z_p\). Put \(v_p(0)=+\infty\).
The local and global formulas are
\[
\begin{split}
C_p(h)&=
\begin{cases}
1+c_p^2\{p-(p+1)p^{-v_p(h)}\},&p\ \mathrm{selected},\\
1+c_p^2\{p\mathbf1_{p\mid h}-1\},&p\ \mathrm{excluded},
\end{cases}\\
\mathbb E_U[S(U)S(U+h)]&=\prod_p C_p(h).
\end{split}                                                          \tag{10}
\]
In particular \(S\in L^2\) almost surely and
\[
 \mathbb E_U S^2
 =\prod_{p\ \mathrm{selected}}(1+p b_p^{2k})
  \prod_{p\ \mathrm{excluded}}(1+(p-1)b_p^{2k}).                         \tag{11}
\]

To verify this, at conductor \(p^\ell\) the Fourier coefficient of \(g_p\)
is \(b_p^k t_p^{k(\ell-1)}\). For a selected prime this equals
\(c_p p^{-(\ell-1)}\). For an excluded prime only conductor \(p\)
contributes besides the constant. Each \(g_p\) is square integrable
because \(v_p(U_p)\) has a geometric tail. Parseval gives the
correlation by summing squared single-character coefficients.
At a selected prime the nontrivial part is
\[
 c_p^2\sum_{\ell\ge1} p^{-2(\ell-1)}
 \bigl(p^\ell\mathbf1_{p^\ell\mid h}
              -p^{\ell-1}\mathbf1_{p^{\ell-1}\mid h}\bigr).
\]
If \(v=v_p(h)<\infty\), the sum without \(c_p^2\) is
\((p-1)\sum_{j=0}^{v-1}p^{-j}-p^{-v}=p-(p+1)p^{-v}\).
At \(h=0\) it is \(p\). This proves the local factors including the
excluded-prime case. The selected local variances are \(O_k(p^{-1})\),
which are summable almost surely over selected primes by their
selection probabilities. The excluded variances are \(O_k(p^{1-2k})\),
which are summable for \(k\ge2\). Independence and the mean-one
property imply the finite products are Cauchy in \(L^2\): the square of
the difference between nested products equals the difference between
their second moments. Their \(L^2\) limit is the same \(S\) as the
previous \(L^1\) product. Translation preserves Haar measure, so
Cauchy–Schwarz gives (10) from the finite products.

## Sharp second-moment route to the random-prime intensity law

On the cutoff and grids of [equation (6)](RANDOM_PRIME_SEMIGROUP_LOCAL_INTENSITIES.md#a-sufficient-condition-for-the-critical-density),
write \(\mu_{n,X}=\mathbb E[r_k(n)\mid\mathcal F_Y]\) and
\[
 E_X=\frac1X\sum_{n\le X}\mu_{n,X}^2,\qquad
 E_\omega=B^2\mathbb E_U S(U)^2 .
\]
Almost surely, the retained fixed-residue first-moment convergence
gives
\[
                  \liminf_X E_X\ \ge\ E_\omega.                       \tag{12}
\]
A matching upper bound \(\limsup_X E_X\le E_\omega\) on those grids
implies the empirical conditional intensity convergence in (6),
and gives convergence of its second moment as well.

Indeed, for every fixed modulus \(m\), the conditional
residue-class mean of \(\mu_{n,X}\), normalized by \(m/X\), tends to
\(B S_m(a)\). This is the residue-restricted concentration and mean
limit used in [the proof of (4)](RANDOM_PRIME_SEMIGROUP_LOCAL_INTENSITIES.md#stronger-almost-sure-critical-gap-bound).
Consequently
\[
 \frac1X\sum_{n\le X}|\mu_{n,X}-B S_m(n)|^2
   = E_X-B^2\mathbb E_U S_m(U)^2+o_m(1).
\]
Nonnegativity gives the lower bound with \(S_m\). Taking a divisible
sequence of moduli generating all profinite coordinates gives
\(S_m=\mathbb E[S\mid U\bmod m]\to S\) in \(L^2\) and proves (12).
If the sharp upper bound holds, the iterated limit of the displayed
mean-square discrepancy is zero. For each bounded Lipschitz test,
Cauchy–Schwarz bounds the discrepancy of empirical averages by its
Lipschitz constant times the root mean-square discrepancy. The
empirical law of \(S_m(n)\) tends to its Haar law by periodicity.
Passing to the divisible sequence proves weak convergence to
\(\operatorname{Law}_U(BS)\). The second-moment convergence follows
from (12) and the sharp upper bound.

The statistic \(E_X\) has an exact additive-energy interpretation.
Complete the random prime selection above \(Y\) in two conditionally
independent ways, retaining the same primes through \(Y\). If
\(r_k^{(i)}\) is the representation count for completion \(i\), then
\[
 E_X=\frac1X\mathbb E\!\left[\sum_{n\le X}r_k^{(1)}(n)r_k^{(2)}(n)
                                   \,\middle|\,\mathcal F_Y\right].
                                                                         \tag{13}
\]
For \(k=2\) this is a four-summand additive energy with different
completions on the two sides of the sum constraint. The weight of
one pair is the selection probability of the union of its
unexposed prime supports. Shared primes within a pair prevent
replacement by a product of two one-summand weights.
The single-completion overlap condition in (6) is still needed.
The sharp upper bound in (12) is a sufficient route for the
intensity law, which could also follow from a weaker argument that
does not control the empirical second moment.

## Mixed shifted binary counts for prescribed-prime semigroups

Let \(q_1,q_2\ge3\) be fixed. For \(i=1,2\) put
\(L_i=\operatorname{lcm}(2,q_i)\),
\(\mathcal B_i=\{p:p\nmid L_i,\ p\not\equiv1\bmod q_i\}\), and write \(S_i\)
for the binary singular series in [(V1)](PRESCRIBED_PRIME_BINARY_VARIANCE.md#definitions-and-the-proved-fixed-width-theorem).
For an odd prime set \(d_p(h)=p\mathbf1_{p\mid h}-1\) and
\(\sigma_p(h)=1+d_p(h)/(p-1)^2\).
The exact mixed kernel at shift \(h\) is
\[
\begin{split}
 K_{12}(h)
  ={}&\gcd(L_1,L_2)\mathbf1_{\gcd(L_1,L_2)\mid h}
     \prod_{p\in\mathcal B_1\cap\mathcal B_2}
                   \left(1+\frac{d_p(h)}{(p-1)^4}\right)\\
  &{}\times
     \prod_{\substack{p\in\mathcal B_1\\p\mid L_2}}\sigma_p(2-h)
     \prod_{\substack{p\in\mathcal B_2\\p\mid L_1}}\sigma_p(2+h).
\end{split}                                                            \tag{14}
\]
In particular if \(q_1=q_2=q\), this reduces to
\(L\mathbf1_{L\mid h}\prod_{p\in\mathcal B}(1+d_p(h)/(p-1)^4)\).
Its diagonal is the known four-form constant \(K_q\).

To prove (14), restrict each singular series to finitely many of its
bad primes, including any bad prime dividing the other modulus.
The two compatible progressions require \(N\equiv2\bmod L_1\)
and \(N+h\equiv2\bmod L_2\). Their joint Haar mass, multiplied by
\(L_1L_2\), is \(\gcd(L_1,L_2)\) if the congruences are compatible
and zero otherwise. For \(p\in\mathcal B_1\) dividing \(L_2\),
the second congruence supplies \(N\equiv2-h\bmod p\), and the
reverse case supplies \(N+h\equiv2+h\bmod p\).
Every remaining exclusive bad prime contributes its mean \(1\).
For a shared bad prime the factor \(1+d_p/(p-1)^2\) is precisely
the \(k=2,t_p=0\) local law in (10). Its shifted correlation
is \(1+d_p(h)/(p-1)^4\). CRT multiplies the factors.

For an arithmetic average truncate the bad-prime products at
\(z=(\log X)^{1/2}\). The dyadic divisor-sum argument of
[(6) in the variance reduction](PRESCRIBED_PRIME_BINARY_VARIANCE.md#the-exact-power-window-variance-reduction)
bounds each squared \(L^2\) truncation error by \(O_{q_i}(1/z)\).
The same bound applies to \(S_i(N+h)\) uniformly for \(|h|\le X^\eta\),
\(\eta<1\), since the shifted ranges lie in \([1,3X]\)
and the divisor count in that interval is at most \(3X/d\).
Each truncated series is bounded by
\(L_i\prod_{3\le n\le z}(1+1/(n-1))=O_{q_i}(z)\).
The common period of the two truncated factors is
\(\exp(O_{q_1,q_2}(z\log z))\), so its average in the dyadic sample
has negligible endpoint error. The products in (14) converge
uniformly in \(h\) with tail \(O(z^{-2})\) by absolute summability.
Cauchy–Schwarz on the \(L^2\) tails yields
\[
 \frac1X\sum_{X<N\le2X} S_1(N)S_2(N+h)
      =K_{12}(h)+O_{q_1,q_2,\eta}((\log X)^{-1/4}).                     \tag{15}
\]

Let \(5/6<\theta_i<1\) be fixed, \(H_i=X^{\theta_i}\), and define
\(R_{i,X}=r_{q_i,H_i}/(H_i\mu_{i,X}^2)\) using the normalization of
[the power-window variance theorem](PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md).
In these widths one obtains the corresponding **actual count**
correlation
\[
 \frac1X\sum_{X<N\le2X}R_{1,X}(N)R_{2,X}(N+h)
     =K_{12}(h)+O_{q_1,q_2,\theta_1,\theta_2,\eta}((\log X)^{-1/8}),
 \qquad |h|\le X^\eta .                                               \tag{16}
\]
The power-window estimate gives normalized squared \(L^2\) errors
\(O((\log X)^{-1/4})\) on the original interval. The shifted interval
adds at most \(2|h|\) boundary terms. Pointwise
\(r_{q_i,H_i}(M)\le H_i+2\), and \(S_i(M)\ll_{q_i,\epsilon}M^\epsilon\)
for each \(\epsilon>0\) by bounding each large-prime factor by \(p^\epsilon\).
Thus the squared \(L^2\) cost of the added boundary terms is negligible,
taking \(2\epsilon<1-\eta\). The normalized counts and model have bounded
\(L^2\) norms, and Cauchy–Schwarz transfers (15) to (16).

The model correlations are available at every width. The actual-count
conclusion uses the proved range of the power-window variance. In the
shorter range the absolute diagonal at \(q_1=q_2,h=0\) still requires
the sharp energy bound in
\(\texttt{GAP-ARITH-FBS-SHARP-FOUR-FORM-ENERGY}\).
