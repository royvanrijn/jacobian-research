# A larger growing window for prime-discriminant core gaps

Fix \(r\in\{3,7\}\) and sample \(p\) uniformly from the primes \(p\le X\)
with \(p\equiv r\bmod8\). Write \(E_p\) for the eligible gaps in the
split-prime core. Let \(E^{(r)}\) be the independent fair odd-prime model
with 2 admitted precisely when \(r=7\).
For every fixed \(0<\theta<2/19\), the law of \(E_p\cap[2,X^\theta]\)
converges in total variation to the law of the entire \(E^{(r)}\).

## Improving the monomial error

Use the disjoint supports and the normalized centered tenth-moment
polynomial of the [original transfer](RANDOM_SPLIT_PRIME_SUMSETS.md#transfer-to-a-power-sized-window-of-discriminants).
The number of disjoint supports satisfies
\(K\gg n^{9/40}/(\log n)^2\). The comparison moment is \(O(K^{-5})\).
Every sign monomial has odd squarefree modulus \(Q\le n^{20}\), with
a bounded number of factors, and their coefficient \(\ell^1\) norm is
bounded. Discard \(p\le X^\theta\), an event of probability \(o(1)\).
For the remaining primes and \(n\le X^\theta\), the original signs equal
the retained extensions
\(t_q(a)=(-1/q)(a/q)+\mathbf 1_{q\mid a}\) evaluated at \(a=p\).
A gap forces every disjoint support to fail.

For \(n\le L=(\log X)^{16}\), the original fair-bit comparison by
Siegel–Walfisz is unchanged: the total error over these targets tends to
zero. For \(L<n\le X^\theta\), use the original majorization of the
nonnegative moment on primes by \(O(\log X)\) times the moment over all
integers \(a\le X\) in the fixed class. Its reference bits have means
\(1/q\). Expand the divisibility terms in each sign monomial.
The all-divisibility term has mean \(1/Q+O(1/X)\). Each other term
has \(d\mid a\) for a divisor \(d\mid Q\) and a primitive quadratic
character of squarefree conductor \(e=Q/d>1\). On writing
\(a=d(8v+b)\) for the appropriate residue \(b\), this character is a
constant times a primitive interval character sum in \(v\) shifted
modulo \(e\). The length is \(O(X/d)\), with an empty interval when
\(d>X\).

Apply Burgess with moment parameter \(2\), as stated in
[Heath-Brown, equation (1)](https://arxiv.org/pdf/1203.5219).
For any fixed \(\nu>0\), after normalization by the number of integers
in the class, the error of the term is
\[
 \ll_\nu X^{-1}(X/d)^{1/2} e^{3/16+\nu}
 \le C_\nu X^{-1/2} Q^{3/16+\nu}.
\]
The character is primitive since each odd prime in the squarefree
conductor contributes a nontrivial quadratic character. The change of
residue modulo 8 is invertible modulo the conductor.

The union bound on progression-average errors through \(X^\theta\) is
\[
       O_\nu\!\left(\log X\,X^{-1/2+\theta(19/4+20\nu)}\right)=o(1)
\]
on choosing \(\nu\) small enough. The comparable sum of
\(O(\log X\,K^{-5})\) terms above \(L\) tends to zero since
\(K^{-5}\ll n^{-9/8}(\log n)^{10}\).
For \(B<n\le L\), use the original per-target estimate
\(O(K^{-5})+O_A((\log X)^{-A})\) and let \(B\to\infty\).
Together these bounds give
\[
 \lim_{B\to\infty}\limsup_{X\to\infty}
     \Pr(E_p\cap(B,X^\theta]\ne\varnothing)=0 .
\]
The fixed-window law and the almost surely finite independent-model
gap set now give the asserted total variation convergence by truncation.

## Equivalent remaining tail condition

For any such fixed \(\theta\) and any fixed \(\epsilon>0\), the
[uniform polynomial bound](PRIME_DISCRIMINANT_CORE_UNIFORM.md#an-effective-polynomial-bound)
confines the full arithmetic gaps to \(n\le C_\epsilon X^{3+\epsilon}\).
Full total variation convergence of \(E_p\) is equivalent to
\[
 \Pr\big(E_p\cap(X^\theta,C_\epsilon X^{3+\epsilon}]
                             \ne\varnothing\big)\longrightarrow0 .
\]
For sufficiency, couple the full and truncated arithmetic sets on the
same sampled prime. They differ only on the displayed tail event.
For necessity, under full total variation convergence the probability
of any gap beyond \(X^\theta\) differs from its independent-model
counterpart by at most the total variation distance. The latter
counterpart tends to zero by almost sure finiteness. This leaves the
stated intermediate range as the exact remaining obligation.
