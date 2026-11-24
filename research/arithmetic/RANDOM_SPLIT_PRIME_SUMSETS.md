# Random prime semigroups and growing windows for quadratic splitting cores

Choose odd primes independently, and let \(C\) be the multiplicative
semigroup they generate, including 1. The status of 2 is deterministic.
An eligible target is \(n\ge2\), restricted to even targets if 2 is
excluded. Write \(E=\{n\text{ eligible}:n\notin C+C\}\) and
\(M=\sup(E\cup\{0\})\).

**Theorem.** If every odd prime is selected with the same probability
\(\rho>0\), then for every fixed \(0<\delta<1/4.267\) there are
positive constants \(c,C_0\) such that for every eligible \(n\),
\[
      \Pr(n\in E)\le C_0\exp(-c n^\delta).
\]
The random set \(E\) is almost surely finite and both \(|E|\) and \(M\)
have finite moments of every order.

**Theorem.** Suppose odd prime \(q\) is selected with probability
\(q^{-\gamma}\), for \(\gamma>0\). For \(0<\gamma<1/2\), the same
stretched-exponential bound holds for every
\[
 0<\delta<\min\left(1-2\gamma,\frac{1-\gamma}{4.267}\right).
\]
The selected primes have relative density zero almost surely.
With 2 admitted, \(C\) is almost surely an asymptotic additive basis
of order two. At \(\gamma=1/2\) the eligible gap set almost surely
has positive upper density. For \(\gamma>1/2\), \(C+C\) has density
zero almost surely. These negative conclusions also hold with 2 excluded.

## Rough decompositions with controlled supports

Put \(b_2=4.267\), fix \(0<\theta<1/b_2\), and set \(z=n^\theta\).
For even \(n\), choose \(a\) odd in \([n/3,n/2]\).
For odd eligible \(n\), choose \(a\equiv n-2\bmod4\) in that
interval, so \(a\) is odd and \(b=n-a\) has dyadic valuation one.
In both branches \(a,b\) have size comparable to \(n\). Initially
retain the pairs for which \(ab\) is free of odd prime factors
less than \(z\).

For an odd prime \(q\), the number \(\omega_n(q)\) of omitted
residues for \(a\) is two unless \(q\mid n\), in which case it is
one. For every odd squarefree \(d\) the number of progression
points with \(d\mid a(n-a)\) is
\[
        N\,\frac{\omega_n(d)}d+O(\omega_n(d)),
\]
where \(N\asymp n\) is the number of points and \(\omega_n\)
is extended multiplicatively. This follows by counting each root
in the progression, whose step is coprime to \(d\).
The two-dimensional sieve hypothesis holds with an absolute
constant. Indeed \(g(q)=\omega_n(q)/q\le2/q<1\), and for odd \(q\)
\((1-2/q)^{-1}=(1-1/q)^{-2}(1+O(q^{-2}))\).
Mertens' product estimate in a ratio of two intervals then gives
\[
 \prod_{w\le q<z}(1-g(q))^{-1}
 \le(\log z/\log w)^2(1+O(1/\log w)),
 \qquad
 V(z):=\prod_{q<z}(1-g(q))\gg(\log z)^{-2}
\]
over odd primes. The latter bound uses the convergent product of
\((1-2/q)/(1-1/q)^2\).
[Rosser–Schoenfeld, Theorems 7–8](sources/rosser_schoenfeld_1962.pdf)

Use the two-dimensional Diamond–Halberstam–Richert lower bound
with \(D=n^\xi\), where \(b_2\theta<\xi<1\).
The lower sieve function is positive at \(s=\xi/\theta>4.267\).
The remainder weights are at most \(8^{\omega(d)}\), and
\(\omega_n(d)\le2^{\omega(d)}\) for squarefree \(d\). Their
product is bounded by \(\tau_{16}(d)\).
For any fixed integer \(k\), \(\sum_{d<D}\tau_k(d)\le
D(1+\log D)^{k-1}\), by counting factorizations and dropping
the product constraint on the first \(k-1\) factors. Thus the
remainder is \(o(n/(\log n)^2)\), uniformly in \(n\).
The lower bound gives \(\gg_\theta n/(\log n)^2\) pairs.
[Johnston, v3, (1.12) and (6.11)–(6.12), applying Diamond–Halberstam Theorem 9.1](sources/johnston-2510.10853v3.pdf)

Delete pairs with an odd square factor in either summand or
an odd common prime factor. Square-factor deletions cost
\(O(n/z+\sqrt n)\). A common odd prime must divide \(n\) and
be at least \(z\), so there are at most \(1/\theta\) such
primes and their deletions cost \(O_\theta(n/z+1)\). The resulting
family \(\mathcal D_n\) still has size \(\gg_\theta n/(\log n)^2\).
The odd prime support of a pair is the support of the odd part
\(R_i\) of \(ab\), which is squarefree and satisfies
\(n^2/18\le R_i\le n^2\). Its number of factors is at most a
constant \(J\) depending on \(\theta\).

For each odd prime \(q\) in these supports, at most \(O(n/q+1)\)
pairs contain it, by the two residue classes for \(a\).
A greedy choice of one pair removes at most \(O_\theta(n/z)\)
pairs with overlapping support. Consequently there are
\[
                   K\gg_\theta \frac{n^\theta}{(\log n)^2}    \tag{1}
\]
pairs with pairwise disjoint odd-prime supports. In the
constant-probability model, each succeeds with probability at
least \(\rho^J\), independently across this family. The gap
probability is at most \(\exp(-\rho^J K)\).
Taking \(\delta<\theta<1/4.267\) proves the first tail bound.

## Weighted overlaps in the sparse model

For the sparse model, use all pairs in \(\mathcal D_n\).
Write \(m=|\mathcal D_n|\), \(I_i\) for success of pair \(i\),
\(W=\sum_i I_i\), and \(\mu=\mathbb E W\).
The probability of success is \(R_i^{-\gamma}\asymp n^{-2\gamma}\),
so \(\mu\asymp m n^{-2\gamma}\). The joint probability for pairs
\(i,j\) with overlapping support is
\[
          (R_iR_j)^{-\gamma}\gcd(R_i,R_j)^\gamma .
\]
If the gcd exceeds one, its \(\gamma\)-th power is at most the
sum of \(d^\gamma\) over its nontrivial divisors \(d\).
Let \(D(d)=\#\{i:d\mid R_i\}\). For any \(d\) that occurs here,
\(\omega(d)\le J\) and the CRT gives
\(D(d)\le 2^J(n/d+1)\). In particular \(D(d)\ll_\theta 1\)
if \(d>n\), and \(\sum_{d>n}D(d)\le 2^Jm\).

The ordered distinct-pair overlap sum \(\Delta\) in Janson's
inequality is therefore bounded by
\[
 \begin{aligned}
 \Delta
 &\ll n^{-4\gamma}\sum_{d>1}d^\gamma D(d)^2\\
 &\ll n^{-4\gamma}
       \left(n^2\sum_{z\le d\le n}d^{\gamma-2}+m n^{2\gamma}\right)\\
 &\ll_\gamma n^{2-4\gamma}z^{\gamma-1}+\mu .
 \end{aligned}
\]
The first sum starts at \(z\), since all nontrivial divisors
that occur have prime factors at least \(z\).
The \(d>n\) term uses the bounded multiplicities and \(d\le n^2\).

Janson's lower-tail inequality on independent Bernoulli prime
variables gives \(\Pr(W=0)\le\exp(-\mu^2/(2(\mu+\Delta)))\).
The exponent is bounded below by a constant times
\[
       \min\left(\frac{n^{1-2\gamma}}{(\log n)^2},
                 \frac{n^{\theta(1-\gamma)}}{(\log n)^4}\right).
\]
This proves the second tail bound for \(\gamma<1/2\).
Summing the failure probabilities gives almost sure finiteness.
The union bound for \(M\ge x\) preserves a stretched-exponential
tail after decreasing \(c\), since
\(\sum_{n\ge x} e^{-c n^\delta}\le
e^{-c x^\delta/2}\sum_{n\ge1}e^{-c n^\delta/2}\).
The moment claim follows from \(|E|\le M\).
[Janson–Warnke, introductory equation (2) with \(\varepsilon=1\)](sources/janson-warnke-1406.1248v1.pdf)

For any \(\gamma>0\), put \(\eta=\min(\gamma,1/2)\). The
expected number of selected primes through \(x\) is
\(O_\gamma(x^{1-\eta})\) by summing over integers. Relative
to \(\pi(x)\asymp x/\log x\), Markov's inequality gives
\(O_{\epsilon,\gamma}(x^{-\eta}\log x)\)
for exceeding any fixed relative fraction \(\epsilon>0\).
This is summable along dyadic \(x\), and monotonicity extends
relative density zero to all \(x\).

## Failure at and above the critical value

For \(\gamma>1/2\), take \(0<s<1/2\) with \(s+\gamma>1\).
The expected value of \(\sum_{a\in C}a^{-s}\) is finite.
The odd-prime Euler factors have the form
\(1+q^{-\gamma-s}/(1-q^{-s})\), and the dyadic factor is finite.
Hence \(|C\cap[1,x]|=O_\omega(x^s)\) almost surely.
The number of sumset values through \(x\) is
\(O_\omega(x^{2s})=o(x)\).

For \(\gamma=1/2\), condition on all odd primes through \(y\ge2\),
and write \(\mathcal F_y\) for the resulting sigma field.
Let \(A(x)=|C\cap[1,x]|\), and put \(f_y(a)=\Pr(a\in C\mid
\mathcal F_y)\). This function is multiplicative. For an exposed
odd prime its positive-exponent values are either all 1 or all 0.
For an unexposed prime they all equal \(q^{-1/2}\).

Write \(f_y=(n\mapsto n^{-1/2})*g_y\).
At an unexposed \(q\), \(g_y(q)=0\) and
\(g_y(q^a)=q^{-1/2}(1-q^{-1/2})\) for \(a\ge2\).
Thus \(\sum_d |g_y(d)|d^{-1/2}<\infty\) for each fixed \(y\).
The elementary estimate \(\sum_{m\le u}m^{-1/2}\sim2\sqrt u\)
and dominated convergence give
\[
  x^{-1/2}\mathbb E(A(x)\mid\mathcal F_y)\ \longrightarrow\
  c_y=2L_2
        \prod_{\substack{3\le q\le y\\q\ {\rm prime}}} L_q
        \prod_{q>y}(1+q^{-3/2}),
\]
where \(L_2=(1-1/2)\sum_{a\ge0}f_y(2^a)2^{-a/2}\) is deterministic,
and \(L_q=1+q^{-1/2}\) at a selected exposed prime and \(1-1/q\)
otherwise. The same convolution gives a bound \(O_y(\sqrt x)\)
for all \(x\ge1\) in each exposed configuration.
Since \(\mathbb E L_q=1+q^{-3/2}\), the \(c_y\) form a
nonnegative martingale along the primes. They are bounded above
almost surely by the nonnegative martingale maximal inequality.

If \(d\) is squarefree and all its primes exceed \(y\), then
\[
 \sum_{\substack{a\le x\\d\mid a}} f_y(a)
 = d^{-1/2}\sum_{\substack{e\mid d^\infty\\d\mid e}}
                   \ \sum_{\substack{u\le x/e\\(u,d)=1}}f_y(u).
\]
The limsup after division by \(\sqrt x\) is at most
\(c_y d^{-1}\prod_{q\mid d}(1-q^{-1/2})^{-1}\), by summing the
geometric series and using the fixed-configuration \(O_y\) bound.
The joint membership probability for \(a,b\) is \(f_y(a)f_y(b)\)
times \(\prod\sqrt q\) over unexposed common primes.
Expand this product using local increments \(\sqrt q-1\).
The same fixed-configuration bound permits dominated summation over
squarefree divisors, because the product below converges. It gives
\[
 \limsup_{x\to\infty}\frac{\mathbb E(A(x)^2\mid\mathcal F_y)}x
 \le C_* c_y^2,\qquad
 C_*=\prod_{q\ge3}
     \left(1+\frac{\sqrt q-1}{q^2(1-q^{-1/2})^2}\right)<\infty .
\]

In a fixed exposed configuration, the total conditional expected
number of ordered representations of \(n\le x\) is at most
\(\mathbb E(A(x)^2\mid\mathcal F_y)\).
For all sufficiently large \(x\), at least a fixed positive
fraction of the eligible integers in \((x/2,x]\) therefore
have conditional expected representation count at most
\(16(1+C_* c_y^2)\), by Markov's inequality.
Exclude targets with a representation using only the finitely
many deterministic selected primes. There are only polylogarithmically
many such targets in \(x\).
For each remaining \(n\), every possible success has conditional
probability \(w_i\le1/\sqrt3\).
The failures are decreasing events in the remaining independent
bits. Harris' inequality gives
\[
  \Pr(n\in E\mid\mathcal F_y)\ge \prod_i(1-w_i)
    \ge \exp\left(-\frac{16(1+C_*c_y^2)}{1-1/\sqrt3}\right).
\]
For all large \(x\), the conditional expected density of gaps
through \(x\) is thus at least a positive constant depending
only on \(c_y\). Boundedness of the density and reverse Fatou
give the same lower bound for the conditional expectation of
its limsup, hence for the conditional probability of positive
upper density. As \(y\) grows, this lower bound stays positive
almost surely because \(\sup_y c_y<\infty\). The conditional
probabilities converge almost surely to the indicator of the
event of positive upper density, by bounded martingale convergence
and the fact that the increasing sigma fields generate the model.
This proves the critical conclusion.
[Janson–Warnke, introductory equation (3)](sources/janson-warnke-1406.1248v1.pdf)

## Transfer to a power-sized window of discriminants

For a prime \(p\equiv r\bmod8\), \(r=3,7\), use \(E_p\) as in
the [prime-discriminant result](PRIME_DISCRIMINANT_CORE_UNIFORM.md).
Let \(E^{(r)}\) be the independent model with \(\rho=1/2\), with
2 admitted precisely for \(r=7\).

**Theorem.** Sample \(p\) uniformly from primes \(p\le X\) in either
fixed class \(r=3,7\bmod8\). Then \(E_p\cap[2,X^{1/12}]\)
converges in total variation to \(E^{(r)}\) on the space of finite
sets. In particular,
\[
 \lim_{B\to\infty}\limsup_{X\to\infty}
 \Pr\big(E_p\cap(B,X^{1/12}]\ne\varnothing\big)=0.              \tag{2}
\]

Discard \(p\le X^{1/12}\), of relative probability \(o(1)\).
Fix an eligible \(n\le X^{1/12}\) and choose the disjoint supports
in (1) with \(\theta=9/40\). They give \(K\gg n^{9/40}/(\log n)^2\).
For a prime \(p>n\), the splitting sign at odd \(q\) is
\(\epsilon_q(p)=(-1/q)(p/q)\) by quadratic reciprocity.
For a general integer \(a\equiv r\bmod8\), extend the sign by
\[
       t_q(a)=(-1/q)(a/q)+\mathbf 1_{q\mid a}\in\{-1,1\}.
\]
The extension has the same value on the retained primes. It
has mean \(1/q\) under a uniform residue modulo \(q\).

Write the success indicator of support \(S_i\) as
\(W_i=\prod_{q\in S_i}(1+t_q)/2\). For independent reference
bits, either fair bits or the bits of mean \(1/q\), write
\(\mu_i=\mathbb EW_i\). In both models \(\mu_i\ge2^{-J}\).
The disjoint supports make the \(W_i\) independent, so the
normalized centered tenth moment is \(O(K^{-5})\). In its
expansion, any term with a singleton support index has mean
zero, leaving \(O(K^5)\) terms of absolute value at most one
before division by \(K^{10}\).
On a gap at \(p>n\) every \(W_i=0\). Thus the tenth power of the
normalized centered sum is at least \(2^{-10J}\).
The sign-monomial expansion of that polynomial has coefficient
\(\ell^1\) norm at most \(2^{10}\). Since \(t_q^2=1\), each monomial
is supported on a squarefree \(Q\) with \(Q\le n^{20}\) and at
most \(10J\) prime factors.

Put \(L=(\log X)^{16}\). For \(n\le L\), compare the fair-bit
moment to the average over the retained primes. Every nonempty
sign monomial averages to zero on the reduced residues modulo
\(8Q\) in the given class modulo 8, because the character in
the odd part is nonprincipal. Siegel–Walfisz applies uniformly
to \(8Q\le8(\log X)^{320}\). Partial summation gives the
unweighted progression counts from the logarithmically weighted
statement. Summing its error over the reduced classes and
dividing by the number of primes gives
\(O_A((\log X)^{-A})\) for any fixed \(A\).
Omitting the primes at most \(X^{1/12}\) also costs less than
any such inverse logarithmic power. The polynomial coefficient
bound therefore gives a per-target gap probability
\(O(K^{-5})+O_A((\log X)^{-A})\) in this range.
[Green, Proposition 4.2 and (4.4)](sources/green-roth-primes-2005.pdf)

For \(L<n\le X^{1/12}\), majorize the prime average of the
nonnegative polynomial by \(O(\log X)\) times its average over
all positive integers \(a\le X\) with \(a\equiv r\bmod8\).
Use independent biased bits of mean \(1/q\) as the comparison.
For a monomial \(\prod_{q\mid Q}t_q(a)\), expand the factors
\(\mathbf1_{q\mid a}\). The term with all divisibility indicators
has progression mean \(1/Q+O(1/X)\), which is the CRT reference
mean up to the error. Each other term has a divisor \(d\mid Q\)
with \(d\mid a\) and a primitive quadratic character modulo
\(e=Q/d>1\). On \(a=d u\), the condition on \(u\) modulo 8
turns the sum into a constant times a character sum over an
interval, since \(8\) is invertible modulo \(e\).
Pólya–Vinogradov gives \(O(\sqrt e\log(2e))\) for that sum.
This bound follows from the Fourier expansion of the primitive
character, whose Gauss sums have absolute value \(\sqrt e\),
and the geometric-series bound summed over the nonzero
frequencies. The number of terms in the expansion is at most
\(2^{10J}\). The monomial mean error is therefore
\(O_J(\sqrt Q\log(2Q)/X)\). It follows that the per-target prime
gap probability in this range, apart from the globally discarded
primes, is
\[
      O\left(\log X\left[K^{-5}
                       +\frac{n^{10}\log(2n)}{X}\right]\right).
\]
The first term sums for \(n>L\) to
\(O(\log X\, L^{-1/8}(\log L)^{10})=o(1)\), since
\(K^{-5}\ll n^{-9/8}(\log n)^{10}\).
The second sums through \(X^{1/12}\) to
\(O(X^{-1/12}(\log X)^2)=o(1)\).
For \(B<n\le L\), the small-window bound sums to at most
\(O(\sum_{n>B}n^{-9/8}(\log n)^{10})+o(1)\).
This proves (2).

At fixed \(B\), the earlier finite-window law gives convergence
in total variation on the finite set of possible truncated gaps.
The independent gap set almost surely has finite maximum by the
first theorem. Truncating both distributions at \(B\), then using
(2) and letting \(B\) increase, proves the stated total-variation
claim.

**Open question.** To obtain convergence of the full arithmetic
exception set on the space of finite sets, it remains to control
gaps with \(X^{1/12}<n\le C_\varepsilon X^{3+\varepsilon}\)
with probability tending to zero as \(X\) grows.
The upper end uses the uniform polynomial theorem. The present
transfer argument uses monomial moduli up to \(n^{20}\) and
does not give a summable progression error on this remaining range.
