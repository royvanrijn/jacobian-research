# The squarefree twin correlation near zero

Write \(\omega(n)\) for the number of distinct prime divisors and set
\[
 T_X(s,t)=\sum_{\substack{3\le n\le X\\ n\ {\rm odd}}}
       \mu^2(n)\mu^2(n+2)s^{\omega(n)}t^{\omega(n+2)},\qquad
 T_X(t)=T_X(t,t).
\]
Every pair has coprime entries. Thus for \(0\le t\le1\), \(T_X(t)\)
is the expected number of these squarefree pairs whose entries both
belong to the random semigroup with odd primes independently selected
with probability \(t\). The lower cutoff 3 removes the pair \((1,3)\),
which would add \(t\). The polynomial satisfies
\[
             T_X(t)=t^2\pi_2(X)+\sum_{j\ge3} a_j(X)t^j,
             \qquad a_j(X)\ge0,                       \tag{1}
\]
where \(\pi_2(X)=|\{p\le X:p,p+2\text{ prime}\}|\).
Indeed each entry has at least one prime factor and each squarefree
entry has exactly one if and only if it is prime.

**Theorem.** There are absolute constants \(C,c>0\) and an integer
\(J\ge3\) such that for all large \(X\) and \(|s|,|t|\le1\),
\[
           |T_X(s,t)|\le C|st|\,\frac{X}{(\log X)^{2-|s|-|t|}}.   \tag{2}
\]
For real \(0<t\le1\),
\[
            T_X(t)\ge c\frac{X}{(\log X)^2}
                     \max\{t^3,\ t^J(\log X)^{2t}\}.             \tag{3}
\]
In particular \(T_X(t)\asymp_t X(\log X)^{2t-2}\) at every fixed
positive \(t\le1\). The first term of (3) is a consequence of Chen's
theorem on a prime with a prime or semiprime at distance two.

**Theorem.** Write \(L_X=\log\log X\). Uniformly for
\(0<|t|\le 1/L_X\),
\[
       \frac{T_X(t)}{t^2}=\pi_2(X)+
                  O\!\left(|t|L_X\,\frac{X}{(\log X)^2}\right).
                                                               \tag{4}
\]
After removing the singularity at \(u=0\), the polynomials
\[
      H_X(u)=\frac{(\log X)^2}{X}\,
                     \frac{T_X(u/L_X)}{(u/L_X)^2}
\]
form a normal family on compact sets as \(X\to\infty\).
Any subsequential limit is entire, has nonnegative Taylor
coefficients, and satisfies \(|H(u)|\le C e^{2|u|}\).
In particular, if \(H_X(u)\) has pointwise limits on a fixed
real interval of positive length, then it converges locally
uniformly on the whole complex plane, and the
normalized twin-prime count converges to \(H(0)\).

## A sieve bound with forced prime factors

It suffices first to work on a dyadic interval \(I=(x/2,x]\)
at large \(x\), and with real \(s,t\in[0,1]\).
Set \(y=x^\eta\) for a small absolute \(\eta>0\).
For \(q\) odd below \(y\), independently declare \(0\bmod q\)
forbidden with probability \(1-s\) and \(-2\bmod q\)
with probability \(1-t\). For an odd \(n\), the probability
of avoiding all these classes is
\[
 W(n)=s^{\omega_{<y}(n)}t^{\omega_{<y}(n+2)}.
\]
There is no conflict between the two classes at an odd prime.

Take \(a,b\) each either 1 or an odd prime less than \(y\).
For \(a>1\) require \(a\mid n\), and for \(b>1\) require \(b\mid n+2\).
The case \(a=b>1\) is empty. Otherwise these conditions and oddness
define one progression modulo \(2ab\).
On that progression the corresponding \(a,b\) residue classes
must be allowed, an event of probability
\(s^{\mathbf1_{a>1}}t^{\mathbf1_{b>1}}\).
The other class at each of these primes is irrelevant. Conditional on
their allowance, sieve the progression over the remaining odd primes
below \(y\). For a realized collection of forbidden classes, let
\(g(q)\) be its number of forbidden residues divided by \(q\).
At active primes \(0<g(q)\le2/q<1\).
The number of members hitting all forbidden sets at a squarefree
\(d\) supported on active primes is
\[
                  \frac{|I|}{2ab}g(d)+O(2^{\omega(d)}),
\]
by the Chinese remainder theorem and interval counting.

Here is a uniform Selberg bound for the number of survivors.
Put \(R=x^{1/8}\), and \(h(q)=g(q)/(1-g(q))\) on the active primes.
Extend \(g,h\) multiplicatively to their squarefree products, and put
\(G(R)=\sum_{d\le R}h(d)\) over these products.
The Selberg weights \(\lambda_d\) are supported on \(d\le R\)
among these products and defined by
\[
  \sum_{\substack{d\le R\\ r\mid d}}\lambda_d g(d)
             =\frac{\mu(r)h(r)}{G(R)}.
\]
Möbius inversion gives
\(\lambda_d=\mu(d)G_d(R/d)\prod_{q\mid d}(1-g(q))^{-1}/G(R)\),
where \(G_d\) restricts the sum to products coprime to \(d\).
Thus \(\lambda_1=1\), and \(|\lambda_d|\le1\), since
\(G_d(R/d)\prod_{q\mid d}(1+h(q))\le G(R)\).
The identity \(1/g(\gcd(d,e))=\sum_{r\mid\gcd(d,e)}1/h(r)\)
gives
\(\sum_{d,e}\lambda_d\lambda_e g(\operatorname{lcm}(d,e))=1/G(R)\).
On a survivor, the squared weight
\((\sum_{d:\ n\text{ hits the sets for }d}\lambda_d)^2\) is 1,
so summing it bounds the survivor count above by
\[
                O\!\left(\frac{x}{abG(R)}+R^2(1+\log R)^2\right).
\]
The remainder follows from
\(2^{\omega(\operatorname{lcm}(d,e))}\le2^{\omega(d)}2^{\omega(e)}\)
and \(\sum_{d\le R}2^{\omega(d)}\le\sum_{d\le R}\tau(d)\).

Under the full product weight \(\prod_q(1+h(q))\), a squarefree \(d\)
is obtained by independent inclusion of active primes with probabilities
\(g(q)\). Its expected logarithm is at most
\(2\sum_{q<y}(\log q)/q\ll\log y\).
Choose \(\eta\) so that this is at most \((\log R)/2\).
Markov gives \(G(R)\ge\frac12\prod_q(1+h(q))\), hence
\(1/G(R)\ll\prod_q(1-g(q))\). The prime sum used here follows
from Mertens' product bound, by partitioning successively at
square roots.
[Rosser–Schoenfeld, Theorems 7–8](sources/rosser_schoenfeld_1962.pdf)

Now average over the forbidden classes. By independence across primes
the mean product is
\[
 \prod_{\substack{3\le q<y\\q\nmid ab}}
       (1-(2-s-t)/q)\ll(\log x)^{s+t-2},
\]
uniformly for \(s,t\in[0,1]\). Comparing the factors to
\((1-1/q)^{2-s-t}\) has a uniformly absolutely convergent ratio.
Omitting at most two primes changes the bound by a bounded factor.
Including the allowance probability gives
\[
 \sum_{\substack{n\in I,\ n\ {\rm odd}\\a\mid n,\ b\mid n+2}}W(n)
 \ll s^{\mathbf1_{a>1}}t^{\mathbf1_{b>1}}
       \left(\frac{x(\log x)^{s+t-2}}{ab}
                      +x^{1/4}(1+\log x)^2\right).             \tag{5}
\]

For squarefree \(m\) equal to \(n\) or \(n+2\) in this interval,
and \(0<v\le1\),
\[
 v^{\omega(m)}\le v^{\omega_{<y}(m)}
       \left(v+\frac{2}{\log x}\sum_{\substack{q<y\\q\mid m}}\log q\right).
\]
A prime factor at least \(y\) gives the first term. If all prime
factors are below \(y\), the prime logarithms sum to \(\log m\).
Apply the inequality to both entries, expand the product and use
(5), dropping squarefreeness from its right side.
Every term carries \(st\), either from the first terms in the brackets
or from the allowance probabilities when a divisor is specified.
The main term is \(O(st\,x(\log x)^{s+t-2})\), since
\(\sum_{q<y}(\log q)/q\ll\log x\).
The combined remainder is
\(O(st\,x^{1/4}(1+y)^2(1+\log x)^2)\), which is smaller for
fixed sufficiently small \(\eta\).
Sum the intervals with \(x\ge\sqrt X\) and bound the remaining
terms by \(st\sqrt X\). Nonnegative coefficients extend the estimate
to complex \(s,t\) by absolute values, and zero parameters by continuity.
This proves (2).

## Lower bounds for the weighted count

For the diagonal and \(n\in(x/2,x]\) odd, exclude each odd prime
below \(y\) with probability \(1-t\), independently across primes.
At an excluded prime forbid both \(n\equiv0,-2\bmod q\).
At an allowed prime forbid both \(n\equiv0,-2\bmod q^2\).
For a fixed configuration, this is a sieve with local proportions
\(g(q)=2/q\) or \(2/q^2\). Label each \(n\) by the product of
its offending primes below \(y\). For squarefree \(d\), the count
of labels divisible by \(d\) is the main term with density \(g(d)\)
and remainder \(O(2^{\omega(d)}d)\). This follows by counting the
offending classes modulo \(d^2\), at most \(2^{\omega(d)}d\) in all.
The dimension product bound is uniform because \(g(q)\le2/q\).
Apply the lower-bound fundamental lemma at level \(x^{1/4}\),
taking \(\eta\) smaller if needed. The total remainder
\(O(x^{1/2}(1+\log x))\) is negligible against
\(x\prod_q(1-g(q))\gg x/(\log x)^2\).
Delete also the entries with a prime square divisor at least \(y^2\).
At most \(O(x/y+\sqrt x)\) pairs are deleted, uniformly over
configurations. Each surviving pair has at most
\(J=\lceil 3/\eta\rceil\) further prime factors in total.
[Cao–Huang–Zhang, Theorem 3.2](sources/cao-huang-zhang-2509.23320v1.pdf)

For a squarefree pair the probability of surviving the random
conditions below \(y\) is \(t^{\omega_{<y}(n(n+2))}\), so multiplying
the expected survivor count by \(t^J\) lower bounds its full
weighted total. The average of the sieve product is
\[
 \prod_{3\le q<y}\bigl(1-2/q+2t(1/q-1/q^2)\bigr)
                        \asymp (\log x)^{2t-2}
\]
uniformly for \(t\in[0,1]\), by the same Mertens comparison.
This gives the second term of (3) through \(X=x\).
For the first term, Chen's theorem gives \(\gg X/(\log X)^2\)
primes \(p\le X\) for which \(p+2\) is prime or a product of two
primes counted with multiplicity. The prime-square case costs
\(O(\sqrt X)\) choices. Each remaining pair has degree two or three
and thus weight at least \(t^3\). Equivalently,
\(a_2(X)+a_3(X)\gg X/(\log X)^2\), with \(a_2(X)=\pi_2(X)\).
[Halberstam's proof of Chen's theorem, p. 281](sources/halberstam-chen-twin-1975.pdf)

## The scale that isolates twin primes

Take \(r=1/L_X\) in (2). Then
\(T_X(r)/r^2\ll X/(\log X)^2\), since
\((\log X)^{2r}=e^2\).
By (1), for \(|t|\le r\) the contribution of degree at least three
after division by \(t^2\) has absolute value at most
\((|t|/r)T_X(r)/r^2\). This proves (4).
The same upper estimate at \(t=u/L_X\) bounds \(H_X(u)\) by
\(C e^{2|u|}\) on each compact set for all sufficiently large \(X\);
continuity covers zero. Montel's theorem gives subsequential limits
and Cauchy's formula preserves nonnegative coefficients.
If the real-interval limits exist, all subsequential limits agree
there and hence everywhere by the identity theorem. In particular
\(H_X(0)=(\log X)^2\pi_2(X)/X\) converges.

If \(t_X=o(1/\log\log X)\) is positive, then
\[
 \frac{T_X(t_X)}{t_X^2 X/(\log X)^2}
             =\frac{\pi_2(X)}{X/(\log X)^2}+o(1).
\]
Thus a positive limiting lower bound on the left gives a
lower bound of twin-prime scale. Convergence on a fixed real
\(u\)-interval to \(S e^{2u}\), where \(S\) is the usual
Hardy–Littlewood twin-prime constant, would give the full
Hardy–Littlewood asymptotic. The fixed-\(t\) lower sieve and
Chen bound give smaller powers of \(t\) in this normalization.
The existing [independent-linear-form correlation theorem](MULTIPLICATIVE_HARDY_LITTLEWOOD.md)
requires pairwise linearly independent coefficient rows. The
shifted forms \(n,n+2\) have the same coefficient row.

**Open question.** Obtain a positive lower bound at the scale
\(t_X^2 X/(\log X)^2\) for any positive
\(t_X=o(1/\log\log X)\), or identify a limit with \(H(0)>0\)
on the rescaled parameter window. Equation (4) makes the former
bound equivalent to a lower bound of twin-prime scale.
