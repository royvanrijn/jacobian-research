# Doubling and small cofactors in the squarefree twin correlation

Use the diagonal \(T_X(t)\), the coefficients \(a_j(X)\), and the
uniform estimates from the [squarefree correlation note](SQUAREFREE_TWIN_CORRELATION_NEAR_ZERO.md).
Write \(B_X=X/(\log X)^2\) and \(L=\log\log X\).
The normalized polynomials \(F_X(u)=T_X(u/L)/((u/L)^2B_X)\)
form a normal family, with \(F_X(0)=\pi_2(X)/B_X\).
Pointwise convergence on any fixed real interval would determine
the limit at zero by the identity theorem.

**Theorem.** Put \(D_X(t)=8T_X(t)-T_X(2t)\).
For \(t>0\) one has \(D_X(t)\le4t^2\pi_2(X)\).
If for fixed \(c,\varepsilon>0\) the inequality
\(T_X(2c/L)\le(8-\varepsilon)T_X(c/L)\) holds along an
unbounded sequence of \(X\), then on that sequence
\(\pi_2(X)\gg \varepsilon c B_X/L\).
Conversely, a lower bound \(\pi_2(X)\gg B_X\) on any sequence
would give a fixed doubling margin below 8 there for all
sufficiently small fixed \(c>0\).

Indeed the exact identity is
\[
 D_X(t)=4t^2 a_2(X)
       -\sum_{j\ge4}(2^j-8)a_j(X)t^j.                        \tag{1}
\]
The cubic coefficient drops out. A margin \(\varepsilon T_X(t)\)
therefore gives
\(\pi_2(X)\ge \varepsilon T_X(t)/(4t^2)\gg\varepsilon t B_X\)
by Chen's theorem in the form \(a_2+a_3\gg B_X\).
For the converse, evaluate the uniform upper bound at \(r=1/L\).
For \(t=c/L\) with \(0<c\le1/2\), the sum in (1) is at most
\((2c)^4T_X(r)\ll c^2 t^2 B_X\).
A full-scale lower bound for \(a_2\), with sufficiently small
fixed \(c\), gives \(D_X(t)\gg t^2 B_X\). Use
\(T_X(t)\ll t^2 B_X e^{2c}\) to convert this to a fixed margin.

The estimates used here are compatible, as polynomial inequalities,
with a missing quadratic coefficient. In the formal polynomial
\(T_X^\circ(t)=T_X(t)-a_2t^2+a_2t^3\), positivity is preserved,
\(a_2^\circ+a_3^\circ=a_2+a_3\), and for real \(0<t\le1\),
\(tT_X(t)\le T_X^\circ(t)\le T_X(t)\).
This preserves the uniform upper bound, the Chen-scale lower bound,
and the order of growth at each fixed positive \(t\).
All nonzero degrees of \(T_X^\circ\) are at least three, so
\(T_X^\circ(2t)\ge8T_X^\circ(t)\).

## Higher-order one-sided extraction

For an odd positive integer \(r\), put \(m=r+1\) and define
\[
 \mathcal D_r=-\frac1{r!}\prod_{k=3}^{r+2}
                   \left(t\frac{d}{dt}-k\right).
\]

**Theorem.** For \(t>0\),
\[
 \begin{split}
 \mathcal D_r T_X(t)&=t^2\pi_2(X)-E_r(t),\\
 E_r(t)&=\sum_{j\ge r+3}\binom{j-3}{r}a_j(X)t^j\ge0 .
 \end{split}                                                   \tag{A1}
\]
At \(t=c/L\), for fixed \(c>0\) and \(c+m/2\le L\), one has
\[
 E_r(t)\ll e^{2c}t^2B_X\left(\frac{2ec}{m}\right)^m.             \tag{A2}
\]
In particular, one may take \(r\) to be the least odd integer at least
\((1+\varepsilon)\log L/\log\log L\) for fixed \(\varepsilon>0\).
Then \(E_r(c/L)=o(t^3B_X)\). Along any unbounded sequence at this
choice of \(r\), a positive lower bound \(\mathcal D_r T_X(t)\gg t^3B_X\)
is equivalent to \(\pi_2(X)\gg B_X/L\), with constants allowed to
depend on \(c\) and \(\varepsilon\).

On \(t^j\) the multiplier is \(-\prod_{k=3}^{r+2}(j-k)/r!\).
It is \(1\) at \(j=2\), zero at \(3\le j\le r+2\), and
\(-\binom{j-3}{r}\) at \(j\ge r+3\). The falling factorial also gives
\[
 \mathcal D_r T_X(t)
   =-\frac{t^{r+3}}{r!}\left(\frac{d}{dt}\right)^r
                         \left(\frac{T_X(t)}{t^3}\right).
\]
The bound \(\binom{j-3}{r}\le 2^{j-3}\) gives
\(E_r(t)\le2^{-r-6}T_X(4t)\). For a growing evaluation radius
\(\lambda>1\), the binomial series gives
\[
 \binom{j-3}{r}\lambda^{-j}
 \le \sum_{\ell\ge r}\binom{\ell}{r}\lambda^{-\ell-3}
 =\lambda^{-2}(\lambda-1)^{-m}.
\]
Consequently \(E_r(t)\le\lambda^{-2}(\lambda-1)^{-m}T_X(\lambda t)\).
The uniform upper sieve gives \(T_X(\lambda t)\ll
(\lambda t)^2 B_X e^{2\lambda c}\) for \(\lambda t\le1\).
Taking \(\lambda=1+m/(2c)\) proves (A2). The relevant decay condition is
\(m\log(m/(2ec))-\log L\to+\infty\), which holds at the stated order.

The leading constant \(1\) in the order is sharp for a bound from the
polynomial envelope alone. A formal addition
\(B_X t^2(2eLt/m)^m\) has nonnegative coefficients and obeys that
envelope, since \((ex/m)^m\le e^x\) for \(x\ge0\).
It leaves every fixed coefficient unchanged eventually as \(m\to\infty\).
Its first surviving degree contributes exactly
\(t^2 B_X(2ec/m)^m\) to \(E_r(c/L)\).

The same error control applies to any subpolynomial with nonnegative
coefficientwise discrepancy from \(T_X\). If it omits \(o(B_X/L)\)
twins, the difference of the two signed statistics is \(o(t^3 B_X)\)
at the stated order. Indeed (A1) applies to the discrepancy, and its
nonnegative tail is bounded by \(E_r(t)\). This includes \(U_{X,\delta}\)
below, whose omitted quadratic coefficient is at most \(X^{1/3}\).

Strict positivity of the exact statistic already suffices for
infinitude when \(r\le K_0\log L\) for fixed \(K_0\).
To see this, let \(v=A\log L/L\)
with \(A\) a sufficiently large fixed constant. The weighted lower
sieve gives \(T_X(v)\gg B_X v^J L^{2A}\) for the fixed \(J\) in
the squarefree correlation note. The degrees at most \(r+2\) contribute
at most \(2^{r+2}T_X(v/2)\ll 2^{r+2}B_Xv^2 L^A\). Degrees
above \(M=\lceil K\log L\rceil\) contribute at most
\(2^{-M}T_X(2v)\ll2^{-M}B_Xv^2 L^{4A}\). Choosing \(A\) large
relative to \(J\) and the bound on \(r/\log L\), then \(K\) large
relative to \(A\), makes both contributions \(o(T_X(v))\).
On the intervening degrees the binomial multiplier is at least one,
and replacing \(v\) by \(t\) costs at most the factor \((t/v)^M\)
in the lower bound. Thus
\[
 \frac{E_r(t)}{t^2}
 \ \gg\ B_X\exp\!\big(-O_{c,J,K_0}(\log L\log\log L)\big)
 \ \longrightarrow\ \infty.
\]
If \(\mathcal D_r T_X(t)>0\), (A1) gives
\(\pi_2(X)>E_r(t)/t^2\).

### Finite differences on a positive interval

In the normalized variable, the statistic in (A1) is the Taylor
polynomial at \(u=c\) evaluated at zero:
\[
 \frac{\mathcal D_r T_X(c/L)}{(c/L)^2 B_X}
     =\sum_{k=0}^r\frac{(-c)^k}{k!}F_X^{(k)}(c).
\]
This follows by applying Leibniz's rule to the derivative formula.
A finite-difference version uses distinct nodes \(u_0,\dots,u_r\)
in a fixed interval \([a,b]\) with \(0<a<b\). Set
\[
 w_i=\prod_{j\ne i}\frac{-u_j}{u_i-u_j},
 \qquad I_{X,r}=\sum_{i=0}^r w_i F_X(u_i).
\]
The interpolating polynomial of degree at most \(r\) has value
\(I_{X,r}\) at zero. The interpolation remainder, obtained by repeated
Rolle's theorem, equals
\(F_X^{(m)}(\xi)\prod_i(-u_i)/m!\) for some \(0<\xi<b\).
Even \(m\) and coefficient positivity give its sign.
Cauchy's estimate on the circle of radius \(m/2\) about \(\xi\),
using \(|F_X(z)|\ll e^{2|z|}\) for \(|z|\le L\), gives, when
\(b+m/2\le L\),
\[
 0\le F_X(0)-I_{X,r}
      \ll e^{2b}\left(\frac{2eb}{m}\right)^m.                    \tag{A3}
\]
Thus \(I_{X,r}\le\pi_2(X)/B_X\) with error \(o(1/L)\) at the
same order.

The Chebyshev extremal nodes
\(u_i=(a+b)/2+((b-a)/2)\cos(i\pi/r)\) give
\[
 \sum_i |w_i|=\frac{\rho^r+\rho^{-r}}2,
 \qquad \rho=\frac{\sqrt b+\sqrt a}{\sqrt b-\sqrt a}.
\]
Indeed the weights alternate, as do the values of the degree-\(r\)
Chebyshev polynomial mapped to \([a,b]\). Applying the interpolation
formula to that polynomial gives the stated absolute sum. Since the
polynomial has supremum \(1\) on the interval, its value at zero also
lower bounds the absolute weight sum of any linear formula on that
interval which is exact through degree \(r\). At the chosen order,
\(\rho^r=L^{o(1)}\). Additive errors \(O(L^{-1-\eta})\) in the
normalized values at all nodes, for any fixed \(\eta>0\), therefore
have total effect \(o(1/L)\).

**Open question.** Prove a positive lower bound of order \(1/L\) for
the normalized Taylor statistic or its finite-difference version
along an unbounded sequence. The uniform sieve bounds control the
remainder. The required sign and size of the arithmetic statistic
remain open.

## The large-divisor sum in the twin-prime estimate

Let \(I=(X,2X]\), let \(N_2(I)=\#\{p\in I:p,\ p+2\text{ prime}\}\),
and let \(\mathfrak S_2=2\prod_{q>2}(1-(q-1)^{-2})\), with \(q\) prime.
Use the classical expansion \(\Lambda(n)=\sum_{d\mid n}\mu(d)\log(n/d)\).
[Tao, equations (1)–(3) and the Euler product](https://terrytao.wordpress.com/2016/07/17/notes-on-the-bombieri-asymptotic-sieve/)

**Theorem.** For \(D=X^{1/2}/(\log X)^A\) with a sufficiently large
fixed \(A\), define
\[
 R_X(D)=
 \sum_{2\le m\le (2X+2)/D}\log m
 \sum_{\substack{d>D\\ X<md-2\le 2X\\ md-2\ {\rm prime}}}
             \mu(d)\log(md-2).
\]
Then
\[
 (\log X)^2 N_2(I)=\mathfrak S_2 X+R_X(D)+o(X/L).               \tag{B1}
\]
In particular the quantitative target on \(I\) is equivalent to
\(R_X(D)\ge-\mathfrak S_2 X+\kappa X/L\) for some fixed \(\kappa>0\)
along an unbounded sequence, with \(\kappa\) allowed to change.
Every summand in \(R_X(D)\) has \(md\) composite.

To prove (B1), expand \(\Lambda(p+2)\) in
\(\sum_{p\in I}\log p\,\Lambda(p+2)\).
Divisors \(d>D\) give \(R_X(D)\) on putting \(m=(p+2)/d\), since
\(\log m=0\) at \(m=1\). All divisors here are odd for large \(X\).
For \(d\le D\), apply the classical Bombieri–Vinogradov theorem
at each prefix length between \(X\) and \(2X\) to the class
\(-2\bmod d\). Partial summation with \(\log p\log((p+2)/d)\)
costs \(O((\log X)^2)\) in the total error, by integrating the
summed progression error. With \(A\) large the small-divisor
contribution is
\[
 \left(\int_X^{2X}\log(x+2)\,dx\right)\!
       \sum_{\substack{d\le D\\d\ {\rm odd}}}\frac{\mu(d)}{\varphi(d)}
 -X\sum_{\substack{d\le D\\d\ {\rm odd}}}
                         \frac{\mu(d)\log d}{\varphi(d)}
 +O\big(X(\log X)^{-2}\big).
\]
[Johnston, Theorem 1.1 for the classical distribution estimate](sources/johnston-2510.10853v3.pdf)

The Dirichlet series of \(g(d)=\mathbf 1_{d\ {\rm odd}}\mu(d)/\varphi(d)\)
is \(H(s)/\zeta(1+s)\), where
\[
 H(s)=\frac1{1-2^{-1-s}}
 \prod_{q>2}\frac{1-q^{-s}/(q-1)}{1-q^{-1-s}},
 \qquad H(0)=\mathfrak S_2.
\]
The Euler factors give coefficients \(h(n)\) of \(H\) with
\(\sum_n |h(n)|n^\eta<\infty\) for fixed \(0<\eta<1/2\).
The classical zero-free estimate for the Mertens function and
partial summation give, for any fixed \(K\),
\(\sum_{n\le Y}\mu(n)/n=O_K((\log Y)^{-K})\) and
\(\sum_{n\le Y}\mu(n)\log n/n=-1+O_K((\log Y)^{-K})\).
The limits are the value and negative derivative of
\(1/\zeta(1+s)\) at \(s=0\). Convolving with \(h\), and splitting
its variable at \(\sqrt D\), gives the same estimates for the two
sums over \(d\), with limits \(0\) and \(-\mathfrak S_2\).
Thus the small-divisor contribution is
\(\mathfrak S_2 X+O(X(\log X)^{-2})\).

Proper prime powers \(p+2\) contribute \(O(\sqrt X(\log X)^3)\).
For the prime terms the logarithmic weight is
\((\log X)^2+O(\log X)\) on \(I\). The uniform twin upper sieve bounds
their number by \(O(B_X)\), so replacing their weights by
\((\log X)^2\) costs \(o(X/L)\). This proves (B1).
The higher-order extraction applies on \(I\) as well, by coefficient
domination by \(T_{2X}\). Its required positive lower bound on \(I\)
therefore reduces to the one-sided estimate for \(R_X(D)\) above.
Cancellation for the unweighted Möbius function on shifted primes
does not estimate the sum over the growing family \(m\ge2\).

## Reweighting test for transfer from other shifts

Put \(I_h=\{n:n,\ n+h\in I,\ n\ {\rm odd}\}\) for even \(h\).
Write \(T_{X,h}^{I}(v)=\sum_{n\in I_h}\mu^2(n)\mu^2(n+h)
v^{\omega(n)+\omega(n+h)}\) and normalize by
\(F_{X,h}^{I}(u)=T_{X,h}^{I}(u/L)/((u/L)^2B_X)\), with the
removable value at zero understood.

**Theorem.** For each sufficiently large \(X\) there are nonnegative
prime weights \(w_X\), equal to \(1\) outside \(I\), such that the
multiplicative function
\[
 f_{v,X}(n)=\mu^2(n)v^{\omega(n)}\prod_{q\mid n}w_X(q)
 \qquad(q\ {\rm prime})
\]
has gap-two correlation on \(I\) with quadratic coefficient zero.
For every fixed \(C\), its normalized correlations on \(I_h\) satisfy
\[
 \big|F_{X,h}^{I,*}(u)-F_{X,h}^{I}(u)\big|
       \ll_C \frac{L^3}{\log X}
 \qquad (|u|\le C,\quad h\ {\rm even},\quad h\ne0,\pm2).
                                                                    \tag{B2}
\]
Moreover, for \(q\le\lfloor L\rfloor\), every residue \(a\bmod q\)
and every interval \(J\subset I\),
\[
 \sum_{\substack{p\in J\\p\equiv a\bmod q}}\big(w_X(p)-1\big)
   \ll B_X\exp(-(\log X)^{1/3}).                                  \tag{B3}
\]
Thus macroscopic prime marginals in these progressions and the other
shift correlations can agree to finer than the precision in (A3)
while the gap-two quadratic coefficient vanishes in the reweighted
model.

To construct the weights, set \(y=\lfloor L\rfloor\),
\(W=\operatorname{lcm}(1,\dots,y)\), and partition \(I\) into blocks of
length comparable to \(H=X\exp(-(\log X)^{1/3})\). Chebyshev's bound
gives \(\log W=O(y)\). In each block and reduced residue modulo \(W\),
let \(P\) count its primes and \(Q\) count the lower endpoints \(p\)
of twin pairs there, including those with \(p+2\) just outside \(I\).
Siegel–Walfisz gives \(P\asymp H/(\varphi(W)\log X)\) uniformly.
The dimension-two upper sieve in the residue class gives
\[
 Q\ll \frac H{W(\log X)^2}
             \prod_{\substack{q\mid W\\q>2}}\!(1-2/q)^{-1}.
\]
The product is over primes. One may use the Selberg argument in the
squarefree correlation note up to a small fixed power of \(X\),
omitting primes dividing \(W\). Its CRT error is \(O(X^{1/4}(\log X)^2)\)
with sieve weight cutoff \(X^{1/8}\), negligible here. Mertens'
product bound gives \(Q/P\ll(1+\log y)/\log X\).
Set the weight of every removed lower endpoint to zero and that
of every other prime in the block and residue to \(P/(P-Q)\).
No odd composite in \(I\) has a prime factor in \(I\), so the only
altered endpoint weights in \(I_h\) are prime weights. In particular
every twin pair in \(I\) has zero weight. Exact balance in full
blocks and \(q\mid W\) imply (B3), since each of at most two partial
blocks has total deleted mass \(O(H/(\log X)^2)\) by the same sieve
without a progression.
[Green, Proposition 4.2](sources/green-roth-primes-2005.pdf)

For (B2), put \(v=|u|/L\). First delete the lower twin endpoints
without rescaling. An affected term at shift \(h\ne0,\pm2\) has
two prime forms \(p,p+2\) and a squarefree form \(p+h\) or \(p-h\).
These forms are distinct and contribute a term of absolute value
at most \(v^{1+\omega(p\pm h)}\). For each sign, omit primes at most
three and primes dividing \(h(h\mp 2)\) from the Selberg sieve.
The other primes below \(X^\gamma\), for small fixed \(\gamma>0\),
have two forbidden prime classes and an additional class forbidden
with probability \(1-v\). The expected sieve product is bounded by
\(O(L^3(\log X)^{v-3})\). Indeed the omitted divisors of the nonzero
integer \(h(h\mp2)\) cost at most \(O(L^3)\) in the product, by
Mertens and splitting at \(\log X\). The prime-factor logarithm
bound in the squarefree correlation note for the third form gives
another factor \(v\), either directly or from the allowance of one
prescribed prime divisor; the logarithms of prescribed primes sum
against \(1/q\) to \(O(\log X)\). A prescribed prime in an omitted
class costs the same factor \(v\) in the original weight.
The CRT errors are negligible after summing the prescriptions.
The total deleted weight is therefore
\(O(v^2 X L^3(\log X)^{v-3})\).

Rescaling the remaining weights changes terms by a relative
\(O((1+\log y)/\log X)\). The shifted version of the two-class
upper sieve bounds their unscaled mass by
\(O(v^2 X L^2(\log X)^{2v-2})\), uniformly in nonzero \(h\).
Here one omits primes dividing \(h\), at Euler product cost \(O(L^2)\).
If both prescribed primes coincide at a prime dividing \(h\),
the progression has that prime only once in its modulus, and the
original weight pays \(v^2\) for its two occurrences.
The corresponding sum uses
\(\sum_{q<X^\gamma}(\log q)^2/q\ll(\log X)^2\).
These are the forced-prime estimates from the squarefree
correlation note with the indicated local changes.
Divide the two loss estimates by \(v^2B_X\) to obtain (B2).
The assertion at \(u=0\) follows by continuity.

The construction changes the arithmetic weights. It demonstrates
that a transfer based only on these progression marginals and
off-target correlations requires additional information about
gap two, even at the extraction precision.

## Reduction to a two-prime equation

For \(0<\delta\le1/6\), let \(U_{X,\delta}(t)\) be the subpolynomial
of \(T_X(t)\) from pairs
\[
 n=ap,\qquad n+2=bq,\qquad
 a,b\le X^\delta,\qquad p,q>X^{1/3}\text{ prime}.               \tag{2}
\]
The odd cofactors \(a,b\) are squarefree and coprime, and their
factorizations in (2) are unique because \(X^\delta<X^{1/3}\).
Equivalently \(U\) counts solutions to \(bq-ap=2\) in this range,
weighted by \(t^{2+\omega(a)+\omega(b)}\).

**Theorem.** Uniformly for \(0<\delta\le1/6\), \(0<t\le1\),
\[
 0\le T_X(t)-U_{X,\delta}(t)
       \ll \frac{t^3}{\delta}\frac{X}{(\log X)^{2-2t}}
                 +t^2\sqrt X.                               \tag{3}
\]
For \(t=c/L\), \(0<c\le1/2\), the error in the doubling combination
has the stronger bound
\[
 \left|D_X(t)-\big(8U_{X,\delta}(t)-U_{X,\delta}(2t)\big)\right|
 \ll \frac{c^4}{\delta}\frac{B_X}{L^3}
                       +\frac{c^2}{L^2}\sqrt X.               \tag{4}
\]
Thus for normalized critical-scale limits one may take
\(\delta=L^{-1/2}\) and restrict to subpower cofactors, with an
\(O_A(L^{-1/2})\) error on \(|u|\le A\) after normalization.
For the finer \(B_X/L^3\) scale in the doubling combination,
keep \(\delta\) fixed. At \(t=c/L\) the error in (4) is
\(O(c/\delta)\) times \(t^3 B_X\), apart from the negligible
\(\sqrt X\) term.

To prove (3), consider odd squarefree pairs on \(n\in(x/2,x]\),
with \(x\ge\sqrt X\). Put \(y=x^\eta\) for small fixed \(\eta>0\),
and let \(B_y(m)=\sum_{p<y,\ p\mid m}\log p\).
Call an endpoint good in this block if it has exactly one
prime factor at least \(y\) and the product of the other factors
is at most \(x^\delta\).
If both endpoints are good, they lie in (2). Indeed each singled-out
prime is greater than \(\frac12 x^{1-\delta}\), hence greater
than \(X^{1/3}\) for large \(X\).

For either endpoint \(m\), the general weight bound from the
uniform sieve proof is
\[
      t^{\omega(m)}
           \le t^{\omega_{<y}(m)}
                     \left(t+\frac{2 B_y(m)}{\log x}\right).
                                                               \tag{5}
\]
For a bad endpoint, make the following replacements for the
factor in parentheses, according to the reason it is bad.

- With at least two primes at least \(y\), use \(t^2\).
- With no prime at least \(y\), use
  \(C(\log x)^{-2}\sum_{p\ne q<y,\ pq\mid m}\log p\,\log q\).
  In this case \(B_y(m)=\log m\gg\log x\), and
  \(\sum_{p\mid m}(\log p)^2\le\eta\log x\log m\).
  Thus the sum over distinct primes is at least a constant
  times \((\log x)^2\).
- With one prime at least \(y\) and remaining product greater
  than \(x^\delta\), use \(t B_y(m)/(\delta\log x)\).
  Here \(B_y(m)>\delta\log x\).

For the other endpoint still use (5). Apply the same random
forbidden-class upper sieve as in (5) of the earlier note, permitting
up to two distinct prescribed primes on the bad side and one on
the other side. For prescribed squarefree products \(d,e\) of those
primes, \((d,e)>1\) is impossible for odd endpoints. When
\((d,e)=1\), the progression has modulus \(2de\) and the forced
classes give allowance factor \(t^{\omega(d)+\omega(e)}\).
The proof of that sieve bound is unchanged. Removing up to three
primes changes the Euler-product bound by a bounded factor, and
the same CRT remainder is uniform in \(d,e\).
Every term gives \(t^3\) from explicit powers or allowance factors.
The logarithmic weights sum against \(1/(de)\) by
\(\sum_{p<y}(\log p)/p\ll\log x\).
The total sieve remainder is
\(O(\delta^{-1}t^3 x^{1/4}(1+y)^3(1+\log x)^2)\),
which is smaller than the block main bound for fixed
sufficiently small \(\eta\).
The omitted mass on the block is therefore
\(O(\delta^{-1}t^3 x(\log x)^{2t-2})\).
Sum the blocks and use \(t^2\sqrt X\) for small \(n\) to get (3).

On \(n>\sqrt X\), the omitted polynomial has no quadratic term
because each twin prime pair is in (2) with \(a=b=1\).
Its cubic term disappears from the doubling combination.
By coefficient positivity its degrees four and up at \(2c/L\)
are bounded by \((2c)^4\) times the omitted polynomial at \(1/L\).
By (3) restricted to that large-\(n\) part, this is
\(O(c^4 B_X/(\delta L^3))\). Small \(n\) contribute
\(O(t^2\sqrt X)\) to the signed combination. This proves (4).
The assertion on complex compact sets uses coefficient positivity
in (3) with \(t=|u|/L\).

## What the lifted dilation theorem estimates

The fixed shift is allowed in Theorem 5.1 of the weighted-dilation
preprint. Its sum includes \(W_\ell(n)W_\ell(n+k)\), where
(5.1) assigns the falling factorial \((\omega_g(n))_\ell\)
at every prime group \(g\). Thus the support requires at least
\(\ell\) distinct primes in every group at both endpoints, for
positive fixed \(\ell\). The group count \(s\asymp\log\log x\)
tends to infinity. The theorem also requires the endpoint form
(5.3) and discrepancy condition (5.4), with its damping parameter
fixed in advance.
In (5.3) each endpoint also has a divisor in a fixed
positive-power interval \([x^\tau,x^\eta]\), \(0<\tau<\eta<1/4\),
coming from a long rough-integer factor. For the subpower
cofactors in (2) and \(n>\sqrt X\), every divisor of \(n=ap\)
either divides \(a=X^{o(1)}\) or is at least \(p\ge X^{1/2-o(1)}\).
There is no divisor in that middle interval on the corresponding
scale, even after stripping group primes.
[Weighted-dilation preprint, Theorem 3.1 and §5.1](sources/weighted-dilation-graphs-sep24-2026.pdf)

For \(t=c/L\), any unweighted restriction to pairs with total degree
at least \(M\) contributes at most
\[
            2^{-M}T_X(2t)\ll_c 2^{-M}t^2 B_X .
\]
This follows term by term. In particular the support of the
marked theorem, with \(M\ge2\ell s\) on its corresponding scale,
makes a negligible contribution to the original \(T_X(t)\) at
the scale \(t^3 B_X\) for fixed \(c>0\).
The bound is for the original weights restricted to that support.
The marked estimates have their own weights and normalization.
A transfer to the signed sum in (4) would require a new step
applicable to the unmarked or sparsely marked configurations,
including \(a=b=1\) in (2).

Bombieri's asymptotic sieve under near-level-one distribution
relates certain almost-prime factorization statistics to a single
parity parameter. The relevant factorization tests have fixed
support in the simplex or prescribed boundary behavior. The
critical-scale weights in (2) concentrate on small cofactors
and require separate estimates at the boundary. The doubling
combination selects the even degree two, kills degree three,
and subtracts the remaining degrees. A bound on the
prime-or-semiprime total leaves the signed combination undetermined.
[Bombieri, Main Theorem and §I](sources/bombieri-asymptotic-sieve-1976.pdf)
[Tao, equations (7)–(9) and Remark 4](https://terrytao.wordpress.com/2016/07/17/notes-on-the-bombieri-asymptotic-sieve/)

**Open question.** Prove a positive signed margin at the scale
\(t^3B_X\), with \(t=c/L\) for some fixed \(c>0\), for the
two-prime sum (2). In view of (4), a margin exceeding the
\(O(c/\delta)\) error in units of \(t^3B_X\) suffices.
The current lifted dilation theorem does not provide that estimate.
