# Fixed-cofactor Frobenius reduction for mirrored unit-gap binomials

Consider an ordinary binomial whose unequal rows are \(r\ge2\) mirror
pairs \((a_i,a_i+1),(a_i+1,a_i)\), where \(a_i\ge1\) and \(\min_i a_i=1\).
Let \(a_*=\max_i a_i\). Equal rows give scalar factors. For a positive
integer cofactor \(m\), take \(Q=p^q\), \(q\ge1\), and \(N=mQ\), where
\[
                     p>\max\{r,2a_*,m(a_*+1)\}.
\]
The adjacent pair at terminal \(N\) has a coprimality test over
\(\mathbb F_p\) whose polynomial degrees are bounded by \(m+r-2\),
independently of \(q\). In particular, if \(r=2\) or \(3\) and
\(m=2\) or \(3\), the actual pair is coprime for all powers of every
sufficiently large prime, for every list of bases as above.
The \(m>1\) test separates the upper reduction into a cofactor part
and the Eulerian part from the
[prime-power case](MULTIVARIATE_BALANCED_UNIT_GAP_EULERIAN_REDUCTION.md).

## The residue-field test

Put \(\ell=r-1\), \(\epsilon=(-1)^\ell\), and define
\[
\begin{aligned}
 f_u&=\binom{m-1}{u}^{-\ell}
       \prod_{i=1}^r
         \binom{ma_i+u}{ma_i}\binom{m(a_i+1)-1-u}{ma_i},\\
 L(z)&=\sum_{u=0}^{m-1} f_u z^u,\qquad
 A(z)=\sum_{u=0}^{m-1}f_u
                   \prod_i\frac{ma_i}{ma_i+u}\,z^u,\qquad
 B(z)=z^{m-1}A(1/z).
\end{aligned}
\]
All the displayed rational coefficients are \(p\)-adic units.
Let \(R_\ell\) be the Eulerian polynomial in the normalization
\(\sum_{k\ge1} k^\ell x^k=xR_\ell(x)/(1-x)^{\ell+1}\).
Set \(D(z)=(1-\epsilon z)R_\ell(\epsilon z)L(z)\). If \(m\)
is even and \(A(-1)\) is a \(p\)-adic unit, remove all factors \(z+1\)
from \(D\) over \(\mathbb F_p\).
Compute \(Y=z^Q\bmod D(z)\) by modular exponentiation. The test is
\[
                  \gcd_{\mathbb F_p[z]}(D,\,z A(Y)+Y B(Y))=1 .       \tag{1}
\]
A constant \(D\) makes the test vacuous. Failed coprimality in (1)
records a common residue root only. The test has degree bound
\(\deg D\le m+r-2\) and uses only the remainder of its second input.

## Exact residue formulas

At degree \(j\), divide the ordinary coefficients by
\(j!\prod_i(ja_i)!^2\). The resulting factorial coefficient is
\[
 V_j(k)=(k!(j-k)!)^\ell
      \prod_i
        \binom{ja_i+k}{ja_i}\binom{j(a_i+1)-k}{ja_i}.
\]
For \(0\le k<N\) write \(k=uQ+v\), \(0\le v<Q\).
At order \(N\), divide \(V_N\) by \((N!/N)^\ell\).
For \(k\ge1\) the first quotient is
\((k/\binom{N-1}{k-1})^\ell\). The denominator is a unit because
the lower \(q\) digits of \(N-1\) are all \(p-1\) and \(m<p\).
The terms at \(k=0,N\) and every index \(p\mid k\) reduce to zero.
For \(p\nmid k\) the denominator reduces to
\(\binom{m-1}{u}(-1)^{v-1}\). The remaining row binomials
reduce by digit separation to the two top binomials in \(f_u\).
Thus the upper reduction is
\[
 L(z^Q)\,\epsilon\sum_{v=1}^{Q-1}v^\ell(\epsilon z)^v
                 =z L(z)^Q R_\ell(\epsilon z)(1-\epsilon z)^{Q-r}.
                                                               \tag{2}
\]
For the factorization, apply \((x\,d/dx)^\ell\) to
\((1-x^Q)/(1-x)\) in \(\mathbb F_p[x]\), then put \(x=\epsilon z\).
The numerator has zero derivative, and \(Q>r\).

For the lower order \(n=N-1\), divide \(V_n\) by
\((n!)^\ell(-mQ)^r\).
The ratio of the first factorial factor to \((n!)^\ell\)
is \(\binom nk^{-\ell}\), always a unit. For a level
\(P=p^b\) with \(1\le b\le q\), the increasing row binomial
makes a carry exactly when \(k\bmod P\ge a_i\).
Its reflected row makes a carry exactly when
\(k\bmod P<P-a_i\). The two tests follow by taking factorial
floors, since \(na_i\bmod P=P-a_i\).
The inequality \(P>2a_*\) ensures at least one carry in each pair.
At level \(Q\), the pair with \(a_i=1\) makes two carries
unless \(v=0\) or \(v=Q-1\). At those two endpoints all pairs make
exactly one carry at every level. There are no higher-level carries
because every factorial argument is less than \(pQ\).
This proves integrality and confines the lower reduction to
\(v=0,Q-1\).

At \(v=0\) the increasing unit is
\(\binom{ma_i+u-1}{ma_i-1}\) by digit separation.
For the decreasing binomial, the factorial arguments upstairs and
downstairs have high parts \(m(a_i+1)-1-u,\ ma_i-1,\ m-1-u\).
Their low parts are \(Q-a_i-1,\ Q-a_i,\ Q-1\).
The top factorial quotient is
\(ma_i\binom{m(a_i+1)-1-u}{ma_i}\).
The ratio of low-digit factorial products is \((-1)^{q-1}/a_i\):
the least digit gives \((p-a_i-1)!/((p-a_i)!(p-1)!)=1/a_i\),
and each of the other \(q-1\) positions gives \(-1\).
Stripping the valuation \(q\) supplies an additional \((-1)^q\)
from Wilson's theorem. Hence the binomial divided by \(-mQ\)
reduces to \(\binom{m(a_i+1)-1-u}{ma_i}\).
Here the digit factorial formula for the first unit is the
\(p\)-level consequence of [Granville's prime-power congruence](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node2.html).
The factor \(\binom nk^{-\ell}\) reduces to
\(\binom{m-1}{u}^{-\ell}\) in this band. Reflection supplies the
other band, so the lower reduction is exactly
\[
                         A(z^Q)+z^{Q-1}B(z^Q).                    \tag{3}
\]

The lower constant and leading coefficients are units. Every
common root in characteristic zero therefore reduces to a nonzero
common root of (2) and (3). When \(m\) is even, \(n\) is odd
and reciprocity gives the exact lower factor \(z+1\).
The derivative of (3) at \(-1\) is a unit if \(A(-1)\) is a unit,
since \(Q=0\) in the residue field and \(B(-1)=-A(-1)\).
No root of the quotient by that exact lower factor can then
reduce to \(-1\). The upper even moment is strictly positive at
the real point \(-1\), by the integral of an even power of the
difference of distinct monomials against independent exponentials.
This justifies removal of the residue \(-1\) in (1).

For even \(m\), the rational number \(A(-1)\) is strictly positive.
Indeed \(A(z)\) is a positive scalar times the degree-\((m-1)\)
moment with Gamma shapes \(a_i\) on the increasing rows and \(a_i+1\)
on the decreasing rows, as follows by expanding its factorial
coefficients. Start instead with shape \(a_i\) on both rows of
every mirror pair, writing \(X_i,Y_i\) for its independent Gamma
variables. The shape increase adds the density factor
\(\prod_i Y_i\) up to a positive scalar. The monomials satisfy
\(U-V=\prod_i(X_iY_i)^{a_i}(\prod_iY_i-\prod_iX_i)\).
Average with the simultaneous swap \(X_i\leftrightarrow Y_i\).
The odd power integral becomes a positive scalar times
\(\mathbb E[(\prod_iY_i-\prod_iX_i)(U-V)^{m-1}]\), which is positive.
Thus the simple-residue removal is available for all sufficiently
large primes in each fixed profile and cofactor.

## Cofactors two and three

For \(r=2,3\) the Eulerian factors are \(R_1=1\) and \(R_2=1+x\).
Their roots together with the factor \(1-\epsilon z\) only give
candidates \(z=\pm1\). All coefficients of \(A\) and \(B\) are
positive, so \(z=1\) is excluded outside finitely many primes.

At \(m=2\), divide \(L\) by its constant term to obtain \(1+z\).
Its only residue root and the possible point \(z=-1\) are
excluded by the even-terminal argument, after excluding the prime
divisors of the numerator of \(A(-1)\). This proves the result
for cofactor two.

At \(m=3\), divide \(L,A,B\) by \(f_0\) and put
\[
 g=\prod_i\frac{3a_i+1}{3a_i+2},\qquad
 b=\prod_i\frac{3a_i}{3a_i+1},\qquad d=bg,\qquad t=2g .
\]
Then \(L=1+t z+z^2,\ A=1+tb z+d z^2,\ B=d+tb z+z^2\).
The point \(z=-1\) gives \(A(-1)+B(-1)=2(1-d)>0\).
At a root of \(L\) over \(\overline{\mathbb F}_p\), put \(y=z^Q\).
Frobenius puts \(y\) in the same two-root set, so
\(z=y\) or \(z=1/y\). If \(z=y\), equation (3) tests
\(A(y)+B(y)=t(2b-1-d)y\). Strict convexity on the positive
axis of \(\prod_i3a_i/(3a_i+x)\) gives \(2b<1+d\).
If \(z=1/y\), multiply (3) by a nonzero power of \(y\).
The required vanishing is \(A(y)+y^2B(y)=0\), while
\[
 y^{-2}(A(y)+y^2B(y))
             =t^2(1-b)-2(1-d).
\]
This rational number is negative. When \(g\le1/2\), half its
value is \(2g^2-1+b(g-2g^2)\le g-1<0\).
For \(g>1/2\), \(b\ge2-1/g\). To verify this, each factor
of \(b\) equals \(f(x)=2-1/x\) at the corresponding factor of
\(g\), and \(f(x)f(y)-f(xy)=2(1-x)(1-y)/(xy)\ge0\).
All partial products of the factors of \(g\) exceed \(1/2\),
so induction applies. Substitution bounds half the value above
by \(-2(1-g)^2<0\).
All residue obstructions at cofactor three are therefore excluded
by finitely many explicit prime divisors of nonzero rational
numerators, uniformly in the power \(q\).

For arbitrary \(r\), the same conclusions at \(m=2,3\) hold
when \(q\) is a multiple of \(\operatorname{lcm}(1,\ldots,r-2)\),
with the empty least common multiple set to one. They also hold
for every \(q\) at admissible primes where \(R_{r-1}\) splits.
Indeed under either condition all roots of \(R_{r-1}(\epsilon z)\)
satisfy \(z^Q=z\), by their residue field degrees. On these
candidates equation (3) reduces to \(A(z)+B(z)=0\).
For \(m=2\) this polynomial is a nonzero rational multiple of
\(1+z\), and the possible root \(-1\) was already excluded.
For \(m=3\), the coefficients of \(A+B\) are proportional to
\((1+d,2tb,1+d)\). Since \(2b<1+d\) and \(0<t<2\), its
discriminant is negative. The Eulerian polynomial has only real
negative roots, by the [Eulerian recurrence and sign argument](MULTIVARIATE_BALANCED_UNIT_GAP_EULERIAN_REDUCTION.md#3-explicit-families-with-finitely-many-exceptional-primes).
Thus a nonzero rational resultant excludes finitely many primes
in this remaining modular test.

For instance, with \(r=2\) and both \(a_i=1\), cofactor two
works for \(p>5\). At cofactor three the two nontrivial normalized
tests are \(-188y/625\) and \(-352y^2/625\), so \(p>6\),
\(p\notin\{11,47\}\), suffices.
