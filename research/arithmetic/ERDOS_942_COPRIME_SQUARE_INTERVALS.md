# Erdős #942 and coprime occupants of square intervals

Write \(h(q)=\#\{t:q^2\le t<(q+1)^2,\ t\text{ powerful}\}\).
The left endpoint \(q^2\) is powerful, so the count in the
open interval is \(h(q)-1\).
[Erdős #942](https://www.erdosproblems.com/942), from
[Erdős, printed p. 35](https://www.renyi.hu/~p_erdos/1976-39.pdf#page=11),
asks for matching logarithmic power upper and infinitely often lower
bounds. The upper bound remains open.

The [short power window theorem CPW1](ABC_COPRIME_SHORT_POWER_WINDOWS.md#a-positive-density-primitive-construction)
gives positive density of intervals with any fixed number of pairwise
coprime nonsquare powerful occupants. Its constants depend on that number.
Selecting a larger centre for each successive fixed multiplicity
also gives qualitative unboundedness.
[Narumi–Tachiya, Theorems 1–2 and Lemmas 1–2](https://arxiv.org/html/2611.20656v1)
use simultaneous approximation to obtain ordinary occupancy of order
at least \(\log q/\log\log q\) infinitely often. Their approximation
conditions do not impose pairwise coprimality.

## Growing coprime multiplicity

**Theorem (Growing coprime occupancy).**
For each integer \(m\ge2\) there is an interval \((q^2,(q+1)^2)\)
containing at least \(\lceil m/2\rceil\) pairwise coprime nonsquare
powerful integers, with
\[
 m\le q,\qquad \log q\ll m^4
\]
with an absolute constant. Consequently along an unbounded sequence
of \(q\) the coprime occupancy is
\(\gg(\log q)^{1/4}\).
In particular the lower bound \((\log q)^{1/5}\) holds for all
sufficiently large \(q\) in that sequence.

We seek integers of the form \((B^2+d_i)^3L_i(B^2)^2\), with one
common integer square root approximation \(B^3P(B^2)\).
Truncating a binomial expansion alone leaves an error of order
\(B\) in the square root. We choose \(P\) so the first discarded
coefficient vanishes for all \(i\). This reduces the error to
\(O(1/B)\) with the construction data fixed. Resultants then let
us choose \(B\) with all the powerful values pairwise coprime.

### Cancellation at all selected shifts

Let \(c_j=\binom{-3/2}{j}\), so \(c_0=1\) and
\[
 c_j=(-1)^j(2j+1)\binom{2j}{j}/4^j,\qquad
 c_{j+1}/c_j=-(2j+3)/(2j+2).
\]
For a vector of shifts \(\mathbf d=(d_1,\ldots,d_m)\), define
\[
 R(T)=T\prod_i(T-d_i)=\sum_{j=1}^{m+1}r_jT^j,\qquad
 P(Z)=\sum_{k=0}^m p_k Z^k,\qquad p_k=c_{m+1}r_{k+1}/c_{k+1}.
\]
In particular \(P\) is monic and \(p_k\) is homogeneous of degree
\(m-k\) in the shifts. Define the monic polynomial
\[
 L_i(Z)=\sum_{k=0}^m p_k\sum_{j=0}^k c_j d_i^j Z^{k-j},\qquad
 F_i(Z)=(Z+d_i)L_i(Z).
\]
These polynomials depend on \(\mathbf d\). In the Laurent expansion
of \(P(Z)(1+d_i/Z)^{-3/2}\) at infinity, \(L_i\) is the polynomial
part. The coefficient of \(Z^{-1}\) is
\(\sum_k p_k c_{k+1}d_i^{k+1}=c_{m+1}R(d_i)=0\).

For example, with \(m=2\) and \((d_1,d_2)=(24,48)\), this yields
\[
 P=Z^2+84Z+1680,\quad L_1=Z^2+48Z-264,\quad L_2=Z^2+12Z-48 .
\]
At \(B=90725\) the two resulting powerful integers are coprime
and both exceed \((B^3P(B^2)-1)^2\) while lying below
\((B^3P(B^2))^2\). These claims can be checked by evaluating the
three integer polynomials and taking a gcd. The proof for arbitrary
\(m\) follows.

### Making the resultants nonzero

For any pair \(i\ne j\), the resultant in \(Z\) of \(F_i,F_j\)
is a nonzero polynomial in the shifts. To prove this, specialize
\(d_i=a\ne0\) and all other shifts to zero. Put
\(\lambda=c_{m+1}/c_m\). Then
\[
 P(Z)=Z^{m-1}(Z-\lambda a),\qquad
 F_j(Z)=Z^m(Z-\lambda a),
\]
so it suffices to check \(F_i\) at \(0\) and at \(\lambda a\).
The polynomial part formula gives
\[
 L_i(0)=a^m(c_m-\lambda c_{m-1}),\qquad
 L_i(\lambda a)=c_m a^m .
\]
For the second equality the two truncated sums from the two
terms of \(P\) cancel except for the term of index \(m\).
The first quantity is nonzero because \(c_m/c_{m-1}\ne\lambda\).
The factors \(Z+a\) are nonzero at both points because
\(a\ne0\) and \(\lambda\ne-1\). This proves the assertion about
the resultant.

The coefficient of \(Z^k\) in \(F_i\) has degree at most
\(m+1\) in the shifts. Its Sylvester determinant with \(F_j\)
thus has degree at most \(2(m+1)^2\). The product over pairs
has degree at most \(m(m-1)(m+1)^2\). A nonzero polynomial of
total degree \(D\) cannot vanish on a Cartesian grid with more
than \(D\) entries in each coordinate. This follows by expanding
in one coordinate, using the univariate root bound, and inducting
on the number of coordinates. Thus we can choose positive integers
\(t_i\le 1+(m+1)^4\) for which all the resultants in the vector
\(\mathbf t\) are nonzero. In particular the \(t_i\) are distinct.

### Choosing the integer parameter

Put
\[
 A=4^{2m+1}(2m+3)!,\qquad d=2m(m+1),\qquad
 M=A\prod_{p\le 2d}p,\qquad d_i=Mt_i.
\]
For the vector \(\mathbf t\), all coefficients of \(P,L_i,F_i\)
have denominator dividing \(A\). Indeed
\((2k+3)\binom{2k+2}{k+1}\) divides \((2m+3)!\) for
\(0\le k\le m\), and the powers of \(4\) from
\(c_{m+1}/c_{k+1}\) and an additional \(c_j\) divide
\(4^{2m+1}\). This also clears the denominators in the sums.
Scaling \(\mathbf t\) by \(M\) multiplies the coefficient of
\(Z^k\) by \(M^{m-k}\) in \(P,L_i\), and \(M^{m+1-k}\) in
\(F_i\). Thus the polynomials for \(\mathbf d\) have integer
coefficients, and all the coefficients below the leading one are
divisible by every prime at most \(2d\).

Use \(F_i(Z;\mathbf t)\) to denote the unscaled polynomials, and put
\[
 R_0=\prod_{i<j}
 \left|\operatorname{Res}_Z(A F_i(Z;\mathbf t), A F_j(Z;\mathbf t))\right|.
\]
This is a positive integer. For a prime \(p\) not dividing \(M R_0\),
no two of the integers \(F_i(B^2;\mathbf d)\) can be zero modulo
\(p\). Otherwise \(B^2/M\) modulo \(p\) would be a common root of
the corresponding unscaled polynomials, by homogeneity and
\(p\nmid M\). Since \(A\) is also a unit, the prime would divide
\(R_0\).

Every prime dividing \(M\) is at most \(2d\). Impose
\(B\equiv1\bmod \operatorname{rad}(M)\). Then all
\(F_i(B^2;\mathbf d)\) are one modulo the primes dividing \(M\).
For every other prime dividing \(R_0\), exclude the zeros
of the monic polynomial \(\prod_i F_i(B^2;\mathbf d)\) in \(B\).
It has degree \(d\), so at most \(d\) classes are forbidden
modulo any such prime \(p>2d\).

For the unscaled vector, each coefficient of \(R(T)\) has modulus
at most \(\prod_i(1+t_i)\), and \(1\le |c_j|\le2j+1\).
The formulas for \(P,L_i,F_i\) give
\(\log\max(2,H(A F_i(Z;\mathbf t)))=O(m\log(m+1))\),
where \(H\) is maximum absolute coefficient.
The Sylvester determinants have size \(2(m+1)\), so
\(\log R_0=O(m^4\log(m+1))\). This estimate uses the unscaled
polynomials, which keeps \(M\) out of the product of resultants.
Also \(\log M=O(m^2)\). To see the primorial bound, the product of
primes in \((l,2l]\) divides \(\binom{2l}{l}\le4^l\).
Sum the logarithmic bound over intervals between powers of two.

Write \(B=1+\operatorname{rad}(M)b\) and let \(b\) run through
\(T\) consecutive integers. Each \(p\mid R_0\), \(p\nmid M\)
excludes \(b_p\le d\) residue classes of \(b\). Put
\[
 J=1+\lceil\log_2 R_0\rceil,\qquad Y=100(d+1)J .
\]
There are at most \(J\) primes to consider. For those at most \(Y\)
put \(g_p=b_p/p\), \(S=\sum g_p\) and \(V=\prod(1-g_p)\).
We have \(g_p<1/2\), \(S=O(d\log Y)\) by the harmonic sum,
and \(V\ge e^{-2S}\). By the Chinese remainder theorem any
\(j\) prescribed prime exclusions intersect in a fraction
\(\prod g_p\) of classes modulo their product. The same
fraction holds in a progression whose step is coprime to all
those primes. In \(L\) successive terms of such a progression
the count differs from \(L\prod g_p\) by at most \(\prod b_p\).

Use inclusion–exclusion truncated at consecutive even and
odd orders of size at least \(10(S+1)\) and at most \(10(S+1)+3\).
The main term tails are at most
\(\sum_{j\ge10(S+1)}S^j/j!\le V/8\), by \(j!\ge(j/e)^j\)
and \(V\ge e^{-2S}\). The total counting error is at most
\(E=(1+dY)^{O(S+1)}\), with an absolute constant in the exponent,
by summing \(\prod b_p\) up to those orders. Thus the count
passing the primes at most \(Y\) is between \((7/8)VL-E\)
and \((9/8)VL+E\), uniformly over those progressions.
Take \(E\ge1\). Its logarithm is \(O(m^2(\log(m+1))^2)\).

For \(Y<p\le\sqrt T\), apply the upper small prime count to
each forbidden step \(p\) progression with \(L\le T/p+1\).
Summing the main terms costs at most \(2VTdJ/Y+2VdJ\),
and the errors cost at most \(EdJ\). The primes exceeding
\(\sqrt T\) cost at most \(dJ(\sqrt T+1)\), even without
a small prime restriction. Take
\[
 T=\left\lceil\big(100(d+1)(J+1)(E+1)/V\big)^2\right\rceil .
\]
The total exclusions are smaller than the lower small prime
count at \(L=T\), since \(dJ/Y\le1/100\) and
\(V\sqrt T\ge100(d+1)(J+1)(E+1)\). An admissible \(B\)
therefore occurs within \(O(\operatorname{rad}(M)T)\) beyond
any prescribed positive threshold. The logarithm of this
length is \(O(m^2(\log(m+1))^2)\).

### Passing to one open interval

Set \(D=\max_i d_i\) and let \(H\ge1\) be the maximum coefficient
modulus over \(P,L_i\) for \(\mathbf d\). The same coefficient
bounds give \(\log(2D)=O(m^2)\) and
\(\log(2H)=O(m^3)\). For \(Z=B^2\ge2D\) the binomial
series converges absolutely. The contribution to the remainder
from the \(k\)-th term after omitting its polynomial part
and the \(Z^{-1}\) term is at most
\[
 |p_k| Z^k\sum_{j\ge k+2}(2j+1)(d_i/Z)^j
 \ll (k+2)|p_k|D^{k+2}Z^{-2}.
\]
The \(Z^{-1}\) terms sum to zero by the cancellation identity.
Consequently for positive \(B\), as soon as \(L_i(B^2)>0\),
\[
 \left| B^3P(B^2)-(B^2+d_i)^{3/2}L_i(B^2)\right|
 \ll (m+2)^2 H D^{m+2}/B .
\]
Choose a common threshold for \(B\) larger than an absolute
constant times \((m+2)^2 H D^{m+2}\). This also ensures
\(d_i<2B+1\), \(B>m+2\), and positivity of \(P(B^2)\) and
all \(L_i(B^2)\), by monicity. The logarithm of the threshold
is \(O(m^3)\). Apply the preceding progression count beyond this
threshold. Its length bound gives \(\log B=O(m^3)\).

Each \(T_i=(B^2+d_i)^3L_i(B^2)^2\) is powerful. The base \(B^2+d_i\)
lies strictly between \(B^2\) and \((B+1)^2\), so \(T_i\) is
nonsquare. Its prime support equals that of \(F_i(B^2;\mathbf d)\).
Thus the \(T_i>1\) are pairwise coprime and distinct. Their
square roots lie at distance strictly less than one from the
integer \(Q=B^3P(B^2)\), and none equals \(Q\). At least half
lie on the same side of \(Q\). The two possible intervals
have lower square roots \(q=Q-1\) and \(q=Q\).
Since \(P\) is monic and the threshold makes its leading term
dominate, \(q\ge m\). Also \(\log Q=O(m\log B+\log H)=O(m^4)\).
This proves the theorem.

The argument is an elementary polynomial construction with an elementary
coprimality sieve. It gives no logarithmic upper bound on \(h(q)\).
