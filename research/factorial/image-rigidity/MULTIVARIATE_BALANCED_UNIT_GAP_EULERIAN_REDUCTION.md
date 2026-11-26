# Eulerian reduction for balanced unit-gap factorial binomials

Suppose a factorial binomial has \(r\ge2\) increasing rows of
gap one with smaller exponents \(A_1,\ldots,A_r\ge1\) and \(r\)
decreasing rows of gap one with smaller exponents
\(B_1,\ldots,B_r\ge1\). Equal rows contribute only a scalar.
The order-\(Q-1,Q\) pair has a modular coprimality test of degree
bounded independently of \(Q\) whenever \(Q=p^q\), \(q\ge1\), and
\[
                 p>\max(r,\max_i A_i+\max_i B_i).
 \tag{1}
\]
In particular, for profiles with a smaller exponent one on each
side, every such prime outside an explicit finite set depending
only on \(r\) certifies coprimality for all its powers. For \(r=2\)
the conclusion for all sufficiently large primes holds without
the smaller-exponent-one hypothesis.

This treats balanced orientations without requiring identical lists
of smaller exponents. The preceding
[linked beta reduction](LINKED_BETA_FROBENIUS_REDUCTION.md)
covers one row on each side by a different upper-order expression.

## 1. The residue-field test

Write \(A_*=\min_i A_i,\ B_*=\min_i B_i\). Define
\[
\begin{aligned}
 L_A(w)&=\sum_{k=0}^{A_*-1}(-1)^{rk}
           \prod_{i=1}^r\binom{A_i-1}{k}
           \prod_{i=1}^r\frac{B_i!\,k!}{(B_i+k)!}\,w^k,\\
 L_B(w)&=\sum_{k=0}^{B_*-1}(-1)^{rk}
           \prod_{i=1}^r\binom{B_i-1}{k}
           \prod_{i=1}^r\frac{A_i!\,k!}{(A_i+k)!}\,w^k,\\
 F(w,Y)&=w^{B_*}L_A(w)+Yw^{B_*-1}L_B(1/w),\\
 \psi&=L_A(1)+L_B(1).
\end{aligned}
\]
The denominators are units modulo \(p\), because the factors in
each quotient after cancellation have arguments at most
\(\max B_i+A_*-1\) or \(\max A_i+B_*-1\). For \(\ell\ge1\), use the
Eulerian polynomial \(R_\ell\) in the normalization
\[
                  \sum_{k\ge1} k^\ell w^k
                    =\frac{w R_\ell(w)}{(1-w)^{\ell+1}}.
\]
Put \(\ell=r-1\). The two tests
\[
                 \psi\not\equiv0\pmod p,\qquad
       \gcd_{\mathbb F_p[w]}\bigl(R_\ell(w),F(w,w^Q)\bigr)=1
 \tag{2}
\]
together certify coprimality over characteristic zero. The second
test requires just the class of \(w^Q\) in degree \(r-2\). If
\(r=2\), \(R_1=1\) and the second test is vacuous.

## 2. Reduction of the moments

Divide the factorial coefficient of \(z^k\) by the scalar
\(m!\prod_i(mA_i)!(mB_i)!\). This gives the integral coefficient
\[
 V_m(k)=(k!(m-k)!)^{r-1}
       \prod_{i=1}^r
       \binom{mA_i+k}{mA_i}\binom{m(B_i+1)-k}{mB_i}.
 \tag{3}
\]
We make further scalar normalizations in the localization
\(\mathbb Z_{(p)}\) and use \(w=(-1)^{r-1}z\).

At order \(Q\), divide by \((Q!/Q)^{r-1}\). For
\(1\le k\le Q-1\) the factorial ratio for a single factor is
\[
 \frac{k!(Q-k)!}{Q!/Q}=
        \frac{k}{\binom{Q-1}{k-1}}\equiv(-1)^{k-1}k\pmod p .
\]
The denominator is a unit by the digit formula, since every
lower digit of \(Q-1\) is \(p-1\). Both endpoint factorial ratios
are \(Q\). All row binomial factors at interior indices reduce
to one, using the digit formula with \(Q\mid mA_i,mB_i\).
Since \(r-1\ge1\), neither endpoint survives. Apart from
the scalar \((-1)^{r-1}\), the upper reduction is therefore
\(\sum_{k=1}^{Q-1} k^\ell w^k\). In characteristic \(p\) the
operator \(w\,d/dw\) annihilates \((1-w)^Q\). Applying its
\(\ell\)-th power to the geometric series
\((1-w)^Q/(1-w)\) yields the exact reduction
\[
                    w R_\ell(w)(1-w)^{Q-r}.
 \tag{4}
\]
Here we used \((1-w)^Q=1-w^Q\), and \(Q>r\).

At order \(n=Q-1\), divide by \((n!)^{r-1}Q^r\). The first
factor in (3) divided by \((n!)^{r-1}\) is a \(p\)-unit and
reduces to \((-1)^{k(r-1)}\). This is absorbed by \(w\).
At each level \(x=p^j\) for \(1\le j\le q\), an increasing
binomial contributes a carry exactly when \(k\bmod x\ge A_i\).
A decreasing binomial contributes one exactly when
\(k\bmod x<x-B_i\). This is the factorial-floor computation
from the linked beta reduction applied row by row.
Because \(x>\max A_i+\max B_i\), at least one orientation
makes all \(r\) carries. The sum equals \(r\) only in the two
bands
\[
                    k\bmod x<A_*,\qquad k\bmod x\ge x-B_* .
\]
Higher levels give zero because all factorial arguments are
less than \(pQ\). At level \(Q\) these two bands are the
two end intervals, and in those intervals the carry count is
exactly \(r\) at every smaller level as well. Thus the
normalization is \(p\)-integral and only the two end bands
survive.

At \(k=0\) each decreasing row has unit \(-1\) after division
by \(Q\). This is formula (4) of the linked beta reduction
with \(h=1,j=0\), using the decreasing base in place of \(a\).
At \(k<A_*\), that unit changes by the factor
\[
 \frac{(n)^{\underline k}}{(n(B_i+1))^{\underline k}}
       \equiv \frac{B_i!\,k!}{(B_i+k)!},
\]
where the underline denotes the falling factorial.
Each increasing row reduces to
\((-1)^k\binom{A_i-1}{k}\). The opposite end band is
the reflection exchanging \(A_i,B_i\). Since \(Q-1\) is even,
the variable change contributes the same sign to both bands.
Thus, apart from the scalar \((-1)^r\), the lower reduction is
\[
                     L_A(w)+w^{Q-1}L_B(1/w).                  \tag{5}
\]
Its leading and constant coefficients are units. In particular
zero is never a common root of (4) and (5). At \(w=1\) equation
(5) gives \(\psi\). At a different nonzero root of (4) it
vanishes precisely when \(F(w,w^Q)=0\). This proves (2).
The degree drop in the upper reduction causes no problem:
the lower degree is unchanged and its leading coefficient is a
unit, so the Sylvester resultant modulo \(p\) may be evaluated
as a product of upper values over the lower roots, with the
unit power of the lower leading coefficient taken at the
original upper degree.

## 3. Explicit families with finitely many exceptional primes

The Eulerian recurrence
\[
 R_1(w)=1,\qquad
 R_\ell(w)=(1+(\ell-1)w)R_{\ell-1}(w)
                  +w(1-w)R_{\ell-1}'(w)
\]
follows by applying \(w\,d/dw\) to the defining series. It shows
by induction that \(R_\ell\) is monic of degree \(\ell-1\), has
constant term one and has \(\ell-1\) distinct negative roots.
Indeed, at each negative zero of \(R_{\ell-1}\), the sign
of \(R_\ell\) is opposite to that of \(R_{\ell-1}'\). These
alternating signs and the signs at zero and negative infinity
give one root in each successive interval. This is the
classical [Eulerian root theorem of Frobenius](https://www.e-rara.ch/zut/content/titleinfo/5929554).

If \(A_*=B_*=1\), expression (5) is \(1+w^{Q-1}\).
A nonzero root of \(R_\ell\) in common with it satisfies
\(w^Q=-w\), hence satisfies \(R_\ell(-w)=0\) as well by
Frobenius. Over the rationals \(R_\ell(w)\) and \(R_\ell(-w)\)
have disjoint roots by the preceding sign theorem.
In this family all powers of every prime satisfying (1) and
\[
             p\nmid \operatorname{Res}(R_\ell(w),R_\ell(-w))
\]
therefore give characteristic-zero coprimality. The root \(w=1\)
is excluded because \(\psi=2\).

If \(r=2\), the coefficients of both \(L_A,L_B\) are
nonnegative, with constant coefficient one. Thus \(\psi>0\)
and (2) proves coprimality for every power of each prime
satisfying (1) and \(p\nmid\operatorname{num}(\psi)\).

More generally, if \(r\) is odd, \(L_A(w)\) and \(L_B(w)\)
are positive at every negative real \(w\). For any two roots
\(w,Y\) of \(R_\ell\) over \(\mathbb C\), the ratio \(Y/w\)
is positive. Hence
\(L_A(w)+(Y/w)L_B(1/w)>0\), and \(F(w,Y)\ne0\).
It follows that the rational number
\[
             \operatorname{Res}_w\!
                \left(R_\ell(w),
                   \operatorname{Res}_Y(R_\ell(Y),F(w,Y))\right)
\]
is nonzero. For odd \(r\) and \(\psi\ne0\), excluding the prime
divisors of the numerators of this number and \(\psi\) suffices
for all powers of the remaining primes in (1).
A failed modular test leaves characteristic-zero coprimality open.

