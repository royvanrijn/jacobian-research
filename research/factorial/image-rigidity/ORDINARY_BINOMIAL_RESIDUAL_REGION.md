# Ordinary binomial Strong Factorial

The [original polynomial triangle](POLYNOMIAL_TRIANGLE_AND_TERMINAL_SIX.md)
defines the actual moment polynomials in (T1) and the gap cutoff in (T2).
Work throughout with ordinary factorial moments and their rational polynomials.

**Theorem.** For every characteristic-zero field, every nonzero ordinary
binomial \(f=X^a(A+B X^d)\) with \(a\ge0,\ d\ge1\), and every integer
\(N\ge2\), at least one of the factorial moments
\(\mathcal L(f^{N-1}),\mathcal L(f^N)\) is nonzero.

For \(A B\ne0\), the moment at order \(m\) equals
\(A^m(am)!\,P_{m,a,d}(B/A)\). Rational coprimality of
\(P_{N-1,a,d},P_{N,a,d}\) therefore gives the assertion by Bézout
over any characteristic-zero field. For a nonzero monomial,
every positive factorial moment is nonzero.
The [earlier terminal theorem](ORDINARY_BINOMIAL_TERMINALS_THROUGH_TWENTY_FOUR.md)
covers \(2\le N\le24\).
The [finite certificate](#the-finite-certificate)
covers \(25\le N<3150\), and the
[Padé cutoff](#coprimality-from-terminal-3150) covers \(N\ge3150\).

## Residual region

Use the explicit common cutoff \(N_*=3150\) from
[the sharper Padé comparison](#coprimality-from-terminal-3150).
The [earlier prime-gap cutoff](#a-lower-terminal-cutoff-from-verified-prime-gaps) covers \(N\ge20{,}000\).
The [analytic cutoff](#an-explicit-terminal-cutoff) already covers
\(N\ge\lceil e^{15}\rceil=3{,}269{,}018\). There is also an independent
[support-uniform asymptotic proof](SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SYNTHESIS.md#3-a-finite-covering-of-every-positive-lower-exponent-support)
of the existence of a common cutoff. For each \(N\), keep \(D_N\)
from (T2). Below the cutoff it suffices to prove rational coprimality at integer triples
\((a,d,N)\) satisfying
\[
\begin{gathered}
25\le N<N_*,\qquad 2\le d<\min\{D_N,128(N-1)\},\qquad
5\le a,\qquad d<3a,\\
d^2<128(a+d),\qquad
a<\min\{d(N-2),\,10^{24}d^2,\,10^{20}d^4\},\\
(a,d)\notin\{(p^w,p^s):p\text{ prime},\ 0\le w\le s,\ s\ge1\},\\
(a,d)\notin\{(2p^w,p^s):p\text{ prime},\ 1\le w<s\}.
\end{gathered}
\tag{T5}
\]
The endpoint coefficient comparison below imposes the additional bound
\[
 s=\frac{d^2}{a+d}<\sigma_N,\qquad
 \sigma_N=
 \begin{cases}
  2\log N+3,&25\le N<200,\\
  2\log N+2,&N\ge200.
 \end{cases}
 \tag{T5a}
\]
In particular \(d<\sigma_N(N-1)\) in (T5).
The prime-gap bounds (T33) and (T40) give \(d\le140\) throughout
this remaining terminal range.
The sharper Padé comparison requires a full Kummer coset for
\(N\ge275\) and excludes \(\alpha=d/(a+d)<.03,\ a+d\le N\)
for \(N\ge200\). The [near-terminal prime tests](ORDINARY_BINOMIAL_COMMON_SLOPE_BOUND.md#prime-factors-near-the-terminal-order) also constrain the candidate degrees. For \(N\ge2200\) only
\(\alpha<.03,\ a+d>N\) remains.

The [strengthened Gamma inner annulus](ORDINARY_BINOMIAL_STRICT_GAMMA_INNER_ANNULUS.md#the-product-bound-and-its-ordinary-consequences)
also imposes \(6a<(5N-11)d\) at these orders. Its nonlinear
condition (5) gives a further exact modulus comparison.
At prime gap \(d=p\), impose the two further necessary conditions
\[
 \gamma=a(N-1)+1,\quad C=\gamma+pN,\qquad
 v_p\!\left(\frac{pN}{(C)_a}\right)\le-2,\qquad
 [Y^2](\gamma+Y)_a\equiv0\pmod p.
 \tag{T6}
\]
The two polynomials are the actual pair \(P_{N-1,a,d},P_{N,a,d}\)
from (T1). Equivalence holds for the conjecture over \(\mathbb C\) and
for the conjecture over all characteristic-zero fields.

The [terminal-through-twenty-four theorem](ORDINARY_BINOMIAL_TERMINALS_THROUGH_TWENTY_FOUR.md)
covers \(2\le N\le24\). The common eventual cutoff covers \(N\ge N_*\).
The triangle note's Section 1 gives \(d<D_N\), and
[sharp shape separation](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#3-shape-separation-and-a-sharper-soft-complement)
excludes \(a\ge d(N-2)\). The
[quadratic](UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md) and
[quartic](UNIFORM_QUARTIC_STRONG_FACTORIAL_TAIL.md) all-terminal tails
give the other two upper bounds on \(a\).
The [effective small-ratio theorem](EFFECTIVE_SMALL_RATIO_ORDINARY_BINOMIAL_SFC.md)
excludes \(d^2\ge128(a+d)\) at every order. Its remaining
inequality gives \(a>d^2/128-d\). Together with \(a<d(N-2)\)
it gives \(d<128(N-1)\). Thus for fixed \(N\) there are
\(O(N^3)\) candidate support pairs, using at most \(d(N-2)\)
lower exponents per gap.

The [constant-plus-power](ALL_UNIVARIATE_TWO_TERM_SFC.tex),
[adjacent](UNIFORM_PADE_PARABOLA_SFC.tex),
[lower-one](ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md),
[lower-two](ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md),
[lower-three](STEPANOV_SFC_PROGRESS_AND_LOWER_THREE_CLOSURE.md) and
[lower-four](#lower-exponent-four) results give \(d\ge2\) and \(a\ge5\).
The [relative-gap theorem](RELATIVE_GAP_THREE_STRONG_FACTORIAL.md)
covers \(d\ge3a\) for \(a\ge2\).
The [same-prime support theorem](SFC_PRIME_POWER_DIVIDED_DIFFERENCE_CLOSURE.md#4-every-lower-power-of-the-same-prime-is-closed)
covers \(a=p^w,d=p^s\) for \(0\le w\le s\).
The [upper-response argument at offset \(2p^w\)](PRIME_POWER_HASSE_COPRIMALITY.md#offset-twice-a-prime-power)
covers \(a=2p^w,d=p^s\) for \(1\le w<s\).
The [arbitrary-offset native-prime criterion](NATIVE_PRIME_ARBITRARY_OFFSET_COPRIMALITY.md)
gives (T6), including the exclusion of \(a=2p,d=p\) at every order.

Each restriction excludes a proved domain. Coprimality on (T5)–(T6)
and on the complement gives the full result by (T1) and rational Bézout.
Conversely, a nonconstant rational gcd on this region gives a common
complex root, hence a complex counterexample. The finite certificate below proves coprimality throughout this region.

### The finite certificate

A possible common factor for \(25\le N<3150\) must have
\(5\le a\), \(2\le d<3a\), \(a<d(N-2)\), and satisfy
\[
 6a<(5N-11)d,\qquad d^2<19(a+d).
\]
The first inequality is the strengthened Gamma annulus.
The second follows from (T5a), since \(\sigma_N<19\) throughout this
range. By (T33) and (T40) the gap caps are \(35,56,76,110,140\)
on terminals \(25\)–\(50\), \(51\)–\(100\), \(101\)–\(300\),
\(301\)–\(1000\), and \(1001\)–\(3149\), respectively.
Below \(275\) the caps use only the first part of (T33).
For \(2200\le N<3150\), the Padé exclusions also require
\(b>N\) and \(97d<3a\).
These inequalities give a containing region, so the
check can omit the other support exclusions in (T5)–(T6).

For each \(N,a\) the bounds \(6a<(5N-11)d\), \(d<3a\),
the gap cap, and, when \(N\ge2200\), \(a+d>N\) and
\(97d<3a\), give an integer interval \([d_0,d_1]\).
The lower endpoint is at least two. For any candidate
degree \(h\) and either \(m=N-1,N\), (T7) makes the
first \(dh\) factorial arguments \(am+1,\ldots,am+dh\)
free of prime factors above \(dN\). If \(N\ge275\) or \(a\ge2d\),
(T10) is available and lowers that cap to
\[
 \min\{dN,\ \max(a+2d-1,\lfloor\sqrt{(a+d)N}\rfloor)\}.
\]
Use \(dN\) when these sufficient conditions for (T10) are
unavailable. This choice of cap is nondecreasing in \(d\),
including the switch to \(dN\) for \(a<2d\) at a small terminal.
In either order, a prime above the cap at \(d_1\) in the
first \(2d_0\) arguments therefore excludes the whole interval.
The verifier then checks the remaining \(d\) against the
individual cap in the first \(2d\) arguments, after the two
inequalities \(a<d(N-2)\) and \(d^2<19(a+d)\).

For each support left by these tests, try the even degrees
\(2\le h<N-1\). Check the starting-point inequality for (T10)
directly, or use the uniform extension for \(N\ge275\).
If (T10) is available, impose its degree bound. At each
remaining degree impose the full (T7) prime bound on the
first \(dh\) arguments in both orders.

The Newton test is the degree allocation (T12) at odd primes
\(aN<p\le(a+d)N\).
For these primes \(p>N\), \(p>d\), and \(p^2>(a+d)N\)
because \(d<3a\). Thus all binomial coefficients are units,
\(am<p\), and the ordinate at index \(k\) is
\(\lfloor(am+dk)/p\rfloor\).
Writing \(q_m=\lfloor(a+d)m/p\rfloor\le3\), the lower
hull uses at most five points
\[
 (0,0),\qquad
 \left(\left\lfloor\frac{jp-am-1}{d}\right\rfloor,j-1\right)
 \quad(1\le j\le q_m),\qquad (m,q_m).
\]
The intermediate points listed are the right ends of the
constant-ordinate intervals. The lower hull of these points is
nondecreasing. Every omitted point is to the left of a retained
point of the same ordinate, so it is on or above this hull.
For each common slope with reduced denominator \(v\)
admit the multiples of \(v\) up to the smaller of its
two widths, including zero. The set of attainable total
degrees is the sumset over slopes. Intersect this set
with the candidate degrees after each prime, stopping
when it is empty.

The [exact verifier](verify_ordinary_binomial_residual.cpp)
computes an odd smallest-prime-factor table through
\(1{,}156{,}520{,}634\), the largest factorial argument
allowed by its bounds. Segmented Eratosthenes marking
fills the table and prime divisors are obtained by division.
All prime, cap and hull comparisons are integer comparisons.
The square root used in a prime cap is corrected by integer squaring.
The first prime-cap check leaves \(21{,}583\) groups \((N,a)\).
After the individual \(d\) checks, \(74{,}362\) supports remain.
The degree and prime-block tests exclude \(62{,}593\).
The Newton allocations exclude \(11{,}705\), leaving only \(64\)
supports, all of gap two. For those the verifier forms the two
monic reciprocals modulo \(1{,}000{,}000{,}007\) and computes their
gcd by the Euclidean algorithm. The modulus is checked by trial
division and the coefficients are computed by the successive
binomial and factorial ratios. All gcds are one. Any nonconstant
common factor over \(\mathbb Q\) would be monic integral by
Gauss's lemma and would remain a nonconstant common divisor
after reduction, so the entire finite region is coprime.

From the repository root, replay with a C++20 compiler and
OpenMP. The prime-factor table takes about \(1.1\) GiB.

```sh
g++ -std=c++20 -O3 -fopenmp research/factorial/image-rigidity/verify_ordinary_binomial_residual.cpp -o /tmp/verify-ordinary-binomial-residual
OMP_NUM_THREADS=12 /tmp/verify-ordinary-binomial-residual
```



### A smaller high-gap complement

Use the endpoint normalization of \(P_{m,a,d}\) at \(m=N-1,N\)
from [the center and coefficient estimates](lower_exponent_initial_and_boundary.tex),
Lemma \(\mathrm{lem:centers}\) and equation (7).
With \(G_m=(((a+d)m)!/(am)!)^{1/m}\), that normalization is
\[
 P_{m,a,d}(w/G_m)=1+w^m+\sum_{k=1}^{m-1}c_{m,k}w^k,\qquad
 c_{m,k}\le\binom{m}{k}\exp\!\left(-\frac{s k(m-k)}{2m}\right).
\]
The majorant is symmetric in \(k,m-k\). Suppose \(s\ge\sigma_N\),
and write \(C=3\) or \(2\) for the constant in \(\sigma_N\).
For \(1\le k\le4\) the majorant is at most \((AE^k)^k/k!\),
where \(A\ge e^{-C/2}\) and
\(E\ge\exp((\log N+C/2)/(N-1))\). Indeed \(m\le N\) and
\(\binom{m}{k}\le m^k/k!\).
For \(5\le k\le m/2\), write \(u=k/m\).
By \(\binom{m}{k}\le(e/u)^k\), the logarithm of the majorant
divided by \(k\) is at most \(1-\log u-\sigma_N(1-u)/2\).
This is convex in \(u\), so the majorant is at most \(q^k\) if
\[
 q\ge e\max\{AE^5/5,\ 2e^{-\sigma_N/4}\}.
\]
The constants in the table suffice. The upper bounds for \(E\)
and \(2e^{-\sigma_N/4}\) are checked at the smallest \(N\) of
each row with its own \(C\), since they decrease with \(N\).
\[
 \begin{array}{c|c|c|c|c}
   N&A&E&q&\eta\\\hline
   25\le N<200&0.224&1.218&0.52&7/8\\
   N\ge200&0.368&1.033&0.24&19/20
 \end{array}
\]
In each row \(2\sum_{k=1}^4(AE^k)^k/k!+2q^5/(1-q)<\eta\)
by rational arithmetic. This bounds \(\sum_k c_{m,k}\) by \(\eta\).

Comparison with the constant or leading term puts the normalized
root moduli strictly between \((1-\eta)^{1/m}\) and
\((1-\eta)^{-1/m}\). The center ratio satisfies
\(G_N/G_{N-1}>(N/(N-1))^{d/2}\) by the cited lemma.
Since \(d>s\ge\sigma_N\), we have \(d\ge10\) in the first row
and \(d\ge13\) in the second. In both rows,
\[
 \frac{d}{2}\log\!\left(\frac{N}{N-1}\right)
 \frac{N(N-1)}{2N-1}
 \ \ge\frac{d(N-1)}{2(2N-1)}\ >-\log(1-\eta).
\]
This separates the two annuli after undoing the normalization and
proves (T5a). It uses only the two actual ordinary shapes.

## Arithmetic restrictions on every common factor

Let \(b=a+d\), \(M=bN\), and use the monic reciprocals
\(R_m(x)=x^mP_{m,a,d}(1/x)\), \(m=N-1,N\). Any monic rational common
factor \(F\) lies in \(\mathbb Z[x]\) by Gauss's lemma, with nonzero
constant term. Every nonconstant such factor has even degree \(h\).
Indeed the even-order reciprocal is a positive Gamma integral at
every real \(x\), so \(F\) has no real zero. The
[consecutive nondivisibility theorem](BINOMIAL_IRREDUCIBILITY_REDUCTION.tex)
also gives \(h<N-1\).

### Large primes and short products

For every triple satisfying (T5) and every nonconstant common factor
of degree \(h\), set \(H=dh\) and
\[
 B_h=\min\{dN,\ \max(2H,\ a+2d-1,\ \sqrt{M})\}.
\]
Then the integer
\[
 F(0)\,(a(N-1)+1)_H\,(aN+1)_H
 \tag{T7}
\]
has all its prime divisors at most \(B_h\). The bound applies at both
orders and for both prime and composite gaps. Real square-root
comparisons in this criterion can be implemented by integer squaring.

Write \(y_m(k)=v_p([t^k]P_{m,a,d})\). The slopes of its lower Newton
polygon are the valuations of the reciprocal roots, with
multiplicities given by edge lengths.

First suppose \(p>dN\). The bound \(a<d(N-2)\) gives \(p^2>M\),
and \(p>N\) makes all binomial coefficients units. Each interval
\((am,bm]\) contains at most one multiple \(u_m\) of \(p\).
If it has one, its positive valuation slope is \(1/w_m\), where
\(w_m=1+\lfloor(bm-u_m)/d\rfloor\). Indeed the coefficient valuations
are zero until the single jump to one. Equal positive slopes at both
orders would imply equal \(w_m\), hence
\[
 u_N-u_{N-1}\in[b-d+1,b+d-1]=[a+1,a+2d-1].
 \tag{T8}
\]
This is a positive interval below \(p\) and contains no multiple of \(p\).
Thus all common reciprocal roots are \(p\)-adic units.

Next suppose \(p>\max(2H,a+2d-1,\sqrt M)\). Here \(p\) need not
exceed \(N\). Because \(p^2>M\ge N\), coefficient valuations are
\[
 y_m(k)=\lfloor(am+dk)/p\rfloor-\lfloor am/p\rfloor
            +\epsilon_m(k),\qquad
 \epsilon_m(k)=\mathbf 1_{\{k\bmod p>m\bmod p\}}.
 \tag{T9}
\]
The binomial term follows by the single-floor formula for
\(v_p(m!),v_p(k!),v_p((m-k)!)\) and is zero at \(k=m\).

A positive valuation \(s\) of a root of \(F\) has reduced
denominator at most \(h\). The corresponding edge of the Newton
polygon of \(F\) has integral height, with positive width at most
\(h\), so \(s\ge1/h\). If an edge of slope \(s\) for \(P_m\) begins
at \(k\), convexity and (T9) give
\[
 s(m-k)\le y_m(m)-y_m(k)<1+\frac{d(m-k)}p
                         <1+\frac{m-k}{2h}.
\]
Thus \(m-k<2h\). The middle difference is a positive integer
less than two, so it is one. This edge is the last edge: any
intermediate vertex would have integer ordinate strictly between
\(y_m(k)\) and \(y_m(m)\), as every subsequent slope is positive.
Consequently \(s=1/\ell\), \(\ell=m-k\le h\).
The factorial floor difference from \(k\) to \(m\) is at most one.
Hence \(\epsilon_m(k)=0\) by (T9). Since \(\ell<p\), this means
\(k\) is in the same block of \(p\) consecutive integers as \(m\),
and \(\epsilon_m(k+1)=0\). The point at \(k+1\) lies on or above
the edge, so \(y_m(k+1)\ge y_m(k)+1\). Therefore a multiple \(u_m\)
of \(p\) lies in \((bm-d\ell,bm-d(\ell-1)]\). The common slope
uses the same \(\ell\) at both orders and gives (T8), again
impossible. Every prime \(p>B_h\) satisfies at least one of the
two sufficient bounds above, so \(F(0)\) is a \(p\)-adic unit.

If such a prime divided \((am+1)_H\) at either order, the
rising-factorial coefficients at \(k\ge h\) would all vanish modulo
\(p\). Then \(R_m\bmod p\) would be \(x^{m-h+1}\) times a polynomial
of degree \(h-1\). Its monic divisor \(F\bmod p\) has degree \(h\)
and nonzero constant term, a contradiction. This proves (T7).
In particular, a prime \(p>B_2\) in either first block of \(2d\)
integers excludes a degree-two factor. A prime \(p>dN\) in either
such block excludes every common factor by even parity, since
\(B_h\le dN\) and \(dh\ge2d\).

The first block in the upper order also gives a degree bound without
factoring that block. Its logarithm exceeds \(H\log(aN)\).
Its prime support is at most \(B_h\), and the valuation at \(p\)
is at most \(H/(p-1)+\lfloor\log_p M\rfloor\) by counting
multiples of each power of \(p\) in the block. Thus
\[
 H\log(aN) <
 H\sum_{p\le B_h}\frac{\log p}{p-1}
   +\sum_{p\le B_h}\lfloor\log_p M\rfloor\log p .
 \tag{T7a}
\]
In particular, one can use the bound with any larger prime threshold.

A uniform consequence uses Laishram and Shorey's Theorem 1(a).
A block of \(H\ge2\) integers beginning at
\(x>\max\{H+13,279H/262\}\) has a prime divisor greater than \(2H\).
Apply this at \(x=aN+1\). If the starting-point inequality holds,
a candidate degree \(h\) must satisfy
\[
 2dh<\max\{a+2d-1,\sqrt{(a+d)N}\}.
 \tag{T10}
\]
For example the starting-point inequality holds throughout the part
of (T5) with \(a\ge2d\), since \(N\ge25\), \(d\ge2\), and \(h<N-1\).
[Laishram and Shorey, *The greatest prime divisor of a product of consecutive integers*, Theorem 1(a)](https://www.impan.pl/shop/en/publication/transaction/download/product/82021).

In particular, a whole part of the containing region is coprime:
\[
 \frac{279}{262}d\le a\le2d+1,\qquad (a+d)N\le16d^2.
\]
Here \(a\ge d+1\), so \(aN+1-dh>13\) at every possible degree
\(h<N-1\), and \(aN+1>279dh/262\). The right side of (T10) is at
most \(4d\), which excludes every even \(h\ge2\).

### A uniform upper bound on the remaining degrees

If \(N\ge20{,}000\), the hypotheses of (T5) alone force any common
factor to satisfy (T10), even for \(a<2d\). In fact, two primes just
above \(X=aN\) first give
\[
 h<\frac{3aN}{500d}.
 \tag{T11}
\]
It follows that \(aN+1>\max(dh+13,279dh/262)\), as required for
(T10). Thus at these terminals the only candidate degrees are
below \(\max(a+2d-1,\sqrt{(a+d)N})/(2d)\).

Put \(\eta=1/1000\). Use the three windows
\[
 X(1+\eta,1+2\eta],\quad X(1+5\eta,1+6\eta],\quad
 X(1+13\eta,1+14\eta].
\]
Each contains a prime by Dusart's Proposition 5.4, which gives
a prime in \((x,x+x/\log^3 x]\) for \(x\ge89{,}693\).
Here the left endpoints satisfy \(x\ge X\ge10^5\), so
\(\log x>11\) and \(x/\log^3 x<\eta X\).
[Dusart, *Explicit estimates of some functions over primes*, Proposition 5.4, printed p. 242](https://doi.org/10.1007/s11139-016-9839-4).

Discard any window meeting an interval
\[
 [(M-b)/j,(M+2d)/j],\qquad 1\le j\le4.
\]
At most one of these intervals can meet the envelope of all three
windows. Their centers \(M/j\) have spacing at least \(X/12\), their
endpoints lie within \(6a\) of their centers, and the envelope has
width \(13\eta X\), with \(13\eta X+12a<X/12\). Each forbidden
interval has width at most \(b+2d<10a\), smaller than the minimum
window separation \(3\eta X\). Thus there are primes \(p_1<p_2\) in
two distinct windows avoiding all four forbidden intervals.

For either prime \(p\), we have \(p>X\ge m\), \(p>d\), \(p^2>M\),
and \(am<p\) at both orders. The entire Newton ordinate is
\(\lfloor(am+dk)/p\rfloor\). Both orders have the same
\(q=\lfloor bm/p\rfloor\le3\), and \((q+1)p-bm>2d\), by the
forbidden intervals. If \(q=0\), all slopes are zero and
\(h\le N<(p-X)/d\) from \(p>M\).

If \(q\ge1\), let \(k_j(m)=\lfloor(jp-am-1)/d\rfloor\) for
\(1\le j\le q\), the last index before jump \(j\). The slope-zero
width is \(k_1(m)\). Every positive edge of the lower hull on
\([0,k_q(m)]\) goes between two of the \(k_j(m)\). At either
vertex of a positive interior edge, a further point on the same
plateau would lie below the supporting line. Consequently its slope
is at most \(1/\lfloor p/d\rfloor\), by the separation of the
\(k_j(m)\). The final edge from \(k_q(m)\) to \(m\) has height
one and width
\[
 w_m=\left\lceil (bm-qp+1)/d\right\rceil
     \le\lceil p/d\rceil-2<\lfloor p/d\rfloor.
\]
Its slope is strictly greater than the interior slopes, so
appending it preserves the hull, with the intervening ordinates
all equal to \(q\) and lying above the appended segment. Since \(b>d\), we have
\(w_N>w_{N-1}\), and neither final slope is shared.

For an interior edge of the upper hull, write \(J_e\) for its
positive integer height and \(w_e=J_e p/d+e\) for its width,
with \(|e|<1\) by the formula for \(k_j(N)\).
If \(F\) uses \(r_e\) roots from that edge, then
\(0\le r_e\le w_e\) and \(r_eJ_e/w_e\) is integral, by the
integral ordinates of its own Newton polygon. Summing the widths
used gives
\[
 h=h_0+Jp/d+E,\qquad 0\le h_0\le k_1(N)<(p-X)/d,\qquad
 J=\sum_e r_eJ_e/w_e\in\{0,1,2\}.
\]
The interior edges have total height \(q-1\le2\), so \(|E|<2\)
when \(J>0\). When \(J=0\), every \(r_e=0\) and \(E=0\).
Thus for \(J>0\) a prime imposes an interval from \(Jp/d-2\) to
\(Jp/d+2+(p-X)/d\). At \(p_1,p_2\) these intervals with the
same positive \(J\) are disjoint since \(p_2-p_1>p_1-X+4d\).
This follows from the window endpoints and \(\eta X>4d\).
The intervals for different \(J\), including \(J=0\), are
also disjoint across the primes, since \(p\in(X,(1+14\eta)X]\)
and \((1-28\eta)X>4d\). Therefore both local conditions must
have \(J=0\). The smaller usable prime lies in one of the first
two windows, and hence \(h<(p_1-X)/d\le6\eta X/d\), proving (T11).

### All-prime degree and constant-term tests

The Newton common-degree argument in
[Residual lifting and parity for binomial factorial moments](BINOMIAL_RESIDUAL_LIFTING.tex)
applies to the actual coefficients
\(c_{m,k}=\binom mk(am+1)_{dk}\) with arbitrary \(a,d\). Here it
can also constrain the constant term of the monic reciprocal factor.
For each prime \(p\le B_h\) and each slope \(\lambda=u/v\) in
lowest terms, let \(w_m(\lambda)\) be its edge width for
\(P_{m,a,d}\), using zero if the slope is absent. For this prime,
a common factor of degree \(h\) requires an integer allocation
\[
 \sum_\lambda v j_\lambda=h,\qquad
 0\le j_\lambda\le
 \frac{\min(w_{N-1}(\lambda),w_N(\lambda))}{v}.
 \tag{T12}
\]
Only slopes present in both polygons are needed, with \(v=1\) at
zero. The corresponding exponent of \(p\) in \(F(0)\) is
\(\sum_\lambda u j_\lambda\). This is because the root valuations
of a factor form a submultiset at each prime and each edge width is
a multiple of its reduced slope denominator. Multiplying its
reciprocal roots gives its constant term up to sign.

One may strengthen (T12) by requiring \(j_\lambda\) to be the
degree of a divisor of the gcd of the two edge residuals in
\(\mathbb F_p[Y]\). For an edge starting at index \(i\), its residual
uses the coefficients at \(i+vj\), divided by
\(p^{y_m(i)+uj}\) before reduction, as coefficients of \(Y^j\).
To see the multiplicativity, adjoin \(\pi\) with \(\pi^v=p\),
substitute \(t=\pi^{-u}z\), and divide each polynomial by the power
of \(\pi\) given by its minimum coefficient valuation. The reductions
are nonzero and multiply. Each is a monomial \(z^i\) times its edge
residual at \(Y=z^v\), so each common factor contributes a residual
dividing both, up to a nonzero scalar. For the factor of the
reciprocals, apply this argument to its own reciprocal in \(t\).

The [sharp Gamma-product annulus](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#2-every-positive-shape-a-sharper-product-annulus)
gives the simultaneous Archimedean requirement
\[
 \bigl((aN+1)_d\bigr)^h
       < |F(0)|
       < \bigl((a(N-1)+1+d(N-2))_d\bigr)^h .
 \tag{T13}
\]
For the lower bound use the upper-order roots. For the upper bound
use the lower-order roots. Both degrees are at least two. Thus the
products of the primes at their respective minimum and maximum
attainable exponents from (T12) bound \(|F(0)|\) on each side, with
integral comparisons against (T13). If an allocation is empty or
the resulting interval is disjoint from (T13), the degree is
excluded. Primes above \(B_h\) require no allocation for this step
by (T7).

The [strengthened inner annulus](ORDINARY_BINOMIAL_STRICT_GAMMA_INNER_ANNULUS.md)
can replace the lower endpoint in (T13) by its formula (7).
The [joint valuation bounds](ORDINARY_BINOMIAL_COMMON_SLOPE_BOUND.md)
also constrain the slopes and prime support in (T12). They use
both orders' terminal intervals and apply at primes with higher
factorial powers.

A [Chinese remainder obstruction](BINOMIAL_NEWTON_PRIME_SET_OBSTRUCTION.md)
shows why a fixed finite list of primes cannot certify the slope-width
test at all orders in the unbounded containing region. The
obstruction concerns slope and width alone.

### An exponential Padé bound on small common factors

The constant-term test has a second constraint from a nonzero integral
resultant. This adapts the Laguerre Padé norm argument to all ordinary
gaps by using all \(d\) rotations of the exponential. The integral
normalization used below keeps the \(p\)-adic loss proportional to the
moment order. [Nguyen-Dang, *Uniform irreducibility and Galois groups of generalized Laguerre polynomials*, v1, Section 3.2.1](https://hexagonmath.org/pdf/2611.00087v1).

Let \(F\) be an irreducible monic common factor of degree \(h\), put
\(m=N-1\), \(r=(b-1)m\), \(c=(d-1)/2\), \(s=d^2/b\) and
\(\alpha=d/b\). Put \(K_\pm=1/2\pm\alpha/14\),
\(r_0=b(1-\alpha K_+)=a+dK_-\), and define
\[
 g_+(v)=b\log(b/v)-b+v,\qquad
 g_-=r_0-a-a\log(r_0/a),\qquad g=\max\{g_+(r_0),g_-\}.
\]
\[
\begin{aligned}
 L_0&=\frac{m\sqrt a}{d\bigl(1+c/(a+c)+(a+2c)^{-1/2}\bigr)}
       \exp\!\left(-\frac{(1+c/\sqrt a)^2}{2}\right),\\
 L&=\max\{L_0,(m!e^{-mg-1})^{1/(m+1)}\}.
\end{aligned}
\]
For some \(k\in\{1,\ldots,d\}\), there is a monic factor of
\(F(-z^d)\) of degree \(D=hk\) whose absolute constant term \(K\)
satisfies \(K^d=|F(0)|^k\), hence \(K\mid |F(0)|\), and
\[
\begin{aligned}
 L_k&=\left(\prod_{j=0}^{k-1}
          \max\left\{L,\frac{m}{d}
                        \left(\frac{2aj^2}{d^2}-d\right)\right\}\right)^{1/k},
 \qquad V=\log(M+1),\\
 \log L_k
 &<
 \sum_{p\mid K}\min\left\{\frac{v_p(K)}D,\frac1{p-1}\right\}\log p
          +\frac{\omega(K)V}{N}.
\end{aligned}
\tag{T14}
\]
Here \(\omega(K)\) counts the distinct prime divisors of \(K\).
Either (T12)–(T13) or a prime bound for the larger norm \(|F(0)|\)
can be used to bound the right side.

**Construction.** Set \(\widehat P(z)=(am)!R_m(-z^d)/m!\),
and \(Q(z)=[e^z\widehat P(z)]_{\le r}\), where the brackets
truncate in \(z\). For \(0\le j\le bm\) write \(s=bm-j\).
Expansion of the factorial quotient gives
\[
 [z^j]e^z\widehat P(z)
   = \frac{s!}{m!}[w^s](1+w)^{am}((1+w)^d-1)^m.
 \tag{T15}
\]
Indeed the term of \(R_m\) with index \(i\) contributes
\((-1)^{m-i}\binom mi (am+di)!/(m!(am+di-s)!)\),
interpreted as zero if \(s>am+di\). The binomial expansion gives
(T15), including these zero terms. In particular \(Q\in\mathbb Z[z]\)
has degree \(r\), leading coefficient \(d^m\), and the terms of
degrees \(r+1,\ldots,bm\) vanish. If \(\omega_j\) ranges over the
\(d\) distinct \(d\)-th roots of unity, then
\(e^{-\omega_j z}Q(\omega_j z)-\widehat P(z)\) has order at least
\(bm+1\) at zero.

Use the Wronskian of \(\widehat P\) and the \(d\) functions
\(e^{-\omega_j z}Q(\omega_j z)\), with derivative orders \(0\) through
\(d\) as rows. After removing the exponential factors it is a
polynomial of degree at most \(dm+dr=dbm\). Its top coefficient
is nonzero. It uses row zero of the \(\widehat P\) column and, up
to nonzero column factors, the determinant with entries
\((-\omega_j)^i\) for \(1\le i\le d\). This is a Vandermonde times
\(\prod_j(-\omega_j)\). To see its order at zero, subtract the
\(\widehat P\) column from the other columns and take constant
linear combinations of the resulting \(d\) independent series so
that their leading orders are distinct. These are at least
\(bm+1,\ldots,bm+d\). The determinant vanishes to order at least
\(d\,bm\), after subtracting the sum of the row derivative orders.
The polynomial determinant is therefore a nonzero multiple of
\(z^{dbm}\). At any nonzero root of \(\widehat P\) at least one
\(Q(\omega_j z)\) is nonzero, by its zeroth row. Choose a rotated
root of \(F(-z^d)\) where \(Q\ne0\) and let \(G\) be its monic
minimal polynomial over \(\mathbb Q\). Then \(G\in\mathbb Z[z]\)
and \(\operatorname{Res}(G,Q)\) is a nonzero integer. The tower
through \(-z^d\) gives \(D=hk\), with exactly \(k\) roots of \(G\)
above each root of \(F\). Multiplying these roots in absolute
value gives \(K^d=|F(0)|^k\).

**Local bound.** Write \(c_s=[x^s](1+x+\cdots+x^{d-1})^m\).
For \(\ell\ge1\), the exponential tail coefficient is
\[
 [z^{bm+\ell}]e^z\widehat P(z)
  =(-1)^m\binom{m+\ell-1}{m}
       \sum_s c_s\frac{(am+s)!}{(am+s+m+\ell)!}.
 \tag{T16}
\]
This follows by replacing the factorial ratio
\((am+di)!/(am+di+\ell)!\) by
\(\int_0^1 t^{am+di}(1-t)^{\ell-1}\,dt/(\ell-1)!\)
and summing the binomial coefficients. Factor
\(1-t^d=(1-t)(1+t+\cdots+t^{d-1})\) and evaluate the remaining
beta integrals.

Fix a prime \(p\) and a root \(z\) of \(G\) with \(e=v_p(z)>0\).
Put \(\mu=1/(p-1)\). Legendre's factorial formula and
\(s_p(n)/(p-1)\le\log_p(n+1)\), where \(s_p\) denotes the digit sum,
give
\[
 v_p(Q(z))\ \ge\ (bm+1)e-(m+1)\min\{e,\mu\}-\log_p(M+1).
 \tag{T17}
\]
For the digit inequality, each digit \(u\) satisfies
\(p^{u/(p-1)}\le u+1\) by concavity of \(\log\).
The product of \(u+1\) over the digits counts integers whose
digits do not exceed those of \(n\), so it is at most \(n+1\).
If \(e\le\mu\), the coefficient of \(z^{r-j}\) in (T15)
is divisible by \((m+j)!/m!\), of valuation at least
\(j\mu-\log_p(M+1)\). Each term of \(Q(z)\) therefore has valuation
at least \(re-\log_p(M+1)\), giving (T17).
If \(e>\mu\), the \(p\)-adic exponential converges at \(z\)
and \(Q(z)\) equals the negative tail, since \(\widehat P(z)=0\).
The factorial denominator in (T16) has valuation at most
\((m+\ell)\mu+\log_p(M+1)\). Indeed its starting argument
\(am+s\) is at most \(r\), and Legendre's formula bounds the
digit error by the digit sum at this starting argument, uniformly
in \(\ell\). The valuation lower bound for the terms at \(z\)
increases with \(\ell\), so the first tail index proves (T17).
For \(e=0\) integrality gives \(v_p(Q(z))\ge0\).

A second formula for the same tail makes its local losses
depend on arithmetic progressions. For \(\ell\ge1\),
\[
 [z^{bm+\ell}]e^z\widehat P(z)=
 \frac{(-1)^m d^m}{(\ell-1)!}
 \sum_{j=1}^\ell
 \frac{(-1)^{j-1}\binom{\ell-1}{j-1}}
      {\prod_{i=0}^{m}(am+di+j)} .
 \tag{T17a}
\]
To obtain it, expand \(1/(am+di+1)_\ell\) into simple fractions
and use
\(\sum_{i=0}^m(-1)^i\binom mi/(x+di)=d^m m!/\prod_{i=0}^m(x+di)\),
which follows by comparing residues at all the simple poles.
In particular the first tail coefficient has absolute value
\(d^m/\prod_{i=0}^m(am+di+1)\).

For \(e>\mu\), (T17a) and the vanishing of \(\widehat P(z)\)
give the additional lower bound
\[
 v_p(Q(z))\ \ge\ (bm+1)e+m v_p(d)-
 \sup_{j\ge1}\left\{
 v_p\!\left(\prod_{i=0}^m(am+di+j)\right)-
                       (j-1)(e-\mu)\right\}.
 \tag{T17b}
\]
Indeed \(v_p((\ell-1)!)\le(\ell-1)\mu\), and each index
\(j\le\ell\) benefits from \((\ell-1)(e-\mu)\ge(j-1)(e-\mu)\).
This supremum is a maximum over a finite range: the valuation
of the product grows at most as \((m+1)\log_p(bm+j)\).
For \(p\nmid d\), its exact evaluation amounts to counting
indices \(i\in[0,m]\) in the single residue class
\(-d^{-1}(am+j)\bmod p^w\) for each power \(p^w\).
One can use the larger of (T17) and (T17b) in the local norm.

**Complex bound.** Write \(R=|z|\); the root annulus gives
\(am<R\le bm\). Equations (T15)–(T16) give the two estimates
\[
\begin{aligned}
 |Q(z)|&\le \frac{d^m R^r}{m!}
           \int_0^\infty t^m e^{-t}(1+t/R)^r\,dt,\\
 |Q(z)|&\le \frac{d^m R^{bm+1}}{m!}
           \int_0^1 s^m(1-s)^{am}e^{s\Re z}\,ds.
\end{aligned}
\tag{T18}
\]
The first follows from the nonnegative coefficients of the
polynomial in (T15) and integration of each factorial. The second
sums the entire exponential tail in its beta integral before using
\(1-(1-s)^d\le ds\); it applies at a root of \(\widehat P\).
In particular if \(am-\Re z>0\), the second estimate gives
\(|Q(z)|< R^{bm+1}/(((am-\Re z)/d)^{m+1})\).

For \(R\le(a+c)m\), use
\(\log(1-s)\le-s-s^2/2\) in the second estimate.
The square bound \(s^2/2\ge us-u^2/2\) with
\(u=(c+\sqrt a)/a\) leaves a decaying exponential of rate
\(m\sqrt a\), at cost \(\exp(m(1+c/\sqrt a)^2/2)\).
For \(R\ge(a+c)m\), put
\(v=(r/R-1)_++(a+2c)^{-1/2}\) in the tangent bound
\(\log(1+t/R)\le\varphi(v)+(t/R)/(1+v)\), where
\(\varphi(v)=\log(1+v)-v/(1+v)\le v^2/2\).
Then \(R-r/(1+v)\ge m\sqrt a/(1+v)\), using
\(a(a+2c)\le(a+c)^2\). Also
\((a+2c)v^2\le(1+c/\sqrt a)^2\) and
\(1+v\le1+c/(a+c)+(a+2c)^{-1/2}\).
Integrating the resulting gamma bound in (T18) in both cases
gives \(|Q(z)|<R^{bm+1}/L_0^{m+1}\), where strictness also follows
from \(d^m<d^{m+1}\).

To obtain the other gain in \(L\), keep \(((1+u)^d-1)^m\)
in the first integral of (T18) and \((1-(1-u)^d)^m\) in the second,
before the estimates by \(d^m\). More precisely,
\[
 \frac{|Q(z)|}{R^{bm+1}}\le\frac{1}{m!}\min\left\{
 \int_0^\infty e^{-Ru}(1+u)^{am}((1+u)^d-1)^m\,du,\
 \int_0^1 e^{u\Re z}(1-u)^{am}(1-(1-u)^d)^m\,du
 \right\}.
 \tag{T18a}
\]
The first bound follows directly from (T15) by integrating each
factorial after scaling the variable by \(R\). The second follows
from the beta integral in (T16).

If \(R\le mr_0\), discard the last factor in the second integral
and maximize its exponent to bound the integral by \(e^{mg_-}\).
If \(R\ge mr_0\), bound \((1+u)^d-1\) by \((1+u)^d\) in the
first integral and reserve \(e^{-u}\). The remaining exponent is
at most \(m g_+(r_0-1/m)\).
Since \(s=\alpha d<9a/4\),
\[
 r_0-b/2=a/2-s/14>19a/56>1.
\]
Thus \(|g_+'(v)|<1\) on \([r_0-1/m,r_0]\), and integrating
\(e^{-u}\) gives the bound \(e^{mg_+(r_0)+1}\).

Both normalized losses are increasing in \(\alpha\) for
\(0<\alpha\le3/4\), as follows from
\[
 \frac{g_+(r_0)}s=K_+^2\int_0^1\frac{u\,du}{1-\alpha K_+u},
 \qquad
 \frac{g_-}s=K_-^2\int_0^1\frac{u\,du}{1-\alpha+\alpha K_-u}.
\]
The first integrand increases directly. For the second, the
derivative of \(K_-^2/(1-\alpha+\alpha K_-u)\) has the sign of
\(K_-(1-u/2)-(1-\alpha)/7\), which is positive. The strict
integral bounds above give \(|Q(z)|/R^{bm+1}<e^{mg+1}/m!\).
In particular, if \(U=\log N\ge23\), then
\[
 \log L>U-1.02-g.
\]
Indeed \(m!\ge(m/e)^m\), \(0.01(m+1)>\log m\) and
\(\log m>U-0.01\).

Among \(k\) distinct rotations above one root of \(F\), the
\(j\)-th smallest absolute argument, indexed from \(j=0\), is at
least \(\pi j/d\), by the spacing \(2\pi/d\). For such a rotation
with argument \(\theta\),
\(am-\Re z\ge am(1-\cos\theta)-dm\ge 2amj^2/d^2-dm\).
Therefore the product of the complex gains over all roots of \(G\)
is at least \(L_k^D\). Summing (T17) over roots with positive
valuation and primes dividing \(K\), and comparing with the
complex bound on the same integral resultant, proves (T14).

### Large primes and the degree of the Padé lift

Use \(U=\log N\), \(B=H^4/M\), and \(P=2H+d\).
Every prime dividing either initial block \((am+1)_H\) for
\(m=N-1,N\) divides \(F(0)\). Indeed, reduction at that prime
leaves \(R_m\) equal to a power of \(x\) times a polynomial of
degree less than \(h\), as in (T7). It follows that the source
of the prime in the proof of (T10) also constrains the lift.

For any \(p\), the sum of the largest \(h\) valuations of
roots of \(R_m\) is at most \(H/(p-1)+\lfloor\log_p M\rfloor\).
To see this, use the upper concave hull of the backward
ordinates
\[
 U_m(j)=v_p(c_{m,m})-v_p(c_{m,m-j})
       =v_p\left(\prod_{i=bm-dj+1}^{bm}i\right)-v_p\binom mj .
\]
Counting multiples of each prime power in the product gives
\(U_m(j)\le dj/(p-1)+\lfloor\log_p M\rfloor\), which also
bounds the hull at \(j=h\). Thus for \(p>2H\),
\[
 v_p(F(0))\le\lfloor\log_p M\rfloor .
\]
If the starting-point inequality of (T10) holds, Laishram
and Shorey give such a \(p\) in the upper initial block.
For any irreducible lift of degree \(hk\), the identity
\(K^d=|F(0)|^k\) then gives
\[
 \ell=\left\lfloor\frac{\log M}{\log(2H+1)}\right\rfloor,\qquad
 d\mid k v_p(F(0)),\qquad k\ge d/\ell .
 \tag{T19}
\]
All primes dividing \(F(0)\) divide \(K\) as well. In
particular a prime dividing \(F(0)\) to exponent one forces
\(k=d\). This makes \(F(-z^d)\) irreducible up to its leading
sign and gives all \(d\) rotations in (T14). A
[local ramification test](ORDINARY_BINOMIAL_COMMON_SLOPE_BOUND.md#a-prime-support-test-for-lifts-of-smaller-degree)
also limits the prime support of a lift of degree below \(dh\).

**Rotations in a proper lift.** Fix a root \(\xi\) of \(F\)
and set \(E=\mathbb Q(\xi)\) and \(E_{\rm cyc}=E(\mu_d)\), where
\(\mu_d\) is the group of \(d\)-th roots of unity. For \(z^d=-\xi\)
let \(k_0=[E_{\rm cyc}(z):E_{\rm cyc}]\). Over \(E_{\rm cyc}\)
the extension is Galois and the ratio \(\sigma(z)/z\) identifies
its group with \(\mu_{k_0}\subseteq\mu_d\). Thus the \(k\) lifts
above \(\xi\) belonging to \(G\) form a union of cosets of
\(\mu_{k_0}\) acting by multiplication. All these cosets have
the same size because all the lifts generate the same extension
of \(E_{\rm cyc}\).

For the prime \(p>2H>d\) supplied in (T19), some prime of \(E\)
over \(p\) has integer-normalized valuation \(u\) of \(\xi\) with
\(1\le u\le v_p(F(0))\le\ell\), by the norm formula. Adjoining
\(\mu_d\) is unramified at \(p\): the \(d\)-th roots of unity in a
finite extension of the residue field lift by Hensel to the
corresponding unramified extension of the completion.
Since \(z^{k_0}\in E_{\rm cyc}\), its \(d/k_0\)-th power equals
\(-\xi\), and \(d/k_0\) divides \(u\). Consequently
\[
 k_0\ge d/\ell.
 \tag{T19c}
\]
Every coset consists of evenly spaced points on a circle. It
supplies at least \(\lfloor k_0/2\rfloor\) rotations with
\(\Re z\le0\), at each of which (T18) gives gain at least
\(mA\), where \(A=a/d\). The same \(k_0\) applies above each
conjugate of \(\xi\). These gains and the baseline \(L\) can
be used together in (T14).

If \(q=d/k_0\ge2\), then the prime support is smaller.
At a prime \(p>H+d\), the \(q\)-th power relation in \(E_{\rm cyc}\)
and unramified base change make every positive valuation of
\(\xi\) in its completion over \(\mathbb Q_p\), with normalization
\(v_p(p)=1\), at least \(q/h\). Indeed its integer-normalized
valuation in that completion is a multiple of \(q\) and the
ramification index is at most \(h\). The norm valuation at such
\(p\) is also a multiple of \(q\).

Suppose \(p>H+d\) and \(p^{q+1}>M\). The affine bound on \(U_m(j)\)
with error \(q\) puts the maximizing indices for a selected positive
slope below \(qh/(q-1)\le2h\). In this range \(j!\) is a unit and
the factorial jumps are spaced at least \(h+1\) apart with height
at most \(q\). Thus the only possible maximizers are zero and the
first jump, with a unique such slope at each order. The norm
valuation is at most \(H/(p-1)+q<q+1\). If it is positive, it is
therefore \(q\), and the edge at each order must have height \(q\).
The common slope gives equal widths, so the difference between
the two reached multiples of \(p^q\) belongs to \([a+1,a+2d-1]\),
as in (T8). Consequently the norm and the two initial blocks
have prime support at most
\(\max\{H+d,M^{1/(q+1)},(a+2d-1)^{1/q}\}\).

Under the starting-point hypothesis of (T10), the first-block
prime \(p>2H\) forces \(M^{1/(q+1)}>2H\).
Indeed \(H+d<2H\). If \((a+2d-1)^{1/q}>2H\), then
\(s<(1+\alpha)/(4h^2)\). Using \(b<dN\), we have
\[
 \frac{M}{(a+2d-1)^{3/2}}>
 \frac{1}{\sqrt s(1+\alpha)^{3/2}}>
 \frac{2h}{(1+\alpha)^2}>1.
\]
Thus \(M^{1/(q+1)}>(a+2d-1)^{1/q}\)
when the latter could exceed \(2H\). In all cases with \(q\ge2\),
\[
 2H<M^{1/(q+1)},\qquad
 \text{each prime in the norm or either initial block is at most }
 M^{1/(q+1)}.
 \tag{T19d}
\]

**Short products for unequal edges.** Suppose \(B\ge10\). For \(p>P\) we have \(p^4>M\).
A positive valuation \(\lambda\) of a root of \(F\) is at
least \(1/h\), and at each order \(U_m(j)-\lambda j\) has at
least two maximizing indices. Its maximum is at least zero.
By the preceding affine bound all maximizing indices are
less than \(6h\). Call \(p\) exceptional if either terminal
interval \([bm-6H,bm]\) contains a positive multiple of \(p^3\).
There are at most \(2M/P^3\) exceptional primes. For each
order and each positive multiplier \(q\le M/P^3\),
the interval can contain at most one \(qp^3\) with \(p>P\),
since successive cubes at that scale differ by more than
\(3P^2>6H\). Each exceptional prime contributes at most
\(\log M\) to \(\log |F(0)|\) by the integer valuation bound.

At a nonexceptional prime the sharper estimate
\(U_m(j)\le dj/(p-1)+2\) on \(j<6h\) puts all maximizing
indices below \(4h<p\). In this range the valuation of the
binomial is nondecreasing since \(j!\) is a unit. Between
factorial jumps \(U_m(j)-\lambda j\) therefore decreases
strictly. Consecutive factorial jumps are separated by at
least \(\lfloor p/d\rfloor\ge2h+1\) indices, and each jump has
height at most two. The values of \(U_m(j)-\lambda j\) at
successive jump indices decrease strictly as well.
Consequently the maxima are at \(0\) and the first jump,
at index \(\ell_m\), where the factorial interval first
contains a multiple \(u_m\) of \(p\). The slope is
\(\lambda=e_m/\ell_m\), with height \(e_m\in\{1,2\}\).
This is the unique possible shared positive slope of
denominator at most \(h\). Its contribution to \(v_p(F(0))\)
is at most \(\min(e_{N-1},e_N)\).

Write \(u_m=bm-d\ell_m+q_m\) with \(1\le q_m\le d\).
If the two heights agree, their widths agree and
\(u_N-u_{N-1}\in[a+1,a+2d-1]\). This difference is divisible
by \(p^{v_p(F(0))}\) since each \(u_m\) is. As \(p>2d-1\),
the product of these prime powers divides
\(\binom{a+2d-1}{2d-1}\). Its log is less than
\(2d(1+\log(a/d+2)-\log2)\), by \(\binom{n}{k}\le(en/k)^k\)
and monotonicity of \(k(1+\log(n/k))\) up to \(k=n\).

If the heights are unequal, then \(v_p(F(0))=1\) and the
width with height two is twice the other width. The integer
\(2u_{N-1}-u_N\), or \(2u_N-u_{N-1}\) according to the
order of the heights, is a multiple of \(p\) in one of
\[
 b(N-2)+[2-d,2d-1],\qquad b(N+1)+[2-d,2d-1].
 \tag{T19a}
\]
The prime \(p^2\) divides the terminal integer at the
height-two edge. The reduced slope denominator is at most
\(h\), so this edge has width at most \(2h\) and puts the
square multiple in \([bm-2H,bm]\). In particular
\(p\le\sqrt M\).

Write \(S\) for the sum of \(\log p\) over these unequal
edges for primes dividing \(F(0)\). The two intervals
(T19a) have fewer than \(3d\) integers each, all positive
and less than \(2M\). There is also at most one prime \(p>P\)
per multiplier \(q\le M/P^2\) in a terminal interval of width
\(2H\), since consecutive odd squares differ there by more
than \(4P\). The bounds are
\[
 S\le\min\{6d\log(2M),\ (M/P^2)\log M\}.
 \tag{T19b}
\]
The first bound can be refined by taking the product of
the two binomial coefficients of width \(3d-2\) whose
numerator intervals are (T19a), since \(p>3d-2\).

**Consequences for common factors.** The prime estimates used
here and below are
\[
 \pi(x)\le \frac{x}{\log x}\left(1+\frac{1.2762}{\log x}\right)
 \quad(x>1),\qquad
 \sum_{p\le x}\frac{\log p}{p-1}<\log x-0.54
 \quad(\log x\ge8).
\]
Dusart's Theorem 6.11 bounds the sum with denominator \(p\)
by \(\log x+E+0.2/\log x+0.2/\log^2x\) for \(x>2974\).
Here \(E=-\gamma_{\rm E}-\sum_p\log p/[p(p-1)]\) and
\(\gamma_{\rm E}>0.57\). Adding the prime-power correction
gives the second estimate. The first is Theorem 6.9.
[Dusart, *Estimates of some functions over primes without R.H.*, v1](https://arxiv.org/abs/1002.0442v1).

For \(B\ge10,\log P\ge8\), put \(A=a/d\) and
\[
 \Delta=\log(aN/H)-\frac{2}{h}(1+\log(A+2)).
\]
The strict norm bound (T13), the affine estimate at
\(p\le P\), and the three contributions at \(p>P\) above give
\[
 \Delta<
 \log(P/H)+\frac{\pi(P)\log M}{H}
 +\frac{S}{H}+\frac{\log M}{4B}.
 \tag{T20}
\]
Here the last term bounds the charge from the cubes.
The stronger equal-height bound also gives (T20) with
\(\Delta_*=\Delta+(2/h)\log2\) on the left.
Since \(P/H\le2.5\), \(\log P>8\) and \(\log M<4\log P\),
the first two terms total less than \(13\).
Also \(\log M<3U+5\) by (T5).

When (T10) applies, we have \(\Delta>U/2-1/2\).
For \(h=2\), use \(A+2<7A\) to obtain \(\Delta>U-4\).
For \(h\ge4\), if \(H<(a+2d-1)/2\le7a/2\) in (T10),
then \(A>2h/7>1\). Using \(A+2<3A\) gives
\[
 \Delta>U-(2/h)\log h-(1-2/h)\log(7/2)
                 -(2/h)(1+\log3)>U-5.
\]
In the other case, \(H<\sqrt M/2\) implies
\(aN/H>\sqrt{aN}\) by \(b<4a\). Then \(A+2<a\) and
\(h\ge4\) give \(\Delta>U/2-1/2\).
We also have the second estimate
\(\Delta>U-\log h-3\) directly from \(A>1/3\) and \(A+2<7A\).

If \(U\ge34\), it follows that \(S/H>U/2-13.5-(3U+5)/(4B)>0\).
Thus **some large prime has norm exponent one**, so \(k=d\)
in this region. Also (T19b) first gives \(h<800\):
bound \(h\) by \(6(3U+6)/(U/2-13.5-(3U+5)/40)\),
which decreases for \(U\ge34\). Then use
\(\Delta>U-\log h-3\), substituting \(\log 800\) and
then \(\log 76\), to obtain successively \(h<76\) and \(h<60\).
Using this in (T20) with the second bound in (T19b) gives
\[
 h\le58,\qquad
 \frac{H^3}{M}<
 \frac{3U+5}{4(U-16-\log60-(3U+5)/40)}
 <\frac52 .
 \tag{T21}
\]
The last expression before \(5/2\) decreases for \(U\ge34\).
In particular, since \(b<d(N-1)\), we obtain
\(d/N<s\) and \(H/N<\min(sh,5/(2sh^2))\).

**Counting the unequal edges.** Suppose \(U\ge15\) and \(H<b\)
with \(b>N\). For either order \(m\) put \(X=bm\) and count
primes \(p>P\) with \(X-2H\le qp^2\le X\) for an integer \(q>0\).
Work on intervals \([D,\eta D)\) up to \(\sqrt X\), where
\(\eta=9/8\) and \(D=\eta^j P\). These lattice points are within
vertical error \(\epsilon=2H/D^2\) below \(f(x)=X/x^2\).

Each prime permits at most one multiplier because \(p^2>2H\).
Each multiplier permits at most one prime in the terminal interval,
because successive odd squares above \(P\) differ by more than \(4P\).
Thus for distinct points
\[
 1\le |q_2-q_1|\le 2X|p_2-p_1|/D^3+\epsilon .
\]
Since \(\epsilon<1/(2H)\le1/8\), the prime abscissas have
spacing at least \(g(D)=\max(2,7D^3/(16X))\).

Partition the interval into blocks of width \(\rho=(D^4/X)^{1/3}\),
clipping the last block if needed. Their number is at most
\(0.132(X/D)^{1/3}\) since \(\eta-1=0.125\) and
\(\rho/D\le X^{-1/6}<0.007\).
For three points in one block the integer area determinant
is the determinant for \(f\) plus that of the vertical errors.
The first has magnitude at most \((3X/D^4)\rho^3/4=3/4\),
by the second divided difference and \(f''(x)\le6X/D^4\).
The error determinant is at most \(\epsilon\rho\le(2H/X)^{1/3}<0.01\),
using \(H<b\), \(X\ge b(N-1)\), \(U\ge15\).
Thus all lattice points of a block are collinear.

For a block of \(n\ge3\) points, the chord of \(f\) through
the extreme abscissas lies at most \(\epsilon\) above \(f\)
at a middle abscissa. Both distances from the middle point
to the extremes are at least \((n-2)g(D)/2\).
The lower curvature \(6X/(\eta D)^4\) therefore gives
\[
 n\le2+2\eta^2\sqrt{2/3}\frac{D\sqrt{H/X}}{g(D)}.
\]
The same bound holds when \(n\le2\). The first terms from
the blocks sum to less than \(7(X/P)^{1/3}\), since
\(0.264/(1-\eta^{-1/3})<7\).
For the second terms, split the geometric sums at
\(D_*=(32X/7)^{1/3}\), where \(g(D)\) switches.
Their total is less than \(3.4X^{1/18}\sqrt H\), because
\[
 0.132\,\eta^2\sqrt{2/3}\,(32/7)^{2/9}
 \left(\frac1{1-\eta^{-2/3}}+\frac1{1-\eta^{-7/3}}\right)<3.4 .
\]
Both orders satisfy \(X\le M\) and each prime has
\(\log p\le\tfrac12\log M\). Hence
\[
 \frac{S}{H\log M}<
 \frac{7}{(2B)^{1/3}}+\frac{3.4M^{1/18}}{\sqrt H}.
 \tag{T22}
\]

For a sharper log charge, weight each strip by \(\log(\eta D)\).
This replaces the first bound by
\(7(X/P)^{1/3}(\log P+3.06)\) because
\(\log\eta/(1-\eta^{-1/3})<3.06\).
In the second bound the terms with \(D\le D_*\) have
weight at most \(\log(\eta D_*)\). Above \(D_*\), the
function \(D^{-7/3}\log(\eta D)\) is decreasing. The geometric
sum there is bounded by the one starting at \(D_*\), with
mean log weight \(\log D_*+\log\eta/(1-\eta^{-7/3})
<(\log X)/3+1\). Both parts therefore give a combined
weighted charge below \(3.4X^{1/18}\sqrt H((\log X)/3+1)\).
Summing the two orders gives
\[
 \frac{S}{H}<
 \frac{14(\log P+3.06)}{(P/H)^{1/3}B^{1/3}}
 +\frac{6.8 M^{1/18}}{\sqrt H}\bigl((\log M)/3+1\bigr).
 \tag{T22a}
\]

### Full Kummer cosets for \(\log N\ge15\)

Take \(U=\log N\ge15\), \(m=N-1\), and assume the support bounds
in (T5) and (T5a) apart from the upper bound on \(N\).
In particular \(s<\sigma=2U+2\), \(N\ge20{,}000\), and (T10)
applies. We prove that every possible irreducible common factor
requires \(k_0=k=d\) in its Padé lift.

If the subgroup of order \(k_0\) is proper, put \(q=d/k_0\ge2\)
and \(P=M^{1/(q+1)}\). By (T19d), \(2H<P\) and all norm
primes are at most \(P\). Also \(P^{q+1}=M<dN^2=HN^2/h\), so
\[
 \frac PN<(2h)^{-1/q}e^{-(1-2/q)U}.
\]
Dusart's prime count preceding (T20), \(P>2H\ge8\), and
\(V/\log M<1.001\) give
\[
 \frac{\omega(K)V}{N}
 <1.616(q+1)4^{-1/q}e^{-(1-2/q)U}.
\]
For \(q\ge3\) this decreases with \(q\) on the required range.
In fact \(q<\log M/\log8<2U-\log4\), and the logarithmic
derivative in \(q\) is \(1/(q+1)-(2U-\log4)/q^2<0\).
At \(q=3,U=15\) the bound is below \(0.03\).
For \(q=2\) it is below \(2.43\). In the latter case
\(\omega(K)V/N<4.848P/N\), and we can use
\[
 P/N=(b/N^2)^{1/3}<(d/N)^{1/3}.
\]
Thus the digit loss is below \(0.121\) if \(d<50\), or below
\(0.01\) if \(\alpha\ge0.03\) by \(b\le \sigma/(0.03)^2\).

Use the prime sum preceding (T20) with \(x=\max(e^8,P)\).
For the geometric mean \(\mathcal L\) of gains from (T14), this
gives
\[
 \log\mathcal L<\max(8,\log P)-0.54+\frac{\omega(K)V}{N}.
 \tag{T23}
\]
In (T18a) the baseline gain \(m e^{-1.01-g}\) applies for
\(U\ge15\), since \(0.01(m+1)>\log m\) and \(\log m>U-0.001\).

**Moderate and large \(\alpha\).** For \(k_0\ge2\) at least a
third of every coset has \(\Re z\le0\). There (T18) gives
gain at least \(m(1-\alpha)/\alpha\).
On a row \([\alpha_-,\alpha_+]\) of the table we use
\(g\le\kappa s\) by the formulas in (T18a).
It follows that
\[
 \log\mathcal L>
 U-0.001-\tfrac23(1.01+\kappa\sigma)
          +\tfrac13\log\bigl((1-\alpha_+)/\alpha_+\bigr).
\]
The upper bound in (T23) is at most
\(\max(8,(U+\log(\sigma/\alpha_-^2))/3)-0.51\).
The right-hand columns compare these expressions at \(U=15\).
\[
 \begin{array}{c|c|c|c}
 [\alpha_-,\alpha_+]&\kappa&\text{lower bound}&\text{upper bound}\\\hline
 [0.03,0.1]&0.134&11.9&8.1\\
 [0.1,0.5]&0.176&10.5&7.5\\
 [0.5,0.75]&7/32&9.2&7.5
 \end{array}
\]
The lower expressions increase more quickly with \(U\), so
these cosets are impossible.

**Small \(\alpha\).** Put \(t=\log(1/\alpha)\), so \(3.5<t<U\),
and \(g_0=1.01+0.128s\). The baseline gain is \(m e^{-g_0}\).
If \(s\le1\), we may replace \(g_0\) by its minimum with
\(4+\tfrac12\log s\), since (T14) gives \(L_0>m e^{-4}/\sqrt s\).
In that formula \(d/\sqrt a<2\sqrt s\), \(c/\sqrt a<1\), and
the parenthesized denominator is less than \(3\).
A gain uses all rotations in each coset. Index them by increasing absolute argument, so the angle in position \(j\) is at least \(\pi j/k_0\). Put \(J=\lceil k_0\sqrt\alpha\rceil\). For \(j\ge J\), the second integral in \(\mathrm{(T18)}\) gives individual gain at least
\[
 m(1-\alpha)\alpha^{-1}(j/k_0)^2 .
\]
For \(j\le k_0/2\), use \(1-\cos(\pi j/k_0)\ge2(j/k_0)^2\) and \(\Re z\le bm\cos(\pi j/k_0)\), which holds since the threshold cosine is nonnegative, irrespective of the actual angle. For \(j>k_0/2\), use \(\Re z\le0\). The other positions have the bare gain. Since \(\sum_{j=1}^{k_0-1}\log(j/k_0)\ge-k_0\),
\[
 \log\mathcal L
 \ge \log m+t+\log(1-\alpha)-2
                -(\sqrt\alpha+1/k_0)(t+g_0).                  \tag{T24}
\]
Suppose \(k_0\ge2\) and \(q=2\). For \(d\ge50\),
\(g_0/k_0\le 2.02/d+0.256\alpha<0.049\), and
\[
 \sqrt\alpha(t+g_0)\le e^{-1.75}(4.766+0.256U)
\]
by monotonicity in \(t\ge3.5\). Thus (T24) gives
\(\log\mathcal L>0.955U+0.96t-2.92\).
This exceeds \(\max(8,(U+2t+\log\sigma)/3)+1.89\) from
(T23) at \(U=15,t=3.5\), with both differences increasing
in \(U,t\). For \(d<50\), the baseline exceeds \(U-1.22\)
in logarithm since \(s<1.5\). Use \(M<50N^2\) and digit
loss below \(0.121\) in (T23).

For \(q\ge3,k_0\ge2\), use the left third with gain
\(m(1-\alpha)/\alpha\) and the baseline \(m e^{-g_0}\)
on the rest. This gives \(\log\mathcal L>0.829U+t/3-0.86\).
The upper bound is \(\max(8,(U+2t+\log\sigma)/4)-0.51\);
subtracting it gives a positive difference using \(t<U\)
and \(U\ge15\).

**One rotation.** If \(k_0=1\), then \(q=d\) and (T19d) gives
\(d\log(4d)<2U+\log d\), by \(M<dN^2\).
This forces \(d<U\), since the inequality fails at \(d=U\)
and the left side minus \(\log d\) increases for larger \(d\).
If \(d\ge3\), the baseline with \(s<0.75U\) and \(g\le7s/32\)
gives \(\log\mathcal L>0.835U-1.011\). This beats
\(\max(8,(2U+\log U)/4)-0.51\) for \(U\ge15\).
If \(d=2\), use \(s<1.5\) and \(M<2N^2\). The baseline
\(\log\mathcal L>U-1.34\) beats
\(\max(8,(2U+\log2)/3)+1.89\).
This completes the full-coset proof.

### Support restrictions below the terminal cutoff

Continue with \(U\ge15\), \(\sigma=2U+2\), and full cosets.
First suppose \(\alpha\ge0.03\).
Then \(b\le\sigma/(0.03)^2<N/4\), so (T10) and (T7)
give prime support at most \(\sqrt M\).
Use \(V/\log\sqrt M<2.001\) and \(\log\sqrt M>U/2\)
in the prime count. The right side of (T14) is less than
\[
 \tfrac12(U+\log b)-0.54+2.35\sqrt{b/N}.
\]
Here \(\log\sqrt M>8\), since \(\alpha<0.75,\ a\ge5,\ d\ge2\).
If \(d\ge10\), the left-half-plane fraction is at least
\(0.45\). In the first and third rows above, the log gains
are greater than
\[
 U-0.001-0.55(1.01+\kappa\sigma)
            +0.45\log\bigl((1-\alpha_+)/\alpha_+\bigr).
\]
In the middle row, the fraction with \(|\arg z|\ge\pi/3\)
is at least \(2/3-1/d\ge17/30\). Those rotations have
\(\Re z\le am\); the second integral in (T18a) is at most
one, giving gain at least \(m e^{-1.01}\). The log gain
is greater than \(U-1.011-(13/30)\kappa\sigma\).
At \(U=15\) the three lower expressions exceed
\(13.07,\ 11.54,\ 10.09\), respectively.
Using \(b\le\sigma/\alpha_-^2\), the upper expressions are
less than \(12.45,\ 11.08,\ 9.41\), respectively.
The lower expressions increase more quickly. If \(d\le9\),
use \(b\le300\), \(s<6.75\) and the baseline
\(U-1.011-7s/32\) in logarithm, which also beats the upper
bound. Therefore a surviving support has \(\alpha<0.03\).

For small \(\alpha\) we use a product over the rotation angles.
On a coset of size \(n=k_0\), let \(\theta_0\) be the nearest
argument to zero in absolute value and \(u_j=|1-e^{i\theta_j}|^2\).
The cyclotomic product gives
\[
 \prod_{j\ne0}u_j\ge (2n/\pi)^2
\]
because \(|\theta_0|\le\pi/n\) and
\(|\sin(n\theta_0/2)/\sin(\theta_0/2)|\ge2n/\pi\), with
the same inequality at \(\theta_0=0\) by continuity.
Among the other arguments let \(r\) be the number with
\(u_j<4\alpha\). The count in that arc gives
\(r/n\le2\arcsin(\sqrt\alpha)/\pi<0.111\).
For the remaining rotations (T18) gives gain at least
\(m(1-\alpha)u_j/(4\alpha)\). If \(\cos\theta_j\ge0\) use
\(am-\Re z\ge bm(u_j/2-\alpha)\ge bm u_j/4\); otherwise
use \(am-\Re z\ge am\).
Apply the baseline to the \(r+1\) other rotations and bound
their omitted \(u_j\) by \(4\alpha\) for \(j\ne0\). With
\(C=(1-\alpha)/(4\alpha)\), this gives
\[
 \log\mathcal L\ge \log m+\log C
    +\frac{2\log(2n/\pi)-g_0-\log C}{n}
    -\frac{r}{n}(g_0+\log(1-\alpha)).
 \tag{T24a}
\]
Use \(n=d\) from now on.

If \(b\le N\), (T10) and (T7) bound the prime support by
\(X=1.03\sqrt M\), since \(a+2d-1<1.03b\le1.03\sqrt M\).
Here \(\log X>8\) and \(V/\log X<2.001\). The digit term
is below \(2.42\) since \(X/N\le1.03\) and \(\log X>U/2\).
The right side of (T14) is thus less than
\(\tfrac12(U+\log b)+1.91\).
For \(d\ge10\), (T24a) and \(g_0\le1.266+0.256U\) give
\[
 \log\mathcal L>0.945U+0.9t-1.70.
\]
Indeed \(\log C>t-1.417\), \(\log C<t\), and the total
coefficient on \(g_0\) in the loss is at most \(0.211\).
Since \(\log b=2t+\log s\le U\), the surplus over the upper
bound exceeds \(0.395U-0.45\log\sigma-3.61>0\).
For \(d\le9\), \(s<0.27\) and the left third yield
\(\log\mathcal L>U+t/3-0.71\). Subtract the upper bound
using \(\log b=t+\log d\) and \(t<U\); the surplus exceeds
\(U/3-(\log9)/2-2.62>0\). Thus \(b>N\) is necessary.

Next suppose \(B=H^4/M<2048\). In (T14) split the prime
sum at \(x=\max(e^8,HV)\), which is at least \(\log K\).
The sum is below \(\log x+0.461\), by the bound before
(T20) below \(x\) and \(\log K/(x-1)<1.001\) above it.

The identities \(B=h^4s^2b/N\) and \(B=s^3h^4e^{2t-U}\)
give \(s<8\sqrt2\) and \(s^3<128e^{U-2t}\), using \(b>N\).
They also give
\[
 V<\tfrac83 U+1.62,\qquad
 H=(B/h)^{1/3}e^{(U+t)/3}<10.08e^{(U+t)/3}.
\]
For the bound on \(V\), use \(M=se^{U+2t}\) and \(t<U\).
Splitting the digit count at \(y=N/16\) gives
\[
 \frac{\omega(K)V}{N}
 \le \frac{\pi(y)V}{N}+\frac{HV^2}{N\log y}
 <0.236+9.46\,U\frac HN.
 \tag{T24b}
\]
The last constants follow at \(U=15\) by the prime count
and the bound on \(V\). The two quantities bounding
\(\pi(y)V/N\) and \(V^2/(U\log y)\) decrease thereafter.

For \(d\ge50\), (T24a) gives \(\log\mathcal L>U+0.98t-1.72\),
since \(g_0<2.459\) and \(g_0/d<0.025\).
Subtracting the prime and digit bounds leaves a surplus greater than
\[
 \tfrac23 U+(0.98-\tfrac13)t
 -\log(10.08(\tfrac83 U+1.62))-2.417
 -95.3568\,Ue^{(t-2U)/3}.
\]
This is concave in \(t\); at the endpoints \(3.5,U\)
the values exceed \(3.59,1.60\) at \(U=15\) and increase
thereafter. In each case the linear term grows faster
than the logarithm and the last cost decreases.
For \(d<50\), since \(b>N\) we have \(s<2500/N<0.001\)
and \(t>U-\log50\). The baseline and left third give
\(\log\mathcal L>U+t/3-0.71\).
Here \(M<50N^2\), \(H<(102400)^{1/4}N^{1/2}\),
and \(V<2U+\log51\). The digit bound in (T24b) is less
than \(1.65\). The total upper bound
\[
 U/2+\tfrac14\log102400+\log(2U+\log51)+0.461+1.65
\]
is less than \(4U/3-(\log50)/3-0.71\) at \(U=15\)
and increases more slowly. Hence \(B\ge2048\) is necessary.

Finally, if \(b>N\) and \(\alpha<0.03\), the prime support
in (T7) is at most \(1.03b\) by (T10).
The left side of (T7a) divided by \(H\), after subtracting
its prime sum, exceeds
\[
 \log(aN)-\log(1.03b)+0.54>U+0.47.
\]
The remaining sum in (T7a) is at most
\(\pi(1.03b)\log M<2.24b\), since \(\log M<2\log b\)
and \(\log b>U\ge15\). In summary, at \(U\ge15\)
a possible common factor requires
\[
 k_0=k=d,\qquad b>N,\qquad \alpha<0.03,\qquad
 B\ge2048,\qquad
 \frac Hb<\frac{2.24}{U+0.47}.
 \tag{T24c}
\]

### Joined degree exclusions from Newton polygons

Assume the ordinary residual support constraints (T5) and \(N\ge20000\).
For a monic irreducible common factor \(F\) of even degree \(h\ge2\), put
\(H=dh\), \(b=a+d\), \(M=bN\), \(L=\log M\), and \(t=\log H\).
Then the following two exclusions apply.

The bounds are
\[
\begin{aligned}
  H &< \max\{b,(6M)^{1/4}\} && (L\ge18),\\
  H &< M^{1/4}/3 && (L\ge24,\ b\le H).
\end{aligned}
\tag{T27}
\]
The second bound gives coprimality for \(a\le d\), \(e^{24}\le bN\),
and \(20{,}000\le N\le648d^3\), since \(b\le 2d\le H\) and
\(M\le1296d^4\). Supports outside (T5) already satisfy coprimality.

**Terminal prime support.**

For each order \(m\in\{N-1,N\}\), use the backwards ordinates
\[
 U_m(j)=v_p\left(\prod_{i=bm-dj+1}^{bm}i\right)-v_p\binom mj .
\]
The factorial count gives \(U_m(j)<dj/(p-1)+3\) for \(j>0\)
when \(p^4>M\). A positive slope belonging to \(F\) satisfies
\(\lambda\ge1/h\), and its edge uses at least two maximizing
indices of \(U_m(j)-\lambda j\), with maximum at least zero.
If \(p>2H+d\), all these indices are less than \(6h\).
If neither \([bm-6H,bm]\) contains a positive multiple of
\(p^3\), the affine error improves from 3 to 2 on these indices,
so they are less than \(4h<p\). The same \(j<p\) conclusion
holds with cubes allowed if \(p>3H+d\).

In either case, \(j!\) is a unit, so the binomial valuation
is nondecreasing in the relevant indices. Consecutive factorial
multiples are reached at indices spaced at least
\(\lfloor p/d\rfloor\) apart. This spacing exceeds \(2h\)
in the noncubic case, or \(3h\) in the cubic case with
\(p>3H+d\). The weighted values at successive jumps decrease
strictly because the jump heights are at most two or three,
respectively. They also decrease strictly between jumps.
Thus any selected edge ties zero to the first jump at \(w_m\),
reaching a multiple \(u_m\) of \(p\). Its height \(r_m\) is
at most two or three, respectively, and \(\lambda=r_m/w_m\).
The terminal offsets satisfy
\[
 e_m=bm-u_m\in[d(w_m-1),dw_m-1].
\]

The multiples \(u_N,u_{N-1}\) are distinct. If they were the
same, \(b>d\) would give \(w_N>w_{N-1}\). The longer upper
binomial numerator interval would contain the lower one, and
both factorial intervals would contain only this first multiple.
This gives \(r_N\le r_{N-1}\), contrary to equality of their
positive slopes. Also \(|w_N-w_{N-1}|=|r_N-r_{N-1}|/\lambda\).
Hence
\[
 0<|u_N-u_{N-1}|
   \le b+d-1+
   \begin{cases} H&\text{in the noncubic case},\\
                 2H&\text{with cubes allowed}.
   \end{cases}
\]
For \(b\le H\), this contradicts \(p\mid u_N-u_{N-1}\) once
\(p>2H+d\) in the noncubic case or \(p>3H+d\) in the cubic case.

**Norm comparison.**

Assume \(b\le H\). Equation (T10) gives \(2H<\sqrt M\),
because \(a+2d-1<2H\). Thus \(t<L/2-\log2\).
At any prime the usual hull bound gives
\(v_p(F(0))\le H/(p-1)+\lfloor\log_p M\rfloor\).
At \(p>2H\) the integer bound is
\(v_p(F(0))\le\lfloor\log_p M\rfloor\).
The strict lower norm estimate (T13) gives
\(\log|F(0)|/H>\log(aN)>L-\log4\), using \(b<4a\).

First suppose \(H^4/M\ge6\), and put \(P=2H+d\).
All \(p>P\) have \(p^4>M\), and the terminal support argument
leaves only primes with a positive multiple \(q p^3\) in one
of the two six-\(H\) terminal intervals. Each interval has
at most one such prime per \(q\le M/P^3\), because the spacing
of consecutive cubes above this threshold exceeds \(6H\).
Their normalized norm charge is at most \(2ML/(HP^3)\le L/24\).
For the other case, \(H\ge M^{1/4}/3\), take \(P=3H+d\).
Then \(p^4>M\) for \(p>P\) and there are no norm primes above
\(P\).

In either case put \(r=P/H\) and \(c=1.2762\). The Dusart prime count
used in (T20) gives
\[
 \frac{\pi(P)L}{H}\le
 \frac{rL}{t+\log r}+\frac{crL}{(t+\log r)^2}.
\]
For \(t+\log r>2\), both terms and \(\log r\) increase with \(r\).
Use \(r_0=5/2\) in the first case and \(r_0=7/2\) in the second.
Define
\[
 f(t)=L-t-\frac{r_0L}{t+\log r_0}
             -\frac{cr_0L}{(t+\log r_0)^2}.
\]
This is concave in \(t\), so its minimum on the relevant
interval is at an endpoint.

In the first case
\(t\in[(L+\log6)/4,L/2-\log2]\).
The prime sum bound
\(\sum_{p\le P}\log p/(p-1)\le\log P+0.2\)
would require \(f(t)<\log10+0.2+L/24\).
At \(L=18\) the two endpoint margins exceed \(0.45\) and \(0.88\).
In the second case \(t\in[L/4-\log3,L/2-\log2]\) and
\(P\ge e^{L/4}>350\). The bound
\(\sum_{p\le P}\log p/(p-1)\le\log P-0.49\)
would require \(f(t)<\log14-0.49\).
At \(L=24\) the two margins exceed \(0.46\) and \(3.17\).
All four endpoint margins increase with \(L\). At an endpoint,
\(t+\log r_0=\beta L+C\) for \(\beta\in\{1/4,1/2\}\)
and \(C>0\). The derivative of \(r_0L/(\beta L+C)\) is
\(r_0C/(\beta L+C)^2<0.2\) on the stated ranges, and that
of \(cr_0L/(\beta L+C)^2\) is negative. The affine
contribution to the derivative is at least \(1/2-1/24\).
This proves both contradictions.

**Prime-sum bounds for small arguments.**

For \(P\ge3000\) use the Dusart bound preceding (T20).
For smaller \(P\), the excess of the prime sum over \(\log P\)
decreases between primes. The following exact-rational check
uses normalized logarithms \(r\in[1,2]\). The remainder of
the ten-term sum with \(z=(r-1)/(r+1)\) is at most
\(2z^{21}/[21(1-z^2)]\).

    from fractions import Fraction as Q
    from math import isqrt

    D = 10**8
    def log_interval(r):
        z = (r - 1) / (r + 1)
        lo = 2 * sum(z**(2*j+1) / Q(2*j+1) for j in range(10))
        return lo, lo + 2*z**21 / (21*(1-z*z))

    lower_2, upper_2 = log_interval(Q(2))
    total = 0
    for n in range(2, 3001):
        prime = all(n % q for q in range(2, isqrt(n)+1))
        if not prime and n != 350:
            continue
        k = n.bit_length() - 1
        lo, hi = log_interval(Q(n, 2**k))
        lower = int(D*(lo + k*lower_2))
        upper = int(D*(hi + k*upper_2)) + 1
        if prime:
            total += (upper + n - 2) // (n - 1)
        assert total - lower <= D//5
        if n >= 350:
            assert total - lower <= -49*D//100

### A joined degree bound when the terminal shift grows

Assume the ordinary residual exponent inequalities. Let \(F\)
be a monic irreducible common factor of even degree \(h\ge2\). Put
\(H=dh\), \(b=a+d\), \(M=bN\), \(B=H^4/M\),
\(L=\log M\), \(U=\log N\), and write \(\gamma=b/H\).
First suppose \(B\ge6\). Use the reverse edges from (T27).

**Shared reverse edges.**

Take \(P\ge2H+d\) and \(p>P\). Then \(p^4>M\).
Call \(p\) exceptional if either \([bm-6H,bm]\),
\(m=N-1,N\), contains a positive multiple of \(p^3\).
As in (T27), the affine bound on the reverse ordinate gives
maximizing indices \(j<6h\) for every selected positive slope.
At nonexceptional primes the bound improves to \(j<4h<p\).
The binomial valuation is nondecreasing there. Successive
factorial jumps are spaced at least \(2h+1\) indices apart
with heights at most two. Hence the only selected positive
slope takes the edge from zero to the first reached multiple
\(u_m\) of \(p\), at width \(w_m\), with height \(r_m\in\{1,2\}\).

Put \(e_m=bm-u_m\in[d(w_m-1),dw_m-1]\).
If \(u_N=u_{N-1}\), then \(w_N>w_{N-1}\) since \(b>d\).
The upper binomial numerator then contains the lower one,
whereas the factorial valuation is the same, contradicting
the common positive slope. Also
\[
 |w_N-w_{N-1}|=|r_N-r_{N-1}|/\lambda\le h .
\]
Thus \(\Delta=u_N-u_{N-1}\ne0\) belongs to
\([b-H-d+1,b+H+d-1]\).
Both terminal multiples are divisible by
\(p^{v_p(F(0))}\), as the norm exponent is at most
\(\min(r_N,r_{N-1})\). The same divisibility holds for \(\Delta\).

The exceptional primes have normalized log norm charge at
most \(2ML/(HP^3)\), by the cube multiplier count. For primes
at most \(P\), use the affine hull estimate and the prime sum
and count from (T27). With \(r=P/H\), this gives charge at
most
\[
 \log H+\log r+0.2+
 \frac{rL}{\log P}\left(1+\frac{1.2762}{\log P}\right).
\]
The strict lower norm bound is
\(\log|F(0)|/H>U+\log(\gamma-1/h)+\log H\).

**Ratios from one to three.**

Use \(P=2H+d\). If \(1\le\gamma\le c\le3\) then every
nonexceptional norm prime above \(P\) equals \(\Delta\)
itself. Indeed the bound on \(\Delta\) gives \(|\Delta|<2p\),
and its negative range has absolute value below \(P\).
The possible primes lie in \((P,P+(c-1)H]\).
Their log sum divided by \(H\) is at most
\[
 E(c)=(c+1.5)\log(c+1.5)-2.5\log2.5-(c-1)\log(c-1)
\]
with the last term zero at \(c=1\).
This is the entropy upper bound for the binomial coefficient
with numerator interval \((P,P+\lfloor(c-1)H\rfloor]\).
Its denominator is prime to these primes.

The prime-count term is less than 10. Increasing \(r\) to
\(5/2\) increases the prime-count bound and
\(\log P\ge L/4+\log(5/2)+(\log6)/4\) at that value.
The constant after \(L/4\) exceeds \(1.2762\),
so the bound is less than \(4(5/2)\).
The cube charge is at most \(L/24<(3U+5)/24\).
Using \(\gamma-1/h\ge1/2\), a common factor would require
\[
 U-\log2<
 \log2.5+0.2+10+E(c)+(3U+5)/24 .
\]
This excludes \(c\le1.5\) at \(U\ge16\), \(c\le2.5\) at
\(U\ge17\), and \(c\le3\) at \(U\ge18\). The respective
entropy constants are less than \(1.352,2.647,3.092\).

**Ratios at least three.**

Use \(r=2+2/h\) and \(P=rH\). The nonexceptional norm over
\(p>P\) divides
\[
 \binom{b+H+d-1}{2H+2d-1}.
\]
The numerator covers \(\Delta\), it has positive entries
because \(\gamma\ge3\), and the denominator is prime to
all these primes. Its logarithm divided by \(H\) is at most
\[
 r\left(1+\log\frac{\gamma+1+1/h}{r}\right).
\]

The prime-count term is at most \(4r+0.04\). To see this
uniformly, write \(C_0=\log r+(\log6)/4\),
\(w=4C_0/(L+4C_0)\), and \(c_0=1.2762\).
Relative to \(4r\) the term is at most
\((1-w)(1+(c_0/C_0)w)\), whose maximum on \([0,1]\) is bounded
by \(1+(c_0-C_0)_+^2/(4c_0C_0)\).
Since \(r\in(2,3]\) and \(C_0>1.14108\), the resulting
additive error is less than \(0.04\).
The cube charge is at most \(2L/(Br^3)<(U+5/3)/r^3\).

Combining the bounds gives, with
\[
 \begin{split}
 A(r)&=\frac1r+\frac1{r^2}+\frac1{r^3},\\
 C(r)&=5r-(r-1)\log r+0.24+\frac5{3r^3}
       +r\log(1+r/6)-\log(1-(r-2)/6),
 \end{split}
\]
the necessary condition
\[
 \log\gamma>A(r)U-\frac{C(r)}{r-1}.
\]
Here the two ratio-dependent terms in the norm rearrangement
are \(r\log(1+(1+1/h)/\gamma)-\log(1-1/(h\gamma))\),
which are bounded by their values at \(\gamma=3\).

For \(U\ge13\) the displayed lower bound decreases with \(r\).
Write \(C(r)/(r-1)=5r/(r-1)-\log r+E_0(r)/(r-1)\).
Direct differentiation gives \(E_0'>0\) and \(0<E_0<1.701\)
for \(2\le r\le3\). Since
\[
 (r-1)^2|A'(r)|=1-4/r^3+3/r^4\ge11/16,
\]
the derivative of the lower bound, times \((r-1)^2\),
is at most \(-11U/16+6.701+2<0\).
Consequently the following rounded lower bounds hold.

| Degree | Necessary \(\log(b/H)\) exceeds |
|---:|:---|
| \(h=2\) | \(13U/27-7.252\) |
| \(h\ge4\) | \(78U/125-8.287\) |
| \(h\ge6\) | \(237U/343-8.799\) |
| \(h\ge8\) | \(532U/729-9.106\) |

For \(U\ge18\), (T27) covers \(\gamma\le1\) and the first
part of this note covers \(1\le\gamma\le3\).
Thus the displayed bound with \(A(r),C(r)\) applies to all
surviving factors with \(B\ge6\). Without assuming \(B\ge6\)
in advance this yields
\[
 H<\max\!\left\{(6M)^{1/4},\
 b\exp\!\left[-A(r)U+\frac{C(r)}{r-1}\right]\right\}.
 \tag{T28}
\]

For example, if \(A_0=a/d\),
\(d^3\ge(3/8)(A_0+1)N\) and
\(A_0+1\le2\exp(13U/27-7.252)\), then all degrees are
excluded. If \(d^3\ge(3/128)(A_0+1)N\) and
\(A_0+1\le4\exp(78U/125-8.287)\), then all degrees
at least four are excluded.

### Common-factor degree and scale

Assume the reduced conditions (T24c) at \(U\ge15\).
In particular \(B\ge2048\), \(b>2H\), and
\(\log P\ge(\log M)/4+\log2+\tfrac14\log2048>8\),
so (T20) applies with \(P=2H+d\).
Write \(A=a/d\) and
\(W=e^{U/2}/\sqrt{2U+2}>308\), so \(A+1>W\) by \(b>N\).
For \(h=2\), \(\Delta-U=\log(A/(A+2))-1-\log 2>-1.70\).

For \(h\ge4\), (T24c) gives \(h< h_0=2.24(A+1)/(U+0.47)\).
At fixed \(A\), \(\Delta-U=\log(A/h)-(2/h)(1+\log(A+2))\)
has its minimum on \([4,h_0]\) at an endpoint. Its derivative
has the sign of \(2(1+\log(A+2))-h\).
The endpoint \(h=4\) increases with \(A\) and exceeds \(0.97\)
for \(A>307\). At \(h=h_0\) the expression equals
\[
 \log\frac{U+0.47}{2.24}
 +\log\left(1-\frac1{A+1}\right)
 -\frac{2(U+0.47)}{2.24}\frac{1+\log(A+2)}{A+1}.
\]
It increases with \(A\), so substitute \(A+1=W\).
The negative term decreases in magnitude with \(U\ge15\),
by logarithmic differentiation using \(0.46<W'/W<0.5\).
The resulting lower bound exceeds \(1.6\) at \(U=15\)
and increases. Hence \(\Delta>U+0.95\) for \(h\ge4\).

The prime count in (T20) gives \(\pi(P)\log M/H<10\),
since \(P/H\le2.5\) and \(\log P>(\log M)/4+1.2762\).
Using \(S\le6d\log(2M)\) and \(V<3U+5\) in (T20),
for \(h\ge4\) we obtain
\[
 h<\frac{6(3U+6)}{U+0.95-\log2.5-10-(3U+5)/8192}<61 .
\]
The upper expression is decreasing for \(U\ge15\). Thus
every possible factor has \(h\le60\), including \(h=2\).

We next use the arithmetic progression bound for the
Padé digits. At any cutoff \(P'\ge2H+d\), primes
\(p>P'\) have \(p^4>M\) and \(p\nmid d\).
For a positive lifted valuation \(e\), we have
\(e\ge1/H\) and \(e-\mu\ge1/(2H)\).
Counting in each residue class modulo a prime power
bounds the progression valuation in (T17b) by
\(N\mu+\lfloor\log_p(bm+j)\rfloor\).
With \(L_p=\lfloor\log_p M\rfloor\), the first increase
of the floor beyond \(L_p\) needs \(j-1\ge b>2H\).
Subsequent increases are more than \(2H\) steps apart.
The penalty in (T17b) therefore absorbs every increase.
The supremum there is at most \(N\mu+L_p\), which keeps the
same main term as (T17) and gives only \(L_p\le3\) digits.
Using (T17) below \(P'\), the digit term can be replaced by
\[
 \frac{\pi(P')V}{N}+
 \frac3N\sum_{p>P'}v_p(F(0))\log p .
 \tag{T26}
\]
In fact, at each large norm prime the loss is at most
\(L_p\log p/N\) after averaging over the lifted roots.

Put \(\delta=H/N\) and take \(r=2+2/h\), \(P'=rH\).
By (T24c), \(c=b/H>6.9\). The nonexceptional part of
the last sum in (T26) has log mass at most
\(Hr(1+\log((c+r/2)/r))\) by the binomial argument of
(T28). This expression increases with \(r\le3\).
Also \(c=e^t/h\), so
\[
 r(1+\log((c+r/2)/r))
 \le 3(t+1-\log6+\log(1+1.5/6.9))<3(t-0.595).
\]
The primes below \(P'\) give at most \(12\delta\) digits,
since \(V<\log M+0.001\) and
\(\log P'>(\log M)/4+2.59\).
The exceptional primes give at most \(6\delta\log M/(Br^3)\).
Consequently the total digit loss is less than
\(\delta(9t+J(U))\), with \(J(U)=6.65+6(3U+5)/16384\),
since \(\log M<3U+5\).

For \(d\ge50\), (T24a) gives \(\log\mathcal L>U+0.98t-1.64\).
Use \(\log C>t-1.417\), \(g_0/d<0.025\), and
\((2/\pi)\arcsin(e^{-t/2})g_0<0.19\).
For the last bound use \(t>\tfrac12(U-\log(2U+2))\)
from \(b>N\), \(g_0<1.266+0.256U\), and
\(\arcsin v\le v/\sqrt{1-v^2}\). The product of this
last bound and the bound on \(g_0\) decreases from \(U=15\).
The prime sum in (T14) is below \(\log(HV)+0.461\)
by the split at \(x=HV\), since \(\log(HV)>8\).
If \(\delta\le0.045\), the lower bound on the gain exceeds
these two upper bounds by more than
\[
 (.98-9\delta)t-\log\delta-\log(3U+5)-2.101-J(U)\delta .
\]
Use \(\delta=0.045\) and
\(t=(U-\log(2U+2))/2\) to lower bound the expression.
It is greater than \(0.10\) at \(U=15\) and increases
thereafter.

For \(d<50\), \(h\le60\) gives \(H<3000\),
so \(\delta<3000e^{-U}\). The digit bound
\(\delta(9U+J(U))\) is below \(0.14\).
As \(b>N\), the lower bound \(\log\mathcal L>U+t/3-0.71\)
from the small-\(d\) case of (T24c) applies.
It beats \(\log(3000(3U+5))+0.461+0.14\) at \(U=15\)
and increases more slowly than the gain thereafter. We have proved
\[
 h\le60,\qquad H>0.045N
 \qquad (U\ge15).
 \tag{T29}
\]

### An explicit terminal cutoff

Suppose \(U=\log N\ge15\) and an irreducible common factor
satisfies (T24c) and (T29). Put
\[
 \delta=H/N,\qquad v=sh/\delta=e^{U-t},\qquad L=\log M .
\]
The strengthened Gamma annulus following (T5a) gives \(b/d<(5/6)N\),
so \(v>6/5\).
The identities
\[
 M=\frac{\delta N^3}{hv},\qquad B=hv\delta^3 N
 \tag{T30}
\]
keep the size of the terminal intervals and the Padé scale
in the same estimate.

For \(h=2\), the improved surplus in (T20) is
\(\Delta_*=U-1+\log(A/(A+2))>U-1.01\), since \(A>307\)
by (T29)'s degree argument. For \(h\ge4\) the same lower
bound follows from \(\Delta>U+0.95\).
Using \(\log P>U+\log(0.09)\) and \(V<3.24U\) gives
\[
 \frac{S}{HL}>
 \frac{U-1.01-\log2.5
   -\frac{2.5(3.24U)}{U+\log(0.09)}
      \left(1+\frac{1.2762}{U+\log(0.09)}\right)
   -(3.24U)/8192}{3.24U}>0.05 .
\]
The first comparison uses \(B\ge2048\) in (T20).
For the second, the displayed ratio increases from \(U=15\),
where it exceeds \(0.0502\). The multiplier bound (T19b) also
gives \(S/(HL)\le1/(s(2h+1)^2\delta)\). Together with \(v>6/5\)
and \(h(2+1/h)^2\ge12.5\), this gives \(\delta^2<4/3\).
In particular \(\delta<1.4\) and \(L<3U\).

**A bound from prime powers.** Use \(\vartheta(x)=\sum_{p\le x}\log p\)
and Rosser and Schoenfeld's Theorems 9 and 10, namely
\(\vartheta(x)<1.017x\) for \(x>0\) and
\(\vartheta(x)>0.98x\) for \(x\ge7481\).
[Rosser and Schoenfeld, *Approximate formulas for some functions of prime numbers*](https://doi.org/10.1215/ijm/1255631807).

At \(p\le P\), the affine hull bound before (T19) contributes
at most \(H/(p-1)+\lfloor\log_p M\rfloor\) to \(v_p(F(0))\).
The extra term summed over these primes is exactly
\(\sum_{k\le L/\log 2}\vartheta(\min(P,M^{1/k}))\).
For constants \(\ell\le\delta\le u\) and \(v\ge w\), define
\[
 \begin{aligned}
 V_0&=3U+\log(u/(2w))+0.001,\\
 Z&=\frac1{\ell}\sum_{k=4}^{8}
         (u/(2w))^{1/k}e^{(3/k-1)U}
       +\frac{3U}{\ell\log2}(u/(2w))^{1/9}e^{-2U/3},\\
 C_0&=1.017\bigl(5+\min(2.5,(2w\ell^2)^{-1/3})+Z\bigr),\\
 G(U,\ell,u,w)&=U-1.01-\log2.5-C_0-V_0/8192 .
 \end{aligned}
\]
Here \(V<V_0\) follows from (T30), \(h\ge2\), and
\(V-L<0.001\). The terms at \(k=1,2\) cost at most
\(1.017(2P/H)\le1.017\cdot5\) after division by \(H\).
The term at \(k=3\) uses \(M^{1/3}/H\le(2w\ell^2)^{-1/3}\).
The remaining terms cost at most \(1.017 Z\) by (T30),
using \(L/\log2<3U/\log2\) for the terms \(k\ge9\).
Using \(C_0\) for the extra terms in the proof of (T20),
and keeping its other bounds with \(\Delta_*\), gives
\[
 S/H>G(U,\ell,u,w).
 \tag{T31}
\]

By the weighted count (T22a) and (T30), on the same domain
\(S/H<W(U,\ell,w)\), where
\[
 W(U,x,w)=e^{-U/3}\left[
  \frac{14(U+\log(2.5x)+3.06)}{(5w)^{1/3}x}+
  \frac{6.8(U+\tfrac13\log(x/(2w))+1)}{(2w)^{1/18}x^{4/9}}
 \right].
\]
For the first term use \(P/H\le2.5\) in the logarithm and
\((P/H)h=2h+1\ge5\) in the denominator.
Before substituting \(\ell\), the right side is bounded
by \(W(U,\delta,w)\). Both terms decrease in \(x\) on
the stated ranges, while \(xW(U,x,w)\) increases.
This follows by differentiating the products of the
logarithmic factors and the indicated powers of \(x\).
Thus \(\delta S/H<uW(U,u,w)\) as well.
At fixed \(\ell,u,w\) used below, \(W(U,\ell,w)\) decreases
and \(G(U,\ell,u,w)\) increases for \(U\ge15\).
In \(W\) each affine factor multiplying \(e^{-U/3}\)
is greater than \(3\). In \(G\), \(Z\) decreases and
\(U-V_0/8192\) increases.

At \((\ell,u,w)=(0.21,1.4,6/5)\) the difference \(G-W\)
exceeds \(0.074\) at \(U=15\). This forces \(\delta<0.21\).
Also \(G(U,0.045,0.21,6/5)/V_0\) increases with \(U\).
Indeed \(C_0\) decreases and
\(3(1.01+\log2.5+C_0)+\log(0.21/2.4)+0.001>0\).
The multiplier bound and (T31) now give \(sh<0.675/\delta\):
the upper constant \((2/25)V_0/G\) is below \(0.675\)
at \(U=15\). For \(\delta\ge0.10\) we can instead bound
\(sh=\delta v\). In each corresponding row of the table
below, \(G>W\) at the displayed \(w\) would contradict
\(v\ge w\). Hence \(sh<uw\).

**The Padé comparison.** By (T24a), at \(U\ge15\)
\[
 \log\mathcal L>U+t-1.60 .
\]
Indeed \(t>(U-\log(2U+2))/2\) and \(d>0.045N/60\).
The term \((t+g_0)/d\) is below \(0.009\), and the
angle term is below \(0.19\), by the decreasing bound
used for (T29). Also \(\log C>t-1.39\) and
\(\log m>U-0.001\).

For a row \([\ell,u]\), now use \(V_0\) at \(w=6/5\).
The total log norm above \(P\), divided by \(H\),
is less than \(U+0.127+W(U,\ell,6/5)+V_0/8192\).
Indeed \(A+2<(5/6)N+1/6\) gives \(\log(A+2)<U-0.18\).
So the equal-height contribution, divided by \(H\),
is below \((2/h)(U+0.127)\) by the sharpened binomial bound.
The unequal and exceptional contributions use (T22a)
and the cube bound in (T20).
Set
\[
 Q=\min\{V_0,\ U+0.127+W(U,\ell,6/5)+V_0/8192\}.
\]
Both entries bound the log norm above \(P\) divided by \(H\).
Also \(HQ>P\) and \(\log(HQ)>8\).
Splitting the prime sum in (T14) at \(HQ\), the primes
above it contribute at most the log norm above \(P\)
divided by \(HQ-1\). Thus the prime sum is less than
\(\log(HQ)+0.461\).

In (T26), the exact large-prime loss has \(L_p\le2\)
unless \(p\le M^{1/3}\), where \(L_p\le3\). The log mass
charged a third time is at most
\(\max(0,\vartheta(M^{1/3})-\vartheta(P))\).
The upper and lower \(\vartheta\) bounds apply since
\(P>0.09N>7481\). By (T30), this extra loss after division
by \(N\) is at most
\(\max(0,1.017(u/2.4)^{1/3}-1.96\ell)\).
Combining the three norm charges above \(P\) and the
prime count for the digits below it, the digit loss is bounded by
\[
 D(U,\ell,u)=
 u\left[\frac{2.5V_0}{z}\left(1+\frac{1.2762}{z}\right)
          +2(U+0.127)+V_0/4096\right]
 +2uW(U,u,6/5)+\max(0,1.017(u/2.4)^{1/3}-1.96\ell),
 \quad z=U+\log(2.5\ell).
 \tag{T32}
\]
For the head count \(r/\log(rH)(1+1.2762/\log(rH))\)
increases with \(r\le2.5\), and we use \(\delta\ge\ell\)
in its denominator and \(\delta\le u\) in its numerator.

With \(H=she^t\), the surplus of the log gain over
the two bounds exceeds \(U-\log(sh Q)-2.061-D(U,\ell,u)\).
The table uses \(J\) as an upper bound for \(sh\).
The last three rows bound \(v\) by (T31), using the margins
in the fourth column. The first row uses \(0.675/\delta\le15\).
\[
\begin{array}{c|c|r|c|c}
 [\ell,u]&w&J &G-W\text{ exceeds at }U=15
   &U-\log(JQ)-2.061-D\text{ exceeds at }U=15\\\hline
 [0.045,0.10]&\text{--}&15    &\text{--}&0.36\\
 [0.10,0.14] &9        &1.26  &0.109&1.63\\
 [0.14,0.18] &3.6      &0.648 &0.117&0.83\\
 [0.18,0.21] &1.8      &0.378 &0.093&0.30
\end{array}
\]
The surplus lower bounds increase for \(U\ge15\).
For either branch of \(Q\), its logarithmic growth rate
is below \(0.08\), using \(W'<0\), \(V_0>3U-4\),
and \(U+0.127+W+V_0/8192>U\).
In \(D\), the rational factor \(V_0/z\,(1+1.2762/z)\)
decreases. In fact it is \(3(1-k/z)(1+1.2762/z)\) with
\(k=\log(2.5\ell)-(\log(u/2.4)+0.001)/3<0.02\) in every row.
The term \(W(U,u,6/5)\) decreases as well. The remaining
\(U\)-dependent terms in \(D\) grow at rate at most
\(0.21(2+3/4096)<0.43\). The positive margins contradict
(T14) for all the rows.

Thus \(N\ge\lceil e^{15}\rceil=3{,}269{,}018\) is sufficient.
The integer rounding follows by summing the exponential series
through index \(100\) and bounding its tail by
\((15^{101}/101!)/(1-15/102)\).

### A lower terminal cutoff from verified prime gaps

Use the analytic cutoff above to restrict this argument to
\(25\le N\le3{,}269{,}017\). Write \(G\) for a bound on gaps
between primes in the intervals used below.

For an even record prime gap \(G\), let \(X_G\) be the lower endpoint of the first larger record gap. Every real interval \((u,v]\) with \(u\ge2\), \(v<X_G\), and \(v-u\ge G\) contains a prime. Indeed the consecutive primes straddling \(u\) would otherwise have gap greater than \(v-u\), and the first of them is below \(X_G\). The values below come from Table 8 of Oliveira e Silva, Herzog, and Pardi. The table lists lower prime endpoints and computes the gaps through \(4\cdot10^{18}\). [Oliveira e Silva, Herzog, and Pardi, Table 8, printed page 2050](https://doi.org/10.1090/S0025-5718-2013-02787-1).

**Bounds on the gap and common degree.**

Suppose
\[
 G\ge4\sigma_N,\qquad G\left(\frac{N(5N-11)}6+2\right)<X_G.
\]
Then
\[
 d<G/2,\qquad H<G\quad\text{when (T10) applies}.
 \tag{T33}
\]

For the first claim let \(g=G/2\). If \(d\ge g\), take \(r=\lfloor d/g\rfloor\). Then \(g\le d/r<2g\) and the interval
\[
 (aN/r,(aN+2d)/r]
\]
has length at least \(G\) and upper endpoint less than \(G(N(5N-11)/6+2)<X_G\). Its left endpoint exceeds \(dN\): the bound \(d^2/(a+d)<\sigma_N\) gives \(a>d^2/\sigma_N-d\ge d^2/(2\sigma_N)\ge d^2/g\ge rd\). Thus a prime \(p>dN\) divides \(rp\) in the first \(2d\) factorial terms, contrary to (T7). The exact upper endpoint may be used in place of its uniform estimate for given \(a,d,N\).

For the second claim, let \(C=\max\{a+2d-1,\sqrt M\}\); by (T10), \(2H<C\). If \(H\ge G\), take \(r=\lfloor H/G\rfloor\). Now \(G\le H/r<2G\) and
\[
 (aN+H)/r<2G\left(\frac{N(5N-11)}{6h}+1\right)
 \le G\left(\frac{N(5N-11)}6+2\right).
\]
The left endpoint exceeds \(aG/d>2\), since \(r\le H/G<dN/G\).
Take a prime \(p\in(aN/r,(aN+H)/r]\). If \(C=\sqrt M\), then
\[
 p>aGN/H>2aG\sqrt M/(a+d)>\sqrt M,
\]
using \(d<3a\) and \(G\ge52\) for the table. If \(C=a+2d-1\), the inequality \(2H<C\) and \(h\ge2\) imply \(a>2d\). By the first claim, \(r\le H/G<dN/G<N/2\), so \(p>2a>C\). In both cases \(p>\max\{2H,a+2d-1,\sqrt M\}\) and \(p\) divides \(rp\) in the first \(H\) factorial terms, contrary to (T7).

Use the first applicable row. Both displayed scalar inequalities at the row endpoint imply the inequalities throughout that row, using \(2\log N+3\) as an increasing majorant for \(\sigma_N\).

| Inclusive upper endpoint for \(N\) | \(G\) | \(X_G\) |
|---:|---:|---:|
| 50 | 72 | 155,921 |
| 100 | 114 | 1,349,533 |
| 300 | 154 | 17,051,707 |
| 1,000 | 222 | 189,695,659 |
| 2,000 | 282 | 1,294,268,491 |
| 6,000 | 354 | 10,726,904,659 |
| 10,000 | 456 | 42,652,618,343 |
| 19,999 | 468 | 182,226,896,239 |
| 32,767 | 516 | 461,690,510,011 |
| 65,535 | 602 | 2,614,941,710,599 |
| 3,269,017 | 1132 | 43,841,547,845,541,059 |

The row ending at \(19{,}999\) gives \(d\le233\) below the
new cutoff. All tested intervals, including their upper endpoints,
lie below their indicated \(X_G\).

**Extending (T10) to \(N\ge275\).**

It suffices to handle \(275\le N<20{,}000\) and \(a<2d\); the proof of (T10) already handles \(a\ge2d\), and (T11) handles all larger \(N\). Let \(X=aN\), so \(1375\le X<2{,}640{,}000\), since \(a<2d\), \(d^2/(a+d)<2\log N+2\), and \(d<66\). Set
\[
 w=
 \begin{cases}
 34,&X\le7160,\\
 36,&7160<X\le11750,\\
 148,&11750<X.
 \end{cases}
\]
Then \(w<X/40\). Each of the windows
\[
 (X,X+w],\quad
 (X+2w+4d,X+3w+4d],\quad
 (X+6w+12d,X+7w+12d]
\]
has length \(w\) and contains a prime. Indeed the envelope ends at \(X+E\), \(E=7w+12d\), and \(E/X<7/40+36/275<1/3\). Its endpoint lies below the first record gap exceeding \(w\) in each range: \(9551\), \(15683\), and \(4{,}652{,}353\), respectively.

As in (T11), discard any window meeting a forbidden interval
\[
 I_j=[(M-b)/j,(M+2d)/j]\quad (1\le j\le4).
\]
The envelope has no intersection with \(I_1\), since \(a<2d\) gives \(M-b>(3/2)(1-1/N)X>4X/3\). If \(I_j\) with \(2\le j\le4\) meets the envelope, \(M+2d>jX\). For \(j'<j\) this gives
\[
 \frac{M-b}{j'}>\frac{jX-(a+3d)}{j-1}
 \ge X+\frac{X-(a+3d)}3>X+E,
\]
where \(E/X<7/40+36/275<(1-10/275)/3\) and \(d<3a\). Thus at most one forbidden interval meets the envelope. Such an interval has length at most \((a+3d)/2<5d/2\), less than the separation between windows. Two windows remain.

Apply the same Newton analysis as (T11) to a prime \(p\) in either remaining window. The part of that proof using the floor ordinates and the forbidden intervals is independent of the window endpoints. In particular a candidate \(h\) satisfies
\[
 h=h_0+Jp/d+\epsilon,\quad
 0\le h_0<(p-X)/d,\quad J\in\{0,1,2\},
\]
with \(\epsilon=0\) for \(J=0\) and \(|\epsilon|<2\) otherwise. For two remaining primes \(p_1<p_2\), the windows give \(p_2-p_1>p_1-X+4d\). The intervals for equal positive \(J\) are therefore disjoint across the primes. Those for distinct \(J\) are also disjoint, since \(X>2E+4d\), by \(2(7/40+36/275)+12/275<1\). Thus \(J=0\) at both primes. The smaller usable prime lies in one of the first two windows, so
\[
 H<p_1-X\le3w+4d<\left(\frac3{40}+\frac{12}{275}\right)X<0.12X.
\]
This implies \(X+1>\max\{H+13,279H/262\}\), the starting-point hypothesis of (T10).

Thus (T33) gives both degree bounds at every \(N\ge275\) in
the table. It also gives an order-only test below \(20{,}000\).
If a prime \(p\ge3G\) divides \(N(N-1)\) at \(275\le N<20{,}000\), then
\(p>2H\) and \(p^4>M\). The latter follows in each row from
\((3G)^4>(5/6)(G/2)N(N-1)\), by checking the upper row endpoint.
The [terminal-prime test](ORDINARY_BINOMIAL_COMMON_SLOPE_BOUND.md#a-prime-divisor-of-either-terminal-order),
equation (6), would require \(p<3H\), a contradiction.

**Full cosets at the lower terminals.** We now take \(N\ge20{,}000\)
within the prime-gap table. Put \(U=\log N\), \(\sigma=2U+2\).
Write \(H_{\max}=G-1\) and \(D=G/2-1\), so \(H\le H_{\max}\),
\(d\le D\). The three applicable rows have
\[
 (N_{\min},H_{\max},D)=(20{,}000,515,257),\
 (32{,}768,601,300),\ (65{,}536,1131,565).
 \tag{T34}
\]
In particular \(d/N<0.015\).

Suppose \(q=d/k_0\ge2\). Equation (T19d) applies, with
prime support at most \(P=M^{1/(q+1)}\).
Use the estimates leading to (T23), but take
\(\log P+0.2\) for the prime sum by the bound in (T27).
We can take \(\log P-0.49\) if \(P\ge350\).
The digit loss is at most \(0.151\) for \(q\ge3\),
from the same decreasing bound in \(q\) as in (T23).
At \(q=2\) it is less than \(4.848(d/N)^{1/3}<1.196\),
or less than \(0.659\) when \(d<50\). In the case \(d=2\)
this last estimate is less than \(0.226\).

The baseline \(me^{-1.01-g}\) following (T23) applies
since \(m=N-1\ge19{,}999\), and \(\log m>U-0.001\).
If \(k_0\ge2\) and \(\alpha\ge0.03\), use the three
rows and left-third argument from (T23). Here for \(q=2\)
the digit loss is less than \(0.191\) by
\(b\le\sigma/(0.03)^2\) and \(P/N=(b/N^2)^{1/3}\).
The three log gain bounds exceed
\[
 (U+\log(\sigma/\alpha_-^2))/3+0.391
\]
by more than \(0.95,0.41,0.50\), respectively, at
\(U=\log20{,}000\). The margins increase thereafter.

For \(\alpha<0.03\) and \(k_0\ge2\), the lower bounds
following (T24) still apply with \(t=\log(1/\alpha)\).
For \(q=2,\ d\ge50\), the bound \(0.955U+0.96t-2.92\)
beats \((U+2t+\log\sigma)/3+1.396\) at \(U=\log20{,}000\),
\(t=3.5\), and increases more quickly in both variables.
If \(q=2,\ d<50\), the baseline \(U-1.22\) exceeds
\((2U+\log50)/3+0.169\) for \(P\ge350\), and exceeds
\(\log350+0.859\) otherwise. If \(q\ge3\), the lower bound
\(0.829U+t/3-0.86\) beats \((U+2t+\log\sigma)/4+0.351\)
using \(t<U\).

If \(k_0=1\) and \(d\ge3\), then \(d<U\) by the same
\(d\log(4d)<2U+\log d\) comparison as above.
The lower bound \(0.835U-1.011\) exceeds
\((2U+\log U)/4+0.351\). If \(k_0=1,\ d=2\),
the lower bound \(U-1.34\) exceeds \((2U+\log2)/3+0.426\).
Thus \(k_0=k=d\) throughout the range (T34).

**Moderate and large \(\alpha\).** Suppose \(\alpha\ge0.03\).
Then \((1+\alpha)^2 b\le(1+\alpha)^2D/\alpha<N\),
since \((1+\alpha)^2/\alpha\) decreases for \(\alpha<1\)
and \(D/N<0.015\). In (T7) and (T10) the prime support is
at most \(X=\sqrt M>350\).
Use \(V/\log X<2.001\) and \(X\ge\sqrt{7N}\).
The right side of (T14) is less than
\[
 \tfrac12(U+\log b)-0.49+2.46\sqrt{b/N}.
 \tag{T35}
\]
For a rotation with \(\cos\theta\le1-\alpha\), we have
\(\Re z\le am\), so (T18a) gives gain at least \(m e^{-1.01}\).
Among the equally spaced rotations the fraction violating
the cosine condition is at most \(\arccos(1-\alpha)/\pi+1/d\).
Combining this with \(s/d=\alpha\) improves the log gain
to \(\log m-1.01-\kappa(\beta s+\alpha)\),
where \(g\le\kappa s\) and \(\beta\ge\arccos(1-\alpha)/\pi\).
If \(d\le9\), use \(\kappa=7/32,\ \beta=0.42,\ s\le6.75\).
The lower bound exceeds \(U-1.80\), while \(b\le300\)
makes (T35) smaller already at \(N=20{,}000\).

For \(d\ge10\) the left-half-plane fraction is at least
\(0.45\) and belongs to the cosine range just considered.
We can add \(0.45\max(0,1.01+\log((1-\alpha)/\alpha))\)
to the log gain bound by using (T18) there.
On an interval \([\ell,u]\) of values of \(\alpha\), take
\(S_0=\min(\sigma,Du)\) and \(b_0=\min(\sigma/\ell^2,D/\ell)\).
The following table gives certified upper bounds
for \(\kappa,\beta\) and a lower bound on the surplus at each
of the three starting \(N_{\min}\) in (T34). Use
\[
 U-1.011-\kappa(\beta S_0+u)
       +0.45\max(0,1.01+\log((1-u)/u))
 -\left(\tfrac12(U+\log b_0)-0.49+2.46\sqrt{b_0/N}\right).
\]
The \(\kappa\) values follow from (T18a) at the upper endpoints.
\[
 \begin{array}{c|c|c|c}
 [\ell,u]&\kappa&\beta&\text{minimum surplus exceeds}\\\hline
 [0.03,0.036]&.128&.086&.12\\
 [0.036,0.045]&.129&.096&.20\\
 [0.045,0.07]&.131&.120&.11\\
 [0.07,0.11]&.134&.151&.21\\
 [0.11,0.17]&.140&.189&.51\\
 [0.17,0.28]&.150&.245&.67\\
 [0.28,0.45]&.170&.315&.62\\
 [0.45,0.60]&.190&.370&.53\\
 [0.60,0.75]&.218&.420&.08
 \end{array}
\]
All of these differences increase in each band of \(N\).
The affine gain grows at rate at least \(1-2\kappa\beta\).
The log term in (T35) grows at rate at most
\(1/2+1/\sigma\), and its square-root term decreases.
This excludes \(\alpha\ge0.03\).

**Small \(\alpha\) and \(b\le N\).** Use \(X=1.03\sqrt M\) for
the support, by (T7) and (T10). If \(\log X<8\), the right
side of (T14) is below \(7.81\) using the prime sum at \(e^8\),
\(V<16\), \(N\ge20{,}000\) and the prime count.
Otherwise use the prime sum at \(X\), and
\(\pi(X)V/N<2.392\) from \(V/\log X<2.001\), \(X/N\le1.03\).
This gives an upper bound
\[
 \tfrac12(U+t+\log d)+1.89
 \tag{T36}
\]
on the right side, since \(\log b=t+\log d\).
For \(d\ge20\), (T24a) gives \(\log\mathcal L>U+0.95t-A_D\),
where
\[
 A_D=1.418+0.111(1.01+0.00384D)+0.05434.
\]
Use \(g_0/d\le0.05434\), \(g_0\le1.01+0.00384D\),
and \(\log C>t-1.417\) in that inequality.
For each band in (T34) the surplus over (T36)
at \(N=N_{\min},t=3.5,d=D\) exceeds \(0.16\).
For \(10\le d\le19\), the same argument with \(g_0<1.083\)
gives \(\log\mathcal L>U+0.9t-1.65\), which also
exceeds (T36) at \(N=20{,}000,\ t=3.5,\ d=19\).
If \(d\le9\), use \(g_0<1.045\) and the left-half-plane
fraction. It is at least \(0.4\) except at \(d=3\),
where it is at least \(1/3\).
The log gains exceed \(U+0.4t-0.641\) and \(U+t/3-0.71\),
respectively. They exceed (T36) by \(t\le U-\log d\)
and \(U\ge\log20{,}000\).
All these lower bounds also exceed \(7.81\), which
handles the smaller \(\log X\) case.

**Small \(\alpha\) and \(b>N\).** Use \(V_0=2U+\log(D+1)\)
and \(x=H_{\max}V_0\ge\log K\).
The right side of (T14) has prime sum less than
\(\log x+0.461\) by splitting at \(x\), as in (T24b).
For the digit term, for \(y>1\) we have
\[
 \omega(K)\le\frac{\log K}{\log y}+\pi(y)-\frac{\vartheta(y)}{\log y}.
\]
Each norm prime above \(y\) contributes at most its log
divided by \(\log y\). Below \(y\) the excess is at most
\(1-\log p/\log y\).
At \(y=2048\), the prime count gives \(\pi(y)<314\)
and Rosser and Schoenfeld's Theorem 10 gives \(\vartheta(y)>0.96y\).
Hence the total upper bound in (T14) is
\[
 \log(H_{\max}V_0)+0.461+
    e^{-U}\left(60V_0+\frac{H_{\max}V_0^2}{\log2048}\right).
 \tag{T37}
\]
It grows at rate at most \(2/V_0\) within each band,
because its last two terms decrease with \(U\).

Here \(\alpha<D/N<0.013\) and \(s<D^2/N<4.88\),
so \(g_0<1.635\).
Equation (T24a) gives \(\log C>t-1.40\) and angular
fraction \(r/d<0.073\), with \(t>U-\log D\).
If \(d\ge50\), it gives \(\log\mathcal L>U+0.98t-1.56\).
If \(2\le d<50\), then \(g_0<1.026\).
Use \(1-1/d\) at least \(1/2,0.8,0.9\)
on \(2\le d\le4,\ 5\le d\le9,\ 10\le d\le49\),
respectively, in the coefficient of \(t-1.40\).
In each range, use \(t>U-\log d\) and the lower endpoint
of the range for \(g_0/d\).
All three gains exceed \(1.5U-2\) for \(U\ge\log20{,}000\).

At the three starting \(N_{\min}\) in (T34),
the differences between (T37) and the respective smaller
gain bounds are strict: the gains exceed (T37) by
more than \(0.42,\ 1.55,\ 1.47\).
The lower bounds have slopes at least \(1.5\), so this
excludes the last case throughout the bands.
Thus every consecutive pair with terminal \(N\ge20{,}000\)
is coprime over \(\mathbb Q\). Rational Bézout transfers the
result to every characteristic-zero field.

### Sharper Padé estimates at smaller terminals

Continue to write \(U=\log N\), \(m=N-1\), and \(t=\log(1/\alpha)\)
for \(\alpha=d/b\). By Stirling, the baseline gain from (T18a)
exceeds \(me^{-1.01-g}\) for \(N\ge200\). In fact the needed
inequality \(0.01(m+1)+\tfrac12\log(2\pi/m)>0\) holds at \(m=199\)
and its left side increases. Also \(\log m>U-0.006\).

Write \(S(x)=\sum_{p\le x}\log p/(p-1)\). For each pair \((x_0,c)\) below, \(S(x)<\log x-c\) holds for all \(x\ge x_0\).
\[
 (x_0,c)=(21,.22),(34,.3),(43,.33),(117,.425),
          (137,.449),(343,.49),(100,.39),(150,.45).
\]
We will also use \(\pi(x)\log x/x<1.224\) for \(x\ge117\).
Both bounds below \(3000\) follow by checking the left endpoints
and the primes, using the rational logarithm enclosures in (T27).
For the count it suffices to check through \(300\), since
\(\log x/x\) decreases for \(x>e\) and the Dusart bound preceding
(T20) then gives \(1+1.2762/\log x<1.224\). For the sum above
\(3000\) the bound preceding (T20) implies all the displayed
inequalities. We use \(\pi(x)\log x/x<1.256\) for \(x>1\) as well.
[Rosser and Schoenfeld, Corollary 1, equation (3.6)](../../arithmetic/sources/rosser_schoenfeld_1962.pdf).

**A product bound for small angles.**

Suppose \(\alpha<0.03\) and put \(g_0=1.01+0.128s\) and
\(C'=(1-\alpha)/(2\alpha)\). Let \(n\) be the size of a
Kummer coset and order its angles by absolute value,
starting at \(\theta_0\). Put \(u_j=|1-e^{i\theta_j}|^2\).
At every point with \(u_j\ge4\alpha\), (T18) gives gain at
least \(mC'u_j(1-2\alpha/u_j)\). If the cosine is nonnegative,
use \(|z|\le bm\) in \(am-\Re z\). If the cosine is negative,
use \(|z|>am\). Use the baseline \(me^{-g_0}\) at \(j=0\)
and at all other points with \(u_j<4\alpha\).

Write \(r\) for the number of these other points. The identity
\(\prod_{j=0}^{n-1}u_j=|1-e^{in\theta_0}|^2\) gives
\(\prod_{j>0}u_j\ge(2n/\pi)^2\), also when \(\theta_0=0\)
by continuity. This uses \(|\theta_0|\le\pi/n\) and
\(\sin v\ge2v/\pi\) on \([0,\pi/2]\). The removed factors
\(C'u_j\) are less than \(2(1-\alpha)\).
For the retained factors, \(\log(1-y)\ge-2(\log2)y\)
for \(0\le y\le1/2\). The angle bound \(|\theta_j|\ge\pi j/n\)
gives \(r/n\le(2/\pi)\arcsin\sqrt\alpha<.111\). It also bounds
the correction from those logarithms by the integral of the
decreasing function \(2(\log2)\min\{1/2,\alpha/(2\sin^2(\pi x/2))\}\)
over \(0\le x\le1\). Thus the geometric mean \(\mathcal L\)
of gains in (T14) satisfies
\[
 \log\mathcal L\ge
 \log m+(1-1/n)\log C'
 +\frac{2\log(2n/\pi)-g_0}{n}
 -\frac{r}{n}\bigl(g_0+\log(2(1-\alpha))\bigr)-A,\qquad
 A=\frac{2\log2}{\pi}\bigl(\arcsin\sqrt\alpha+\sqrt{\alpha(1-\alpha)}\bigr)
 <.153.
 \tag{T38}
\]
This applies uniformly over all cosets of size \(n\).
For \(2\le n\le9\), \(r=0\) because
\(\sin^2(\pi/(2n))>.03\). In that range one can use
\(A=2(\log2)\alpha(n^2-1)/(3n)\), from
\(\sum_{j=1}^{n-1}\csc^2(\pi j/(2n))=2(n^2-1)/3\).

**Proper Kummer cosets.**

Assume \(N<20{,}000\) so \(d\le233\) by (T33). If (T10)
applies and \(q=d/k_0\ge2\), (T19d) gives \(4d\le2H<P\),
where \(P=M^{1/(q+1)}\) bounds all norm primes. For any
\(X\ge P\), the upper bound in (T14) is at most
\[
  S(X)+\frac{\pi(X)\log(M+1)}{N}.
  \tag{T39}
\]
If \(X^j\ge M\) in the applications below, then
\(\log(M+1)/\log X<j+.001\) since \(M\ge7N\) and \(X\ge9\).

For \(\alpha\ge.03\), \(M=dN/\alpha<PN/(4\alpha)\)
gives \(P^q<N/(4\alpha)\). On an interval \([\ell,u]\)
take \(X=\sqrt{N/(4\ell)}\). Then \(X>P\), \(X^3>M\),
and \(s=d\alpha<\sqrt{Nu/64}\).
The differences between the baseline from (T18a) and
(T39) exceed
\[
 \tfrac12\log(4\ell N)-1.016-\kappa\sqrt{Nu/64}
       -c-\frac{3.001\cdot1.256}{\sqrt{4\ell N}},
\]
using \(g\le\kappa s\) and \(S(X)\le\log X+c\).
This is concave in \(U\), so the endpoint margins suffice.
The first two rows use (T10) from \(a\ge2d\). The last row
uses the extension of (T10) proved above.
\[
 \begin{array}{c|r|c|r|rr}
 [\ell,u]&N_{\min}&\kappa&c&
  \text{margin at }N_{\min}&\text{margin at }20{,}000\\\hline
 [.03,.1]&200&.134&-.30&>.028&>2.34\\
 [.1,1/3]&200&.160&.20&>.39&>1.60\\
 [1/3,.75]&275&.218&.20&>1.14&>.51
 \end{array}
\]
Here \(\kappa\) follows from (T18a) at the upper endpoints.

For \(\alpha<.03\), \(a>32d\) makes (T10) applicable.
If \(q\ge3\), \(M<dN^2<PN^2/4\) gives
\(P<X=N^{2/3}/4^{1/3}\), \(M<X^4\), and
\(d<N^{2/3}/4^{4/3}\). Use \(S(X)<\log X-.22\) and
the Dusart count as in (T20). The surplus of the baseline is
greater than
\[
 U-1.016-.00384D-\log X+.22
       -4.001(X/N)(1+1.2762/\log X).
\]
Use \(D=26\) on \(200\le N<2000\) and \(D=233\) on
\(2000\le N<20{,}000\). The surplus increases on both
intervals and exceeds \(0.72\) and \(1.04\) at their left endpoints.

If \(q=2\), the same inequality \(M<PN^2/4\) gives
\(P<X=N/2\) and \(M<X^3\). Formula (T39) is below
\(U+.862\) using \(S(X)<\log X-.33\) and the global count.
For \(d=4\) the coset of size two gives gain at least \(m(1-\alpha)/\alpha\)
at one point in the left half-plane. Its average log gain with
the baseline exceeds \(U+1.21\). For \(6\le d\le18\), use (T38)
with \(n=d/2\) and the small-\(n\) bound on \(A\). With
\(t>\log(100/3)\), \(\log C'>t-.7237\) and \(g_0<1.01+.00384d\),
the resulting gain exceeds \(U+1.89\) for each \(n=3,\ldots,9\).
For \(d\ge20\) use \(1-1/n\ge.9\), \(g_0\le1.90472\),
\(g_0/n\le.10868\) and discard \(2\log(2n/\pi)/n>0\).
The gain exceeds \(U+1.94\).
If \(d=q=2\), use \(X=(2N^2)^{1/3}\) instead. The
surplus of the baseline over (T39) exceeds
\[
 U-1.016-.00384\cdot2-\log X+.33-3.001\cdot1.256\,X/N.
\]
This exceeds \(0.029\) at \(N=200\) and increases thereafter.
Thus \(k_0=k=d\) is necessary at every \(N\ge275\). For
\(\alpha\le1/3\) this holds already at every \(N\ge200\).

**Small \(\alpha\) and \(b\le N\).**

Take \(N\ge200\), \(\alpha<.03\), \(b\le N\), and a full
coset of size \(n=d\). By (T7) and (T10), we can take
\(X=1.03\sqrt M\) in (T39), so \(X/N\le1.03\).
For \(d\ge20\), (T38) gives \(\log\mathcal L>U+2.14\),
using \(g_0/d\le.05434\), \(g_0\le1.90472\),
\(\log m>U-.002\) since \(N\ge b>d/.03\), and dropping the
positive \(2\log(2d/\pi)/d\). For \(10\le d\le19\),
(T38) gives \(\log\mathcal L>U+2.41\) by substituting each
\(d\), with \(g_0\le1.01+.00384d\). In both ranges,
\(X>343\) and \(V/\log X<2.001\), so (T39) is below
\[
 U+\log(1.03)-.49+2.001\cdot1.03(1+1.2762/\log343)
 <U+2.052.
\]
For \(4\le d\le9\), the small-\(n\) formula in (T38) gives
gain greater than \(U+2.24\) for each \(d\). Here \(X>137\),
and (T39) is below \(U+2.177\), replacing \(.49\) and \(343\)
in the displayed bound by \(.449\) and \(137\).

For \(d=3\), write \(\mu=\tfrac12\log(b/N)\). Formula (T38)
gives gain greater than \(U+2t/3-.435\). The upper bound
in (T39) is less than \(U+\mu-.419+2.523e^\mu\) since
\(X>137\). With \(t=U-\log3+2\mu\), the surplus is concave
in \(\mu\in[(\log100-U)/2,0]\).
At \(\mu=0\) it is greater than \(0.26\) for \(N\ge200\).
At the other endpoint the lower gain minus \(U\) is greater
than \(1.90\) while the upper bound minus \(U\) is less than
\(1.02\) and decreases with \(N\).

For \(d=2\) we can take \(X=1.015\sqrt M>117\) since \(b\ge67\)
and \(a+2d-1=b+1\). At antipodal angles put \(c=|\cos\theta|\),
\(A_0=(1-\alpha)/\alpha\) and \(B_0=e^{-g_0}\).
The two gains have product at least
\(m^2 A_0(1+c)\max(B_0,A_0-c/\alpha)\) by (T18).
This expression is concave on \(0\le c\le1-\alpha(1+B_0)\)
and increasing on the remaining interval. The smaller
endpoint gives the lower bound \(m^2 A_0B_0(2-\alpha(1+B_0))\).
Consequently \(\log\mathcal L>U+t/2-.2\). Formula (T39)
is less than \(U+\mu+\log(1.015)-.425+2.001\cdot1.015\cdot1.224\,e^\mu\).
Since \(t=U-\log2+2\mu\), the surplus is minimized at \(\mu=0\),
and exceeds \(0.026\) already at \(N=200\).
This excludes \(\alpha<.03,\ b\le N\) from every possible
common factor at \(N\ge200\).



**Local refinement of the gap bound.**

For any record gap \(G\ge52\), the second part of (T33)
also applies when the actual upper endpoint \(aN+H\) is
below \(X_G\), and an available bound \(H\le H_0\) gives
\(H_0/G<N/2\). In fact \(r=\lfloor H/G\rfloor\), if positive,
then satisfies \(r<N/2\) and \((aN+H)/r<X_G\). The same
two cases for \(C\) give \(p>C\) and contradict (T7).
In particular, once (T10) applies, a bound \(aN+H_0<X_G\)
under these hypotheses gives
\[
 H<G,\qquad d<G/2 .
 \tag{T40}
\]

**Moderate and large \(\alpha\) in full cosets.**

Take \(2200\le N<20{,}000\) and \(\alpha\ge.03\).
The global bounds \(d\le233\), \(H\le467\) and
\(a\le(97/3)d\) give \(aN+H<189{,}695{,}659=X_{222}\).
Thus (T40) gives \(H\le221\) and \(d\le110\).
For \(N\le3500\) these bounds in turn give
\(aN+H<17{,}051{,}707=X_{154}\), so \(d\le76\).

The beta tail factor in (T18a) gives an additional gain at
rotations away from the positive axis. Fix \(c\ge\cos\theta\)
and \(u\ge\alpha\) with \(u\le.75\). For \(c\ge0\) put
\(W=(1-u-c)/u\). Make the substitution \(v=1-e^{-x/d}\)
in the second integrand of (T18a). Its logarithm divided by
\(m\) is bounded by
\[
 \frac{\Re z}{m}(1-e^{-x/d})-\frac{a}{d}x+\log(1-e^{-x})
 \le -Wx+\log(1-e^{-x}).
\]
Here \(\Re z/m\le bc\), \(d(1-e^{-x/d})\le x\),
and \((1-\alpha-c)/\alpha\ge W\).
For \(W>0\), the right side is at most \(-J(W)\), where
\[
 J(W)=(1+W)\log(1+W)-W\log W.
\]
Consequently the log gain exceeds
\(\log m+\max\{-1.01+J(W),\log W\}\).
The second term uses (T18) directly. For the first,
the bound \(J(W)<5\), which follows from
\(J(W)\le1+\log(1+W)\) and \(W\le 97/3\),
absorbs the difference between \(m\) and \(m+1\) in
the Stirling estimate at \(m\ge2199\).
If \(c<0\), the same argument uses \(W=(1-u)/u\) in the
beta integral and \(\log(W(1-c))\) in (T18).

At \(c\ge0,\ W\le0\), use the baseline
\(\log m-1.01-\kappa(u)du\), where \(\kappa(u)=g(u)/s\)
in (T18a) evaluated at \(\alpha=u\).
For the smaller terminal band we can also drop
the beta tail factor and maximize the remaining integrand.
The log gain then exceeds
\[
 \log m-1.01-\frac d u
 \left[c-(1-u)-(1-u)\log\frac{c}{1-u}\right].
\]
Indeed the loss at \(\alpha\) is zero for \(c\le1-\alpha\),
and otherwise is the displayed expression with \(u=\alpha\).
It increases with \(\alpha\), since the derivative has the sign of
\(\log(c/(1-\alpha))-c+(1-\alpha)\ge0\).
Use the better of this bound and the baseline. Denote the resulting
lower bound minus \(\log m\) by \(\Phi_u(c)\).

For an interval \([\ell,u]\) of values of \(\alpha\) and fixed \(d\),
take \(c_j=\cos(\pi j/d)\) for \(0\le j<d\), using the angle
ordering in (T38). This gives
\(\log\mathcal L>\log(N-1)+d^{-1}\sum_{j=0}^{d-1}\Phi_u(c_j)\).
By (T7) and (T10), a valid prime threshold is
\(X=\max(a+2d-1,\sqrt M)\). Put
\[
 X_0=\max\{d(1+1/\ell),\sqrt{dN/\ell}\},\qquad
 X_1=\sqrt{\max(7,d/u)N}.
\]
Then \(X_1\le X\le X_0\) and \(V/\log X<2.001\).
The right side of (T39) is less than
\[
 \log X_0-c_0+2.001(1+1.2762/\log X_1)X_0/N .
 \tag{T41}
\]
Use \(c_0=.39\) for \(2200\le N\le3500\), where
\(X_1>100\), and \(c_0=.45\) for \(3500\le N<20{,}000\),
where \(X_1>150\). The term \(\log X_0\) grows with slope
at most \(1/2\) in \(U\), and the final term decreases.
The lower bound on the gain grows with slope greater than one.

It suffices to compare at \(N_0=\max\{N_{\min},
\exp((d\ell-2)/2)\}\), since \(d\ell\le s<2\log N+2\).
For the first band use intervals of width \(0.005\) from
\(0.03\) to \(0.05\) and width \(0.01\) thereafter to \(0.75\).
For the second band use width \(0.02\) throughout.
Intervals with \(d(1/\ell-1)<5\) have no eligible support.
The directed bounds give surplus more than \(0.070\) and
\(0.089\), respectively. The [endpoint verifier](verify_ordinary_binomial_pade_intervals.py)
checks these comparisons and the small-argument prime estimates.
Thus full cosets with \(\alpha\ge.03\) are excluded from
\(N=2200\) onward.



### Coprimality from terminal 3150

The remaining case for \(3150\le N<20{,}000\) is a full
coset with \(\alpha<.03\) and \(b>N\).
The strengthened Gamma annulus gives \(b<(5/6)(N-1)d\).
In this range the row \(G=354\) of (T33) applies through
\(N=6000\) by its uniform endpoint inequality.
To improve the caps, successively apply (T40) with
\((G,N_{\max})=(320,5000),(292,4150),(288,3467),(282,3295)\).
At each step the previous caps \(d\le D_{\rm old}\), \(H\le H_{\rm old}\)
give \(aN+H\le(5/6)D_{\rm old}N_{\max}^2+H_{\rm old}<X_G\).
The respective \(X_G\) are \(3{,}842{,}610{,}773\),
\(2{,}300{,}942{,}549\), \(1{,}453{,}168{,}141\) and
\(1{,}294{,}268{,}491\). This gives the first four rows below.
The last three rows follow directly from (T33).

| \(N\) | \(G\) | \(D=G/2-1\) | \(A_v\) | \(A_c\) | Lower bound for surplus |
|---:|---:|---:|---:|---:|---:|
| 3150–3295 | 282 | 140 | 2.504 | 2.885 | .023 |
| 3296–3467 | 288 | 143 | 2.506 | 2.886 | .146 |
| 3468–4150 | 292 | 145 | 2.514 | 2.889 | .31 |
| 4151–5000 | 320 | 159 | 2.522 | 2.892 | .44 |
| 5001–6000 | 354 | 176 | 2.531 | 2.896 | .86 |
| 6001–10000 | 456 | 227 | 2.58 | 3.01 | .59 |
| 10001–19999 | 468 | 233 | 2.58 | 2.92 | 2.37 |

**An upper bound for the norm charge.**

Let \(w=\log(b/N)\), \(v=\log M=2U+w\),
\(t_D=U+w-\log D\), \(\xi=t_D+.337\),
\(P=4M^{1/3}\), and \(z=\log P\). Then \(0\le w\le\log(5D/6)\).
In each row \(P>830\) and \(P>5G/2>2H+d\).
At any \(p>P\), we have \(p^3>M\). The affine Newton
error is at most two, so the argument of (T19a) now
applies without a condition on \(H^4/M\). More precisely,
all maximizing indices for a selected positive slope lie below
\(4h<p\). The factorial jumps are spaced more than \(2h\)
indices apart with heights at most two. The binomial valuation is
nondecreasing in this range, so any selected positive
slope joins zero to the first jump at each order.

If the two edge heights are equal, their widths are equal.
By the short binomial product from (T19a), these primes
have total logarithmic weight at most
\[
 2d\bigl(1+\log(a/d+2)-\log 2\bigr)< 2d(t+.337)\le 2D\xi .
\]
The last inequality uses that \(x(\log(b/x)+.337)\) increases
for \(x\le D<b\), since \(D\le233\) and \(b>N\).
Their contribution to (T14) is at most
\(2D\xi(1/(P-1)+V/(Nz))\).

If the heights differ, at least one of the two reached terminal
integers is a multiple of \(p^2\). The smaller width is at most
\(h\), and the other is twice it. Hence a window
\([bm-2H,bm]\), \(m=N-1,N\), contains \(qp^2\), with
positive integer \(q\le Q=M/P^2=P/64\). For each \(q\) and
window there is at most one such \(p>P\), since the spacing of
its squares exceeds \(2H\). Their digit charge is at most
\(2QV/N\).
Put \(Z=b(N-1)-2H>.9995M\).
Since \(\sqrt{Z/q}\le p\le\sqrt M\), their prime sum is bounded by
\[
 v\,\frac{\sum_{q\le Q}\sqrt q}{\sqrt Z-\sqrt Q}
 \le \frac{1.002v}{64}\left(\frac23+\frac{1}{2Q}\right)
 <\frac{.71v}{64}.
\]
The middle bound uses
\(\sqrt Z-\sqrt Q\ge\sqrt M(.9995-1/P)>\sqrt M/1.002\)
and \(\sum_{q=1}^n\sqrt q\le 2n^{3/2}/3+\sqrt n/2\)
by concavity and the trapezoid sum. Here \(Q>12\).

At \(p\le P\), use \(S(P)<\log P-.49\) and the Dusart
prime count. Since \(V<v+.001\) and \(v=3(z-\log4)\),
the ratio \(V/z\) is bounded above by
\(3-(3\log4-.001)/z\).
Both this expression and its product with \(1+1.2762/z\)
increase in \(z\). At the upper corner of each row they are
less than \(A_v,A_c\), respectively. Thus the right side of
(T14) is less than
\[
 R=z-.49+\frac PN\left(A_c+\frac{2(v+.001)}{64}\right)
 +\frac{.71v}{64}+\frac{2.004D\xi}{P}+\frac{2A_vD\xi}{N}.
 \tag{T42}
\]

**The angular gain and endpoint comparison.**

Use (T38) and \(\log m>U-.001\).
Its error terms satisfy \(A<.883\sqrt\alpha\) and
\(r/d<.641\sqrt\alpha\).
For the first use
\(\arcsin x+x\sqrt{1-x^2}\le2x\).
For the second use \(\arcsin x\le1.006x\) at \(x^2\le.03\),
from integrating
\((1-y^2)^{-1/2}\le1+y^2/(2(1-x^2)^{3/2})\) on \([0,x]\).
If \(d\ge20\), put \(E=D/b\ge\alpha\).
Then \(\log C'>t-.6932-1.031E\),
\(g_0\le1.01+.128DE\), and \(g_0/d\le1.01/d+.128E\).
For fixed \(b\), the functions
\((1-1/d)(\log(b/d)-.6932)\) and
\([2\log(2d/\pi)-1.01]/d\) decrease from \(d=20\) through
\(D\) as used here. For the first, \(t<U<10\) at the starting
\(d\) and \(t\) decreases as \(d\) increases.
Consequently the log gain is greater than
\[
 L=U-.001+(1-1/D)(t_D-.6932)
     +\frac{2\log(2D/\pi)-1.01}{D}
     -1.159E-\bigl[.641(1.7032+.128DE)+.883\bigr]\sqrt E.
 \tag{T43}
\]
For smaller \(d\), use \(D=2,3,9,19\) on the respective
ranges \(d=2,\ d=3,\ 4\le d\le9,\ 10\le d\le19\)
in (T42) and in \(t_D\). From (T38) and the same antipodal
bound as for \(b\le N\), the log gain exceeds \(L=U+f t_D-c\)
with \((f,c)=(1/2,.2),(2/3,.435),(.75,.68),(.9,.925)\),
respectively. For \(4\le d\le9\) use \(r=0\) and the
small-coset error in (T38). For \(10\le d\le19\) use
\(g_0<1.083,\ r/d<.112,\ A<.153\), and
\(2\log(2d/\pi)/d\ge2\log(38/\pi)/19\).
The same \(A_v,A_c\) apply with the smaller \(D\).

For fixed \(N\), \(L\) is concave in \(w\) since its nonaffine
terms are negatives of decreasing exponential functions.
The bound \(R\) is convex in \(w\). The only potentially
negative second derivative term is from \(2.004D\xi/P\).
Together with the derivative from the \(P/N\) term it is
nonnegative if
\[
 (P^2/N)\left(A_c+\frac{2(v+.001+6)}{64}\right)
 \ge 2.004D(6-\xi).
\]
It suffices to check at the lower \(N\) endpoint, \(w=0\),
and the largest \(D\) in each row. The left side increases
with \(U,w\). The right side decreases, and a smaller \(D\)
can only decrease its positive part. The \(U\) derivative
of \(L\) is at least \(1.5\), while that of \(R\) is at most
\(2/3+2\cdot.71/64<.69\).
In fact the \(P/N\) term decreases with \(U\), and so do
the two equal-height terms since \(\xi>3\).

The surplus \(L-R\) exceeds the margins in the table at
the lower \(N\) endpoint and both endpoints \(w=0,\log(5D/6)\).
For the four smaller \(D\) ranges it exceeds \(3.8\).
The [exact-rational check](verify_ordinary_binomial_hard_pade_endpoints.py)
certifies these margins, the convexity bounds, the gap endpoints
and the constants for the small-\(d\) gains. Concavity in \(w\)
and growth in \(U\) imply \(L>R\) throughout each row,
contrary to (T14).

The proper-coset exclusion and the two other angular bounds
cover every remaining branch through \(N=19{,}999\).
The earlier prime-gap cutoff covers \(N\ge20{,}000\).
Therefore every ordinary binomial consecutive pair
with terminal \(N\ge3150\) is coprime over \(\mathbb Q\)
and over every characteristic-zero field by rational Bézout.



## Lower exponent four

Every nonzero ordinary binomial \(X^4(A+BX^d)\), for every \(d\ge1\),
has a nonzero ordinary factorial moment in each consecutive positive
pair over every characteristic-zero field.

Adjacency covers \(d=1\) and the relative-gap theorem covers \(d\ge12\).
The [all-prime-gap theorem](STEPANOV_LOWER_FOUR_ELEVENTH_GAP_CLOSURE.md#7-all-prime-gaps-at-lower-exponent-four)
covers all prime gaps. The [dyadic theorem](DYADIC_BINOMIAL_SFC.tex)
covers \(d=4\), and the [gap-eight theorem](STEPANOV_LOWER_FOUR_EIGHTH_GAP_CLOSURE.md)
covers \(d=8\). These complements leave \(6,9,10\).
At all three gaps the [terminal-through-twelve theorem](STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md)
and the eventual theorem already covered the terminals outside \(13\le N<N_*\).
Within that prefix the [prime-power Hasse criterion](PRIME_POWER_HASSE_COPRIMALITY.md#lower-four-gap-nine-and-the-offset-twice-the-prime)
closes gap nine. The [gap-six proof](ORDINARY_LOWER_FOUR_GAP_SIX.md)
and [gap-ten proof](LOWER_FOUR_GAP_TEN_STRONG_FACTORIAL.md) now cover
the two other gaps at every order. Hence the whole lower-four
statement follows without an evaluation of \(N_*\).
