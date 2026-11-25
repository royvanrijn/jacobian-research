# A joint valuation bound for ordinary binomial common roots

Put \(b=a+d\), \(M=bN\), and let \(a,d,N\) be integers with \(a\ge1\),
\(d\ge2\), and \(N\ge3\). The two ordinary moment polynomials are
\[
 P_m(t)=\sum_{k=0}^m\binom mk(am+1)_{dk}t^k,\qquad m=N-1,N.
\]
Use their monic reciprocals \(R_m(x)=x^mP_m(1/x)\). The bounds below
use the relation between the two orders and hold at every prime. They
give a necessary norm inequality and a new prime support bound for a
common factor in the [ordinary residual region](ORDINARY_BINOMIAL_RESIDUAL_REGION.md).

## A bound for every common slope

For a prime \(p\), write \(L_p=\lfloor\log_p M\rfloor\),
\(\mu_p=d/(p-1)\), and \(E_p(J)=\lfloor\log_p(b+dJ)\rfloor\) for integers
\(J\ge0\). Every common reciprocal root, in \(\overline{\mathbb Q}_p\),
has valuation \(\lambda\) satisfying
\[
 \lambda\ \le\
 \mu_p+\max\left\{\frac{L_p}{J+1},\,3E_p(J)\right\}
 \qquad(J\ge0).                                                    \tag{1}
\]

To prove this, let \(c_{m,k}=[t^k]P_m\) and set
\[
 U_m(j)=v_p(c_{m,m})-v_p(c_{m,m-j})
       =v_p\!\left(\prod_{u=bm-dj+1}^{bm}u\right)-v_p\binom mj ,
       \qquad 0\le j\le m.
\]
An edge of slope \(\lambda\) in the lower Newton polygon of \(P_m\)
requires at least two minimizers of \(\lambda j-U_m(j)\). The value at
\(j=0\) is zero. Counting multiples of each \(p^s\) in the factorial
interval gives \(U_m(j)\le \mu_p j+L_p\). Thus, if
\(\lambda>\mu_p+L_p/(J+1)\), every minimizer is at \(j\le J\).

Write \(E=E_p(J)\). All factorial intervals under consideration lie in a
single interval
\[
 [\,b(N-1)-dJ+1,\ bN\,]\cap\mathbb Z_{>0}
\]
of length at most \(b+dJ<p^{E+1}\). There is at most one integer
\(u_*\) in it divisible by \(p^{E+1}\). Write \(V_m^E(j)\) for the
factorial valuation restricted to the powers \(p^s\) with \(s\le E\),
and \(B_m^E(j)\) for the corresponding truncated binomial valuation.
For \(j'>j\), counting the added multiples gives
\[
 V_m^E(j')-V_m^E(j)\le\mu_p(j'-j)+E,\qquad 0\le B_m^E(j)\le E.
\]
The binomial bound follows because each level in Legendre's formula
contributes either zero or one. Its remaining contribution
\(B_m^{>E}(j)\) is nondecreasing in \(j\le J\): the denominator \(j!\)
has no multiples of \(p^{E+1}\), so only the growing interval
\([m-j+1,m]\) contributes at these levels.

Suppose also \(\lambda>\mu_p+3E\). Then
\(\lambda j-V_m^E(j)+B_m^E(j)+B_m^{>E}(j)\) is strictly increasing
for \(0\le j\le\min(J,m)\). The only remaining contribution to
\(\lambda j-U_m(j)\) is a single downward jump of size \(v_p(u_*)-E\)
when the factorial interval first contains \(u_*\). If that occurs
at \(\alpha_m\), the only possible minimizing indices are \(0,\alpha_m\).
An edge therefore requires this jump and equality of those two values.
A common slope uses the same \(u_*\) at both orders, with
\[
 \alpha_m=1+\left\lfloor\frac{bm-u_*}{d}\right\rfloor,\qquad
 \lambda\alpha_m=U_m(\alpha_m).
\]
In particular \(\Delta=\alpha_N-\alpha_{N-1}\ge1\), since \(b>d\).

In the difference \(U_N(\alpha_N)-U_{N-1}(\alpha_{N-1})\), the factorial
contribution above level \(E\) cancels. The binomial numerator interval
\([N-\alpha_N+1,N]\) contains \([N-\alpha_{N-1},N-1]\). Hence the
binomial contribution above \(E\) can only decrease the difference.
The two factorial intervals intersect at \(u_*\). Their difference is
described by at most two added or removed intervals, of net length
\(d\Delta\). At each level \(p^s\), the count of multiples in an
interval differs from its length divided by \(p^s\) by at most one.
The truncated factorial difference is consequently at most
\(\mu_p\Delta+2E\). The truncated binomial difference can add at most
\(E\). Thus \(\lambda\Delta\le\mu_p\Delta+3E\), contradicting the
assumed slope. This proves (1), including \(E=0\).

## An exact norm and prime support test

Let \(F\in\mathbb Z[x]\) be a common monic factor of degree \(h>0\).
Its constant term is nonzero. Any positive valuation of one of its
roots is at least \(1/h\), since the reduced denominator of a Newton
slope divides the width of its edge. For a prime \(p>2dh+1\), take
\(J=2hL_p\) in (1). If \(p>b+2dhL_p\) then \(E_p(J)=0\) and the
right side of (1) is strictly less than \(1/h\). It follows that every
prime divisor of \(F(0)\) belongs to
\[
 \mathcal P_h=
 \{p\le M:\ p\text{ prime},\
             p\le2dh+1\ \text{or}\ p\le b+2dhL_p\}.                 \tag{2}
\]
Primes above \(M\) see only coefficient units. The same prime support
bound applies to \((a(N-1)+1)_{dh}\) and \((aN+1)_{dh}\).
Indeed, if a prime outside \(\mathcal P_h\) divided either block,
the corresponding \(R_m\bmod p\) would be \(x^{m-h+1}\) times a
polynomial of degree \(h-1\). The reduction of \(F\) is monic of
degree \(h\) with nonzero constant term, which precludes divisibility.
In the residual region one can intersect (2) with the prime support
in (T7).

Taking \(J=L_p\) in (1) gives the slope excess bound
\[
 A_p=\max\{1,3E_p(L_p)\}.
\]
If \(p>b\) one can additionally take
\(J=\lfloor(p-b-1)/d\rfloor\), where \(E_p(J)=0\). Define
\[
 C_p=\mu_p+
 \begin{cases}
 A_p & p\le b,\\[2pt]
 \min\left\{A_p,\dfrac{L_p}{1+\lfloor(p-b-1)/d\rfloor}\right\}
       & p>b.
 \end{cases}
\]
Then \(v_p(F(0))\le\lfloor h C_p\rfloor\). Combining this with the
[sharp Gamma-product inner annulus](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#2-every-positive-shape-a-sharper-product-annulus)
gives the necessary inequality
\[
 \bigl((aN+1)_d\bigr)^h
       <|F(0)|
       \le\prod_{p\in\mathcal P_h}p^{\lfloor h C_p\rfloor}.          \tag{3}
\]
All tests in (2)–(3) use integer or rational arithmetic. In particular,
the logarithm floors are comparisons of integer powers. To enumerate
\(\mathcal P_h\) it suffices to use primes at most
\(\max(2dh+1,b+2dh\lfloor\log_2 M\rfloor)\).
Equation (3) excludes a candidate degree without constructing the
order-\(N\) polynomials.

## A prime support test for lifts of smaller degree

Let \(F\) be irreducible over \(\mathbb Q\) of degree \(h\) and let
\(hk\) be the degree of an irreducible factor of \(F(-z^d)\).
If \(k<d\), every prime dividing \(F(0)\) or either initial
factorial block of length \(H=dh\) is at most
\[
 \max\{2H,\sqrt{a+2d-1},M^{1/4}\}.
 \tag{4}
\]
More precisely, if \(p>2H\), \(p^4>M\) and \(p\mid F(0)\), then
\[
 u=v_p(F(0))\in\{2,3\},\qquad
 p^u\mid(a+1)_{2d-1},\qquad d/\gcd(d,u)\mid k .
 \tag{5}
\]
Thus a prime of type \(u=2\) forces \(k=d/2\) and one of type
\(u=3\) forces \(k\in\{d/3,2d/3\}\). Both types cannot occur
for the same lift.

Fix \(p>2H\) with \(p^4>M\) and a positive valuation
\(\lambda=u/v\) of a root \(\xi\) of \(F\), in lowest terms.
If \(\lambda\ge2/h\), the strict count
\(U_m(j)<dj/(p-1)+3\le j/(2h)+3\) for \(j>0\) puts all
maximizers of \(U_m(j)-\lambda j\) below \(2h\).
In this range \(j!\) is a unit and the whole factorial interval
has length less than \(p\). The binomial valuation is nondecreasing,
so the only possible maximizing indices are \(0\) and the first
index \(w_m\) reaching a multiple \(q_m\) of \(p\). Both must be
maximizers, so the edge has height \(r_m=\lambda w_m\) at most
\(v_p(q_m)\le3\). Each order has at most one such slope.

Let \(E=\mathbb Q_p(\xi)\). Its ramification index is divisible
by \(v\) and at most \(h\). If \(\lambda<2/h\), this forces
\(u=1,\ v>h/2\) and the ramification index is \(v\).
If \(\lambda\ge2/h\), adjoin \(\pi\) with \(\pi^v=p\)
to an unramified extension with sufficiently large residue field.
Since the scaled polynomial \(R_m(\pi^u Z)/R_m(0)\) has
integral coefficients and reduction \(1+cZ^{w_m}\), \(c\ne0\),
all its residue roots are simple by \(w_m<p\).
Hensel lifting gives the roots of valuation \(\lambda\) in
the extension of ramification index \(v\). Thus \(E\) again
has ramification index \(v\), and the integer valuation of
\(\xi\) in \(E\) is \(u\).

After completing at the embedding that sends a root of \(F\)
to \(\xi\), the degree-\(k\) factor of \(Z^d+\xi\) has Newton
slope of reduced denominator \(d/\gcd(d,u)\). Its degree is
therefore a multiple of this denominator. If \(k<d\), this
excludes \(u=1\), so every positive valuation of \(F\) at \(p\)
is the unique slope \(\lambda\ge2/h\). Its numerator \(u\ge2\)
divides both edge heights \(r_m\le3\), hence \(r_m=u\)
and \(w_m=v\). The width of this slope for \(F\) is exactly
\(v\), so \(v_p(F(0))=u\). The common width places
\(q_N-q_{N-1}\) in \([a+1,a+2d-1]\), and both \(q_m\)
are divisible by \(p^u\). This proves (5). As every prime
dividing either initial block must divide \(F(0)\), the bound
(4) follows.

## A prime divisor of either terminal order

Let \(F\) be a nonconstant monic common factor of the actual ordinary
reciprocals, of degree \(h\ge2\), and put \(H=dh\), \(b=a+d\) and \(M=bN\).
Suppose a prime \(p>2H\) divides \(m\in\{N-1,N\}\) and \(p^4>M\).
Then a necessary condition for \(F\) is
\[
 p\nmid b,\qquad 3H-d\le p<3H,\qquad
 p^2\mid b(m/p)-1. \tag{6}
\]

To prove this, first \(R_m(x)\equiv x^m\pmod p\).
For \(1\le k<p\), the binomial coefficient \(\binom mk\) is
divisible by \(p\). For \(k\ge p\), the rising factorial of length
\(dk\) contains a multiple of \(p\). Hence all roots of \(F\)
have positive \(p\)-adic valuation.

Use the reverse ordinates
\[
 U_m(j)=v_p((bm-dj+1)_{dj})-v_p\binom mj .
\]
For a slope \(\lambda\ge1/h\), the strict bound
\(U_m(j)<dj/(p-1)+3\) puts every maximizing index of
\(U_m(j)-\lambda j\) below \(6h\). In this range the factorial
multiples of \(p\) are \(bm-tp\), entered at indices
\(j_t=1+\lfloor tp/d\rfloor\). Only \(t\in\{0,1,2\}\) can occur.
For \(j<p\), the binomial valuation is \(v_p(m)\) at positive
indices. If the range reaches \(p\), then necessarily \(d=2\).
At the isolated index \(j=p\) the binomial valuation falls to
\(v_p(m)-1\), but
\[
 U_m(p)=v_p(b)+v_p(bm-p)+1\le4<p/h.
\]
Indeed if \(p\mid b\) then \(v_p(b)\le2\) and \(v_p(bm-p)=1\).
For \(p<j<6h<2p\), the binomial numerator also contains \(m-p\),
so the valuation is at least \(v_p(m)\) and is nondecreasing.
The isolated index \(p\) cannot maximize the weighted ordinate.
Away from this index the binomial valuation does not decrease,
including the comparison from \(j=p-1\) to \(j=p+1\) when both
indices occur. Thus all candidates are zero or factorial jump indices.

If \(p\mid b\), the first jump at \(j_0=1\) has height
\(v_p(b)>0\). The other factorial multiples, if present, are
simple because \(p^2\mid bm\) and \(0<t\le2<p\). Successive
jump indices have gaps at least \(\lfloor p/d\rfloor\ge2h\).
The weighted value at any subsequent jump is strictly smaller
than at the first. Thus the only slope \(\ge1/h\) has width
one, inconsistent with all \(h\) roots of \(F\) having positive
valuation.

If \(p\nmid b\), the first jump cancels against the binomial
valuation. The high-slope argument using \(j<2h\) therefore
excludes \(\lambda\ge2/h\). Every chosen slope must now have
numerator one and reduced denominator greater than \(h/2\).
All \(h\) roots of \(F\) have positive valuation, so the only
allocation uses a single edge of width \(h\), at slope \(1/h\).

At the next jump, \(U_m(j_1)\le3\). At the third jump,
\(U_m(j_2)\le4\): the two additional factorial multiples include
at most one multiple of \(p^2\), and the original terminal
multiplicity \(v_p(m)\) has canceled against the binomial term.
Also \(j_2=1+\lfloor2p/d\rfloor\ge4h+1\), so this third
jump cannot maximize \(U_m(j)-j/h\). The first jump gives a
negative weighted value. Consequently an edge at \(1/h\)
requires \(U_m(j_1)=j_1/h\), with \(j_1\ge2h+1\).
The only possibility is \(U_m(j_1)=3\) and \(j_1=3h\), so
\(p^3\mid bm-p\) and \(3H-d\le p<3H\), as asserted.

In particular the criterion excludes any candidate degree with
such a \(p\) dividing \(m\) to exponent at least two. If \(m\)
itself is prime and \(m>2H\), the residual inequality
\(b<d(N-1)\) gives \(b<m^2\) and \(M<m^4\).
The necessary congruence would give \(b\equiv1\pmod{m^2}\),
which is impossible for \(1<b<m^2\). In particular, if (T10)
applies and \(a+2d\le N\), the upper bound \(\max(a+2d-1,\sqrt M)\)
from (T10) is at most the prime terminal order \(m\).
Indeed \(a+2d-1\le N-1\) and \(M\le N(N-2)<(N-1)^2\).
Thus ordinary coprimality holds throughout this part of (T5).

**Arbitrary valuation depths.** There is a separate test that
does not require \(p^4>M\). If \(p>2H\) divides a terminal
order \(m\), put \(r=v_p(m)\) and \(L_p=\lfloor\log_p M\rfloor\).
The condition
\[
 p^r\ge 2H L_p
 \tag{7}
\]
excludes the candidate degree \(h\).
Indeed, for any candidate \(\lambda\ge1/h\), the strict counting
bound \(U_m(j)<dj/(p-1)+L_p\) at positive \(j\) puts all
maximizing indices at \(dj<2H L_p\le p^r\).
For \(1\le j\le m\) in this range, \(v_p(bm-i)=v_p(i)\)
for \(1\le i<dj\), and the binomial numerator similarly gives
\(v_p\binom mj=r-v_p(j)\).
Therefore
\[
 U_m(j)=v_p(b)+v_p((dj-1)!)+v_p(j).
\]
In particular \(U_m(1)=v_p(b)\).
For \(j\ge2\), the last two terms sum to strictly less than
\((j-1)/h\). If \(p\nmid j\), use
\(v_p((dj-1)!)<dj/(p-1)\le j/(2h)\).
If \(p\mid j\), add \(v_p(j)\le j/p<j/(4h)\) and use
\(j\ge p\ge5\). Thus all maximizing indices are in \(\{0,1\}\).
As \(R_m\bmod p=x^m\), this contradicts the \(h\ge2\) positive
valuations required by \(F\).

## Prime factors near the terminal order

For either order \(m\), every prime dividing
\(\binom mh(am+1)_{dh}\) must divide \(F(0)\). The factorial
factor gives the usual degree drop (T7). If it is a unit at
a prime dividing \(\binom mh\), then \(dh<p\), so \(h<p\).
Writing \(j=m\bmod p\), Lucas reduction gives \(j<h\).
It kills every binomial coefficient with \(h\le k<p\)
and the factorial coefficient kills the terms with \(k\ge p\).
Again \(R_m\bmod p\) is a power of \(x\) times a polynomial
of degree less than \(h\), proving the assertion.

There is a sharper exclusion at the lower end of this window.
Suppose a prime \(p\) divides \(N-1-j\), where \(0\le j<h/2\),
and \(p^4>M\). Then either of the following conditions excludes
degree \(h\):
\[
\begin{split}
 &p>2H+d\quad\text{and neither }[bm-6H,bm]\ (m=N-1,N)
                      \text{ contains a positive multiple of }p^3,\\
 &p>\max\{3H+d,\ j(H+d-1)+2H-1\}.
\end{split}
\tag{8}
\]

Indeed \(N-1\bmod p=j\) and \(N\bmod p=j+1\) in both cases.
At the lower order, the same Lucas argument puts at least
\(h-j\) roots of \(F\) at positive valuation. In the first
case, the first-jump argument in
[the unequal-edge bound](ORDINARY_BINOMIAL_RESIDUAL_REGION.md#large-primes-and-the-degree-of-the-padé-lift)
applies. It uses \(p^4>M\), the absence of terminal cubes and
the spacing \(p>2H+d\). Thus the shared positive slope is
unique, and equals the ratio of the height to the width of the
first jump in each reverse polygon. Let \(u_m\) be the multiple
of \(p\) reached at that jump, \(w_m\) its width, and
\(z_m=bm-u_m\in[d(w_m-1),dw_m-1]\).
Both widths are at least \(h-j\ge j+2\) and less than \(p\).
The binomial valuation at the first jump therefore equals
\(r_p=v_p(N-1-j)\ge1\) at both orders. Each remaining edge height
\(v_p(u_m)-r_p\) must be one since \(v_p(u_m)\le2\).
The common slope is \(1/w\) with equal widths \(w\le h\).
Writing \(\delta=z_N-z_{N-1}\), we have \(|\delta|\le d-1\),
\(b\equiv\delta\pmod p\) and \(z_{N-1}\equiv j\delta\pmod p\).
But \(j\delta-z_{N-1}\) is negative because
\(w-1\ge j+1\), and its absolute value is less than \(2H<p\).
This proves the first exclusion.

In the second case the same first-jump description holds
with cubes allowed. The bound \(U_m(i)<di/(p-1)+3\)
puts every maximizing index for a slope at least \(1/h\)
below \(9h/2<p\). Jumps are spaced at least \(3h+1\) indices
apart and have height at most three, so each later jump has
strictly smaller weighted value than the first.
The width comparison and the binomial valuation \(r_p\)
hold as before. Now both \(u_m\) are divisible by \(p^2\)
and the edge heights are at most two. Equal heights give the
same contradiction, as the equal widths are at most \(2h\)
and \(|j\delta-z_{N-1}|<p\) follows from the second threshold.
If the heights differ, they are one and two. The slope has
numerator one, the larger width is twice the smaller, and
the smaller width is at most \(h\). Thus
\(|\delta|\le H+d-1\) and \(0<z_{N-1}\le2H-1\).
The congruence forces \(z_{N-1}=j\delta\) by the threshold on
\(p\), and \(b\) is a unit since \(0<z_{N-1}<p\).
But \(p^2\mid u_m\) at both orders now gives \(b\equiv\delta\)
modulo \(p^2\), hence \(p^2\mid b(N-1-j)\). This contradicts
\(r_p=1\), which follows from the height-two edge.

### A cubic-depth exclusion

The threshold in (8) can be made uniform when cubes occur.
Suppose \(p\) is prime, \(p>2H+d\), \(p^4>M\), and
\[
 p\mid N-1-j,\qquad 0\le j<h/2.
 \tag{9}
\]
Then there is no common factor of degree \(h\).

Write \(r_-=j\), \(r_+=j+1\) for the residues of the two orders modulo \(p\), and \(v=v_p(N-1-j)\ge1\). Lucas reduction kills \(\binom{N-1}{k}\) for \(j<k<p\). The factorial part kills all terms of \(P_{N-1}\) at \(k\ge p\). Thus at most \(j\) of the common roots are \(p\)-adic units. At least \(h-j\ge j+2\) have positive valuation, using the evenness of \(h\).

For a selected slope, the strict affine bound \(U_m(k)<dk/(p-1)+3\) at positive \(k\) forces every maximizing index to satisfy \(k<6h\). Such indices lie below \(p\) if \(d\ge3\). At indices \(k<p\), writing \(r=m\bmod p\), the binomial valuation is zero through \(k=r\) and equals \(v\) from \(r+1\) on. In this range positive maximizing indices must be factorial jump indices. The first three possible jumps \(k_1,k_2,k_3\) reach consecutive multiples of \(p\) backwards from \(bm\). They are spaced at least \(\lfloor p/d\rfloor\ge2h+1\) apart, with at most one multiple divisible by \(p^2\), and no fourth jump can enter \(k<6h\).

The case \(d=2\) gives the same restriction to maximizing indices \(k<p\). In fact \(p>4h+2\) and \(k<6h<2p\). Beyond \(p+r\), the binomial valuation is at least \(v\), and the factorial valuation is at most five, so \(U_m(k)\le4<\lambda k\). On \([p,p+r]\) the binomial valuation is \(v-1\). An index can attain \(U_m(k)>4\) only if \(v=1\), the third jump has already entered, and one of the first three multiples is cubic. The earliest and best weighted value in this range is then at \(k_3=p+k_1\) with \(k_1\le r\). If the cubic is at the first jump, compare \(U_m(k_3)=5\) with \(U_m(k_1)=3\); the distance is \(p>4h\). If the cubic is at the third jump, compare with \(U_m(k_1)=1\) at the same distance. If it is at the second jump, \(k_2<p\) because \(k_1\le r\) and \(r\le h/2<p/2\), and \(U_m(k_2)=3\) at a distance at least \(2h+1\). In every case \(U_m(k_3)-\lambda k_3\) is strictly below the earlier value. This excludes all \(k\ge p\).

Among the three remaining jump indices, any simple jump after the first has smaller weighted value than its predecessor, since its height is at most one and the binomial valuation does not decrease. If the unique nonsimple jump is third, it loses to the first: the rise is at most four across a distance at least \(4h+2\). If it is second and \(k_1\le r\), it loses to the first because the binomial carry gives a rise at most \(3-v\le2\) over at least \(2h+1\) places. If \(k_1>r\) and is simple, its weighted value is negative by the same carry. Thus any selected edge must tie zero to just one jump, either first or second. There is only one possible slope of this form, since the tie requires \(\lambda=\max_{k>0}U_m(k)/k\). The selected jump has width at least \(h-j\), so it lies after the residue \(r\) in both orders. Its reached integer \(u_m\) is divisible by \(p^2\), either because it is the nonsimple first jump or because it is the nonsimple second jump. Its height is at most three after the binomial carry, and its width is at most \(3h\).

Set \(z_-=b(N-1)-u_{N-1}\), \(z_+=bN-u_N\), and \(\delta=z_+-z_-\). Then \(0<z_\pm<3H\), since each selected width is at least \(h-j\ge2\). The two square multiples give
\[
 b\equiv\delta\pmod{p^2},\qquad z_-\equiv bj\equiv j\delta\pmod{p^2}.
\]
The absolute value of \(z_--j\delta\) is less than \(3H(j+1)\le3H^2/4<p^2\), so \(z_-=j\delta\). If \(p\mid b\), both offsets are multiples of \(p\) strictly between zero and \(2p\). Both equal \(p\), giving \(\delta=0\), a contradiction. Otherwise the exact offset relation gives \(p^2\mid b(N-1-j)\), so \(v\ge2\).

This rules out the second jump because its height is at most \(1+3-v\le2\), while its width is at least \(2h+2\). At both first jumps the height is therefore one, and the widths equal some \(w\ge h-j\) by the common slope. Then \(|\delta|\le d-1\) and \(z_-\ge d(w-1)>j|\delta|\), contradicting \(z_-=j\delta\).

### An exact degree congruence at depth one

Let \(p\le m\), \(p^2>M\), and \(r=m\bmod p\) for either actual order \(m\). Then
\[
 h\bmod p\le2r
\]
is necessary for a common degree \(h\).

Indeed \(p>d\) follows from \(p^2>M\ge dN\ge dp\). Since \(m<p^2\), the binomial valuation is \(\epsilon(k)=\mathbf 1_{\{k\bmod p>r\}}\). The forward ordinates
\[
 y_m(k)=\lfloor(am+dk)/p\rfloor-\lfloor am/p\rfloor+\epsilon(k)
\]
and the reverse ordinates both increase by exactly \(d\) when their index increases by \(p\) within their domains. Put \(\mu=d/p\). For \(\lambda<\mu\), all minimizing indices of \(y_m(k)-\lambda k\) are in the first \(p\) positions. If \(r<k<p\), counting the multiples of \(p\) and adding \(\epsilon(k)=1\) gives \(y_m(k)>dk/p\). No such index minimizes because index zero gives value zero. All slopes below \(\mu\) therefore have total width at most \(r\). For \(\lambda>\mu\), the periodicity of the reverse ordinates likewise confines their maximizing indices to \(k<p\). At \(r<k<p\), subtracting the carry gives \(U_m(k)<dk/p\), so the total width above \(\mu\) is also at most \(r\). All remaining slopes equal \(\mu\), whose reduced denominator is \(p\), and the total width chosen by a factor from those slopes is a multiple of \(p\). This proves the degree congruence.

## A uniform consequence for low degrees

Let \(0<\kappa<e-1\) be fixed. On the part of the ordinary containing
region with
\[
 b\le \frac{d\log M}{(\log\log M)^2},
\]
all common factors have degree \(h>(\log M)^\kappa\) once \(M\) is
sufficiently large in terms of \(\kappa\) alone. Here the containing
region uses the explicit constraints of (T5) without the
upper terminal cutoff. In particular this applies uniformly to many
varying supports and to each fixed support satisfying the explicit
support constraints of (T5) for large orders.

Write \(\ell=\log M\) and \(P_0=d\ell/\log\ell\). The proof uses
\(d^2<128b\) and \(d<3a\) from (T5). The displayed restriction on \(b\)
gives \(d<128\ell/(\log\ell)^2\) and \(b=o(P_0)\), uniformly.
If \(h\le\ell^\kappa\), all primes in (2) satisfy
\[
 p\le B=\max\{2dh+1,2b,4dh\ell\}
\]
for large \(\ell\): for \(p>2b\) with \(p\le b+2dhL_p\), one has
\(p\log p\le4dh\ell\). The contribution of the \(\mu_p\) terms
to \(\log|F(0)|/h\) is at most
\(d\log B(1+\log B)=o(d\ell)\), by comparison with the harmonic sum.
For \(p\le P_0\), the excess term satisfies
\[
 (\log p)A_p\le
 \max\{\log P_0,\,3\log(b+d\ell/\log 2)\}=O(\log\ell).
\]
Using \(\pi(P_0)=O(P_0/\log P_0)\) gives \(O(d\ell/\log\ell)=o(d\ell)\)
for their total contribution.

For \(p>P_0\), the second cap in \(C_p\) bounds the excess contribution
per prime by \(d\ell/(p-b-1)\). Since \(b=o(P_0)\), the total over
\(P_0<p\le B\) is at most
\[
 (1+o(1))d\ell\sum_{P_0<p\le B}\frac1p
 \ \le\ d\ell\bigl(\log(1+\kappa)+o(1)\bigr).
\]
For the last inequality use the prime reciprocal estimate
\(\sum_{p\le x}p^{-1}=\log\log x+C+O((\log x)^{-2})\) and
\(\log B/\log P_0\le1+\kappa+o(1)\). Both the needed prime count and
the prime reciprocal estimate follow from Dusart's Theorems 6.9 and
6.10. [Dusart, *Estimates of some functions over primes without R.H.*](https://arxiv.org/pdf/1002.0442).

On the other hand \(aN>M/4\) and the strict lower bound in (3)
gives \(\log|F(0)|/h>d(\ell-\log 4)\). This is impossible because
\(\log(1+\kappa)<1\). The result excludes a range of low degrees
while leaving the coprimality test for any degrees that survive
both these bounds and the residual degree restrictions.
