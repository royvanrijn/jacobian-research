# Critical boundary polynomials of an even binary quadratic

Let \(K\) have characteristic zero and let
\(\Lambda_2(x^i y^j)=i!j!\). Write
\[
 f=a x^2+bxy+c y^2+k,\qquad M_j=\Lambda_2(f^j),
 \qquad abck\ne0 .
\]
The following bounds hold at every positive start index.

**Theorem.** Every such phase has a nonzero moment in every seven
consecutive orders. On the plane \(a+c=b\) six orders suffice.
On the rank-one stratum \(4ac=b^2\) five orders suffice. The three
[known coefficient planes](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md)
retain their four-window. For any remaining possible four-zero
block, after normalizing \(k=1\), the following restrictions hold
in every realization over \(\overline{\mathbb Q}_2\).

- Put \(r=n+2\), and let \(l_*\) be the largest power of two
  at most \(r+1\). The two axis valuations cannot both exceed
  \(-3/2\), and cannot both be below \(-2+1/l_*\).
- At an even starting order they cannot straddle those two
  regions. At an odd starting order, if \(v_2(a)<-2+1/l_*\)
  and \(v_2(c)>-3/2\), then either
  \(v_2(a)+v_2(c)=-2-v_2(n+3)\) when \(v_2(c)<-1\), or
  \(v_2(a)\le-1-v_2(n+3)\) when \(v_2(c)\ge-1\).
  In the latter case, a strict inequality forces \(2c\)
  to have residue one. These rules also hold with \(a,c\)
  exchanged.
- If \(v_2(a),v_2(c)>-5/3\) and \(v_2(b)\ge-1\), the moment
  at the index divisible by four in the block is a unit.

The longer window bound controls arbitrary complex coefficients.
The valuation restrictions use the actual boundary polynomials.
The four-window outside the known planes remains open.

## Squarefreeness of the once-shifted boundary moments

Scale \(f\) by \(1/k\) to take \(k=1\). Define the axis polynomials
\[
 g_j(A)=\sum_{i=0}^j\binom ji(2i)!\,A^i,\qquad
 h_r(A)=\sum_{i=0}^r\binom ri\frac{(2i+2)!}{2!}\,A^i .
\]
They satisfy
\[
 g'_j(A)=2j h_{j-1}(A),\qquad
 A g'_j(A)=j\bigl(g_j(A)-g_{j-1}(A)\bigr),\qquad
 A h'_r(A)=r\bigl(h_r(A)-h_{r-1}(A)\bigr).
 \tag{1}
\]
The first and the other two identities follow by differentiating
the corresponding factorial moments and replacing \(A x^2\) by
\((1+A x^2)-1\) in the latter two.

We prove that every \(h_r\) is squarefree. For \(r\ge1\), take the
reciprocal of \(h_r(B/4)\). Number its coefficients by \(d=0,\ldots,r\)
and shift all valuations by \(s_2(r)\), where \(s_2\) is binary
digit sum. The shifted ordinates are exactly
\[
                  s_2(d)+v_2(r+1-d).                      \tag{2}
\]
In fact \((2i+2)!/2!=2^i(i+1)!(2i+1)!!\), so the
coefficient at \(i=r-d\) of \(h_r(B/4)\) has valuation
\(v_2(r!)-i-v_2(d!)+v_2(i+1)\). This is (2) minus
\(s_2(r)\) by the factorial valuation formula.

For even \(r\), the binary prefixes \(P\) obtained by adding the
set bits of \(r\) from largest to smallest are all even. At these
indices the second term in (2) is zero. The lower hull joins
successive prefixes, with an edge of width \(l\) and rise one
for each set bit \(l\). Indeed between two consecutive prefixes
\(P,P+l\), one has \(s_2(P+u)=s_2(P)+s_2(u)\) for \(0<u<l\),
and the extra term in (2) is nonnegative.

For odd \(r\), the hull first joins
\((0,v_2(r+1))\) to \((1,1)\), and then uses the binary
prefixes of \(r-1\) translated by \((1,1)\).
The translated prefix indices \(1+P\) are odd and have
ordinate \(1+s_2(P)\). At an odd interior index \(1+P+u\)
between \(1+P\) and \(1+P+l\), the digit sum is at least
\(2+s_2(P)\). An even interior index instead has the form
\(P+w\) with \(2\le w\le l\) even. For \(w<l\) its digit
sum is at least \(s_2(P)+1\) by disjoint digits, and for
\(w=l\) it is \(s_2(P)+1\) since \(P\) omits that bit.
In either case \(v_2(r+1-P-w)\ge1\). Hence all those interior
points lie above the indicated chords. The first slope is
nonpositive and the subsequent slopes are positive.

Each edge has width equal to the denominator of its slope
in lowest terms. This includes the first edge of width one
when \(r\) is odd. By the Newton polygon product rule, an
irreducible factor on that slope has degree a positive multiple
of the denominator. There is room for exactly one such factor,
with multiplicity one. Each irreducible factor is separable in
characteristic zero, and distinct slopes have disjoint root
valuations. This proves squarefreeness over \(\mathbb Q_2\)
and over \(\mathbb Q\). [Newton polygons, Theorem 2.1](https://www.impan.pl/shop/publication/transaction/download/product/91656)

By (1), \(h_r,h_{r-1}\) have no simultaneous zero: \(A=0\)
makes each value one, and at a nonzero common zero \(h_r'\)
would vanish. Thus \(g_m'\) and \(g_{m-1}'\) are coprime in
every characteristic-zero field, including \(m=2\) where
\(h_0=1\). In particular a nonzero \(A\) cannot have
\(g_{m-2}(A)=g_{m-1}(A)=g_m(A)\).

## Five weighted boundary zeros are impossible

Let \(a,c\) be distinct and nonzero and let the two weights
be nonzero. Suppose the two boundary polynomials have
\(g_j(c)=\lambda g_j(a)\), with a fixed \(\lambda\ne0\),
for five consecutive \(j\). The univariate factorial
Stein identity gives
\[
 g_j(A)=1+A L_j(A),\qquad
 L_j(A)=2 j[(2j-1)g_{j-1}(A)-2(j-1)g_{j-2}(A)] .
 \tag{3}
\]
For \(j=1\) omit the term with zero multiplier.
Indeed \(\Lambda_1((1+A x^2)^j)-1=2jA\Lambda_1(x(1+A x^2)^{j-1})\),
and applying the same monomial identity to \(x(1+A x^2)^{j-1}\)
evaluates the last functional.

For each of the last three indices \(j\) of the block, both
inputs to \(L_j\) satisfy the common ratio. Therefore
\[
 (\lambda a-\lambda c)L_j(a)=1-\lambda,\qquad
                 g_j(a)=1+aL_j(a).
\]
Thus the last three values of \(g_j(a)\) are equal, which
the squarefreeness result excludes.

For the binary phase, put
\[
 D=a+c-b,\quad \Delta=ac-b^2/4,\quad
 T_j=(c-b/2)g_j(c)+(a-b/2)g_j(a).
\]
The [exact weighted face recurrence](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md#4-the-exact-reduction-for-general-b-and-its-limitation)
reads
\[
       D M_j=4\Delta j^2 M_{j-1}
                    -4\Delta j(j-1)M_{j-2}+T_j.            \tag{4}
\]
At orders outside the three known coefficient planes the
weights are nonzero and \(a\ne c\). Seven zero moments make
the last five \(T_j\) zero. If \(D=0\), a block of six zero
moments makes five consecutive \(T_j\) zero, including the
order immediately after the block. If \(\Delta=0\), every
vanishing \(M_j\) already gives \(T_j=0\). This proves the
three unconditional window bounds in the theorem.

## Exclusion using the boundary ratios at the two final orders

For a fixed \(r\ge1\), the two boundary cancellations from a
four-zero block require equality of the projective values
\[
             [g_{r+1}(a):g_r(a)]=[g_{r+1}(c):g_r(c)],
               \qquad r=n+2 .                              \tag{5}
\]
The [axis two-window](UNIVARIATE_QUADRATIC_ALL_SUPPORTS_STRONG_FACTORIAL.md#2-constant-plus-square-a-squarefree-moment-polynomial)
ensures each pair defines a projective point, even when a
single moment vanishes. Both weighted relations yield (5).

For the inner disk put \(U=2A\). The coefficient at index
\(i\) of \(g_j(U/2)\) has valuation
\(i+s_2(j-i)-s_2(j)\). For even \(j\) this is
at least \(i/2\). Indeed for even \(i\) one has
\(s_2(i)\le i/2\); for odd \(i\) and even \(j\),
adding \(i\) to \(j-i\) involves at least one carry
and improves the bound \(i-s_2(i)\) by one.
For odd \(j\) the valuation is at least \((i-1)/2\)
at every \(i\ge1\), since \(i-s_2(i)\ge(i-1)/2\).
Use the even index of the pair in the denominator.
Its higher terms have positive valuation for
\(v_2(U)>-1/2\), so its inverse is defined by a
geometric series with coefficients at least \(i/2\)
in valuation at degree \(i\). This follows by multiplying
the degreewise bounds in each product of positive-degree terms.

The ratio of the odd-order polynomial to the even-order
polynomial now has the form
\[
                1\ \pm U+\sum_{i\ge2} e_i U^i,\qquad
                         v_2(e_i)\ge(i-1)/2 .
\]
The linear coefficient is the difference of the two
successive order indices, which is \(\pm1\), and the
higher bounds follow by convolution. For any two arguments
with valuation above \(-1/2\), the difference of the
higher terms has valuation strictly greater than that of
the difference of the arguments. Factor the difference
out of each monomial and use
\(v_2(e_i)+(i-1)\min(v_2(U),v_2(U'))>0\).
Thus this chart for (5) is injective. Its ratio at \(A\)
has valuation \(1+v_2(A)\) if \(v_2(A)<-1\). For
\(v_2(A)\ge-1\) it is integral and reduces to \(1+2A\)
in the residue field.

For the outer region put \(Z=1/(4A)\) and
\(Q_j(Z)=g_j(A)/((2j)!A^j)\). Its constant coefficient
is one and the coefficient at \(d\ge1\) has valuation
\(s_2(d)\), by \(v_2((2j)!)=2j-s_2(j)\).
For \(j=r,r+1\), \(s_2(d)\ge d/l_*\), because the
set bits that sum to \(d\le r+1\) are all at most \(l_*\).
For any two points in \(v_2(Z)>-1/l_*\), choose \(t\ne0\)
whose valuation is the smaller of their two \(Z\)
valuations. Writing \(Z=tX\) places the two points
in \(|X|_2\le1\). The coefficients of positive degree
of both \(Q_j(tX)\) have valuation at least some
\(\delta>0\). The same bound holds for the positive
degrees of their ratio \(R(X)=Q_{r+1}(tX)/Q_r(tX)\),
by geometric expansion. The boundary map, up to a fixed
nonzero factor depending only on \(r\), is \(A R(X)\).
Its values at \(a,c\) differ by \(a-c\) plus the
difference of
\((4t)^{-1}\sum_{d\ge1}e_d X^{d-1}\), where
\(e_d=[X^d]R\). The latter difference has magnitude at
most \(2^{-\delta}|a-c|_2\), since
\(a-c=(4t)^{-1}(X_c-X_a)/(X_c X_a)\).
This proves injectivity on the whole outer region.

For an outer \(a\) and inner \(c\), the ratio with odd order
in the numerator has valuation at \(a\) equal to
\[
     1+v_2(a) \quad(r\ \text{even}),\qquad
    -1-v_2(r+1)-v_2(a) \quad(r\ \text{odd}),
\]
by the unit values of \(Q_j\). In the first case
\(r\ge4\), so the valuation is below \(-1/2\);
the inner ratio always has valuation greater than \(-1/2\).
For odd \(r\), compare with the two inner cases just
proved to obtain the stated restrictions. Their equality
and residue conditions concern the original boundary values
at \(r=n+2,n+3\).

## The interior unit at an index divisible by four

This test uses a full binary moment rather than (5).
Write \(v=v_2\). At an integer \(j>0\) divisible by four,
let \(d_0=2^{v(j)}\ge4\). For counts \(p,q,r,s\) of the
monomials \(a x^2,c y^2,bxy,1\), with sum \(j\), put
\(l=p+q+r\). The constant term with \(l=0\) is one.
For \(l>0\), using the factorial valuation formula and
\(s_2(x+y)\le s_2(x)+s_2(y)\) gives the lower bound
\[
   p v(4a)+q v(4c)+r v(4b)-s_2(r)+s_2(j-l)-s_2(j)        \tag{6}
\]
on the valuation of the corresponding multinomial term.
In detail the factorial indices are \(2p+r,2q+r\), and
their digit sums are at most \(s_2(p)+s_2(r)\) and
\(s_2(q)+s_2(r)\). Subtract these two bounds from the
multinomial digit sums to get (6).

Under the stated coefficient assumptions choose
\(\gamma=\min(v(4a),v(4c),1/2)>1/3\).
For \(r=1\) the \(b\) contribution in (6) is at least
zero. For \(r\ge2\) it is at least \(r-s_2(r)\ge(r-1)/2\),
as each set bit other than one contributes at least two
to \(r\). The same weaker bound with \(\gamma(l-1)\)
for all non-digit-sum terms holds also at \(r=0\).
When \(l<d_0\), the digit difference is
\(v(j)-1-s_2(l-1)\ge0\). At \(l=1\) it is strictly
positive, and for \(l>1\) the coefficient terms contribute
at least \((l-1)\gamma>0\).
When \(l\ge d_0\), the lower binary polygon starting at
\(l=0\) gives \(s_2(j-l)-s_2(j)\ge-l/d_0\). To see
this bound, read the binary prefix polygon of \(j\) from
its final endpoint, so its first width is \(d_0\) and
its subsequent widths are larger. Consequently (6) is at
least
\((l-1)\gamma-l/d_0>0\), since \(d_0\ge4\), \(l\ge d_0\)
and \(\gamma>1/3\). All nonconstant terms are strictly
positive in valuation, and \(M_j\) is a unit.

For applying the local statements to algebraic coefficients,
they hold in each chosen \(2\)-adic embedding. For the
unrestricted characteristic-zero existence question at a
fixed start index, a solution to the polynomial equalities
and inequalities over any characteristic-zero field would
give an algebraic solution by the weak Nullstellensatz
after adjoining inverses for the required nonzero
expressions. Thus the valuation restrictions give necessary
conditions on an algebraic counterexample to the full
four-window, without assuming that every original solution
already has algebraic coefficients.
