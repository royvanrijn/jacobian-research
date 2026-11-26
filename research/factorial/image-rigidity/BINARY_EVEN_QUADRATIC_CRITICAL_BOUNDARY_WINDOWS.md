# Critical boundary polynomials of an even binary quadratic

Let \(K\) have characteristic zero and let
\(\Lambda_2(x^i y^j)=i!j!\). Write
\[
 f=a x^2+bxy+c y^2+k,\qquad M_j=\Lambda_2(f^j),
 \qquad abck\ne0 .
\]
The following bounds hold at every positive start index.

**Theorem.** Every such phase has a nonzero moment in every six
consecutive orders. The exact four-monomial window holds on the
plane \(a+c=b\) and on the rank-one stratum \(4ac=b^2\), as well as
on the three [previously covered coefficient planes](BINARY_EVEN_QUADRATIC_SINGLE_FACE_STRONG_FACTORIAL.md).
For any remaining possible four-zero
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

## From four weighted boundary zeros to two critical points

Let \(a,c\) be distinct and nonzero and let the two weights
be nonzero. Suppose the two boundary polynomials have
\(g_j(c)=\lambda g_j(a)\), with a fixed \(\lambda\ne0\),
for four consecutive \(j\). The univariate factorial
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

For each of the last two indices \(j\) of the block, both
inputs to \(L_j\) satisfy the common ratio. Therefore
\[
 (\lambda a-\lambda c)L_j(a)=1-\lambda,\qquad
                 g_j(a)=1+aL_j(a).
\]
Thus the last two values of \(g_j(a)\) are equal, and so
are the corresponding two values for \(c\). At the last
order \(m\), (1) makes \(a,c\) critical points of \(g_m\).
Their values are nonzero because \(g_m\) is squarefree.
Writing \(g_{m-1}(c)=\lambda g_{m-1}(a)\) in (3) also
gives the same number at the two critical points:
\[
       F_m(A)=\frac{g_m(A)-1}{A g_m(A)},\qquad F_m(a)=F_m(c).
                                                           \tag{4a}
\]
In fact \(F_m=L_m/g_m\) there, and both numerator and
denominator satisfy the same boundary ratio.

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
weights are nonzero and \(a\ne c\). Six zero moments make
the last four \(T_j\) zero. The injectivity proof below rules
out (4a), and proves this six-window.

## The backward axis recurrence

For \(A\ne0\), put \(w=1/A\), \(e_i=1/(2^i i!)\) and
\(U_i=e_i g_i(A)\). Equation (3), read backwards, gives
\[
           U_{i-2}=(2i-1)U_{i-1}-wU_i+w e_i,\qquad i\ge2.
                                                           \tag{7}
\]
The starting pair \(U_{m-2},U_{m-1}\) determines \(U_0\)
using the forcing terms at indices \(2,\ldots,m-1\).
The coefficient of \(U_{m-2}\) is a polynomial \(Q(w)\)
over \(\mathbb Z\) with odd constant, since at \(w=0\)
one follows the adjacent steps, multiplying only by odd
integers. The coefficient of \(U_{m-1}\) is a polynomial
in \(w\mathbb Z[w]\), since it can enter only by a skip
with coefficient \(-w\). Each forcing term \(w e_i\)
has a coefficient polynomial with odd constant, by the
same adjacent-step argument. These observations include
the initial case \(m=3\), where the remaining recurrence
has just one step. They imply precise valuation comparisons
when \(v_2(w)>0\).

## Injectivity of the critical point boundary statistic

For \(i\ge2\), write \(R_i(w)\) for the polynomial
multiplier carrying a forcing term at index \(i\) in (7)
to \(U_0\). Set \(R_1=0\). Then \(R_2=1\), and
\[
        R_i=(2i-3)R_{i-1}-wR_{i-2}\quad(i\ge3).
\]
The two terms carry a value at \(U_{i-2}\) along its next
adjacent step and its next skip, respectively. In particular
\(R_i(0)\) is odd for \(i\ge2\).

Let \(m\ge4\), let \(A\) be a critical point of \(g_m\),
and put \(y=g_m(A)=g_{m-1}(A)\), \(w=1/A\),
\(Q=R_m\) and
\[
        C=\sum_{i=2}^{m-1}(e_i/e_{m-1})R_i,\qquad d=1/e_m .
\]
The equation for \(U_0\) made explicit by (7) is
\[
  1=e_{m-1}y\bigl((2m-1)Q-wR_{m-1}\bigr)
                      +w e_m Q(1-y)+w e_{m-1} C .
\]
Indeed \(U_{m-2}=(2m-1)e_{m-1}y+w e_m(1-y)\),
and the contribution entering from \(U_{m-1}\) after
this initial step is \(-w e_{m-1}y R_{m-1}\).
Substituting \(1/y=1-F_m(A)/w\) gives
\[
 \frac{F_m(A)}{2m}
    =\frac{N(w)-d/(2m)}{D_1(w)-d/w},\quad
 N=(2m-1)Q+w(C-R_{m-1}),\quad D_1=Q+2mC.                  \tag{8a}
\]
These denominators are units at the critical points, as
shown below.

In the integral power series \(N/D_1\) at zero,
all nonconstant coefficients are even and its linear
coefficient is \(2\) modulo \(4\). To see this,
\(Q(0)\) is odd, \(C\equiv R_{m-1}\pmod2\), and
\(C-R_{m-1}\equiv 2(m-1)R_{m-2}\pmod4\).
Divide numerator and denominator by \(Q\). The numerator
contributes \(2(m-1)\pmod4\) to the derivative at zero.
The denominator contributes \(-2m\pmod4\), since
\((R_{m-1}/R_m)'(0)\equiv1\pmod2\) by the recurrence for
\(R_m\) and the odd constant coefficients. Modulo two the
quotient is constant. Thus at two distinct points
\(w,w'\) of positive valuation its difference of values
has exactly valuation \(1+v_2(w-w')\):
every term of degree at least two has strictly higher
difference valuation after factoring \(w-w'\).

At a critical point the positive Newton slopes for \(h_{m-1}\)
give \(v_2(w)=2-1/l\) for powers of two \(l\ge2\).
For even \(m\) the sole additional point has
\(v_2(w)=1+v_2(m)\). This is the first, width-one
edge in (2). In particular both \(Q(w)\) and \(D_1(w)\)
are units. At any pair of distinct critical points,
\[
             v_2(d)-v_2(w)-v_2(w')>1 .
\]
For two positive slopes, use \(v_2(w),v_2(w')<2\)
and \(v_2(d)=2m-s_2(m)\ge6\) for \(m\ge4\).
If one point comes from the additional edge the sum
of the two valuations is less than \(3+v_2(m)\).
Here \(2m-s_2(m)\ge m+1\) and
\(m+1\ge4+v_2(m)\) for even \(m\ge6\); the remaining
case \(m=4\) follows from \(v_2(d)=7\).
Two points cannot both come from the single width-one
edge. These estimates also give \(v_2(d/w)>1\).

For such a pair the correction between (8a) and \(N/D_1\)
has a difference of values of valuation strictly above
\(1+v_2(w-w')\). In detail, that correction is
\[
           \frac{N(d/w)-D_1(d/(2m))}
                         {D_1(D_1-d/w)} .
\]
The Laurent term \(d/w\) has difference valuation
\(v_2(w-w')+v_2(d)-v_2(w)-v_2(w')\), which is strictly
above the required bound. Its values and the constant
\(d/(2m)\) have valuation strictly above one. Integral
polynomials and their unit-valued inverses have differences
divisible by \(w-w'\) in the valuation ring. These facts
give the bound for the correction by the product and
quotient rules for differences. Thus (8a) is injective
on distinct critical points. In (4a) the last index is at
least four, so this proves that four weighted cancellations
are impossible.

## The plane \(a+c=b\)

Assume a four-zero block starts at \(n\) and \(D=0\).
If \(a=c\), the previous coefficient planes already apply.
Otherwise \(\alpha=c-b/2=(c-a)/2\) and \(\beta=-\alpha\),
while \(\Delta=-(a-c)^2/4\) is nonzero. Equation (4),
including at the order immediately after the four-zero
block, gives
\[
               g_j(a)=g_j(c),\qquad j=n+2,n+3,n+4 .
\]
Set \(m=n+4\). Subtracting the two recurrences (3) at \(m\)
shows \(L_m(a)=L_m(c)=0\), and hence
\(g_m(a)=g_m(c)=1\).

We show that \(g_{m-1}\) is injective on the nonzero
\(2\)-adic solutions of \(g_m(A)=1\). The identity
\[
        (g_m(A)-1)/(2mA)
          =\sum_{i=0}^{m-1}\binom{m-1}{i}(2i+1)!\,A^i
\]
comes from (3) and the monomial expansion of
\(\Lambda_1(x(1+Ax^2)^{m-1})\).
Its constant is one and the coefficient at \(i>0\)
has valuation at least \(i\), so every nonzero solution
has \(v_2(A)\le-1\). Thus \(|w|_2\le1/2\).
Write \(y=g_{m-1}(A)\). Since \(g_m(A)=1\),
the initial step (7) gives
\[
          U_{m-2}=(2m-1)e_{m-1}y,\qquad U_{m-1}=e_{m-1}y.
\]
After division by \(e_{m-1}\), the backward equation
for \(U_0=1\) therefore reads
\[
        1/e_{m-1}=P(w)y+
                    w\sum_{i=2}^{m-1}(e_i/e_{m-1})R_i(w)
\]
with \(P(0)\) and \(R_{m-1}(0)\) odd and all polynomials
integral. The factors \(e_i/e_{m-1}\) with \(i<m-1\)
and the constant \(1/e_{m-1}\) are even integers.
Expanding \(1/P(w)\) about zero after division by its
unit constant shows that \(y\) is an integral power
series in \(w\), with an odd linear coefficient.
On \(|w|_2\le1/2\), every term of index at least two
has difference of values bounded in magnitude by
\(\tfrac12|w-w'|_2\) times its coefficient magnitude.
The linear term has difference of magnitude
\(|w-w'|_2\). Thus \(a\ne c\) cannot have the same
\(g_{m-1}\) value. This excludes the four-zero block.

## The rank-one stratum

A further consequence of (7) is that two distinct critical
points of \(g_m\) cannot have critical values with sum one.
Indeed let \(A\) be one such point and write \(y=g_m(A)\).
Equation (1) gives \(g_{m-1}(A)=y\). By the Newton
polygon (2), each critical point has
\(\lambda=v_2(4A)=1/l\) for a power of two \(l\ge2\),
except for a single point when \(m\) is even.
That point comes from the first edge of the odd
\(h_{m-1}\) and has \(\lambda=1-v_2(m)\le0\).

At an ordinary positive slope \(\lambda\) in this list,
\(v_2(w)=2-\lambda>0\). Start (7) with
\[
        U_{m-2}=(2m-1)e_{m-1}y+w e_m(1-y),\qquad
        U_{m-1}=e_{m-1}y .
\]
The backward equation for \(U_0=1\) has the form
\[
              1=Pe_{m-1}y+Qw e_m(1-y)+E .                 \tag{8}
\]
Here \(P,Q\) are \(2\)-adic units, by the adjacent-step
coefficients proved above, and \(v_2(E)>v_2(w e_m)\)
because the remaining forcing terms use \(e_i\) at
\(i<m\) and integral polynomial multipliers. Put
\(C=v_2(w e_m)<0\) and \(B=v_2(e_{m-1})\). Then
\(C-B=1-\lambda-v_2(m)\). If \(m\) is odd, \(B<C<0\);
collecting \(y\) in (8) gives \(v_2(y)=C-B=1-\lambda>0\).
If \(m\) is even, \(C<B\); instead collect \(y-1\)
to find \(v_2(y-1)>0\), since the coefficient of \(y-1\)
has valuation \(C\) and all other terms have valuation
strictly above \(C\). At the single exceptional point
with \(\lambda\le0\), the normalized reciprocal of
\(g_m(A)\) has all nonconstant coefficients of positive
valuation, so the leading term of \(g_m(A)\) dominates.
In this case \(v_2(y)=m\lambda-s_2(m)<0\).
For odd \(m\) the sum of any two critical values has
positive valuation. For even \(m\) the sum of two
values from the positive slopes has residue zero,
and a sum including the single remaining critical
point has negative valuation. In either case the
sum of values at distinct critical points cannot be one.

Finally suppose \(\Delta=0\), and work in an algebraic
extension containing \(u,v\) with \(a=u^2,\ c=v^2,\
b=2uv\). Outside the known planes we have \(u\ne v\)
and \(u\ne -v\). At four consecutive zero moments
(4) reduces to \(u g_j(a)=v g_j(c)\).
At the last two indices \(j=m-1,m\), where \(m=n+3\),
the two inputs to (3) also satisfy this ratio. Subtraction
therefore gives \(L_j(a)=-1/[u(u+v)]\). Consequently
\[
               g_j(a)=\frac{v}{u+v},\qquad
               g_j(c)=\frac{u}{u+v},\qquad j=m-1,m .
\]
By (1), both \(a,c\) are distinct critical points of
\(g_m\), and the displayed values sum to one. This
contradicts (8) and the critical value calculation.

For these two four-window proofs the reduction to a
\(2\)-adic algebraic solution is legitimate. At fixed
order the equations are polynomial over \(\mathbb Q\).
Adjoin inverses of the required nonzero expressions.
A characteristic-zero solution makes the ideal proper,
and by base extension and the weak Nullstellensatz it
has a solution over \(\overline{\mathbb Q}\), hence
over \(\overline{\mathbb Q}_2\).

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
