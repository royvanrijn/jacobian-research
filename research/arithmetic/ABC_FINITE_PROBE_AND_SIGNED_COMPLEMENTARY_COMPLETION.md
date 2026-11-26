# Finite valuation probes and mixed-sign complementary products

Positive quotient masks on finitely many scaled-degree intervals suffice
for exact word recovery at an arbitrarily small radical-scale modulus.
The probe count stays bounded as the budget decreases. For words with
two positive weights and one negative weight, the two exact
complementary products also give a power saving on distinct-row
fibres. This second step uses the original radical caps on all three
rows.

## Finite probes by Möbius inversion

Use the height-first balanced valuation boxes, the prescribed prime
pool \(P\), and the full-divisor chart of
[the small-modulus word theorem](ABC_SMALL_MODULUS_TRIPLE_PRODUCT_COMPLETION.md#signed-words-and-separation-of-weighted-slot-sums).
Write \(n=X^L\). Fix the integer weights of the word, the chart
allowance, and \(\delta,\eta>0\) before taking sufficiently late
budgets. The small modulus is the first power of two exceeding
\(n^\delta\). For integer \(0<q<q'\), allow the step \(q'-q\) to
depend on the probe and set
\[
 \phi_{q,q'}(e)=\lfloor e/q\rfloor-\lfloor e/q'\rfloor,
 \qquad M_{q,q'}(W)=\prod_{p\notin P}p^{\phi_{q,q'}(v_p(W))}.
\]
Each pair of residues at a probe is
\((M_{q,q'}(A)/M_{q,q'}(C),\,M_{q,q'}(B)/M_{q,q'}(C))\).

**FP1.** The fixed modular word permits
\(O(n^\eta X^\varepsilon)\) pairs of exact complementary
products using a number of probes bounded in terms of the
limiting measures, fixed allowances and fixed weights, uniformly
as \(L\to0\). In particular the positive three-row completion
and the signed six-slot energy estimate in (SC1) and (SC3)
hold with a bounded number of probes.

For any integer \(e\ge1\) one has the exact identities
\[
 \phi_{q,q'}(e)=\sum_{d\ge1}1_{[dq,dq')}(e),\qquad
 1_{[q,q')}(e)=\sum_{k\ge1}\mu(k)\phi_{kq,kq'}(e).
 \tag{FP2}
\]
Expanding the second sum gives coefficient
\(\sum_{k\mid d}\mu(k)=1_{d=1}\) on each interval at
multiple \(d\). The sums are finite for each \(e\). Also
\(0\le\phi_{kq,kq'}(e)\le\phi_{q,q'}(e)\), since the
dilated mask selects a subset of the base mask's intervals.

Here are the choices that keep these masks short. Fix the
allowed total exceptional raw mass \(\tau>0\) used in (SC3)
and let \(W_*=\sum_i|w_i|\). Choose a mask allowance
\(0<\zeta<\delta/(4W_*)\). Let \(\nu\) be the sum of the
three limiting measures outside \(P\). Their first moments
are bounded by the uniform bound on the sum of the three
letter heights. Denote a bound for this sum by \(H_0\). The
same bound controls the first moment in the inner limit at
each fixed \(L\). Choose \(R\) sufficiently large that
\(H_0/R\) is small compared with both \(\tau\) and \(\zeta\).
Choose \(0<r<R\) small enough that the nonatomic mass below
\(r\) is less than \(\tau/10\). Choose the positive essential
atoms in \([0,R]\) to leave atomic mass less than
\(\min(\tau,\zeta)/10\) on that interval, apart from the
atom at zero, which has its own essential window. These are
finitely many positive atoms. Continuity endpoints may be
used for \(r,R\).

For every essential positive atom \(s\) remove small
neighborhoods of each \(s/d\) in \([r,R]\), for
\(1\le d\le 2+\lceil R/r\rceil\). Their total nonatomic
mass can be made less than \(\tau/10\). Remove slightly
larger neighborhoods if needed so that the remaining
intervals have positive distance from these aliases. For
any fixed \(\omega>0\), cover the remaining compact
nonatomic mass, up to a loss less than \(\tau/10\), by
finitely many disjoint half-open intervals \(I_j=[a_j,b_j)\)
inside \((r,R)\). Take their widths at most \(\omega/W_*\).
Their endpoints can be chosen at continuity points and
their widths can be made uniformly as small as required.
If a dilated interval \([d a_j,d b_j)\) intersected an
essential atom in \([0,R]\), the corresponding alias would
lie in \(I_j\). The positive distance also excludes every
shrinking essential window at sufficiently late budgets.
The zero window lies below the first interval. Positive
atoms outside \([0,R]\) can be left exceptional because the
whole tail has small raw mass.

Set \(q_j=\lceil a_j/L\rceil,\ q'_j=\lceil b_j/L\rceil\).
The intervals \([q_j,q'_j)\) are disjoint and lie strictly
below \(R/L\) at late budgets. As the mask list use
\((kq_j,kq'_j)\) with \(1\le k\le\lceil R/r\rceil\) and
\(\mu(k)\ne0\). The number of base intervals and this bound
for the dilations are fixed before \(L\) is taken small.

For \(Le\le R\) the height of a base mask is controlled by
the sum of the masses on a bounded number of harmonic
intervals. Indeed \(Lq_j\ge r\), so only at most
\(\lceil R/r\rceil\) intervals can occur. On the nonatomic
part the mass of an interval tends uniformly to zero with
its width. Every harmonic width is at most
\(\lceil R/r\rceil(b_j-a_j+O(L))\). The essential atoms do
not meet these intervals, and the mass of the remaining
atomic part is bounded by its total mass on \([0,R]\).
To keep its contribution below \(\zeta/10\), one may shrink
the base intervals to have \(d(b_j-a_j+O(L))<r/2\) for all
the indicated \(d\), so the harmonic intervals of a given
base mask are disjoint. These bounds pass to the actual
boxed slots by weak convergence, using fixed slightly
enlarged intervals with continuity endpoints and the strict
margins in the mass bounds. For \(Le>R\) use
\[
 \phi_{q_j,q'_j}(e)
       \le 1+(Le)\frac{b_j-a_j+O(L)}{r^2}.
\]
That part of the normalized logarithmic mask height is at most
\(H_0/R+H_0(b_j-a_j+O(L))/r^2\). Choosing the intervals
narrow enough gives total normalized mask height less than
\(\zeta\), uniformly over the rows at sufficiently late
budgets, with inner rounding margins. Each dilated mask
has the same bound by (FP2). Exceptional slots are those
outside the essential windows and the union of the base
intervals. Their normalized raw mass is less than \(\tau\)
by the choices of \(r,R\), the atoms and the intervals.

Now lift the two word coordinates at each probe by rational
reconstruction, as in (SC3). Its determinant test applies
to the full list because each mask has height at most
\(n^\zeta\). Include the supports of these lifted rationals,
the prescribed pool and the recorded exceptional primes in
a known pool. On any other prime, the two signed differences
of mask exponents vanish at each probe. Essential-window
terms have mask zero. Apply (FP2) separately to each base
interval to conclude that the signed sum of membership
indicators in each good interval vanishes for both differences.
Subtract the left endpoint of the interval from every scaled
degree in it. The signed scaled degree sum on the good slots
then has absolute value at most \(\omega+o(1)\) for either
difference. This gives exactly the short rational corrections
used in (SC3). The finite atom allocations, prescribed
congruences at the two full divisors, amplified wide-cofactor
counts and exact positive product completion apply with their
stated margins. Take \(\omega,\tau\) according to those
margins. The probe refinement changes neither the full-divisor
price nor the order of the limits.
[Small-modulus signed recovery (SC3)](ABC_SMALL_MODULUS_TRIPLE_PRODUCT_COMPLETION.md#signed-words-and-separation-of-weighted-slot-sums)
· [Amplified wide-cofactor count](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)
· [Möbius inversion](https://dlmf.nist.gov/27.5)

## Mixed-sign completion using the three radical caps

Let \(a,b,c\) be fixed positive integers and let
\(\rho,\eta>0\). For sufficiently small fixed \(L>0\),
take any family of positive primitive rows with all three
denominators in \(X/2<C\le X\) and with each original letter's
radical at most \(n^\rho\), where \(n=X^L\). A row is
represented by \(x=A/C\), so the complementary fraction is
\(1-x=B/C\). Consider three pairwise distinct rows \(x,y,z\)
with prescribed exact word values
\[
 x^a y^b z^{-c}=u,\qquad
 (1-x)^a(1-y)^b(1-z)^{-c}=v.                                 \tag{SP1}
\]
Put \(s=(-1)^{a+b-c}\). The count is
\(O_{a,b,c,\rho,\eta,L,\varepsilon}(n^{\rho/2+\eta}X^\varepsilon)\)
uniformly in the prescribed positive rationals \(u,v\).
The same estimate holds when every sign is reversed, by
reciprocating the right-hand sides.

A nonempty fibre bounds the heights of \(u,v\) by
\(X^{a+b+c}\). Let \(P_*\) be the pool of primes dividing
their reduced numerators and denominators and the nonzero
reduced numerators of \(u-1,\ v-1,\ u/v-s\). Its prime
product is polynomial in \(X\) with exponent depending only
on the weights. For any prime outside the pool, the two
valuation coordinates of a row are a nonnegative multiple of
one of the vectors \((1,0),(0,1),(-1,-1)\), for its
\(A,B,C\) roles respectively. The multiple is the original
prime valuation. Equation (SP1) says that the weighted
vectors from \(x,y\) sum to the positive weighted vector
from \(z\). The sum of positive multiples of two different
rays is on none of the three positive rays. It is also
nonzero. Consequently every such prime of \(x\) or \(y\)
occurs on \(z\) in the same role, and an overlap of \(x,y\)
can only use the same role in all three rows.

An overlap in \(A\) makes the complementary fractions units
congruent to \(1\) modulo the prime, forcing \(v\equiv1\).
An overlap in \(B\) gives \(u\equiv1\). On \(C\) use the
row ratios \(x/(1-x)=A/B\) and their analogues. They
are units congruent to \(-1\) because \(A+B=C\).
The corresponding word unit is \(u/v\), so an overlap
forces \(u/v\equiv s\). By construction of the pool, an
overlap is therefore possible only in the respective cases
\(v=1,\ u=1,\ u/v=s\).

The pair \(u=v=1\) gives no distinct solution in the shell
for sufficiently large \(X\). For \(c\le a+b\), weighted
Hölder gives
\[
 z^{c/(a+b)}+(1-z)^{c/(a+b)}
      = (x^a y^b)^{1/(a+b)}
         +((1-x)^a(1-y)^b)^{1/(a+b)}\le1.
\]
If \(c<a+b\), the left side is greater than one. If
\(c=a+b\), equality forces \(x=y\). For \(c>a+b\), the
height \(H\) of a reduced fraction satisfies
\[
 (X/2)^c \le H(z^c)=H(x^a y^b)
                          \le H(x)^a H(y)^b\le X^{a+b},
\]
which is impossible for large \(X\). Thus in a nonempty
fibre at most one of the above three exceptional equalities can
hold. This uses positivity to rule out \(u/v=-1\).

Choose two nonexceptional roles. In each role the supports
of \(x,y\) outside the pool are disjoint and their product
is bounded by the radical of that role of \(z\), at most
\(n^\rho\). Over the two chosen roles the product of all
four outside radicals is at most \(n^{2\rho}\). Thus at
least one of \(x,y\) has pair outside radical product at
most \(n^\rho\). Bin the two individual normalized outside
costs for that row in a fixed mesh sufficiently fine in
terms of \(\eta\). Apply the two-letter dictionary bound
(MK1) with pool \(P_*\) on both roles and the given raw
caps on all three roles. The two binned costs sum to at
most \(\rho\) plus the mesh allowance. Counting over the
bins and over the choice of \(x\) or \(y\) gives the
exponent \(\rho/2+\eta\).

For completeness, after fixing \(x\) put \(U=u/x^a\)
and \(V=v/(1-x)^a\). Then \(z=(y^b/U)^{1/c}\).
On the interval where \(0<y,z<1\),
\[
 \frac{d}{dy}\log\!\left(\frac{(1-y)^b}{(1-z)^c}\right)
       =\frac{b(z-y)}{y(1-y)(1-z)}.
\]
If \(b\ne c\), the ratio \(z/y\) is strictly monotone, so
the derivative vanishes at most once and the level \(V\)
has at most two preimages. If \(b=c\), the ratio \(z/y\)
is constant. A nonzero derivative gives at most one preimage;
the remaining case gives \(y=z\) and is excluded. The same
argument applies after fixing \(y\). Each rational fraction
determines its primitive row, so this finishes the count.
[Two-letter dictionary bound with prescribed pools (MK1)](ABC_MOVING_TARGET_KERNEL_EXTRACTION_BARRIER.md#a-two-letter-count-with-prescribed-prime-pools)
· [Fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For a fixed \(m\)-slot modular word in (FP1) containing
a mixed-sign three-slot pattern, restrict to tuples whose rows
are distinct on those three slots. On a chart of \(T\) rows
satisfying the above caps, fix the other rows, recover the
remaining exact word values through (FP1), and apply (SP1).
The count is
\[
            O\!\left(T^{m-3}n^{\rho/2+\eta}X^\varepsilon\right).
                                                               \tag{SP2}
\]
Split the normalized allowance between the two bounds.
On a balanced chart with \(\rho=1/3+o(1)\) and
\(T=n^{1/3-o(1)}\), this saves \(n^{1/6-o(1)}\) over
fixing \(m-2\) rows and using only the two-row bound.
Here and in this comparison the fixed allowance can be
made arbitrarily small. Positive three-slot patterns already
have the stronger bound in (SC3). The estimate (SP2) also
applies when no three weights agree.
[Positive triple completion](ABC_EXACT_COMPLEMENTARY_PRODUCT_TRIPLE_COMPLETION.md#completion-from-two-exact-products)

For least target repair, multiplication by the prescribed common
factor leaves both projective fractions unchanged. The same
signed count therefore holds for distinct minimally repaired
primitive quotients in each original height shell and chart,
with the original radical caps on their primitive letters.
[Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## Height obstruction and the signed rational-curve example

For weights \((1,1,-1)\), the relation \(u+v=1\) forces
a repeated row. Indeed (SP1) gives
\[
 y\bigl(1-(1+v/u)x\bigr)=1-v-x,\qquad z=xy/u .
\]
When \(u+v=1\), the first equality is equivalent to
\((1-x/u)(y-u)=0\). Either \(x=u\) and \(z=y\), or \(y=u\)
and \(z=x\).

For \(u+v\ne1\), the exact map to the third row is
\[
                   z=\frac{x(1-v-x)}{u-(u+v)x}.
                                                               \tag{SP3}
\]
After homogenizing in \(x=A/C\), the forms
\(F(A,C)=A((1-v)C-A)\) and \(G(A,C)=C(uC-(u+v)A)\)
have degree two and no common projective zero. At \(A/C=0\)
the denominator is nonzero. At \(A/C=1-v\) it equals
\(v(u+v-1)\), and at \(C=0\) the numerator is nonzero.

Let \(Q=\max(H(u),H(v))\). Clearing both coefficient
denominators gives integral binary forms of coefficient
size at most \(2Q^2\). Write \(K=2Q^2\). In the chart
\(C=1\), the denominator is linear with nonzero leading
coefficient. If \(g_0+g_1 t\) denotes this linear form
and \(f(t)\) is the numerator, then
\(R_1=g_1^2 f(-g_0/g_1)\) is a nonzero integer of
absolute value \(O(K^3)\). The identity
\(g_1^2 f(t)-R_1=(g_0+g_1t)h(t)\) has a linear
cofactor of coefficient size \(O(K^2)\). On
\(0<t<1\) it gives
\(\max(|f(t)|,|g_0+g_1t|)\gg K^{-2}\). Homogenize and
multiply the identity by \(C\) to recover the full degree-two
denominator form. The gcd of the two integral values at
\((A,C)\) therefore divides \(R_1 C^3\). In the chart
\(A=1\) the numerator is linear, or a nonzero constant
if \(v=1\). The same argument, or the constant form,
gives divisibility of the gcd by \(R_2 A^3\) for a
nonzero integer \(R_2=O(K^3)\). Since \((A,C)=1\), the
gcd divides \(R_1R_2\). Thus
\[
                         H(z)\gg H(x)^2 Q^{-16}.
\]
In the common shell \(X/2<C\le X\) this implies
\(\max(H(u),H(v))\gg X^{1/16}\) in every distinct
\((1,1,-1)\) fibre. More generally \(c>a+b\) in (SP1)
implies \(H(u),H(v)\ge 2^{-c}X^{c-a-b}\), by applying
the height inequality to \(z^c=x^a y^b/u\) and to its
complementary counterpart.

There is a rational-curve obstruction to a signed extension
of the unrestricted-height, subpower completion. Fix
\(u=v=1/3\) and take any rational \(1/20<x<1/10\).
Then \(y=(2/3-x)/(1-2x)\) and \(z=3xy\) lie in \((0,1)\)
and satisfy \(x<z<y\) and both equations in (SP1).
Their heights are \(O(H(x)^2)\). Taking odd numerators
in the given interval over one power-of-two denominator
of size \(\asymp\sqrt X\), with a sufficiently small
fixed constant, gives \(\gg\sqrt X\) distinct triples
of primitive rows of height at most \(X\). This example
uses neither a shared height shell nor a vanishing radical
budget.

## The remaining image bound

Finite probes and mixed-sign completion control how many
original tuples can have a given record. They leave open an
upper bound on the number of occupied records. For example,
if per-role normalized mask height bounds at probe \(j\)
are \(t_{A,j},t_{B,j},t_{C,j}\), the generic determinant
reconstruction for a threefold positive product uses
modulus exponent greater than both \(3(t_{A,j}+t_{C,j})\)
and \(3(t_{B,j}+t_{C,j})\). Covering the diffuse mass to
within a small loss requires
\(\sum_j(t_{U,j}+t_{C,j})\ge c_U+c_C-o(1)\) for
\(U=A,B\), where \(c_U\) is the nonatomic mass on that
role. Each such pair has mass at least \(1/3\) at the
endpoint. Consequently the product of the two-coordinate
ambient group sizes from this sufficient reconstruction
rule has normalized exponent at least \(2-o(1)\), while
a saturated family requires only \(n^{1-o(1)}\) distinct
threefold products. This comparison concerns the ambient
group estimate. A power saving from this approach requires
additional restrictions on the actual occupied word image
or a separate way to exploit its expansion.
