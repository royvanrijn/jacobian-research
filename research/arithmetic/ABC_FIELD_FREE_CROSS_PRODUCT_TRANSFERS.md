# Field-free product and comparison transfers for primitive pairs

The determinant of two products of actual scaled operands
is the product of two reduced shifts. Its primitive
additive triple can recover the original pair by a
split quartic count. The only additional record required
in this recovery is the common parity kernel that
disappears from the auxiliary triple. This record is
cheap on typical pairs of an original near-full family.

## An ordinary primitive triple from the product difference

**FP1.** Take two distinct positive primitive rows
\(a_i+b_i=c_i\le X\). Exclude pairs obtained by swapping
the two summands of the same row. Set
\[
 G=(c_1,c_2),\quad u=c_2/G,\quad v=c_1/G,\quad
 (A,B)=(ua_1,ub_1),\quad (C,D)=(va_2,vb_2).
\]
The common sum is \(W=A+B=C+D\).
Write \(h=A-C,\ j=A-D\). Neither is zero, because
an equality would identify or interchange the reduced
primitive proportions. The identity
\[
                 AB-CD=-h j
\]
gives the auxiliary triple with positive values
\[
                  P=AB/d,\quad Q=CD/d,\quad |P-Q|,
                     \qquad d=(AB,CD).
\]
The three absolute values are put in additive order.
Here \(d=(a_1b_1,a_2b_2)\), since \(u,v\) are coprime
and every prime of \(u\) is prime to \(va_2b_2\),
with the reverse assertion for \(v\). Thus \(d\)
is prime to \(uv\).

Assume an actual common divisor \(M\mid G\) has
\(M\ge X^\nu\) for fixed \(\nu>0\) and
\(a_1\equiv a_2\bmod M\). This is an ordinary absolute
residue cell with a shared target lock. The auxiliary
has fixed positive numerical height unless
\(a_1b_1=a_2b_2\). In fact \(u^2\mid P,\ v^2\mid Q\),
\(d\) is a unit at \(M\), and
\[
   M\mid v^2P-u^2Q,\qquad
              |v^2P-u^2Q|\le P Q .
\]
The divisibility follows from both original products
\(a_i b_i\equiv-a_i^2\bmod M\).
If they are unequal the integer on the left is
nonzero, and \(\max(P,Q)\ge\sqrt{M}\).
The equality case consists of at most
\(N X^\varepsilon\) ordered pairs in any original
\(N\)-row family, because fixing \(a_1b_1\) leaves
only divisor-many possibilities for \(a_2,b_2\).
The normalized-lock alternatives \(M\mid h\) or
\(M\mid j\) instead make \(P-Q\) divisible by \(M\)
when the scaled operands are units there. For a
common target lock, this unit condition follows
from that shift divisibility and primitivity.

The raw radical mass of the product auxiliary
can be read exactly, including possible prescribed
zero-price support. Write \(\ell(Z)\) for the
normalized radical of the absolute nonzero integer
\(Z\) outside one fixed pool. Put
\[
 \delta=\ell(d),\qquad
 \delta_0=\log_n\operatorname{rad}_{S^c}(d),
 \quad S=\operatorname{supp}(P Q(P-Q))\cup P_0,
\]
where \(P_0\) is the chosen prime pool and
\(\operatorname{rad}_{S^c}\) counts only primes
outside \(S\). Then
\[
 \ell(P Q(P-Q))=
     \ell(AB)+\ell(CD)+\ell(hj)-2\delta-\delta_0 .
                                                               \tag{FP1}
\]
For each prime of \(d\), all three raw values
\(AB,CD,hj\) contain the prime. Normalization
leaves it on exactly one auxiliary role or
removes it altogether. This gives the two
discounts. All other primes have the same
incidence before and after division.

## Subpower recovery from one auxiliary and a parity overlap

**FP2.** With \(P,Q\) prescribed as actual exact auxiliary
values, write \(D_i=\operatorname{sf}_2(a_i b_i)\), the
squarefree kernel of the original product, and
\(D_0=(D_1,D_2)\). Fix a polynomial-height prime pool
\(P_0\). If
\[
     \operatorname{rad}_{(\operatorname{supp}(P Q(P-Q))
                       \,\cup\,P_0)^c}(D_0)\le n^\zeta ,
\]
the original pair count is
\(O(n^\zeta X^\varepsilon)\) at any fixed budget.
The two complete original targets may both vary.
All original caps can be retained in the recovered pairs.

Indeed \(\operatorname{sf}_2(P Q)=\operatorname{sf}_2(D_1D_2)\)
because the scaling and normalization factors in the
product form a rational square. Hence \(P,Q\) give the
symmetric difference of the two parity supports.
The primes on that difference can be allocated to
the two fields at subpower cost. Choose \(D_0\)
outside the prescribed support at cost \(n^\zeta\),
then choose its remaining prime factors on that
support at height-subpower cost. This fixes \(D_1\).

Write \(r=A/C=s/t\) in lowest terms with \(s,t>0\),
and \(x_1=A/W,\ x_2=C/W\). From the auxiliary ratio,
\[
 x_2=\frac{P-Qr}{P-Qr^2},\qquad
 x_1=\frac{s(Pt-Qs)}{Pt^2-Qs^2},
 \qquad
 x_1(1-x_1)=\frac{Pst(Pt-Qs)(t-s)}
                         {(Pt^2-Qs^2)^2}.
\]
The denominator in the first expression cannot
vanish on a genuine pair: if it vanishes, the
numerator vanishes as well, and the two identities
would force \(r=1\). This is excluded by the
nonzero displacement \(h\). All four linear
factors of the numerator have distinct rational
projective zeros since \(P,Q\) are positive
and \(P\ne Q\). The squareclass of
\(x_1(1-x_1)\) is \(D_1\). Thus the
projective direction \((s:t)\) occurs on the
split quartic double cover of the line
\[
                  D_1 z^2 = P s t(Pt-Qs)(t-s).
\]
The affine equation is understood projectively,
with the square coordinate of weight two.
Move a branch point to infinity to obtain a
cubic with three distinct rational roots,
then scale the leading coefficient. It is an
elliptic equation with rational \(2\)-torsion.
The transformation uses only a fixed number of
rational operations and the input degrees are
fixed. Both the equation and the mapped rational
points therefore have polynomial height in \(X\),
uniformly in \(P,Q,D_1,s,t\).
Dujella's rational-torsion bounded-height theorem,
as used for MTK1, now gives \(X^{o_X(1)}\) directions.
Each exact direction determines \(x_1,x_2\)
and hence both original primitive rows. Summing
over the chosen parity record gives the bound.
[Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/html/2312.03655v1)
· [Four distinct linear forms (MTK1)](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#the-mixed-kernel-count)

In particular an outside radical cap \(n^q\) on
the actual product auxiliary in FP1, at fixed
numerical depth, gives original pair exponent
\[
                             q/3+\zeta+\eta .
                                                               \tag{FP2}
\]
For the auxiliary count, fix its complete powers
on the prescribed pool at height-subpower cost.
These powers give arbitrary fixed coefficients.
Use exact-order singleton layers and the
two-cheapest-letter bound MCB2 on the remaining
radicals. Their three normalized masses sum to
at most \(q\) because the primitive auxiliary
roles are pairwise coprime. Half of the two
cheapest costs at most \(q/3\). FP1 gives
numerical height at least \(X^{\nu/2}\), above
the MCB2 threshold at small fixed budgets.
This argument only uses an outside radical cap
and permits prescribed coefficients of
polynomial height. A positive additive ordering
of the auxiliary leaves boundedly many possible
roles for the prescribed \(P,Q\).
[Arbitrary-coefficient two-cheapest count (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Effective radical budgets (EB4–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

The parity record is cheap on typical pairs.
Suppose all three original cap exponents are at
most \(1/3+\tau\). For any fixed first row, a
second row sharing at least \(\sigma\) of
normalized radical with the first row's
fixed complete product support has count
\(O(n^{1/3+\tau-\sigma/3+\eta}X^\varepsilon)\)
by EB4 and EB5. Thus power-full pairs of a
near-full family can be restricted to
\(\zeta\le\sigma\) provided their pair loss
is below this saving, after the actual
chart and cap tolerances. The omitted common
parity mass can be much smaller than either
individual parity field price. If the pairs
have equal targets and the fixed-target
multiplicity approaches \(1/6\), the same
conclusion follows fibre by fibre: sharing
\(\sigma\) on the two free operands forces
one operand to spend at least \(\sigma/2\)
on the first row's pool, so FL1 bounds such
second rows by \(n^{1/6+\tau/2-\sigma/4+\eta}\).
For a chart of exponent \(d_0\) and pair
exponent \(2d_0-\rho\), the sufficient
deletion margin is
\(\sigma/3>1/3-d_0+\tau+\rho\).
For an equal-target selected exponent
\(d_0+u_0-\rho\), it is
\(\sigma/4>1/6-u_0+\tau/2+\rho\).
The strict margins absorb the height allowances.
[Fixed-letter count with prescribed pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Thus on an exponent-full pair sector in a
balanced chart of exponent \(d_0\) near
\(1/3\), at actual lock price and chosen
pair loss \(\rho\), the cross-product auxiliary
needs outside radical cap at least
\(3(2d_0-\rho-\zeta)-o(1)\) if the inequalities
are to coexist. For full-scale generic pairs
on arbitrarily cheap charts this approaches
the two-row cap \(2\). At equal target and
complete multiplicity \(u_0\), a selected
pair exponent \(d_0+u_0-\rho\) instead
requires \(3(d_0+u_0-\rho-\zeta)\).
These are necessary conditions on the
actual product differences, not assumptions
about the radical of an arbitrary new shift.

## Two other field-free recovery records

**FP3.** If the two original complete targets
are independently prescribed, one cross-product
auxiliary can replace the quartic recovery.
Normalize the pair \((a_1b_2,a_2b_1)\) by its
gcd and add the positive difference as the
third entry. Its exact ratio together with
the prescribed \(c_1,c_2\) recovers both
original rows at divisor cost by TR1. With
moving targets, use both \((a_1b_2,a_2b_1)\)
and \((a_1a_2,b_1b_2)\) as the two
normalized comparisons. Their exact ratios
\(\rho,\sigma\) recover the original odds
from the positive square roots of
\(\rho\sigma\) and \(\sigma/\rho\).
Every such comparison is numerically
tall in a common absolute residue cell:
the cross-products agree modulo \(M\),
and their gcd is a unit there. Their
nonzero reduced difference is therefore
a multiple of \(M\).
[Exact pair recovery from a target ratio (TR1)](ABC_DIFFUSE_FACTOR_GRAPH.md#recovery-from-a-matched-target-ratio)

A different second auxiliary is the primitive
comparison of two scaled endpoints, for
example \(A,C,|h|\) divided by
\(\gcd(A,C)\). Its ratio is the \(r\) used in
FP2, so fixing this comparison and the
product auxiliary recovers the two rows
directly without a parity record. An
opposite-role comparison uses \(A,D\) and
then interchanges \(a_2,b_2\) during the
recovery. An absolute same-residue lock
suffices for numerical depth even when
the scale factors are not units mod \(M\).
For example, for normalized endpoints
\(S=A/g_-,T=C/g_-\) with \(g_-=(A,C)\)
one has \(u\mid S,v\mid T\), and
\(M\mid vS-uT\). If \(a_1\ne a_2\),
this nonzero difference is at most \(ST\)
in absolute value. Use the other
same-role comparison if \(a_1=a_2\)
and \(b_1\ne b_2\).
For opposite roles, use instead the
positive \(vS+uT\) divisible by \(M\)
and at most \(2ST\).
The two ratios here can be recovered
from their exact auxiliary additive
triples up to a bounded choice of order.
The product auxiliary must retain the
height floor of FP1. Under the absolute
lock handle equal original operand products
at the separate divisor bound of FP1.
Apply the same MCB2 count used in FP2 to
both actual auxiliary triples with their
own caps and prescribed pools, giving
one-third of the sum of the two
actual auxiliary radical prices.
This alternative trades the parity
overlap price in FP2 for an additional
ordinary radical price. Both recoveries
allow all original caps and genuine
least-repair eligibility.

## A boundary for three-term identities of the same records

**FP4.** In three independent variables
\(B,h,j\), let the seven linear forms be
\[
 B,\quad A=B+h+j,\quad C=B+j,\quad D=B+h,\quad
 h,\quad j,\quad W=2B+h+j .
\]
Among homogeneous polynomial identities given by
a nonzero linear combination of three monomials in
these records, with no common factor and using
both \(h,j\), the only one is
\(AB-CD=-hj\), up to scaling and reordering.
All monomials have nonnegative integer exponents.

Indeed a factor in two monomials must divide
the third, so the supports are pairwise
disjoint, by irreducibility of the seven
distinct linear forms. If \(h,j\) are in
different terms, reduce modulo \(h\).
The remaining term with \(j\) cannot equal
the one without it: at \(h=0\) none of the
other remaining forms is proportional to \(j\).
Hence both are in one term. Reducing at
\(h=0\) requires equal exponent totals in
the other two terms on each class
\(\{A,C\},\{B,D\},\{W\}\). Reducing at
\(j=0\) gives equal totals on
\(\{A,D\},\{B,C\},\{W\}\). Their exponent
difference is consequently
\((k,k,-k,-k,0)\) in the order \(A,B,C,D,W\).
Using disjoint supports and interchanging
terms if needed, the two monomials are
\((AB)^k\) and \((CD)^k\) for \(k>0\),
with opposite coefficients. Their difference is
\[
         -h j\sum_{r=0}^{k-1}(AB)^r(CD)^{k-1-r}.
\]
The sum is nonzero on both \(h=0\) and \(j=0\),
so the two shift exponents of the third
monomial must both be one. Its only other
available factor is \(W\), to exponent
\(2k-2\). At \(h=0\) the sum is
\(k[B(B+j)]^{k-1}\), which is proportional to
\((2B+j)^{2k-2}\) only for \(k=1\).
This classifies the proposed monomial
transfers. The common-field square identities
use square roots of the cross-products, an
additional algebraic structure outside this
record algebra.
[Common-field auxiliary triples (CP4)](ABC_CROSS_PHASE_AND_AUXILIARY_TRIPLES.md#two-auxiliary-abc-triples-in-the-common-field-equal-target-case)
