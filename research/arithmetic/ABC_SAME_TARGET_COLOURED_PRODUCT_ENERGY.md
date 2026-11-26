# Same-target product energy and a coloured collision

The two complementary modular records have few nontrivial
twofold collisions when either input pair shares a target.
For common-target threefold collisions, the detected free
primes give an explicit overlap discount near the largest
target multiplicities. A six-row example with all summands
of full numerical height shows that exact role-coloured
valuation records can still have nontrivial threefold
collisions on targets of arbitrarily large fullness order.

## Twofold collisions with one same-target pair

Use a balanced full-divisor chart as in (FP1) with \(T\)
distinct primitive rows \(A+B=C\le X\), where \(n=X^L\).
Give each row the two-coordinate complementary mask record
\(R\) at the probes of that theorem. Fix the probe allowances
before the small-budget limit. A collision between two
ordered lists is nontrivial if their multisets of primitive
rows differ.

**CE1.** For every fixed \(\eta>0\), at sufficiently small
fixed budgets the number of nontrivial collisions
\[
                       R(x_1)R(x_2)=R(y_1)R(y_2)
\]
in which at least one of the two input pairs shares a
complete target is \(O(Tn^\eta X^\varepsilon)\).
An exact collision with prescribed quotient values \(U,V\)
in place of the modular identity has the corresponding
\(O(TX^\varepsilon)\) bound without any radical caps,
provided the heights of the original rows are at most \(X\).
If the shared target is required to lie in a prescribed set,
\(T\) can be replaced by the number of chart rows on those
targets. The other pair can still use arbitrary chart rows.

To prove the exact assertion, represent rows by their
fractions \(A/C\) and take \(U=x_1x_2/(y_1y_2)\) and
\(V=(1-x_1)(1-x_2)/((1-y_1)(1-y_2))\). If a row
occurs on both sides, cancel it. The two remaining
fractions are determined by their exact ratios \(U,V\).
Indeed if they are \(x,y\), the equations
\(x=Uy\), \(1-x=V(1-y)\) fix them unless \(U=V=1\).
That case is trivial. The cancelled row costs at most
\(T\) choices and lies on the shared target. Orient the
collision so that the positive pair shares the target.
If that pair repeats a row, fix it. The two complementary
products of the other pair determine its sum and product,
hence its unordered multiset. If only the negative pair
repeats, fix one positive row \(x\). The repeated fraction
\(y\) on the negative side must satisfy
\[
          U y^2/x+V(1-y)^2/(1-x)=1 .
\]
The left side is a genuine quadratic, since its leading
coefficient is positive. Each of its at most two roots
determines the remaining positive row. In either case the
freely fixed row lies on the shared target.

Suppose now that all four rows are distinct. Exchange the
two lists if necessary so that the positive pair shares a
target \(C\), and fix its first row. If a prime \(p\nmid C\)
divides the denominator of \(y_1\) or \(y_2\), the contribution
from that negative slot to both word valuations is strictly
positive. At that prime the other negative slot can have
a negative contribution in at most one coordinate, by
primitivity. The positive slots give nonnegative contributions
in both coordinates. Thus at least one of \(v_p(U),v_p(V)\)
is positive. This shows that the two remaining denominator
supports lie in the fixed polynomial-height pool supplied by
\(C,U,V\). The two complete denominators are accordingly
enumerated at subpower cost. Every nonempty exact fibre has
\(H(U),H(V)\le X^4\).

It remains to complete three distinct fractions at prescribed
denominators, possibly different. They obey two equations of
the type
\[
                      x y/z=u,\qquad
               (1-x)(1-y)/(1-z)=v .
\]
The split factorization (CT2) applies here as well. The
denominators only enter it through their known product,
which clears \(z\) and \(x-y\). More explicitly, for
\(t=u+v\) and \(K=1-v\), multiply both linear terms in
\[
             (t^2z+Kt-2u)^2-\bigl(t(x-y)\bigr)^2
                               =4uv(t-1)
\]
by a common multiple of the required denominators of
polynomial height. If \(t\ne1\), factor the prescribed
nonzero integer on the right after this scaling. The
factors recover \(z,x-y\), and the exact sum gives
\(x,y\). If \(t=1\), (CT1) forces \(z=x\) or \(z=y\),
which is excluded. Thus after fixing the first row the
completion costs a subpower. The modular assertion follows
by first applying (FP1) to the four-slot identity and then
using the exact count uniformly for each recovered pair.
[Split factorization (CT1–CT2)](ABC_SAME_TARGET_COMPLEMENTARY_PRODUCTS.md#subpower-completion-for-the-basic-mixed-sign-word)
· [Finite signed word recovery (FP1)](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md#finite-probes-by-möbius-inversion)
· [Uniform fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For instance, let \(D_C\) denote the number of chart rows
with target \(C\). In any subcollection of targets with
\(D_C\ge n^u\) for fixed \(u>0\), all but a power-small
proportion of its ordered same-target pairs have their
modular product unshared by any other unordered pair in
the chart. Indeed \(\sum_C D_C^2\ge n^u\sum_C D_C\)
over the subcollection. Take the exponent allowances in
(CE1) smaller than \(u\).
This conclusion counts the complementary record with two
coordinates. It need not hold for a single role-blind mask
product. For \(C=39\) the two pairs
\((4,35),(17,22)\) and \((5,34),(11,28)\) already
have exactly the same combined multiset of prime and valuation
slots across their free letters.

## A free-support discount for common-target threefold energy

Use the same chart with the original free-letter caps
\(\operatorname{rad}(A),\operatorname{rad}(B)\le n^\rho\).
At the balanced endpoint one may take \(\rho\) arbitrarily
close to \(1/3\). Choose finite probes as in (FP1) so that
each row has detected normalized free radical mass at least
\(1/3-\xi\), for any prescribed \(\xi>0\) and sufficiently
late budgets. Here detection means that the prime occurs
in a nontrivial mask at some probe. This coverage follows
from the lower bound \(1/3\) for the sum of the two free
nonatomic masses and the exceptional-slot allowance in that
construction.

The full-divisor allocation argument gives the stronger
bound \(n^\eta X^\varepsilon\sum_C D_C^{m-1}\) for
any fixed fold \(m\), at every positive multiplicity.
It also allows the two lists to have different
target radicals.
[Full-divisor completion FC1](ABC_FULL_DIVISOR_COMPLEMENTARY_WORD_ENERGY.md#fc1-the-complete-common-target-and-its-radical)

**CE2.** Fix
\(0<\sigma<(1/3-\xi)/6\) with \(\sigma\le\rho\).
In a complete-target bin, the number of nontrivial threefold
collisions of the complementary modular records in which
all six rows have the same complete target is
\[
             O\!\left(n^{(\rho-\sigma)/2+\eta}X^\varepsilon
                                \sum_C D_C^2\right).
                                                               \tag{CE2}
\]
If all bin multiplicities are \(D_C=n^{u+o(1)}\),
this saves on \(\sum_C D_C^3\) whenever
\(u>(\rho-\sigma)/2\) by a fixed margin. In the balanced
limit the saving applies for any fixed \(u>5/36\) with
small enough allowances.

Consider first six distinct rows in such a collision.
Since all six share \(C\), the three masks on \(C\) on
each side cancel in each modular coordinate. At every
probe the products of the three \(A\)-masks agree exactly
as integers, and likewise for \(B\). This follows from
the mask height cap and the prescribed reconstruction
modulus of (FP1). Each detected free prime of a row must
therefore occur in the same role on some row of the
opposite side. For a given row the two roles and the
three partners give six possible classes for this outside
overlap. Its total detected logarithmic mass is at least
\((1/3-\xi)\log n\). One opposite row shares a radical
packet of weight at least \(\sigma\log n\) in one role.

Fix the first member of this pair and one additional
slot. With their target \(C\) now fixed, use (FL1) to
count the other member. Its shared role has outside
radical at most \(n^{\rho-\sigma}\) after pricing the
anchor's prime pool at zero. Its count is at most
\(n^{(\rho-\sigma)/2+\eta}X^\varepsilon\).
The two freely fixed rows cost \(D_C^2\). The three
unfixed slots carry the pattern \((1,1,-1)\) up to
reversing signs, so (CT1) completes them at subpower
cost for each exact pair from (FP1). There are finitely
many choices of overlap roles and positions. Split the
prescribed losses between (FP1), (FL1) and the completion.

For repetitions, first fix the exact output pair supplied
by (FP1). If the two lists share a row, cancel one
copy from each and apply the exact part of (CE1) to the
two remaining pairs on their common target. The residual
collision is nontrivial and the cancelled row costs \(D_C\).
With no common row between signs and exactly five active
rows, fix two slots so that the other three have two equal
positive magnitudes and one negative magnitude, reversing
signs if necessary. Then (CW1) gives a subpower completion.
If there are at most four active rows with no common row
between signs, fix all but two distinct rows on the same
side when such a pair exists. Their exact positive
weighted complementary products have at most two
representations. Indeed \(x^a y^b=u\) makes \(y\) strictly
decreasing in \(x\). The derivative of the logarithm of
the complementary product has the sign of \(y-x\) and
changes sign at most once. If no same-side pair exists,
the active exponents are \(3,-3\); take cube roots of
the two ratios and solve the linear equations as in
(CE1). The repetitions thus cost at most
\(D_C^2 n^\eta X^\varepsilon\) after summing the exact
pairs. This proves (CE2).
[Fixed-letter count with prescribed pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Equal positive signed weights (CW1)](ABC_SAME_TARGET_EQUAL_WEIGHT_WORDS.md#equal-positive-weights-and-arbitrary-negative-magnitude)

## An all-tall coloured threefold collision

**CE3.** There are infinitely many complete targets \(C\)
of arbitrarily high fixed fullness order with six distinct
positive primitive rows \(A+B=C\), all with \(A,B\asymp C\),
that give a nontrivial threefold collision of every
role-coloured valuation mask. The collision holds for the
complete \((p,e)\)-slot records, with no modular reduction.
Within each side of the collision the free prime supports
of different rows are disjoint outside one fixed finite
prime pool.

Take integers \(j\ge0\) with \(j\equiv0\bmod210\) and put
\[
\begin{aligned}
 w&=14+90j,& h&=1+9j,& s&=2w,\\
 v&=43+275j,& k&=32+204j,& d&=41+261j,& q&=29+185j.
\end{aligned}
\]
The complete target is \(C=289+1845j\). The left and right
lists are
\[
\begin{array}{c|cc|cc}
  &A_{\rm left}&B_{\rm left}&A_{\rm right}&B_{\rm right}\\ \hline
  1&h&9k&2h&7d\\
  2&3s&5d&s&9q\\
  3&2v&7q&3v&5k .
\end{array}
\]
Every row sums to \(C\), by the linear formulas. The six
summands in role \(A\) have distinct slopes
\(9,540,550,18,180,825\), all strictly less than \(1845/2\).
Their constant terms have the same order as these slopes.
Each is therefore the smaller summand of its row, and the
rows are distinct even after allowing summand exchange.
The positive linear coefficients of the smaller summands
give \(A,B\asymp C\) uniformly as \(j\) grows.

For primitivity, if \(A=a_0+a_1j\), any prime common to
\(A,C\) divides \(1845a_0-289a_1\). In the left and right
row orders the six determinants are respectively
\[
             -756,\ -1080,\ -280,\ -1512,\ -360,\ -420 .
\]
All their prime divisors are in \(\{2,3,5,7\}\). The base
rows at \(j=0\) are primitive, and \(j\equiv0\bmod210\)
preserves their residues at all these primes. This proves
primitivity of every displayed row.

Role \(A\) uses the identical list of multiplicative atoms
\(h,2,3,s,v\) on both sides. Role \(B\) uses the identical
list \(9,k,5,d,7,q\). Each displayed pairing of atoms
within a letter is coprime. Every such pairing has a fixed
factor supported on \(2,3,5,7\), and the required
coprimalities follow at \(j=0\) and persist modulo \(210\).
Consequently each role has exactly the same multiset of
prime and valuation slots on the two sides. The target role
has three copies of \(C\) on each side. In particular
separate masks on all three roles, or their complementary
ratios, give this exact collision for any choice of probes.
The nonconstant atoms \(h,w,v,k,d,q\) are nonproportional
linear polynomials in \(j\). For two distinct such atoms,
a prime of their common support divides their fixed
nonzero coefficient determinant. Within either side
each nonconstant atom appears on only one row. Adjoining
the primes of all these determinants and \(2,3,5,7\)
gives the required fixed overlap pool.

Finally set \(M=1845\cdot210=387450\).
Since \(17\) is a unit modulo \(M\), for every fixed
\(r\ge2\) with \(r\equiv2\bmod\varphi(M)\), every positive
integer \(K\equiv17\bmod M\) has
\(K^r\equiv289\bmod M\). Hence \(C=K^r\) lies in the
displayed target progression. One may take \(K\) squarefree
and obtain \(\gg_r X^{1/r}\) such examples in each
sufficiently large shell \(X/2<C\le X\). In fact Möbius
inversion counts squarefree integers in this fixed unit
progression with positive density, by summing over square
divisors and using the Chinese remainder theorem. For these
choices \(\operatorname{rad}(C)=C^{1/r}\).

The overlap condition in this example is within each
side of a collision. Cross-side repetition of detected
primes is necessary and occurs here. The example makes
no claim for the two small free-letter radical caps.
It shows why internal disjointness, common exact-power
targets and numerical height alone do not force diagonal
coloured threefold energy. Conversely (CE1) and (CE2)
bound the multiplicity of occupied complementary records.
An upper bound for the actual same-target scalar image
with both free-letter radical restrictions is still needed
for the global exceptional-set problem.
[Same-target scalar mask image (SR3)](ABC_SCALAR_MASK_TARGET_MOMENTS.md#a-scalar-image-for-weighted-target-occurrence)
