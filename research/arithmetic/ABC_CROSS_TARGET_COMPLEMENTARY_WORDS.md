# Exact complementary data across different primitive targets

Two locked original roles give complementary congruences whose
exact outputs recover both primitive rows, including when the
targets differ. With only one lock, the reduced odds quotient
identifies which operand support can be charged to the pair.
At higher orders, exact complementary words force common
prime packets of specified signs and roles. Reserving two
or three rows for exact completion saves a fixed power.

## Close coefficient pairs across distinct complete targets

**CTO1.** Take \(n=X^L\), with the height limit at each fixed budget.
Consider an arbitrary family of genuine primitive rows
\(A+B=C\) of height between \(X^{h_0}\) and \(X^B\),
with \(h_0,B>0\) fixed and all three original radical caps
bounded in base \(n\). Fix a prescribed divisor \(D\)
of one original role, say \(C\), on the chart. Assume
\(D\ge X^w\) and \(Le\ge\sigma\) at every \(p^e\parallel D\)
for fixed \(w,\sigma>0\). A counted pair of distinct rows
has a full relative odds quotient
\[
             E=\frac{A_i/B_i}{A_j/B_j}= r^k z,\qquad
   H(r)\le X^{K_r/k},\qquad H(z)\le X^{c\min(1,kL)} ,
\]
where \(r,z\) are positive rational units at \(D\).
Here \(K_r>0\) is fixed and the order \(k\ge k_0\) may be
chosen by the pair. For instance this holds when the two
rows have an eligible common-order full coefficient ratio
with small relative height and each original root ratio
has height at most \(X^{K_r/(2k)}\). Original prime support
on \(A_i B_i A_j B_j\) ensures the needed unit property.

For every fixed \(\eta>0\), choose \(c>0\) sufficiently small
and \(k_0\) sufficiently large in terms of these fixed data.
At all sufficiently late fixed budgets the number of
distinct exact outputs \(E\) is \(O(n^\eta X^\varepsilon)\),
uniformly in the chart and in \(D\). In fact \(E\equiv1\bmod
D\) by primitivity and the two original equations. Use the
same rational output count as FRP3 and FRP4. For the finite
list \(k\le J/L\) the amplified congruence proof applies.
For larger orders record \(r\) at cost \(O(n^{2K_r/J})\).
Fixing \(r,k\) gives a single residue for \(z\), which
determines it at height \(X^c\) if \(2c<w\). If \(r\ne1\)
then \(k\le K_r\log_2 X\). If \(r=1\) this height bound
recovers \(z=1\) directly at sufficiently large \(X\),
and \(E=1\) cannot occur for distinct primitive rows.
Choose \(J\) so large that the enumeration is within the
desired allowance.
[Amplified exact quotient count](ABC_SIGNED_VALUATION_PHASES.md#a-rational-congruence-bound-for-distinct-exact-outputs)
· [Deep-modulus pair recovery (FRP4)](ABC_FIXED_FIBRE_RELATIVE_ODDS_PAIRS.md#a-doubled-full-modulus-for-two-moving-queries)

### A one-lock overlap bound

Fix an exact \(E=P/Q\) in lowest terms and let \(\mathcal P\)
be the support of \(PQ\). Each prime outside \(\mathcal P\)
has the same signed valuation on \(A_i/B_i\) as on \(A_j/B_j\).
Original pairwise coprimality forces the prime to occur on
the same one of the two free roles with the same exponent,
or to be absent from both. Thus prescribing the pool
\(\mathcal P\) leaves only the same-role shared support
outside it. Equivalently, put \(A_i=g_A a_i\), \(A_j=g_A a_j\)
and \(B_i=g_B b_i\), \(B_j=g_B b_j\) with the two same-role
gcds removed. The cross products give \(P=a_i b_j\),
\(Q=a_j b_i\) exactly. Indeed every prime of a factor of
\(a_i b_j\) is excluded from both factors of \(a_j b_i\) by
one of the original coprimality conditions or one of the two
gcd divisions. The four reduced multipliers can also be
recovered at divisor cost.

Suppose the logarithmic radical caps on the outside common
supports on \(A\) and \(B\) are \(s_A,s_B\) in base \(n\).
Let \(\ell_C\) be any original outside cap after the known
\(D\)-support has been removed from \(C\). If \(D\) has
radical exponent \(f\) and the original cap on \(C\) is \(t_C\),
one may take \(\ell_C=t_C-f\). Fix the supported powers at
height-subpower cost. The original two-cheapest-letter
bound counts the first row with exponent at most
\[
                  \frac12\min(s_A+s_B,\ s_A+\ell_C,\ s_B+\ell_C)
\]
up to any fixed positive normalized allowance. The exact
odds \(A_i/B_i\) and \(E\) recover the odds of the second
row, and primitivity recovers that complete row. Summing
the subpower exact quotient allowance gives the same bound
for all the close pairs with these outside overlaps. It is
uniform in the full targets of the two rows, with no
target dictionary charge. If only the sum \(s_A+s_B\le\rho\)
is capped, the optimum upper bound from these comparisons is
\[
             \min(\rho/2,\ \rho/4+\ell_C/2).
\]
In particular, at any small prescribed outside-overlap
threshold the close pairs have arbitrarily small exponent.
[Arbitrary-coefficient two-cheapest bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

### A reduced-odds radical cost at finite free valuation depth

There is a second one-lock count that does not
require the outside overlaps to be small. For a
fixed output \(E=P/Q\), write
\[
           t=\log_n\operatorname{rad}(PQ).
\]
Let \(r_A,r_B\) be the original two free-role
caps in base \(n\). Each prime counted by \(t\)
occurs on a free operand in at least one of the
two rows. At least one row therefore has
combined free radical mass at least \(t/2\)
on the prescribed support of \(PQ\). Fix the
supported powers in that row at height-subpower
cost and use these two roles in MCB2. This
counts the row at exponent at most
\((r_A+r_B)/2-t/4\). Its exact primitive ratio
and \(E\) recover the other row. Thus for
outputs with \(t\ge t_0\), the total
one-lock close-pair exponent is at most
\[
                      (r_A+r_B)/2-t_0/4+\eta .
\]
This estimate retains the target's original
cap and applies in either direction.

In particular suppose every prime on the two
free roles has scaled valuation at most
\(S_0\), with \(S_0>0\) fixed before the budget
limit. Then the height of the reduced \(E\)
is at most \(X^{2S_0t}\), since at most two
free operands contribute to its positive or
negative valuation at one prime. On the other
hand \(E\ne1\) and \(E\equiv1\bmod D\)
give \(|P-Q|\ge D\). If \(D\ge X^w\), then
\(t\ge w/(2S_0)\), so all these close pairs
have exponent at most
\[
                       (r_A+r_B)/2-w/(8S_0)+\eta .
\]
If a common role coloring gives disjoint
prime supports on \(A,B\) across the chart,
each prime contributes to at most one free
role, and \(w/(4S_0)\) is available in place
of \(w/(8S_0)\). This last count is within
the given role coloring. In the presence
of unrestricted high-order tails the
exponent estimate with \(t_0\) still applies,
but the depth of \(D\) need not give a fixed
lower bound on \(t\).

### Two complementary locks recover all eligible pairs

**CTO2.** Suppose the chart instead fixes deep divisors \(D_A\mid A\)
and \(D_B\mid B\) with positive numerical weight and positive
scaled valuation floors on those divisors. For both relative
quotients
\[
             E_1=\frac{A_i/C_i}{A_j/C_j},\qquad
             E_2=\frac{B_i/C_i}{B_j/C_j},
\]
assume the same kind of full narrow-root and relative
coefficient-height representations. The two orders need
not equal each other, and can vary by pair. Require the
variables in the first representation to be units at \(D_B\)
and those in the second at \(D_A\), as automatic when they
use only their original two free roles. Then \(E_1\equiv1\bmod
D_B\) and \(E_2\equiv1\bmod D_A\). The exact output count
applies to both. For each pair \((E_1,E_2)\) solve
\[
           \frac{A_i}{C_i}=E_1\frac{A_j}{C_j},
           \qquad 1-\frac{A_i}{C_i}
                      =E_2\left(1-\frac{A_j}{C_j}\right).
\]
If \(E_1=E_2\), both equal one by addition, so the rows
coincide. Otherwise these two linear equations determine both
reduced fractions and hence both original rows. The total
count of eligible distinct ordered pairs is \(O(n^\eta
X^\varepsilon)\) for any fixed \(\eta,\varepsilon>0\), after
using sufficiently small width and budgets and sufficiently
large fixed order cutoff. Permuting the roles gives the
same conclusion for any two different locks: the two
quotients then use the remaining role as denominator and the
two recovered ratios are related by the actual signed
additive equation. A wholly atomic profile on all but one
role is excluded by the balanced endpoint's pair nonatomic
lower bounds, so there are always two roles eligible for
cheap positive-depth locks in the usual limiting chart
construction.

This tests a proposed occurrence extraction without either
a complete-target catalogue or any assumption of target
multiplicity. To use it to exclude a limiting full family,
the open task is to force a positive exponent of pairs
with jointly close full coefficients on the two
complementary quotients, after fixing the actual two locks.
[Saturation profile and affordable full divisors](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#a-one-grid-estimate-in-a-genuine-fixed-target-fibre)
· [Two complete complementary products on probe transcripts](ABC_SMALL_MODULUS_TRIPLE_PRODUCT_COMPLETION.md#recovering-the-two-exact-product-ratios)

## Exact signed complementary words with all three caps

**SWP1.** Work at \(n=X^L\) on genuine primitive triples \(A+B=C\)
in \(X/2<C\le X\) and with
\(\operatorname{rad}(ABC)\le n^{B_*}\), for fixed \(B_*>0\).
The height limit precedes the decreasing-budget limit.
Write \(x_i=A_i/C_i\). Fix \(m\ge2\), nonzero integer
weights \(w_1,\ldots,w_m\), and positive rational targets
\(u,v\) for the signed complementary products
\[
                \prod_{i=1}^m x_i^{w_i}=u,\qquad
                \prod_{i=1}^m(1-x_i)^{w_i}=v .
\]
Assume all \(m\) rows are distinct. There are at most two
ordered solutions for \(m=2\). Let \(a,b\) count the
positive and negative weights, respectively, so \(a+b=m\).
For \(m\ge3\), the count at sufficiently late fixed
budgets is
\[
 O\left(n^{mB_*/6-mB_*/(6\max(3,ab))+\eta}
                         X^\varepsilon\right)
\]
for any fixed \(\eta,\varepsilon>0\), uniformly in the
two exact targets. If \(ab=0\) the exponent is \(mB_*/9\).

To see the two-row recovery with arbitrary weights
\(p,q\ne0\), prescribe
\(x^p y^q=U\) and \((1-x)^p(1-y)^q=V\) for \(0<x,y<1\).
Then \(y=(U/x^p)^{1/q}\) on its admissible \(x\)-interval,
and the derivative of the logarithm of the second product is
\[
                 \frac{p(y-x)}{x(1-x)(1-y)} .
\]
The equation \(y=x\) has at most one solution unless
\(p+q=0,\ U=1\), when every solution has \(y=x\).
In the first case there are at most two solutions,
and the exceptional case is excluded by distinctness.
The exact fractions recover the primitive rows.

For the price bound, use a prescribed pool \(P\) of
all primes in the numerators and denominators of \(u,v\).
Its product has polynomial height whenever a tuple exists.
For any prime outside \(P\), both target valuations
are zero. A row's contribution has one of the three
vectors \((e,0),(0,e),(-e,-e)\) times its nonzero
signed weight, according as the prime lies on \(A,B,C\).
Original primitivity guarantees only one role per row.
A prime outside \(P\) cannot occur in just one row,
and one occurring in exactly two rows occurs on the same
role and with opposite signs of weights. This uses the
linear independence of every two of the three directions.
This constraint distinguishes the signed complementary
words from a single multiplicative relation.

Let \(b_i=\log_n\operatorname{rad}(A_iB_iC_i)\) and
\(D_0=mB_*-\sum_i b_i\). Put
\[
 S=\sum_{p\in P}m_p\log_n p+
   \sum_{p\notin P,\;m_p\ge3}(m_p-2)\log_n p ,
\]
where \(m_p\) is the number of rows containing \(p\).
Both \(S,D_0\) are nonnegative. In a chosen exposure order,
charge each outside prime on its first row. The total
new mass is \((mB_*-D_0-S)/2\). For a fixed initial
part of the exposure, the previous supports and \(P\)
give a prescribed polynomial-height pool. By setting
its prime prices to zero in EB4, the next row with
new mass at most \(t\) costs at most \(n^{t/3+\eta'}\)
with any positive counting allowance, after the
arbitrary inner allowance. One may pad the budget \(B_*L\)
by a fixed small multiple of \(L\) to use the height shell.
Use fixed finite meshes for all masses before counting.

The incidence of primes outside \(P\) of multiplicity
at least three, together with the pool incidence, is
at most \(3S\). The remaining primes occur on one of
the \(ab\) opposite-sign pairs. If \(Q_*\) denotes the
largest exclusive two-row packet, then
\[
                        mB_*\le D_0+3S+2abQ_* .
\]
Place a maximizing pair last and recover it with the
two-row argument after fixing the earlier rows. This
saves a priced mass of at least \(Q_*/3\). If there is
no two-row packet, any pair can be placed last and
take \(Q_*=0\). Thus the count exponent is at most
\[
        \frac{mB_*-D_0-S}{6}-\frac{Q_*}{3}+\eta .
\]
Since
\(D_0+3S+2abQ_*\le\max(3,ab)(D_0+S+2Q_*)\),
this is the asserted exponent. All fixed pools and
exposure orders are uniform in the outputs; a bin
for the finitely many parameters and a choice of
maximizing pair do not change the exponent.
[Prescribed-support pricing (EB4–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Additive reciprocal exposure and pair recovery (JR3)](ABC_JOINT_RECIPROCAL_PAIR_RECOVERY.md#sequential-pricing-and-a-saving-at-every-fixed-order)

### Equal positive words and the exclusive triple packet

**SWP2.** When all weights are equal and have the same sign, take
the corresponding root of \(u,v\) and use weight one.
A prime outside the target pool then occurs on at least
one row on each of the three roles. Indeed the signed
valuation sums on the three roles are equal, and if
some row contains the prime these positive sums are
nonzero. Consequently the multiplicity is at least
three. With \(S\) defined as above using \(m_p-3\) in
place of \(m_p-2\), total new mass is
\((mB_*-D_0-S)/3\). For \(m\ge4\), reserve a triple of
rows with maximum exclusive outside packet \(Q_*\).
The incidence of all higher-multiplicity primes and
the pool is at most \(4S\), and
\[
                       mB_*\le D_0+4S+3\binom m3 Q_*.
\]
Recover the three reserved rows from their two
complementary products by EC1, at height-subpower
cost. The exponent becomes at most
\[
        \frac{mB_*}{9}-\frac{mB_*}{9\binom m3}+\eta ,
                          \qquad m\ge4,
\]
since \(\binom m3\ge4\). For \(m=3\) the three-row
completion is already subpower. The new exponent
improves the bound \((m-3)B_*/3\) from fixing all but
three rows when \(m\ge5\).
[Exact complementary product completion (EC1)](ABC_EXACT_COMPLEMENTARY_PRODUCT_TRIPLE_COMPLETION.md#completion-from-two-exact-products)

### Consequences for mixed words and the lock test

For nontrivial equal positive-product collisions of
order \(k\), merge repeated rows and cancel rows which
occur equally often on the two sides. Each canceled
row costs at most \(n^{B_*/3+\eta}\), and uses at
least two slots. The nonzero effective weights are
fixed after one of finitely many multiplicity patterns.
The total surviving signed weight is zero. One active
row is therefore impossible, and two active distinct
rows cannot have identity output. For at least three
active rows the signed bound, including the canceled
rows' cost, saves a fixed power from \(kB_*/3\) at
original half-order \(k\). Each nonempty
active pattern saves at least
\(\delta_k=B_*/(3\max(3,k^2))\) from this exponent,
since the number \(m\) of active rows is at least
two and \(ab\le k^2\). The factors from freely
canceled rows do not consume this margin. If the
original family has
exponent \(B_*/3\), the equal positive-product energy
is thus asymptotically diagonal at the exponent scale.
The same estimates for distinct rows apply to all
fixed signed weights in a polynomial-height
two-coordinate target catalogue, with its actual
dictionary exponent added.

In particular, on a genuine two-lock chart with
\(D_A\mid A\), \(D_B\mid B\), the full words in the two
complete ratios have residue one at \(D_B,D_A\)
respectively by original addition and primitivity.
If their narrow-root and short-remainder hypotheses
are those of the exact quotient count, that count
makes the two-output catalogue exponent arbitrarily
small. For \(m\) distinct slots the preceding
saving then applies to all tuples meeting the two
word hypotheses, whether or not the complete
targets agree. This is a quantitative obstruction
to proposed compressed signed word configurations.
It retains the full three-role radical caps.
[Two-base exact unit output count](ABC_SIGNED_VALUATION_PHASES.md#a-rational-congruence-bound-for-distinct-exact-outputs)
