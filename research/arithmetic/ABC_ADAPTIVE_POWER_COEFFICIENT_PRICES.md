# Conditional power coefficients and an adaptive two-cap obstruction

A fixed target and a full divisor on one free role make an
adaptive power on the other role inexpensive after its
coefficient is known. With two variable powered roles, a
small separating prime pool for their joint coefficient
ratio can replace a small coefficient dictionary. A
primitive two-cap matching shows why an economical joint
dictionary cannot be extracted from only the two powered
roles, even when orders may adapt to the row.

## A one-grid estimate in a genuine fixed-target fibre

Use \(n=X^L\) with \(L\) fixed before the height limit.
Consider actual positive primitive rows of height at
most \(X\) with a fixed complete value \(T\) in any
original role, written
\[
          \epsilon_U U+\epsilon_V V+\epsilon_T T=0,
                       \qquad \epsilon_U,\epsilon_V,\epsilon_T\in\{-1,1\}.
\]
Retain the original radical caps and any least-repair
eligibility condition. Suppose \(V\) admits a divisor
\(D\ge X^w\), for a fixed \(w>0\), with
\(\operatorname{rad}(D)\le n^\mu\). The divisor may
vary from row to row. Suppose
\[
                             U=A y^k
\]
with integral positive \(A,y\), with \(k\ge K(L)\)
an integer chosen per row, and with \(K(L)\to\infty\).
The coefficient dictionary in the original fibre
has size at most \(n^{b+o_X(1)}\).

**APC1.** For every \(\eta>0\), at sufficiently small
fixed \(L\) the count is
\[
                        O(n^{b+\mu+\eta}X^\varepsilon)
                           \qquad(\varepsilon>0).
\]
The constants are uniform in \(T\), the divisor
and coefficient dictionaries, and the permitted
order choices. One can use a fixed lower bound
on the orders in place of \(K(L)\to\infty\) if
\(k w>1\) uniformly. An actual divisor dictionary
may be priced at its cardinality instead of by
the radical cap.

To prove the estimate, fix an integer \(J\) with
\(1/J<\eta/4\). In the range \(k\le J/L\) there
are finitely many orders at each fixed budget.
Fix \(D\) by recording its squarefree support
and its full exponents. This costs at most
\(n^\mu X^{o_X(1)}\), uniformly at polynomial height.
By primitivity the values \(A,y,T\) are units at
\(D\). The signed equation prescribes the residue
of \(y^k\) there. For each fixed \(k\) the number
of possible unit residues of \(y\) is
\(X^{o_X(1)}\). Indeed at each prime power there
are at most \(2k\) such roots, by the structure
of the group of units of a prime power. The product
of these bounds is height-subpower since \(D\le X\).
Since \(y\le X^{1/k}<D\), each residue determines
the root. Its coefficient then fixes \(U\), and
the prescribed \(T\) fixes \(V\).
For \(k>J/L\) and \(y>1\), enumerate
\(y\le n^{1/J}\) and \(k\le\log_2 X\).
When \(y=1\), the coefficient already determines
\(U\). This branch needs no enumeration of \(D\).
Choosing the inner allowances inside \(\eta\)
proves the assertion.
[Unit power roots in exact-layer moments](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#joint-moments-in-a-slot)

In a positive-exponent fixed-letter fibre at
normalized exponent \(u>0\), FL2 gives nonatomic
mass at least \(u\) on each of the two free roles
in limiting valuation profiles. Here is the
resulting lock choice on the role opposite the
grid. Take a fixed compact positive interval
\([\alpha,\beta]\) carrying positive nonatomic
mass. Exclude primes of individual normalized
size above \(\tau\), with \(\tau\) fixed and
small. There are boundedly many such primes at
each fixed cutoff. Their weak limiting mass is
supported on finitely many atoms, so the selected
band still contains positive mass on the small
primes. Choose continuity endpoints of the
original measure. Thus for any sufficiently
small fixed \(\mu>0\), one can take
\(\tau<\mu/2\) and greedily select powers in
the band with total radical mass in
\([\mu/2,\mu]\). Their whole product has height
at least \(X^{\alpha\mu/2}\). The choices are
made in this order before taking the height
limit. The same argument is uniform in late
boxes of the fixed profile by a contradicting
subsequence of actual rows.
[Positive nonatomic mass on both free roles (FL2)](ABC_FIXED_LETTER_RADICAL_CAPS.md#strict-positive-scale-inequalities-on-the-two-free-letters)
· [Large-prime cutoff and the nonatomic remainder (JRM5)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#retaining-the-nonatomic-mass-after-large-prime-fixing)

It follows that on an exponent-full positive
fibre an adaptive grid on either one free role,
with \(K(L)\to\infty\), needs a conditional
coefficient dictionary of exponent at least
\(u-o(1)\) when the complementary lock is chosen
at arbitrarily small cost. In particular an
adaptive order whose off-grid radical has
normalized mass less than \(u\) by a fixed
margin gives a power-small subset of that
fibre. One can use the canonical coefficient
\(\prod_p p^{v_p(U)\bmod k}\), since it is the
coordinatewise remainder of any valid integral
coefficient. The counting statement applies
after summing the prescribed complete targets
at their actual dictionary price, where the
uniform hypotheses and margins hold.

## A separating pool for the joint coefficient

Now use actual primitive \(abc\) rows in the
shell \(X/2<C\le X\) and at total radical
budget \(L\). Let \(U,V\) be the proposed
powered roles and \(W\) the remaining role.
Suppose
\[
              U=A x^k,\qquad V=B y^k,\qquad
       A=\prod_p p^{v_p(U)\bmod k},\quad
       B=\prod_p p^{v_p(V)\bmod k}.
\]
The common \(k\) may depend on the row, with
the same lower bound \(K(L)\to\infty\) as
above. For the pool catalogue, treat all
\(k>\log_2 X\) as a single state. The
canonical coefficients and roots are the
same throughout that range. Let
\(Q_k=A/B\) be the canonical
coefficient ratio. Its reduced numerator
and denominator are \(A,B\), by primitivity.
Suppose \(W\) has a divisor \(D\ge X^w\)
with radical at most \(n^\mu\), where \(w>0\)
is fixed. At a fixed \(k\) allow a prescribed
catalogue of at most \(n^{b_0+o_X(1)}\)
prime pools. In each selected pool state
\(S\), assume that the valuations on \(S\)
separate the actual canonical coefficient
ratios of that state. For some fixed
\(a\in[0,1)\) require
\[
                   \sum_{p\in S,\ p\le X}p^{-a}
                              \le g(X)=o(\log X)
\]
uniformly over the pool catalogue, the
orders and the permitted orientations.
A pool's prime product may have polynomial
height; the displayed summability is the
condition used in the estimate.

**APC2.** For any fixed \(0<t\le1\), any
fixed counting margin \(\eta>0\), and all
sufficiently small fixed \(L\), the count is
\[
 O\!\left(
 n^{\,b_0+\max\{(1-(1-a)t)/3,\ at+\mu\}+\eta}
                                      X^\varepsilon\right).
\]
The constants are uniform in the chosen
pool states under the displayed uniform
bound. If a common \(k\) is fixed below
\(K(L)\), the same estimate holds whenever
\(kw>2\) by a fixed margin.

Fix \(k,S\). Let \(R_S\) be the product of
the primes on which the reduced coefficient
ratio has nonzero valuation in \(S\).
If \(R_S\ge n^t\), these primes occur on
the original \(U,V\), so they already
account for at least \(t\) of the true
triple radical in base \(n\). Price them
at \(a\) times their natural cost using
the prescribed entire \(S\). By EB5 the
count is at most
\(n^{(1-(1-a)t)/3+\eta'}X^\varepsilon\).
This branch uses the original total
radical constraint including the third
letter. In the other branch, Rankin's
trick bounds the number of allowed
projected supports by
\[
       n^{at}\prod_{p\in S,\ p\le X}(1+p^{-a})
                             \le n^{at}X^{o_X(1)}.
\]
The signed valuation vectors on any fixed
support of a polynomial-height rational
have height-subpower count. The separation
assumption therefore gives at most
\(n^{at+o_X(1)}\) exact coefficient ratios,
and this is also the bound for the
coefficient pairs in this state.
[Effective original radical pricing (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Separating prime pools (VS2–VS4)](ABC_VALUATION_SPAN_COSETS.md#strict-saving-from-a-separating-prime-pool)

For each such coefficient pair and each
choice of \(D\), the original equation
prescribes one value of \((x/y)^k\bmod D\).
All four quantities \(A,B,x,y\) are units
there. For \(k\le J/L\), with \(J\) fixed
as in the preceding proof, there are
height-subpower many unit-root classes.
Two pairs in one slope class satisfy
\(x_1 y_2\equiv x_2 y_1\bmod D\).
Both cross products are at most
\(X^{2/k}<D\) at sufficiently small fixed
\(L\), or whenever the asserted fixed
order inequality holds. Their equality
fixes the reduced positive pair, since
\(\gcd(x,y)=1\). Thus the light branch
costs \(n^{at+\mu+o(1)}\) in this range.
For \(k>J/L\) the two roots cost
\(n^{2/J}\), without enumerating \(D\).
There are at most \(\log_2 X\) orders
with either root larger than one. In
the stipulated single larger-order
state the coefficients specify the
original row. Choose \(J\) large enough
for the counting allowance, sum the
actual pool catalogue and absorb the
order count.
This proves APC2 without a uniform
subspace-theorem constant in \(k\).

For \(0\le\mu<1/3\), the heavy branch decreases
in \(t\), and the light branch is nondecreasing.
An optimal split and the resulting exponent are
\[
 t_*=\frac{1-3\mu}{1+2a},\qquad
 E=b_0+\frac{a+(1-a)\mu}{1+2a}.
\]
At \(a=0\) the light branch is constant, and
this split is still optimal. For \(a=1/2\) the
exponent is \(1/4+b_0+\mu/4\). The available
saving below one third, before count allowances,
is
\[
        \frac{(1-a)(1-3\mu)}{3(1+2a)}-b_0 .
\]
The formula includes the actual third-divisor
price in the choice of the split. A fixed positive
value in the last display gives a strict saving.
In a balanced chart with individual caps
\(1/3+o(1)\), one may also use the light
bound \(2a/3+\mu+o(1)\) when every
coefficient projected support is at most
\(n^{2/3+o(1)}\). If \(W\) has positive
nonatomic mass, the preceding compact
band construction permits arbitrarily
small fixed lock price, with a positive
\(w\) chosen before the budget limit.

For example, in each prescribed coset of
coefficient ratios of multiplicative
rank \(O((\log X)^\kappa)\), with a
fixed \(\kappa>0\), choose a basis of
prime valuation coordinates on the
difference span. Its size is
\(O((\log X)^\kappa)\), and its valuations
separate the coset. The partition sum
is sublogarithmic for every fixed
\(a>\max(0,1-1/\kappa)\) with \(a<1\).
A varying family of cosets or pools
requires the explicit catalogue charge
\(b_0\). A separating pool chosen
independently for each row without a
bounded catalogue does not give this
count.
[Valuation coordinate bases in rational cosets (VS1)](ABC_VALUATION_SPAN_COSETS.md#counting-low-rank-cosets-by-separating-valuations)

## Price of a cover by small relative coefficient heights

Let \(n=X^L\) and consider genuine primitive shell rows with
the original total radical budget \(L\), retaining any separate
letter caps and repair eligibility. Fix an orientation of two
letters \(U,V\) and put \(W\) for the third. Suppose a row
admits an assigned representation
\[
            U/V=Q r^k,\qquad H(Q)\le X^{C_Q},\qquad
                            H(r)\le X^{C_r/k},
\]
with \(Q,r\) positive rationals, fixed \(C_Q,C_r>0\), and
\(k\) a positive integer. All prime divisors in the numerator
or denominator of \(Q,r\) lie in the support of \(UV\).
For instance \(r=x/y\) for two exact weighted powers of
common order \(k\), and \(Q=A/B\). Canonical coefficients
give \(C_r=1\). Signed coefficients from rounding each
valuation \(e/k\) to a nearest nonnegative integer also
fit, with \(C_r=2\). The rounded root degree \(q\)
vanishes if \(e<k/2\), and otherwise \(kq\le2e\).
Unequal literal orders can first be regrouped at any
common divisor, using the corresponding effective roots.

For each row suppose that \(W\) has a divisor \(D\) with
\[
 X^w\le D\le X,\qquad \operatorname{rad}(D)\le n^\mu,\qquad
         p^e\parallel D\ \Longrightarrow\ Le\ge\sigma
\]
for fixed \(w,\sigma>0\). At each order and orientation
take a prescribed catalogue of at most \(n^{b_0+o_X(1)}\)
prime pools \(S\) with
\(\sum_{p\in S,\ p\le X}p^{-a}=o(\log X)\) uniformly
for fixed \(0\le a<1\). As in APC2 the support of the
projection of \(Q\) to \(S\) is the squarefree product
\(R_S\) of primes with nonzero valuation there. At
orders \(k>C_r\log_2 X\) the root is one and the slope
is \(Q\) itself. Combine these orders in one state,
with the same uniform pool and cover requirements.

Fix \(0<t\le1\). On the light branch \(R_S<n^t\),
after fixing \(k,S,D\) and the exact signed projection
\((v_p(Q))_{p\in S}\), assume the actual ratios in that
class are covered by at most \(n^{b_1+o_X(1)}\) subsets
with pairwise relative height
\[
                 H(Q_1/Q_2)\le X^{c\min(1,kL)}.
                 \tag{APC4a}
\]
Use right side \(X^c\) in the combined root-one state.
The cover can be chosen conditionally in each such
class. In particular an absolute height bound for the
class representative is already supplied by membership.
A ball around an arbitrary rational center with half
the exponent in (APC4a) satisfies the diameter bound,
by \(H(z_1/z_2)\le H(z_1)H(z_2)\).

**APC4.** For every fixed \(\eta>0\) there are
\(c>0\), an integer \(k_0\), and \(L_0>0\) depending
on the fixed numerical data and allowance such that,
for \(0<L<L_0\) fixed and \(k\ge k_0\), the count under
the preceding conditions is
\[
 O\left(n^{\,b_0+\max\{(1-(1-a)t)/3,\ at+\mu+b_1\}+\eta}
                                   X^\varepsilon\right).
 \tag{APC4b}
\]
In particular the theorem allows every lower threshold
\(K(L)\to\infty\). Put \(u=\mu+b_1\). For \(u<1/3\),
if the light-branch cover is available at
\(t=(1-3u)/(1+2a)\), the optimized exponent is
\[
 b_0+\frac{a+(1-a)u}{1+2a}.
\]
When \(b_1=0\) this is exactly APC2's optimized
exponent, with the separation of exact ratios
replaced by a subpower cover at the stated
relative width. The additional lower-depth
condition on \(D\) is used in the proof.

The heavy branch uses the EB5 proof in APC2 verbatim.
It prices actual primes of \(UV\) at discount
\(1-a\) in the original three-letter budget. On
the light branch there are at most \(n^{at+o_X(1)}\)
projected supports and \(X^{o_X(1)}\) signed exponent
vectors per support. There are \(n^{\mu+o_X(1)}\)
possible \(D\)'s. Count after fixing these data and
one member of the prescribed cover. Choose one
actual coefficient ratio \(Q_0\) from that subset
as a reference. It can depend on these fixed
data, with no further catalogue charge. Both
\(Q_0\) and all \(Q,r\) in the subset are units
at \(D\). Put \(z=Q/Q_0\).
The primitive equation yields
\[
                 r^k z=(U/V)/Q_0\equiv\pm Q_0^{-1}\pmod D .
\]
The sign is fixed by the orientation of the addition.
For fixed \(Q_0\), each distinct output \(r^kz\)
determines \(U/V\) and hence the original coprime
pair and row. Outputs have polynomial height
uniformly in the class.

Let \(K=\lfloor\sigma/L\rfloor\), so \(D\) is \(K\)-full.
Choose a fixed output amplification order \(\ell\)
with \(2/(\ell\sigma)<\eta/4\). For \(k\le J/L\)
with \(J\) fixed, use the rational two-base congruence
count with exponents \(k,1\) and numerical bounds
\(\xi_r=C_r/k,\ \xi_z=c\min(1,kL)\). The amplified
height constants \(C_i=\max(1,\ell K\xi_i/w)\)
satisfy
\[
 \frac{C_1C_2}{K}\ \le\
 \frac{1}{K}+\frac{\ell C_r}{wk}+\frac{\ell c}{w}
                             +\frac{\ell^2\sigma C_r c}{w^2}.
\]
Choose \(k_0\) large, \(c\) small and then \(L_0\)
small so this lies below the absolute reciprocal
constant of that count. The output count per class
is \(O(n^{\eta/4}X^\varepsilon)\), using
\(\log_X D\le1\) and \(K L\ge\sigma/2\). Its
constant is uniform over \(Q_0\) and \(D\). The
orders used in this step form a finite list
at each fixed budget. [Rational distinct-output congruences and amplification](ABC_SIGNED_VALUATION_PHASES.md#a-rational-congruence-bound-for-distinct-exact-outputs)

For \(k>J/L\), enumerate rationals \(r\) of height
at most \(n^{C_r/J}\) at cost \(O(n^{2C_r/J})\).
The prescribed unit congruence then determines \(z\):
two such rationals congruent modulo \(D\) have a
cross determinant of absolute value at most \(2X^{2c}\),
which is smaller than \(D\) if \(2c<w\) and the height
is large. Use \(2C_r/J<\eta/4\). The total number
of noncollapsed orders is logarithmic in \(X\).
In the collapsed state \(r=1\) the same uniqueness
applies, with no order enumeration. Adding the
actual cover size, catalogue costs and allowances
proves (APC4b). No separate lower bound on the
radical of \(D\) or restriction on its primes is used.

### Charging a centre's shared high-degree remainders

There is a second count for such covers that uses
the original radical budget directly. Let the row
have raw radical mass at least \(\rho\) on primes
satisfying
\[
                   |v_p(Q)|\ \ge \theta \min(k,1/L)
\]
for fixed \(\theta>0\). In a subset satisfying
(APC4a), the total absolute valuation weight of
\(Q/Q_0\) is at most \(2c\min(1,kL)\) in base \(X\).
At any such prime outside \(\operatorname{supp}(Q_0)\)
the contribution is at least
\(\theta\min(k,1/L)\log_X p\). Hence the row shares
at least \((\rho-2c/\theta)_+\) of its original raw
radical with the center. The same proof in the
root-one state uses cutoff \(\theta/L\). This
argument can use any fixed \(c>0\), without the
numerical restriction of the local congruence count.

Suppose a catalogue of pairs \((D,Q_0)\) with center
height \(X^{O(1)}\) and count \(n^{b+o_X(1)}\) is
prescribed, and each assigned \(D\mid W\) has raw
radical at least \(d_0\). For this count only, one
may also take \(D=1,d_0=0\), so no modulus depth
or height is required. The reference can be an
actual ratio in the same divisor state, so its
support is disjoint from \(D\). In any case the
row's shared center primes occur on \(UV\), and
are disjoint from \(D\) by primitivity. The
prescribed polynomial-height pool
\(\operatorname{supp}(D)\cup\operatorname{supp}(Q_0)\)
can be priced at zero by EB5. Its captured raw
mass is at least \(d_0+(\rho-2c/\theta)_+\).
Summing the catalogue gives coefficient at most
\[
                   b+\frac{1-d_0-(\rho-2c/\theta)_+}{3}
\]
up to any fixed positive counting allowance at
small fixed budgets. Thus a wide relative cover
can still save where its centers capture more
than three times their raw catalogue price,
with the divisor capture included. This uses
the actual supports and all three original
radical caps. [Zero-price polynomial-height prime pools (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Full coefficient price on a primitive two-cap matching

Use the existing complete-label primitive
matching construction with \(L=1/m\),
\(J=2m-1\), \(d=1/3\), \(\beta=d/J\),
and degree indices
\(\mathcal E=\{2m+1,\ldots,4m-1\}\)
on each of the two free roles \(U,V\).
It assigns one prime from a separate pool
to each slot. The matching has
\(n^{d-o_X(1)}\) rows at fixed large
integer \(m\), and its selected slot
fibres satisfy the simultaneous size
bound of that construction. Both free
roles have tall values and radical caps
\(n^d\), with diffuse limiting degree
profiles as \(m\to\infty\). The third
role \(T=U+V\) is an actual complete
integer with radical above the cap.
[Simultaneous slot-fibre construction](ABC_OVERLAPPING_FULL_RECORD_COUNTEREXAMPLES.md#simultaneous-entropy-saturation)

**APC3.** On that matching, allow every
row to choose exact representations
\(U=A x^{k_1}\), \(V=B y^{k_2}\) with
arbitrary positive integral factors and
arbitrary orders \(k_1,k_2\ge2\).
The resulting ordered coefficient
dictionary needs \(n^{d-o_X(1)}\)
entries for the full matching, even for
unequal orders and noncanonical
coefficients. Inside any fixed slot
record of mass \(t<d\) carrying
\(n^{d-t-o_X(1)}\) rows, it needs
\(n^{d-t-o_X(1)}\) entries. Thus
the joint record of the fixed slots
and the adaptive coefficients cannot
have an exponent below \(d\) by a fixed
margin across the matching.

Indeed, fix the two coefficient values
and the orders. In any slot with
\(k_1\nmid e\), its prime divides \(A\),
since \(v_p(A)\equiv e\bmod k_1\).
It is the unique prime of \(A\) in that
slot's prime pool. The same holds for
\(B\) and the slots with \(k_2\nmid e\).
There are at most \(m-1\) divisible
indices in \(\mathcal E\) for each
order at least two. For order two this
follows from the odd endpoints; for
larger orders the interval bound is
\(\lceil(2m-1)/3\rceil\le m-1\)
for \(m\ge2\). Thus the two recorded
complements name together at least
\(2m\) slots of total mass
\(2m\beta>d\). The simultaneous
fibre bound leaves only \(O_m(\log X)\)
rows per assignment, independently of
the chosen orders. An order attached
to a root larger than one is at most
\(\log_2 X\). If the root is one, the
coefficient gives that whole letter
and no order record is needed. Therefore
each coefficient pair has only
polylogarithmically many preimages
under any adaptive order selection.
This remains true after any row
deletions or fixed slot records.

The same lower bound permits rational
coefficients and roots when the coefficient
heights are at most a fixed power \(X^H\).
In an off-grid slot its prime must occur in
the reduced numerator or denominator of
the coefficient. This leaves at most
\(O_H(\log X)\) candidates in each fixed
slot pool. The number of slots is fixed at
each \(L\), so this adds only a polylogarithmic
factor. For a root \(x\ne1\), its rational
height satisfies
\(H(x)^k=H(U/A)\le X^{1+H}\).
Thus its adaptive order still has only
\(O_H(\log X)\) choices.

The construction supplies a barrier to
recovering an affordable two-role
coefficient dictionary from primitive
addition, tall numerical heights and
the two original operand caps. It
leaves the small radical of the
remaining letter, and the existence of
an economical full lock there, as
additional arithmetic inputs.

### Prime degree indices and several adaptive grids

The same two-cap construction permits prime
degree indices and consequently several adaptive
grids per role. For each integer \(m\) take the
set \(\mathcal E_m\) of prime indices in
\([2m,4m]\), put \(J=|\mathcal E_m|\),
\(\beta=d/J\), and choose
\[
                   L=\frac{3}{J^{-1}\sum_{e\in\mathcal E_m}e}.
\]
Then \(L\sim1/m\), and
\(L\beta\sum_e e=1\) exactly. Choose the
disjoint slot prime pools and their values
as in the same two-cap matching construction.
The prime number theorem gives the same
limiting uniform measure on scaled degrees
\([2,4]\). In particular the high slots
\(Le\ge79/25\) on both roles have summed mass
tending to \(7/25\) and summed numerical
weight tending to \(1253/1250\). Thus their
fixed-target CRT count, the deletion of
low-radical targets, and the simultaneous
sampling of slot fibres apply with the
same strict margins.

Each order at least two divides at most one
of the chosen degree indices. If primes in
the moving part of a letter are assigned
to disjoint exact power roots, at most one
nonexceptional slot per root can be removed
from the coefficient record. For any fixed
total of fewer than \(J\) grid orders across
the two roles, the combined off-grid slots
therefore have mass strictly greater than
\(d\). Fixing the coefficients and orders
leaves only polylogarithmically many matching
rows, as in APC3. The number of orders per
role can grow with the small-budget limit
subject to this bound, since the slot lists
and the number of recorded orders remain
fixed during each height limit. The claim
covers the exact divisor-grid assignments
in FV1.
[Several adaptive divisor grids (FV1)](ABC_FINE_VALUATION_DIVISIBILITY.md#a-bounded-number-of-adaptive-valuation-grids)


## The extraction obligation with all three original caps

Put \(n=X^L\). Suppose that a sequence of fixed
budgets \(L_i\downarrow0\) has genuine original
primitive shell families of exponent
\(1/3-o_i(1)\) at arbitrarily large inner
heights. A sufficient extraction keeps at
least \(n^{1/3-\varepsilon_0}\) rows for fixed
\(\varepsilon_0\ge0\), with a count exponent
at most \(1/3-\Delta\) before arbitrarily
small allowances, for some fixed
\(\Delta>\varepsilon_0\). The outer budgets
and inner heights remain a subsequence of
the given ones. Use two genuine powered roles
and a divisor of fixed positive numerical
weight on the third role. Allow the power
order to tend to infinity with the outer
budget or use a qualifying fixed lower order.

By CPW5 a direct joint coefficient and
divisor catalogue suffices if its total
price is at most \(1/3-\Delta\). By APC2
a prescribed separating-pool catalogue
suffices if
\[
       b_0\le
       \frac{(1-a)(1-3\mu)}{3(1+2a)}-\Delta,
           \qquad 0\le a<1,\quad 0\le\mu<1/3 .
\]
Here \(\mu\) charges the true third divisor and
\(b_0\) the true number of pool states. If the
records split into finitely many regimes before
the inner limits, sum their actual prices or
retain one regime at the stated extracted
exponent. Choose the count allowances and
\(\varepsilon/L_i\) at each fixed budget
inside \(\Delta-\varepsilon_0\). Either
count then contradicts the extracted lower
bound at arbitrarily large heights. This
would eliminate all near-saturation sequences
to which the extraction applies.

APC3 shows why retaining the third original
cap matters. Two capped operands can
carry the entire coefficient dictionary
even after adaptive unequal orders and with
large heights and primitive addition.
Their third letters can have large radicals.
The open input is a numerical or arithmetic
restriction on this joint coefficient and
divisor data in actual three-cap families.
[Elementary common-power reconstruction (CPW5)](ABC_COPRIME_SHORT_POWER_WINDOWS.md#an-adaptive-common-power-sector-with-a-third-letter-lock)
