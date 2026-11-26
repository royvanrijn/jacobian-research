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
The constants are uniform in the chosen pool states
under the displayed uniform bound. More generally it
suffices to have \(D\ge X^{2/k}\) on each row, so a fixed
order with \(kw\ge2\) qualifies under the exact height
bounds. The small fixed budget is used in the heavy branch.

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

For the light branch apply CPW5 after fixing the
coefficient pair and the actual divisor. Its slope
condition is \(D\ge X^{2/k}\), which follows from
the stated hypotheses on \(k\) and \(w\).
The joint recording cost in that branch is
\(n^{at+\mu+o_X(1)}\). The order states number
\(O(\log X)\) under the root-one convention, and
CPW5 includes primes dividing the order and the large-order
branch. [Primitive slope reconstruction (CPW5)](ABC_COPRIME_SHORT_POWER_WINDOWS.md#an-adaptive-common-power-sector-with-a-third-letter-lock)

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

### Comparable unequal orders by rational approximation

The same count and prices apply to two independently
chosen literal orders \(k_1,k_2\) with
\(\kappa_-\le Lk_i\le\kappa_+\) for fixed positive
\(\kappa_-,\kappa_+\). Assume
\(U/V=Qx^{k_1}/y^{k_2}\),
\(H(x)\le X^{C/k_1}\), \(H(y)\le X^{C/k_2}\),
and the same support condition on \(Q,x,y\).
Take \(S\) at the actual pair of orders and require
the light cover for original \(Q\) with diameter
\(H(Q_1/Q_2)\le X^{c'}\) for a sufficiently small
fixed \(c'>0\). This includes coprime literal orders.

Fix a large integer \(M_0\). The pigeonhole
principle applied to the fractional parts of
\(j k_2/k_1\) for \(0\le j\le M_0\) gives
\(1\le a_1\le M_0\) and an integer \(a_2\) with
\[
             |a_1 k_2/k_1-a_2|\le 1/M_0 .
\]
For \(M_0>\kappa_+/\kappa_-\) we have \(a_2>0\).
Let \(g\) be a nearest integer to \(k_1/a_1\) and
\(E_i=k_i-a_i g\). Then
\[
 |E_1|\le a_1/2,\qquad |E_2|\le k_1/(a_1 M_0)+a_2/2,\qquad
 \frac{\kappa_-}{2M_0}\le gL\le \kappa_++1
\]
at small fixed budgets. For
\(r=x^{a_1}/y^{a_2}\) and
\(z=(Q/Q_0)x^{E_1}/y^{E_2}\), we have the exact
identity \(r^g z=(U/V)/Q_0\).
The bound \(a_i/k_i\le2/g\) for small \(L\) gives
\(\log_X H(r)\le4C/g\). The error powers contribute
at most
\[
 \frac{C(|E_1|/k_1+|E_2|/k_2)}{gL}
       \le \frac{2C}{M_0\kappa_-}+O_{\kappa_-,\kappa_+,C}(L M_0^2).
\]
First set the local width \(c\) from the congruence
argument with root height constant \(4C\). Choose
\(M_0\) large, then \(c'\) small using the positive
lower bound on \(gL\), and then \(L\) small, so the
sum of the two error powers and the original
relative-height charge \(c'\) is at most
\(c\min(1,gL)\). Now use the same distinct-output
count with order \(g\). This regrouping occurs
inside the relative-height completion step.
The prime projection and the heavy branch still
use the original \(Q\), and the number of
literal order pairs is finite at each fixed
budget before the height limit.

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

### Weighted overlap of two coefficient ratios

There is a sharper support charge when coefficient
valuations satisfy \(|v_p(Q)|\le B k\) for fixed \(B>0\).
For a projection on a prescribed pool \(S\), define
\[
 \mathcal W_S(Q)=\sum_{p\notin S}\frac{|v_p(Q)|}{k}\log_n p .
\]
Let \(Q_0\) be a reference in the same class and
suppose \(H(Q/Q_0)\le X^{c kL}\). Outside the support
intersection, the absolute valuation contributions
from \(Q\) and \(Q_0\) are disjoint and sum to at
most \(2c\) in the units of \(\mathcal W_S\).
On the shared support outside \(S\), the two
weights together cost at most \(2B\) per unit of
raw radical mass. Thus
\[
 \log_n\left(\prod_{\substack{p\notin S\\p\in\operatorname{supp}(Q)
                                  \cap\operatorname{supp}(Q_0)}}p\right)
     \ \ge\ \frac{(\mathcal W_S(Q)+\mathcal W_S(Q_0)-2c)_+}{2B}.
 \tag{APC5}
\]
In particular if both weights are at least \(A\)
then a fixed polynomial-height center captures
at least \((A-c)_+/B\) outside the projected primes.
That captured mass can be added to the projected
and actual third-divisor masses in the zero-price
EB5 count. This argument also applies with the
center record alone, independently of the
congruence height threshold in APC4.

### Fine-grid distances at balanced saturation

For a consequence at the balanced endpoint, use
a near-saturating sequence of genuine rows with
stabilized scaled valuation profiles, with the
inner height limit at each fixed budget as usual.
On two chosen roles let \(M\) be the combined
nonatomic radical mass. One can restrict the
prime slots of each row on those roles to
submeasures converging to the nonatomic parts
of the two profiles. For each fixed precision,
remove small neighborhoods of finitely many
atoms covering all but that precision of the
atomic mass, with continuity endpoints and with
nonatomic mass below that precision. Cut off the
low and high ends with the same negligible
nonatomic loss, again at continuity endpoints.
Uniform weak convergence then gives the
restricted profile at that precision. Let the
precision tend to zero sufficiently slowly
across the outer budgets. The degree-index
selectors are fixed at each retained budget
before using arbitrarily large inner heights.
In particular the selection never fixes a
prime value or charges a record of such values.

Choose common fine-grid orders \(k\in[K,2K)\)
with \(K\asymp L^{-\alpha_0}\), \(1/2<\alpha_0<1\).
For each row let \(Q\) be its centered signed
remainder ratio on the selected slots of the
two roles. This observable is supported on the
original letters and has
\(|v_p(Q)|\le k/2\). The two signed moment
calculations in MPC2 give
\[
 \mathcal W_{\varnothing}(Q)=M/4+o(1),\qquad
             \frac{\log_X H(Q)}{kL}=M/8+o(1)
\]
for all but a proportion tending to zero of
the orders, uniformly over the stabilized rows.
Uniformity follows by choosing violating rows
along a subsequence and applying the same
compact-band and resonance argument to their
convergent restricted profiles. Averaging the
bad-order proportions over the rows shows that
for all but a proportion tending to zero of the
orders these conclusions hold on all but a
proportion tending to zero of the rows.
[Centered remainder moments on nonatomic bands (MPC2)](ABC_MOVING_SIGNED_POWER_COEFFICIENTS.md#centred-remainder-benchmark-on-diffuse-spectra)

For the quantitative bound keep the full empirical
distribution of the normalized absolute remainders.
The MPC2 second-moment proof also applies to the
bounded-variation indicators
\(\mathbf 1_{\{\|e/k\|\le v\}}\) at finitely many
fixed \(v\in[0,1/2]\). Their circle means are
\(2v\). Choose a finite mesh in \(v\) before
the outer cutoff, then use monotonicity between
mesh points. Thus for most \(k\) and most rows
the distribution function of the selected raw
mass on normalized magnitudes at most \(v\)
is uniformly close to \(2Mv\). To obtain
uniform \(o(1)\) error choose the mesh to
refine sufficiently slowly with the budgets,
using the same stabilization and averaging.
Call rows with these estimates and the signed
height estimates good at the given order.

Fix \(0\le c<M/4\) and fix an actual good center
\(Q_0\). In either row delete the selected primes
that belong to both coefficient supports,
of shared raw radical mass \(I\). Each remaining
measure has mass at least \(M-I-o(1)\). The two
remaining signed coefficient contributions to
\(\log_X(\operatorname{num}(Q/Q_0)
          \operatorname{den}(Q/Q_0))\) are disjoint.
A submeasure of raw mass \(v\) in the uniform
law of density \(2M\) on \([0,1/2]\) has first
moment at least \(v^2/(4M)\), by taking its mass
from the smallest magnitudes. With a uniform
distribution-function error \(o(1)\), this
lower bound has error \(o(1)\) uniformly for
\(0\le v\le M\): integrate the lower tails
\((v-2Mt-o(1))_+\) over \(0\le t\le1/2\).
The two rows therefore contribute total
absolute weight at least
\((M-I)^2/(2M)-o(1)\) in units of \(kL\).
Since \(H(Q/Q_0)\le X^{c kL}\) allows total
absolute weight at most \(2c\), this proves
\[
       I\ge M-2\sqrt{cM}-o(1).
\]
In particular (APC5) already gives \(I\ge
2(M/4-c)-o(1)\) using only the signed means.
The distribution improves that estimate at
the smaller radii.

**APC6.** For any fixed counting allowance \(\eta>0\)
and sufficiently late fixed budgets, uniformly
in any such good reference \(Q_0\), the number
of good original rows with
\(H(Q/Q_0)\le X^{c kL}\) is at most
\[
 O\left(n^{(1-M+2\sqrt{cM})/3+\eta}X^\varepsilon\right),
                  \qquad 0\le c<M/4.
 \tag{APC6}
\]
Prescribe the center's polynomial-height support
in EB5 and use the shared raw mass just proved.
The original three caps remain in the count.
In particular for any fixed \(0<\delta<M/4\),
take \(c=M/4-\delta\) and choose the allowance
small enough to obtain \(O(n^{1/3-\delta/3}
X^\varepsilon)\). For a near-saturating family
the exceptional fraction among pairs of good
rows is thus power-small. The height product
bound in the other direction gives normalized
relative height \(M/4+o(1)\) for the remaining
pairs of good rows.

In fact the reduced numerator and denominator
each attain this scale on typical pairs. For
any prescribed reference support, EB5 also
makes overlap of positive fixed raw mass
power-small on a saturating family. Centered
rounding caps the magnitude of each prime
valuation by \(k/2\), so the canceled
numerator–denominator factor of \(Q/Q_0\)
has numerical height exponent at most \(I kL/2\).
Each unreduced side has exponent
\((M/4+o(1))kL\) by the signed moment calculation.
The same cancellation argument applies to
any fixed-length product of independently
chosen good coefficients, each raised to
a fixed signed integer power. A union bound
over the finitely many pairs makes the
total shared support negligible on typical
tuples. For fixed exponents \(a_i\), both
unreduced sides, and hence both reduced
sides, have normalized exponent
\((M/8)\sum_i|a_i|+o(1)\) in units of \(kL\).

Balance and the one-third diffusion theorem
give total nonatomic mass at least \(1/2\).
Choose the two most diffuse roles, which gives
\(M\ge1/3\). The conclusion therefore reaches
every balanced limiting profile, including
mixed atomic profiles after the stated
restriction. It requires no assumption on
complete dictionary exponent, multiplicity
or existence of full divisors in a cell.
It shows why the individual numerical
height bound from fine-grid centering
does not by itself provide a low-cost
relative cluster. The pair law concerns coefficients on the selected slots
after removing the atomic parts. [One-third balance and diffusion](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13)

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


### Numerical relative-height covers on the same matching

Keep \(L=1/m\), \(n=X^L\), \(d=1/3\), \(J=2m-1\) and
\(\beta=d/J\) from APC3. Write \(H\) for the usual reduced
rational height. There is an absolute \(c_*>0\) such that a
centre \(Q_0\in\mathbb Q^\times\) of height at most \(X^B\)
with \(B\) fixed covers only \(O_{m,B}((\log X)^2)\) rows
of that matching under the following rule. On a covered
row choose any integer \(k\ge16\) and positive rational
representations \(U=A x^k,\ V=B' y^k\). With \(Q=A/B'\),
require
\[
                 H(Q/Q_0)\le X^{c_*\min(k,m)L}.
\]
The choices and the centre may depend on the rows covered.
The bound is uniform in the centre and applies to any
subfamily. Thus a cover of the full matching costs
\(n^{d-o_X(1)}\) centres. Within a fixed slot record of
price \(t<d\) carrying \(n^{d-t-o_X(1)}\) rows, the
same rule costs \(n^{d-t-o_X(1)}\) centres.

To prove this, put \(r=x/y\). The identity \(r^k=(U/V)/Q\)
and \(H(Q)\le H(Q_0)H(Q/Q_0)\le X^{B+c_*}\) give
\(k\le (1+B+c_*)\log_2 X\) when \(r\ne1\), since then
\(H(r)\ge2\). Group those rows by their common chosen
\(k\). If \(r=1\) then \(Q=U/V\); keep all such rows
in one separate group.

For a fixed \(k\ge16\) put \(u=\min(k,m)\) and
\(\theta=1/64\). Call a slot expensive if its fixed
degree \(e\) satisfies \(\operatorname{dist}(e,k\mathbb Z)
\ge\theta u\). For \(k>8m\) all slots are expensive.
For \(16\le k\le8m\), the number of remaining degrees
on one role is at most
\[
          (2m/k+2)(2\theta u+1)
              \le m/4+2\le J/4
\]
for large \(m\). The first inequality counts the
multiples of \(k\) meeting the widened degree
interval, then counts the nearby integers per
multiple. Thus at most \(J/2\) slots across the
two roles are outside the expensive set.

Compare two rows in the same order group, with
ratios \(Q,Q'\). In a slot where their primes
differ, the two primes are absent from every
other slot in both rows. Their two valuations in
\(Q/Q'\) are congruent to \(e,-e\) modulo \(k\),
up to the orientation of the role, because
\(U/V=Q(x/y)^k\). For each expensive differing
slot both absolute valuations are at least
\(\theta u\). At sufficiently large \(X\) for fixed \(m\)
each of the primes has \(\log_n p\ge\beta/2\).
For any rational \(z\),
\(\log H(z)\ge \tfrac12\sum_p |v_p(z)|\log p\).
Consequently each expensive differing slot
contributes at least \(\beta\theta u/2\) to
this lower bound in base \(n\). On the other
hand \(H(Q/Q')\le H(Q/Q_0)H(Q'/Q_0)\le
X^{2c_*uL}\). Taking \(c_*=\theta d/16\), at most
\(4c_*/(\beta\theta)=J/4\) expensive slots can
differ. There are at most \(3J/4\) differing
slots in total, so the two rows agree on more
than \(J\) slots.

In the group \(r=1\), every differing slot has
the two valuation magnitudes equal to its
degree \(e\ge2m+1\). The same argument with
\(u=m\) and every slot expensive works because
each relative height to the centre is at most
\(X^{c_*mL}\), regardless of the chosen order.

Fix one row in any nonempty group. Every other
row agrees with it on some subset of more than
\(J\) slots. For each of the finitely many such
subsets at fixed \(m\), the matching's simultaneous
slot-fibre bound allows \(O_m(\log X)\) rows with
the prescribed primes, since their total cost
exceeds \(d=\beta J\). Summing over subsets and
groups proves the bound. Each source row in this
construction may also be required to have
\(\operatorname{rad}(T)>n^{d+\vartheta}\) for any fixed
\(0<\vartheta<4/75=d-7/25\). Indeed its candidate
bound for low-radical targets becomes
\(n^{d+\vartheta+t_m+o_X(1)}\). After sampling with
the stated probability the expected number of
these edges is \(n^{\vartheta+t_m+o_X(1)}\), below
the expected full sample exponent by a fixed
margin at large \(m\), since \(t_m\to7/25\).
Delete these edges together with the collisions.
Every simultaneous slot-fibre bound survives
deletion, and the sample still has exponent \(d\).
As the two original operand radicals each have
size \(n^{d+o_X(1)}\), this version of the
obstruction exceeds the original total radical
cap by a fixed power. It retains both operand caps.

## Conditional separation with a joint divisor pool

Keep the genuine primitive shell and total budget \(L\)
of APC2. Fix \(\beta,\mu_0\ge0\), and at each orientation
and order allow a prescribed catalogue of at most
\(n^{b_0+o_X(1)}\) states \((S,T)\). Require
\[
        \sum_{p\in S\cup T,\ p\le X}p^{-a}=o(\log X)
                     \quad\text{uniformly},\qquad 0\le a<1 .
\]
Each row is assigned a state and an actual divisor \(D\mid W\)
with \(D\ge X^{2/k}\) and
\(\operatorname{rad}_{T^c}(D)\le n^{\mu_0}\). Here and below
\(\operatorname{rad}_{T^c}\) retains the primes outside \(T\).
Fix the actual divisor and the signed coefficient valuations
on \(S\). Require at most \(n^{\beta+o_X(1)}\) possible
exact coefficient ratios in this conditional class; \(\beta=0\)
includes exact separation. Use canonical positive coefficients
and the single state for orders above \(\log_2 X\), as in APC2.

**APC7.** For any fixed \(0<t\le1\) and any \(\eta>0\), at
sufficiently small fixed \(L\) the count is
\[
 O\!\left(n^{\,b_0+\max\{(1-(1-a)t)/3,\ at+\mu_0+\beta\}
                   +\eta}X^\varepsilon\right).
\]
In particular if \(u=\mu_0+\beta<1/3\), the optimum is
\(b_0+[a+(1-a)u]/(1+2a)\).

For a fixed state measure the combined support used in the
records, namely the primes of \(Q_k\) on \(S\) and of \(D\)
on \(T\). These sets are disjoint by primitivity and both
occur on the original triple. If the product is at least
\(n^t\), price the union \(S\cup T\) by EB5 to obtain the
heavy count. Here the actual primes of \(D\) use the original
third letter's radical budget. Otherwise Rankin enumeration
gives at most \(n^{at+o_X(1)}\) combined supports. Assigning
each prime to one of the two records and choosing exact
signed valuations costs \(X^{o_X(1)}\), by the uniform
fixed-support estimate. The outside support of \(D\) then
costs \(n^{\mu_0}\), with complete exponents again at subpower
cost. The remaining ratio choices cost \(n^{\beta+o_X(1)}\)
conditionally. CPW5 completes the rows. Summing the order
states and balancing the two counts proves the formula.
[Effective radical pricing (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Slope reconstruction (CPW5)](ABC_COPRIME_SHORT_POWER_WINDOWS.md#an-adaptive-common-power-sector-with-a-third-letter-lock)

With \(T=\varnothing\) and \(\beta=0\), this already allows
the exact separation to depend on \(D\). If part of the
third divisor is supplied by the prescribed sparse alphabet,
the correct input price is its outside radical, together
with the conditional dictionary. For example, at \(a=1/2\)
the values \(b_0=3/50\), \(\mu_0=1/50\), \(\beta=0\) give
\(63/200<1/3\) before allowances. With full divisor price
\(\mu=9/50\), the APC2 formula at the same \(a,b_0\) would
give \(71/200\). These are comparisons of sufficient record
budgets.

## Conditional recovery from the root ratio

**APC8.** Use a centered signed decomposition
\(U/V=Q r^k\) at a fixed order in a genuine primitive shell
with total radical budget \(L\). In each prime position choose
a nearest nonnegative integer to \(v_p(U)/k\) or \(v_p(V)/k\)
on the respective role. Fix a tie convention. If \(D\mid W\)
and \(D\ge X^{kL}\), the actual root ratio \(r\) and \(D\)
determine the original row at that order. For every divisor
weight \(w>0\), at \(D\ge X^w\) this applies throughout
\(kL\le w\), including orders of size \(L^{-\alpha}\) for
fixed \(0<\alpha<1\) at sufficiently small fixed \(L\).

Indeed the reduced numerator and denominator of \(Q\) are
each at most \(\operatorname{rad}(UV)^{k/2}\le X^{kL/2}\).
The equation fixes \(Q\equiv\pm r^{-k}\pmod D\).
Two possible positive coefficient ratios therefore have
congruent cross products, both at most \(X^{kL}\). The
difference has magnitude strictly less than \(D\), so they
are equal. The exact complete odds \(Q r^k\) recover the
original coprime positive pair. All quantities are units
at \(D\), by primitivity.

Consequently the APC7 estimate also holds with the signed
root valuations on \(S\) and at most \(n^{\beta+o_X(1)}\)
conditional exact root ratios. Use \(D\ge X^{kL}\) in this
route. The heavy cost has the same formula because the
projected root primes occur on the original powered roles.
A direct joint root-ratio and divisor dictionary also works
at its actual price. Where both \(D\ge X^{kL}\) and
\(D\ge X^{4/k}\) hold, the joint centered coefficient and
divisor dictionary and the joint root-ratio and divisor
dictionary have the same exponent at each fixed order.
For the converse direction the centered root height is at
most \(X^{2/k}\) and the same unit slope argument as in CPW5
applies to the signed decomposition. These statements retain
all three caps and do not require any valuation limit.
[Centered decomposition](ABC_MOVING_SIGNED_POWER_COEFFICIENTS.md#a-count-with-unrecorded-coefficients)

A small relative cover for the root ratio can replace
conditional exact root enumeration. After fixing the actual
divisor and projection, choose an actual center \(r_0\)
in a class of diameter \(H(r/r_0)\le X^{c/k}\).
The positive rational output \(Q(r/r_0)^k\) has
height at most \(X^{kL/2+c}\) and a prescribed residue
\(\pm r_0^{-k}\bmod D\). Thus it is unique when
\(D\ge X^{kL+2c}\). This gives the same count with \(\beta\)
the conditional cover price, without a lower bound on the
individual scaled valuations of \(D\).

## A mixed fourth-power defect replaces exact power records

**APC9.** Fix which original role is \(W\). Define the
signed fourth-power-free part \(\kappa\) of \(U V W^2\),
taking a minus sign when \(W=U+V\). The support of this
record consists of the nonzero valuation residues modulo
four on \(U,V\) and the odd valuations on \(W\).
There are \(O_\varepsilon(X^\varepsilon)\) original primitive
rows per exact \(\kappa\), uniformly at height at most \(X\).
If \(4\mid k\) in a weighted power representation of \(U,V\),
only the remainders modulo four of their coefficients enter.

To prove the count, take signed \(P,Q\) of magnitudes \(U,V\)
with \(r=P-Q\) of magnitude \(W\), with \(P>0\). Then
\(P Q r^2=\kappa m^4\) for a positive integer \(m\).
Set \(x=Pr/m^2,\ y=Pr^2/m^3\). Direct substitution gives
\(y^2=x^3-\kappa x\). This curve is nonsingular and \((0,0)\)
has exact order two. The coefficients and the point have
polynomial height uniformly, so Dujella's Corollary 1.3
gives \(X^{o_X(1)}\) possible points. Since
\(x^3/y^2=P/r\) and \(\gcd(P,r)=1\), each point recovers
the signed primitive pair and the original row, up to the
bounded number of orientations.
[Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/html/2312.03655v1)

In the original total radical budget, suppose a prescribed
catalogue of \(n^{b_0+o_X(1)}\) pools \(S\) has
\(\sum_{p\in S,\ p\le X}p^{-a}=o(\log X)\) uniformly,
with \(0\le a<1\). If each assigned row has
\(\operatorname{rad}_{S^c}(|\kappa|)\le n^\mu\), then
\[
              O\!\left(n^{\,b_0+[a+(1-a)\mu]/(1+2a)+\eta}
                                    X^\varepsilon\right)
\]
holds for \(\mu<1/3\), any \(\eta>0\), and sufficiently
small fixed \(L\). Split at \(n^t\) on the projected defect
support. In the heavy branch these primes all divide the
original \(UVW\), so use EB5. In the light branch record
both parts of the defect support at \(n^{at+\mu+o_X(1)}\),
and choose residues in \(\{1,2,3\}\) at subpower cost.
Apply the fixed-\(\kappa\) count, then take the APC2 split.
A direct actual defect dictionary can also be used at its
own price. If the valuations of the defect on \(S\)
separate the actual records, the same bound holds with
\(\mu=0\) independently of the outside defect mass.
In the light branch the projection then determines the
whole record. This criterion uses the third letter's actual
parity and needs no numerical divisor.

For example, view the defect exponent vectors in a state
as elements of \((\mathbb Z/4\mathbb Z)^{\{p\le X\}}\)
and fix the sign. Let \(G\) be their difference subgroup.
If \(G\) has at most \(d\) generators, then
\(G[2]\) has dimension at most \(d\) over \(\mathbb F_2\).
Indeed its cardinality is \(|G/2G|\le 2^d\).
Choose at most \(d\) prime coordinates which separate
\(G[2]\), by ordinary linear algebra on the coordinates
\(0,2\). Projection on these coordinates is injective on
\(G\). Any order-four element in its kernel would have a
nonzero double in the kernel on \(G[2]\).
Thus a prescribed catalogue of such states with
\(d=O((\log X)^s)\) for fixed \(s>0\) supplies
separating pools for every fixed \(a>\max(0,1-1/s)\)
with \(a<1\). Use the actual state catalogue price and the
\(\mu=0\) formula. This only asks for finite-module rank
of the joint residue record.

For a necessary condition at saturation, write \(x_0\) for
the combined raw mass at odd valuations and \(y_T\) for the
mass at valuations \(2\bmod4\) on role \(T\), in base \(n\).
For any fixed tolerance \(\xi>0\) on an exponent-full sequence,
apart from a power-small fraction, simultaneously for each
choice of \(W\) one has
\(x_0+\sum_{T\ne W} y_T\ge1/3-\xi\).
Indeed the left side is exactly the radical price of the
mixed record and a cheaper support has a power-small
count by fixed-\(\kappa\) recovery. In particular
\(x_0+(2/3)\sum_T y_T\ge1/3-\xi\).
The inequalities allow broad residue laws and do not furnish
the needed defect extraction in arbitrary profiles.

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
budget or require the exact CPW5 condition
\(D\ge X^{2/k}\). In particular \(kw\ge2\)
qualifies under exact numerical bounds.

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
\(b_0\) the true number of pool states. For the
high-depth divisors of APC4, it suffices instead
to give a cheap conditional relative-height
cover on the light branch. Replace \(\mu\)
by \(\mu+b_1\) in the formula, where \(b_1\)
prices the cover after both the actual divisor
and the projected valuation vector are fixed.
Comparable unequal orders are also permitted under the
regrouping conditions of APC4. APC7–APC8 allow the same
price condition with \(u=\mu_0+\beta\) for conditional
exact ratio records, using the appropriate reconstruction
threshold on the actual divisor. In the centered root route
the ratio recorded is \(r\) and the threshold is \(D\ge
X^{kL}\). Alternatively, APC9 permits the same price condition
with \(u=\mu\) for the outside mixed-defect support and
requires no divisor. Each route must have fixed numerical
and counting margins along the retained outer budgets and
uniform subpower catalogue errors at arbitrarily large
inner heights.

This weaker extraction still needs genuine
three-cap arithmetic. The matching obstruction
forces the full small-radius cover price with
only two caps, and APC5–APC6 show what overlaps a
low-radius cover would impose on the diffuse
parts of actual near-saturating families. If the
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
