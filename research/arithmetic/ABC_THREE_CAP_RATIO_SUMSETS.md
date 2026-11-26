# Higher sums from all three radical caps

A rational observable whose reduced denominator contains the primes
of all three complete letters has near-diagonal additive moments
at the one-third counting scale. The estimate applies to arbitrary
fixed rational weights and right-hand sides of bounded polynomial
height. For ordinary complete-letter ratios, a separate recovery
argument gives near-maximal higher sumsets in one reciprocal
orientation, even when numerical weight escapes to high orders.

## A pointwise bound when the denominator records all three letters

Fix \(B,K_0,K_1,\eta>0\) and an integer \(j\ge1\).
For sufficiently small fixed \(L>0\) put \(n=X^L\) and let
\(\mathcal E\) be any family of positive primitive \(a+b=c\) with
\(X/2<c\le X\) and
\[
                         \operatorname{rad}(abc)\le n^B .
\]
Attach a rational \(f(x)\) to each row \(x=(a,b,c)\).
Assume \(H(f(x))\le X^{K_0}\) and that the prime support of its
reduced denominator is exactly the support of \(abc\).
For any fixed nonzero rational weights \(w_1,\ldots,w_j\)
of reduced height at most \(X^{K_1}\), and any rational \(t\),
one has
\[
 \#\{(x_1,\ldots,x_j)\in\mathcal E^j:
                   \textstyle\sum_{i=1}^j w_i f(x_i)=t\}
                  \ll_{B,K_0,K_1,j,\eta,L,\varepsilon}
                      X^{(Bj/6+\eta)L+\varepsilon}.
 \tag{RS1}
\]
Here \(H\) is the maximum of the absolute reduced numerator
and the positive reduced denominator. The coefficient sizes,
the attached values and the right-hand side can vary with \(X\).
The implied constant is uniform in them and in the family.
The small-budget threshold may depend on the parameters
fixed before the height limit.

Group equal rows in a solution tuple, add their weights, and
separate the \(k\) distinct active rows with nonzero resulting
weights from the \(s\) cancelled row groups. Thus \(k+2s\le j\).
There are only a bounded number of equality patterns at fixed
\(j\). Put into a prime pool \(P_0\) the supports of the numerators
and denominators of all the nonzero subset sums of the \(w_i\),
and the support of the reduced denominator of \(t\).
When a solution exists, the height of \(t\) is bounded by a
fixed power of \(X\), by the bounds on the original weights
and values. The pool is therefore supported on a polynomial-height
integer uniformly, including \(t=0\), whose denominator is one.

At a prime outside \(P_0\), a prime occurring on just one of the
\(k\) active rows would give a unique term of negative valuation
in the weighted sum. The valuation of the sum would then be
negative, in contradiction to the choice of the pool. Thus all
primes outside \(P_0\) occurring in the union of the active
rows' radicals occur on at least two active rows. Their union
has normalized radical mass at most \(kB/2\) in base \(n\).

Count the active rows in their order in the equality pattern.
Once the preceding rows are fixed, prescribe the pool consisting
of \(P_0\) and the prime supports of those rows. If the new
radical mass of the current row outside this pool is at most
\(v\) in base \(n\), the number of choices is
\[
                             n^{v/3+\zeta+o_X(1)}
 \tag{RS2}
\]
for arbitrarily small fixed \(\zeta>0\) and sufficiently small
fixed budgets. This is the effective radical count (EB4)
with zero prices on the prescribed pool. The pool has polynomial
height with exponent depending only on the parameters fixed before
the height limit, so it has \(O(\log X/\log\log X)\) primes
and satisfies the required sublogarithmic partition sum. One can
use a fixed padded relative budget \((B+\xi)L\) with \(\xi>0\)
small, since \(c>X/2\). The effective cap at this step is \(vL\)
in base \(X\), with an arbitrarily small fixed allowance in
(EB4). Bin the successive new masses on a sufficiently fine
finite mesh at fixed \(j,B,\eta\). Their sum is the union
mass, at most \(kB/2\), up to the rounding losses. Thus the
active rows cost at most \(n^{kB/6+o(1)}\) with arbitrarily
small fixed exponent allowances. Each cancelled group costs
at most the unrestricted shell bound \(n^{B/3+o(1)}\),
as also follows from (RS2) with the empty pool. The total
exponent is at most \((k+2s)B/6\) up to the allowances,
which proves (RS1).
[Effective radical counting (EB3–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)
· [Limiting small-budget coefficient](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

The same proof allows the denominator supports to agree with
\(abc\) outside any one prescribed polynomial-height prime pool.
Terms with distinct rows but identical attached values cause
no problem, since the valuation test is on the active rows
after merging equal rows. The count thus applies with or
without injectivity of the attachment.

### Radical surplus and separation by letter

For a fixed equality pattern with \(k\) active rows and \(s\)
cancelled groups, let \(b_i=\log_n\operatorname{rad}(a_ib_ic_i)\)
on the active rows. Write \(m_p\) for the number of active row
radicals containing \(p\), and define the surplus
\[
 S=\sum_{p\in P_0} m_p\log_n p+
   \sum_{\substack{p\notin P_0\\m_p\ge2}}(m_p-2)\log_n p .
\]
For a chosen order of the active rows, assign the mass of each
prime outside \(P_0\) to the first row containing it and to the role
\(T\in\{a,b,c\}\) where it occurs on that row. Let \(v_{i,T}\)
be these new masses, \(v_i=\sum_T v_{i,T}\), and put
\[
             A=\sum_{i=1}^k(3\max_T v_{i,T}-v_i).
\]
All quantities are nonnegative, and the prime multiplicity test
gives the exact identity \(2\sum_i v_i=\sum_i b_i-S\).

In the same uniform fixed-budget sense as (RS1), the subfamily
where \(\sum_i b_i\le b_*\), \(S\ge S_*\), and \(A\ge A_*\)
has count
\[
       \ll X^{\{sB/3+(b_*-S_*-A_*)/6+\eta\}L+\varepsilon}.
 \tag{RS5}
\]
Here \(b_*,S_*,A_*\) are fixed normalized bounds, and the order
of active rows can be chosen separately in the finitely many
permutation classes. To see the stronger rowwise exponent,
perform the anatomic decomposition and prime pricing of (EB4),
but apply (GD3) to singleton root dictionaries in place of
(EB2). The numerical caps and small initial coefficient are
those of the original anatomic decomposition. In particular,
the powers at primes of price zero retain their full numerical
weights. Round the actual prices of each root on sufficiently
fine meshes at fixed truncation. Their totals on the three
letters are bounded by \(Lv_{i,T}\) up to the rounding losses.
The contribution of (GD3) is half the sum of the two smallest
dictionary costs. This function is coordinatewise nondecreasing,
so the exponent is at most
\[
 \frac{v_i-\max_T v_{i,T}}2
       =\frac{v_i}{3}-\frac{3\max_T v_{i,T}-v_i}{6}
\]
in base \(n\), with an arbitrarily small fixed allowance.
Summing, rounding and counting cancelled groups as in (RS1)
give (RS5).
[Grouped dictionary two-letter bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)
· [Prime-pricing decomposition (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

For example, at fixed order \(j\), tuples with a fixed positive
deficit \(kB-\sum_i b_i\), a fixed positive surplus \(S\), or a
fixed positive asymmetry \(A\) in any row order are each
power smaller than \(n^{Bj/6}\) in the small-budget counting
scale. The residual active tuples must have almost all their
radical incidence on primes occurring exactly twice outside the
pool. Indeed on primes with \(m_p\ge3\) the total incidence
mass is at most three times their surplus, and pool incidence
itself contributes to \(S\). Their pair packets must also be
balanced in the three roles at either end, up to the tolerance
on \(A\). To check one pair of active rows \(i,l\), put
\(i\) second last and \(l\) last. The new mass at \(i\) consists
exactly of the primes occurring outside the pool on that pair
and no other active row. Small asymmetry balances their three
role masses on \(i\), and the reversed choice balances them on
\(l\). If a prime occurs on exactly two active rows outside
the pool, its exponents in the two reduced denominators also
agree, since otherwise the sum has a unique term of smallest
negative valuation. For \(F\) and \(G\) below these are the full
exponents on the attached complete letters.

### Improvement on a fixed complete-letter fibre

Fix a complete letter \(T\), and suppose the other two roles
\(U,V\) have raw normalized radical caps \(\ell_U,\ell_V\)
bounded independently of \(L\). Add the prime support of \(T\)
to \(P_0\). The active union outside the pool now has mass at
most \(k(\ell_U+\ell_V)/2\). For a new row with outside masses
\(u,v\) on the two free roles, (FL1) bounds its contribution
by \(n^{\min(u,v)/2}\), with fixed allowances as above. This
is at most \(n^{(u+v)/4}\). A cancelled group has at most
\(n^{(\ell_U+\ell_V)/4}\) choices up to allowances by (FL1)
with empty pools. Consequently the fixed-weight pointwise
count in this fibre is
\[
       \ll X^{\{j(\ell_U+\ell_V)/8+\eta\}L+\varepsilon}.
 \tag{RS6}
\]
The bound is uniform in \(T\) under a bounded raw cap for
its radical, which is automatic under the original cap.
In a balanced sector, the \(2m\)-term fibre energy is thus
at most \(n^{m/6+\eta}\) before arbitrarily small inner
allowances and with sufficiently fine balance tolerances.
[Uniform fixed-letter count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

### A single ratio and an observable separable by letters

For \(z=a/c\in(0,1)\), define
\[
                   F(z)=z-\frac1z+\frac1{1-z}.
\]
It is strictly increasing on \((0,1)\), since
\(F'(z)=1+z^{-2}+(1-z)^{-2}>0\).
Its reduced denominator has support and exponents exactly those
of \(abc\). In fact
\[
                 F(a/c)=\frac{a^2b-c^2b+c^2a}{abc}.
\]
At a prime dividing \(a,b,c\), respectively, the numerator is
congruent to \(-c^2b,c^2a,a^2b\), a unit by primitivity.
Its reduced height is \(O(X^3)\), and \(F(a/c)\) determines
\(a/c\) and hence the primitive row.

A second observable is
\[
                    G(a,b,c)=\frac1a+\frac1b-\frac1c.
\]
Its reduced denominator is also exactly \(abc\).
At each prime dividing a complete letter the corresponding term
has the unique negative valuation. This observable is
separable into the three reciprocal-letter functions. Its height
is \(O(X^3)\). In fact \(G=N/D\) in reduced terms with
\(N=c^2-ab\) and \(D=abc\). Since \(c>\sqrt N\) and
\(c^3-Nc=D\), strict increase of \(u^3-Nu\) on
\([\sqrt N,\infty)\) recovers \(c\). The sum \(a+b=c\) and
product \(ab=D/c\) recover the unordered pair \(\{a,b\}\).
Thus \(G\) is injective with the convention \(a\le b\).

For a finite set \(A\) of rationals write \(mA\) for its set of
\(m\)-fold sums. Applying (RS1) with \(j=2m\) and opposite
unit weights bounds the \(2m\)-term additive energy of the
row multiset of attached values by
\(n^{Bm/3+\eta}X^\varepsilon\). This is the diagonal exponent
when the row family has exponent \(B/3\).
For weights \(w_1,\ldots,w_m\), apply it with the two lists
\((w_1,\ldots,w_m)\) and \((-w_1,\ldots,-w_m)\).
Cauchy--Schwarz gives, for \(A=f(\mathcal E)\),
\[
                  |w_1 A+\cdots+w_m A|
              \gg |\mathcal E|^{2m}
                        n^{-Bm/3-\eta}X^{-\varepsilon}.
 \tag{RS3}
\]
The counts use original rows, so this assertion also holds
if \(A\) is an image with repetitions discarded. The same attachment is used in every slot of (RS3).

## Higher sums of ordinary ratios

For an ordered pair of distinct complete roles \(U,V\), put
\[
                R=\{U(x)/V(x):x\in\mathcal E\},
                 \qquad R^{-1}=\{1/z:z\in R\}.
\]
Every such ratio determines the primitive row. For any fixed
\(m\ge1\), \(K>0\) and \(\eta>0\), there are
\(\kappa,L_0>0\) such that the following holds at every fixed
\(0<L<L_0\). If \(|\mathcal E|\ge n^{B/3-\kappa}\), then one
of \(R,R^{-1}\), say \(R_*\), satisfies
\[
                   |w_1 R_*+\cdots+w_m R_*|
                              \ge n^{Bm/3-\eta}
 \tag{RS4}
\]
for every list of fixed nonzero rational weights of height at
most \(X^K\), for sufficiently large \(X\) at the fixed budget.
The orientation can be chosen independently of the list
of weights. The power statement is uniform in the weights.
This proof uses ordinary ratios rather than the
denominator-support hypothesis of (RS1).

Choose a large fixed \(D\). For a complete letter strip its
full prime powers with scaled valuation \(Lv>D\).
The stripped factor has radical at most \(n^{1/D}\)
and height at most \(X\), so it has at most
\(n^{1/D+o_X(1)}\) possible values. Write the two complete
values as \(A P,B_0 Q\) before deciding which role is the
numerator. On at least half the rows one common orientation
makes \(P\le Q\). Choose that orientation now, before the
coefficient pool and all subsequent tuple selections. It
corresponds to one of the two original ratios
\(\frac{A P}{B_0 Q}\).

Fix the weights and let \(S_0\) be the prime support of the
product of their numerators and denominators. For a fixed
\(\rho>0\), discard rows having normalized radical mass
at least \(\rho\) on \(S_0\). By (EB5) this is a
power-smaller set of rows if \(\kappa\) is sufficiently
small. Similarly, all but a power-smaller fraction of
distinct ordered \(m\)-tuples from the selected family have
total shared radical at most \(n^\rho\). For each pair of
slots condition on the first row and use its complete
prime support as the prescribed pool for the second, at
overlap threshold \(\rho/m^2\). The union of the exceptional
pair classes covers every tuple failing the desired bound.
Tuples with repeated rows have a power-smaller count as well.
All these estimates apply before fixing the stripped factors.
[Shared-support discounts and order of conditioning](ABC_SECANT_ENDPOINT_TRANSFER.md#small-intersections-and-output-multiplicities)

Now fix the \(2m\) stripped factors \(A_i,B_i\) in the tuple,
losing at most \(n^{2m/D+o_X(1)}\). Thus we retain at least
\(\lvert\mathcal E\rvert^m n^{-2m/D-o_X(1)}\) tuples up to
a fixed constant. Consider one exact value
\[
                         t=\sum_i w_i\frac{A_i P_i}{B_i Q_i}
\]
on this retained class. Record the squarefree common support
\(J\) consisting of primes shared between at least two entire
original rows, with product at most \(n^\rho\). There are
at most \(n^\rho\) such records. Remove from the reduced
denominator of \(t\) the primes of \(S_0,J\), and of all
the fixed factors \(A_i,B_i\). The remaining prime powers
are exactly the private parts \(Q_i^*\) of the \(Q_i\),
namely their full parts outside \(S_0\) and \(J\).
Indeed a prime in such a part divides no other row,
is absent from the coefficients, and neither divides
the other letter of the original pair on its own row nor
the stripped tail of the same letter. Its exact negative
valuation therefore survives the sum.

One can allocate the private prime powers of the denominator
among the slots in \(X^{o_X(1)}\) ways, since the denominator
has polynomial height. The remaining parts \(Q_i/Q_i^*\)
are supported on \(S_0\cup J\); their complete values can be
recorded at subpower cost by fixed-support enumeration.
Each has normalized radical at most \(2\rho\), by the first
discard and the bound on \(J\). Its numerical height is
therefore at most \(X^{2D\rho}\), since retained scaled
valuations are at most \(D\). At this point all \(Q_i\) and
all original denominators \(B_iQ_i\) are fixed.

The exact sum also fixes \(P_i\) modulo \(Q_i^*\).
Multiply it by \(Q_i^*\); every other term then vanishes
modulo that integer, with denominators units at its
primes, and the coefficient of \(P_i\) in the \(i\)-th
term is a unit. This argument includes the empty private
factor by taking its congruence to be void. Since
\(0<P_i\le Q_i\), the possible \(P_i\) lie in one progression
whose integer index ranges over an interval of length at
most \(1+X^{2D\rho}\). The same index presents the original
numerator \(A_iP_i\) as an affine function with arbitrary
real coefficients, and the original denominator is fixed.
The rank-one small-feature estimate on a complete-letter
fibre bounds the number of such original rows by \(n^{\tau}\),
up to inner allowances, for any prescribed fixed \(\tau>0\),
as soon as \(2D\rho\) is sufficiently small in terms of that
estimate. Apply this separately in each slot. For one exact
sum the retained preimage count is therefore at most
\(n^{\rho+m\tau+o_X(1)}\).
[Small-feature saving on a single fixed-letter fibre](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)
· [Uniform enumeration on fixed prime support](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

Choose \(D\) large enough to make \(2m/D\) small compared
with \(\eta\). Choose \(\tau\) small compared with \(\eta/m\),
then choose \(\rho\) small compared with \(\eta\) and with
the rank-one numerical threshold divided by \(D\).
Finally take \(\kappa\) small compared with the overlap gaps
and \(\eta/m\), and shrink the budget threshold. The loss
\(m\kappa+2m/D+\rho+m\tau\) and the inner allowances
fit inside \(\eta\). This proves (RS4). In particular
the ambient heights of the stripped tails never enter the
short-feature threshold, and no control on the escaped moment
is needed.

## Conditional compression and minimal repair

At total budget \(B=1\), any bound
\[
                         |m f(\mathcal E)|
                               \le n^{m/3-\sigma}
\]
for a denominator-recording attachment as above gives,
by (RS3), a count exponent at most
\(1/3-\sigma/(2m)\) with arbitrarily small fixed allowances.
For ordinary ratios, a bound of this shape for both orientations
excludes a near-full family by (RS4). These are sufficient
compression hypotheses. No such compression follows from
the cap on each original radical alone by direct rational
denominator enumeration, which allows much larger sumsets.

For \(G(a,b,c)\), a trade matching the complete multisets of
all three letters produces an exact additive equality in its
image. For each left ordered \(m\)-tuple, there are at most
\((m!)^2\) ways to permute its \(a\)- and \(b\)-lists on the
right. These lists already determine the ordered \(c\)-list.
Thus fixed-size trades supply at most
\((m!)^2|\mathcal E|^m\) energy identities, even in a dense
additive graph on all three complete dictionaries. A power
compression from such an energy argument requires additional
coincidences between different multisets.

The least target repair takes \((a,b,c)\) to
\((Z_M(c)a,Z_M(c)b,H_M(c))\) with \(H_M(c)=Z_M(c)c\).
The forward ratio \(a/c\) and its reciprocal \(c/a\) are
unchanged, as is \(F(a/c)\). Apply (RS4) to the two
orientations of the repaired ratio, or (RS1)--(RS3) to its
single image under \(F\). The reciprocal-letter expression
on the repaired triple is \(G/Z_M(c)\). Its reduced denominator
is exactly \(Z_M(c)abc\), because \(Z_M(c)\) is supported on
\(c\) and the original numerator is coprime to \(abc\). This
has exactly the original radical support. Its height is
polynomial for bounded \(ML\), since
\(Z_M(c)\le \operatorname{rad}(c)^M\); (RS1)--(RS3) and
(RS5)--(RS6) apply directly to this repaired expression.
In the central normalization
\(Z=Y^{1/M}=n^{1/3+o(1)}\), a critical primitive-quotient
family has \(Z^{1-o(1)}\) rows. Both its \(F\)-image and its repaired reciprocal-letter
image have \(m\)-fold sumset \(Z^{m-o(1)}\) and \(2m\)-term
energy \(Z^{m+o(1)}\), with the moment counted on original rows.
At least one ordinary ratio orientation has \(m\)-fold
sumset \(Z^{m-o(1)}\) as well. These are expansion
consequences of criticality. The missing implication is an
arithmetic reason for compression in every candidate endpoint,
or a count that uses expansion and the third radical constraint
to exclude it.
[Least target repair and its primitive quotient](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

### A power saving for nontrivial complete-letter trades

Fix \(B>0\), \(m\ge 2\) and \(\eta>0\). At all sufficiently
small fixed \(L\), a family \(\mathcal E\) in the shell of (RS1)
has
\[
 \#\{\text{nontrivial ordered }m\text{-against-}m
       \text{ complete-letter trades in }\mathcal E\}
 \ \ll_{B,m,\eta,L,\varepsilon}\
       X^{\{Bm/3-B/9+\eta\}L+\varepsilon}.
 \tag{RS-trades}
\]
A trade matches the multiset of complete values separately
in every role on its two sides. It is nontrivial if the two
row multisets differ. The same estimate holds for the
complete repaired letters with bounded \(ML\), always
counting trades of minimally repaired primitive quotients.

For ordinary trades apply (RS5) to the signed equality of
the \(G\)-values. After merging equal rows, let \(k\) and
\(s\) be the number of active rows and cancelled groups.
A nontrivial trade leaves \(k\ge2\), with \(k+2s\le2m\).
For each active row and each role, its complete value must
occur on another active row in that role, since the signed
multiplicity of that value is zero. Write \(D=kB-\sum_i b_i\)
for the total raw-budget deficit of the active rows, and
retain the definitions of \(S\) and \(A\) from (RS5).

Let \(q_{i,l}\) denote the mass of primes outside the fixed
pool occurring on exactly the active pair \(\{i,l\}\).
If such a prime occurs in role \(W\) of row \(i\), then the
other active row with the same complete \(W\)-value must be
\(l\), since that value contains the prime. If a pair had
such primes from two roles at row \(i\), then two complete
values would agree, making the rows equal. Thus each nonzero
packet \(q_{i,l}\) lies in a single role of row \(i\).
Also row \(i\) has at most one such neighbor in each role.
If its value in that role occurred on three distinct active
rows, every prime of the value would have multiplicity at
least three. The graph of nonempty packets therefore has
at most \(3k/2\) edges.

Put \(P=\sum_{i<l}q_{i,l}\), and let \(A_{\max}\) be the
largest \(A\) over the orders of the active rows. In the order
putting \(i\) penultimate and \(l\) last, its entire new mass
is \(q_{i,l}\) and lies in a single role. Hence
\(A_{\max}\ge2q_{i,l}\), so \(2P\le (3k/2)A_{\max}\).
The total incidence of pool primes and primes of multiplicity
at least three is at most \(3S\). Consequently
\[
 kB\le D+3S+2P
          \le D+3S+(3k/2)A_{\max}
          \le (3k/2)(D+S+A_{\max}).
\]
In particular \(D+S+A_{\max}\ge 2B/3\).
Choose a maximizing order from the finitely many active
orders, bin the three nonnegative scores on a fine fixed
mesh, and apply (RS5). The exponent is bounded by
\((k+2s)B/6-B/9\) up to arbitrarily small allowances,
as required.

For minimally repaired trades use \(G/Z_M(c)\) in (RS5)
and balance complete repaired values in their designated roles.
Every original prime in a role divides the complete repaired
value of that role, and the repair introduces no prime outside
the original row's radical. For a prime on exactly two active
rows, its designated complete repaired value can thus repeat
only at the other row of the pair. Agreement of two complete
repaired values implies agreement of the third and of the
primitive quotients. Three occurrences of a designated
complete repaired value again force multiplicity at least
three of every original prime in that value. These are the
only properties of the complete values used in the packet
argument.

## Relation to overlap bounds and sum-product estimates

The general full-denominator bound of Bloom and Maynard,
Theorem 2, controls the \(2m\)-term energy of a set of reduced
fractions in \([0,1]\) with denominators at most \(Q\) and at most
\(d\) fractions for each denominator by
\((\log Q)^{Cm}(Qd)^m\). For fixed \(m\), their discussion
following Theorem 2 also points to a stronger
cardinality-sensitive estimate when denominators are distinct.
The proof of Lemma 5 of the cited STOC paper by Bourgain,
Dilworth, Ford, Konyagin and Kutzarova already bounds the
denominator condition \(q_i\mid\prod_{\ell\ne i}q_\ell\);
it gives \(X^\varepsilon\prod_i|\mathcal D_i|^{1/2}\) for
polynomial-height denominator dictionaries \(\mathcal D_i\).
General reduced numerators in a rational relation imply the
same denominator condition. A fixed inhomogeneous term
contributes a singleton dictionary.
[Bloom and Maynard, Theorem 2 and subsequent discussion](https://arxiv.org/pdf/2011.13266v2)
· [Bourgain, Dilworth, Ford, Konyagin and Kutzarova, STOC Lemma 5](https://denka.web.illinois.edu/rip_stoc.pdf)

In the setting of (RS1), a prescribed reduced denominator of
\(f(x)\) determines the original row's radical. Even with the
stated fixed-pool extension it determines the support of the
row outside that pool, so all candidate complete letters have
polynomial-height prescribed prime support. Fixed-support
enumeration gives subpower row multiplicity. For a reduced
weight \(w=\alpha/\beta\), a prescribed reduced denominator
\(q'\) of \(w f(x)\) implies that the reduced denominator of
\(f(x)\) divides \(|\alpha|q'\), since
\(f(x)=\beta\,w f(x)/\alpha\). This introduces only subpower
divisor multiplicity. The STOC denominator result therefore
gives the unstratified bound \(|\mathcal E|^{j/2}X^\varepsilon\)
directly, including odd and weighted pointwise counts, before
any use of the one-third bound. The arithmetic refinement here
is (RS5), which gives explicit savings for incidence surplus
and for separation failure among the three role costs. The
ordinary-ratio result (RS4) also allows large denominator
fibres, where the distinct-denominator argument need not apply.

Bloom and Kuperberg's Theorem 3 uses another set of parameters.
For odd \(k\), it bounds sums of \(k\) reduced fractions that
land in the integers, with numerator bound \(H\le Q\) and
denominator bound \(Q\), by
\(H^{(k+1)/2}Q^{(k-1)/2}(\log Q)^{O_k(1)}\). That
bound addresses numeric numerator and denominator sizes.
The stratification in (RS5) instead uses radical overlap and
individual role costs.
[Bloom and Kuperberg, Odd moments and adding fractions, Theorem 3](https://doi.org/10.1112/plms.70068)

The pole argument uses the full prime union across a fixed
linear pattern. Its count comes from charging each prime just
once, when it is new relative to the preceding rows. The
original radical caps enter through the bound on total
incidence mass and through the effective one-third count.
Weighted GCD sums measure a different quantity. For instance,
Bloom and Walker's Theorem 2 bounds the pairwise kernel
\((a,b)^{2\alpha}/(ab)^\alpha\) under an assumption on
multiplicative energy. Such a bound alone does not impose
the absence of a prime of multiplicity one in a rational
sum. The condition of repeated denominator primes follows
here from the exact equation.
[Bloom and Walker, GCD sums and sum-product estimates, Theorem 2](https://arxiv.org/html/1806.07849)

The few-prime-factor assumption of the integer sum-product
literature is available after a further regularization at
the balanced endpoint. For a fixed large \(J\), price
primes below \((\log X)^J\) at \(a\log_X p\) with
\(1-1/J<a<1\). Rows placing a fixed positive normalized
mass there are power smaller by (EB5). Keeping mass at most
\(\rho\log n\) on these primes bounds their number by
\((\rho+o_X(1))\log n/\log\log X\), by the primorial bound.
The larger primes contribute at most
\(B\log n/(J\log\log X)\). Letting \(\rho\) decrease and
\(J\) increase after their fixed-budget discards gives
\(\omega(abc)=o(\log n/\log\log X)\) on an exponent-preserving
sequence. The sum-product alternatives of Agrawal, Bloom and
Petridis still do not yield (RS1) or (RS4): at polynomial
height each fixed product of complete integers already has
divisor-many representations, and thus the large product
branch is automatic on an exponent scale. Recent product
expansion for translates has the same limitation on these
integer dictionaries.
[Agrawal, Bloom and Petridis, Theorems 8 and 9](https://arxiv.org/html/2512.04931v2)
· [Fang, He and Zhu, Products of translates, Theorems 2 and 3](https://arxiv.org/html/2611.30174v1)

## Full congruences from the last occurrence of a denominator prime

Use the projective coordinate \(F(a/c)\) above. Fix an active order
in a weighted rational equation and let \(y=(a,b,c)\) be the next row.
The earlier rows have been exposed. Put
\[
 S_{\rm rem}=t-\sum_{\text{earlier }i}w_i F(a_i/c_i).
\]
Here the \(w_i\) are the resulting nonzero weights after merging
equal rows. For the current row, mark primes outside \(P_0\) which
divide its denominator and are absent on every later row.
The later terms are integral at such a prime \(p\).
If \(p^e\parallel W\) on the current role \(W\), then
\(e=-v_p(S_{\rm rem})\). The local expression of the projective
coordinate is
\[
 F(a/c)\equiv \epsilon_W\frac{\ell_W}{W}\pmod{\mathbb Z_p},
 \qquad
 (\epsilon_a,\epsilon_b,\epsilon_c)=(-1,1,1),\qquad
 (\ell_a,\ell_b,\ell_c)=(c,c,a).
\]
In particular, with \(w\) the current weight,
\[
              w\epsilon_W p^e\ell_W-(p^e S_{\rm rem})W
                       \equiv0\pmod {p^{2e}}.                  \tag{RS7}
\]
Both coefficients are \(p\)-integral and the coefficient of \(W\)
is a unit. The factor \(p^{2e}\) comes from multiplying the
integral difference by \(p^e W\). The direction of the
congruence is fixed by the earlier rows. One can also use a
prime whose later weighted suffix has valuation at least
\(-r_p\) for \(0\le r_p<e\), with modulus \(p^{2e-r_p}\).
The equalities of the current and remainder pole orders still hold.

The marked prime set and its allocation to roles have only
\(X^{o_X(1)}\) possibilities, since they come from the reduced
denominator of the known \(S_{\rm rem}\), of polynomial height.
In the variant with \(r_p\), its choices over all primes are
at most the divisor count of that denominator. After recording
the primes and their roles, the inequality or divisibility
that selected them may be ignored for an upper bound.
The map from \((a,b)\) to each \((W,\ell_W)\) is unimodular.
Consequently the congruence forces the determinant of any two
current candidate rows to vanish modulo \(p^{2e}\), in any
pair of complete roles. The resulting modulus can be
multiplied by full divisors already fixed on the current row
at primes disjoint from the marked ones. At such divisors one
coordinate is zero, which gives the ordinary unsquared
determinant divisibility.

In a common dyadic height box, write \(h_W\) for the upper
numerical exponent of role \(W\) and \(H_*=\min_{U\ne V}(h_U+h_V)\).
If the product of all these moduli exceeds \(2X^{H_*}\), the
primitive current row is unique for the recorded data. In fact
the determinant on the least-height pair is zero, and
primitivity fixes the scale of the two proportional triples.
This step uses actual reduced pole orders and the homogeneous
linear principal parts of \(F\). General polynomial-height
rational coordinates need not have these principal parts.
The same coordinate \(F\) applies to the minimally repaired
quotient since \(a/c=(Z_M(c)a)/(Z_M(c)c)\).

## Near-maximal relations pair the repeated prime supports

Fix \(k\ge2\), \(B>0\) and the polynomial-height bound on weights
from (RS1). Choose \(0<\rho<B/(4(k-1))\). For sufficiently
small fixed \(L>0\), the count of solutions with \(k\) distinct
active rows to a fixed weighted equation in \(F(a/c)\)
has a power saving from \(n^{kB/6}\), except possibly on
tuples having the following property. Join two active rows
if the normalized radical mass of their exclusive common
primes outside \(P_0\) is at least \(\rho\). Exclusive here
means occurring on exactly those two active rows. Then the
resulting graph is a perfect matching. More precisely, there
exist \(\delta,L_0>0\), depending only on the fixed parameters
and \(\rho\), such that tuples outside this matching condition
are
\[
                         \ll X^{(kB/6-\delta)L+\varepsilon}
                   \qquad (0<L<L_0).                          \tag{RS8}
\]
The estimate is uniform in the right-hand side and the weights.
In particular it holds for the entire distinct-row solution
set if \(k\) is odd. It also holds when there is just one active
row, by injectivity. When estimating a tuple with cancelled
groups, sum over their choices at the usual cost \(B/3\)
per group. The saving then applies to each pattern with an
odd number of active rows.

Here is the proof. Let \(D_{\rm raw}=kB-\sum_i b_i\), and let
\(S,A\) have the meanings in (RS5). Fix \(R>36/\rho\) and fix
\(s_0>0\) with
\[
 s_0<\min\{\rho/2,\ B/16,\ 1/(48R)\}.
\]
By (RS5), discard all tuples with \(D_{\rm raw}\), \(S\), or
the asymmetry \(A\) in any active row order exceeding \(s_0\).
All these discarded sets have a uniform power saving, since
there are only a fixed number of orders. In what remains,
each row carries at least \(B-4s_0\) of exclusive pair mass.
In fact its total incidence is at least \(B-D_{\rm raw}\) and
the total pool incidence and higher-multiplicity incidence
over all rows is bounded by \(3S\). Thus each row has some
neighbor in the stated graph.

Suppose there is a row \(i\) with two different neighbors.
For its primes outside \(P_0\) and with \(Lv_p\le R\),
compare the numerical weights of the two exclusive pair
packets. Let \(j\) be the neighbor with the smaller such
weight \(q\). Order \(i\) penultimate and \(j\) last. The
new raw weight \(v_i\) at step \(i\) is exactly the exclusive
pair mass with \(j\) outside \(P_0\), so \(v_i\ge\rho\).
The exponent charged for this row by (GD3) in (RS5) is
\[
                  v_i/3-(3\max_W v_{i,W}-v_i)/6
                               \ \ge\ \rho/4 .
\]
We will replace it by at most \(3/R\).

Fix the full prime powers on row \(i\) in \(P_0\) and those
with \(Lv_p>R\), enumerating their three corresponding
products. The former have subpower cost by supported-part
enumeration. All latter tails together cost at most
\(3/R+o_X(1)\), since the sum of the three numerical
weights is at most three. Let \(p_0\) be the numerical
weight of the union of the fixed powers, each counted once.
Let \(H_{\rm tot}=\sum_W\log_X W\) on the current row.
The two exclusive compact packets used to select \(j\)
are disjoint and avoid the fixed powers, so
\(q\le(H_{\rm tot}-p_0)/2\). The other unfixed primes
of row \(i\) which also occur on row \(j\) have multiplicity
at least three among the active rows. They have normalized
raw weight at most \(3S\) and numerical weight at most
\(3RS\). All remaining unfixed prime powers on \(i\) are
last occurrences. Thus the determinant modulus from (RS7)
and the ordinary fixed powers has numerical weight at least
\[
 p_0+2(H_{\rm tot}-p_0-q-3RS)
                    \ \ge\ H_{\rm tot}-6RS.
\]
This exceeds \(H_*\) by a fixed margin. Indeed, in the current
dyadic box \(H_{\rm tot}=H_*+\max_W h_W+o_X(1)\) and
\(\max_W h_W=1+o_X(1)\) because \(X/2<c\le X\).
Our choice of \(s_0\) gives \(6RS<1/8\).
Hence after the cheap fixing and marking there is at most
one current row. Prescribing the dyadic height boxes costs
only a height-subpower factor.

Count as in (RS5), binning the effective raw masses on a
fixed finite mesh and taking the minimum of the ordinary
row bound and this reconstructed row bound for step \(i\).
The designated row loses at least \(\rho/4-3/R>\rho/6\)
relative to the ordinary exponent before rounding errors.
All choices of \(i\), \(j\), and their order range over a
fixed finite set. Choose the mesh error and the small-budget
allowances below both this improvement and the discard
margin. This proves (RS8). Every valuation cutoff and
positive tolerance was fixed before the inner height limit.

The matching statement cannot require that the paired rows
are identical in an ordered family. The identity
\(F(z)+F(1-z)=1\) realizes a match between a primitive row
and its swapped summands. In a near-maximal even distinct-row
relation the large exclusive pair packets are forced, but
their matched rows can still differ. The existing
full-denominator odd-moment estimates use numerical
numerator and denominator bounds. The argument here gets its
gap from exact comparison modulo the square of the last
prime power, the three radical prices, and the primitive
height box.
[Bloom and Kuperberg, Theorem 3](https://arxiv.org/abs/2312.09021)

For sorted rows and one fixed exact right-hand side, pair
recovery on a split genus-two curve improves this result
at every fixed order. The variable-output extension below
continues to apply when only the outside radical of the
right-hand denominator is small.
[Uniform affine recovery of a distinct pair](ABC_PROJECTIVE_RECIPROCAL_GENUS_TWO_RECOVERY.md#a-saving-at-every-order-of-distinct-active-rows)

## Dense primitive graphs with no additional reciprocal identities

The private-prime realization of a finite additive pattern
already supplies the required tags by prescribing depth one
and choosing all tags above the prescribed moment order.
The complete-grid construction below gives a direct example.
[Private prime powers and exact additivity](ABC_DIFFUSE_FACTOR_GRAPH.md#private-prime-powers-and-exact-additivity)

Fix an integer \(N\ge2\) and a positive stride \(D\) divisible
by every prime up to \(N\). Choose distinct tag primes, all
greater than \(4N\) and coprime to \(D\), for the \(N\)
labels \(a_i=A+Di\), the \(N\) labels \(b_j=B+Dj\), and the
\(2N-1\) labels \(c_k=A+B+Dk\), with \(1\le i,j\le N\) and
\(2\le k\le2N\). We arrange that a tag divides exactly its
own label, exactly to the first power.

At an \(a_i\) tag \(p\), set \(A\equiv -Di+p\bmod p^2\)
and choose \(B\bmod p\) avoiding every zero on the \(b\)
and \(c\) labels. There are at most \(3N-1<p\) forbidden
residues. At a \(b_j\) tag, exchange the roles of \(A,B\).
At a \(c_k\) tag impose \(A+B\equiv-Dk+p\bmod p^2\).
Choose \(A\bmod p\) to avoid zeros on the \(a\) labels
and, with \(B\equiv-A-Dk\bmod p\), the \(b\) labels.
At most \(2N\) residues are forbidden. The congruences at
each tag now prescribe one admissible pair of residues
modulo \(p^2\). In each tagged role the other indices
differ by nonzero multiples of \(D\) modulo \(p\).

First use CRT to choose \(A>0\) meeting its tag congruences
and \(A\equiv1\) modulo the primes up to \(N\). Factor the
product of its \(a_i\). Any non-tag prime of this product
is greater than \(N\). For each one, choose \(B\) modulo
that prime avoiding the at most \(N\) residues that make
some \(b_j\) vanish. Choose \(B\) by CRT with these
conditions and all the prescribed tag congruences. It can
be arbitrarily large, so we may impose \(a_i<b_j\).
Then every \((a_i,b_j,c_{i+j})\) is positive and primitive.
There are exactly \(N^2\) rows using dictionaries of size
\(N,N,2N-1\).

Consider the \(m\)-versus-\(m\) additive energy of \(G\)
with \(m\) smaller than the least tag prime. At the tag of
any one label, every term in the reciprocal equation is
integral except for the occurrences of that label. After
multiplication by the tag and reduction modulo it, the
net multiplicity of the label times a unit is zero.
The net multiplicity has absolute value at most \(m\),
so it is zero. This holds for all three dictionaries.
Conversely, matching all three role multisets gives the
identity. Therefore this complete primitive graph has
only trade identities in every fixed reciprocal moment
once \(N\) is sufficiently large. For example, taking
\(B\ge\exp(N^2)\) makes \(\log N/\log(\max c_k)\to0\). It does not furnish a family saturating a
small radical budget at positive exponent.

## Odd relations with variable output denominator

For a real parameter \(Q\) write \(\operatorname{rad}_{P^c}(Q)\)
for the product of the primes dividing the positive integer \(Q\)
outside a prescribed prime pool \(P\). The parameter \(Q\) below is
the reduced denominator of a rational. Assume \(P\) is supported
on an integer at most \(X^K\), with a fixed \(K\). Fix \(k\),
\(B>0\), and nonzero rational weights \(w_1,\dots,w_k\) of polynomial
height as in (RS1). Enlarge \(P\) by the prime supports of the
weight numerators and denominators. Consider distinct primitive
rows \(y_i=(a_i,b_i,c_i)\) in the shell \(X/2<c_i\le X\), with
\(\operatorname{rad}(a_ib_ic_i)\le n^B\), where \(n=X^L\).
For each tuple put
\[
 t=\sum_{i=1}^k w_i F(a_i/c_i),\qquad Q=\operatorname{den}(t).
\]
In particular \(t\) varies with the tuple. Write \(b_i\) for the
normalized raw radical mass of row \(i\). Let \(m_p\) denote the
number of rows containing \(p\) in their radicals, and define
\[
 H=\sum_{\substack{p\notin P\\m_p=1}}\log_n p,\qquad
 S=\sum_{p\in P}m_p\log_n p+
     \sum_{\substack{p\notin P\\m_p\ge3}}(m_p-2)\log_n p.
\]
For any ordering the new mass \(v_{i,W}\) in role \(W\) is charged
at the first row using a prime outside \(P\), and
\(A=\sum_i(3\max_W v_{i,W}-\sum_Wv_{i,W})\).

The condition \(\operatorname{rad}_{P^c}(Q)\le n^\mu\)
implies \(H\le\mu\). A prime outside \(P\) with a unique
denominator occurrence has a unique term of negative valuation
in the sum and therefore survives in \(Q\). If \(U=\sum_{i,W}
v_{i,W}\) is the outside union mass, the exact accounting is
\[
                      2U=\sum_i b_i+H-S.
\]
Fix normalized bounds \(b_*,S_*,A_*,\mu\). The number of tuples
satisfying \(\sum_i b_i\le b_*\), \(S\ge S_*\),
\(A\ge A_*\) in a prescribed order, and the denominator restriction
above is
\[
                 \ll X^{\{(b_*+\mu-S_*-A_*)/6+\eta\}L+\varepsilon}
                 \qquad (0<L<L_0).                       \tag{RS9}
\]
Here \(\eta,\varepsilon>0\) are arbitrary, and \(L_0>0\) depends
on \(\eta\) and the fixed parameters. For each order and finite
mesh class apply the same sequential (GD3) count as in (RS5).
The step costs are \(\frac13v_i-\frac16(3\max_W v_{i,W}-v_i)\)
up to the allowances. Their sum gives (RS9).
The cost includes no enumeration of the support of \(Q\).

Now fix \(k\ge2\) and \(0<\rho<B/(4(k-1))\).
There are fixed \(\mu_0,\delta,L_0>0\) such that the tuples with
\(\operatorname{rad}_{P^c}(Q)\le n^\mu\), \(0\le\mu\le\mu_0\),
in which the exclusive-pair graph of threshold \(\rho\) fails
to be a perfect matching are
\[
                 \ll X^{\{(kB+\mu)/6-\delta\}L+\varepsilon}
                 \qquad (0<L<L_0).                       \tag{RS10}
\]
The graph uses primes outside \(P\), and a prime is exclusive to
a pair when it occurs on those two rows alone. In particular
(RS10) applies to all such tuples for odd \(k\). Reducing
\(\mu_0\) below \(3\delta\) gives an exponent at most
\(kB/6-\delta/2\). For \(k=1\) one can instead bound the
count by \(X^{\mu L/3+\eta L+\varepsilon}\) using the effective
radical inequality, since the outside primes of the unique row
all survive in \(Q\).

To prove (RS10), fix \(R>36/\rho\), and positive \(s_0,\mu_0\)
so small that
\[
 s_0<\min\{\rho/2,B/16,1/(96R)\},\qquad
 \mu_0<\min\{B/16,1/(16R)\}.
\]
Use (RS9) to discard \(D=kB-\sum_i b_i>s_0\),
\(S>s_0\), or \(A>s_0\) in any of the finitely many
row orders, with a fixed loss from \((kB+\mu)/6\).
Every remaining row has at least \(B-4s_0-\mu\) mass on
exclusive pair primes. The incidence on outside singleton primes
is at most \(\mu\), and that on pool primes and higher-multiplicity
primes is at most \(3S\). Each row thus has a neighbor of
mass at least \(\rho\). If the graph fails to match, choose a row
\(i\) with two neighbors and let \(j\) be the neighbor with
smaller compact exclusive numerical weight \(q\), for cutoff
\(L v_p\le R\). Put \(i\) penultimate and \(j\) last.
The new raw mass at \(i\) is at least \(\rho\), so its
(GD3) step costs at least \(\rho/4\).

Record the full powers on \(i\) carried by the pool or by
\(L v_p>R\), at cost at most \(3/R+o_X(1)\). Their
combined numerical weight is \(p_0\).
Let \(H_{\rm tot}=\sum_W\log_X W_i\). Then
\(q\le(H_{\rm tot}-p_0)/2\) by disjointness of the two
exclusive pair packets. Among the unfixed primes remove those
occurring on row \(j\). Their total numerical weight is at most
\(q+3RS\). Further remove the primes dividing \(Q\). This
costs numerical weight at most \(R\mu\), because the remaining
powers have scaled order at most \(R\). Mark all other
unfixed primes of the current row.

The right-hand side is unknown, but at a marked prime it is
integral. Put \(S_0=-\sum_{\text{earlier }l}w_l F(a_l/c_l)\),
which is determined by the earlier rows. Since the later terms
are integral at the prime, its order \(e\) on the current row
satisfies \(e=-v_p(S_0)\). The congruence (RS7) holds with
\(S_0\) in place of \(S_{\rm rem}\), because their difference
is \(t\) and \(p^e tW\in p^{2e}\mathbb Z_p\).
All the mark records come from the known denominator of \(S_0\),
with a subpower number of choices. Thus the product of the
squared marked moduli and the ordinary fixed factors has
numerical weight at least
\[
 p_0+2(H_{\rm tot}-p_0-q-3RS-R\mu)
          \ \ge\ H_{\rm tot}-6RS-2R\mu .
\]
It exceeds the least pair-height threshold by a fixed margin
because the tallest letter lies in the shell \(X/2<c\le X\).
The current row is therefore reconstructed after the records.
Its charged exponent is at most \(3/R\), giving a loss
at least \(\rho/4-3/R>\rho/6\) before the small mesh and
height allowances. Choose the allowances below this margin
and the discard margin. This proves (RS10) uniformly in
the unknown output numerator and denominator.

The finite-order odd-fraction result of Bloom and Kuperberg
counts integer sums using bounds on the numerical numerators
and denominators. It does not supply (RS10), which uses the
three individual radical prices and squared projective locking
while allowing all outputs with a small outside denominator
radical. The same proof applies to the coordinate of any least
repair, since \(a/c\) is unchanged.
[Bloom and Kuperberg, Theorem 3](https://arxiv.org/abs/2312.09021)
