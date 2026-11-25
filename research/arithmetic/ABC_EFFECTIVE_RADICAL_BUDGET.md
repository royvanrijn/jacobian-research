# The one-third bound for an effective radical budget

The numerical height of an anatomic root and the number of
its possible values enter the small-budget argument in
different ways. The numerical size is needed to make a
multiplicative residue moment diagonal. Once it is diagonal,
the number of root tuples in the moment is a dictionary
cardinality. Keeping these quantities separate gives the
coefficient \(1/3\) for the sum of the actual dictionary
exponents. This permits discounts on prime sets even when
their full powers retain their original numerical weights.

## A diagonal moment with two scales

Consider an anatomic box with fixed initial coefficients,
fixed exponents at each height limit and pairwise coprime
complete terms \(T\). Let \(D\) be the count exponent to
base \(X\). Coordinate \(x_{T,i}\) is in a dyadic interval
at exponent \(\alpha_{T,i}\) and ranges over an arbitrary
dictionary of exponent at most \(\beta_{T,i}\), where
\(0\le\beta_{T,i}\le\alpha_{T,i}\). Here the exponents can
be taken as inner height limits. Put
\[
 H_T=\sum_i\alpha_{T,i},\qquad
 B_T=\sum_i\beta_{T,i},\qquad B=\sum_T B_T .
\]
Use the size convention \(\sum_T H_T\le2\lambda\), with
\(c\) in a height shell comparable to \(X\). The two
tall terms then have numerical weight \(1\) in the inner
limit. The usual decomposition from radical budget
\(\lambda\) satisfies this size convention at small
truncation loss.
[Bernert, Browning, Lichtman and Teräväinen, Proposition 2.1 and Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

Freeze a set of coordinates on one letter \(T\) at
dictionary cost \(f=\sum_{i\ \mathrm{frozen}}\beta_{T,i}\).
If their full product has numerical weight at least
\(5\lambda\), let \(M\) be that product. The two free
complete letters are units modulo \(M\) and are equal
up to sign there. The root dictionaries can be restricted
to units for this record. For a free coordinate write
\[
 F_i(\chi)=\sum_{x\ \mathrm{in\ its\ dictionary}}\chi(x)^i
\]
on the character group of the units modulo \(M\).
For an integer moment order \(m_i\) with
\(m_i\alpha_{T',i}<5\lambda\), fix the first product in
the congruence counted by \(\mathbb E_\chi|F_i|^{2m_i}\).
There are at most \(X^{m_i\beta_{T',i}+o_X(1)}\) choices
of its ordered factors. The other product lies in a
power-root residue class, and the number of such classes
is \(X^{o_X(1)}\) at fixed \(i\). Each class contains at
most one product in the height range for large \(X\),
and each product has \(X^{o_X(1)}\) ordered factorizations.
Thus the moment is \(X^{m_i\beta_{T',i}+o_X(1)}\).
[Multiplicative moments for residue locks, equation (4)](ABC_OMITTED_COORDINATE_GRAPH.md#higher-multiplicative-moments-for-a-residue-lock)

For each free letter \(T'\) with \(H_{T'}>0\), take
\(m_i=\lceil H_{T'}/\alpha_{T',i}\rceil\) at the positive
inner height exponents. Then \(\sum_i 1/m_i\le1\) and
\(m_i\alpha_{T',i}\le2H_{T'}\le4\lambda\). The coordinates
of zero inner height exponent can first be fixed at inner
subpower cost in each copy of the energy. They only
change unit twists. Hölder now bounds the residue energy
of the free letter by \(X^{B_{T'}+o_X(1)}\). When
\(H_{T'}=0\) the same count follows from the entire
dictionary cost being zero. Cauchy–Schwarz between the
two free letters gives
\[
                       D\le f+\frac{B-B_T}{2}.                 \tag{EB1}
\]
Both complete free letters recover the third letter.
All moment orders and the number of coordinates are
fixed before the inner height limit. No bound on their
growth in the subsequent budget limit is needed here.

## The one-third coefficient for coordinate dictionaries

For every fixed \(\eta>0\), at all sufficiently small
fixed \(\lambda\) the dictionary-box count satisfies
\[
                              D\le B/3+\eta\lambda .            \tag{EB2}
\]
The assertion is in the height-limit convention with
\(\sum_T H_T\le2\lambda\), the two tall heights tending
to one, and any fixed truncation at the height limit.
The small-budget threshold is independent of dictionary
cardinalities, coefficient values and the selected
truncation family. The usual \(X^\varepsilon\) allowance
absorbs all fixed-height-limit errors.

Here are the limiting details. If (EB2) fails, pass to a
sequence \(\lambda\to0\) of inner limits with
\(D/\lambda\to d\), \(B/\lambda\to b\),
\(B_T/\lambda\to b_T\) and \(d\ge b/3+\eta\).
All three complete heights have limiting values. Let
\(u=b-2d\), so \(u<d\) and at least one letter has
\(b_T>u\). Call any such letter supported. The limit
\(d\) is positive, since the hypothetical excess is a
fixed positive multiple of \(\lambda\).
The following proof uses the same cases as the
[original one-third argument](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13).

A supported letter reduces to boundedly many high
coordinates and one short low cofactor. To see this,
fix \(\theta=(b_T-u)/2>0\). By (EB1), any collection
of that letter's coordinates of numerical weight
at least \(5\lambda\) has normalized dictionary
cost at least \(\theta-o(1)\). List the coordinates
with \(\beta_{T,i}/\lambda\ge\theta/2\).
There are boundedly many since \(B_T\le2\lambda\).
Pass to a subsequence where each listed degree is
fixed or tends to infinity. Retain the latter as
individual high coordinates. Each has actual root
exponent at least \(\theta\lambda/2\), so its full
power eventually has numerical weight at least
\(5\lambda\), and its normalized dictionary cost
is at least \(\theta-o(1)\).

Choose a sufficiently large fixed \(J\) greater than the
bounded listed degrees. Every unretained coordinate of
degree at least \(J\) has normalized dictionary cost
less than \(\theta/2\). If their total numerical weight
reaches \(5\lambda\), take an initial segment stopping
when it first does so. The preceding part costs at
most \(5\lambda/J\), since \(\beta_{T,i}\le\alpha_{T,i}\).
Including the last coordinate still costs less than
\((\theta-o(1))\lambda\), a contradiction for large
enough \(J\). All unretained factors therefore have
total numerical weight \(O_J(\lambda)\), including
the factors of degree below \(J\) by the raw height
sum. Group them as a single low cofactor. It and the
finitely many individual roots are of height \(X^{O_J(\lambda)}\).
At least one retained root on a supported tall letter
has positive limiting numerical weight.

A full modulus of positive limiting numerical weight
and minimum prime valuation \(\gg1/\lambda\) permits
the bounded-root modular release (FC2) on any two
supported letters at normalized cost tending to zero.
This uses their actual root and low-cofactor heights.
It also permits a supported letter and a single
unreleased root of height \(X^{O(\lambda)}\).
The complete unit ratio determines both free letters by
primitivity. When two tall letters of these forms have
ratio \(1+O(X^{-c_0})\) for fixed \(c_0>0\), the
short-strip release (SR3) has the same normalized cost.
Individual roots on distinct letters can instead be
released by (30) after all other coordinates are fixed,
provided their degrees tend to infinity. Grouped
cofactors stay coprime to retained roots since the
anatomic factors are pairwise coprime; the fixed initial
coefficients give fixed twists.
[Modular release (FC2) and real release (SR3)](ABC_OMITTED_COORDINATE_GRAPH.md#signed-monomial-congruences-over-full-moduli)
· [Bounded high-root release (30)](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

**All three letters supported.** Choose a positive-weight
retained root on each of the two tall letters. The degree
of such a root is \(\gg1/\lambda\) because its raw height
exponent is at most \(2\lambda\). Each of the two roots
must cost at least \(d\) in normalized dictionary exponent:
fix it as the modulus and use (FC2) on the two other
supported letters. Conversely, release both chosen roots
from their exact additive equation using (30), fixing all
other coordinates. This gives \(d\le b-2d\), contradicting
\(u<d\).

**Two letters supported.** Write \(S\) for the third, so
\(b_S\le u<d\). If \(S\) has positive limiting height \(h_S\),
fix its coordinates of degrees at least \(h_S/(4\lambda)\).
Their full product has positive limiting weight, as the
lower degrees carry at most \(h_S/2\) by the raw height
sum. This costs at most \(b_S\). Apply (FC2) to the two
supported letters to get \(d\le b_S\).
If \(S\) has limiting height zero, the supported letters
are the two tall ones. Their signed ratio has a real
short-strip precision of fixed positive numerical weight.
Equation (SR3) gives \(d=0\). Both cases contradict
the choice of limits.

**One letter supported.** Write \(U,V\) for the two
unsupported letters with \(b_U,b_V\le u\).
For a positive-height \(U\), consider its high
coordinates of degrees at least \(h_U/(4\lambda)\).
They have positive collective weight. There must be
one such coordinate with positive limiting weight
and positive normalized dictionary cost that can be
selected as follows. If every high coordinate has
cost at most \((d-b_V)/2+o(1)\), group the coordinates
in order into blocks of cost at most
\((d-b_V)/2+\xi\) for a sufficiently small fixed
\(\xi>0\). Consecutive completed blocks have combined
cost above this cap. The number of blocks is bounded,
and one block has positive limiting weight. Fixing
\(V\) and that block gives cost strictly less than
\(d\), and (FC2) releases the supported letter modulo
the block. This is a contradiction. Consequently one
high coordinate has normalized cost bounded away
from zero. Since \(\alpha_{U,i}\ge\beta_{U,i}\) and its
degree is \(\gg1/\lambda\), it has positive limiting
numerical weight. Write its normalized cost as
\(\gamma_U\), and similarly choose one with cost
\(\gamma_V\) on \(V\) when that letter has positive height.

If both heights are positive, lock the selected root
of \(U\) and freeze \(V\) apart from its selected root.
Equation (FC2) counts the supported letter and that
root at cost tending to zero. Interchanging \(U,V\)
gives
\[
 d\le\gamma_U+b_V-\gamma_V,\qquad
 d\le\gamma_V+b_U-\gamma_U .
\]
The sum contradicts \(2d>2u\ge b_U+b_V\).
If \(V\) has height zero, the supported letter and \(U\)
are tall. Freeze the complement of the selected root of
\(U\) and use (SR3), giving \(d\le b_U-\gamma_U\).
Fixing \(V\) and using the selected root of \(U\) as
modulus in (FC2) gives \(d\le b_V+\gamma_U\).
Their sum is the same contradiction. These are all
possibilities since \(3u<b\).

This proves (EB2). The moments were taken at fixed
truncation before the height limit. The bounded number
of released coordinates, the positive numerical margins,
and the comparisons with \(\lambda\) were chosen in
the subsequent budget argument.

## A strict effective saving when one letter is free to record

Fix \(b_->0\). There is a saving depending only on
\(b_-\) in the dictionary-box exponent when one letter
has all its dictionary costs zero in the inner limit
and \(B\ge b_-\lambda\). In this range,
\[
                         D\le B/3-\delta(b_-)\lambda
                                                                    \tag{EB9}
\]
holds at sufficiently small fixed budgets, before the
inner exponent allowance. The zero-cost letter can
vary through a subpower dictionary.

Suppose there is a sequence with the gap tending to
zero. Normalize by \(\lambda\), use the notation in the
proof of (EB2), and pass to \(d=b/3>0\). The support
cutoff is \(u=d\). Some letter has cost above \(d\),
since one has cost zero. The supported decomposition
and its releases apply for every strict \(b_T>u\), so
they also apply at this equality. If two letters are
supported, the zero-cost third letter has cost strictly
below \(d\), and the same modulus or short-strip argument
gives \(d=0\).

In the remaining case name the supported letter \(T\),
the zero-cost letter \(V\), and the other letter \(U\),
with \(b_U\le d\). If \(U\) has positive limiting
height, the high-block argument in (EB2) produces a
positive-weight high root of cost \(\gamma_U>0\)
on \(U\). In that argument one first fixes all of \(V\),
whose cost is zero, so its strict margin is \(d>0\).
If \(V\) has positive height, its high part supplies
a positive-weight full modulus at zero cost. Fix the
complement of the selected root on \(U\). By (FC2),
\(T\) and that root can be released together, yielding
\(d\le b_U-\gamma_U<d\). If \(V\) has zero height and
\(U\) positive height, \(U\) and \(T\) are tall and the
same cost inequality follows by (SR3). If \(U\) has
height zero, \(V\) and \(T\) are tall. Fix \(V\) at zero
cost and release \(T\) in the real short strip, giving
\(d=0\). This excludes the equality sequence and yields
a uniform positive gap by the fixed lower bound \(b_-\).

One can use an upper budget \(B_0\) in place of \(B\)
whenever \(B_0\ge b_-\lambda\). If \(B<(b_-/2)\lambda\),
(EB2) supplies the gap from the budget comparison, with
a smaller allowance. Otherwise apply (EB9) with \(b_-/2\).
This also applies to the charged-prime count below when
one actual letter is supported on primes of price zero.
Every root on that letter then has a subpower dictionary.
The coefficient and box losses can be chosen smaller
than the gap.

## Counting by a sublogarithmic prime partition sum

For each \(X\) let primes \(p\le X\) have prices
\(0\le w_X(p)\le\log_X p\). Suppose the partition sum
\[
                       \sum_{p\le X}X^{-w_X(p)}
                                  \le g(X)=o(\log X)            \tag{EB3}
\]
is uniform over the allowed choices of prices. Put
\(\mathfrak b_X(N)=\sum_{p\mid N} w_X(p)\).
An integer coordinate of height at most \(X\) and
price at most \(s\) ranges over at most
\(X^{s+o_X(1)}\) values, uniformly in the prices and \(s\).
Indeed its possible squarefree supports number at most
\[
 X^s\prod_{p\le X}(1+X^{-w_X(p)})
                  \le X^s\exp(g(X)).
\]
On any one support of a value at most \(X\),
the exponent vectors cost \(X^{o_X(1)}\),
by the uniform fixed-support enumeration. Equivalently,
the support has \(O(\log X/\log\log X)\) primes and
stars and bars counts the possible exponents whose sum
is \(O(\log X)\).
[Uniform fixed-support enumeration](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

Consequently for every \(\eta>0\) and all sufficiently small
fixed \(\lambda>0\), the number of primitive triples in
\(X/2<c\le X\) with \(\operatorname{rad}(abc)<c^\lambda\)
and \(\mathfrak b_X(abc)\le b_0\) satisfies
\[
                    \ll_{\lambda,\eta,\varepsilon,g}
                         X^{b_0/3+\eta\lambda+\varepsilon},
                    \qquad 0\le b_0\le\lambda.                  \tag{EB4}
\]
The bound is uniform in the prices and \(b_0\) under
(EB3). To see this, use the BBLT anatomic factorization
with initial coefficient height at most \(X^\rho\) and
raw root sum at most \(\lambda+O(\rho)\), with
\(\rho=o(\lambda)\) in the outer limit. Fix the
three coefficients at cost \(3\rho\) and use dyadic
root intervals at subpower cost. All positive-index
coordinates across the three letters are pairwise
coprime. Their total price is at most
\(\mathfrak b_X(abc)\), regardless of possible
overlap with the initial coefficients. Box the coordinate
prices into intervals of sufficiently small length,
using a total upper endpoint error \(o(\lambda)\);
at fixed truncation this introduces only a bounded
number of price bins independent of \(X\).
The root enumeration gives dictionary exponents
at most these endpoints, while the raw root caps
remain the same. One may take their coordinatewise
minimum. Apply (EB2), with a smaller allowance, and
sum over the boxes. The truncation error, binning
and inner subpower errors fit in (EB4).
[Anatomic decomposition](https://arxiv.org/html/2410.12234v2#S2)

## Fixed supports and polynomial prime prices

Let \(\mathcal P_X\) be prescribed, and suppose
\(\sum_{p\in\mathcal P_X}p^{-a}=o(\log X)\) for some
fixed \(0\le a<1\), with a uniform bound. Price those
primes at \(a\log_X p\) and the others at \(\log_X p\).
Condition (EB3) follows from this hypothesis and
\(\sum_{p\le X}1/p=O(\log\log X)\). If a fixed
fraction \(\chi>0\) of the radical budget is carried
on the pool, in the sense that
\(\sum_{p\mid abc,\ p\in\mathcal P_X}\log_X p\ge\chi\lambda\),
then (EB4) gives
\[
            X^{\{1-(1-a)\chi\}\lambda/3+\eta\lambda+\varepsilon}.
                                                                    \tag{EB5}
\]
This is a strict saving for \(\eta<(1-a)\chi/3\).
Every prescribed pool of size \(O((\log X)^\kappa)\)
for fixed \(\kappa>0\) satisfies the condition with
\(a>\max(0,1-1/\kappa)\), \(a<1\), by ordering its primes
and bounding their sum by \(\sum_{j\le O((\log X)^\kappa)}
(j+1)^{-a}\). Thus a positive fraction of the radical
budget in any fixed polylogarithmic prime alphabet is
sufficient, even when most of the radical is outside it.
This includes the earlier
[three-divisor pool sector](ABC_SMOOTH_SUPPORT_DIVISOR_PACKING.md).

In the polylogarithmic case one may take the limiting
price multiplier \(a=1-1/\max(1,\kappa)\) in the
exponent, at the cost of an arbitrarily small fixed
addition to \(\eta\). More generally suppose
\(\#\{p\in\mathcal P_X:p\le y\}\le C y^\alpha\) for all
\(y\ge2\) and some fixed \(0\le\alpha<1\).
Summation by parts gives \(\sum_{p\in\mathcal P_X}p^{-a}
=O_{C,a,\alpha}(1)\) for every \(a>\alpha\). Thus the
pool can be polynomially large. A proportion \(\chi\)
of the budget on that pool gives coefficient
\((1-(1-\alpha)\chi)/3\) up to arbitrarily small
fixed allowance at sufficiently small budgets.

If the primes in the pool belong to one prescribed
integer of height \(X^{O(1)}\), then \(a=0\) is allowed.
For a specified complete letter \(T=t\) in the shell
with \(\operatorname{rad}(t)\ge X^{\chi\lambda}\)
this gives the uniform fixed-value exponent
\[
                         (1-\chi)\lambda/3+\eta\lambda.
                                                                    \tag{EB6}
\]
The same conclusion applies to any fixed polynomial-height
list of divisors of bounded length whose *union* of prime
supports contributes logarithmic radical weight
\(\chi\lambda\) to every triple. At the balanced endpoint
a complete letter has radical exponent
\((1/3+o(1))\lambda\). Thus every fixed complete-value
fibre there has exponent at most \((2/9+o(1))\lambda\).
When the zero-price pool contains all primes of the fixed
letter, (EB9) also gives an absolute gap below the \(2/9\)
bound on each fixed-value fibre at balance, by taking
\(b_-=1/2\) and sufficiently small balance allowances.
A saturated balanced family therefore has dictionary
exponent at least \(1/9+\delta_1\) in base \(X^\lambda\)
for each letter, for some absolute \(\delta_1>0\)
at sufficiently late budgets. On the uniform family
the probability of any individual marginal value is
at most \(X^{-(1/9+\delta_1)\lambda}\) after reducing
\(\delta_1\) to absorb the height allowances.

For a different application put
\(K=\log X/(\log\log X)^2\) and choose
\[
                    w_X(p)=\log_X\max(1,p/K).
\]
The price partition sum is \(o(\log X)\):
primes at most \(K\) cost at most \(K\) terms and the
others contribute at most \(K\sum_{p\le X}1/p\).
For a support of \(k\) primes the discount is at least
\(k\log_X K-K\log_X K=k\log_X\log X-o_X(1)\),
using \(k=O(\log X/\log\log X)\) for \(abc\le X^3\).
It follows that for each fixed \(0<u\le1\) the additional
constraint
\[
                   \omega(abc)\frac{\log\log X}{\log X}
                                      \ge u\lambda
\]
gives exponent \((1-u)\lambda/3+\eta\lambda+\varepsilon\)
at sufficiently small fixed budgets, for arbitrarily small
fixed positive \(\eta\). There is no upper bound on the
largest prime in this sector.

## Occurrence frequencies and anatomic entropy at saturation

For any given shell family \(\mathcal F\), put
\(\pi_X(p)=\Pr_{\mathcal F}(p\mid abc)\) on its uniform
measure and take
\[
                       w_X(p)=\min(\log_X p,-\log_X\pi_X(p)).
\]
If \(\pi_X(p)=0\), use the first term. Then
\(\sum X^{-w_X(p)}\le \sum p^{-1}+\mathbb E_{\mathcal F}
\omega(abc)=O(\log X/\log\log X)\) uniformly. Thus the
prices can depend on the observed family in (EB4).
For every fixed \(\zeta>0\), the subfamily of that family
with
\[
                \sum_{p\mid abc}\bigl(\log_X(p\pi_X(p))\bigr)_+
                                               \ge\zeta\lambda
                                                                    \tag{EB7}
\]
satisfies (EB4) with \(b_0=(1-\zeta)\lambda\) when
\(0<\zeta\le1\). The discount event is empty
if \(\zeta>1\). In particular, at a limiting saturation
sequence, for each fixed positive \(\zeta\) such rows
form a power-smaller fraction at the fixed sufficiently
late budgets, with the allowance in (EB4) taken smaller
than \(\zeta/3\). For every fixed \(a<1\), the primes
occurring with frequency at least \(p^{-a}\) therefore
carry \(o(\lambda)\) typical logarithmic radical mass
in the iterated endpoint limit. The same assertions apply
with occurrence frequencies recomputed on any subfamily
that is itself near saturation.

The frequency statement is two-sided in mean. Write
\(Z=\sum_p\pi_X(p)\), \(J=\sum_{p\le X}p^{-1}\), and sum
only over the primes of positive frequency when logarithms
of \(\pi_X(p)\) are taken. The log-sum inequality, applied
to the weights \(\pi_X(p)/Z\) and \((1/p)/J\), gives
\[
  \sum_p\pi_X(p)\log\frac1{p\pi_X(p)}
                   \le Z\log(J/Z)\le J/e .
\]
By (EB7) and the bound \(\sum_{p\mid abc}\log p\le
\lambda\log X\), the mean positive part
\(\sum_p\pi_X(p)(\log(p\pi_X(p)))_+\) is
\(o(\lambda\log X)\) at saturation. The displayed log-sum
bound implies the same for the mean negative part.
Consequently
\[
             \sum_p\pi_X(p)|\log(p\pi_X(p))|
                                     =o(\lambda\log X).
                                                                    \tag{EB8}
\]
In particular, for any fixed \(\mu>0\), primes of frequency
outside \([p^{-1-\mu},p^{-1+\mu}]\) carry \(o(\lambda)\)
mean logarithmic radical mass in base \(X\). Markov's
inequality gives the corresponding typical statement
with loss tending to zero. The one-sided power saving
in (EB7) is stronger for its own exceptional rows.

There is also a coordinate entropy statement. In an
anatomic box write \(P_i\) for the marginal distribution
of coordinate \(i\) on a uniform family of size \(X^{D+o_X(1)}\).
On each row set \(J_i=-\log_X P_i(x_i)\) and consider
\(S=\sum_i\min(J_i,\alpha_i)\). At a fixed truncation
bin the \(J_i\) in intervals with total upper endpoint
error below any fixed positive multiple of \(\lambda\).
A bin with upper endpoint \(s_i\) contains at most
\(X^{s_i}\) values of \(x_i\). Its dictionary exponent is
therefore at most \(\min(s_i,\alpha_i)\) up to inner errors.
By (EB2), for every fixed \(\zeta>0\), the part of the
box with \(S<3D-\zeta\lambda\) is power-smaller at
sufficiently small budgets, using allowance below
\(\zeta/3\). Since \(S\le\sum_i J_i\), taking means gives
\[
                       \sum_i H(x_i)\ge
                           (3D-o(\lambda))\log X
\]
in the limit sense. At saturation the raw root-height
caps sum to \((1+o(1))\lambda\), and \(3D=(1+o(1))\lambda\).
Thus a typical row has
\(\sum_i(\alpha_i-J_i)_+=o(\lambda)\), and the total
marginal entropy deficit from the raw caps tends to
zero on the \(\lambda\log X\) scale. The joint entropy
is the logarithm of the actual family size.

## Consequences for the least target repair

Let \(M=\lceil K_0/\lambda\rceil\) for a fixed \(K_0>0\)
and use the central radical sector of the least target
repair. Thus \(R_c\ge c^{(1/3-\xi)\lambda}\),
\(R_c\le c^{(1/3+\xi)\lambda}\), and
\(R_aR_b\le c^{(2/3+\xi)\lambda}\). In this sector
\[
                  H_M(c)\ge R_c^M,\qquad
                  c\le Y^{1/[M(1/3-\xi)\lambda]}
                  \quad\hbox{if }H_M(c)\le Y.
\]
If a prescribed prime pool with the partition condition
above carries at least a fixed positive fraction
\(\chi\lambda\) of the logarithmic radical in base \(c\),
(EB5) gives a primitive quotient count
\(O(Y^{(1-\sigma)/M+\varepsilon})\) for some fixed
\(\sigma>0\) at sufficiently small \(\xi\) and \(\lambda\).
Here use budget \((1+3\xi)\lambda\) on the original triples.
The saving survives the height change by choosing
\(\xi\) and the exponent allowance sufficiently small
relative to \((1-a)\chi\). Apply the estimate shell by
shell when the partition bound holds on each original
height shell. For a fixed polylogarithmic pool at scale
\(Y\), split off \(c\le Y^{1/(4M(1+3\xi)\lambda)}\)
by the elementary bound of exponent \((1+3\xi)\lambda\);
on larger shells \(\log Y\ll_{\lambda,M}\log c\)
and the same pool is admissible.

For each given primitive target \(c\) in the central sector,
(EB9) with that same budget gives at most
\[
                      Y^{(2/3-\sigma_0)/M+\varepsilon}
\]
distinct quotients with \(H_M(c)\le Y\), uniformly in \(c\),
for an absolute \(\sigma_0>0\) and sufficiently small
absolute \(\xi\). The price cap at height \(c\) is at
most \((2/3+4\xi)\lambda\), so the zero-price target
saves a fixed positive multiple of \(\lambda\) from
one third of this cap. Division by
\((1/3-\xi)\lambda\) gives the claimed repair exponent
for small enough \(\xi\). Both counts use the primitive
quotient recovered by the least multiplier.
[The least target repair and its height](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

All these conditions allow a family whose original
radical primes are sufficiently rare in the uniform
edge measure. To settle the unrestricted bound one
still needs a saving for the configurations without a
positive effective discount, or a separate argument
forcing such a discount in every saturated family.
