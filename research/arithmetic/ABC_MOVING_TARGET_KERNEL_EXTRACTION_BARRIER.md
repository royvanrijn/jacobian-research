# A power-class kernel barrier for moving-target codegrees

The moving-target graph estimates pay for the relative power class
of the target, which can be much cheaper than the two marginal
classes. The three original radical caps also give a lower bound
on the price of an eligible-pair population. When almost all
radical mass lies on valuations not divisible by the chosen
power order, this lower bound reaches the unsaved \(1/3\)
exponent. The bound also applies to the target-exclusive record
and to counting all distinct pairs within a column.

## A two-letter count with prescribed prime pools

Work in a height shell \(X/2<c\le X\) of primitive positive
triples \(a+b=c\), and put \(n=X^L\). The height limit is taken
at fixed \(L>0\) before the small-budget limit. For \(U=A,B,C\)
let \(\ell_U\) be a bounded normalized upper cap for the
complete radical of the letter, and let \(S_U\) be a prescribed
prime pool with product at most \(X^K\). Suppose the radical
outside \(S_U\) of this letter is at most \(n^{d_U}\), where
\(d_U\ge0\). Write
\[
        \Psi(d_A,d_B,d_C)
          =\tfrac12\min\{d_A+d_B,\ d_A+d_C,\ d_B+d_C\}.
\]
For each \(\eta>0\), all sufficiently small fixed \(L\), and
each \(\varepsilon>0\), the count is
\[
                  O\bigl(n^{\Psi(d_A,d_B,d_C)+\eta}
                          X^\varepsilon\bigr).                \tag{MK1}
\]
The constants are uniform in the pools and the bounded caps.
This is the two-letter dictionary inequality with fixed
factors priced at zero.

For completeness, use the exact-slot refinement of the
anatomic decomposition with parameter sufficiently small
relative to \(L\eta\). Fix its small exceptional coefficients
and its one packed high root on each letter at a total cost
less than any prescribed fraction of \(\eta\) in base \(n\).
The packed roots remain fixed original coordinates, with
zero dictionary cost in (GD3). For each remaining squarefree root
position split the root into its part supported on \(S_U\)
and its part outside \(S_U\). The supported parts of the finite
list of roots can be fixed by uniform supported-value
enumeration, at inner subpower cost. They remain portions of
the original roots for the application of (GD3), so possibly
large in-pool complete powers do not become exceptional
coefficients. Dyadically bin the outside factors. Their
total root-height cap on \(U\) is at most \(d_U L\), up to
inner bin errors, since the slots are squarefree with
disjoint support. The original root in each slot now has
dictionary cardinality at most the range size of its outside
factor, since its in-pool factor is fixed. The packed root
and any other fixed positions have zero dictionary cost.
Apply (GD3) to the singleton dictionaries. Sum over the
truncation and support-fixing choices and choose their
allowances small relative to \(\eta\). This proves (MK1).
[Two-letter dictionary count (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)
· [Exact valuations in the height-first convention](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)
· [Polynomial-height support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## An anchored upper bound for eligible queries

Fix an integer \(q\ge2\). For a positive integer \(z\) set
\[
        R_q^0(z)=\prod_{\substack{p\mid z\\q\mid v_p(z)}}p.
\]
This records only the primes at positive orders divisible by
\(q\). Besides the complete caps \(n^{\ell_U}\), impose
\(R_q^0(U)\le n^{s_U}\) on the relevant letters, with
\(0\le s_U\le\ell_U\). The order may be fixed arbitrarily in
each box. The estimates below do not use any uniform counting
theorem for a weighted degree-\(q\) equation.

Take an anchor row \(i\) and count query rows \(j\) in its
group. A prescribed pool \(\mathcal P\), of polynomial prime
product, may depend on the anchor and the group. It can also
include primes of the query's own summand kernels for tests
on its target, since those primes are absent from that
target. Alternatively any additional pool choices from a
dictionary may be summed at their explicit cost.

**MK2.** Let the paid union radical be at most \(n^\beta\)
under either of the following hypotheses.

- The group is a column of a fixed right \(q\)-power class,
  so \([B_j]=[B_i]\) in
  \(\mathbb Q_{>0}^{\times}/(\mathbb Q_{>0}^{\times})^q\).
  The record prices either the outside relative target
  defect, or the target defect on primes exclusive to one
  of the two complete targets. It suffices to price just
  the query-exclusive part for this upper bound.
- The right class is allowed to move. The record prices
  the union of its outside relative defect and the outside
  target-exclusive defect, as in (MT1).

Then the number of eligible query rows is at most
\[
             O\bigl(n^{\beta+\Psi(\ell_A,s_B,s_C)+\eta}
                        X^\varepsilon\bigr).                 \tag{MK2}
\]
The query may have the same left label as the anchor. No
packing condition on the roots is needed for this bound.
If the query-dependent pool in the second case contains
its own left-kernel primes, these can be omitted since
they divide neither the query's right term nor its target.
A query-dependent right-kernel pool cannot be omitted
from the relative right test in this manner.

In the first case, every prime on the query's right letter
whose valuation is not divisible by \(q\) belongs to the
support of the fixed right class. In the second case such
a prime outside the actual anchor's right support and the
fixed pool is necessarily in the relative right record.
Similarly any prime on the query target whose positive
valuation is not divisible by \(q\) lies in the anchor
target support, the fixed pool, or the paid target record.
This remains true if only the exclusive part is priced,
irrespective of the differences of depths at common
target primes.

Enumerate the paid union radical at cost \(O(n^\beta)\).
For that record take the known right support and the paid
primes and fixed pool as the right pool. Take the anchor's
complete target support and the paid primes and fixed
pool as the target pool. Outside these pools the two
query letters have only orders divisible by \(q\), so their
outside caps are \(s_B,s_C\), respectively. Apply (MK1)
to these two pools and the empty pool on the other letter.
All pool products are of polynomial height uniformly in
the data. This proves the bound. The argument applies
with either orientation of the signed additive equation
among the original letters. If the right class is allowed
to move but its relative defect is not included in the
record, the right-support step is lost; in general this
proof then has \(\ell_B\) in place of \(s_B\).

## A necessary aligned radical mass for a strict graph bound

Let \(W\) be the number of actual rows and \(B_1\) the number
of groups in the graph. Write their normalized exponents
as \(d=\log_n W\), \(b=\log_n B_1\), and write \(a\) for
the normalized exponent of the left \(q\)-class dictionary.
Suppose a fraction \(\rho\ge n^{-\theta}\) of the ordered
pairs under consideration satisfies the record cap
\(n^\beta\) uniformly. For (MP6) the pairs have distinct
left labels in the same column. For (MT2) they can have
equal left labels and include all distinct row pairs in
each group. With \(M_0\) bounding the actual row
multiplicity per kernel pair, assume \(W\ge2M_0B_1\) in
the first case; otherwise the linear term of (MP6)
already reaches the observed row count. In the second
case assume \(W\ge2B_1\).

**MK3.** Under the pool conditions in (MK2),
\[
           \beta+\theta
                \ge d-b-\Psi(\ell_A,s_B,s_C)-o(1),             \tag{MK3}
\]
where the small error represents arbitrarily small fixed
normalized allowances followed by the height-first limit.
In the moving-right application with all distinct pairs
allowed in one group, one may set \(b=0\).

Indeed the pair population before eligibility is at least
\(W^2/B_1-M_0W\) in the first case or \(W^2/B_1-W\) in the
second. The respective assumptions make it \(\gg W^2/B_1\).
On the other hand the sum of (MK2) over anchors is at
most \(W n^{\beta+\Psi(\ell_A,s_B,s_C)+\eta}
X^\varepsilon\). Compare these bounds and take
\(\varepsilon\) sufficiently small after fixing \(L\).
[Moving-target graph count (MP6)](ABC_MOVING_TARGET_POWERFUL_CODEGREES.md#counting-all-original-rows)
· [Graph count with equal left labels (MT2)](ABC_MOVING_RIGHT_AND_TARGET_POWER_CLASS_CODEGREES.md#all-original-pairs-in-a-column)

To express the consequence at balance, also assume the
\(s_A\) cap on the left letter. On fixing a left
\(q\)-class its nondivisible support is prescribed, so (MK1)
on that pool gives
\[
                 a\ge d-\Psi(s_A,\ell_B,\ell_C)-o(1).
\]
This bound applies when the left labels are distinct
power-class dictionary entries, including the canonical
full kernels. Even in the favourable regime where the
angular-spacing bracket in (MP2) or (MT1) is bounded and
there is no additional diagonal equation exponent, the
graph square-root term has exponent at least
\[
 \begin{split}
 a+\tfrac12 b+\tfrac12(\beta+\theta)
      \ \ge\ \tfrac32d-\Psi(s_A,\ell_B,\ell_C)
                         -\tfrac12\Psi(\ell_A,s_B,s_C)-o(1).
                                                                    \tag{MK4}
 \end{split}
\]
If \(d\to1/3\) and \(\ell_A,\ell_B,\ell_C\to1/3\), its
lower bound is
\[
              \tfrac13-\tfrac12 s_A-\tfrac14(s_B+s_C)-o(1)
\]
using the upper bounds \(s_U\le\ell_U\). In particular if
all three \(q\)-divisible radical masses vanish, such an
eligible-pair estimate cannot lower the exponent below
\(1/3\) by a fixed amount. Conversely a fixed margin
\(\delta>0\) in the square-root bound requires
\(2s_A+s_B+s_C\ge4\delta-o(1)\). If the left dictionary
is already at exponent \(1/6+o(1)\), the same comparison
requires \(s_B+s_C\ge4\delta-o(1)\), without an assumption
on \(s_A\).

This is a necessary condition on the record price in these
codegree arguments. The estimate retains the three
original radical caps but does not rule out a stronger
method using the numerical weights or more information
about the eligible pairs. It also does not give a
positive lower bound on the alignment mass for some
common order in every balanced family.

## An anchored test for many-row stars

Fix an integer \(m\ge2\) before the height limit. In a group
with one right \(q\)-power class, take an ordered star of \(m\)
distinct actual rows. Fix its first \(m-1\) rows as anchors.
Suppose the paid record includes the outside defect radical
of
\[
                         G^{w_0}\prod_{j=1}^m T_j^{w_j},
\]
where \(G\) divides the first anchor's actual right letter
and \(\gcd(w_m,q)=1\). For instance \(G\) may be a common
right divisor that depends on the whole star. The weights
can be prescribed by roles; their magnitudes do not enter
the proof.

A nondivisible prime on the last target outside the anchor
target supports, the first anchor's right support, and the
prescribed pool must occur in this defect. Indeed the other
factors are units at that prime and \(w_m\) is invertible
modulo \(q\). Thus the last target outside those fixed
supports and the paid primes has cap \(n^{s_C}\).
The right letter outside the column's kernel support
has cap \(n^{s_B}\). Applying (MK1) after enumeration
of the paid union gives the same bound (MK2) for the last
query. This argument also covers a moving-right group
if the paid union reveals the last nondivisible right primes
outside the anchors' right supports. A relative right
class with an anchor suffices, as does a product of right
classes with a unit weight on the query at order \(q\).

**MK5.** With \(d_j\) the actual row count in group \(j\),
suppose a fraction at least \(n^{-\theta}\) of its total
ordered distinct \(m\)-stars, summed over groups, is
eligible for the above test with paid union cap \(n^\beta\).
Then the same necessary inequality (MK3) holds, provided
\(W/B_1\) is large in terms of \(m\).

Write \(\mathcal S\) for the total star count and set
\(\psi=\Psi(\ell_A,s_B,s_C)\). The anchored query bound
gives at most \(O(n^{\beta+\psi+\eta}X^\varepsilon
\sum_j d_j^{m-1})\) eligible ordered stars.
The falling factorial sum \(\mathcal S=\sum_j(d_j)_m\)
is comparable to \(\sum_j d_j^m\) in the stated range.
In fact for \(d_j\ge2m\) the falling factorial is at least
\(2^{-m}d_j^m\), whereas the other groups contribute
\(O_m(B_1)\) to the ordinary moment. The power mean
bounds that moment below by \(W^m/B_1^{m-1}\).
Hölder now gives
\(\sum_jd_j^{m-1}\ll_m B_1^{1/m}\mathcal S^{(m-1)/m}\).
Dividing the eligible count by \(\mathcal S\) yields
\[
             \rho\ll_m n^{\beta+\psi+\eta}X^\varepsilon
                                   B_1/W .
\]
This is (MK3) at the exponent level. The inequality
concerns the local record price for a large star
population and does not assume any particular
completion estimate for those stars.

## The exponent cost of larger stars

**MK6.** In (MS3), allow two orders. Let \(e\) be the
omitted root order and \(q\mid e\) the prime target class
order. Take \(s_A,s_B\) to bound the normalized \(e\)-divisible
radicals and \(s_C\) to bound the normalized \(q\)-divisible
radical. The anchored argument for (MK5) still applies.
It uses the right class at order \(e\) and the target record
at order \(q\). In a moving-right group the relative right
charge must also be taken at order \(e\).

For the actual right-column dictionaries of (MS3), apply
(MK1) separately to each fixed left and right \(e\)-class.
Writing \(a,b\) for their dictionary exponents gives
\[
 \begin{split}
 a&\ge d-\Psi(s_A,\ell_B,\ell_C)-o(1),\\
 b&\ge d-\Psi(\ell_A,s_B,\ell_C)-o(1).
 \end{split}
\]
Assume the star bracket is bounded. With the defect cap
\(n^\beta\) and eligible fraction \(n^{-\theta}\), the
exponent of the second term in (MS3) is bounded below by
\[
 \begin{split}
 a+(1-1/m)b+(\beta+\theta)/m
 &\ge (2-1/m)d-\Psi(s_A,\ell_B,\ell_C)\\
 &\quad -(1-2/m)\Psi(\ell_A,s_B,\ell_C)
                 -\Psi(\ell_A,s_B,s_C)/m-o(1).
 \end{split}
\]
This follows by substituting (MK5) and the two dictionary
lower bounds. At balanced saturation the right side is
\[
       \tfrac13-\tfrac12s_A
                     -\bigl(\tfrac12-\tfrac1{2m}\bigr)s_B
                     -\tfrac{s_C}{2m}-o(1).
\]
If the dictionary exponents are already \(1/6+o(1)\), then
(MK5) alone gives
\(\tfrac13-(s_B+s_C)/(2m)-o(1)\). Thus increasing the star
size still requires alignment mass for a fixed margin, using
the full three radical caps. The right dictionary lower
bound here refers to actual right \(e\)-classes. The local
price (MK5) is the part that applies to arbitrary groupings
with relative right defects.
[Actual all-row star estimate (MS3)](ABC_MOVING_TARGET_POWER_CLASS_STARS.md#an-actual-all-row-graph-estimate)

## Sequential support loss in larger stars

**MK7.** In the fixed right \(e\)-class groups of (MS3), let
\(q\mid e\) be prime and let \(m\ge2\) be fixed. Use complete radical
caps \(n^{\ell_A},n^{\ell_B},n^{\ell_C}\), and divisible-order
caps \(R_e^0(B)\le n^{s_B}\) and \(R_q^0(C)\le n^{s_C}\).
In each group fix a prime pool \(P_0\) of polynomial prime
product, which may depend on the group and its prescribed
modulus. In a star of ordered distinct rows, allow the free
pool to be
\[
          P=P_0\ \cup\ \bigcup_{i=1}^m\{p:p\mid P_i\},
\]
where \(P_i\mid A_i\) is its left label. The estimate also
allows additional free primes if each occurs on \(A_i C_i\)
in at least two rows of the star. Let \(G\) divide every
\(B_i\), and let each \(w_i\), \(1\le i\le m\), be a unit
modulo \(q\). Suppose a fraction at least \(n^{-\theta}\)
of all ordered distinct \(m\)-stars in the groups has
outside defect radical at most \(n^\beta\) for
\(G^{w_0}\prod_i C_i^{w_i}\). As in (MK5), assume
\(W/B_1\ge4m\), where \(W=n^d\) and \(B_1=n^b\).
Then
\[
 \theta\ \ge\
 m\left(d-b-\frac{s_B}{2}-\frac{\ell_A+\ell_C}{4}\right)
 +\frac{[m(\ell_C-s_C)-\beta]_+}{8}-o(1).
 \tag{MK7}
\]
If every prime of the free pool outside \(P_0\) occurs
on the targets in at least two rows, one has the stronger
inequality
\[
 \theta\ \ge\
 m\left(d-b-\frac{s_B+\ell_C}{2}\right)
 +\frac{[m(\ell_C-s_C)-\beta]_+}{4}-o(1).
 \tag{MK8}
\]
Here the errors have the same height-first convention as
(MK3). A cap on a union of records is sufficient whenever
the union includes the indicated target defect.

Fix a group and expose the rows of one star in order.
For step \(i\), use the prespecified pool \(P_0\) together
with the prime supports of \(A_j C_j\) for \(j<i\).
Let \(x_i\) be the base-\(n\) logarithmic radical on \(A_i C_i\)
captured by this pool, and set \(t_i=\log_n\operatorname{rad}(C_i)\).
The selected right class supplies a fixed prime pool containing
every nondivisible-order prime on \(B_i\). Indeed every such prime
lies in the support of its canonical \(e\)-free representative.
At a fixed budget, bin \(t_i,x_i\) with any fixed small mesh,
using upper endpoints for \(t_i\) and lower endpoints for \(x_i\).
Bin the two discounts on \(A_i,C_i\) separately if needed
to determine their individual outside caps. Their sum is \(x_i\)
by primitivity. The two pair choices with \(B_i\) in (MK1)
give a query exponent at most
\[
                 \frac{s_B}{2}+\frac{\ell_A+t_i-x_i}{4}
\]
up to the mesh and counting allowances, by averaging the two
pair bounds. All the pools at this step have polynomial prime
product and are fixed when counting the next row. In particular
the future left labels have no role in the query bound.

For a prime \(p\), let \(t\) be the number of rows in which
it occurs on \(A_i C_i\), and let \(v\le t\) be the number
of occurrences on \(C_i\) at valuations not divisible by \(q\).
Its multiplicity in the captured radicals is \(t\) if \(p\in P_0\),
and \((t-1)_+\) otherwise. If \(t\ge2\), then \(v\le2(t-1)\).
If \(p\notin P_0\) and \(v=t=1\), its only occurrence on
these letters is on one target, so no left label contains it.
The additional free-prime condition cannot remove it in this case.
It cannot divide \(G\), by primitivity at that target, so its
weighted valuation is nonzero modulo \(q\). It therefore occurs
in the paid defect. Thus for each prime, its \(v\) target
occurrences less its indicator in the paid defect are at most
twice its multiplicity in the captured radicals. Multiplying by
\(\log_n p\) and summing gives
\[
               \sum_i x_i\ \ge\
                        \frac{[\sum_i t_i-ms_C-\beta]_+}{2}.
\]
This counts repeated primes across rows and across the two
letters only through their actual multiplicities. The divisor
\(G\) may depend on the whole star, since its support is disjoint
from every target support.

Now multiply the query bounds in the fixed bins and sum the
bins. If \(S=\sum_i t_i\), the variable part of the resulting
exponent is at most \(S/4-[S-ms_C-\beta]_+/8\).
This function is increasing, so its maximum for
\(0\le S\le m\ell_C\) occurs at \(m\ell_C\).
There is no lower target radical hypothesis on the selected
rows. Summing over the \(B_1\) groups and comparing the bound
with \(\gg_m n^{-\theta}W^m/B_1^{m-1}\) proves (MK7).
For (MK8), use previous target supports and the captured
radical on \(C_i\) alone. The same multiplicity argument
now uses \(t\) equal to the number of target occurrences.
The query exponent is at most \((s_B+t_i-x_i)/2\).
The relevant function \(S/2-[S-ms_C-\beta]_+/4\)
is again increasing, proving (MK8).
[Pair population and the original graph estimate](ABC_MOVING_TARGET_POWER_CLASS_STARS.md#an-actual-all-row-graph-estimate)

At balanced saturation, write \(z=1/8\) for (MK7) and \(z=1/4\)
for (MK8). Both bounds imply
\[
          \theta+z\beta\ \ge\
             m\left(\frac16-b-\frac{s_B}{2}
                         +z\left(\frac13-s_C\right)\right)-o(1).
\]
The same pool condition allows (MK5). On fixing the other
rows, the query's own left-label primes are absent from its
target and all other label primes are fixed. An additional
free prime occurring on the query target and on an \(A_j C_j\)
of another row belongs to the fixed anchor pool. The
target-only condition in (MK8) has the same property.
Use (MK5), \(a\ge1/6-s_A/2-o(1)\), and \(\beta\ge0\) in
the exponent \(a+(1-1/m)b+(\beta+\theta)/m\) of (MS3).
Two resulting lower bounds are
\[
 \begin{split}
 &\frac16-\frac{s_A}{2}
       +(1-2/m)b+\frac{1/3-(s_B+s_C)/2}{m}-o(1),\\
 &\frac13-\frac{s_A}{2}-\frac{s_B}{2}
       +z(1/3-s_C)-\frac{b}{m}-o(1).
 \end{split}
\]
Their convex combination with weights \(1/(m-1)\) and
\((m-2)/(m-1)\) removes \(b\). It gives the lower bound
\[
 \frac13-\frac{s_A}{2}-\left(\frac12-\frac1{2m}\right)s_B
       -\frac{s_C}{2m}
       +g_m(1-3s_C)-o(1),
 \qquad
 g_m=\frac{(m-2)(2mz-1)}{6m(m-1)} .
\]
For vanishing aligned masses this is strictly greater than
\(1/3\) when \(m\ge5\) in (MK7) or \(m\ge3\) in (MK8).
In each of these ranges \(g_m\) increases to \(z/3\), since
\[
                g_{m+1}-g_m
                    =\frac{(2z+1)(m+1)-4}{6m(m-1)(m+1)}.
\]
These statements concern the bound obtained from (MS3)
even if its angular-spacing bracket is bounded. A pool chosen
by a more general function of the full left tuple requires
separate control of its support at primes confined to a
single target. In particular a prime of a pair determinant
can be absent from both of that pair's targets and occur
solely on a third target, so the full union of pair determinant
supports need not satisfy the stated additional-pool condition.
For example, take \(e=q=2\), \(A_i=P_i\), \(y_i=B_i=1\),
and \(A_i=29,26,31\) in order. The coefficient determinant
of the last two rows is \(-5\), and only the first target
among \(30,27,32\) is divisible by \(5\).

If the groups instead fix distinct complete right values, (FL1)
bounds each query by half the smaller of the outside caps on
\(A_i,C_i\). Thus (MK7) and (MK8) hold with \(s_B=0\)
regardless of the right valuations. In (MK5) the same argument
replaces \(\Psi(\ell_A,s_B,s_C)\) by
\(\min(\ell_A,s_C)/2\). At balance the complete-value dictionary
has \(b\ge1/6-o(1)\) by (FL1). On a common limiting valuation
profile, (FL2) also gives \(b\ge1/3-\min(c_A,c_C)-o(1)\),
where \(c_A,c_C\) are the nonatomic masses. To apply these
necessary bounds to a graph estimate, the completion of a
left tuple must be uniform over the complete right values.
Any further pool price needed for that completion enters the
graph exponent.
[Fixed-letter query and profile bounds (FL1)–(FL2)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
