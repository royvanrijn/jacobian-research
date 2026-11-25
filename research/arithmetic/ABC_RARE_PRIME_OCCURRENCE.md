# Rare primes and common divisor records at the small-budget abc endpoint

Let \(\mathcal F\) be any set of distinct ordered positive primitive triples
\(a+b=c\) with \(X/2<c\le X\) and \(\operatorname{rad}(abc)<c^\lambda\).
For \(0\le r<1/3\), call a prime rare relative to \(\mathcal F\) if it
divides \(abc\) in at most \(X^{r\lambda}\) rows of \(\mathcal F\).
All occurrences refer to distinct primitive triples.

## Uniform bounds for rare-prime sectors

Fix \(r<1/3\) and \(\eta>0\). There are \(\delta,\lambda_0>0\)
depending only on \(r,\eta\) such that, for every fixed
\(0<\lambda<\lambda_0\), the rows of \(\mathcal F\)
with
\[
 \sum_{\substack{p\mid abc\\p\text{ rare relative to }\mathcal F}}
       \log p\ \ge\ (1/2+\eta)\lambda\log X
 \tag{RP1}
\]
number \(O_{\lambda,r,\eta,\varepsilon}
(X^{(1/3-\delta)\lambda+\varepsilon})\) for every \(\varepsilon>0\).
The constant is uniform in \(\mathcal F\).

For \(0\le r<1/6\), the same bound, with
\(\delta,\lambda_0>0\) depending only on \(r\), holds for rows with at
least two distinct rare primes anywhere in \(abc\). Consequently,
for this range of \(r\) the same bound also holds when (RP1) is weakened
to a threshold of \((1/3+\eta)\lambda\log X\), with the constants then
allowed to depend on \(\eta\). If the graph on \(\mathcal F\) joining
two rows that share a prime has maximum degree at most
\(X^{r\lambda}-1\), every prime that occurs is rare and the bound
applies to all of \(\mathcal F\) for \(r<1/3\). The same conclusion with
a degree bound \(X^{r\lambda}\) follows by increasing \(r\) slightly.

There is also a saving for the coexistence of a rare prime and
macroscopic primes on a different letter. Fix \(0\le r<1/6\) and
\(\xi,\eta>0\). Impose that a prime rare relative to \(\mathcal F\)
occurs on a letter \(T\), while some other letter \(U\) satisfies
\[
       \sum_{\substack{p\mid U\\p\ge X^{\xi\lambda}}}\log p
                        \ \ge\ (r+\eta)\lambda\log X .
 \tag{RP6}
\]
The same bound holds, with \(\delta,\lambda_0>0\) depending on
\(r,\xi,\eta\).

### Joint recurrence of divisor records

At fixed truncation, take an anatomic box with its initial coefficients
fixed. Choose one root representation per row of a subfamily
\(\mathcal B\subseteq\mathcal F\) in this box. Write \(\mathcal I\)
for the positive-index coordinates across all three letters
and \(x_i(z)\) for the root of a row \(z\) at coordinate \(i\).
A divisor record at \(z\) is a tuple \(g=(g_i)_{i\in\mathcal I}\)
of positive integers with \(g_i\mid x_i(z)\). Put
\[
 G=\prod_{i\in\mathcal I}g_i,\qquad
 b(g)=\frac{\log G}{\lambda\log X},\qquad
 \mathcal B_g=\{z'\in\mathcal B:g_i\mid x_i(z')\text{ for all }i\}.
\]
For any fixed \(\kappa>0\), outside an \(O(X^{-\kappa\lambda/2})\)
fraction of \(\mathcal B\), every divisor record at the row satisfies
\[
                 |\mathcal B_g|\ \ge\
                   |\mathcal B| X^{-(b(g)+\kappa)\lambda}.
 \tag{RP2}
\]
Here and below the number of root coordinates and \(\lambda\) are
fixed at the height limit. The box convention is
\(\sum_i\log_X x_i\le 2\lambda+o_X(1)\).

At a fixed bound \(G\le X^{u\lambda}\) the number of records is at
most \(X^{u\lambda+o_X(1)}\). For a specified \(G\) the number of
ordered factorizations into \(|\mathcal I|\) factors is \(X^{o_X(1)}\),
uniformly up to polynomial height. Bin \(b(g)\) on a fixed grid of
mesh at most \(\kappa/4\). In a bin with upper endpoint \(u\), a
violation of (RP2) gives
\(|\mathcal B_g|<|\mathcal B| X^{-(u+3\kappa/4)\lambda}\).
Sum this bound over the records and bins. This estimates the number
of rows belonging to even one violating record.

In particular, the record \(g_i=x_i(z)\) and \(g_j=1\) for \(j\ne i\)
makes every prime in \(x_i(z)\) recur, with its prescribed power on
the same letter, throughout \(\mathcal B_g\). Taking all the exact
roots in a chosen packet as \(g_i\) gives the same conclusion for the
packet jointly. Exact equality of those roots can also be imposed in
(RP2) by the same proof using records of exact values. Recurrence of
a positive-index prime here does not require that its valuation in the
small initial coefficient vanish.

### Consequence for all balanced limiting profiles

Write \(L=\sum_i\alpha_i\), where \(X^{\alpha_i}\) is the upper
endpoint for \(x_i\) in the box. Use the limiting measures
\(\nu_T=\lim\sum_{i\text{ on }T}(\alpha_i/L)\delta_{e_iL}\)
and write \(c_T\) for the nonatomic mass of \(\nu_T\). Along any
sequence of limiting boxes at \(d=1/3\), with ambient families
containing the boxes, put \(a_r=1/3-r>0\).
On a typical row, the rare logarithmic radical mass on \(T\),
normalized by \(\lambda\log X\), satisfies
\[
 \limsup \frac{1}{\lambda\log X}
       \sum_{\substack{p\mid T\\p\text{ rare relative to }\mathcal F}}
           \log p
 \ \le\
 \sum_{\substack{t\ge0\\ \nu_T(\{t\})\ge a_r}}\nu_T(\{t\}).
 \tag{RP3}
\]
The inequality is in probability in the iterated endpoint limit.
More precisely, every fixed positive excess on the right occurs
on a vanishing fraction of the rows. If \(c_T>r\), a typical row
has no rare primes on \(T\), including primes with negligible
individual logarithmic mass.

To prove this, the box count is \(X^{D+o_X(1)}\) with
\(D/\lambda\to1/3\) and \(L/\lambda\to1\). A prime dividing an
initial coefficient occurs in every row of the box and is not rare.
For other primes, (RP2) implies that a rare prime on a typical row
can only occur in a coordinate with
\(\alpha_i/\lambda\ge D/\lambda-r-\kappa+o_X(1)\).
Indeed its absolute multiplicity in the ambient family is at least
the multiplicity guaranteed by the root record in the box.

Let the tolerance tend to zero in the outer limit. One may do
so slowly after fixing each truncation and then taking its height
limit. There are boundedly many candidate positions of normalized
cap at least \(a_r-o(1)\). Their scaled degrees are bounded because
\(\sum_i e_i L(\alpha_i/L)\) is bounded. Pass to a subsequence on
which their caps and scales converge. For each limit scale, the sum
of these limiting caps is at most the mass of the atom there, by
weak convergence in shrinking neighborhoods, or by first taking
fixed neighborhoods and then shrinking them. Each selected cap
has limit at least \(a_r\). The candidate caps bound the rare prime
mass on every row to which (RP2) applies, proving (RP3).
The same argument shows that a selected coordinate on \(T\) requires
\(1/3-c_T\ge a_r\), proving the absence statement when \(c_T>r\).
[Limiting degree measures and balance](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

For every fixed \(\theta>0\), a typical row has every rare prime
\(p\ge X^{(a_r-\theta)\lambda}\). To see this, use the frequency
\(\pi_{\mathcal B}(p)\) on the selected box family. If \(p\) is
rare relative to the ambient family then
\(\pi_{\mathcal B}(p)\le X^{r\lambda}/|\mathcal B|\).
A rare prime below the asserted size has
\[
   |\log(p\pi_{\mathcal B}(p))|
                  \ \ge\ (\theta-o(1))\lambda\log X .
\]
The mean over rows of the sum of these absolute deviations over their
prime factors is \(o(\lambda\log X)\) by (EB8). Markov's inequality
excludes even one such prime from a typical row.
[Two-sided frequency estimate (EB8)](ABC_EFFECTIVE_RADICAL_BUDGET.md#occurrence-frequencies-and-anatomic-entropy-at-saturation)

For \(r<1/6\), suppose a rare prime occurs on \(T\) on a nonvanishing
fraction of the rows. By (RP3), \(c_T\le r\). Both opposite letters
have positive nonatomic mass, since \(\Gamma\ge1/2\) and their
individual total masses are \(1/3\). Each can supply arbitrarily
cheap positive-weight full moduli. The pair inequalities with these
moduli imply
\[
          c_U,c_V\ \ge\ 1/3-c_T\ \ge\ 1/3-r,
             \qquad \Gamma\ge 2/3-r .
 \tag{RP7}
\]
On either opposite letter the mass from primes \(p\ge X^{\xi\lambda}\)
is at most \(r\) in the limit, for each fixed \(\xi>0\). Indeed their
coordinates have caps at least \(\xi\lambda\) and the same
finite-position argument as for (RP3) charges their total to the
atomic mass. Here the initial coefficient is too small to contain
any such prime at sufficiently late budgets.

Two distinct rare primes on the same letter also cannot occur on a
nonvanishing fraction. Choose \(\theta>0\) with \(a_r-\theta>1/6\).
Each such prime costs at least \(a_r-\theta\) in the normalized
radical budget by the frequency estimate, while the whole letter
has mass \(1/3+o(1)\). Thus there is at most one rare prime in a
typical row for \(r<1/6\), including within the same letter.
At \(r=0\), if the private prime occurs on \(T\) on a nonvanishing
fraction, its size on those typical rows is \(X^{(1/3+o(1))\lambda}\).
The remaining radical of \(T\) is \(X^{o(\lambda)}\), \(\nu_T\)
is one atom of mass \(1/3\), and the other two measures are fully
nonatomic by (RP7).

### Passage to a fixed exponent loss

If any asserted uniform sector bound fails, take budgets
\(\lambda_j\downarrow0\) for which the retained rows at arbitrarily large
heights have count exponent at least \((1/3-o(1))\lambda_j\).
Take a fixed truncation at each \(\lambda_j\) with coefficient-height
loss \(o(\lambda_j)\). Fix the coefficients and a dyadic root box.
These choices cost \(o(\lambda_j)\) in the exponent, and all the
retained rows inherit the original occurrence threshold. In fact
restriction can only decrease absolute multiplicities. Take inner
height limits for the counts and caps. The dictionary bound with raw
root caps gives \(D\le L/3+o(\lambda_j)\), while the radical bound
gives \(L\le(1+o(1))\lambda_j\). It follows that
\(D/\lambda_j\to1/3\), \(L/\lambda_j\to1\), and the selected boxes
give a limiting balanced profile.
[Dictionary bound (EB2)](ABC_EFFECTIVE_RADICAL_BUDGET.md#the-one-third-coefficient-for-coordinate-dictionaries)

At this profile \(\Gamma=\sum_T c_T\ge1/2\). By (RP3) the sum
of rare masses is at most \(\sum_T(1/3-c_T)=1-\Gamma\le1/2\)
on the normalized scale. This contradicts (RP1).
For rare primes on different letters, pass to a subsequence with the same
specified letters \(U,V\). The absence statement requires
\(c_U,c_V\le r\). The third letter \(T\) has
\(c_T\ge1/2-2r>0\). It is therefore positive-height and not
finite-scale, and it supplies arbitrarily cheap positive-weight
full moduli. The pair bound (GE1) then gives
\(c_U+c_V\ge1/3\), contrary to \(2r<1/3\).
For the smaller radical threshold, sparsity on at most one
letter means the normalized rare mass is at most \(1/3+o(1)\)
because every letter has normalized cap \(1/3+o(1)\) and its
fixed initial coefficient has cap \(o(1)\).
Two rare primes within one letter and the sector in (RP6)
are excluded by the preceding profile arguments.
[Balanced endpoint inequalities (GE), (GE1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13)

For the overlap-graph consequence, it remains to exclude a count
at the endpoint with total radical mass below the threshold in (RP1).
Fix any \(\eta<1/2\). The remaining rows have radical budget in base
\(X\) at most \((1/2+\eta)\lambda\) and the effective-budget estimate
gives a strict loss there, with prices \(w_X(p)=\log_Xp\).
[Effective-budget bound (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

## Common powerful moduli with controlled extra support

For fixed \(b_0,\rho>0\) with \(2b_0<\rho<1-b_0\), consider a
saturation sequence of boxes as above. For each fixed \(\sigma>0\),
at sufficiently late fixed budgets and then at large heights,
a typical row \(z\) has the following property simultaneously for
every squarefree divisor record \(g\) at \(z\) with \(b(g)\le b_0\).
The fibre \(\mathcal B_g\) contains a family \(\mathcal H_g\) of at
least
\[
 X^{((\rho-2b_0)/3-\sigma)\lambda}
 \tag{RP4}
\]
rows, including \(z\), in which distinct rows \(z_1,z_2\) satisfy
\[
 \log\operatorname{rad}\gcd(a_1b_1c_1,a_2b_2c_2)-\log G
                      <\rho\lambda\log X .
 \tag{RP5}
\]
Here a squarefree record means that its product \(G\) is squarefree.

Indeed (RP2) gives a fibre of size
\(X^{(1/3-b(g)-o(1))\lambda}\) with arbitrarily small fixed tolerance.
For any prescribed \(y\) of the fibre, assign zero price to the
prime divisors of its product \(abc\) and use usual prices on the
other primes. The partition hypothesis (EB3) holds uniformly:
the prescribed product is at most \(X^3\), so its support has
\(O(\log X/\log\log X)\) primes. A row failing (RP5) with \(y\)
has priced budget at most \((1-b(g)-\rho)\lambda\), because it
already shares every prime of \(G\) with \(y\). By (EB4) the count
of such rows is at most \(X^{((1-b(g)-\rho)/3+\xi)\lambda+\varepsilon}\)
with an arbitrarily small fixed \(\xi>0\) at sufficiently small fixed
budgets, uniformly in \(y\) and \(g\). Starting with \(z\), select
one row at a time and discard it and the rows that fail (RP5) with
it. The lower exponent from the fibre size divided by the maximum
discard size is \((\rho-2b(g))/3\) up to arbitrarily small losses.
This proves (RP4). It uses the radical cap for each complete triple,
including all three letters.

A shorter sequence can control overlap with the entire previous
union. For any fixed \(0<\tau<\rho/(2b_0)-1\), one can choose
\(\lfloor(\log X)^\tau\rfloor\) ordered rows of the fibre, starting
with \(z\), such that each row shares outside \(G\) radical weight
less than \(\rho\lambda\log X\) with the union of all preceding
supports. Choose \(a\in[0,1)\) so that
\[
             2b_0<(1-a)\rho,\qquad a>1-\frac1{1+\tau}.
\]
At a given step the previous union has
\(O((\log X)^{1+\tau}/\log\log X)\) primes. Assign zero price to
the primes of \(G\), price \(a\log_X p\) to the remaining primes
in that union, and the ordinary price to all others. The polynomial
logarithmic size bound on the union gives (EB3) by summing its
\(j\)-th possible contribution at most \((j+1)^{-a}\).
The separately zero-priced support of \(G\) has \(o(\log X)\)
primes. A bad row has effective budget at most
\((1-b(g)-(1-a)\rho)\lambda\). By (EB4) its count is power smaller
than the whole fibre, uniformly at every step. The already chosen
rows also have smaller count. Thus the selection can continue.
In a bounded scaled-degree band of squarefree coordinates,
stripping the previous supports and \(G\) gives pairwise coprime
powered products. Outside \(G\) this removes at most
\(\rho\lambda\log X\) of each row's radical weight. For an upper
scaled-degree bound \(v\) the associated numerical-height loss
is at most \(v\rho+o(1)\) in base \(X\).

In a balanced limiting profile, fix a positive compact range
\([u,v]\), \(u>0\), on a letter \(T\) with positive nonatomic
mass there. Given any small bound on individual coordinate caps,
delete neighborhoods of the finitely many atoms above a smaller
bound. Along the sequence, each remaining cap on the range is below
the original bound. Otherwise an individual coordinate would
produce another atom of at least that mass in the range.
The deleted neighborhoods can be chosen to remove arbitrarily
little nonatomic mass. The one packed coordinate per letter in the
anatomic decomposition can also be discarded: its raw cap is
\(o(\lambda)\) at the chosen truncations. The remaining positive-index
coordinates are squarefree. Thus the range still supplies a packet of
roots of total normalized cap between \(b_0/2\) and \(b_0\) for all
sufficiently small fixed \(b_0>0\), by choosing the individual cap
bound below \(b_0/4\) and collecting roots up to the first crossing
of \(b_0/2\). The dyadic lower endpoints give the same cost and
numerical weight in the inner limit. The exact packet record on a
typical row can be used as \(g\), and every row in (RP4) has the
same divisor \(\prod_i g_i^{e_i}\) on \(T\), of numerical height
at least \(X^{u b_0/2-o(1)}\) and with fullness of order at least
\(1/L\). Each selected position has scaled degree in \([u,v]\).
This construction also gives joint recurrence for any small packet
outside the large-atom neighborhoods, by (RP2), without choosing
a particular range or applying (RP4).
[Anatomic decomposition, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

The threshold \(2b_0\) in this intersection argument has a
support-only obstruction. Fix \(0<b<1/3\) and
\(s=1/3-b\). For each letter take two lists of distinct primes,
with all six lists disjoint. The first lists have common length
\(X^{b\lambda+o_X(1)}\) and primes in
\([X^{b\lambda}/4,X^{b\lambda}/2]\). The second lists have common
length \(X^{s\lambda+o_X(1)}\) and primes in the analogous interval
with \(s\) in place of \(b\). Such lists follow from the prime
number theorem, including when \(b=s\).
An abstract row is a pair of indices \((i,j)\), and each letter
uses its \(i\)-th prime of the first list and its \(j\)-th prime
of the second list. The total radical is at most \(X^\lambda\),
and there are \(X^{\lambda/3+o_X(1)}\) rows with near-reciprocal
individual prime frequencies on the logarithmic scale. Fixing an
\(i\)-prime on one letter fixes the \(i\)-primes on both other
letters. Every pair of distinct rows in that fibre therefore
shares outside radical of weight \((2b+o_X(1))\lambda\log X\).
This construction concerns prime supports and does not impose
the additive equation.

## Minimally repaired quotients with rare primes

Fix \(\ell>0\), \(0\le r<1/6\), a sufficiently small fixed
\(\lambda>0\), and a fixed positive integer repair order \(M\).
Let \(\mathcal Q\) be any family of distinct primitive quotients
\(a+b=c\) in the sector
\[
  \operatorname{rad}(abc)<c^\lambda,\qquad
  \operatorname{rad}(c)\ge c^{\ell\lambda},\qquad
  H_M(c)=\prod_{p\mid c}p^{\max(v_p(c),M)}\le Y .
\]
Here \(H_M\) is the target height after the least repair with orders
\((1,1,M)\). Put \(t=\min\{1,1/(M\ell\lambda)\}\). Count rows in
\(\mathcal Q\) in which at least two distinct primes each occur in
at most \(Y^{rt\lambda}\) members of \(\mathcal Q\). There are
\(\delta,\lambda_0>0\), depending only on \(r\), such that for
\(0<\lambda<\lambda_0\) this count is
\[
       O_{\lambda,M,\ell,r,\varepsilon}
                     \bigl(Y^{(1/3-\delta)t\lambda+\varepsilon}\bigr)
 \tag{RP8}
\]
for every \(\varepsilon>0\), uniformly in \(\mathcal Q\).
Occurrence may equivalently be measured on the minimally repaired
representatives, because the least repair multiplier is supported
on primes of the quotient.
[Target minimal repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

Indeed \(H_M(c)\ge c\) and \(H_M(c)\ge\operatorname{rad}(c)^M
\ge c^{M\ell\lambda}\), so every quotient has \(c\le Y^t\).
Choose \(0<\kappa<1\) with \(r/(1-\kappa)<1/6\). The quotients
with \(c\le Y^{(1-\kappa)t}\) satisfy (RP8) by the unrestricted
one-third estimate with an allowance sufficiently small in terms
of \(\kappa\). For the others use original height shells \(X/2<c\le X\)
with \(Y^{(1-\kappa)t}<X\le 2Y^t\).
In each shell the absolute occurrence threshold is at most
\(X^{(r/(1-\kappa))\lambda}\). Apply the two-rare-prime theorem
and sum over the shells. The \(O(\log Y)\) factor is absorbed into
the height allowance.

In the range \(M\ell\lambda\ge1\), the exponent in (RP8) is
\((1/3-\delta)/(M\ell)\) before the allowance. The argument counts
primitive quotients with a single least repair each, so it incurs
no contribution from arbitrary common multipliers.
[Unrestricted one-third estimate](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

## Remaining arithmetic problem

These results impose occurrence constraints across the entire
balanced endpoint. They leave open the case in which the
nonatomic mass is carried by many individually recurring roots
and each recurring packet fixes only a small part of the radical.
An additional arithmetic inequality would need to bound stars
like (RP4) using the primitive equation and the numerical weight
of their common powered divisor.
