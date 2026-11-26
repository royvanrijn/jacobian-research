# Primitive secants in near-full abc families

Nearly every secant of a near-full balanced family gives a primitive
additive triple of positive numerical height and almost the full
original radical budget. This applies across the balanced valuation
profiles. It supplies an exponent saving whenever a sufficiently large
fraction of the secants compresses the radical by a fixed factor in
the exponent. The argument uses two-direction recovery and a count of
rich lines in the actual plane of complete letters.

## A uniform bound for queries short in two directions

Let \(n=X^L\), with \(L>0\) fixed before \(X\to\infty\). In the
shell \(X/2<c\le X\), take positive primitive triples \(a+b=c\)
with \(\operatorname{rad}(abc)\le c^L\) and
\(\operatorname{rad}(W)\le n^{1/3+\kappa}\) for \(W=a,b,c\).
Identify a triple with \(x=(a,b)\). For distinct integer points
\(x,y\) in this plane put
\[
 u(x,y)=\frac{x-y}{\gcd(|x_1-y_1|,|x_2-y_2|)},\qquad
 H(x,y)=\max(|u_1|,|u_2|,|u_1+u_2|).
\]
The height depends only on the unoriented direction. If none of the
three entries vanishes, their absolute values, with the largest last,
give a positive primitive additive triple. Its radical is denoted
\(R(x,y)=\operatorname{rad}(|u_1u_2(u_1+u_2)|)\).

There exist absolute \(\theta,\kappa_0,\Delta,L_0>0\) such that,
for \(0<L<L_0\), \(0\le\kappa\le\kappa_0\) and every
\(\varepsilon>0\), uniformly in distinct anchors \(y,z\in[0,X]^2
\cap\mathbb Z^2\),
\[
 \#\{x\ne y,z:H(x,y),H(x,z)\le X^\theta\}
       \ll_{L,\varepsilon} X^{(1/3-\Delta)L+\varepsilon}.
 \tag{SC1}
\]
Only the queries are required to have the original three caps.
When an imported estimate uses a strict radical inequality,
increase the shell budget by an arbitrarily small fixed multiple
of \(L\) before the height limit.

First, a line through \(y,z\) has at most
\(n^{2/9+O(\kappa)+\zeta}\) such queries for arbitrary fixed
\(\zeta>0\) at sufficiently small fixed \(L\). A noncoordinate
line with nonzero intercept uses (AS1), a fixed-letter line uses
(FL1), and a line through the origin has at most one primitive
positive query. The coefficients of a line through two integer
points of height at most \(X\) are of polynomial height. The same
bound applies to every line meeting the original shell in at least
two queries, using those queries as the two points.
[Affine-slice count (AS1)](ABC_AFFINE_SLICE_RADICAL_COUNTS.md#a-fixed-weighted-shift-inside-the-three-letter-count)
· [Fixed-letter count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Off this line put \(u=u(x,y)\), \(v=u(x,z)\),
\(d=[u,v]=u_1v_2-u_2v_1\ne0\). The intersection identity is
\[
                  dx=dy+[z-y,v]u.
\]
Thus each of \(da,db,dc\) is a fixed homogeneous linear function
of the integer vector
\[
                  Z=(d,u_1v_1,u_1v_2,u_2v_1,u_2v_2).
\]
Its entries have absolute value at most \(2X^{2\theta}\). The
coefficients of the linear functions may depend on the anchors.
There is no modular division by \(d\).

Here is the rank argument with the uniformity in \(\theta\).
If (SC1) fails for every choice of the four constants, take
\(\theta_j,\kappa_j,\Delta_j,L_j\to0\) and arbitrarily large
violating heights at each fixed \(L_j\). When extracting exponents,
use an inner exponent allowance smaller than the prospective gap
at that fixed budget. Refine by the original anatomic boxes and
take the weak degree limits after the inner height refinements.
The unconditional one-third limiting bound makes the limiting count
exponent \(1/3\). Individual masses tend to \(1/3\): they are
bounded above by \(1/3+o(1)\), and a total mass bounded below one
by a fixed gap would have a smaller limiting count. By (GE1) at
least one letter \(W\) has positive nonatomic mass in some compact
band of positive scaled degrees. Choose this role along a
subsequence. All refinements and truncations are made before the
inner height limit. [One-third count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [Nonatomic mass (GE1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

Fix \(\xi>0\) small compared with the nonatomic mass in that band.
At any stage with a nonempty query class, select a feature witness
\(Z_0\) from it, with its nonzero integer \(d_0W_0\).
The witness prime pool has polynomial height. Discard queries
with normalized radical mass at least \(\xi\) on this pool in \(W\).
By (EB5) their count exponent is at most \(1/3-\xi/4\) for
sufficiently late budgets, allowing for the inner height losses.
This bound is uniform in the witness. On the rest, select a full
divisor \(F\mid W\) coprime to \(d_0W_0\) with normalized
radical cost at most \(p\) and \(F\ge X^{c_*(p)}\) for
\(c_*(p)>0\). Here \(p>0\) is fixed sufficiently small compared
with \(\xi/5\). The selection uses full powers in the compact
band after removing small neighborhoods of its finitely many
atoms of mass at least \(p/4\). These neighborhoods can exclude
arbitrarily little nonatomic mass. The available individual powers
then have cost at most \(p/2\) for sufficiently late boxes.
A greedy packet has cost in \([p/2,p]\) and the band's lower
degree supplies \(c_*(p)>0\). This selection depends only on the
original profiles, which remain valid on the retained rows.
There are at most \(n^{p+o_X(1)}\) different locks by radical
enumeration at height \(X\). Fix one, losing at most that factor.
[Full-divisor selection from a compact band](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#limiting-setup-and-cost-of-fixed-prime-supports)
· [Pool discount (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

If the new feature rows still span the old linear space of rank
\(r\le5\), take \(r\) spanning retained rows and a nonzero
minor \(\Delta_0\) of \(r\) coordinate columns. Integer Cramer's
rule expresses \(\Delta_0 Z_0\) as an integer combination of those
rows. Apply the linear function giving \(dW\).
As \(F\mid W\) on them, \(F\mid\Delta_0 d_0W_0\), hence
\(F\mid\Delta_0\). But
\(\lvert\Delta_0\rvert\le r!2^r X^{2r\theta_j}<F\)
for late \(j\) and large \(X\).
The rank must therefore drop every time a lock is fixed. At most
five locks have total cost at most \(5p\), smaller than
\(\xi/4\). Thus the discard and lock steps continue to leave
a positive exponent class at every rank descent. Rank zero
is impossible because \(d\ne0\). This contradicts the choice
of the sequence and proves (SC1). The linear calculation is
related to the earlier rational feature argument inside fixed
divisor records. Here it is applied before any affordable record
is chosen and needs positive nonatomic mass on only one of the
three original letters.
[Rational feature rank descent in a fixed divisor record](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#lifting-short-rational-features-from-two-congruences)

As a separate consequence, fix a rank bound \(r\) and an exponent
\(K\). Suppose all three complete letters, multiplied by the same
nonzero integer \(d(x)\) of magnitude at most \(X^K\), are
linear functions of assigned integer features in \([-X^\theta,
X^\theta]^r\). The same argument gives a strict saving below
\(1/3\) for sufficiently small fixed \(\theta,\kappa,L\) depending
on \(r,K\). The witness pool is still of polynomial height,
the number of drops is at most \(r\), and the bound on a witness
minor is \(r!X^{r\theta}\). Specifying two such linear numerators
suffices, because \(a+b=c\) gives the third.

## A pair inequality on the full balanced family

Let \(\mathcal E\) be any family in the original shell, with
\(E=|\mathcal E|\). Keep the constants of (SC1), shrinking
\(\kappa_0,L_0\) if necessary. Fix \(0<b\le B\) and a
polynomial-height prime pool \(\mathcal P=\{p:p\mid P\}\),
where \(1\le P\le X^K\) and \(B,K\) are fixed. A pair of distinct
queries is counted by \(I_{b,B,\mathcal P}\) when all three
direction entries are nonzero and
\[
 R(x,y)\le n^B,\qquad
 \prod_{\substack{p\mid R(x,y)\\ p\notin\mathcal P}}p\le n^b.
\]
Allow arbitrary fixed \(\zeta>0\) in the normalized exponent.
At all sufficiently small fixed \(L>0\), the following bound
holds uniformly in the family and in the pool:
\[
 I_{b,B,\mathcal P}
 \ll_{B,K,b,L,\zeta}
 E^{3/2} n^{\max(1/3-\Delta,b/3)/2+\zeta}
       + E n^{\,2/9+2\kappa/3+\zeta}.
 \tag{SC2}
\]
For \(P=1\) one can set \(B=b\). Counts here are of ordered pairs.
Pairs with a zero direction entry can also be included in (SC2),
as can all short pairs without any condition on their radicals.

For the proof let \(S\) be the number of ordered short pairs
\(x\ne y\) with \(H(x,y)\le X^\theta\). Write \(s(x)\) for the
number of their anchors at query \(x\). By (SC1), for any prescribed
fixed \(\epsilon_0>0\),
\[
 S^2/E\le \sum_x s(x)^2
       \le S+O_{L,\epsilon_0}\big(E^2 n^{1/3-\Delta+\epsilon_0}\big).
\]
Solving gives \(S\ll E+E^{3/2}n^{(1/3-\Delta+\epsilon_0)/2}\).
This step has no condition on the radicals of the short directions.

For the tall directions, let \(M\) be the number of unoriented
primitive directions satisfying the two radical restrictions.
Each derived positive primitive triple gives at most a bounded
number of such directions. The direction heights lie in
\([X^\theta,X]\): all original letters are at most \(X\), so their
differences before normalization are also at most \(X\).
Group this range in dyadic height shells with upper endpoint \(Y\).
In each shell use (EB4), pricing the primes of \(P\) at zero and
the rest ordinarily. The pool has \(O_{K,\theta}(\log Y/\log\log Y)\)
primes uniformly. The raw relative budget at shell height is
at most \(B L\log X/\log(Y/2)=O_{B,\theta}(L)\). Round this
relative budget up on a fixed mesh of spacing \(\epsilon_1 L\)
at the fixed original budget \(L\). The mesh has finitely many
possible budgets with bound \(O_{B,\theta}(L)\) independently
of \(X\). The effective cap in base \(Y\) is \(b L\log X/\log Y\).
The count from (EB4) is \(n^{b/3}\) up to arbitrarily small
fixed normalized allowances, because this effective exponent
times \(\log Y/\log X\) is \(bL\). Thus
\(M\ll_{B,K,b,L,\zeta} n^{b/3+\zeta}\), using part of the
allowance and summing the height shells. The same argument for
\(P=1\) uses the unrestricted one-third limiting bound.
[Effective radical count (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

Let \(J\) be the largest number of input points on one line.
From the first paragraph of the proof of (SC1), one has
\(J\ll n^{2/9+2\kappa/3+\zeta}\). By the rich-line consequence
of the Szemerédi–Trotter theorem, the number of lines containing
at least \(k\ge2\) points of \(\mathcal E\) is
\(O(E^2/k^3+E/k)\). Summing the number of pairs on lines of
dyadic richness \(k,2k,\ldots,J\) gives
\[
                       O(E^2/k+EJ).
\]
On lines of richness below \(k\), every point has at most one
line through it in each permissible direction, so at most
\(EMk\) ordered pairs use these lines. If \(1\le M\le E\), take
\(k=\max(2,\lceil\sqrt{E/M}\rceil)\), giving
\(O(E^{3/2}M^{1/2}+EJ)\). If \(M>E\) the first term already
dominates the trivial pair bound. If \(M=0\) this branch is empty.
Apply the count on \(M\) and add the short pairs, splitting the
exponent allowances as needed. A primitive direction with a
zero entry has height one, so was already included on the
short branch. This proves (SC2).
[Szemerédi and Trotter, Extremal problems in discrete geometry](https://trotter.math.gatech.edu/papers/38.pdf)

## Consequences at the one-third scale

Fix \(0<\eta<1\) and \(B,K>0\). For some \(\sigma,\kappa_1,L_1>0\)
depending only on these parameters, if \(0<L<L_1\),
\(0\le\kappa\le\kappa_1\) and \(E\ge n^{1/3-\kappa_1}\), at most
\(O_{L,\eta,B,K}(E^2 n^{-\sigma})\) ordered pairs are short or have
\(R(x,y)\le n^B\) with outside-\(\mathcal P\) radical at most
\(n^{1-\eta}\). This follows by setting \(b=1-\eta\) in
(SC2), or \(b=\min(B,1-\eta)\) when \(B<1-\eta\). Choose
\(\kappa_1,\zeta\) small compared with \(\Delta,\eta/3,1/9\).
In particular, taking \(P=1\) shows that all but a power-small
fraction of the secants have both positive numerical height
and radical larger than \(n^{1-\eta}\), without a raw secant
cap hypothesis.

Conversely, if at least \(E^2 n^{-\rho}\) ordered pairs are short
or pass the radical test in (SC2), then
\[
 E\ll n^{\,\max\{1/3-\Delta+2\rho,\ b/3+2\rho,\
                               2/9+2\kappa/3+\rho\}+2\zeta}.
 \tag{SC3}
\]
For a fixed \(b<1\), sufficiently small fixed \(\rho,\kappa,\zeta\)
give a strict saving in the original row count. This counts
arbitrary pairs in the shell. It does not force a power proportion
of such pairs. The explicit missing estimate on this route is
a lower bound of this size for the number of short or low-radical
primitive differences in every near-full family. All secant heights,
coefficients and caps are on the primitive normalized vectors.
Representing the input by least repairs gives the very same
test on their primitive quotients.
[Least-repair equivalence](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md)

## Sequentially fresh radical outside earlier secants

The fixed-pool estimate can be used with a pool chosen from
one endpoint and an arbitrary bounded list of previous endpoints.
Fix a positive integer \(k\), \(B,K>0\), and \(0<\eta<1\).
Allow a prescribed pool \(\mathcal P_0\) supported on an integer
at most \(X^K\), which may depend on the central query \(x\).
At the near-full scale following (SC2), all but
\(O(E^{k+1} n^{-\sigma})\) ordered tuples
\((x,y_1,\ldots,y_k)\) have all directions tall with nonzero
entries. For each \(i\) with \(R(x,y_i)\le n^B\), they satisfy
\[
 \log_n\prod_{\substack{p\mid R(x,y_i)\\
                  p\notin \mathcal P_0\ \cup\ \operatorname{supp}(abc(x))
                         \ \cup\ \bigcup_{j<i}\operatorname{supp}(R(x,y_j))}}p
                    >1-\eta .
 \tag{SC4}
\]
Here \(\operatorname{supp}\) denotes prime support. In particular, on
the branch with all \(k\) raw caps, their joint fresh radical outside
the central query and \(\mathcal P_0\) is greater than
\(n^{k(1-\eta)}\). The raw joint radical passes the same lower bound
on all but a power-small fraction of tuples, without any raw cap
hypothesis, by using \(B=k\). If the joint raw bound fails, all
individual raw caps then hold.

To prove this, choose a sufficiently small \(\omega>0\).
The short-pair bound and the rich-line estimate at
\(k_{\rm line}=\lceil n^\omega\rceil\) first discard a
power-small fraction of tuples without using any prime pools.
In the remaining pairs, for fixed \(x\) a specified direction
leads to fewer than \(k_{\rm line}\) different anchors.
Condition on \(x,y_1,\ldots,y_{i-1}\) and put all the displayed
excluded primes into the pool. This is supported on an integer
of height \(X^{O_{K,k}(1)}\), independently of the next anchor.
The proof of the slope count for (SC2) then gives at most
\(n^{(1-\eta)/3+\zeta}\) candidate tall slopes with the stated
raw cap and with outside radical at most \(n^{1-\eta}\).
Consequently there are at most
\(O(n^{(1-\eta)/3+\zeta+\omega})\) next anchors that pass
the cheap test on geometrically good pairs. At the near-full
input scale, choosing the losses small compared with \(\eta\)
makes this a power-smaller choice. Sum over the other endpoints
and the \(k\) stages. If \(B<1-\eta\), use \(b=B\) in the slope
bound, which only improves the count. Arbitrary later
endpoints play no role in the pool. This order of conditioning is
essential.

For example, if a near-full pair relation has raw direction
radical at most \(n^{1+o(1)}\), then its directions form a
second family of exponent \(1/3-o(1)\) at the same absolute
radical scale. To see this, discard short and rich pairs and
use the rich-line bound for any smaller direction dictionary.
For an upper bound of \(n^{1/3-\xi}\) directions the number
of possible ordered pairs is power smaller than \(n^{2/3}\)
by the line count used in (SC2). The shell adjustment in the
proof gives positive relative budget tending to zero at the
derived heights. Input rank and arbitrary derived valuation
profiles remain allowed. By the conditioning argument, a
near-threshold direction is nearly radical-disjoint from each
of its input rows outside a power-small exceptional pair set.
For one input use (SC4) with \(k=1\), and then interchange
the inputs. At the threshold, any already established output
sector with a fixed exponent gap contributes only a power-small
pair class, by applying the same direction-dictionary bound.
This includes the existing sufficiently high-powerful output
sector when its hypotheses hold at the derived height.
[High-powerful triple transfer](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md)

## The radical of a physical coordinate difference

There is a separate constraint on the complete dictionaries
when one retains the actual coordinate difference before
normalizing the direction. Work with two input subfamilies
\(\mathcal A,\mathcal B\) of the same balanced shell, and fix
a role \(W\). Write \(F_{\mathcal A},F_{\mathcal B}\) for the
largest sizes of their respective fixed-\(W\) fibres.
For a fixed bound \(B>0\) and any \(\zeta>0\), there exist
\(\sigma,\kappa_2,L_2>0\) such that if \(0<L<L_2\),
\(0\le\kappa\le\kappa_2\) and \(0\le b\le B\), then
\[
 \begin{split}
 &\#\{(x,y)\in\mathcal A\times\mathcal B:
            W(x)\ne W(y),\quad
            \operatorname{rad}(|W(x)-W(y)|)\le n^b\}\\
 &\qquad\ll_{B,\zeta,L}
 (|\mathcal A|+|\mathcal B|)n^{1/3-\sigma}
            +F_{\mathcal A}F_{\mathcal B}\,
          n^{\min\{(2/3+2\kappa+b)/3,\ 1/3+\kappa\}+\zeta}.
 \end{split}
 \tag{SC5}
\]
This bound is uniform in \(b\) in the displayed range.

Use the positive-scale consequence of (SH1). For some fixed
\(0<r_0<R_0\), \(\gamma>0\), the shell where any one letter
has less than \(\gamma\) normalized radical in exact valuations
\(v_p(W)L\in[r_0,R_0]\) has a uniform exponent gap from \(1/3\).
To see this finite-parameter version, choose a slightly wider
band and a smaller \(\gamma\) than the compact limiting bound
in (SH1). A violating sequence would have a balanced limiting
profile with the smaller mass on the widened band. Test against
a continuous function supported in the widened band and equal
to one on the compact limiting band. The anatomic truncation
errors can be made arbitrarily small before taking the weak
limit. This contradicts the limiting bound. Discard pairs with
an input in that sector.
[Uniform positive scaled-degree mass (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)

Fix any first input \(x\). With \(\xi>0\) smaller than both
\(\gamma/2\) and \(\zeta/4\), (EB5) discards at a uniform
power saving the second inputs \(y\) with at least \(\xi\)
normalized radical mass of \(W(y)\) on the primes of \(W(x)\).
Together the discarded pairs are covered by the first term of
(SC5), after reducing \(\sigma\) if needed. On the remaining
pairs put \(g=(W(x),W(y))\). Its radical is at most \(n^\xi\).
Also \(W(y)/g\ge X^{r_0(\gamma-\xi)}\), since the compact-band
prime powers outside \(W(x)\) survive fully. With \(s=\min(W(x),W(y))\)
and \(t=\max(W(x),W(y))\), the triple
\[
                 (s/g,\ (t-s)/g,\ t/g)
\]
is positive primitive and has radical at most \(n^{2/3+2\kappa+b}\).
It has height between \(X^{r_0(\gamma-\xi)}\) and \(X\).
The height and budget mesh used for (SC2), with no pool
discount, bounds the number of output triples by
\(n^{(2/3+2\kappa+b)/3}\) up to an arbitrarily small fixed
normalized loss, uniformly in \(b\). The two pure quotients have radical
caps \(n^{1/3+\kappa}\) each. Apply the two-letter bound (GD3)
with singleton root groups to the output boxes. The bounded
raw cap on the third quotient keeps the relative input budget
at \(O_B(L)\) even when its cost is omitted. This yields the
second bound \(n^{1/3+\kappa}\) up to the same type of loss.
There are \(n^{\xi+o_X(1)}\) choices for \(g\). The triple, \(g\) and
orientation recover both complete letters exactly. Each
has at most the stated number of input preimages. This proves
(SC5). The third original radical constraint enters the overlap
discount, the positive-scale bound and the fixed-letter caps
used below. [Two-letter count (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)

### Complete-letter dictionaries forced by bounded differences

Consider a height-first endpoint sequence. Suppose a selected
relation of pairs of exponent \(2/3-o(1)\) has bounded physical
difference radical on role \(W\), meaning that
\(\operatorname{rad}(|W(x)-W(y)|)\le n^B\) for a fixed \(B\).
The equal-letter pairs are negligible, since (FL1) gives at
most \(n^{1/6+o(1)}\) second inputs for each first input.
Split both sides into fibre-degree bins in \(W\), using the
sizes of the fibres in the original two families. Also bin the
radical of the difference in normalized bins of sufficiently
small fixed width. Arbitrarily fine meshes can be taken
diagonally after the inner height limit. Keep an exponent
\(2/3-o(1)\) class. Let \(f,f'\le1/6\) be the two limiting
fibre exponents and \(b\le B\) the limiting normalized radical
of the difference. Equation (SC5) gives
\[
         \frac23\le f+f'+
                   \min\left\{\frac{2/3+b}{3},\frac13\right\}.
 \tag{SC6}
\]
In particular \(f=f'=1/6\) and \(b\ge1/3\).
Dividing the ambient family sizes by the lower fibre sizes in
these bins bounds both resulting \(W\) dictionaries by
\(n^{1/6+o(1)}\). The two projections to individual input
rows retain exponent \(1/3-o(1)\), since either row has at
most \(n^{1/3+o(1)}\) partners. This uses the two pure-letter
output caps to handle every fixed \(B\), including when the
original difference radical substantially exceeds \(n^{1/3}\).

Apply the same refinement on all three roles. If all three
physical difference radicals are bounded by \(n^B\) on a pair
relation of exponent \(2/3-o(1)\), there are two projected
families at the balanced endpoint in which all three complete
label dictionaries have the minimum exponent \(1/6\).
This conclusion needs no minimum nonatomic-mass hypothesis
on either input. In particular it applies if the primitive
secant radical and the radical of its common scaling
gcd each have a fixed bound \(n^{O(1)}\).
If the product of the three individual physical difference
radicals is at most \(n^{1+o(1)}\), each has exponent
\(1/3+o(1)\) on the retained relation, by (SC6).
The same conclusion follows when the primitive secant
radical is at most \(n^{1+o(1)}\) and that of the scaling
gcd is \(n^{o(1)}\): their product with the latter raised
to the third power bounds the product of the three physical
radicals.

In the central least-repair normalization \(Z=Y^{1/M}\) and
\(Z=n^{1/3+o(1)}\), consider a class whose occurring target
radicals each have multiplicity at most \(Z^{u+o(1)}\) for
fixed \(u<1/2\). The same cap holds for each complete target.
If the class attains the critical row exponent, for every
fixed \(B\) its pairs of distinct complete targets with
\(\operatorname{rad}(|c-c'|)\le n^B\) are power smaller
than the full pair count. Use (SC5) with
\(F_{\mathcal A},F_{\mathcal B}\le n^{u/3+o(1)}\)
and the second limb of its minimum. This tests pairs of
distinct targets across the actual low- or
intermediate-multiplicity part of the repaired family.
It gives no upper bound on multiplicity or on the number
of occurring targets by itself.
[Repaired target-fibre and multiplicity range](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)
· [Extension to adaptive pools without a raw difference cap](ABC_FIXED_ENDPOINT_SHIFT_RADICAL_COUNTS.md#dense-pair-relations-and-target-fibres)

The later moving-curve and short-cofactor estimates remove
the bounded-radical hypothesis on the common scaling gcd.
They put its prime factors on a pool chosen from the first
row up to a numerically short part, then use a two-letter
count which accepts that unrecorded factor.
[Primitive coordinate transfer without a scaling gcd cap](ABC_UNRECORDED_COFACTOR_SECANT_TRANSFER.md#primitive-secant-coordinates-with-an-arbitrary-scaling-gcd)

## Scope of the reduction

The short-direction proof uses the simultaneous input
constraints through the positive nonatomic band and the
prime-price discount. The direction count then uses the
numerical height of the recovered primitive triple. The
coordinate inequality uses instead a height lower bound
for the normalized pair of complete values. These are
different recoveries. The earlier product secants use
projective cross-products, while the present direction is
an affine difference of complete triples.
[Earlier projective secant transfer](ABC_SECANT_ENDPOINT_TRANSFER.md)
· [Earlier rich common-target shifts](ABC_POPULAR_SHIFT_RADICAL_AMPLIFICATION.md)

Neither the scalar normalization nor the height-first
estimates supply an upper bound on the radical of the
scaling gcd. They also do not force the abundance of
compressed secants needed in (SC3). The residual arithmetic
question on this route is whether the original three caps
force that abundance, or otherwise force a bounded
physical-difference relation on a role whose fixed-letter
multiplicity is below the half-radical exponent. The latter
option covers the intermediate multiplicities and leaves the
minimum-dictionary case to a separate argument. Any use of
only low primitive secant radicals must account explicitly
for the scaling gcd before invoking the complete-letter
conclusion (SC6). All these counts use primitive quotients.

The external incidence input is the original planar point-line
theorem. Bernert–Browning–Lichtman–Teräväinen supply the
prime-power decomposition used in the underlying repository
counts. Their exceptional-set theorem bounds individual
triples. The new claim here is a restriction on pairs of
triples near the repository's small-budget endpoint.
[Original planar incidence theorem](https://trotter.math.gatech.edu/papers/38.pdf)
· [Bernert–Browning–Lichtman–Teräväinen, v2, Section 2](https://arxiv.org/html/2410.12234v2#S2)

## A two-cap obstruction with actual diffuse least repairs

The target cutoff and two of the three original radical caps allow
counterexamples to a bounded physical-difference extraction. The
construction has actual primitive sums, comparable heights and
diffuse valuation profiles on the two capped letters.

Fix \(m\ge3\), set \(q=2^m\), \(e_i=q(m+i)\) for \(1\le i\le m\),
\(D=\sum_i e_i\), and \(L=3m/D\). Fix a positive \(\vartheta\) so
small that \((1+\vartheta)^D<11/10\) and
\(2^{-1/D}(1+\vartheta)<1\). Take an inner parameter \(Q\to\infty\).
For each slot \(i\), take a prime palette in its own nonempty
subinterval of \([Q,(1+\vartheta)Q]\). Take \(m\) more palettes
in disjoint subintervals of
\([2^{-1/D}Q,2^{-1/D}(1+\vartheta)Q]\). In both types,
order the disjoint subintervals increasingly with the slot \(i\).
The prime number theorem lets all \(2m\) palettes have the same
size \(J=Q^{1+o_Q(1)}\) at fixed \(m,\vartheta\).
A word in \(\{1,\ldots,J\}^m\) chooses a prime \(p_i\) and
an \(s_i\) in the corresponding palettes of each type. Set
\[
      c=\prod_i p_i^{e_i},\qquad a=\prod_i s_i^{e_i},
      \qquad b=c-a,\qquad X=((1+\vartheta)Q)^D .
\]
Every row is positive primitive, since the palettes of the
two types are disjoint. We have
\(Q^D\le c\le X\), \(Q^D/2\le a\le X/2\), and
\(b\ge Q^D-X/2\), so all three heights are comparable.
The increasing orders and the increasing logarithms of the chosen
primes give \(\log c\ge (D/m)\log\operatorname{rad}(c)\), by
the rearrangement inequality for these two lists. The same holds
with \(a\) in place of \(c\). Thus their radicals are each at
most \(c^{L/3}\), and at most \(n^{1/3}\) for \(n=X^L\).
Put \(M=e_m\) and \(Y=n^{M/3}\). Each target
uses the least repair \(H_M(c)=\operatorname{rad}(c)^M\le Y\).
There are \(n^{1/3+o_Q(1)}=Y^{1/M+o_Q(1)}\) distinct primitive
quotients and target radicals. Each cap has the wholly
nonatomic limiting law \((1/3)\operatorname{Unif}[2,4]\)
when \(m\to\infty\) after the inner height limit, since
\[
                          Le_i=\frac{6(m+i)}{3m+1}.
\]

For any fixed \(B>0\), the number of ordered distinct target
pairs with \(\operatorname{rad}(|c-c'|)\le n^B\) is at most
\[
                         n^{1/3+B/(2(m-1))+o_Q(1)} .
 \tag{SC7}
\]
To prove this, use the exact \(q\)-th roots of the two targets
and bin by the number \(j\) of slots shared. There are at most
\(Q^{2m-j+o_Q(1)}=n^{2/3-j/(3m)+o_Q(1)}\) ordered pairs
in the bin. Cancel the common gcd of the roots to get
coprime \(A_0,B_0\), ordered so \(A_0<B_0\). Their difference
of \(q\)-th powers divides the original difference. Its
\(m-1\) nonlinear factors
\(\Phi_{2^t}(B_0,A_0)\), \(2\le t\le m\), have disjoint prime
support outside \(2\). At least one has radical at most
\(O_m(n^{B/(m-1)})\). The fixed-support count enumerates its
values up to \(n^{B/(m-1)+o_Q(1)}\), and there are subpower
many pairs \(A_0,B_0\) per such value by the fixed cyclotomic
norm argument. The normalized roots identify which slots
were removed because all palettes are disjoint, and the
removed primes have at most \(Q^{j+o_Q(1)}\) choices.
The second bound on this bin is therefore
\(n^{B/(m-1)+j/(3m)+o_Q(1)}\). The minimum of the two
exponents is at most their mean, which proves (SC7) after
summing over \(j\). [Cyclotomic norm and fixed-value count](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#primitive-cyclotomic-factors-and-their-fixed-value-counts)

The same proof bounds the pairs with bounded physical difference
radical on the capped summand. The remaining summand meets its
required cap on only a power-small part of the construction. The roots of \(a,c\)
are already coprime, so if \(\operatorname{rad}(b)\le n^{1/3}\)
the same factor-count without a common-gcd choice gives at
most \(n^{1/(3(m-1))+o_Q(1)}\) rows. For late \(m\), (SC7)
excludes bounded-radical difference abundance while even the
two actual original caps and the central least-repair target
cutoff are retained. A three-cap input is necessary in any
such forcing argument.
