# Real paired Laurent covers and matching walks

Real existence is uniform in every even degree at least eight. The construction
can be made odd, so its numerator is an even polynomial. A matching-walk
description classifies the corresponding real chambers and determines exactly
which real-root-count spectra can occur. Coefficients can be chosen real
algebraic; rational coefficients in degree twenty-six remain unconstructed.

These are statements about the passport
((2,2,1^(n-4))^(n/2),(n-1,1)). The real-cover machinery is established
Hurwitz theory. In particular El Hilany and Rau's
[2019 primary paper](https://arxiv.org/abs/1712.05639v2) treats precisely
rational functions with one simple finite pole, and Theorem 2.4 gives their
real-dessin realization. Here the elementary matching reduction specializes
that framework to paired simple critical points. No priority claim is made
for real existence or the matching description.

## Theorem 1. A finite matching-walk classification

Let n=2b>=8. Fix b distinct ordered real numbers s_1<...<s_b.
Consider connected degree-n real covers of the target projective line with
poles of orders (n-1,1), precisely these b finite branch values, precisely
two simple critical points over every s_i, and all those n critical points
real. Quotient by real source isomorphisms,
while keeping the target fixed.

Put q=b-1 and define matchings on {0,...,n-1} by

\[
 A=\{A_i=(2i-1,2i):1\le i\le q\},\qquad
 B=\{B_i=(2i,2i+1):1\le i\le q\}.
\]

The covers are classified by the following finite objects. Choose one
auxiliary edge e=(0,j), 1<=j<=n-1. Starting at A and ending at B, perform
b steps, each toggling two disjoint edges of a matching. Every path edge
A_i,B_i is toggled once, and e is added once and removed once.
An absent edge may be added only when both endpoints are free.
There are no further edges or relabeling identifications in this normal form.

All these walks are transitive. For n>=10 their monodromy is A_n.
For n=8 all 72 walks have monodromy A8, as certified in Proposition 5.

**Proof.** Trivialize the cover over the upper open hemisphere, and label
its n sheets. Reflection across the target real axis identifies its lower
sheets with the same labels. Over each of the b+1 real intervals the
gluing is an involution c_i; its transpositions form a matching, and its
fixed labels are exactly the real points of a regular fibre.

At a real branch level, the two simple critical points correspond to a
double transposition sigma_i. A critical transposition (u,v) is real
exactly when c_{i-1} either fixes u,v individually or swaps them.
Consequently a level with two real critical points toggles two disjoint
matching edges, with

\[
                       c_i=c_{i-1}\sigma_i.
\]

The two outer intervals each have two real points: one near the simple
pole and one near the odd-order pole. Thus c_0 and c_b each have b-1
edges. Their product is the finite branch product, a cycle of length n-1
with one fixed label.

The union of two matchings has components that are alternating paths,
alternating cycles, or common edges. On an alternating path the product
is one cycle on its vertices. On an alternating cycle it gives two
cycles; a common edge gives two fixed labels. Therefore a product of
type (n-1,1) forces their union to be exactly an alternating path on
n-1 labels and one isolated label. Put the isolated label at 0 and
orient the path from the endpoint fixed by B to the endpoint fixed by A.
This gives precisely the displayed A,B, with no residual automorphism.

Every one of its n-2 path edges has different initial and final status,
so must be toggled at least once. Transitivity requires an edge incident
to isolated label 0. This edge is absent at both endpoints, so must be
toggled at least twice. The b steps permit exactly n toggles. They force
each path edge once and just one auxiliary edge twice; no other event
is possible. Conversely their union is a connected tree, so the
generated permutation group is transitive.

For completeness, the converse construction is explicit. Take n copies of
the upper closed hemisphere and n of the lower one, and glue their open
real intervals using the matchings c_i. At each s_i the transition is
the stated double transposition, giving two simple ramification points.
The two toggled edges are fixed by the induced reflection at their
ramification points, so both points are real. At infinity the transition
is AB or its inverse and has pole profile (n-1,1).
The total ramification is n+(n-2)=2n-2. Connectivity and Riemann--Hurwitz
give genus zero. The pulled-back complex structure makes reflection
antiholomorphic; the complex realization is the usual Riemann existence
theorem. Its real locus is nonempty, so the genus-zero real source is
the real projective line. Both poles are real and have distinct orders.
Real source coordinates put them at infinity and zero, yielding

\[
                         f(z)=P_{n-1}(z)+a/z,\qquad a\ne0.
\]

This also proves the asserted bijection: a labeled upper-hemisphere
trivialization gives exactly the gluing involutions, and source
isomorphism changes all labels simultaneously. The canonical path has
removed precisely that freedom.

For n>=10 the (n-1)-cycle forces primitivity: a block containing its
fixed sheet is invariant under this cycle and hence is either that
singleton or the whole set. The primitive double-transposition theorem
in Jones, Proposition 6.1(3), gives A_n. Every branch permutation is even,
so this is the full group. Its exact statement and primary-source
inspection are retained in the
[existing priority ledger](alternating_hurwitz_priority_ledger.json). ∎

The ordered real branch-configuration chamber is contractible. The usual
Hurwitz covering therefore has one real chamber for each walk, with no
branch collision or ramification change allowed. This is a classification
of this all-real locus with fixed target and ordered real branch values.
It does not classify the complex Hurwitz components or rational points.

There is also a unique monic Laurent representative for each such real
cover after its poles are put at zero and infinity: the remaining source
scaling is fixed by its leading coefficient, since n-1 is odd.
Thus the chamber statement equally applies to that normalized real
coefficient locus.

## Theorem 2. Odd covers are reversible walks

Suppose the branch set is stable under s->-s. Let r fix 0 and send
j to n-j on the path labels. A cover in Theorem 1 has an odd Laurent
representative if and only if its matching walk satisfies

\[
                         c_{b-i}=r c_i r
                         \quad(0\le i\le b).                   \tag{1}
\]

Its auxiliary edge is necessarily (0,b). In every even degree n>=8,
such a walk exists with maximum real-fibre size exactly six.

**Proof of the criterion.** The target involution s->-bar(s) preserves
the upper hemisphere and reverses its real intervals. A permutation r
lifts it exactly when (1) respects every gluing. The resulting
antiholomorphic involution commutes with ordinary reflection, so their
product is a holomorphic involution over s->-s. It fixes the two poles,
which have different orders. On real genus-zero source coordinates
with these poles at zero and infinity, this nontrivial involution is
z->-z. Hence f(-z)=-f(z).

Conversely an odd real f gives just this lift and its induced sheet
permutation. It interchanges A and B, so on their canonical path it
must be the displayed reversal. The isolated label stays fixed.
The unique auxiliary edge must be stable under r. Since r's only
fixed path label is b, that edge is (0,b).

**A peak-six construction.** Write b=2h or 2h+1, with h>=2,
and retain q=b-1. All operations below toggle the two indicated disjoint
edges; a path A-edge is removed and a path B-edge is added.

If b=2h, the first half removes the A-edges with indices h-1,...,q
and adds the auxiliary edge (0,b). First remove A_h and A_q together.
Next add the auxiliary edge while removing the largest-index remaining
edge in that range. Pair the other removals, in decreasing index order,
with additions B_q,B_(q-1),... in decreasing order.
There are h steps. Its final matching is

\[
 \{A_1,\ldots,A_{h-2}\}\ \cup\
 \{B_{h+2},\ldots,B_q\}\ \cup\ \{(0,b)\},
\]

which is r-invariant. Complete the walk by reversing those h steps and
reflecting every edge under r.

If b=2h+1, instead remove the A-edges with indices h,...,q in
the first half. Begin by removing A_(h+1),A_q, then add (0,b) while
removing the largest remaining index; pair the remaining removals with
B_q,B_(q-1),... as above. The h-step midpoint is

\[
 \{A_1,\ldots,A_{h-1}\}\ \cup\
 \{B_{h+3},\ldots,B_q\}\ \cup\ \{(0,b)\}.
\]

Insert the central step removing A_(h-1) and adding B_(h+2).
These are disjoint, and this step takes the midpoint to its reflection.
Finish with the reflected reversed first half.

All additions are allowed: the central A-edge incident to b is removed
before the auxiliary edge is added, and each added B_i has its incident
A-edges already removed. Apart from the first two special removals, a
decreasing A-index is at least two less than its paired B-index.
Hence their supports are disjoint. These descriptions also cover h=2,
where the indicated residual ranges are empty.

The first step removes two edges and raises the real-fibre size from
two to six. Every subsequent first-half step is neutral, as is the
central step when present. Reflection completes the pattern.
The walk has (1), uses the auxiliary edge, and is transitive. ∎

In particular, writing f=(H+a)/z gives an even H with H(0)=0.
The construction supplies the symmetric degree-26 coefficient problem
over the reals; it gives no rational coefficients.

## Theorem 3. Complete real-root-count spectra

For the covers of Theorem 1, the possible sets of sizes of regular real
fibres are exactly

\[
 \{2,6,10,\ldots,M\},\qquad
 M\equiv2\pmod4,\quad 6\le M\le
 \begin{cases}
 n-2,&n\equiv0\pmod4,\\
 n,&n\equiv2\pmod4.
 \end{cases}                                                  \tag{2}
\]

For odd covers, they are exactly the same sets with the sharper bound

\[
                         M\le
 \begin{cases}
 n-2,&n\equiv0\pmod4,\\
 n-4,&n\equiv2\pmod4.
 \end{cases}                                                  \tag{3}
\]

All spectra can be realized for any ordered real branch set in the
general case and any negation-symmetric one in the odd case.

**Necessity.** In a matching walk the real-fibre size is
R_i=n-2|c_i|. It starts and ends at two and changes by -4,0,or4.
Its spectrum is therefore exactly the initial segment ending at its
maximum M. A maximum of two is impossible. Initially no B-edge can be
added. A neutral first step would have to add the auxiliary edge between
the two initially free vertices 0,n-1 and remove one A-edge. The only
free vertices would then be the endpoints of that removed A-edge,
which do not form a B-edge. Thus no B addition is yet possible, and
the next allowed step must remove two edges, giving six real points.
Otherwise the first step already removes two edges. This proves M>=6.
Congruence and R_i<=n give (2).

For an odd cover with n=4h+2, zero is a branch value and its critical
points are a pair x,-x. If they are real, their folds have opposite
directions; near zero at most n-2 roots are real and congruence improves
this to n-4. The number of positive branch values is h. Reaching n
from the outer value two would require all h positive crossings to
increase the count by four. No crossing would remain to reduce the
count before zero, a contradiction. Hence M<=n-4. For our all-real
locus the critical pair is real; the
[earlier uniform argument](ALTERNATING_HURWITZ_UNIFORMITY_AND_REAL_NODES.md)
also treats nonreal critical pairs over real values.

**Every general spectrum occurs.** Put e=(0,1) and M=4h+2.
If q>=2h+1, remove A_1,...,A_(2h) in consecutive pairs; this gives
h steps and reaches M. Then add e while removing A_(2h+1).
For i=1,...,q-2h-1, add B_i while removing A_(2h+1+i).
Remove e while adding B_(q-2h). Finally add the remaining 2h B-edges
in consecutive pairs. Every middle step is neutral, and the last steps
reduce R by four. These are b valid steps and the maximum is M.

The remaining possibility is q=2h, so n=4h+2 and M=n.
Remove all A-edges in consecutive pairs, add e and B_1, remove e
and add B_2, then add the remaining B-edges in pairs. This again has
b steps, all allowed, and maximum n.

**Every odd spectrum occurs.** Theorem 2 gives maximum six in every
degree. There is a symmetric extension raising both degree and maximum
by four. From a walk of degree n, add two new path labels at each end,
send old path label j to j+2 and keep isolated label 0 fixed.
Start by removing the two new outer A-edges, perform the shifted old
walk, and finish by adding the two new outer B-edges.
This is a reversible degree-(n+4) walk. During the shifted old walk,
the real-fibre size is the old size plus four; the two outer intervals
have size two, so its maximum is M+4. Starting with the peak-six
construction in degree n-4k realizes M=6+4k whenever that base degree
is at least eight in the divisible-by-four case or ten otherwise.
These are exactly all maxima in (3). ∎

General degree-ten covers can consequently have a totally real regular
fibre. The odd degree-ten pencil cannot: its spectrum is {2,6}.
The distinction is the imposed source/target involution.

## Theorem 4. Uniform real algebraic maximal-node models

For every even n>=8 and every negation-symmetric ordered set of b=n/2
distinct real numbers, there is a monic even real polynomial H of degree
n, with H(0)=0 and a real nonzero, such that:

* f=(H+a)/z has full geometric monodromy A_n;
* its n finite critical points are simple and real, paired two at each
  chosen branch value;
* the line t=a passes through b real ordinary nodes of the reduced
  degree-n discriminant curve, with both branches real and transverse.

The maximum b of collinear ordinary nodes is attained. If the branch
values are rational, H and a may be chosen in a finite real number field.
For admissible n=m^2+1, m odd, the discriminant is a literal square over
that field. No rational model over Q is asserted.

**Proof.** Theorem 2 and the gluing argument provide the real odd cover.
Its two real poles normalize it to the stated form. The monic
normalization is a real source scaling. Its group is A_n by Theorem 1,
including its degree-eight finite control. The normalization of the
discriminant curve is

\[
                 z\longmapsto(H'(z),zH'(z)-H(z)).
\]

At a simple critical point r of f, H''(r) is nonzero, and its tangent
has slope r. Paired distinct real critical points therefore give two
smooth real transverse branches; r!=0 makes the horizontal line
transverse to each. The b nodes consume all n intersections, so attain
the line bound. These are the normalization and node arguments of the
[all-degree discriminant bridge](ALTERNATING_HURWITZ_DISCRIMINANT_BRIDGE.md).

Here is the real-algebraic coefficient assertion. For rational branch
values, arithmetic Riemann existence first gives an isomorphic cover
over Qbar. Its real reflection is algebraic over Qbar: in genus-zero
coordinates it is M composed with conjugation for a Möbius map M;
reflection sends algebraic points to algebraic points because it lies over
target conjugation and the map is finite over Qbar. Thus three algebraic
points have algebraic images, so M has algebraic coefficients. It gives descent from Qbar to the real algebraic numbers.
The unique simple pole is rational over this real closed field, hence
its genus-zero source is a projective line over that field.
The odd involution also descends, since it commutes with reflection.
Normalizing the two poles and leading coefficient over the real
algebraic field gives the desired even polynomial. Finitely many
coefficients lie in some finite real number field.

Finally the leading discriminant factor for monic degree n=2b is
(-1)^(b+1)(n-1)^(n-1). If n=m^2+1 with odd m, b is odd and this
equals the square m^(2n-2). Every chosen finite value contributes
multiplicity two. Thus the complete discriminant is a square over
the coefficient field. ∎

For n=26 one may, for example, prescribe the thirteen branch values
-6,-5,...,6. This proves existence over a real number field with all
critical points real. It does not provide the rational degree-26
coefficients requested by the arithmetic moonshot, or a bound on that
coefficient field.

## Proposition 5. Exact real chamber counts

For a fixed ordered real branch set, the counts of the covers in
Theorem 1 are:

| Degree | All-real covers | Peak 6 | Peak 10 |
|---:|---:|---:|---:|
| 8 | 72 | 72 | 0 |
| 10 | 1,848 | 840 | 1,008 |
| 12 | 63,576 | 7,848 | 55,728 |

For a fixed negation-symmetric ordered branch set, the odd-cover counts
are 4, 8, 96, and 269,668,672 in degrees 8, 10, 12, and 26 respectively.
The degree-twelve odd counts split into 36 with peak six and 60 with
peak ten. In degree twenty-six the exact split is:

| Maximum real roots | Odd covers |
|---:|---:|
| 6 | 32,296 |
| 10 | 9,108,696 |
| 14 | 83,956,320 |
| 18 | 128,943,360 |
| 22 | 47,628,000 |

These are finite exact counts, not rational-point counts.

**Finite reduction and computation.** A state records the subset of
removed A-edges, the subset of added B-edges, and the auxiliary phase
absent/active/removed. Every step chooses an unordered pair of eligible
events with disjoint supports. A B-edge is eligible after all incident
A-edges have been removed; if it meets the auxiliary edge, that edge
must also have been removed. Auxiliary addition is eligible after its
incident A-edge has been removed. This enumerates every matching walk
once. Unbounded integer dynamic programming also records its largest
R_i. Summing over j=1,...,n-1 gives the first table.

For the odd counts fix j=b and run only the first floor(b/2) steps.
If b is even, retain exactly the r-invariant middle matchings.
If b is odd, retain middle matchings differing from their r-reflection
by two disjoint edges. Insert that forced central step and reflect the
first half. The first-half maximum is already the full maximum.
These tests give exactly the reversible walks, rather than a subset.

A separate exhaustive degree-eight calculation supplies an independent
control. It enumerates all 93,296 DT quadruples with fixed product AB.
For each it tests all seven involutory conjugators of AB, retaining the
all-real transitive gluings. There are 504 labeled pairs, hence 72
absolute covers after division by the free centralizer of order seven.
The canonical 72 tuples all generate a group of order 20,160, hence A8;
four are reversible. This both checks the walk count and handles the
degree-eight exception to the imported primitive-group theorem. ∎

## Exact evidence and limits

The [standard-library checker](verify_real_paired_laurent_matching.py) and
[retained receipt](real_paired_laurent_matching_verification.json) record
the counts and degree-eight independent enumeration. The checker also
verifies every constructed spectrum in even degrees eight through sixty.
Those finite controls accompany the written uniform constructions;
they do not replace their all-degree proofs. An optional
--include-degree26-all run performs the approximately one-minute
unrestricted degree-twenty-six walk count as an additional control.
The retained default receipt uses the smaller all-real counts and the
complete odd degree-twenty-six count.

Run from the repository root:

~~~sh
.venv/bin/python research/arithmetic/verify_real_paired_laurent_matching.py
~~~

The [source ledger](real_paired_laurent_matching_source_ledger.json)
pins the direct real-Hurwitz comparison and existing indispensable
inputs. The general real-dessin correspondence is established prior
work; the matching classification, explicit schedules and counts are
proved here for the stated passport without asserting novelty.

The arithmetic gap is unchanged: no degree-twenty-six rational model,
uniform rational parametrization, rational-point theorem for higher
Hurwitz spaces, or new polynomial Keller realization follows.
The results give a uniform real-algebraic existence theorem, an exact
real combinatorial classification and complete real-root-count spectra.
The checks are internal and exact; no formal verification, independent
external replay or external human review is claimed.

