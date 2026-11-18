# Prior covering results, uniform equality, and optimal Hasse ranks

The numerical minimum total coset degree used in the alternating Hasse
construction is a consequence of Maroti's 2005 covering theorem.
That theorem applies in every degree congruent to two modulo four,
without a prime-half restriction. The Arala-based Hasse family
therefore attains the globally minimum rank in every admissible degree.

This note records the exact prior-work overlap and proves the equality
classification in every such degree at least fourteen. The inspected
sources do not explicitly state that weighted uniqueness result.
Its priority remains unresolved; no novelty claim is made.

## The known numerical optimum

For a finite noncyclic group G, let sigma(G) be the minimum cardinality
of a covering by individual proper subgroups. Define

\[
 d_{\rm wt}(G)=\min\left\{\sum_i[G:H_i]:
 G=\bigcup_i\bigcup_{g\in G}gH_i g^{-1},\quad H_i<G\right\}.
\]

A family in this minimum gives an ordinary covering with at most
sum_i[G:N_G(H_i)] individual subgroups. Since H_i<=N_G(H_i),

\[
                         \sigma(G)\le d_{\rm wt}(G).          \tag{1}
\]

**Theorem A (consequence of Maroti).** For n>=6, n=2 mod4,

\[
                         d_{\rm wt}(A_n)=2^{n-2}.             \tag{2}
\]

**Proof and precise import.** Maroti,
*Covering the symmetric groups with proper subgroups*,
J. Combin. Theory Ser. A110 (2005), 97--111,
[doi:10.1016/j.jcta.2004.10.003](https://doi.org/10.1016/j.jcta.2004.10.003),
Theorem4.1, gives sigma(A_n)=2^(n-2) when n>2 is even and
not divisible by four. See the
[retained author manuscript dated15October2004](sources/maroti-covering-symmetric-2004.pdf),
printedp9, with proof and explicit covering pp10--11.
The n6 case is expressly included by the earlier Bryce--Fedri--Serena
result cited there as [4]; the remaining proof begins at n10.
The odd-degree exceptions seven and nine elsewhere in the paper do not apply.

Put b=n/2. The same proof uses precisely

\[
 U_k=(S_k\times S_{n-k})\cap A_n\quad(1\le k<b,\ k\ {\rm odd}),
 \qquad W_b=(S_b\wr S_2)\cap A_n.                            \tag{3}
\]

The latter stabilizes an unordered partition into two b-element sets.
An odd cycle shorter than b supplies U_k; a longer one supplies its
odd-sized complement. A b-cycle supplies W_b. If every cycle is even,
alternating colors along each cycle give two b-element sets exchanged
by the permutation. These conjugate subgroups therefore cover A_n.
Their total index is

\[
 \sum_{\substack{k<b\\k\ {\rm odd}}}\binom nk+
                             \frac12\binom nb=2^{n-2},       \tag{4}
\]

by pairing complementary odd subsets. Equation(1) and Maroti's theorem
give the reverse inequality. ∎

The numerical minima256 forA10,16777216 forA26, and the entire
prime-half numerical formula are thus known-theorem consequences
in a stronger degree range. The retained finite and elementary proofs
remain valid alternative proofs and controls.

## Exact comparison with the other primary literature

| Source and locator | Exact overlap |
| --- | --- |
| Maroti(2005), Theorem4.1 and proof, printedpp9--11 | Supplies the ordinary minimum and exactly the family(3); equation(1) yields the full weighted optimum(2). |
| Bubboloni--Praeger, *Normal coverings of finite symmetric and alternating groups*, arXiv:1011.4368v1, Proposition7.6(a),(c), printedpp23--24 | Gives gamma(A_(2p))=(p+1)/2 and exhibits both W_p and W_2 with the same odd intransitive classes. It minimizes the number of classes, rather than the sum of indices. |
| Bubboloni--Sonn, *Intersective S_n polynomials with few irreducible factors*, arXiv:1507.08593v1, Proposition1.1, printedpp1--2 | Gives the fixed-subgroup/decomposition-group criterion. Its r(G),s(G) count irreducible factors or subgroup classes, rather than total degree. |
| Sonn, *Two remarks on the inverse Galois problem for intersective polynomials*, JTNB21(2009), Proposition,p436 | Every realizable nonsolvable group occurs for an intersective polynomial. Its proof allows all proper subgroups and asserts no minimum rank. |

The inspected versions are retained as
[Bubboloni--Praeger](sources/bubboloni-praeger-normal-1011.4368v1.pdf),
[Bubboloni--Sonn](sources/bubboloni-sonn-intersective-1507.08593v1.pdf),
and [Sonn](sources/sonn-intersective-two-remarks-2009.pdf).
The [prime-half proof](MINIMUM_TOTAL_DEGREES_ALTERNATING_PRIME_HALVES.md)
independently proves weighted uniqueness in that subfamily.
Neither the inspected statement of Maroti's Theorem4.1 nor the stated
Bubboloni--Praeger Proposition7.6 asserts the weighted equality
classification below. This observation does not establish absence
from other literature or from a detailed equality analysis.

## A uniform cutoff for all other imprimitive actions

Let I(a,r)=(ar)!/((a!)^r r!). This is the index in A_(ar) of
the stabilizer of a partition into r blocks of size a, intersected
with A_(ar); the symmetric stabilizer contains odd permutations.

**Lemma B.** If r>=3,a>=2 and ar>=8, then

\[
                            I(a,r)>2^{ar-2}.                 \tag{5}
\]

**Proof.** For a>=3 start with I(3,3)=280>128. At fixed r3,

\[
 \frac{I(a+1,3)}{I(a,3)}
 =\frac{3(3a+1)(3a+2)}{(a+1)^2}>8,
\]

since the numerator minus eight times the denominator is
19a^2+11a-2>0 for a>=1. This proves(5) at r3 for all a>=3.

At fixed a>=2 and r>=3,

\[
 \frac{I(a,r+1)}{I(a,r)}
 =\frac{\binom{(r+1)a}{a}}{r+1}
 \ge(r+1)^{a-1}\ge2^a.                                     \tag{6}
\]

For the binomial bound divide the labeled letters into a prescribed
groups of size r+1 and choose one letter from each.
The final inequality uses r+1>=4 and a>=2.
Induction on r proves(5) for all a>=3.
For a2 start with I(2,4)=105>64; increasing r multiplies I
by2r+1>4. The omitted pair(a,r)=(2,3) has n6 and I=15<16,
so its omission is necessary. ∎

## Uniform classification of equality

**Theorem C.** For n>=14,n=2 mod4, every family attaining
d_wt(A_n)=2^(n-2) consists, up to conjugation and reordering,
of exactly the subgroups(3).

**Proof.** Write R=2^(n-2). For every even n>=14,

\[
                           n!>2^{n-1}3^n.                   \tag{7}
\]

At n14 the two sides are87178291200 and39182082048.
Increasing n by two grows the left side by(n+1)(n+2)>36
and the right side by36.

Maroti's
[*On the Orders of Primitive Groups*, Corollary1.2, printedp3](sources/maroti-primitive-orders-2002.pdf)
gives order less than3^n for a primitive subgroup of S_n not
containing A_n. This bound imports the classification of finite
simple groups. Every proper primitive H<A_n therefore has
[A_n:H]>n!/(2*3^n)>R.

LemmaB also excludes every imprimitive subgroup preserving three
or more blocks from a family of total cost at most R.
Every remaining constituent is intransitive or contained in W_b,
b=n/2, up to conjugation.

For each odd k<b, an element g_k of type(k,n-k) has odd order.
It cannot exchange two b-element blocks and has no invariant
b-element subset. It is therefore absent from W_b.
An intransitive subgroup covering g_k preserves its two supports,
so is contained in U_k and has index at least binomial(n,k).
The same constituent cannot cover g_l for a different odd l<b:
a permutation of type(l,n-l) has no invariant k-element subset
unless k=l.

The type(b,b) occurs in none of these U_k and requires another
constituent. An intransitive one has index at least binomial(n,b);
a subgroup of W_b has index at least binomial(n,b)/2.
Together with the forced odd costs, equation(4) gives a sum at
least R. Equality forces every containment to be equality, the
last subgroup to be W_b, and no additional constituent. ∎

The numerical part is prior work by TheoremA. This gives a complete
uniform proof of the equality classification; its priority remains
unresolved after the stated search. No novelty claim is made.

## Minimum Hasse ranks in every admissible degree

**Corollary D.** Let m>=3 be odd,n=m^2+1 and R=2^(n-2).
There are infinitely many pairwise linearly disjoint A_n extensions
K_j/Q for which

\[
 \mathcal A_j=
 \left(\prod_{\substack{k<n/2\\k\ {\rm odd}}}K_j^{U_k}\right)
                           \times K_j^{W_{n/2}}
\]

has rank R, a point over every completion of Q, and no Q-point.
Every such finite etale Q-algebra split by an A_n extension has
rank at least R. For n>=26, equality forces the displayed family.

For each admissible n, one fixed polynomial Keller map of A4_Q
contains infinitely many of these minimum-rank algebras as complete
fibres. Its generic arithmetic and geometric groups are S_R.

**Proof.** Here n=2 mod8, hence n=2 mod4.
The [Arala Hasse theorem](ARALA_LAURENT_PENCILS_AND_KELLER_FIBRES.md),
Theorem3, supplies the infinite disjoint family of rank R.
No change to its inertia or exceptional-prime argument is needed.

For any algebra split by K/Q with group A_n, write its field
factors as K^H. Absence of a Q-point makes every H proper.
Every conjugacy class occurs as an unramified Frobenius by
[Chebotarev's theorem, Lenstra's notes, Theorem3.1,p3](sources/lenstra-chebotarev-density-notes.pdf).
Local solubility therefore makes those subgroup classes cover A_n.
TheoremA proves the minimum rank, and TheoremC classifies equality
when n>=26. Corollary4 of the same Arala note supplies the single
A4 map, full-fibre isomorphism and generic symmetric group. ∎

Thus the Arala Hasse upper bound is globally minimum also in
admissible composite-half degrees, such as n50. The numerical
lower bound is prior work. The specific optimal arithmetic and
fixed-Keller transfer has no completed global priority assessment.

This does not construct a rational paired-simple degree26 pencil,
change the obstructions for that different double-transposition
passport, or make the ambient Keller group alternating.

## Evidence and priority limits

The [checker](verify_alternating_weighted_cover_priority.py) checks
the exact induction bases/ratios, finite imprimitive and primitive
cutoffs, rank identities and all even cycle types in selected degrees.
It also checks admissible composite-half ranks. The
[receipt](alternating_weighted_cover_priority_verification.json)
pins theorem locators and inputs.
Maroti's covering theorem and order theorem are imported, not
reproved by the finite computation.

The [source ledger](alternating_weighted_cover_priority_source_ledger.json)
records the precise overlap and inspected versions.
This is a source comparison, not an exhaustive novelty assessment.
No formal verification or external human review is claimed.

