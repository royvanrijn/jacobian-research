# Minimum total degrees of alternating covers with prime half-degree

For n=2p with p>=7 prime, the least total degree of a union of
transitive permutation actions of A_n in which every group element
fixes a point is exactly 2^(n-2). The optimal action is unique up to
isomorphism: it consists of the odd subset actions of sizes less
than p and the action on unordered partitions into two p-element
sets.

The proof excludes costly primitive and matching actions using
Maroti's order bound, then forces every remaining constituent with
two-cycle witnesses. It imports the classification of finite simple
groups through that order bound; it does not use a classified list
of maximal subgroups of A_n.

## Theorem 1. The minimum and its equality case

Let p>=7 be prime and n=2p. For a finite family of proper subgroups
H_1,...,H_r of A_n suppose

\[
                A_n=\bigcup_{i=1}^r\ \bigcup_{g\in A_n}gH_i g^{-1}.
 \tag{1}
\]

Then

\[
                \sum_i[A_n:H_i]\ \ge\ 2^{n-2}.              \tag{2}
\]

Equality is attained precisely, up to conjugation and reordering, by

\[
 U_k=(S_k\times S_{n-k})\cap A_n
             \quad(1\le k<p,\ k\ {\rm odd}),\qquad
 W_p=(S_p\wr S_2)\cap A_n.                                 \tag{3}
\]

Here W_p stabilizes an unordered partition into two p-element sets.
In particular an optimal family has (p+1)/2 constituents.

**Upper bound.** If an even permutation has an odd cycle of length
d<p, its support is stabilized by U_d. If d>p its complement gives
U_(n-d), and if d=p its support gives a partition fixed by W_p.
If every cycle has even length, color each cycle alternately. The
two color classes have size p and are exchanged by the permutation,
which again fixes an unordered partition. Thus (3) covers A_n.

The indices are binomial(n,k) and binomial(n,p)/2. Complementation
pairs all odd subsets, including the middle layer since p is odd.
Consequently

\[
 \sum_{\substack{k<p\\k\ {\rm odd}}}\binom nk+
                              \frac12\binom np=2^{n-2}.    \tag{4}
\]

**Cost cutoffs.** Write R=2^(n-2). For every even n>=14,

\[
 n!>2^{n-1}3^n,\qquad (n-1)!!>2^{n-2}.                       \tag{5}
\]

Both strict inequalities hold at n=14:

\[
 14!=87178291200>39182082048=2^{13}3^{14},\qquad
             13!!=135135>4096=2^{12}.
\]

Increasing n by two multiplies their left sides by (n+1)(n+2)
and n+1, respectively, while their right sides grow by 36 and
four. These factors prove (5) inductively.

By [A. Maroti, *On the Orders of Primitive Groups*, Corollary 1.2,
printed page 3](https://renyi.hu/~maroti/primitive.pdf),
a primitive subgroup H of S_n not containing A_n has |H|<3^n
in every degree. The [retained primary source](sources/maroti-primitive-orders-2002.pdf)
has been inspected at that corollary; its proof uses the classification
of finite simple groups. Therefore every proper primitive H<A_n
has

\[
             [A_n:H]>\frac{n!}{2\cdot3^n}>R.                 \tag{6}
\]

A transitive imprimitive subgroup in degree 2p preserves a block
system of block size two or p, since p is prime. It lies in a
conjugate of

\[
 W_2=(S_2\wr S_p)\cap A_n
       \quad\hbox{or}\quad W_p.
\]

The index of W_2 is

\[
        [A_n:W_2]=\frac{n!}{2^p p!}=(n-1)!!>R.             \tag{7}
\]

All these index formulas hold in A_n since the corresponding
symmetric stabilizers contain odd permutations.

**Forced constituents.** It suffices to consider a covering family
whose total index is at most R. Equations (6)-(7) exclude primitive
subgroups and subgroups contained in W_2 from that family.
Every constituent is consequently intransitive or contained in W_p.

For each odd k<p choose a permutation g_k of cycle type (k,n-k),
including its one fixed letter when k=1. Both cycles have odd
length, so g_k belongs to A_n and has odd order. It cannot lie in a
conjugate of W_p: an odd-order permutation cannot exchange the two
blocks, and it has no invariant p-element subset.

It follows that g_k must be covered by an intransitive H. The
H-orbits are unions of the two g_k-orbits. Since H is intransitive,
these two orbits remain separate, so H is contained in a conjugate
of U_k. That constituent has cost at least binomial(n,k).

One constituent cannot serve two different k,l in this range.
Indeed a permutation of type (l,n-l) preserves a k-element subset
only when k=l, because its only proper nonempty invariant-subset
sizes are l and n-l. Thus distinct odd k force distinct
constituents, with total cost at least the sum in (4).

Finally g_p has type (p,p). None of these forced constituents
covers it, since it has no invariant subset of odd size less than p.
It needs another constituent. If that subgroup is intransitive,
its two orbits must be the two p-element supports, giving cost at
least binomial(n,p). Otherwise it is contained in W_p, giving cost
at least binomial(n,p)/2. This proves (2) by (4).

If the total cost is exactly R, every lower bound just used is
equality. Each forced subgroup therefore equals a conjugate of U_k,
and the final subgroup equals a conjugate of W_p. There is no
room for any additional constituent. Conversely (3) attains R,
which proves the equality statement. ∎

This theorem concerns total coset degree. It does not assert that
(p+1)/2 is the unrestricted minimum number of subgroup classes:
a family with greater total degree is outside the cost cutoff
argument used for the equality classification.

## Corollary 2. Minimum ranks of alternating Hasse algebras

Let K/Q be a Galois extension with group A_(2p), p>=7 prime.
Every finite etale Q-algebra split by K that has a point over
every completion of Q and no Q-point has rank at least 2^(2p-2).

If equality holds, its field-factor multiset, relative to an
identification Gal(K/Q)=A_(2p), is exactly

\[
 \left\{K^{U_k}:1\le k<p,\ k\ {\rm odd}\right\}
                                  \ \cup\ \{K^{W_p}\}.     \tag{8}
\]

**Proof.** Write the algebra as a product of fields K^H.
Absence of a rational point makes every H proper. At an unramified
prime, a local point is equivalent to its cyclic decomposition
group fixing a coset in at least one of the actions A_n/H.
Every conjugacy class occurs as an unramified Frobenius by
[Chebotarev's theorem, Lenstra's notes, Theorem 3.1, page 3](sources/lenstra-chebotarev-density-notes.pdf).
Thus these subgroup classes cover A_n. Theorem 1 gives both
the lower bound and its equality classification. ∎

When additionally n=2p=m^2+1 for an odd integer m, the
[Arala Hasse construction](ARALA_LAURENT_PENCILS_AND_KELLER_FIBRES.md),
Theorem 3, attains this minimum for infinitely many pairwise
linearly disjoint K/Q. Its Corollary 4 places the entire
minimum-rank family inside one fixed polynomial Keller map of
four-dimensional affine space.

For n=14 the unique optimal factor degrees are 14,364,2002,1716,
which sum to 4096=2^12. For n=22 the minimum is 1048576=2^20.
For n=26 the unique optimal factor degrees are

\[
       26,\ 2600,\ 65780,\ 657800,\ 3124550,\ 7726160,\ 5200300,
\]

which sum to 16777216=2^24.
The earlier degree-ten minimum 256 is established by its separate
classified finite calculation; degree ten is outside the hypothesis
p>=7 here. Degrees 26,82,122,226 also lie in the admissible
prime-half intersection. No assertion that infinitely many primes
of the form (m^2+1)/2 exist is required or made.

## Exact controls and scope

The [checker](verify_minimum_alternating_prime_halves.py) verifies the
integer bases and induction ratios in (5), all index and sum formulas,
the two-cycle witness incidence and equality costs for finite prime
controls. It independently checks every even cycle type in degrees
14,22,26 and 34 against the proposed cover. These finite controls support
the written proof; they do not prove its all-prime quantifier or
Maroti's classification-dependent bound.

The [receipt](minimum_alternating_prime_halves_verification.json)
and [source ledger](minimum_alternating_prime_halves_source_ledger.json)
record exact arithmetic, inputs and imported scope. No formal
verification or external human review is claimed.
