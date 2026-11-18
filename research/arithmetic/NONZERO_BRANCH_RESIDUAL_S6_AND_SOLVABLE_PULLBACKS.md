# Nonzero branch residue fields obstruct solvable Hasse parametrizations

Every rational base change with solvable arithmetic monodromy leaves a
density-zero set of everywhere locally soluble rank-256 parameters in
every integer progression. This includes all substitutions
\(s=b+ct^d\), with rational \(b,c\), \(c\ne0\), and \(d\ge1\).
In particular, squaring the parameter removes the zero-branch obstruction
but does not repair the smaller subgroup cover.

The reason is structural: the residual sextic at a nonzero branch has
full symmetric group \(S_6\). Its \(A_6\) subgroup survives every solvable
residue-field extension and supplies an obstructing Frobenius class.
Riemann--Hurwitz prevents a rational map from killing the whole orbit
of four nonzero branch values.

The conclusion concerns positive-density parametrizations of the fixed
rank-256 subgroup product. It leaves thinner infinite families open.
The rank-301 everywhere locally soluble family is unaffected.

## Setting

Use the [regular degree-ten pencil](RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md)

\[
 P_s(z)=9z^{10}-864z^8+28161z^6-350150z^4
                 +1127500z^2-sz+250000.
\tag{1}
\]

Its normal closure over \(\mathbb Q(s)\) is regular with group \(A_{10}\).
For every full-group specialization with normal closure \(M_s\), put

\[
 A_s^{256}=M_s^{H_1}\times M_s^{H_3}\times M_s^{W_5}.
\tag{2}
\]

Here \(H_1,H_3\) stabilize a letter and a three-subset, and \(W_5\)
stabilizes an unordered partition into two five-element blocks. The
[subgroup-cover theorem](ALTERNATING_SUBGROUP_COVERS_AND_HASSE.md)
gives total rank \(10+120+126=256\). A completion has a point exactly
when its decomposition group fixes a coset of one of these subgroups.

A nonconstant \(\phi\in\mathbb Q(t)\) has **solvable arithmetic
monodromy** here when the normal closure of
\(\mathbb Q(t)/\mathbb Q(s)\), \(s=\phi(t)\), has solvable Galois group.
The arithmetic qualification includes its constant-field extension.

## Lemma 1. Exact nonzero branch data

Let \(a\) satisfy the irreducible quartic

\[
                   J(a)=a^4-\frac{109}3a^2+\frac{1690}9=0,
\]

and set \(K=\mathbb Q(a)\) and

\[
 \alpha=\frac{-29800a^3+1453400a}{3},\qquad
                    q_a(z)=z^2+az-\frac{10}3.
\tag{3}
\]

Then \(K=\mathbb Q(\alpha)\), and its four conjugates are exactly the
four nonzero finite branch values of (1). Moreover

\[
                       P_\alpha(z)=9q_a(z)^2G_a(z),
\tag{4}
\]

where

\[
\begin{split}
G_a(z)={}&z^6-2az^5+(3a^2-268/3)z^4+(-4a^3+172a)z^3\\
 &+(14251/9-199a^2/3)z^2+(298a^3/3-10034a/3)z+2500.
\end{split}
\tag{5}
\]

The quadratic and sextic are squarefree and coprime over \(K\).
Thus (4) has exactly two double roots and six simple roots.

**Proof.** The reduction \(J\bmod11=a^4+4a^2+2\) is irreducible,
as certified by modular powers and gcds in the retained checker.
Polynomial division of \(P_s\) by \((z^2+az-10/3)^2\) yields (3)--(5)
after reducing modulo \(J\). The powers \(1,\alpha,\alpha^2,\alpha^3\)
have nonzero determinant in the basis \(1,a,a^2,a^3\):

\[
                 \frac{608148060849659200000000000000000}{729}.
\]

Their characteristic polynomial is

\[
 W(S)=S^4-\frac{522379320400000}{243}S^2
                 +\frac{68130645548641000000000000000}{59049}.
\tag{6}
\]

It is the nonzero branch quartic already present in the discriminant
of (1). The determinant proves \(K=\mathbb Q(\alpha)\) and hence
irreducibility of \(W\). The same exact calculation checks the
quadratic discriminant and the required coprimality. The squarefree
sextic reductions below also certify its characteristic-zero
squarefreeness. ∎

## Theorem 2. The residual sextic has full group S6

The splitting field \(E/K\) of \(G_a\) has group \(S_6\).

**Proof.** At the following degree-one primes of \(K\), the exact
squarefree factor degrees of \(G_a\) are:

| Rational prime | Residue of a | Factor degrees |
|---|---:|---|
| 7 | 1 | 1,5 |
| 107 | 39 | 6 |
| 1151 | 291 | 1,1,1,1,2 |

These are roots of \(J\bmod p\), and \(J\)'s discriminant is a unit
at all three primes. All coefficients are integral there. The retained
input gives the literal factors; modular power/gcd tests certify each
factor's irreducibility and their product, without relying on a
factorization command. The Frobenius factorization theorem supplies a
5-cycle, a 6-cycle and a transposition in \(\operatorname{Gal}(E/K)\).

The 6-cycle makes this action transitive. A nontrivial block system
on six letters has block size two or three. A 5-cycle cannot act
nontrivially on its two or three blocks, nor nontrivially within such
small blocks. Hence the group is primitive. Conjugates of its
transposition form an invariant graph on the six letters; the
connected components are blocks. Primitivity makes that graph
connected, and its edge transpositions generate \(S_6\). ∎

## Lemma 3. The obstructing class survives solvable extensions

Let \(L/K\) be any finite extension whose normal closure is solvable.
After adjoining \(L\) and the two roots of \(q_a\), the residual sextic
group still contains \(A_6\), and therefore an element of type \((2,4)\).

**Proof.** Choose a finite Galois solvable extension \(F/K\) containing
\(L\) and the quadratic splitting field of \(q_a\). The intersection
\(E\cap F\) is Galois over \(K\), with solvable group. Every solvable
quotient of \(S_6\) kills its perfect subgroup \(A_6\). Consequently

\[
     \operatorname{Gal}(EF/F)
         \cong\operatorname{Gal}(E/(E\cap F))\ \supseteq A_6.
\tag{7}
\]

Passing to a smaller base field can only enlarge the residual group.
The permutation \((12)(3456)\) is even and belongs to \(A_6\).
For completeness, perfectness is elementary here: all 3-cycles are
conjugates of commutators of two double transpositions, for example
\([(12)(45),(23)(45)]=(132)\), and 3-cycles generate \(A_6\). ∎

## Lemma 4. Every solvable rational pullback has a surviving odd contact

For every nonconstant rational \(\phi\) with solvable arithmetic
monodromy, and for any chosen nonzero branch value \(\alpha\), there
is a finite point \(\beta\) with \(\phi(\beta)=\alpha\), odd
ramification index \(e\), and residue extension
\(K(\beta)/K\) having solvable normal closure.

**Proof.** For a degree-\(d\) rational map, a target fibre whose every
index is even contributes at least \(d/2\) to its ramification sum.
Four such fibres would contribute at least \(2d\), contradicting
Riemann--Hurwitz, whose total is \(2d-2\). The property of being such
a fibre is Galois invariant. Since the four nonzero values form one
irreducible orbit, none of their fibres can have every index even.
Thus the chosen fibre has an odd-index point. It is finite because
\(\phi(\infty)\) is rational or infinity, whereas \(\alpha\) has
degree four.

Base change the solvable Galois normal closure of \(\phi\) to \(K\).
Normalize the source and take a point over the target \(\alpha\).
Its residue Galois group is a decomposition group modulo its inertia,
hence solvable. The residue extension containing \(\beta\) is a
subfield of such an extension; its normal closure is solvable. This
argument also applies when \(\alpha\) is a branch value of \(\phi\). ∎

## Theorem 5. Density zero after every solvable rational base change

Fix such a \(\phi\), and integers \(b,m\) with \(m\ne0\). Among integer
indices \(n\) for which \(s=\phi(b+mn)\) is finite, unbranched and has
full group \(A_{10}\), those satisfying

\[
                \operatorname{Spec}(A_s^{256})(\mathbb Q_v)
                         \ne\varnothing\quad\text{for every }v
\]

have natural density zero. The progression also contains infinitely
many full-group parameters failing this local condition, with
pairwise linearly disjoint normal closures.

**Proof.** Choose \(\alpha,\beta,e\) from Lemma4. Apply Lemma3 to
\(L=K(\beta)\). In the compositum with a solvable Galois closure of
\(L\) and the critical quadratic, there is an automorphism fixing
\(\beta\) and both critical roots, and acting on the six residual
roots as \((2,4)\). Extend it to a finite normal closure \(\Omega/\mathbb Q\).
Chebotarev gives a positive-density conjugacy class of rational primes
with this Frobenius. Choose embeddings corresponding to that class.
Then \(\beta\), \(\alpha\), and the two critical roots lie in
\(\mathbb Q_p\); the residual sextic has factor degrees \((2,4)\).

Exclude the finitely many primes of ramification, denominators,
nonunit contact coefficients, root collisions, and divisors of \(m\).
All relevant algebraic data are fixed, so these exclusions are finite.
At every remaining prime the local expansion is

\[
           \phi(t)-\alpha=(t-\beta)^e u(t),\qquad
                               u(t)\in\mathbb Z_p^*,
\tag{8}
\]

on the residue disk of \(\beta\). If \(v_p(t-\beta)=1\), the right
side has odd valuation \(e\). Coprime Hensel factorization of (4)
gives two quadratic factors near its separate critical roots and
unramified factors of degrees two and four. The translated quadratics
have linear coefficient divisible by \(p^e\) and constant coefficient
of valuation \(e\): compare with their double-root factors modulo
\(p^e\), or evaluate \(P_s\) at the critical roots. Their discriminants
therefore have odd valuation, so both quadratic factors are ramified
and irreducible.

The decomposition group has root orbits \((2,2,2,4)\), with central
double-transposition inertia. Since it is a subgroup of \(A_{10}\),
a Frobenius lift on the two ramified pairs flips either both or
neither; multiplying by inertia chooses neither. Hence the local
group is \(C_2\times C_4\). It lies in none of the conjugates of
\(H_1,H_3,W_5\). There are no invariant odd subsets. Inertia has
fixed sheets and therefore preserves each block of any invariant
unordered \(5+5\) partition; its two ramified pairs must each be
monochromatic. The entire decomposition group preserves those pairs,
so it cannot exchange the two blocks. An invariant block would
then be a union of the even-sized orbits, impossible for size five.
This proves local failure of (2).

For each of these primes, \(v_p(b+mn-\beta)=1\) occupies exactly
\(p-1\) residue classes modulo \(p^2\). For any finite set of them,
the Chinese remainder theorem bounds the upper density of locally
soluble indices by

\[
                        \prod_p\left(1-\frac1p+\frac1{p^2}\right).
\tag{9}
\]

Positive prime density makes \(\sum_p1/p\) diverge. Exhausting the
prime set in (9) proves density zero.

Finally, the pullback retains regular arithmetic and geometric
\(A_{10}\). Indeed the original normal closure has no nontrivial
intersection with the solvable normal closure of \(\phi\): that
intersection would give a solvable quotient of the perfect group
\(A_{10}\). The same argument over \(\overline{\mathbb Q}\) proves
geometric connectedness. Select one of the forbidden congruence
classes at one prime. The integral Hilbert irreducibility theorem
gives infinitely many full-group integers in this refined progression.
Apply it successively to the compositum of the regular generic
extension with previously chosen normal closures. Its full group
is the direct product, giving linear disjointness while retaining
the local obstruction. ∎

## Theorem 6. Conditional transfer in every degree two modulo four

Let \(n=2b\ge10\), with \(b\) odd. Suppose a rational pencil
\(H(z)-sz+a_0\), \(a_0\ne0\), has regular normal closure \(A_n\).
Suppose one orbit of finite branch values has at least four members,
and at one of them, \(\alpha\), its factorization over
\(K=\mathbb Q(\alpha)\) is

\[
                   h_n q(z)^2G(z),
\]

where \(q\) is monic quadratic, \(G\) is monic of degree \(n-4\),
both are squarefree and coprime, and \(G\) has group \(S_{n-4}\) over
\(K\). For every solvable arithmetic rational base change, the
subgroup product

\[
       \prod_{\substack{1\le k<b\\k\text{ odd}}} M_s^{H_k}
                         \ \times\ M_s^{W_b}
\tag{10}
\]

has rank \(2^{n-2}\), and its everywhere locally soluble full-group
integer parameters have density zero in every progression. Each
progression contains infinitely many linearly disjoint full-group
failures of local solubility.

Here \(H_k\) stabilizes a \(k\)-subset, and \(W_b\) an unordered
\(b+b\) partition. These classes cover every element of \(A_n\).
The assertion is conditional on the specified rational cover and
its residual symmetric group; it constructs no new pencil.

**Proof.** An odd cycle gives an invariant odd subset of size at most
\(b\), taking its complement if necessary. Below \(b\) use \(H_k\);
at \(b\) use \(W_b\). If all cycles are even, alternate two colors
around each cycle to obtain a \(b+b\) partition exchanged by the
permutation. The symmetric-group orbits of these subsets and partitions
are also alternating-group orbits, since their stabilizers contain odd
permutations. The total coset degree is

\[
 \sum_{\substack{k<b\\k\text{ odd}}}\binom nk+\frac12\binom nb
             =\frac12\sum_{k\text{ odd}}\binom nk=2^{n-2}.
\]

The proofs of Lemmas3--4 apply with \(S_{n-4}\) and its perfect
subgroup \(A_{n-4}\). That subgroup contains the even permutation
with two cycles of lengths \(2,n-6\). After fixing an odd-contact
preimage and the critical roots, Chebotarev therefore supplies local
groups \(C_2\times C_{n-6}\), with root orbits
\((2,2,2,n-6)\). Every orbit is even, so no odd subset is fixed.
The same monochromatic-pair argument excludes an unordered \(b+b\)
partition, since \(b\) is odd. The Hensel, CRT and Hilbert arguments
of Theorem5 now give all the conclusions. ∎

This conditional theorem includes every admissible degree \(m^2+1\),
\(m\ge3\) odd, whenever its additional cover hypotheses hold.
The exact degree-ten calculation verifies those hypotheses for one
cover; no degree-26 existence or residual-group assertion is imported.

## Consequences and sources

There is a small exact control for the most immediate attempted repair,
\(s=t^2\). At \(p=347\), take \(a\equiv337\). Then
\(\alpha\equiv143\), the two critical roots are \(17,340\), and the
residual sextic has irreducible factor degrees \((2,4)\). The value
143 is a square, with square root 166. Hensel lifting modulo \(347^2\)
gives \(a=84311\), \(\alpha=10553\), and \(\beta=105654\).
For every

\[
                         t\equiv106001\pmod{347^2},
\]

the contact \(t^2-\alpha\) has valuation one: its residue modulo
\(347^2\) is \(115204=332\cdot347\). Thus every full-group
specialization in this progression has the displayed local obstruction.
The checker verifies the lifts, factors and contact exactly. This
control does not replace the all-solvable-map argument.

The substitutions \(\phi(t)=b+ct^d\) have solvable normal closures:
adjoin the \(d\)-th roots and \(d\)-th roots of unity. Thus Theorem5
applies to every such substitution, including \(s=t^2\), which erases
the original rational zero-branch inertia. It also applies to any
composition whose arithmetic normal closure is solvable.

The Chebotarev factorization and density statements are imported from
Lenstra, *The Chebotarev Density Theorem*,
[Fact2.1 and Theorem3.1, page3](sources/lenstra-chebotarev-density-notes.pdf).
Integral Hilbert irreducibility is Zywina,
[*Hilbert's irreducibility theorem and the larger sieve*, Theorem1.2, page2](sources/zywina-quantitative-hit-2010.pdf).
Their exact statements and use boundaries remain those of the
[earlier ramification source ledger](alternating_ramification_density_source_ledger.json).
The local Hensel argument is written above and extends the zero-branch
calculation in the [density note](ALTERNATING_PENCIL_RAMIFICATION_DENSITY.md).

The [checker](verify_alternating_solvable_pullback_obstruction.py),
[input](alternating_nonzero_branch_residual_inputs.json),
[receipt](alternating_solvable_pullback_verification.json), and
[source ledger](alternating_solvable_pullback_source_ledger.json)
retain the exact branch identities, all finite-field factor certificates,
the field-generation determinant and the local subgroup test. The
quantified pullback and density claims use the written arguments and
the named primary theorem imports, not finite sampling.

This does not exclude a thinner infinite set of locally soluble
parameters, arbitrary nonsolvable rational base changes, or different
subgroup products. It neither changes the rank-301 construction nor
claims an obstruction to generic alternating Keller maps. No novelty,
formal verification or external human review is asserted.
