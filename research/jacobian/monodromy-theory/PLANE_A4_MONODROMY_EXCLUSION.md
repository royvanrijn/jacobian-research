# Every plane Keller action of A4 is excluded

Let k be an algebraically closed field of characteristic zero. A polynomial
Keller map F:A^2_k -> A^2_k has nonzero constant Jacobian. Its geometric
monodromy acts faithfully and transitively on its generic fibre, of degree

\[
 N=[k(x,y):k(F_1,F_2)].
\]

**Theorem.** The generic geometric monodromy group of such a plane Keller
map cannot be abstractly isomorphic to A4, in any faithful transitive action.
In particular the natural four-sheet action is impossible.

This is a consequence of Domrina's published degree-four theorem together
with the retained [point-stabilizer normal-generation theorem
KMFI1](THEORY.md#4-the-fixed-sheet-inertia-theorem). It makes no novelty
claim. The conclusion concerns generic field monodromy in dimension two;
it does not exclude A4 on affine charts, as a special fibre group, or for
polynomial maps in higher dimensions.

## The group-theoretic reduction

Write G=A4 and let H be a point stabilizer of a faithful transitive action.
Then the action is the coset action G/H, its degree is [G:H], and its kernel
is the core, the intersection of all conjugates of H. Thus H is core-free.
KMFI1 requires

\[
 \langle gHg^{-1}:g\in G\rangle=G. \tag{1}
\]

Here (1) is a necessary theorem for actual polynomial Keller maps, not a
sufficiency assertion about candidate inertia permutations.

The proper core-free subgroups of A4, up to conjugacy, are exactly 1, C2
and C3. To see that this list is exhaustive, Lagrange's theorem leaves
orders 1,2,3,4,6. There is no subgroup of order six: it would have index
two and give a nontrivial homomorphism A4 -> C2, whereas every three-cycle
must map to the identity and three-cycles generate A4. A subgroup of order
four consists of the identity and the three double transpositions, since
A4 has no elements of order four. It is the normal Klein subgroup V4,
so is not core-free. The subgroups of orders two and three are cyclic;
their distinct conjugates have trivial intersection.

| Stabilizer H | Degree [A4:H] | Normal closure of H in A4 | Consequence |
|---|---:|---|---|
| 1 | 12 | 1 | Violates (1). |
| C2 | 6 | V4 | Violates (1). |
| C3 | 4 | A4 | Requires the degree-four geometric exclusion. |

The middle row follows because conjugation permutes the three double
transpositions and they generate V4. The last row follows because the
conjugates of a subgroup C3 contain all eight three-cycles, which generate
A4. Therefore every possible plane Keller action of A4 would have N=4.
This reduction in fact applies in every dimension; its final exclusion
uses a theorem specific to dimension two.

## The unconditional four-sheet theorem

A. V. Domrina, *On four-sheeted polynomial mappings of C^2. II. The general
case*, Izvestiya: Mathematics 64:1 (2000), 1-33,
[DOI 10.1070/IM2000v064n01ABEH000273](https://doi.org/10.1070/IM2000v064n01ABEH000273),
proves that a polynomial map C^2 -> C^2 with four generic preimages cannot
have nonzero constant Jacobian. The exact locator is the **unnumbered main
theorem in the Introduction, printed page 1 / PDF page 1**. The preceding
paragraph defines topological degree to be generic fibre cardinality.
Thus the hypothesis is N=4, not a bound on the degrees of F_1 or F_2.
The [retained primary PDF](sources/domrina-four-sheeted-general-2000.pdf)
contains the inspected statement and proof.

The paragraph immediately following the theorem says that Part I treated
one dicritical component and that Part II treats the remaining cases.
Part II therefore supplies the general exclusion. Its Proposition 1.2,
printed/PDF page 2, reduces the remaining hypothetical map to two
dicritical components. Lemma 3.14 and Figure 11, printed/PDF pages 14-15,
list the final diagram types; Sections 4-7 prohibit them, ending with
Lemma 7.11 and the explicit conclusion of the main theorem on page 32.
We adopt this published theorem; we do not claim an independent replay
of its splice-diagram proof.

For completeness, the complex theorem excludes degree-four Keller maps
over every characteristic-zero field. The finitely many coefficients
of a hypothetical map lie in a field k0 finitely generated over Q, which
embeds into C. Generic degree is invariant under these field extensions:
after localizing the target at one nonzero polynomial, the source
coordinate algebra is a finite free module of rank N over the target,
and this rank remains N under field extension. Such a localization exists
by generic finiteness and generic freeness; the localized source remains
an integral domain after a field extension because it is a localization
of a polynomial ring. Nonzero constant Jacobian also persists under the
embedding. An N=4 map would therefore give the complex map prohibited
by Domrina. Apply this argument after the group-theoretic reduction to
finish the proof. No embedding of the entire original field into C is
required.

## The separate degree-five source obligation

The natural degree-five A5 action passes (1): its stabilizer A4 normally
generates the simple group A5. Its double transpositions have a fixed
sheet and normally generate A5 as well. These tests give no plane
realization.

The publisher's [primary article page for H. Zoladek,
*An application of Newton-Puiseux charts to the Jacobian problem*,
Topology 47:6 (2008), 431-469](https://www.sciencedirect.com/science/article/pii/S0040938308000207)
announces a result on topological degree. Secondary sources report the
cutoff five and point to Theorem 6.12. In this audit the full original
Theorem 6.12 and its reduction hypotheses have **not yet been inspected**;
the publisher PDF and author-hosted routes did not return the article.
The all-group degree-five statement remains a primary-source
verification obligation. The natural five-sheet \(A_5\) action has a
separate [written exclusion via the Orevkov budget](SIX_SHEET_BOUNDARY_BUDGET_AND_MONODROMY.md),
which forces \(S_5\) and leaves that symmetric case to be excluded.
The reported polynomial-degree gcd statements are outside these proofs.

## Exact finite controls and assurance

The [checker](verify_plane_a4_exclusion.py) enumerates every subgroup of
A4 from literal permutations, computes cores and normal closures, and
checks all three faithful transitive types. It also controls the necessary
normal-generation tests for the natural A4 and A5 actions and pins the
retained primary PDF. The [receipt](plane_a4_exclusion_verification.json)
records its actual execution. The general exclusion rests on the proof
above, KMFI1 and Domrina's theorem; finite group enumeration does not
verify the external complex-geometric proof. The
[source ledger](plane_a4_exclusion_source_ledger.json) records exact
locators and read/import limits. No formal verification, external review
or independent replay is claimed.
