# Subset-orbit forms and inverse-cover rationality

Work over an algebraically closed field k of characteristic zero.
Let M_c and M_s denote the A12 Galois-closure surface fields of
the level-two cone and alternating shear used in
[KMFI19](SURFACE_FORMS.md). Independent identity coordinates are
omitted from the notation.

**KMFI20.** If J is any subgroup of A12, and r_m(J) is its number
of orbits on m-element subsets of the twelve inverse labels, then
smooth projective models satisfy

\[
 p_g(M_c^J)\ \ge r_6(J)-r_1(J),\qquad
 p_g(M_s^J)\ \ge r_6(J)-r_2(J).
 \tag{1}
\]

Here p_g means the dimension of the space of regular two-forms.
These are lower bounds, not claimed exact genera.

**KMFI21.** For either family, every proper subgroup J<A12
transitive on the twelve inverse sheets has a target field M^J
which is not stably rational. Thus no connected degree-twelve
subgroup reduction of either audited family can produce a
polynomial map between ordinary affine spaces, even after
common independent stabilization and arbitrary birational recharts.

**KMFI22.** These same fields obstruct many other actions of the
original group A12 when formed by root configurations.
For the cone, the fields of unordered k-subsets or ordered
k-tuples of distinct inverse roots are not stably rational for
2<=k<=6. For the shear this holds for 3<=k<=6.
Their generic group is still A12, but these source fields do
not admit affine-space polynomial models.

These are statements about the two specified normal closures.
They are not exclusions of the abstract subgroup actions or of
the configuration permutation actions for every possible Keller
map. The full spectrum and the specified PSL2(11) occurrence
problems remain open.

The [checker](verify_orbit_forms.py), [model](orbit_forms_model.json),
[receipt](orbit_forms_receipt.json) and
[source manifest](sources/orbit-forms-manifest.json) retain exact controls.
The finite-group classification inputs in Section 5 are attributed
established results; the audit does not independently reprove their
completeness. No formal verification or independent whole-proof
review is claimed.

## 1. Whole families of global cut forms

For an m-subset A of the twelve root labels let

\[
 f_A(w)=\frac1{\prod_{i\in A,j\notin A}(w_i-w_j)},\qquad
 h_m=m(12-m).
\]

It is a nonzero rational function in M since the generic roots
are distinct. On the affine charts consider

\[
 \omega_{c,A}=f_A\,dP\wedge dQ,\qquad
 \omega_{s,A}=Y^4 f_A\,dP\wedge dY.
 \tag{2}
\]

The complete branch geometry and finite quotient descent have
already been proved in KMFI19. Here we vary the cut size.
If r roots form a cluster and k of them belong to A,
the denominator has k(r-k) differences within that cluster.
The permissible range is
max(0,r-(12-m))<=k<=min(r,m).

For all m=2,...,6, the cone forms in (2) extend globally.
Their projective expressions in weights (s,Q,P)=(1,6,11) are

\[
 s^{h_m-18}f_A\,\iota_E(ds\wedge dP\wedge dQ).
 \tag{3}
\]

The exponents are 2,9,14,17,18. They are all nonnegative.
The finite three-, four- and five-root cluster bounds are
respectively 2,4,6, no larger than the pulled-back volume orders
in KMFI19. At the weighted boundary the branch is smooth and
transverse, with at most a three-root collision; the same bound
applies. The only other charts are finite quotient charts.

For the shear, all m=3,...,6 work. In weights (s,Y,P)=(1,2,9)
their projective expressions are

\[
 s^{h_m-20}Y^4f_A\,\iota_E(ds\wedge dP\wedge dY),
 \tag{4}
\]

with exponents 7,12,15,16.
All finite bounds from KMFI19 remain valid:
the Y^4 factor cancels the unramified lost triple, and the
four-root origin has volume-plus-numerator order twelve,
matching its maximum cut order twelve.
At the boundary seven-root cluster, differences have order four
under the weighted (3,28) blowup.
The maximum cut is twelve, so a denominator order at most 48
is opposed by volume-plus-boundary order

\[
 30+3(h_m-20)=51,66,75,78.
\]

Thus all four values of m extend, including all exceptional
divisors and all weighted axis charts. The same resolutions and
finite quotient argument of KMFI19 apply term by term.

The raw m=2 shear form does not satisfy this boundary estimate:
it would give 30 against a possible denominator order 40.
No assertion about additional m=2 forms is made here.

Each form in (2) is nonzero: its scalar factor is nonzero and
the base coordinates are algebraically independent.
This observation will give a module-level nonvanishing proof
without relying on a special root evaluation.

## 2. Incidence modules: an elementary irreducibility proof

Let M_m be the permutation vector space on all m-subsets,
for 0<=m<=6. Define

\[
 U:M_{m-1}\to M_m,\quad
 U[B]=\sum_{A\supset B,\ |A|=m}[A],
 \qquad D=U^{\mathsf T}.
\]

These maps are equivariant for S12 and A12.
On level k of the subset lattice, the up/down operators satisfy

\[
 DU-UD=(12-2k)I.
\]

For k<6 this proves U injective over R:
||Uv||^2=||Dv||^2+(12-2k)||v||^2.
Since the matrices have rational coefficients, injectivity holds
in every characteristic-zero field.

There is a multiplicity-free decomposition

\[
 M_m=M_{m-1}\oplus V_m,\quad
 V_m=\ker D,\quad
 \dim V_m=\binom{12}{m}-\binom{12}{m-1}.
 \tag{5}
\]

Here each V_m is irreducible even for A12, and the V_m are
pairwise nonisomorphic. For completeness, the commutant dimension
of M_m for S12 is m+1: ordered pairs of m-subsets are classified
by their intersection sizes 0,...,m.
The same is true for A12. The stabilizer of any such ordered pair
in S12 contains an odd permutation, since its at most four
membership cells divide twelve points and one cell has size at
least three. Hence the S12 orbit does not split.

Now use induction and Maschke's theorem. M_(m-1) has m distinct
irreducible constituents, each once. Its injected image in M_m
accounts for m of the m+1 commutant dimensions. Semisimplicity
and the sum-of-squared-multiplicities formula force exactly one
new irreducible constituent, again once. This proves (5) and
irreducibility without needing an irreducible-character table.
These V_m are the usual two-row modules, but that terminology
is not required for the argument.

## 3. Why the cut forms contain the whole new constituent

For every (m-1)-subset B,

\[
 \sum_{\substack{A\supset B\\|A|=m}} f_A=0.
 \tag{6}
\]

Here is the exact rational identity.
Put R equal to the remaining 13-m labels and
D_B=product_(b in B,j notin B)(w_b-w_j).
For A=B union {i}, rewrite

\[
 f_A=D_B^{-1}
      \frac{\prod_{b\in B}(w_b-w_i)}
           {\prod_{j\in R,\ j\ne i}(w_i-w_j)}.
\]

The numerator is a polynomial of degree m-1 in w_i.
The Lagrange leading-coefficient identity makes its sum over i
zero, because m-1<|R|-1=12-m for m<=6.
Thus (6) holds universally, before specializing the root family.

Map [A] to the corresponding global form (2). This is an
A12-equivariant linear map. Equation (6) annihilates U(M_(m-1)),
so the map factors through V_m.
It is nonzero because each individual cut form is nonzero.
Irreducibility in Section 2 makes this map injective.
The different V_m are nonisomorphic, so all their images form
a direct sum.

Consequently the cone Galois-closure model contains copies of
V_2,...,V_6 in its regular two-forms, of dimensions

54,154,275,297,132, with total 912.

The shear contains V_3,...,V_6, of total dimension 858.
These are guaranteed subspaces; no equality with the full
space of forms is asserted.

For a subgroup J, the invariant dimension of M_m is r_m(J).
Taking J-invariants in (5) is exact in characteristic zero, hence

\[
 \dim V_m^J=r_m(J)-r_{m-1}(J).
\]

The invariant forms descend to a smooth projective model of M^J
by KMFI19's finite quotient argument. Summing these dimensions
telescopes to (1), proving KMFI20.
The resulting regular forms are a stable-rationality obstruction.

## 4. Exact subgroup controls

The checker obtains the following full orbit counts directly
from the retained permutation generators:

| subgroup on twelve labels | r1 | r2 | r3 | r4 | r5 | r6 | lower bound for either target p_g |
|---|---:|---:|---:|---:|---:|---:|---:|
| PSL2(11) | 1 | 1 | 1 | 2 | 2 | 6 | 5 |
| transitive M11 | 1 | 1 | 1 | 2 | 2 | 3 | 2 |
| M12 | 1 | 1 | 1 | 1 | 1 | 2 | 1 |

Thus the new method also excludes subgroup restriction to both
Mathieu groups for these two constructions.
For PSL2(11) it strengthens KMFI19's nonzero-form certificate to
a lower bound of five linearly independent forms.

This is compatible with rationality of the original source and
target. For G=A12 all r_m are one, so both lower bounds are zero.
For the point stabilizer A11, which fixes one label and is
transitive on the rest, r1=r2=r6=2, and the bounds are again zero.
A zero lower bound is not a general rationality criterion.

## 5. Every proper transitive subgroup is excluded

We use two established classification inputs.

1. Beaumont--Peterson, *Set-transitive permutation groups*,
   Canadian Journal of Mathematics 7 (1955), Theorem 10,
   [primary paper](https://doi.org/10.4153/CJM-1955-005-x):
   at degree twelve a group transitive on the six-subsets
   contains A12. The exceptional degrees in that theorem
   are 5,6,9.
2. The complete degree-twelve primitive-group list is
   M11, M12, PSL2(11), PGL2(11), A12, S12.
   The [official PrimGrp library](https://gap-packages.github.io/primgrp/)
   and its degree-twelve record are retained in the source manifest.
   Every group in this list is two-homogeneous.
   The M11 and M12 orbit assertions are checked from the stored
   generators, and the projective linear assertions follow from
   their usual two-transitive action.

The catalogue is an external completeness input. The checker
does not claim to classify primitive groups from a generator sample.
The retained file is pinned by its Git blob SHA and content hash.

For the cone, let J<A12 be proper and transitive.
Then r1=1 and r6>=2 by Beaumont--Peterson, so (1) gives p_g>=1.

For the shear, if J is primitive the complete list makes r2=1,
and Beaumont--Peterson again gives p_g>=1.

If J is imprimitive, a nontrivial block system has block size
b=2,3,4 or 6. Put W_b=S_b wr S_(12/b), the full setwise stabilizer
of that system. Then J is a subgroup of W_b and of A12.
The orbits of W_b on m-subsets are classified by the multisets
of their occupancies in the blocks. They give

| block size b | r2(W_b) | r6(W_b) | difference |
|---:|---:|---:|---:|
| 2 | 2 | 4 | 2 |
| 3 | 2 | 5 | 3 |
| 4 | 2 | 5 | 3 |
| 6 | 2 | 4 | 2 |

The formal S12 module V_3+...+V_6 therefore has a nonzero
W_b-fixed subspace. Its image in the global form space is
A12-equivariant, so this subspace is fixed by J as well.
It does not matter that W_b itself contains odd permutations:
only its even subgroup acts on the Galois field, and J lies there.
Thus the shear quotient again has p_g>=1.

This proves KMFI21 for all proper transitive J, not just the
three named almost-simple subgroups in Section 4.
Since J is transitive, the base-changed degree-twelve point field
stays connected and has normal-closure group J.
Its target M^J is not stably rational, so no independent
stabilization or birational rechart can make it an affine-space
polynomial map.

An intransitive J can split the twelve-sheet base change into
several components. Selecting such a component is a different
construction and is not covered by the transitive-subgroup
conclusion. The orbit inequality (1) itself still applies.

## 6. Other actions of A12: source obstructions for configurations

Let H_k=(S_k times S_(12-k)) intersection A12 be the stabilizer
of a k-subset, 1<=k<=6. In its action on m-subsets, for m<=6,

\[
 r_m(H_k)=\min(k,m)+1.
\]

The orbit is determined by intersection size with the distinguished
k-subset. Passing to the even subgroup does not split these
orbits: the four membership cells contain an odd stabilizer
permutation, just as in Section 2.
Consequently

\[
 p_g(M_c^{H_k})\ge k-1,\qquad
 p_g(M_s^{H_k})\ge\max(0,k-2).
 \tag{7}
\]

These fields are the connected unordered-k-root configuration
covers, of degree binomial(12,k) over the original rational
target, with generic group A12 in its k-subset action.
For k=6 the action has complementary-pair blocks; no primitivity
claim is being made there.

For ordered distinct k-tuples, the stabilizer fixes k labels
pointwise and is A_(12-k) on the rest. For k<=6,

\[
 r_6=2^k,\quad r_1=k+1,\quad
 r_2=1+k+\binom{k}{2}.
\]

Thus the source-genus bounds are

\[
 p_g(M_c^{A_{12-k}})\ge2^k-k-1,\qquad
 p_g(M_s^{A_{12-k}})\ge2^k-1-k-\binom{k}{2}.
 \tag{8}
\]

The degree is 12!/(12-k)!, and the ordered configuration is
connected since A12 is k-transitive in this range.

Every positive lower bound in (7) or (8) excludes a polynomial
affine-space model of that configuration field pair after
common independent stabilization and birational changes:
an affine source field must be rational.
This proves KMFI22.

For example, the cone's unordered-pair and ordered-pair covers,
of degrees 66 and 132, already have p_g>=1.
For the shear, the unordered-triple and ordered-triple covers,
of degrees 220 and 1320, have p_g>=1.
Therefore adjoining several inverse answers of an actual Keller
map does not automatically give another affine-space Keller map,
even with the original target fixed and generic group unchanged.

The statements do not rule out another normal closure realizing
one of these permutation actions. They also do not realize or
exclude a genuinely new PSL2(11), Mathieu, A4, or odd-dihedral
Keller action. The full objective remains open.
