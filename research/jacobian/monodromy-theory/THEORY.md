# Symmetry groups of polynomial Keller maps: affine inertia and realization

This note contains the fixed-sheet inertia theorem, its group-theoretic
consequences and the degree-twelve PSL2(11) local sieve. The complete
monodromy classification remains open. Mathematical scope and assurance
are recorded in [MATH_STATUS.json](../../MATH_STATUS.json).

## 1. The object being classified

Work first over an algebraically closed field k of characteristic zero.
A Keller map is a polynomial morphism F:X=A^m -> Y=A^m with
det DF in k*. It is dominant, generically finite and etale. Write

\[
 K=k(Y),\quad L=k(X),\quad M=\text{Galois closure of }L/K,\qquad
 G=\operatorname{Gal}(M/K),\quad H=\operatorname{Gal}(M/L).
\]

The degree is n=[L:K]=[G:H]. Analytic continuation over C, or the geometric
etale fundamental group in general, acts on the cosets G/H. This is the
monodromy action. It is faithful because M is the Galois closure.

The action must be recorded together with H. The deck group of the generic
cover is N_G(H)/H, equivalently the centralizer of G in S_n. It need not
equal G. Polynomial source symmetries covering target symmetries form yet
another object. These three meanings of "symmetry group" must be kept separate.

For a nonclosed ground field, geometric and arithmetic monodromy are
different groups, G_geom normal in G_arith <= S_n. Their quotient describes
the algebraic constant field of the Galois closure. The F11 example below
makes this distinction necessary, even with an explicit rational formula.

The absolute spectrum consists of actions realized by polynomial Keller
self-maps of affine space, allowing dimension to vary. The chart spectrum
allows smooth affine source and target varieties and a chosen volume form.
A field-versal spectrum additionally asks to realize all torsors of a
specified action. These are different problems.

## 2. Three distinct realization steps

A finite etale algebra in a special fiber does not prescribe the generic
monodromy group. A special target slice can have a smaller group than the
ambient cover. The universal pencil H(W)-sW+t has group S_n for every
characteristic-zero degree-n polynomial H; choosing a Chebyshev or a
monomial seed alone does not preserve its dihedral or cyclic slice group.
See the retained [universal-pencil proof](../../foundations/verified/UNIVERSAL_SYMMETRIC_MONODROMY.md).

There is, however, a general cover-preserving chart construction. If
f:V -> U is a finite etale map of smooth affine curves with coordinates x,t
on open subsets of the projective line, write t=f(x). On a chart where x
is a coordinate and f' is a unit,

\[
 (x,z)\longmapsto (f(x),z/f'(x))
 \tag{2.1}
\]

has determinant one and preserves the function-field extension after adding
one independent transcendental variable. A monic root presentation
P(t,x)=0 gives the analogous derivative-unit construction.

Such charts can have many nonconstant units. Since A^m has only constant
units, adjoining affine coordinates does not make these charts affine
spaces. Changing the boundary through affine modifications is a separate
problem. A candidate absolute realization must exhibit polynomial
coordinates everywhere, a nonzero constant full Jacobian, an exact generic
field identification, and actual affine source and target spaces.

## 3. Every affine divisor is met without ramification

**Lemma.** Let F:A^m -> A^m be Keller over k. For every prime divisor
D=V(q) in the target there is a prime divisor E in the source dominating D,
and its ramification index in L/K is one.

Indeed q(F) is nonconstant, since F is dominant. A nonconstant polynomial
over an algebraically closed field is not a unit and has a prime factor.
Every component E of V(q(F)) has dimension m-1. The restriction of a
quasi-finite morphism is quasi-finite, so its image has dimension m-1;
its image closure is D. At the generic points the local rings are DVRs.
The local map is etale, hence has ramification index one and separable
residue extension. This works for m=1 as well.

In fact q(F) is squarefree: every component just considered has index one.
Distinct target prime divisors have disjoint sets of source prime-divisor
factors, since such a factor cannot dominate two different divisors.
Consequently pullback on the free divisor groups is a primitive injection:
each target prime maps to a nonempty sum of distinct source primes, and
the supports of these sums are disjoint. Choosing one source prime over
each target prime gives a splitting. The cokernel is free.

One elementary consequence, useful for quadratic resolvents, is

\[
 b(F)\in k(X)^{*d}\quad\Longrightarrow\quad
 b\in k(Y)^{*d}\qquad(d\ge2).
 \tag{3.1}
\]

Write the divisor of b in target prime divisors and compare valuations on
one source prime over each. Every exponent must be divisible by d; the
remaining constant is a d-th power in algebraically closed k. This
argument needs neither a prescribed primitive polynomial nor a Jacobian
ansatz. It cannot be applied to a selected target slice instead of the
whole affine target.

## 4. The fixed-sheet inertia theorem

**KMFI1.** For every prime divisor D of the affine target, a geometric
inertia group I_D of M/K fixes at least one sheet of G/H. The normal
subgroup generated by all these inertia groups is G. Consequently

\[
 \boxed{\langle gHg^{-1}:g\in G\rangle=G.}
 \tag{4.1}
\]

Here "geometric inertia" means inertia after strict henselization at D;
it is cyclic in residue characteristic zero. The group is determined up
to conjugacy. The theorem concerns divisors of A^m, not arbitrary
divisors at infinity on a compactification.

To prove the first assertion, use E from Section 3. In the standard
inertia action on embeddings of L into M, the orbits have sizes equal to
ramification indices of geometric branches. Index one gives a fixed
sheet. Equivalently, for a suitable conjugate of H one has I_D <= gHg^-1.
An unramified branch may have a nontrivial residue-field extension before
strict henselization; this does not remove its inertia fixed points.

For the second assertion, let N be the normal subgroup generated by all
the I_D, and normalize Y in M^N. This normalization is finite, since Y
is excellent and M^N/K is finite. All codimension-one inertia is killed,
so the finite morphism is unramified in codimension one. Its source is
normal and its target regular. Purity of the branch locus makes it
finite etale everywhere. Affine space over an algebraically closed
characteristic-zero field has no nontrivial connected finite etale cover.
Thus M^N=K and N=G. Since every I_D lies in a conjugate of H, (4.1) follows.

The purity input is [Stacks, Lemma 58.21.4](https://stacks.math.columbia.edu/tag/0BMB).
For simple connectedness, descend the finite etale cover to a finitely
generated characteristic-zero field, embed it in C and use comparison
with the contractible analytic space C^m. Equivalence between algebraic
finite etale covers and finite analytic covers is the Riemann-existence
comparison of SGA 1, Expose XII, Theorem 5.1 and Corollary 5.2.
Isomorphism and geometric connectedness descend/extend with algebraically
closed constants. These are the only global inputs in this proof.

An equivalent formulation of (4.1) is that there is **no nontrivial
intermediate Galois field** K < E <= L with E/K Galois. Such a field
corresponds to a normal subgroup N of G containing H. Its normalization
over Y would be the finite etale cover just excluded. No polynomial
factorization through E needs to be assumed.

This obstruction also applies to an arithmetic Galois intermediate
field over a characteristic-zero k: finite etale covers of A^m_k are
constant covers after geometric triviality, whereas k is relatively
algebraically closed in k(x_1,...,x_m). We use the geometric formulation
for the action classifications below.

The assumptions matter. Kummer covers of G_m can have regular monodromy,
because an affine divisor has been deleted and the target is not simply
connected. In characteristic p, x -> x^p-x is finite etale on A^1 with
regular C_p monodromy. These are controls against dropping the target
or characteristic hypotheses.

## 5. Absolute exclusions beyond regular actions

**KMFI2.** The following are necessary consequences of KMFI1.

1. No nontrivial regular action occurs: H=1 cannot normally generate G.
   This recovers the classical Galois-case obstruction.
2. No nontrivial nilpotent group occurs in any faithful transitive action.
   For a finite nilpotent group every maximal subgroup is normal. A proper
   H lies in such a subgroup, so its normal closure is proper.
3. No proper normal block stabilizer is allowed. Under the
   block/intermediate-field correspondence, a normal block stabilizer
   would give a nontrivial Galois intermediate field. This does not
   exclude all imprimitive actions or algebraize arbitrary intermediate
   fields into polynomial factors.
4. H must surject onto every quotient in which its normal closure would
   otherwise be proper; in particular H -> G_ab is surjective. This last
   test alone is weaker than (4.1), for example for perfect groups.

For the natural degree-n action of
D_n=<r,s | r^n=s^2=1, srs=r^-1>, of order 2n, the point stabilizer is <s>.
Conjugating s by r gives r^2s, so

\[
 \langle\!\langle s\rangle\!\rangle
   =\langle s,r^2\rangle
   =\begin{cases}D_n& n\text{ odd},\\
                 \text{an index-two subgroup}&n\text{ even}.
     \end{cases}
 \tag{5.1}
\]

Thus **every even-degree natural dihedral action is impossible absolutely,
in every affine dimension**, irrespective of masks, blowdowns or
stabilization. For even n both reflection classes fail. In fact these and
the regular action exhaust faithful transitive dihedral actions: the
rotation part of any stabilizer is a normal subgroup contained in its
core, so a core-free stabilizer is trivial or a reflection. Therefore
D_n for even n has no absolute faithful transitive realization at all.

There is an independent elementary check of both exclusions. Every
quadratic subextension has the form K(sqrt(b)); equation (3.1) forbids
such a subextension inside L. More generally, over algebraically closed
k every cyclic prime-degree extension is Kummer, so (3.1) forbids it.
A finite nilpotent group has a normal maximal subgroup of prime index
containing any given proper stabilizer. The even dihedral action has
the explicit quadratic quotient (5.1). Thus these exclusions do not
depend on the global purity proof of KMFI1. On the even Dickson cover
itself, the identity P_(2m)(a,u)+2u^m=P_m(a,u)^2 displays the quadratic
resolvent directly.

Odd natural dihedral actions pass (4.1); this is not a construction.
D_3=S_3 is already realized by the retained symmetric family. D_5 and the
larger odd dihedral actions still require absolute realizations or new
obstructions. The Dickson derivative-unit charts remain valid in both
parities. Their even-degree absolute search is superseded by (5.1).

The same computation gives a useful affine-group test. For an abelian
group V with a faithful complement A acting on it, in the transitive
action of V semidirect A on V,

\[
 \langle\!\langle A\rangle\!\rangle=[V,A]\rtimes A.
 \tag{5.2}
\]

Indeed a normal subgroup containing A contains every commutator [v,a],
and quotienting by [V,A] makes A normal. Hence the necessary condition is
V/[V,A]=0. It excludes generalized dihedral actions when V/2V is nonzero.
For F_20=C_5 semidirect C_4 this quotient vanishes; its natural action
survives this test. These universal group statements have group-theoretic
proofs; the bounded finite replay illustrates and checks them.

## 6. The PSL2(11) local sieve

**KMFI3.** For the natural degree-twelve action of PSL2(F_11), the complete
element-order and cycle dictionary is

| order | number of elements | cycles | allowed nontrivial affine-divisor inertia? |
|---:|---:|---|---|
| 1 | 1 | 1^12 | unramified |
| 2 | 55 | 2^6 | no |
| 3 | 110 | 3^4 | no |
| 5 | 264 | 5^2 1^2 | yes |
| 6 | 110 | 6^2 | no |
| 11 | 120 | 11 1 | yes |

The dictionary follows either by the split/nonsplit/unipotent
classification or by enumerating SL2(F_11)/{+/-I}; the replay does both
the projective action and its fixed-point test. In general the natural
PSL2(q) action permits only split semisimple and unipotent nontrivial
divisor inertia. Nonsplit semisimple inertia is a derangement.

The modular (2,3,11) cover on X_0(11) has two derangement branch values.
A Mobius rechoice of the target point at infinity can move only one of
them out of the affine line. Thus that product presentation cannot
satisfy the affine-divisor test. This is a statement about that cover and
its affine branch placement, not an exclusion of the natural action.

Suppose now that a degree-twelve cover of P^1 with this natural group
has genus zero and every *finite* branch inertia has a fixed sheet.
Finite inertia can have only orders 5 and 11, of indices 8 and 10.
The index at infinity belongs to {0,6,8,10}. Riemann--Hurwitz gives

\[
 8a+10b+\operatorname{ind}(\sigma_\infty)=22.
 \tag{6.1}
\]

For nonnegative integers a,b the unique solution is a=2, b=0,
ind(sigma_infty)=6. Thus the unique possible genus-zero passport with
this branch-placement condition is

\[
 5^2 1^2\mid5^2 1^2\mid2^6.
 \tag{6.2}
\]

This is a necessary passport classification under the genus-zero
hypothesis, not a uniqueness theorem for dessins or a genus-zero theorem
for arbitrary inverse curves of Keller maps.

The passport really occurs. In PSL2(F_11), matrices modulo scalar sign

\[
 A=\begin{pmatrix}0&1\\10&3\end{pmatrix},\quad
 B=\begin{pmatrix}0&1\\10&8\end{pmatrix},\quad
 C=AB=\begin{pmatrix}1&3\\3&10\end{pmatrix}
 \tag{6.3}
\]

have orders 5,5,2, satisfy ABC=1 projectively and generate all 660
elements. Their natural cycle indices are 8,8,6, giving genus zero.
The checker certifies generation rather than inferring it from a
passport. Riemann existence gives a rational cover over C.

The passport occurs in van Hoeij--Vidunas's F11 example. Its explicit
rational formula and arithmetic coefficient field are retained in
[model.py](model.py) and [model.json](model.json). The
[primary permutation data](sources/Dessins_and_ModuliFields.txt) and
[finite checker](verify.py) compare the literal dessins with the projective
group action. This rational cover does not itself give a polynomial
self-map of affine space.

## Current scope

The fixed-sheet inertia theorem, its affine-group exclusions and the
degree-twelve local sieve above are the retained results. A complete
classification of monodromy groups of polynomial Keller maps remains
open. The [canonical receiver and profile results](../universal-keller-structure/PROOF.md)
give complementary positive constructions under their stated hypotheses.

The geometric arguments are written proofs using branch purity and the
characteristic-zero covering theorem cited above. Finite computations
check the explicit group actions, formulas and passports; they do not
supply an independent verification of those general geometric proofs.
No formal verification or external review is claimed.
