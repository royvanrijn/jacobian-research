# Global oriented normal forms and the homological Seifert barrier

Theorems HC4MYKB0/1/2, over C. There is a complete polynomial source
normal form for the following actual class: potentials linear in two
fibre variables, with a finite simple coefficient zero set and one
common nonzero full critical Hessian determinant. Within this class
the number and multiset of critical values determine the full source
orbit.

This gives explicit source-alignment gates for K2 and S_rho, separate
from existence of two- or three-point Keller fibres. A primary-source
audit also closes the homological Seifert-form route without assuming
regularity at infinity. Neither gate for the original models is solved.

## 1. The global polynomial normal form (HC4MYKB0)

Let b=(b1,b2) be actual polynomial base coordinates, and let

    h=d(b)+U f1(b)+V f2(b), R=C[b1,b2,U,V].

Assume I=(f1,f2) is the radical ideal of N>=1 distinct points b_i.
Then j_i=Jac_b(f1,f2)(b_i) is nonzero. The full critical scheme of h
has exactly one point over each b_i, with value c_i=d(b_i), and
Hessian determinant j_i^2. Suppose these squares are one common m!=0.

Define, for all finite N at once,

    C_N(z)=product_(i=1)^N(z-i),
    a_N(z) C_N'(z)=1 modulo C_N(z), deg a_N<N,
    D_c(i)=c_i, deg D_c<N,
    R_(N,c)=a_N(z) X P+C_N(z) Q+D_c(z).                 (1)

The inverse a_N exists because C_N is squarefree. Then h is
polynomially source-equivalent to (1), with actual polynomial inverses.

Here is the construction. Use the finite-point plane automorphism
from JV2 to put the b_i at (X,z)=(0,i). Its constant Jacobian multiplies
all j_i by the same scalar, so their squared equality persists.
Write these transformed values again as j_i. Interpolate a(z) with

    a(i)=j_i/C_N'(i).

The pair (a(z)X,C_N(z)) generates (X,C_N). Its point Jacobians are
exactly j_i. JV2's determinant-repair lemma therefore relates it to
the transformed f by a GL2(C[X,z]) matrix with determinant one.
The corresponding transpose or inverse-transpose fibre map, as
dictated by the coefficient identity, converts h to

    d_0(X,z)+a(z)X P+C_N(z) Q.

Since d_0-D_c vanishes on every base point, it lies in (aX,C_N).
A polynomial translation of P,Q removes this difference. Every step
so far has a polynomial inverse and fixes the critical values.

It remains to equalize the unsquared signs. Choose j with j^2=m,
using the transformed m, and interpolate p(z) with p(i)=0 or 1 so

    (1-2p(i)) j_i=j.

Put e=1-2p and

    M(p)=[[1-p^2,2p-p^3],[-p,1-p^2]], det M=1.

For (X_old,P_old)^t=M(p)(X,P)^t, the product difference

    q=X_old P_old-e X P

is divisible by p(p-1), hence by C_N. Replace the old Q by

    Q_old=Q-a(z)q/C_N(z).

This is an actual polynomial source automorphism retaining z. Its
inverse first applies M^-1 to X_old,P_old, then adds a q/C_N to Q_old.
It changes the potential to

    a(z)e(z)X P+C_N(z)Q+D_c(z).

All its coefficient point Jacobians are now j. Comparing this row
with (a_N X,C_N), whose point Jacobians are one, JV2 gives a global
invertible fibre matrix and hence (1). This is the weighted version
of JX0's explicit root sign twist. It changes the embedded base X
when needed; unsquared fixed-base orientation is not assumed invariant
under arbitrary source maps.

Conversely N and the multiset {c_i}, with repetitions, are preserved
by every exact polynomial source equivalence. Any permutation can be
chosen in the finite-point interpolation. Thus two members of the
stated class are source-equivalent if and only if these data agree.
The common Hessian scalar m is absorbed by the constant source
Jacobian; no value scaling or translation is included in this
classification.

For N=1, (1) is a nondegenerate quadratic plus its critical value.
If I is the unit ideal instead, row completion makes h an actual
polynomial coordinate; this is the separate zero-critical-point
endpoint. Nonradical coefficient zero schemes and arbitrary
four-variable polynomials are outside this normal-form theorem.

## 2. Exact existence and alignment gates (HC4MYKB1)

Let E_N mean existence of a plane Keller map with one actual fibre
of exactly N points, not a statement about generic degree. For
every N>=2 and every critical-value list c, the following are equivalent:

- R_(N,c) has a unit-Hessian representative under arbitrary polynomial
  source automorphisms;
- it has one through a base-preserving affine fibre automorphism;
- E_N holds.

The forward implication is HG0/Tang: its N finite Morse critical
points survive source change and force a nonautomorphic gradient;
the coefficient Keller zero fibre has exactly N points. Conversely
an E_N pair can be moved to the prescribed N base points. Its ideal
is (X,C_N), and its point Jacobians have the same constant ratio to
those of (a_N X,C_N). JV2 gives the required invertible polynomial
fibre matrix. The resulting potential D_c+P K1+Q K2 has literal
Hessian determinant Jac(K)^2, for arbitrary c.

Thus the normal form in section 1 has unit existence exactly E_N,
in every source degree. This does not construct or exclude any
E_N. For N=1 the normal form already has an invertible unit gradient;
its *universal gradient-inversion* question remains the old singleton
endpoint. Unit existence must not be substituted for that endpoint.

More generally, if a four-variable polynomial f has exactly N>=2
finite Morse critical points, one common nonzero critical Hessian
determinant, and critical-value multiset c, then

    f has a full-source unit member
       iff E_N and f is source-equivalent to R_(N,c).   (2)

Necessity uses Tang and section 1 on its actual Keller cotangent
representation. Sufficiency composes the source equivalence with
the conditional fibre map above. The two conditions in (2) are
separate; matching finite critical data alone is not an alignment.

For the actual affine frontier this gives two concrete statements.

First, K2 has a unit member if and only if E2 holds AND K2 is
polynomially source-equivalent to h1. Indeed both have critical values
{-1,1}; h1 already belongs to section 1's class, so its normal form
is R_(2,(-1,1)). In particular K2 unit existence would *force* a
full K2/h1 source equivalence. Neither E2 nor this equivalence is
proved. K2 remains redundant in the universal test list by JX2.

Second, write the three-point target entirely explicitly:

    C_3=(z-1)(z-2)(z-3),
    a_3=(3z^2-12z+10)/2,
    D_rho=[(rho+3)z^2-(3rho+13)z+2rho+12]/2,
    R_rho=a_3(z)X P+C_3(z)Q+D_rho(z).                  (3)

Here a_3 C_3'=1 modulo C_3 and D_rho has values 1,-1,rho at
1,2,3. For every rho, including rho=+-1,

    S_rho unit existence
       iff E3 and S_rho is source-equivalent to R_rho. (4)

This isolates the actual source/action alignment from the global
three-point Keller question. No converse from E3 alone is asserted.
For the universal family requirement, (4) is equivalently

    no unit S_rho for every rho
       iff not E3 OR no rho admits this alignment.     (5)

The complete admitted affine frontier is therefore the two original
endpoint inversion requirements, absence of E2, and the disjunction
(5). These reductions do not realize arbitrary auxiliary N-point
models as original Meng--Yang graphs or cover arbitrary HC4 embeddings.

## 3. The homological Seifert form is already forced (HC4MYKB2)

We inspected the full primary
[Neumann--Norbury, math/9910054v1](https://arxiv.org/html/math/9910054v1),
including the statements and proofs of Theorems 2.2,2.3,3.2,3.3.
Their setup allows isolated finite singularities with unrestricted
singularities at infinity. In particular their regular-link Seifert
construction is available here; no unproved tameness assumption
is inserted.

The precise imported facts are: the global positive homology splits
into the variation images of all finite and infinite Milnor pieces;
these images inject; the regular-link Seifert form is block triangular,
with local blocks and zero cross terms at the same value; its cross
terms at different values are tied to off-diagonal homological
monodromy by Theorem 3.3.

All normalized models and required Keller cotangent models have
reduced integral homology only in degree three, of rank N. By KA,
the N finite Morse vanishing cycles already form an integral basis,
and each local monodromy on H3 is the identity. The same torsor/link
argument applies directly to (1), using only its simple coefficient
zeros, not a global Keller condition.

The imported injective direct sum leaves no positive-degree relative
homology for any Milnor piece at infinity: in degree three the finite
pieces exhaust the lattice; in every other positive degree the
entire target homology is zero. This does not prove that those
noncompact pieces are empty, contractible or topologically regular.

In degree three the graded transpose in the source is the ordinary
transpose. Each finite Morse piece has the same rank-one Seifert
form epsilon, epsilon=+1 or -1 according to the single global sign
convention. Finite variation is an isomorphism, so the transported
variation generators are the primitive local cycles already used
in KA, up to sign.

For different critical values, Theorem 3.3 gives L_i h_ij=-L_ij
or its transpose. Since h_ij=0, all such Seifert cross terms vanish.
At equal values they vanish by Theorem 3.2. Hence

    L=epsilon I_N, intersection form L-L^t=0.           (6)

This holds also at every value collision and for the explicit
universal targets (1),(3). N=0 is the empty form.

For S_rho, in KA's positive puncture basis rather than its local
vanishing basis, this is explicitly

    L_puncture=epsilon [[2,1,-1],[1,2,-1],[-1,-1,1]].    (7)

Indeed the vanishing basis has row matrix
M=[[1,0,1],[0,1,1],[0,0,1]], so L_puncture=M^-1 L M^-t.
It is integrally congruent to epsilon I_3 for every rho.

The form can be defined by a compact-core push-off inside the
regular-link sphere, as in the source's proof. Exact polynomial
source equivalences transport that oriented construction; changing
an integral basis gives congruence. In any event the actual computed
forms on both sides of the necessary comparison are already equal.
Neither a finer homological Seifert calculation nor a nonzero
positive-degree infinity-variation block can supply an obstruction.

The diagonal calculation does not identify the regular link itself,
the embeddings of Milnor pieces, nonhomological thimble data, or
the topology of homologically invisible pieces at infinity.
The source explicitly distinguishes this global embedding data
from its homological constraints. Actual critical gluing, a proof
or disproof of the alignments (2)--(4), globally Keller coefficient
existence, the old endpoints and unrestricted HC4/JC2 remain open.

## Reproducibility and boundary

[The source context](sources/global-normal-forms-seifert-context.json)
pins the new primary, its precise proved imports and checked
hypotheses, and retained JV,JX,JQ,KA inputs.
[The checker](verify_global_normal_forms_seifert.py) checks the weighted
sign-twist inverse, coefficient ideal and orientation identities,
the explicit three-point target, its entire critical scheme and
Hessian controls, and the block-matrix consequences.
[Receipt](receipt-global-normal-forms-seifert.json).
The all-N normal-form construction and imported topological theorems
are not inferred from finite symbolic samples. Independent replay,
formal verification and external review are not claimed.

