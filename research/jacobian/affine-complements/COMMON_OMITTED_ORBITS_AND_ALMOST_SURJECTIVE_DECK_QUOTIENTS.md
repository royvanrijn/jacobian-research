# Common omitted orbits give an almost-surjective deck quotient

PF2RULED136. This is an all-degree reduction using the ACTUAL whole plane
chart and target, rather than a rational replacement of their fields.
It does not prove full finite-order multiplication, trace regularity or
Wright/HC4/JC2. Finite point gaps and the mixed-fibre profile remain OPEN.

For any subgroup of the generic deck group, remove exactly the omitted
curve orbits missed by every translated plane chart. The resulting actual
affine quotient is almost covered by the original plane. Each removed
orbit costs one unit of quotient Euler characteristic. Noncyclic groups
leave quotient Euler at least3, giving a sharper degree budget and an
almost-surjectivity condition at the unresolved index-five frontier.

## 1. Actual setup and the curves missed by every translated chart

Retain an ACTUAL complex polynomial plane Keller map, with original
X=A2, field L/F of degree d, entire etale completion U=U_N, and all
omitted primes. By HC4MYLA0,

    U minus X = disjoint union of b0 smooth A1 curves E_i,
    pi1(U)=1, units(U)=C*, chi(U)=1+b0.

Let G<=Aut_F(L) be ANY subgroup, g=|G|. Its action on U is regular and
FREE. Put T=U/G, q:U->T, c=chi(T)=(1+b0)/g. The quotient is smooth
affine and etale over the ORIGINAL target A2, with its exact coframe.
We do not assume G preserves X or permutes the original boundary basis.

Consider the closed invariant set

    Z=intersection_(gamma in G) (U minus gamma(X)).

Its curve components are among the original E_i, because the identity
is included. The union D of those components is invariant. A stabilizer
of an E_i in D is trivial: any finite-order automorphism of A1 has a
fixed point, while even an element acting trivially on that A1 would
fix points of U. Either case contradicts the free action unless the
element is the identity. Thus D is a union of r FREE curve orbits,
each containing g pairwise disjoint A1s:

    #components(D)=g*r.

The quotient C=D/G is consequently a disjoint union of r embedded smooth
A1 curves in T. Its full inverse image is D. Equivalently, C is EXACTLY
the divisorial part of T minus q(X): a point of T is missed by q(X)
precisely when its entire G orbit upstairs lies outside X.

## 2. Delete the divisorial gaps, retaining the original plane and target

Define

    T'=T minus C,   U'=q^-1(T')=U minus D.

Both are affine: on the smooth affine surfaces the deleted curve unions
are effective Cartier divisors, whose complements are affine open
immersions. The finite etale quotient restricts to U'->T' with group G.
The WHOLE original X lies in U'. Its complement there is still a disjoint
union of b0-g*r original unramified A1s; no ramified prime or relative
defect has been discarded from the ORIGINAL finite normalization.

Moving loops off that divisor complement gives pi1(U')=1, just as for U.
Hence U' is the universal finite etale cover of T' and pi1(T')=G. Its
exact volume, scalar units and negative log Kodaira dimension are retained,
the latter from domination by X. The original target functions restrict
to an actual ETALE map T'->A2, and the composite X->T'->A2 is still the
SAME original polynomial Keller map.

The map

    X=A2 -> T'

is etale, has Galois generic field extension L/L^G of degree g, and is
ALMOST SURJECTIVE: T' minus its image is a finite set. Indeed every
divisorial missed component was deleted, and etale images are open.
Thus it is a Galois affine pseudo-covering in Miyanishi's precise sense.
Surjectivity is NOT claimed. Upstairs the union of the gamma(X) misses
only a finite G-stable set in U'; since the action is free, its cardinality
is g times the number of missed points downstairs.

The open T' is the LARGEST open subset of T over which this original X
is almost surjective: including any nonempty part of a deleted irreducible
curve in an open would include its generic point, again missing a divisor.
This construction is unique; no ruling or chart choice was added.

T' is an auxiliary open quotient, NOT generally the entire etale locus
of the finite intermediate in L^G. Its function field and target map
remain actual. The full finite normalization and its original discriminant
are unchanged, so this deletion is not a trace-discriminant vanishing
argument and does not replace that normalization by a selected open.

## 3. Exact Euler and degree budgets for the removed orbits

Every removed curve has compactly supported Euler characteristic1.
Consequently, with c'=chi(T'),

    chi(U')=1+b0-g*r,   c'=c-r,
    rank Pic(T')=c-r-1.                               (1)

The rank statement follows by the same whole-plane Gysin and finite-cover
transfer argument as HC4MYLA0: U' has free boundary Picard group and only
H0/H2 rational cohomology. In particular c'>=1. No integral permutation
structure is assumed for the remaining boundary basis.

If G is NONCYCLIC, the low-Euler exact-volume lemma PF132 applies to
T' and forces c'>=3. Hence

    r<=c-3.                                          (2)

Use PF133's ACTUAL deck-normalized ledger, b=d/g,

    b=c+S+J,
    S=sum_C(w_C/g)>=2,   J=TOTAL Eta/g>=0,

where S,J are integers. The lower bound on S is for noncyclic G: then
b>1 and the genuine fixed-field intermediate has a visible branch.
Combining (1)--(2) gives the sharper exact ledger and bound

    b=c'+r+S+J,   b>=r+5.                            (3)

All original branch supports, invisible inertia, movement in fixed blocks
and excess remain in S,J. The new integer r counts common OMITTED
UNRAMIFIED curve orbits, not ramification or residual support.

Equality b=r+5 forces c'=3, S=2 and J=0. As in PF133, there is exactly
one original branch image, its quotient inertia is a transposition,
w_C=2g, and residual defect is zero. Normal generation and transitivity
give a natural S_b quotient. Eta=0 and PF130 give smooth original and
intermediate normalizations, residue-one critical curves, immersive image
normalization, and

    T=V (entire intermediate etale locus),
    chi(V)=r+3, U=h^-1(V), U->V finite etale,
    b0=g*(r+3)-1, chi(T')=3.                         (4)

This is a necessary profile, not an exclusion or construction. In
particular the first index-five case has r=0, so T'=T=V: its ORIGINAL
plane map A2->V is already a Galois affine pseudo-covering, with only
possible finite missed points. PF134 still leaves one mixed two-component
bad ruling fibre there. Its Picard rank is2, so it is NOT a Platonic
A1-fibre space with every fibre irreducible and Picard rank1. The published
irreducible-fibre context does not close this profile.

## 4. A prior whole-plane cyclic control makes the deletion nonvacuous

Reuse the known Danielewski surface

    W_k={x^2*y=z^k-1}, k=g*h, g>=2,

with free cyclic subgroup C_g acting by
(x,y,z)->(epsilon*x,epsilon^-2*y,epsilon*z). The original plane chart
at the root1 omits exactly D_beta={x=0,z=beta}, beta!=1. Its g translates
retain the g roots in the orbit of1. All other root components are
common omitted curves. Thus

    #common curves=k-g, r=h-1,
    chi(W_k/C_g)=h, chi(T')=1.

The translated charts cover U' in this control, so its plane map to T'
is even surjective. The noncyclic inequality c'>=3 does not apply to C_g.
This reuses the retained atlas and free action; their checkers are not
rerun. No etale map from this control quotient to the PLANE target has
been supplied, so it is not an actual Keller example. It does show why
the divisorial deletion is needed even with a whole plane, free finite
symmetry and exact unit volume.

## 5. Source boundary and the remaining obligation

Miyanishi's primary 2005 report, printed1110/PDF38, supplies the exact
definition of affine pseudo-covering and the finite Galois-normalization
context. The complete relevant statements were read; the construction
above is proved directly and imports no new classification. The full
2005 Journal of Algebra article was not obtained from the blocked publisher;
its abstract is not used as a theorem eliminating general mixed fibres.

Shared HC4MYMO0's complete-ruling plane chart may omit a WHOLE fibre and
need not be almost surjective. The almost-surjectivity here instead uses
the G-translates of the ACTUAL original chart and the explicit deletion
of their common omitted curves. Neither assertion is substituted for
the other or credited as a new classification of pseudo-coverings.

The source-base-change audit HC4MYMA0 retains the residual discriminant;
splitting off a known inverse branch cannot remove it. Likewise T' still
has an actual finite-target/coframe problem, and almost surjectivity does
not make its curve fibres irreducible, erase finite missed points, or
force trace-pairing regularity. The next genuine gate is this almost-
surjective mixed-fibre case, or a global finite-order/trace/energy theorem.
No claim of full Wright/HC4/JC2 closure, formal verification, independent
replay, external review or novelty is made. The written results are all-
degree; finite orbit and Euler controls support rather than replace them.
