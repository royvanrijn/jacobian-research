# Nonspherical reductive pairs and the two exceptional quotients

15 November 2026. Internal written deduction from Long's
[pinned spherical-pair draft](https://github.com/long-mathematics/mathieu-second-conjecture-spherical-pairs/blob/e9edac262b4b64eb68feb0e6e80496cefe4a1dee/mathieu_second_conjecture_spherical_pairs.tex),
using only its parts that do not assume sphericity, and the classical finite
subgroup classification in
[Flicker, printed page 215](https://doi.org/10.1007/s00229-021-01274-x).
The local [source audit](exceptional_sl2_source_audit.json) binds the sources.
The [checker](verify_exceptional_sl2.py) certifies the finite groups and
all-degree rational character identities below. The general reduction is
a written proof; independent specialist review and formal verification
remain outstanding.

Sphericity can be removed completely for connected stabilizers and for a
larger class controlled by their finite component groups. For arbitrary
stabilizers, two explicit three-dimensional homogeneous spaces suffice to
test the proposed classification. Their moment problem remains open here.

**Later result.** The [Haar-shear theorem](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md)
now closes moment-nullcone detection for all reductive homogeneous pairs,
without sphericity. Its witnesses are eventually marked-null. The two
exceptional global Mathieu--Zhao questions remain open.

## 1. The unconditional extension

Let G be connected complex reductive and H a closed reductive subgroup,
possibly disconnected. Put

    A=C[G/H],     D=[G,G],     F=H/H^0.

Let h be normalized Reynolds, Z the elements with every positive pure
moment zero, and N the G-nullcone in the coordinate algebra, as in the
[component theorem](SPHERICAL_PAIRS_AND_COMPONENT_NULLCONES.md).

**Theorem 1.** Assume F has no quotient isomorphic to A_4, S_4 or A_5.
Then, with no sphericity assumption, the following are equivalent:

* D is contained in H;
* ker h is Mathieu--Zhao;
* Z=N.

If D is not contained in H, there are fixed regular P,Q with

    h(P^m)=0,             h(Q P^m)>0        for every m>=1.       (1)

In particular the conclusion holds when H is connected, when F is abelian,
or whenever 12 does not divide |F|. The numerical sufficient condition
uses the orders 12,24,60 of the projective exceptional groups, rather than
24,48,120 of their binary lifts.

For example G=SL_2^r and H=diagonal SL_2, with r>=4, has connected
stabilizer and is nonspherical: dim(G/H)=3r-3>2r=dim(Borel).
It nevertheless has the fixed witnesses (1).

The hypothesis is sufficient. An exceptional quotient of F is not itself
a counterexample to the classification: many such stabilizers are already
covered by the witness constructions below.

## 2. Reduction without a spherical projection hypothesis

Take a finite central isogeny

    pi : Gtilde=Z_0 times product_j S_j -> G,
    Htilde=pi^(-1)(H),

where the S_j are simply connected simple groups. This does not change
the homogeneous space, normalized functional or coordinate-ring group
orbits.

If an active factor S_j has rank at least two, Long's universal-root-doublet
lemma and line--doublet proposition give (1). Neither step uses sphericity.
The detection lemma obtains an irreducible module and Htilde-fixed vector
from evaluation at the base point. The actual defining root summand,
equivariant orthogonal projection and radial identity then give the
nonzero positive marked moments.

Otherwise every higher-rank factor is contained in Htilde and acts trivially.
Quotienting them out leaves an isomorphic homogeneous space with acting group

    G'=Z_0 times SL_2^r

and reductive stabilizer H'. Suppose its derived group is not contained
in H'. Each factor projection H_i=p_i(H') is a closed reductive subgroup
of SL_2.

The proper positive-dimensional possibilities are a torus or its normalizer.
For finite H_i, the finite subgroup classification gives cyclic, binary
dihedral, binary tetrahedral 2T, binary octahedral 2O, or binary icosahedral 2I.
Include the trivial finite subgroup among the cyclic cases.

Cyclic and binary dihedral groups are contained, after conjugation, in the
normalizer N(T). Thus whenever a proper projection is not exceptional,
there is a genuine equivariant surjection

    G'/H' -> SL_2/N(T).

The regular Weyl-even circuit in Long's draft pulls back to give (1).
This includes nonspherical spaces such as SL_2/Gamma for every finite
cyclic or binary dihedral Gamma; no Borel orbit dimension bound is used.

If every projection is surjective onto SL_2, Long's subdirect-stabilizer
argument still applies without sphericity. Its connected derived subgroup
is a product of diagonal SL_2 blocks. Singleton blocks would put the whole
derived group in H'. A nonsingleton block supplies an H'-stable alternating
line in the tensor product of two standard modules and a defining doublet.
The line--doublet proposition again gives (1), including disconnected H'.

The only remaining branch has a finite exceptional factor projection.
This proves the exhaustive reduction. The alternatives need not be disjoint:
a pair can have an exceptional projection and another projection that
already supplies witnesses.

## 3. Why projective exceptional groups must be quotients of F

Return to the original central cover Gtilde -> G and retain the exceptional
SL_2 factor of the preceding argument. The composite

    Gtilde -> SL_2 -> PSL_2

kills ker(pi), since the kernel is central and the SL_2 center is {+1,-1}.
It therefore descends to a surjective algebraic homomorphism

    alpha : G -> PSL_2.                                    (2)

The image alpha(H) is respectively A_4, S_4 or A_5. It is finite, so the
connected group H^0 maps trivially to it. Consequently (2) induces a
surjection from the original component group F onto alpha(H).

This tracks the component group through the central cover. Counting
components of Htilde instead of H can give a false hypothesis; the
central kernel has been removed by projectivization.

If F has none of these quotients, the exceptional branch is impossible.
Every non-toral pair therefore has witnesses (1). Hilbert--Mumford and
invariance make all unstable elements eventually marked-null, so these
witnesses are outside N and violate the Mathieu property.

Conversely D contained in H makes G/H a torus. Duistermaat--van der Kallen
and strict separation give Z=N and the Mathieu property. This completes
the proof of Theorem 1.

## 4. Exactly two model questions

Choose binary polyhedral subgroups 2T,2O,2I of SL_2. The binary tetrahedral
group is contained in 2O (and also in 2I). Thus the projective A_4 is
contained in an S_4 subgroup of PSL_2. Put

    Y_O = SL_2/2O = PSL_2/S_4,
    Y_I = SL_2/2I = PSL_2/A_5.

**Theorem 2.** For every connected reductive pair (G,H) with D not contained
in H, either the proof of Section 2 supplies (1), or there is a surjective
equivariant morphism G/H -> Y_O or G/H -> Y_I, with a surjective acting-group
homomorphism G -> PSL_2.

Indeed, use alpha from (2). If alpha(H)=A_4 or S_4, enlarge it to S_4;
if alpha(H)=A_5, keep it. This gives the asserted quotient map on the
original G/H, not just on a covering space.

For any such quotient q:X -> Y,

    h_X(q^*f)=h_Y(f).                                     (3)

Composition with q^* is normalized and invariant, so uniqueness proves
(3). Pullback is injective. On each finite-dimensional submodule the
acting-group orbits agree through the surjective group homomorphism and
the linear embedding q^*. Therefore membership in the respective nullcones
is also preserved and reflected for pulled-back elements.

It follows that the following two reductions are exact:

* The universal connected reductive-pair statement
  "ker Reynolds is Mathieu--Zhao iff G_der is contained in H" is equivalent
  to failure of the Mathieu property on **both** Y_O and Y_I.
* The universal connected reductive-pair statement
  "Z=N iff G_der is contained in H" is equivalent to Z differing from N
  on **both** Y_O and Y_I.

Necessity follows by specializing to the two non-toral model pairs.
For sufficiency, use (3) and the exhaustive reduction. Failure of the
Mathieu property on either model supplies one pure-null element and a fixed
marker nonzero at infinitely many indices; it need not supply strict
positivity at every index. Such a marked witness is enough to violate
nullcone equality, by Hilbert--Mumford. Thus failures of the Mathieu property
on both models would prove both universal classifications. Conversely,
failure of nullcone equality alone is not asserted to imply failure of the
Mathieu property.

If both models admit the stronger witnesses (1), that stronger all-positive
marked conclusion also transfers to every non-toral connected reductive pair.

Both models have dimension three and finite stabilizer, so their Borel
orbits have dimension at most two. They are genuinely outside the spherical
classification. Neither model is solved by this reduction.

## 5. Exact character series and degree bounds

Let V_n=Sym^n(C^2), let Gamma be a finite subgroup of SU_2, and set

    m_n=dim(V_n^Gamma).

Finite-group averaging of characters gives the rational identity

    sum_(n>=0) m_n t^n
       = (1/|Gamma|) sum_(g in Gamma) 1/(1-tr(g)t+t^2).       (4)

The denominator is the symmetric-power character generating function for
a determinant-one two-dimensional matrix. Thus (4) is an all-degree
identity, not a fit to computed coefficients.

The checker generates all group elements exactly as unit quaternions in
Q(sqrt(2)) or Q(sqrt(5)). It verifies group orders, inverses, projective
element orders, and the inclusions of 2T in 2O and 2I. Pairing conjugate
quadratic traces evaluates the right side of (4) over Q(t). Exact polynomial
cancellation proves the following formulas.

| Gamma | Character series | First n>0 with m_n>0 |
|---|---|---:|
| 2T | (1-t^24)/((1-t^6)(1-t^8)(1-t^12)) | 6 |
| 2O | (1-t^36)/((1-t^8)(1-t^12)(1-t^18)) | 8 |
| 2I | (1-t^60)/((1-t^12)(1-t^20)(1-t^30)) | 12 |

These familiar invariant-series identities are used as certificates, with
no novelty claim. A separate symmetric-character recurrence checks the
coefficients through n=32.

Algebraic Peter--Weyl identifies the left SL_2 types in C[SL_2/Gamma] as
V_n^* with multiplicity m_n. Therefore no nonconstant function on Y_O has
left highest-weight degree below 8, and none on Y_I has degree below 12.
The degree is the maximum n among its nonzero left irreducible components.

The same lower bounds apply to the least total degree of a representative
polynomial in the four matrix entries of SL_2. A polynomial of total degree
at most d has left types of highest weight at most d; averaging over Gamma
preserves this degree bound. Conversely a nonzero Gamma-invariant binary
form of degree 8 or 12 supplies a nonconstant matrix coefficient of that
degree. Thus 8 and 12 are the exact minimum degrees of nonconstant functions
in these two coordinate rings.

Any nonzero pure-null counterexample P must be nonconstant, so these are
rigorous lower bounds on its possible degree. Any marker giving a nonzero
marked moment must also be nonconstant. The bounds do not assert that
degree 8 or degree 12 already contains a counterexample. Highest-weight
vectors at the first allowed degrees are pure-null but unstable and hence
do not solve either exceptional problem.

## 6. Two barriers to reusing the existing constructions

For G=SL_2 an irreducible module is V_n. A nontrivial algebraic homomorphism
SL_2 -> SL_2 is an automorphism, so a defining two-dimensional summand in
an irreducible module forces n=1. The exceptional groups contain the
quaternion subgroup Q_8. Its standard generators i,j have commutator -1.
They cannot preserve a common complex line: their scalar restrictions
would commute on it, whereas -1 acts as -Id. Thus the line--doublet
criterion cannot by itself handle either exceptional model.

There is a separate obstruction to powers of the normalizer circuit.
For any exceptional Gamma and any maximal torus T,

    C[SL_2]^(right N(T)) intersect C[SL_2]^(right Gamma) = C.   (5)

To prove (5), let L be the Zariski closure of the group generated by N(T)
and Gamma. Its Lie algebra contains Lie(T), is T-stable, and is invariant
under the Weyl element. The root decomposition forces Lie(L) to be either
Lie(T) or sl_2. In the first case L^0=T, so L normalizes T, which would put
Gamma in N(T). This is impossible for projective A_4,S_4,A_5: every finite
subgroup of a torus normalizer is cyclic or dihedral projectively. Hence
L=SL_2. A regular function invariant under both groups is constant.

If a positive power of a nonconstant right-N(T)-invariant function descended
to SL_2/Gamma, that power would contradict (5). Consequently simply raising
the verified normalizer witness to a power cannot solve either exceptional
case. Averaging or taking a norm over Gamma creates different functions;
their pure moment identities require a new proof.

## 7. Averaging the normalizer witness fails in two different ways

Let P_N be the regular normalizer circuit from the earlier note, and average
it on the right over Gamma:

    E_Gamma(P_N)(g) = (1/|Gamma|) sum_(gamma in Gamma) P_N(g gamma).

This average is regular and right-Gamma-invariant. Write U=u^2 and S=tau^2
as before. The circuit has the exact expansion

    P_N = V + (1-6S+5S^2)
          + U(-6S+10S^2) + U^2(-2S+10S^2)
          + 5U^3 S^2 + U^4 S^2.                         (6)

Under the substitution g -> diag(t,t^(-1))g, the six summands have
exponents -4,0,4,8,12,16 respectively. These are substitution exponents;
the representation on functions uses the inverse substitution.

The -4 part V belongs to the right irreducible type V_4. The zero part
1-6tau^2+5tau^4 has right types contained in V_0,V_4,V_8, and its Haar mean
is zero, so its V_0 part vanishes. One way to see the degree bound is to
view it as a degree-four even polynomial on the adjoint sphere: its
spherical harmonics have degrees 0,2,4, corresponding to SL_2 types
V_0,V_4,V_8. Thus the character gaps in Section 5 decide all nonpositive
parts of the averaged circuit.

For **Gamma=2I**, neither V_4 nor V_8 has a fixed vector. The entire
nonpositive part of E_(2I)(P_N) vanishes. The checker verifies it is
nonzero by its exact restriction

    E_(2I)(P_N)(diag(t,t^(-1)))
      = -(sqrt(5)/24)t^4 -(1/24)t^8
        +(sqrt(5)/192)t^12 -(31/5120)t^16.                 (7)

All its substitution exponents are positive. The inverse cocharacter
therefore contracts it for the action on functions. It is pure-null and
lies in the nullcone. More quantitatively, any fixed marker Q of
matrix-entry degree at most e has

    h(Q E_(2I)(P_N)^m)=0 whenever 4m>e.                   (8)

The exponents of Q are between -e and e, so no zero torus exponent remains
in (8). This is an all-order safe example, not an exceptional counterexample.

For **Gamma=2O**, the negative part again vanishes. In the quaternion
orientation used by the checker, its adjoint orbit of the first coordinate
axis is precisely the six signed coordinate axes. Writing (x,y,z) for
the corresponding unit-sphere coordinates, the zero-exponent part is

    C(x,y,z) = (5/3)(x^4+y^4+z^4)-1.

It is real and nonzero. All remaining exponents are positive, so the first
two pure Haar moments of the full average are determined by C:

    h(E_(2O)(P_N))=0,
    h(E_(2O)(P_N)^2)=integral_(S^2) C^2 = 16/189.        (9)

Indeed E(x^4)=1/5, E(x^8)=1/9 and E(x^4 y^4)=1/105, which give (9)
directly. The checker confirms the axis orbit, (7) and this rational
moment value. Thus direct finite averaging fails in both standard
exceptional models: the octahedral average loses pure-nullity, while the
icosahedral average retains pure-nullity only by becoming unstable.

These conclusions concern E_Gamma(P_N). They do not decide averages of
higher powers E_Gamma(P_N^k), k>1, or other invariant families. Those
constructions require new exact moment calculations and an all-order proof.

## 8. Disconnected acting groups

Let now G have reductive identity component K and let B=K intersect H.
The connected quotient K/B has the stabilizer component group B/B^0.
If this group has no A_4,S_4,A_5 quotient, Theorem 1 applies without
sphericity to K/B. The component proof in
[the earlier theorem](SPHERICAL_PAIRS_AND_COMPONENT_NULLCONES.md) then gives:

    ker Reynolds on C[G/H] is Mathieu--Zhao iff [K,K] subset H;

    Z=N iff [K,K] subset H and the component character map is injective.

The negative case is algebraic zero extension from one component. The
positive case uses the equal finite component trace and GVC2FTS; moment
nullcone equality additionally needs one simultaneous character separator.

Use B/B^0 in this hypothesis. It is a subgroup of H/H^0, rather than
necessarily the whole component group. A prohibition on quotients of H/H^0
alone need not pass to its subgroups. Abelian H/H^0, connected H, and
12 not dividing |H/H^0| are sufficient for the needed subgroup condition.

## 9. The next mathematical question

The universal extension now asks for exact moment information on Y_O and
Y_I, in function spaces starting at highest weights 8 and 12. The finite
group and character-series calculations close the reduction and degree
obligations; they supply no exceptional moment witness.

A productive next construction must respect right binary-polyhedral
invariance from the start and establish all-order pure moment cancellation.
The line--doublet criterion, powers of the normalizer circuit, and direct
finite averaging have the explicit obstructions above. Any finite-prefix experiment must retain the
distinction between a pure-null unstable function and a counterexample
outside the nullcone; a nonzero marked tail would certify the latter.

The [spectral and orbit-span theorem](EXCEPTIONAL_SPECTRAL_BARRIER_AND_POWER_AVERAGES.md)
now excludes the full span of averaged complex right translates of P_N,
computes its ordinary moment ideals, and supplies an exact all-k spectral
second-moment test. Standard power averages k=2,...,6 fail that test on
both remaining models; untested powers and larger families remain open.
