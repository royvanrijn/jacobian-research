# Quantitative Lefschetz bounds from meridian centralizers

14 November 2026. PF2RULED47 gives a componentwise finite-monodromy
obstruction in every field degree. It does not require every ramified
prime to have zero excess or symmetric monodromy. The quartic two-prime
case and arbitrary higher-degree exclusions remain open.

Let C be a chosen ramified image component with immersed polynomial
normalization, and suppose every ramified prime over THIS support has
residue degree one and zero excess. Other branch-image supports and
their ramification are unrestricted. Use the original finite normalization
and its actual monodromy group G. Let sigma be the generic C meridian,
and b=[G:C_G(sigma)] its conjugacy-class size. On a smooth projective
target compactification X resolving and separating the entire branch
union at infinity, put

    eta=deg N_(P1->X)=-K_X.Cbar-2,
    delta=sum of finite delta invariants of C.

If eta>0, then

    b <= Cbar^2/eta = 1+2*delta/eta,
    (b-1)*eta <= 2*delta.                              (1)

All finite self-contacts are included through their delta invariants;
they are not assumed nodal. No positivity condition is needed on the
other image components. The target compactification, not a chosen source
completion or the complement obtained by deleting a prime, enters (1).

## 1. The exact neighborhood and the imported index theorem

Here is a more general representation form of the argument. Let X be a
smooth projective complex surface, E a removed boundary, and B a reduced
curve in X minus E. Fix an irreducible component C of B. Assume its
normalization h:H->Cbar is an immersion everywhere on X, and that it
meets E transversely at smooth points of E, away from the other strict
components of B. Fix a finite quotient

    rho:Pi=pi1(X minus(E union Bbar)) ->> G.

At every finite singular value of B on C, assume the image of the small
local complement group centralizes the meridian of each selected local
branch of C. This is a LOCAL representation hypothesis, not a claim that
the local group itself is abelian. It is enough for Section2. Let sigma
be a based generic C meridian.

Nori's neighborhood construction (Definitions1.8-1.10 and paragraphs
1.11-1.12, printed313-314) gives a connected complex surface V, a
holomorphic local isomorphism q:V->X and an embedding H->V whose
composite with q is h. H is compact; distinct normalization branches
stay distinct in V even when their images meet in X. Its normal bundle
has degree

    eta=H_V^2=deg(h^*T_X)-deg(T_H)
       =-K_X.Cbar+2g(H)-2=Cbar^2-2*delta.               (2)

The last equality is adjunction with p_a(Cbar)=g(H)+delta. The chosen
infinity resolution leaves no singularity of Cbar at infinity, so delta
in (2) is its finite affine delta. For the polynomial normalization H=P1,
(2) is exactly the eta in the statement. It is the normal bundle of the
immersed normalization, not the self-intersection of the singular image.

Put V'=V minus q^(-1)(E union Bbar). If eta>0, the line bundle
O_V(H)|H is ample. Nori's Weak Lefschetz Theorem parts A and C then give

    K=image(pi1(V')->Pi),
    [Pi:K] <= (Div h)^2/H_V^2 = Cbar^2/eta.             (3)

The input divisor Div h is Cbar because h is birational onto its image.
There is no residue-degree factor here. The primary statement is printed
305-306; the proof of part C and Definition3.25 are on printed330. This
is WLT(C), not Proposition3.26 with a deleted set disjoint from Cbar:
Cbar meets the boundary at infinity, so that latter substitution would
be invalid. The removed set in WLT is allowed to contain Cbar itself.

## 2. The monodromy of that whole neighborhood centralizes sigma

We prove

    rho(K) subset C_G(sigma).                           (4)

Shrink V around the compact H. Let S be the finite subset of H above
the singular values of B on C and above the boundary intersections.
Off disjoint small disks around S, the reduced pulled-back removed
divisor is just H. Its complement is an oriented punctured-normal-disk
bundle over H minus those disks. Its fiber meridian is central in its
fundamental group: complex normal-bundle transition functions preserve
the oriented generator of the normal circle. Thus its entire image
centralizes sigma. This includes base loops; they are not assumed
trivial or discarded.

Near a point of S lying over a finite singular value a, q is a local
isomorphism. The extra complement pieces map into the local complement
of B at a. By hypothesis their monodromy images centralize the selected
branch meridian. The paths joining these pieces to the normal-disk
bundle transport that meridian to the same sigma, so the based images
also lie in C_G(sigma). A singular value with several normalization
preimages gives several separate pieces; the selected meridians are
transported along their respective copies of H.

At a point over infinity the strict Cbar and E are transverse, with no
other strict branch component through that point. The local complement
is a product of two punctured disks. Its image is abelian and commutes
with the Cbar meridian. This step retains the meridian around infinity;
it does not assume identity preferred longitude.

These pieces cover a sufficiently small V'. Van Kampen gives generators
from their fundamental groups and the punctured normal bundle. The
underlying gluing graph is a central piece with one local piece for
each point of S; its attaching overlaps are connected annular collars.
No additional graph loop is introduced. Every generator has image in
C_G(sigma), proving (4). Tangential residual branches change the local
groups, but not this conclusion when their representation satisfies
the stated centralization condition.

Write L=rho(K). Surjectivity of rho gives the exact subgroup inequalities

    [G:C_G(sigma)] <= [G:L] <= [Pi:K].                  (5)

Combining (2), (3) and (5) proves (1). In particular the estimate applies
to any finite quotient satisfying the LOCAL centralization condition;
it does not require an abelian Orevkov quotient or connectedness of a
moduli space of polynomially parametrized curves.

## 3. Check the local condition on an actual Keller normalization

Return to the full finite normalization pi:N->A2 of a plane Keller map.
Choose a branch-image support C such that all ramified primes over C
have residue degree one and zero excess. PF2OMIT6 implies that C has
immersive normalization and these primes map birationally onto it.
PF2OMIT2's local cyclic lemma identifies their local finite factors
with (s,t)->(s^e,t), with the individual generic e retained.

At a finite singular value on C, each selected normalized branch has
one such cyclic factor for each ramified prime over C. Its branch
meridian sigma is the product of those cycles. Other branches of C
occupy disjoint local finite factors. Other branch-image supports also
occupy different factors, since omitted primes in N are disjoint and
the selected factor ramifies over only its own smooth branch. Those
other factors may have nonabelian monodromy, for example a cusp S3;
they still act on sheet clusters disjoint from the selected cycles.
Consequently their permutations commute with the selected meridian.
This proves precisely the local hypothesis of Section1.

The argument retains all preimages of the target value and the entire
finite normalization. A local factor is not obtained by restricting to
the nonfinite open U union one prime. Images of unramified omitted
primes add no branch component to the space used for rho.

Take any X which resolves and separates the FULL branch union at
infinity and gives each chosen Cbar a transverse intersection with E.
Apply (1) separately to every image support satisfying the selected
zero-excess condition. Normal bundles and canonical numbers must be
computed on that simultaneous X: separating another branch at infinity
may decrease eta. The canonical number of an isolated component alone
cannot be reused without the relevant mutual infinity intersections.

No statement here makes G symmetric in general. A central sigma has
b=1 and gives no inequality beyond delta>=0, even in a nonabelian G.
If selected primes have positive excess or ramified residue maps, their
image normalization or their local factor may fail the hypotheses.
No such case is silently treated as cyclic.

## 4. All-degree symmetric consequences and the remaining quartic case

There is an explicit degree bound even BEFORE the global monodromy
group is identified. Suppose sigma is a single transposition and G is
any transitive group on d sheets. Let T be the normal subgroup generated
by all G-conjugates of sigma. The graph of these transpositions has
equally sized connected components of size m>=2, because G is transitive
and permutes them. Edge transpositions generate the full symmetric group
of each connected component. Thus

    m divides d, T=product_(d/m copies) S_m,
    b=d*(m-1)/2.                                       (6a)

For the class-size equality, conjugation inside each S_m gives every
transposition of that component, and conjugation by G reaches every
component. No conjugate crosses these components by their definition.
So all and exactly (d/m)*binomial(m,2) transpositions occur in the class.
This is a proof for arbitrary finite degree, not a subgroup census.

It follows from (1) that

    d*(m-1)<=2+4*delta/eta,
    d<=2+4*delta/eta, for eta>0.                         (6b)

If ell(d) is the least prime divisor of d, then m>=ell(d), giving the
stronger d*(ell(d)-1)<=2+4*delta/eta. In prime field degree m=d and T=S_d,
so G=S_d and the sharper symmetric inequality below applies. These
components are the orbits of the normal transposition subgroup; their
size is not the ramification index or a source-coordinate pole. No
intermediate polynomial plane or target-invariant source action is
inferred. This supplies an explicit field-degree bound from the chosen
image geometry even when other ramified images have arbitrary inertia.

If the actual group is S_d and sigma has m_l cycles of length l,
including fixed letters, then its conjugacy class has size

    b=d! / product_l(l^(m_l)*m_l!).                     (6)

This is an explicit all-degree use of the actual generic ramification
profile. It does not assert that every formal profile is realizable.
For a sole index-two prime over C, sigma is a single transposition, so
b=binomial(d,2). Therefore for eta>0

    (binomial(d,2)-1)*eta <= 2*delta,
    d*(d-1) <= 2+4*delta/eta.                           (7)

This sharper symmetric estimate bounds the field degree by the target image geometry whenever its
immersed-normalization bundle has positive degree. At eta=0 or eta<0
the ampleness hypothesis fails, and (7) is not asserted.

In PF2RULED46's remaining quartic packet, the sole ramified image C has
S4 monodromy and eta=a_n-1 in that note's target infinity convention.
Every finite singularity is two smooth branches. If k are ordinary
nodes and T sums the intersection multiplicities of all tangential
double values, then delta=k+T. Formula (7) gives

    5*(a_n-1) <= 2*(k+T) when a_n>1.                    (8)

For a_n<=1 the inequality also holds automatically since its right side
is nonnegative. Combine (8) with the separate all-contact Nori bound
of PF2RULED46 to get the uniform necessary inequality

    a_n <= 1+min(2*T, floor(2*(k+T)/5)).                (9)

At eta>0 this adds a field-degree/singularity budget without discarding
tangencies. It is an exclusion test, not a bound forcing a finite list
of target image equations. Arbitrarily negative eta and large delta
are not excluded. The whole quartic two-prime case remains open.

## 5. A stronger uniform application to the retained degree-six gate

Retain the COMPLETE degree(6,10), gap21 first-image family C of
PF2R7FAM1, and suppose its distinct second ramified image B has
same-target-coordinate polynomial normalization degrees (3,5), as in
PF2R7B5CONTACT1. This second-degree assumption remains explicit and is
not derived for every unknown B. Consider the oriented residual rows
A,C,N of PF2OMIT13, in which the second image is immersive, of generic
index two, three, two respectively. Row N's first image can have its
ordinary cusp; it does not invalidate the selected B condition.

These covers have S6 monodromy. For A,N all generic branch meridians
are single transpositions, so transitivity gives S6 by the sheet graph.
For C the generic cycles are one transposition and one 3-cycle, so
PF2OMIT8's small-cycle normal-generation lemma gives S6. The second
image's conjugacy size b is15 for index two and40 for index three.

Let I=I_infinity(C,B). Reuse ONLY the infinity calculation of
PF2R7NORI1, which did not use affine nodality or finite incidence:
four initial shared centres have B multiplicities 2,2,1,1, and
n=(I-20)/2 subsequent shared centres have B multiplicity one. On the
simultaneous infinity model this gives

    Bbar^2=25-10-n=15-n.

Every birational degree(3,5) polynomial quintic has projective arithmetic
genus6 and its one infinity branch has delta2. Hence its TOTAL affine
delta is4, including all tangencies or multiple image values. An
immersive B therefore has

    eta_B=Bbar^2-2*delta_B=7-(I-20)/2.                  (10)

If I<34, the retained contact spectrum makes eta_B a positive integer.
But (1) would require b<=1+8/eta_B<=9, whereas b is15 or40. Contradiction.
Consequently

    I_infinity(C,B)>=34                                (11)

in EVERY such row A,C,N, without a nodality assumption on B, without
mutual finite singular-locus avoidance, and without a bound or penalty
for finite C/B tangencies. The complete previously proved polynomial
quintic contact spectrum then leaves only

    I in {34,38,40}.

At I34 the normal degree is zero; at I38,40 it is negative. Those cases
are not excluded by WLT(C). The earlier conditional bound I+2T_CB>=34
is preserved in its sealed PF2R7NORI1 note and receipt; (11) removes its
finite-tangency escape for the stated immersive second-image rows.
It does not assert that any remaining quintic is an actual Keller image.
The nonimmersive second-image row B is outside this application.

## 6. Source audit, controls and next gate

The essential new primary import is Nori WLT(A),(C), with its compact
normalization neighborhood rather than the image curve's neighborhood.
Complete PDF pages2,3,10,27,28 were rendered and visually inspected,
covering the theorem, neighborhood construction and quantitative proof.
The retained Neumann--Norbury convention identifies a_n=-K_X.Cbar-1;
the formula is used on the correct simultaneous infinity model.

The [source manifest](sources/lefschetz-centralizer/manifest.json) pins
these two PDFs, the predecessor PF2RULED46
receipt, and the unchanged PF2R7NORI1 infinity data. The checker verifies
the index/adjunction algebra, exact conjugacy classes and centralizers,
normal-transposition block controls, local cyclic-plus-cusp controls,
selected contact values and old source
seals. It does not certify Nori's theorem or turn finite group examples
into realized covers. The proof of neighborhood centralization and its
application to the actual finite local factors is written above.

[verify_lefschetz_centralizer.py](verify_lefschetz_centralizer.py) writes
[lefschetz_centralizer_receipt.json](lefschetz_centralizer_receipt.json).

The next general task is to extract more finite-algebra or peripheral
information when eta<=0, or compare several selected components on a
single simultaneous target model. The sole remaining quartic packet,
the remaining quintic contacts34,38,40 and unrestricted higher-degree
cases remain open. Formal verification, independent full reproof,
external review, novelty priority and full-goal completion remain false.
