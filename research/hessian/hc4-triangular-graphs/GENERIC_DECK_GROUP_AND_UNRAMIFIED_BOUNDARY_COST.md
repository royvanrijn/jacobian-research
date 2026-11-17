# Generic deck groups have a compulsory unramified boundary cost

Date: 2026-11-16. Claims HC4MYLA0 and HC4MYLA1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: the ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the alternative full marked-energy integrality implication is also OPEN.
The results below do not supply either implication. External review/outreach
remains DEFERRED.

The [KT weight-two filter](KULIKOV_EXACT_COFRAME_AND_TATE_BOUNDARY_FILTER.md)
alone does not make the ordinary field-trace pairing regular. We instead
retain its actual affine-plane chart and the canonical
[KU maximal etale completion](CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md),
and study their GENERIC rational deck group. This gives an all-degree
boundary restriction, plus a sharp family showing why cohomology and exact
volume cannot replace the missing polynomial target pair.

## 1. Actual setup and the canonical place where generic decks are regular

Let K:X=A^2->Y=A^2 be complex Keller, R=C[X], B=C[K], L=Frac(R),
F=Frac(B), d=[L:F]. Let Z be the full finite normalization over Y and U
its FULL etale locus over Y. KU0/1 proves that U is smooth AFFINE and
contains X openly. Its complement over X consists of precisely b0 omitted
UNRAMIFIED prime curves E_i, smooth disjoint copies of A^1. The omitted
ramified primes are outside U and are NOT included in b0.

KT1 and the integral smooth-pair Gysin sequence give

    O(U)^*=C^*, Pic(U)=Z^b0,
    H^0(U,Q)=Q, H^2(U,Q)=Q(-1)^b0,
    H^j(U,Q)=0 for j!=0,2,  chi(U)=1+b0.                  (1)

The first Chern class on the actual boundary divisor basis identifies
Pic(U) tensor Q with H^2(U,Q(1)), naturally under automorphisms.
Furthermore pi1(U)=1. Indeed any loop in U can be moved off the finitely
many smooth complex divisors E_i by ordinary real transversality: their
real codimension is two and the loop has dimension one. Thus pi1(X)->pi1(U)
is onto, and pi1(X)=1. This is a statement about the FULL actual U, not
about an arbitrary intermediate etale surface.

Set G=Aut_F(L), g=|G|. This is the finite GENERIC rational deck group;
in permutation notation it is N_M(H)/H for generic monodromy M and its
sheet stabilizer H. It is NOT the monodromy group M. KMFI47 already proves
that its subgroup extending polynomially on X is trivial. Nontrivial
generic decks can still have source poles, as its higher-dimensional
examples show; no generic-deck-triviality theorem is assumed.

Every sigma in G preserves the integral closure S=O(Z), hence acts
regularly on Z and preserves its etale locus U. Its action on U is FREE.
If sigma fixes u in U, the local inverse for U->Y makes sigma the identity
on an analytic neighbourhood of u. Since U is irreducible, sigma is the
identity everywhere. This argument works even above a target value over
which other points of Z are ramified.

The affine quotient

    q:U->V=Spec(O(U)^G)                                   (2)

is consequently finite etale of degree g; V is smooth affine and
pi1(V)=G. The map U->Y descends to an ETALE map e:V->Y, since this can
be checked after the finite etale cover q. Its field is L^G and generic
degree is d/g. V is an OPEN SUBSET of the full finite normalization in
L^G, contained in that field's etale locus. We do NOT identify it with
that ENTIRE etale locus: further points may be etale in the intermediate
field but ramified upstairs. The image of X in V need not be an open
embedding. Neither q|X nor e is asserted finite.

## 2. LA0: divisibility, the whole rational character, and quotient Picard groups

The same free action on the d generic sheets gives g|d. Euler
multiplicativity for the ACTUAL finite etale cover (2) gives

    c=(1+b0)/g=chi(V) is a positive integer.               (3)

It is positive because transfer gives H^j(V,Q)=H^j(U,Q)^G, so V has
only H^0 and H^2, and chi(V)=1+dim H^2(V,Q). Therefore

    g divides gcd(d,1+b0), and b0>=g-1.                    (4)

In particular no omitted unramified prime, or gcd(d,1+b0)=1, forces the
GENERIC rational deck group to be trivial. This does not force d=1:
many transitive nonregular monodromy actions have trivial generic decks.

There is a stronger representation statement, without an integral
permutation-lattice assumption:

    Q direct_sum H^2(U,Q(1)) is isomorphic to Q[G]^c.      (5)

Here is a direct character proof that avoids a noncompact Lefschetz
shortcut. For every subgroup A<=G, the free cover U->U/A has degree |A|.
Transfer and Euler multiplicativity give

    dim H^2(U,Q)^A=(1+b0)/|A|-1.                          (6)

For an element s of order n, the rational character has the same value
on all generators of <s>: its eigenvalues are roots of unity and their
power substitutions by units modulo n are Galois conjugate. Induct on n.
Using (6) for <s>, the sum of the traces on that cyclic subgroup is
1+b0-n. The identity contributes b0; every element of smaller order
contributes -1 by induction. The remaining phi(n) equal generator traces
are therefore all -1. At the identity the trace is b0. Adding the trivial
character gives the value gc at the identity and zero elsewhere, exactly
c copies of the regular character. Semisimplicity in characteristic zero
proves (5), including rational representation isomorphism. The Tate twist
in (5) is supplied by (1), not by a new purity assertion for arbitrary covers.

Naturality of c1 transfers this to Pic(U) tensor Q. Finite-cover norm and
pullback identify Pic(V) tensor Q with (Pic(U) tensor Q)^G. In particular

    rank Pic(V)=dim H^2(V,Q)=c-1.                         (7)

The torsion is also exact:

    Pic(V)_tors = Hom(G,C^*).                             (8)

For completeness, a line bundle on V with trivial pullback is the trivial
bundle on U equipped with a descent cocycle in O(U)^*=C^*. Such cocycles
are precisely characters of G, and changing a trivialization changes none
of them. Conversely a character descends to a line bundle. Any torsion
bundle on V pulls to zero in the free Pic(U), and every character has
finite order. This proves (8), not just an inclusion. Norm shows that
the remaining pullback image has the same rational rank as the invariant
classes, proving (7) without a Brauer-group vanishing assumption.

CRUCIAL MARKING: G need NOT preserve X or permute its chosen E_i. Formula
(5) is a RATIONAL REPRESENTATION statement, not a claim that the integral
boundary basis is a G-permutation basis or that its lattice is stably free.
Replacing (3) by g|b0 by an assumed boundary permutation is WRONG.

## 3. LA1: a noncyclic generic deck group costs at least twice its order minus one

Suppose equality b0=g-1 holds in (4). Then c=1, so (6) says that V has
the rational cohomology of A^2: it is a smooth Q-homology plane.
Its logarithmic Kodaira dimension is -infinity. We record the needed
argument to avoid importing an unjustified source-recognition assertion.

The dominant generically finite morphism X=A^2->V is q restricted to X.
Choose smooth SNC completions and resolve its graph only over the source
boundary. The inverse image of the target boundary lies in the source
boundary. Pullback of logarithmic differentials therefore injects any
nonzero m-th log-canonical section on V into one on A^2: locally dlog of
a boundary equation is a sum of integral multiples of source dlog equations
plus a regular form. Generic finiteness in characteristic zero makes
the top differential nonzero. But on P^2 with its infinity line,
K+L_infinity=-2H, so A^2 has no positive log pluricanonical sections.
The same is true after boundary resolution. Thus all positive log
plurigenera of V vanish. This uses the log-pullback foundation already
retained in UALOGVOL1, not a claim that V itself is an affine plane.

Now adopt ONE NEW PROVED EXTERNAL INPUT: Dubouloz--Palka,
The Jacobian Conjecture fails for pseudo-planes, Adv.Math.339(2018),248-284,
arXiv:1701.01425v2, Lemma2.9, printed/PDFp6. It states that a smooth
Q-homology plane of negative logarithmic Kodaira dimension has an
A^1-fibration over A^1, with reduced fibres A^1, and fundamental group
the free product of the cyclic multiple-fibre groups. In particular a
FINITE such fundamental group is CYCLIC. The statement is explicitly
credited there to Miyanishi/Fujita. The primary statement and its exact
hypotheses were read; the classification itself is adopted, not independently
reproved. The previously cached HTML and the newly cached actual 26-page
v2 PDF are pinned. No statement about arbitrary rational affine surfaces,
arbitrary log dimensions, or arbitrary c>=2 is imported.

Since pi1(V)=G is finite, equality in (4) therefore forces G CYCLIC.
Combining this with the divisibility (3), we obtain the all-degree bound

    G NONCYCLIC ==> c>=2 ==> b0>=2*g-1.                   (9)

This applies to EVERY subgroup of the generic deck group as well, by
applying the same quotient argument to that subgroup. For instance a
Klein-four DECK subgroup requires at least7 omitted unramified primes;
an S3 DECK subgroup requires at least11. These are group-order boundary
costs, not statements about natural degree-four or degree-three MONODROMY.
In a noninvertible Keller map at least one ramified omitted prime is
also necessary by retained purity, so the total omitted count is at least
2*g when G itself is noncyclic. No new field-degree bound follows.

The quotient V has an exact nonvanishing volume and exact global
coordinate coframe e^*dP,e^*dQ, because e is etale. Its canonical bundle
is trivial. This concerns the OPEN quotient V, not the full finite
normalization in L^G; it does not contradict the retained nontorsion
canonical class of a nontrivial finite intermediate normalization.

## 4. Sharp cyclic controls, including the ACTUAL plane chart

The lower bound b0>=g-1 and the weight/volume data cannot be improved by
discarding the actual target-plane requirement. For EVERY k>=2 consider
the KNOWN smooth Danielewski surface

    W_k={x^2*y=z^k-1}.

Dubouloz--Palka Examples2.11/2.12 describe its simply connectedness and
the following FREE cyclic action (their parameter a=k-1):

    (x,y,z)->(epsilon*x,epsilon^-2*y,epsilon*z), epsilon^k=1. (10)

The action preserves the defining equation. If a nonidentity element fixes
a point, x=z=0, contrary to x^2*y=z^k-1. It cycles the k components of
x=0 transitively.

There is an explicit open A^2 chart, with coordinates (x,t),

    z=1+x^2*t,
    y=((1+x^2*t)^k-1)/x^2.                               (11)

The quotient in (11) is a POLYNOMIAL. The chart omits precisely the k-1
components x=0,z!=1. At the retained component z=1 the inverse t is
regular: if R_k(z)=(z^k-1)/(z-1), then t=y/R_k(z), with R_k(1)=k!=0.
Together with the inverse t=(z-1)/x^2 off x=0 this proves the WHOLE open
embedding, not just a birational parametrization. Its complement is k-1
disjoint A^1s. Units are scalar and Pic(W_k)=Z^(k-1), by the same actual
chart localization argument as in (1).

The nowhere-zero volume

    Omega=dx wedge dz/x^2

extends regularly on W_k (on x=0 use dx wedge dy/(k*z^(k-1))). It is
EXACT, with the explicit invariant regular primitive

    theta=-(z/k)*((k+1)*x*dy+(k+2)*y*dx),  dtheta=Omega.  (12)

This is the m=2 specialization of the EXISTING UADANVOL1 formula; exact
volume on higher-exponent Danielewski surfaces is NOT claimed new.
Both terms in theta have cyclic weight0 under (10), so theta and Omega
descend to the quotient T_k=W_k/C_k. This is a smooth rational affine
Q-homology plane with pi1(T_k)=C_k, scalar units, exact nonvanishing
volume, and Pic(T_k)=C_k. It is NOT A^2. Its covering surface has the
same pure weight-two completion cohomology as (1), with b0=k-1, and the
Picard representation is the augmentation of the cyclic regular module.

The restriction A^2->T_k of the finite quotient map is even SURJECTIVE
ETALE of generic degree k: the k cyclic translates of the chart (11)
cover W_k. Its fibres have k points away from the image of x=0 and ONE
point on that reduced curve. Thus it is genuinely nonproper, although
both source and target have Euler characteristic1 and the target has
exact volume. None of this is a plane Keller counterexample: T_k is a
pseudo-plane target. No etale map T_k->A^2, or polynomial exact coordinate
coframe on T_k, is constructed. Those missing functions would be essential
to using it in the HC4/JC2 problem. This control keeps the actual plane
source and loses exactly the required plane target.

## 5. What this changes and what it does not

LA0/1 gives an unbounded, actual-source restriction on generic rational
deck symmetry, and determines its full rational Picard character. It
does not prove a nontrivial deck group exists: symmetric and alternating
natural monodromy packets usually have TRIVIAL generic deck group and
survive this test. It does not prove all generic deck groups cyclic:
noncyclic groups with b0>=2*g-1 remain open. It supplies no bound on b0
independent of g, no integral permutation-lattice theorem, and no regular
ordinary trace pairing from Hodge purity or a polarization.

The sharp cyclic pseudo-plane controls have actual open-plane sources,
pure weight-two completion cohomology, and descended exact volume. They
show why another unmarked cohomology or volume calculation cannot by
itself supply the required target coordinate pair. The universal
trace-dual and marked-energy integrality implications remain OPEN.

The checker checks exact identities, rational character arithmetic,
sharp chart/volume controls, source locks and scope. It is not a proof
of the published Q-homology-plane classification, algebraic descent, or
HC4. Independent replay, formal verification, external review and full
goal completion remain false.
