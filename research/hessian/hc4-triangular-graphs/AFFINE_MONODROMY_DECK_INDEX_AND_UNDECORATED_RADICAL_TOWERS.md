# Affine monodromy bounds generic decks and excludes the full radical towers

HC4MYMN0. 16 November 2026. Actual complex polynomial PLANE Keller maps,
their FULL finite normalizations and ENTIRE étale completions.

**Theorem.** Suppose the natural sheet monodromy of a nonautomorphic
plane Keller map is a subgroup of AGL(n,F_Q), acting on ALL d=Q^n
affine points, where Q is any prime power. For EVERY noncyclic subgroup
C of its generic rational deck group, with h=|C|,

    d/h >= 3*Q.                                         (1)

In particular for Q=2 the integer d/h is a power of two and is at least
EIGHT. For Q=p prime>=3 it is at least p^2. The deck group is NOT the
monodromy group. No condition on the generic inertia classes is imposed.

**Corollary.** For EVERY q>=2, the natural degree 4*2^q action of

    H_q=Aff(F2^2,F2^q) semidirect S4

is impossible as plane Keller monodromy, with ANY inertia classes.
This excludes degrees16,32,64,... in this natural family, without
degree-minimality or a degree cutoff. It strengthens the coordinator's
specified lifted-transposition exclusion by a different quotient argument.
It is not an exclusion of all binary affine monodromy, arbitrary actions
of the same abstract groups, or the q=1 degree-eight case. Full HC4/JC2
and the trace-pairing/marked-energy implications remain OPEN.

## 1. Credit and the low-Euler surface lemma

The coordinator's [inward coupling packet, Sections1--2](../../jacobian/cross-programme-connections/INWARD_COUPLED_RESTRICTIONS_2026_11_16.tex)
supplies the following NEW combination: a smooth connected affine complex
surface X with negative logarithmic Kodaira dimension, finite pi1(X),
chi(X)<=2, and a regular nowhere-zero EXACT two-form has CYCLIC pi1(X).
That combination is credited to the packet, not claimed discovered here.
The shared geometry draft PF2RULED132 audits it and its decorated tower
application. The stronger undecorated conclusion below does not require
the relative cover of the quartic intermediate to be étale.

Here is our audit of the lemma, including its two classical inputs.

* Adopt the standard Miyanishi–Sugie A1-ruling existence theorem, as
  stated in Dubouloz–Kishimoto,
  [1408.1328v1](https://arxiv.org/html/1408.1328v1), Section1.2,
  Theorem2, printed/PDFpage3. Over C no finite constant-field extension
  is needed. The complete statement, proof and adjacent discussion were
  read; the original deep Miyanishi–Sugie classification is NOT
  independently reproduced. The handoff's Section1.1.2 locator is
  corrected to Section1.2 for this version.
* Adopt the classical good/bad orientable two-orbifold classification,
  Thurston, [Chapter13](https://people.math.osu.edu/davis.12/courses/8800/Thurston13.pdf),
  Theorem13.3.6 and its table, printed312/PDF18. The full pertinent
  proof pages18--24 and the bad-orbifold discussion were read; page18
  was visually checked. The underlying geometric classification is an
  explicit standard input, not proved by our finite controls.

Use ML0's ALREADY ADOPTED local structure: fibre reductions are disjoint
unions of A1s, with positive multiplicities; one reduced irreducible A1
in every fibre gives a Zariski A1-bundle. No second import of that package
is made. Let rho:X->B be the ruling supplied above. Restoring special
fibre divisors to the ordinary A1-bundle gives the orbifold presentation
of pi1(X): each base meridian gamma_i is subject to gamma_i^m_i=1, where
m_i is the gcd of ALL component multiplicities of that fibre. This follows
by van Kampen, since a component meridian projects to gamma_i raised to
its own multiplicity. It also agrees with the retained general ruling
presentation in the Miyanishi2005 report used by PF130.

The group surjects onto pi1(B), so finite pi1 leaves B=A1 or P1. For
A1 the presentation is a free product of cyclic groups, finite only
when at most one nontrivial factor occurs, hence cyclic.

For P1, Euler integration gives

    chi(X)=2+sum_i(number of reduced fibre components -1).

Thus chi(X)<=2 makes EVERY fibre irreducible. If the finite orbifold
group is noncyclic, it is GOOD and spherical, with signature (2,2,n),
(2,3,3), (2,3,4) or (2,3,5). Its universal orbifold cover is a finite
holomorphic map P1->P1 with exact local ramification orders m_i.
Bad teardrops and unequal spindles have cyclic groups and are NOT
assigned such a cover. The cover can be taken algebraic on these
projective complex curves.

Normalize the corresponding base change X'->X. Its generic field is the
base extension of C(B)(z), so it is connected. At a fibre m_i*C_i, the
normalized root extension has ramification index m_i/gcd(m_i,m_i)=1.
It is unramified in codimension one; purity over smooth X makes it
finite ÉTALE everywhere, including possible codimension-two points.
Consequently X' is smooth affine and its pulled-back form remains exact
and nowhere zero. Its ruling over the covering P1 has REDUCED fibres.
They can SPLIT into several A1s: finite étale covers of each C_i=A1
are trivial, but the fibre need not stay irreducible. Generic reducedness
and the Cartier divisor's Cohen–Macaulay property ensure reducedness
of the whole fibre, with no embedded components.

Now choose a multiplicity-one component in every fibre and apply ML0's
FULL component-selection theorem, not just a bundle assertion for the
unselected surface. It prohibits precisely this complete-type ruling
on an exact symplectic affine surface. This proves the low-Euler lemma.
Neither ruling existence nor the good orbifold cover has been silently
assumed in the exceptional bad-orbifold cases.

## 2. Apply the lemma to the free deck quotient, not a block intermediate

Let K:X=A2->Y=A2 be an actual plane Keller map, d>1. Write N=Spec S
for its full finite normalization and U for the ENTIRE étale locus
of N->Y, including ALL omitted unramified primes. Retained KU/LA0 give

    U smooth affine, X open in U, pi1(U)=1,
    chi(U)=1+b0>0.

Here b0 counts only omitted UNRAMIFIED primes. Put
D=Aut_(C(Y))(C(X)), the generic rational deck group. LA0 proves that
every subgroup C<=D acts REGULARLY and FREELY on U, with

    U -> V_C=U/C finite étale of degree h,
    pi1(V_C)=C,       chi(V_C)=(1+b0)/h a positive integer. (2)

It is crucial to use this free quotient. It need NOT equal the entire
étale locus of its intermediate finite normalization, and the original
plane chart need not be C-invariant. Neither identification is needed.

All hypotheses of the surface lemma hold on V_C except possibly the
Euler upper bound. It is smooth affine by finite étale descent. The
dominant generically finite map X=A2->V_C makes its logarithmic Kodaira
dimension negative: the log-canonical pullback argument is already
proved in LA0/1 Section3. The target functions P,Q descend, the map
V_C->Y is étale by descent, and

    omega=dP wedge dQ=d(P*dQ)

is regular, nowhere zero and EXACT. In particular the primitive P*dQ is
C-invariant, so exactness is not merely guessed from trivial canonical
bundle. It is not transferred to the ramified part of N, where the
pulled-back target volume can have zeros.

The low-Euler lemma and (2) now give, for EVERY noncyclic C<=D,

    chi(V_C)>=3,      chi(U)>=3*h,      b0>=3*h-1.        (3)

This cost is also written in the shared PF132 draft. It strengthens
LA1's earlier 2*h-1 cost through the coordinator's lemma. We credit
that strengthening; the new use here is its combination with ALL
nontrivial affine inertia, without a relative-etaleness hypothesis.

## 3. The complete affine support estimate

A nonidentity affine permutation of F_Q^n has either no fixed point or
a translate of ker(A-I) as its fixed set. If it has a fixed point and
is nonidentity, A!=I and that kernel has dimension at most n-1. Thus

    #fixed <= d/Q,
    moved support >= d*(1-1/Q).                          (4)

This holds for EVERY nonidentity element, including order-four lifts,
vertical translations, and all other inertia types in the natural
affine action. The source field remains C; F_Q describes ONLY a finite
permutation representation, not the characteristic of the Keller map.

Apply retained PF129 to M=N itself. If A_tot is the sum of generic branch
MOVED SUPPORTS, and Eta>=0 its full finite excess, then

    chi(U)=d-A_tot-Eta.                                 (5)

At least one branch divisor exists: otherwise purity and connectedness
would make N a finite étale cover of simply connected A2, of degree1.
Every generic branch permutation is nonidentity. Equations (4)--(5)
therefore yield

    0<chi(U)<=d/Q.                                     (6)

No critical point, omitted divisor or residual ramification is dropped
in (5). A support is the number of moved letters, not the permutation
index or e-1. In fact two branch images would cost at least d, so the
same argument forces ONE branch image, though the proof of (1) needs
only the nonempty branch and (6).

Combining (3) and (6) proves 3*h<=d/Q, hence (1). LA0 gives h|d, so
d/h is a power of the characteristic prime of F_Q. In the binary case
it is a power of two at least6, hence at least8. For Q=p>=3 prime it
is a power of p at least3p, hence at leastp^2. These are necessary
inequalities; no profile at equality is asserted realized.

## 4. Entire natural radical towers, with no inertia decoration

Use the [knot construction, Section2](../../meridional-rank/representations/LIFT_DIMENSION_AND_RADICAL_TOWERS.tex)
for the definition and credit of H_q, and the coupling packet for its
deck-centralizer observation. Set W=F2^2, R=F2^q and t=2^q. Its natural
action on W direct_sum R has degree d=4*t and consists of

    (x,z) -> (g(x), z+b(x)),
    g in S4=AGL(2,2), b:W->R an AFFINE function.

Every such map is affine on F2^(q+2), so (1) applies with Q=2.
The CONSTANT radical translations form a central subgroup C=R of
order t in the permutation centralizer. They therefore supply generic
deck transformations. One can also check that this is the full centralizer:
a commuting permutation first has form (x,z)->(f(x),z+c(x)); affine
functions separate x-points and force f(x)=x; the base S4 then forces
c constant. No transvection or branch decoration enters this argument.

For q>=2, C=(C2)^q is noncyclic but d/|C|=4, contradicting the binary
bound8. Equivalently, (6) makes chi(U/C)<=2 while pi1(U/C)=C is
noncyclic, directly contradicting the audited surface lemma.

This argument applies EVEN WHEN the quartic block intermediate has
residual ramification. It uses the étale free-deck quotient U/C, not
the potentially ramified relative map over that intermediate's entire
étale locus. It therefore legitimately removes the coordinator packet's
inertia restriction for the FULL NATURAL H_q ACTION. The older decorated
statement remains valid and is strictly weaker here. A claim about
different transitive actions of the same abstract H_q would need a new
permutation/deck calculation and is NOT made.

For q=1 the subgroup C2 is cyclic, so this proof makes no exclusion.
General binary affine monodromy with trivial or cyclic generic deck
group, or with noncyclic deck index>=8, also remains outside this
exclusion. The published generic-degree-four PLANE exclusion is retained;
q=0 is not reopened. The knot representations themselves remain valid:
their source groups have none of the affine-surface hypotheses used here.

## 5. Evidence, adoption and the original goal

The proof audits the coordinator's low-Euler combination and its genuine
surface hypotheses. Its classical ruling/orbifold inputs retain the
explicit foundation limits in Section1. The shared source packet and
primary PDFs are pinned; its checker and every prior registered checker
are NOT run. Its Mathieu and Factorial assertions are not adopted here.
Finite new controls check the all-affine fixed-point estimate, generic
deck arithmetic, order-four/vertical inertia countercontrols and the
distinction from a free quotient. The all-Q,n,q proofs are (1)--(6),
not an inference from an enumeration.

FULL THEOREM TARGET: prove or break full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.                 OPEN

Full marked-energy integrality, cyclic finite-target realizations and
the remaining monodromy classes are OPEN. No actual noninverse, blanket
solvable-monodromy theorem, formal verification, independent replay,
external review or novelty certificate is asserted. Review and outreach
remain DEFERRED.
