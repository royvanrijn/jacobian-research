# A quartic colouring obstruction and an etale triple cover of the sign surface

15 November 2026. PF2RULED54 gives an actual exclusion, uniform in source
pole orders and target image degrees: in the surviving quartic packet,
the target knot at infinity must have determinant divisible by3. It also
gives a finite linear lifting test stronger than three-colourability and
an intrinsic etale-cover gate on the sign double surface. None of these
conditions is asserted sufficient for a Keller map.

## 1. The whole quartic packet and the knot at infinity

Keep the original finite normalization pi:N->A2 of a degree-four plane
Keller map on its marked source U=A2. PF2RULED45/46 leave exactly one
ramified omitted prime R of index2 and one unramified omitted prime D.
Their residue degrees are one and their image normalizations are immersive.
The whole branch curve is C=pi(R), and the actual monodromy is S4.
At a double C value the two local branch meridians act by disjoint
transpositions. Other finite factors are unramified. The D image and all
its fibres are retained. No compatible ruling or special source chart
is required anywhere in this note.

The normalization of C is A1. Hence its link at sufficiently large
affine infinity is a knot K. The primary infinity epimorphism is

    pi1(S3 minus K) ->> pi1(A2 minus C).                 (1)

We use the argument in Leidy--Maxim, proof of Theorem4.7, complete PDF11.
It applies to an affine plane curve without a transversality assumption
at projective infinity. It was already retained in PF2RULED34; its
use here does not assume that N is the marked Wright ring of that claim.

Composing (1) with actual monodromy gives a surjection onto S4 sending
a knot meridian to a single transposition. Let V4 be its normal Klein
four subgroup. The action on the three partitions of four letters into
two pairs identifies S4/V4 with S3, and sends a transposition to a
transposition. Therefore

    pi1(S3 minus K) ->> S3, meridians -> transpositions. (2)

Identify the three transpositions with colours in F3. Conjugation of
one by another gives the Wirtinger rule

    c_under1+c_under2=2*c_over.

Because (2) is onto, the colouring is nonconstant. Conversely any
nonconstant Fox3 colouring gives such an S3 epimorphism. The correspondence
is also Lemma3.10 of Przytycki's primary colouring paper, PDF18.

The colouring group modulo constant colours is H1 of the double branched
cover of the knot, with F3 coefficients (Przytycki Theorem3.3, PDF16).
In a Seifert matrix V, that integral homology has presentation V+V^t:
cutting along the Seifert surface and gluing the two sheets imposes the
sum of its two normal push-offs. Its order is

    det(K)=|det(V+V^t)|=|Delta_K(-1)|.

Equivalently the Fox crossing matrix is the Alexander crossing matrix
at -1; delete a constant-colour column and a redundant relation. A
nonzero F3 nullvector means its determinant vanishes modulo3. Thus

    3 divides det(K).                                  (3)

Every knot-at-infinity candidate with determinant prime to3 is excluded
from the ENTIRE surviving quartic packet, at every source pole order
and without a target degree cutoff. This is a degree-four consequence
of the S4 quotient. It is NOT asserted for S_d when d>=5: those symmetric
groups have no corresponding S3 quotient.

## 2. A closed arithmetic test through any number of cable stages

Use ordinary cable notation C_(w,p)(J), where w>=2 is the winding,
p!=0, gcd(w,p)=1. Neumann--Norbury's notation has w=q_i and p=p_i.
Keep the actual knot cabling parameters, not the Puiseux-pair or
relative-twist parameters from another convention.

The satellite Alexander formula, in Livingston--Melvin Theorem2 and
its Seifert corollary, PDF6, is

    Delta_Cable(t)=Delta_T(w,p)(t)*Delta_J(t^w).         (4)

Put D=det(J). Evaluation at -1 gives the exact recurrence

    det(C_(w,p)(J)) = |p|,          w even;
                      |w|*D,       w odd and p even;
                      D,           w and p odd.         (5)

Indeed Delta_J(1)=1, while for a coprime torus knot the determinant
is the absolute odd parameter when exactly one parameter is even,
and1 when both are odd. Mirroring does not change these absolute values.
Starting with the unknot determinant1 computes the determinant of any
finite iterated-torus cabling sequence. Polynomially parametrized knots
have such descriptions, as retained in PF2RULED35/46.

Combining (3),(5) excludes whole infinite families, not a bounded list.
For example ANY final even-winding stage with 3 not dividing p is
impossible, regardless of all inner stages. Such a stage can erase
three-colourability of its companion. For odd winding and even p, the
factor w can instead create it. Passing this arithmetic test supplies
no polynomial curve and no finite Keller map.

## 3. A stronger linear test: can a colouring lift to S4?

Three-colourability alone is insufficient. The exact extra group condition
can be computed by linear algebra over F2 from ANY finite presentation.
It can be applied to the knot group as an infinity obstruction, or to
the full Orevkov quotient O(C) as a global obstruction. A knot-group
solution need not extend to the affine complement or kill the required
finite local commutators.

Fix an epimorphism phi:Gamma->S3 taking specified meridians to
transpositions. Identify S3=GL(2,F2), acting on V=F2^2. The affine group
V semidirect S3 acts on the four points of V and is S4. A lift of phi
is precisely a crossed homomorphism z:Gamma->V:

    z(gh)=z(g)+phi(g)z(h).

Each presented relation is a homogeneous linear condition on the two
coordinates of z at each generator. A meridian lifts to a transposition
exactly when its affine lift has order2, namely

    (I+phi(mu))*z(mu)=0.                               (6)

In that coset the other two lifts are four-cycles. Conjugating a lift
by a translation changes z by the coboundary z_v(g)=v-phi(g)v. Define

    H1_par(Gamma,V_phi) = {relation cocycles satisfying(6)}/B1.

Then

    a surjective transposition S4 lift exists
       iff H1_par(Gamma,V_phi) is nonzero.              (7)

Here is the group proof, not an inference from finite examples. If H is
the lifted image, H intersect V is invariant under all GL(2,F2), so it
is either0 orV, since this two-dimensional module is irreducible. In
the second case H=S4. In the first H is a complement isomorphic to S3.
Every such complement is translation-conjugate to the standard one:
H1(S3,V)=0. To prove that last fact, average over the normal subgroup
C3 (its order is invertible in F2) to remove a cocycle there; V^C3=0
then forces its value on a remaining transposition to vanish. Thus
complements are exactly the coboundary class, proving (7).

For the actual quartic packet, at least one infinity three-colouring
MUST have a nonzero space in (7). If all of them have zero space, the
knot is excluded even when (3) passes. The test is linear once phi is
chosen; it has no coefficient-depth or Laurent-order cutoff.

Complete small controls distinguish the implications: the trefoil has
determinant3 and nonconstant three-colourings, but its parabolic lift
space is zero, so it has no required S4 representation. The torus knot
T(3,4) has such a representation and a nonzero lift space. That is
only a group-at-infinity control; the older positive-braid and marked
geometry restrictions still apply to an actual candidate.

## 4. The sign double surface must have an etale cyclic triple cover

Let G define C in the actual affine target. The normal sign surface is

    S=Spec C[P,Q,zeta]/(zeta^2-G).

It is the sign subfield of the actual S4 Galois closure. Take the quotient
of that Galois closure by V4. This gives a connected S3 Galois cover of
the target, and factors as a connected cyclic cover T->S of degree3.
It is not a cubic subfield of the original quartic source field.

In this zero-excess packet the map T->S is etale at EVERY affine point,
including the singular points of S. Outside C this is immediate. Near
a smooth branch value the local S3 image is C2, killed after the sign
base change. At a double value the two S4 branch transpositions are
disjoint. Their product lies in V4, so their S3 images are the SAME
transposition. Again the entire local image is C2; its sign kernel is
trivial. Locally T is three disjoint copies of the normal sign germ,
including a tangency germ zeta^2=g1*g2. This proves actual etaleness,
rather than merely absence of divisorial ramification.

Consequently

    H1_et(S,mu_3) has a nonzero class alpha,
    tau^*(alpha)=-alpha,                               (8)

where tau is the sign involution. The inversion follows from conjugation
in S3. Thus vanishing of this anti-invariant group excludes the quartic
packet. The primary Kummer sequence (Stacks Tag03PK) gives

    0 -> O(S)^*/O(S)^(*)3 -> H1_et(S,mu_3)
      -> Pic(S)[3] -> 0.                               (9)

If all units are cubes, (8) requires nonzero Pic(S)[3]. We do NOT assume
that units are constant or that the Picard group is torsion-free without
evidence for the actual S. Passing (8) is still weaker than the full
S4 lift and the marked source A2 realization.

## 5. A necessary negative control for the sign-surface route

Smoothness alone does not make an affine sign surface simply connected
or exclude an etale cyclic triple cover. In coordinates (x,u), use

    f=(xu-2)^2-u^3,
    W(s)=s^3-3u*s-2(xu-2), Disc_s(W)=-108*f.            (10)

The reduced branch curve f=0 is smooth and isomorphic to C*, with

    u=t^2, x=t+2/t^2, t=(xu-2)/u, t!=0.

It is Oka's example as recorded by Neumann--Norbury, PDF4, after the
polynomial coordinate change u=x^2-y. The cubic is irreducible: as a
polynomial in x it is primitive linear with coefficients -2u and
s^3-3u*s+4, which are coprime. Its nonsquare discriminant gives S3.
Since the branch curve is smooth, the normalized S3 Galois cover gives
an actual connected etale cyclic triple cover of the smooth sign surface
zeta^2=f. This curve has TWO normalization ends, so it is not an allowed
ramified image for our quartic packet. The control nevertheless disproves
the proposed general shortcut from smooth sign surface to trivial covers.
No assertion that the same cover exists for an A1-normalized candidate
is made.

## 6. Verification and the remaining global gate

The [checker](verify_quartic_colouring.py) writes
[the receipt](quartic_colouring_receipt.json). Primary PDFs and exact
page audits are pinned in [the manifest](sources/quartic-colouring/manifest.json).
It checks the S4 quotient, transposition colour rule, closed cable
determinants, full finite F2 cocycle systems and their affine lifts,
and the smooth C* sign-cover counter-control. These support the written
unbounded arguments. The older finite cover, infinity epimorphism and
Kummer facts are imported, not independently reproved by these controls.

The result excludes determinant-prime-to3 infinity types in the entire
quartic problem, without assuming common or complementary source pencils.
It leaves the determinant-divisible-by3 types, their nonzero lift classes,
the global etale-cover gate, marked realization and higher field degrees
open. PF2RULED53's simultaneous finite-algebra conditions remain available
in their r2q1 scope. No full quartic/higher-degree exclusion, finite Keller
pair, independent full reproof, formalization, external review or novelty
priority is claimed.
