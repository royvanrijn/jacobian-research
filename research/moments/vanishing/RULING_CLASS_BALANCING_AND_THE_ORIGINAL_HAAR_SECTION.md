# Ruling-class balancing and the original Haar section

17 November 2026. Written all-degree geometric construction from a supplied
proper resolved pencil. The full exceptional Mathieu--Zhao statements remain
**OPEN** for BOTH quotients. Independent replay, formal verification and
external review are not claimed. Review remains deferred.

The previous genuine-boundary theorem identifies the group but leaves the
original section as input. Here that section is obtained from a FIXED ruling
line bundle and a rational boundary intersection system. This is an explicit
NS/Picard interface from a resolved original phase, retaining the original
ambient class. It does not prove that purity forces its local component zero.

## 1. Proper resolved original data

Let G=O or I and X=PGL2/G=SL2/Gamma. The fixed projective compactification is

    Z=P3/G, B=Q/G,
    Q={det=0}=P1_left times P1_right,
    B=P1_left times (P1_right/G) ~= P1 times P1.

On Q the right action preserves the left coordinate. Therefore

    A_B=O_B(1,0)                                      (1)

is an actual line bundle downstairs. Its pullback is O_Q(1,0), not a
phase-dependent weighting of the volume.

For ANY nonconstant original regular f take a smooth projective W with
pi:W->Z and a resolved proper phase p:W->P1, isomorphic to X over the
original open. Put D=(W minus X)_red. Resolve the boundary and pencil in
boundary centres. Over a common full regular open S of the phase line,
Y_s=p^-1(s) is smooth projective and

    F_s=Y_s intersect X,
    C_s=Y_s intersect D

is the ENTIRE reduced SNC fibre boundary. All geometric connected components
and all node branches remain. D->B is a morphism, since its reduced support
maps into B. Define the fixed line bundle on the whole boundary

    A=(pi|D)^* A_B.                                  (2)

The original Haar unit is represented, up to a nonzero rational lattice
normalization, by c1(A) in

    H_c^3(X)(1)=H^3(W,D)(1).

Indeed upstairs the relative sequence for (P3,Q) gives

    H_c^3(PGL2)(1)=coker(Pic(P3)_Q -> Pic(Q)_Q)
                  =Q^2/Q*(1,1),                    (3)

and O_Q(1,0) maps nontrivially. Proper modifications in the boundary leave
the open unchanged. Rational finite transfer proves that (2) is nonzero
downstairs as well. The actual nonzero factor is fixed by pairing with the
original Omega; it must be retained as the lattice multiplier below.

## 2. A boundary intersection system always has a solution

Fix a generic fibre Y and label the geometric boundary components C_i. Put

    S_ij=C_i dot C_j on Y,
    a_i=deg(A|C_i).                                  (4)

This S is the intersection matrix on the PROJECTIVE FIBRE SURFACE. It is
not a negative Laplacian, and is different from the semistable normalization
matrices used later for local component calculations.

Then

    S c=a has a rational solution.                   (5)

Moreover beta=sum_i c_i[C_i] is a UNIQUE class in the rational span V of
the boundary classes in NS(Y)_Q. Only the coefficient vector can be nonunique.

Here is why solvability is an actual geometric theorem, rather than an
extra compatibility assumption. The compact real Haar cycle lies in
U_s=X minus F_s for sufficiently large |s|, and Omega has nonzero integral
on it. Thus H^3(X)->H^3(U_s) is injective there, and its constant line remains
injective on the regular locus. Compact-support duality gives a surjection

    H^3(W,D union Y) -> H^3(W,D).                     (6)

The target is pure weight 2. In weight 2, the relative sequence and the
proper Mayer--Vietoris sequence for D union Y show that the class of A has
a compatible rational degree-two lift on Y, modulo restriction of classes
from W. Polarizable pure Hodge structures are semisimple, so choose this
lift of type (1,1). The degree-two Lefschetz theorem makes its Y component
an NS class beta_0. If its D component differs from c1(A) by a restriction
from W, lift that difference by a rational (1,1) class on W and subtract
its restrictions on BOTH D and Y. Consequently

    beta_0 dot C_i=a_i for every i.                   (7)

Only the weight-two pieces are used: H^2(C) is pure of type (1,1), so lower
weights in H^2(D) or H^2(D union Y) introduce no missing degree condition.
This is a lift in the pure GRADED piece. It is not a Hodge lift through
the full mixed extension; the latter failure is exactly what (10) records.

The boundary B on Z is an ample Q-Cartier divisor: its finite pullback is
the ample determinant quadric. Restricting its pullback to Y gives a class
h in V with h^2>0, because Y maps birationally to its projective fibre
closure in Z. By the Hodge index theorem, the intersection form is
nondegenerate on V. Equivalently, the kernel of S consists precisely of
rational numerical relations among the C_i. Orthogonally project beta_0
onto V. Its projection beta has (7), proving (5) and uniqueness in V.
This reasoning applies separately to every connected generic branch.

For clarity about the classical inputs: the surface Hodge index statement
here follows from the primitive (1,1) Hodge--Riemann form and the positive
ample direction. It is not an assertion that the raw intersection form
polarizes all of H^2. Deligne, Theorie de Hodge II, Sections 2.2.4--2.2.6,
gives the Chern-class comparison and primitive polarization used; the
degree-two Lefschetz theorem and strictness are named classical imports.

## 3. The resulting line bundle is the actual original section

The data S and a are invariant under the finite permutation of generic
boundary components. Average a solution of (5) over this permutation group.
Its coefficients are constant on each generic component orbit, hence come
from coefficients on the horizontal irreducible boundary divisors D_j of W.
Clear denominators, choosing m>0, and set

    B_W=sum_j n_j D_j, n_j=m c_j,
    L_K=(A^m tensor O_W(-B_W)|D)|C_K.                (8)

Vertical boundary divisors may be assigned coefficient zero in this choice;
they remain part of the model and of every subsequent special-fibre
calculation. Formula (8) has degree zero on EVERY generic normalization,
so ell=[L_K] lies in J=Pic^0(C_K).

For the actual map

    psi:J -> G_f ~= Pic^0(C_K)/im Pic^0(Y_K),

the geometric comparison gives

    psi(ell)=[a_0]sigma_f, a_0 != 0.                 (9)

After rational lattice clearing and choosing the common orientation, one
can take a_0>0. This is the multiplier determined by the ORIGINAL ambient
class in (3), not an assumed generic degree or an unweighted substitution.

To verify the EXTENSION in (9), use the relative Picard construction, not
only the cohomology groups. In the proper normalization complex of
(W,D union Y), a compatible NS pair on D,Y maps to the mismatch line bundle

    A^m|C tensor O_Y(-B_W|Y)|C.                      (10)

This is exactly (8). It is the boundary morphism defining the 1-motive,
with all node trivializations and triple-intersection descent retained.
The published input is Barbieri-Viale--Rosenschon--Saito, Construction 3.1
and full comparison 5.4--5.6, already used for the actual boundary motive.
Thus (10) computes the extension class, not just its graded dimensions.

The negative Picard group in this Mayer--Vietoris presentation initially
quotients by both Pic^0(D) and Pic^0(Y). There is no extra discarded kernel:
H_c^2(X)=0, so H^1(W)->H^1(D) is onto. Full Picard realization then implies
that the image of Pic^0(D) on C is already contained, up to isogeny, in
the image of Pic^0(Y), by restriction from W. The quotient is precisely
the genuine G_f proved in the preceding note.

Finally the compatible NS pair in (10) is a GLOBAL invariant lattice class
on the original regular parameter locus. The retained canonical theorem
states that the compact-middle summand C_f has no global invariants there.
Hence this invariant lattice class, which maps to a nonzero original
ambient generator, belongs to its canonical rank-one ambient summand.
Pulling the 1-motive back along this rank-one class is exactly M_f. This
proves (9); no arbitrary choice of a compact-middle lattice is hidden.

Different solutions c change beta by a numerically trivial class, hence
change (8), after a tensor power if needed, by a restriction from Pic^0(Y).
They give the SAME point of G_f up to isogeny. This construction removes
the need to supply sigma_f or an arbitrary lift independently once the
resolved pencil and its actual line-bundle restrictions are supplied.
It is a written producer, not an implemented resolution/Picard algorithm.

## 4. An actual original NS representative localizes at the chosen root

Choose a left line ell in the determinant quadric; for a pure exceptional
phase choose the unique strict dominant left-root line. Upstairs blow it up:

    W_0=Bl_ell P3, D_0=Q' union E.

The restriction O(E)|Q' is O_Q(ell)=O_Q(1,0). Therefore the line bundle

    A_root=(pullback of O_Q(1,0) to D_0)
                       tensor O_W0(-E)|D_0          (11)

is TRIVIAL on Q'. On E=P1_bottom times P1_normal-row it is

    A_root|E=O_E(-1,1),
    (Q' intersect E)=diagonal.                       (12)

Indeed N_(ell/P3)=O(1) direct_sum O(1) and
O(E)|E=O_E(1,-1), while the left ruling bundle pulls back trivially from
the fixed left point. Restriction of (12) to the diagonal is trivial,
with the gluing inherited from (11). Its class on E minus diagonal is
nonzero, in agreement with the original quadric residue calculation.

The right G action preserves ell and all these data. Take a common positive
power to descend the equivariant line bundle through the finite quotient
(killing finite stabilizer characters), then pull it to a smooth proper
phase resolution dominating W_0/G. Equivalently clear the Cartier index
of the quotient exceptional divisor. This changes only the recorded
nonzero lattice multiplier. The ORIGINAL Omega is retained by rational
finite transfer and the comparison (3).

Consequently on the final boundary D there is an actual representative
of a multiple of the ambient Haar lattice which is trivial on the entire
reduced union of components mapping into Q'. Possible nontrivial
restrictions lie above E. This is a line-bundle support statement with
its gluing, stronger than merely locating poles of a differential form.
It is still a statement about the UNBALANCED ambient representative.

Replacing A by A_root does not change the balanced point (9). If an ambient
line bundle O_W(T) is subtracted from A, subtract T|Y from beta as well.
Since T is a boundary divisor here, its restriction belongs to V, and
the uniqueness of the boundary projection gives exactly that change.
The two terms cancel in (8). This explicitly records the compensation
which must not be dropped when localizing the original class.

## 5. An exact original-phase control recovers its known Kummer section

On PGL2 itself take f=xw/Delta, the original SL2 phase xw. This control is
not right-O/I invariant and is not pure. Its generic compact fibre is

    Y_s={(1-s)xw+s yz=0}, s!=0,1.

The boundary has four lines: top-row zero and bottom-row zero, with ruling
degree a_i=0, and the two column-zero lines, with ruling degree a_i=1.
Opposite lines have the same ruling class on Y_s. Thus the boundary
projection beta in (5) is the class of the top-row-zero line ell. Blowing
up ell realizes this class as E|Y_s. The balanced line bundle (8) is exactly
A_root|C_s; no further NS class is guessed.

The other three boundary lines lie in Q', where A_root is trivial. On E,
the remaining boundary line is the graph

    [z:w] -> [s z:(s-1)w]

from the bottom-row P1 to the normal-row P1. Substitution gives
(1-s)*(s z)*w+s*((s-1)w)*z=0, as required. The line bundle
O_E(-1,1)=Hom(normal tautological line,bottom tautological line) is
trivialized on this graph by diag(s,s-1)^(-1). At the two intersections
with the diagonal its values relative to the diagonal identity
trivialization are 1/s and 1/(s-1). Orienting the actual square boundary
cycle to take the second value divided by the first gives holonomy

    g(s)=s/(s-1),

with the cycle oriented accordingly. It has order 1 at phase zero.

This is exactly the previously recorded original ordinary Cauchy function

    C_1(s)=integral_0^1 dt/(s-t)=log(s/(s-1))

on the physical determination near infinity. The geometric mismatch
construction thus recovers an ACTUAL section and its zero-value winding,
including the original period normalization up to the recorded lattice
factor. It is not the assignment of an arbitrary obstructing Picard point.
This control provides no exceptional pure counterexample and no inference
from exponential generating functions.

## 6. The exact corrected coupling now has specified geometric input

Apply the [genuine-boundary quotient](GENUINE_BOUNDARY_PICARD_QUOTIENT_AND_LOCAL_CYCLE_REFLECTION.md)
and [raw Picard construction](RAW_PICARD_DEGREES_AND_A_CONSTRUCTED_IDENTITY_CORRECTION.md)
to (8). For s=u^e, take marked regular strict semistable models of all
generic boundary normalizations, and retain every special component E_(v,k).
Extend (8) using its actual restrictions, with marked blowups if needed.
More specifically these models can dominate the corresponding actual boundary
surfaces D_j after base change. Their two paired node sections have the same
map to D: they agree generically, and proper separatedness gives agreement
at the closed point. Thus the glued nodal model maps to D, and the fixed
line bundle A^m tensor O_W(-B_W)|D pulls back to an INVERTIBLE raw extension.
No separately guessed node valuation is needed. Any additional component
contracted to a point of D has degree zero for this pulled-back extension.
This does not erase special curves already present on D, which can carry
nonzero degree, as the Long control shows. All components are retained.
Let

    delta_(v,k)=deg(L|E_(v,k)),
    M_v d_v=delta_v.                                 (13)

Here M_v IS the negative Laplacian intersection matrix on a normalization
model, unlike the projective-surface matrix S in (4). The actual abelian
identity correction is constructed by the preceding degree-balancing
lemma. For every intrinsic character chi,

    nu_f(chi)=1/(a_0 e) sum_edges (psi^*chi)_b
                         (d_tailAttachment-d_headAttachment). (14)

All rational cycles of the genuine generic boundary graph occur, because
psi has proper kernel and induces a generic torus isogeny. Under eventual
purity, ALL original markers are safe if and only if every value in (14)
is zero. Neither S c=a nor the mere existence of an invertible extension
of (8) implies that further vanishing.

On any whole normalization model mapping into the good union Q', A_root
is trivial. There the raw degree input is entirely the balancing term:

    delta_(v,k)=-deg(O_W(B_W)|E_(v,k)),               (15)

with the appropriate already-cleared tensor powers. Thus the support
calculation does not erase good components: it identifies exactly what
their remaining correction is. On models meeting the root region keep
the A_root contribution as well. All exceptional degrees and both
attachments at every node remain in (13)--(15).

## 7. What is still required before contraction or full safety

To contract a proposed good subgraph, first show that (14) is zero on
EVERY cycle of that good subgraph. Formula (15) now expresses this as
an intersection/attachment identity for the actual boundary-balancing
divisor B_W. Triviality of A_root on the good union alone proves no such
identity: the compensation B_W is part of the original section.
After justified contraction retain all loops, parallel edges and remaining
attachments, and test their periods by (14).

The ambient resolution graph, the generic fibre-boundary graph, the
projective surface matrix S, and the semistable normalization matrices
M_v are distinct data. Neither an ambient tree nor a tree among one
normalization's special components forces (14) zero. The already recorded
actual Long control has nonzero period 1 with raw degrees (1,-1); it
prevents any unrestricted geometric zero theorem covering all phases.
It is not right-2O/2I invariant and supplies no exceptional counterexample.

The outstanding exceptional theorem is now a specific identity: use purity,
the strict dominant root and the retained capacity/automorphism constraints
to prove that the degrees of THIS balanced ruling representative annihilate
every actual attachment current. Alternatively construct an original
exceptional invariant pure phase for which (14) is nonzero and supply one
fixed original regular bad marker. No such vanishing theorem or phase is
proved here. All family-transfer and primitive-branch results remain intact;
no ordinary-resolvent algebraicity is inferred from Gaussian EGFs.

Both full cases remain OPEN, and external review remains deferred.
