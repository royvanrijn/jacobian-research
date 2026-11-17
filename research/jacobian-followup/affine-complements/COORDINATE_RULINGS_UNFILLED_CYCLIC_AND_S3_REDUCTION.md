# Coordinate rulings, unfilled cyclic periphery and the S3 reduction

15 November 2026. PF2RULED90 gives three degree-free tools and reduces
ALL incoming S3 state2 trajectories, in BOTH source pencils and at
arbitrary companion depth/later tail length, to one explicit solid
trefoil history. No prefix eta hypothesis is used in this reduction.
It does not close that remaining history or the full goal.

## 1. Coordinate rulings on the WHOLE affine source

Let N be smooth affine over C, U=A2_(x,y) a dense open, with complement
finitely many disjoint closed affine lines R_1,...,R_b, with b>=1.
The empty-boundary case is simply N=U and is separate. No volume or
cover-degree assumption is needed here. EVERY polynomial coordinate
u of U extends to a morphism N->P1. All omitted lines are vertical.
For a generic LINEAR coordinate u=a*x+b*y, in fact

    pi^-1(A1)=U,
    (pi^-1(infinity))_red=R_1 disjoint-union ... R_b.     (1)

The only possibly degenerate fibre is at infinity. Its multiplicities
are the positive integers -v_Ri(u). Its reduced components are exactly
the omitted lines, not new unnamed components.

First, every nonconstant map A1->an affine variety is finite onto its
image: choose a nonconstant coordinate polynomial p(t), and C[t] is
finite over C[p(t)] and hence over the image coordinate algebra.
Consequently each coordinate line {u=c} in U is CLOSED in N.

At a hypothetical indeterminacy point of u:N-->P1, write u=a/b in
the regular local UFD with no common divisor. Both a,b vanish there.
For general c, the local curve a-cb=0 has a branch outside the finitely
many boundary components. It would put the point in the closure of
{u=c}, contradicting that this affine line is already closed. Thus
the map extends. The same argument, or purity of fibres of a map
from a smooth surface to a curve, rules out a horizontal R_i: a
general fibre through it would acquire an extra point on {u=c}.
Each finite fibre is its original line {u=c}, with multiplicity1,
plus any vertical omitted components over c. The infinity fibre
consists only of omitted components.

At least one of x,y has a pole at EACH R_i. Otherwise shrink N around
U and the generic point of R_i, removing all remaining pole divisors
and isolated exceptions. Both coordinates then define a morphism
N'->A2 restricting to the identity on U. Its section U->N' is closed
because the morphism is separated; U is also dense, a contradiction.
A generic linear combination therefore has a pole at every R_i:
each valuation excludes at most one projective coefficient choice.
This proves(1), including simultaneous treatment of every component.

For ANY coordinate-induced ruling, its image base is A1 or P1, and
Euler additivity gives

    sum_t(number of reduced fibre components over t -1)
       = b if the base is A1; b-1 if it is P1.           (2)

Indeed chi(N)=1+b. The supports were identified above, so this does
not require a classification of degenerate fibres. In the generic
linear ruling the budget b-1 occurs entirely at infinity. The
complete-type base P1 is essential; no affine-base or Ga-action
conclusion follows. Dubouloz--Kishimoto1212.0521v1 Section1.2.2
independently records extension of cylinder projections on smooth
affine surfaces; their Section1.1 warns against the Ga inference.
The stronger coordinate/whole-fibre statement here has the preceding
direct proof. It supplies a global source constraint, not an immediate
exclusion of the remaining finite packets.

## 2. ALL-DEGREE unfilled cyclic-pattern peripheral transfer

Let K be ANY oriented knot pattern of signed winding w in an
unknotted solid torus V, and gcd(n,w)=1, n>=2. Restrict the STANDARD
cyclic n-cover of S3 minus K to V minus K, and call it P. There is
one OUTER and one INNER boundary torus. Keep the ACTUAL preferred
bases

    outer: m_o=lift(M^n), l_o=lift(L);
    inner: m_i=lift(mu^n), l_i=lift(lambda).

There are integral splittings, with T=H1(Sigma_n(K)),

    H1(P)=Z*m_i direct-sum Z*l_o direct-sum T,
    m_o=w*m_i, l_i=w*l_o.                              (3)

Every boundary coordinate in T is ZERO. T need not be finite, and
no Seifert, cable or graph-manifold hypothesis is needed.

For m_o, lift a meridian disk of V, punctured at K. Its outer boundary
is m_o, and its inner boundaries sum to w copies of m_i. For l_i,
glue the outside lifted solid torus to P. The result is the full
unbranched cyclic exterior E_n(K), and Mayer--Vietoris kills l_o.
The preferred longitude bounds a lifted Seifert surface in E_n(K),
so l_i is an INTEGRAL multiple of l_o in H1(P). Its image in the
actual two-component link-exterior abelianization fixes that multiple
to w. Projection to the axis-meridian coordinate sends l_o to1 and
m_i to0; splitting it and using PF65's actual cyclic-exterior
splitting proves the displayed decomposition. This proof never
extends an exterior map through a meridian filling.

Filling m_i recovers PF87's filled-pattern transfer. More generally,
if a one-boundary companion E has actual maps

    H1(E)=Z*m_o plus T_E, l_o=-A*m_o,

attaching P gives

    H1(E union P)=Z*m_i plus T_E plus T,
    l_i=-A*w^2*m_i.                                   (4)

For an ACTUAL solid companion with primitive meridian A*m_o+b*l_o,
an ACTUAL solid continuation therefore has new meridian

    A*w^2*m_i+b*l_i.                                  (5)

Here the new solid topology must be established separately. The
quotient by the old cap has relation A*w*m_i+b*l_o; actual solidness
forces gcd(A*w,b)=1 and T=0, so(5) is primitive. In the cable
applications below the cap distance1 proves solidness. Homology alone does not recognize a
solid torus. No arbitrary-pattern oriented linking splitting is asserted.

## 3. A null-homologous filling core cannot be an interior L-space slope

Let E have one incompressible torus boundary and m a filling slope.
If its core K_m is ZERO in H1(E(m);Z), then

    m is outside the INTERIOR of the L-space interval. (6)

If the filling is not an L-space this is automatic. Otherwise an
interior slope makes K_m Floer simple, by Rasmussen--Rasmussen
1508.05900v2 Section3.1, in particular Lemma3.5 and its following
two-sided argument. Ni--Wu1205.7053v2 Proposition5.1/Theorem1.2
then forces a null-homologous Floer-simple knot in that L-space to
have genus0. Equivalently their correction-term expression has
PD[K_m]=0, so its right side vanishes and the primitive Seifert
class has2g=0. The resulting disk compresses boundary E, contradiction.
These Floer results are NEW primary imports, not independently
reproved here. Complete relevant HTML statements and proof paragraphs
are pinned in the manifest. The hypothesis is integral zero, not
merely rational torsion or zero self-linking. No graph or field-degree
bound occurs in(6); graph hypotheses enter the gluing application.

## 4. The capped S3 companion and its zero core at ALL lengths

For state2 let J be the actual companion knot. The S3 action has
one three-sheet orbit and one global fixed sheet. Let E^sharp(J)
be the connected THREE-sheet unbranched exterior cover with its
internal fixed-meridian boundary capped, leaving the transposition
boundary (m=M^2,L). Its normal meridian fill is W_3(J), with
index2 core R_2. This is an actual irregular D3 branched cover;
no full S4/A4 cover is being called dihedral.

The complete PF82 graph says that the first S3 image is on the
initial torus, or is a birth from state1(C2) or6(C3). Afterward
only cyclic2 self loops occur before the final2->3 arrival.
Initially write e=2b,o=3t, with gcd(2b,3t)=1. PF73 gives

    E^sharp: disk(b,2b,t), H1=Z*m plus Z/b, L=-bt*m;
    W_3: H1=Z/b, [R_2]=0.                             (7)

One can also get the zero core from the exact fibre-order formula:
the closed cones(b,2b,t) have lcm2bt, Euler-1/(2bt), fibre order1.
At a first S3 birth, all older companion restrictions are cyclic.
PF65 replaces them by the unknotted torus pattern and adds unmarked
cyclic-branched homology summands. The torus S3 model again has
zero R_2. Every state2 loop uses(3) on its active double orbit and
the genuine fixed-sheet meridian fill, so R_2 stays zero at all depths.

Thus if W_3 has finite H1, E^sharp has the primitive meridian
splitting Z*m plus H1(W_3), with L an integral multiple of m and
NO extra summand coordinates. If H1(W_3) has positive rank, this
unmarked free summand persists in the final attachment and later
tails, contradicting PF61's finite source-end homology. It can be
discarded at this gate, not silently divided out.

## 5. Nonsolid E^sharp excludes BOTH pencils at every later length

The final2->3 cable image is S3 on a DIFFERENT three-sheet set,
so the whole image is S4. Its actual equations give r=3s odd,
k even, s odd. After all true branch/fixed fills its main piece P
has two companion boundaries, annulus base and cone s. It attaches
to E^sharp(J) and the ORDINARY exterior E(J); the latter comes from
the old global fixed sheet and MUST be retained. The integral
meridian balance is m+n=0 (PF73/81). Drilling R gives a pair of
pants with cone s, including the negative-Euler product when s=1.

If E^sharp is nonsolid, both companions are irreducible and boundary-
incompressible graph manifolds. Their meridians lie outside their
L-space INTERIORS by(6): R_2=0 for E^sharp; the ordinary J is
null-homologous in its S3 meridian fill. For s=1, P is a product
matching those meridians, so the actual gluing is not an L-space.
For s>=3, P(n) has disk cones(s,r) and is nonsolid; m+n=0 makes
m its rational longitude. The finite unmarked summands above give
b1(E^sharp union P)=1. PF81's TWO-step meridian-balance argument,
using the retained graph-manifold L-space gluing theorem, applies
with ALL its nonsolid/incompressibility conditions checked. The
closed end is an irreducible graph non-L-space, hence has a nontrivial
left-orderable group by the retained published theorem. PF74 excludes
the common source. Its ACTUAL drilled R exterior is also irreducible,
incompressible and graph, so PF82's non-L-space tail propagation
excludes ANY later simple tail. No hypothetical prefix source is used.

For complementary sources the drilled root has TWO genuine companion
edges, contradicting PF77/79's rooted endpoint. Original adjacent
fibre determinants r*B-k, B the previous ordinary fibre coefficient,
are nonzero by coprimality, and finite-index lifts preserve distinct
fibre lines. If a later tail follows, that vertex instead has degree3,
still forbidden by the source path. All relevant pieces have negative
Euler bases and the full retained peripheral recognition applies.
This is not a comparison based only on an unfilled orbit graph.

## 6. Exactly which histories can have SOLID E^sharp?

Initial disk(b,2b,t) is solid ONLY at b=t=1: the trefoil T(2,3).
Its genuine cap is m+L, so A_0=1. First S3 births cannot be solid:
from C2 the fixed cap leaves a three-boundary piece attached to
TWO nonsolid unbranched old companions; from C3 it leaves an annulus
with cones(r/2,r), omitting1, attached to one nonsolid unbranched
cyclic companion. These are incompressible gluings. In particular
they are not solid merely because an ordinary branch cycle was filled.

A subsequent state2 loop has odd winding w>=3. If k is even, its
unfilled active annulus has cones(w,w); capping a solid companion
adds |Aw-k/2|>0, so it is nonsolid. If k is odd, its active fibre
is k*m+2wL and the new exterior is disk(w,|2Aw-k|).
It is solid exactly when |2Aw-k|=1. Equation(5) gives A_new=A*w^2.
Once nonsolid, incompressible gluing with further genuine annuli
keeps it nonsolid, with the actual nonzero fibre determinants.

Consequently EVERY surviving incoming S3 trajectory in EITHER
pencil has initial T(2,3), followed ONLY by the solid cyclic2 chain

    A_0=1,
    w_i odd>=3, k_i=2*A_(i-1)*w_i +/-1,
    A_i=A_(i-1)*w_i^2=(product_(j<=i) w_j)^2.           (8)

The chain may be empty. Its last S3 arrival has r=3s odd,k even.
Capping the solid component leaves an annulus with cones
(s,|Ar-k/2|), attached to the ACTUAL ordinary exterior E(J).
Orders1 are omitted; if both are1 the annulus is a product and
must be transported with the full peripheral map, not a fictitious
canonical edge. Further simple tails remain possible only outside
the previous PF81/83 exclusions. This is an exact necessary solid
history, not a claim that every such integer list is genuine infinity
data or realizes an affine/finite map. The remaining cases stay OPEN.

## 7. Independent dihedral linking control and evidence

Cahn--Catania--Chimgee--Del Guercio--Kendrick2112.14790v2,
Section6.1, computes mutual linking2 for the index1/index2 lifts
of the trefoil. Its displayed two bounding-chain linear systems
are replayed in the checker. An INDEPENDENT full unfilled three-sheet
CW calculation in our actual preferred bases gives

    L_2=-m_2+2m_1, L_1=2m_2-4m_1.

The off-diagonal coefficients agree with the published mutual value2
and with the integral balance2L_2+L_1=0. The paper's chosen diagonal
push-off framings are not identified with our preferred lifted
longitudes. This is a small independent mark/sign check, not an
implementation of the complete dihedral algorithm or a new exclusion.
The method of1611.10330 and Section7.1 of2112.14790 can be used for
further ACTUAL D3 components; they do not automatically apply to
the full S4/A4 end or prove a ribbon/cobordism hypothesis.

The [checker](verify_global_ruling_unfilled_cyclic.py),
[receipt](global_ruling_unfilled_cyclic_receipt.json) and
[manifest](sources/global-ruling-unfilled-cyclic/manifest.json) retain
the precise new primary HTML imports, prior complete topology audits,
full unfilled peripheral relations, capped S3 controls, all-length
state classifications and ordinary target-cover examples. Finite
controls support the written universal arguments; no formal
verification, external review or independently reproved deep Floer
theorem is asserted.

Remaining: the solid trefoil-square S3 history(8), in its still-open
common/complementary sectors; and PF89's later A4 terminal types
(2,4),(2,8),(4,10) with their order/tail restrictions. Earlier seed
closures remain. Full longer-stage/quartic/higher-degree/degenerate
goals and finite Keller-map realization remain OPEN. No novelty or
priority claim is made.
