# Cyclic pattern mark transfer and the A4 incoming reduction

15 November 2026. PF2RULED87 proves an integral marked replacement
for arbitrary knot patterns and all cyclic cover degrees with coprime
winding. It closes ALL incoming A4 state7 COMMON-source trajectories,
at any companion depth and later tail length, without prefix eta.
In the COMPLEMENTARY pencil, every first A4 birth at a cable and every
initial non-trefoil A4 torus trajectory also closes. The remaining A4
cases start with T(2,3), up to interchanging exponents, and its actual
fixed-sheet cap; only trajectories outside earlier closed seeds remain.
The full quartic/higher-degree goal is OPEN.

## 1. Integral marked replacement for any cyclic knot pattern

Let K be any oriented knot pattern in an unknotted solid torus V in
S3. Its signed winding is w. Take the n-fold cyclic branched cover
of S3 along K, with n>=2 and gcd(n,w)=1. Write P for the inverse
image of V. The inverse image of the complementary solid torus is
one unbranched solid torus; hence P has ONE boundary. In the actual
outer basis write

    m=lift(M^n), l=lift(L),

where M,L are the meridian/preferred longitude of V. L has zero
linking with K; M has linking w and generates the cyclic action.
Let R be the ramification knot in P. Then INTEGRALLY

    H1(P)=Z*l direct-sum H1(Sigma_n(K)),
    m=0, [R]=w*l.                                      (1)

No finiteness, graph-manifold, Seifert or cable assumption is needed.
The second summand may have positive rank. In (1) all boundary and
marked coordinates in that summand are ZERO.

Here is a direct proof retaining the markings. Lift a meridian disk
of V transverse to K. Its full branched preimage is an oriented
surface with ONE boundary, the closed lift M^n, since w is a unit
modulo n. Thus m=0. The genuine branched map P->V induces a homology
map to Z sending l to1: L lifts with degree1. It splits off Z*l.
Attaching the complementary lifted solid torus kills l. The H0 map
in Mayer--Vietoris is injective, so the quotient is exactly the
homology of the closed cyclic branched cover Sigma_n(K).

The ramification knot in that closed cover is nullhomologous. Indeed
an oriented Seifert surface for K has all loops of linking number
zero; a copy lifts, and its boundary can be joined to R in the
branched filling solid torus. Therefore the image of R in H1(P)/<l>
is zero. Thus R is an INTEGRAL multiple of l. The actual map P->V
sends R to K with degree1, determining that multiple to be w.
This proves (1), including the absence of a hidden torsion mark.
No map to a knot EXTERIOR is extended after filling its meridian.

If an arbitrary connected companion piece E is attached along this
boundary, with the true m,l identifications, Mayer--Vietoris gives

    H1(E union P)=H1(E(m)) direct-sum H1(Sigma_n(K)),
    [R]=(w*[l],0).                                     (2)

Other boundary fillings on E can be retained throughout. This is a
homology and mark theorem, not a homeomorphism or unconditional
oriented linking splitting. Coprime winding and the standard cyclic
restriction are essential: multiple outer boundaries need a different
argument. Nor may one replace a noncyclic companion by an unknot.

For a cable (r,k), gcd(r,k)=1, the active cyclic n-sheet pattern has
w=r, provided gcd(n,r)=1 and L has identity monodromy. Formula (2)
proves the mark transfer r*l for ALL k, including gcd(n,k)>1.
This avoids trying to extend a two-cone lemma to a many-cone pattern.

## 2. All-length A4 histories and a two-torsion companion mark

Retain PF80/82/84's actual preferred words, full branch/fixed meridian
cycle fills and complete period12/all-length state graph. State7 is
A4 with M a3-cycle plus one fixed sheet and L=1. Its nonself
predecessors are exactly states6(C3) and8(double-transposition C2).
After its first A4 image it has only cyclic3 self loops until final
arrival7->3. Each such loop has gcd(r,3)=1, with the active cyclic3
orbit and one degree1 fixed sheet. These assertions do not bound
the number of loops or assume that a prefix has a source packet.

Fill ALL companion meridian cycles, and denote by R3 the core of
the index3 filling. At the first A4 level its order divides2:

* If the image was already A4 on the initial torus, its generator
  orders are2 and3, because A4 has trivial centre and the exponents
  are coprime. Write e even, o=3t odd. The closed cover has base S2,
  cones(e/2,e/2,t,3t), generic fibre degree1 and base degree4.
  PF72's exact fibre order is

      |lcm(e/2,t,3t)*(-4/(eo))|=2.

  The actual R3 longitude is the fibre after M^3 is filled, because
  3 divides eo. Thus R3 has exact order2.
* If A4 first appears at a cable, its older companion action is
  cyclic on every orbit. PF65's actual cyclic companion replacement
  gives the A4 cover of the unknotted torus pattern, plus explicit
  cyclic-branched companion summands. The R3 mark lies wholly in
  the first summand. The preceding fibre-order calculation applies
  to that torus pattern, so it again has order2. The added summands
  can have positive rank; they contribute no marked coordinates.

At any state7 self loop, apply (2) on the active three sheets and
retain the degree1 fixed piece, which performs the other companion
meridian fill. Hence R3_new=r*R3_old and extra homology is unmarked.
Inductively the order of R3 divides2 at EVERY companion depth.

The final7->3 cable has a three-sheet S3 orbit and one fixed sheet.
The actual simple-meridian/trivial-longitude equations force

    r=2u even, 3 divides k, cable image S3.

PF72's integral two-cone transfer applies to its true filled disk
with cones(u,2u): m=M^3 is zero, the ramification mark is r*L,
and the added summand is Z/u. The degree1 piece fills the fixed
companion boundary. Consequently the final simple R mark is ZERO.
Every later state3 tail preserves zero R by PF82 (or (2) with n=2).
PF67's common-source mark is nonzero under its actual smooth exact
whole-line/common-pencil hypotheses, at all finite field degrees
and boundary counts. Thus ALL incoming7 COMMON sources are excluded.

## 3. Complementary first A4 births are permanently branched

Now retain PF77/79's actual smooth exact complementary source with
TWO entire omitted lines and their whole charts. Its marked exterior
is a canonical Seifert PATH, with one-cone internal annuli and a
two-cone terminal disk. Full peripheral recognition, fibre uniqueness,
irreducibility and incompressibility are required as in PF79.

At a first A4 cable birth from state6 or8 the cable action is
transitive, its fibre permutation is identity, and the new periphery
has two boundary orbits(3,1). The older companion has TWO orbits,
(3,1) or(2,2). Cap the new fixed-sheet boundary by its TRUE meridian.
It is a primitive base section, so no exceptional fibre is added.
The remaining piece has THREE boundary circles: two toward genuine
unbranched companion exteriors, one toward the next level/R piece.
All those companion exteriors are nonsolid and boundary-incompressible.
This is a negative-Euler canonical vertex of degree3.

Later cyclic3 levels retain the active boundary annulus, while their
degree1 fixed pieces propagate exactly the fixed meridian cap. They
do not cap either old companion torus at the birth. Original adjacent
fibre determinants are r*B-k with integer B; they are nonzero since
gcd(r,k)=1 and r>=2. Finite-index boundary lifts preserve distinct
rational fibre lines. Thus the branching survives every cyclic3
level, final arrival and later simple tail. It contradicts the
source PATH. Earlier graph cycles/genus already fail PF62/69.

This closes all complementary first A4 births, regardless of how
many cyclic levels came before or after them.

## 4. Initial non-trefoil A4 companions also close

For initial T(e,3t), capping its fixed meridian leaves the actual
one-boundary piece E' with disk cones

    e/2,e/2,t,3t, omitting orders1.                      (3)

If e>2 there are at least THREE genuine cones. This terminal piece
is nonsolid with incompressible boundary and negative Euler base.
The same genuine fibre determinants retain it through every cyclic3
loop and the final arrival/tail. A complementary source terminal
has exactly TWO cones, so these cases close.

Suppose e=2 and t>=3 (necessarily odd). The terminal is disk(t,3t),
again nonsolid and negative Euler. A cyclic3 self loop with 3|k has
an active annulus with THREE cones(r,r,r). It is internal after
the final arrival; this already contradicts the source one-cone
internal budget. All remaining loops have 3 not dividing k.
Their cyclic3 pattern has NO added homology: PF65's torus Alexander
module is trivial modulo t_deck^3-1 when both r,k are units modulo3.
Indeed its Alexander polynomial reduces to1 if either exponent is1
modulo3, and to t_deck^2 if both are2 modulo3. This follows by
evaluating the torus formula at1 and the two primitive cube roots;
the residue degree is at most2, so these evaluations determine it.
Both residues are units. The checker retains independent full CW
controls. Thus these loops
preserve H1(W)=Z/(2t), and keep the terminal disk(t,3t).

If the final arrival has r=2u>=4, its drilled root annulus has two
cones(u,2u). These are not consecutive, unlike the only source
two-cone corner merger. With a later tail this becomes an internal
two-cone annulus, also impossible. The only remaining arrival is
r=2; its root is an annulus with one cone2, R-fibre distance1.
With no later tail the source possibilities have n_R=3 (collar
case) or7 (rho=3,q=2). With any later tail, all even-k tails fail
the two-cone budget, so tails have odd k and add no homology; their
R-fibre distance2 forces the collar case n_R=3. In all cases

    delta=|H1(end)|=2t, n_R in {3,7}.                    (4)

The normalized coefficients at BOTH terminal cones in (3) are
+/-inverse(2), modulo their respective cone orders. The original
T(2,3t) odd-exponent cone has that coefficient. Its lifts have
base degrees3 and1 and fibre degree1, so the residue stays the same
modulo t and3t. Fixed meridian capping does not change it. Use PF85's
retained local covering convention and allow BOTH orientation signs.

Every source free-tail cone has normalized coefficient +/-1 by the
actual final negative-two chain. If its order is kappa, matching
either target cone requires 2=+/-1 mod kappa. Hence kappa=3.
For t>=3 the only possibility is t=3, with source terminal winding
a=9 and free-tail kappa=3. This argument uses normalized residues,
not just cone orders or homology.

The following weaker terminal congruence works at EVERY source depth,
including a single characteristic node:

    n_D=kappa mod a,
    delta=n_R*kappa mod a.                             (5)

For at least two nodes PF83 gives
n_D=a^2 Z+a-(a-1)kappa; a divides p_D, so the full cross square
vanishes modulo a. For one node a=p_D, S_D=(sigma-1)p_D and
n_D=(sigma+1)p_D-sigma+1, while kappa=p_D+1-sigma.
This proves (5) in the previously separate single-node case too.
It requires actual characteristic winding/free-tail assignments.

For a=9,kappa=3, (4)-(5) would say6 is congruent to0 or3 modulo9,
which is impossible. Therefore every initial NON-TREFOIL A4
trajectory closes in the complementary pencil, at every depth and
later tail length, with no prefix eta condition.

## 5. Scope and retained residual

The [checker](verify_cyclic_pattern_a4.py),
[receipt](cyclic_pattern_a4_receipt.json) and
[manifest](sources/cyclic-pattern-a4/manifest.json) retain arbitrary-
degree cyclic-pattern CW controls, actual A4 first-level/loop/arrival
markings, complete state predecessors and incoming root choices,
branching examples, and terminal residue/congruence checks. Finite
controls support the written universal proofs, rather than prove
unbounded statements by extrapolation. Deep primary topology inputs
are imported from the complete retained audits; no fresh visual PDF
audit or independent full reproof is asserted. Related primary web
records were consulted without importing an additional theorem.

The complementary A4 residual starts at T(2,3), with any subsequent
cyclic3 levels and its actual fixed-sheet cap. That capped companion
is a SOLID TORUS. Its cap slope can later create a new nonsolid
terminal, so the retained-terminal argument above cannot be applied
at its starting point. Earlier PF81/83 seed exclusions still hold;
this note does not reopen them or assert that every residual is viable.

Remaining: incoming2(S3), complementary incoming7 trajectories of
this trefoil origin outside prior closed seeds, and initialstate3
torus trajectories outside the closed intrinsic ranges. PF85/86's
cyclic6/9 and S4 incoming4/5 closures remain. No full quartic,
higher-degree, degenerate-stage, finite Keller-map, formal verification,
external-review or novelty/priority claim is made. Full goal OPEN.
