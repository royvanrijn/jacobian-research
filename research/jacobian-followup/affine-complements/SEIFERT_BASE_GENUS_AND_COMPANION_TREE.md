# Seifert base genus and the companion permutation tree

15 November 2026. PF2RULED69 strengthens the all-degree orbit-tree
obstruction of PF62. In a finite branched end assembled from lifted
Seifert pieces, handles in the bases of those pieces survive EVERY
boundary-torus gluing and filling. A smooth affine plane extension
with the PF61 disjoint-line boundary hypotheses therefore requires
both the entire lifted-piece graph to be a tree AND the base of
every lifted Seifert piece to have genus zero.

This excludes 32 further periodic representation classes in PF63's
complete quartic census. All 32 have tree gluing graphs but a genus-one
S4 companion piece, hence end first Betti number at least2. They are
excluded at every actual coprime exponent lift and winding in their
residue classes. Other representations of the same target keys remain;
no whole residue cell or remaining noncyclic sector is declared closed.

## 1. An integral quotient that no boundary gluing removes

Let a compact connected orientable 3-manifold be obtained by gluing
connected pieces along boundary tori and filling other boundary tori
by solid tori. Some pieces P_v are Seifert fibred with ORIENTABLE
base surface B_v of genus g_v. Other pieces may be arbitrary. Let G
be the connected graph of the pieces and cutting tori, before the
internal solid-torus fillings; such fillings do not change G.

For each selected Seifert piece the projection to its underlying
base, followed by capping all base boundary circles, induces

    H1(P_v;Z) ->> H1(B_v capped;Z)=Z^(2g_v).                 (1)

To check surjectivity, choose standard handle loops avoiding the
exceptional values. Each lifts to a path whose endpoints can be
joined along its connected circle fibre. This gives closed lifts
of all the handle generators. Every boundary torus maps through a
boundary circle, which is zero after capping. Thus (1) kills the
image of the ENTIRE H1 of every boundary torus, not just its meridian.
Every filling slope is killed as well, even a fibre slope. No claim
that such an arbitrary filling preserves a Seifert fibration is needed.

The first-homology Mayer--Vietoris sequence has the form

    C -> H1(M;Z) ->> H1(G;Z) -> 0,
    C=coker(direct-sum H1(edge tori) -> direct-sum H1(pieces)).

More precisely C injects into H1(M). Fillings can be included as
additional relations from their boundary tori; their core generators
can be eliminated using a complementary boundary slope. By (1), C
surjects onto the direct sum of the capped-base homology groups. The
last displayed exact sequence splits as a sequence of abelian groups
because H1(G) is free. Therefore there is an integral surjection

    H1(M;Z) ->> Z^(b1(G)+2*sum_v g_v).                      (2)

The splitting is not canonical; the rank inequality and the base
quotients are independent of that choice. In particular

    b1(M) >= b1(G)+2*sum_v g_v.                             (3)

No incompressibility, minimal JSJ decomposition, regular covering,
cyclic companion restriction, simple branching, or prescribed
gluing slopes are assumed in this lemma. Orientability of the
selected bases IS assumed. The torus-knot and cable bases used below
and all their lifted bases are orientable.

PF61 gives finite H1 for the actual end of a SMOOTH affine N containing
A2 with omitted DISJOINT entire A1s. Under those explicit hypotheses,
(3) forces G to be a tree and all g_v=0, at ALL finite field degrees
and cabling lengths. Smoothness in higher degrees is not automatic.
The claim is necessary, and does not supply an affine cover or source.

Hatcher's standard construction from a circle bundle over a punctured
base and solid-torus attachments, complete PDF33-34/printed32-33, is
the primary model used here. The Mayer--Vietoris quotient and its
survival under arbitrary boundary operations are proved above.

## 2. Compute the base on fibre orbits, not individual sheets

Take one connected sheet orbit Omega_v of a finite cover of an
orientably Seifert-fibred target piece. Let H be the permutation of
an oriented regular fibre. It commutes with the piece action. Its
cycles on Omega_v all have the same length m: transitivity carries
any one cycle equivariantly to any other. These are the lifts of
one generic target fibre. Set

    F_v = Omega_v / <H>, n=|F_v|=|Omega_v|/m.                (4)

The induced action on F_v is the monodromy of the cover of base
orbifolds, of degree n; the lifted fibre degree is m. This follows
also by lifting the local product S1 times a disk away from cone
values. Over a cone value, completing the punctured local cover
gives a branched disk cover on bases, with cycle lengths from the
induced cone permutation. This proves the needed covering assertion
directly from the fibre models, including nontrivial H.

Suppose the target underlying base has genus h, b boundary circles
and j cone points. Write sigma_i for the induced monodromies around
cone points and delta_l for those around boundary circles. For a
permutation pi, c(pi) counts ALL its cycles, including fixed letters.
Riemann--Hurwitz for surfaces with boundary gives

    chi(B_v)=n*(2-2h-b-j)+sum_i c(sigma_i),
    |boundary B_v|=sum_l c(delta_l),
    2g_v=2-n*(2-2h-b-j)-sum_i c(sigma_i)-sum_l c(delta_l).    (5)

One can derive the first equality without a branched formula: remove
small disks about the j cone values, use multiplicativity of Euler
characteristic in the resulting unbranched cover, then attach one
disk for each cycle over each removed value. Formula (5) is valid
at arbitrary degree. Powers and fibre factors in the ACTUAL peripheral
words must be retained before passing to F_v.

For a torus companion <a,b|a^p=b^q>, take H=a^p=b^q,
M=a^u b^v with qu+pv=1 and L=H M^(-pq). On F_v, a and b have
orders dividing p,q. The base cone loops are a^u,b^v; their cycle
counts equal those of a,b because gcd(u,p)=gcd(v,q)=1. The unique
base boundary loop is M. Therefore, using induced permutations,

    2g_v = n+2-c(a)-c(b)-c(M).                              (6)

For a cable piece with base an annulus and one cone value,
the three terms are the two base boundary loops and its cone loop;
the right side is again n+2 minus their three cycle counts.

The distinction (4) is essential. In the cyclic4 representation of
T(2,3), let M=(0123), a=M^3, b=M^2. Then H=M^6=M^2 has two cycles,
so n=2. The induced a,M are transpositions and b is the identity:
(6) gives genus0. Applying (6) to the FOUR individual sheets instead
would give genus1, falsely excluding PF68's finite-H1 C4 control at
p=2,q=3,r=3. Nontrivial fibre monodromy cannot be discarded.

## 3. An all-degree companion permutation tree

Consider a two-piece torus-companion/cable decomposition at arbitrary
degree d. Suppose BOTH piece actions are transitive on all d sheets
and the companion image has trivial centre. The regular fibre H is
central in the companion group, so its image is the identity. The
peripheral image is consequently generated just by M, since
L=M^(-pq). PF62's entire orbit-tree condition now implies that the
edge action is transitive, so M is a d-cycle.

There are n=d base sheets and c(M)=1. The new genus-zero requirement
is the unrestricted identity

    c(a)+c(b)=d+1,
    index(a)+index(b)=d-1, index(pi)=d-c(pi).                (7)

Equivalently form a bipartite graph with one black vertex for each
cycle of A=a^u, one white vertex for each cycle of B=b^v, and one
edge per sheet joining its two cycle vertices. Unit powers imply
<A,B>=<a,b>, so this graph is connected. It has d edges and
c(a)+c(b) vertices. Thus (7) says precisely that it is a TREE.
This is the usual permutation graph underlying a dessin; the
elementary vertex/edge argument is all that is needed here.

The result applies, for example, to a full Sd companion image for
d>=3 and a transitive adjacent cable image. The hypothesis is on
the ACTUAL companion image; full global monodromy alone does not
give it. In longer decompositions, apply (5) separately on each
actual orbit of each vertex, not an independently conjugated piece.

For quartic S4/S4 pieces, (7) leaves only a transposition and a
3-cycle among noncyclic torus generators of coprime orders. If one
generator is a 4-cycle and the other a 3-cycle, c(a)+c(b)=3 rather
than5, so the base has genus1 and the end has b1>=2. This excludes
the WHOLE generator-cycle-type sector at arbitrary p,q,r, not all
S4/S4 representations. The alternative involving a double
transposition and a 3-cycle generates A4, not S4.

## 4. Complete quartic residue census and independent controls

PF63's exhaustive certificate has 596 periodic conjugacy classes in
1152 residue cells: 112 already fail the global graph-cycle gate and
484 have tree graphs. Applying (4)-(6) to that retained census gives

| Companion/cable images | Tree classes, base genus0 | Tree classes, base genus1 |
|---|---:|---:|
| C3/S4 | 84 | 0 |
| C4/S4 | 96 | 0 |
| A4/S3 | 96 | 0 |
| S3/S3 | 64 | 0 |
| S4/S4 | 64 | 32 |
| S4/C2 | 48 | 0 |
| Total | 452 | 32 |

Every cable base has genus0. The new32 classes occupy eight residue
cells, four classes per cell:

    (p,q,r) mod12 = (3,4,3), (3,8,3), (4,3,3), (4,9,9),
                   (8,3,3), (8,9,9), (9,4,9), (9,8,9).

No entire formerly nonempty tree cell becomes empty. These counts
describe the genus gate, before combining other PF65-68 exclusions;
they are not counts of actual affine covers or of surviving Keller
maps. Keep the sealed PF63 certificate unchanged.

There is no exponent bound behind the gate: permutations in S4 have
exponent12, so H=a^p, M=beta^r and every induced quotient action are
fixed in the corresponding residue cell, for every actual coprime
lift. Coprimality of a residue pair alone is not asserted to imply
coprimality of every integer lift.

As a full group-presentation control, p=3,q=4,r=3 gives H1=Z^2
for a new excluded class. The same representation has H1=Z^2 at
(3,16,3), Z^2 direct-sum Z/5 at (15,4,3) and (3,4,15), consistent
with (3). These are controls; the unbounded exclusion comes from
the proved quotient. The genus-zero p=2,q=3,r=9 S4/S4 control
still has H1=Z/24 and nonzero R mark (quotient Z/6); it is retained
as an open-sector control, not excluded by this result.

[verify_seifert_base_genus.py](verify_seifert_base_genus.py) checks
all induced fibre quotients, companion and cable genera in the
complete residue census, the all-degree symbolic genus formula,
boundary-insensitive handle columns in a Seifert presentation, and
full satellite presentations with independent spanning trees.
[seifert_base_genus_receipt.json](seifert_base_genus_receipt.json)
pins this evidence and [the source manifest](sources/seifert-base-genus/manifest.json).

The remaining even cyclic4 complementary sector and genus-zero
noncyclic representations still require genuine integral peripheral
homology, marked linking and finite affine-cover compatibility.
General keys and longer sequences are not assumed binomial, and
higher-degree smoothness is not assumed automatically. PF63's
polynomial immersions still exist. No full quartic/higher-degree
exclusion, finite Keller map, independent full reproof, formal
verification, external review, or novelty/priority claim is made.
