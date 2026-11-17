# The monodromy orbit tree and an unbounded eta0 winding gate

15 November 2026. PF2RULED62. Programme archived; full goal active.

Every finite branched end cover supplied by a smooth affine A2 extension
has finite H1, by PF2RULED61. Cutting the target along cabling tori gives
a cheaper all-degree obstruction: the graph of lifted pieces must be a
TREE. Its vertices and edges are computed solely from permutation-group
orbits. This applies at every finite degree and every cabling length,
without coefficient bounds or a classification of symplectic fillings.

For the entire primitive-binomial two-stage class with first torus pair
(2,3), the remaining quartic eta0 condition fixes the last cable parameter
to w+1. Exact S4 exponent periodicity reduces ALL windings to twelve
finite cases. Combined with the sign-surface determinant-square theorem,
a surviving winding must satisfy

    w=3 modulo6,  OR  w=8 modulo18.                    (1)

Even in these classes the orbit-tree condition removes additional
representations. It does not prove an affine realization for any survivor.

## 1. A graph quotient of a branched cover

Let a link exterior in S3 be cut along disjoint tori into connected
pieces P_v, with base decomposition graph a finite tree Gamma. Iterated
torus cabling is one case. Take a connected degree-d branched cover M
with transitive sheet action rho on Omega={1,...,d}. No regular-cover,
simple-branch or symmetric-group hypothesis is needed in this section.

Choose base paths through the target tree. Let H_v be the image of
pi1(P_v), and H_e the image of the peripheral torus group at an edge.
The inclusions H_e into the adjacent H_v use these SAME base paths;
arbitrary independently conjugated permutation lists are not incidence
data. Connected covering theory gives

    vertices above v = H_v-orbits on Omega,
    edges above e    = H_e-orbits on Omega.             (2)

An edge orbit lies in one orbit of each adjacent vertex group, which
specifies its incidence. Let G_rho denote the resulting finite multigraph.
An ordinary covering of the base tree is NOT being asserted: different
vertex and edge orbit sizes allow parallel edges and cycles upstairs.

Fill all boundary components over the branch link. Each filling attaches
a solid torus to a boundary torus of ONE lifted piece. It does not join
different pieces, so it leaves G_rho unchanged. This includes every
meridian cycle, including fixed unramified sheets. The lifted cutting
tori are disjoint closed tori in M. They need not be incompressible or
form a minimal JSJ decomposition for the following argument.

Collapsing the connected piece spaces in the graph-of-spaces description
gives a surjection of fundamental groups onto pi1(G_rho). Equivalently,
the H0 part of Mayer--Vietoris has the exact sequence

    H1(M;Z) -> Z^(edges G_rho) -> Z^(vertices G_rho) -> Z -> 0.

The middle map is the incidence matrix. Thus G_rho is connected and

    H1(M;Z) ->> Z^g,
    g=sum_e orb(H_e)-sum_v orb(H_v)+1.                 (3)

Here orb(H) is the number of its orbits on ALL d letters. A cycle
cannot be removed by meridian fillings inside the individual pieces.
No surgery slope or normal framing is discarded in making this claim.

Now suppose this is the actual boundary cover of a smooth affine finite
normalization containing A2 with omitted disjoint affine lines. The
whole integral source lattice in PF61 gives finite H1(M). Consequently

    g=0: the complete monodromy orbit graph MUST be a tree.     (4)

This is an all-degree and all-length necessary condition with smoothness
explicit. It does not assert smoothness for arbitrary higher-index
normalizations. Passing (4) leaves all vertex homology, gluing slopes,
the full linking form and the canonical knot marking to be checked.

## 2. Useful local consequences and a multistage warning

For a single cabling torus, write A for the companion image, B for the
cable-space image and P for the common peripheral image. Then

    g=orb(P)-orb(A)-orb(B)+1.                           (5)

If both A and B act transitively, (4) requires P to be transitive too.
But P is abelian. A transitive abelian permutation group acts regularly:
an element fixing a letter fixes every letter, hence is the identity.
Therefore every element of P is semiregular, with all cycles of the
same length. In particular a companion meridian with MIXED cycle
lengths forbids simultaneous transitivity of the two pieces in every
degree. This is stronger than global transitivity of the whole cover.

In degree four, simultaneous transitivity requires the companion
meridian to be a 4-cycle. The mixed types (2,1,1) and (3,1) are excluded
by the preceding argument. The identity cannot normally generate a
nontrivial companion image. For type (2,2), normal generation by the
meridian makes the companion image a subgroup of A4. Its transitive
possibilities are V4 and A4. In V4 a double transposition normally
generates only C2, and in A4 it normally generates at most V4. Neither
can be the whole image of a knot group whose meridian normally generates
it. Thus type (2,2) is impossible for a transitive companion image.
No claim that every higher-degree semiregular meridian is a d-cycle is made.

At multiple stages, check the WHOLE orbit graph, not each edge separately.
A group-theoretic control over the base path A--B--C has transitive C2
images at A and C and trivial images at B and both edges, on two letters.
Each one-edge graph is a tree: one endpoint joined to two B lifts.
Together they form the four-cycle A--B1--C--B2--A. The total g is1.
This is an orbit-incidence control, not a claimed knot presentation.
It shows why separately passing one-edge counts cannot replace (3).

## 3. All windings in the primitive-binomial trefoil class

Suppose an actual surviving quartic branch has the immersive polynomial
normalization of PF2RULED60 with degrees 3w,2w, primitive next/final key
F=X^2-Y^3 of degree k<6w and gcd(k,w)=1, w>=2. These are explicit
actual-binomial hypotheses. Lower-weight corrections to a general key
are not dropped. PF60 gives eta=k-w-1>=0, while PF59 gives eta<=0 for
the whole quartic. Hence

    k=w+1, eta=0, K=C_(w,w+1)(T2,3).                   (6)

Keep the positive peripheral convention already proved in PF61:

    <a,b | a^2=b^3>, M=a*b^-1, L=a^2*M^-6,
    h=M^(w+1)*L^w, z^w=h, beta=M^-1*L^-1,
    mu=z*beta, lambda=h*mu^(-w*(w+1)).                (7)

The actual monodromy has image S4, mu a single transposition and lambda1.
Since w*(w+1) is even, rho(h)=1. Thus all such representations are
exactly the permutations a,b,beta,tau satisfying

    a^2=b^3, M=a*b^-1=beta^w,
    L=a^2*M^-6=beta^(-w-1),
    (tau*beta)^w=1, tau a transposition,
    <a,b,beta,tau>=S4.                                (8)

S4 has exponent12. Therefore the ENTIRE representation set (8), not
merely its cardinality, depends only on w modulo12. The companion,
cable and edge subgroups are respectively

    A=<a,b>, B=<beta,tau>, P=<M,L>.                   (9)

Their actual permutations also depend only on that residue, so the
complete orbit graphs do too. This exact finite-period reduction is
the unbounded proof; no inference is made from testing finitely many
large windings.

Exhausting all permutations in (8), and simultaneous conjugation, gives
the table below. Every orbit has24 representations; the trivial
centralizer follows also from the full S4 image.

| w mod12 | representations | conjugacy classes | tree classes | cycle classes |
|---|---:|---:|---:|---:|
| 2 | 72 | 3 | 3 | 0 |
| 3 | 72 | 3 | 3 | 0 |
| 8 | 168 | 7 | 4 | 3 |
| 9 | 144 | 6 | 6 | 0 |
| all other residues | 0 | 0 | 0 | 0 |

The residue2 row is essential: its smallest winding2 fails the sign
determinant-square test, but larger windings in that SAME representation
row need not fail. Discarding an entire residue from its smallest
example would be an invalid extrapolation.

For w=8 modulo12, each rejected class has A=A4 and B=S4, both transitive,
while P=C3 has a three-letter orbit and a fixed letter. Thus g=2-1-1+1=1.
It is excluded for EVERY w in this residue and EVERY smooth affine source
in Section1. At w8 this explains the three Z+Z/26 computations in PF61;
their torsion26 is not asserted to persist at larger winding. The four
tree classes are retained, including those with an intransitive cable
image. The argument does not replace those subgroups by S4.

## 4. Combine the infinite residue test with the actual sign surface

The PF2RULED55 sign-surface theorem requires 9 to divide the target knot
determinant for every actual surviving quartic. The cabling Alexander
formula of PF54 gives, in (6),

    D_K=w+1 if w is even; D_K=3w if w is odd.          (10)

For odd w the representation table leaves w=3 or9 modulo12, exactly
w=3 modulo6, which already makes 9 divide3w. For even w the table leaves
w=2 modulo6. Requiring 9|(w+1) then gives w=8 modulo18. This proves (1).

More precisely, the even surviving classes split as

    w=8 modulo36:  seven peripheral classes, only four pass the tree test;
    w=26 modulo36: three peripheral classes, all three pass the tree test.

The odd classes have three or six tree representations according as
w=3 or9 modulo12. These are NECESSARY infinity gates. No polynomial
realization for arbitrary w, sign Picard certificate, affine monodromy
extension, finite quartic normalization or marked plane source is supplied.
The w3 six-node curve and its remaining source-sector comparisons in
PF61 are retained, not generalized to all w by numerical analogy.

## Remaining obligation and rejected shortcuts

The orbit graph is only the graph part of end homology. Tree candidates
must still match the integral source lattice, full linking form, canonical
ramification marking, affine sign Picard torsion, and the whole finite
target algebra. The theorem applies to higher degrees after actual
smoothness and finite-cover data are supplied; it does not manufacture
an S3 quotient for S_d at d>=5.

Literature searches on rational affine surfaces and symplectic fillings
did not justify assuming a unique boundary contact structure or identifying
a boundary manifold from H1. No such classification is used. Exploratory
coset enumerations reaching a cap are not evidence that a boundary group
is infinite. The graph proof needs neither those enumerations nor named
three-manifold identifications.

[Checker](verify_orbit_tree.py) and
[manifest](sources/orbit-tree/manifest.json) distinguish the written
all-degree graph proof, exact exponent-period reduction, finite exhaustive
S4 certificate, and retained primary cabling conventions. No full quartic
or higher-degree exclusion, finite Keller pair, formal verification,
independent full reproof, external review, novelty priority or goal
completion is claimed.
