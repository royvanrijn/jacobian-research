# Global attachment changes the local normal-r7 excess gate

13 November 2026. PF2BATT1 is a consequence of the external plane
Stein-boundary theorem and branch purity. PF2R7BND1 is its application
to the two retained nonimmersion classes, with explicit ramified-divisor
hypotheses. Work over C. No Tang structural theorem is used.

## 1. The external theorem and its precise quantifiers

Use Alexander Borisov, *On the Stein factorization of resolutions of
two-dimensional Keller maps*, Beitr. Algebra Geom.56(2015),299-312,
[DOI10.1007/s13366-014-0208-4](https://link.springer.com/article/10.1007/s13366-014-0208-4),
published online 7 June 2014. The pinned
[author manuscript](https://people.math.binghamton.edu/borisov/documents/papers/ampleramification.pdf)
has a different title and uses Theorem3.4 and Proposition3.3 for the
boundary structure. Its Theorem3.1 asserts **existence** of a ramified
dicritical, not ramification of every dicritical.

The imported structure is this: in the projective finite Stein surface,
the source-boundary curves meet only at one point over target infinity.
Each affine boundary normalization is A1. Each such curve is unibranch;
it is smooth wherever the surface is regular. A boundary curve has at
most one further exceptional point, a cyclic quotient surface singularity.
The manuscript proves this via a minimal source boundary tree and
contraction of its side chains. We use this stated theorem; the checker
does not reprove that global classification.

Let F:A2 -> A2 be Keller and pi:N -> A2 its finite normalization in the
source function field. This is exactly the affine restriction of that
projective Stein surface, by uniqueness of normalization. Write
U=A2_source inside N and B=N\U. The imported theorem therefore applies
to all of B, including generically unramified boundary primes.

For clarity about unibranchness, a finite point with two preimages on
one boundary normalization would require a connected exceptional fibre
joining two points of its strict transform. Together with that transform
this is a cycle in the smooth A2 boundary tree. Different boundary
curves cannot meet in N either: the imported intersection point lies
over infinity and has been removed. Normalization of each boundary
curve is a bijection A1 -> D, hence a homeomorphism in the complex
topology, so chi(D)=1 even if its unique possible surface-singular
point is present. Consequently chi(B) is the number of boundary primes.

## 2. What a weakly unramified omitted prime can do

Suppose a boundary prime D has generic normal ramification index e=1.
All ramification primes of pi lie in B, because pi restricts to the
etale source map on U. Boundary primes are disjoint. Thus at any point
q of D, every codimension-one generization is unramified: its generic
point is either outside B, on D, or on no other boundary prime.

Apply [purity of branch locus, Stacks tag0BMB](https://stacks.math.columbia.edu/tag/0BMB)
at q. The source local ring is normal, the target local ring is regular,
pi is finite and has equal local dimensions, and all height-one
generizations are unramified. It follows that pi is etale at q.
In particular N is regular along D; the imported boundary theorem then
makes D smooth. Its inclusion in N is an immersion, so its composite
map into A2 is immersive as well.

Let C_D=pi(D). This is a component of the nonproper-value set. Its
normalization is A1: the nonproperness theorem gives a polynomial
parametrization; a nonconstant map A1 -> C_D lifts to its rational
normalization, which can have at most one puncture. The finite map
D -> C_D lifts to a finite map D -> C_D^nu. It is unramified because
its composite into A2 has nonzero differential. Since D is A1, the
resulting polynomial self-map of A1 has nowhere vanishing derivative,
and in characteristic zero it is linear. Thus

\[
 \boxed{e_D=1\ \Longrightarrow\ f_D=1,
 \quad D\simeq\mathbb A^1,\quad D\longrightarrow\mathbb A^2
 \text{ is immersive}.}\tag{U}
\]

Every local surface fibre at q on D has length one. Such a boundary
prime cannot have as image a curve whose normalization is nonimmersive
at a finite parameter. This conclusion is specific to a globally
realized plane Keller boundary. The previously retained local cubic
opens that delete an unramified companion do not have source A2 and
are not counterexamples to (U).

The polynomial parametrization input is Jelonek's nonproperness theorem,
as pinned in [the pointed-pencil audit](../../hessian/hc4-pointed-pencil/sources/manifest.json).
We have not assumed every nonproper component is a branch component
of the finite map pi.

## 3. The collided tame-node model is excluded globally

For the retained normal local cover

\[
 \{uv=w^n\}\longrightarrow\mathbb A^2_{u,v},\qquad n\ge2,
\]

there are two different analytic ramification branches through the
origin: (u,w)=0 and (v,w)=0. At the general point of either branch
the normal ramification index is n. Both must be omitted from U in
a Keller realization. If their two germs belong to different global
boundary primes, those primes intersect in N, which is forbidden.
If they belong to one global boundary prime, that prime has two
analytic branches there, also forbidden. Hence this finite map germ
cannot occur in the finite normalization of a plane Keller map.

This closes the potential global attachment of every collided node
control from [HC4MYRB2](../../hessian/hc4-radial-intersection/PROOF.md), including
the odd-n extra -1 radial horizontal. It does not exclude the abstract
local normal covers. It also does not exclude an A_(n-1) surface
singularity under another projection: the finite map, not merely
the surface's analytic type, is the excluded object.

More generally, at a finite normalization point at most one analytic
ramification branch can be present. Thus the S3 local assignment for
two different smooth target branches cannot identify their ramified
normalization points on the same global boundary prime. This will
remove the old purely local A5 excess option below.

## 4. Application hypotheses for the two r7 cusp classes

Assume now that F has geometric degree six, and that the **entire
reduced ramification divisor of pi** is one boundary prime R with
generic (e,f)=(2,1). Let C=pi(R) be one of the two birational
polynomial curves classified by
[PF2R7CUSP1](../jc2-normal-r7/README.md):

\[
 \mathrm{I}: A_2+A_5+8A_1,\qquad
 \mathrm{II}: A_2+11A_1.
\]

The retained normal-r7 ramified row motivates this hypothesis. The
statement here explicitly audits a single **ramified** row; it does
not infer a complete missing-boundary ledger from it. All additional
boundary primes, if present, are generically unramified and are allowed.
Their images are not assumed to equal C or to share its infinity point.
This distinction is essential when applying the result to a proposed
global F2 cover.

Each C has normalization A1 and just one nonimmersion parameter, the
ordinary cusp. Row I has nine two-branch singular values (eight nodes
and the A5); row II has eleven. Both branches of each two-branch
singularity are smooth and correspond to distinct points of C^nu.

The map R^nu -> C^nu has degree one. Its points above those two
normalization parameters cannot identify in N, since R is unibranch.
Their local ramified orbits are therefore separate. Since R is the
only ramification prime, each such orbit consists of exactly two
sheets: at its point just one smooth target branch carries the simple
ramification. The other target branch acts trivially on that orbit.
Thus every two-branch value has two double points plus two simple
points in the complete normalization fibre: four distinct N-points.
In particular the A5 cannot contribute the locally permissible S3
collision from the earlier excess-only calculation.

At the ordinary cusp, two meridians have transposition images a,b
satisfying aba=bab. Equal transpositions give a two-sheet local orbit;
distinct allowed ones share one letter and give a three-sheet orbit.
The cusp has only one point on R^nu, so only one such ramified orbit.
Put epsilon=0 in the first case and epsilon=1 in the second. The
complete normalization fibre has 5-epsilon points. No fibre is being
identified here with the fibre of the original affine source.

## 5. A complete Euler calculation forces the missing primes

Let s be the number of two-branch values: s=9 or 11. Normalization
identifies two points at each such value and one point at the cusp,
so

\[
 \chi(C)=1-s,\quad
 \chi(C\setminus\mathrm{Sing}C)=-2s,\quad
 \chi(\mathbb A^2\setminus C)=s.
\]

Away from C the **finite normalization map** is etale of degree six.
Over its smooth points it has five distinct points. The complete
Euler stratification, counting the actual finite fibres, is

\[
 \chi(N)=6s+5(-2s)+4s+(5-\epsilon)=5-\epsilon.\tag{E}
\]

If there are t additional unramified boundary primes, all t+1 boundary
curves are disjoint and have Euler characteristic one. Since U=A2,

\[
 1=\chi(U)=\chi(N)-(t+1),\qquad
 \boxed{t=3-\epsilon\in\{2,3\}.}\tag{M}
\]

By (U), every one has residue degree one and an immersed image. Since
C's normalization is nonimmersive at its cusp, none can have image C.
At least one additional nonproper-value component is therefore forced.
Two or three **boundary primes** are forced; their images need not be
distinct from each other. We do not count them as two or three distinct
target components.

If the full nonproper-value locus is asserted to be just C, (M) and
(U) are a contradiction. This excludes that global stratum for both
rows, while allowing arbitrary unramified rows over C at the outset.
It strengthens the old conditional one-row excess test: row I's A5
option fails globally, and omitted unramified primes over the same
cusp curve cannot repair either row.

This is also consistent with the complete Orevkov budget when that
formula is used: e*f contributes 2+t and the sole remaining excess
is epsilon, giving 5=2+t+epsilon. The proof above derives (M) directly
from finite-normalization Euler characteristic and does not import an
unproved complete one-row deletion ledger.

## 6. What has not been proved

The argument does not exclude unramified boundary primes whose immersed
images are other curves. In particular, knowing the branch locus of
the finite map does not determine the full nonproper-value set of its
distinguished affine open. The terminal cubic's ramification
classification does not by itself exclude those additional images.
Their global infinity positions and compatibility with the carrier
remain a required audit, not a completed consequence of the passport.

The immersed normal-r7 family, normal r5, nonnormal conductor slices,
the whole (75,125) row and JC2 remain open. This note does not claim a
new coordinate-degree bound or a global counterexample.

The finite checker enumerates all S6 transposition pairs for the cusp
and the two smooth-branch collision test, verifies the Euler ledger for
both exact singularity types, and checks the two monomial ramification
branches of n=2,...,12 node controls. It reads the retained classification
data without regenerating the large r7 elimination. The global theorem,
purity application and topological arguments are written mathematical
inputs; finite calculations do not certify them. No formal verification,
external review of this application, independent reproof of Borisov's
full theorem or novelty priority is claimed.
