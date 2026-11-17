# Rational source ends and the entire trefoil common-source exclusion

15 November 2026. PF2RULED74 proves a SOURCE theorem with no field-degree,
ramification-order, omitted-line-count or cluster-depth bound, under the
explicit common-pencil hypotheses below. Its end is, up to orientation,
a connected sum of rational negative-definite plumbing links. In particular,

    EVERY homomorphism pi1(end N) -> L is trivial
    for every left-orderable group L.                         (1)

This is stronger than saying the source group is not left-orderable;
that statement alone would not rule out its having a left-orderable quotient.
We do NOT assert that all these source groups are finite or cyclic.

Every remaining trefoil S3/S3 target admits a NONTRIVIAL left-orderable
Seifert quotient. Thus its common source is impossible. With PF68/71/72,
the ENTIRE actual primitive-binomial TWO-STAGE QUARTIC class with first
pair (2,3), in either order, is now excluded in the COMMON source pencil,
at all windings, source orders and poles. Complementary sources remain
open where previous gates permit them. This is not a full target-key or
quartic exclusion. PF63 polynomial immersions still exist.

## 1. The all-degree source hypotheses and reduction

Let N be a smooth affine complex surface containing U=A2_(x,f), with
N minus U a finite disjoint union of ENTIRE affine lines E_i. Assume
f is regular on N and constant, with finite value lambda_i, on each E_i.
Assume x:N->P1 is regular and has a pole along every E_i. These are
PF67's common-coordinate hypotheses, without imposing any volume order.
They are not automatic for arbitrary higher-degree packets.

Use PF67's completion and factorization outside U:

    X -> P1_x times P1_f

is a sequence of boundary point blowups, Gamma=X minus U is the full
SNC boundary, and each retained closure E_i is a free LEAF. The initial
sections H_x={x=infinity}, H_f={f=infinity} belong to T=X minus N.
All retained valuations are centred at points (infinity,lambda_i).

Prune unnecessary blowups in reverse order: retain exactly the centres
on the ancestry of the E_i. Every removed exceptional curve belongs to T;
its contraction is the inverse of a boundary blowup and changes neither
N nor its end. The image is still outside N: over f=infinity it lies on
the remaining infinity fibre, and over a finite value it lies on the
remaining boundary ancestry or H_x. All blowups over f=infinity can be
removed since no E_i is centred there. Thus H_f has weight0 and is a
leaf attached to H_x. All other vertices are exceptional over finitely
many smooth points of H_x. No original finite f-fibre is a boundary vertex.

Apply Neumann Proposition2.1 R6 to this zero end and its H_x neighbour.
Their genera and puncture counts are0. The graph is a TREE, so each
remaining connected branch has exactly ONE attaching edge; the extra
S1xS2 count 2g+sum(k_j-1) is0. No unmarked summand is dropped.
The resulting graph is the disjoint union of the exceptional trees
with the retained E_i leaves removed. Neumann's convention for a disjoint
plumbing graph is CONNECTED SUM (printed303), not disjoint union of
three-manifolds. Consequently, with M(G) the oriented plumbing boundary,

    end(N) = orientation reversal of #_j M(G_j).             (2)

Each G_j is a subgraph of the exceptional graph of point blowups over
a smooth surface point. Its intersection matrix is negative definite,
and its vertices are rational. An empty graph contributes S3.
These statements use the full boundary, not only the omitted-line Gram.

## 2. Rationality at arbitrary depth, without a bounded graph census

For an effective exceptional cycle Z put chi(Z)=-(Z^2+K.Z)/2, with
K given by adjunction. We prove chi(Z)>=1 whenever Z is nonzero,
on ANY sequence of point blowups over a smooth point.

For one further blowup pi with exceptional F, write

    Z=pi^*(Z_old)+aF, a an INTEGER,
    K_new=pi^*(K_old)+F.

Here a can be NEGATIVE, since total transforms and strict transforms
are different bases. Orthogonality and F^2=-1 give

    chi_new(Z)=chi_old(Z_old)+a(a+1)/2.                    (3)

The last term is nonnegative for EVERY integer a. If Z_old is nonzero,
induction applies, since pushforward preserves effectivity. If Z_old=0,
effectivity and Z nonzero force a>0, and the triangular term is at least1.
This proves the assertion by induction. A cycle supported on a subgraph
has the same form and adjunction pairings, so the assertion restricts to
every G_j. In particular it holds for its nonzero fundamental cycle.
Artin's rationality criterion, in Nemethi Section2.1, therefore says
that every G_j is rational. This is an unbounded proof, not a finite
Laufer-algorithm search or an affineness claim for arbitrary clusters.

Import Nemethi Theorem1.6: for a connected negative-definite rational
plumbing graph, a NONTRIVIAL link group is not left-orderable. The link
is irreducible; as an orientable irreducible closed three-manifold it
is P2-irreducible. Boyer--Rolfsen--Wiest Theorem1.1(1) then shows that
EVERY homomorphism from that group to any left-orderable group is
trivial. Otherwise that theorem would make the source group itself
left-orderable. The S3/trivial-group case has the same conclusion directly.

By van Kampen, (2) gives a FREE PRODUCT of the factor groups. A homomorphism
from a free product which is trivial on every factor is trivial. This
proves (1), including reducible source ends. No irreducibility of the
whole source end or target end is assumed. Orientation reversal does
not change fundamental groups. Equivalently the factors are rational
L-space links, but an L-space label alone is not used as a quotient theorem.

## 3. A Seifert group quotient for EVERY S3/S3 target

Keep PF73's actual e=2b,o=3t,r=3s, A=bt,k=3s(6bt-2b-3t)+1,
D=k/2-As. The actual main companion E' has

    H1(E')=Z*m plus Z/b, L=-Am.

Its abelianization followed by projection onto Z maps m->1,L->-A.
This is compatible with a solid torus whose meridian is eta=L+A*m.
The degree1 companion knot exterior likewise maps onto a solid-torus
group with its preferred longitude killed. PF73's pushout lemma gives
a SURJECTION of ACTUAL end groups onto the closed Seifert model

    pi1(actual M) ->> pi1(W), W=P'(eta).                    (4)

This simultaneously abelianizes both companion leaves; it is not a
homeomorphism or group equality. No incompressibility is needed.
P' has base D and cones s,k. The eta filling has third order

    |B|, B=k/2-rA,                                       (5)

since h=(k/2)m+rL and det(h,eta)=k/2-rA. It is nonzero:
2B=3s(4bt-2b-3t)+1 cannot vanish. Also H1(W)=Z/D by PF73.
Thus a nontrivial left-orderable pi1(W) excludes ANY common source
with the hypotheses in Section1. This is a new reusable three-cone
test for other first pairs. We do not classify its rationality for
all b,t,s here, or infer cone twists from homology when orders are
not pairwise coprime.

## 4. Every remaining trefoil model is left-orderable

For the trefoil b=t=A=1, PF73 already forces

    D odd, D=3 mod8, D>=11, s=2D-1, k=6D-2.

In (5), B=-(3D-2). Write J=3D-2=(3s-1)/2. The actual quotient
W has base S2 and three cone orders s,k,J, with H1(W)=Z/D.
The orders are PAIRWISE COPRIME: gcd(s,k)=gcd(s,J)=1, and
gcd(k,J)=gcd(2,J)=1 because D is odd.

Choose the orientation of W so its rational Euler number is negative.
The determinant formula fixes it as -D/(skJ). Its normalized Seifert
data are then uniquely determined: for each cone alpha_i,

    omega_i*(skJ/alpha_i)=-D mod alpha_i,
    0<omega_i<alpha_i.

Pairwise coprimality is essential for this uniqueness. The solution is

    e0=-2; (alpha_i,omega_i)=(s,1),(k,k-1),(J,J-1),         (6)
    -2+1/s+(k-1)/k+(J-1)/J=-D/(skJ)<0.

Thus, up to orientation (irrelevant for the group), W has a connected
negative-definite star plumbing with central weight-2, a single weight-s
leg, and two all-negative-two legs of lengths k-1,J-1. This conclusion
uses the actual Seifert structure, its cone orders, Euler number and
normalized invariants; it is not a homeomorphism guessed from H1.

Nemethi Example1.12, equation(1.13), gives a direct non-rationality test:
for some integer l>=0 require sum floor(-l*omega_i/alpha_i)<=l*e0-2.
Take l=s+1. Since s>=21 and J=(3s-1)/2, one has

    1<l/s<2, 0<l/k<1, 0<l/J<1,
    sum floor(-l*omega_i/alpha_i)=-2-l-l=-2l-2.            (7)

The criterion holds with EQUALITY, for EVERY surviving D, without a
growing denominator search. Hence the graph is non-rational, and
Nemethi Theorem1.6 makes pi1(W) LEFT-ORDERABLE. It is nontrivial:
its abelianization has order D>=11. By (4), the actual target group
has this nontrivial left-orderable quotient, contradicting (1).

Therefore ALL remaining trefoil S3/S3 COMMON sources are excluded.
PF68/71/72 had closed every other common representation sector in
this actual-binomial two-stage quartic class. This closes the entire
trefoil COMMON-source case at all windings, including all source
configurations for the target keys paired with PF73's odd-j controls.

## Evidence, imported results and remaining frontier

[The checker](verify_rational_end_orderability.py) pins complete primary
page audits, PF73 and its predecessors; checks the unbounded cycle and
Seifert identities, exact R6 attachment counts, rational source cluster
controls, both trefoil periodic classes, normalized stars, the universal
floor witness and the retained non-rationality certificate. The written
arguments establish universality; numerical controls do not.
[The receipt](rational_end_orderability_receipt.json) and
[manifest](sources/rational-end-orderability/manifest.json) identify
the imported Neumann R6, Artin/Nemethi and Boyer--Rolfsen--Wiest theorems.
Their deep topology is imported, not independently reproved here.

The all-degree source theorem retains its explicit smooth, entire-line,
global common-coordinate hypotheses. The target application is the
actual-binomial TWO-STAGE QUARTIC trefoil class only. Other first pairs
in S3/S3 common sectors now have the three-cone quotient test (4)-(5),
but are not declared excluded. Passing COMPLEMENTARY sectors, other
keys, longer sequences, higher degree and the full goal remain OPEN.
No finite Keller map, general higher-degree smoothness, independent
full reproof, formal verification, external review or novelty/priority
claim is asserted. The Neumann splitting and graph-rationality/LO
imports are conditional on their precise published hypotheses as checked.
