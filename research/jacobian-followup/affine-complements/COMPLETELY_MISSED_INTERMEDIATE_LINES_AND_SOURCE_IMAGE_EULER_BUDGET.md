# Completely missed intermediate curves are affine lines with exact Euler cost

PF2RULED140. 16 November 2026. ALL actual complex plane Keller maps,
ALL embedded intermediate fields and ALL original boundary primes.
This supplies the omitted-curve geometry left unspecified in HC4MYKV0,
then couples the actual replica source image of MR0/PF139 to Wright's
ENTIRE-intermediate Euler budget. Full Wright/HC4/JC2 remains OPEN.

## 1. Setup and statement

Let B=C[P,Q] subset R=C[x,y] be an actual plane Keller inclusion,
F=Frac B, L=Frac R. For F subset E subset L retain the FULL finite
normalizations and maps

    N --h--> M --p--> Y=A2, deg h=m, deg p=b, d=m*b.

Write V for the ENTIRE etale locus of p, I for the open image of the
WHOLE original X=A2 in V, and D for the reduced union of the curve
components of V minus I. Put r=#components of D. Then:

1. D is a disjoint union of r SMOOTH EMBEDDED A1s, closed even in M.
   Every component maps finitely and immersively to Y with residue
   degree1 onto the normalization of its image.
2. Above each D_j the entire inverse image is a disjoint union of
   original omitted A1s B_(j,a). N is SMOOTH along them, each maps
   ISOMORPHICALLY to D_j, and h is locally (u,v)->(u^e,v).
   Thus sum_a e_(j,a)=m and there is NO finite local excess there.
3. T=V minus D is AFFINE and canonically

       T=Spec(R intersection E), chi(T)=chi(V)-r.           (1)

   The original X->T is etale and almost surjective. Its possible
   FINITE missed points are not eliminated in general.

With PF129 notation A_M for total intermediate moved critical support
and Eta_M for finite excess, the exact source-image ledger is

    b=chi(T)+r+A_M+Eta_M.                                 (2)

For b>1, A_M>=2 and Eta_M>=0, hence

    b>=chi(T)+r+2.                                       (3)

This is an actual arbitrary-intermediate statement, not restricted to
deck quotients, relative Galois covers or equivariant factors.

## 2. Entire boundary images cannot cross a lower critical axis

We reuse the retained original boundary theorem PF2BATT1 and the
uniformized-axis geometry in PF123: the original omitted primes are
pairwise disjoint, have A1 normalizations, and their reduced inverse
images under small cyclic quotient uniformizations are smooth axes.
PF123 also makes the lower critical curves quotient axes, hence locally
Q-Cartier. The original X maps into V, as already used in KU/KV.

Let B_i be an original omitted prime whose generic image under h is
not a lower critical curve. It cannot meet h^-1(Crit p). Indeed, a
local positive power of a lower critical-axis equation is Cartier.
Its pullback has pure divisorial zero set contained in N minus X.
At a hypothetical point on B_i, disjointness leaves only B_i as a
possible component of that zero set. This would force B_i generically
into the lower critical curve, a contradiction. This is PF130's
Q-Cartier/disjoint-axis argument, now used also for unramified B_i.

Consequently h(B_i) lies ENTIRELY in V. Its image is closed in M
because h and the normalization of B_i are finite. If its generic
image is one of the completely missed curves D_j, its whole image
equals that D_j. In particular these are WHOLE affine boundary images,
not punctured images obtained by intersecting V with a larger curve.
The closure of a missed curve cannot acquire points of I: M minus I
is closed. Therefore h^-1(D) is supported wholly on original omitted
primes, at every point and not just generically.

## 3. Smooth inverse-image descent and the local power map

Fix q in D and z in N over q. Take a small cyclic uniformization
rho:C2->(N,z); if N is smooth take the identity. The germ

    g=h*rho:(C2,0)->(V,q)=(C2,0)

is FINITE. The reduced inverse image of the ENTIRE local union D is
contained in the one original boundary axis through z, and is nonempty;
it is therefore exactly that smooth axis. The previously imported
Giraldo--Roeder Theorem A, used in PF123, makes D smooth at q.
This treats the UNION, so distinct components cannot meet and a
single component cannot have a cusp or self-intersection in V.

The exact reduced/set-theoretic convention was rechecked in the primary
[1911.00628v1, Theorem A and Section2.2](https://arxiv.org/html/1911.00628v1).
We reuse its retained bytes and import boundary; no independent reproof
of that theorem or new hypersurface theorem is asserted.

There is a stronger elementary local conclusion. Choose target
coordinates (s,t) with D={s=0}, and source coordinates (u,v) with
g^-1(D)_red={u=0}. Absorb a nonvanishing analytic kth root to write

    s(g(u,v))=u^k, t(g(u,v))=G(u,v), k>=1.                (4)

The only possible critical divisor of g is u=0: off the original
boundary the original Keller map and the uniformization are etale,
and p is etale here. Thus G_v is a unit times a nonnegative power
of u. A positive power would make G(0,v) constant, contradicting
finiteness of g. Hence G_v is a UNIT, and (u,G) are source coordinates:

    g=(u^k,v).                                           (5)

Every element of the small uniformizing group fixes v and u^k, so
acts by (u,v)->(zeta*u,v). A nontrivial such element fixes a whole
line, contradicting smallness/free action off0. The group is trivial.
Thus N itself is smooth at z, and h has local form (5). This proves
the local-power and zero-excess assertions without assuming them from
generic residue information or ignoring singular original points.

Each D_j is the finite image of an original A1 normalization. It is
smooth, rational, affine and has only one point at infinity, so is A1
by the same finite-image argument already used in PF123. Each B_(j,a)
above it is now smooth, has normalization A1, and maps unramifiedly
at EVERY point by (5). A finite unramified polynomial A1->A1 has degree1.
This gives the asserted isomorphisms and sum e_(j,a)=m.

Finally p is etale along D_j, so its restriction is immersive and finite
onto its curve image. The normalization of that image is A1, and the
finite polynomial map from D_j to it is unramified, hence also degree1.
Global image curves of different D_j may COINCIDE, and an image may
self-identify distinct parameters. Embeddedness is in M, not necessarily
in the FINAL target Y. No identification of those image points is erased.

## 4. The actual maximal image and the Euler ledger

On smooth affine V, D is effective Cartier; its complement T is affine
(locally principal on an affine cover, so the open immersion is affine).
The original image I is contained in T with only finitely many missed
points. Hartogs on normal T, or the faithfully flat local descent along
X->I followed by Hartogs, gives

    Gamma(T,O)=Gamma(I,O)=R intersection E.

This is the canonical maximal-ring identification of KV0, now with the
NEW disjoint-line geometry. It is NOT a surjectivity claim for X->T:
finite point gaps inside T remain possible without MR0/MQ0 hypotheses.

Euler additivity gives chi(V)=chi(T)+r. PF129's exact ENTIRE-intermediate
identity chi(V)=b-A_M-Eta_M gives (2). If b>1 the finite normal cover
M->Y must branch; otherwise purity and affine-plane finite-cover
triviality force b=1. Every nontrivial inertia moves at least two sheets,
so A_M>=2. This proves (3). No finite-map hypothesis is imposed on T->Y
or X->T in using this calculation; the finite map is the original M->Y.

For completeness, let U be the entire ORIGINAL etale completion in N,
and let u_j be the number of e_(j,a)=1 components above D_j. Exactly
those u_j whole axes lie in U, so

    w_j=sum_(e_(j,a)>1)e_(j,a)=m-u_j,
    chi(U intersection h^-1(T))=chi(U)-sum_j u_j.          (6)

These w_j are PACKET moved supports, not necessarily the whole support
over their final target curve: other intermediate sheets may lie over
the same image. Removing them leaves the relative PF130 Euler defect

    m*chi(T)-chi(U intersection h^-1(T))
      =Delta+Eta_res-sum_j w_j.                          (7)

The removed packets have zero excess by Section3. The original FULL
normalization and its discriminant are not replaced by these opens.

## 5. A uniform quotient-degree cost for the replica mixed alternative

Now assume MR0's replica equivariance, and write its source ruling as
u=h(y), deg h=m, delta=p_0(y)partial_x and p_0=h'*rho(h). MR0 removes
the finite source-image gaps, so X SURJECTS etale onto the same T.
For a fibre let s_a be its number of components and let nu_a be the
number of DISTINCT root multiplicities occurring in h-a. Its exact
source-line multiplicity rule implies s_a>=nu_a. Therefore

    chi(T)=1+sum_a(s_a-1)>=1+sum_a(nu_a-1),
    b>=r+3+sum_a(nu_a-1).                                (8)

The sums are finite; ordinary fibres have nu_a=s_a=1. This bound uses
the actual target quotient degree b, not the relative polynomial degree m.

In PF139's mixed alternative, the wholly critical fibre has normalized
parts>=2 and gcd1. They cannot all be equal, so nu_(a_*)>=2.
It has q>=2 distinct roots; its contribution m-q to h' leaves q-1>=1
critical zeros elsewhere. Hence there is ANOTHER critical value a!=a_*.
PF139's derivative budget makes that other fibre have a simple root;
being critical it also has a multiple root, so nu_a>=2. Thus

    chi(T)>=3, b>=r+5.                                   (9)

Together with PF139's m>=5*g, a minimal-degree mixed equivariant factor
requires d=m*b>=5*g*(r+5). Here g=|pi1(T)| is NOT silently a generic
deck order. Equation(9) itself needs no minimality if the stated mixed
polynomial fibre exists. For Morse h of degree m>=3, each of its m-1
distinct critical values has multiplicities2 and1, giving the all-degree
comparison b>=m+r+2 before any minimal-degree exclusion.

These are necessary conditions for ACTUAL equivariant factors. They do
not construct or eliminate every mixed profile, do not reopen noncyclic
generic decks closed by MP0, and are not a growing partition campaign.

## 6. Source and completion boundaries

The complete KV0, PF123, PF129, PF130 and MR0 proofs are reused at the
stated scopes. The original boundary and smooth-inverse-image theorems
remain explicit adopted inputs, not independently re-proved. The new
local power-map argument is written in Section3. No new external
classification/theorem import, predecessor mathematical replay, novelty,
formal verification, independent whole proof or external review is claimed.

The shared MS0 first-source-derivative theorem is also credited: its
maximal additive open concerns one directional saturation in the full N.
Here the deleted curves concern EVERY intermediate source image, without
an additive-action hypothesis. Neither result makes MS0's first-derivative
integrality tests or full finite trace-lattice stability automatic.

General finite source-image POINT gaps, the mixed and cyclic actual
target compatibility, trivial-deck and non-equivariant cases, the full
finite trace-dual multiplication and marked-energy implication remain
OPEN. The full user goal is ACTIVE. The original ramified primes and
all finite normalization data are retained throughout.
