# All-index peripheral shear, including cyclic singularities and residue joining

PF2RULED127, 16 November 2026. The exact ramified shear is determined for
EVERY critical prime of EVERY genuine finite intermediate of an actual
complex plane Keller normalization. No index, residue, degree, knot type
or iteration-length bound is imposed. Retain the ORIGINAL finite N, its
WHOLE marked source A2 and ALL omitted primes throughout.

For a critical prime E of an intermediate M, let e>=2 be its transverse
index, f its degree onto the image normalization, n its cyclic quotient
order (n=1 when M is smooth along E), and delta_C the sum of the finite
delta invariants of its image component C. With positive meridian mu and
preferred Seifert longitude lambda, normalize the actual torus lattice as

    (e,0), (a,f), 0<=a<e.

Then

    a = f*(2*delta_C-1) + n^(-1)  modulo e.             (1)

The inverse is taken modulo e; gcd(n,e)=1 by PF106/PF123. In particular:

* every individual peripheral orbit is CYCLIC as an abelian action;
* EVERY zero-shear orbit with gcd(e,f)>1 is excluded;
* if e is even, a+f is ODD, so every even-index/even-residue orbit with
  even shear is excluded, including every such zero-shear orbit.

These are complete all-degree marked peripheral exclusions, including
singular sources and imprimitive original monodromy. They do not exclude
all remaining monodromy or prove the full finite algebra. The shared HC4
full-sheet empty-fibre obstruction remains OPEN.

## 1. The normalized power map and its location descend to intermediates

PF123 already proves that every intermediate critical curve is a disjoint
smooth A1, with at most cyclic quotient-axis germs of type1/n(1,-e),
gcd(n,e)=1 and n<=e-1. We add the GLOBAL single-point/power-map location
needed here, rather than assuming the full original Borisov theorem
applies to a surface without a whole-plane open embedding.

Fix E and any original critical D mapping onto it. The induced maps of
normalizations are finite maps of A1s. Their composition is the original
normalized image map. The retained Borisov input (PF105) identifies that
composition, in affine coordinates, with z->z^F; it is unramified at every
regular point of N along D, and has at most the one cyclic quotient point.

First, if a point q on E is singular in M, every preimage p on D is
singular in N. Otherwise lift the finite smooth source germ at p through
the small cyclic uniformization of M at q, as in PF123. On the axis the
uniformization has degree n>1. The local normalized degree D->C^nu would
therefore be n times a positive integer, contradicting its degree one at
a regular original point. Since D has at most one singular point and
maps onto all of E, M has at most ONE singular point on E.

Second, suppose q is smooth in M. The inverse image of E in the local
factor at any p over q is the single original critical axis: another
component would be a second critical branch of the original finite map.
Off E the map N->M is etale, by the Jacobian chain rule. A connected
normal finite germ over a smooth surface branched only along a smooth
axis is the cyclic completion of a cover of its complement; pi1 of that
complement is Z. Its completion is (s,z)->(s^k,z). Hence p is smooth and
D->E is locally unramified. Since D->C^nu is unramified there, E->C^nu
is unramified at every SMOOTH point of M on E.

Finally write the two polynomial maps D->E->C^nu as H,F. From
F(H(z))=z^r_D for a positive integer r_D, F has only one root over0:
two different roots would both have nonempty disjoint inverse images
under H. Thus F is a power after affine changes, as is H. Its only
possible finite ramification point is the single quotient point of E
just identified. Consequently the inherited data are exactly

    E=A1_u, C^nu=A1_t,  t=u^f,
    n=1 ==> f=1,

and for n>1 the quotient point is u=0. This proves the required
hereditary power-map LOCATION, not an intermediate-plane recognition.

## 2. A cyclic cover of a neighborhood of the ENTIRE affine prime

For n>1, nE is Cartier on M: E meets no other singular point, and its
local Cartier index at its quotient point is exactly n. The restriction
of O(nE) to E=A1 is trivial. Affineness lifts a nowhere-zero section on E
to a section tau of O(nE) on M. Its zero set is disjoint from E. By the
Chinese remainder theorem choose a principal affine neighborhood V of E
avoiding that zero set and the other critical curves. On V, O(nE) is
trivial and its canonical section becomes a function s with div(s)=nE.

Normalize V in the extension obtained by adjoining an nth root of s.
Near a smooth point of E, s is a unit times an nth power of a local
parameter, so the normalized root is etale there. At the quotient point,
the marked local type has faithful normal character; the normalized
root is exactly its small smooth cyclic uniformization. It is connected
of degree n because the local Cartier index is n. Its inverse axis E'
is a connected degree-n cyclic cover of A1 branched only at0. Therefore

    E'=A1_y,      u=y^n,      t=y^(n*f).               (2)

The cover is smooth along ALL of E', not asserted smooth everywhere.
It is unramified in codimension one, NOT a normal ramification cover of
index n along E. Normalizing the root is essential. For n=1 take V itself.

This is the standard normalized-root construction. Ambro,
*Cyclic covers and toroidal embeddings*, arXiv1310.3951v2, Section3,
Proposition3.1 and Lemmas3.2,3.4--3.5, supplies a primary comparison for
finiteness, the normalized root of a unit times a parameter power, and
ramification index n/gcd(n,ord(s)). Section4 gives the irreducibility
criterion. The local uniformization and whole-axis argument above are
written explicitly. No flatness or etaleness at the quotient point of
the BASE is assumed. The cited primary HTML bytes are pinned.

Since E' is an affine line in a smooth surface neighborhood, its conormal
line bundle is trivial. Affineness lifts its generator to a function x
vanishing simply on E', and lifts y to a regular function. The pair(x,y)
is an etale coordinate frame along E'. Shrink away from E' if needed.
For the composite pi' to the target, the full Jacobian divisor near E'
is (e-1)E'. Thus, in these coordinates,

    Jac(pi')=x^(e-1)*U(x,y),       U(0,y)=c!=0.         (3)

The leading coefficient is constant because it is a unit on A1_y.
This global whole-axis constancy cannot be inferred merely from one
isolated analytic germ with an unspecified unit.

## 3. The correct target normal framing for arbitrary plane branches

Let h(t)=(P_C(t),Q_C(t)) be the polynomial normalization of C and put
g=gcd(P_C',Q_C'), M_C=deg(g), G=h'/g. Choose a polynomial w with
det(G,w)=1. Its holomorphic normal push-off nu_w of K at infinity has

    L=lk(K,nu_w)=2*delta_C+M_C,
    nu_w=mu^L*lambda.                                 (4)

PF126 proved this for the simple-axis image types. For an ARBITRARY
irreducible local plane branch, take a generic first coordinate of
multiplicity m and a reduced local equation b. A vertical displacement
has local self-intersection ord(b_Q(h(t))). The retained characteristic-
zero polar/conductor formula gives

    ord(b_Q(h(t)))=2*delta_branch+m-1.

Indeed the coefficient of epsilon in b(P_C(t),Q_C(t)+epsilon) is b_Q(h(t));
after removing the diagonal epsilon factor, conservation of local
intersection gives that order. Different branches contribute twice their
pairwise intersection multiplicities. Since ord(g)=m-1 at each chosen
normalization parameter, their sum proves(4), including all Puiseux pairs,
cusps, multiple values and arbitrary contacts. The fixed-ball disk-chain
argument of PF126 avoids a global deformation assumption at infinity.

This reuses the polar/conductor import of PF115, not a new theorem:
Garcia Barroso--Ploski, arXiv1812.06512v1, Proposition1.1(iii),(iv) and
its proof (the retained author/book PDF numbers the corresponding
proposition2.1). Its characteristic-zero specialization is essential.
The primary HTML was re-read and pinned; no fresh PDF audit is claimed.

## 4. Leading normal coefficient and exact monodromy

Set r=n*f. At a generic point of E', the target tubular coordinates are

    (t,zeta) -> h(t)+zeta*w(t),

whose Jacobian along zeta=0 is g(t). After adjusting the source tangent
coordinate to hold t fixed, generic ramification index e gives

    zeta=A(y)*x^e+higher x-order,
    t=y^r along E'.

Taking Jacobians and using(3) gives, up to a nonzero constant,

    A(y)=1/[g(y^r)*r*y^(r-1)].                         (5)

This is an identity in the rational function field of the axis; poles at
the exceptional normalization parameter are permitted. It is derived on
the generic smooth branch and used on a LARGE circle where all factors
are units. It does not assume rank one at the finite quotient point or
at any nonimmersive image point.

For a small target normal displacement zeta=epsilon the e roots satisfy

    x^e=constant*epsilon*g(y^r)*y^(r-1) to leading order.

One full y-circle therefore shifts the chosen positive meridian e-cycle
by r*M_C+r-1. Analytic root continuation for sufficiently small epsilon
gives the same permutation as the leading equation. That y-circle makes
r=n*f circuits of the target normalization; it is nu_w^(n*f).

Let lambda^f act on the selected meridian cycle as mu^k. The peripheral
elements commute, and lambda joins exactly f such cycles. Equation(4)
and the preceding root continuation give

    n*k+n*f*L = n*f*M_C+n*f-1 modulo e.

Using L=2delta_C+M_C and invertibility of n modulo e yields

    k=f*(1-2delta_C)-n^(-1) modulo e.

The normalized stabilizer vector(a,f) means mu^a lambda^f fixes a sheet,
so a=-k modulo e. This proves(1), with its sign checked in the genuine
finite examples below. All longitudinal joining is retained.

## 5. Complete consequences for marked covers

The joint peripheral action on a torus orbit of size e*f is the regular
abelian group Z²/<(e,0),(a,f)>. Its invariant factors are

    gcd(e,a,f),       e*f/gcd(e,a,f).

If a prime divided e,a,f, equation(1) would make n^(-1)=0 modulo that
prime, a contradiction. Therefore gcd(e,a,f)=1 and the orbit action is
CYCLIC. This excludes every NONCYCLIC JOINT PERIPHERAL ORBIT of every
genuine intermediate, not only an old knot family or a bounded degree.
For e=1 the orbit is automatically cyclic; every ramified orbit is a
critical prime covered by the theorem.

The ENTIRE peripheral image need not be cyclic: its actions on different
orbits can give a product of cyclic groups. Nothing here makes the global
monodromy cyclic or the generic rational deck group equal to monodromy.

For even e, n and its inverse are odd and 2delta_C-1 is odd. Thus

    a=f+1 modulo2.                                    (6)

Zero shear requires gcd(e,f)=1; even-index/even-residue zero-shear
joining is one excluded subclass. Nonzero shear cannot be dropped in
this test: an odd shear can pass where zero shear fails. Existing e<=6
residue-one results remain credited; these new exclusions apply at all
higher indices as well, without a local Morse condition.

For fixed(e,f,a,delta_C), the quotient order is RECOVERED by

    n = (a-f*(2delta_C-1))^(-1) modulo e,

as its unique representative1,...,e-1. A nonunit is forbidden. If this
recovers n=1 but f>1, the hereditary power-map location excludes it.
For several critical primes over the SAME C, the single odd integer
T=2delta_C-1 must satisfy all congruences

    f_i*T=a_i-n_i^(-1) modulo e_i,    T=1 modulo2.       (7)

Linear-congruence reduction and generalized CRT test this exactly, with
no parameter sweep. A solution is only peripheral compatibility; the
actual curve still must have the required delta and finite cover.

There is also an EXACT singular-order law along a genuine critical
factorization. If D_i above E_j has residue ratio l=f_i/f_j, then

    n_j=l*n_i,   n_i*f_i=n_j*f_j,
    gcd(l,e_j)=1.                                     (8)

To prove this, uniformize both source pairs at the points over the
quotient point and lift the finite factor map between the two smooth
uniformizations. The inverse image of the target axis is just the source
axis; the Jacobian chain rule leaves ramification only on that axis.
The smooth-axis complement has pi1=Z, so this lifted normal completion
is cyclic and maps the axes ISOMORPHICALLY. Its equivariance induces an
injection of the source deck group into the target deck group (both act
faithfully on the axes). The coarse axis map consequently has local degree
n_j/n_i. The global power-map location makes this degree exactly f_i/f_j,
proving the equality. The smooth cases have n=1 and give the same result.
The last assertion follows from gcd(n_j,e_j)=1.

Thus the uniformized normalized degree n*f is preserved while the curve
remains CRITICAL through a finite tower. Every original quotient order
above E_j is at most n_j<=e_j-1, even if its original index is much larger.
This is not asserted at the final target image C, which can be singular
and is not such a smooth critical axis. As a consistency check, lattice
containment gives a_i=l*a_j modulo e_j; formula(1) then gives the weaker
congruence n_j=l*n_i modulo e_j, exactly as the equality requires.

## 6. Genuine finite covers check the singular correction and joining

The retained PF123 smooth covers

    (x,y)->(x^k-3y², x^k*y-2y³), k>=1,

have e=2k,n=f=1, ordinary cusp delta_C=1. Equation(1) gives a=2 modulo2k,
which is NONZERO for k>=2. Weighted rotation gives the same answer:
the trefoil central fibre mu^6 lambda acts on the critical cycle as mu^4.
The open is C* times C, not a Keller whole plane. Thus zero shear is not
automatic beyond simple ramification, even on a smooth finite surface.

PF116's ALL-m dessin covers have e=5m-2,n=3m-1,f=1,delta_C=1. The identity
5n-3e=1 gives a=6 modulo e. Their actual trefoil regular fibre lifts
with order1, so mu^6 lambda=1, independently verifying the same shear.
The smooth-source formula a=2 would be wrong for these singular covers.

PF105's ALL-s quotient covers have e=2s+1,n=2,f=s and image branch
(2,4s-1), hence delta_C=2s-1. Equation(1) gives a=3 modulo e. Their
weighted action makes the torus-knot central fibre mu^(8s-2) lambda
have order dividing s, so lambda^s=mu^(-s*(8s-2)); its normalized shear
is again3 modulo2s+1. On the smooth upstairs cover n=1,f=2s the shear
is6 modulo e. This tests genuine nontrivial RESIDUE JOINING, not formal
passports. These covers have the wrong marked open and are not Keller
realizations or exclusions of the permitted nonzero-shear data.

The same PF105 polynomials, with source weights(4s-1,1), are invariant
under mu_(2s): their weighted degrees are 2s(4s-1) and4s. Hence N_s
also factors through M_s=C2/mu_(2s), of quotient order2s and residue1,
with unchanged index2s+1. This genuine degree-s further quotient checks
the EXACT law n_M=s*n_N and the invariant n*f=2s. Its regular-fibre
lift closes after one circuit, giving a_M=8s-2 modulo2s+1, while
a_N=s*a_M=3 modulo2s+1, as required by the lattice inclusion.

## 7. Relation to linking recursions and remaining goal

PF101's actual filling basis is unchanged. The vector of normalized
shears in F*e=a is now fixed by(1) at every ramified orbit; unramified
orbits have e=1,a=0. The equation still needs rationally null cores,
and integral torsion defects remain unless PF111's original complete
coverage is established. The linking ENTRIES are not reconstructed
from these congruences. No dihedral representation is substituted for
a different actual permutation cover, and the sealed authors-code
comparison is not rerun.

PF126 is recovered when e=2: n=1,f=1 gives a=0. Its binary affine
residual-cover and common-form results remain intact. The new scope
is all indices, all normal cyclic sources, residue joining and the
hereditary global location lemma. The proof is written without a bounded
degree or coefficient campaign. Cheap exact checks supply transcription,
source seals, lattice and genuine weighted-family controls, not the
unbounded proof itself.

Natural A/S and remaining binary affine compatibility, the shared HC4
EMPTY ORIGINAL full-sheet even-radial-zero profile and the full finite
algebra remain OPEN. Shared trace/marked-energy priorities and the
published quartic PLANE retirement remain. No intermediate is assumed
to be a plane. No historical priority, independent replay, formal
verification, fresh PDF visual audit or external review is claimed.
The full research goal remains ACTIVE.
