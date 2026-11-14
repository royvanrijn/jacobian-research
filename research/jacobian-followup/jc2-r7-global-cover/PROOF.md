# A global Nori obstruction for the remaining two-branch-image cover

14 November 2026. PF2R7NORI1 gives a new global necessary condition,
without assuming high contact for the unknown second image. Suppose
the degree-six finite normalization has the S6 monodromy of PF2OMIT8,
with first image C in the complete degree(6,10) gap21 family. Suppose
the second image B has polynomial normalization degrees(3,5), only
ordinary affine nodes, and mutual finite singular-locus avoidance with C.
Let I be their projective infinity intersection and let

    T=sum I_y(C,B), summed over finite nontransverse intersections.

Then

                         I+2T >=34.                 (1)

In particular all finite-transverse pairs with I<=32 are impossible
as these branch images. Contact34 is the equality case;38 and40 are
below the positivity threshold and are not excluded. The earlier high-
contact quintic cases remain open in the immersed rows A and C.

PF2R7LINE1 separately gives the correct general-target-line Hurwitz
budgets. They have genus3--5 for two index2 branch images of degrees10,5,
and genus5--7 when the second index is3. The old single-prime genus2
adjoint pencil cannot be transferred to these rows.

This is an authored written application of Nori's published theorem,
with exact blowup and permutation checks. The imported topology theorem
is not reproved by the checker. No formalization or independent external
review, global cover construction or JC2 conclusion is claimed.

## 1. A usable primary theorem, and the exact kernel to which it applies

We import [Nori, Zariski's conjecture and related problems, Proposition3.27,
journal page331](https://www.numdam.org/article/ASENS_1983_4_16_2_305_0.pdf).
For a nodal divisor D on a smooth projective surface, transverse to
another curve E, his strict inequality C_i^2>2r(C_i) for each component
of D implies that the kernel of

    pi1(X minus(D union E)) -> pi1(X minus E)

is abelian. Here r counts nodes of that component. We apply it only
with D the single irreducible strict transform of B. E contains C and
the whole added boundary. We do not assert positivity for C or
abelianness of the full complement group.

For the actual residual finite normalization, PF2OMIT8 supplies a
surjective monodromy representation

    pi1(A2 minus(C union B)) -> S6.

Its generic B meridian is a single transposition in row A, or a single
3-cycle in row C. The optional unramified omitted H in row A remains
in the complete source ledger. Its identity inertia does not add a
branch curve to this finite normalization, so H does not change the
monodromy space used here.

Filling B kills its meridian, so that meridian lies in the kernel N
above. The image of N in S6 is normal and contains the B inertia.
Such a normal subgroup is nonabelian: a transposition has conjugate
transpositions(12),(23) which do not commute; a 3-cycle has conjugates
(123),(124) which do not commute. Consequently Nori's positivity
hypotheses would contradict the actual required monodromy. A curve
meeting those hypotheses cannot be the second ramified image here.
This argument needs no presentation of the unknown complement group.

## 2. The precise hypotheses on the image pair

The first image has a single infinity branch of characteristic(4;10,31)
and infinity-line contact10, uniformly by PF2R7TOP1. A polynomial
degree(3,5) normalization has projective degree5, one infinity branch
of characteristic(2;5), and infinity-line contact5. Its normalization
is rational. With only ordinary affine nodes, arithmetic genus6 and
infinity delta2 force exactly four such nodes.

PF2R7B5CONTACT1 and PF2R7B5RAT1 give the complete necessary contact
spectrum

    I in{20,22,24,26,28,30,32,34,38,40}.

We retain the explicit assumption that B avoids Sing(C) and C avoids
Sing(B) at finite points. It is forced by the pointwise ledger in row C
and in row A when H=C or H=B. For example, c1,b2 would require at least
2+6 units in row C, and at least2+4+1 in either shared-image row A.
The symmetric singular-C condition follows likewise from PF2OMIT11.
For row A with H distinct this avoidance is an additional hypothesis;
we do not silently infer it from the degree ledger.

Every finite C/B intersection is then between two smooth branches.
A nontransverse one has integer multiplicity m>=2 and contributes m
to T. The first image can have other singularities away from B; they
will lie in E and do not need to be resolved for this application.

## 3. Exact loss of B self-intersection at infinity

Resolve only far enough to make B smooth and transverse to the added
boundary, separated from C there. The two branch characteristics and
their pairwise intersection numbers, including those with the infinity
line, determine the embedded resolution type. This classical local
input is recalled in [Campillo--Greuel--Lossen, Section1](https://arxiv.org/pdf/math/0511407v2)
and was already pinned for PF2R7TOP1.

The first four common centers have B multiplicities

    2,2,1,1,

and C multiplicities4,4,2,2. Thus they account for20 intersection units
and cost10 in B self-intersection. A concrete representative uses

    B: X=s^2,Y=s^5,
    C: X=r^4,Y=r^10+r^(I-10)+r^31 for I>20,
       X=r^4,Y=2r^10+r^31 for I20.

The infinity line is Y0. The four blowup charts are successively

    (X,Y)->(X,Y/X)->(X,Y/X^2)
      ->(X^3/Y,Y/X^2)->(X^3/Y,Y^2/X^5-1).

After them B is(s,0), transverse to the exceptional coordinate X0.
For I20, C is already at a different point. Otherwise the remaining
contact is I-20. It is even and at most20, while the remaining C cusp
has first odd exponent21. Hence exactly n=(I-20)/2 additional shared
centers are required, with multiplicities2 on C and1 on B at every
center. Each removes two contact units and costs one B self-intersection.
After the last blowup C meets a different point and B is transverse
to exactly one boundary component. Further resolution of C does not
touch B. Therefore

              B_infinity^2=25-10-n=15-(I-20)/2.      (2)

The checker verifies valuations, separation, all multiplicities and
Noether's contact sum for every possible I. It also checks36 as a
local-model control, while not admitting36 into the actual polynomial
contact spectrum.

## 4. Finite tangencies, including the last boundary crossing

A transverse finite C/B crossing already meets Nori's hypothesis and
is left unchanged. For a smooth tangency of multiplicity m>=2, use
local branches B:Y0 and C:Y=X^m. Blow up at their common point along B.
After m-1 blowups the two strict branches are transverse, but B passes
through the node of E consisting of C and the last exceptional curve.
One final blowup is required so that B meets E transversally at a
smooth point of E. Thus there are m centers of B multiplicity1,
not m-1. The B self-intersection decreases by exactly m.

Perform this at every finite tangency. All B nodes remain untouched,
since C avoids them. On the resulting smooth projective X,

    B_tilde^2=15-(I-20)/2-T, r(B_tilde)=4.

Let E consist of the strict C, the total infinity boundary and these
finite exceptional curves. B_tilde is nodal and intersects E transversally,
away from singularities of E. Both complements in Nori's theorem agree
with the original complements for removing C and optionally B: every
blowup center lies in the removed C or infinity boundary, and all
exceptional curves are included in E.

Nori would give an abelian N if B_tilde^2>8. Section1 excludes that.
Therefore15-(I-20)/2-T<=8, which is exactly(1). This proves PF2R7NORI1.
It is a global covering obstruction with a named topological import,
not an inference from local length nonnegativity.

With T0 it excludes I20,22,24,26,28,30,32. At34 the square is8 and the
strict positivity hypothesis fails. At38 and40 it is6 and5. Neither
the equality case nor the two high-contact cases is closed by this
argument. Tangencies or finite singular-locus coincidences can also
escape the stated clean-pair version; their hypotheses must be checked.

## 5. The correct target-line genus budgets

Let a general target line pull back to the smooth connected curve of
genus gamma on the actual minimal resolution. Put s=sum f_i over all
primes above generic target infinity. It is the number of geometric
points over infinity, and1<=s<=6, since sum e_i*f_i=6.

More generally, if B has polynomial normalization degrees(3k,5k) in
the same coordinates, its image degree is5k. The affine ramification
degree on a general line is10+5k*(e_B-1), and Riemann--Hurwitz gives

    2gamma-2=-12+[10+5k*(e_B-1)]+[6-s],
    s=6+5k*(e_B-1)-2gamma.                          (3)

Unramified omitted primes do not contribute to this ramification sum.
No assertion that individual f_i are one is needed. For k1 this gives

| Other index | General-line genus | Sum of infinity residue degrees |
|---|---:|---:|
| 2 | 3,4,5 | 5,3,1 respectively |
| 3 | 5,6,7 | 6,4,2 respectively |

Also K.H=2gamma-8 on this degree6 model. JC2POS1 gives
Q^2=24gamma-4rho-2>0, so rho<=6gamma-1. Hodge index gives
rho>=ceil(10-(2gamma-8)^2/6), as well as rho>=2. The exact table is in
budgets.json. These are necessary numerical windows, not realized
surfaces or complete boundary attachments. The old genus2 adjoint
pencil used affine branch degree10 alone. That hypothesis is false
in these two-prime rows, so its rational-conic reduction cannot simply
be reused. This proves PF2R7LINE1.

## 6. Why the nodal-complement shortcut does not apply

The fresh literature check found [Cogolludo-Agustin--Elduque, published
8 April2026](https://link.springer.com/article/10.1007/s13398-026-01853-1),
with affine nodal pairs having the same local infinity data and different
complement groups. This motivates retaining the actual global embedding;
their examples do not compute the group of our specific pair.

Orevkov's [negativity theorem, Section1](https://www.math.univ-toulouse.fr/~orevkov/fg.pdf)
does require a separation condition on distinct projection branches
at infinity. Our C fails that condition for every affine linear
coordinate system. In p=T^6 coordinates, the first odd term of q(T)
is a*T^-11 with a nonzero. The branches T and -T have q difference
tending to zero. For any other linear projection, normalize it to
Z=q+b*p. Its two corresponding roots have

    T'=-T-(a/5)*T^-20+lower terms,
    p(T')-p(T)=(6a/5)*T^-15+lower terms.

This is again nonzero and tends to zero. Every independent transverse
linear coordinate differs by a nonzero multiple of this p difference
on equal-Z fibres. Thus the needed negativity condition is unavailable
in the retained affine coordinates. The paper's own2003 correction
also withdraws its stated Jacobian reduction; no such reduction is
used here. Nori's precise component-wise application above replaces
this shortcut, without claiming the whole group is abelian.

The remaining tasks are genuine global ones: the contact34 equality
case, contact38/40, special finite incidences, other B charts, actual
joint complement relations with S6 inertia and the normal source A2.
An exploratory key-polynomial identity is retained in discovery.json
for future pencil work; it is not a premise of these results, and a
field factorization on an image curve is not an intermediate field
of the degree-six source extension.
