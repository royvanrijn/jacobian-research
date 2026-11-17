# Original contact primitive and polar-divisor residues

16 November 2026. Written all-degree differential calculation on the
ORIGINAL exceptional quotients. Both FULL Mathieu--Zhao cases remain
**OPEN**. External review remains deferred. No new degree search or
holonomic computation is used.

The preceding [dominant-root calculation](DOMINANT_ROOT_RESIDUE_AND_FIBRE_GRAPH_GLUING.md)
identifies the whole ambient Haar line with a nonzero quadric residue.
There is now an explicit primitive for that residue OFF the leading
polar curve, and an explicit primitive for the original volume OFF its
charged lift. Restriction to these opens therefore KILLS the unit class.
The relative support and gluing data must be retained to compute the
intrinsic extension. These primitives do not prove that extension zero.

The general fact that slice-pole complements kill the Haar class was
ALREADY proved in GHAARPOLERESIDUE20261115, Section 6 of the retained
[quadratic trace proof](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md).
The present additions are the right-invariant contact primitive for
ANY charged h, its explicit dominant-root and fixed-cylinder residues,
and the exact family-extension presentation (20)--(23). The old pole
obstruction and interval transgression are preserved, not rediscovered
as a new full Mathieu argument.

## 1. A global original contact form

Let X=PGL2/G, G=O or I. Work first upstairs on PGL2. For a local matrix
representative g of a projective matrix put

    lambda=(dg g^-1)_(11)-(dg g^-1)_(22).             (1)

Rescaling g by a local invertible scalar adds a scalar matrix to
dg g^-1, so its diagonal difference is unchanged. Thus lambda is a
GLOBAL REGULAR algebraic one-form. Right multiplication by a constant
matrix preserves dg g^-1. Therefore lambda descends to the ORIGINAL
quotient X=SL2/Gamma. No averaging or substituted volume is required.

In the dominant-root chart

    g=[[r v,r(v u-1)],[1,u]], det g=r,

direct matrix multiplication gives

    lambda=dlog r-2v du,
    d lambda=-2 dv wedge du,
    Omega=dlog r wedge dv wedge du
         =-(1/2) lambda wedge d lambda.               (2)

Omega has the same original normalization as the preceding note.
All formulas rescale with a different fixed Haar normalization.

Let E be the effective left-torus vector field, so E=r partial_r in
this chart. Then

    lambda(E)=1, i_E d lambda=0.                     (3)

In particular E never vanishes, even on the original finite quotient.
The kernel of contraction with E on cotangent spaces has dimension 2.

## 2. An explicit primitive from ANY nonzero charged original function

Choose ANY nonzero original regular torus eigenfunction h of charge
k!=0, so E(h)=k h. Let D_h be its zero divisor and set

    B_h=-(1/(2k)) lambda wedge dh/h
         on X minus D_h.                             (4)

This is an ORIGINAL regular two-form on that open and a rational
two-form on X. Define theta_h=(1/k) dh/h. By (3), theta_h(E)=1, so
theta_h-lambda is horizontal. A horizontal one-form wedges to zero
with the horizontal two-form d lambda in dimension 3. Hence

    d lambda wedge dh/h
       =k d lambda wedge lambda=-2k Omega,
    d B_h=Omega.                                     (5)

This is a global identity, with no phase purity hypothesis. It applies
in every degree and to any chosen nonzero charged component of a phase.
Also dh wedge B_h=0, so for every positive integer j,

    d(h^j B_h)=h^j Omega.                            (6)

The two-form h^j B_h is regular on ALL of X. Formula (6) is consistent
with the retained one-sided Haar vanishing; it supplies no new mixed
pure-to-safe implication.

At a smooth point of the reduced support of D_h write h=z^m a, with
a invertible. With residue defined by the coefficient of dz/z FIRST,

    Res_(D_h) B_h=(m/(2k)) lambda|_(D_h).             (7)

The restriction on the right is closed. Indeed E is tangent to the
reduced divisor (differentiate z^m a and use E(h)=kh), and never zero.
On the smooth two-dimensional divisor d lambda is horizontal with
respect to that tangent E, so it vanishes. Thus (7) is an explicit
closed ORIGINAL one-form on every smooth component, with its actual
multiplicity. Singular intersections require the relative support
complex; no assertion that normalization alone captures them is made.

Since the original nonzero H^3(X;Q) is one-dimensional and represented
by Omega, (5) proves

    H^3(X;Q) -> H^3(X minus D_h;Q) is ZERO.            (8)

Exactness over C suffices to prove the rational assertion by extension
of scalars. This is ordinary cohomological localization of the ambient
unit on D_h. It is NOT localization or splitting of its family extension.

## 3. The primitive of the dominant-root quadric residue

Use the dominant-root blowup and its zero section E^o=Q2. A charged
function has the exact form

    h=r^k A(u,v).

Here A is a section of a power of the normal line bundle, not in
general a global scalar function on Q2. Let C be its zero divisor.
On the two charts of the preceding calculation,

    u'=1/u, v'=u^2 v-u, r'=r/u^2,
    A'=u^(2k) A,
    alpha=v du, alpha'=-v' du'=alpha-dlog u.

On Q2 minus C define

    beta_h=alpha+(1/(2k)) dlog A.                    (9)

The two transition terms cancel EXACTLY, so beta_h is global, and

    d beta_h=omega=dv wedge du,
    Res_C beta_h=m/(2k)                              (10)

at a smooth component of multiplicity m. Equivalently, its pullback
to X minus D_h is

    beta_h=(theta_h-lambda)/2,
    B_h=beta_h wedge theta_h.                        (11)

The right side of the first equality is basic for the left torus:
its contraction with E is zero, and it is invariant. It also shows
that beta_h is right-G invariant, so the calculation descends on the
free quadric quotient used in the preceding note. Finite quotient
singularities are not silently treated as smooth boundary charts.

For the leading polar component h of a surviving mixed phase in the
dominant-root gauge, write k=-a<0. Then

    f=r^(-a) A + terms of strictly larger charge,
    beta_h=alpha-(1/(2a))dlog A,
    Res_C beta_h=-m/(2a).                            (12)

Thus the nonzero quadric class becomes EXACT off C, with explicit
logarithmic residues on C. This also proves C is nonempty: otherwise
(9) would be a global primitive for the nonzero class computed on the
compact sphere. The zero curve cannot be omitted as irrelevant data.

## 4. What survives in the phase family, and what does not follow

Let F_s be a regular phase fibre and U_s=X minus F_s. Formula (5)
continues to hold on

    U_s minus D_h.

Consequently the ambient unit in H^3(U_s) restricts to zero there.
An argument that discards D_h, computes a constant residue system on
its complement, and invokes a UNIT-RETAINING quotient theorem is
invalid: this restriction does not retain the unit.

When D_h is smooth, the divisor support/Gysin sequence explicitly
places the ambient unit in the image

    H^1(D_h minus F_s)(-1) -> H^3(U_s),               (13)

with representative (7), up to the standard comparison orientation
and 2*pi*i normalization. Formula dB_h=Omega proves the image statement.
For singular D_h use cohomology WITH SUPPORT and its residue complex,
retaining component multiplicities and intersections. The image
statement follows from (8) and localization, without assuming smoothness.

Crucially (13) is a map INTO the complement system. A subsystem whose
image contains the unit does not give a retraction from the whole
system: classes outside that image can still have monodromy into the
unit. Thus even a split source in (13) would not prove the intrinsic
extension split. The actual relative support maps, their kernels and
the ensuing boundary gluing have to be computed.

Similarly the Cauchy form has an additional term. On X minus D_h,F_s,

    d(B_h/(s-f))
       =Omega/(s-f)+df wedge B_h/(s-f)^2.             (14)

Neither the second term nor the residues (7) may be dropped. Multiplying
by an arbitrary original marker produces dQ wedge B_h as well. Therefore
this rational primitive proves neither ordinary marked resolvent
algebraicity nor marked logarithm vanishing. It must not be promoted
from an ambient de Rham identity to a full actual period identity.

## 5. One fixed smooth original cylinder works for EVERY phase

There is no need to assume the leading polar divisor smooth in order
to obtain a smooth localization model for the unit. Choose once and
for all the minimum invariant binary form in the BOTTOM row:

    h_O=F_8(z,w)=z^8+14z^4 w^4+w^8,
    h_I=F_12(z,w)=z w(z^10+11z^5 w^5-w^10).          (15)

These are the retained Klein forms, in the normalizations of
[the quotient model](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md).
They have distinct projective roots. Directly, F_8 has nonzero end
coefficients and its affine polynomial has no common root with its
derivative. For F_12 the roots 0,infinity are simple; the other roots
are given by the two distinct nonzero roots of t^2+11t-1, with t=z^5/w^5.
The root-orbit lower bounds 6 on O and 12 on I show the 8 or 12 roots
form ONE orbit: two orbits would require at least 12 or 24 roots.

Upstairs in SL2 the zero locus of h is a DISJOINT union of n=8 or 12
smooth surfaces, one for each bottom-row root line. The bottom row
cannot vanish. On one component choose an eigenbasis for its cyclic
stabilizer and write matrices as

    [[b a,-a^-1],[a,0]], a!=0, b in A1.              (16)

The stabilizer in Gamma has order ell=|Gamma|/n, namely 6 or 10.
It is diagonal in this basis and acts by a->zeta*a, b->b, with zeta
a primitive ell-th root. The entire ORIGINAL quotient divisor is
therefore the smooth connected cylinder

    D=D_h ~= Gm_c times A1_b, c=a^ell.                (17)

This is a fixed original divisor, independent of the phase, degree,
torus support and purity. Its effective charge is k=-n/2. On (16)
formula (1) restricts to lambda=-2dlog a, so (7) descends to

    Res_D B_h=(1/N)dlog c,
    N=|Gamma|/2=24 or 60.                            (18)

The scalar in (18) uses exactly the original Omega normalization in
(2); rescaling Omega rescales it too. This is not an arbitrary marked
form or a weighted pullback of Haar measure.

## 6. The intrinsic global extension is a fixed-divisor localization extension

For ANY regular phase f, choose a common regular algebraic parameter
base for U_s, D minus F_s, and their complements. Extra exceptional
values of f|D may be removed; this changes no global splitting test.
Write

    B=H^3(U_s), P=H^1(D minus F_s)(-1),
    I=image(P -> B),
    R=image(B -> H^3(U_s minus D)).                  (19)

The smooth-divisor localization sequence gives an EXACT sequence of
geometric variations

    0 -> I -> B -> R -> 0.                           (20)

By (5),(18) its kernel I CONTAINS the entire ambient line L=Q(-2).
Moreover I is pure of weight 4. Indeed P has weights 3,4, whereas
the retained full complement B has weights 4,5,6; strictness gives
I subset W4 B. In fact P itself is pure Tate of weight 4 here: deletion
of curves from the cylinder creates only weight-2 H^1, since H^1(D)
already has weight 2. The weaker weight bound suffices for (20).

Use the retained canonical decomposition

    W4 B=L direct_sum C_f,

where C_f is the entire globally nontrivial compact-middle summand.
Then I=L direct_sum C', where C'=I intersect C_f. The unique flat
map r_I:I->L fixing L kills C'. Push out (20) by r_I. The resulting
extension is represented by

    0 -> L -> B/C' -> R -> 0,
    epsilon_D in Ext^1(R,L).                         (21)

This is an EXACT global presentation of the original unit-splitting
problem. A global retraction B->L restricts to r_I, because C' has no
trivial global constituent. Conversely a splitting of (21) supplies
that retraction. Thus, before purity,

    epsilon_D=0 iff the original full-source unit
                        splits GLOBALLY.             (22)

For precision about the intrinsic boundary source, set

    K=C_f/C' subset R, V=R/K,
    E_boundary=B/C_f.

The canonical lift of K in B/C' is C_f/C'. Quotienting it in (21)
gives exactly 0->L->E_boundary->V->0. Hence

    epsilon_D=q^* epsilon_boundary, q:R->V.           (23)

Pullback in (23) is injective, since Hom(K,L)=0, by the elementary
Hom/Ext sequence. This also proves equivalence with the intrinsic
global extension, without mistaking a unit-containing subsystem for
a quotient. Under eventual purity the retained intrinsic theorem now
identifies (22) with safety for EVERY original marker.

At a SINGLE puncture one must still take the quotient by K, or the
canonical reduced local source, BEFORE testing its connecting class.
Global Hom(K,L)=0 does not give local Hom=0. Local invariants of R
need not lift all invariants of V, and a naive connecting test on
R^T alone can miss the intrinsic obstruction. No invalid local
reflection through an unreduced compact-middle quotient is asserted.

Thus (20)--(23) are a uniform concrete support presentation of the
ACTUAL global unit extension, using one fixed original cylinder and
the explicit residue (18). They do not prove epsilon_D=0 on every
pure phase, nor identify the dominant-root polar curve's source with
this fixed-divisor presentation. Those are further mathematical steps.

## 7. The next all-degree residue/extension calculation

There is now explicit input on both sides of the localization sequence:
the nonzero original Haar class, a global primitive off its leading
charged divisor, and its closed residues m lambda/(2k). Track these
THROUGH the actual finite-phase support and Picard/NS comparison maps,
including singular intersections, then make the actual abelian
identity-component correction. The desired conclusion remains zero
pairing of the resulting corrected cochain with EVERY canonical
surviving intrinsic graph character at phase zero.

The localization extension (20)--(23) is now available for EVERY phase.
Compute its actual gluing or its intrinsic reduced local class, using
the explicit fixed-cylinder residue (18). The desired dominant-root
pole-to-extension localization and the full vanishing are still
unproved. The preceding relative-forest statement is usable
only after actual support localization, on the fibre graph with all
good components contracted. An ambient tree supplies neither step.

All original capacity, automorphism, strict-leading and primitive
constraints remain inputs, as does the completed family-transfer
theorem. No new pure phase with a fixed original bad marker is
constructed. The ordinary/exponential generating-function distinction
is unchanged. Both FULL exceptional Mathieu--Zhao cases remain OPEN.
