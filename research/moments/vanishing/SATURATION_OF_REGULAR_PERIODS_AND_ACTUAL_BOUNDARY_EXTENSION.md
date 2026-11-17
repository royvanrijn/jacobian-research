# Saturation of regular periods and the actual boundary extension

16 November 2026. Internal written proof using the retained Hodge and
actual-period theorems. Independent replay, formal verification and
external review are not claimed. Both FULL exceptional MZ cases remain
**OPEN**. Review is deferred and is not the next task.

For every original phase, ALL smooth models projective over the phase
line give the SAME regular-period image after adjoining the literal unit:
the canonical weight-four sector of the FULL original marked module.
No holomorphic extension of the original volume is needed for this
unit-adjoined equality. Without adjoining the unit, its presence is
determined exactly by whether the ambient cohomology class lifts.
An original-volume proper model is the unit-containing special case.
Thus the missing marked sector is intrinsically model-independent for
EVERY phase, not just phases admitting a holomorphic-volume model.

The resulting extension test uses the actual global period kernel, not
an identification of a scalar Cauchy function with the full source. The
known exceptional uniform marker gives an explicit constant line in the
missing sector; pulling back the extension to that line is precisely the
Kummer extension of log((s+1)/s), for EVERY regular Bezout parameter b.
The theorem does not prove purity or supply a pure exceptional failure.

See the [checker](verify_haar_saturation.py),
[receipt](haar_saturation_verification.json), and
[source audit](haar_saturation_source_audit.json).

## 1. Canonical full-source and actual-function sectors

Let X=SL2/Gamma for ANY finite Gamma, and let f be ANY nonconstant
regular phase. On a common full regular parameter base S put

    B_f=H3(X minus X_s),
    0 -> L=Q(-2) -> B_f -> H2(X_s)(-1) -> 0.          (1)

Superscripts are suppressed in the displayed H3,H2: these denote
COHOMOLOGY. The retained
[canonical-boundary theorem](CANONICAL_BOUNDARY_REDUCTION_AND_LOWEST_WEIGHT_INVARIANTS.md)
proves, on the ORIGINAL global base,

    W4 B_f=L direct_sum C_f,
    C_f ~= W2 H2(X_s)(-1),
    C_f has NO global trivial subquotient.             (2)

The last assertion uses global pure semisimplicity and the nonproper
Leray/lowest-weight argument, not local monodromy semisimplicity.
Deleting extra regular values changes no global monodromy image.

Let q:B_f -> M_f be the FULL actual marked-period function quotient,
with literal unit u=q(L) retained BEFORE quotienting by C[s]. Rational
connections are over C(s); the actual kernel need not be Hodge. Define

    I_f=q(C_f), R_f=q(W4 B_f)=C(s)u direct_sum I_f,
    D_f=M_f/R_f, E_f^act=M_f/I_f.                      (3)

All these are finite-rank global differential modules. In particular
R_f is the image of a geometric weight piece, NOT an assumed Hodge
filtration on M_f. Since C_f is semisimple with no trivial constituent,
its actual image I_f has the same properties and does not contain u.
Consequently

    Hom_D(I_f,C(s))=0,
    Hom_D(R_f,C(s))=C*r_f, r_f(u)=1, r_f(I_f)=0.        (4)

Thus R_f has a UNIQUE global unit retraction. This statement applies to
EVERY original phase and ALL original markers, without a degree, weight
support, genus, graph or proper-model hypothesis.

## 2. Universal saturation for EVERY proper phase model

Let f extend through an OPEN inclusion X subset Y to a smooth Y
PROJECTIVE over A1_f. Such models exist by graph compactification and
resolution for EVERY phase. We do NOT assume Omega extends. Put

    B_Y=H3(Y minus Y_s),
    K=ker(H3(Y)->H3(X)),
    J_Y=image(H3(Y)->H3(X)) subset L.                  (5)

The rank of J_Y is either0 or1. The first Gysin image in localization
is in K, by the same argument as the retained
[proper-unit proof](PROPER_VOLUME_UNIT_PERIODS_AND_COMPLEX_HOPF_CLOSURE.md):
for a sufficiently large physical s, the original compact middle cycle
avoids X_s, so L injects into B_f. Naturality kills the restricted Gysin
image; flatness makes this true on the common regular base. Consequently

    E_Y=B_Y/image(K),
    0 -> J_Y -> E_Y -> P_Y -> 0,
    P_Y=ker[H2(Y_s)(-1)->H4(Y)].                       (6)

J_Y is either zero or pure weight4, and P_Y is pure polarizable weight4.
Thus E_Y is PURE weight4 and globally semisimple in BOTH cases. This
uses the dimension3 ambient weights and H3(X)=Q(-2); it does not use
holomorphic extension of the original form. Restriction factors through
E_Y and gives

    image(B_Y -> B_f) subset W4 B_f.                   (7)

We next prove that the image contains ALL C_f. Projective smooth fibres
Y_s compactify X_s. The lowest-weight image theorem therefore gives

    H2(Y_s)(-1) -> W2 H2(X_s)(-1) ~= C_f ONTO.         (8)

Pure global semisimplicity lifts C_f to an isomorphic subvariation C'_f
of the left side. Since C'_f has NO trivial constituent, its Gysin map
to the CONSTANT H4(Y) is zero. Localization lifts it through B_Y.
Naturality then shows that (7) surjects onto C_f. Since (7) is a pure
semisimple subobject of L direct_sum C_f, its nontrivial isotypic part
is exactly the canonical C_f. Only its unit line remains to be decided.

If J_Y=L, (7) contains that ambient line and hence is ALL W4 B_f.
If J_Y=0, E_Y=P_Y. We show that every invariant vector of P_Y maps to
zero in B_f, ruling out a unit line created by fibre residues.

Here is the required projective-fibre invariant-cycle lift, with its
proof. On S the classical Leray sequence for Y_S->S has only columns0
and1 because S is a graph up to homotopy. Its edge map
H2(Y_S)->H0(S,R2 f_*Q) is onto. For a fixed smooth projective fibre,
restriction is a strict mixed-Hodge morphism to pure weight2. Its image
therefore equals the image of W2 H2(Y_S). The lowest-weight image theorem
identifies W2 H2(Y_S) with the image of H2(bar_Y), for a smooth proper
compactification bar_Y containing Y. Thus EVERY global invariant fibre
class lifts to H2(Y) via restriction from bar_Y. This is the curve-base
fixed-part argument, using the already retained compactification image
and MHS strictness. No nonproper invariant-cycle theorem with missing
hypotheses or new general fixed-part theorem is imported.

Apply this to beta in the invariant part of P_Y. Choose alpha in H2(Y)
restricting to beta (with its Tate twist understood). The product

    alpha cup eta_s in H3(Y minus Y_s),
    eta_s=(f-s)^*[dlog z/(2*pi*i)]                     (9)

has residue beta. But its restriction to X minus X_s is zero, since
H2(X)=0 makes alpha|X=0. When J_Y=0, residue identifies E_Y with P_Y,
so (9) proves that the image of EVERY trivial constituent of P_Y is
zero. Its globally semisimple image in W4 B_f thus contains no L.
We have proved the exact ALL-PHASE dichotomy

    image(H3(Y minus Y_s)->B_f)
       = C_f                 if J_Y=0,
       = L direct_sum C_f    if J_Y=L.                 (10)

In particular adding the SAME original ambient line always saturates
EXACTLY W4 B_f. Disconnected fibres and their every component are
included throughout. No degree, support or genus restriction is used.

## 3. Universal model-independence after the ACTUAL period kernel

Define Reg_Y(f) to be the image in the FULL original marked-function
module of the ENTIRE algebraic de Rham cohomology of Y minus Y_s,
under the actual original compact-cycle orbit period map. This means
full de Rham hypercohomology, not only global regular 3-forms when the
complement is nonaffine. Restriction and integration over that same
cycle orbit commute. Applying q to (10) gives

    Reg_Y(f) = I_f               if J_Y=0,
    Reg_Y(f) = C(s)u + I_f=R_f   if J_Y=L,
    Reg_Y(f) + C(s)u = R_f       ALWAYS.               (11)

The direct sum in the second case is canonical by (2)--(4). These are
identities AFTER exact-form relations and ALL actual GLOBAL period
relations. The actual kernel need not be Hodge. The unit-adjoined image
is the SAME for EVERY proper phase model Y of EVERY original phase.
Thus choosing another model cannot recover any missing marked direction
outside R_f. The missing unit-adjoined sector is EXACTLY D_f=M_f/R_f.

If the SAME original normalized volume does extend holomorphically to
Y, its nonzero period proves J_Y=L. This recovers and strengthens the
proper-volume result: its entire regular image is R_f, and the actual
unit Cauchy function belongs to it because Omega_Y/(s-f) is regular.
Without that volume hypothesis we do NOT assert that the unit Cauchy
function belongs to Reg_Y(f). Adjoining the literal1 is a cohomological
sector operation, not a claim that the original integrand extends.
In particular mixed pure SAFE phases with intrinsic volume poles are
fully covered by the unit-adjoined equality, without pretending that
their original top form becomes holomorphic.

For any phase and any model Y, v in M_f belongs to the unit-adjoined
regular image iff its class in M_f/I_f belongs to C(s)u. The unit Cauchy
submodule need not generate R_f; R_f need not equal the full M_f.

## 4. The exact boundary-lifting obstruction, for EVERY phase

Apply Hom_D(-,C(s)) to 0 -> R_f -> M_f -> D_f -> 0. The connecting map is

    delta_f: Hom_D(R_f,C(s)) -> Ext1_D(D_f,C(s)).

By (4), define the canonical class

    epsilon_f=delta_f(r_f).                            (12)

The pushout along r_f is exactly

    0 -> C(s)u -> E_f^act=M_f/I_f -> D_f -> 0.          (13)

This gives a useful necessary-and-sufficient statement, BEFORE PURITY:

    u splits GLOBALLY in FULL M_f iff epsilon_f=0.      (14)

Indeed a full unit retraction restricts to r_f, since (4) leaves no
choice. Conversely vanishing of the connecting class extends r_f to
M_f. This is the explicit lifting obstruction, not the assumption that
one can extend a functional from a submodule.

Under EVENTUAL PURITY, the retained intrinsic/global safety theorem
identifies (14) with safety of ALL original markers. Thus the FULL
exceptional goal is epsilon_f=0 on the entire exceptional eventual-pure
locus, or an actual invariant pure phase with nonzero epsilon_f and a
fixed original bad marker. The boundary source controlling (13) has
weights5 and6 after removal of C_f, but D_f and E_f^act themselves need
not inherit a MHS from their actual kernels.

This formulation is global. Although r_f is the unique GLOBAL
retraction of R_f, local retractions on a punctured disk can be more
numerous, because C_f may have local fixed vectors. We do NOT identify
the local Ext class in (13) with an arbitrary choice of a local full-
source retraction. For local corrected-component tests, keep the retained
canonical geometric source and its actual genus/torus corrections.

## 5. An explicit constant boundary line for EVERY Bezout parameter

On both exceptional quotients retain the actual phases and fixed markers
from the [uniform-marker theorem](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md):

    g_b=(1+Y_b)(1+F0+Y_b F0),
    C_b=W_b(1+Y_b),
    A_C(s)=h(C_b/(s-g_b))=log((s+1)/s).                (15)

Here b is ANY regular invariant parameter, not just a torus-equivariant
one. All functions in (15) are original regular invariant functions.
Let v=A_C. Then

    partial_s v=a(s)u, a(s)=1/(s+1)-1/s.               (16)

The rank2 module J=C(s)<u,v> is a nonsplit extension of the trivial
line by itself: no rational meromorphic primitive of a(s) exists,
as its residues are -1 at0 and +1 at-1. Its only rank1 global submodule
is C(s)u. Since R_g is globally semisimple, J cannot be a submodule of
R_g; hence

    J intersection R_g=C(s)u,
    e_b=[v] in D_g is NONZERO and horizontal.          (17)

Thus e_b spans an ACTUAL global constant line in the missing boundary
sector. Pulling (13) back along C(s)e_b -> D_g gives EXACTLY J:

    epsilon_g restricted to C(s)e_b
       =[ (1/(s+1)-1/s) ds ] in H1_dR(S).              (18)

The residue vector is (-1,+1), in the order (0,-1). This computes a
rank-one restriction of the canonical ACTUAL boundary extension class,
not only a raw marker pole on a chosen compactification. The equality
survives all exact-form and actual global function relations. Equivalently
I_g intersection J=0, so the Kummer module injects into E_g^act.

For the torus-equivariant parameters of weight n+m, the
[original-source model theorem](CONSTANT_VOLUME_MODELS_AND_ACTUAL_MARKER_DEFECTS.md)
supplies proper models, so (17) describes exactly their missing regular-
period direction. These equivariant phases and every nonconstant
polynomial recomposition are already known NOT eventually pure.
For other b, (12)--(18) still apply canonically, without presuming a
proper model or purity. No new pure exceptional phase is produced.
The known uniform-marker all-H detection is retained; no claim that
all marked orders are nonzero or that purity reflects from H(g_b) is made.

## 6. What this settles and what it leaves

There is no remaining freedom to recover marked functions outside R_f
by choosing a different proper phase model: its unit-adjoined ENTIRE
regular period image is exactly R_f for EVERY original phase. Its
unique global unit functional is intrinsic, in all degrees. Original-
volume extension is only a sufficient condition to include the unit
without adjoining it. The boundary lift epsilon_f is the precise additional
obligation; being able to retract the regular sector is insufficient.

The finite-base question for the universal bundle is more delicate.
Its unit volume extends over finite-base boundary, but that partial
model is NOT projective over the phase line. We cannot apply (6) or
replace its unit-adjoined image by R_f. At a finite-base crossing, residues on
individual affine-line components can be exact while intersection
gluing still contributes. Neither that raw calculation nor (9) proves
finite-base elimination or localizes (18) exclusively to finite base
or infinity. That localization remains OPEN.

NEXT FULL-CLASS target: use actual unweighted eventual-purity identities
and full original marked/global relations to force the corrected local
unit component, equivalently the established intrinsic safety test, to
vanish on EVERY original exceptional phase. A proof only for another
support class does not finish this task. No degree/parameter campaign,
Gaussian Bernstein duplication, universal software proof or new review
gate has been introduced. Both FULL SL2/2O and SL2/2I cases remain OPEN.
