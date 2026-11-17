# Global torus bundle and unit-invisible infinity poles

16 November 2026. Internal written proofs with retained geometric inputs.
Independent replay, formal verification and external review are not claimed.
Both COMPLETE exceptional MZ cases remain **OPEN**; review is deferred.

This treats the ENTIRE original phase algebra and EVERY original marker,
without a support or degree restriction. A global line-bundle model removes
the need to trivialize the torus torsor rationally. Its graph compactification
preserves the FULL actual marked-period function module. For every phase
with positive AND negative torus weights, its ORIGINAL unit volume has no
finite-base boundary pole, regardless of ramification or critical collisions.

The remaining volume poles lie above the base infinity conic. An explicit
all-degree family of actual mixed-weight pure SAFE phases has a permanent
logarithmic pole there, while its corrected local unit component is zero.
Thus neither exceptional invariance nor mere pole presence distinguishes
the obstruction. No conclusion about vanishing of the full component vector
is inferred from finite-base volume regularity.

See the [checker](verify_haar_bundle.py), [receipt](haar_bundle_verification.json),
[source audit](haar_bundle_source_audit.json) and retained
[full ownership obligation](exceptional_mz_ownership_20261116.json).

## 1. The complete original algebra on one global bundle

Let Z0=P1_top times P1_bottom parameterize the two row lines ell1,ell2
of a projective 2-by-2 matrix. Let B be the diagonal, and X0=Z0 minus B.
Its points are ordered distinct lines; X0 is the affine sphere quadric Q2.
Put

    L=Hom(ell2,ell1)=O(-1,1),
    Z=P(O direct_sum L), pi:Z->Z0.

Let D0,Dinf be the two disjoint sections. Then there is an algebraic
isomorphism

    PGL2 = L^times over X0
          = Z minus (pi^-1 B union D0 union Dinf).     (1)

Indeed a projective matrix with rows u,v determines their lines and
the nonzero map t:ell2->ell1 with t(v)=u. Conversely choose any nonzero
v in ell2 and form the matrix with rows t(v),v. Its projective class is
independent of the choice, and its determinant is nonzero exactly off B.
These constructions are regular in the usual line-bundle charts.

The effective left diagonal torus scales t by its character rho=lambda^2.
There is no global fibre coordinate on the nontrivial bundle. Every full
regular function on PGL2 nevertheless has a finite Laurent expansion in
local fibre coordinates; its coefficients glue as sections of powers of L.
Thus (1) retains the charged coefficients globally, instead of introducing
poles by a rational frame. Every phase and marker on SL2/Gamma for
Gamma=2O or2I is a right-barGamma-invariant regular function on this
PGL2 open, since the binary centre acts trivially. This is full-algebra
coverage, including zero-free three-weight and arbitrary larger supports.

The right action of barGamma=O,I extends to Z, acting on both row lines
and on their Hom line. On B the line is End(ell), canonically trivialized
by the identity map. Hence barGamma acts TRIVIALLY on that fibre scalar.
It acts nontrivially on the base B=P1. The orbit bounds6/12 therefore
spread right-base roots but do not spread a fixed fibre root t=t0.
This distinction will account for the control in Section5.

## 2. The original normalized volume and EVERY marked period

In affine row-line coordinates xi,eta choose row representatives (1,xi),
(1,eta) and fibre ratio t. A projective matrix is

    [x:y:z:w]=[t:t*xi:1:eta], Delta=t(eta-xi).

The invariant PGL2 volume on this chart is, up to its fixed nonzero
normalizing constant,

    dx wedge dy wedge dw/Delta^2
       = dxi wedge deta/(eta-xi)^2 wedge dt/t.         (2)

Write omega=dxi wedge deta/(eta-xi)^2 and Theta=omega wedge dlog t.
A change of fibre frame adds a base 1-form to dlog t, whose wedge with
the base top form is zero. Hence Theta is GLOBAL. On Z its exact divisor is

    div(Theta)=-2*pi^*B-D0-Dinf.                       (3)

There are no hidden frame divisors or substituted dualizing generators.
Formula (2) is an explicit computation of the ORIGINAL invariant volume.

The compact core K=PSU2=SO3 is the circle bundle over the anti-diagonal
real sphere M in X0. Normalize omega on M and dlog t on each effective
circle so the fibre integral is1. Fubini for this oriented circle bundle
then gives integral_K Theta=1. This does NOT replace K by M times S1:
the bundle has Chern number of absolute value2 on that sphere. Indeed
the first tautological line has degree-1, while the second is pulled
back through the orientation-reversing antipodal map and has degree+1;
their Hom line has degree-2.

For every original phase f, every original marker Q and every polynomial R,

    h_Gamma(Q R(f))=integral_K Q R(f) Theta.            (4)

This follows from the retained normalized Haar-volume identity; passing
from SU2 to its central quotient and integrating right-invariant functions
does not insert a Jacobian. Conversely any PGL2 marker can be averaged
over barGamma, without changing its marked period for invariant f. Thus
using the cover adds NO new actual marked-period functions.

Resolve the projective graph of f:Z -> P1, as an isomorphism on PGL2,
and restrict to finite phase values. Denote the smooth model by Y and
the phase by fbar:Y->A1. It is projective; pi extends to Y. Keep the
original open PGL2 and its full boundary, including polar divisors of
all marker forms. Then (4) holds on this ONE model for every Q,R.
In particular all actual Cauchy germs

    A_Q(s)=integral_K Q Theta/(s-f)

are exactly the original ones. Their derivatives, all analytic
continuations and ALL global function relations agree. Therefore the
FULL actual marked-period FUNCTION module, with literal1 retained
BEFORE quotienting by C[s], is preserved. This is not an assertion
that every Q Theta extends holomorphically to Y, nor an identification
of the entire geometric cohomology of Y with the original source.

Every polynomial H(f) is covered as well. Nonconstant H has a pole only
at infinity, so finite H(f) valuations are exactly finite f valuations;
compose the same projective model with H. No primitive-purity reflection
or scalar branch cancellation is used.

## 3. ALL finite-base unit-volume poles disappear for every mixed phase

Suppose f has at least one positive and one negative effective torus
weight. At the generic point of D0 its lowest Laurent term gives a
pole of f; at the generic point of Dinf its highest term gives a pole.
Thus neither strict transform of these sections remains in Y over a
finite f value.

For ANY boundary divisorial valuation E over Z put

    b=ord_E(pi^*B), a0=ord_E D0, ai=ord_E Dinf,
    A=A_Z(E),
    ord_E Theta=A-1-2b-a0-ai.                          (5)

This is the exact discrepancy formula for (3). The sections are disjoint,
so a0 and ai cannot both be positive.

If E lies over the FINITE BASE X0, then b=0. A remaining boundary E
cannot be D0 or Dinf itself, since those map to infinite phase. Its
centre is a PROPER subvariety of one of the sections, say D0. Choose
local parameters x1=0 for D0 and a second parameter x2 vanishing at
the generic centre. Then

    A>=ord_E x1+ord_E x2>ord_E x1=a0.

For completeness the elementary inequality follows by pulling back
dx1 wedge dx2 wedge dx3: each factor contributes its valuation, and
in their wedge only one dt/t can occur at a uniformizer of E. The
Jacobian order is at least the sum of these valuations minus1. Adding1
gives the displayed bound. Generic codimension-one centres give only
the original divisor itself. Integrality now gives A-1-a0>=0 in (5).

Consequently

    for EVERY mixed-sign f in the FULL algebra,
    Theta is HOLOMORPHIC at EVERY finite-base
    boundary divisor over EVERY finite phase value.   (6)

This covers all degrees, arbitrary coefficient vanishing and arbitrary
critical-point collisions. Inside the original open the graph is an
isomorphism and Theta is already regular; further smooth modifications
cannot create a pole of a regular form. No isolated-critical-point,
simple-branch or finite support hypothesis has been imposed. The same
conclusion holds simultaneously for every nonconstant H(f).

For a one-sided phase, the opposite generic section instead has finite
phase value and order-1 for Theta. This is compatible with known pure
safe phases and explains why the mixed-sign assumption in (6) is needed.
When purity removes the neutral part, one-sided safety is already known;
we do not re-register that safe class as new.

## 4. The exact remaining infinity condition and its limitations

For mixed-sign f all possible unit-volume poles on Y are now above B.
If every E with b>0 and finite f satisfies

    A>2b+a0+ai,                                       (7)

then (5) is nonnegative everywhere, and the retained
[projective original-volume period theorem](SEMISTABLE_LOG_PERIODS_AND_SINGLE_ISOTYPIC_SAFETY.md)
excludes eventual purity of EVERY nonzero H(f). Thus every nonconstant
eventually pure mixed phase must have an actual finite-phase valuation
above B violating (7). This is a necessary pole statement, not an
existence or failure certificate.

The remaining proper-volume obstruction is not rational-frame
preservation or finite-base critical collisions: those parts are now
settled on the ENTIRE algebra. It is holomorphic volume extension at
infinity. But even complete knowledge of that polar set is weaker than
the corrected LOCAL UNIT component. The following actual controls prove
that one cannot identify the two.

## 5. Actual pure SAFE mixed phases with a permanent logarithmic pole

Let F be ANY nonzero right-Gamma-invariant homogeneous binary form of
even degree 2D>0. Such forms exist in every multiple of8 for2O and12
for2I by taking powers of a minimum-degree form. Put on SL2

    f_F=F(x+z,y+w).                                    (8)

This is an ORIGINAL regular invariant phase. It is the left-SL2
translate of the known highest-weight phase F(x,y), so it is pure in
EVERY positive order and SAFE against EVERY original marker, with
sufficient bound 2D*m>deg_matrix Q after the linear left transport.
This safe class and its automorphism transport are RETAINED results,
not a new pure-family discovery.

Relative to the standard torus, both extreme terms F(x,y) and F(z,w)
are nonzero, so f_F is mixed-sign. In fact its left V_(2D) orbit has
all 2D+1 weights in this unipotent translate; neither original degree
nor support size is bounded in the family.

In the projective bundle chart, with R=eta-xi and W=t+1, (8) reads

    f_F=F(W,W*xi+R)/[(W-1)^D R^D].                    (9)

Blow up the smooth curve R=W=0, lying above B at fibre scalar t=-1.
In the chart W=R v the exceptional divisor E is R=0 and

    f_F=R^D F(v,v*xi+1)/(Rv-1)^D,
    Theta=unit* dxi wedge dR wedge dv/[R(Rv-1)].       (10)

The numerator F(v,v*xi+1) is not identically zero: its arguments
parametrize a dense open set of the binary plane. Hence

    ord_E f_F=D>0,  f_F(E)=0,
    ord_E Theta=-1                                   (11)

for EVERY such form and EVERY degree. Equivalently A=2,b=1,a0=ai=0
in (5). Right invariance cannot spread this fibre root: by Section1
barGamma acts trivially on the scalar t over B. The leading term has
the fixed factor (t+1)^(2D), however widely its right-base roots spread.

The pole is INTRINSIC to this finite-f valuation. It cannot disappear
on another smooth model proper over A1_f carrying the same rational
Theta. The valuative criterion gives E a centre on any such proper
model; a common smooth resolution would pull a holomorphic form to a
form of nonnegative E order, contradicting (11).

Nevertheless the retained
[full intrinsic local criterion](FULL_INTRINSIC_LOCAL_TATE_CRITERION.md)
gives ZERO corrected local unit component at0 for f_F, since it is
pure and safe. The same pole survives for every nonzero polynomial
H(f_F) with H(0)=0. Each such phase is still one-sided under the
transported torus and pure safe, and E remains over value0. Thus the
countercontrol includes all such recompositions, not just primitive f_F.

This disproves BOTH proposed shortcuts:

* Exceptional mixed-sign phases automatically have a proper ORIGINAL-
  volume model over their phase.
* A finite-phase logarithmic volume pole detects a nonzero corrected
  local unit/torus-gluing component.

It does not disprove MZ or give a new failure: these are actual SAFE
phases with unit-invisible poles. The family lies in a transported
nullcone; no automorphically closed residual failure is constructed.

## 6. The FULL theorem obligation after this universal reduction

Original unweighted and FULL marked-function preservation now has a
single model for EVERY actual phase, with no rational frame. Every
possible full failure is mixed-sign, so (6) covers its finite-base unit
volume. The remaining positive target is still: unweighted purity
forces ZERO corrected local unit component/torus-valuation vector on
the ENTIRE residual class. The negative target is still an actual
invariant pure phase with a fixed original bad marker.

Two implications have NOT been proved and must not be hidden:

1. Finite-base regularity of Theta has NOT been shown to remove the
   corresponding directions functorially from the FULL actual unit
   extension. Marker forms may still have boundary poles. A reduction
   of the full component calculation to the infinity boundary requires
   an actual source/period-kernel preservation proof.
2. Infinity volume poles have NOT been separated into unit-invisible
   and unit-obstructing directions in general. The all-degree control
   (11) proves that bare pole order or orbit size cannot make this
   distinction. One must retain corrected gluing data, global relations
   and all original markers, then use unweighted purity.

The earlier safety/proper-model dichotomy remains a sufficient route
but its FULL coverage is unproved. A blanket proper-model alternative
is explicitly false by (8)--(11). The known strict leading-root and
automorphically closed residual constraints remain in force; this
global model does not prove that every residual phase contracts under
an automorphism. General Gaussian/formal Bernstein work remains owned
by the separate thread and is used only as input here. Both COMPLETE
exceptional MZ cases remain OPEN. Review is DEFERRED, not the next gate.
