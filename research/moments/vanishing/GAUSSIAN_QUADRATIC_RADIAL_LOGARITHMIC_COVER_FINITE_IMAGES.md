# Quadratic radial Gaussian phases on logarithmic covers

Let K be a characteristic-zero field. Use an original Witt Gaussian
pair Z,W and an independent ordinary Gaussian Y, with E[ZW]=E[Y^2]=1
and the other basic pairings zero. Put S=ZW. Consider

    P=P0+R, P0=A(Y)S^2+H(Y)S+G0(Y),                (L1)

where H,G0 are arbitrary polynomials and every monomial of R has
strictly positive Z-minus-W weight, or strictly negative weight after
reversing the pair. R=0 is allowed, with no degree restriction on R.
For ORIGINAL polynomial markers define

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P.

The new scalar obstruction is the following.

**Logarithmic-cover scalar obstruction.** Suppose A has degree at least
two and every root of A has multiplicity at most two. Then, for ANY
H,G0 in K[Y], the ACTUAL original unit cannot have the form

    Phi_P[1]=p(t)exp(psi(t)),
    p in K[t], p(0)=1, psi in tK[t].               (L2)

The scalar exponential psi is unrestricted. In particular these phases
are neither pure nor safe and their full actual images are not finite.
This conclusion permits arbitrary strictly one-sided additions R.

Combining this obstruction with the retained affine-coefficient theorem
gives a full finite-image classification. Allow A=0 or a nonzero A
whose every root has multiplicity at most two, including constants and
affine polynomials. Throughout this class,

    M_P finite iff A=H=0 and G0=c+bY.               (L3)

Every finite original vacuum is primitive, with

    Phi_P[1]=exp(ct+b^2t^2/2),
    M_P=exp(ct+b^2t^2/2)K[t].                       (L4)

For nonconstant P, V_P=M_P; for P=c constant, V_P=K exp(ct).
Safety holds exactly when P0=0 and therefore implies purity. In the
ZERO-boundary subfamily G0=0, polynomiality of the original unit,
purity, safety and zero core are equivalent. Arbitrary-boundary
purity in the constant/affine-A corner is not classified by (L2).
Triple or higher roots of A and cubic radial terms are not covered.

The proof retains the original Gamma endpoint, the full relative
extension, all puncture cycles and Gaussian Stokes corrections. The
key estimate is on absolute chain measure on a fixed logarithmic
curve; a bound obtained only after cancelling scalar periods would
not suffice. The written proof and its comparison imports have their
stated internal-review scope. Independent replay, external review
and formal verification are not asserted.

## 1. Original marker reduction and the retained branches

For a core in K[S,Y], every unbalanced original monomial Z^iW^jY^ell,
i!=j, has zero period. Every balanced one is precisely S^iY^ell, with

    E[S^iY^(2j)]=i!(2j-1)!!,
    E[S^iY^(2j+1)]=0.

Thus K[S,Y] represents the entire original core marker family. The
[angular-core theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md)
preserves both the embedded full module and the original unit:

    M_(P0+R)=M_P0, Phi_(P0+R)[1]=Phi_P0[1].         (L5)

It suffices to prove the scalar obstruction for P0.

For the unified full-image theorem (L3), the
[general actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
supplies a scalar identity p exp(psi) whenever M is finite. The
[boundary-image theorem](GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md)
then makes G0=c+bY and fixes the finite exponent to ct+b^2t^2/2.
Neither theorem assumes that p=1. For degree A<=1 the complete
finite criterion and the zero-boundary scalar criterion are the
[affine-coefficient quadratic-radial theorem](GAUSSIAN_QUADRATIC_RADIAL_AFFINE_COEFFICIENT_FINITE_IMAGES.md).
For A=0 they also follow from the retained
[S-affine theorem](GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md).
The new work proves (L2) for degree A>=2 without imposing any boundary
or exponential normalization beyond the displayed scalar hypothesis.

## 2. Exact relative source with arbitrary G0 and psi

Work first over C. Put B=C(t)[S,Y] and write

    D_S=partial_S-1+t(2AS+H),
    D_Y=partial_Y-Y+t(A'S^2+H'S+G0').

The actual original Wick functional kills D_Y B and D_S(SB), exactly
and coefficientwise. The S factor in the relative primitive preserves
the Gamma endpoint S=0. The exact top relative de Rham source is

    H_original=B/(D_Y B+D_S(SB)),                   (L6)

with complex SB -> B dS+SB dY -> B dS wedge dY. Its time connection
after the scalar gauge is partial_t+P0-psi'(t). Rational-time
amplitudes have finite Laurent principal parts at t=0, so the formal
functional is defined on this entire source. No convergence of other
original marked periods is assumed.

On A!=0 put q=1-tH-2tAS. The exact relation

    [qS^n]=n[S^(n-1)], n>=1,                        (L7)

follows from D_S(S^n)=nS^(n-1)-qS^n. It reduces the relative S source
to the independent classes v=[dS] and e_B=[q dS]: a nonzero D_S(SF)
has highest S degree at least two. The class e_B is the actual
endpoint class, not an exact class in this relative complex.

Put r=1/t, consistently choose sqrt(r)=1/sqrt(t), and on the
normalized square-root cover define

    u=sqrt(tA)v, x=(H-r)/(2sqrt(rA)),
    G=G0(Y)/r-Y^2/2-psi(1/r), F=G-x^2.             (L8)

The exact connection in BOTH Y and t is

    nabla e_B=dG e_B,
    nabla u=dF u-dx e_B.                            (L9)

One can verify (L9) entirely by (L7). The Y connection on v is obtained
by reducing (-Y+tG0'+tH'S+tA'S^2)v; the time connection is the
reduction of (AS^2+HS+G0-psi')v. Scaling by sqrt(tA) gives the diagonal
coefficient dF and off-diagonal coefficient -dx. In particular no
assumption A''=0, affine G0, or prescribed psi enters this identity.
The original seed density after S integration is

    eta=sqrt(r) omega, omega=dY/sqrt(A).            (L10)

## 3. Localization, descent and the actual unit quotient

Let D be the REDUCED support of {A(Y)=0}. Its components are smooth
S-lines, each transverse to the original boundary S=0. Repeated
factors of A do not change this support. Relative Gysin gives

    H^0(D,D intersect {S=0};exp(phi|D))
      -> H_original -> H_localized,                (L11)

where phi=-S-Y^2/2+tP0-psi(t). The left term is zero: a horizontal
function on a connected S-line vanishing at S=0 is identically zero.
Hence the original source injects into the localized source. Its
localization cokernel and residue classes are retained.

This is not evaluation of inverse spatial markers by the Wick rule.
At an ordinary time point, a horizontal functional on a differential
submodule of the original source extends linearly to the localized
source fibre and then horizontally by analytic parallel transport.
The extension need not be rational.

On the punctured Y-line the normalized cover V^2=A is finite etale,
possibly disconnected. Pullback and trace divided by two identify
the original source with its invariant descent summand. The deck
operation sends u and x to their negatives and fixes e_B; it respects
(L9). All projectors have constant coefficients. Extra covering-curve
boundary cycles are not substituted for original-line cycles.

The ordinary-point period input is Hien--Roucairol, *Integral
representations for solutions of exponential Gauss--Manin systems*,
Bull. SMF 136(4) (2008), 505--532, Theorems 2.7 and 3.5
([published primary text](https://www.numdam.org/item/10.24033/bsmf.2564.pdf)).
They give perfect absolute exponential pairings and representation
of local holomorphic horizontal functionals away from finitely many
parameter values. Apply the natural pairings and Stokes' formula to
the relative de Rham and rapid-decay chain cones of the original
smooth affine pair. This induces the relative perfect pairing used
here. It is not a direct citation about an arbitrary rank-two system.
After deleting finite parameter poles and the fixed spatial poles,
the phase is regular on a smooth affine total space. Time and phase
are algebraically independent because of the nonzero S^2 term; on
the endpoint, -Y^2/2 gives the corresponding independence.

Under the scalar hypothesis (L2), the gauged original formal functional
takes w=v/p to one. Its restriction to the rational cyclic differential
span C(t)<partial_t>w takes the successive covariant iterates to
1,0,0,... and therefore has a rational horizontal unit quotient.
Every rational relation is respected by the actual formal functional.
The [seed-functional lemma](../../../results/gaussian-boundary/seed_period_comparison.tex),
Lemma lem:seedfunctional, now permits ordinary-point extension from
that cyclic source, through (L11), and relative period representation.
Its seed period would be identically one. Sections 4-8 exclude this.

## 4. A fixed logarithmic normalization, including positive genus

Write A=D0(Y)^2 B0(Y), with B0 squarefree and D0 containing once each
double root. The multiplicity hypothesis implies that D0 and B0
have disjoint roots. Normalize V^2=A by T^2=B0 and V=D0 T.
This is a FIXED smooth projective curve C after adding infinity.
If B0 is constant, its normalization is two disjoint smooth components.
Otherwise it may have positive genus; no rational parametrization is
assumed.

The fixed differential omega=dY/V is logarithmic on every component:

* At a simple A-root use the branch coordinate T. Since Y-y0 is a
  unit times T^2, omega is regular.
* At a double A-root the normalization is unramified and omega has
  a simple pole, with constant residue on each sheet.
* If L=degree A is even, use Y^-1 at each infinity point. The order
  of omega is L/2-2, logarithmic for L=2 and regular for L>=4.
* If L is odd, use Y^-1/2 at its infinity point. The order is L-3>=0.

All poles of F,G,x belong to the fixed preimages of {A=0,infinity}.
Regular punctures are retained even when their phase has no pole.
F is spatially nonconstant for ANY H,G0,psi: its coefficient of r is
the nonconstant function -1/(4A). G is also spatially nonconstant,
since its parameter-constant spatial term is -Y^2/2. A scalar psi
cannot remove either fact.

For clarity, if A=a2((Y-m)^2-a^2), a2*a!=0, the normalization has

    Y=m+(a/2)(z+z^-1),
    V=sqrt(a2)(a/2)(z-z^-1),
    omega=dz/(sqrt(a2)z).                           (L12)

For a double quadratic root, V=+/-sqrt(a2)(Y-y0) on the two components
and omega=+/-dlog(Y-y0)/sqrt(a2). These exact examples illustrate
the logarithmic geometry; higher genus uses the same local orders.

## 5. The complete fixed-curve logarithmic chain estimate

Here is the fixed-curve version of the retained
[fixed-puncture logarithmic seed bound](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#fixed-puncture-logarithmic-seed-bound).
Let C be any fixed smooth projective complex curve, D a fixed finite
set, omega a fixed meromorphic differential with only simple poles
in D, and f(r,z) a rational parameter family with spatial poles in D.
The zero-exponential seed block has a complete geometric chain frame
whose absolute integrals are O((1+|log r|)^N) on a generic sector.
Multiplication by a scalar of positive r-order makes them tend to zero.
The same statement holds componentwise for a disconnected fixed curve.

The three exact comparison inputs are those retained in the
[primary manifest](../../../results/gaussian-boundary/sources/laurent_primary_manifest.json):

1. Roucairol, Section 1.2 and Lemma 1.4: the regular formal summand is
   recovered by nearby cycles with proper-image compatibility, and
   a vertical exponential phase pole has no moderate nearby cycles.
2. Sabbah, Section 1.b, formula (1.4) and Proposition 1.5: the localized
   finite/pole and finite/finite node complexes retain the gluing maps
   and monodromy.
3. Hien--Roucairol, Theorems 2.7 and 3.5: perfect exponential period
   pairing and ordinary-point period representation.

These inputs are local on curves or compatible with proper image,
and do not impose genus zero. They do not themselves provide the
quantitative logarithmic estimate; the following argument does.

Mark D and all branches of the fibres f=0,1,infinity. Split their
convergent Puiseux sections after a finite parameter cover and separate
their collisions by blowups. The contact orders are finite and
decrease, so this terminates. The resulting proper marked model
resolves the phase: an unresolved zero/pole crossing would also
contain a value-one branch. Auxiliary sections delete no new points
of the original connection.

Near a fixed pole a write omega=c_a dw/w plus a holomorphic differential,
with c_a constant. The pullback of w is a monomial in the normal-crossing
coordinates times a holomorphic unit. Its logarithmic differential
therefore stays logarithmic, with no negative r-order from a moving
residue. The holomorphic part pulls back regularly. On a fibre through
uv=r its relative coefficient is bounded by a constant times |du/u|.
This local argument replaces the P1 partial-fractions expansion. The
holomorphic differentials and handles of positive-genus components
remain in compact finite-phase cells with bounded coefficients.

The full nearby-cycle complex consists of all finite-phase components,
their original puncture circles, finite/finite annuli, and descending
attachments at finite/pole nodes, with all gluing maps and monodromy.
Triangulate the compact pieces and keep the entire annulus complexes.
Proper comparison gives a finite chain frame of the regular block,
with constant coefficients on a small parameter sector. This retains
all genus cycles and all punctures, not a selected set of stationary
points.

The estimates hold on ABSOLUTE chain measure. On a compact finite-phase
piece the phase and differential coefficients are bounded. At a node,
the logarithmic coefficient stays bounded and the bridge length is
O(log(1/|r|)). Choose fixed coordinate discs at all punctures. A single
power |r|^M contains the finitely many coefficient and Puiseux transition
scales. In the annuli between |r|^M and a fixed radius, descending
inverse-phase paths are real algebraic curves of fixed degree in
fixed algebraic charts. The
[curve-length estimate](../../../results/gaussian-boundary/analytic_core.tex),
Lemma lem:curvelength, bounds each dyadic annulus by a constant in
logarithmic measure; there are O(log(1/|r|)) such annuli. Finite
algebraic covers give the same local bound. There are finitely many
branches and no unbounded retracing on a strictly descending segment.

Beyond the last transition at an ACTUAL phase pole, its leading polar
monomial dominates and |omega/df|<=C/|f|. The scalar leading phase
coefficient cancels in this ratio. A descending tail is bounded by
an integral exp(-s) ds/s, with truncation after the phase has become
large and negative. A puncture where f is regular is instead represented
by its puncture-circle and adjacent finite/finite complexes; it is not
treated as an irregular tail. These estimates complete the full-chain
bound. Finite logarithmic factors from monodromy do not defeat a
strictly positive scalar r-order.

## 6. Incomplete descending chains and regular puncture values

The same estimate applies to an open inverse-f path starting at an
algebraic/Puiseux point z0(r) with Re f(z0)<=C. Choose a leftward
phase-plane path to an actual phase pole, detouring around both the
finite critical phase values AND the finite f-values of all regular
fixed punctures. Small detours stay in Re f<=C+C'. If the start is
critical, its ramified inverse has an integrable algebraic endpoint.

Choose the cutoff exponent M larger than the orders of all moving
start/critical points and all inverse contact orders with regular
punctures, as well as coefficient transition scales. All are algebraic
convergent Puiseux data, so one finite M works. Phase-plane detour
radii may also be chosen as Puiseux powers. No exponentially small
spatial scale is introduced.

The compact and dyadic-annulus estimates of Section 5 apply. The finite
value detours prevent an open chain from ending at a regular puncture:
its inverse lift stays polynomially away and costs at most logarithmic
length. Regular puncture circles remain separate frame members. True
tails go to actual phase poles, where |omega/df|<=C/|f| and exponential
decay controls the continuation. Thus

    integral_open omega exp(f)=O((1+|log r|)^N).     (L13)

This is needed when a Gaussian sector jump begins at bounded phase
rather than at exponentially negative phase. It includes incomplete
paths and never assumes that an eliminated S-coordinate stays bounded.

The distinction between regular and irregular punctures matters even
in this application. If H=A Y^2 and G0=A Y^4/4, then exactly

    F=-r/(4A)-psi(1/r).                             (L14)

For A=Y^2, omega=dY/Y has a logarithmic pole at infinity while F is
regular there. Its finite F-value must be included among the phase
detours above; no false polar-tail estimate is applied at infinity.

## 7. Two entire Gaussian sectors with their exact jump

Let E(x)=integral_0^x exp(z^2)dz and define

    J_+(x)=exp(-x^2)(i sqrt(pi)/2-E(x)),
    J_-(x)=exp(-x^2)(-i sqrt(pi)/2-E(x)).

Then J'=-2xJ-1 and

    J_+-J_-=i sqrt(pi)exp(-x^2).                   (L15)

On the fixed closed sectors

    -pi/6<=arg x<=7pi/6, and
    -7pi/6<=arg x<=pi/6,

respectively, |J_+/-|<=C/(1+|x|). These sectors cover the whole x
circle and include the closed upper/lower half-planes strictly inside
their asymptotic domains. To prove the uniform bound, use

    I(z)=integral_0^infinity exp(-s^2-2zs)ds,
    J_+(x)=iI(-ix), J_-(x)=-iI(ix).

For |arg z|<=3pi/4-delta rotate the s ray to
theta=-sign(arg z)(pi/4-delta/2), taking theta=0 at arg z=0.
The integrand is entire and quadratic decay permits the rotation.
On s=exp(i theta)v,

    Re(s^2)>=sin(delta)v^2,
    Re(zs)>=sin(delta/2)|z|v.

The resulting integral is bounded by the minimum of a fixed Gaussian
constant and 1/(2sin(delta/2)|z|). Take delta=pi/12. In particular
the bound is uniform at x=0; no large-x assumption is needed.

The exact local horizontal values of (L9) are

    h_B=C_B exp(G),
    h_u=exp(F)(C_A-C_B E(x)).

Hence the endpoint density of v is eta J_+/- exp(G), while its absolute
Gaussian density is eta exp(F). Changing endpoint sectors adds
EXACTLY a constant multiple of the absolute density. The identity
J_-(-x)=-J_+(x) compensates the sign of 1/sqrt(A) and makes the endpoint
density descend to the original Y-line.

## 8. The complete relative frame for arbitrary G0 and psi

The Y long exact sequence has the form

    0 -> H^1(boundary G) -> H^1(relative source)
      -> H^1(absolute Gaussian F) ->0.              (L16)

The H^0 terms vanish because F,G are spatially nonconstant and hence
have no algebraic horizontal exponential coefficient. There is no
H^2 on the affine curve. This sequence retains the full endpoint
extension. For nonaffine G0 it is not restricted to the one ordinary
Gaussian contour of an affine boundary.

Take ONE common proper marked model of F,G,x: mark the fixed punctures,
all branches of the zero, one and pole fibres of F and G, and the
zero, one and pole sections of x. Split them after a finite parameter
cover and separate their collisions. Extra markings delete no new
spatial points. This common refinement retains the full G nearby-cycle
complex by proper-image compatibility; omega remains logarithmic on
the same refinement.

Use the complete zero-exponential boundary-G chain frame supplied by
Section 5. Compact finite-phase pieces and finite/finite necks have
bounded phase; descending G tails can be chosen at fixed imaginary
phase. The absolute measure |omega|exp(Re G) on their algebraic chain
pieces has O((1+|log r|)^N) mass. This includes shrinking regular-
puncture loops and critical pieces escaping to infinity. The bound
is on chain measure, not a cancellation of a scalar omega-period.

Choose J_+ on upper-half-plane x pieces and J_- on lower pieces. The
pointwise bound in Section 7 gives endpoint seed contributions
O(sqrt(r)*(1+|log r|)^N). There are uniformly finitely many switches:
the real-x equations on the algebraic chain pieces have fixed degree.
Their points have Puiseux coordinates after a generic parameter ray
and finite cover. At a switch,

    Re F=Re G-x^2<=Re G<=C.

Close the exact Gaussian fibre jump by transporting a closed Gaussian
thimble over an open inverse-F chain to an actual phase pole. The
critical and regular-puncture phase detours of Section 6 are retained.
By (L13) each complete closing surface costs only
O(sqrt(r)*(1+|log r|)^N). These corrections need not be exponentially
small and are not discarded. Continuously chosen cycles on a small
generic sector are flat relative homology classes.

The construction may be performed on the base Y boundary source,
lifting its chain pieces to the normalized cover. A simple-root loop
need not lift to a closed single-sheet path; the J_-(-x) identity and
the half-density sign give the exact descent. The invariant base
boundary classes are retained before taking the Gaussian deck summand.
Estimating all lifted chains does not replace the original source
with extra boundary classes of a higher-genus cover.

Now take the zero-exponential part of (L16). Formal exponential
decomposition is exact. Its boundary image is exactly the boundary
zero block just lifted. The absolute zero block has the bound of
Section 5. The zero part of any lift of a nonzero-exponential boundary
factor lies in the absolute block. A full absolute frame together
with all boundary lifts has an invertible block-triangular period
matrix; all factors, gluing maps and monodromy have been retained
before taking the zero summand. No regular semisimplicity is assumed.

Consequently every zero-exponential original seed period has a bound

    O(|r|^(1/2)*(1+|log r|)^N).                     (L17)

Division by p adds deg p>=0 to this order. The constant seed period
one furnished by Section 3 has zero-exponential component one, which
contradicts (L17). This proves the scalar obstruction (L2) for arbitrary
H,G0,p,psi.

The argument is first over C. The finitely generated field of the
phase and scalar-unit coefficients embeds into C, preserving degree,
nonzero coefficients and root multiplicities in characteristic zero.
The original formal identity transfers coefficientwise. Thus (L2)
holds over every characteristic-zero field.

## 9. Unified finite images and all original markers

For degree A>=2 the new scalar obstruction excludes full finite image
by the actual-image theorem. The retained affine-coefficient theorem
handles the remaining degrees, including A=0. This proves necessity
in (L3). Sufficiency is explicit: when P0=c+bY, let gamma>=1 be the
least positive angular weight of a nonzero R after choosing its sign.
For every original marker Q of spatial degree D,

    Phi_P[Q]=exp(ct+b^2t^2/2)
       *sum_(n=0)^floor(D/gamma) t^n/n!
           E[Q(Z,W,Y+bt)R(Z,W,Y+bt)^n].             (L18)

For R=0 only n=0 is used. All summands are polynomial, and Q=1
leaves only n=0. Thus the vacuum is primitive and M is exactly the
module in (L4). The literal fixed-marker theorem supplies V=M for
nonconstant P, with the constant exception already stated.

Safety makes (L4) polynomial and therefore c=b=0 over a field.
Conversely (L18) then proves safety. In the zero-boundary subfamily,
a polynomial unit has scalar normalization psi=0. The new obstruction
and the retained zero-boundary affine-A criterion force A=H=0.
Purity and safety are consequently equivalent there, with zero core.
For nonaffine G0 in the constant/affine-A corner, no purity converse
is inferred from this theorem.

The safe-recomposition ideal follows from the
[general safe-recomposition theorem](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md):
it is (U-c) exactly when A=H=0 and G0=c is constant, and is the zero
ideal otherwise. The safe original cutoff is floor(D/gamma), with
the zero-phase R=0 case separate. With b!=0 a valid normalized degree
bound is D+(1+deg_Y R)floor(D/gamma); these are upper bounds, not an
arbitrary-R filtered-surjectivity claim. The
[mixed-cutoff theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md) gives
the corresponding finite-jet and certified-tail consequences.

For every fixed coefficient support with degree A<=2, G0=0 and bounded
degree H, all A-root multiplicities are automatically at most two.
The field pure locus is therefore the origin. Its rational moment
ideal has radical the coefficient maximal ideal and finite colength;
some finite prefix contains a power m^nu. No effective new prefix,
nilpotence index or length is claimed for the quadratic-A slices.
Any actual containment certificate gives coefficient nilpotence over
Q-algebras and the marked cutoff nu-1+floor(D/gamma). It does not make
coefficients zero or vacua primitive over nonreduced rings. The same
global finite-colength claim is not made for all higher-degree A
coefficient spaces: their triple-root strata are outside this theorem.

## 10. Remaining scope and verification

After the proved boundary reduction, every invariant core of original
spatial degree at most six has support

    delta S^3+S^2 A2(Y)+S H4(Y)+bY+c,
    deg A2<=2, deg H4<=4.

The entire delta=0 branch is now classified for full finiteness.
Outside it and the
[retained quadratic-Y family](GAUSSIAN_VARIABLE_RADIAL_QUADRATIC_FIBRE_FINITE_IMAGES.md), the remaining support
requires delta!=0 and deg H4>=3, with known parity subfamilies still
available independently. A triple root of a higher-degree A gives
a volume pole of order greater than one on the normalized cover;
the logarithmic argument does not classify that case.

Run the
[focused checker](verify_gaussian_quadratic_radial_logarithmic_cover_finite_images.py)
with `--output gaussian_quadratic_radial_logarithmic_cover_finite_images_verification.json`
from this directory. Its
[receipt](gaussian_quadratic_radial_logarithmic_cover_finite_images_verification.json)
pins this proof, exact source controls and all indispensable retained
imports, including the three cached primary texts. The checker verifies
the full arbitrary-function source connection, bounded original Wick
controls with nonaffine boundaries, logarithmic normalizations of
positive-genus and disconnected covers, exact Gaussian jumps and
the regular-infinity cancellation (L14). Triple roots are rejected
by the scoped geometric backend.

The unrestricted logarithmic chain comparison, joint-model completeness,
open-chain filling and scalar exclusion are the written proof above.
Finite tests do not certify those all-order statements. The existing
primary-source and internal-review limits are unchanged; independent
replay, formal verification and external review are not asserted.
