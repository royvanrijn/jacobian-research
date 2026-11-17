# Nonhomogeneous leading divisors and purity-adapted Ga globalization

17 November 2026. Written uniform regularity and safety theorem in all
degrees. Both FULL exceptional Mathieu--Zhao statements remain **OPEN**.
External review is deferred; independent replay and formal verification
are not claimed.

The previous rational-Ga theorem extends from single positive-weight
seeds to the ENTIRE strict positive-weight sector, with arbitrary mixed
Peter--Weyl degrees and right multiplicities. A constant-volume birational
comparison preserving the Ga fibration up to a regular base automorphism
cannot have poles if its image of such a seed is regular. It is an actual
global volume-preserving automorphism, so every original marker and the
full corrected extension transport.

The new input is a divisor multiplicity estimate for arbitrary invariant
polynomials, not only homogeneous binary forms. The resulting pure phases
remain in the known automorphic-nullcone safe class. Coverage of all
surviving mixed pure phases is not established.
Section3a also treats EVERY pure phase: all but one Ga direction force
a global comparison. This transports the corrected class; it does not
make that class zero without a safe endpoint or another argument.

## 1. Every invariant affine curve divisor obeys the orbit-degree bound

Let Gamma=2O or2I act linearly on C^2, and put B=C[z,w]^Gamma. Let

    s0=6 on2O,       s0=12 on2I,

the minimum projective root-orbit sizes. Suppose 0!=ell in B has total
ordinary polynomial degree N. It need NOT be homogeneous. For EVERY
prime divisor E of the punctured Klein base U, one has

    ord_E(ell)<=N/s0.                                  (1)

Here ell is the actual polynomial on the bottom plane, so N has no
ambiguity from determinant representatives or rational coordinates.

**Proof.** Pull back to C^2 minus0. The quotient is etale, so choose an
irreducible polynomial factor q of ell over E and write M=ord_q(ell)
=ord_E(ell). Invariance makes this multiplicity the same on every
distinct Gamma-translate of q. Let Q be the product of those distinct
translates, normalized up to constants. Then Q^M divides ell, and Q
is semi-invariant: Gamma permutes its irreducible factors and changes
their product only by a nonzero scalar.

Its highest homogeneous part Q_L is also semi-invariant, of degree
L=deg Q>0. Its projective zero set is nonempty and Gamma-stable. It
therefore contains an orbit of at least s0 points. Hence L>=s0, and

    N>=M L>=M s0.

This proves (1). No classification of the affine factor q, its genus,
or its singularities was used. The factor need not pass through the
origin. All its translates and their multiplicities are retained. QED.

In particular, applying homogeneous root bounds directly to a general
affine factor would be unjustified; it is the semi-invariant orbit
PRODUCT and its homogeneous leading part that provide the estimate.

## 2. A numerical globalization criterion for arbitrary regular endpoints

Keep the original Ga torsor X=SL2/Gamma -> U and

    D=z partial_x+w partial_y,
    d=deg_D f>=1,       ell=D^d f/d!,       N=deg ell.

Let Phi be a birational self-map preserving this fibration up to ANY
REGULAR automorphism of U, and suppose

    Phi^*Omega=lambda Omega,       lambda in C^*,
    f and g=Phi^*f are ORIGINAL regular phases.

The [sharp two-translate theorem](TWO_INTEGRAL_GA_TRANSLATES_AND_EXCEPTIONAL_RATIONAL_SHEAR_CLOSURE.md)
reduces the base automorphism by a global linear-normalizer lift. That
lift preserves Omega, Ga degree and total polynomial degree of ell.
In a regular local slice, a pole of order r of the remaining affine
fibre term costs leading multiplicity at least r ceil(d/2). By (1),

    s0 r ceil(d/2)<=N.                                (2)

Consequently the sufficient condition

    N<s0 ceil(d/2)                                    (3)

forces every local fibre term regular. The map globalizes across ALL
base divisors. The same theorem uses the nonzero original Haar line
to force the fibre linear part to be one: any other global affine
linear part would produce a global Ga slice and make Omega exact.
Thus Phi is a GLOBAL Omega-preserving automorphism.

This criterion allows nonhomogeneous and nonreduced ell. It is a new
usable degree form of the divisorial comparison, not a numerical bound
on all possible phase degrees. Failure of (3) does not produce an
obstruction or a counterexample. It only leaves rational poles possible.

All divisor multiplicities in (1)--(2) are measured on the ORIGINAL
punctured base. No finite-fibre boundary graph, special-fibre component
or corrected Neron vector is replaced by this curve-divisor data.

## 3. Every strict positive-weight seed satisfies the criterion

Suppose f is any nonzero original phase with ONLY strict positive
left-torus weights. Its support size and matrix degree are unrestricted.
In a Peter--Weyl summand of highest weight n and actual torus weight w,
the Ga degree of that component is

    j=(n+w)/2.

Put d=deg_D f. A component contributes to D^d f only if j=d. Its
bottom polynomial degree is therefore

    n=2d-w.

Every contributing w is positive, and all weights are even on these
quotients. Hence the complete, possibly NONHOMOGENEOUS, ell satisfies

    N=deg ell<=2d-2<2d.                               (4)

Distinct summands and right multiplicities cause no problem: D^d
retains only the components with j=d, and their degrees all obey (4).
At least one survives by the definition of d.

Since s0>=4,

    2d<=s0 ceil(d/2),

so (4) implies (3). We have proved:

**Positive-cone globalization theorem.** On BOTH exceptional quotients,
let f have arbitrary strict positive weight support. Every constant-
volume birational self-map preserving the Ga fibration up to a regular
base automorphism whose phase image Phi^*f is regular is an ACTUAL
global Omega-preserving automorphism.

This conclusion requires no purity hypothesis: f is already pure SAFE
by weight selection, and regularity of the image forces globalization.
The rational fibre term can be fully inhomogeneous. The base
automorphism need not be fixed in advance.

## 3a. Purity globalizes every comparison except possibly one Ga direction

There is an ALL-PURE-LOCUS consequence, independent of positive support
for the full phase. Let f be any nonzero EVENTUALLY PURE original phase
on either exceptional quotient, and let n be its highest matrix degree.
The [strict-leading theorem](STRICT_LEADING_NULLCONES_AND_FULL_THIN_CLOSURE.md)
gives a UNIQUE common leading left factor L^m with m>n/2. Lower
isotypic components and opposite weights are retained.

Choose ANY left Ga root derivation D such that D(L)!=0 on linear left
binary forms. On each nonzero leading binary coefficient L^m h,
the Ga degree is at least m: degree for a locally nilpotent derivation
is additive on products in a characteristic-zero domain, and L has
degree one in this direction. Explicitly the highest derivative of
a product is a nonzero binomial factor times the product of the two
highest derivatives. Thus

    d=deg_D f>=m>n/2.

The original D preserves the matrix filtration, and ell=D^d f/d! is
a bottom polynomial of degree N<=n. Hence

    N<=n<2d<=s0 ceil(d/2).                             (4a)

The argument of Sections1--2 applies in this left-conjugated Ga
fibration: right Gamma and its orbit sizes are unchanged. Every
constant-volume birational comparison preserving THIS fibration up to
a regular base automorphism, and taking f to a regular phase g, is
therefore a GLOBAL Omega-preserving automorphism. This covers EVERY
eventually pure f and EVERY such nonannihilating Ga direction, with
no degree or support bound.

There is precisely ONE exceptional left Ga direction, up to scalar:
the nilpotent direction with kernel line L, namely D(L)=0. Thus any
NONGLOBAL comparison of the stated type with a pure phase endpoint
must use this direction. When its Ga degree is positive it must also
satisfy N>=s0 ceil(d/2), by (2). These are necessary conditions for a
rational comparison obstruction, not a failure producer.

The exceptional direction cannot simply be deleted from the statement.
For the original pure SAFE bottom invariant f=F_n(z,w), n=8 or12,
D=z partial_x+w partial_y annihilates f. The nonglobal rational shift
with parameter 1/F_n has constant volume Jacobian one and fixes f.
It has a regular phase endpoint despite its poles. The phase remains
SAFE; this is a control for the comparison theorem, not a Mathieu
counterexample.

For a nonexceptional direction the global comparison transports ALL
original fixed markers and the full corrected component. It does NOT
by itself prove that the component is zero: both endpoints could carry
the same unknown component. If one endpoint is independently known
safe, then (5) and the intrinsic criterion prove safety of the other.
Thus a rational safe normal form in ANY such direction now suffices
without a separate marker-regularity or unit-retention gap.

This directional localization concerns poles and nonglobality of
Ga comparisons. It is NOT support localization of the intrinsic
corrected Neron vector, and it does not identify an ambient tree with
the finite-fibre graph or justify contraction of good components.
Existence of a suitable normal form and vanishing of the actual
coupling in the remaining direction are still unproved.

## 4. Full original-marker safety and all polynomial recompositions

Let g=Phi^*f in the preceding theorem. For ONE FIXED original regular
marker Q and EVERY m,

    h(Q g^m)=h(((Phi^-1)^*Q) f^m).                     (5)

The inverse marker is regular and independent of m. If a is the minimum
positive weight in f and B_Q the largest absolute weight of this inverse
marker, then a*m>B_Q is a sufficient cutoff. Thus g is pure and safe
for all original markers, regardless of mixed support in its original
coordinates or genera of the phase boundary components.

For every polynomial H,

    H(g) eventually pure iff H(0)=0;
    when H(0)=0, H(g) is safe for all original markers. (6)

The unit moments equal H(0)^m by the same weight calculation transported
through Phi. If H(0)=0, H(f) has only strict positive weights, giving
safety by (5). This includes the zero polynomial; a nonzero constant
fails purity. No arbitrary-ring or marker-degree invariance under a
nonlinear automorphism is assumed.

Because Phi is a genuine isomorphism of the ORIGINAL phase complements
over the SAME phase parameter and preserves Omega, it transports the
full source, canonical quotient and corrected boundary class. On this
pure class the intrinsic extension splits and the corrected cycle
coupling vanishes. This conclusion follows with every original
component retained. It does not rely on contracting an ambient tree
or discarding the leading-divisor localization data.

The actual ordinary marked resolvents on this class are finite Laurent
polynomials, by the proven tail vanishing. No algebraic exponential
generating function is substituted for an ordinary resolvent.

## 5. Genuine multiweight seeds and the exact coverage boundary

The retained original phases

    F_8(top)+G_12(top) on2O,
    F_12(top)+G_20(top) on2I

are examples with different positive weights and different isotypic
degrees; see the [primitive safe controls](ACTUAL_BRANCH_SUPPORT_AND_PRIMITIVE_SAFE_BIFURCATIONS.md).
The theorem globalizes every regular comparison of the stated kind
from these entire seed phases, not merely separate comparisons of
each summand. It allows all other positive components as well.

Their original one-sided safety was already closed. The new result is
the regularity/globalization bridge for their full rational fibration
comparisons. Thus the new construction coverage is strictly broader
than the previous single-weight seed statement, while the resulting
safe points remain in the known automorphic-nullcone class.

A seed with several NEGATIVE weights is not covered by Section3.
If it is eventually pure, Section3a still applies in its nonexceptional
Ga directions, provided the proposed map respects that fibration.
Its leading coefficient can have arbitrarily large degree relative
to d, so (4) fails. The previous negative SINGLE-weight theorem is
retained; it cannot be extended by simply reversing the torus while
leaving the map's Ga fibration fixed. Likewise arbitrary mixed seeds,
nonregular birational base changes and nonconstant volume Jacobians
remain outside the present coverage.

For the full exceptional theorem one still needs a proof covering
every surviving pure phase, or an actual pure phase with ONE FIXED
original bad marker. A rational normal form to a positive seed in
this compatible fibration would now automatically be a regular
unit-retaining normal form. No such universal normal-form existence
or universal corrected-coupling vanishing is proved here.

All capacity, automorphism, primitive-fibre, relative Ga and full
family-transfer inputs remain in force. No finite degree campaign,
extra decidability result or external review is introduced. Both
FULL exceptional Mathieu--Zhao cases remain OPEN.
