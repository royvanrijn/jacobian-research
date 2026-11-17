# Complex homogeneous Gaussian growth and the complete Borel trichotomy

16 November 2026. Internal written proofs using a named classical regularity
theorem. No novelty, independent replay, external review or formal verification
is claimed. The [checker](verify_complex_homogeneous_growth.py),
[receipt](complex_homogeneous_growth_verification.json) and
[source audit](complex_homogeneous_growth_source_audit.json) delimit the evidence.

The [previous arithmetic theorem](HOMOGENEOUS_MOMENT_RIGIDITY_AND_EXCEPTIONAL_MARKER_OBSTRUCTIONS.md)
gives a lower bound at EVERY nonzero sufficiently late index for algebraic
coefficients. Here a different proof establishes an ALL-COMPLEX theorem in
every dimension and homogeneous degree. Its lower conclusion is a LIMSUP,
not a pointwise lower bound or a statement about every nonzero subsequence.

For a homogeneous complex phase P of degree D>=1 and an arbitrary fixed
complex polynomial Q, write A_n=E[Q P^n] for a standard real N-Gaussian.
Either A_n is eventually zero, or its factorial growth order is exactly D/2.
More precisely the generalized Borel series

    B_s(z)=sum_(n>=0) A_n z^n/(n!)^s,    s real,

has radius zero for s<D/2, finite positive radius for s=D/2, and is entire
of exact order 1/(s-D/2) for s>D/2. In particular, when D>2, convergence
of the ordinary exponential moment series at any positive radius is already
equivalent to eventual vanishing of this particular marked tail.

This is a global growth theorem, not a full Gaussian or exceptional MZ
classification. It supplies no algorithm deciding whether a given tail is
eventually zero, and it does not remove the Jacobian weight in the exceptional
Bezout transfer. The genuine mixed-isotypic exceptional problem remains open.

## 1. Compact algebraic resolvents cannot be transcendental entire functions

**MVCOMPACTREGULAR20261116.** Let X be a smooth complex algebraic variety of
dimension k, P a regular function, omega a regular algebraic top form, and
gamma a compact oriented real k-cycle on X(C). Finite complex linear
combinations of such cycles are allowed. For any regular Q, put

    H(z)=integral_gamma Q omega/(1-zP).

This defines a convergent germ at zero. The germ belongs to a finite-rank
algebraic differential system regular singular at every point of P1,
including infinity. If this selected germ extends to an entire function on C,
it is a polynomial. Thus its ordinary coefficient sequence is either eventually
zero or has finite positive convergence radius.

The same conclusion holds for finite sums of such germs, their derivatives,
and polynomial differential operators applied to them. Cancellation between
the summands is allowed; no conclusion about each separate summand is needed.

Proof. Work over K=C(z) with the smooth open family

    U={(x,z) in X times A1:1-zP(x)!=0}.

It is smooth over A1 because it is open in the product of a smooth variety
and A1. On the generic fibre Q omega/(1-zP) is a closed algebraic top form.
Its class belongs to finite-dimensional algebraic de Rham cohomology. On a
nonempty open subset of the parameter line the integration cycle determines
a local horizontal homology class: near zero it is the SAME compact cycle,
since the denominator is uniformly nonzero there. If zero is an exceptional
parameter for the cohomology presentation, use a nearby ordinary parameter;
the original germ already extends across zero by its compact integral.

The classical input is Deligne's open local monodromy theorem, stated as
Theorem14.3 in Katz, *Nilpotent connections and the monodromy theorem*,
Publ. Math. IHES39(1970), printed pp230-231. It applies to the smooth generic
fibre without a properness hypothesis and makes its Gauss-Manin connection
regular at every place of C(z). Only this regularity assertion is used.
Katz Remark12.3, printed p221, gives the analytic local fundamental matrix
as a power matrix times powers of logarithms. A rational de Rham section
adds only finite pole orders. The retained primary PDF is
[here](sources/katz_nilpotent_connections_1970_regular_periods_20261116.pdf).

Consequently each continuation of H has at most power-logarithmic growth
in a sector at infinity. If H is entire, finitely many sectors cover a large
circle, using fixed continuations from its one entire branch. There are
constants M,C with |H(z)|<=C(1+|z|)^M there after absorbing logarithms into
one extra power. Cauchy's estimate makes every Taylor coefficient of index
greater than M zero. This proves polynomiality. The initial compact integral
has positive radius; a nonpolynomial germ therefore has finite positive
radius, by the contrapositive. Direct sums of the regular systems handle
sums, and differentiating a period pairs the SAME flat cycle with a
covariant derivative of its rational section. Hence all the stated operations
retain the same regularity argument. For k=0 the integral is a finite sum
of rational functions and the conclusion is immediate.

Neither generic D-finiteness alone nor convergence near zero suffices for
this argument. For instance exp(z) is entire and D-finite but is irregular
at infinity. We are using the actual compact algebraic resolvent, not an
arbitrary exponential period or a formal source presentation.

## 2. Every homogeneous marked normalization is such a regular period

Let P be homogeneous of degree D>=1 in N>=1 variables, with arbitrary complex
coefficients. Decompose Q=sum_l Q_l into homogeneous pieces and set

    k_n=ceil(Dn/2), R_n=2^k_n (N/2)_(k_n), b_n=A_n/R_n.

On the unit real sphere use normalized measure sigma. Gaussian polar
integration gives, when Dn+l is even,

    E[Q_l P^n]
      =2^((Dn+l)/2) (N/2)_((Dn+l)/2)
          *integral_sphere Q_l P^n d sigma;                         (1)

the term is zero when Dn+l is odd. This formula also holds for the two-point
sphere N=1. No reality of P or Q is used.

For D even, only l=2a contributes. Dividing (1) by R_n multiplies its sphere
moment by the polynomial

    r_l(n)=product_(j=0)^(a-1) (N+Dn+2j).                           (2)

For D odd, put epsilon=l mod2 and a=(l-epsilon)/2. Only n congruent to
epsilon modulo2 contributes, and its multiplier is

    r_l(n)=product_(j=0)^(a-1) (N+Dn+epsilon+2j).                   (3)

Empty products are1. Let H_l(z)=integral_sphere Q_l/(1-zP) d sigma and
theta=z*d/dz. For D even the normalized generating function is

    F(z)=sum_n b_n z^n=sum_(l even) r_l(theta) H_l(z).              (4)

For D odd it is

    F(z)=sum_l r_l(theta)
             *(H_l(z)+(-1)^l H_l(-z))/2.                          (5)

These identities keep all inhomogeneous marker pieces and both index
parities; replacing the marker by its top homogeneous part is unjustified.

For N>=2 the complex affine quadric sum x_i^2=1 is smooth and the
restriction of sum_i (-1)^(i-1) x_i dx_1 wedge ... omit dx_i ... wedge dx_N
is a regular top form whose real sphere restriction is its volume form.
Thus every H_l is a compact algebraic resolvent from section1, multiplied
by a harmless fixed normalization constant. For N=1 it is rational.
Equations(4)-(5) and section1 make F a regular period germ at all points
of P1. Its initial convergence also follows directly by bounding |P| on
the compact sphere and absorbing the fixed polynomials r_l(n) into an
exponential bound.

We have proved the exact dichotomy

    F polynomial <=> F entire <=> A_n eventually zero.            (6)

Otherwise 0<rho(F)<infinity. By Cauchy-Hadamard there are C>1 and infinitely
many n with |b_n|>=C^(-n), and there is an upper bound |b_n|<=C^n for all
sufficiently large n. The lower estimate is only on SOME infinite subsequence.

## 3. Complete generalized Borel growth in every real normalization order

**GMCHOMBOREL20261116.** In the notation above, if the marked tail is not
eventually zero, then

    limsup_(n->infinity) log|A_n|/(n log n)=D/2,                   (7)

where log0=-infinity. Stirling's formula gives

    log R_n=(D/2)n log n+O(n),
    log(R_n/(n!)^(D/2))=O(n),

and the nth-root of the latter ratio has a finite strictly positive limit.
Together with the finite positive radius of F this proves the claimed
critical Borel radius and (7). For s<D/2, the infinite lower-bound
subsequence gives nth roots tending to infinity in B_s, hence radius zero.
For s>D/2, the upper bound gives an entire series. The coefficient formula
for its order,

    order(B_s)=limsup n log n / log(1/|A_n/(n!)^s|),

with zero coefficients contributing0, equals 1/(s-D/2): the upper estimate
bounds it from above and the SAME infinite lower-bound subsequence attains
the bound. This formula follows from Cauchy's coefficient estimate and
optimizing the maximal term, or the usual power-series characterization of
entire order. Thus the order is exact, even if infinitely many coefficients
are zero or exceptionally small.

If the tail is eventually zero, every B_s is simply a polynomial. For D>2,
take s=1: the ordinary exponential moment series has positive radius iff the
marked tail is eventually zero. At D=2 that series instead always converges
near zero, and a nonzero infinite tail has finite positive radius. At D=1
an infinite tail gives an entire exponential moment function of exact order2.

An immediate extension of exact-law rigidity is useful: if an exact prescribed
non-eventually-zero sequence t_n has finite factorial LIMSUP order tau and

    E[QP^n]=c^n t_n, c!=0, for all sufficiently large n,

then D=2tau. The target coefficients need not be rational or algebraic.
This follows directly from(7); it does not specialize infinitely many
inequalities or remove transcendental coefficients by a finite prefix.

## 4. Scope, controls and the exceptional target

The checker compares direct Wick expectations with the parity-correct
sphere formulas for selected odd/even degrees and inhomogeneous markers.
It checks exact rational examples at D=1,2,3,4,6, polynomial differential
operators in(4)-(5), regular and irregular scalar equations including
infinity, and the Borel threshold/order calculator. These are normalization
and boundary controls, not a finite computational proof of regularity or
asymptotic growth. No asymptotic fit or degree/support search is used.

The earlier arithmetic result remains stronger at each individual nonzero
index for algebraic inputs and remains unchanged. Neither proof licenses
a complex pointwise lower bound. Homogeneity is substantive: for example
P=Z^4+T in a circular isotropic coordinate Z independent of a real standard
Gaussian T has exactly the moments of T, hence order1/2 despite degree4.
This also rules out extending the theorem to arbitrary inhomogeneous phases
by reading only their highest degree.

Both exceptional quotients have the actual homogeneous four-Gaussian lift.
The result applies to it and every fixed polynomial marker, with no new
variable or source substitution. It gives a precise convergence criterion
for each marked tail. It does NOT prove that pure vanishing forces all those
tails finite, or that a weighted Bezout identity is an unweighted one.
The leading-nullcone/lower-opposite-weight cancellation implication remains
the full exceptional MZ gate. External review and outreach remain deferred
under the user's closure-first priority; assurance flags remain false.
