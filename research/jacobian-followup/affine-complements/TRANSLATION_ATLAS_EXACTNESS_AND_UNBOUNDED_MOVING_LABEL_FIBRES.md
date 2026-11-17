# Whole-plane translation atlases determine exactness, not target realization

PF2RULED138. 16 November 2026. Characteristic zero, geometrically.
This completes PF137's transition-residue test under precise affine
atlas hypotheses. An explicit all-degree family then shows why the
single-base-root conclusion cannot extend to moving root labels.
The full finite-normalization/trace-dual and Wright/HC4/JC2 goal is OPEN.

## 1. An exact criterion on a whole affine atlas

Let S be a smooth integral AFFINE surface with a finite open cover
U_i=A2_(x,w_i). Require the same global base x on every chart, and

    U_i intersection U_j = D(a_ij(x)) times A1_(w_i),
    w_j=w_i+h_ij(x), h_ij in k(x).                            (1)

The intersections are the ACTUAL chart intersections, not just generic
fibre identifications. Their rational transition cocycle satisfies
h_ij+h_jk=h_ik. The form Omega=dx wedge dw_i glues to a unit volume.
Then

    Omega has a GLOBAL REGULAR primitive
      iff every h_ij(x)*dx has zero finite residues.             (2)

An etale map S->A2 is therefore impossible if ANY such residue is
nonzero: units on S are constant, since they restrict to units on a
dense whole plane chart. As in PF137, an etale target would supply
an exact scalar multiple of Omega. Vanishing residues do NOT construct
a target pair, its closed coframe, a finite order or trace regularity.

Necessity is PF137 Section2's direct argument. On U_i let
theta_i=-w_i*dx. If Omega=dtheta, the polynomial Poincare lemma gives
theta-theta_i=dH_i with H_i polynomial. Consequently

    d(H_j-H_i)=theta_i-theta_j=h_ij(x)*dx.                     (3)

The left function has zero fibre derivative, so is rational in x.
Rational derivatives have zero residues at every finite pole.

For sufficiency, fix a reference chart0. Partial fractions and the
residue hypothesis give rational R_i with R_i'=h_0i; set R_0=0.
Put R_ij=R_j-R_i. Then R_ij'=h_ij, and the R_ij satisfy the EXACT
additive cocycle identity on all triple intersections. A pole of a
rational primitive always gives a pole of its derivative. Thus R_ij
has no pole where h_ij is regular, and is regular on the intersection
in (1). No integration constants on triples are left unresolved.

The R_ij form an O_S-valued Cech1-cocycle. Since S is affine,
H1(S,O_S)=0, and this affine-intersection cover computes it. Hence
there are REGULAR H_i on U_i with R_ij=H_j-H_i. Equations(3) now
show that theta_i+dH_i glue to the required global regular theta.
Their derivative is Omega on every chart, proving (2).

The classical cohomology inputs are explicitly imported from
[Stacks Lemma30.2.2, tag01XB](https://stacks.math.columbia.edu/tag/01XB)
and Lemma30.2.6 in [tag01X8](https://stacks.math.columbia.edu/tag/01X8).
Their complete statements and proof passages were read and cached;
the affine vanishing theorem is not independently re-proved here.
No classification of affine-line surfaces is used. AFFINENESS of S
and the whole cover in (1) are essential retained hypotheses.

## 2. Application to the actual cyclic source-image universal cover

Use HC4MYMP0/MQ0's notation: G=C_g is any nontrivial ACTUAL generic
deck subgroup, W=union_(gamma in G) gamma(X), X the ORIGINAL plane,
and V=q(X). W and V are affine and W->V is finite etale. There is a
global affine source-coordinate ruling y on W with reduced fibres;
each translated chart selects one whole unit A1 in each fibre.

Normalize its chart coordinates to (y,w_gamma) so that the ACTUAL
volume is dy wedge dw_gamma. On the generic A1 fibre the transition
is affine linear. Equality of the volumes makes its slope1; hence
it is a rational translation in y. Two selected whole components
of a fibre are equal or disjoint, by MQ0. Accordingly each actual
pair intersection is a base complement of finitely many values times
A1. If a rational translation is regular at a base value, separatedness
identifies both chart maps there: they agree on the dense generic open
and therefore on the whole local plane. Conversely a pole prevents
intersection there. Thus this is exactly an atlas of form (1).

The ACTUAL target already makes Omega exact, so all these residues
vanish. Conversely, on such an affine translated-chart cover those
residue conditions account for the WHOLE de Rham class of this Omega:
there is no further obstruction to its regular primitive. This closes
an exact-volume-only continuation, not the target realization problem.
It does not identify V with the entire intermediate etale locus or
replace any omitted ramified prime of the original full normalization.
It also does not infer that arbitrary A1-fibred surfaces have (1).

## 3. Exact free-additive surfaces with any number of degenerate fibres

Let b(x) be ANY square-free nonconstant polynomial, r=deg b, and
m>=2 ANY integer. Define

    a=b^m, sigma=(1-m)*b', R=b^(1-m), R'=sigma/a,
    S_(b,m)={a(x)*y=z*(z-sigma(x))}.                         (4)

At a root beta of b, sigma(beta)=(1-m)*b'(beta)!=0. Thus the
hypersurface is smooth: away from a=0 use its y derivative; over
a=0 its two distinct z roots have nonzero z derivative. It is integral
by the primitive linear-in-y equation, as in PF137.

There are two WHOLE plane charts:

    U_0: z=a*w,       y=w*(a*w-sigma),
    U_1: z=sigma+a*v, y=v*(sigma+a*v).                       (5)

They cover S. Their actual intersection is D(b) times A1, with
w=v+sigma/a=v+R'. The common unit form Omega=dx wedge dz/a is
dx wedge dw and dx wedge dv, respectively. Thus (2) already proves
exactness. Here is an explicit REGULAR primitive certificate as well.

Since gcd(sigma,b)=1, choose a polynomial B with

    B*sigma=1 modulo b^(m-1).

Set

    H_0=b*B*w,
    H_1=b*B*v+(B*sigma-1)/b^(m-1).                          (6)

Both are polynomials on their respective WHOLE charts. On the overlap,

    H_1-H_0=-R,
    (-w*dx)-(-v*dx)=-R'*dx=d(H_1-H_0).                    (7)

Consequently -w*dx+dH_0 and -v*dx+dH_1 glue, with derivative Omega.
This avoids extending a merely rational primitive through its poles.
It gives exactness for every r and m by one Bezout identity, not by
a bounded coefficient experiment.

The following derivation is globally defined and locally nilpotent:

    delta(x)=0, delta(z)=a, delta(y)=2z-sigma,
    delta^2(z)=0, delta^2(y)=2a, delta^3(y)=0.               (8)

It preserves the equation and has no zero: use delta(z) off a=0
and delta(y)=+/-sigma!=0 over a=0. Its flow is explicitly

    z -> z+t*a, y -> y+t*(2z-sigma)+t^2*a, x -> x.

It is a FREE Ga action. On U_0 it is the original source translation
partial_w and preserves the entire chosen plane chart.

The x-ruling has exactly r degenerate fibres, at the ARBITRARILY
prescribed distinct roots of b. Each is TWO disjoint REDUCED A1s;
every other fibre is A1. U_0 omits exactly one closed A1 at each
of those r values. Hence units are constant and

    chi(S_(b,m))=1+r, pi1(S_(b,m))=1, Pic(S_(b,m))=Z^r.     (9)

Euler additivity proves the first formula. For pi1, general position
moves every loop off the added smooth complex divisors, so pi1(U_0)
surjects onto pi1(S), and the former is trivial. For Pic, use the
normal divisor localization sequence: the r omitted divisors generate,
there are no relations because units on U_0 are constants, and
Pic(U_0)=0. These are the same whole-plane topology/divisor arguments
retained in HC4MYLA0; smoothness identifies Pic with the class group.

For r>=2 this disproves the inference that a smooth affine surface
with an open whole plane, scalar units, exact unit volume and a FREE
source-preserving additive action has at most one degenerate affine
ruling fibre. It also disproves any bounded fibre-count conclusion from
these conditions alone. The root labels in (4) MOVE with x; PF137's
separated-variable single-root theorem is not contradicted.

## 4. The retained target and assurance boundaries

There is NO etale map S_(b,m)->A2 constructed here, no claimed generic
deck group and no actual full finite normalization for these controls.
They are NOT Keller counterexamples or realized intermediates. Passing
(2),(6)--(9) does not give the special primitive P*dQ or two global
functions with Jacobian one. No target invariance follows from (8).

PF137 still excludes every separated a(x)y=p(z) model with multiple
distinct base roots. PF130's other exact x-dependent mixed(2,3) control
is credited and remains outside that separated theorem. MP0 eliminates
actual noncyclic decks; MQ0 eliminates actual source-image point gaps.
The present criterion and all-r family show why the next step must use
the ACTUAL finite target order, its trace-dual multiplication or full
marked-energy integrality, rather than a further exact-volume or fibre
count constraint alone. Trivial decks also remain in scope.

The checker gives small exact cocycle/Bezout/chart/flow controls and
pins the primary cohomology source. Written Sections1--3 prove the
unbounded statements. No classification, novelty, formal verification,
independent whole-proof replay, external review or full-goal closure
is asserted; predecessor mathematical checkers are not executed.
