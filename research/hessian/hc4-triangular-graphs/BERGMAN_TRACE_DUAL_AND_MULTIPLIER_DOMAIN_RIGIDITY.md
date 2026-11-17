# Bergman trace duality and the multiplication-domain gap

2026-11-16. HC4MYLW0/1. ALL DEGREES, complex polynomial PLANE Keller maps.

HC4 and the plane Jacobian conjecture remain open.
The remaining unproved implication is: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
This note addresses multiplication of the ACTUAL finite trace dual directly.
It proves that ordinary square integrability is automatic, identifies the
missing fourth-power estimate, and shows exactly how a polynomial that is
not integral destroys density of a Bergman multiplication domain. The last
density assertion needed for actual Keller coordinates is NOT supplied.

## 1. Setup and one bounded target set seeing the full boundary

Write K=(P,Q):X=A2_(x,y)->Y=A2, Jac(K)=J in C*, B=C[P,Q], R=C[x,y],
d=[Frac(R):Frac(B)], and pi:Z=Spec S->Y for the FULL finite normalization.
Let Lambda=Hom_B(S,B) be embedded in R by field trace. Retained KU1 gives
finite local freeness and

    omega_Z=Lambda*Omega, Omega=dP wedge dQ=J dx wedge dy.

KN0 gives the WHOLE-plane open X in Z. The
[quotient-point application of Flenner--Zaidenberg Theorem 2.1](FULL_NORMALIZATION_QUOTIENT_POINTS_AND_CANONICAL_L2_EXTENSION.md)
gives FINITE QUOTIENT singularities at every omitted singular point.
That point theorem is sufficient for the global square-integrability
proof below; cyclicity and the separate unverified KQ0 differential
source-layer theorem are not needed here. No higher-dimensional version
is inferred. KN0 gives

    div_Z(Omega)=sum (e_i-1)E_i.

Choose a bounded open target BALL V whose closure is compact and which
meets pi(E_i) for EVERY omitted prime. Such a ball exists because there
are finitely many primes. The intersections contain generic smooth points:
exceptional algebraic subsets of each image curve are finite. If the
boundary is empty any nonempty bounded ball is allowed. Put U=K^-1(V).
The open U may be disconnected; its Bergman Hilbert space

    H=A^2(U)={f holomorphic on U : integral_U |f|^2 dV_X<infinity}

is defined componentwise. All omitted primes are retained in this choice,
including unramified ones. For any nonnegative measurable F, change of
variables off the algebraic exceptional set gives

    integral_U F dV_X
       =|J|^-2 integral_V sum_(K(z)=w) F(z) dV_Y.    (1)

There are d sheets almost everywhere, hence Vol(U)=d|J|^-2 Vol(V)<infinity.
This identity uses finite generic degree and local biholomorphy, NOT
properness of K. The finite pi ensures pi^-1(closure V) is compact in Z.

## 2. LW0: the polynomial L2 vectors are EXACTLY the finite trace dual

We have the GLOBAL embedded identities

    R intersect H = Lambda,
    {r in R : integral_U |r|^(2d)<infinity} = S.     (2)

At a generic omitted prime choose smooth analytic coordinates (u,v), with
E={u=0}. If r has pole order q>0 there, (1), or the source volume itself,
gives the local radial test

    integral |r|^p dV_X
       comparable to integral_0^epsilon rho^(2e-1-pq) d rho,
    finite iff p*q<2e.                             (3)

Units are bounded above and below after shrinking the chart. At p=2 this
is exactly q<=e-1, the codifferent bound. Consequently r in R intersect H
satisfies every height-one canonical-module condition and lies in Lambda.
The ball sees every omitted prime; source divisors cause no polynomial pole.

For the converse Lambda*Omega is a regular reflexive canonical form on Z.
At smooth points its coefficient is holomorphic and its canonical L2 density
is locally integrable. At a finite quotient point take a finite small
quotient chart C2/G->Z. Away from the origin the chart is etale, so the
pullback is a holomorphic top form there. Hartogs extends its coefficient
over the origin. Its canonical L2 integral is finite on a smaller cover
chart; division by |G| gives the integral downstairs. Compactness of the
finite preimage of closure V gives the global norm. On X the canonical
density is |J|^2 |r|^2 dV_X. Thus Lambda is contained in H, proving the first
identity including ALL finite singular points. No Grauert--Riemenschneider
or Lp extension theorem is newly imported.

For the second identity, elements of S are bounded on pi^-1(closure V)
and belong to every finite Lp(U), since U has finite volume. If r is not
in S it has a pole q>=1 at an omitted prime. Its ramification index e<=d,
so p=2d makes p*q>=2e and (3) diverges. This proves (2). The cutoff is
uniform and sufficient; no optimality is claimed.

Formula (3) is a GENERIC-divisor statement for arbitrary p. For p>2 it is
NOT asserted to exhaust possible extra integrability restrictions at cyclic
quotient points. The global p=2 statement has the separate quotient proof;
the p=2d detection uses boundedness for S and generic poles for its converse.

## 3. A bounded Hilbert trace retraction exists AUTOMATICALLY

For f in H define Tf off the exceptional target curve by summing its values
over the d inverse branches. It is single-valued and holomorphic there.
Cauchy--Schwarz and (1) give

    ||Tf||_(A2(V)) <= sqrt(d)*|J|*||f||_H.          (4)

A holomorphic function with finite L2 norm off a complex hypersurface
extends across it: Laurent coefficients exclude negative normal powers at
generic points, and Hartogs handles the remaining finite set. Thus T is a
bounded map H->A2(V). For Cg=g composed with K,

    ||Cg||_H^2=d|J|^-2 ||g||_(A2(V))^2,
    TC=d Id, C^*=|J|^-2 T,
    Pi=(1/d)CT is the ORTHOGONAL projection onto C(A2(V)). (5)

These are direct change-of-variable identities; they do not assume that
the original map U->V is proper, and are not a new import of Bell's proper
mapping theorem. In particular an analytic trace retraction is automatic
even in the putative noninverse case. It is NOT the algebraic retraction
R->B forbidden by KO2: the polynomial vectors on which this Hilbert trace
acts are exactly Lambda, rather than all R. No polynomial density in H or
extension of T to all source polynomials is asserted.

## 4. The missing multiplication estimate is EXACTLY L4

For ANY finite B-generating family lambda_1,...,lambda_m of Lambda,
the following conditions are equivalent:

1. every lambda_i belongs to L4(U);
2. every element of Lambda belongs to L4(U);
3. Lambda*Lambda subset Lambda;
4. the ordinary pairing Tr(Lambda*Lambda) is contained in B;
5. the Keller map is invertible.

Coefficients in B are bounded on V, so 1 implies 2. If 2 holds, every
lambda^2 belongs to R intersect H=Lambda. Polarization gives 3. A finite
B-module closed under multiplication is integral over B by the determinant
trick. Since S subset Lambda and S is the full integral closure, 3 implies
Lambda=S; hence its ordinary trace pairing is regular. The retained KU/KN
inverse criterion gives 4 iff 5. Conversely an inverse makes U relatively
compact in the source, so every polynomial belongs to L4(U). This also
establishes the equivalence to 3 without assuming a primitive presentation.

Equivalently define the POSITIVE Hermitian energy, off the exceptional set,

    E(w)=sum_i sum_(K(z)=w) |lambda_i(z)|^2.

Then E is AUTOMATICALLY L1(V), whereas E^2 is L1(V) if and only if the inverse
exists. Indeed the square of a finite sum of nonnegative numbers is bounded
below by their squares and above by their number times the sum of squares;
(1) identifies the resulting integrals with the fourth powers in 1.

The ordinary bilinear traces Tr(lambda_i lambda_j) are already L1(V), by
the same inequality and Cauchy--Schwarz. Meromorphic L1 does NOT imply
regularity: a simple complex pole is L1. L2 on the BASE would remove all
such poles, but ordinary source L2 only yields L1 for products. In a supplied
finite free basis this is precisely the inverse-Gram regularity gap. Positivity
removes cancellation from E; it does not automatically improve its exponent.

## 5. LW1: natural polynomial multiplication domains are either total or NOT DENSE

For r in R let M_r have its NATURAL maximal domain

    Dom(M_r)={f in H : r*f in H}.

For the single ball V chosen in Section 1,

    r in S
      iff M_r is bounded and defined on ALL H
      iff Dom(M_r) is DENSE in H.                  (6)

The forward direction follows from boundedness of r on the compact finite
preimage of closure V. For the converse suppose r has a pole of order q>0
along E. Choose a generic smooth p in E above V, avoiding all other boundary
components and singularities, and a polydisk W in Z about p with W intersect
X=W minus E. For any f in H, the form f*Omega on W minus E has finite
canonical L2 norm and therefore extends holomorphically across E. If its
coefficient in regular coordinates is F_f, mean value gives a BOUNDED
functional

    ell_p(f)=F_f(p), |ell_p(f)| <= C_p ||f||_H.      (7)

This functional is NONZERO: the finite canonical module has a global
polynomial section lambda*Omega not vanishing at p, and lambda belongs to H
by (2). If f belongs to Dom(M_r), r*f*Omega must also extend holomorphically.
The pole of r forces F_f to vanish to order at least q on E. In particular

    Dom(M_r) subset ker(ell_p) proper closed subset H. (8)

Thus its domain is NOT DENSE. This is stronger than failure of a bound.
It involves the actual finite normalization and an interior canonical
evaluation at a missing divisor, not an estimate at the outer boundary of V.

Consequently dense natural domains for multiplication by BOTH source
coordinates would force x,y in S, hence the complete inverse. Density for
all finite Lambda generators would likewise force Lambda=S and inversion.
These density statements are UNPROVED for actual Keller sources. They do
not follow from source polynomiality or the Hilbert trace retraction.

## 6. Sharp connected constant-Jacobian LAURENT controls, every e>=2

On X_e=Spec C[u,u^-1,v]=Gm times A1 take

    P=u^e, Q=v/u^(e-1), Jac_(u,v)(P,Q)=e.

Its full finite normalization is S_e=C[u,Q], finite free of degree e over
B_e=C[P,Q], and X_e is its open D(u), with v=u^(e-1)Q. The finite codifferent
is Lambda_e=u^(1-e)S_e and embeds in the actual LAURENT source. It has the
B_e-basis lambda_a=u^-a, 0<=a<=e-1. For a unit target bidisk its source
preimage is the generalized Hartogs triangle

    0<|u|<1, |v|<|u|^(e-1).

The source Jacobian is constant, the target volume has a GLOBAL primitive
P*dQ on the entire finite normalization, and every lambda_a is L2. But

    integral |u|^(-p*a) dV_(u,v)
       =2*pi^2 / (2e-p*a) when p*a<2e,

and it diverges at equality or above. Thus lambda_(e-1) fails L4 for EVERY
e>=2; e=2 gives the logarithmic threshold 4, so the uniform exponent is
sharp in local ramification geometry. Natural multiplication by u^-1 has
a NON-DENSE domain: it excludes the allowed top pole mode u^(1-e).

For a,b>0, a+b<=2e-2, the ordinary trace is

    Tr(u^(-a-b))=e/P if a+b=e, and 0 otherwise.

The simple pole is L1 on the base and not L2. A single square trace may
vanish by cyclic cancellation while a mixed trace has that pole. The
positive energy detects the full family, not one self-pairing.

This is a connected ACTUAL constant-Jacobian Laurent map and an actual
finite flat normalization/codifferent, not a hypothetical algebra packet.
It is NOT a polynomial map on the WHOLE plane: u is a nonconstant unit,
the formula for Q has a denominator, and the target divisor P0 has NO
retained source sheet. It is not an HC4/JC2 counterexample. Whole-plane
input is needed to improve the multiplication estimate.

## 7. Entire coordinate multiplication need not have dense domain, even on Runge sets

For EVERY integer a>=3, consider the connected open polynomial polyhedron

    W_a={|y|<1, |x^a*y|<1} subset C2.

It is complete Reinhardt and logarithmically convex, hence pseudoconvex.
Its Taylor series converge on compact subpolydisks; finite subcovers of a
compact set give uniform polynomial approximation, so it is Runge. It
contains the entire line y=0 and has finite volume. Direct integration gives

    ||x^m y^n||^2
      =pi^2 / ((m+1)*(n+1-(m+1)/a))
        iff m+1<a*(n+1),                           (9)

and divergence otherwise. Distinct monomials are orthogonal. In particular
Vol(W_a)=pi^2*a/(a-1), x is L2 for all a>=3, and x^(a-2) is an allowed
Bergman vector whose product with x is NOT L2. The coefficient functional
of x^(a-2)y^0 is continuous and nonzero. Every function in the natural
domain of M_x has that coefficient zero, since multiplication shifts it
to the forbidden monomial x^(a-1). Thus Dom(M_x) is NOT DENSE.

This whole-ambient-plane, Runge, finite-volume example rules out a generic
functional-analysis assertion of dense coordinate multiplication. It is
NOT an actual Keller pullback: the defining polynomial map (x^a*y,y) has
Jacobian a*x^(a-1)*y, which vanishes. No inference from this control to the
missing actual Keller density assertion is made.

## 8. Evidence and the remaining actual-source problem

The L2 theorem uses the audited finite-quotient point theorem and the
retained canonical-module input. Change of variables, Laurent removability, Hartogs, finite
matrix/module arguments and Hilbert mean value estimates are written out;
the point theorem is attributed in its linked application note, which
retains the exact primary PDF and scope of that import.
Bell-type trace identities and Bergman operator issues were checked against
primary literature for context, not imported as properness or density
theorems for this setting. No historical-priority claim is made.

New exact controls check the all-parameter Laurent Jacobian, radial and
monomial norms, threshold strictness, cyclic trace pairings and Hilbert
retraction scalings. They do not replay ANY registered mathematical checker.
Finite controls do not certify universal analytic arguments or the retained
boundary theorem; independent replay, formalization and external review
remain absent.

The new distinction is precise: the ACTUAL finite trace dual already is
the polynomial L2 space and has a bounded Hilbert trace, but algebraic
multiplication requires L4, or equivalently dense natural multiplication
domains for its generators. Neither improvement is proved automatic. A
successful proof must use whole-plane KELLER structure beyond Runge geometry,
finite volume, constant Jacobian on a punctured source, or bounded trace.
The ordinary trace-pairing/full marked-energy implication and unrestricted
HC4/JC2 remain OPEN.
