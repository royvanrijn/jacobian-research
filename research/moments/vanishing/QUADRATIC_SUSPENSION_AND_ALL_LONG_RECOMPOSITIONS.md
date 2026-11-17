# Quadratic suspension and all Long recompositions

16 November 2026. Internal written proof with retained resolution,
log-canonical-threshold and projective-period inputs. Independent replay,
formal verification and external review are not claimed. Both COMPLETE
exceptional MZ cases remain **OPEN**; review is deferred.

The entire nonlinear-composite gap for the opposite-weight Long
architecture is closed on both exceptional quotients, in every degree.
The proof constructs a separate three-dimensional compact-period model
with EXACTLY the same unweighted moments. Its original volume extends to
a projective model. The projective-period theorem then applies directly
to EVERY polynomial recomposition; no primitive-purity reflection is used.

For u,v of opposite nonzero torus weights and P=(1+u)(1+v+uv), the exact
classification is: H(P) is eventually pure iff H is the zero polynomial,
or uv=0 and H(1)=0. Every surviving phase is safe. Thus for nonzero
u and v the entire polynomial algebra C[P] has eventual radical {0}.

See the [checker](verify_haar_suspension.py),
[receipt](haar_suspension_verification.json) and
[source audit](haar_suspension_source_audit.json).

## 1. A suspension preserving every unweighted moment

Let X be a smooth complex surface, tau a nonconstant regular function,
Omega an algebraic 2-form and M a compact oriented real 2-cycle in X
with integral_M Omega=1. Form the affine-over-X threefold

    W: ab=(tau^2-1)/4,
    S=tau+a+b.                                         (1)

It has the rational 3-form

    Theta=Omega wedge da/a=-Omega wedge db/b
          =Omega wedge dS/(a-b).                       (2)

These formulas agree on overlaps: d(ab) is a multiple of d tau, which
wedges to zero with Omega. On a!=0 or b!=0 they are holomorphic, including
at tau=1 and tau=-1. In particular they define the SAME form across
the compact cycle below, not a different dualizing form on a singular
compactification.

Use z on Cstar to parametrize a compact cycle by

    (x,z) -> (x, a=(tau+1)z/2, b=(tau-1)/(2z)), |z|=1.

Its image is in W's smooth locus: a and b cannot both vanish. The
pullback of Theta is Omega wedge dz/z. Normalize by 1/(2pi i).
The resulting compact oriented singular 3-cycle C has mass1; possible
self-intersections of its image do not change its cycle or its integral.
Constant-term integration gives, for EVERY m>=0,

    integral_C S^m Theta/(2pi i)
      = integral_M L_m(tau) Omega,                     (3)

since

    CT_z [tau+(tau+1)z/2+(tau-1)/(2z)]^m=L_m(tau).

This is the same Legendre identity proved in the
[scalar-centre theorem](LEGENDRE_POISSON_JETS_AND_ALL_SCALAR_LONG_CENTRES.md).
Thus (3) holds for EVERY polynomial R(S), with no purity assumption,
Jacobian weight, marker insertion, branch norm or subsequence argument.

## 2. An explicit compactification criterion for the original 3-form

Take a smooth projective compactification Z of X on which tau extends
to a morphism Z->P1. Resolve it so the infinity fibre is SNC and the
divisor of Omega is SNC there. Assume:

* Omega is holomorphic over every FINITE height value.
* At each infinity component E, write b_E=-ord_E(tau)>0 and
  k_E=ord_E(Omega). Require k_E+ceil(b_E/2)>=0.          (4)

It suffices that these properties hold on a resolution; the application
below verifies the second for ALL divisorial valuations over infinity.

The function field of W is the quadratic extension of C(Z)(S) given by

    y^2=S^2-2S tau+1, y=a-b.                           (5)

Normalize Z times A1_S in this field. This is finite over that product,
hence projective over A1_S. Its projective resolution will give the
desired smooth model. We prove the actual form Theta is holomorphic
on such a resolution, including the exceptional divisors.

Over FINITE tau the model is locally (1). Resolve the nonzero function
f=(tau^2-1)/4 on the smooth base surface by point blowups. Omega remains
holomorphic. Locally f is a unit times r^b t^c; absorb the unit into a.
The equation ab=r^b t^c has a resolution preserving regularity of
Omega wedge da/a. Here is an explicit induction, rather than an appeal
to an unspecified singular dualizing generator. Blow up (a,r), b>=1.
In the chart a=r a', the strict transform is

    a'b=r^(b-1)t^c,
    dr wedge dt wedge da/a=dr wedge dt wedge da'/a'.

In the other chart r=a r', the strict transform solves for b and is
smooth; the same form becomes, up to sign, dr' wedge dt wedge da.
Decrease b to zero and then c. The ending equation ab=unit is smooth.
Regular coefficients of the pulled-back Omega remain regular throughout.
These proper blowups are isomorphisms on the original smooth locus.
Equivalently they show the local residue form has no negative
discrepancy; any further smooth resolution retains its holomorphicity.

Over tau=infinity set T=1/tau. Locally T is a unit times
product r_i^b_i (there are at most two r_i). With Y=T y, (5) becomes

    Y^2=T[(S^2+1)T-2S],
    Theta=T Omega wedge dS/Y.                          (6)

Divide out the square part of T and locally take a square root of its
unit. Writing epsilon_i=b_i mod2, the normalized equation has the form

    Y'^2=(product r_i^epsilon_i) F,
    F=(S^2+1)T-2S,
    Theta=unit*(product r_i^(k_i+ceil(b_i/2)))
                   dr_1 wedge dr_2 wedge dS/Y'.        (7)

A missing r_i means an ordinary second base coordinate. At a finite
S!=0, F is a unit on T=0. At S=0, F has derivative -2 in the S
direction, so F=0 is smooth and transverse to ALL r_i=0. The branch
divisor in (7) is therefore SNC near infinity, including the corner
T=S=0 where the two generic branch loci collide.

The normalized double cover Y'^2=x_1...x_l of a smooth 3-fold, l<=3,
has a resolution on which the residue volume dx_1 dx_2 dx_3/Y'
is holomorphic. Explicitly blow up intersections of two branch
components. In a chart x_2=x_1 x_2', set Y''=Y'/x_1; both numerator
Jacobian and denominator gain x_1, so the residue form is unchanged.
There are fewer branch components meeting above that chart. Continuing
resolves the SNC double cover into smooth or etale charts with no pole.
This is the usual elementary toric resolution, with its volume
calculation included. By (4) the prefactor in (7) is regular. Thus the
ORIGINAL Theta extends holomorphically at every finite S boundary.

The local resolutions prove extension on a global projective resolution:
take common resolutions with these local models, and use that a
meromorphic top form on a smooth variety is regular iff it has no
divisorial pole. Regular top forms remain regular under smooth blowups;
a common resolution also reflects a pole at a generic divisor. Choose
the global resolution to be an isomorphism over W's smooth locus, so
C and its period are retained. We obtain a smooth Y projective over
A1_S with the ORIGINAL Theta holomorphic and compact period1.
No assertion about every marker form extending is needed.

The retained
[projective original-volume period theorem](SEMISTABLE_LOG_PERIODS_AND_SINGLE_ISOTYPIC_SAFETY.md)
now gives, for EVERY nonconstant polynomial H,

    integral_C H(S)^m Theta/(2pi i)
        is NOT eventually zero.                       (8)

The morphism S is nonconstant: for a generic fixed x the conic in (1)
has a nonconstant a+b. Nonzero constant H is immediate from mass1.

## 3. The strict leading sphere bound gives BOTH parts of the criterion

Let X=Q2 be the affine sphere quadric, Z0 its smooth projective
completion and B its smooth infinity conic. Write a minimal-degree
nonconstant height as tau=H/u^D, where u=0 defines B and D>=3. Suppose
every root of H|B has multiplicity m with

    D>2m.                                              (9)

The original sphere volume has pole order2 along B. For ANY divisorial
valuation E above B put a=ord_E(u)>0, h=ord_E(H), and let A be its log
discrepancy over Z0. The exact volume order is

    k=A-1-2a.                                          (10)

At a FINITE height boundary, h>=Da. At a leading base point of
multiplicity m, the local total multiplicity of H is <=m, so the
retained elementary threshold bound gives A>=h/m. Consequently
A>2a, by (9), and k>=0. This is exactly the finite-volume argument
in the [torus volume theorem](VALUATIVE_BOUNDARY_VOLUMES_AND_TORUS_INVARIANT_CLOSURE.md).

At an INFINITE height divisor, b=Da-h>0. Criterion (4) is equivalent,
by integrality, to

    2(k+1)+b=2A+(D-4)a-h>0.                            (11)

If the centre is outside H|B then h=0 and A>=a, because (Z0,B) is
smooth with coefficient1 boundary. The expression in (11) is at least
(D-2)a>0. The same covers the generic point of B.

At a root with m>=2, use BOTH A>=a and A>=h/m. If h<=ma, the
expression is at least (D-m-2)a>0, since D>=2m+1. If h>=ma, it is at
least (D-4)a+(2/m-1)h. For m=2 this is (D-4)a>0. For m>2 the
coefficient 2/m-1 is negative, so h<Da gives

    (D-4)a+(2/m-1)h > (2D/m-4)a >0.

At a SIMPLE root, H is smooth transverse to B, so the SNC pair gives
A>=a+h. The expression is then at least (D-2)a+h>0. Thus (11)
holds for EVERY infinity valuation, including all exceptional divisors
of any subsequent resolution. This proves (4).

The only classical valuative inputs here are the retained
log-discrepancy definition and lct(H)>=1/mult(H), together with the
elementary SNC inequalities A>=a and A>=a+h. Inversion of adjunction
or a classification of canonical threefold singularities is not needed
for this construction. Embedded resolution and normalization finiteness
in characteristic zero are standard retained geometric inputs.

## 4. Exceptional invariance makes the criterion automatic

For every nonconstant tau in C[Q2]^O or C[Q2]^I, the degree-2D leading
zero divisor on B=P1 is invariant. The retained minimum orbit sizes
6 and12 give D>=3m and D>=6m respectively at EVERY root. Since this
positive-degree section has a root, D>=3. Hence (9) and the degree
condition of Section3 hold for EVERY such height, with all complex
coefficients and lower terms. No degree or support search is used.

Now let A_Gamma=C[SL2]^Gamma for Gamma=2O,2I and u,v have opposite
nonzero left-torus weights. Put

    P=(1+u)(1+v+uv), tau=1+2uv.

If u and v are both nonzero, the retained domain/constant-unit argument
makes tau nonconstant. Its neutral moment functional is exactly the
normalized sphere period on Q2. By the universal Legendre identity,

    h_Gamma(R(P))=integral_C R(S) Theta/(2pi i)
        for EVERY polynomial R.                       (12)

Sections2--3 give the proper ORIGINAL-volume model of this SEPARATE
suspension. Apply (8) with R=H^m. Every nonzero polynomial H has
H(P) NOT EVEN EVENTUALLY PURE. This directly handles every H; it does
not infer purity of any primitive phase from purity of a composite.
Nor does it assume a proper model or an unweighted-Haar-preserving map
for the full original mixed-weight phase P itself.

If uv=0, one of u,v is zero and P=1+w with w in a single nonzero
weight, possibly w=0. Then h_Gamma(H(P)^m)=H(1)^m. Thus eventual
purity is equivalent to H(1)=0. In that case H(1+w) has only strictly
one-sided weights, or is zero, and is safe against ALL original markers.
We have the exact ALL-RECOMPOSITION classification

    H(P) eventually pure iff
       H=0, or [uv=0 AND H(1)=0].                      (13)

Every eventually pure phase in this entire recomposition class is SAFE.
This holds in all degrees on BOTH quotients and under retained Haar
transport by algebraic automorphisms and conjugate tori. It is a
complex-field theorem; no arbitrary-ring promotion is claimed.

## 5. The entire equivariant Bezout Long construction is closed

For EVERY regular Bezout b of weight L=n+m, (n,m)=(8,12) or(12,20),
the original functions u=Y_b=V+G0 b and v=F0 have opposite nonzero
weights and are both nonzero. The determinant identity excludes Y_b=0.
Consequently EVERY nonzero polynomial H satisfies

    H(g_b), g_b=(1+Y_b)(1+F0+Y_bF0),
        NOT EVEN EVENTUALLY PURE, in EVERY degree.     (14)

This closes ALL nonlinear recompositions of the entire torus-equivariant
primitive Long branch. The older scalar-centre theorem and its Poisson
jets remain valid, now subsumed for this application. The prior
all-b nonzero-component and uniform-marker statements remain: the
[pullback theorem](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md)
does not by itself assert purity, and (14) supplies a genuine purity
exclusion on this weight branch. Other homogeneous b weights and mixed b
are NOT excluded. The distinct Hopf phase's general nonlinear composites
are NOT settled by this quadratic suspension.

For a boundary control, the noninvariant linear sphere height has
D=1, k=-2 and b=1 at generic infinity. Then k+ceil(b/2)=-1 and the
volume condition fails. This allows the known Long pure example and
its pure composites; one must not omit the infinity estimate or replace
the original volume by a generator on a singular graph. The new proof
also expressly retains known exceptional one-weight SAFE phases through
the degenerate locus in (13).

The WHOLE exceptional target is still purity forcing all corrected
actual local torus orders to vanish, or an actual invariant pure failure
outside the closed constructions. Both COMPLETE exceptional MZ cases
remain OPEN. No generic proper-model assumption, marker-to-unweighted
transfer, nonlinear purity reflection, compute campaign or review gate
has been introduced.
