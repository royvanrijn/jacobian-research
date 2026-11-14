# Determinant1020: explicit MW17 equation, marked basis and source inverse

The fixed source-degree20 fibration now has its rational pencil, polynomial
Weierstrass equation, all seventeen marked section coordinate pairs, and
explicit rational maps in both directions to the original source. The height
Gram is exactly the retained Gram, has determinant1020, and is saturated.
The inverse is given by fixed determinant ratios, rather than expanded
polynomials. Mathematical status is recorded in [MATH_STATUS.json](../MATH_STATUS.json).

## Use and exact outputs

The equation is `v^2=w^3+A(u)*w+B(u)`, with degrees8,12. Its parameter is
exactly the original certified `u=N1/N0`; the origin is source I4_theta1.
All coefficient arrays below use ascending powers, with exact rational strings.

- [Pencil numerator functions](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-rational-pencil-1024.json)
- [Weierstrass coefficients](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-weierstrass-invariants-QQ.json)
- [Seventeen polynomial section pairs and full intersection/height checks](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-seventeen-sections-QQ.json)
- [Exact identification with the particular physical basis](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-seventeen-marking-identification.json)
- [Integral trace-to-basis transport](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-trace-basis-transport.json)
- [Source inverse divisors and rational calibration](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-inverse-divisors.json)
- [Pinned inverse minors and exact nonvanishing witnesses](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-source-inverse-certificate.json)
- [Determinantal formula implementation](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/determinantal_source_inverse.py)
- [Sage loader and map evaluation API](scripts/load_det1020_weierstrass.py)

Run the loader with Sage, or import `load_model()` and `load_marked_basis()`.
The returned list has the exact order of the retained physical17 basis.
`load_known_sections()` preserves the earlier four-component interface.
`source_to_weierstrass(t,x,y,field=QQ)` and
`weierstrass_to_source(u,w,v,field=QQ)` evaluate the rational maps on their
stated affine domains; finite fields are also supported. Denominator zeros
raise arithmetic errors. Rational evaluation can be expensive. No coefficient
minimization or specialization-performance claim is made.

## Proof dependencies

The [physical fibration proof](DET1020_PHYSICAL_MW17_FIBRATION_2026-09-15.md)
supplies the rational integral marking, primitive nef rootless fiber D,
origin, Picard19 bound, torsion-free MW group and saturated determinant1020
height lattice. The [exact pencil](DET1020_EXACT_RATIONAL_PENCIL_2026-09-15.md)
supplies L(D), its source functions and plane identity. The
[cubic conversion](DET1020_CUBIC_WEIERSTRASS_2026-09-15.md) supplies the complete
H=D+O_source system, cubic embedding, polynomial Weierstrass model and exact
cubic/Weierstrass maps. Those claims retain their original scopes.

## Seventeen sections and their physical identification

Exact source resolution charts give parametrizations for the sixteen
nonidentity fiber components and the identity component of the old I4 fiber.
Their four numerator restrictions satisfy the quartic identity over Q. The
covering degrees agree with D.C, and the recorded integral MW coordinates
come from the inverse unimodular physical transport. The sixteen nonzero
component traces, together with the old zero trisection trace, form an
integral basis: their17x17 transition matrix V has determinant plus or minus1.
Its inverse specifies each particular retained generator.

For a component with restricted numerators n0,...,n3, intersect with the new
fiber in the finite algebra defined by `n1(z)-u*n0(z)=0`, after removing the
fixed common factor of n0,n1. Substitute `[1:n2/n0:n3/n0]` into the exact
Weierstrass map. The trace checker verifies every unit inverse, separability,
and that powers of the resulting w coordinate form a basis. Its monic
minimal polynomial U and the polynomial V representing v form a divisor
pair. Repeated exact reduction

    U <- monic((w^3+A*w+B-V^2)/U),  V <- -V mod U

computes its elliptic sum; zero remainder and strict degree decrease are
checked at each step. The old zero trace is the residual tangent point at
the cubic origin, since its degree3 hyperplane class is equivalent to
`2*origin + residual point`. Apply the integral inverse transition matrix.

The retained height4 classes imply polynomial coordinate bounds4,6 on this
minimal rootless K3 equation. Lifting and interpolation produced candidates;
all seventeen resulting identities `v_i^2=w_i^3+A*w_i+B` are then verified
exactly over Q[u]. The numerical lifts alone are not proof of a section.

For distinct polynomial sections the finite intersection is the degree of
the gcd of coordinate differences. At infinity it is
`min(4-deg(w_i-w_j),6-deg(v_i-v_j))`, with a zero difference assigned infinite
order. Polynomial coordinates mean each section is disjoint from the zero.
The rootless K3 height formula therefore gives diagonal4 and off-diagonal
`2-intersection`. The entire matrix equals the retained positive definite
Gram and has determinant1020. Since the full MW determinant is1020, the
index is1. The first point is exactly minus the old I10_theta2 component.

To identify all the particular basis labels, a separate checker recomputes
traces over GF(1013), without reading any p-adic checkpoint. It uses the
seven values0,1,2,4,6,9,10. Every component algebra is finite etale and every
necessary inverse is a unit. Trace specialization therefore commutes with
reduction at these values. All resulting marked points are affine.
Seven distinct residue values and the unit Vandermonde determinant imply
that the unknown characteristic-zero marked section polynomials are
1013-integral and that their reductions equal the candidate polynomials.

The coefficients A,B are1013-integral, and their reduced discriminant has
24 distinct roots and degree24. The global Weierstrass K3 thus has smooth
proper reduction, including the infinity chart. Neron-Severi specialization
tensored with Q is injective, and K3 NS is torsion-free with Pic^0=0; see
[van Luijk, Lemma2.1 and Proposition2.2](https://msp.org/ant/2007/1-1/ant-v1-n1-p01-s.pdf).
Equal reduced section maps give equal specialized divisor classes, hence
equal characteristic-zero classes. A section has self-intersection minus2,
so two distinct effective sections cannot represent that same class. This
identifies every coordinate pair with its actual retained physical section.

## Explicit inverse: old t by a39x39 determinant system

Write mu for the old I4 fiber location. On the generic new elliptic curve,
the full old I10 fiber is a degree20 divisor A0, whereas the old I4 fiber
is B0+Z, where B0 has degree19 and Z is the new origin. Consequently

    div(t/(t-mu)) = A0-B0-Z.

The inverse-divisor packet supplies every component parametrization: A0
uses all ten I10 components, B0 uses I4_theta0,theta2,theta3. For each
component, use its finite algebra and exact Weierstrass coordinate pair.
Let U19(w) be the product of the monic w-minimal polynomials of B0.

The nonvanishing witnesses check that these algebras are etale, the w
projections are primitive, their pole polynomials are pairwise coprime,
v is a unit on B0, and U19 is a unit on A0. These are generic properties
because all witnesses are reductions of the specified rational functions.
Thus B0 has no two points exchanged by elliptic negation and no2-torsion
point, and A0 avoids both B0 and its negative.

Multiplying `t/(t-mu)` by U19 cancels its finite poles. The result is
regular off Z and has pole order exactly39 at Z. Since the affine elliptic
coordinate ring is Q(u)[w,v]/(v^2-w^3-Aw-B), it has a unique form
`A19(w)+v*B18(w)`, with the top coefficient of B18 nonzero.
Normalize that coefficient to1.

The39 unknown coefficients, ordered A0,...,A19,B0,...,B18, satisfy39 linear
conditions: vanish at A0 and at the negative of B0. Rows are the power-basis
coordinates of

    [1,w,...,w^19, sign*v,sign*v*w,...,sign*v*w^18],

with sign minus on the finite-pole components and plus on the zero
components, in the packet order. This defines an explicit matrix M over
Q(u). Existence of the preceding rational function proves its rank at most38.
The retained finite-field nonzero minor proves rank at least38. Select the
pinned38 rows, let C be their first38 columns and b their final column.
The normalized coefficients are exactly

    c_j = -det(C with column j replaced by b)/det(C),  0<=j<38,
    c_38 = 1.

This is a finite rational formula; no nullspace choice or equation solving
is left unspecified. Linear solves in the evaluator implement these same
Cramer ratios. Set `R0=(A19(w)+v*B18(w))/U19(w)`.

The exact I4_theta1 resolution chart supplies
`Y_cubic=N3/N0=beta(u)*(t-mu)+O((t-mu)^2)`, with beta nonzero.
The pointed-quartic conversion gives
`v/w=6*l1/Y_cubic+O(1)` at Z. For example, using
`T=l1*X_cubic/Y_cubic+l2`, the cubic equation gives
`T=q*Y_cubic/l1^2+O(Y_cubic^2)`, while the quartic map gives
`v/w=6*q/(l1*T)+O(1)`. The normalized R0 has the same leading term as v/w.
Comparison with `t/(t-mu)=mu/(t-mu)+1` fixes the scalar:

    kappa = mu*beta/(6*l1),  r=kappa*R0,  old t=mu*r/(r-1).

Both beta and kappa are supplied as exact rational functions of u. Divisors
fix the function up to a Q(u) scalar; the exact leading term fixes that scalar.
This proves recovery of old t in characteristic zero.

## Explicit inverse: old x,y by Sylvester cofactors

First recover `[Zc:Xc:Yc]` using the certified Weierstrass-to-cubic formulas,
and set h=Xc/Zc. Write the retained source numerator functions as
`Ni=Ai(t,x)+y*Bi(t,x)`. Let F(t,xi,u) be the exact degree20 old-plane
polynomial and define

    G = ((A2-h*A0)*(B1-u*B0)-(B2-h*B0)*(A1-u*A0))/(xi-a(t)).

This division is exact: at xi=a(t), all Ai=C(t)*Bi, by the already checked
cancellation at the negative old section. The numerator has degree at most20,
so G has degree at most19. On the source, F and G both vanish at xi=old x
by their exact defining identities. The retained witness gives degrees20,19.

Form the39x39 matrix S with ascending xi columns0 through38 and rows
`xi^i*F` for i=0,...,18, followed by `xi^i*G` for i=0,...,19.
Its right kernel contains `(1,x,...,x^38)`, so its generic rank is at most38.
A retained nonzero38x38 minor proves its rank is exactly38. Let R consist
of the pinned38 rows. Alternating maximal minors give its kernel, hence

    old x = -det(R without column1)/det(R without column0),
    old y = (u*A0(t,x)-A1(t,x))/(B1(t,x)-u*B0(t,x)).

The second denominator is the one in the previously certified plane inverse;
it is also nonzero at the controls. This proves the entire inverse on a dense
open set of the source. The earlier cubic embedding is birational, so these
are inverse rational maps of the same generic elliptic curve and K3 surface.
No inference from three sampled compositions to a rational identity is used.

## Replay, audit and limitations

    sage -python research/elkies-k3/scripts/certify_det1020_complete_equation.py --check
    sage -python research/elkies-k3/scripts/load_det1020_weierstrass.py

The final checker replays the exact rational component restrictions, integral
transport, all17 characteristic-zero section identities and intersections,
good-reduction marking identification, inverse divisors and calibration,
and fixed nonzero inverse minors. The original-source points2P at t=2,3,4
provide three exact GF(1013) composition controls. These are fixed map tests,
not a specialization or point-search campaign. The inverse proof combines
these nonvanishing witnesses with the written divisor and kernel identities.
The prior pencil and cubic checkers remain separately replayable dependencies.

The requested fixed-target compilation is complete. Independent whole-proof
implementation, formal verification, external review, literature priority,
coefficient minimization and rank-jump experiments are not claimed.
