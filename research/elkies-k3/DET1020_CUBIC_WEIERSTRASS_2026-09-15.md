# Determinant1020: exact cubic, MW17 Weierstrass model and four sections

**Completed endpoint:** [full17 marked basis and explicit source inverse](DET1020_COMPLETE_MW17_EQUATION_2026-09-15.md).
This dependency note retains its earlier proof scope; its remaining-work
statements describe that stage, not the current construction status.

The fixed degree20 pencil now has a polynomial Weierstrass equation over
Q(u), with coefficient degrees8 and12 and nonzero discriminant of degree24.
Its parameter remains exactly the certified u=N1/N0. Exact birational maps
between the new cubic and Weierstrass models are supplied. Four independent
sections have explicit polynomial coordinates of degrees4 and6.

The full seventeen section coordinate pairs and an explicit inverse to the
original source coordinates t,x,y remain open. The rank17 and saturated
seventeen-class marking are inherited from the existing fixed-fibration
proof; four computed points alone do not prove rank17.

## Reproducible outputs

- [Exact four-dimensional linear system](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-rational-trisection-1024.json)
- [Quartic and cubic, retaining the original pencil](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-trisection-quartic-QQ.json)
- [Polynomial Weierstrass coefficients](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-weierstrass-invariants-QQ.json)
- [Exact compositional forward and inverse cubic maps](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-weierstrass-compositional-QQ.json)
- [Four component sections and integral marking coordinates](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-component-sections-QQ.json)

Coefficient lists are in ascending powers of u. The equation is

    v^2=w^3+A(u)*w+B(u),  A=-27I, B=-27J.

The coefficients are large; no constant-height minimization is claimed.
The [Sage loader](scripts/load_det1020_weierstrass.py) constructs the curve
and its four known points without a search or rank computation.

## A small extension of the completed pencil calculation

Let O be the old source zero section and set H=D+O. The physical marking
gives D.O=3, O^2=-2, H^2=4 and H.O=1. Since D is nef and every other
irreducible curve meets O nonnegatively, H is nef and big. Characteristic-
zero K3 vanishing gives h0(H)=4.

Its source ambient divisor is E=40f_infinity+20O+P. This changes the
previous459-function presentation to461 functions: add x^10 and the next
normalized phi term, now with phi indices0 through58. The normalized phi
expansion includes even x-powers through10 and y*x-powers through8.
The same sixteen component conditions apply.

The four reconstructed numerator polynomials are checked exactly over Q.
Each is linear in y, has weighted degree at most46, even x-degree at most11
and y-coefficient x-degree at most9, and satisfies N(a,-C)=0. Division by
x-a gives poles bounded by E. Exact component tests give membership in
L(H); coefficient-vector independence and h0(H)=4 give completeness.
The numerical lift is discovery only, not a proof input replacing these
rational checks.

On a generic new fiber, H has degree3. The restriction sequence has kernel
H0(O), of dimension1, and H1(O)=0. Thus the four global sections restrict
to the complete three-dimensional degree3 system, embedding the smooth
genus-one fiber as a plane cubic. This does not change the elliptic
fibration or seek a different divisor D.

Keep the original N0,N1 and select two complementary numerator functions.
Their restrictions to the rational component Z=I4_theta1 lie in the span
of the restrictions of N0,N1. Subtracting those exact constant linear
combinations centers the extra functions N2,N3 at Z. The cubic origin is
then [1:0:0] in coordinates [Z:X:Y]=[N0:N2:N3].

The35 quartic monomials in the four numerator functions have exact rank34
in the source coordinate ring. Their unique relation vanishes identically
after reducing y^2=x^3+fx+g. It contains the line N0=N1=0; substituting
N1=u*N0 and dividing by N0 gives the cubic. The original-source forward
map is therefore explicit. Birationality follows from the complete
degree3 system on the generic fiber, independently of an uncomputed
formula for the inverse old coordinates.

## Exact Weierstrass conversion without large intermediate expansions

Write the centered cubic as

    F=Z^2 L(X,Y)+Z Q(X,Y)+C(X,Y),
    W=2LZ+Q,  W^2-(Q^2-4LC)=4LF.

For a binary quartic aT^4+bT^3+cT^2+dT+e, the invariants used are

    I=12ae-3bd+c^2,
    J=72ace+9bcd-27ad^2-27b^2e-2c^3.

The checker also proves the conversion directly. For the pointed quartic

    V^2=aT^4+bT^3+cT^2+dT+q^2

define

    x0=(2q(V+q)+dT)/T^2,
    y0=(4q^2(V+q)+d(V+3q)T+2qcT^2+qbT^3)/T^3.

Exact universal polynomial identities give

    y0^2=x0^3+c*x0^2+(bd-4q^2a)*x0+q^2b^2+ad^2-4q^2ac,
    T=(2q*y0+d*x0+2q^2*b)/(x0^2-4q^2*a),
    V=(x0*T^2-d*T-2q^2)/(2q).

The checker verifies both curve substitutions and the inverse recovering
T,V before instantiating these identities with the exact Q[u] coefficients.
It checks the linear binary-coordinate change and its invariant scaling,
then the shift and scaling to A=-27I, B=-27J. This is exact compositional
verification, not numerical sampling or finite-order expansion.

The sign q=-Qnew(0) is the branch of the specified cubic origin: locally
F=0 gives W=-Q-2C/Z at [1:0:0]. The maps therefore carry the retained
I4_theta1 origin to Weierstrass infinity. Their denominators are nonzero
as functions; the formulas describe birational maps on dense open sets.

## Four explicit sections and the remaining seventeen-basis endpoint

The source components I10_theta2, theta4, theta7 and theta9 have new
fiber degree1. Their exact nodal restrictions give cubic coordinates
[1:alpha+beta*u:gamma+delta*u]. Substituting the certified map cancels
the denominators exactly and produces x-polynomials of degree4 and
y-polynomials of degree6. Every Weierstrass identity is checked over Q[u].

Each point is identified with its physical component class. The output
records its integral coordinates in the retained seventeen-dimensional
MW basis. These four vectors are independent and each has height4;
their full mutual height Gram is computed from the fixed marking.
They are not a replacement for the requested full basis.
The first retained basis point is exactly minus I10_theta2, as its
recorded marking vector is (-1,0,...,0).

The next stage is to recover the remaining marked generators, for example
by transporting divisor traces of the known source components and section,
and to construct the inverse rational map to old t,x,y. All seventeen
retained basis heights equal4; on this minimal K3 model their section
coordinates have the corresponding polynomial degree bounds4 and6.

## Replay and proof boundary

    sage -python research/elkies-k3/scripts/certify_det1020_cubic_weierstrass.py --check

This replays the exact four-function membership checks, the quartic/cubic
identity, polynomial invariants, compositional birational maps and four
section identities. Prior K3, nefness, rational marking, Picard and
saturation results remain named inputs. Independent whole-proof
implementation, formal verification, external review and novelty are
unclaimed. The fixed-target construction goal remains open.
