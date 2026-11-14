# Determinant1020: the exact rational degree20 pencil

**Completed endpoint:** [full17 marked basis and explicit source inverse](DET1020_COMPLETE_MW17_EQUATION_2026-09-15.md).
This dependency note retains its earlier proof scope; its remaining-work
statements describe that stage, not the current construction status.

The complete two-dimensional Riemann–Roch space of the
[certified degree20 fiber](DET1020_PHYSICAL_MW17_FIBRATION_2026-09-15.md)
is now given by explicit rational functions. The
[coefficient certificate](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-rational-pencil-1024.json)
supplies integer polynomials N0,N1 with

    h0=N0/(x-a(t)),   h1=N1/(x-a(t)),   u=h1/h0=N1/N0,

on the source y²=x³+f(t)x+g(t), with P=(a(t),C(t)). A term is stored as
`[[t exponent, x exponent, y exponent], integer coefficient]`. The two
numerators have526 and524 terms and are linear in y. This is an exact
rational pencil, not a finite-field fit or finite-precision assertion.
The [subsequent cubic and Weierstrass conversion](DET1020_CUBIC_WEIERSTRASS_2026-09-15.md)
is now certified. The full seventeen section coordinate pairs remain uncomputed.

## Exact membership and completeness

In source divisor coordinates,

    D=40f_infinity+19O+P-sum(c_i Theta_i),
    c=(7,10;6,11;5,9,10;9,16,22,25,27,26,22,17,9).

The chains are I3 at lambda, I3 at1, I4, I10. Put
E=40f_infinity+19O+P. Each numerator has x-degree at most10 in its even
y part and at most9 in its y coefficient. Its horizontal pole at O is
therefore at most21; division by x-a lowers this by2, giving at most19.
The exact identity N_i(a,-C)=0 cancels the other possible pole at -P.
The remaining pole at P has order at most1.

Every monomial t^j x^k y^e in either numerator satisfies
j+4k+6e<=46. Since deg(a)=6 with leading coefficient1, x-a has vertical
pole order6 at infinity, giving pole order at most40 for h_i. There
are no other ordinary vertical poles.

The checker verifies all sixteen nonidentity component valuations over
Q, using split nodal charts. At an I_n fiber let r(s) be the critical
branch3r²+f=0, let d(s)=sqrt(3r(s)), and put H(s)=r³+fr+g. Then
ord(H)=n. On component j use

    U=s^j z, V=H/(s^j z), y=(U+V)/2,
    x=r+X, X*sqrt(3r+X)=(U-V)/2.

The inverse series has rational coefficients

    X=sum_(k>=1) binom(-k/2,k-1)/k * ((U-V)/2)^k/d^(3k-2).

The tangent sign at the I3 at1 is chosen so P meets component1, matching
the physical marking. The other chain orientations are recorded by the
same fixed tangent roots throughout. Each denominator valuation is
checked in the chart, and both cleared numerators vanish to at least
the required component coefficient plus that denominator valuation.
The largest required order is32; finite exact series suffice.

The local computation uses a faithful polynomial encoding. Terms s^a z^b
satisfy |b|<=a. For a required precision N, encode them by
q^((2N+1)a+b), and truncate at q^(2N²). Terms of s-order below N are
retained uniquely, while all terms of order at least N are removed.
The checker also verifies the source equation in each encoded chart.
This avoids arithmetic in a rational-function coefficient field.

Thus h0,h1 belong to L(D). Their numerator coefficient vectors are
independent over Q. Since the existing primitive nef certificate gives
h0(D)=2, these functions span the complete pencil. Its old-fiber degree
is therefore the certified intersection D.f=20. No modular rank bound
is needed for this final exact membership proof.

## Plane presentation and rational origin

Write N_i=A_i+yB_i. Eliminating y gives

    (u*A0-A1)²-(x³+fx+g)*(B1-u*B0)²=0,
    y=(u*A0-A1)/(B1-u*B0).

The checker divides out the fixed denominator component x-a and the
common vertical coefficient factor. The resulting exact
[plane presentation](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-rational-plane.json)
has x-degree20 and2958 terms in(t,x,u). It presents the generic genus-one
fiber birationally; it is not a Weierstrass equation.

The existing source component I4_theta1 has D-intersection1. Restricting
the exact pencil to its nodal chart gives a nonconstant fractional linear
function of z over Q. Its numerator and denominator coefficients are
in the [origin certificate](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-rational-origin.json).
This gives a convenient actual rational section as origin.

The same certificate translates the seventeen section classes
consistently. If v_i are the previous frame basis rows and Z is the new
component origin, put

    v_i'=v_i-(v_i.Z)D,      S_i'=Z+(G_ii/2)D+v_i'.

The rows D,Z+D,v_i' form an integral unimodular transport, preserving
the full height Gram G of determinant1020. Each S_i' has square-2 and
fiber degree1. Thus no saturation or marking information is lost by
this choice of origin. Their rational coordinate pairs remain uncomputed.

## Replay and discovery boundary

    sage -python research/elkies-k3/scripts/certify_det1020_rational_pencil.py --check

This replays exact rational pole/cancellation identities, all16 component
checks, independence, the plane identity, the rational component origin
and the translated integral height marking. K3 geometry and the prior
nef/divisor certificate remain named written inputs. Independent whole-
proof implementation, formal verification and external review are unclaimed.

Discovery used the459-dimensional ambient presentation in the retained
[compiler notes](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/rr-compiler-progress.md).
The67 replay had dimension2; the59 dimension11 result was retained as
a mismatch, without a rational conclusion. Unit-pivot lifting produced
compatible residue kernels. Reconstructing the auxiliary basis directly
was poorly conditioned: at1024 digits,344 entries failed. Normalizing
by the leading pole coefficients reduced this to17 failures. Reconstructing
the actual numerator polynomials instead recovered every coefficient at
the same precision. All final claims above depend on exact Q checks,
not on reconstruction heuristics. Higher-precision partial checkpoints
are preserved and are no longer needed for this pencil.

The [current conversion proof](DET1020_CUBIC_WEIERSTRASS_2026-09-15.md)
supplies the cubic, Weierstrass equation, their exact maps and four sections.
The remaining gates are the full seventeen basis coordinates and the
explicit inverse to the old source coordinates.
