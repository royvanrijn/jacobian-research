# Every nontrivial constant-endpoint profile has a primitive nonpower lift

Internal proof addendum, 20 November 2026. Fix a characteristic-zero field
K and an integer a>=0. Consider a NONTRIVIAL normalized constant-endpoint
profile p=t^a(C+t²D)/x. The retained marked Davenport--Zannier classification
gives

    C=H, D=-HP², P(0)=H(0)=1,
    P squarefree, r=degP>=1, degH=ar+1,

and every geometric H-root is a P-root. Purity and the fixed original
marker Z detecting infinitely many EVEN powers are already proved there.

**Theorem ConstantEndpointPrimitiveNonpowerLifts.** The minimal balanced
original lift F of every such profile is radial-primitive and is not a
nonzero scalar times a proper polynomial power over ANY field extension
of K. Its balanced degree is d=(a+2)(r+1). Thus the retained geometric
counts enumerate primitive NONPOWER Image-counterexample families in
this exact profile class and its stated x-scaling equivalence. They do
not count or classify arbitrary Image polynomials modulo arbitrary
changes of variables.

## 1. The literal minimal lift and its radial residue

Use the original bilinears

    R=xi1*z1+xi2*z2, Z=xi1*z2,
    W=2xi2*z1, T=xi1*z1-xi2*z2,
    T²+2ZW=R².

Put h=ar+1 and form the homogeneous polynomials

    P_h=R^r P(Z/R), H_h=R^h H(Z/R),
    B=2W R^(2r)-T²*(P_h²-R^(2r))/Z,
    F=T^a H_h B.                                     (N1)

The division by Z is literal polynomial division, since P_h(R,0)=R^r.
The degree of B in bilinears is2r+1. Therefore F has balanced degree

    a+h+2r+1=(a+2)(r+1)=d.

At R1 its angular polynomial is exactly

    t^a H[2y-t²(P²-1)/x]=t^a H(1-t²P²)/x=p.

The retained exact least-degree formula makes d minimal, rather than
merely a degree of a possible lift. If p_r=lcP and h_h=lcH, then

    F modR = -h_h*p_r²*T^(a+2)*Z^(h+2r-1) !=0.       (N2)

Here R is a prime bilinear and neither T nor Z vanishes in its integral
quotient. The two leading coefficients are nonzero. Thus R does not
divide F, over K or any extension. This proves radial primitivity.

## 2. Geometric factor orders exclude every proper power

Work temporarily over an algebraic closure k of an arbitrary extension
field of K. The geometric roots alpha of P are nonzero and distinct.
For each define

    L_alpha=R-Z/alpha.

Because P(0)=H(0)=1 and every H-root belongs to P,

    P_h=product_alpha L_alpha,
    H_h=product_alpha L_alpha^(h_alpha),
    sum_alpha h_alpha=ar+1.                          (N3)

Zero multiplicities h_alpha are allowed. Every L_alpha is a prime
bilinear over k: its coefficient matrix is upper triangular with diagonal
entries1,1, hence has rank2. A bihomogeneous polynomial of bidegree(1,1)
can factor nontrivially only into a xi-linear and a z-linear factor,
which would give matrix rank1. This also proves absolute irreducibility
of R and T. The L_alpha are mutually nonassociate and none is associate
to T or R; they also do not divide Z. Characteristic zero ensures that
T, with diagonal entries1,-1, is distinct from every L_alpha.

Modulo T the polynomial B is2W R^(2r), nonzero in that prime quotient.
No L_alpha divides T. Consequently

    ord_T(F)=a.                                     (N4)

The quadric identity gives the exact polynomial equation

    ZB=R^(2r+2)-T² P_h².

Modulo each L_alpha its right side is R^(2r+2), nonzero in the prime
quotient. Hence L_alpha does not divide B, and

    ord_(L_alpha)(F)=h_alpha.                        (N5)

Suppose F=cG^k for a nonzero scalar c and integer k>=2, over any extension
field. Extend that field to an algebraic closure. Irreducible factor
orders in the polynomial UFD are multiples of k, so(N4)--(N5) imply

    k|a, k|h_alpha for all alpha,
    hence k|(ar+1) and k|a.

Their gcd is1, a contradiction. When a0, k|a supplies no information,
but sum h_alpha=1 already forces k1. Thus the argument includes a0
without invoking a nonexistent positive T-order. The contradiction over
the larger algebraic closure excludes proper powers over the original
extension field as well. QED.

## 3. Exact primitive-nonpower family count and field scope

**Corollary ConstantEndpointPrimitiveNonpowerFamilyCount.** Fix a>=1,r>=1.
Over an algebraically closed characteristic-zero field, the normalized
constant-endpoint profiles with degP=r give exactly

    (1/r)*binomial(ar,r-1)

families of radial-primitive NONPOWER Image counterexamples at balanced
degree(a+2)(r+1), under their specified SOURCE-scaling equivalence.
For a0 there is one such class at r1 and none at other r.

This is the retained CE-TreeCount applied to the theorem just proved.
For precision, x->s*x on normalized coefficients is C_s(x)=C(s*x),
D_s(x)=D(s*x). On original variables let phi_s scale xi1 bys andz1
bys^(-1), fixing xi2,z2. Then R,T are fixed, Z maps tosZ and W toW/s;
the literal minimal lifts obey

    F_s=s*phi_s(F).

Dehomogenization recovers the profile, so this specified equivalence
introduces no additional identifications of profile classes. Arbitrary
linear or polynomial changes of the original variables are NOT quotiented.
The number is geometric: over a prescribed nonclosed K, one retains only
the marked classes fixed by Gal(Kbar/K). The retained canonicalq=1
section provides a K-model exactly for those classes, and each K-model
family is its free K* source-scaling orbit. No assertion that all the
geometric families have rational or real coefficients is made.

## 4. Exact controls and existing antecedents

The standalone checker reconstructs literal ORIGINAL four-variable
polynomials for three normalized profiles:

* a0,r1: P=1+x,H=1+x, balanced degree4;
* a1,r2: P=(1+x)(2+x)/2,H=(1+x)^3, balanced degree9;
* a2,r2: P=(1+x)(1+2x/3),H=(1+x)^5, balanced degree12.

It verifies the bilinear identity, the division in(N1), all balanced
monomial degrees, radial nondivisibility and a nonzero R0 specialization,
and the ACTUAL T/L_alpha factor orders by exact polynomial division.
It checks original pure and Z-marked contractions at m1,2 against direct
angular expansion. Squaring a fixture doubles the factor-order gcd,
providing a negative control for the nonpower criterion. These finite
controls are separate from the unrestricted proof(N2)--(N5).

All-order purity, infinitely many fixed-Z survivals and geometric
classification/counts are retained antecedents. No new witness degree,
minimum pair dimension or historical-priority claim is made here.

The [checker](verify_constant_endpoint_primitive_nonpower_lifts.py),
[inputs](constant_endpoint_primitive_nonpower_lifts_inputs.json),
[source ledger](constant_endpoint_primitive_nonpower_lifts_source_ledger.json)
and [receipt](constant_endpoint_primitive_nonpower_lifts_verification.json)
retain the literal original-variable controls.
