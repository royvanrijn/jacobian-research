# The elliptic contact34 family and explicit key models of the first image

14 November 2026. PF2R7ELL1 identifies every polynomial degree(3,5)
image of infinity contact at least34 with a coefficient point of an
explicit cubic E_C. Its smooth projective completion has genus1 except
at39 normalized first-image coefficient points, where it is nodal.
The seven higher-contact choices are the zero fibre of a degree7
function on E_C. Its critical coefficient scheme is finite free of
rank8, and the universal evaluation Jacobian is exactly the derivative
of that degree7 function. The critical resultant of the quintic is
minus15 times the same derivative.

Exact contact34 witnesses give genuine nodal image pairs with sixteen
transverse affine intersections and mutual singular-locus avoidance.
Thus Nori's equality case is populated at the curve level; these local
tests do not close it. A rational one-parameter ansatz cannot cover
the smooth elliptic coefficient curve for a fixed C.

PF2R7KEY1 gives the exact implicit equation of C in nine key positions,
and a finite polynomial target map(U,V)->(U,G_C). In L this restricts
to a degree3 map onto a rational cubic. In K, away from the common
class-I points, it is finite birational onto a degree9 curve, with an
explicit length6 partial-normalization defect supported over one U value.

These are authored curve/function-field and coefficient-family proofs
with exact identities, monic Groebner and good-prime rank checks, without
formalization or independent external review. The polynomial target
map is not an automorphism or an intermediate field of the original
degree6 Keller extension. No global Keller covering or JC2 closure follows.

## 1. All quintics of contact at least34

Use the complete first-image charts from PF2R7FAM1, the closest quintic
G_C and its degree9 restriction g from PF2R7B5CONTACT1. Retain all
normalization and target freedoms as in PF2R7B5REC1:

    P=s^3+A*s+M,
    Q=s^5+b*s^4+c*s^3+e*s^2+f*s+n.

Five constant-linear equations reconstruct b,c,e,f,n polynomially from
A,M and the first-image parameter. Their actual implicit quintic F_B
satisfies the exact identity already reconstructed there,

    G_C-F_B=r9*U^3+r8*U*V+k6*U^2+k5*V+k3*U+k0.      (1)

All higher implicit coefficients have thereby been matched. The first
image pole degrees of the displayed monomials are18,16,12,10,6,0.
By the complete thirteen-element restriction basis, contact at least34
is equivalent to r9=0. Put f0=(27/5)*r9. Contact is exactly34 when
r8 is nonzero, and at least38 when r8=0. Thus the full coefficient
curve E_C is f0=0, not just an arbitrary subfamily. Distinct coefficient
points give distinct B images by the unique centered monic normalization.

## 2. The explicit coefficient cubic and its degeneration set

In either chart f0 has the form

    f0=A^3+B1*A*M+P1*A-9*M^2+B0*M+R0,

where the coefficient polynomials depend only on d or z. Set

    Y=18*M-(B1*A+B0).

Then f0=0 is isomorphic to the Weierstrass cubic

    Y^2=(B1*A+B0)^2+36*(A^3+P1*A+R0).             (2)

The leading cubic coefficient is36. The projective point at infinity
is unique and smooth. Its parameter poles are ord_infinity(A)=-2
and ord_infinity(M)=-3. The complete cubics are in models.json.

Their discriminants are

    Disc_L=-(3048192/78125)*F4(d^3),
    Disc_K=-(62208/78125)*z^3*H9(z),

where

    F4(w)=53248w^4+142848w^3+144720w^2+65340w+10125,

and

    H9(z)=390625z^9+937500z^8+1031250z^7+1062500z^6
          +853125z^5+402000z^4+94800z^3+7680z^2-4800z+1024.

Both primitive factors are squarefree and nonzero at0. Exact gcds with
the Weierstrass c4 show that each vanishing fibre is an irreducible
nodal cubic, not a cusp. Thus E_C is smooth of genus1 except at12 L
parameters and27 K x parameters. The shared L/K class-I coefficient
points, with d^3=x^3=-2/3, meet neither degeneration set, so there are
39 distinct exceptional normalized first-image coefficient points.
These are degenerations of the auxiliary coefficient curve E_C; they
are not a new count of singular first-image curves C.

For a fixed smooth E_C, any rational one-parameter family of these
quintics defines a rational map P1->E_C, extending to its smooth
projective completion. Riemann--Hurwitz excludes a nonconstant such
map. This does not exclude individual coefficient points or families
in which C itself varies. At the39 nodal auxiliary fibres a rational
parametrization can become possible.

## 3. Seven high-contact points and eight critical coefficients

The function r8 on E_C has pole order7 at its single infinity point:
its unique highest weighted term is(5/9)A^2*M. Therefore

    r8:E_C -> A1

is finite of degree7. Its zero scheme is exactly the finite-free
seven-choice interface of PF2R7B5REC1. All seven points there count
coefficient multiplicity, with generic reducedness already witnessed.

Write D=f0_M*partial_A-f0_A*partial_M. For a smooth E_C this is dual
to its everywhere nonzero differential dA/f0_M, using the alternate
expression where f0_M=0. It extends over infinity. The exact identities
are

    J:=P_s*D(Q)-Q_s*D(P)=D(r8) modulo f0,
    Res_s(P_s,Q_s)=-15*D(r8) modulo f0.             (3)

Thus the evaluation map E_C times A1_s -> A2_(U,V), given by(P,Q),
is ramified on exactly the vertical fibres corresponding to nonimmersive
quintics. It is not a Keller map. When E_C is smooth, D(r8) has a pole
of order8 at infinity, hence eight finite zeros counted with multiplicity.
Equivalently the degree7 map on the projective elliptic curve has
total ramification14, of which6 is at infinity and8 is finite.

There is a stronger coefficient-scheme statement, valid even at all
39 nodal auxiliary fibres. Order A,M by weights2,3, with M preferred
at equal weight. A monic Groebner basis of the ideal(f0,Res(P_s,Q_s))
has leading monomials M^2,A^4, with constant leading units. The second
generator is the first remainder of the resultant, multiplied by3/175.
The two leading monomials are coprime and their S-pair reduces to zero.
The eight standard monomials A^i*M^j,0<=i<=3,0<=j<=1 give a finite
free rank8 scheme over each entire first-image parameter base.

This counts nonimmersive coefficient points with multiplicity. At a
nodal E_C, D vanishes at the node and is not a local tangent-frame;
the smooth-surface ramification interpretation is not asserted there.
The rank8 coefficient statement remains valid without that interpretation.

## 4. The Nori equality case has locally admissible images

At d1 and separately z1, set A0 in f0. The two monic coefficient
equations are

    L: M^2-(34/135)M-80/729=0,
    K: M^2-(28861/1440)M+99270355/995328=0.

Both are reduced. For every root in each pair, the exact matrices
modulo101 certify: r8 nonzero; four simple off-diagonal collision pairs;
no triple value of B; B avoids Sing(C); C avoids Sing(B); and all16
affine C/B intersection roots are simple. The finite-free dimensions
are2,8,48 and32. All relations are monic in successive variables, and
the leading intersection coefficient is first proved a unit, exactly
as in PF2R7B5GEN1. Hence the ranks lift to characteristic zero and
test both coefficient roots simultaneously, without numerical embeddings.

Each B has four ordinary nodes, C has twelve nodes, and there are
sixteen distinct transverse intersections, all smooth on both curves.
Thus I_infinity=34 and the prepared B self-intersection is8, exactly
Nori's strict-positivity boundary. Openness gives a nonempty open of
the elliptic coefficient curve of locally compatible equality34 images
at these first-image parameters, and on nonempty first-image base opens.
This supplies image pairs, not a connected S6 covering or source A2.
If a row C or shared-image row A were realized with such a pair, its
sixteen intersections would be singleton values by the existing ledger.
These statements prove PF2R7ELL1.

## 5. Exact implicit equations in nine positions

The pole basis of PF2OMIT16 shows that g^2, of pole18, has an exact
relation with p^3,pq,pg,p^2,q,g,p,1. Direct cancellation in the original
normalizations gives

    R_C(U,V,W)=W^2+ell(U)*W+H(U)+L(U)*V,            (4)

where deg ell<=1,deg H<=3,deg L<=1, and

    R_C(p,q,g)=0, g=G_C(p,q).

The displayed nine coefficient positions and their complete rational
functions are retained in models.json. Substitution W=G_C makes a
monic V-degree6 polynomial of total degree10. Since the image C is
irreducible of degree10 and the equation vanishes identically, this
is its exact reduced implicit equation, not a multiple or partial
initial form. The target polynomial map(U,V)->(U,G_C) is finite of
degree3 because G_C is monic cubic in V.

## 6. The polynomial component compresses to a rational cubic

For L, L(U)=0 identically. With u=t^3+d*t the literal identities are

    p=u^2+u,
    g=-u^3/6561+(73/531441+49d^3/1417176)*u.

The image in(U,W) is a rational cubic, with normalization u, since
the degrees2,3 force birationality. The original C normalization maps
to it by the degree3 polynomial u(t). Completing the square in(4) gives

    (W+ell(U)/2)^2
      =(4U+1)*(-648U+147d^3-64)^2/72301961339136.

It is nodal for3d^3+2 nonzero and cuspidal for3d^3+2=0.
The distinction concerns this cubic image, not all singularities of C.

## 7. The Laurent component is a partial normalization of a degree9 curve

For K with z nonzero,

    coefficient_U(L)=z^12*(3z+2)^2/129140163.

At z=-2/3 the constant coefficient vanishes too, and the restriction
has degree3 as in the shared cubic case. Otherwise L(U)=aU+b has a
nonzero a. Equation(4) gives

    L(p)*q=-g^2-ell(p)*g-H(p),

so C(p,g)=C(p,q)=C(t). The restriction is finite birational. The image
C' has polynomial normalization degrees(6,9) and projective degree9.

Put A=C[p,q],A'=C[p,g]. PF2OMIT16 gives the free C[p]-basis of A,
reordered as1,g,q,qg,q^2,q^2g. Then A' has the explicit basis

    1,g,Lq,Lqg,L^2q^2,L^2q^2g.                    (5)

Each element in(5) is polynomial in p,g by(4). Conversely their span
contains1 and is closed under g, using(4) and the monic cubic equation
g=G_C(p,q). The checker replays these six closure identities with
independent formal coefficients. Thus the span is precisely A'.
Consequently, as C[p]-modules,

    A/A'=(C[p]/L)^2 direct_sum(C[p]/L^2)^2.         (6)

Its vector-space dimension is6, all supported over p=-b/a. This is a
normalization defect length, not a number of distinct singular points.
The pole degrees of(5) are0,9,16,25,32,41, with distinct residues mod6.
The affine pole semigroup is<6,9,16>, with conductor36 and18 gaps.
Thus C' has affine normalization defect18, while C has12.

Finally g^2+c3*p^3 has pole16, with nonzero coefficient from L(p)q;
the other terms of(4) have smaller pole. In infinity coordinates
X=p/g,Y=1/g, subtracting the leading analytic cubic in X leaves order11,
while X has order3. Hence C' has infinity characteristic(3;11) and
delta10, consistent with arithmetic genus28 minus affine delta18.
At the common z=-2/3 points, L vanishes identically and the field
degree is3 rather than1, as the quadratic relation and pole9 exclude
g being a function of p alone. These statements prove PF2R7KEY1.

None of these finite target transformations preserves the Keller
Jacobian condition automatically. Neither a curve-field compression nor
the elliptic parameter family supplies an intermediate field in the
primitive degree6 source extension. Those missing global bridges,
other image-degree charts and the complete covering/source construction
remain part of the original research goal.
