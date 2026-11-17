# Complete complex fibres and nonproper locus of Gao F6

14 November 2026. Claims **GAOF6FIB1**, **GAOF6NP1**. All fibres and
strata are over C, and strata carry their reduced structures.

## Prior results and scope

Use precisely F6 in [Gao, arXiv:2608.00222v2](https://arxiv.org/html/2608.00222v2),
Sections 4.5.1 and 9. Gao already proves the simple-root correspondence,
the omitted surface, the x=0 boundary isomorphism, and the explicit
identity factor. Those results are reproduced or imported explicitly,
not claimed as additions. The same revision leaves full F4--F7 fibre
classifications open in Sections 11.1 and A.3.

The [Korea Superintelligence Labs exact-fibres note](https://ideosphere.ai/papers/keller-fibers)
already gives an exact rational six-point F6 fibre and real-root constraints.
Its selected fibres are credited prior results. The inspected manuscript
does not contain the complete boundary counts or all multiplicity strata
below. No exhaustive novelty priority is claimed.

The additions here are the complete complex cardinality stratification,
explicit rational descriptions of all curve strata and exceptional points,
and the entire nonproper-value locus. This does not classify F5 or F7,
and gives no proof of the planar Jacobian conjecture or unrestricted HC4.

## Definition, polynomiality and etaleness

In source coordinates (x,y,z1,z2,z3), set

    v=(xy,x^2*z1,x^3*z2,x^4*z3),
    gamma=1-29v1+999v1^2+355v1*v2-41553v1^3+v4,
    u1=1+27v1-5v2,
    u2=v1-12v2+(2128/5)v1^2+v3,
    u3=-4v1-10v2,
    s=gamma*u1, a=gamma^3*u2, b=gamma^4*u3, C=gamma*x.

Put

    G2=a+s*b+s^2-s^3, G3=2s*G2+a,
    psi=2s^4-2s^5, G4=b-G2^2+psi,
    X1=-det d(G2,G3,G4)/d(s,a,b)
      =2a+2bs-b-10s^5+8s^4-2s^3+5s^2-2s,
    X=(X1,G2+sX1,G3+s^2X1,G4+aX1),
    S=X+gamma*(1,s,s^2,a).

Then F6=(C,S1/C,S2/C^2,S3/C^3,S4/C^4), with polynomial extensions.
One can certify those extensions cheaply without expanding every high
degree component. Each

    E_i(gamma,u)=S_i(gamma,gamma*u1,gamma^3*u2,gamma^4*u3)/gamma^i

is polynomial. After the polynomial stage, it is polynomial in v1..v4.
Assign weights 1,2,3,4 to these variables. A monomial of weight j becomes
x^j times a source monomial, so output i is polynomial exactly when all
terms of weight below i vanish. The checker substitutes
v=(e*y,e^2*z1,e^3*z2,e^4*z3) and verifies these finitely many coefficients
through weight four. Since all weights are positive, no unexpanded
higher-weight term can affect that test. The same exact jets give the
complete boundary restriction below.

For etaleness, the v-stage has determinant -290, and
det dS/d(gamma,s,a,b)=gamma, both checked symbolically. The other
determinant factors are x^10, gamma^8 from scaling the three sweep
coordinates, gamma from padding with C=x*gamma, and C^-10 from the
target twist. Their product is -290 on a dense open and hence globally.
Thus all fibres are reduced and zero-dimensional.

The stage has a polynomial inverse in (gamma,u1,u2,u3):

    v1=(u1-1-u3/2)/29, v2=-u3/10-2v1/5,
    v3=u2-v1+12v2-(2128/5)v1^2,
    v4=gamma-1+29v1-999v1^2-355v1*v2+41553v1^3.

When x!=0, divide vi by x^i to recover y,z1,z2,z3 uniquely.

## GAOF6FIB1: the complete boundary count

Write the original target as (c,A,B,D,Z). The preimage of c=0 is the
disjoint union x=0 and gamma=0,x!=0, since gamma=1 at x=0.
The exact weighted jets give F6(0,y,z1,z2,z3)=(0,H1,H2,H3,H4), where

    H1=y,
    H2=29z1-(4636/5)y^2,
    H3=(6039y^3-5205y*z1+25z2)/5,
    H4=(74000249y^4-463585y^2*z1+125y*z2
        -26025z1^2-50z3)/25.

The diagonal coefficients are (1,29,5,-2), so this triangular map is an
isomorphism. It contributes exactly one point to every boundary fibre.

On gamma=0,x!=0, cancellation in the four expressions E_i gives

    A=(1-2u1)/x,
    B=(u1-u1^2)/x^2,
    D=(u1^2+u2)/x^3,
    Z=(u1^4-2u1*u2+u2+u3)/x^4.

Set delta=A^2+4B. The first two equations are equivalent to

    u1=(1-xA)/2, delta*x^2=1.

The remaining two determine uniquely

    u2=x^3*D-u1^2,
    u3=x^4*Z-u1^4+(2u1-1)*u2.

The stage inverse then reconstructs a unique source for each nonzero x.
Therefore every c=0 fibre has three points if delta!=0 and one point
if delta=0. There is no omitted boundary target, and no other source
piece has been left out. The full source fibre, not just its x=0 part,
has now been counted.

## Nonzero chart and three free sextic coefficients

For c!=0 put Y=(cA,c^2 B,c^3 D,c^4 Z). Gao's exact inverse is

    a(s)=Y3-2sY2+s^2Y1,
    b(s)=Y4+Y2^2-Y1Y3-psi(s),
    R(s)=2s^6-2s^5-s^3+(Y1+1)s^2
         +(Y4+Y2^2-Y1Y3-2Y2+Y1)s+Y3-Y2,
    gamma(s)=Y1-X1(s,a(s),b(s))=R'(s)-2R(s).

All sweep residuals vanish when R=0. At a root gamma=R', hence precisely
the simple roots lift. Their reconstruction is x=c/gamma, followed by
u1=s/gamma, u2=a/gamma^3, u3=b/gamma^4 and the stage inverse.
Every actual source over c!=0 occurs in this way, uniquely.

Write

    P=R/2=s^6-s^5-s^3/2+A0*s^2+B0*s+C0,
    A0=(Y1+1)/2,
    B0=(Y4+Y2^2-Y1Y3-2Y2+Y1)/2,
    C0=(Y3-Y2)/2, h=Y2.

This is a polynomial coordinate change on Y-space. Its inverse is

    Y1=2A0-1, Y2=h, Y3=h+2C0,
    Y4=2B0-h^2+(2A0-1)*(h+2C0)+2h-(2A0-1).

Thus the multiplicity classification is entirely in coefficient space
A3_(A0,B0,C0); all its nonzero-chart target strata are products with
C* times A1. This also follows from Gao's identity factor, but does not
require importing its longer source coordinate transformation.

Let a_g=deg gcd(P,P') and b_g=deg gcd(P,P',P''). The number of simple
roots, and hence of source points, is 6-2a_g+b_g. The complete table is:

| Partition | (a_g,b_g) | Fibre size | Coefficient stratum |
| --- | --- | ---: | --- |
| 1+1+1+1+1+1 | (0,0) | 6 | nonempty open |
| 2+1+1+1+1 | (1,0) | 4 | nonempty surface stratum |
| 2+2+1+1 | (2,0) | 2 | A1 minus 13 points |
| 2+2+2 | (3,0) | 0 | one point, already Gao's omitted case |
| 3+1+1+1 | (2,1) | 3 | A1 minus 9 points |
| 3+2+1 | (3,1) | 1 | six points |
| 3+3 | (4,2) | 0 | empty |
| 4+1+1 | (3,2) | 2 | three points |
| 4+2 | (4,2) | 0 | empty |
| 5+1 | (4,3) | 1 | empty |
| 6 | (5,4) | 0 | empty |

There are eleven degree-six partitions, all listed. The gcd tests,
boundary quadratic, and formulas below classify every target; the two
partitions sharing (4,2) are both impossible in this family.

## Excluded partitions and the omitted coefficient point

For 6,4+2,2+2+2 the polynomial is a square. Its leading coefficients
force it to be q(s)^2, where

    q=s^3-s^2/2-s/8-5/16.

This cubic is squarefree, so only 2+2+2 occurs. Its coefficient point is
(A0,B0,C0)=(21/64,5/64,25/256). This is Gao's prior omitted locus after
restoring c and h. A partition 3+3 would be (s^2+u*s+v)^3. The first
two constraints force u=-1/3,v=-1/9, whose s^3 coefficient is 5/27,
contradicting -1/2.

For 5+1 write P=(s-r)^5*(s-(1-5r)). Its s^4 coefficient forces
r=0 or 1/3; neither gives the required s^3 coefficient. Thus it is
also excluded. These arguments exhaust all impossible rows.

## Triple and quadruple roots: all exceptional coefficient points

A triple root r forces

    P=(s-r)^3*T3(s),
    T3=s^3+(3r-1)s^2+(6r^2-3r)s+10r^3-6r^2-1/2.

There is no denominator or exceptional parameter in this family. Set

    K4=40r^3-20r^2-1,
    K32=4320r^6-4320r^5+540r^4-448r^3+384r^2+36r+35.

Exact identities are T3(r)=K4/2 and disc(T3)=-K32/4. The two polynomials
are coprime, with nonzero discriminants

    disc(K4)=-75200,
    disc(K32)=-21557457167069587580544614400.

Therefore K4=0 gives exactly three 4+1+1 points, K32=0 gives exactly
six 3+2+1 points, and their complement is 3+1+1+1, isomorphic to A1
minus nine points. The excluded partitions above and the coprimality
ensure that none of these degenerations hides a higher multiplicity.
The distinguished triple/quadruple root is unique; the reverse maps
are regular on constant-partition strata by monic gcds/subresultants.

## Two double roots: a rational curve stratum

Let sigma,pi be their sum and product. Then

    P=Q2(s)^2*T2(s),
    Q2=s^2-sigma*s+pi,
    T2=s^2+(2sigma-1)s+3sigma^2-2sigma-2pi.

The remaining leading-coefficient constraint gives

    pi=(8sigma^3-6sigma^2-1)/(12sigma-4).

There is no point at sigma=1/3, since the numerator is nonzero there.
Define

    J4=5sigma^3-5sigma^2-1,
    J222=8sigma^3-8sigma^2+sigma+3,
    J321=1440sigma^6-2880sigma^5+2220sigma^4-728sigma^3
         +48sigma^2+24sigma+13.

Then

    disc(Q2)=-J4/(3sigma-1),
    disc(T2)=-J222/(3sigma-1),
    Res(Q2,T2)=J321/[16*(3sigma-1)^2].

The polynomials 3sigma-1,J4,J222,J321 are pairwise coprime, and each
is squarefree. Consequently the exact 2+2+1+1 stratum is A1_sigma
minus 1+3+3+6=13 points. Its coefficient parameterization is an
isomorphism: its regular inverse is obtained from the coefficients of
the monic gcd(P,P')=Q2.

The three J4 points give the three 4+1+1 coefficient points. The three
J222 points all map to the same 2+2+2 coefficient point, corresponding
to the three choices of two among its three double roots. The six J321
points give the six 3+2+1 points. None of these punctures is counted
as part of the exact 2+2+1+1 stratum.

## Occurrence and all cardinalities

The prior rational six-point witness has

    R=2s*(s+1)*(s-1/2)*(s-3/2)*(s^2+5/4).

For four points use Y=0, where R=s^2*(s-1)*(2s^3-1), also prior data.
For three points use r=0 in the triple-root family. For two points use
sigma=0 in the two-double-root family, with Q2=s^2+1/4 and T2=s^2-1/2.
The six K32 points give one-point fibres, and q^2 gives empty fibres.
All satisfy the coefficient constraints and the stated exact partitions.
Hence the global possible cardinalities are exactly {0,1,2,3,4,6}.
The boundary contributes sizes one and three only. Five is impossible.

## GAOF6NP1: the full nonproper-value locus

Let D0(A0,B0,C0)=disc_s(P). Its multiple-root incidence has a chosen root
r and the equations P(r)=P'(r)=0. Writing f=s^6-s^5-s^3/2, they solve as

    B0=-f'(r)-2A0*r,
    C0=r*f'(r)-f(r)+A0*r^2.

Thus the incidence is A2_(r,A0). Its map to coefficient space is finite,
since r is integral by the monic equation P(r)=0 and all other variables
are polynomial in r and the coefficients. It is generically one-to-one
over the single-double-root stratum. This A2 is the normalization of
the irreducible discriminant surface. The discriminant polynomial is
reduced, because its C0 derivative is nonzero at (1/2,0,0).

Use Gao's polynomial target coordinates

    q=Z+B^2-A*D, k=cD-B, l=cq-2D,
    (c,A,B,D,Z) -> (c,A,k,l,q).

Their inverse is D=(cq-l)/2, B=cD-k, Z=q-B^2+A*D. In these coordinates
the monic sextic is

    P_c=s^6-s^5-s^3/2+(1+cA)*s^2/2
        +(cA+2c^2*k+c^3*l)*s/2+c^2*k/2.

Its discriminant has exact c-valuation two. The quotient

    E(c,A,k,l)=disc_s(P_c)/c^2

is polynomial and satisfies

    E(0,A,k,l)=-(27/256)*(A^2-4k).

The full nonproper-value locus is V(c*E) times the free q-line. On c!=0,
off the discriminant the six simple-root inverse branches exist locally
and remain bounded, excluding escape. At a multiple root, approach by
simple roots; gamma=2P' tends to zero, so x=c/gamma diverges. Such
approximations exist at every root-incidence point: this incidence is
itself A3_(s,A0,B0), by solving for C0, and its simple-root open is
dense. Fixing h and using the invertible coefficient change realizes
the approximating coefficient targets as actual sweep targets.

Every c=0 target is nonproper as well. Approach it along (c,A,B,D,Z)
with the last four original coordinates fixed. The limiting polynomial is

    P0=s^2*(s-1)*(s^3-1/2).

Choose a root alpha of 2s^3-1. It is a simple root of P0, is not 1,
and has gamma=2P0'(alpha)!=0. Its nearby analytic root branch therefore
has x=c/gamma tending to zero. If the reconstructed source coordinates
were bounded, all vi=x^i*(source coordinate) would tend to zero. The
stage would then give gamma->1,u1->1, hence s->1, contrary to s->alpha.
Thus the inverse branch escapes to infinity while its target approaches
the specified boundary point.

Finally E is irreducible: after inverting c, the transformed target
coordinates identify the three coefficients with an affine A3 and E
with its irreducible discriminant up to a unit. The nonzero restriction
at c=0 excludes a factor c. The complete nonproper locus consequently
has exactly two irreducible components. They intersect on c=0 along
A^2-4k=0, equivalently A^2+4B=0 in original coordinates, precisely the
boundary one-point locus.

## Reproduction and assurance

Run from the repository root:

```sh
timeout 60s .venv/bin/python research/jacobian/gao-f6-fibres/verify.py
```

The checker verifies exact identities, weighted jets, determinant factors,
inverses, all partitions, nonvanishing discriminants, puncture disjointness,
witnesses and source hashes. Geometric exhaustiveness, regular gcd maps,
normalization and analytic escape are written arguments, not finite-checker
outputs. This is not an independent proof review or formal verification.
No external review, author endorsement, publication or novelty priority
is claimed. MATH_STATUS.json is authoritative; unrelated missing historical
artifacts still prevent a fully green repository-wide validation.
