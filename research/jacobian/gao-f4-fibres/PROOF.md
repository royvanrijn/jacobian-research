# Complete complex fibres of Gao's F4, and an elliptic multiplicity stratum

13 November 2026. Claims **GAOF4FIB1**, **GAOF4ELL1**, **GAOF4NP1**.
All varieties and fibres below are over C, with reduced structures on
strata. Fibre counts mean distinct source points. All finite fibres are
reduced because the map is etale.

## Source and contribution boundary

The map is exactly F4 in [Gao, arXiv:2608.00222v2](https://arxiv.org/html/2608.00222v2),
Sections 4.4.1 and 7. Gao already proves the simple-root correspondence on
the nonzero first-coordinate chart, the five omitted curves, and that the
restriction x=0 is an isomorphism onto the target boundary. Sections 11.1
and A.3 explicitly leave the complete F4--F7 fibre classifications open.

The earlier [Korea Superintelligence Labs note](https://ideosphere.ai/papers/keller-fibers)
studies selected rational fibres and reality constraints. Its KF1/KF13/KF14
include a rational five-point F4 fibre and a root correspondence on a
restricted chart. Those are prior results, not claimed additions here.
The pinned 11-page PDF and its claim list were examined; the complete
boundary count and the elliptic stratum below were not found there.
This is a comparison with inspected sources, not an exhaustive priority
claim. Gao's companion scripts are described in Appendix B but are not
present in the downloaded two-file source archive; no replay of those
unavailable scripts is claimed.

The additions are: all boundary fibres, every multiplicity stratum and
its cardinality, a smooth elliptic completion of the two-double-root
stratum, and the full nonproper-value hypersurface with rational
normalization of its discriminant on the nonzero chart. This does not
solve F5 or F7, the planar Jacobian conjecture, or unrestricted HC4.
F6 has since received a [separate classification](../gao-f6-fibres/PROOF.md).

## Definition and coordinate verification

In source coordinates (x,y,z,t), put

    gamma=1+xy+x^2*t, U=1+xz,
    V=-7/9+(22/9)xy+(8/3)xz, s=gamma*U, a=gamma*V,
    C=gamma*x.

Define the three-component vectors

    A(s)=(s/3-4s^2/3, 2s^2/3-8s^3/9,
          2s+7s^3/9-2s^4/3),
    B(s)=(s^2-1,s^3,s^4+s^2+2), Delta(s)=(1,s,s^2),
    S=A(s)+a*B(s)+gamma*Delta(s).

Then F4=(C,S1/C,S2/C^2,S3/C), using the polynomial extensions.
For direct reproduction after cancelling powers of gamma, the last three
coordinates are N1/x, N2/x^2, N3/x, where

    N1=1+U/3-V-(4/3)gamma*U^2+gamma^2*V*U^2,
    N2=U+(2/3)U^2-(8/9)gamma*U^3+gamma^2*V*U^3,
    N3=2U+2V+gamma^2*(U^2+7U^3/9+V*U^2)
       -(2/3)gamma^3*U^4+gamma^4*V*U^4.

The checker expands these numerators and verifies divisibility by x,x^2,x.
The resulting degrees are (4,11,12,21). To verify etaleness without a large
expanded determinant, the four successive Jacobians are

    (x,y,z,t) -> (x,gamma,U,V): -22*x^4/9;
    (x,gamma,U,V) -> (x,gamma,s,a): gamma^2;
    (x,gamma,s,a) -> (C,S1,S2,S3): 2*gamma^2;
    (C,S1,S2,S3) -> (C,S1/C,S2/C^2,S3/C): C^-4.

Their product is -44/9 on a dense open, hence everywhere as a polynomial
identity. In particular every fibre is zero-dimensional and reduced.

## GAOF4FIB1: boundary fibres

Write a target as (c,A,B,D). When c=0, C=x*gamma=0 splits its preimage
into the disjoint pieces x=0 and gamma=0,x!=0: at x=0 gamma is 1.
There can be no third piece.

Gao's x=0 restriction is

    (0,-26y/9-35z/9,
       37y^2/9+16yz/3+11z^2/3-22t/9,
       20y/3+29z/3).

It has precisely one inverse at every boundary target:

    y0=-(87A+35D)/18, z0=(30A+13D)/9,
    t0=(9/22)*(37y0^2/9+16y0*z0/3+11z0^2/3-B).

On gamma=0,x!=0, the other target equations become

    A=(1+U/3-V)/x, B=(U+2U^2/3)/x^2, D=2(U+V)/x.

Set K=2A+D and d0=3K^2-32B. The first and third equations give uniquely

    U=3Kx/8-3/4, V=Dx/2-U.

Substitution into the second gives d0*x^2=12. If d0 is nonzero this has
exactly two nonzero complex roots; if d0=0 it has no solution. For either
root, the original source is recovered uniquely by

    xy=(9/22)*(V+7/9-(8/3)*(U-1)),
    z=(U-1)/x, y=(xy)/x, t=(-1-xy)/x^2.

These identities also verify gamma=0. Therefore

| Boundary condition | Fibre size |
| --- | ---: |
| c=0 and 3(2A+D)^2-32B != 0 | 3 |
| c=0 and 3(2A+D)^2-32B = 0 | 1 |

For example the origin has its single source point at the origin.
The target (0,0,0,1) has three points; one of the two with x!=0 is
(2,10/11,-1/2,-31/44). Counting only the x=0 isomorphism would miss two
points on the dense open of the boundary.

## Nonzero chart and the complete multiplicity table

For c!=0 write Y=(cA,c^2 B,cD) and define

    P(s)=s^5+s^4+A0*s^3+B0*s^2+D0*s+(B0+6),
    A0=4+(9/2)Y1, B0=-6-9Y2, D0=9Y1+(9/2)Y3.

The affine coefficient map Y -> (A0,B0,D0) is invertible. For clarity,
Gao's lifting proof can be reproduced without exceptional s-values. Put

    phi=2s-2s^3/9-2s^4/9,
    a(s)=(Y3-2sY2+s^2Y1-phi)/2,
    gamma(s)=Y1-A1(s)-a(s)*B1(s), R=(2/9)P.

The sweep residuals are (0,R/2,sR), and

    R'-2*gamma=sR.

Thus at a root, gamma=P'/9. Precisely the simple roots lift. For each,
x=c/gamma and the stage inverse above (with t=(gamma-1-xy)/x^2)
recovers one and only one source point. All divisions are nonzero.

Let a_g=deg gcd(P,P') and b_g=deg gcd(P,P',P''). A root of multiplicity
m contributes m-1 to a_g and max(m-2,0) to b_g, so the number of simple
roots is 5-2*a_g+b_g. All degree-five partitions are exhausted by:

| Partition | (a_g,b_g) | Fibre size | Occurrence |
| --- | --- | ---: | --- |
| 1+1+1+1+1 | (0,0) | 5 | nonempty open |
| 2+1+1+1 | (1,0) | 3 | nonempty surface stratum |
| 2+2+1 | (2,0) | 1 | elliptic curve minus 13 points |
| 3+1+1 | (2,1) | 2 | A1 minus 13 points |
| 3+2 | (3,1) | 0 | 5 points |
| 4+1 | (3,2) | 1 | 5 points |
| 5 | (4,3) | 0 | empty |

These are statements in Y-space. Their corresponding target strata on
c!=0 are products with C*. The table and the boundary equation give
constructible tests for every target, including all exceptional fibres.
No numerical root tests or generic-witness extrapolation are required.

To certify nonemptiness of the positive-cardinality cases, suitable P are

    5 points: s^5+s^4+4s^3-15s^2-9;
    3 points: s^2*(s-1)*(s^2+2s+6);
    2 points: (s-1)^3*(s^2+4s+17/2);
    1 point:  s*(s^2+s/2-6)^2.

The one-point example has Y=(-7/2,0,15). At c=1 its unique source is
(1/4,383/44,-4,145/11), directly checked in the expanded polynomials.
The triple-double and quadruple strata below establish nonemptiness of
the remaining occurring partitions. Globally the possible fibre sizes
are exactly {0,1,2,3,5}; four is impossible.

## The triple-root strata

Write

    P=(s-r)^3*(s^2+(1+3r)s+b).

The constant-minus-quadratic coefficient constraint is

    -b*r^3+3b*r-8r^3-3r^2-6=0.

The denominator e=r(r^2-3) never vanishes on this incidence curve:
8r^3+3r^2+6 is coprime to e. Consequently

    b=-(8r^3+3r^2+6)/e.

If Q is the residual quadratic, direct identities give

    Q(r)=H4(r)/e, disc(Q)=3h(r)/e,
    H4=4r^5+r^4-20r^3-6r^2-6,
    h=3r^5+2r^4+2r^3-2r^2-r+8.

Here disc(H4)=-43861243392 and disc(h)=1292134656. Both polynomials
are coprime to e and to each other. Thus H4=0 gives five distinct 4+1
points, h=0 gives five distinct 3+2 points (already Gao's omitted targets),
and the complement is the 3+1+1 stratum. Its parameter r excludes three
denominator roots and ten other roots: A1 minus 13 points, or P1 minus
14 points. Each polynomial has a unique root of multiplicity at least
three, making the parameterization injective and regular in reverse on
each constant-partition stratum via the monic gcd/subresultants.

The potential partition 5 would force P=(s+1/5)^5 from its s^4
coefficient. Its constant coefficient minus its s^2 coefficient minus 6
is -18999/3125, so it is impossible. In particular none of the listed
quadruple points conceals a quintuple root.

## GAOF4ELL1: the two-double-root elliptic curve

Let sigma and pi be the sum and product of the two double roots. Then

    P=(s^2-sigma*s+pi)^2*(s+1+2*sigma).

The family constraint becomes the affine cubic

    H22=(1+2*sigma)*(pi^2-sigma^2)-2*(1+sigma)*pi-6=0.       (E)

Its coefficient map is

    A0=-3*sigma^2-2*sigma+2*pi,
    B0=sigma^2*(1+2*sigma)+2*pi*(1+sigma),
    D0=pi^2-2*sigma*pi*(1+2*sigma).

The partition is exactly 2+2+1 when

    sigma^2-4*pi != 0,
    pi+(1+2*sigma)*(1+3*sigma) != 0.

On this stratum gcd(P,P')=s^2-sigma*s+pi. Its coefficients are regular
functions on the constant-gcd-degree locus, computed by a nonzero
subresultant normalized to be monic. This gives an algebraic inverse to
the coefficient map, proving an isomorphism rather than only a pointwise
bijection.

The projective closure of (E) is

    (T+2*sigma)*(pi^2-sigma^2)-2*(T+sigma)*pi*T-6*T^3=0.

It has exactly three points at infinity: [0:1:0], [1:1:0], [1:-1:0],
all smooth by its partial derivatives. Away from sigma=-1/2, complete
the square using W=(1+2*sigma)*pi-(1+sigma):

    W^2=4*sigma^4+4*sigma^3+2*sigma^2+14*sigma+7.

The quartic has discriminant -12854016, hence is squarefree. This affine
chart is irreducible and smooth. At sigma=-1/2 the sole affine point is
pi=-6, with H22_pi=-1, also smooth. There is no component contained in
that vertical line or the line at infinity. The projective cubic is
therefore smooth and irreducible, and has genus one. Its rational points
at infinity supply an origin, so it is an elliptic curve over Q as well.

The forbidden divisor sigma^2-4*pi=0 is parameterized by
sigma=2r, pi=r^2. Equation (E) then becomes H4(r)=0, giving five distinct
4+1 points. For the other forbidden divisor, let r be the root occurring
three times and u=(-1-3r)/2 the other double root. Substituting
sigma=r+u, pi=r*u in (E) gives -3h(r)/4. These are the five 3+2 points.
They are disjoint from the first five, since a common point would have
partition 5. Thus the exact 2+2+1 stratum is the elliptic projective curve
minus 3+5+5=13 distinct points. The rational point (-1/2,-6) survives.

## GAOF4NP1: discriminant normalization and the nonproper locus

Let D(Y)=disc_s(P). Its zero set consists of polynomials with a multiple
root. The incidence of a chosen double root has the expression

    P=(s-r)^2*(s^3+(1+2r)s^2+b*s+d),
    (r^2-1)d+2r*b-r^2*(1+2r)-6=0.

Put T_r=r^2*(1+2r)+6. An explicit polynomial parameterization is

    b=r*T_r/2+(r^2-1)w, d=-T_r-2r*w,
    w=-(b-r*T_r/2)-r*(d+T_r)/2.

It identifies the incidence variety with A2_(r,w). The map to Y-space
is finite: r satisfies the monic polynomial P and b,d are polynomial
in r and its coefficients by division by (s-r)^2. It is generically
one-to-one, because a generic point has exactly one double root. The
three-point example above witnesses that nonempty open. Hence this A2
is the normalization of the irreducible discriminant surface.

The discriminant polynomial is reduced: at Y=0, its derivative in Y2 is
nonzero, as directly checked. Together with irreducibility of its zero
set, this proves D itself is irreducible up to scalar. Define

    E(c,A,B,D)=disc_s(P(s;cA,c^2 B,cD))/c^2.

Here the final D in the argument denotes the fourth target coordinate,
not the discriminant polynomial D(Y). The quotient is polynomial; the
exact c-valuation of its numerator is two, and

    E(0,A,B,D)=-393660*[3(2A+D)^2-32B].

The nonproper-value locus of F4 is exactly V(c*E). In particular it has
two irreducible hypersurface components, c=0 and E=0. To see this on
c!=0, away from the discriminant all five roots and their inverse
branches vary analytically and remain finite, excluding an escaping
sequence. At a multiple root, approach by simple roots; gamma=P'/9
tends to zero and x=c/gamma tends to infinity, giving an escaping
sequence. Such approximations exist by the nonempty simple-root open
in the finite root incidence, or by perturbing coefficients.

Every target (0,A,B,D) is also nonproper. Along targets (c,A,B,D), the
polynomial tends as c->0 to

    P0=s^2*(s-1)*(s^2+2s+6).

Choose either root alpha of q0=s^2+2s+6. It is simple, so a nearby
analytic root s(c) exists, with gamma(c)=P'(s(c))/9 tending to a nonzero
value. In the inverse formula z=(s-gamma)/c. Modulo q0,

    s-P0'/9 = -3s-4,

which is nonzero at either alpha. Thus z diverges while its target
tends to the specified boundary point. This proves the entire c=0
hyperplane is nonproper, including the three-point locus.

Finally, on c!=0 the weighted target change is an isomorphism, so E
is irreducible after inverting c by irreducibility of D(Y). Its nonzero
boundary restriction excludes a factor c, hence it is irreducible in
the full polynomial ring. The second component intersects the first
exactly along 3(2A+D)^2-32B=0, the boundary one-point locus.

## Verification and limits

Run `timeout 60s .venv/bin/python research/jacobian/gao-f4-fibres/verify.py`
from the repository root. The checker expands the polynomial map, checks
the determinant chain, boundary inverses and the other boundary equation,
Gao's lifting identity, all partition cases, the parameterizations,
discriminants, coprimalities, infinity points, and explicit rational
examples. It also verifies the source manifest hashes.

The geometric implications (etale fibres, regular gcd maps, genus,
normalization and analytic escape) are written arguments. This is not
formal verification or an independent end-to-end proof review. No external
review, source-author endorsement, publication, or novelty priority is
claimed. The authoritative claim scope is in MATH_STATUS.json. Local
checks do not resolve unrelated missing historical artifacts in the
repository-wide validation.
