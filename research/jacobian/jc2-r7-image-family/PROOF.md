# All gap-21 image coefficients lie on two rational curves

14 November 2026. PF2R7FAM1 gives a complete set-theoretic classification
of the normalized coefficient family used in the normal-r7 research.
PF2R7FAM2 identifies the exact nonimmersion loci in these two charts.
The result applies to the image curve even when additional ramified
primes are present. It does not assume the excluded entire-single-prime
stratum and does not construct a finite cover or a Keller map.

The conclusion is about complex coefficient points. No equality of the
possibly nonreduced coefficient schemes is asserted. The three small
graph relations below have literal memberships in the original five
gap equations; no printed Groebner basis is a proof premise. Exact
algebra supports the written completeness proof. There is no formal
verification, independent external review or novelty-priority claim.

## 1. Normalization and explicit answer

Use the retained PF2D6O1 polynomial normalization

    p=t^6+a*t^4+b*t^3+c*t^2+d*t,
    q=[p^(5/3)+B*p^(4/3)+D*p^(2/3)+E*p^(1/3)]_+.

The free terms C*p+F were removed by affine target changes. The five
equations rho1=rho3=rho5=rho7=rho9=0 prescribe the gap. Birationality
forces b nonzero by the earlier parity identity. A complex parameter
scaling and diagonal target scaling put b=1. A residual cube-root
choice remains; this is not a claim of unique representatives modulo
all affine transformations.

Every solution with b=1 belongs to one of the following two families.
The remaining coefficients c,D,E are the polynomial graph functions in
Section2; [family.py](family.py) gives all coefficients explicitly.

| Family | Parameter | a | d | B |
|---|---|---|---|---|
| L | d in C | 2d | d | -7d/6 |
| K | x in C, x nonzero | -(15x^6+8x^3+4)/(4x^2) | -(15x^6+4x^3+4)/(12x^2) | -(5x^3+1)/(3x^2) |

On L the other coefficients are particularly small:

    c=d^2,
    D=(18-35*d^3)/81,
    E=7*d*(5*d^3-9)/243.

Both substitutions satisfy all five equations. Their next obstructions are

    rho11(L)=1/19683,
    rho11(K)=-2/(59049*x^3).

Thus every point has gap exactly21. Every point is birational: the field
degree [C(t):C(p,q)] divides both6 and10, hence is one or two. If it were
two, its nontrivial automorphism would be a Mobius transformation fixing
p and q. It fixes infinity, the unique pole of p, so is affine. A
nontrivial affine involution has form t maps to -t+beta. The missing t^5
coefficient of p forces beta=0. Then p is even, contradicting b=1.

The images L and K are closed irreducible rational coefficient curves.
L has normalization A1. K has normalization Gm. Their intersection
consists exactly of the three points x^3=-2/3 on K, with d^3=-2/3 on L.
These are coefficient curves, not the target plane image C itself.

## 2. Three certified graph relations

Write f1,...,f5 for the integer-primitive normalizations of the five
rho equations with b=1. Frozen rational polynomial multipliers give
three identities g_j=sum_i M_ij*f_i, where

    g0=18Bd+45d^2-6Ba-30da+7a^2-6c,
    g1=108B^2a+405d^2a-108Ba^2-378da^2+a^3
       +216Bc-216cd+252ca+324D-72,
    g2=810d^3-810d^2a-36Ba^2+189da^2+7a^3
       +108Bc+243cd-126ca+81D-18.

The checker reconstructs the five f_i and three g_j, then replays all
three memberships with separate sparse Fraction arithmetic. Discovery
used the retained old basis only to suggest these targets; the literal
memberships, not that basis, justify their use here.

Solve g0 for c and g2 for D. Solve rho1 for E; its coefficient is1/3.
Every denominator used in these three graph functions is a nonzero
rational constant. In particular

    c=-Ba+3Bd+7a^2/6-5ad+15d^2/2.

After this substitution g1 is exactly

    9*(a-2d)*Q,
    Q=36B^2-108Ba+288Bd+95a^2-510ad+675d^2.

Also rho5=-(2a/3)*rho3 on this graph. Denote the integer-primitive
residual rho3,rho7,rho9 equations by F3,F7,F9. Their fixed normalization
is part of [family.py](family.py). Every original solution therefore has
the graph values and satisfies (a-2d)*Q=F3=F7=F9=0.

## 3. The exceptional plane a=2d contributes only L or Q=0

Set u=6B+7d. On a=2d the first two residuals are

    F3=d*u*f, F7=u*g,
    f=16d^3+d*u^2+8,
    g=416d^6-16d^5*u+26d^4*u^2-d^3*u^3
      +216d^3-7d^2*u+4.

If u=0 we are on L. Suppose u nonzero. If d=0, F7=4u is nonzero.
Thus d is nonzero and f=g=0. Put l=d^2*u+8d^3+4. Exact identities give

    g=(-d^2*u+26d^3)*f+l,
    d^3*f-(d^2*u-8d^3-4)*l=8*(2d^3+1)*(5d^3+2).

If 2d^3+1=0, l=0 forces u=0, a contradiction. If 5d^3+2=0,
l=0 forces u=2d, hence 6B+5d=0. But on this plane

    Q=(6B+5d)*(6B+7d).

Consequently every solution belongs to L or Q=0. No exceptional d or
rank-drop branch has been removed by a generic division.

## 4. The conic branch and its complete Laurent parameter

On Q=0 introduce the homogeneous linear expressions

    X=6B-7a+21d, Y=6B-13a+33d, Z=5X-3Y,
    N=5X^2-2XY+3Y^2.

Their exact relation is Q=(N-6Zd)/6. If Z=0, Q=0 gives
Y=5X/3 and N=10X^2=0, hence X=Y=0. This is L. Off L, Z is nonzero and

    d=N/(6Z), a=2d+(X-Y)/6, B=(X+7a-21d)/6.

Set R=5X^4Y^2+2Z^3. Substitution into the first two residuals gives

    F3=X*(5X^2-6XY+5Y^2)*R/Z^4,
    F7=R*P7/(2Z^6).

The checker reconstructs the polynomial P7. Suppose R were nonzero.
If X=0, its cofactor P7=-54Y^4 is nonzero since Z is nonzero. Otherwise
put r=Y/X. The first residual forces k(r)=5r^2-6r+5=0. Modulo k, the
degree-three-in-X coefficient of P7/X^4 is zero and its constant is

    -256*(139r-95)/125.

Thus r=95/139, but k(95/139)=62500/19321 is nonzero. This is impossible.
Therefore R=0. In particular X and Y are both nonzero off L.

Now put x=XY/Z, which is nonzero. Divide R=0 by Z^2 to obtain

    5*x^2*X^2+2Z=0.

Combining this with xZ=XY and Z=5X-3Y gives successively

    Y=-(5/2)*x^3*X,
    X=-(3x^3+2)/x^2,
    Y=(5/2)*x*(3x^3+2).

Substituting these into d,a,B gives exactly K in Section1. This covers
every off-L solution. Its three missing values x^3=-2/3 already lie on
L, and direct substitution verifies K there too. Conversely both L and
K satisfy every original equation, so the two images exhaust the complex
coefficient points. The proof includes d=0, a=2d, Z=0 and every excluded
denominator stratum explicitly.

For the normalization statement, the inverse x=XY/Z is valid on a dense
open of K. Both x and x^-1 are integral over the coefficient ring: from
B=-(5x^3+1)/(3x^2) one has

    5x^3+3B*x^2+1=0, (x^-1)^3+3B*x^-1+5=0.

Thus Gm maps finitely and birationally onto its closed coefficient image
and is its normalization. The simultaneous equations a=2d,6B+7d=0
on K have exactly x^3=-2/3. Their d values have d^3=-2/3 and are distinct.
This proves PF2R7FAM1 at the set-theoretic scope stated above.

## 5. Exact immersion loci

Let A(t),B(t) be the numerator polynomials of p'(t),q'(t) after the
displayed family substitution and rational simplification. Clearing
these denominators never loses a finite parameter of L or a nonzero
parameter of K. The checker computes the exact resultants in t:

    Res_L=-1011316948992*(3d^3+2),
    Res_K=5263136739342894229381323813814272
          *x^67*(3x^3+2)*(15x^3-4).

The derivative degrees are fixed five and nine, with nonzero leading
coefficients on their parameter domains. The resultant criterion
therefore says exactly:

- L is immersed except at d^3=-2/3;
- K is immersed except at x^3=-2/3 or x^3=4/15.

The first K locus is exactly the intersection with L. This agrees with
the retained two-class nonimmersion classification PF2R7CUSP1. No
generic-node or complete equisingularity assertion follows from these
derivative resultants. Multibranch collisions and tangencies require
their separate divided-difference equations. This is PF2R7FAM2.

## 6. What changes in the global research problem

The known ramified image C in a residual degree-six r7 candidate is now
an explicit member of L or K. If C is immersed, its parameter avoids
the finite loci in Section5. The second ramified image, any unramified
omitted image, and their intersection parameters remain unknown. The
PF2OMIT6--8 length, radial and S6 matching constraints remain necessary;
they are not supplied by this coefficient classification.

In particular the generic member of one coefficient component cannot
stand in for its entire component in a braid calculation. Additional
discriminant strata of the divided differences must be computed before
transporting complement presentations. The earlier capped braid and
finite-slope probes remain historical experiments. This result gives
two explicit one-parameter inputs for that next global calculation;
it does not close whole normal r7, other terminal types, (75,125) or JC2.
