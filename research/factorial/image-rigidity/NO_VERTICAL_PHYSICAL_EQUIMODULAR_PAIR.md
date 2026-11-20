# No vertical physical equimodular pair

This note proves an elementary obstruction for physical stationary pairs. It is independent of the square-root-gap eventual theorem. No global labeling, activation or Strong Factorial theorem is inferred.

Fix real a>1 and integer d>=2, D=a+d, and let
Omega(S)=S^d(D-S)/(S-a), lambda(S)=d S^D exp(-S)/(S-a),
using the physical wedge -pi/d<arg S<0 and its positive boundary where the expressions are defined. Powers with noninteger exponent in the displayed wedge use its logarithm.

Claim: two DISTINCT points S1,S2 of this wedge cannot satisfy simultaneously Omega(S1)=Omega(S2), |lambda(S1)|=|lambda(S2)|, and Re S1=Re S2.

A common value Omega=0 cannot give a distinct pair: the only permitted zero is S=D; S=0 is outside the stated domain. Thus the common value is nonzero.

Suppose otherwise. Write S_j=X-i y_j with y1>y2>=0. The wedge gives X>0. Set R=X^2+y2^2>0 and t=(X^2+y1^2)/R>1. The equal-modulus conditions, first for lambda and then for Omega divided by lambda, give
(tR+C)/(R+C)=t^D, C=a^2-2aX;
(tR+E)/(R+E)=t^a, E=D^2-2DX.
For r>1 let K_r=(t^r-t)/(t^r-1)>0. Then C=-R K_D, E=-R K_a, whence
R=a D d/(D K_D-a K_a),
X=(D^2 K_D-a^2 K_a)/(2(D K_D-a K_a)).
Indeed K_r=1-(t-1)/(t^r-1) is strictly increasing in r. Since K_D>K_a, the denominator is positive.

Here is the needed sharp location of X. Put v=(log t)/2>0. Then
K_r=e^v sinh((r-1)v)/sinh(rv).
The ratio z(v)=K_D/K_a strictly decreases from a(D-1)/(D(a-1)) to1. Indeed its logarithmic derivative is
-[F(D)-F(D-1)]+[F(a)-F(a-1)]<0,
where F(x)=x coth(xv) and
F''(x)=2v csch^2(xv)[xv coth(xv)-1]>0 for x>0. The strict last sign follows from z*cosh(z)-sinh(z)>0 for z>0: this function vanishes at0 and has derivative z*sinh(z)>0.
The function X(z)=(D^2 z-a^2)/(2(Dz-a)) has derivative -a D d/[2(Dz-a)^2]<0. At the two endpoint ratios it equals a+(d-1)/2 and (D+a)/2. Consequently
  a < a+(d-1)/2 < X < (D+a)/2 < D.

Now consider y between y2 and y1, and define
 F1(y)=D/(X^2+y^2)-1/((X-a)^2+y^2),
 F2(y)=1/((D-X)^2+y^2)-a/(X^2+y^2).
The derivatives with respect to y of log|lambda(X-i*y)|^2 and log|Omega(X-i*y)/lambda(X-i*y)|^2 are respectively2*y*F1 and2*y*F2. Hence the same two equal-modulus conditions say int y F1 dy=int y F2 dy=0. The numerator of F1 has positive coefficient D-1 at y^2, so it crosses sign once from negative to positive; the numerator of F2 has negative coefficient1-a and crosses once in the opposite direction. Their vanishing weighted integrals force both crossings strictly inside the interval. With c1,c2 denoting their separate crossings, int (y-c1) F1 dy>0 and int (y-c2) F2 dy<0. Therefore int F1 dy<0 and int F2 dy>0.

The continuous physical argument of Omega(X-i y) is
 theta(y)=-d atan(y/X)+atan(y/(D-X))+atan(y/(X-a)).
Differentiating and using D=a+d gives exactly
 theta'(y)=(D-X) F2(y)-(X-a) F1(y).
It follows that theta(y1)-theta(y2)>0. But both physical arguments lie in(-pi,pi): the first summand lies in(-pi,0] and the other two have sum in[0,pi), since a<X<D. Thus their strictly positive difference is less than2pi. Equal nonzero Omega values require an integer multiple of2pi, a contradiction.

This proves nonvanishing of Re(S_R-S_L) at every distinct physical equimodular pair for ALL real a>1 and integer d>=2, without an everywhere fixed vertical divider. It does not by itself prove actual-cycle activation, globally uniform contour estimates, turning-point neighborhood exclusion, any cutoff, or all-order SFC. The above-cusp local radial Airy identity is already canonical in COALESCING_BINOMIAL_SADDLE_BRIDGE.md Section4 and is not new.

## Endpoint a=1

**Corollary NoVerticalPhysicalEquimodularPair.** The same exclusion holds for every real a>=1 and integer d>=2.

Only a=1 remains to be checked. The second equal-modulus equation has K_1=0, so E=D^2-2DX=0 and X=D/2. Because D=1+d and d>=2, one has 1<X<D. The function F2 is identically zero: its numerator is X^2-(D-X)^2=0. The first equal-modulus equation still gives integral_(y2)^(y1) y F1(y)dy=0. Its numerator has strictly positive y^2 coefficient D-1=d, so its unique crossing c1 is strictly inside the interval and positive. Therefore integral F1<0 exactly as in the main proof. The argument identity now gives

    theta(y1)-theta(y2)=-(X-1)*integral F1>0.

The same physical argument bounds put this difference strictly below2pi. This again contradicts equality of nonzero Omega values. No division by K_1 or a-1 occurs in this endpoint proof.

## Exact support and reproduction

Run `python3 verify_no_vertical_physical_equimodular_pair.py --repo-root /path/to/jacobian-research --output review.json`. The standard-library guard verifies full denominator-cleared identities, a differential-ring derivation of the hyperbolic convexity formula, exact rational candidates satisfying both modulus equations, their two strict crossings, and sufficient rational certificates placing candidates in the physical wedge. It also authenticates this note, the checker, inputs, and the retained stationary context. The stationary-context source is not a logical premise of this self-contained proof.

Finite controls support the written universal argument and do not replace its integration, convexity or strict-argument proofs. The earlier high-precision constructed-candidate search is optional numerical research evidence and is not a proof premise or runtime dependency. No historical priority, formal verification, external review, global fixed branch order, cycle activation, cutoff or Strong Factorial assertion is claimed.
