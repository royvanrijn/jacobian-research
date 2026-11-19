# Higher starting heights need not reduce polynomially to cubic height

20 November 2026.

Work over any characteristic-zero field K. All phases are regular
unit-center phases

 p=t^a(C(x)+t²D(x))/x,
 C,D in K[x], C(0)=-D(0)=1, a>=3 odd.

Write the selected rational square root in REDUCED form X=P/J,
P(0)=J(0)=1, gcd(P,J)=1, and C=HJ²,D=-HP². Here H need not be
squarefree. The retained all-height endpoint criterion uses

 u=xX^a/C=x+O(x²), omega=W dx=X^(a-1)dx/C,
 T=W/u', Phi(u)=T(v(u)),

where v is the selected formal inverse. Purity is exactly
Phi(u)+Phi(-u)=2.

The proposed universal divisibility

 P^(a-1) divides H                                      (R0)

is FALSE, including on nonconstant-Phi families of prime degree.
It is equivalent to polynomiality of Ctilde=C/X^(a-1).
The formal formulas

 Ctilde=C/X^(a-1), Dtilde=D/X^(a-1)

always preserve u and Phi when regarded as rational functions, but
need not produce polynomials. This is not a counterexample to the
retained existential compensating-pole law.

## 1. A two-root Phi=1 counterexample already at degree15

**Theorem TwoRootHigherStartReductionObstruction.** Put

 U=1+x, V=1+2x/3, X=UV,
 C=V^(a-2) U^(a+3), D=-V^a U^(a+5).                  (R1)

For EVERY odd a>=3, the phase is all-pure and unsafe, of exact least
degree3a+6. Its selected coordinate and differential are

 u=x V²/U³,       u'=V/U^4=W,       Phi=1.             (R2)

To verify the derivative, multiply the logarithmic derivative by xUV:

 UV+(4/3)xU-3xV=1.

Thus u'/u=1/(xUV), and(R2) follows. Since Phi=1, the endpoint
criterion proves ALL positive pure moments vanish, without a
finite-prefix inference.

The reduced denominator J is1 and H=C. At V=0 (x=-3/2),

 ord(P)=1, ord(H)=a-2 < a-1,
 ord(u)=2, ord(omega)=1.

In particular

 Ctilde=U^4/V,       Dtilde=-V U^6,                    (R3)

and Ctilde has an actual simple pole. The desired polynomial
reduction fails. At U=0 (x=-1),

 ord(P)=1, ord(H)=a+3,
 pole(u)=3, pole(R)=3, R=u.

This is precisely a compensating pole of the kind forced somewhere
by the retained existential law. There is no contradiction with it.

The degrees are degC=2a+1,degD=2a+5. The least-degree normal-form
formula max(degC+a-1,degD+a+1) therefore gives3a+6.
For a3 this is15, below17.

### One fixed unsafe marker and original lifts

Marker x is nonzero at infinitely many positive even powers.
Indeed Phi=1. If its even marked series had finite support, then

 v(u)+v(-u)=A(u²)

for a polynomial A. After substituting u(x), the other formal branch
would be j(x)=A(u(x)²)-x, a rational involution with u(jx)=-u(x).
Rational-map degrees force j to be Mobius. It would fix the unique
pole x=-1 of u, and its two zeros x0 and x=-3/2, whose different
multiplicities1 and2 prevent their exchange. A Mobius transformation
fixing three distinct points is the identity, contradicting u(jx)=-u.
Hence the even marked series is not polynomial.

This also proves that the selected formal return on the original
x-line is not rational. The stable field K(u,R,T)=K(u) nevertheless
has its involution u -> -u.

An explicit original balanced lift has degree d=3a+6. With the
standard bilinears R,Z,W,T satisfying T²+2ZW=R², put

 U_h=R+Z, V_h=R+2Z/3, P_h=U_h V_h,
 C_h=V_h^(a-2)U_h^(a+3),
 B=2WR^4-T²*(P_h²-R^4)/Z,
 F=T^a C_h B.                                         (R4)

The quotient is a literal polynomial because P_h(R,0)=R².
Dehomogenization gives exactly(R1), without a scalar clearing.
At R0, B=-(4/9)T²Z³ and C_h=(2/3)^(a-2)Z^(2a+1), so F is nonzero
there: it is radial-primitive.

The prime bilinear T has exact ordera in F, since B modT=2WR^4
and C_h is not divisible by T. The prime bilinear V_h has exact
ordera-2: multiplying B by Z and reducing modulo V_h givesR^6,
which is nonzero in that domain. The bilinears have rank-two
coefficient matrices and remain irreducible over every field
extension. Since gcd(a,a-2)=1, F is not a scalar times a proper
power over any extension. Original marker Z detects infinitely
many even powers by the exact factorial/angular normalization.

There is also NO polynomial generator of K(u) inside K[x].
The map u has degree3. At infinity its value is4/9 and

 u-4/9=-(4+3x)/(9U³).

Thus the infinity fiber also contains the finite point-4/3, and
infinity has ramification index2. Every generator of K(u) differs
from u by a Mobius transformation. A polynomial generator would
have a singleton pole fiber at infinity of multiplicity3, impossible.
Consequently the cubic-height polynomial-quotient mechanism does
not extend automatically to higher starting heights.

## 2. Exact local patterns at every numerator root

At a geometric P-root alpha let p=ord_alpha P>=1,h=ord_alpha H.
The point is nonzero and J is a unit there. Then, without any
additional purity assumption,

 ord_alpha u=ap-h,
 ord_alpha omega=(a-1)p-h.                            (R5)

The rational primitive part R is regular if the second order is
nonnegative or is-1; a simple differential pole contributes only
a logarithm. Otherwise its pole order is

 h-(a-1)p-1.

Thus h<(a-1)p is exactly a finite ZERO of omega, with u having
zero order at leastp+1. Such finite zeros are permitted at higher
height and are realized by(R1). If h>ap, u has pole orderh-ap
and R has pole order(h-ap)+p-1, which is at least as large.
The retained existential law forces at least one such compensating
P-root in a nontrivial pure phase; it does not impose(R0) on all
other roots.

In the special Phi=1 subclass, omega=du. If u has nonzero orderk
at a P-root, its derivative has orderk-1. Comparing with(R5)
forces p1. Orderk0 is impossible, since then du is regular but
omega has negative order-p. Hence every P-root in this subclass
is simple and has

 h in {0,...,a-1} union {a+1,a+2,...}.

A J-root or an H-root outside P would give different u' and omega
pole orders, so neither is possible: J1, and all H-roots lie in P.
The logarithmic residues of u are

 k_alpha=a-h_alpha=1/(alpha*P'(alpha)) in Z-{0},
 1+sum_alpha k_alpha=0

when P is nonconstant. The residue at0 is1 and u has a finite
nonzero limit at infinity, giving the last identity. These are
necessary patterns in this subclass, not a claimed full existence
classification.

## 3. Admissible polynomial pullbacks at every fixed height

**Lemma SameHeightAdmissiblePullback.** Fix any integer a>=0 and a
normalized unit-center pair Cbar(z),Dbar(z). Let w in K[x] satisfy
w(0)=0,w'(0)=1. Define the rational functions C,D and the selected
formal square root X by

 C=Cbar(w)*(w/x)^(a-1)/(w')^a,
 D=Dbar(w)*(w/x)^(a+1)/(w')^(a+2),
 X=(w/(xw')) Xbar(w).                                (R6)

For a0 a negative exponent is interpreted rationally. Xbar and X
may be algebraic for a nonpure target; their selected formal branches
suffice here. In the pure odd-height applications they are rational
by the retained descent theorem. Assume the
displayed C,D are polynomials. Their center is normalized and the
selected square-root identity holds. Their exact coordinate and
Jacobian satisfy

 u=ubar(w),       T=Tbar(w),       Phi=Phibar.          (R7)

Indeed substitute(R6) into u=xX^a/C; the powers cancel to
w Xbar(w)^a/Cbar(w). Differentiating this identity and substituting
into X^(a-1)/(C u') proves the T formula. The formal inverses at0
are compatible because w'(0)=1.

Therefore purity of the source and target is equivalent at this starting height.
Indeed(R8) with marker1 identifies every positive pure moment in BOTH
directions; no safety equivalence for arbitrary source markers is inferred.
More strongly, for EVERY polynomial marker Gbar(z) and every m>=0,

 Lambda_x(Gbar(w)*p^m)=Lambda_z(Gbar*pbar^m).           (R8)

When am is odd both sides are zero. When am is even the endpoint
formula has the SAME kappa_(a,m), and(R7) identifies every marked
coefficient. The identity does not require purity of the target.
Polynomiality of C,D is an explicit premise; arbitrary w need not
be admissible.

This is a rational-differential change-of-variable construction.
The chain rule and formal endpoint identity prove it directly;
no rational-return, polynomial-descent or composition-classification
theorem is being imported.

## 4. Nonconstant Phi and a prime59 obstruction

Take w=x(1+x)², and now write

 A=1+x, B=1+3x, Q=1+3x/4, lambda=27/4,
 Z0=1+lambda*w=B²Q.

(The symbol Z0 here is a univariate polynomial, distinct from the
original bilinear Z.) Let n>=1 and a>=3 odd. Use the RETAINED
hyperbolic target at the SAME heighta:

 S=S_(2n)(lambda*w),
 Cbar=Z0^(n+a) S^(a-1)/(2n)^(a-1),
 Dbar=-Z0^(n+a+2) S^(a+1)/(2n)^(a+1).

Its selected Xbar=Z0*S/(2n), coordinate
ubar=(Z0^n-Z0^(-n))/(2n*lambda), and nonconstant Tbar=2/(1+Z0^(2n))
give all-order purity by the retained hyperbolic proof.

Since w'=AB and Z0=B²Q, the pullback(R6) is polynomial:

 X=A B Q S/(2n),
 C=A^(a-2) B^(2n+a) Q^(n+a) S^(a-1)/(2n)^(a-1),
 D=-A^a B^(2n+a+2) Q^(n+a+2) S^(a+1)/(2n)^(a+1).     (R9)

At A=0, w=0 and S=2n, while B,Q are nonzero. Therefore X has
a simple zero and C has exact multiplicitya-2, again violating(R0).
The lowered cubic pair has

 Ctilde=B^(2n+1) Q^(n+1)/A,

an actual simple pole. The source target has nonconstant Phi,
so(R7) preserves a genuinely nonconstant Phi in this counterexample.

The exact degrees are

 degC=3(2a-1)n+1, degD=3(2a+3)n+1,
 d=3(2a+3)n+a+2.                                    (R10)

For a3,n2 this is d59, with degC31,degD55 and degX12.
Prime59 already occurs in the retained hyperbolic degree spectrum
(the corresponding old index is3n); no new prime-degree occurrence
or historical-priority claim is made.

By(R8), the fixed marker w(x) has the target marker z's nonzero
mixed moment at EVERY positive even exponent. Explicitly m=2r gives

 Lambda(w p^(2r))=kappa_(a,2r)*mu_(n,lambda,r) !=0,

where mu is the retained nonzero hypergeometric coefficient. Thus
the homogeneous original marker Z(R+Z)² of fixed degree3 also
detects every even exponent, with the exact factorial factor
(2rd+4)!. Independently the all-height two-marker theorem ensures
one of the original Z,Z² detects infinitely many even exponents.
This is a purity/unsafety statement, not a finite moment inference.

The general lift formula is(R4) with the homogeneous P_h and C_h
of degreesp=degX andc=degC, replacing R^4 by R^(2p).
Its exact T-order isa. At A_h=R+Z the exact C_h-order isa-2 and
P_h vanishes, so ZB=R^(2p+2) modulo A_h is nonzero.
The same coprime-order argument proves failure of scalar proper
powers. Its nonzero top D term atR0 proves radial primitivity.
These statements apply to all allowed parameters in(R9), and the
checker controls the prime59 coefficients directly.

## 5. Sources and verification limits

The source ledger pins canonical editions of the all-height endpoint
and two-marker safety theorem, the existential compensating-pole law,
the retained hyperbolic target and the original factorial normalization.
The present counterexamples, local orders and admissible pullback
are proved explicitly.

Classical rational differential integration is the relevant context:
Rosenlicht, *Liouville's theorem on functions with elementary
integrals*, Pacific J. Math.24(1968),153--161, proposition
printedpp154--155 and the trace/norm arguments, was inspected
through the author journal source. The retained exact/log separation
proof already records that antecedent. No part of its classification
is required for the elementary identities above; no novelty of
rational differential pullbacks or logarithmic derivatives is claimed.

The standard-library verifier reconstructs the ACTUAL polynomial
coefficients, proves rational derivative and square-root identities
by clearing denominators, computes geometric multiplicities,
checks the full nonconstant-Phi source pullback and exact degrees,
and compares raw angular moments with the preserved target marker.
It separately checks original finite contraction controls and the
homogeneous divisor residues used by the written nonpower argument.
Bounded controls supplement the unrestricted proof and do not
establish purity from a finite prefix or constitute external review.


The [checker](verify_higher_start_cubic_reduction_counterexamples.py),
[inputs](higher_start_cubic_reduction_counterexamples_inputs.json),
[source ledger](higher_start_cubic_reduction_counterexamples_source_ledger.json)
and [receipt](higher_start_cubic_reduction_counterexamples_verification.json)
retain the exact polynomial and original-contraction controls.
