# The exact quarter-critical residue obstruction for a binary six-window

This is a refinement of the [critical six-window theorem](BINARY_QUADRATIC_CRITICAL_DYADIC_SIX_MOMENT_WINDOWS.md). It proves a necessary residue profile for an actual six-zero block at delta=3/4. It does not construct such a native zero block.

Use the same natural-integral binary phase
 q=alpha X^2/2+beta XY+gamma Y^2/2+dX+eY,
 f=q+k, v(k)=3/4,
 Rbar=N/D, N=1+(dbar+ebar+betabar)T+dbar ebar T^2,
 D=(1+alphabar T+dbar^2T^2)(1+gammabar T+ebar^2T^2).
Assume Rbar is NONPOLYNOMIAL. Put mu=residue(k^4/8), a nonzero element of the residue field, and
                    B(T)=(1+mu T^4)N(T)/D(T).

**Proposition.** If SIX actual factorial moments at n,...,n+5 all vanish, then n is2 modulo8 and, up to swapping X and Y, the normalized residues obey
 alphabar=ebar=0,
 gammabar=dbar=betabar=lambda !=0,
 mu=lambda^4.                                               (1)
Conversely, precisely this profile makes the critical residue TEST vanish at the six positions2,...,7 modulo8:
 Rbar=1/(1+lambda T)^3,
 B(T)=1+lambda T,
 (B_0,...,B_7)=(1,lambda,0,0,0,0,0,0).                       (2)
The converse in this sentence concerns the finite residue test ONLY, not exact characteristic-zero moment zeros.

## Proof

If B is nonpolynomial, the independently proved critical six-window theorem already excludes any six-zero actual block. Thus B must be a polynomial. Since degD>=1 and degN<=2,
                         degB=4+degN-degD<=5.              (3)
The critical endpoint calculation identifies a necessary eight-periodic test vector with the coefficient prefix B_0,...,B_7. Its index0 equals1. A circular six-zero block cannot include0, hence it can start ONLY at residue1 or2.

If it starts at1, then B_1,...,B_6 vanish. By(3), B=1. The identity D=(1+mu T^4)N forces degD=4 and degN=0. But degD=4 implies dbar ebar!=0 and therefore degN=2. Contradiction.

If it starts at2, then B_2,...,B_7 vanish. By(3), degB<=1. The constant case was just excluded. Thus degB=1 and
                         degD-degN=3.                     (4)
The possibility degD=4,degN=1 is impossible by the same native degree-four implication. The only remaining possibility is
                         degD=3, N=1.                     (5)
Over an algebraic closure, 1+mu T^4=(1+lambda T)^4 for a unique lambda with lambda^4=mu. Since D divides this fourth power, (5) forces
 D=(1+lambda T)^3, B=1+lambda T.                            (6)

The native denominator has degree3 only when exactly one of dbar,ebar is nonzero and the opposite diagonal residue is nonzero. Swap variables if necessary so dbar!=0 and ebar=0. Then
 D=(1+alphabar T+dbar^2T^2)(1+gammabar T), gammabar!=0.
Both factors divide the same single-root polynomial in(6). Their constant coefficients are1, so
 1+gammabar T=1+lambda T,
 1+alphabar T+dbar^2T^2=(1+lambda T)^2.
This gives gammabar=lambda,alphabar=0,dbar^2=lambda^2. Frobenius is injective on a field, hence dbar=lambda. Finally N=1 gives dbar+betabar=0, so betabar=lambda. These equalities prove(1). In particular lambda belongs to the ORIGINAL residue field because it equals dbar; passing to an algebraic closure introduced no extra profile.

Direct substitution gives(2). The test has a six-zero circular block only at start2 modulo8 because its first two entries are nonzero. The critical endpoint moment implication now proves the stated necessary index and profile.

## Exact local translation and limits

Normalize the WHOLE phase further by its unit d, if the displayed orientation is chosen. The remaining profile says
 v(alpha)>0, v(e)>0,
 v(gamma-d)>0, v(beta-d)>0,
 v(k)=3/4, v(k^4-8d^4)>3.                                  (7)
These are valuations of the natural coefficients alpha=2A,gamma=2C,beta=B. They are not the exact identities alpha=e=0 or gamma=beta=d. No coordinate rescaling or translation preserving Lambda is silently invoked.

Thus with a nonpolynomial matching residue, every quarter-critical dense six-term phase outside (7), and every positive start outside n=2 mod8 even within (7), has a nonzero moment in its six-window. Polynomial Rbar remains a separate exception. The actual one-axis/mixed critical profile (7) still requires a higher native arithmetic argument. The test cancellation alone is not an SFC counterexample.

No finite field enumeration is a premise, no smaller-support window is inferred, and no canonical/DAG mutation, formal/external assurance or novelty claim is made.