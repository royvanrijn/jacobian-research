# Conditional closure of the half-valuation quadratic forcing branch

Let Lambda(X^d)=d!, h_j(r,z)=Lambda((rX^2+X+z)^j), and
T_j=h_j+jh_(j-1), j>=1. Normalize v_2(2)=1.

The missing hypothesis below concerns the ordinary coefficients of a MONIC
Newton factor. It is not asserted to hold. Under exactly that hypothesis,
this note excludes every half-valuation triple and therefore gives the precise
last implication needed for the universal three-forcing-window.

## 1. The exact coefficient gate and conclusion

Let M>=2 be 1 or2 modulo4 and let z be a root with v(z)=1/2 of

 Pi_M(z)=2(z+M)^3+3(z+M)^2+M=0.                       (1)

Put p=M if M is odd and p=M-1 if M is even; then p=1 modulo4. Define

 P_p(V)=2^{s_2(p)} z^{-p} T_p(zV/4,z).

Its leading coefficient is the rational unit 2^{s_2(p)}(2p)!/4^p.
Divide by that leading coefficient to make it monic.
The divided-power native factorization recalled below gives a monic factor
g_e(V) for each higher binary degree e=2^rho>=4 occurring in p. Put
S=p mod e and ell=s_2(S). If Q_e is the corresponding unreversed
ordinary divided-power factor, the normalization is EXACTLY

             g_e(V)=V^e Q_e(1/V)/Q_e(0).                         (2)

In particular the leading coefficient of g_e is1. Multiplying Q_e by a
scalar unit changes neither (2) nor the gate.

The field K=Q_2(z) is ramified quadratic. Let sigma be its nontrivial
Q_2 automorphism and write Anti(a)=(a-sigma(a))/2.

**Hypothesis H (monic ordinary coefficient gate).** For every such M,z,p,e,
writing g_e(V)=sum_(j=0)^e a_j V^j with a_e=1, one has

                   v(Anti(a_j))>=ell+2, 0<=j<e.                 (H)

This is a gate on the actual native monic factor. It is not a claim for
an arbitrary divided-power polynomial or for arbitrarily scaled factors.
It is weaker than some bounded diagnostic minima, and no such diagnostic
is a proof of(H).

**Theorem 1 (conditional half-branch exclusion).** Assume(H). For every
k>=1 there is no triple T_k(r,z)=T_(k+1)(r,z)=T_(k+2)(r,z)=0 in Q_2bar
with z satisfying(1), M=k+1, and v(z)=1/2.

The native cubic implies M=1 or2 modulo4 in this valuation branch, so
these are all the half-valuation possibilities, not a parity subfamily.

## 2. Exact local inputs, the field, and the degree-one exception

The one imported native-factor input is
[the ramified divided-power pattern](QUADRATIC_FACE_FORCING_HALF_VALUATION_DIVIDED_POWER_PATTERN.md),
SHA256 43c7815c2c91a2dac6c5cc4ce89584a8520d263a4191af3156db2181bbfc1f60, Lemmas1--2. Its exact consequences used here are these:

* If z^2/2=1 modulo2O_K, P_n has factors of every higher binary degree
  e>=4 and a final factor of degree n mod4, when n=1 or2 modulo4.
  Every root of the degree-e factor has valuation1/e.
* For n=1 modulo4, the final degree-one root has valuation1. For
  n=2 modulo4, the final degree-two roots have valuation1/2.
* For n=1 modulo4, each monic higher g_e has constant valuation1 and

 v(a_j)>=s_2(e-j),
 v(a_j)>=s_2(e-j)+1 when j=2 or3 modulo4.                         (3)

The fixed-degree divided-power Hensel proof is valid over O_K with the
ideal2O_K, even though O_K/2O_K is nonreduced. No irreducibility of g_e
over K is required: the whole degree-e Newton factor is used.

Here the hypothesis on z is automatic. For W=M(M+1)/2 (odd) and
u=z^2/2, direct rearrangement of(1) gives

 (2M+1)(3u+W)+2z(u+3W)=0,
 u+W/3=-16Wz/[3(3(2M+1)+2z)].                                  (4)

Thus v(u+W/3)=9/2, so u=1 modulo2O_K.
For completeness, after substituting w=2z the cubic has the simple
unit root1 modulo2. Its Hensel lift supplies the unique negative-valuation
Q_2 root; the other two roots constitute the quadratic factor containing z.
Write t=Tr(z), D=-Norm(z), so z^2=t z+D and v(D)=1. In

 u+W/3=(D/2+W/3)+(t/2)z,

the two summands have respectively integer and half-integer valuations.
Equation(4) therefore forces v(t)=5. Consequently
pi=z-t/2 has valuation1/2, sigma(pi)=-pi, and pi^2 belongs to Q_2 with
valuation1. Every coefficient a in K has the unique decomposition
A+B pi with A,B in Q_2; its rational and anti parts have distinct valuation
classes (integers and half-integers) when nonzero.

Suppose a half-valuation triple exists. Its native equation is Pi_M=0.
Both T_p and T_(p+1) are zero, whichever parity M has. Put beta=4r/z.
The roots of the normalized P_p are all nonzero, so r!=0. If beta came
from its final degree-one factor, it would have valuation1. But T_(p+1)=0
makes the same beta a root of the even normalized polynomial P_(p+1),
whose final factor has root valuations1/2 and whose higher factors have
root valuations1/e for e>=4. None has root valuation1. This excludes the
final degree-one factor. Thus beta lies in a unique higher factor g_e,
with e>=4 and

 v(beta)=1/e,  v(r)=-3/2+1/e.                                    (5)

This also settles p=1 immediately, since it has no higher factor. Below
one may therefore assume p>=5. Distinct root-valuation classes, not
irreducibility over the ramified field, give uniqueness of the chosen
whole Newton factor.

## 3. Native derivative identities and BOTH middle parities

The elementary factorial identities give, for j>=2,

 h_j=2rj(2j-1)h_(j-1)+j(j-1)(1-4rz)h_(j-2)+z^(j-1)(z+j),        (6)
 2r^2 partial_r h_j=(2rj-1)h_j+j(1-2rz)h_(j-1)+z^j.             (7)

For(7), set U=Lambda(X phi^(j-1)), phi=rX^2+X+z. Then
h_j-z^j=j(2rU+h_(j-1)) and
r partial_r h_j=j(h_j-U-z h_(j-1)); elimination proves(7).
One may derive(6) by applying Lambda(R')=Lambda(R)-R(0) to
phi' phi^(j-1) and phi^j, using (phi')^2=4r phi+1-4rz.
These are polynomial identities, without hidden nonzero pivots.

For adjacent lower zeros T_(j-1)=T_j=0, direct elimination in(6)--(7)
gives

 r^2(2j-1+2z) partial_r T_j=z^(j-1)(z+j)^2.                      (8)

For adjacent upper zeros T_p=T_(p+1)=0, put D_+=2p+1+2z. The exact
other orientation is

 r^2 partial_r T_p = z^p [(p+z)(p+z+1)
                +rD_+(p(p+1)-z^2)]/[(4rz-1)(p+1)D_+].           (9)

Here is its full elimination, to make the orientation and sign reviewable.
At those two zeros, (6) gives

 2rp(p+1)D_+ h_(p-1)=z^p(z+p+1).

At T_p=0, using(6) to replace T_(p-1) in the derivative expression gives

 r^2(4rz-1)partial_r T_p
   =rp(D_+-1-2rzD_+)h_(p-1)+r z^p(z+p).

Substituting the first equation, noting D_+-1=2(p+z), gives precisely(9).
After(5), all its denominators are nonzero: D_+ is a unit,
4rz-1 is a unit since v(4rz)=1+1/e>0, p+1 has valuation1, and r!=0.
Thus this use of rational expressions loses no hypothetical root.

If M is odd, p=M and the triple supplies the lower zeros T_(p-1),T_p.
In(8), the denominator 2p-1+2z is a unit and p+z is a unit. Hence

                 v(r^2 partial_r T_p)=(p-1)/2.                  (10)

If M is even, p=M-1 and the triple instead supplies the upper pair
T_p,T_(p+1), so(9) is used. Equation(4) implies

 E=z^2+M(M+1)/3,  v(E)=11/2,
 p(p+1)-z^2=2M(2M-1)/3-E,
 v(p(p+1)-z^2)=2.

The first numerator term in(9) has valuation1/2. The second has valuation
v(r)+2=1/2+1/e, strictly larger. The denominator has valuation1, entirely
from p+1=M. Thus(10) holds for this parity as well. In particular the
derivative is nonzero. This is the part that cannot be replaced by the
lower-pair formula at a nonzero preceding T_(p-1).

Scaling the derivative to P_p, using(5), now gives in BOTH cases

 v(P_p'(beta))=s_2(p)+1-2/e.                                     (11)

Indeed P_p'=2^{s_2(p)-2}z^{1-p}partial_r T_p and
(10) gives v(partial_r T_p)=(p-1)/2+3-2/e.
The rational leading-unit normalization of P_p does not change(11).

Let h be the number of binary bits of p above e. The lower factors have
total degree S and their values at beta have total valuation S/e;
each of the h higher factors has value of valuation1. This follows
termwise from their root valuations: if d<e each difference has valuation
1/e, and if d>e each difference has valuation1/d. Therefore

 v(g_e'(beta))=ell+2-(S+2)/e=:tau,                               (12)

using s_2(p)=ell+h+1. No cancellation among other-factor values is assumed.

## 4. The monic anti gate removes the ramified cancellation problem

Write a_j=A_j+B_j pi, A_j,B_j in Q_2. Hypothesis(H) says
v(B_j pi)>=ell+2. The native bounds(3) also pass to A_j: integer-valued
A_j and half-integer-valued B_j pi cannot cancel, so

 v(A_j)>=s_2(e-j),
 v(A_j)>=s_2(e-j)+1 for j=2 or3 modulo4.                         (13)

Put R(V)=sum A_j V^j and A(V)=sum B_j pi V^j, so g_e=R+A and R is monic.
Every nonzero anti derivative summand at beta has valuation

 v(j B_j pi beta^(j-1))>=ell+2+v_2(j)+(j-1)/e>=ell+2>tau.

Their sum also has valuation greater than tau. By(12), R'(beta) must
therefore have valuation exactly tau.

The nonzero rational derivative summands j A_j beta^(j-1) have pairwise
distinct fractional parts (j-1)/e modulo1. They have a unique least
valuation; it cannot cancel. Since S=p mod e is1 modulo4,
1<=S<=e-3. The only derivative index with the fractional part in(12) is

                 j_0=e-S-1,  j_0=2 modulo4.

It has v_2(j_0)=1 and (j_0-1)/e=1-(S+2)/e. Thus equality with tau forces

                 v(A_(j_0))=ell.                               (14)

But s_2(S+1)=s_2(S)=ell because S is1 modulo4. The extra bound in(13)
at j_0=2 modulo4 gives v(A_(j_0))>=ell+1, contradicting(14).
Every possible higher factor is excluded, proving Theorem1.

This argument bounds the ENTIRE anti derivative above tau before using
rational fractional uniqueness. It never assumes that two low derivative
summands with half-integer coefficient valuations cannot cancel.

## 5. Global and symmetric consequences remain CONDITIONAL

[The exact residual-family theorem](QUADRATIC_FACE_FORCING_THREE_ZERO_RESIDUAL_FRONTIER.md),
SHA256 736e2675513962010892597610563b2d00a1e8047e20e322efd8c7d665404be2,
Theorem1, proves that any characteristic-zero triple must have
M=q(4q+3)^2 and z=-8q(q+1)(2q+1)+/-sqrt(-3q(q+1)), with q=1 or2 modulo4.
Every embedding has v_2(z)=1/2. Its algebraicity argument embeds the
number field generated by a hypothetical r,z into Q_2bar; it does not
assume that the whole original coefficient field embeds there.
Consequently(H), together with that proved residual classification,
would rule out every three-zero forcing block over every characteristic-zero
field. This is a one-way consequence of the unproved coefficient gate.

As a separate application, [the full symmetric quadratic recurrence](FULL_SYMMETRIC_BINARY_QUADRATIC_STRONG_FACTORIAL.md), Section2 Equation10,
SHA256 c29933ff9bfee00d55703e049cde90c7ceb1cc577ada03a34dce4b9413530c00,
would then exclude five consecutive moments of
f=a(x^2+y^2)+bxy+e(x+y)+c with ec!=0: a five-zero block at n forces
T_(n+2),T_(n+3),T_(n+4)=0 after division by e.
In particular the actual support-five b=0,aec!=0 family would satisfy
its Strong Factorial window CONDITIONALLY on(H). Neither this note nor
the existing six-window theorem proves(H), that five-window assertion,
full binary quadratic SFC, SFC(2), or SFC(3).

## 6. Diagnostic normalization and status

The private diagnostic scripts reporting factor minima first formed
normg_k=G_PD,k/G_PD,0 using the full quadratic-field inverse, and then
weighted by2^k/k!. Thus their ordinary coefficient at V^(e-k) was
(2^k/k!)normg_k, exactly the monic normalization(2).
For a coordinate a+bz, the reported anti valuation is v_2(b)+1/2 because
Anti(a+bz)=b pi. No arbitrary scalar unit was left in those reported
minima. This confirms the meaning of those diagnostics; it does not
promote their bounded observations into a proof of(H).

The [native high-block congruence](QUADRATIC_ODD_NATIVE_HIGH_BLOCK_ANTI_CONGRUENCE.md)
and the proposed formal-quotient estimate are separate statements; neither
is an input to the conditional implication proved here. A proof that
transfers the necessary weighted bound to the actual monic g_e is still
required. Hypothesis(H) is not proved by the native pattern, the native
high-block theorem, or bounded factor diagnostics. There is no converse
from the universal three-window to(H) asserted in this proof.
