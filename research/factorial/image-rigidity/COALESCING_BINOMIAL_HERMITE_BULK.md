# Uniform local Hermite bulk for the actual binomial pair

## Scope and exact normalization

This note proves a LOCAL uniform stationary-branch and Gaussian-packet lemma, and an exact two-order action-separation identity. It does not assert that all zeros of an actual moment polynomial lie in the bulk tube or that the original half-line integral has only the two controlled packets. Those global hypotheses remain explicit in Section5.

Let d>=2 and a>=10000*d^2 be integers. Put epsilon=1/sqrt(a), delta=d/sqrt(a)<=1/100, D=a+d and

s=a+sqrt(a)*v, w=a^d*(-1+delta*z),
h(v)=[(1+epsilon*v)^d-1]/delta.

The common-scale phase Phi=log[s^a(w+s^d)]-s differs by a constant from

F(v,z)=epsilon^(-2)[log(1+epsilon*v)-epsilon*v]+log(z+h(v)).

The exact stationary polynomial is

P(v,z)=z*v+(v-delta)h(v)-1,

which tends to P0(v,z)=v^2+z*v-1 as delta tends to zero. No fixed-degree Gaussian limit is being used as a uniform theorem.

## 1. Two uniformly isolated stationary branches

**Lemma UniformHermiteBulkStationaryBranches.** On the complex rectangle |Re z|<=1/10, |Im z|<=3/2, the polynomial P has exactly two roots in |v|<3. They are simple analytic branches V_plus,V_minus. Each lies within3delta of the corresponding root v_plus,v_minus of P0, where the square root in v_plus/minus=(-z plus/minus sqrt(z^2+4))/2 has positive real part. In particular

Re(V_plus-V_minus)>6/5.

Proof. Write z=x+iy. Re(z^2+4)>=7/4, so Re sqrt(z^2+4)>13/10. Also |Im sqrt(z^2+4)|<=|xy|/(13/10)<=3/26. The P0 roots obey |v_plus|,|v_minus|<6/5. To prove the last bound directly, if a root has modulus r>=6/5, the identity z=1/v-v implies

[Re z/(r-1/r)]^2+[Im z/(r+1/r)]^2=1.

But its left side is at most(3/11)^2+(45/61)^2<1. The other root is the negative reciprocal, so both roots have moduli between5/6 and6/5.

For |v|<=R, Taylor's formula and |h''(v)|<=delta*exp(delta*R) give

|h(v)-v|<=delta*R^2*exp(delta*R)/2.

Therefore P-P0=(v-delta)(h-v)-delta*v has modulus at most

delta*[R+(R+delta)R^2 exp(delta*R)/2].

On each circle |v-v_plus/minus|=3delta, take R=123/100; this bound is<11delta/5, using exp(123/10000)<=10000/9877. Meanwhile |P0|>3delta*(13/10-3/100)=381delta/100. Rouche gives one root, counting multiplicity, in each disk. The disks are disjoint, hence the roots are simple and analytic in z. Their real-part separation exceeds13/10-6/100>6/5.

To exclude additional roots in the local region, use |v|=3. Since |z|<8/5, |P0|>16/5; the same remainder bound, with exp(3/100)<=100/97, is<1/5. Thus exactly two P roots lie in |v|<3. This root count does not assume that the degree-(d+1) stationary polynomial has only two roots globally.

## 2. Exact curvature and local remainder constants

**Lemma UniformHermiteBulkGaussianPacket.** At both selected branches let A=-F_vv(V,z). Then Re A>1/2 and |A|<3. On |v-V|<=1/20, F is analytic, |F_vvv|<6 and, on the horizontal local segment,

Re[F(V+u)-F(V)]<=-u^2/5.

For every real m>=9999 the local integral over that horizontal segment has relative error<1/5 from the Gaussian leading term exp(mF(V))*sqrt(2pi/(mA)). This is an assertion about the LOCAL integral, not the complete half-line period.

Proof. The exact stationary curvature is

A=[1+V^2+(epsilon-delta)V]/(1+epsilon*V)^2.

The root bounds in Section1 imply |V|<=123/100, |Re V|>=57/100 and |Im V|<=1089/1300. Thus the real part of the numerator exceeds61/100 and its modulus is at most6313/2500. Put t=123/10000. Since |epsilon*V|<=t,

|(1+epsilon*V)^(-2)-1|<=(2t+t^2)/(1-t)^2.

The rational inequalities

61/100-(6313/2500)*(2t+t^2)/(1-t)^2>1/2,
(6313/2500)/(1-t)^2<3

prove the curvature claims.

On the local disk |v|<=32/25 put H=625/617, so |h'|<=H,|h''|<=H/100,|h'''|<=H/10000. At the stationary point z+h(V)=(1+epsilon*V)^d/V. Bernoulli's inequality gives |z+h(V)|>=9877/12300. Moving by at most1/20 changes h by at most H/20. Hence throughout the disk

|z+h(v)|>=9877/12300-125/2468>3/4.

The logarithms are therefore analytic there. The elementary third-derivative bound is

|F_vvv|<= (2/100)/(1-32/2500)^3
             +(H/10000)/(3/4)+3(H/100)H/(3/4)^2+2H^3/(3/4)^3 <6.

Consequently the cubic Taylor remainder is at most|u|^3. Re A>1/2 and |u|<=1/20 give the stated decay1/5.

A direct Gaussian comparison bounds the absolute error by

[25/m+100/m*exp(-m/2000)]*exp(m ReF(V)).

The leading modulus is greater than sqrt(2/m), because pi>3 and |A|<3. Thus the relative error is at most[25+100exp(-m/2000)]/sqrt(2m), which decreases for m>=9999 and is<1/5 there. No numerical root approximation enters these estimates.

## 3. The actual radial shift stays in the same local tube

Let N>=2 be an integer in the complement of the retained large-lower-exponent theorem, namely

N>a/(2d)+3/2.

Then N,n=N-1>=9999 under the present a>=10000d^2,d>=2 hypothesis. Suppose additionally |Re z|<=1/20 and |Im z|<=29/20. For theta in[0,1] put

m_theta=N-theta, kappa_theta=m_theta/N,
w_theta=kappa_theta^(-d)*w,
z_theta=[1+w_theta/a^d]/delta.

Because d/(N-1)<2delta^2,

|z_theta-z|<=[2delta+3delta^2]/(1-2delta^2)<=203/9998<1/40.

Indeed kappa_theta^(-d)-1<=exp[d/(N-1)]-1<=2delta^2/(1-2delta^2), while |z|<3/2. Hence every z_theta stays in the outer rectangle of Sections1-2. Both ACTUAL orders, and the stationary interpolation between them, have the two controlled local branches. The scaling of the original stationary points is

S_theta,plus/minus=kappa_theta*[a+sqrt(a)*V_plus/minus(z_theta)].

It follows that Re(S_theta,plus-S_theta,minus)>(6/5)kappa_theta*sqrt(a).

## 4. An exact weighted-action drop

**Lemma ActualPairBulkActionDrop.** Set L_theta(s)=(N-theta)log g(s,w)-Ns, and let

Delta(theta)=Re[L_theta(S_theta,plus)-L_theta(S_theta,minus)].

Then

d/dtheta[Delta(theta)/(N-theta)]
 =-N*Re(S_theta,plus-S_theta,minus)/(N-theta)^2.

Consequently

Delta(1)/(N-1)<Delta(0)/N-(6/5)sqrt(a)*log[N/(N-1)],

and, since N>=6,

Delta(1)<[(N-1)/N]Delta(0)-sqrt(a).

Proof. Stationarity cancels derivatives of S_theta, so Delta' is minus the real log-g difference. Also Delta=(N-theta)Re(logg_plus-logg_minus)-N Re(S_plus-S_minus). Substitute this into Delta' and differentiate the quotient. Integrate, using the real-part separation in Section3. Finally (N-1)log[N/(N-1)]>(N-1)/N and(6/5)(N-1)/N>=1 for N>=6.

This identity treats the weighted ACTUAL lower order without assuming that exp(-s) or1/g is a slowly varying amplitude on a local Gaussian/Airy scale. The interpolation uses only local analytic branches, not nonexistent globally entire noninteger powers.

## 5. What the local theorem would imply, and the remaining global gap

**Corollary ConditionalHermiteBulkNoncollision.**

Assume one has proved for the ORIGINAL actual integrals a complete exact two-packet decomposition at both orders, with each packet within relative1/5 of its selected Gaussian leading term, and all endpoint/other-saddle contributions absorbed in those error bounds. The within-order Gaussian prefactor ratio is between1/sqrt6 andsqrt6 by Section2. At an upper-order zero, |Delta(0)|<log[(3/2)sqrt6]<log4. Section4 then gives Delta(1)<log4-sqrt(a). Therefore the lower-order packet ratio is less than16exp(-sqrt(a))<1, and the lower order cannot vanish.

This is a CONDITIONAL two-packet noncollision consequence. Sections1-4 have NOT proved that every possible upper moment zero satisfies the inner z-rectangle hypothesis, or that the complete original contour has the assumed two packets and globally small remainders. Those are the next mathematical obligations. The root-localizing curve must be established on the whole relevant w-plane, including the Airy endpoint zones and matching regions. The exact local Airy transversality is retained separately in [the coalescing-saddle note](COALESCING_BINOMIAL_SADDLE_BRIDGE.md).

No all-binomial SFC theorem, no all-gap-two theorem and no all-parameter CFU remainder claim follows from the present local results alone.


## Exact evidence

Run from a checkout:

    python3 research/factorial/image-rigidity/verify_coalescing_binomial_hermite_bulk.py --repo-root . --output /tmp/hermite_bulk_review.json

The checker and [receipt](coalescing_binomial_hermite_bulk_verification.json)
support the stated exact constants. The [shared source ledger](coalescing_binomial_local_source_ledger.json)
pins the current accounts and retained sources. The whole parameter bounds
and conditional action argument are proved above; finite controls do not
prove global contour or root localization. All review is internal, with no
external human review, formal verification or priority assertion.
Retained20November2026.
