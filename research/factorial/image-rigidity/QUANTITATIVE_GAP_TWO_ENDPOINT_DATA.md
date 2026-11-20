# Quantitative matched-cubic data at the gap-two endpoint

This note proves explicit bounds for the matched critical-value and stationary-amplitude data. It does not prove a full inverse or divided-amplitude bound uniform on a nonzero parameter window, or a numerical Strong Factorial threshold. In particular it does not establish a_0<=4*10^16. A [separate basepoint lemma](HERMITE_CRITICAL_CUBIC_INVERSE_DISK.md) proves a full inverse disk only at delta=0,z=2i; Uniform inverse persistence is not established by this note; the full thin window is treated in [the uniform cubic inverse theorem](GAP_TWO_FULL_ENDPOINT_UNIFORM_CUBIC_INVERSE.md). That theorem does not assert the whole broader complex disk of the matched-data bounds.

Let

F_delta(v,z)=4*delta^(-2)*[log(1+delta*v/2)-delta*v/2]
             +log(z+v+delta*v^2/4),

with its removable delta0 value. The matched cubic normalization has phase Phi+t^3/3-zeta*t and limiting v_t=i, zeta_z=-i at(v,z,delta)=(-i,2i,0). Denote the even/odd interpolation of the two stationary Jacobians by A,B, so the leading canonical integral is A*Ai(xi)-m^(-1/3)*B*Ai-prime(xi), up to one common nonzero factor.

**Lemma ExplicitMatchedEndpointData.** These symmetric normalization data extend holomorphically through the turning curve on

|delta|<=10^(-3), |z-2i|<=11/100.

On the smaller polydisc |delta|<=10^(-8), |z-2i|<=101/1000 they satisfy

|A|>49/50, |B|<13/25, |zeta|<1/8, |zeta_z+i|<1/4.          (Q1)

The statement concerns the matched data. Where a regular cubic coordinate exists, these are exactly its leading interpolation coefficients and unfolding parameter. No off-stationary inverse-Jacobian bound is inferred from(Q1).

## 1. Exact Hermite formulas

At delta0 write the stationary pair as

v_1=-i*exp(-u), v_2=-i*exp(u), z=2i*cosh(u).

Define b by

b^3=(3/4)*(sinh(2u)-2u), b/u tends to1.

The two critical actions give zeta=b^2. The stationary Jacobian relation F_vv*v_t^2=2t, with v_1 corresponding to t=b, yields

A^2=-(b/2)*coth(u/2),
B^2=-tanh(u/2)/(2b), AB=-1/2,
B/A=tanh(u/2)/b.                                        (Q2)

These extend regularly through u0, with A=i and B=i/2.

Every z in |z-2i|<=101/1000 has its selected pair parametrized by |u|<1/3. On |u|=1/3, Rouche compares u^2/2 with (z-2i)/(2i) and the remaining cosh series:

101/2000+[cosh(1/3)-1-1/18]<1/18.

Thus cosh(u)-1=(z-2i)/(2i) has exactly two roots, opposite in sign, in that disk.

Put S(u)=(3/(4u^3))*(sinh(2u)-2u), so b/u=S(u)^(1/3). The positive Taylor coefficients show, on |u|<=1/3,

|S-1|<1/40,
|b/u-1|<1/117,
|tanh(u/2)/(u/2)-1|<1/100,
|sinh(u)/u-1|<1/50.                                     (Q3)

For the cube-root inequality integrate its derivative along the segment from1 to S; using (40/39)^(2/3)<40/39 gives1/117. For the tanh inequality, with r=1/6, bound its numerator directly by cosh(r)-sinh(r)/r and its denominator below by2-cosh(r).

Equations(Q2)--(Q3) imply

|A^2+1|<= (1/117+1/100)/(1-1/100)<1/50.

Consequently |A|>989/1000, |B|<500/989 and

|zeta|< (1/9)*(118/117)^2.

Differentiating the critical-value formula gives

zeta_z=-i*sinh(u)/b,

so |zeta_z+i|<(1/50+1/117)/(1-1/117)<3/100.                (Q4)

## 2. A complex outer polydisc for the actual family

The exact stationary polynomial is

P_delta=v^2+z*v-1+delta*(v^3/4-v)-delta^2*v^2/4.

On |v+i|=1/2 compare it with(v+i)^2. For |z-2i|<=11/100 and |delta|<=1/1000, its other terms have total modulus at most

(11/100)*(3/2)
 +(1/1000)*[(3/2)^3/4+3/2]+(1/1000000)*(3/2)^2/4<1/4.

Thus exactly two selected roots v_1,v_2 remain in |v+i|<1/2. Their symmetric sums and product are holomorphic in delta,z, including delta0, by the contour-trace formulas. Write c=(v_1+v_2)/2, q=(v_1-v_2)/2. The exact factorization gives

F_v=[(v-c)^2-q^2]*K(v),
K(v)=-L(v)/[(1+delta*v/2)*(z+v+delta*v^2/4)],
L(v)=1-delta^2/4+delta*(v+v_1+v_2)/4.                    (Q5)

For every v in the same physical disk,

|z+v+delta*v^2/4-i|<=9769/16000,
|L-1|<=4501/4000000,
|L/(1+delta*v/2)-1|<1/500.

The reciprocal-disk bound therefore gives the literal inequalities

Re(K/i)>=1/(1+9769/16000)
          -(1/500)/(1-9769/16000)>3/5,
|K|< (501/500)/(1-9769/16000)<13/5.                       (Q6)

Also |K'/K|<3, directly from

K'/K=L'/L-(delta/2)/(1+delta*v/2)
              -(1+delta*v/2)/(z+v+delta*v^2/4).

## 3. Matched data without a full inverse-map estimate

Define the normalized weighted average

J=(3/4)*integral_-1^1 (1-t^2)*K(c+q*t) dt.

It is holomorphic in c,q^2 and the parameters; the odd powers of q integrate to zero. The weight is positive and has total integral1. Hence Re(J/i)>3/5 and |J|<13/5 by(Q6). Let mu=J/i and take its principal fractional powers in the right half-plane. Put

b=q*(-i)*mu^(1/3), zeta=b^2=-q^2*mu^(2/3).

The exact integral of F_v between the two stationary points is -4*b^3/3, so this is precisely the matched critical-value unfolding. It satisfies |zeta|<49/100, since |q|<1/2 and |mu|^(1/3)<7/5.

For bounding only the stationary Jacobians, introduce the AUXILIARY analytic function

H(v)=i*mu^(1/6)*(K(v)/i)^(-1/2).

At v_1,v_2 its square is b/[q*K(v_j)], exactly the stationary cubic Jacobian square. The branch is fixed by H=i at the Hermite merger. Define

A=(H(v_1)+H(v_2))/2,
B=(H(v_1)-H(v_2))/(2b).

The divided difference is regular at q0. Symmetry makes A,B holomorphic in delta,z across the turning curve. On the outer polydisc,

|H|^2< (7/5)/(3/5)=7/3<4,
|H'|=(|H|/2)*|K'/K|<3,
|q/b|<5/4.

Thus |A|<2 and |B|<15/4<4. Importantly H(v) is NOT claimed to equal the inverse Jacobian away from the two stationary points.

## 4. Cauchy transfer to delta<=10^(-8)

Apply the Cauchy coefficient bound in the complex delta disk of radius1/1000, with fixed |z-2i|<=101/1000. If a holomorphic function has bound M there, its change from delta0 is at most

M*|delta|/(1/1000-|delta|)<=M/99999.

Using M=2,4 and49/100 for A,B,zeta gives respectively changes below21*10^(-6),41*10^(-6), and1/200000. Combining these with(Q3)--(Q4) proves the first three inequalities in(Q1).

For zeta_z, use the additional z margin from101/1000 to11/100. Cauchy's derivative bound gives |zeta_z|<=490/9 on the outer delta disk, so its delta change is at most490/(9*99999)<3/5000. Together with(Q4) this proves |zeta_z+i|<3/100+3/5000<1/4.

## 5. A quantitative Robin budget, conditional on the actual remainder bound

**Proposition ConditionalRobinBudget.**

The effective zero estimate should retain its m dependence. If a genuine complete-period representation satisfies

P_m=A*Ai(xi)+m^(-1/3)*B_*Ai-prime(xi)+R_m,
|R_m|<=C_R/m*(|Ai(xi)|+|Ai-prime(xi)|),

then at a zero

|Ai/Ai-prime|<=K(m)*m^(-1/3),
K(m)=[B_0+C_R*m^(-2/3)]/[c_0-C_R/m].                    (Q7)

For example, the much weaker geometric bounds c_0>=9/10, B_0<=7/10, |zeta|<=1/8, together with C_R<=10^10 and m>=10^16, imply K(m)<11/10. Indeed m^(2/3)>=4*10^10 and C_R/m<=10^(-6). The Robin product guard K(m)^2*|zeta|<1/4 and cubic guard K(m)^3/m<1/16 then have large exact margins. At a=4*10^16 the actual radial drift overwhelms both zero bands.

This is a CONDITIONAL budget. The present note does not supply C_R<=10^10, a regular inverse map on a fixed cubic disk throughout the stated endpoint window, or all nonlocal/vector comparison constants needed for that value. Those are the precise remaining numerical obligations. Qualitative compactness from the retained theorem cannot be substituted for these quantitative bounds.

## 6. An exact negative control for the remaining gap

The auxiliary H is not the full Jacobian even at the Hermite merger. With v=-i+u and z=2i, the canonical coordinate has

t=-i*u+u^2/4+(11i/80)*u^3+O(u^4).

Hence the true inverse Jacobian, expressed as a function of v, is

dv/dt=i+u/2+(13i/80)*u^2+O(u^3).

The auxiliary function is H(v)=i*sqrt(1-i*u), whose expansion is i+u/2+(i/8)*u^2+O(u^3). The difference is (3i/80)*u^2, already nonzero at quadratic order. This rules out turning the bounds in Section3 into an unjustified full inverse-map or divided-amplitude bound.


## Exact reproduction and classical context

    python3 research/factorial/image-rigidity/verify_quantitative_gap_two_endpoint_data.py \
      --repo-root . --output /tmp/quantitative_gap_two_endpoint_review.json

The standard-library verifier checks the rational series and polydisc
margins, universal normalization identities, true and auxiliary Taylor
series, three exact critical-inverse bounds, and the explicitly
conditional Robin calculation. Its source ledger requires all SIX
listed proof/checker/context inputs.

The cached primary source
[Olde Daalhuis and Temme, Uniform Airy-type expansions of integrals](sources/olde-daalhuis-temme-1994-uniform-airy.pdf),
SIAM Journal on Mathematical Analysis25(1994),304--321, supplies
classical method context. Its Theorem4.1 and Cauchy formulas on
printed308--309 bound remainder amplitudes only after a bound on the
full transformed analytic amplitude is supplied. Theorem5.2 on
printed312 has its own stated large-saddle-parameter assumptions.
Neither result is imported as a certificate for the missing actual
C_R bound or for uniform inverse persistence here.

Updated 20 November 2026.
