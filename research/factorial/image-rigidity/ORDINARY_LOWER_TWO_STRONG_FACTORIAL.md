# Ordinary lower-two Strong Factorial theorem

## 1. Statement and exact scope

**Theorem OrdinaryLowerTwoStrongFactorial.**

Let K be any characteristic-zero field and let L_K:K[X]->K be the factorial functional L_K(X^j)=j!. For every d>=1, N>=2 and A,B in K, the two equations

L_K((X^2(A+BX^d))^(N-1))=L_K((X^2(A+BX^d))^N)=0

force A=B=0. Thus the ordinary Strong Factorial condition holds for every binomial whose least exponent is two. This statement concerns ordinary factorial moments. It does not assert the corresponding theorem for every lower exponent or every positive Gamma parameter.

The proof has three parts: known small gaps, an explicit weighted two-saddle tail, and a finite computation of actual consecutive moments. **Theorem LowerTwoPracticalWeightedTail.**

The actual normalized consecutive moment polynomials for X^2(A+BX^d) are coprime over Q whenever

(T) d>=6, D=d+2, N>=max(4000,128D).

The retained sparse theorem overlaps this tail for every d>=46. The remaining finite domain is exactly

(F) 6<=d<=45, 2<=N<max(4000,128(d+2)).

It contains 176816 pairs, with N<=6015. Both complete computations described in Section 8 have finished successfully. Analytic inequalities are certified by a separate exact rational ledger; modular finite computations certify only their stated finite domain.

## 2. Actual moment integrals and packets

Assume AB!=0; a missing monomial has a nonzero factorial moment and is immediate. Put D=d+2, H=D^D exp(-D), W=D^D exp(-d)/4 and c=5D/8. Write A/B=N^d w and g(s)=s^2(w+s^d). On a local logarithm set Phi(s)=log(g(s))-s. Define entire polynomial integrals

I_N(w)=int_0^infinity g(s)^N exp(-Ns) ds,
J_N(w)=int_0^infinity g(s)^(N-1) exp(-Ns) ds.

Direct substitution X=Ns gives the ACTUAL consecutive moments

L(f^N)=B^N N^(DN+1) I_N(w),
L(f^(N-1))=B^(N-1) N^(D(N-1)+1) J_N(w).

The same N scale occurs in both integrals. Split each at c into its upper and lower packets I_U,I_L,J_U,J_L. All integrands used to move contours are polynomials times exp(-Ns); no contour crosses an integrand pole.

Upper estimates hold on |w|<=6W/5, and lower estimates on |w|>=4W/5. The stationary equation and curvature at a stationary point are

w(s-2)=s^d(D-s),
a(s)=-Phi''(s)=[s^2-(D+1)s+2D]/s^2.

There are analytic selected stationary points satisfying

|s_U-D|<1/2, |s_L-2|<1/D.

They are unique in these disks. For the upper disk Rouche's theorem follows, using e>163/60, from

(3/5)D^2(D-3/2)[60D/(163(D-1/2))]^(D-2)<1.

At D=8 this is 85761906966528/93776847890045<1. The ratio of successive bounds is at most (9/8)^2(15/13)(64/163)=1215/2119<1: the base and both displayed polynomial growth ratios decrease. For the lower disk the maximum of |s|^d|D-s| on its boundary is at its right endpoint. Differentiating its squared modulus with respect to the boundary cosine gives the positive lower bound 2d^2-4+D^-2-4(d+1)/D. Consequently the Rouche ratio is at most

5(1-2/D-1/D^2)[(87/32)(2D+1)/D^2]^d.

Its D=8 value is (235/64)(1479/2048)^6<1. Its successive quotient is at most (64/47)(87/32)(19/81)<1. These arguments also fix the branches at w=0 and w=infinity respectively.

The curvature disks satisfy

Re a_U>1/(2D), |a_U|<3/D;
Re a_L>1/4, |a_L|<1.

For completeness map the disk |s-c0|<=r by v=1/s to the disk with center c0/(c0^2-r^2) and radius r/(c0^2-r^2). The real part of 1-(D+1)v+2Dv^2 is a quadratic in the boundary cosine. For the upper disk it decreases, and the minimum is (D/2+3/4)/(D-1/2)^2>1/(2D). For the lower disk it increases, and the minimum is at s=2+1/D and is >1/4. The derivative signs reduce respectively to

-(D+1)+(4D^2+4D)/(D^2-1/4)<0,
-(D+1)+8D(1-1/D)/(4-1/D^2)>0.

The norm bounds follow directly from

a-1/D=(s-D)[(D-1)s-2D]/(D s^2),
a-1/2=(s-2)(s/2-D)/s^2.

Write S_U=Re Phi(s_U), S_L=Re Phi(s_L), x_U=Re s_U and x_L=Re s_L. Then

|S_U-log H|<7/80,
|S_L-log|w|-(log4-2)|<1/32.

Indeed radial integration of the stationary value has derivative 1/(w+s_U^d); the upper q bound at D=8 is 6144/78125, giving error <=6144/71981<7/80. For the lower point the Taylor error of 2log(s)-s is <=1/240 and |s_L^d/w|<1/(Dd-1)<=1/47; hence the log error is <1/46 and their sum is <1/32.

## 3. Local Gaussian estimates, including the weighted integral

The selected square root of a is nonzero on each packet domain. Put n=N-1. The leading terms are

L_I*=exp(N Phi(s_*)) sqrt(2pi/(N a_*)),
L_J*=exp(n Phi(s_*)-s_*) sqrt(2pi/(n a_*)).

In particular L_J*/L_I*=sqrt(N/n)/g(s_*). The SAME square-root factor appears at both saddles. This identity is essential for the pair determinant.

Here are local Taylor constants (rho,B,b), meaning a local half-length rho, remainder |Phi(s_*+u)-Phi(s_*)+a_*u^2/2|<=B|u|^3, and quadratic decay at least b u^2:

upper: rho=1/2, B=3/D^2, b=1/(7D);
lower: rho=1/8, B=1/3, b=1/12.

We give explicit uniform checks, rather than an asymptotic assertion. On the upper local disk |s-D|<=1 set q_D=(3/10)D^2[3D/(8(D-1))]^d and y=q_D/(1-q_D). At D=8 q_D=69984/588245, and

2[D/(D-1)]^3+D^2(y+3y^2+2y^3)<18.

This bounds D^2|Phi'''|. The base decreases; q_(D+1)/q_D<=(3/8)(9/8)^3 and D^2 y^j decreases for j=1,2,3. The decay inequality reduces after clearing positive denominators to

6(D-8)^3+89(D-8)^2+354(D-8)+187>0.

On the lower local disk |s-2|<=1/4 use q_D=5D^-2[(87/32)(9/4)/D]^d. The third derivative is at most

(4/7)^3[4+d^3(y+3y^2+2y^3)]<2.

Here q_(D+1)/q_D<=(783/128)/18 and d^3q_D decreases since (7/6)^3(783/128)/18<1. The factor (D/(D+1))^(D-2)<1/2 used in this step follows at D=8 from (9/8)^6>2 and then from its increasing logarithm. The same bound gives the lower decay constant.

An elementary real Gaussian comparison yields absolute local error, for exponent m and constant amplitude, at most

[B/(m b^2)+exp(-m b rho^2)/(m b rho)] exp(m S).

To verify this bound integrate mB|u|^3 exp(-mbu^2) over the line and complete the Gaussian beyond rho. For J the normalized amplitude is exactly exp(-u), since g^(N-1)exp(-Ns)=exp(n Phi(s))exp(-s). Thus multiply the first local term by exp(rho), add exp(rho)/(m b), and keep the Gaussian completion. This avoids dividing a nonlocal integral by g. We use exp(1/2)<5/3 and exp(1/8)<8/7.

The I leading modulus lies between 1 and 5 times sqrt(D/N) at the upper saddle, and between 2 and 6 times 1/sqrt(N) at the lower one. These follow from the curvature bounds and 3<pi<22/7.

## 4. Complex contour certificates and their uniformity

For the upper contour keep the vertical segment at c and move the half-line c..infinity to height Im(s_U). For the lower contour keep 0..1 real, move 1..c to height Im(s_L), and use the two vertical segments. Every deformation is valid by the entire integrand and exponential decay at positive infinity. Logs in estimates are local notation for log modulus; a possible zero of g away from a saddle does not invalidate the contour.

Three exact Bernstein certificates establish the parts where a coarse real-envelope estimate would be insufficient. If P(l+(h-l)t)=sum a_i t^i has degree m, its Bernstein coefficients are b_k=sum_(i<=k) a_i binom(k,i)/binom(m,i). Positivity of every b_k proves positivity on the ENTIRE interval, not only on sampled points.

Put A/B=(3/10)8^8(60/163)^6=234827336908800000/18755369578009.

(UC) On 7D/8<=x<=5D/4, |Im(s)|<=1/2, Re Phi''<-1/(4D). After taking x=beta D and bounding r_D=(3/10)D^2[60/(163 beta)]^d, the D=8 certificate on [7,10] is

P_U(x)=-5B^2 x^14+1248B^2 x^12+10AB x^8-9216AB x^6-5A^2 x^2+2208A^2.

All fifteen degree-fourteen Bernstein coefficients are positive. To check the reduction for every D, r_(D+1)/r_D<4/7, the ratios of d(d+1)/D and d^2/D at successive D are at most 32/27 and 49/36. Thus both correction terms d(d+1)y/D and d^2y^2/D decrease. The base term satisfies Re(D/s^2)>(39/40)D/x^2 (the weaker 39/40 is deliberate). The resulting inequality is precisely P_U>0 after clearing positive denominators.

(UL) On 5D/8<=x<=7D/8, |Im(s)|<=1/2, Re Phi'>1/100. For |q|<=r<1 the exact disk-Mobius minimum gives

Re[(1/s)/(1+q)] >= [Re(1/s)-r/|s|]/(1-r^2).

Using Re(1/s)>=1/x-1/(4x^3), 1/|s|<=1/x gives

Re Phi'>=D/x-1-d r/[x(1+r)]-D/[4x^3(1-r^2)].

For beta in [5/8,7/8], r_(D+1)/r_D<3/4. The quantity (1-2/D)r/(1+r) decreases, with ratio at most (28/27)(6/7)=8/9; the last displayed error also decreases. The D=8 certificate on [5,7] is

P_left(x)=-101B^2 x^15+800B^2 x^14-200B^2 x^12-600AB x^8+101A^2 x^3-200A^2 x^2.

All sixteen Bernstein coefficients are positive. This is a COMPLEX-strip assertion; merely maximizing w on the positive real axis would not prove it.

(LR) On 5/2<=x<=7/2, |Im(s)|<=1/D, Re Phi'<-1/20. The corresponding Mobius maximum gives

Re Phi'<=2/x-1+d r/[x(1+r)]+d r^2/[D^2 x^3(1-r^2)].

Here r_D=(100/99)5D^-2[(87/32)x/D]^d; the factor 100/99 absorbs |s|^d/x^d<=exp(1/(2Dx^2))<100/99. The ratio r_(D+1)/r_D<17/32, and r_8<1/4. The two correction terms decrease since

(7/6)(17/32)(5/4)/(1+17/128)<1,
(7/6)(17/32)^2<1.

Set a/b=6022586125125/49539595901075456. The D=8 certificate on [5/2,7/2] is

P_right(x)=-152a^2 x^15+1280a^2 x^14-15a^2 x^12-960ab x^8+152b^2 x^3-320b^2 x^2.

All sixteen Bernstein coefficients are positive.

The remaining contour bounds follow by elementary monotone envelopes, recorded explicitly here. On the upper central interval UC integrates twice to Re Phi<=S_U-u^2/(8D). The left endpoint 7D/8 loses at least D/5120, and the right endpoint 5D/4 loses at least D/1280. On x>=5D/4, |q|< (3/10)D^2(3/10)^d<1/50 and Re Phi'<-1/10. On the vertical segment at c, the modulus relative to H is bounded by

[15/32+3/8] exp(1/25)=(27/32)exp(1/25).

The high and low contributions and their imaginary corrections decrease with D. Since

(27/32)exp(1/25+7/80+1/50)<1,

the vertical loses 1/50 relative to S_U. For example the inequality follows from exp(59/400)<=400/341. The integrated amplitude is at most 1/2.

For the lower central interval [3/2,5/2], Re Phi''<-3/20. Use |s|<251/100 and q_D=5D^-2[(87/32)(251/100)/D]^d. With y=q_8/(1-q_8),

39/20-30y-36y^2>(3/20)(251/100)^2=189003/200000.

All correction terms decrease with D, by the same ratio argument as in the lower local disk. Thus Re Phi<=S_L-3u^2/40. Its right endpoint loses at least 9/1280 (a weaker bound than the direct integrated curvature).

Keep [0,1] real. On [1,3/2] and its vertical endpoint, |s|^2<=145/64 and |s^d/w|<1/500. Their envelope divided by 4exp(-2) is at most (145/256)(501/500)exp(1/2). Hence

(145/256)(501/500)exp(1/2+1/32+1/100)<1,

using exp(11/20)<7/4. The envelope on [0,1] is smaller by monotonicity, and the total length/amplitude bound is 3. These pieces lose 1/100 relative to S_L.

For the lower far part [7/2,c], including its right vertical segment, use the real envelope

E_D(x)=(x/2)^2 exp(2-x)+(5/4)(x/D)^D exp(D-x).

For fixed x<D the high term decreases with D, while the low term decreases with x>=2. Divide [7/2,5] into [7/2,4], [4,9/2], [9/2,19/4], [19/4,5]; on each use the low term at its left endpoint and the high term at its right endpoint, at D=8. On [5,c] use the low term at 5 and the high term at c at D=8. Each of the five resulting exact exponential expressions E satisfies

E exp(1/32+1/196+1/200)<1.

The imaginary correction is at most exp(1/196), because D(Im s_L)^2/(2x^2)<=1/(2D x^2). The whole far contour has length <2D and loses 1/200. The certificate ledger checks these five expressions with rational Taylor bounds, not floating-point rounding.

## 5. Complete packet errors

The following explicit relative-error bounds include every segment, in the order local phase/amplitude, local Gaussian completion, central tails, left pieces, vertical/far pieces. They are stated separately for the two ACTUAL exponents. Put n=N-1. For N>=max(4000,128D),

E_IU =147/sqrt(ND)+14sqrt(D/N)exp(-N/(28D))+16sqrt(D/N)exp(-N/(32D))
       +100/sqrt(ND)exp(-ND/5120)+(1/2)sqrt(N/D)exp(-N/50)
       +10/sqrt(ND)exp(-ND/1280) <5/6;

E_JU =exp(1/2)[147/sqrt(nD)+7sqrt(D/n)]+14sqrt(D/n)exp(-n/(28D))
       +32sqrt(D/N)exp(-N/(64D))+200/sqrt(ND)exp(-ND/10240)
       +(1/2)sqrt(N/D)exp(-N/100)+10/sqrt(ND)exp(-ND/2560) <3;

E_IL =24/sqrt(N)+48/sqrt(N)exp(-N/768)+160/(3sqrt(N))exp(-3N/2560)
       +(3/2)sqrt(N)exp(-N/100)+10/sqrt(N)exp(-9N/1280)
       +D sqrt(N)exp(-N/200) <2/5;

E_JL =240/(7sqrt(n))+48/sqrt(n)exp(-n/768)+320/(3sqrt(N))exp(-3N/5120)
       +(3/2)sqrt(N)exp(17/8+1/100-N/100)+10/sqrt(n)exp(-9n/1280)
       +D sqrt(N)exp(-N/200) <18/25.

Here is the derivation of the weighted nonlocal terms, which cannot be obtained simply by dividing an integral by g. For the upper central tail, n u^2/(8D)-|u|>=N u^2/(16D) when N>=80D+2; direct integration outside |u|>=1/2 gives absolute bound (32D/N)exp(nS_U-x_U-N/(64D)). On the upper left the derivative of n Re Phi-x is >=n/100-1>=N/200. Since x_U-7D/8<=3D/16 and N>=1922, the endpoint exponential is <=exp(nS_U-x_U-ND/10240). On the upper vertical x_U-c<=7D/16, and N>=50D gives loss N/100. On the upper far part x_U-5D/4<0, so the N-1 exponent gives at least ND/2560 loss.

For the lower central weighted tail, n(3u^2/40)-|u|>=3N u^2/80 for |u|>=1/8 and N>=216; integration gives (640/(3N))exp(nS_L-x_L-3N/5120). The weighted left pieces introduce at most exp(17/8+1/100); the right strip has x_L<5/2 and so introduces no adverse factor. On the far part x_L<=17/8<7/2, its extra amplitude absorbs the replacement n by N. The displayed errors follow by division by the leading modulus bounds from Section 3.

For exact budgeting, n>999N/1000. In E_IU, for D<=17 the local term is <33/40 and its two Gaussian tails are individually <1/1000. For D>=18 the local term is <29/40, while those tails are <3/200 and <3/100. Uniformly its last three terms are <1/500,1/10000,1/10000. Both sums are <5/6.

For E_JU the first bracket is less than (5/3)(1000/999)(33/40+5/8)<5/2. Its remaining terms are respectively <1/50,2/5,1/20,1/1000,1/1000, whose sum with 5/2 is <3. The six E_IL terms are respectively <38/100,5/1000,1/100,1/1000,1/1000,1/1000; the six E_JL terms are respectively <545/1000,5/1000,165/1000,1/1000,1/1000,1/1000. Their sums are <2/5 and <18/25. These are strict rational certificates. For terms involving sqrt(N)exp(-cN), the derivative is negative on N>=4000. For the last lower term use D<=N/128 and decreasing N^(3/2)exp(-N/200). The two upper Gaussian tails depend only on N/D>=128. Thus no finite sampling of unbounded D or N is used.

## 6. No common zero of the actual pair

The bounds E_IU<5/6 and E_IL<2/5 first locate every zero of I in the annulus 4W/5<|w|<6W/5. On |w|=4W/5 the upper stationary height exceeds the lower by more than 1/10:

log(5/4)-19/160>1/10,

since exp(7/32)<5/4. The ratio |I_L/I_U| is at most 54exp(-N/10)<1. The same estimate holds inside this circle: the coefficients of I_L as a polynomial in w are positive, so its modulus is bounded by its value on the positive real boundary, whereas the upper saddle estimate holds throughout its disk. Thus I has no zero there.

On |w|=6W/5 the lower height exceeds the upper by more than 1/20, since log(6/5)-19/160>1/20 follows from exp(27/160)<6/5. The upper-to-lower ratio is at most 10sqrt(D)exp(-N/20)<1. Outside the circle apply the positive-coefficient comparison to I_U/w^N, whose coefficients in 1/w are nonnegative, and the lower stationary estimate, whose height tracks log|w|. This proves the outer zero-free region. The endpoint inequalities are far stronger than needed for N>=4000 and D<=N/128.

At a zero of I, both saddle domains therefore apply and I_U=-I_L!=0. The packet errors and leading ratios give

|J_U/I_U| <24sqrt(N/n)/|g(s_U)|,
|J_L/I_L| >(1/5)sqrt(N/n)/|g(s_L)|.

The second constant is valid because (1-18/25)/(1+2/5)=1/5. For the first, (1+3)/(1-5/6)=24. Finally

log|g(s_U)/g(s_L)| >= d-119/160-log(6/5),

using |x_U-D|<=1/2, |x_L-2|<=1/D, D>=8 and the two height errors. At d=6 the exponential of the right side exceeds 150, as certified by a rational exponential lower bound; it increases with d. Consequently |J_U/J_L|<120/150<1 and J!=0. This proves (T). In particular the proof uses the weighted actual consecutive pair; it does not replace the lower shape 2N-1 by the upper shape 2N+1.

## 7. Sparse overlap and the finite domain

The retained sparse theorem in lower_exponent_initial_and_boundary.tex states that, for a>=1,d>=2, the inequality

d^2>=4(a+d)log(8N)

makes all first N actual moment polynomials simple and pairwise coprime. Its proof is the cited Rouche argument; here a=2. At d=46, N0=128(d+2)=6144,

8N0=49152<exp(2116/192).

The difference d^2/[4(d+2)]-log(1024(d+2)) has derivative

1/4-1/(d+2)^2-1/(d+2)>0.

Thus the sparse theorem covers every 2<=N<=N0 for all d>=46, while Section 6 covers N>=N0. The fixed prime gaps d=2,3,5 are covered at every terminal by Theorem thm:prime of LOWER_TWO_PRIME_GAP_SFC.tex, whose beta=1 is the ordinary factorial functional. Theorem thm:fourth-gap of LOWER_TWO_FOURTH_POWER_GAP_SFC.tex covers d=4. The adjacent-exponent theorem in UNIFORM_PADE_PARABOLA_SFC.tex covers d=1. These are exact all-terminal antecedents, not inferences from finite scans.

It remains precisely the domain (F). Summing max(4000,128(d+2))-2 over d=6..45 gives 176816, with N<=6015. For every finite pair the actual normalized moment polynomial is

A_(d,N)(z)=sum_(k=0)^N binom(N,k) [(2N+dk)!/(2N)!] z^(N-k).

It is monic, has degree N and nonzero constant coefficient. Its consecutive partner uses N-1 and base shape 2(N-1)+1=2N-1; this differs from the upper shape 2N+1. The primary finite checker and both complete runs use the direct reciprocal orientation T_(d,N)(t)=t^N A_(d,N)(1/t), whose constant coefficient is one. Both endpoints are nonzero; reciprocal reversal preserves degree and coprimality. Thus this computation proves the coprimality of the displayed monic A polynomials as well.

**Corollary LowerTwoFiniteEquivalence.**

The full ordinary lower-two Strong Factorial statement over C, equivalently over every characteristic-zero field, holds if and only if ALL the actual rational consecutive pairs in (F) are coprime over Q. The forward direction restricts the full statement to (F): a nonconstant rational gcd would have a complex root, whose nonzero coordinate yields a nonzero binomial with two vanishing actual moments. The reverse direction combines rational coprimality in (F), the weighted tail, the sparse overlap and the retained small-gap theorems. Rational Bezout identities transfer this conclusion to every characteristic-zero field. This equivalence is between the characteristic-zero statements; it does not claim that rational coprimality forces success at either prescribed finite prime.

## 8. Exact finite computation and characteristic-zero transfer

The primary full finite computation checked all 176816 actual pairs in (F), at p=1000003, and found gcd one in every case, with no fallback prime. Every factorial argument is at most 47*6015=282705<p; all modular factorials used in the normalization are therefore units and every actual polynomial retains its full degree. The checker proves p prime by exact integer trial division and verifies the domain and deterministic per-gap case streams, literal coefficient fixtures, and malformed-shape/repeated-factor rejection controls.

A separately written complete replay checked the same 176816 pairs at p=1000033. It derives coefficients from binomial and sliding-rising recurrences rather than factorial tables and checks every full extended-gcd Bezout identity. It also succeeded in every case, with no fallback. This is additional independent evidence for the same finite rational coprimality statement; neither fixed prime is asserted to work for any domain outside (F).

If two rational polynomials with degrees preserved modulo p had a nonconstant common rational factor, their resultant would be zero over Q and hence modulo p after clearing unit denominators. A modular gcd one therefore proves rational coprimality. Conversely, rational coprimality is not claimed to force success at either specified prime. The resulting rational Bezout identity remains valid in every characteristic-zero field. The same applies to the analytic complex coprimality theorem in Section 6, because the two actual moment polynomials have rational coefficients. This proves the theorem in Section 1 for every characteristic-zero K.

## 9. Supporting scope and classical context

The numerical ledger uses exact rational Taylor bounds for exponential functions, rational square-root enclosures, and the three full Bernstein coefficient lists. A positive ledger result supplements the written contour argument; it does not turn sampled saddle values or numerical polynomial roots into proof. The finite modular checker certifies the finite domain (F), and the analytical ledger certifies the stated constants on their unbounded domains. Neither is substituted for the other's hypothesis.

Closest primary comparisons inspected were Bihun--Calogero, arXiv:1507.07808, Proposition 2.1 (generic hypergeometric zero relations), DLMF 18.34 (the adjacent-gap Bessel-polynomial representation and its specific orthogonality hypotheses), and Wolfs, arXiv:2401.08312, Section 4 and Theorem 4.6 (multiple Bessel-like hypergeometric formulas). Those statements do not supply actual changing-shape consecutive coprimality for arbitrary gap d. No claim of priority or novelty relative to all literature follows from this comparison. The proof above is self-contained apart from the explicitly named retained small-gap and sparse antecedents.
