# Uniform relative-gap weighted saddle tail

## 1. Statement and hypotheses

**Theorem RelativeGapThreeWeightedTail.** Let a>=2 and d>=3a be integers, D=a+d and rho=D/a>=4. For every integer N>=max(4000,128D), the two ACTUAL normalized ordinary factorial moment polynomials at exponents N-1 and N for X^a(A+BX^d) are coprime over Q. Consequently no nonzero binomial in this domain has both consecutive moments zero over any characteristic-zero field.

The proof uses the exact three full Bernstein certificates and elementary contour estimates in [the ordinary lower-two theorem](ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md), Sections2-6. It does not assert that arbitrary a has the same bounds without rescaling. This note proves the new parameter monotonicity, rescaled contours and weighted Gaussian estimates needed for that extension. The all-terminal finite obligation in Section8 is separate from the analytical theorem.

Classical comparison: [DLMF2.4(iv)-(v)](https://dlmf.nist.gov/2.4#iv) describes saddle expansions and the need for uniform treatment at coalescing saddles. Those asymptotic statements do not give the explicit finite onset here. The inequality rho>=4 keeps the two selected stationary disks separated; none of the following argument covers the coalescing regime a/d large.

## 2. Exact actual pair and parameter majorants

Put

H=D^D exp(-D), W=D^D exp(-d)/a^a, c=5D/8,
g(s)=s^a(w+s^d), Phi=log(g)-s.

With A/B=N^d w the actual moments are B^N N^(DN+1)I and B^(N-1)N^(D(N-1)+1)J, where

I=int_0^infinity g(s)^N exp(-Ns)ds,
J=int_0^infinity g(s)^(N-1) exp(-Ns)ds.

Both use the SAME N scale. Let n=N-1; the lower integrand is exp(nPhi(s))exp(-s). Split at c into U and L packets. The upper packet domain is |w|<=6W/5 and the lower packet domain is |w|>=4W/5.

For continuous a>=2,rho>=4 define

F_U(rho,beta)=log rho-(rho-1)(1+log beta),
F_L(rho,beta)=(rho-1)(1+log beta)-rho log rho.

The upper and lower perturbation majorants at x=beta D and x=beta a are respectively

q_U=(6/5)exp(aF_U),
q_L=(5/4)exp(aF_L).

The following monotonicity facts are used only for positive majorants; the analytic exponents a,d themselves remain integers.

(M1) For beta>=5/8 and beta<=7/8, q_U decreases in a and rho. Indeed F_U(4,5/8)<0, because exp(3)(5/8)^3>4 already follows from e>8/3; F_U decreases with rho since 1+log(5/8)>1/2 and 1/rho<=1/4. The correction h=(1-1/rho)q_U/(1+q_U) also decreases with rho: its logarithmic derivative is 1/[rho(rho-1)]+a(1/rho-1-log beta)/(1+q_U)<1/12-1/4<0, using q_U<1. The error 1/[D^2(1-q_U^2)] decreases in both parameters.

(M2) For beta>=7/8, q_U and the products [d(d+1)/D]q_U and [d^2/D]q_U^2 decrease in a,rho. For the first correction, the logarithmic a-derivative of its polynomial prefactor is at most1/a<=1/2, whereas F_U<=log4-18/7<-1/2, using log4<3/2 and log(7/8)>-1/7. Its rho-derivative is at most2/(rho-1)<=2/3, whereas a partial_rho F_U<=2(1/4-6/7)=-17/14. The second correction is even smaller; replacing q by y=q/(1-q) preserves these signs.

(M3) For beta<=7/4, q_L decreases in both parameters. At rho4,beta7/4, F_L=log[(7/4)^3 exp(3)/4^4]<0; the exact rational bound exp(3)<(87/32)^3 verifies this. Also partial_rho F_L=log(beta/rho)<=log(7/16)<0. The product (rho-1)q_L/(1+q_L) decreases: its logarithmic rho derivative is at most1/3+2log(7/16)/(1+q_L)<0 when q_L<=1/4. Multiplication by the harmless constant100/99 retains this conclusion.

(M4) The products (d^3/a)q_L at beta9/8, and (d(d-1)/a)q_L and (d^2/a)q_L^2 at beta251/200, decrease in a,rho. Their polynomial logarithmic derivatives are bounded by2/a<=1 and3/(rho-1)<=1 in the first case, and by6/(5a)<=3/5 and11/[5(rho-1)]<=11/15 for the d(d-1)/a correction; the d^2/a correction has the smaller bounds1/a and2/(rho-1). At beta9/8, F_L(4,beta)<-47/24 and partial_rho F_L<=-29/24, using log4>4/3,log(9/8)<1/8. At beta251/200, the weaker log beta<51/200 gives F_L<-941/600 and partial_rho F_L<-647/600. These strict bounds dominate the displayed polynomial derivatives. Higher powers of y=q/(1-q) only improve monotonicity.

## 3. Selected stationary points and uniform local data

Stationarity and curvature are

w(s-a)=s^d(D-s),
A(s)=-Phi''(s)=[s^2-(D+a-1)s+aD]/s^2.

There are unique analytic branches |s_U-D|<1/2 and |s_L-a|<a/(2D)=1/(2rho). For the upper Rouche circle, the ratio is bounded by

R_U=(12/5)(d+1/2)rho^a exp(-d)(1-1/(2D))^(-d).

It decreases in a,rho. To see this set L=-log(1-1/(2D)). Since L<1/(2D-1), partial_a(dL)=(rho-1)[L-1/(2D-1)]<0. Thus partial_a log R_U<=1/a+log rho-(rho-1)<0. For rho, dropping a negative term gives partial_rho log R_U<=1/(rho-1)+a/rho-a+a/(2D-1)<=1/3-3a/4+2/15<0. At a2,rho4 the exact e>163/60 bound is85761906966528/93776847890045<1.

For the lower circle use the triangle bound, without a boundary-correlation assertion:

R_L<=(5/2)rho(d+1/(2rho))rho^(-rho a)exp(d+1/2).

The factor before exp(1/2) decreases in a,rho: its logarithmic derivatives are at most1/a+(rho-1)-rho log rho<0 and1/rho+1/(rho-1)-a log rho<0. At2,4, use exp(1/2)<5/3 and e<87/32; the resulting exact rational bound is<1 (approximately0.6291). These disks select the branches from w0 and winfinity.

The curvature bounds remain

Re A_U>1/(2D), |A_U|<3/D;
Re A_L>1/(2a), |A_L|<2/a.

For a general disk the v=1/s transformation turns Re A into a quadratic in the boundary cosine. Its derivative is negative on the upper disk (use a<=D/4) and positive on the lower disk. At s_U=D-1/2 the minimum is [(D+a)/2-1/4]/(D-1/2)^2, at least the a2 minimum used in the retained proof. At s_L=a+a/(2D), multiplying the minimum by a gives

[1/2+(a+1)/(2D)+a/(4D^2)]/[1+1/(2D)]^2>1/2.

The norm estimates follow from A-1/D=(s-D)[(D-1)s-aD]/(Ds^2) and A-1/a=(s-a)[(a-1)s-aD]/(as^2). For example the latter difference has modulus at most18/(25a)<1/a, using D>=8,a/D<=1/4.

The stationary heights obey

|RePhi(s_U)-log H|<7/80,
|RePhi(s_L)-log|w|-a log a+a|<1/32.

For the upper height radial stationary-parameter integration gives q/(1-q), with q<32(60/163)^6 and q/(1-q)<7/80. Here (1-1/(2D))^(-d)<exp(1/2), rho^a exp(-d) is maximized at2,4 and exp(1/2)<5/3. For the lower height the logarithm Taylor error is at most a/[4D(2D-1)]<=1/240; stationarity gives |s_L^d/w|<1/[2rho d-1]<=1/47, hence its log error<1/46.

Use local constants

U: radius1/2, B_U=3/D^2, b_U=1/(7D);
L: radius a/16, B_L=4/(3a^2), b_L=1/(6a).

The whole upper local disk lies in |s-D|<=1. Its q bound is (6e/5)rho^a exp(-d). Both q and D^2q decrease in a,rho; the logarithmic derivatives for D^2q are at most1+log4-3<0 and1/2+1/2-2<0. At2,4 use q=(96/5)(60/163)^5. The exact inequality

2(8/7)^3+64[y+3y^2+2y^3]<18, y=q/(1-q),

bounds |Phi'''| by18/D^2. The retained positive polynomial6(D-8)^3+89(D-8)^2+354(D-8)+187 proves decay b_U. The whole lower local disk lies in |s-a|<=a/8. There

a^2|Phi'''| <=(8/7)^3[2+(d^3/a)(y+3y^2+2y^3)]<8,

by M4 and the retained a2 inequality. Curvature minus the cubic remainder gives b_L=1/(6a).

## 4. The whole complex contours reduce to the certified corner

Move the upper half-line c..infinity to height Im(s_U), adding the vertical at c. For the lower contour retain [0,a/2] on the real line and move [a/2,c] to height Im(s_L), adding the two endpoint verticals. All integrands are entire polynomials times exp(-Ns). One never divides a nonlocal integrand by g(s).

The upper central interval is [7D/8,5D/4]. M2 and the same estimate Re(D/s^2)>(39/40)D/x^2 show that the full complex curvature inequality RePhi''<-1/(4D) has its worst case at a2,rho4. It is therefore certified by the retained degree-fourteen P_U with all fifteen positive Bernstein coefficients. The upper left strip [5D/8,7D/8] uses the exact disk-Mobius bound

RePhi'>=D/x-1-d r/[x(1+r)]-D/[4x^3(1-r^2)].

After x=beta D, M1 shows its negative correction and error have their largest magnitude at2,4. The retained P_left and all sixteen positive Bernstein coefficients give RePhi'>1/100 on the WHOLE complex strip. The upper far strip x>=5D/4 obeys RePhi'<-1/10 by the same q bound and |q|<1/50. These arguments prove, respectively, the same upper central Gaussian tails and the left/right endpoint losses D/5120 and D/1280.

On the vertical at c the high contribution normalized by H is at most [(5/8)exp(3/8)]^D, decreasing in D. The low contribution is at most (6/5)[(5rho/8)exp(1-5rho/8)]^a, decreasing in both a and rho. Their maxima are again at2,4, and are respectively<15/32 and<3/8. The imaginary correction is exp[8/(25D)]<=exp(1/25). Thus the same vertical loss1/50 and integrated amplitude<=1/2 apply. Every upper weighted nonlocal term is UNCHANGED: s_U-D still has modulus<=1/2, n/N is unchanged and N>=128D.

The lower central interval is [3a/4,5a/4]. Its imaginary displacement is at most a/(2D), and |s|<(251/200)a. M4 shows that its curvature numerator divided by a is smallest at2,4. The exact retained inequality

39/20-30y_8-36y_8^2>189003/200000=(3/20)(251/100)^2

therefore proves RePhi''<-3/(10a). The lower phase consequently loses at least3u^2/(20a) from its saddle.

For the lower right strip [5a/4,7a/4], put x=a beta. The bound r<=(100/99)(5/4)exp(aF_L(rho,beta)) absorbs |s|^d/x^d<exp(1/(8D beta^2))<100/99. The exact Mobius estimate is

RePhi'<=a/x-1+d r/[x(1+r)]+d r^2/[4rho^2 x^3(1-r^2)].

Its main correction is [(rho-1)/beta]r/(1+r). M3 shows this decreases in both parameters. The last correction has prefactor (rho-1)/(4rho^2a^2 beta^3), also decreasing. Thus the retained P_right, with all sixteen positive Bernstein coefficients, gives RePhi'<-1/20 throughout this strip. The central endpoint5a/4 loses at least9a/2560. No assertion about maximization at real positive w is used.

For the lower left pieces |s|^2<=(145/256)a^2, and |s^d/w|<1/500 by M3 with beta151/200. Their normalized modulus is bounded by

[(sqrt(145)/16)exp(1/4)]^a(501/500).

The base is less than one. Its maximum a2 is exactly the retained lower-left cap, including height error1/32 and loss1/100. The real origin segment has a smaller envelope by monotonicity. Their total length is less than3a/2.

For the lower far pieces, x=a v with7/4<=v<=5rho/8, the real envelope divided by |w|a^a exp(-a) is

E(a,rho,v)=[v exp(1-v)]^a+(5/4)[(v/rho)^rho exp(rho-v)]^a.

For fixed v<rho both bases are below one, and the high base decreases in rho. Hence on7/4<=v<=5/2 the envelope is at most its a2,rho4 value. The four retained interval caps apply after x_old=2v on [7/2,4],[4,9/2],[9/2,19/4],[19/4,5]. For v>=5/2 bound the low base at5/2 and the high term at c=5D/8. Both then decrease with D or a, so the retained fifth cap applies. The imaginary factor is at most exp(1/196), since D(Im s_L)^2/(2x^2)<=1/(8D v^2). Including the1/32 height error, every far part loses1/200. Its entire length is <2D. This verifies the ENTIRE contour for every a,rho, including the domain expansion beyond the base corner.

## 5. Uniform weighted local Gaussian error

The lower local interval has radius a/16, so the estimate exp(|u|)<=exp(a/16) is NOT uniform and is not used. For any m,c>0 and epsilon=1/24, Young's inequality gives

exp(|u|)<=exp(mc u^2/24+6/(mc)).

Consequently the local phase-comparison integral has multiplier

C_phase=(24/23)^2 exp(6/(mc)),

and the amplitude-variation integral has multiplier

C_amp=(24/23)exp(6/(mc)).

These follow by integrating |u|^3 and |u| against exp[-(23/24)mc u^2]. The unweighted Gaussian completion is unchanged. This proves a uniform estimate for the ACTUAL lower exponent, whose amplitude is exp(-u), while keeping exactly the same Gaussian leading term as in the retained proof.

The leading moduli are >sqrt(2D/N) at U and >sqrt(3a/N) at L, because pi>3 and the curvature norm bounds are3/D and2/a. For the lower J local comparison, n/a>511, so exp(36a/n)<511/475. Define the rational constants

C_p=(24/23)^2(511/475), C_a=(24/23)(511/475).

The complete local relative error is at most

E_local(a,N)=48C_p/sqrt(3na)+6C_a sqrt(a/(3n))<9/20.

Here N>=max4000,512a. For2<=a<=125/16, put n>=3999 and t=sqrt(a); the expression is A/t+Bt, convex in t, so its maximum is at a2 or a125/16. Exact rational square-root enclosures give respectively less than0.44994533 and0.35553366, both<9/20. For a>=125/16 put n>=512a-1; both 1/sqrt[a(512a-1)] and sqrt[a/(512a-1)] decrease, reducing to the same endpoint. This is a proof on an unbounded parameter region, not a finite a scan.

For the upper weighted local comparison n/D>127968/1000. Therefore exp(42D/n)<=127968/85968, which gives C_phase<5/3 and C_amp<8/5. Keeping even the weaker leading lower bound sqrt(D/n), its local relative error is at most

(5/3)(1000/999)(33/40)+(8/5)(5/8)<5/2.

Thus the retained upper weighted local budget remains valid.

## 6. Full packet budgets and actual pair separation

All upper packet error terms are precisely those in ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md Section5: the local constants/radius, whole contour constants and N/D onset are unchanged. The new Young bound replaces its local factor exp(1/2) and still gives the same upper weighted local cap5/2. Hence

E_IU<5/6, E_JU<3.

For the lower unweighted local and Gaussian completion, the leading sqrt(3a/N) is larger than the old2/sqrtN at every a>=2. The local comparison is bounded by the old24/sqrtN. The completion exponent is Na/1536>=N/768. The central absolute tail is (320/(3N))exp(NS_L-3Na/5120), and after division by the leading term it is bounded by the old lower central error160/(3sqrtN)exp(-3N/2560). All right/far terms are similarly bounded by the old terms, since the endpoint losses scale by a/2>=1.

The lower left term needs an explicit new length estimate. Its relative unweighted error is <=(sqrt3/2)sqrt(aN)exp(-N/100)<=N/16 exp(-N/100), because a<=N/512. It is<1/1000 for N>=4000. For J, the gain from exp(-s) is at most exp(a+1/8), and n=N-1. The exponent satisfies

-n/100+a+1/8 <=-N/128

for N>=4000,a<=N/512. Thus its weighted relative error is <=N/16 exp(-N/128)<1/1000. This replaces the nonuniform fixed-a left gain rather than silently reusing exp(17/8).

For the weighted lower central tail, n*3u^2/(20a)-|u|>=3Nu^2/(40a) for |u|>=a/16 and N>=216, so its absolute bound is (640/(3N))exp(nS_L-x_L-3Na/10240). This again is at most the old lower weighted central relative term. On the right strip the amplitude only helps, since x_L<5a/4. On the far contour x>=7a/4 and x_L<=a+a/(2D), so x_L-x+1/200<0; n may be replaced by N in its exponent. The whole far relative error is bounded by DsqrtN exp(-N/200), exactly as before.

Summing, E_IL<2/5. For E_JL the six budgets are bounded by9/20,5/1000,165/1000,1/1000,1/1000,1/1000, whose sum is<18/25. Therefore all four retained packet inequalities hold uniformly for a>=2,rho>=4.

The zero-free comparison carries over without alteration of the annulus4W/5<|w|<6W/5. The lower height now tracks log|w|+a log a-a, exactly balanced with H by W. The lower leading grows like sqrt(a/N), so its ratio to the upper leading has factor sqrt(a/D)=1/sqrt(rho)<=1/2; the old coarse constants54 and10sqrtD remain valid. Positive polynomial coefficients of the original packets justify extension of the endpoint comparisons to the inner/outer regions, as in the retained proof.

At an I zero, I_U=-I_L!=0. For the ACTUAL J leading term the ratio is still sqrt(N/n)/g(s_*), with the SAME square-root factor at both saddles. The two height errors and stationary displacements give

log|g(s_U)/g(s_L)|>=d-119/160-log(6/5).

Indeed x_U-x_L>=d-1/2-1/(2rho)>=d-5/8. Since d>=6, the ratio is>150. The four packet budgets give |J_U/I_U|<24sqrt(N/n)/|g_U| and |J_L/I_L|>(1/5)sqrt(N/n)/|g_L|, hence |J_U/J_L|<120/150<1. Thus the actual pair has no common complex zero. Its rational coefficients give the asserted Q-coprimality and characteristic-zero transfer.

## 7. All-terminal sparse overlap

**Corollary RelativeGapThreeLargeGapAllTerminals.**

For every characteristic-zero field K, every integers a>=2 and d>=max(64,3a), every N>=2 and every A,B in K not both zero, at least one of the ordinary factorial moments of X^a(A+BX^d) at orders N-1,N is nonzero.

The older actual sparse theorem gives pairwise coprimality through terminal N whenever d^2>=4Dlog(8N). In this region D=a+d<=4d/3, hence d^2/(4D)>=3d/16. At d64 the exact inequality exp(12)>4096*64/3 follows already from (8/3)^12>4096*64/3. Also

d/d_d[3d/16-log((4096/3)d)]=3/16-1/d>0.

Consequently every d>=64 is covered for all2<=N<=128D by the sparse theorem, and for N>=128D by the new tail. The number4000 is irrelevant at d>=64 because128D>=128(d+2)>4000. This establishes all-terminal SFC on the unbounded relative-gap class d>=64,d>=3a without a finite scan.

## 8. Genuine finite equivalence, not an executed closure

**Corollary RelativeGapThreeFiniteEquivalence.**

Since least exponents1 and2 are already proved at every gap, the all-terminal relative-gap statement a>=2,d>=3a is equivalent over C, equivalently over EVERY characteristic-zero field, to actual rational consecutive coprimality in

9<=d<=63, 3<=a<=floor(d/3), 2<=N<max4000,128(a+d).

There are532 support pairs(a,d) and3694040 actual terminal pairs, maximum terminal10751 and maximum factorial903084. Every parameter occurs in the ACTUAL shapes a(N-1)+1 and aN+1, separately; no fixed-shape replacement is licensed. The characteristic-zero forward implication is restriction, and the reverse combines all known least-exponent1/2 results, Section7 and the uniform tail. A modular computation at a degree-preserving prime would be one possible way to prove this finite rational obligation, but no finite success is asserted by this analytical note.
