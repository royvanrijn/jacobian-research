# Native binary quadratic series modulo four and a finite quarter-critical reduction

This note proves finite identities for the remaining quarter-critical six-zero question. It does not exclude such a block. The general native rational formula in Sections1--2 has no rank or residue-profile assumption. Section3 imports the explicit quarter-gap profile from `BINARY_QUADRATIC_QUARTER_GRID_LIFTED_NATIVE_PROFILE.md`; this hypothesis is stronger than merely saying the coefficient-grid denominator has 2-adic valuation two.

Work over a characteristic-zero real-valued nonarchimedean extension of v2, normalized by v(2)=1. All congruences are coefficientwise in O/4O, or in O/2O when explicitly stated. They concern formal series, without convergence or evaluation at a nonintegral formal variable.

## 1. Shifted univariate edges over the universal integer coefficient ring

For q_x=alpha*X^2/2+d*X put

 F_x(T)=sum_n Lambda(q_x^n)T^n/n!,
 C_(j,x)(T)=sum_n Lambda(X^j*q_x^n)T^n/(j!*n!).

Every coefficient belongs to Z[alpha,d]. For i+l=n its scalar coefficient is

 (2i+l+j)!/(j!*2^i*i!*l!)
 =binom(2i+l+j,j)*(2i-1)!!*binom(2i+l,2i),                 (1)

with (-1)!!=1. Thus the following derivations take place in a universal integer polynomial ring before specialization; reduction and differentiation in alpha,d preserve the relevant congruences.

Write D_x=1-alpha*T-d^2*T^2 and U_x=1+d*T. The exact native recurrence gives

 D_x F_x=U_x+2alpha*T^2*F_x'.                              (2)

Indeed the coefficients u_n satisfy u0=1,u1=alpha+d and u_n=alpha(2n-1)u_(n-1)+d^2u_(n-2) for n>=2. Substitution of F0=U_x/D_x once on the derivative side of(2) gives, modulo4,

 F_x=U_x/D_x+
       2alpha*T^2*(alpha+d+d^3*T^2)/D_x^3.                 (3)

For the substitution, F_x-F0 is divisible by2 coefficientwise. Its derivative is also divisible by2. The omitted term2d^2T in the numerator of F0' is multiplied by another2 and vanishes modulo4.

Let A_x=C_(1,x). On each monomial of F_x the extra factor from inserting X is 2i+l+1=n+i+1. Hence the EXACT operator identity is

 A_x=F_x+T*partial_T F_x+alpha*partial_alpha F_x.

Applying this identity to(3), and reducing only the already twice-multiplied correction modulo2, gives

 A_x=[1+(alpha+2d)T+d^2T^2]/D_x^2
                  +2alpha^2*T^2/D_x^3.                  (4)

This proof divides by no specialized coefficient alpha. In the correction, in characteristic two the T and alpha derivatives of D_x^(-3) cancel, leaving alpha^2*T^2/D_x^3.

For the two remaining shifted edges needed below use Hasse derivatives in d. Since the coefficient of s^j in F_x(alpha,d+s;T) is T^j*C_(j,x), the modulo-two identity

 F_x(alpha,d+s;T)=(1+(d+s)T)/(D_x+s^2T^2)

gives

 C_(2,x)=U_x/D_x^2,  C_(3,x)=1/D_x^2 modulo2.             (5)

All denominators are formal units with constant coefficient one. The division by T^j in extracting these Hasse derivatives is legitimate because their coefficient series starts at T^j. Define D_y,U_y,F_y,A_y analogously with gamma,e replacing alpha,d.

## 2. A finite rational formula for every native binary quadratic modulo four

For

 q=alpha*X^2/2+beta*XY+gamma*Y^2/2+dX+eY,
 R(T)=sum_r Lambda(q^r)T^r/r!,

formal expansion of the cross term yields the EXACT series identity

 R(T)=sum_(j>=0) beta^j*j!*T^j*C_(j,x)(T)*C_(j,y)(T).     (6)

Equation(1) proves integrality of every factor in this expression. For j>=4, j! is divisible by4 (in fact by8), so all such summands vanish modulo4. For j2,3 use(5), and note 3!=2 modulo4. This gives the universal native congruence

 R=F_x F_y+beta*T*A_x A_y
       +2beta^2*T^2*(U_x U_y+beta*T)/(D_x^2 D_y^2) mod4.  (7)

Here F_x,A_x are the explicit expressions(3),(4), not independent states. In particular, with

 L_x=U_x D_x^2+2alpha*T^2*(alpha+d+d^3*T^2),
 J_x=[1+(alpha+2d)T+d^2T^2]D_x+2alpha^2*T^2,
 N=U_x U_y+beta*T,

one has

 R=[L_x L_y+beta*T*J_x J_y+2beta^2*T^2*N*D_x D_y]
                                          /(D_x D_y)^3 mod4. (8)

The numerator has degree at most ten and the denominator degree at most twelve before specialization. Thus the next layer of the native problem is a finite rational identity, not an arbitrary holonomic or recurrence state. The cross multiplicity cannot simply be discarded modulo4: j2 and j3 contribute the last term in(7).

## 3. Arbitrary starts n=8H+2 reduce to two small native systems

For this section let k!=0,v(k)=3/4, U=k^4/8, and a_j=k^j/j!. Let the natural coefficients of q be integral and normalize d=1. Assume the necessary full modulo-two quarter profile already proved in the imported note: with C=R1 and R_r=Lambda(q^r)/r!,

 U=1 modulo2O,
 R(T)=[1+C*T+(lambda+k*C)T^2+lambda*C*T^3]/(1+U*T^4)
                                               modulo2O, (9)

where lambda=k^2/2. In particular R_(r+4)=U R_r modulo2O for every r>=0. Also, writing

 E_s=sum_(j=0)^s a_j R_(s-j)=Lambda((q+k)^s)/s!,

the same profile gives E0=1,E1=C+k,E2=E3=0 modulo2O. These are necessary data from the actual six-zero block, not freely assignable recurrence values.

Suppose the actual six moments at n,...,n+5 vanish with n=8H+2. Then for s2,...,7 divide M_(8H+s) by its nonzero r=s coefficient (8H+s)_s*k^(8H). The normalized coefficient of R_r is

 c_(s,r)=(8H)!/(8H+s-r)! * k^(s-r).                      (10)

For r<=s, put j=s-r. Its ratio to a_j is binom(8H+j,j)^(-1). For 0<=j<4 this ratio is one modulo4. For 4<=j<=7 it is (-1)^H modulo4: among factors 8H+i relative to i, only i4 contributes a nontrivial odd part modulo4, namely 1+2H. The coefficients a_j in this range are integral, so

 c_(s,s-j)=a_j (j<4), (-1)^H*a_j (4<=j<=7) modulo4.       (11)

For r>s+8 every normalized base cost is at least two. One way to see this is to put m-r=8u+w,0<=w<8. Relative to u=H,w0, the high cost is

 s2(u)-s2(H)+2(H-u)>=H-u,

because decreasing an integer by one decreases its binary population by at most one. The low cost s2(w)-w/4 is nonnegative. When r>s+8, H-u>=2. Thus all these tails vanish modulo4 without any information beyond integrality of R_r.

For the first tail put r=s+t,1<=t<=8. If H is even, v((8H)_t)-3t/4>=2, so again every term vanishes. For odd H, let w=8-t. Since

 v((8H)_t)+v(w!)=v(8!)=7,

the integer (8H)_t*w!/128 is odd. As k^8=64U^2,

 c_(s,s+t)=2U^(-2)*a_(8-t) modulo4.                      (12)

The error has valuation at least2 because every a_w,w0..7, is integral. Terms w3,7 already have valuation greater than one before the factor2; they can equivalently be omitted.

The entire first tail in(12) cancels modulo4. In O/2O, for j0..3,

 a_(j+4)=U a_j,
 R_(s+8-j)=U R_(s+4-j).

The two halves of sum_(w0..7) a_w R_(s+8-w) are therefore equal, so their sum is zero modulo2O. Multiplication by 2U^(-2) gives zero modulo4. This cancellation uses the established full native profile(9), not just the residue-field pole count.

Only the low-head sign change remains. For odd H it is

 -2 sum_(j4..s) a_j R_(s-j)=-2U E_(s-4) modulo4.

Since E0=1,E1=C+k,E2=E3=0 modulo2O, the normalized six-zero equations are exactly the following NECESSARY finite systems:

 H even: E2=E3=E4=E5=E6=E7=0 modulo4;
 H odd:  E2=E3=E6=E7=0, E4=2U, E5=2U(C+k) modulo4.       (13)

All E_s here are the actual small native polynomials Lambda((q+k)^s)/s!, not arbitrary symbols. Thus no R8,...,R15 remain in the final equations, and H enters only through its parity. Conversely, these modulo-four conditions alone do not imply actual large-index zero moments.

Equivalently the prefix of exp(kT)R(T) through degree seven is

 [1+(C+k)T]*(1+2(H mod2)*U*T^4) modulo (4,T^8).          (14)

Since exp(-kT) has integral coefficients through degree seven, one can multiply its finite prefix to obtain an explicit target prefix for R(T). This is an entirely finite operation; no integrality is asserted for higher coefficients of exp(kT) (its eighth coefficient has negative valuation here).

## 4. The remaining problem

Combining(8) with(13) or(14), and the native quarter profile, gives two finite polynomial systems for the next layer. They may have compatible solutions; this note makes no contradiction or all-start SFC claim. In the profile, v(e)=1/2,v(alpha)>=1/2,gamma is a unit, so the edge expression 2gamma*k-e^2 has valuation exactly one. Any later use of that pivot must retain this proof of its nonvanishing rather than assume a generic minor is nonzero.

