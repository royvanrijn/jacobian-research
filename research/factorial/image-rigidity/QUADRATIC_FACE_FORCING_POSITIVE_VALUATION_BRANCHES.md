# Positive 2-adic branches of the symmetric quadratic face forcing

Let h_j(r,z)=Lambda((r X^2+X+z)^j), where Lambda(X^d)=d!, and put
T_j(r,z)=h_j(r,z)+j h_(j-1)(r,z), j>=1. A valuation below is normalized
by v(2)=1 on an algebraic extension of Q_2.

**Proposition 1.** There is no algebraic pair (r,z), z!=0, with v(z)>=1
and T_k=T_(k+1)=T_(k+2)=0 for any integer k>=1. No condition on r is
needed. Thus the positive branches v(z)>=1 of the necessary cubic below
cannot support a three-zero forcing block.

This proposition does not settle the branch v(z)=1/2, and does not by
itself prove the universal three-window. The remaining half-valuation
family is isolated in the [residual-frontier proof](QUADRATIC_FACE_FORCING_THREE_ZERO_RESIDUAL_FRONTIER.md).

## 1. Native cubic and normalized coefficient reductions

The native factorial recurrence is

 h_j=2rj(2j-1)h_(j-1)+j(j-1)(1-4rz)h_(j-2)+z^(j-1)(z+j).

For completeness, apply Lambda(P')=Lambda(P)-P(0) to phi^j and
phi' phi^(j-1), using (phi')^2=4r phi+1-4rz. Eliminating
Lambda(phi' phi^(j-1)) gives the displayed recurrence without a
coefficient-dependent division.

Two consecutive T zeros give

 2ri(2i-1+2z)h_(i-1)+z^(i-1)(z+i)=0.

Combining this at consecutive i gives, for a three-zero block starting k,
with i=k+1,

 P_i(z)=2(z+i)^3+3(z+i)^2+i
      =2z^3+(6i+3)z^2+6i(i+1)z+i(i+1)(2i+1)=0.       (1)

These identities use no coefficient-dependent division. Since P_i has
nonzero constant coefficient, z=0 is not a root. Write S=v_2(i(i+1))>=1.
The lower Newton polygon of (1) has vertices (0,S),(2,0),(3,1).
Consequently its roots have valuation S/2 or -1.

Put R=2r and U_j(R)=h_j(R/2,z)/j!, V_j(R)=T_j(R/2,z)/j!
=U_j(R)+U_(j-1)(R). A direct multinomial expansion gives

 [R^a]U_j=(2a-1)!! sum_(d=0)^(j-a) binom(j+a-d,2a) z^d/d!, 0<=a<=j, (2)

with (-1)!!=1. If v(z)>=1, every d>=1 term is in the maximal ideal,
since v(z^d/d!)>=s_2(d)>=1. Thus U_j and V_j are integral polynomials
with odd leading coefficient (2j-1)!!. Their reductions are independent
of z. The formal generating series of those reductions are

 sum Ubar_j t^j=(1+t)/(1+Rt+t^2),
 sum Vbar_j t^j=(1+t^2)/(1+Rt+t^2).                  (3)

For example, Ubar_j=sum_a binom(j+a,2a)R^a; summing first in j proves
(3) directly in F_2[R][[t]]. For j>=1 write v_j=Vbar_j. Then
v_1=R, v_2=R^2, and v_j=R v_(j-1)+v_(j-2) for j>=3.
Equivalently v_j=R W_(j-1), where W_0=1,W_1=R and W_j=R W_(j-1)+W_(j-2).
Consecutive W's are coprime by this recurrence. Hence two consecutive
v_j have no common root except R=0, over any extension of F_2.
Also the coefficient of R in v_j is 1 for odd j and 0 for even j.

## 2. Constant terms and the exclusion except for two small starts

Let E_j(z)=sum_(d=0)^j z^d/d!. The constant term is

 V_j(0)=2E_j(z)-z^j/j!.                              (4)

For v(z)>=1, E_j(z) is a unit, and v(V_j(0))>=1. If v(z)>1 this valuation
is exactly 1 for every j>=1. If v(z)=1, it is exactly 1 whenever j is not
a power of two, because v(z^j/j!)=s_2(j)>=2 in that case.

Suppose two consecutive V_j have a common root R. Their unit leading
coefficients imply v(R)>=0. By (3) and coprimality of the W's, the common
root reduces to R=0, so v(R)>0. The odd member of the consecutive pair
has unit linear coefficient. If v(R)<1, its linear term is uniquely of
least valuation: the constant has valuation >=1, all terms of degree >=2
have valuation >=2v(R)>v(R). Therefore v(R)>=1.

For an even index J whose constant term in (4) has valuation 1, evaluation
at v(R)>=1 is impossible: its linear coefficient is even by (3), so every
positive-degree term has valuation >=2, while its constant has valuation 1.
This excludes v(z)>1 immediately, and excludes v(z)=1 whenever the triple
contains an even index which is not a power of two.

In the remaining v(z)=1 case, the cubic (1) forces S=2. If k is odd, the
middle i=k+1 is even, with v_2(i)=2; if it is a power of two then i=4,
hence k=3. If k is even, both even indices k and k+2 would need to be
powers of two; the only positive pair of powers of two differing by 2
is 2,4. Hence k=2. It remains to exclude exactly starts 2 and 3.

## 3. Exact two finite eliminations

For k=2,3 let A_k(z)=Res_r(T_k(r,z),T_(k+1)(r,z)), with the usual
fixed-degree Sylvester resultant. A common finite r forces A_k(z)=0.
Its literal integer inputs are

 T_j(r,z)=sum_(a+b+c=j) j!(2a+b)! r^a z^c/(a!b!c!)
   +j sum_(a+b+c=j-1) (j-1)!(2a+b)! r^a z^c/(a!b!c!).

Thus both resultants are determinants of explicitly specified matrices
of sizes5 and7. Their z-degrees are at most k(k+1): in each Sylvester
term the selected r-exponents sum to k(k+1), while a coefficient of r^a
in T_j has z-degree at most j-a. The column-index sum gives that fixed
r-exponent sum. Therefore the displayed nonzero top modular coefficients
also prove the full characteristic-zero degrees, without an unrecorded
leading-degree hypothesis.
For the primes below the full leading degrees of A_k and P_(k+1) are
preserved; the displayed polynomial Bezout identities prove coprimality.
All polynomials in each row are reduced modulo its stated prime.

For k=2, p=7:

 A_2=5z^6+z^5+5z^4+4z^3+2z^2+6z+3,
 P_3=2z^3+2z,
 (6z+5)A_2+(6z^4+2z^3+z^2+5z+4)P_3=1.              (5)

For k=3, p=11:

 A_3=9z^12+8z^11+7z^10+7z^9+5z^8+5z^7+8z^6
       +10z^5+2z^3+2z^2+9z+2,
 P_4=2z^3+5z^2+10z+4,
 B_3=8z^2+6z+8,
 C_3=8z^11+9z^10+9z^9+6z^8+7z^7+9z^6+z^5+3z^4
       +z^3+4z^2+9z+10,
 B_3 A_3+C_3 P_4=1.                                (6)

The leading coefficients of A_2,A_3 modulo their primes are 5,9; their
characteristic-zero degrees are 6,12, respectively. The P polynomials
both have degree 3 and leading coefficient 2. Thus modular coprimality
implies coprimality in Q[z], excluding both finite starts. Alternatively,
the nonzero fixed-degree double resultants in (5)--(6) provide the same
necessary-direction contradiction. This proves Proposition 1.

## 4. Algebraic and characteristic-zero scope

If a triple-zero solution is sought over an arbitrary characteristic-zero
field with e,c nonzero, normalize r=a/e,z=c/e. Equation (1) makes z algebraic
over Q; T_k is a nonconstant polynomial in r with leading coefficient
(2k)!, so r is algebraic too. Hence its coefficient field is a number field,
and each embedding into Q_2bar has a well-defined valuation considered
above. Proposition 1 excludes every such embedding with v(z)>=1. It does
not assert that every solution has such an embedding; the -1 and 1/2
branches require separate arguments. No finite test is used to infer an
all-index theorem.
