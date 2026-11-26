# Critical constant shifts and a six-window from the native binary matching residue

This note proves a six-moment window at every positive start in a specified critical local arithmetic sector. Polynomial residues and complete one-pole resonance cancellations are excluded explicitly; no unrestricted binary-quadratic SFC conclusion is asserted.

Let K be a characteristic-zero field with a real-valued nonarchimedean valuation extending v2, with v(2)=1, O its valuation ring and kappa its residue field. Write
 q=alpha X^2/2+beta XY+gamma Y^2/2+dX+eY,
 f=q+k,
with alpha,beta,gamma,d,e in O and k!=0. Put R_j=Lambda(q^j)/j!, R_0=1. The native quadratic coefficient-integrality argument gives R_j in O. Its universal binary matching residue is

 Rbar(T)=sum_(j>=0) Rbar_j T^j=N(T)/D(T),                    (1)
 D(T)=(1+alpha_bar T+d_bar^2 T^2)(1+gamma_bar T+e_bar^2 T^2),
 N(T)=1+(d_bar+e_bar+beta_bar)T+d_bar e_bar T^2.

The complete coefficient-level proofs of integrality and identity(1) are reproduced in Sections1 and5 of the [noncritical companion](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md). They have NO requirement beta^2=alpha gamma or fixed-weight relations. Only these universal native coefficient identities are used here, not a matching-family upper theorem. Thus degD<=4,degN<=2,D(0)=N(0)=1.

**Theorem 1 (critical six-window).** Suppose
 v(k)=1-1/L,  L=2^a, a>=2,
 and let mu be the nonzero residue of k^L/2^(L-1). If
                     B(T)=(1+mu T^L)Rbar(T)               (2)
is NONPOLYNOMIAL, then f has no six consecutive zero factorial moments at any positive start. In particular this proves its actual-support SFC window whenever it has all six binary-quadratic monomials.

A sufficient condition for(2) is that Rbar has at least TWO DISTINCT geometric poles, since 1+mu T^L has just one geometric root. A one-pole residue also works whenever that pole is not canceled completely by this factor. No phase with fewer monomials is assigned its smaller actual-support window by this six-window theorem alone.

For a=1, v(k)=1/2, the separate elementary constant-shift root bound already gives a four-window independently of(1); no pole hypothesis is needed there.

## 1. Exact endpoint minimizers of the native moment expansion

For each positive integer m,
 M_m=Lambda(f^m)=sum_(r=0)^m (m)_r k^(m-r)R_r,              (3)
where (m)_r is the falling factorial. By Legendre's formula, at the critical value delta=v(k)=1-1/L, the coefficient valuation after removing k^m is
 C_r=v((m)_r)-r delta
    =r/L+s_2(m-r)-s_2(m).                                 (4)
Only v(R_r)>=0 remains outside this exact cost.

Put j=m-r. Minimizing(4) is equivalent to minimizing
 F(j)=s_2(j)-j/L over 0<=j<=m.
Write j=2L u+w with 0<=w<2L. Then
 F(j)=s_2(u)-2u+s_2(w)-w/L.
The first part is strictly decreasing in u, since
 s_2(u+1)-s_2(u)<=1<2.
In the second part, every binary place below a has strictly positive weight1-2^i/L, whereas place a has weight0. Therefore its only minima are w=0,L, when these are allowed. A smaller u cannot compensate: each increment of u lowers the first part by at least1, and the low part is nonnegative. Thus use the largest u with2L u<=m, and include both low minima exactly when m has its a-th bit set.

Writing m mod2L=t for0<=t<L, the unique minimizing index in(3) is r=t. Writing m mod2L=L+t, its two minimizing indices are r=t,t+L. All other r have STRICTLY larger coefficient valuation. This includes small m; no asymptotic size condition is used.

## 2. The leading residue is an explicit finite-period test

Set A_m=(m)_t k^(m-t), where t=m modL. This is nonzero and has the minimum coefficient valuation in(3), without any R_t factor included. If m mod2L=t<L, then
                   residue(M_m/A_m)=Rbar_t.               (5)
If m mod2L=L+t, then
 residue(M_m/A_m)=Rbar_t+mu^(-1)Rbar_(L+t).                 (6)

To check the ratio in(6), divide the r=t+L coefficient by the r=t coefficient. It is
       (m-t)_L / k^L.
Here m-t=L(2u+1). Exactly
 v((L(2u+1))_L)=a+v((L-1)!)=L-1.
After dividing by2^(L-1), the remaining rational integer is odd, so its residue is1. The coefficient ratio therefore has residue mu^(-1). Every other term vanishes in the quotient residue by its strict cost inequality and R_r integrality. There is no assumption that R_t or R_(t+L) is a unit.

Consequently, a zero actual moment forces the corresponding entry of the following2L-periodic finite test vector to vanish:
  C_t=Rbar_t                      for0<=t<L,
  C_(L+t)=Rbar_(L+t)+mu Rbar_t     for0<=t<L.               (7)
Multiplication of(6) by the unit mu does not affect whether it is zero. The vector in(7) is exactly the coefficient prefix through degree2L-1 of B(T) in(2). It is not claimed that all later B coefficients are periodic. Only the test vector, indexed by m mod2L, is periodic.

Thus a forbidden block of six native moments would give a circular block of six zero entries of(7). Since C_0=Rbar_0=1, such a block cannot wrap through index0. It remains to exclude an ordinary contiguous six-zero block in the displayed finite prefix.

## 3. Rational coefficient windows needed below

A nonzero proper rational series U(T)/V(T), with V(0)=1 and positive denominator degree r, cannot have r consecutive zero coefficients at ANY nonnegative index. The order-r recurrence holds at all indices>=r and has nonzero last coefficient, so a zero block propagates backward to indices0,...,r-1 and then forward to every index, contradicting nonzero U. Equivalently the recurrence companion matrix is invertible. This statement applies to the reduced denominator; an upper bound for its degree can only enlarge the prohibited zero block.

Let d0=degD in(1). Nonpolynomiality of B implies nonpolynomiality of Rbar. Define the infinite second-half comparison sequence
 S_t=Rbar_(L+t)+mu Rbar_t, t>=0.                           (8)
Its generating series is
 S(T)=[B(T)-sum_(j=0)^(L-1) Rbar_j T^j]/T^L.               (9)
Thus S is nonpolynomial exactly when B is nonpolynomial, and its reduced denominator divides D. Polynomial subtraction and the shift do not create new poles.

If d0>=3, Rbar is proper because degN<=2<d0. Its shifted series and S are also proper with denominator dividing D. Each has no d0 consecutive zero coefficients at any nonnegative start.

If d0<=2, the polynomial part of Rbar has degree at most2-d0. Since L>=4, shifting by L kills that polynomial part, so the polynomial part of S has the same bound. On positive indices at most max(0,2-d0) coefficients can receive a polynomial-part contribution; hence Rbar has no TWO consecutive zero coefficients at a positive start. On nonnegative indices at most3-d0 coefficients can receive that contribution; therefore S has no THREE consecutive zero coefficients, even starting at0. These statements follow by leaving a suffix of at least d0 coefficients in a nonzero proper part. A constant D would make Rbar polynomial and is already excluded, so d0>=1 here.

## 4. The only apparent six-zero crossing forces an impossible native numerator

If a proposed six-zero block lies entirely in the first half, it is a block of Rbar coefficients. The bounds in Section3 exclude it. If it lies entirely in the second half, it is a block of S coefficients, equally excluded.

Suppose instead it crosses index L. Let A be the positive number of its entries before L and Bcount the positive number at or after L; A+Bcount=6.
For d0<=2, A>=2 is impossible by the positive-start Rbar bound. If A=1 then Bcount=5, impossible by the S bound. For d0=3, either A or Bcount is at least3, impossible because both series are proper of denominator degree at most3.

Only d0=4 remains. Both proper series forbid any four-zero block. Therefore the ONLY surviving split would be A=Bcount=3:
 Rbar_(L-3)=Rbar_(L-2)=Rbar_(L-1)=0,
 Rbar_(L+t)+mu Rbar_t=0 for t=0,1,2.                      (10)
Write D=1+d1 T+d2 T^2+d3 T^3+d4 T^4, d4!=0.
The coefficient recurrence is homogeneous from order3 onward because degN<=2; in particular it holds at L,L+1,L+2 since L>=4. Using the first three zeros in(10), it gives
 Rbar_L=mu,
 Rbar_(L+1)=d1 mu,
 Rbar_(L+2)=(d1^2+d2)mu.                                  (11)
Here the first equality uses the second line of(10) at t=0 and Rbar_0=1. In characteristic2 the recurrence signs are plus.
Comparing(11) with the second line of(10) at t=1,2 gives
 Rbar_1=d1,    Rbar_2=d1^2+d2.
The initial coefficient identities D Rbar=N therefore force
                        N_1=N_2=0.                       (12)
But the NATIVE binary matching forms in(1) have
 d4=d_bar^2 e_bar^2 and N_2=d_bar e_bar.
The condition d4!=0 implies N_2!=0, contradicting(12).
This is the extra native relation missing from a mere degree-four recurrence plus degree-two forcing count. The apparent three-plus-three boundary case is impossible.

All six-zero circular test blocks are excluded. Equations(5)--(7) now exclude a six-zero block of ACTUAL factorial moments, proving the theorem.

## 5. Pole criterion and exact remaining exceptions

Over an algebraic closure of kappa, 1+mu T^L=(1+rho T)^L for its unique Frobenius root rho. It has only one geometric root. Thus if Rbar has at least two distinct poles, multiplying by this factor leaves at least one pole, so B is nonpolynomial and the theorem applies. If Rbar has one pole but it is not entirely canceled, the theorem also applies.

The proof does NOT exclude a polynomial Rbar, or a critical constant whose factor1+mu T^L cancels the entire reduced denominator of a one-pole Rbar. Those cases require higher native arithmetic or another argument. No assertion that every phase has a suitable pole structure or noncritical constant valuation is made. A smaller actual support still needs its smaller window.

The complementary noncritical transfer is separately elementary: when0<delta<1 and delta is not1-2^(-a), a unique factorial-cost minimizer gives the period2^a test Rbar_(m mod2^a); an existing w-window for the Rbar coefficients, with Rbar0=1, transfers by circular-prefix reasoning. For delta>=1 the exponential factor exp(kT) reduces to1. Those facts motivate the endpoint theorem but are not used as unproved premises in Sections1--4.
