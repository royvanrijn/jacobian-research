# Exact native reductions at the next critical radius in the cubic residue sector

This note proves necessary conditions at v(k)=7/8. It does not prove a six-window at that radius. In particular two surviving start branches remain; their necessary profiles are not examples of zero moments.

Let E be a characteristic-zero real-valued nonarchimedean field, v(2)=1. Let

    q=alpha X^2/2+beta XY+gamma Y^2/2+X+eY,
    alphabar=ebar=0, betabar=gammabar=1,
    R_j=Lambda(q^j)/j!, R_0=1, f=q+k, v(k)=7/8.           (1)

Natural coefficients belong to O_E; k may belong to a valued extension. Bars denote maximal-ideal residue. The native residue series is Rbar=(1+T)^(-3). All coefficient congruences below, however, are in A=O/2O, which can contain nilpotents.

The inputs are the universal integer-coefficient [native identity and factorial-integrality proof](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md), and the native [critical detector and six-window theorem](BINARY_QUADRATIC_CRITICAL_DYADIC_SIX_MOMENT_WINDOWS.md). The new reductions below use actual normalized moment equations and finite polynomial degrees, not arbitrary recurrence solutions.

**Proposition.** If six moments of f vanish at n,...,n+5, then n modulo16 is either6 or10. Write

    lambda=k^2/2, kappa=k^4/8, U=k^8/128,
    a=1+e+beta, C=R_1=alpha+gamma+a.

Here v(lambda)=3/4, v(kappa)=1/2, and U is a unit with residue one. All assertions below are necessary for a hypothetical block.

* If n=6 modulo16, then in A

      e=lambda, a=k, beta=1+lambda+k,
      gamma=C+alpha+k,
      D(T)P(T)=1+kappa*T^4+U*T^8,                       (2)

  where D=(1+alpha*T+T^2)(1+gamma*T) and P is the actual small-moment prefix of degree at most5. Moreover

      alpha^3*(kappa+U*alpha^4)=0 in A,
      v(alpha)=1/8 or v(alpha)>=1/6.                    (3)

  In this branch v(a)=7/8, so every discrete grid bound (1/e0)Z for the NONCONSTANT coefficient field has 8 dividing e0. In particular this branch is excluded on any grid with 2-part four.

* If n=10 modulo16, then v(e),v(a)>=1/2 and, in A[T],

      D divides 1+U*T^8,
      U=gamma^8,
      alpha^7=0,
      1+U*(1+alpha^4+alpha^6)=0.                        (4)

  Consequently v(alpha)>=1/7 and v(gamma-1)>=1/14. More precisely, if v(alpha)<1/4 then v(gamma-1)=v(alpha)/2; if v(alpha)>=1/4 then v(gamma-1)>=1/8. After adjoining pi with pi^8=2 one may write the TRUE constant as

      k=pi^7*gamma+2chi,  chi integral.                 (5)

No assertion that either branch occurs is made.

## 1. Which start residues are even eligible?

The critical detector for L=8 forces (1+mu*T^8)Rbar to be polynomial. Since the only geometric pole of Rbar is1, mu=1. Thus U has residue one and the detector polynomial is

    (1+T^8)/(1+T)^3=(1+T)^5=1+T+T^4+T^5.

The relevant detector is its prefix of length16, read circularly. Its nonzero entries at0,1,4,5 show that a six-zero block can begin only at residues6,7,8,9,10. This is a necessary leading detector, not a converse to actual zeros.

## 2. Every eligible actual block gives six small native equations

For i=0,...,15 put E_i=Lambda(f^i)/i!. In A, the finite exponential prefix is

    sum_(j=0)^15 (k^j/j!)*T^j
      =(1+kT+lambda*T^2+kappa*T^4)(1+U*T^8).             (6)

Indeed the valuations j*(7/8)-v_2(j!)=s_2(j)-j/8 are nonnegative for j<=15. Exactly the degrees0,1,2,4,8,9,10,12 survive modulo 2O. Their odd integer factors reduce to one. All products of two of k,lambda,kappa vanish in A.

Suppose n=16H+s with s in6,...,10. For each i=s,...,s+5, divide the original zero-moment equation at 16H+i by its nonzero r=i coefficient. For r=i-j<=i the ratio to k^j/j! is the inverse of binom(16H+j,j), an odd integer. It is one modulo2, and all retained k^j/j! are integral. Terms with r>i vanish modulo2: in their polynomial exponent j=16b+c<16H, write

    s_2(j)-j/8=[s_2(b)-2b]+[s_2(c)-c/8], 0<=c<16.

The second bracket is nonnegative. The first bracket decreases by at least one whenever b increases by one, since s_2(b+1)-s_2(b)<=1. Its gap from b=H-1 to H is 1+v_2(H)>=1. Hence every such tail costs at least one valuation unit relative to the r=i coefficient. If H=0 no tail exists. This proves the ACTUAL necessary equations

    E_s=...=E_(s+5)=0 in A.                              (7)

No parity of H enters this first precision.

## 3. The exact polynomial degree constraints

In A[[T]], universal native congruence gives

    R=N/D,
    N=1+aT+eT^2,
    D=(1+alpha*T+T^2)(1+gamma*T+e^2*T^2).

Put B=1+kT+lambda*T^2+kappa*T^4, and let P be the actual prefix sum_(j=0)^(s-1) E_j T^j. Equations(6),(7) give

    D P=B N(1+UT^8) modulo T^(s+6).                      (8)

Initially deg(DP)<=s+3. The coefficients of BN needed below are

    (BN)_1=a+k,
    (BN)_2=e+ka+lambda,
    (BN)_3=ke+lambda*a,
    (BN)_4=lambda*e+kappa,
    (BN)_5=kappa*a,
    (BN)_6=kappa*e.                                    (9)

The coefficient (BN)_4 cannot vanish in A when e has positive valuation: kappa has valuation1/2 and lambda*e has valuation greater than3/4.

If s=7, the T^12 coefficient of(8) must vanish because deg(DP)<=10; it is U(BN)_4, a contradiction. If s=8 the same coefficient contradicts deg(DP)<=11.

If s=9, the coefficients at13 and14 first force kappa*a=kappa*e=0, so v(a),v(e)>=1/2. Then e^2=0 in A and degD<=3; now deg(DP)<=11. Its T^12 coefficient again contradicts the nonzero U(BN)_4. Thus s=9 is also impossible.

## 4. The branch s=6

Here deg(DP)<=9 and (8) is valid through degree11. The coefficients10,11 force

    e+ka+lambda=0, ke+lambda*a=0 in A.

Substitute the first equation into the second. Since k lambda=k^2=0 in A, it follows that lambda*a=0, so v(a)>=1/4. Hence ka=0 and e=lambda in A. Consequently e^2=0, degD<=3, and deg(DP)<=8. Coefficient9 of(8) now forces a+k=0. This proves the first parameter identities in(2). All BN coefficients except0,4 vanish, and (8) becomes the EXACT polynomial identity in(2): its left degree is at most8, while its right agrees through degree11.

Reduce the identity modulo the monic factor A_x=1+alpha*T+T^2. In its polynomial quotient,

    T^4=1+alpha^2+alpha^3*T,
    T^8=1+alpha^4+alpha^6+alpha^7*T.

The coefficient of T in the remainder is alpha^3(kappa+U alpha^4), proving (3). Write x=v(alpha)>0, with x=+infinity allowed. If x<1/8, this product has valuation7x<1, impossible. If 1/8<x<1/6, it has valuation3x+1/2<1, again impossible. The two ranges exhaust the complement of x=1/8 or x>=1/6.

Finally a belongs to the normalized NONCONSTANT coefficient field E. Since a=k modulo2O and v(k)=7/8<1, it has exactly this valuation. Any bound v(E*) subset(1/e0)Z must therefore have8 dividing e0. The constant's field need not itself equal E for this conclusion.

## 5. The branch s=10

Here deg(DP)<=13 and (8) holds through degree15. Coefficient14 forces kappa*e=0, hence v(e)>=1/2. Therefore e^2=0, degD<=3, and deg(DP)<=12. Coefficient13 then forces kappa*a=0, giving v(a)>=1/2.

The product BN now reduces to

    S=1+(a+k)T+(e+lambda)T^2+kappa*T^4,

and S^2=1 in A[T], since all three nonconstant coefficients have valuation at least1/2. Both sides of (8) now have degree at most12, so

    DP=S(1+UT^8)

is exact. Multiplying by S gives D(PS)=1+UT^8, proving its divisibility by D in (4). Because e^2=0, D=A_x(1+gamma*T), with gamma a unit. Evaluation at T=gamma^(-1) gives U=gamma^8. Monic remainder division by A_x, using the displayed T^8 identity, gives alpha^7=0 and the constant identity in (4).

Thus x=v(alpha)>=1/7. Put z=gamma-1, which initially has positive valuation. The constant identity is equivalent in A to

    z^8+alpha^4+alpha^6+z^8 alpha^4+z^8 alpha^6=0.

If x<1/4, the alpha^4 term has valuation4x<1, all the remaining alpha terms are strictly higher, and cancellation requires 8v(z)=4x. If x>=1/4, it simply gives z^8=0 in A, so v(z)>=1/8. These are the claimed bounds, including alpha=0.

To prove(5), adjoin pi^8=2 and write k=pi^7(gamma+y). Since U=gamma^8 modulo2O, y has positive valuation. If 0<v(y)<1/8, its eighth power is uniquely least in (gamma+y)^8-gamma^8: every other nonconstant binomial coefficient is even. The eighth-power difference would have valuation less than one, a contradiction. Hence v(y)>=1/8 and chi=pi^7*y/2 is integral.

The argument stops at these exact necessary profiles. No six-window at radius7/8 is asserted, no higher-radius induction is used, and no numerical values of n are extrapolated.
