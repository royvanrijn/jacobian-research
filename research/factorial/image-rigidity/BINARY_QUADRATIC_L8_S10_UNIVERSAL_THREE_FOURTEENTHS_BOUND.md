# A universal three-fourteenths bound in the L8 s10 branch

This refinement concerns a NECESSARY coefficient valuation for an actual six-zero block. It proves no unrestricted L8 six-window and does not treat other constant radii or polynomial native residue. In particular no stronger bound1/4 is asserted from the single combined native equation used here.

Use the real-valued characteristic-zero setup

    q=alpha X^2/2+beta XY+gamma Y^2/2+X+eY,
    alphabar=ebar=0, betabar=gammabar=1,
    f=q+k, v(2)=1, v(k)=7/8,
    N1=1+e+beta, z=gamma-1.

The constant may lie in a valued extension of the nonconstant coefficient field. Assume an ACTUAL six-zero factorial-moment block starts at n=16H+10. Put R_j=Lambda(q^j)/j!, E_s=Lambda(f^s)/s!, U=k^8/128, h=H modulo2, and

    F_s=E_s-2hU E_(s-8), s=10,...,15,
    G=F_15+gamma F_14.

**Theorem.** Every such hypothetical block satisfies

    v(alpha)>=3/14.                                     (1)

This is universal in the coefficient ramification; it strengthens the earlier necessary bound9/49. The proof uses the exact native inputs and uniform all-H targets in [BINARY_QUADRATIC_L8_NATIVE_BOUND_AND_TWELFTH_GRID_SIX_WINDOW.md](BINARY_QUADRATIC_L8_NATIVE_BOUND_AND_TWELFTH_GRID_SIX_WINDOW.md), and its [necessary-profile input](BINARY_QUADRATIC_L8_EXACT_START_AND_NATIVE_PROFILES.md). Their scopes are recalled before the new finite identity.

## 1. Actual data and the contradiction range

The imported results give

    a=v(alpha)>=9/49,
    v(e),v(N1)>=1/2,
    U=gamma^8 modulo2O,
    1+U(1+alpha^4+alpha^6)=0 modulo2O,
    F_s in4O for s10,...,15.                            (2)

If alpha=0 or a>=1/4, conclusion(1) is immediate. It therefore suffices to exclude

    9/49<=a<3/14<1/4.

In this range the actual s10 remainder also gives v(z)=a/2. Let

    K=alpha^4+z^8, Q=alpha^2+z^4, S=alpha+z^2.

From the constant remainder in(2), using gamma^8=1+z^8 modulo2O,

    K=alpha^6+z^8 alpha^4+z^8 alpha^6 modulo2O.

Because a>=9/49>1/6, all three terms on the right have valuation greater than one. Consequently

    v(K)>=1,
    v(Q)>=1/2,
    v(S)>=1/4.                                         (3)

The last two inequalities follow from Q^2=K+2alpha^2 z^4 and S^2=Q+2alpha z^2, respectively. These are inequalities, not claims that the remainders are zero as field elements.

## 2. A deeper finite native identity

Adjoin pi with pi^56=2. The actual parameters admit the exact integral form

    alpha=pi^8*A, gamma=1+pi^4*tau,
    e=pi^28*D, N1=pi^28*B,
    beta=N1-1-e, k=pi^49*gamma+2chi.                     (4)

All five parameters A,tau,D,B,chi are integral. In the contradiction range their values are

    v(A)=a-1/7>=0, v(tau)=a/2-1/14>=0,
    v(B),v(D),v(chi)>=0.

The exact k-parameterization, which is needed to keep the TRUE U, follows from U=gamma^8 modulo2O by the previous eighth-power least-term argument. It is not an extra coefficient-field or grid assumption.

Write

    L=A^4+tau^8, Q0=A^2+tau^4, S0=A+tau^2.

A universal integer-polynomial calculation in the ring pi^56=2, modulo pi^84, gives, for BOTH h=0 and1,

    G=alpha^7+K*(N1+alpha*e)+R modulo pi^84 O,            (5)

where the COMPLETE low-degree remainder is

    R=pi^72*Q0
      +pi^76*(A^6*B+tau*Q0)
      +pi^80*S0*(tau^2*A+D^2*Q0)
      +pi^81*L
      +pi^82*A*L.                                      (6)

Equation(6) contains15 ordinary monomials. All other native terms have pi degree at least84 in integral parameters. This identity is more precise than the earlier pi^72 five-term calculation; it cannot be inferred merely by extrapolating that earlier precision.

For a transparent finite certificate, (6) expands as follows:

| pi degree | complete coefficient |
|---:|---|
|72| A^2+tau^4 |
|76| A^6 B+tau A^2+tau^5 |
|80| A^3 D^2+tau^2 A^2+tau^2 A^2 D^2+tau^4 A+tau^4 A D^2+tau^6 D^2 |
|81| A^4+tau^8 |
|82| A^5+A tau^8 |

There are no other pi degrees below84 after subtracting alpha^7+K(N1+alpha e).

The universal native inputs are fixed and finite. For r=0,...,15 use

    R_r=sum_(i+j+g+l+m=r)
      [(2i+j+l)!*(2g+j+m)! /
       (2^(i+g)*i!*j!*g!*l!*m!)]
       *alpha^i*beta^j*gamma^g*e^m.                      (7)

There are15504 integer terms. For j<=15, write j!=2^v_j o_j with o_j odd and expand BEFORE truncation:

    k^j/j!=(1/o_j) sum_(t=0)^j binom(j,t)
      pi^(49j+7t-56v_j)*gamma^(j-t)*chi^t.                (8)

All exponents in(8) are nonnegative. In particular the j6/j14 terms of valuation5/4 and the j5/j13 terms of valuation11/8 must be retained at this deeper precision. Invert odd denominators modulo4, form the actual E_s and F_s, then collect G minus the first two terms of(5). Equations(7),(8), the exact parameter substitutions(4) and the five complete coefficient rows above are a finite universal polynomial identity certificate for(5),(6). Only coefficients divisible by4, or pi degrees at least84, may be discarded; 2=pi^56 is retained. No moment-index or coefficient sample is extrapolated.

## 3. Every remainder term is at least three halves

The five rows in(6) have useful physical meanings. Since

    K=pi^32 L, Q=pi^16 Q0, S=pi^8 S0,

they are exactly the following field expressions:

    2Q,
    alpha^6 N1+2zQ,
    e^2 S Q+2alpha z^2 S,
    pi^49 K,
    pi^42 alpha K.                                    (9)

These identities also check the factors of two in the finite table.

Using(2),(3) and a>1/6, their valuations are bounded below by

    3/2,
    min(6a+1/2, 3/2+a/2)>3/2,
    min(7/4, 5/4+2a)>3/2,
    15/8,
    7/4+a,

respectively. The discarded error in pi^84 O has valuation at least84/56=3/2, because all parameters in(4) are integral. The other displayed term in(5) has

    v(K*(N1+alpha e))>=1+1/2=3/2.

It follows that

    G=alpha^7+an element of valuation at least3/2.

But G belongs to4O by the actual all-H moment targets(2). If a<3/14, its alpha^7 term has valuation7a<3/2 and is uniquely least, a contradiction. This proves(1).

The argument stops here. At a=3/14, terms of valuation3/2 can cancel; the single equation G does not justify a further bound1/4. The result is a necessary universal s10 coefficient bound, not an assertion that such a phase exists, nor a proof of every-radius or arbitrary-grid quadratic SFC.