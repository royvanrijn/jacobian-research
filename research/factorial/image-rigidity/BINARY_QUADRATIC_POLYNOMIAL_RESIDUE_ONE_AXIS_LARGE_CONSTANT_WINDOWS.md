# A polynomial-residue one-axis sector with a large constant

This is a complete valued-sector argument. It concerns actual factorial moments, not arbitrary solutions of an interior/edge recurrence. It does not close general binary-quadratic SFC or the polynomial-residue frontier.

Let K be a characteristic-zero field with a real-valued nonarchimedean valuation v extending v_2, normalized by v(2)=1. Write O={x:v(x)>=0}, and set v(0)=+infinity. For polynomials define Lambda(X^i Y^j)=i!j! and its one-variable version Lambda_1(X^i)=i!.

**Theorem.** Suppose

    q=X^2/2+beta*XY+gamma*Y^2/2+d*X+e*Y,
    v(beta-1)>0,  v(d)>0,
    0<sigma=v(gamma)<v(e),
    f=q+k,  delta=v(k)>=1+sigma.                         (1)

Here d or e may be zero; gamma and k are nonzero. Then the factorial moments Lambda(f^n),...,Lambda(f^(n+4)) cannot all vanish for any n>=1. The coefficient field can have arbitrary ramification. The same conclusion holds after interchanging X and Y and applying the corresponding normalization.

The native residue series of q is polynomial: it is exactly 1. Thus this sector is outside the previously proved nonpolynomial-residue constant-shift criteria. If at least one of d,e is nonzero, f has five or six actual monomials, and the theorem implies its actual-support strong factorial window. When d=e=0 the support is four, and this five-window alone makes no four-window claim.

## 1. Native integrality and a congruence over O/2O

Put R_j=Lambda(q^j)/j!, R_0=1. These coefficients are in O. One elementary proof expands q in divided monomials. For a choice of multiplicities, the resulting integer coefficient counts partitions of labeled X- and Y-slots into indistinguishable nonempty blocks of the prescribed multidegrees. Equivalently the coefficient is

    D_X! D_Y! / ( product_nu a_nu! * product_nu (nu_X! nu_Y!)^a_nu ),

where D_i=sum_nu a_nu*nu_i. First allocate the labeled slots to ordered blocks, then divide by the free permutations of blocks of identical nonempty type. This proves integrality without any analytic interpretation.

We will use the following stronger coefficient congruence, in the possibly nonreduced ring O/2O, not merely in the residue field:

    R(T)=sum_j R_j T^j
       == [1+(beta+d+e)T+de*T^2]
          / [(1+T+d^2*T^2)(1+gamma*T+e^2*T^2)] mod 2O.   (2)

Here and below congruences of formal series are coefficientwise. A direct proof is short. For a one-variable zero-constant phase a*X^2/2+b*X, let F(T) be its divided factorial series. Its exact normalized recurrence is

    F_j=a(2j-1)F_(j-1)+b^2 F_(j-2), j>=2,
    F_0=1, F_1=a+b.

Consequently F==(1+bT)/(1+aT+b^2T^2) modulo 2. The series with one extra X in the factorial functional is, modulo 2, (1+b*d/db)F=1/(1+aT+b^2T^2): on a term with u quadratic and w linear factors the extra factorial factor is 2u+w+1, congruent to w+1.

Now expand the cross term beta*XY. With c cross factors, each axis coefficient is divisible by c!, as follows from

    (c+2u+w)! / [c! 2^u u! w!] in Z.

Thus the total cross contribution has a factor c!; every c>=2 term is zero modulo 2. The c=0 and c=1 terms give (2). This is a universal integer-coefficient calculation, so it is valid over O/2O even for arbitrarily small positive valuations. Reducing (2) in the residue field under (1) gives Rbar=(1+T)/(1+T)=1.

For P_j=Lambda(f^j)/j!, formal addition of the constant gives

    P(T)=exp(kT)R(T),
    P_j=sum_(h=0)^j k^h/h! * R_(j-h).                    (3)

For h>=1, Legendre's formula gives

    v(k^h/h!)=h(delta-1)+s_2(h)>=delta>1.

Therefore P_j-R_j has valuation at least delta, and P(T) is congruent to the rational series in (2) modulo 2O at every coefficient.

## 2. The unit axis and its consecutive-ratio slopes

Let

    A_j=Lambda_1((X^2/2+dX+k)^j)/j!,
    A_j^0=Lambda_1((X^2/2+dX)^j)/j!.

All A_j^0 are units: the preceding one-variable residue series is 1/(1+T). Equation (3) on this axis shows v(A_j-A_j^0)>=delta>1, so all A_j are also units.

Define r_i=A_i/A_(i-1), r_i^0=A_i^0/A_(i-1)^0 and

    w_i=r_(i+1)-r_i,   w_i^0=r_(i+1)^0-r_i^0,  i>=1.

The two ratios, and the two slopes, differ by quantities of valuation greater than 1. The exact zero-constant recurrence yields

    r_1^0=1+d,
    w_1^0=2-d/(1+d),
    w_(i+1)^0=2-d^2*w_i^0/(r_i^0*r_(i+1)^0).            (4)

Two consequences will be needed, both including possible cancellations at valuation 1:

* If v(w_i)>=1, then v(w_(i+1))=1.
* If v(w_i)<1, then d!=0 and, writing lambda=v(d)>0,

        v(w_i)=(2i-1)lambda.                            (5)

For the first assertion, the second term in (4) has valuation strictly greater than 1. For the second, trace (4) backwards: a slope below 1 cannot follow a slope at least 1, and while slopes stay below 1 each forward step adds 2lambda. The first slope below 1 starts with valuation lambda. Errors of valuation greater than 1 do not affect either assertion. If d=0 all zero-constant slopes have valuation 1, so the second case cannot occur.

## 3. The small axis and the exact interior forcing

Let

    B_j=Lambda_1((gamma*Y^2/2+eY+k)^j)/j!.

Dividing its phase by gamma gives natural quadratic coefficient 1, linear coefficient of positive valuation v(e)-sigma, and constant of valuation delta-sigma>=1. Its zero-constant divided moments are units, and its constant exponential changes them only by nonunits. Hence

    v(B_j)=j*sigma.

Since gamma-beta is a unit, the element

    rho=(gamma*d-beta*e)/(gamma-beta)

has valuation strictly greater than sigma. Set C_i=B_i+rho*B_(i-1). For every i>=1,

    v(C_i)=i*sigma.                                    (6)

Here is the interior identity, with all signs and index shifts specified. Put

    D=1+gamma-2beta,   Delta=gamma-beta^2,
    V=gamma*d^2-2beta*d*e+e^2.

For i>=2 the native normalized moments satisfy

    D P_i=2Delta*i*P_(i-1)+(V-2kDelta)P_(i-2)
       +(gamma-beta)B_i+(gamma*d-beta*e)B_(i-1)
       +(1-beta)A_i+(e-beta*d)A_(i-1).                  (7)

A formal derivation avoids analytic integrals. Let X(T)=Lambda(X exp(Tf)), Y(T)=Lambda(Y exp(Tf)), with P,A,B the corresponding formal series. The monomial identities Lambda(g_X)=Lambda(g)-Lambda(g(0,Y)) and Lambda(Xg)=Lambda(g)+Lambda(Xg_X) imply

    P-B=T(X+beta*Y+dP),
    P-A=T(beta*X+gamma*Y+eP),
    X+Y=2P+2TP'-T(dX+eY)-2kTP.

Apply the adjugate of [[1,beta],[beta,gamma]] to the first two equations. Substitution into the third cancels the linear-in-T multiple of P and gives

    D P=2Delta*T*P+2Delta*T^2*P'+(V-2kDelta)T^2*P
       +(gamma-beta)B+(1-beta)A
       +T[(gamma*d-beta*e)B+(e-beta*d)A].

Taking the coefficient of T^i proves (7). In particular the forcing sign in V-2kDelta is fixed by this calculation. No division by D or Delta is needed.

Assume now five consecutive P_n,...,P_(n+4) vanish. Put J=n+2>=3. Equation (7), at i=J,J+1,J+2, and division only by the units A_(i-1), gives

    (1-beta)r_i+(e-beta*d)+(gamma-beta)z_i=0,
    z_i=C_i/A_(i-1),  v(z_i)=i*sigma.                  (8)

Subtract the first two equations and then the last two. Since the z_i have strictly increasing valuations, their successive differences have valuations J*sigma and (J+1)*sigma. Thus 1-beta!=0 and

    v(w_(J+1))-v(w_J)=sigma>0,
    v(1-beta)+v(w_J)=J*sigma.                          (9)

The first assertion in Section 2 rules out v(w_J)>=1. Therefore (5) applies and

    d!=0, lambda=v(d)>0,
    (2J-1)lambda=v(w_J)<1,
    v(1-beta)=J*sigma-(2J-1)lambda.                    (10)

In the i=J equation of (8), the last term has valuation J*sigma, strictly larger than v(1-beta). Since r_J is a unit, it follows that

    v(e-beta*d)=v(1-beta).

If sigma>2lambda, then v(e)>sigma>lambda and beta is a unit, so the left side is lambda. Equation (10) would imply sigma=2lambda, a contradiction. We have proved

    sigma<=2lambda,
    (n+1)sigma<1.                                     (11)

Indeed (n+1)sigma=(J-1)sigma<=(2J-2)lambda<(2J-1)lambda<1.

## 4. Removing the unit pole and the final native contradiction

Work in a valued algebraic extension if necessary to factor the two quadratics in (2). Their parameters can be labeled so that

    1+T+d^2T^2=(1+aT)(1+uT),  v(a)=0, v(u)=2lambda,
    1+gamma*T+e^2T^2=(1+bT)(1+cT),
    v(b)=sigma, v(c)>sigma.                            (12)

For example a+u=1, au=d^2; integral roots, their unit sum and their product determine the two valuations. Similarly b+c=gamma and bc=e^2, with v(e)>sigma. If e=0, take c=0. No completeness or convergence of a formal series is required.

Put N_1=beta+d+e and N_2=de. Then N_1 is a unit and v(N_2)=lambda+v(e)>sigma, with the usual convention when e=0. Multiplying the coefficient congruence (2)--(3) by 1+aT gives

    (1+aT)P(T)
      == (1+N_1 T+N_2 T^2)/[(1+uT)(1+bT)(1+cT)] mod2O. (13)

By (11), u has valuation at least sigma; b has valuation sigma; c has strictly greater valuation. Write h_j for the complete homogeneous polynomial of degree j in -u,-b,-c, so these are exactly the denominator inverse coefficients. Every h_j has valuation at least j*sigma. Moreover among every two consecutive h_j,h_(j+1), one attains its lower bound. If v(u)>sigma, only b survives after normalization and every lower bound is attained. If v(u)=sigma, reduce u/b and c/b: the remaining two nonzero unit parameters have complete homogeneous sequence H_j satisfying a second-order recurrence with unit backward coefficient and H_0=1. Two consecutive zero residues would propagate back to H_0=0, impossible.

Apply this assertion to h_n,h_(n+1). For one index i in {n+1,n+2}, the coefficient N_1 h_(i-1) has valuation exactly (i-1)sigma. Both other contributions in the numerator of (13) have larger valuation:

    v(h_i)>=i*sigma>(i-1)sigma,
    v(N_2 h_(i-2))>sigma+(i-2)sigma=(i-1)sigma.

Thus the coefficient of T^i on the right of (13) has valuation exactly (i-1)sigma<1, by (11). On the left it is P_i+aP_(i-1)=0 because the assumed five-zero block contains both indices. This contradicts a congruence modulo 2O. The theorem follows.

## 5. Scope of the result

The argument uses the actual small-axis valuations, the unit-axis consecutive ratios, and the native rational congruence together. It does not use a zero-constant root-preservation assertion. In particular it does not extend automatically to small constants, to v(gamma)>=v(e), to the pure-cross residue profile, or to all finer-grid constant sectors. No finite calculation is a premise of the universal argument. Native finite identity controls may accompany this note only as supplementary checks.