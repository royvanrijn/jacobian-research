# Polynomial fibers of a hyperbolic inverse coordinate

Date: 19 November 2026. This gives an explicit all-order two-pair Image construction and an exact classification within the specified rational inverse-coordinate fiber. The converse is restricted to that fiber.

## 1. Parameters and polynomial phases

Let K be any characteristic-zero field, lambda in K*, N>=1 an integer, a>=1 odd, and ell>=1 odd. Put

    U=1+lambda*x,
    S(U)=(U^(2N)-1)/(U-1)=1+U+...+U^(2N-1),
    X0=U*S(U)/(2N), C1=U^(N+1), A=C1/X0,
    u_N=x/A=(U^N-U^(-N))/(2N*lambda),
    q=2N(ell-1),
    Q(U)=(1+U^(2N*ell))/(1+U^(2N))
        =sum_(j=0)^(ell-1) (-1)^j U^(2Nj).              (F1)

The quotient is a literal polynomial, Q(1)=Q(0)=1 and degQ=q. Define

    X=X0*Q,
    C=C1*X0^(a-1)*Q^a,
    D=-C1*X0^(a+1)*Q^(a+2),
    p=t^a/x*(C+t^2 D)=(A/x)*(Xt)^a*(1-(Xt)^2).         (F2)

C,D are polynomials with C(x=0)=1,D(x=0)=-1, so p is regular in K[x,y,t]/(t^2+2xy-1). Its exact least sphere/balanced degree is

    d=(2ell(a+2)-1)N+a+2.                              (F3)

**Theorem HyperbolicCoordinateFiberImage.** Every phase (F2) is all-pure. The one fixed marker G=xQ has a nonzero marked moment at EVERY positive even exponent. There is an explicit radial-primitive non-proper-power original two-pair polynomial lift F of balanced degree d, with one fixed original marker G_h=ZQ_h of balanced degree q+1, such that

    E2(F^m)=0 for every m>=1,
    E2(G_h F^(2r))=(2rd+q+2)! (2N)^((a+1)2r)
                   *kappa_(a,2r)*mu_(N,lambda,r) !=0   (F4)

for every r>=1. Odd marked moments vanish. Here

    kappa_(a,m)=int_0^1 t^(am)(1-t^2)^m dt
               =2^m*m!/prod_(j=0)^m(am+2j+1)           (m even),
    mu_(N,lambda,r)=lambda^(-1)*(-N^2 lambda^2)^r
      *(1/(2N))_r*(1-1/(2N))_r / ((1/2)_r*r!).         (F5)

Every factor is nonzero over every characteristic-zero field. This theorem includes prime balanced degrees37 (N=2,a=1,ell=3),61 (2,1,5),71 (4,1,3), and89 (2,5,3). The ell=1 subfamily is exactly the already retained hyperbolic family. Degree37 lies outside that earlier general-a degree set: its 2d-1=73 has no factorization into two integers at least5. No minimum-prime or priority claim is made.

## 2. The all-order endpoint criterion

For any regular unit-center odd-a phase p=t^a/x*(C+t^2D), choose the formal square root X=sqrt(-D/C), X(0)=1, and define A=C/X^a and u=x/A. Let v(u) be the inverse branch x and let G=G(x) be a polynomial in K[x]. Formal endpoint contact gives, for even m,

    Lambda(G p^m)=kappa_(a,m)*[u^m]G(v(u))*Phi(u),
    Phi(u)=v'(u)/(A(v(u))*X(v(u))).                    (F6)

The same holds for G=1. Odd m vanish by t parity. To see (F6), substitute w=Xt in the height integral. The numerator w^(am)(1-w^2)^m has zeros of order m at w=+-1; replacing the actual endpoints +-X by +-1 changes the power series only in order at least m+1. Taking the x^m coefficient therefore gives kappa*[x^m]G A^m/X. Formal residue change u=x/A yields (F6). This argument works over any characteristic-zero field and does not assume convergence.

For the old coordinate u_N one has the retained identity

    Phi0=v'(u)/C1(v(u))=2/(1+U(v(u))^(2N)),
    U(v(-u))=1/U(v(u)),
    [u^(2r)]v(u)*Phi0(u)=mu_(N,lambda,r).              (F7)

The last formula and its nonvanishing follow from the self-contained ODE recurrence in [HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md](HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md). Since A is unchanged and X=X0Q, the new factor is

    Phi=Phi0/Q=2/(1+U^(2N*ell)).                        (F8)

It satisfies Phi(u)+Phi(-u)=2, proving every positive pure coefficient zero. For the fixed marker G=xQ, G(v)Phi=v Phi0, so (F7) proves the nonzero every-even marked coefficient, independent of a and ell. This is an all-order proof, not finite vanishing-prefix evidence.

The unique t-degree<=1 normal form is

    p=t*(1-2xy)^((a-1)/2)*((C+D)/x-2yD).

The D summand has top degree degD+a+1, and the C summand has smaller degree; their highest y-degrees are different. Reduction by t^2=1-2xy never increases total degree, so this realizes the exact least degree. Counting gives

    degD=(2a+3)N+1+(a+2)q,
    d=degD+a+1=(2ell(a+2)-1)N+a+2.

## 3. Exact polynomial fiber classification

**Theorem FixedHyperbolicCoordinatePolynomialFiber.** Fix N,a,lambda as above. Let C,D in K[x] be any regular unit-center odd-a phase, put X=sqrt(-D/C), X(0)=1, and assume its formal coordinate x X^a/C equals the specified rational u_N. Then all positive moments vanish if and only if C,D are exactly (F2) for a unique positive odd ell.

Proof. The fixed coordinate is rational. Since a is odd,

    X=(u_N*C/x)/(X^2)^((a-1)/2)

is rational in x. Let Q=X/X0. Equality of the inverse coordinates gives

    C=C0,a Q^a, D=D0,a Q^(a+2),
    C0,a=C1 X0^(a-1)=U^(N+a)S(U)^(a-1)/(2N)^(a-1).

Every nonzero root of S is simple and has multiplicity a-1 in C0,a. A pole of Q there would have multiplicity at least a in Q^a and could not be canceled. At every other nonzero finite U-point, C0,a is a unit. Thus polynomial C forces Q to have no finite poles except U=0: Q lies in K[U,U^(-1)]. Also Q(1)=1.

Let M=(1+U^(2N))Q, a Laurent polynomial with M(1)=2. By (F6), purity is precisely

    1/M(U)+1/M(1/U)=1,
    (1-M(U))*(1-M(1/U))=1.                            (F9)

The units of K[U,U^(-1)] are exactly c U^k. Hence 1-M=cU^k; evaluating at U=1 gives c=-1. It follows that

    Q=(1+U^k)/(1+U^(2N)).                             (F10)

The case k=0 is not Laurent polynomial. For k!=0, Laurent polynomiality requires |k|=2N*ell with ell a positive odd integer. This is seen by evaluating at a primitive 4N-th root of unity, or by elementary binomial divisibility. For k>0 we recover (F1).

For k<0, Q=U^(-2N*ell)Q_positive, so

    ord_(U=0) C=N+a-2Na*ell.

This is negative except the sole case N=a=ell=1. There C=1, but D=-(U+1)^2/(4U^2), not polynomial. Thus the negative case is excluded. The positive cases were already constructed and proved pure. Their distinct degrees make ell unique. This proves the complete polynomial fiber statement without an all-height rational-descent theorem: rationality was an explicit fixed-coordinate assumption.

## 4. Literal original lift, primitivity and nonpower

Use R=xi1*z1+xi2*z2, Z=xi1*z2, W=2xi2*z1, T=xi1*z1-xi2*z2, with T^2+2ZW=R^2. Put U_h=R+lambda Z, and let Y_h and Q_h be the degree2N and degreeq homogenizations of U*S(U) and Q(U) in R,Z. In particular Y_h(R,0)=2N R^(2N), Q_h(R,0)=R^q. Define

    B=8N^2 W R^(4N*ell)
       -T^2*(Y_h^2 Q_h^2-4N^2 R^(4N*ell))/Z,
    F=T^a U_h^(N+1)(Y_h Q_h)^(a-1) Q_h B.            (F11)

The divided difference is literal polynomial. The dehomogenization is exactly (2N)^(a+1) p, which gives the factorial and scalar in (F4) with G_h=ZQ_h. At the original R=0 point (xi1,xi2,z1,z2)=(1,1,1,-1),

    F=(-1)^(N+1) 2^(a+2)
      *lambda^((2ell(a+2)-1)N+1) !=0.                  (F12)

Hence F is radial-primitive and its degree d cannot be reduced by radial homogenization. The exact T multiplicity is a: modulo T the first B term and all displayed factors are nonzero, using the independence of R,Z on the prime hypersurface T=0.

To exclude EVERY scalar proper power also when a>1, specialize the original variables to (xi1,xi2,z1,z2)=(1,0,1,z), and set U=1+lambda z. Then R=T=1,W=0 and

    F_special=-lambda U^(N+1) (Y(U)Q(U))^(a-1) Q(U)
                 *((Y(U)Q(U))^2-4N^2)/(U-1).          (F13)

For U=-v, v>1, Y(-v)=v(v^(2N)-1)/(v+1) is positive and strictly increasing from0 to infinity. Q(-v)=(1+w^ell)/(1+w), w=v^(2N)>1, is positive and nondecreasing. Its derivative in w has numerator (ell-1)w^ell+ell*w^(ell-1)-1, positive for ell>1 and zero for ell=1. Thus Y(-v)Q(-v)=2N has a simple real algebraic root v>1. It avoids all other factors in (F13), and produces a simple zero of F_special. This rules out every proper power over field extensions. All polynomials defining that root have rational coefficients independent of lambda, so the algebraic simple root transports to the algebraic closure of every characteristic-zero field; no ordering assumption on K is used.

## 5. Explicit prime degree37

Take N=2,a=1,ell=3,lambda=1. Let

    Q_h=R^8-(R+Z)^4 R^4+(R+Z)^8,
    Y_h=4R^4+10R^3Z+10R^2Z^2+5RZ^3+Z^4.

Then

    F37=T(R+Z)^3 Q_h *
       [32W R^24-T^2*(Y_h^2 Q_h^2-16R^24)/Z],
    G37=Z Q_h.                                        (F14)

These have balanced degrees37 and9. F37 has746 original monomials, R=0 value-8, and exact T multiplicity1. The exact all-order statements are

    E2(F37^m)=0 (all m>=1),
    E2(G37 F37^(2r))=(74r+10)! 16^(2r)
                        *kappa_(1,2r)*mu_(2,1,r) !=0.

At r=1 the normalized marked value is -1024/35. Direct original four-variable checks at m=1,2 agree; these are controls, while Sections2--4 prove the universal assertions.

## 6. Antecedents and limits

The elementary hyperbolic-coordinate identity and ODE nonvanishing formula are already retained in [HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md](HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md); the endpoint extension to odd starting powers is proved in Section2. Classical DLMF15.4.13 and4.37.16 remain special-function antecedents, not new identities. The Laurent-unit converse is self-contained. No global prime-degree obstruction is asserted; previously unsettled prime degrees such as5,7,11,17,19 remain outside the supplied examples. No minimum prime degree or historical-priority assertion is made.


The [source ledger](hyperbolic_coordinate_fiber_image_source_ledger.json) pins the exact retained coordinate, ODE and original-contraction inputs and records the classical special-function antecedents. The [checker](verify_hyperbolic_coordinate_fiber_image.py) and [receipt](hyperbolic_coordinate_fiber_image_verification.json) bind the finite controls to these proof and source bytes.

```bash
python3 research/foundations/extended-geometry/verify_hyperbolic_coordinate_fiber_image.py \
  --write-receipt research/foundations/extended-geometry/hyperbolic_coordinate_fiber_image_verification.json
```
