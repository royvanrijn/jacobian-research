# A native bound in the L8 cube branch and its twelfth-grid exclusion

This note establishes an all-start six-window at constant valuation7/8 for the cubic native-residue sector on coefficient value grid(1/12)Z. It does not close radius7/8 on arbitrary coefficient fields, polynomial residue, other critical radii, or unrestricted quadratic SFC.

Let E be a characteristic-zero real-valued nonarchimedean coefficient field, v(2)=1, and let

    q=alpha X^2/2+beta XY+gamma Y^2/2+X+eY,
    alphabar=ebar=0, betabar=gammabar=1,
    R_j=Lambda(q^j)/j!, R_0=1, f=q+k, v(k)=7/8.           (1)

The constant may lie in a real-valued extension. Write O for the valuation ring after adjoining algebraic elements as required. The native residue is (1+T)^(-3). Congruences modulo 2O are in its possibly nonreduced characteristic-two quotient.

The first input is the complete [necessary L8 reduction](BINARY_QUADRATIC_L8_EXACT_START_AND_NATIVE_PROFILES.md). Its native integrality and universal modulo-two formula come from the [native coefficient proof](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md). The final width argument uses the [finite least-exponent-width proof](QUADRATIC_CONSTANT_APPELL_MULTIPLICITIES_FROM_DISCRETE_VALUE_GRIDS.md); its needed facts are restated below.

**Theorem 1 (universal s10 native bound).** Under(1), if six factorial moments vanish at n,...,n+5 with n=10 modulo16, then

    v(alpha)>=9/49.                                     (2)

No coefficient-grid hypothesis is imposed in this assertion.

**Theorem 2 (cube-residue twelfth-grid six-window).** If additionally v(E*) is contained in (1/12)Z, f has no six consecutive zero factorial moments at any positive start. The constant need not belong to E. This is not a smaller sparse actual-support assertion.

## 1. The actual s10 data and exact integral parameters

The imported necessary reduction proves that a hypothetical block starts at6 or10 modulo16. At6, a=1+e+beta belongs to the NONCONSTANT coefficient field and has value7/8. That branch is therefore impossible on grid1/12. Starts7,8,9 are impossible without any grid condition.

For the remaining branch write n=16H+10. The necessary s10 data are

    v(e),v(a)>=1/2, a=1+e+beta,
    D(T)=(1+alpha T+T^2)(1+gamma T),
    D divides 1+UT^8 modulo2O, U=k^8/128,
    alpha^7=0, 1+U(1+alpha^4+alpha^6)=0 modulo2O,
    U=gamma^8 modulo2O.                                 (3)

In particular v(alpha)>=1/7 and v(gamma-1)>=1/14. Adjoin pi with pi^56=2. Every actual phase in this branch therefore admits EXACT integral parameters

    alpha=pi^8*A, gamma=1+pi^4*tau,
    e=pi^28*D0, a=pi^28*B,
    beta=a-1-e,
    k=pi^49*gamma+2chi,                                  (4)

where A,tau,D0,B,chi are integral. They need not be units. The last equation is the imported true-constant parameterization: apply its pi8-root statement to pi^7. No leading coefficient, constant error or native equation is silently replaced by a residue-field value.

## 2. Uniform all-H modulo-four small-moment targets

Let

    E_s=Lambda((q+k)^s)/s!, a_j=k^j/j!, h=H modulo2.

The actual s10 divisibility in(3), together with the universal native formula R=N/D where N=1+aT+eT^2, gives

    R_(r+8)=U R_r modulo2O for EVERY r>=0.                (5)

Indeed (1+UT^8)R=N(1+UT^8)/D is a polynomial of degree at most7: the divisor D has degree3 with unit leading coefficient gamma, so the quotient has degree5. This also treats nilpotent coefficients of N correctly.

We prove the necessary targets

    F_s:=E_s-2hU E_(s-8) in4O, s=10,...,15.              (6)

Divide the original zero equation at m=16H+s by its r=s coefficient. For the head r=s-j, j=0,...,s, the coefficient is

    a_j / binom(16H+j,j).

All a_j for j<=15 are integral. The binomial is the product of 1+16H/i for i=1,...,j. Modulo4 only i=8 can contribute a nontrivial factor: the binomial and its inverse are1 for j<8 and (-1)^H for j>=8. Thus the head equals E_s-2hU E_(s-8) modulo4O, using a_(j+8)=U a_j modulo2O for j<=7. This last relation follows directly from their finite factorial coefficients.

For the first tail r=s+t, t=1,...,16, the normalized coefficient is (16H)_t k^(-t). If H>0,

    v_2((16H)_t(16-t)!)=15+v_2(H).                       (7)

To see this, v_2(16H-i)=v_2(i) for1<=i<16, and the binary digit counts of t-1 and16-t sum to4. If H is even the entire first tail vanishes modulo4. If H is odd its coefficient is

    2U^(-2) a_(16-t) modulo4O.

Pair w=16-t with w+8. Relations a_(w+8)=U a_w and(5) make the two terms cancel modulo2O, so the first tail is zero modulo4O in either parity. If H=0 no tail exists.

Every further tail has polynomial exponent j=16b+c with b<=H-2 and0<=c<16. Its valuation cost relative to the r=s coefficient is

    [f(b)-f(H)]+g(c),
    f(b)=s_2(b)-2b, g(c)=s_2(c)-c/8>=0.

The function f drops by at least one at each increment, so this cost is at least2. All further tails vanish modulo4O. This proves(6) at every H, not only for two sampled starts.

## 3. One finite five-term native identity

Substitute the free integral parameters(4) into the actual native coefficients. In the ring with pi^56=2, modulo pi^72, the combination F_15+gamma F_14 is divisible by2 and obeys

    (F_15+gamma F_14)/2
      = A^7+pi^4*(A^4+tau^8)*B
        +pi^12*A*(A^4+tau^8)*D0 modulo pi^16 O.           (8)

The same identity holds for h=0 and1. Individual free F_s are NOT asserted divisible by2; their lower-degree terms must be retained. Only the displayed combination is divided.

Complete finite inputs for (8) are given here to make it an exact polynomial certificate. For r=0,...,15 use the integer native formula

    R_r=sum_(i+j+g+l+m=r)
      [(2i+j+l)!*(2g+j+m)! /
       (2^(i+g)*i!*j!*g!*l!*m!)]
       *alpha^i*beta^j*gamma^g*e^m.                      (9)

There are 15504 integer terms in total. The factorial ratio is integral by divided-block counting. For each j<=15 write j!=2^v_j o_j with o_j odd. Expand the TRUE k before truncation:

    a_j=(1/o_j) sum_(t=0)^j binom(j,t)
         pi^(49j+7t-56v_j)*gamma^(j-t)*chi^t.             (10)

All exponents in(10) are nonnegative. Odd denominators may be inverted modulo4; powers of two must be removed BEFORE discarding high pi degrees. The surviving finite table modulo pi^72 is

    a0=1,
    a1=pi49+pi53*tau+pi56*chi,
    a2=pi42+pi50*tau^2+pi56*chi^2,
    a3=0,
    a4=pi28+pi44*tau^4+pi56*chi^4,
    a5=0, a6=pi70, a7=0,
    a8=1+pi32*tau^8+pi56*(1+chi^8),
    a9=a1, a10=a2, a11=0,
    a12=a4+pi60*tau^8, a13=0, a14=pi70, a15=0.            (11)

Here pi49 means pi^49, and so on. In particular a6 and a14, each of valuation5/4, are retained. Also U=315a8, so the true pi56*chi8 correction to U is present.

Use E_s=sum_j a_j R_(s-j), subtract 2hU E_(s-8), multiply F_14 by gamma=1+pi4*tau, and add F_15. Collecting all pi degrees below72 leaves precisely the following five terms:

    pi56*A7,
    pi60*A4*B, pi60*tau8*B,
    pi68*A5*D0, pi68*A*tau8*D0.

This proves(8). Equations(4),(9),(10) and these five complete residual terms are a finite universal integer-polynomial identity, not a parameter sample or an extrapolation in r or H. A direct independent stdlib reconstruction of all15504 coefficients and both h identities is possible from these formulas; no computer algebra theorem is imported.

## 4. The bound on the actual nonconstant coefficient

Put L0=A^4+tau^8. The constant remainder in(3) gives, modulo2O=pi56 O,

    pi32*L0+pi48*A6=0.

Indeed gamma8=1+pi32*tau8 modulo2O, alpha4=pi32*A4, alpha6=pi48*A6, and the cross term has pi degree64 or more. Dividing by pi32 shows

    L0+pi16*A6 in pi24 O; hence L0 in pi16 O.              (12)

The last two terms in(8) now lie in pi16 O. Since the actual targets(6) give (F_15+gamma F_14)/2 in2O, equation(8) implies A7 in pi16 O. Therefore

    v(A)>=16/(7*56)=2/49,
    v(alpha)=1/7+v(A)>=9/49.

This proves Theorem1, including alpha=0.

On grid1/12, the next permissible value at least9/49 is1/4. Thus v(alpha)>=1/4. The imported s10 remainder then gives v(gamma-1)>=1/8. We also retain v(e),v(a),v(beta-1)>=1/2 from(3).

## 5. Uniform native coefficients at indices three modulo four

These bounds force

    v(R_r)>=1/4 for EVERY positive r=3 modulo4.           (13)

Here is a coefficientwise integer proof, requiring no bound on r. At the cube point alpha=e=0,beta=gamma=d=1, the universal integer native coefficients have residue series (1+T)^(-3), whose r=3mod4 coefficients vanish. Thus the constant Taylor coefficient has valuation at least1. The linear alpha,e,beta-1 contributions have valuation at least1/4. For z=gamma-1, the derivative of the native series at this point, modulo2, is

    T/(1+T)^4.

Its coefficient at r=3mod4 is zero, so the corresponding integer derivative is even; multiplying by z is harmless. All Taylor terms of total degree at least2 have valuation at least2*(1/8)=1/4. All coefficients in this finite polynomial Taylor expansion are integers. This proves(13) for each r uniformly and retains the difference between an even universal coefficient and a merely nonzero small valuation.

## 6. The final finite-width contradiction on grid1/12

Assume the s10 block persists. Put N=n+5=16H+15, and consider the monic Appell polynomial

    P_N(Z)=sum_r (N)_r R_r Z^(N-r).

The root k has multiplicity at least6 by differentiating with respect to Z. Its minimal polynomial F over E has positive least-cost width at delta=7/8. Coefficients of F have values in (1/12)Z, so any two least exponents differ by a multiple of2. Widths add under multiplication by scaling Z=kV and multiplying nonzero polynomials over the maximal-ideal residue FIELD. Characteristic-zero separability gives F^6 dividing P_N. Hence width(P_N)>=12. These are finite exponent-width statements; no completeness or Newton root count is assumed.

Write N=8q+t with q=2H+1 odd and t=7. For the polynomial exponent j=N-r=8b+c, the coefficient cost after subtracting its common base is

    [s_2(b)-b]-[s_2(q)-q]+s_2(c)-c/8+v(R_(N-j)).

The unit R5 at j=8q+2 gives minimum M<=3/4. Only b=q,q-1 can minimize: every earlier b costs at least1. The only possible low indices c are0,2,4; c1 costs7/8 and all other c cost at least5/4. Thus the available exponent span is at most12. The width lower bound forces equality, with least endpoints j=8(q-1) and j=8q+4. At the upper endpoint r=3, so(13) makes its cost at least1/2+1/4=3/4. Therefore M=3/4.

All inequalities are saturated: k has multiplicity exactly6, width(F)=2 and the remaining factor has width zero. Every least exponent of F differs by a multiple of2, so its residue polynomial is a monomial times a binomial with exponent gap2. The residue polynomial of P_N is consequently a monomial times its sixth power. In characteristic two its only nonzero offsets are0,4,8,12.

But j=8(q-1)+2 has r=13, a unit native coefficient in the cubic residue. It has cost3/4=M, and its offset from the smallest least exponent is2. This nonzero least coefficient is forbidden by the sixth-power support, a contradiction.

This excludes s10 on grid1/12. Together with the imported s6 obstruction8|e0 and the grid-free exclusions of starts7,8,9, it proves Theorem2 for all positive starts. It does not assert a cube-residue six-window at radius7/8 on arbitrary finer grids. Nor does this proof treat polynomial residues or radii1-1/L with L>=16. All actual-support consequences must respect the actual number of monomials; no sparse deletion rule is used.
