# A seventh-power cancellation proof for the cubic L8 sector

Let E be a real-valued characteristic-zero nonarchimedean field, v(2)=1,
and let

    q=alpha X^2/2+beta XY+gamma Y^2/2+X+eY,
    alphabar=ebar=0, betabar=gammabar=1,
    R_j=Lambda(q^j)/j!, R_0=1, f=q+k, v(k)=7/8.

The natural coefficients are integral and the constant may belong to a
real-valued extension. The native residue is (1+T)^(-3). Put
a=1+e+beta and U=k^8/128. This note uses the exact
[necessary s10 profile](BINARY_QUADRATIC_L8_EXACT_START_AND_NATIVE_PROFILES.md).
Native integrality and the full O/2O identity are supplied by the
[native coefficient proof](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md).
The proof manifest pins these inputs and the optional stronger bound below.

**Theorem 1 (short universal native bound).** Suppose an actual six-zero
factorial-moment block starts at n=10 modulo16. Write x=v(alpha), with
x=+infinity allowed. Then

    v(R_3)>=min(1/3,2x)>=2/7.                            (1)

Only the earlier necessary x>=1/7 is used. Neither the later9/49 nor
the three-fourteenths native calculation is a premise of this theorem.

**Theorem 2 (short all-four-odd-grid L8 proof).** If the nonconstant
coefficient value group is contained in (1/(4u))Z for an odd positive
integer u, then q+k has no six consecutive zero factorial moments at
any positive start. This is an L8 cubic-residue theorem, not an
unrestricted L8 theorem on every coefficient field.

## 1. Cancel a unit-leading factor in the full modulo-two ring

Assume the hypothetical start-ten block. Work in A=O/2O, a
characteristic-two ring that can contain nilpotents. The necessary
profile gives

    v(a),v(e)>=1/2, x>=1/7,
    X(T)=1+alpha T+T^2, Y(T)=1+gamma T,
    D=X Y divides H=1+U T^8,
    U=gamma^8 modulo2O.                                 (2)

In particular e^2=0 in A. Since the characteristic is two,
H=Y^8. There is therefore a polynomial Q with X Y Q=Y^8.

The factor Y is not asserted to be a polynomial unit. Its LEADING
coefficient gamma is a unit, so multiplication by Y is injective on
A[T]: a nonzero polynomial's highest coefficient remains nonzero
after multiplication by that unit. We may thus cancel Y and obtain

    X Q=Y^7 in A[T].                                    (3)

This cancellation is valid in the nonreduced ring; it is not a division
by an arbitrary nonunit or an appeal to coprime factors.

## 2. A two-dimensional quotient gives a cubic relation

Because X is monic quadratic, A[T]/(X) is free on1,T. Set

    eta=1+gamma T, p=alpha*gamma,
    q0=1+gamma^2+alpha*gamma.

The elements1,eta are another free basis since gamma is a unit.
The relation T^2=1+alpha T gives

    eta^2=p*eta+q0.

Successive multiplication, entirely in this free rank-two algebra, gives

    eta^7=(p*q0^3+p^5*q0)
           +(q0^3+p^4*q0+p^6)*eta.                     (4)

Equation(3) makes eta^7 zero in the quotient. Freeness of the basis
therefore gives, for the actual field lifts of p,q0,

    q0^3+p^4*q0+p^6 in2O.                               (5)

This is a finite polynomial identity over A. No root count, completeness
or division by a nilpotent element has entered.

Since gamma is a unit, v(p)=x. If d0=v(q0) were strictly below
min(1/3,2x), its cube would have value3d0<1, whereas

    v(p^4*q0)=4x+d0>3d0,
    v(p^6)=6x>3d0.

The uniquely least term would contradict(5). Thus

    v(q0)>=min(1/3,2x).                                 (6)

This also covers p=0 or q0=0 with the usual infinite-value convention.

## 3. Relate the quotient invariant to the ACTUAL native R3

The full native identity in A[[T]] is

    R(T)=(1+aT+eT^2)/[(1+alpha T+T^2)(1+gamma T)].

Here e^2 has already vanished by(2). Taking the first three
coefficients, or multiplying this finite recurrence, gives

    R3=C+a(1+alpha^2+alpha*gamma+gamma^2)
          +e(alpha+gamma) modulo2O,
    C=gamma^3+gamma+alpha*gamma^2+alpha^2*gamma+alpha^3.

The last two terms have value at least1/2. The exact field identity

    C=gamma*q0+alpha^2*gamma+alpha^3                     (7)

and(6) imply v(C)>=min(1/3,2x). Since this lower bound is at most1/3,
the error of value at least1/2 cannot spoil it. Therefore

    v(R3)>=min(1/3,2x)>=2/7,

because the necessary profile already has x>=1/7. This proves Theorem1
without either deeper finite native-bound computation.

## 4. A short finite-width contradiction on any four-odd grid

The necessary profile leaves starts6 or10 modulo16. At6, the actual
nonconstant coefficient a belongs to E and has value7/8, impossible
on a grid1/(4u) with u odd. It suffices to treat n=16H+10.
Put N=n+5=16H+15=8q+7, where q=2H+1 is odd.

The monic Appell polynomial P_N(Z)=sum_r(N)_r R_r Z^(N-r)
has k as a root of multiplicity at least6, by differentiating the six
actual zero moments. Let F be k's minimal polynomial over E. It exists
because P_N is monic. At delta7/8, two least exponents of F differ by
an even integer: their value difference times7/8 lies in1/(4u)Z.
Thus F has positive least-cost width at least2. Widths add by scaling
Z=kV and multiplying nonzero polynomials over the maximal-ideal
RESIDUE FIELD. Characteristic-zero separability gives F^6|P_N, so

    width(P_N)>=12.                                    (8)

For a polynomial exponent j=N-r=8b+c,0<=c<8, Legendre's formula
gives normalized term cost

    [f(b)-f(q)]+s_2(c)-c/8+v(R_(N-j)),
    f(b)=s_2(b)-b.

The unit native R5 at j=8q+2 bounds the minimum M by3/4. Any
earlier block b<=q-2 costs at least1; only b=q,q-1 can minimize.
Among low indices, c=1 costs7/8 and every c other than0,2,4
costs at least5/4. Thus least exponents are contained in the span
from j0=8(q-1) to j1=8q+4, of length12.

Equality in(8) forces both endpoints to be least. At the upper
endpoint the native index is3 and its normalized scalar cost is1/2.
Theorem1 therefore forces

    M=1/2+v(R3)>=1/2+2/7=11/14>3/4,

contrary to the unit R5 candidate. This proves Theorem2. No sixth-power
residue-support calculation is needed in this shorter route. The argument
uses actual native coefficients, not arbitrary recurrence states.

## 5. The stronger cusp when a later alpha bound is supplied

**Corollary.** If one also uses the already proved universal s10 bound
v(alpha)>=3/14, then Theorem1 gives

    v(R3)>=1/3.

Indeed2v(alpha)>=3/7>1/3. More generally the formula in Theorem1
has threshold x>=1/6, so the earlier9/49 bound also suffices for this
last sharpening. Neither later bound is used by Theorem2 above.

The optional input is the
[three-fourteenths native bound](BINARY_QUADRATIC_L8_S10_UNIVERSAL_THREE_FOURTEENTHS_BOUND.md).
No stronger universal alpha bound, arbitrary-field L8 closure, polynomial
residue conclusion or unrestricted quadratic SFC is asserted here.
