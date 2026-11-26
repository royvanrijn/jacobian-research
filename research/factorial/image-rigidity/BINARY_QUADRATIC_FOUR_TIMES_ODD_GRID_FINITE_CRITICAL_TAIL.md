# A finite critical-radius frontier for grids with 2-part four

This proof excludes all sufficiently high critical radii in the cubic-residue case
when the coefficient value group is contained in (1/(4u))Z, u odd. It does
not close the remaining finitely many radii or arbitrary quadratic SFC.

Let E have characteristic zero, a real-valued nonarchimedean valuation
v(2)=1, and v(E*) contained in (1/(4u))Z for an odd positive integer u.
Consider the naturally integral binary quadratic

 q=alpha X^2/2+beta XY+gamma Y^2/2+X+eY,
 alphabar=ebar=0, betabar=gammabar=1,
 R_r=Lambda(q^r)/r!, P_m(Z)=sum_r(m)_r R_r Z^(m-r).

Its native residue is Rbar(T)=(1+T)^(-3), with unit coefficients at
indices 0,1 modulo4 and zero residues at indices2,3 modulo4. Set

 J=2u, B=the least power of two with B>=3J=6u.

**Theorem.** If L is a power of two with L>=2B, no constant
k in any real-valued extension with v(k)=1-1/L can produce six
consecutive zero factorial moments of q+k at a positive start.

## 1. Uniform finite-degree periodicity below valuation one half

The coefficients R_r are universal integer polynomials in the natural
coefficients. Reduce those universal polynomials modulo2, and put

 b=beta-1, z=gamma-1, I=(alpha,b,z,e).

The universal native rational identity is R(T)=N(T)/D(T), with

 N=1+(b+e)T+eT^2,
 D=(1+alpha T+T^2)(1+(1+z)T+e^2T^2).

At I=0, N0=1 and D0=(1+T)^3. Write D=D0+deltaD. The polynomial deltaD
has I-order at least1 and T-degree at most4. Modulo I^J the finite
geometric expansion gives

 R(T)=(N/D0) sum_(a=0)^(J-1)(deltaD/D0)^a.               (1)

The denominator of each summand is (1+T)^(3(a+1)). Since B>=3J,
put every summand over (1+T)^B=1+T^B in characteristic two. Its
numerator has degree at most

 2+4a+B-3(a+1)=B+a-1<=B+J-2.

Consequently, in (F2[alpha,b,z,e]/I^J)[[T]],

                  R(T)=A(T)/(1+T^B), deg A<=B+J-2.      (2)

Coefficient comparison implies R_(r+B)=R_r for all r>=J-1 in
this quotient. On specializing to E, every I-monomial of degree
at least J has valuation at least J/(4u)=1/2, and every discarded
even coefficient has valuation at least1. Therefore the ACTUAL
native coefficients obey

        v(R_(r+B)-R_r)>=1/2 for every r>=J-1.            (3)

In particular two indices at least J-1 which agree modulo B have
coefficients congruent modulo valuation at least1/2. This argument
uses finite deformation degree and coefficientwise formal series;
it never evaluates an infinite series at a unit T.

## 2. Six zeros saturate the available least-cost span

Suppose P_n(k)=...=P_(n+5)(k)=0. The critical unit-anchor multiplicity
bound is gcd(L,4u)=4: a unit R_(m modL) makes the least-cost width
of P_m equal L, while any irreducible root factor has positive width
divisible by L/4. Thus n must be2 modulo4; otherwise one of its last
two moments has a unit anchor and root multiplicity at least5.
Put N=n+5, so N=3 modulo4.

For completeness, least-cost width at delta=v(k) is the difference
between the largest and smallest exponents minimizing v(a_i)+i delta
in a finite polynomial sum a_i Z^i. Widths add under products: scale
Z=kU, normalize by a least term, and multiply nonzero residue
polynomials. If F is the minimal polynomial of k over E, the equation
F(k)=0 makes its width positive. Any two least exponents differ by a
multiple of L/4 because

 (i-j)(L-1)/L in (1/(4u))Z, u odd.

Characteristic-zero separability and root multiplicity at least6 in
P_N give F^6|P_N and hence width(P_N)>=6L/4=3L/2. No completeness or
Newton root count is needed.

Write N=Lq+t, 0<=t<L. Then t=3 modulo4 and t>=3. For a coefficient
indexed by r, use its POLYNOMIAL EXPONENT j=N-r=La+c, 0<=c<L.
Legendre's formula, after subtracting the base minimum, gives cost

 [f(a)-f(q)]+g(c)+v(R_(N-j)),
 f(a)=s2(a)-a, g(c)=s2(c)-c/L.                            (4)

The allowed unit coefficient at j=Lq+2, whose R-index is t-2=1
modulo4, bounds the actual minimum M by

                         M<=1-2/L<1.                    (5)

The identity f(a-1)-f(a)=v2(a) shows that only a=q, and also a=q-1
when q is odd, can have base cost less than one. Every low-bit
weight in g is positive, with the two smallest weights1/2 and3/4.
Thus any c with at least two set bits costs at least5/4 and cannot
minimize. The remaining candidates are c=0 or a single power of two
at most L/2 (c=1 is also excluded by (5)). Their total possible span
in j is at most3L/2.

The width lower bound forces equality. Thus q is odd, t>=L/2,
both extremes j=L(q-1), Lq+L/2 are least indices, the root
multiplicity is exactly6, width(F)=L/4, and width(P_N/F^6)=0.

All least exponents of F differ by multiples of L/4. Its initial
residue polynomial, after removing its lowest monomial, is therefore
A+B0 U^(L/4) with A,B0 nonzero. The initial polynomial of P_N is a
monomial times its sixth power. In characteristic two this has FOUR
nonzero offsets, namely

                         0,L/2,L,3L/2.                  (6)

In particular the intermediate exponent j=L(q-1)+L/2 is also a
least index. Its R-index is t+L/2, whereas the smallest exponent
j=L(q-1) has R-index t+L. Formula(4) gives

        v(R_(t+L))=M, v(R_(t+L/2))=M-1/2.                (7)

Both indices are3 modulo4 and hence have positive coefficient
valuation. Thus M>1/2, while (5) makes M<1. The two valuations in
(7) are respectively greater than1/2 and strictly less than1/2.

## 3. Periodicity contradicts the saturated middle coefficient

The two indices in (7) are at least L, so are at least J-1 because
L>=2B>=6J. Their difference is L/2. Since L and B are powers of two
and L>=2B, that difference is a multiple of B. Formula(3) therefore
makes their coefficient difference have valuation at least1/2.
Their distinct valuations in (7) instead make the difference have
valuation M-1/2<1/2, a contradiction. This proves the critical-tail exclusion.

For u=1 the bound is L>=16. For u=3 (a grid1/12), B=32 and the bound
is L>=64; this leaves the critical radii L=4,8,16,32 to other methods.
The one-pole noncube window<=2 plus the old grid transfer already
excludes six zeros for those other one-pole residue types. Multiple
poles are covered by the old two-pole theorem. Thus this result can
reduce a nonpolynomial-residue grid with 2-part four to finitely many
critical radii, without claiming that these remaining radii are closed.

The native universal residue and integrality are proved in
[the noncritical transfer proof](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md).
The exact critical unit-anchor multiplicity and finite-width method are
proved in [the critical detector proof](BINARY_QUADRATIC_CRITICAL_DYADIC_SIX_MOMENT_WINDOWS.md)
and [the grid-width transfer proof](QUADRATIC_CONSTANT_APPELL_MULTIPLICITIES_FROM_DISCRETE_VALUE_GRIDS.md). The new step is finite deformation-degree periodicity combined
with sixth-power saturation of all four anchors. No computational
enumeration, polynomial-residue exclusion or arbitrary-ramification
all-constant theorem is asserted.
