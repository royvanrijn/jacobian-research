# Native periods and all-constant windows on the twelfth coefficient grid

Let E be a characteristic-zero valued field with v(2)=1. Natural binary
quadratic coefficients mean

    q=alpha X^2/2+beta XY+gamma Y^2/2+dX+eY,
    R_r=Lambda(q^r)/r!,  Lambda(X^iY^j)=i!j!, R_0=1.

The natural coefficients are integral. The native series R(T) is formal;
its maximal-ideal residue is denoted Rbar. No series is evaluated at a unit.
For the local arguments, E and the extensions in use are real-valued.
The final theorem explains why an arbitrary valued constant extension causes
no additional case when the coefficient value group is discrete.

The following written inputs are used with their exact scopes. The package
manifest pins every proof artifact; the mathematical links are relative to this
proof directory.

* [Native integer coefficients and the universal congruence over O/2O](QUADRATIC_CONSTANT_NONCRITICAL_DYADIC_WINDOW_TRANSFER.md).
* [The native critical detector](BINARY_QUADRATIC_CRITICAL_DYADIC_SIX_MOMENT_WINDOWS.md), used inside the imported L8 profile proof.
* [Finite least-cost widths, critical unit-anchor multiplicity, and abstract unit-window transfer](QUADRATIC_CONSTANT_APPELL_MULTIPLICITIES_FROM_DISCRETE_VALUE_GRIDS.md).
* [The GENERAL necessary critical single-pole frontier](BINARY_QUADRATIC_TWO_POLE_CONSTANT_SHIFTS_STRONG_FACTORIAL_WINDOWS.md), Section2, Corollary3. This is not merely its conditional grid corollary.
* [Native single-pole classification](BINARY_QUADRATIC_SINGLE_POLE_RESIDUE_UNIT_WINDOWS.md).
* [The gap-free quarter-critical six-window](BINARY_QUADRATIC_QUARTER_CRITICAL_NO_GAP_SIX_WINDOW.md).
* [The necessary L8 reduction](BINARY_QUADRATIC_L8_EXACT_START_AND_NATIVE_PROFILES.md).
* [The universal s10 native bound](BINARY_QUADRATIC_L8_NATIVE_BOUND_AND_TWELFTH_GRID_SIX_WINDOW.md), specifically Theorem1 rather than its grid12-only consequence.

The first new conclusion below strengthens the grid12-only L8 application.
The necessary L8 profiles and the universal native bound are not changed.

**Theorem 1 (L8 on grids with dyadic part four).** Suppose v(E*) is contained
in (1/(4u))Z, where u is an odd positive integer, and the native residue is
the normalized cube described in Section1. At v(k)=7/8, q+k has no six
consecutive zero factorial moments at any positive start. The constant may
belong to a real-valued extension of E.

**Theorem 2 (higher cube radii on the twelfth grid).** Under the same cube
normalization and v(E*) contained in (1/12)Z, no critical radius
v(k)=1-1/L with a power of two L>=16 admits such a six-zero block.

**Theorem 3 (all constants, nonpolynomial residue, twelfth grid).** For every
naturally integral binary q over a valued E with v(E*) contained in
(1/12)Z and nonpolynomial Rbar, q+k has no six consecutive zero factorial
moments at any positive start, for every k in ANY valued extension of E.
No rank or Hessian condition is imposed. Polynomial native residues are
excluded from this statement.

## 1. Two uniform native periods

For this section assume the normalized cubic residue data

    d=1, alphabar=ebar=0, betabar=gammabar=1.

Put z=gamma-1 and a=1+e+beta. Thus alpha,z,e,a have positive
valuation or are zero. In the possibly nonreduced ring O/2O, the universal
integer native identity is

    R(T)=N(T)/D(T),
    N=1+aT+eT^2,
    D=(1+alpha T+T^2)(1+(1+z)T+e^2T^2).                 (1)

On grid1/12, every nonzero positive coefficient has value at least1/12.
Introduce a formal symbol w for e^2. Expanding the two axis factors around
(1+T)^2 and1+T gives, coefficientwise in characteristic two,

    R=N sum_(i,j,k>=0) binom(j+k,j)
      alpha^i z^j w^k T^(i+j+2k)
      /(1+T)^(3+2i+j+k).                               (2)

For a fixed coefficient of T only finitely many triples contribute.
Discard the triples of weight i+j+2k>=6. After specialization, each
discarded coefficient has value at least1/2; all even integer errors
already have value at least1. The remaining denominator exponents are
at most13. Putting the finite sum over (1+T)^16=1+T^16 gives a numerator
of degree at most

    16-3-i+k+deg N <=17,

since k<=2. Hence there is a polynomial A16 of degree at most17 such that
R=A16/(1+T^16) coefficientwise modulo value at least1/2. Consequently

    v(R_(r+16)-R_r)>=1/2 for EVERY r>=2.                 (3)

This is an eventual coefficient identity, not a conclusion obtained by
evaluating a rational expression at a unit.

There is a sharper version when v(alpha)>=1/6. Use weight2i+j+2k
in(2). Again discard weight>=6. Now
3+2i+j+k<=8 for every retained triple. A common denominator1+T^8
has numerator degree at most

    8-3-i+k+deg N <=9.

Thus in this case

    v(R_(r+8)-R_r)>=1/2 for EVERY r>=2.                  (4)

These statements include alpha=0 or e=0. They use only universal native
integer coefficients, the actual coefficient-grid lower bounds, and finite
deformation degree. The numerator coefficients need only be integral;
counting their positive valuations could improve the degree bounds but is
unnecessary.

## 2. A local saturation lemma for every L>=8

Here take any grid1/(4u) with u odd, the normalized cube, and a power
of two L>=8. Suppose a six-zero block starts at n, and v(k)=1-1/L.
Write

    P_m(Z)=sum_r (m)_r R_r Z^(m-r)=Lambda((q+Z)^m).

These are monic polynomials in Z. The imported unit-anchor bound is4:
if R_(m mod L) is a unit, a root of P_m has multiplicity at most
gcd(L,4u)=4. The cube has unit native coefficients at indices0,1
modulo4. If n is0 or1 modulo4, m=n+4 is such an anchor with root
multiplicity at least5. If n is3 modulo4, use m=n+5. Therefore

    n=2 modulo4, N=n+5=3 modulo4.                        (5)

For completeness the remaining finite-width argument is local and works
at all L>=8; it does NOT assume the larger L bound of a previous finite-tail
theorem. At a fixed delta, the width of a polynomial is the span of its
least exponents for v(c_j)+j delta. Widths add under products: scale Z=kV,
normalize by a least term and multiply nonzero polynomials over the
MAXIMAL-IDEAL RESIDUE FIELD. This is distinct from the nonreduced O/2O
used for native congruences.

If F is the minimal polynomial of k over E, it exists because P_N is
monic and P_N(k)=0. Its width is positive. Two least exponents i,j of F
obey (i-j)(L-1)/L in (1/(4u))Z, so their difference is a multiple of
L/4. Characteristic-zero separability and root multiplicity at least6
give F^6 dividing P_N. Thus

    width(P_N)>=6L/4=3L/2.                              (6)

Write N=Lq+t, 0<=t<L. In(5), t=3mod4 and t>=3. For the polynomial
exponent j=N-r=Lb+c, 0<=c<L, subtracting the common base v(N!)+f(q)
from its term value gives

    [f(b)-f(q)]+g(c)+v(R_(N-j)),
    f(b)=s_2(b)-b, g(c)=s_2(c)-c/L.                      (7)

The unit native R_(t-2) at j=Lq+2 gives M<=1-2/L<1 for the minimum
of these normalized costs. Since f(b-1)-f(b)=v_2(b), only b=q and,
when q is odd, b=q-1 can minimize. Earlier blocks cost at least1.
The two smallest distinct low-bit weights in g are1/2 and3/4;
two-bit indices therefore cost at least5/4. The index c=1 costs more
than the displayed unit candidate. All remaining candidates are c=0
or a single power of two from2 through L/2.

Their possible total span is at most3L/2. Equality in(6) forces q odd,
t>=L/2, and least endpoints j0=L(q-1), j1=Lq+L/2. It also forces
the root multiplicity exactly6, width(F)=L/4 and width(P_N/F^6)=0.
After removing a monomial, the initial polynomial of F is a binomial
with gap L/4. Its sixth power in characteristic two has four nonzero
offsets

    0, L/2, L, 3L/2.

In particular j0+L/2 is also least. Formula(7) gives

    v(R_(t+L))=M,
    v(R_(t+L/2))=M-1/2,
    1/2<M<=1-2/L<1.                                    (8)

Both native indices are3 modulo4, hence nonunits, which gives the strict
lower inequality. These indices are at least L. The upper endpoint
also gives v(R_(t-L/2))=M-1/2. No completeness, conjugate-root counting
or generic-rank assertion was used.

On grid1/12, (3) and(8) already exclude every L>=32: the two large
native indices differ by L/2, a multiple of16, while their difference
has value M-1/2<1/2. This proves the high-radius part of Theorem2.

## 3. The broader L8 consequence of the universal native bound

The necessary L8 reduction excludes starts7,8,9 modulo16 and leaves6
or10. At6 it forces the actual coefficient a=1+e+beta in E to have
value7/8. This is impossible on grid1/(4u), u odd: a denominator8
cannot occur in that grid. For s10, its strong O/2O profile and the
universal native bound give, with x=v(alpha), y=v(z),

    v(a),v(e)>=1/2, x>=9/49>1/6,
    y=x/2 if x<1/4; y>=1/8 if x>=1/4,
    z^8+alpha^4+alpha^6+z^8 alpha^4+z^8 alpha^6 in2O.     (9)

The value infinity is allowed. Every term in the last line except
z^8+alpha^4 has value greater than1: if x<1/4, use y=x/2 and
6x>1; if x>=1/4, use y>=1/8. Hence z^8+alpha^4 belongs to2O.
The difference between (z^2+alpha)^4 and z^8+alpha^4 has even integer
coefficients, so

    v(z^2+alpha)>=1/4.                                  (10)

Here no grid rounding was used. With e^2=0 under(9), the exact native
formula modulo2O gives

    R3=gamma^3+gamma+alpha*gamma^2+alpha^2*gamma+alpha^3
       +a(1+alpha^2+alpha*gamma+gamma^2)+e(alpha+gamma).

The last two terms have value at least1/2. Expanding gamma=1+z,
the first line, modulo2O, is

    (z^2+alpha)+z^3+alpha*z^2+alpha^2*(1+z)+alpha^3.

Every term after the first has value greater than1/4 by(9). For example,
when x<1/4, 3y=3x/2>1/4 because x>1/6; the other case is immediate.
Thus v(R3)>=1/4 universally in s10.

Apply Section2 at L=8. Its least upper endpoint has native index3,
so M>=1/2+1/4=3/4. The unit c=2 candidate gives M<=3/4.
Saturation therefore has M=3/4 and initial support offsets0,4,8,12.
But j=8(q-1)+2 has native index13 and is a unit coefficient of cost
3/4=M, at forbidden offset2. This contradiction excludes s10 on
every grid1/(4u), u odd, and proves Theorem1. In particular it supplies
the L8 case needed for the twelfth grid.

## 4. At L16 the actual native equations give polynomial divisibility

From now on use grid1/12 and delta=15/16. Section2 forces q odd and
t>=8. Since t=3mod4, t is11 or15. Therefore the only possible starts are

    n=22 or26 modulo32.                                (11)

Put U=k^16/2^15. Through degree31, the finite exponential is, in O/2O,

    exp(kT)=B(T)(1+UT^16),
    B=1+kT+c2*T^2+c4*T^4+c8*T^8,
    cj=k^j/j!, v(c2)=7/8,v(c4)=3/4,v(c8)=1/2.

Indeed the costs s_2(j)-j/16 below32 leave precisely these degrees
and their translates by16. Products of two nonconstant coefficients
of B vanish, and B^2=1.

Write n=32H+s, s22 or26. For each of its six zero moments, divide by
the nonzero coefficient of the native term indexed by i=s,...,s+5.
The head coefficients are a_j/binom(32H+j,j), j<=31. The binomials
are odd, all a_j are integral, so the head is the actual small E_i
modulo2O. A further tail has exponent j=32b+c, b<=H-1, 0<=c<32;
its relative cost is

    [s_2(b)-2b]-[s_2(H)-2H]+s_2(c)-c/16>=1.

The last low-bit bracket is nonnegative, and the first bracket increases
by at least1 on going backwards. Thus the actual necessary small
equations E_s=...=E_(s+5)=0 hold in O/2O for every H, not just a sampled
prefix. If H=0 the tail is absent.

Let P be the actual prefix sum_(i<s) E_i T^i. The native identity(1)
then gives

    D P=B N(1+UT^16) modulo T^(s+6).                     (12)

For s22, the two sides have degrees at most25 and26, both below28.
For s26, they have degrees at most29 and26, both below32. Thus(12)
is an EXACT polynomial identity in both cases, with the full native
quartic D and without dividing by a possibly nonunit coefficient.

The polynomial N is a unit in (O/2O)[T]. In fact a,e have positive
values at least1/12, so N^16=1 in this ring. Together with B^2=1 this
makes BN a polynomial unit, with inverse B*N^15. Multiplication in(12)
therefore proves

    D divides 1+UT^16 in (O/2O)[T].                      (13)

This uses a POLYNOMIAL inverse, not merely a formal-series inverse.
Nilpotent leading coefficients of D present no problem.

## 5. The L16 remainder forces a refined period

Work in a valued algebraic extension if necessary. Factor the Y quadratic

    1+gamma T+e^2T^2=(1+bT)(1+cT),
    b+c=gamma, bc=e^2, v(b)=0, v(c)=2v(e)>=1/6.

If e=0, take b=gamma,c=0. These values follow from the unit sum and
the positive-value product; no completeness is needed. The small factor
1+cT is a polynomial unit modulo2O. Equation(13) implies that
Ax=1+alpha T+T^2 divides1+UT^16. Evaluating at the actual unit root
T=b^(-1) also gives U=b^16 modulo2O.

Monic remainder division by Ax gives

    T^16=1+alpha^8+alpha^12+alpha^14+alpha^15 T,
    alpha^15=0,
    1+U(1+alpha^8+alpha^12+alpha^14)=0 modulo2O.          (14)

Suppose v(alpha)=1/12, the only positive grid12 value below1/6.
Let y=b-1, which has positive value. In characteristic two modulo2O,
b^16=1+y^16. In(14), alpha^8 has value2/3 while alpha^12 and all
higher displayed terms vanish modulo2O. The remaining relation forces

    16v(y)=8v(alpha)=2/3, hence v(b-1)=1/24.

But v(c)>=1/6, so the ACTUAL coefficient gamma-1=(b-1)+c has the
same value1/24. It belongs to E, contradicting the grid1/12 hypothesis.
No claim that b or c belongs to E is made. Therefore

    v(alpha)>=1/6.                                     (15)

Section1 now gives the native period8 modulo value at least1/2 from
every index r>=2. In the saturated L16 configuration, Section2 gives
native indices t+16 and t+8 with values M and M-1/2 respectively,
where1/2<M<1. Their difference is8 and both indices are large. Period8
requires their coefficient difference to have value at least1/2, whereas
its actual value is M-1/2<1/2. This contradiction excludes L16.
Combined with the period16 exclusion of L>=32, it proves Theorem2.

## 6. All constants and the actual-support consequence

First take a real-valued extension constant and suppose a six-zero block
exists with nonpolynomial native residue on grid1/12. The imported
GENERAL frontier forces a critical radius1-1/L, L>=4, and a single
geometric pole canceled by its critical detector. The native single-pole
classification gives a positive unit window w<=2 except for the exact
cube Rbar=(1+rho T)^(-3). In a noncube case the imported abstract
grid transfer gives a window w+4<=6 for every constant, a contradiction.
In the cube case, its native profile allows a variable swap and division
of the whole phase by an E-unit linear coefficient, preserving grid1/12
and the constant value, to obtain the normalized data used above.

The power-of-two radii are exhausted: L4 is excluded by the gap-free
quarter theorem, L8 by Theorem1, and all L>=16 by Theorem2. This proves
the desired six-window over real-valued extensions.
For the L8 step on this particular grid, Theorem2 of the unchanged
L8 native-bound proof is an alternative to Theorem1 above. That route
does not use the new R3 cusp bridge; both routes use the universal9/49
bound rather than the later stronger three-fourteenths refinement.

Now allow an ANY valued extension for the constant. If k is
transcendental over E, no monic P_m can vanish at k. Otherwise every
nonzero element x of E(k) is algebraic over E. A least-term tie in its
minimal polynomial gives

    (i-j)v(x)=v(a_j)-v(a_i) in (1/12)Z for some i!=j.

An ordered value group is torsion-free. Thus all these values lie in the
rational divisible hull of the base grid, with its usual order. Restricting
to E(k), and making the algebraic adjunctions used above, gives a
real-valued setting to which the proof applies. This proves Theorem3
without an unmentioned extension-rank assumption.

**Corollary (dense actual support).** Under Theorem3, if all five
nonconstant binary quadratic coefficients are nonzero, q+k satisfies its
actual-support strong factorial window for every constant. For k!=0
there are six monomials and Theorem3 applies directly. For k=0 there
are five monomials; the native nonpolynomial unit window is at most4,
so a five-zero block is impossible as well.

**Corollary (necessary remaining frontier).** Any actual six-zero block
with integral natural coefficients and nonpolynomial native residue must
have a NONCONSTANT coefficient value group not contained in (1/12)Z.
If it nevertheless admits a bound (1/(4u))Z with u odd, its constant
radius must be critical1-1/L with a power of two L>=16, and its single
geometric pole must be completely canceled by the critical detector.
In this conditional grid sector, delta>=15/16.

The first assertion is the contrapositive of Theorem3. For the second,
the GENERAL frontier gives the critical single-pole condition; noncube
residues are excluded by their window at most2 and the abstract grid
transfer. The quarter theorem excludes L4 and Theorem1 excludes L8.
There is no arbitrary-coefficient-field delta>=15/16 assertion: grids
whose dyadic part is eight have not been covered by Theorem1. At the specific L8 radius7/8, any NONPOLYNOMIAL native-residue
six-zero block admitting a specified discrete coefficient grid(1/e0)Z
necessarily has8 dividing e0. Indeed, if8 does not divide e0, write
e0=2^s*u with u odd and s<=2. The coefficient grid is then contained
in(1/(4u))Z. The GENERAL frontier and single-pole classification reduce
to the cube case, because the noncube window at most2 and abstract
grid transfer give a window at most6. Theorem1 excludes the cube case.
Thus the hypothesized block cannot occur on such a grid. For these
discrete-grid statements the algebraic-value argument above works with
12 replaced by the specified denominator. No claim is made that a
multiple-of-eight grid actually admits a counterexample.

These are necessary restrictions and do not construct any remaining zero block.

The six-window is not a smaller sparse-support theorem. Polynomial
native residues, coefficient grids beyond the stated classes, and
unrestricted binary quadratic SFC, SFC(2) and SFC(3) remain separate.
No numerical range, generic-rank assumption, canonical update or
dependency-engine evaluation is a premise of this proof.
