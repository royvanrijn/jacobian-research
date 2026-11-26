# Higher critical cube residues on a quarter-integer coefficient grid

This written proof excludes six consecutive native binary quadratic
moment zeros at every critical radius v(k)=1-1/L with L=2^a>=8, in the
specified quarter-integer coefficient grid and cubic residue profile.
It does not assert unrestricted binary quadratic SFC. The separate L=4
proof is not used in the argument below.

Let E be a characteristic-zero real-valued nonarchimedean field, v(2)=1,
with v(E*) contained in (1/4)Z. Let

 q=alpha X^2/2+beta XY+gamma Y^2/2+X+eY,
 R_r=Lambda(q^r)/r!, R_0=1,
 P_m(Z)=sum_(r=0)^m (m)_r R_r Z^(m-r).

All natural coefficients are integral. Assume their residues satisfy

 alphabar=ebar=0, betabar=gammabar=1.                     (1)

Thus the native residue generating series is

 Rbar(T)=1/(1+T)^3,
 Rbar_(4j)=Rbar_(4j+1)=1,
 Rbar_(4j+2)=Rbar_(4j+3)=0.                              (2)

Every nonzero positive valuation of an E-element is at least 1/4.
The constant k may belong to any real-valued extension, with
v(k)=1-1/L, L=2^a>=8. All equalities of moments refer to the original
factorial functional, not arbitrary recurrence states.

**Theorem.** Under these hypotheses, P_n(k),...,P_(n+5)(k) cannot
all be zero for any positive n.

## 1. Width facts and a necessary starting congruence

For a nonzero finite polynomial A(Z)=sum a_i Z^i and slope delta=v(k),
let its least-cost indices minimize v(a_i)+i delta. Its width w(A) is
the largest such index minus the smallest. After replacing Z by kU
and dividing by one least coefficient times its power of k, the
nonzero reduction is a polynomial over the residue field. Its first
and last nonzero exponents are precisely those least-cost indices.
Products of these reductions show

                         w(AB)=w(A)+w(B).                (3)

Indeed a product of nonzero residue polynomials has nonzero leading
and trailing coefficients. This proof uses finite polynomials only.

If A belongs to E[Z], any two least-cost indices i,j satisfy
(i-j)delta in (1/4)Z. At delta=(L-1)/L this implies

              L/4 divides every difference i-j.         (4)

If A(k)=0 then at least two terms attain the minimum, so w(A)>0.
For the irreducible minimal polynomial F of k over E, it follows that
w(F)>=L/4. Existence of any P_m(k)=0 makes k algebraic, because P_m
is monic (R_0=1). Characteristic zero makes F separable. If k is a root
of multiplicity h of P_m, then F^h divides P_m in E[Z], and hence

                       w(P_m)>=h L/4.                   (5)

The exact critical Appell coefficient calculation from
[the discrete-grid multiplicity proof](QUADRATIC_CONSTANT_APPELL_MULTIPLICITIES_FROM_DISCRETE_VALUE_GRIDS.md),
Section2, also gives: if R_(m mod L) is a unit, a root of P_m at this
slope has multiplicity at most gcd(L,4)=4. Its proof follows by the
two base-cost anchors, which give w(P_m)=L when this selected R is a
unit, and then (3)--(4). This statement needs neither completeness
nor a root-counting theorem.

Suppose a six-zero block exists. If n is not 2 modulo4, at least one
of its last two indices m=n+4,n+5 is 0 or1 modulo4. By (2),
R_(m modL) is a unit. The Appell derivative identities
P_m^(j)=(m)_j P_(m-j), j<=m, make k a root of multiplicity at least
five of this P_m. This contradicts the preceding bound4. Therefore

                 n=2 mod4, N=n+5=3 mod4.                 (6)

The six-zero block makes k a root of multiplicity at least six of
P_N. By (5), w(P_N)>=3L/2.

## 2. All possible least coefficients of the terminal Appell polynomial

Write N=Lq+t, 0<=t<L; thus t=3 mod4 and t>=3. Set j=N-r=Lb+c,
0<=c<L. Legendre's formula gives the coefficient cost

 v((N)_r k^(N-r))=N-s2(N)+s2(j)-j/L.

After subtracting the base minimum, the actual term cost is

 [f(b)-f(q)]+g(c)+v(R_(N-j)),
 f(b)=s2(b)-b, g(c)=s2(c)-c/L.                            (7)

All low-bit contributions to g are positive: bit 2^i contributes
1-2^i/L. Also f(b-1)-f(b)=v2(b). If q is even, every b<q has
f(b)-f(q)>=1. If q is odd, exactly b=q,q-1 can have difference zero;
all earlier b have difference at least one. When q=0 or1 the same
assertions apply with absent negative b omitted.

The actual coefficient at j=Lq+2 is allowed because t>=3. Its R-index
is t-2=1 mod4, so it is a unit by (2). Thus the least actual cost M,
relative to the base minimum, satisfies

                         M<=1-2/L<1.                    (8)

Consequently all blocks b other than those just identified are excluded.
Within a permitted block the following low-bit cases are exhaustive.

* If c=0 or L/2, it remains a candidate. In either case its R-index
  is 3 mod4 and its R-coefficient has valuation at least1/4.
* If c is another multiple of4, it contains a bit at most L/4, so
  g(c)>=3/4. Its R-index is3 mod4; its total cost is at least1.
* If c=1 mod4, g(c)>=1-1/L>1-2/L, excluding it by (8).
* If c=2 mod4, g(c)>=1-2/L, with equality only at c=2. Every additional
  bit has strictly positive cost. Thus only c=2 can remain.
* If c=3 mod4, its two lowest bits already give
  g(c)>=2-3/L>1, excluding it.

It follows that every least-cost j belongs to

 {Lq,Lq+2,Lq+L/2}
 union {L(q-1),L(q-1)+2,L(q-1)+L/2},                     (9)

with the second set present only when q is odd, and all candidates
still subject to j<=N. The total possible span is at most3L/2.
Since w(P_N)>=3L/2, equality holds. In particular q is odd, t>=L/2,
and BOTH extreme j=L(q-1),Lq+L/2 must be least-cost indices.

At these extremes the costs are respectively v(R_(t+L)) and
1/2+v(R_(t-L/2)). Both R-indices are3 mod4, so the latter cost is
at least3/4. The former is an integral multiple of1/4. Together with
(8) these force

 M=3/4, v(R_(t+L))=3/4, v(R_(t-L/2))=1/4,
                         w(P_N)=3L/2.                   (10)

## 3. A first-order native coefficient identity excludes L>=16

Put z=gamma-1 and b=beta-1. The four perturbations alpha,b,z,e all
have valuation at least1/4. The integral native polynomial coefficient
identities, reduced modulo2 over universal indeterminates, give

 R(T)=[1+(1+e+beta)T+eT^2]
       /[(1+alpha T+T^2)(1+gamma T+e^2T^2)] mod2.

At alpha=e=b=z=0 this is (1+T)^(-3). Its four first derivatives over
characteristic two are

 partial_alpha R=T/(1+T)^5, partial_e R=T/(1+T)^2,
 partial_b R=T/(1+T)^3, partial_z R=T/(1+T)^4.             (11)

For any r=3 mod4 the last two derivatives have zero coefficient at
T^r, whereas partial_e has coefficient1. The coefficient in
partial_alpha is1 if r=3 mod8 and0 if r=7 mod8: use
(1+T)^(-5)=(1+T)^3/(1+T^8).

The constant native coefficient at such an r is even. Taylor terms
of degree at least two in the four perturbations have valuation at
least1/2. Replacing the integer linear coefficients by their parities
only changes them by even multiples of perturbations. It follows,
coefficientwise and uniformly in all r=3 mod4, that

 R_r=e+alpha mod v>=1/2 if r=3 mod8,
 R_r=e       mod v>=1/2 if r=7 mod8.                     (12)

This uses universal integral coefficients before specialization; no
identity over the residue field is silently promoted to a stronger
congruence. If L>=16, the two indices t+L and t-L/2 in (10) differ by
3L/2, which is divisible by8. Equation (12) makes their difference
have valuation at least1/2. Their different valuations3/4 and1/4 in
(10) instead make that difference have valuation exactly1/4. This
contradiction excludes every L>=16.

## 4. The sixth-power residual polynomial excludes L=8

When L=8, t=3 mod4 and t>=L/2 imply t=7. Equation (10) gives
width(P_N)=12, and (5) forces the root multiplicity h to be exactly6
and width(F)=2. The quotient P_N/F^6 has width zero.

Use the finite residue-polynomial construction from Section1, replacing
Z by kU. By (4), all least exponents of F differ by an even integer.
Since width(F)=2, after removing its lowest monomial its residue is
exactly A+B U^2 with A,B nonzero. Width zero makes the residue of the
quotient a nonzero monomial. Therefore, after removing its own lowest
monomial, the residue of P_N is a nonzero scalar times

 (A+B U^2)^6=A^6+A^4 B^2 U^4+A^2 B^4 U^8+B^6 U^12.      (13)

Here the equality is in characteristic two. In particular its
coefficient at exponent2 must vanish.

But the least-cost exponents of P_N are j=8(q-1) and j=8q+4.
The allowed exponent j=8(q-1)+2 corresponds to R_(t+L-2)=R_13,
which is a unit. Its cost is1-2/L=3/4, the actual minimum by (10).
Its residue is nonzero, and relative to the lowest exponent j=8(q-1)
it occurs at exponent2, contradicting (13). Equivalently the upper
c=2 coefficient corresponds to R_5 and gives the also-forbidden
offset10. This excludes L=8 and proves the
theorem. No first-order Taylor identity is needed for this last case.

## Scope of the possible assembly

The separate quarter-gap L=4 proof can combine with this theorem.
For a naturally integral coefficient field on grid(1/4)Z the necessary
quarter-gap is automatic after the cube-profile normalization.
The noncube single-pole unit-window classification and the two-pole
theorem deal with the remaining nonpolynomial residue types; open
valuation bands and elementary constant ranges are already covered
by the constant-window transfer. Such an assembly would yield an
all-constant SIX-window for nonpolynomial binary residue on the
quarter-integer coefficient grid. That assembly requires its own
checked statement and proof-source links.

No claim is made for a general grid1/(4u) with odd u>1: the positive
valuation step can be smaller than1/4 and the cost/Taylor arguments
above then require changes. Polynomial residues, unrestricted valued
coefficients, smaller sparse support windows and full SFC(2)/SFC(3)
remain separate. This proof does not use an analytic EGF, completeness,
generic rank, numerical experiments or a worldwide-priority assertion.

