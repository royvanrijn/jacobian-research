# Single-face Strong Factorial sectors of an even binary quadratic

Let K be a field of characteristic zero, let
Lambda_2(x^r y^s)=r!s!, and put
f=a x^2+bxy+c y^2+k and M_j=Lambda_2(f^j).
For t in K put B_(t,j)=Lambda_1((t z^2+k)^j).

**Theorem 1 (three coefficient planes).** Assume abck!=0 and at least
one of b=2a, b=2c, or a=c. For every n>=1, at least one of
M_n,M_(n+1),M_(n+2),M_(n+3) is nonzero. These are actual
four-monomial polynomials. The result includes the singular quadratic
matrices on these planes; no discriminant or coefficient-difference
denominator is excluded.

The only imported nonvanishing assertion is the already proved
constant-plus-square two-moment window: for t k!=0 over every
characteristic-zero field, B_(t,r) and B_(t,r+1) cannot both vanish
for any r>=1. Its all-unit divided-power/squarefreeness proof is
Section 2 of [Univariate quadratics with all supports](UNIVARIATE_QUADRATIC_ALL_SUPPORTS_STRONG_FACTORIAL.md),
digest 9254fb9345762aa5a67e36dd3c01a87e5f44821edcde28bd7f197d68e1c0cfa4.
No full binary quadratic or full SFC(2) theorem is imported or proved.

## 1. Polynomial Stein and Euler identities

For every polynomial P, direct evaluation on monomials proves
Lambda_2(partial_x P)=Lambda_2(P)-Lambda_1(P(0,y)),
Lambda_2(partial_y P)=Lambda_2(P)-Lambda_1(P(x,0)),                 (1)
and
Lambda_2(partial_x(xP))=Lambda_2(xP),
Lambda_2(partial_y(yP))=Lambda_2(yP).                             (2)

These are identities of factorial linear functionals, not analytic
integration formulas. Consequently they hold over the given K.

Write N_j=Lambda_2((x+y)f^j). Because f-k is homogeneous of degree two,
x partial_x f+y partial_y f=2(f-k).
Applying (2) to P=f^j and adding gives, for j>=1,
N_j=(2j+2)M_j-2kj M_(j-1).                                      (3)
For j=0 the same identity is N_0=2M_0; no M_(-1) is required.

## 2. The planes b=2a and b=2c

If b=2a, then partial_x f=2a(x+y). Applying (1) to f^(j+1),
and then using (3), gives the exact polynomial recurrence
M_(j+1)=4a(j+1)^2 M_j-4ka j(j+1)M_(j-1)+B_(c,j+1),  j>=1.       (4)
At j=0 it reads M_1=4a M_0+B_(c,1).

No division by c-a occurred. In particular (4) holds when c=a,
where f=a(x+y)^2+k has rank-one quadratic part.

If M_n through M_(n+3) were all zero, substituting j=n+1 and
j=n+2 into (4) would give
B_(c,n+2)=B_(c,n+3)=0.                                          (5)
This contradicts the imported two-moment theorem since ck!=0.

Exchanging x and y preserves Lambda_2 and proves the b=2c plane,
with boundary sequence B_(a,j).

## 3. The plane a=c

Here the two coordinate faces are identical, while
partial_x f+partial_y f=(2a+b)(x+y).
Add the two identities (1) for f^(j+1), divide by 2, and use (3):
M_(j+1)-B_(a,j+1)
 =(2a+b)(j+1)((j+1)M_j-kjM_(j-1)), j>=1.                        (6)
At j=0 it reads M_1-B_(a,1)=2a+b.

Again no division by a discriminant or by 2a+b occurs. Four zero
moments imply two adjacent zero B_(a,j), which is impossible.
In the degenerate direction b=-2a the identity even gives
M_j=B_(a,j) for all j>=1, so the stronger two-moment window holds
despite the four-monomial support. This proves Theorem 1.

## 4. The exact reduction for general b, and its limitation

Set
Delta=ac-b^2/4, D=a+c-b, alpha=c-b/2, beta=a-b/2.
The weighted derivative is
alpha partial_x f+beta partial_y f=2Delta(x+y).
Thus the same polynomial calculation gives
D M_(j+1)
 =4Delta(j+1)^2 M_j-4kDelta j(j+1)M_(j-1)
  +alpha B_(c,j+1)+beta B_(a,j+1), j>=1.                         (7)
This derivation is valid even when D or Delta is zero.

For general abck!=0, any four-zero block forces the two consecutive
weighted face cancellations
alpha B_(c,n+2)+beta B_(a,n+2)=0,
alpha B_(c,n+3)+beta B_(a,n+3)=0.                                (8)
The separate two-moment theorem for each face does not rule out (8).

For example, such a weighted face sequence can itself have adjacent
positive-index zeros without either individual face doing so. In C
take s=sqrt(21)>0, w=(s-6)/360, choose a with a^2=w, and set
c=-a, b=(s-5)/12, k=1.
All four coefficients are nonzero. One has
8640w^2+288w+1=0 and b=8w/(1+24w).
With T_r=alpha B_(c,r)+beta B_(a,r), native expansion gives
T_2=-b(1+24w)+8w=0,
T_3=-b(1+72w)+12w+1440w^2
   =4w(1+288w+8640w^2)/(1+24w)=0.                               (9)
The individual B_(a,2),B_(c,2) are nonzero: w<0, a is nonzero
pure imaginary, and they equal 1+24w plus or minus 4a.
Moreover M_1=1+b is nonzero. This is an obstruction only to an
unconditional adjacent-window assertion for the weighted forcing;
it is not a counterexample to Strong Factorial and does not settle
the actual general four-moment compatibility problem.

A remaining precise target is to exclude the simultaneous full
moment block together with its compatible face cancellations (8)
for coefficients outside the three planes in Theorem 1. No such
general exclusion is proved here.

## 5. Scope

The proof is a polynomial calculation over every characteristic-zero
field, together with the explicitly cited two-moment theorem for a
constant plus a square. No finite-prefix experiment is a universal
proof premise.

Only the stated even four-term sectors are closed. Full inhomogeneous
binary quadratics, arbitrary-support SFC(2), and SFC(3) remain outside
this proof. No novelty or priority assertion is made.
