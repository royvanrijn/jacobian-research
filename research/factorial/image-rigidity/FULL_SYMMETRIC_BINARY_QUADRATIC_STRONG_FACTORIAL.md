# Strong Factorial windows for the full symmetric binary quadratic

Let K be a characteristic-zero field, and define Lambda_2(x^r y^s)=r!s!
and Lambda_1(t^r)=r!.

**Theorem 1 (full symmetric six-term family).** If a,b,e,c in K satisfy
abec!=0, then

    f(x,y)=a(x^2+y^2)+bxy+e(x+y)+c                    (1)

has no six consecutive zero factorial moments at any positive start.
It therefore satisfies the Strong Factorial window for its actual six
nonzero monomials. The conclusion includes all coefficient-dependent
rank and critical-value degeneracies within this stated support.

The argument uses direct factorial identities and a compatible quadratic
face forcing. It does not assume the full univariate quadratic theorem,
any generic observability rank, positivity, or a finite moment cutoff.

## 1. A four-zero obstruction for the quadratic face forcing

For a quadratic P(t)=a t^2+e t+c with ec!=0, put
H_k=Lambda_1(P^k) and

    T_k=H_k+k e H_(k-1),   k>=1.                      (2)

Here a is allowed to be zero in this auxiliary statement. We prove that
T_k has no four consecutive zero values at any positive start.

Scale by e and write r=a/e, z=c/e and h_j=H_j/e^j. Thus z!=0 and
h_j=Lambda_1((r t^2+t+z)^j). Its exact native factorial recurrence is

    h_j=2r j(2j-1)h_(j-1)
       +j(j-1)(1-4rz)h_(j-2)+z^(j-1)(z+j), j>=2.    (3)

For completeness this identity is purely algebraic. Monomial evaluation
gives Lambda_1(R')=Lambda_1(R)-R(0). Put phi=r t^2+t+z and
J=Lambda_1(phi' phi^(j-1)). The identity
(phi')^2=4r phi+(1-4rz), followed by the factorial derivative identity
applied to phi' phi^(j-1), gives
J=2r(2j-1)h_(j-1)+(j-1)(1-4rz)h_(j-2)+z^(j-1).
Applying that identity to phi^j gives h_j=jJ+z^j, proving (3).
No coefficient-dependent denominator is used.

Let t_j=h_j+j h_(j-1)=T_j/e^j. Suppose

    t_k=t_(k+1)=t_(k+2)=t_(k+3)=0,   k>=1.

Set j=k+1, so j>=2. For each index i=j,j+1,j+2, the two zeros
t_(i-1)=t_i=0 give h_i=-i h_(i-1) and
h_(i-2)=-h_(i-1)/(i-1). Substitution in (3) yields

    E_i:=2r i(2i-1+2z)h_(i-1)+z^(i-1)(z+i)=0.        (4)

All divisions so far are by e or a nonzero rational integer.
In particular neither r nor 2i-1+2z is divided out.

For the consecutive equations E_i=E_(i+1)=0, multiply the first by
(i+1)(2i+1+2z) and the second by 2i-1+2z. Use h_i=-i h_(i-1)
and add. The h terms cancel exactly. Dividing only by z^(i-1), which
is nonzero, gives

    (i+1)(2i+1+2z)(z+i)
       +z(2i-1+2z)(z+i+1)=0.

Writing q=z+i, the left side simplifies to

    P_i(z)=2(z+i)^3+3(z+i)^2+i.                      (5)

Thus the three equations (4) imply both P_j(z)=0 and P_(j+1)(z)=0.
But

    P_(j+1)(z)-P_j(z)=6(z+j+1)^2.                    (6)

Since K has characteristic zero and is a field, these equalities force
z+j=-1. Then P_j(z)=j+1, a nonzero rational integer, a contradiction. This proves the four-zero
obstruction for T_k, including every vanishing value of the intermediate
quadratic pivots in (4).

As a direct identity, with B=2i-1+2z, the same cancellation is

    z^(i-1)P_i(z)=B t_(i+1)
      +2(i+1)[1-r(2i+1)B]t_i
      -i(i+1)(B+2)(1-4rz)t_(i-1).                    (6a)

It follows by substituting E_i=t_i-i(1-4rz)t_(i-1) into the elimination
above. It confirms the cancellation before imposing any zero hypotheses.

## 2. Interior recurrence without Hessian divisions

Let M_j=Lambda_2(f^j), R_j=Lambda_2((x+y)f^j), and let H_j be the face
moments in Section 1. Both coordinate faces of f are P(t). Set s=2a+b.
For every polynomial R the factorial identities

    Lambda_2(partial_x R)=Lambda_2(R)-Lambda_1(R(0,y)),
    Lambda_2(partial_x(xR))=Lambda_2(xR)                (7)

follow on monomials, and similarly for y.

The gradient sum f_x+f_y=s(x+y)+2e and the first identities in (7) imply

    j s R_(j-1)+2j e M_(j-1)=2(M_j-H_j),  j>=1.       (8)

The Euler identity x f_x+y f_y=2f-e(x+y)-2c and the second identities give

    R_j=(2j+2)M_j-j e R_(j-1)-2cj M_(j-1).           (9)

Use (8) at j+1, substitute (9), and replace j s R_(j-1) using (8) at j.
The two terms in (j+1)e M_j cancel; division by the rational integer 2
leaves

    M_(j+1)=s(j+1)^2 M_j
       +j(j+1)(e^2-sc)M_(j-1)
       +H_(j+1)+(j+1)e H_j,   j>=1.                 (10)

The coefficient s is never divided out. Neither a Hessian determinant,
2a-b, nor e^2-sc is divided out. Thus (10) holds unchanged when s=0,
when b=2a, when b=-2a, or when e^2=sc.

## 3. Exact support-window indices

Suppose M_n,M_(n+1),...,M_(n+5) were all zero, with n>=1. Apply (10)
at j=n+1,n+2,n+3,n+4. All three interior moments M_(j-1),M_j,M_(j+1)
in those equations lie in the assumed block. Hence

    T_(n+2)=T_(n+3)=T_(n+4)=T_(n+5)=0.

Section 1 excludes this four-zero forcing block. This proves Theorem 1
for every positive starting index. All identities hold over Q in their
coefficients and therefore specialize to every characteristic-zero field;
no embedding of K into the complex numbers or a residue field is assumed.

## 4. Scope and relation to other sectors

Under abec!=0, the six monomials x^2,y^2,xy,x,y,1 are all present. The two
quadratic diagonal coefficients are equal, and the two linear coefficients
are equal. These symmetry requirements are essential to the reduction to
one compatible face forcing; they are not silently removed.

The theorem does not settle the unrestricted five-term family
ax^2+bxy+dx+ey+c, the general six-term binary quadratic, full SFC(2),
or SFC(3). If a,b,e or c vanishes, the actual support is smaller and its
shorter Strong Factorial window is a separate assertion. No such lower
support window is inferred merely from this six-zero exclusion.

The proof includes the rank-one quadratic branches b=2a and b=-2a with
nonzero a, as well as the coefficient-dependent critical-value branch
e^2=(2a+b)c. No unverified rank, resultant, prime-existence, or ordinary
Factorial-to-Strong-Factorial implication is a premise.

The all-range proof is given above. Bounded native identity controls are
supplementary only; no formal-verification or literature-priority claim is made.
