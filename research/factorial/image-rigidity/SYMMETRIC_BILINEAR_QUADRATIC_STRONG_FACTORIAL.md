# Symmetric bilinear quadratic Strong Factorial windows

Let K be a characteristic-zero field and let Lambda_2(x^r y^s)=r!s!.
The assertions below concern the displayed families only, not arbitrary binary quadratics.

**Theorem 1.** For alpha,e,c in K with alpha e c!=0, the polynomial

    g(x,y)=alpha xy+e(x+y)+c                           (1)

has no four consecutive zero factorial moments at any positive start.
Thus g satisfies its actual support-four Strong Factorial window.

## 1. A precise adjacent-zero theorem for affine forcing

Set U_k=Lambda_1((e z+c)^k) and Q_k=U_k+k e U_(k-1), k>=1.
Native factorial evaluation gives

    U_k=k e U_(k-1)+c^k,
    Q_k=2U_k-c^k,
    Q_(k+1)=(k+1)e Q_k+c^k(c+(k+1)e).                (2)

We claim that Q_k=Q_(k+1)=0 is impossible for k>=2. A pair of zeros
would force c/e=-(k+1), since c,e are nonzero. Put z=c/e. Direct expansion is

    Q_k/e^k=P_k(z),
    P_k(z)=z^k+2 sum_(j=0)^(k-1) (k!/j!) z^j in Z[z]. (3)

At z=-(k+1), this integer is nonzero. If k is even, z is odd and P_k(z)
is odd, because all lower coefficients in (3) are even.

If k is odd and k>=3, z is even. Write v for the 2-adic valuation on
nonzero rational integers. The constant term has valuation

    v(2k!)=1+k-s_2(k).

For 1<=j<k, the valuation of the jth summand exceeds that constant valuation
by j v(z)-v(j!)>=s_2(j)>=1. The leading term has valuation k v(z)>=k,
which is strictly larger than 1+k-s_2(k), because an odd k>=3 has at least
two set bits. Hence the constant term is the unique least-valuation summand,
and the integer sum cannot be zero. Its nonzero rational value remains
nonzero in every characteristic-zero field. This proves the claim.

The exception k=1 is real: c=-2e gives Q_1=Q_2=0 and Q_3=4e^3!=0.
We therefore do not assert adjacent forcing coprimality at start one.

## 2. A self-contained binary moment recurrence

Write M_j=Lambda_2(g^j) and R_j=Lambda_2((x+y)g^j). For every polynomial P,
monomial evaluation proves Lambda_2(partial_x P)=Lambda_2(P)-Lambda_1(P(0,y))
and Lambda_2(partial_x(xP))=Lambda_2(xP), with corresponding y identities.

The two boundary polynomials of g are both e z+c. Summing the first identity
for g^j therefore gives

    j alpha R_(j-1)+2j e M_(j-1)=2(M_j-U_j).          (4)

Summing the second identity and using
x partial_x g+y partial_y g=2g-e(x+y)-2c gives

    R_j=(2j+2)M_j-j e R_(j-1)-2cj M_(j-1).           (5)

Use (4) at j+1, substitute (5), and use (4) at j to remove R_(j-1).
The M_j terms involving e cancel, yielding

    M_(j+1)=alpha(j+1)^2 M_j
       +j(j+1)(e^2-alpha c) M_(j-1)+Q_(j+1), j>=1.  (6)

No coefficient-dependent nonzero invariant other than the stated parameters
is divided out. In particular e^2=alpha c is allowed.

If M_n,...,M_(n+3) were zero, (6) at j=n+1,n+2 would give
Q_(n+2)=Q_(n+3)=0. Their start n+2 is at least three, contradicting Section 1.
This proves Theorem 1 for all n>=1 over every K.

## 3. The earlier five-term sector has the same moments

**Corollary 2.** For b,e,c nonzero, the genuine five-term polynomial

    f(x,y)=(b/2)x^2+bxy+e(x+y)+c                      (7)

also has no FOUR consecutive zero factorial moments. This is stronger than
its actual support-five Strong Factorial window.

Set u=x/2+y and v=x/2. Then
f(x,y)=2buv+e(u+v)+c and x+y=u+v. The positive quadrant in x,y maps to
u>=v>=0 with Jacobian dxdy=2 dudv. Since the transformed phase is symmetric
in u and v, integration with factorial weights gives, for every integer m>=0,

    Lambda_2(f^m)=Lambda_2((2buv+e(u+v)+c)^m).         (8)

Indeed twice the integral over u>=v>=0 is the integral over the whole
quadrant; the diagonal has measure zero. This proves a polynomial identity
with rational coefficients in b,e,c, so (8) specializes to every
characteristic-zero field. Theorem 1 with alpha=2b proves the corollary.

Since alpha=2b is a bijection on nonzero parameters in every characteristic-zero
field, identity (8) also proves the converse between the universal assertions:
Theorem 1 and the FOUR-window assertion of Corollary 2 are equivalent.
Both directions use the same window length four. This does not identify a
support-five window with a support-four window.

## 4. Scope

Theorem 1 requires the two linear coefficients to be equal. It does not
settle alpha xy+d x+e y+c when d!=e, the unrestricted five-term family
ax^2+bxy+dx+ey+c, the full six-term binary quadratic problem, SFC(2), or
SFC(3). No assertion about worldwide priority, formal verification, or a
finite-prefix implication is made. The written proofs above are all-range arguments; bounded native controls are supplementary only.