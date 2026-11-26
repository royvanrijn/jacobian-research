# Quadratic polynomials in x+y: exact weighted factorial windows

Let K be a characteristic-zero field, Lambda_1(t^r)=r!, and
Lambda_2(x^r y^s)=r!s!. For g(t)=A t^2+B t+C put
N_j=Lambda_1(g(t)^j), W_j=Lambda_1(t g(t)^j), j>=0.

**Theorem 1 (the aligned one-dimensional family).** Every nonzero
f(x,y)=g(x+y), deg(g)<=2, satisfies Strong Factorial with its actual
number of nonzero binary monomials, over every K. More precisely:

| Nonzero coefficients of g | Actual support of f | Sufficient moment-window length |
| --- | ---: | ---: |
| A,B,C | 6 | 4 |
| A,B; C=0 | 5 | 2 |
| A,C; B=0 | 4 | 3 |
| A only | 3 | 1 |
| B,C; A=0 | 3 | 3 |
| B only | 2 | 1 |
| C only | 1 | 1 |

This is a structured binary family, not arbitrary binary quadratics
or SFC(2). The coefficients of x^2 and y^2 coincide, that of xy is
twice either one, and the two linear coefficients coincide.

## 1. Two exact scalar identities

On monomials, Lambda_2((x+y)^r)=(r+1)!, so
Lambda_2(g(x+y)^j)=W_j.                                           (1)

The univariate factorial identities
Lambda_1(P')=Lambda_1(P)-P(0),
Lambda_1((tP)')=Lambda_1(tP)
give
N_j=j B N_(j-1)+2A j W_(j-1)+C^j, j>=1,                         (2)
W_j=(2j+1)N_j-BjW_(j-1)-2CjN_(j-1), j>=1.                       (3)
For (3), use t g'=2g-Bt-2C. These are polynomial identities in A,B,C,
not analytic integrations. They specialize to every K and every
degenerate coefficient stratum.

## 2. B,C nonzero: four weighted moments cannot vanish

Scaling the phase by B does not change any zero moment, so normalize
B=1 and replace A,C by A/B,C/B. In this section C!=0.

Suppose W_(j-1)=W_j=0. By (2)--(3),
(2C-2j-1)N_j=2C^(j+1).                                         (4)

Now suppose W_n=W_(n+1)=W_(n+2)=0, n>=1, and write j=n+1.
Put d=2C-2j-1. The two adjacent pairs give
d N_j=2C^(j+1),
(d-2)N_(j+1)=2C^(j+2),
while (2) gives N_(j+1)=(j+1)N_j+C^(j+1).
Multiplying the second equality by d, substituting, and cancelling
the nonzero factor C^(j+1) yields
(d-2)(2j+2+d)=2Cd.
Since 2j+2+d=2C+1, this simplifies exactly to
2C+2j+3=0, or C=-n-5/2.                                        (5)

Four consecutive zero W moments would apply (5) at n and n+1,
forcing two different values of C. Thus every four-term window
of W is nonzero. This proof does not require A!=0; when A!=0
the actual binary support is six and the four-window is stronger
than the required six-window.

## 3. B=0 with A,C nonzero: a three-window

If W_n=W_(n+1)=W_(n+2)=0, (2) gives
N_(n+1)=C^(n+1), N_(n+2)=C^(n+2).
Then (3) at n+2 gives
0=W_(n+2)=(2n+5)C^(n+2)-2C(n+2)C^(n+1)=C^(n+2),
a contradiction. Hence no three consecutive W moments vanish.
The actual binary support here has four terms.

## 4. C=0 with A,B nonzero: a two-window

If W_n=W_(n+1)=0, (3) gives N_(n+1)=0 and (2) gives
N_(n+1)=(n+1)B N_n. Thus N_n=N_(n+1)=0.

For completeness, these last adjacent unweighted zeros are impossible
by an elementary recurrence. Define
y_m(z)=sum_(i=0)^m (m+i)! z^i/[2^i(m-i)!i!].
Native expansion gives
N_m=m! B^m y_m(2A/B).
The coefficient identity
(m+i)(m+i+1)-(m-i)(m-i+1)=2i(2m+1)
proves
y_(m+1)=(2m+1)z y_m+y_(m-1), y_0=1, y_1=1+z.
Two adjacent zero values propagate backward to y_0=0, impossible.
This uses no irreducibility theorem or extension of the given field.
The actual binary support is five, so the two-window suffices.

## 5. A=0 with B,C nonzero: the actual three-window

The linear phase has
N_j=jB N_(j-1)+C^j,
W_j=N_j+jB W_(j-1).
Eliminating N_(j-1) yields, for j>=2,
W_j=2jB W_(j-1)-j(j-1)B^2 W_(j-2)+C^j.                          (6)
Three consecutive zeros of W, substituted at their last index,
would imply C^j=0. Thus the actual three-term binary support has
its required three-window.

## 6. Single coefficients and the final scope

If g(t)=U t^d, U!=0 and d in{0,1,2}, then
W_j=U^j(dj+1)! is nonzero for every j>=1.
This includes nonzero constants. The cases exhaust every nonzero
g of degree at most two. Any window longer than a proved nonzero
short window is also nonzero, so the table proves Theorem 1 with
the actual support size.

Every proof step is a finite polynomial identity or a backward
recurrence. No finite-prefix experiment proves a universal assertion,
and no full binary quadratic or higher-dimensional SFC closure follows
from this structured family.
