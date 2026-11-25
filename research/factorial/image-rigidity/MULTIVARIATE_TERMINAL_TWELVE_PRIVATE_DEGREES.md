# Terminal twelve with private degree at least two

Let \(a,b\) be distinct nonnegative integer exponent columns. For
\(j=11,12\) put
\[
 Q_j(z)=\sum_{k=0}^j\binom jk\prod_i((j-k)a_i+kb_i)!\,z^k .
\]
The private degree is the sum of the positive exponents across rows
with a zero in the opposite column. Delete equal rows, and put
\(S_0=\sum d^2/(\ell+d)\) over the shared rows, with smaller exponent
\(\ell\ge1\) and gap \(d>0\).

**Theorem.** If the private degree is at least two, then \(Q_{11}\)
and \(Q_{12}\) have disjoint root arguments under independent positive
rescalings. In particular they are coprime over every
characteristic-zero field.

The [fixed-ray proof](MULTIVARIATE_TERMINAL_TWELVE_HIGH_CURVATURE.md)
covers private degree at least four and the high-strength range of
every remaining configuration. Up to column reversal and row
permutation, the table lists the rest. The displayed cutoff is the
upper endpoint of the range certified below.

| Increasing private gaps | Decreasing private gaps | Cutoff on \(S_0\) | Comparison parameter |
| --- | --- | ---: | --- |
| \(3\) | none | \(3/8\) | Shared curvature |
| \(1,1,1\) | none | \(9/8\) | Shared curvature |
| \(1,2\) | none | \(3/4\) | Shared curvature |
| \(1,1\) | \(1\) | \(7/8\) | Shared curvature |
| \(2\) | \(1\) | \(3/8\) | Shared curvature |
| \(1,1\) | none | \(27/8\) | Outer mean |
| \(1\) | \(1\) | \(13/4\) | Outer mean |
| \(2\) | none | \(23/8\) | Outer mean |

## Two-sided curvature bounds

A positive rescaling and scalar multiplication normalize the endpoints
of each polynomial to one. Its private core \(G_{j,k}\) at index \(k\)
contains \(\binom jk\) times the normalized factorial factors of the
private rows. Its remaining shared multiplier \(F_{j,k}\) satisfies
\[
 e^{-L S_0 k(j-k)/(2j)}\le F_{j,k}\le e^{-S_0 k(j-k)/(2j)}.
\]
The checker uses the upper endpoints \(1/8,1/4,3/8,7/8,5/4\)
with \(L=11/10,6/5,63/50,731/500,3/2\), respectively.
The coefficient for the endpoint \(7/8\) is \(731/500\).
At order twelve, \(291/200\) fails the lower multiplier inequality
for the shared row \((1,2)\) at its first interior coefficient.
The upper multiplier bound is the fixed-ray curvature bound.
For a row of strength \(s=d^2/(\ell+d)\), the second logarithmic
difference of its factorial coefficient is at most \(d^2/(j\ell)\)
by \(\log(1+x)\le x\). The endpoint chord then gives the lower multiplier
bound if \(1+d/\ell\le L\). Otherwise
\(s>\ell(L-1)^2/L\). Under the strength cap \(S_*\), the checker
enumerates only \(\ell\le S_*L/(L-1)^2\) and stops the gap enumeration
once \(d^2/(\ell_{\max}+d)>S_*\). It checks each lower multiplier
inequality at both orders by an integer factorial comparison with
a positive exponential Taylor lower bound. Multiplication gives the
stated bound on \(F_{j,k}\).

The coefficient boxes contain the actual arrays and the reference
\(G_{j,k}\exp[-S_0 k(j-k)/(2j)]\). A straight segment in the box
joins the two. Moving the reference down to \(S_0=0\) ends at the
private cores. At that anchor there is the exact disk certificate
described below.

## Outer-mean bounds

For the other three branches, use the shared reflected outer means
\(u,w\) at orders eleven and twelve. A row contributes
\[
 \prod_{h=1}^{d}\frac{11\ell+h}{11\ell+10d+h}
\]
to \(u^2\). The factor reciprocals are \(x=d/(11\ell+h)\) and
\(x'=11x/(12-xh/d)\). For the floor \(c=1/35\), rows containing any
\(x>c\) have \(\ell<d/(11c)-1/11\) and \(d<S_*(1+1/(11c))\).
The checker enumerates these rows. If a row contributes less than
the squared interval minimum for \(u\), that row cannot occur.
The upper reciprocal on the remaining rows below the floor is at
most \(11c/(12-c)\). This gives simultaneous caps at the two orders.

For \(u<1\) and lower reciprocal cap \(C\), the ratio
\(\rho=\log w/\log u\) satisfies
\[
 \frac{\log(1+(121/12)C)}{\log(1+10C)}
 \le\rho\le
 1-\frac{\log(1-C/12)}{\log(1+10C)} .
\]
Indeed, \(1+(121/12)x\le1+11x'\le(1+10x)/(1-x/12)\).
The lower ratio is decreasing in \(x\) by concavity in
\(\log(1+10x)\) and the upper ratio is increasing by the convex
numerator and concave denominator in its excess above one.
Thus the ratio bounds apply termwise in the logarithmic sums.

The [reflected mean and skew formulas](MULTIVARIATE_TERMINAL_NINE_PRIVATE_ROWS.md#shared-outer-mean-boxes)
give the coefficient exponents in the same \(u\). The normalized
mean loss lies between its value at the reciprocal cap and
\(e_{j,k}=k(j-k)/(j-1)\). The normalized absolute skew is at most
its value at the cap. The checker certifies the skew convexity
polynomials at orders eleven and twelve on \(0\le x\le1/2\), which includes
every retained reciprocal. It then takes outward enclosures of
the logarithms and powers. The coefficient boxes contain the actual
arrays and comparison arrays with shared factors \(u^{e_{11,k}}\)
and \(u^{(121/120)e_{12,k}}\), respectively. This follows from the
checked inclusion of each comparison exponent in its exponent
interval. Connecting by straight coefficient segments and sliding
\(u\) down to zero gives the binomial anchors \(1+z^{11}\), \(1+z^{12}\).
At \(u=1\) all shared unequal rows are absent.

## Exact continuation

Constant and leading coefficient dominance for each private core
confines all roots to \(1/8<|z|<8\). All boxes have the same endpoint
coefficients, with smaller interior upper bounds. The checker uses
ten rays with continuous piecewise affine cosines.
The [Bernstein projection test](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md#moving-rays-and-exact-replay)
excludes roots along the rays. The checker verifies the upper
root counts at the private anchors using disjoint rational Rouché
disks. At the binomial anchors it verifies alternating exact
cosine signs. In increasing argument the eleven upper sectors have
one root each, of orders \(12,11,12,11,12,11,12,11,12,11,12\).
There is one additional negative real root of order eleven.
Join the first and last sectors to their reflected lower sectors.
The single order-eleven root in the last joined sector stays
negative real by conjugation and the annulus bound. Strict
positivity on the real line for the actual order-twelve polynomial
follows from its positive Gamma integral of the twelfth power of the
real binomial. Hence the actual order-twelve pair in that sector is
nonreal. These sector counts give disjoint arguments, and rational
Bézout gives the characteristic-zero assertion.

The [exact checker](verify_multivariate_terminal_twelve_private_degrees.py)
certifies 315 row inequalities, ten skew polynomial inequalities,
4,965 curvature cells, and 55,500 outer-mean cells. Run with SymPy
available.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_private_degrees.py
~~~
