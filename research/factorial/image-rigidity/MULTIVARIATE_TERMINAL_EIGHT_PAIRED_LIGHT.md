# Paired shared terminal-eight separation from outer mean .10 to .48

Let all unequal rows of distinct ordinary exponent columns be shared.
On them write \(\ell_i=\min(a_i,b_i)\), \(d_i=|b_i-a_i|\), and let
\(u\) be the normalized order-seven outer geometric mean. If
\[
        .10\le u\le .48,\qquad\max_i d_i/(7\ell_i+1)\le1/8 ,
\]
the order-seven and order-eight factorial moment polynomials have
disjoint root arguments under independent positive rescalings.
Rational Bézout gives coprimality over every characteristic-zero field.

## Cross-order skew link

For a matched Gauss pair write \(x=d/(7\ell+h)\) and
\(x'=d/(8\ell+h)=7x/(8-xh/d)\), where \(1\le h\le d\).
At \(x\le1/8\), the ratio \(x'/x\) belongs to \([7/8,8/9]\).
The magnitudes of the innermost reflected skews of orders seven
and eight in this pair are
\[
\begin{aligned}
 E(x)&=[6g_3-g_1-g_2-g_4-g_5-g_6]/14,\\
 F(x')&=[3g'_3+3g'_4-g'_1-g'_2-g'_5-g'_6-g'_7]/8,
\end{aligned}
\]
where \(g_h=\log(1+hx)\) and \(g'_h=\log(1+hx')\).
Their orientation signs agree, by the paired log formulas.
Both magnitudes vanish with zero derivative at zero and are
convex in the local strength \(s=\tfrac12\log(1+(j-1)x_j)\)
on the reciprocal interval \([0,1/8]\).
Indeed, if the log weights before their positive normalizing
division are \(a_h\), their second derivative in \(s\) is a
positive factor times
\(\sum_h a_h h(j-1-h)/(1+h x_j)^2\).
The checker certifies the positive numerator by rational
Bernstein coefficients on \([0,1/8]\). Thus each skew
magnitude is increasing and its ratio to \(s\) increases.

The [checker](verify_multivariate_terminal_eight_paired_light.py)
verifies the exact comparisons
\[
       \tfrac{189}{100} E(x)\le F(7x/8),\qquad
       F(8x/9)\le\tfrac{1977}{1000}E(x).
\]
In each derivative, clear the positive linear denominators.
The lower comparison numerator, after its zero constant,
has positive coefficients up to degree three and negative
coefficients thereafter. It has at most one positive zero
by Descartes' rule. The lower comparison difference therefore
takes its minimum at one of the endpoints of \([0,1/8]\).
Its value is zero at zero and strictly positive at \(1/8\)
by the exact logarithm check. For the upper comparison,
the derivative numerator after removal of a factor of \(x\)
has strictly positive Bernstein coefficients on that interval.
This proves both comparisons.

Let \(\gamma_7,\gamma_8\) be the signed total inner skews.
Summing over both orientations gives
\[
 \big|\gamma_8-\tfrac{1933}{1000}\gamma_7\big|
             \le\tfrac{44}{1000}\sum_i E(x_i).
\]
The right side is bounded by the outer logarithm of order seven
times the exact ratio \(E/s\) at the reciprocal cap.

## Models and replay

Write \(U=-\log u\), \(W=-\log w\) for the outer means of orders
seven and eight. When \(x\le t\), the bound on \(x'\) and the
logarithm quotient comparisons used in the
[terminal-seven paired proof](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md)
give
\[
 x'\le 7t/(8-t),\qquad
 \frac{\log(1+(49/8)t)}{\log(1+6t)}
 \le \rho=W/U \le 1-\frac{\log(1-t/8)}{\log(1+6t)}.
\]
This follows factorwise from
\(1+(49/8)x\le1+7x'\le(1+6x)/(1-x/8)\).
The Gaussian ratio \(49/48\) lies in the verified interval.

Use the [paired septic mean and skew groups](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md#coefficient-enclosures)
and the [linked octic mean and skew formulas](MULTIVARIATE_TERMINAL_EIGHT_NEAR_FLAT.md#octic-linked-means-and-skews).
For each local cap \(t_j\), put
\(s_t=\frac12\log(1+(j-1)t_j)\). Terms bounded by single-factor
powers use \(\sum x_j^k\le(t_j^k/s_t)\sum s\).
For the leading mean defects use the convex \(G_a\) ratio
at \(s_t\). The remaining mean defects are fifth-order
differences whose derivative bounds use \(\tanh s\le s\le s_t\).
All log quotient intervals are evaluated by the atanh series
with a geometric remainder.

At a parameter midpoint \(u_0\), the Gaussian normalized
coefficients are \(u_0^{e_{7,k}}\) and \(u_0^{(49/48)e_{8,k}}\),
where \(e_{j,k}=k(j-k)/(j-1)\). Keep the common mean shift
and every linked skew or shape correction as a separate shared
scalar. The correction \(\sigma\) in each coefficient uses
\(e^\sigma=1+\sigma+R\), where
\(0\le R\le A^2/[2(1-A)]\) for \(|\sigma|\le A<1\).
The checker subdivides both \(u\) and the normalized
order-seven innermost skew. It suffices to take
\(\gamma_7\ge0\) after exchanging the columns.
The conditional bound above cuts down the primary order-eight
skew interval. Every linked contribution is projected onto a
rotated imaginary part and converted to Bernstein coefficients
before multiplying by its scalar interval. A strictly signed
Bernstein array excludes that ray on its radial cell.

All interior normalized multipliers are at most one, hence
\((13/12)^8<2\) proves the annulus \(1/12<|t|<12\).
For fixed \(u,\rho\), decrease all skews and defects to zero
by linear interpolation in their logarithmic coefficient
corrections. The skew transfer and the other linear bounds
persist along the path. Its actual and Gaussian endpoint
multipliers are at most one, as are their geometric
interpolations. Then vary \(\rho\) to \(49/48\) and move \(u\)
to .10 using the Gaussian arrays and the rational affine rays.

At the Gaussian anchor put \(X=t+t^{-1}\). The septic
quotient after division by \(t+1\) is
\(X^3-X^2-2X+1+a(X^2-X-1)+b(X-1)+c\), and the octic
quotient is \(X^4-4X^2+2+d(X^3-3X)+e(X^2-2)+fX+g\).
The letters are the first interior coefficients at each order.
Signs at the six ray abscissae and \(X=\pm2\) prove that the
successive upper sectors contain orders \(8,7,8,7,8,7,8\).
The additional septic root is real. The exclusions preserve
these counts. In the final sector joined to its conjugate,
the actual octic pair is nonreal by the positive eighth-power
factorial moment integral, and the unique septic root is real
by conjugacy. This proves disjoint arguments.

Replay from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_eight_paired_light.py
~~~
