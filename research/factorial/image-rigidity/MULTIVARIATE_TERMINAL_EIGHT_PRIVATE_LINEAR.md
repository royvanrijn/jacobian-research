# Terminal-eight separation for one private linear row

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
consist of a single private linear row and any number of shared
rows. Their seventh and eighth moment polynomials have disjoint
root arguments under independent positive rescalings, and are
coprime over every characteristic-zero field.

Together with the
[larger private-degree theorem](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_ROWS.md),
this covers every ordinary profile with a private row at terminal
eight. The [checker](verify_multivariate_terminal_eight_private_linear.py)
uses exact rational and outward dyadic arithmetic.

## Ordinary reciprocal caps

Orient the private row as \((0,1)\). Its normalized multiplier
times the binomial coefficient is
\[
 B_{j,k}=\binom jk\frac{k!}{(j!)^{k/j}} .
\]
Write \(S_0=\sum d_i^2/(\ell_i+d_i)\) over the shared unequal
rows. If \(S_0\ge21/4\), multiply the private coefficients by
the shared upper bounds \(\exp[-S_0 k(j-k)/(2j)]\).
At \(21/4\) these caps pass the fixed rays of the
[terminal-eight high-curvature test](MULTIVARIATE_TERMINAL_EIGHT_HIGH_CURVATURE.md);
the same boxes cover larger \(S_0\).

For \(S_0\le21/4\), let \(u,w\) be the shared outer means at
orders seven and eight. On a shared row with lower base \(\ell\)
and gap \(d\), the Gauss reciprocals are \(x=d/(7\ell+h)\),
\(x'=d/(8\ell+h)\), \(1\le h\le d\). The row's squared
contribution to \(u\) is
\[
 u_{\ell,d}^2=\prod_{h=1}^d\frac{7\ell+h}{7\ell+6d+h}.
\]
All other outer factors are at most one, so \(u\le u_{\ell,d}\).
Set \(t=1/40\). A row with \(d/(7\ell+1)>t\) satisfies
\(\ell<(d/t-1)/7\), and \(S_0\le21/4\) then implies
\(d<(21/4)(1+1/(7t))\). The checker lists these candidate
rows. At lower endpoint \(u_-\) it uses only rows with
\(u_{\ell,d}^2\ge u_-^2\), and takes each order's reciprocal
cap as the larger of \(t\) and the largest reciprocal over
the remaining candidates. Rows omitted from the enumeration
by the floor condition have both reciprocals at most \(t\).
This applies regardless of the orientation, the number
of rows, and the magnitudes of their bases.

Put \(U=-\log u,\ W=-\log w,\ \rho=W/U\) if \(u<1\).
Since \(x'=7x/(8-xh/d)\),
\[
 1+\frac{49}{8}x\le 1+7x'\le \frac{1+6x}{1-x/8}.
\]
If \(c\) is the order-seven cap, the bound on \(\rho\) is
\[
 \frac{\log(1+(49/8)c)}{\log(1+6c)}
 \le \rho\le 1-\frac{\log(1-c/8)}{\log(1+6c)} .
\]
For the lower bound the logarithm in the numerator is a
concave function of \(\log(1+6x)\) and vanishes at zero.
For the upper bound \(-\log(1-x/8)\) is convex and
\(\log(1+6x)\) concave, with both zero at zero.
The corresponding quotients are monotone and the bounds
sum over the shared factors. The interval contains \(49/48\).
If \(u=1\), there are no shared unequal rows and both shared
arrays are entirely one.

## Reflected means and skew links

Use the
[reflected mean-loss and skew formulas](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_LINEAR.md#the-two-normalized-coefficient-arrays)
at orders seven and eight. Specifically set \(m=j-1\),
\(g_h=\log(1+h x_j)\), \(s=g_m/2\), and
\(r=\min(k,j-k)\). Relative to \(s\), the reflected mean
loss is the sum of
\(\operatorname{artanh}((1-2h/m)\tanh s)/s\) for \(0\le h<r\).
Each summand decreases in \(s\). Thus the ratio is between its
value at the cap and \(e_{j,k}=k(j-k)/(j-1)\).
The local skew is the negative of
\[
 \sum_{h=0}^{j-1}\omega_h g_h,\qquad
 \omega_h=\tfrac12\big(\mathbf 1_{\{r\le h<j-r\}}-(j-2r)/j\big)
\]
for the increasing orientation. For every skew magnitude used
below, the checker verifies convexity in \(\log(1+m x_j)\)
by expanding the numerator
\[
 \sum_{h=1}^{m-1}\omega_h h(m-h)
       \prod_{\substack{1\le l<m\\l\ne h}}(1+l x_j)^2
\]
of its second derivative, after removing positive factors.
The same calculation applies with any linked linear combination
of the weights. On \([0,1]\) a negative term \(b x_j^i\),
\(i\ge2\), can be bounded below by \(b x_j^{i-2}\).
After this replacement the remaining coefficients are
nonnegative. The function and its first derivative vanish
at zero. Its magnitude relative to \(s\) therefore increases,
and reversing a row changes only the skew sign.

In the independent coefficient boxes, apply these reflected mean
and skew bounds to the interval of exponents in \(u\).
At order eight multiply the bounds by the corresponding bounds on
\(\rho\). The checker verifies that the symmetric reference
exponents \(e_{7,k}\) and \((49/48)e_{8,k}\) lie in every such
box. In the paired range \(3/10\le u\le17/20\), preserve the
following more precise coefficient links.

Write \(\alpha,\beta,\gamma\) for the signed skews of reflected
pairs at indices \(1,2,3\). At order seven let
\(L=\alpha-\beta\), \(K=\alpha-3\beta+5\gamma\), so that
\[
 \alpha=(5\gamma+3L-K)/2,\qquad
 \beta=(5\gamma+L-K)/2 .
\]
These are the
[septic skew links](MULTIVARIATE_TERMINAL_SEVEN_NEAR_FLAT.md#linked-order-seven-coefficients).
At order eight let \(L=8\alpha-7\beta\), \(K=3\alpha-7\beta+7\gamma\).
Then
\[
 \alpha=(7\gamma+L-K)/5,\qquad
 \beta=(56\gamma+3L-8K)/35 .
\]
For one factor, with \(E\) the shift operator on \(g_h\), the order-eight
links are
\[
 L=(E-1)^4(E+1)(5E^2+4E+5)g_0/4,\qquad
 K=(E-1)^6(E+1)g_0/4 .
\]
These identities follow by expanding the reflected skew weights.
The derivative-numerator test above uses the weights of the
negative signed \(\gamma,L,K\) in each case. The cap evaluations
then bound the absolute values of their sums over the factors.

Write \(G_a(s)=as-\operatorname{artanh}(a\tanh s)\).
At order seven the reflected mean defects relative to the
Gaussian powers at indices \(2,3\) are \(H,(9/5)H-J\), where
\(H=\sum G_{2/3}\) and \(J=\sum((4/5)G_{2/3}-G_{1/3})\).
At order eight the defects at indices \(2,3,4\) are
\(H,2H-J,(12/5)H-J-D\), where
\[
 H=\sum G_{5/7},\qquad J=\sum(G_{5/7}-G_{3/7}),\qquad
 D=\sum((2/5)G_{5/7}-G_{1/7}).
\]
These formulas follow by summing the individual mean differences
\(\operatorname{artanh}((1-2h/m)\tanh s)\).
The individual \(G_a\) are convex, vanish at zero, and have
derivatives \(a(1-a^2)\tanh^2s/(1-a^2\tanh^2s)\).
For \(a>b>0\) and \(c=b(1-b^2)/(a(1-a^2))\), the derivative
of \(cG_a-G_b\) is a positive constant times
\(\tanh^4s/[(1-a^2\tanh^2s)(1-b^2\tanh^2s)]\). It increases
with \(s\). This proves the same convexity and cap-based strength
ratio bound for all the linked defect differences.

## Quadratic projection boxes and continuation

On a subinterval centered at \(u_0\), use the private
coefficients times \(u_0^{e_{7,k}}\) or
\(u_0^{(49/48)e_{8,k}}\) as reference. The shared log correction
at each coefficient is a sum of linear terms whose common
scalars are the signed links, the nonnegative defects, and
one outer-mean scalar. At order seven that last scalar is
\(\log(u/u_0)\); at order eight it is
\((49/48)\log(u/u_0)+(49/48-\rho)U\).
Both have the corresponding common Gaussian vector
\((e_{j,k})_k\). For log correction \(\sigma_k\) with
\(|\sigma_k|\le M_k<1\), use
\[
 \left|\exp(\sigma_k)-1-\sigma_k-\frac{\sigma_k^2}{2}\right|
 \le \frac{M_k^3}{6(1-M_k)} .
\]
The series bound uses \(n!\ge6\) for \(n\ge3\).
The model keeps the linear scalars and their quadratic products
grouped before projecting, and bounds the remainder
coefficientwise.

The six ordered ray cosines interpolate affinely in \(u\) at
rational knots. The private-only upper caps give strict
dominance of the constant and leading terms at radii \(1/10\)
and \(10\), respectively, for the actual and reference arrays
and their straight segments. For each ray the
[paired Bernstein test](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md#projection-replay)
or its independent coefficient version certifies one strict
projection sign on each radial cell. The checker subdivides
the outer-mean interval if needed. Logarithms and exponent
endpoints are enclosed by finite positive series and integer
comparisons.

At fixed \(u\), the coefficient boxes contain the straight
segment from the actual arrays to the symmetric shared
reference with ratio \(49/48\). The paired model
contains the segment since it allows independent intervals
for the linear scalars, products, and remainders. Then decrease
\(u\) to zero on the reference arrays using the continuous
rays. Since every interior reference exponent is positive,
the anchors at zero are \(1+t^7,1+t^8\). The rational
binomial sign test gives one upper root per sector, of orders
\(8,7,8,7,8,7,8\), with a remaining negative real septic
root. The high-strength anchors have the same sector counts.
Join the first and last sectors to their respective conjugates
when continuing the counts. The single septic root in the last
joined sector stays real by conjugation. Its eighth-order pair
is nonreal since the actual eighth moment is strictly positive
for real coefficient ratios by the Gamma integral. This proves
argument separation and rational Bézout gives coprimality.

The replay verifies 82 high-strength radial cells and 4,533
low-strength radial cells. The reciprocal-cap enumeration uses 1,373
candidate rows.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eight_private_linear.py
~~~
