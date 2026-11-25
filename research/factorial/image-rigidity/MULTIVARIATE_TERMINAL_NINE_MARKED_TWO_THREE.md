# Terminal-nine angular separation with a marked (2,3) row

Let two distinct ordinary exponent columns have only shared unequal
rows. Suppose a row has entries \(2,3\) in either orientation. Let
\(u\) be the order-eight outer geometric mean after positive endpoint
normalization and removal of the binomial factors. If \(u\ge19/25\),
the order-eight and order-nine root arguments are disjoint under
independent positive rescalings. The dimension and the remaining
shared exponents are unrestricted. Rational Bézout gives coprimality
in characteristic zero.

## The residual reciprocal bound

The marked row contributes \(u_*^2=17/24\) to the squared
outer mean. Since \(19/25>17/24\), two such rows would give
\(u\le17/24\), so the marked row is unique. Remove it and let
\(v=u/u_*\) be the residual order-eight mean. The inequality
\((19/25)^2>(17/24)(903/1000)^2\) gives \(v>903/1000\).
Swapping both columns reverses both root arguments by reciprocal
substitution up to positive rescaling. Thus take the marked row to
have increasing orientation.

For a residual row with lower exponent \(\ell\) and integer gap
\(d\ge1\), the order-eight reciprocals are
\(x_h=d/(8\ell+h)\), \(1\le h\le d\).
The row's squared outer mean bounds \(v^2\) from above.
\[
 v^2\le\prod_{h=1}^d\frac{8\ell+h}{8\ell+7d+h}.
\]
For \(d=1\), the maximum over \(\ell\le3\) is \(25/32\),
below \((903/1000)^2\). Hence every residual unit row has
\(\ell\ge4\) and \(x_h\le1/33\).
For \(d\ge2\), put \(y=d/(8\ell)\). Each ratio in the product
is at most \((1+y)/(1+8y)\). Thus
\[
  v\le \,\left(\frac{1+y}{1+8y}\right)^{d/2}
    \le \frac{1+y}{1+8y},
  \qquad x_h\le\frac{1-v}{8v-1}<\frac1{33}.
\]
In addition \(v^{-2}\) is the product of the factors \(1+7x_h\)
over all residual rows. At a tested subinterval with left endpoint
\(v_-\), the checker therefore uses
\[
                 c=\min\bigl(1/33,(v_-^{-2}-1)/7\bigr).
\]
Every matched order-nine reciprocal is at most \(8c/(9-c)\).

## The two linked coefficient arrays

The marked row's normalized coefficients, including binomial factors,
are
\[
 C_{j,k}=\binom jk
 \frac{(2j+k)!}{((2j)!)^{1-k/j}((3j)!)^{k/j}},\qquad j=8,9.
\]
Use these as the cores of the grouped logarithmic model.
At a residual midpoint \(v_0\) the references multiply the cores
by \(v_0^{e_{8,k}}\) and \(v_0^{(64/63)e_{9,k}}\),
where \(e_{j,k}=k(j-k)/(j-1)\).
The cap \(c\) bounds the ratio of the adjacent residual outer
logarithms by the matched-factor comparison. The checker keeps the linked outer vector,
the three skew vectors, the two mean-defect vectors, and all
quadratic products as separate interval scalar contributions.
The truncated logarithm and exponential Taylor errors are bounded
per coefficient. The core factorial roots and all reference powers
are bracketed by exact integer powers. The closed interval ending
at \(v=1\) also encloses the limiting core coefficients.
[Grouped logarithmic model](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections).
[Matched-factor comparison](MULTIVARIATE_TERMINAL_NINE_UPPER_BRIDGE.md#bounds-on-ordinary-reciprocals).

These bounds also hold during the following homotopy. Add the same
real parameter \(s\ge0\) to both entries of every residual row,
leaving the marked row fixed. The integer gaps and row orientations
are preserved. The Gauss-factor expressions continue to hold
for the real shifted exponents. In an increasing row their factors
follow from
\[
 \frac{\Gamma(j(\ell+s)+dk+1)}{\Gamma(j(\ell+s)+1)}
    =\prod_{h=1}^d\prod_{m=0}^{k-1}
      \bigl(j(\ell+s)+h+dm\bigr)
\]
and the same identity at \(k=j\) gives the endpoint normalization.
For a decreasing row reverse the coefficient indices.
The order-eight reciprocals decrease with \(s\), so their initial
\(1/33\) cap persists. The product identity for \(v\) still holds
and now makes \(v\) increase to one. All residual normalized
coefficients converge to one. Including the limit at \(s=\infty\)
gives a continuous path to the marked core.

## Cayley exclusions and the anchor

Use \(p=1/8\), \(\lambda=\sqrt p\), and
\(t=(\lambda X-1)/(\lambda X+1)\). Transform each grouped
contribution by the integer Cayley matrix and the powers of
\(\lambda\). The checker proves leading-term dominance on
\(|X|=11/4\). The physical image
of this disk lies in the open left half-plane since
\[
                    p(11/4)^2=121/128<1 .
\]

The seven ray levels \(b>0\) are strictly decreasing in their index
and piecewise affine in \(v\) between rational knots. In physical
coordinates their arguments are \(\pi-2\arctan(\lambda b)\).
Inside the Cayley disk, write \(X=x+iy\) on an upper ray. Its
equation is
\[
                  (y-b)(1+pby)=-pb x^2,qquad 0<y\le b .
\]
The checker narrows the height with three interval updates and
tests the full range \(-11/4\le x\le11/4\).
It excludes zero on each box by a strict real projection of
the transformed polynomial, keeping the scalar groups intact.
It subdivides any unresolved box in \(v\) or \(x\).

At \(v=1\) the checker uses rational disks for the transformed
marked core. With \(D_h\) denoting the Taylor coefficients at a
disk center and \(r\) its radius, the tested inequality is
\[
        |D_0|+\sum_{h=2}^j |D_h|r^h < |D_1|r .
\]
Rouché's theorem gives one root per disk and the disks are
disjoint. Four upper roots of each order lie in alternating
sectors of orders \(9,8,9,8,9,8,9,8\). The remaining nonic root
is real by conjugacy of its disk. The ray exclusions preserve
these counts along the residual homotopy, with the last upper
and lower sectors joined across the negative real axis. The
joined sector has two octic roots and one nonic root. The
nonic root stays real by conjugation. The actual octic has no
real root because the eighth power of a real linear combination
of distinct monomials, with its first coefficient nonzero,
has strictly positive factorial moment. The actual octic roots
in the joined sector therefore have different arguments from
the real nonic root. The other sectors separate all other roots.

The replay certifies 327,949 Cayley projection cells across 350
order, ray, and subband checks. The exact checker also verifies
the two core anchors and the per-ray cell counts. Run the full
certificate or a partition from the repository root.
[Exact checker](verify_multivariate_terminal_nine_marked_two_three.py).

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_marked_two_three.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_marked_two_three.py 903 924 9
~~~
