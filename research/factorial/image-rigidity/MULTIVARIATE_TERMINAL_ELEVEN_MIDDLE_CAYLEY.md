# Terminal-eleven Cayley separation from outer mean .76 to .88

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be the order-ten reflected outer mean after positive
endpoint normalization and removal of the binomial factors. Suppose
\(19/25\le u\le22/25\) and the order-ten Gauss reciprocals
\(x=d/(10\ell+h)\) are at most \(1/55\). Then the tenth and eleventh
factorial moment polynomials have disjoint root arguments under
independent positive rescalings. Their coprimality holds over every
characteristic-zero field.

## Coupled scales and logarithmic moment bounds

Let \(w\) be the order-eleven outer mean, \(p=1-u\),
\(q=1-w\), and \(r=q/p\). Write \(c=\min(1/55,(u^{-2}-1)/9)\)
for a lower-order reciprocal cap. The second cap follows from
\(u^{-2}=\prod_i(1+9x_i)\). The matched upper reciprocal cap is
\(10c/(11-c)\). The [matched outer-ratio inequality](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
gives the lower deficit ratio bound
\[
 r\ge \frac{100}{99}\left(1-\frac{1-u}{198u}\right).
\]
For the upper bound, \(1+10x'\le(1+9x)/(1-x/11)\) termwise.
The convex numerator \(-\log(1-x/11)\) and concave denominator
\(\log(1+9x)\) have value zero at \(x=0\), so their ratio is
nondecreasing. Thus
\[
 \rho=\frac{\log w}{\log u}
 \le 1-\frac{\log(1-c/11)}{\log(1+9c)}.
\]
The inequality \(1-u^\rho\le \rho(1-u)\) gives the same upper
bound for \(r\). The checker verifies the resulting rational
enclosures
\[
              \frac{10084}{10000}\le r\le\frac{10110}{10000},
              \qquad q\le q_0=\frac{243}{1000}.
\]

At either order let \(h=p\) or \(q\), and set
\(T_i=\tanh(\log(1+(j-1)x_i)/2)\) for the order-\(j\)
reciprocals. These satisfy
\[
 \sum_i T_i\le-\log(1-h)\le C h,
 \qquad C=\frac{1-q_0/2}{1-q_0}.
\]
The second inequality bounds \(1/k\) by \(1/2\) in the positive
logarithm series for \(k\ge2\). On each subband the checker uses
\(A=\min(C,T_*/h_-)\), where \(T_*\) comes from the reciprocal
cap and \(h_->0\) is the lower deficit bound for that order.
Hence \(T_i\le A h\) and for every \(d\ge2\),
\(\sum_i T_i^d\le CA^{d-1}h^d\). This also bounds the absolute
values of the signed moments.

These bounds apply to the
[stable weight-seven Cayley expansion](MULTIVARIATE_TERMINAL_TEN_NEAR_FLAT.md#a-grouped-taylor-enclosure-through-deficit-weight-seven).
Each scalar moment gets the bound \(CA^{d-1}/d\) after
division by \(h^d\). The two geometric tails at weight eight
have numerators \(CA^7/8\) and \(Cq_0 A^8/9\), with denominator
\(1-(q_0A)^2\). Products of the scalar bounds give the polynomial
enclosure for the retained exponential terms and the raw
remainder. The Cayley transform cancels the low Taylor degrees
in each retained direction before division by \(\lambda_j^j\),
where \(\lambda_j^2=h\). At the closest eleventh-order ray in
the later subbands, the checker also bounds the derivatives of
the pure Gaussian coefficient vectors in \(\lambda_j\).
It keeps the increment in \(\lambda_j\) as a single common
scalar. The derivative bounds follow from differentiating the
convergent binomial series. In a term
\(\lambda_j^m\sum_{i\ge t}a_i h^i\), with \(|a_i|\le B\),
the corresponding majorants are \(Bh^t/(1-h)\) and its
derivative in \(h\), multiplied by the differentiated powers
of \(\lambda_j\).

## Ray certification and continuation

The checker proves leading coefficient dominance in the
Cayley disk \(|Z|<13/5\). The map is
\(z=(\lambda_j Z-1)/(1+\lambda_j Z)\).
Its pole cannot be a zero of the transformed polynomial,
because the leading coefficient of the original polynomial
is nonzero. The nine rays have positive piecewise affine
levels \(b\) in \(p\), with \(b_j=b\) at order ten and
\(b_j=b/\sqrt r\) at order eleven. Their physical arguments are
\(\pi-2\arctan(\sqrt p\,b)\). The checker verifies \(h b_j^2<1\),
so the rays lie in the physical left half-plane. The whole
Cayley disk may extend beyond that half-plane.

Writing \(Z=X+iY\), a ray in the physical upper half-plane
has \(Y>0\), \(h(X^2+Y^2)<1\), and height equation
\[
 h b_j Y^2+(1-h b_j^2)Y=b_j(1-h X^2).
\]
There is at most one nonnegative height, since the left side
is strictly increasing in \(Y\ge0\). The checker restricts to
\(|X|\le 1/\sqrt h\) and encloses the positive branch by
iterating \(Y=b_j-hb_j X^2/(1+hb_jY)\) on \([0,b_j]\).
It checks \(X\ge0\). Column reversal replaces the physical
polynomial by its reciprocal, so reversal followed by conjugation
replaces \(X\) by \(-X\) in the upper Cayley half-plane.
The coefficient enclosures allow both choices of orientation.

For each actual \(p,r\) the grouped coefficient enclosure is
convex and contains the reference with its logarithmic defects
zero at the same deficits. A straight segment joins them.
The reference then moves in \(p\) at fixed \(r\) to the
anchor \(p=3/25\). At that anchor exact imaginary-axis signs
give one root per upper sector, of orders
\(11,10,11,10,11,10,11,10,11,10\), and an additional negative
real root of order eleven. Throughout the continuation the
checker excludes roots on the rays by centered interval
projection. Join the first and last sectors to their
reflections. The single odd-order root in the joined last
sector stays negative real by conjugation and the nonzero
constant coefficient. The even-order factorial moment is
positive on the real line by its Gamma integral, so the two
even roots in that sector have nonreal arguments. This proves
the claim.

The [exact checker](verify_multivariate_terminal_eleven_middle_cayley.py)
certifies 180,337 projection cells in 432 order, ray, and
subband checks.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_middle_cayley.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_eleven_middle_cayley.py 11 0 135 240
~~~
