# A matched-skew Cayley separator at terminal thirteen

Let two distinct ordinary exponent columns have only shared unequal rows.
Let \(u\) be the endpoint-normalized, binomial-stripped reflected outer
mean at order twelve. If \(93/100\le u\le467/500\) and all order-twelve
Gauss reciprocals \(d/(12\ell+h)\) are at most \(1/73\), the order-twelve
and order-thirteen factorial moment polynomials have disjoint root
arguments under independent positive rescalings. Their coprimality
holds over every characteristic-zero field by rational Bézout.

The first separator depends on the same normalized skew for both
orders. This allows the projected enclosures to follow the linked
coefficients through the tight first angular gap.

## Transport on the scaled Cayley grids

Write \(w\) for the order-thirteen outer mean and set
\(p=1-u,\ q=1-w,\ r=q/p,\ U=-\log u,\ W=-\log w\) and \(R=144/143\).
In the [weight-nine model](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY.md#weighted-cayley-model)
the signed weight-two parameters at the two orders are
\[
 P=\frac{M_2}{2p^2}, \qquad P'=\frac{M'_2}{2q^2},
\]
where \(M_2,M'_2\) are the signed tanh sums with the same row
orientation convention at both orders. The
[excess-loss transport](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
gives
\[
 \left|P'-\left(\frac Rr\right)^2P\right|
 \le \frac{R^2}{q^2}\frac{52}{11}B E^2(W-RU),
 \qquad
 B=\frac{143c}{24+143c},\qquad
 E=\frac{11(1+13c)}{24+143c}.
\]
Here \(c\) is the refined lower reciprocal cap in the
[terminal-thirteen scale enclosure](MULTIVARIATE_TERMINAL_THIRTEEN_CAYLEY.md#scales-and-grouped-coefficient-bounds).
On a subband the right-hand side is bounded with the lower endpoint
of \(q\) and \(W-RU\le(L-R)U\), using an outward rational upper
bound for \(U\). The checker intersects the resulting upper-order
skew interval with the physical bound \([-CA'/2,CA'/2]\), with \(A'\)
the cap on the normalized order-thirteen tanh factors.

Partition the order-twelve interval \([-CA/2,CA/2]\) into sixteen
equal parts for the first ray. In each part use the chosen order-twelve
interval and its transported order-thirteen interval to apply the
[conditional higher moment bounds](MULTIVARIATE_TERMINAL_TWELVE_CAYLEY_LOWER.md#condition-the-logarithmic-moments-on-the-leading-skew).
The first ray uses
\[
                 b=b_0(p)+\frac{49}{2500}P^2
\]
as its positive lower-order height parameter, with \(b_0(p)\)
interpolated between the rational knots in the checker. The
upper-order height parameter is \(b/\sqrt r\). The other ten
ray levels depend only on \(p\). The checker verifies their strict
order and that even the largest admissible first-ray height lies in
the physical left half-plane. Each projection cell can subdivide
the lower skew interval as well as the ratio, deficit scale and
real Cayley coordinate. The transformed coefficients preserve all
weight-nine grouped directions.

## Continuation and exact replay

At fixed \(p,r\), multiply the logarithmic corrections at both
orders by \(s\in[0,1]\). Their skew parameters become \(sP,sP'\),
which satisfy the same transport error bound. The conditional
moment bounds remain valid along this path. The positive heights
of the corresponding rays vary continuously, including the first
ray with \(P\) replaced by \(sP\). The far and near sides of the
angular sectors therefore keep the same root counts, since the
projections exclude roots on the moving rays and the checked leading
dominance keeps all roots in the disk.

At \(s=0\) the coefficients are Gaussian and the first ray uses
\(b_0(p)\). Move the Gaussian ratio to \(R(1-(R-1)p/2)\) and
follow this curve to \(p=.066\). These Gaussian arrays are covered
by \(P=P'=0\). At the anchor, the sign counts place the
upper-half-plane roots in alternating orders thirteen and twelve,
starting at the largest height. As in the
[terminal-thirteen Cayley count](MULTIVARIATE_TERMINAL_THIRTEEN_CAYLEY.md#root-counts-and-replay),
the order-thirteen root in the joined sector at the negative real
axis is real by conjugation, while the even order is strictly
positive on that axis. This gives disjoint arguments in every sector.

The [exact checker](verify_multivariate_terminal_thirteen_linked_cayley.py)
certifies 96,819 projection cells across 88 ray and order checks.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_linked_cayley.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_thirteen_linked_cayley.py 13 0 68 69
~~~
