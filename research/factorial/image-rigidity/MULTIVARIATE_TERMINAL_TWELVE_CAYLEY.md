# Terminal-twelve Cayley separation for light shared rows

Let two distinct ordinary exponent columns have only shared unequal
rows. Let \(u\) be the order-eleven outer mean after positive
endpoint normalization and removal of the binomial coefficients.
At least one row is unequal, so \(0<u<1\). The following conditions
on \(u\) and the largest order-eleven Gauss reciprocal
\(x=d/(11\ell+h)\) suffice for disjoint root arguments of the
order-eleven and order-twelve factorial moment polynomials.

| \(u\) | Reciprocal cap |
|---|---:|
| \([19/20,1)\) | Unrestricted |
| \([23/25,19/20]\) | \(1/90\) |
| \([22/25,23/25]\) | \(1/70\) |
| \([171/200,22/25]\) | \(1/60\) |
| \([19/25,171/200]\) | \(1/55\) |

The claim holds under independent positive rescalings. The two
polynomials are coprime over any characteristic-zero field by
rational Bézout.

## Comparison scales

Let \(w\) be the order-twelve outer mean, \(p=1-u\),
\(q=1-w\), and \(r=q/p\). With \(R=121/120\), the
[matched ratio bound](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
gives
\[
 R\left(1-\frac{(R-1)p}{2u}\right)\le r.
\]
Set \(c\) to the smaller of the specified cap and \((u^{-2}-1)/10\).
For \(c>0\), the upper log ratio \(\rho=\log w/\log u\) is
\[
 \rho\le L=1-\frac{\log(1-c/12)}{\log(1+10c)}.
\]
Indeed \(x'=11x/(12-(h/d)x)\) and \(1+11x'\le
(1+10x)/(1-x/12)\). The ratio of the two logarithms is
nondecreasing in \(x\), by convexity of \(-\log(1-x/12)\)
and concavity of \(\log(1+10x)\). On each subband the checker
takes an outward rational upper bound on \(L\) at the largest
admissible \(c\), reaching \(R\) in the limit \(c=0\).
It verifies \(R\le L<2\). Integrating the tangent bound for
\(s^{L-1}\) on \(u\le s\le1\) gives
\[
  r\le (1-u^L)/(1-u)\le L(1-(L-1)p/2).
\]
In particular, on \(p\in[p_-,p_+]\) both the actual \(r\)
and the reference ratio \(r_*(p)=R(1-(R-1)p/2)\) belong to
the interval
\[
 \big[R(1-(R-1)p_+/(2(1-p_+))),\;
       L(1-(L-1)p_-/2)\big].
\]
For the upper comparison with \(r_*\), the function
\(a\mapsto a(1-(a-1)p_-/2)\) increases on \(1\le a\le2\)
and \(0\le p_-\le.24\).

The closest order-twelve rays in selected subbands use a
shorter subdivision and the following integer refinement of \(c\).
When \(d=1\), keep the largest \(1/(11\ell+1)\) with
\(\ell\ge1\) and \(1/(11\ell+1)\le c\). Otherwise, let
\(x_1=d/(11\ell+1)\). The other reciprocals on the same row
are at least \(x_1/(1+x_1)\); since \(d\ge2\),
\[
 u^{-2}\ge(1+10x_1/(1+x_1))^2.
\]
This gives \(x_1\le(1-u)/(11u-1)\). The maximum of these two
bounds, with the second also capped by \(c\), covers every row.
The other rays permit cap \(1/55\) on all of the middle band.

## Weighted Cayley model

Use the [grouped Cayley expansion](MULTIVARIATE_TERMINAL_ELEVEN_MIDDLE_CAYLEY.md)
with weight nine. At order \(j\) use \(h=p\) or \(q\) and
\(\lambda=\sqrt h\). For an outward deficit bound \(q_0\)
let \(C=(1-q_0/2)/(1-q_0)\). The row variables \(T_i\)
satisfy \(\sum_i T_i\le-\log(1-h)\le C h\).
Writing \(T_*\) for the bound from the reciprocal cap and
\(h_-\) for the lower deficit, put \(A=C\) if \(h_-=0\) and
\(A=\min(C,T_*/h_-)\) otherwise. Then \(T_i\le A h\) and
each normalized logarithmic moment of weight \(d\) has magnitude
at most \(C A^{d-1}/d\). The checker retains weights \(2\)
through \(9\) in the logarithm and every exponential monomial
of weight at most \(9\), with its exact multinomial factor.
An even weight carries a signed moment. The coefficient interval
is nonnegative when every signed factor occurs to an even power.

The omitted even and odd logarithmic series first occur at weights
\(10\) and \(11\). Dividing by \(h^{10}\), the respective bounds
are \(\omega C A^9/[10(1-(q_0A)^2)]\) and
\(k_*q_0 C A^{10}/[11(1-(q_0A)^2)]\), where
\(k_*=\min(k,j-k)\) and \(\omega\) is the sum of the absolute
even logarithmic coefficients at position \(k\).
The exponential slack includes the terms of weight at least ten
among powers one through four. Higher powers are bounded by
\(V^5/[5!(1-V)]\), where \(V\) bounds the retained logarithm at
\(h=q_0\). The omitted logarithmic tail is handled by
\(|e^s-e^t|\le |s-t|/(1-M)\) when \(|s|,|t|\le M<1\).
The resulting bound before the Cayley change is proportional to
\(h^{10}\). After the Cayley normalization the remaining powers
are \(\lambda^{20+m-j}\) for coefficient index \(m\). Each
retained grouped direction cancels its exact initial Taylor
coefficients before division by \(\lambda^j\). On the closest
rays the Gaussian direction also keeps the change in \(\lambda\)
as one common scalar, with the derivative bounded by the
convergent binomial series.

The weight increase matters at degree twelve. The weight-seven
estimate starts with \(\lambda^{16-j}\) in the constant
coefficient, and its interval slack can swamp the near-real ray
margins even when the true arrays remain separated.

## Separation and replay

The exact checker proves leading dominance in the scaled Cayley
disk on each certified ray cell. The Cayley transformation is
\(z=(\lambda Z-1)/(1+\lambda Z)\). The rays have
arguments \(\pi-2\arctan(\sqrt p\,b)\) for positive
piecewise affine levels \(b\). The checker verifies that they
are in the physical left half-plane. The positive height equation
and reflection are those of the
[order-eleven grouped Cayley proof](MULTIVARIATE_TERMINAL_ELEVEN_MIDDLE_CAYLEY.md#ray-certification-and-continuation).
When the whole disk need not map into the physical left half-plane,
the checker clips the real coordinate to \(0\le\Re Z\le1/\lambda\)
on such a ray. Reflection and conjugation cover \(\Re Z<0\).

At fixed actual \(p,r\), the grouped enclosure is convex and
contains its pure Gaussian array as well as the actual array.
Join each actual array to its Gaussian reference at these same
deficits. Then move the Gaussian ratio to \(r_*(p)\) at fixed
\(p\), and follow \(r_*(p)\) to the first deficit of the band.
The anchor has exact alternating signs on the imaginary axis,
with orders \(12,11,12,11,12,11,12,11,12,11,12\) in its
successive upper sectors. Order eleven also has a negative real
root. At the zero deficit endpoint these signs are computed for
the limiting normalized Cayley polynomials. Leading dominance
and the strict projection certificates preserve the counts
along the paths. Across the positive real ray the first upper
and lower sectors contain only order-twelve roots. Across the
negative real ray the odd root remains real by conjugation; the
even real factorial moment is strictly positive by its Gamma
integral, so its two remaining roots are nonreal.

Run the [exact checker](verify_multivariate_terminal_twelve_cayley.py)
to replay the ray and subband certificates. The exact replay checks 404,228 projection cells across 960 order, ray, and subband checks.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_cayley.py
python3 research/factorial/image-rigidity/verify_multivariate_terminal_twelve_cayley.py seam 12 0 50 55
~~~
