# Fully shared light rows at terminal ten

Let two distinct ordinary exponent columns have all unequal rows shared.
Write \(u\) for their normalized order-nine outer reflected mean. In
each band below, suppose every reciprocal \(d/(9\ell+1)\) is at most
the listed cap, where \(\ell\) is the lower row exponent and \(d\)
is the gap.

| Outer mean \(u\) | Reciprocal cap |
| --- | ---: |
| \([0,.52]\) | \(1/8\) |
| \([.52,.62]\) | \(99/1000\) |
| \([.62,.70]\) | \(3/40\) |
| \([.70,.76]\) | \(11/200\) |

Then the order-nine and order-ten factorial moment polynomials have
disjoint root arguments under independent positive rescalings.
Consequently they are coprime over every characteristic-zero field.

## Ordinary caps and exact coefficient enclosures

The order-nine squared outer mean of a shared row is
\[
 u_{\ell,d}^{\,2}
  =\prod_{h=1}^d\frac{9\ell+h}{9\ell+8d+h}
  \le\left(\frac{9\ell+d}{9\ell+9d}\right)^d.
\]
Each row mean is at most one, and their product is \(u\), so each
row has squared mean at least \(u_-^2\) in an interval starting
at \(u_-\). If its reciprocal exceeds \(t>0\), then
\(9\ell+1<d/t\) and
\[
 u_{\ell,d}^{\,2}<\left(1+\frac{8t}{1+t}\right)^{-d}.
\]
With \(t=1/40\) these inequalities give a finite integer search
above the floor whenever \(u_->0\). The checker takes the
maximum of the floor and the exact reciprocals of all candidates
with squared mean at least \(u_-^2\) and within the declared cap.
For example the extra rational knots .641 and .701 allow the candidate
lists to drop a row once its squared mean is too small.

Match the successive-order reciprocals
\(x=d/(9\ell+h)\), \(x'=9x/(10-hx/d)\). For an order-nine cap \(c\),
the order-ten cap is \(9c/(10-c)\). Write
\(s_j=\log(1+(j-1)x_j)/2\), \(T_j=\tanh s_j\), \(U_j=\sum s_j\).
For \(u<1\), the ratio \(\rho=U_{10}/U_9\) lies between
\[
 \frac{\log(1+(81/10)c)}{\log(1+8c)}
 \quad\text{and}\quad
 1+\frac{\log(1/(1-c/10))}{\log(1+8c)}.
\]
Indeed \(1+(81/10)x\le1+9x'\le(1+8x)/(1-x/10)\), and the
logarithm quotient calculation in the
[outer-mean bounds](MULTIVARIATE_TERMINAL_NINE_PRIVATE_ROWS.md#shared-outer-mean-boxes)
applies with the displayed constants. The interval contains \(81/80\).
With \(P_j=\sum_i\epsilon_iT_{j,i}^2/2\), the
[matched tanh lemma](MULTIVARIATE_MATCHED_TANH_MOMENTS.md#matched-tanh-moments-at-consecutive-orders)
transfers \(P_9\) to \(P_{10}\) with center coefficient
\((81/80)^2(1+\eta^2)\) and error at most
\(2(81/80)^2\eta\sum_iT_{9,i}^2/2\), where \(\eta=c/(20+79c)\).

The cap-based independent coefficient bounds use reflected means
and skews. For one factor at order \(j\), put \(m=j-1\),
\(p_h=1-2h/m\), \(T=\tanh s\), and \(r=\min(k,j-k)\).
The mean loss at distance \(r\), divided by \(s\), is
\(\sum_{h<r}\operatorname{artanh}(p_hT)/\operatorname{artanh}T\).
This ratio is nonincreasing in \(T\) for \(p_h\in[0,1]\),
by comparing the positive odd series whose term ratios
\(p_h^{2l+1}\) are nonincreasing in \(l\). Its upper endpoint
as \(T\downarrow0\) is \(e_{j,k}=k(j-k)/(j-1)\).
The absolute skew has even-power coefficients proportional to
\(-\sum_h\omega_h p_h^{2l}\), where
\(\omega_h=(\mathbf1_{\{r\le h<j-r\}}-(j-2r)/j)/2\).
They are nonnegative because the central indices have smaller
absolute \(p_h\) than the exterior indices. For \(q\ge2\),
\(T^q/\operatorname{artanh}T\) increases when \(T^2\le1/2\),
as follows by multiplying its derivative sign expression
by \(1-T^2\) and using \(\operatorname{artanh}T\ge T\).
Consequently the skew-to-\(s\) ratio is at most its cap value.
The checker bounds both quantities by outward rational logarithms
and incorporates the ratio \(\rho\) at order ten.

## Grouped moment and moving-ray comparison

For \(0\le u\le .04\) the checker uses the independent boxes
from the preceding paragraph. Their coefficient rectangles
contain each actual array and the Gaussian reference of
powers \(e_{9,k}\), \((81/80)e_{10,k}\) at its \(u\). On this
interval the binomial anchors \(1+z^9\), \(1+z^{10}\) at \(u=0\)
give the alternating root counts.

In the remaining bands use the
[matched skew and mean groups](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_LIGHT.md#matched-skew-parameter)
with the order-nine and order-ten constants above.
Exchange the columns if necessary so the normalized signed
second moment at order nine is in \([0,1]\). The reference and
error calculations below apply on rectangles in this signed
parameter and \(u\). Up to \(u=.38\) the checker retains three
skew and two mean-defect vectors, with a quadratic exponential
model. Above .38 it retains four skew vectors and three mean
defects and keeps the cubic products grouped as well.

If the first \(a\) skew and the first \(b\) mean-defect vectors
are retained, the corresponding omitted logarithmic corrections
at index \(k\) are bounded by
\[
 E_k=U_j\left[
 \frac{T_*^{\,2a+1}\sum_h|\omega_h|}{(2a+2)(1-T_*^2)}
 +\frac{rT_*^{\,2b+2}}{(2b+3)(1-T_*^2)}
 \right].
\]
Here \(T_*\) is the cap in the corresponding order. This follows
by summing the geometric series and using \(s_j\ge T_j\).
The retained scalars use the same upper skew and lower mean
bounds as the matched-light comparison at terminal nine.
For the cubic exponential model, if \(|\sigma-\tau|\le E_k\)
and \(|\sigma|,|\tau|\le A<1\), then
\[
 \left|e^\sigma-\left(1+\tau+\tfrac12\tau^2+\tfrac16\tau^3\right)\right|
 \le \frac{A^4}{24(1-A)}+E_k(1+A+A^2/2).
\]
This is the geometric bound for the fourth and higher
exponential terms, together with the derivative bound for the
cubic at intermediate points.

The [common coefficient and ray displacement](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_UPPER.md#common-coefficient-and-ray-displacement)
uses the same normalized shift \((u-u_0)/u_0\) for Gaussian
coefficients and the affine ray cosine. Each strict Bernstein
projection sign is checked with outward intervals. The ordinary
cap may vary by outer interval, but the actual profile and
Gaussian reference at a shared endpoint obey both applicable caps. The rays within
each band are ordered, continuous and piecewise affine in \(u\).

The checker verifies constant and leading-term domination
for the annulus in each branch. At fixed \(u,\rho\) it uses the
matched logarithmic homotopy to zero skew and mean defects,
then changes \(\rho\) to \(81/80\) and moves the Gaussian
reference to the lower endpoint of its band. The lower
mean-defect bound remains valid under a common scaling of
all the defects and skews, since its two sides scale as
\(t\) and \(t^{2l}\), respectively. Each nonzero lower endpoint
has independent exact rational Rouché disks for its anchor.

The alternating upper sectors have orders
\(10,9,10,9,10,9,10,9,10\). The sole remaining real anchor root
is negative and belongs to order nine. Along the continuation
it stays real in the last sector joined to its reflection.
Strict positivity of the actual order-ten moment on the real
axis follows from its Gamma integral, so the actual tenth-order
pair there is nonreal. The actual root arguments are disjoint,
and rational Bézout gives the characteristic-zero conclusion.

## Reproduction

Run the [checker](verify_multivariate_terminal_ten_shared_light.py)
from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_ten_shared_light.py --native
~~~

The optional --native flag compiles the [outward-dyadic C++17
implementation](verify_multivariate_terminal_ten_shared_light.cpp) in a temporary directory. The pure Python replay
uses the same interval signs. A native input outside its explicit
safe integer bounds uses the Python implementation. To reduce
native mantissas the checker may rescale independent groups
by powers of two, multiplying the parameter by the same power
and dividing the coefficient vector outward by that power.
The fused base and the shared displacement group are kept at
their common original scale.
