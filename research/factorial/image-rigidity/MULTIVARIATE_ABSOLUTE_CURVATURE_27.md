# All-order factorial angular separation at absolute curvature 27

For distinct nonnegative integer exponent columns \(a,b\), omit
equal rows and put
\[
 S=\sum_i \frac{(b_i-a_i)^2}{\max(a_i,b_i)},\qquad
 V=\sum_{a_ib_i>0}\frac{(b_i-a_i)^2}{\min(a_i,b_i)}
             +\sum_{a_ib_i=0}\max(a_i,b_i),\qquad W=1+V.
\]

**Theorem.** If \(S\ge27\), every consecutive pair of ordinary
factorial moment polynomials has disjoint root arguments under
independent positive rescalings. Therefore they are coprime over
every characteristic-zero field.

The proof refines the
[absolute curvature 64 estimate](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md).
A sharper count of connected hypergraphs reduces the growth
bound on the truncated endpoint logarithm. A fixed absolute
error suffices for the final angular comparison, so the length
of the truncation can depend on \(N e^{-S/2}\) alone.

## Two independent coefficient regions

Use the same endpoint normalization, writing its coefficients
as \(c_{m,k}\) at orders \(m=n,N=n+1\). The endpoint entries are
one and
\[
 0<c_{m,k}\le\binom mk\exp[-S k(m-k)/(2m)] .                    \tag{1}
\]
The case \(N=2\) is immediate: the linear first moment has a
negative real root and the second moment is strictly positive
on the real line by its integral representation as an even power.

Split the coefficients into left and right halves as in
the curvature 64 proof, so \(T_m(z)=L_m(z)+z^m U_m(1/z)\).
For \(3\le N\le256\), \(\binom mk\le(em/k)^k\) and the
convexity of \(1-\log w-(S/2)(1-w)\) on \([1/m,1/2]\)
give
\[
  c_{m,k}^{\,1/k}\le \max\{em e^{-(27/2)(1-1/m)},\,2e e^{-27/4}\}
             <1/128\qquad(1\le k\le m/2).
\]
The first expression as a function of \(m\in[2,256]\) has
no interior maximum, since the derivative of its logarithm
has sign \(m-27/2\). The reflected bound gives the same
estimate for the other half. On \(|w|\le e^{1/m}\) either
half satisfies \(|H-1|<37/200\) and \(|wH'|<1/5\)
by summing a geometric series with ratio at most \(e/128\).

For \(N\ge256\) and \(S\ge2\log N+5\) the
[large-curvature slope theorem](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md)
applies. In the remaining region write
\[
       R=\log N,\qquad D=N e^{-S/2}.
\]
Here \(R>11\), \(n>59000\), and \(S/n<1/2100\), since
\(S/(e^{(S-5)/2}-1)\) decreases for \(S\ge27\).
If \(D\le1/6\), then for \(1\le k\le8\), (1) gives
\[
 c_{m,k} e^{k/m}\le (101/600)^k/k!.
\]
In fact \(c_{m,k} \le D^k e^{S k^2/(2m)}/k!\) and
\(\exp(S k/(2m)+1/m)<101/100\). For \(k\ge8\) up to
\(m/2\), the same convex bound as above has exponentiated
endpoint bounds \(e D e^{4S/m}/8<1/16\) and
\(2e e^{-27/4}<1/16\). Thus on \(|w|\le e^{1/m}\), the terms
with \(k\ge8\) are bounded by \(15^{-k}\). Summing gives the
same half bounds \(|H-1|<37/200\) and \(|wH'|<1/5\)
in this region.

In either region, Rouché on \(|z|=e^{1/m}\) and on the
reciprocal circle confines the roots to
\(|\log|z||<1/m\). The domination is
\(e(1-37/200)>1+37/200\).
The arguments of the half values start at zero on the
positive radius and their derivatives in angle have
modulus at most \((1/5)/(1-37/200)<1/4\).
At roots of the two degrees with a common upper argument
\(\theta>0\), subtract their two phase equations. The half
argument errors differ by less than \(\theta\), so an
integer multiple of \(2\pi\) would be strictly between zero
and \(2\theta\le2\pi\). This proves separation in both regions.

## The endpoint cluster region

Assume now \(N>e^{(S-5)/2}\) and \(D>1/6\), and put
\[
 K=\lceil32D\rceil+12,\qquad \alpha=2L/n,\qquad
              \beta=8eW/n .
\]
Here \(L=\max(1,\max_i |b_i-a_i|/a_i)\le W\) in an
orientation with all starting exponents positive.
The bounds \(S\le V\le S(S+2)\) and \(n>59000\) give
\[
 \frac Dn<\frac1{720000},\qquad
 \frac{WD}{n}<\frac{108}{100000},\qquad
 \frac{WK}{n}<\frac{297}{2500},\qquad K\ge18 .                 \tag{2}
\]
For the second bound use the decrease of
\((S+1)^2e^{-S/2}\) for \(S\ge27\).
The third follows from \(K\le32D+13<110D\).
Both halves have degree at least \(K\), since
\(K/n<297/(2500W)\) and \(W\ge28\).
For \(K\le k\le m/2\), (1) gives
\(c_{m,k}\le8^{-k}\) by testing its convex bound at the
two ends. At \(k=K\) use
\(eD e^{SK/(2m)}/K\le (e/32)e^{297/5000}<1/8\), and at
\(k=m/2\) use \(2e e^{-27/4}<1/8\).

If the starting endpoint has a private row, the
[private coefficient bound](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation)
gives \(c_{m,k}\le e^{-(S-5)k/4}\) for \(k\le m/2\).
Thus on \(|w|\le6/5\),
\[
  |H_m-1|\le\sum_{k\ge1}((6/5)e^{-11/2})^k<1/200.
\]
The logarithm with value zero at zero has norm less than
\(1/190\) in this case.

At a positive starting endpoint, use the
[endpoint interpolation](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation)
with interpolated first factor \(B_m\) at real \(m\in[n,N]\).
Equation (2) ensures all the hypotheses of its polygamma
and finite-difference bounds (5)–(7): the arguments are at
least \(ma_i/2\) or \(m/2\), these lower bounds exceed \(K+1\),
\(\alpha K<1/2\), and \(4W/n<1/2\) follows from \(D>1/6\)
and (2). In particular the edge weights \(w_r(m)\) for
the formal logarithm and their order derivatives satisfy
\[
 |w_r|\le M_r,\qquad |w_r'|\le(2r/n)M_r,\qquad
             e M_r=(r-2)!\alpha^{r-2}\beta .
\]
The coefficient of degree \(k\le K\) in that logarithm is
\(B_m^k\Sigma_k(m)/k!\), with \(\Sigma_k(m)\) the connected
sum on \(k\) labels from the same interpolation.

We use the explicit bound for the first-factor derivative
from the [curvature 64 proof](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md#uniform-half-logarithm-bounds).
Its derivation by products and the Gamma chord only needs
\(m\ge n\) and the same lower bounds on the Gamma arguments.
It gives
\[
 \left|(\log B_m)'-1/m\right|
 \le \frac{5V+S(3+R)}{m^2}<\frac1n .
\]
Indeed the numerator is at most \(40R^2\) by
\(S<2R+5\) and \(R>11\), and \(40R^2/(e^R-1)<1\)
on that range. Integrating the resulting derivative bound
and using (1) at \(k=1,m=n\) gives
\(B_m\le D\exp(S/(2n)+2/n)<101D/100\).

## Sharper connected sum

For \(2\le k\le K\) bound the absolute values of edges by
\(eM_r\). The total of these majorants over all possible
edges on \(k\) labels is at most
\[
 \beta k(k-1)\sum_{r=2}^k
            \frac{(\alpha k)^{r-2}}{r(r-1)}
 \le \frac{\beta k(k-1)}{2(1-\alpha K)} .
\]
Indeed \( \binom kr(r-2)!\le k(k-1)k^{r-2}/(r(r-1))\).
Thus the cost of all edges additional to a discovery
collection, divided by that collection's weight, is at most
the exponential of the displayed bound.

Use the rooted discovery tree with blocks of children as in the
[connected-hypergraph count](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#a-connected-hypergraph-estimate).
For one block of \(s\) children, bound the sum over additional
vertices by
\(\beta(s-1)!\alpha^{s-1}(1-\alpha k)^{-s}\).
The bound follows from \(\binom kt\le k^t/t!\) and the
binomial series for the additional vertices.
Use \(t=\alpha/\beta\le1/(4e)\). At a parent with \(u\)
children the sum over partitions of the product of the block
bounds is at most
\[
 (\beta/(1-\alpha K))^u
       \sum_{\pi\in\mathfrak S_u}t^{u-\#\mathrm{cycles}(\pi)}
  =(\beta/(1-\alpha K))^u \prod_{v=0}^{u-1}(1+tv).
\]
The block factors \((s-1)!\) count cyclic orders in each block.
The product formula also follows by inserting the last label
in the cycle notation, either in a new singleton or at any
of the \(u-1\) positions in existing cycles. To pass from
unordered to ordered rooted labeled trees, divide this weight
at each parent by \(u!\), the number of orders of its children.

Let \(g_u=\prod_{v=0}^{u-1}(1+tv)/u!\). The total weight of
the ordered tree shapes with \(k\) vertices, where a vertex of
outdegree \(u\) has weight \(g_u\), is
\[
               \frac1k [y^{k-1}](1-t y)^{-k/t}.
\]
Indeed the preorder outdegrees \(u_1,\ldots,u_k\) sum to \(k-1\).
They encode a unique ordered tree precisely when the partial
sums of \(u_i-1\) are nonnegative before the final sum \(-1\),
by counting available child positions. Exactly one cyclic shift
of any such sequence gives these partial sums. Start just after
the first minimum among partial sums at indices \(0,\ldots,k-1\).
The earlier partial sums are at least one above that minimum,
so the wrapped prefixes are nonnegative as well. Any good shift
must start at that first minimum, proving uniqueness.
The total \(-1\) prevents a smaller cyclic period. Since
\(\sum_{u\ge0}g_u y^u=(1-ty)^{-1/t}\), the formula follows.

By positivity of the series, evaluate at \(y=1/(1+t)\) to
bound the tree weight by
\[
  \frac{(1+t)^{1/t}}{k}
               \bigl[(1+t)^{(1+t)/t}\bigr]^{k-1}
         \le 3^{k-1}\qquad(k\ge3).
\]
The prefactor is at most \(e/k\) and the base at most
\(e^{1+t}<3\) since \(t\le1/(4e)\). The same \(3^{k-1}\)
bound holds directly for \(k=1,2\).
Each shape has \((k-1)!\) assignments of the nonroot labels.
Combining with the extra-edge factor gives
\[
 |\Sigma_k(m)|\le (k-1)!\,\Xi^{k-1},\qquad
 \Xi=\frac{3\beta}{1-\alpha K}
                   \exp\!\left(\frac{\beta K}{2(1-\alpha K)}\right).
\]
The derivative satisfies the same majorant multiplied by
\(2k/n\). For each hypergraph the product derivative costs
\(2X/n\) times the product of the \(M_r\)'s, with \(X\) the
sum of its edge sizes. The inequality \(X\le k e^{X/k}\)
absorbs that extra factor by changing each \(M_r\) to \(eM_r\).
This proof also bounds terms whose weight is zero.

Consequently the truncated logarithm coefficients satisfy
\[
 |\ell_k(m)|\le B_m p_m^{k-1}/k,\qquad
 |\ell_k'(m)|\le4B_m p_m^{k-1}/n,\qquad p_m=\Xi B_m<13/25 .
\]
For the last bound, (2) gives \(\alpha K<297/1250\) and
\(\beta K/[2(1-\alpha K)]<17/10\) using \(e<68/25\).
Then
\(\Xi B_m\le24e\, (WD/n)(101/100)
       \exp(17/10)/(1-297/1250)<13/25\).
The case \(k=1\) follows directly from the first-factor bound.

Put \(A_m(w)=\sum_{k=1}^K\ell_k(m)w^k\). The geometric and
logarithmic series bound
\[
 \|A_m\|_{|w|\le9/5}\le11B_m/2,\qquad
 \|A_m\|_{|w|\le6/5}\le2B_m,\qquad
 \|A_N-A_n\|_{|w|\le6/5}\le13D/n .
\]
The first two use \(-\log(1-(9/5)(13/25))<143/50\)
and \(-\log(1-(6/5)(13/25))<26/25\), respectively.
At an integer order, \(H_m\) and \(e^{A_m}\) agree in all
coefficients up to and including degree \(K\).
Relative to \(|e^{A_m(w)}|\), the tail of \(e^{A_m}\)
on \(|w|\le6/5\) costs at most \(2(2/3)^K e^{15B_m/2}\).
The polynomial's own tail costs at most
\((3/20)^K e^{2B_m}\) there. Since \(B_m\le101D/100\),
\(K\ge32D+12\), and \(\log(3/2)>2/5\), the first bound is
at most \(2\exp(-24/5-(209/40)D)\le2\exp(-1361/240)\).
The second is below \(e^{-K}\). Their sum is less than
\(1/130\). The resulting analytic logarithm differs from
\(A_m\) by less than \(1/120\) on this disk.
Including the case of a private starting row, for both
orientations we have obtained
\[
 \|\log H_m\|_{|w|\le6/5}\le M:=3(D+1),\qquad
 \|\log H_N-\log H_n\|_{|w|\le6/5}\le\delta:=1/55 .         \tag{3}
\]

## Paired radius and phase comparison

The logarithms in (3) are real on the positive radius. Rouché
at \(|z|=\exp(7(D+1)/m)\), then for the reciprocal, confines
the root log radii \(q_m\) to \(|q_m|\le7(D+1)/m\).
The bound \(7(D+1)/n\le49D/n\) ensures the half arguments
have modulus at most \(21/20\) in this annulus.

For a common upper argument \(0<\theta\le\pi\), let \(J_m(q,\theta)\)
be the difference between the two half logs at
\(e^{q+i\theta}\) and \(e^{-q-i\theta}\). Cauchy's derivative
bounds from the disk \(6/5\) to \(21/20\) give
\[
 |J_N-J_n|\le2\delta,\qquad
 |\partial_\theta(J_N-J_n)|\le14\delta,\qquad
 |\partial_q J_m|\le14M,\qquad
 |\partial_q\partial_\theta J_m|\le210M .
\]
The real equations at \(q_N,q_n\) thus imply
\[
   (n-14M)|q_N-q_n|\le 7(D+1)/n+2\delta,\qquad
                           n|q_N-q_n|<1/25 .
\]
The imaginary part is zero at \(\theta=0\) when \(q\) is
real. Integrating the angular derivative of the order
difference and then using the mixed derivative for the
radius change gives
\[
 \left|\Im(J_N(q_N,\theta)-J_n(q_n,\theta))\right|
     \le (14\delta+210M|q_N-q_n|)\theta<\theta/3 .
\]
The difference of the two root equations again puts a
multiple of \(2\pi\) strictly between zero and \(2\pi\).
This concludes the proof. Integer coefficients give the
Bézout identity over \(\mathbb Q\), and hence over every
characteristic-zero field.

## Exact arithmetic replay

The [rational checker](verify_multivariate_absolute_curvature_27.py)
encloses the exponential constants with directed Taylor sums and
a positive geometric remainder. It also tests the weighted-tree
formula against the ordered-tree recurrence at three rational
values of the weight parameter. The cycle argument above gives
the identity at every tree size.
