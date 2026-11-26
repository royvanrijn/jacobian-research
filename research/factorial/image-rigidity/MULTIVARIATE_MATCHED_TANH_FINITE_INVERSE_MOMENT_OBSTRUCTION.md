# A weighted angular collision with ordinary inverse-tanh moments of any prescribed finite degree

Fix an integer \(m\ge0\). There is a positive real weighting
of mirrored ordinary rows, all of gap \(30\), with the
following properties.

- Its exact factorial coefficient arrays at orders eight and
  nine have a common nonreal unit-circle root. Both
  polynomials have simple unit-circle roots.
- In each orientation, its row sums of every polynomial of
  degree at most \(m\) in the smaller exponent agree exactly
  with those of an unweighted ordinary profile of the same
  gap. The comparator works simultaneously for every
  polynomial in that space.
- In particular, for every moment order, every Gauss shift
  index, and every collection of matched orders, its
  inverse-tanh moments of total degree at most \(m\) agree
  with those of the same ordinary comparator. Both unsigned
  and signed versions have this property.

The equality of aggregate inverse moments holds for any
prescribed finite degree. The exact joint integer multiplicity
criterion still excludes the weighted profile.

## Matching the aggregate data

Use the [weighted matched-row collision](MULTIVARIATE_WEIGHTED_TANH_ANGLE_OBSTRUCTION.md)
with gap \(d=30\), smaller exponent \(1\), and weight
\(\mu=\kappa/30\) on each of its two mirrored rows. The parameter
\(\kappa\) ranges over the same compact interval
\(I=[1021/2000,2553/5000]\).

Choose an integer \(L\ge2\), and for \(1\le r\le m+1\)
let \(\beta_r(L)\) be the Lagrange interpolation weights
for evaluating at \(1\) from the nodes \(L,2L,\ldots,(m+1)L\).
Thus
\[
 \beta_r(L)=\prod_{h\ne r}\frac{1/L-h}{r-h},\qquad
        f(1)=\sum_{r=1}^{m+1}\beta_r(L)f(rL)
                         \quad (\deg f\le m).
\]
For \(m=0\) the product is empty and the identity evaluates
constant functions. Since \(0<1/L<1\),
\[
             |\beta_r(L)|
             \le\prod_{h\ne r}\frac{h}{|r-h|}
             =\binom{m+1}{r}.
\]
At each soft node \(\ell=rL\), append a mirrored pair of gap
\(30\) with the same weight
\[
 v_r(\kappa,L)=W_r-\mu\beta_r(L),\qquad
                   W_r=1+\binom{m+1}{r}
\]
on each row. These are positive and uniformly bounded for
\(\kappa\in I\), since \(0<\mu<1\).
Let the ordinary comparator have \(W_r\) mirrored pairs at
each soft node and no hard row. At either orientation,
\[
              \mu f(1)+\sum_r v_r(\kappa,L)f(rL)
                              =\sum_r W_r f(rL)        \tag{1}
\]
for every polynomial \(f\) of degree at most \(m\).
The identity holds throughout \(I\), not just at the
parameter where the adjusted factorial polynomials collide.

For a fixed Gauss shift \(h\in\{1,\ldots,d\}\) and order
\(j\ge2\), the inverse tanh factor is
\[
           T_{j,h}(\ell)^{-1}
                     =1+\frac{2j\ell+2h}{(j-1)d}.
\]
Thus every polynomial of total degree at most \(m\) in
any finite list of these inverse factors is a polynomial
of degree at most \(m\) in \(\ell\). This remains true if
the factors have different shift indices. Apply (1) before
summing over shifts or orientations. This proves all the
claimed inverse-moment equalities. For \(m\ge1\) they include
integer total row counts, Gauss factor counts, total lower
exponents, and the
[marginal inverse count constraints](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#constraints-from-just-the-inverse-moments).

## Persistence of the factorial collision

For smaller exponent \(\ell>0\), the endpoint-normalized
binomial-stripped coefficient of one mirrored pair at
order \(j\) is
\[
 U_{j,k}(\ell)=\frac{(j\ell+dk)!(j\ell+d(j-k))!}
                     {(j\ell)!(j\ell+dj)!}.
\]
The two zero-endpoint Gamma chord errors are nonpositive.
The trigamma sum bound
\(\psi'(1+x)\le1/x\) for positive \(x\) gives, by integrating
each second derivative against the endpoint Green kernel,
\[
           0\le-\log U_{j,k}(\ell)
                           \le\frac{d^2 k(j-k)}{j\ell}.
\]
For clarity, the bound on the second derivative at the
fractional index \(s=k/j\) of each log-Gamma term is
\(jd^2/\ell\). Each chord error is at most this bound
times \(s(1-s)/2\). The two terms add to the displayed bound.

Consequently the additional multiplier
\(\prod_r U_{j,k}(rL)^{v_r(\kappa,L)}\) tends uniformly to
one on \(I\) as \(L\to\infty\), for the two orders \(j=8,9\).
Endpoint coefficients remain one and all coefficients stay
positive. The two monic reciprocal quartics from the original
collision example therefore converge uniformly in coefficient
space to its quartics. Their resultants vary continuously and
retain the opposite strict signs at the two endpoints of \(I\)
once \(L\) is large enough. In addition, all the strict signs
in the rational root isolation intervals of that example
persist throughout \(I\). Indeed each of those finitely many
signed evaluations has a positive lower bound on the compact
interval. Thus the perturbed quartics have four simple roots
in \((-2,2)\), and only the fourth isolation intervals overlap.
The resultant has a zero in the interior of \(I\), and its
common quartic root gives the required common conjugate
unit-circle pair. The odd moment also retains its simple
structural root \(-1\).

These choices are effective for each \(m\). The log error
above has an explicit upper bound of order \(1/L\) with a
constant depending on \(m\). It can be compared with the
existing strict interval margins and with rational bounds
on the sensitivity of the finite Sylvester determinant.

## Scope of the obstruction

The homogeneous rowwise matched tanh transport remains valid
with these positive weights, using the actual per-row caps.
The hard row has \(u=\ell/d=1/30\) and no soft row has that
ratio. Its joint atom at \(a=h/d=1\) still has noninteger
weight \(\kappa/30\). This violates the exact
[local integer multiplicity reconstruction](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md)
even though all of (1) is enforced.

The earlier sharp unweighted packing restriction on the
positive tanh moments remains a separate constraint.
For sufficiently large \(L\) the strict packing violation in
the hard example persists, since the soft rows contribute
\(O(1/L)\) to every fixed positive tanh moment of order at
least one. The obstruction isolates a limitation of any
fixed finite collection of aggregate inverse moments and
leaves coprimality for ordinary rows open.
