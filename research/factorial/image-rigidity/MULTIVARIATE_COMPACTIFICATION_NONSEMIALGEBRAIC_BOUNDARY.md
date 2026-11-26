# A nonsemialgebraic boundary of the exact factorial compactification

Fix adjacent degrees \(n,N=n+1\) with \(n\ge2\), and restrict
the limiting ordinary curvature to
\[
                     \frac13\le S\le\frac25 .
\]
Consider the real compactified coefficient triples of the
[ordinary fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#exact-compactified-data).
Here the bound is on the liminf and limsup of the curvatures in
the approximating sequence. The exact-data proof in the cited
compactification applies at \(n=2\) as well. Its tail estimates
(3) and (5)–(7), the diagonal row extraction, and the prescribed
gap adjustment only need positive orders, and none uses
\(n\ge3\). The result below holds if the signed total gap is fixed
at any one integer, or is allowed to range over a nonempty finite
set.

**Theorem.** The exact compactified coefficient family is not
semialgebraic, even after projection to the two outer geometric
means. The same holds for the real coefficient pairs expressed
in a common physical variable. An algebraic curve cuts the
projection in infinitely many isolated points accumulating at a
strictly positive curvature profile. This gives an obstruction
to describing the entire exact compactification by one finite
system over the real closed field. Tests using semialgebraic
outer enclosures, or tests restricted to the collision locus,
are separate questions.

## The matched upper boundary

For the individually endpoint-normalized moment of order \(j\)
let \(F_{j,k}\) be the coefficient divided by \(\binom jk\).
Write
\[
  U=-\tfrac12\log(F_{n,1}F_{n,n-1}),\qquad
  W=-\tfrac12\log(F_{N,1}F_{N,N-1}).
\]
Use also \(q=e^{2U}\) and \(p=e^{2W}\). These are rational
observables of the real coefficient arrays. In fact a moment
with positive endpoints \(c_{j,0},c_{j,j}\) and coefficients
\(c_{j,k}\), before normalization, has
\[
         (F_{j,1}F_{j,j-1})^{-1}
              =\frac{j^2c_{j,0}c_{j,j}}{c_{j,1}c_{j,j-1}} .
\]
This makes both observables invariant under scalar and positive
argument normalization. In particular the relative center of
the two physical moments does not affect the projection.

For an ordinary row of smaller exponent \(\ell\) and gap \(d\),
a Gauss slot \(h\) has lower loss
\[
          s=\tfrac12\log\bigl(1+(n-1)d/(n\ell+h)\bigr),
                           \qquad 1\le h\le d .
\]
Its upper loss is \(\operatorname{artanh} T'\) with
\(T=\tanh s\). The
[exact adjacent tanh identity](MULTIVARIATE_MATCHED_TANH_MOMENTS.md#exact-tanh-transport-and-moment-errors-at-every-order)
gives \(T'\le n^2T/(n^2-1-T)\), with equality exactly at
\(h=d\). On \(0\le s<\log n\) the bound is equivalent to
\[
  \operatorname{artanh} T'\le
  g(s),\qquad g(s)=\tfrac12\log
                   \frac{(n^2-1)e^{2s}}{n^2-e^{2s}} .
\]
Here \(g(0)=0\), \(g'(0)=n^2/(n^2-1)\), and \(g\)
is strictly convex since
\(g'(s)=n^2/(n^2-e^{2s})\) is strictly increasing.

In the compactification the retained slots contribute these
same exact pairs of losses. The diffuse tail of variance
\(\tau\) contributes \(\tau(j-1)/(2j)\) to the order-\(j\)
loss, by evaluating its multiplier
\(\exp[-\tau k(j-k)/(2j)]\) at the two outer indices.
Thus if the diffuse lower loss is \(U_{\rm diff}\), its upper
loss is \(g'(0)U_{\rm diff}\). The other diffuse parameter
enters the relative center and disappears from \(q,p\).

## A sharp cap from the ordinary curvature band

A retained row has \(d^2/(\ell+d)\le2/5\).
If \(d=1\), this forces \(\ell\ge2\); if \(d\ge2\), then
\(\ell\ge(5/2)d^2-d\ge4d\).
In the second case \(d/(n\ell+h)<1/(2n+1)\).
It follows that all retained lower losses are bounded by
\[
 s_*=\tfrac12\log q_*,\qquad q_*=1+\frac{n-1}{2n+1}.
\]
The only slots attaining the bound are those of the unit-gap
row of unordered exponents \((2,3)\). Its upper exponential
loss is \(p_*=1+n/(2N+1)\). Also \(2s_*<\log n\),
since \(q_*<n\).

Consider the section \(s_*<U<2s_*\). Strict convexity permits
packing the lower losses under the cap \(s_*\), which gives
\[
                        W\le g(s_*)+g(U-s_*) .             \tag{1}
\]
For details, if two finite masses lie strictly between zero
and \(s_*\), transfer mass from the smaller to the larger until
one reaches a boundary. The sum of their \(g\)-values strictly
increases. At fixed total \(U\) in this section the packed
finite list is exactly \(s_*,U-s_*\). For countably many masses,
truncate the list and fill out the remaining mass by capped
pieces before passing to the limit. The inequality extends,
since \(g(s)\le s\,g(s_*)/s_*\) on the capped interval.
A diffuse contribution can likewise be split into finitely
many capped pieces. Replacing the term \(g'(0)s\) of any
positive such piece by \(g(s)\) strictly increases the bound.
Consequently equality in (1) excludes diffuse mass and
excludes two positive interior slots. It requires exactly the
two packed slots, with equality in the adjacent tanh bound on
both. The latter equality requires \(h=d\) at every retained
slot, which is possible only when each retained row has gap
one.

The equality profiles therefore have the hard row \((2,3)\)
and one soft unit row \((\ell,\ell+1)\), with arbitrary
orientations. Their total curvature is
\(1/3+1/(\ell+1)\), so \(\ell\) ranges exactly over the integers
\(\ell\ge14\). Every such two-row list is allowed by the
compactification, even at fixed signed total gap, by its
vanishing-curvature degree correction with \(\tau=0\).
Equivalently, for \(\ell>14\) the correction can also be made
inside the exact closed curvature band since its total
curvature tends to zero and the two-row curvature is interior.

## The algebraic section

In the variables \(q,p\), the range and equality in (1) are
\[
 q_*<q<q_*^2,\qquad
      \left(n^2-\frac{q}{q_*}\right)\frac{p}{p_*}
                    =(n^2-1)\frac{q}{q_*}.
\]
The projection of this section to \(q\) consists exactly of
\[
                  q_*\left(1+\frac{n-1}{n\ell+1}\right),
                         \qquad \ell=14,15,16,\ldots .
\]
These are distinct isolated points with limit \(q_*>1\).
A semialgebraic subset of the line is a finite union of points
and intervals, and projection preserves semialgebraicity by
real quantifier elimination.
[Perrucci and Roy, Theorem 1](https://arxiv.org/pdf/1609.02879)
Since both the section and the mapping from coefficient arrays
to \(q,p\) use rational operations, the displayed infinite
discrete projection proves the theorem.

The rational equation can also be verified without tanh
coordinates. For a unit-gap row its lower and upper exponentials
of twice the loss are
\(q_\ell=1+(n-1)/(n\ell+1)\) and
\(p_\ell=1+n/(N\ell+1)\). The following exact symbolic checks
verify the relation and its multiplication by the hard row.
~~~python
import sympy as sp

n, l = sp.symbols("n l", integer=True, positive=True)
N = n+1
q_l = 1+(n-1)/(n*l+1)
p_l = 1+n/(N*l+1)
q0, p0 = q_l.subs(l, 2), p_l.subs(l, 2)
assert sp.cancel((n*n-q_l)*p_l-(n*n-1)*q_l) == 0
q, p = q0*q_l, p0*p_l
assert sp.cancel((n*n-q/q0)*(p/p0)-(n*n-1)*(q/q0)) == 0
assert sp.cancel((sp.Rational(2,5)-sp.Rational(1,3)
                  -1/(l+1))-(l-14)/(15*(l+1))) == 0
~~~
