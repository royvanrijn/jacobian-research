# A compact-facet local comparison for Gaussian agreement kernels

Let \(T\) be the centered regular tripod in the Gaussian plane
and fix \(0<\rho<1\). Let \(I\) be a compact subinterval of the
interior of one of its facet rays. In tangential and normal
coordinates put \(E_h=\{(s,t):s\in I,\ |t|<h\}\). Consider a
bounded measurable kernel \(F\) satisfying the range, symmetry,
balance, and triangle conditions
\[
 0\le F\le1,\qquad
 \int F(x,y)\,d\gamma_2(y)=1/3,\qquad
 F(x,y)+F(y,z)-F(x,z)\le1 .
\]
Assume also that its integral operator on \(L^2(\gamma_2)\)
is positive semidefinite and that
\[
       F(x,y)=1_{\{T(x)=T(y)\}}
              \qquad\text{for almost every }(x,y)\in(E_h^c)^2.
\]
Then for all sufficiently small \(h>0\), there is a strict
comparison
\[
      \iint M_\rho F\,d\gamma_2^2 \le S_\rho(T)-cL ,
          \qquad c>0,                                      \tag{LC1}
\]
where \(L\) is the simplex score regret of the soft partition
determined by the unchanged label regions below. Equality of
the objectives holds only for the genuine agreement kernel
of \(T\). The result allows transfer to the third label inside
the strip. The constants can be chosen uniformly for \(\rho\)
on a compact positive-correlation band and a fixed \(I\).

The assumptions fix the agreement on pairs of unchanged points
before rotation averaging. The proof controls arbitrary additional
positive-semidefinite feature directions within the strip by
triangles with the unchanged labels.
[The agreement triangle relaxation](GAUSSIAN_AGREEMENT_TRIANGLE_RELAXATION.md)

## Features and triangles with unchanged labels

There is a measurable real Hilbert-space feature map \(V\)
with \(F(x,y)=\langle V(x),V(y)\rangle\) almost everywhere
and \(\|V(x)\|\le1\) almost everywhere. To see the norm bound,
let \(P_n\) be the conditional projections for a refining
sequence of finite partitions. The trace of the positive
compression \(P_n F P_n\) is the sum over atoms \(B\) of
\(\gamma_2(B)^{-1}\int_{B^2}F\,d\gamma_2^2\), which is at most
one. Hence the positive kernel operator has trace at most one
and its spectral feature map is square integrable. The same
atom calculation shows that the feature mean on each atom has
norm at most one. Conditional expectation convergence gives
the almost-everywhere bound.

On each unchanged label region, pairs have agreement one.
The feature mean on that region has norm one, so all its
features coincide almost everywhere with that mean \(e_i\).
Cross-region agreement zero makes the three \(e_i\) orthonormal.
The total feature mean \(m=\int V\,d\gamma_2\) satisfies
\(\langle m,e_i\rangle=1/3\) by the row constraints, and
\(\|m\|^2=1/3\) by their integral. Thus \(m=\sum_i e_i/3\).
Put \(p_i(x)=\langle V(x),e_i\rangle\). Pairing with the
unchanged regions shows \(0\le p_i\le1\), and the mean identity
gives \(\sum_i p_i=1\) and \(\int p_i\,d\gamma_2=1/3\).
This is the induced soft partition.

Put \(w(x)=V(x)-e_{T(x)}\) and \(q(x)=1-p_{T(x)}(x)\).
The vector \(w\) vanishes outside the strip and has integral zero.
The simplex scores are \(U_i=T_\rho1_{T_i}\), and
\[
 L=\int \sum_i p_i(x)\big(U_{T(x)}(x)-U_i(x)\big)\,d\gamma_2(x).
\]
Expanding the feature interaction gives the exact objective
difference
\[
   \iint M_\rho F\,d\gamma_2^2-S_\rho(T)
        = -2L+\iint M_\rho(x,y)
                    \langle w(x),w(y)\rangle\,d\gamma_2(x)d\gamma_2(y).
                                                               \tag{LC2}
\]

For any unchanged label \(i\), the triangle rows with a
point in that label give
\[
 p_i(x)+p_i(y)-1 \ \le\ F(x,y)
             \ \le\ 1-|p_i(x)-p_i(y)|.                       \tag{LC3}
\]
These rows hold almost everywhere by integrating over its
positive-measure unchanged region. Suppose the two reference
labels on the strip are \(i,j\), and let \(k\) be the spectator.
Write \(z(x)=p_k(x)\). For \(x,y\) on the same reference side,
(LC3) gives
\[
     0\le\langle w(x),w(y)\rangle\le2\min(q(x),q(y)).
\]
For opposite sides, say \(T(x)=i,T(y)=j\), the inner product is
\(F(x,y)-q(x)-q(y)+z(x)+z(y)\). The two lower bounds in (LC3)
for \(i,j\) give \(\langle w(x),w(y)\rangle\ge-2\min(q(x),q(y))\).
Their upper bounds give
\(\langle w(x),w(y)\rangle\le\min(z(x),z(y))\).
In particular the absolute inner product is at most
\(2\min(q(x),q(y))\) on all pairs in the strip.

## The normal flux and the one-pair spectral gap

Write \(d\nu(s)=\phi_2(s,0)\,ds\) for Gaussian arclength on the
reference ray and \(\omega(t)=e^{-t^2/2}\), so that
\(d\gamma_2=d\nu(s)\omega(t)\,dt\) in these coordinates.
Let \(H(s)=\int w(s,t)\omega(t)\,dt\), with \(H=0\) off \(I\).
It has zero \(\nu\)-mean. The preceding bounds yield, for
almost every \(s\),
\[
 \|H(s)\|^2\le
      4\int |t|q(s,t)\omega(t)\,dt
      +h\int z(s,t)\omega(t)\,dt.                           \tag{LC4}
\]
For the same-side pair integral use the layer-cake formula
for \(\min(q(s,t),q(s,t'))\). On either half-line, a level
set of weighted mass \(u\) has weighted first moment at least
\(u^2/2\), since its density is bounded by one and filling the
interval \([0,u]\) minimizes the moment. The factor \(2\) in
the same-side inner-product bound gives the first term in
(LC4). For the two opposite-side integrals use
\(2\min(z(s,t),z(s,t'))\le z(s,t)+z(s,t')\) and the fact that
each half-strip has transverse weighted length at most \(h\).

Put \(a(s)=|\partial_t(U_i-U_j)(s,0)|\). It is bounded above
and separated from zero on \(I\). The score difference of
the two incident labels has size
\(a(s)|t|+O(t^2)\), and the spectator score gap has a uniform
positive lower bound on the compact strip. It follows for
small \(h\) and constants \(C,\delta>0\) that
\[
 L\ge (1-Ch)\int_{E_h}a(s)|t|q(s,t)\,d\gamma_2
                       +\delta\int_{E_h} z(s,t)\,d\gamma_2.
\]
Indeed the incident-label transfer has mass \(q-z\) and the
spectator transfer has mass \(z\). The incident cost that was
substituted for \(z\) is \(O(h)z\), absorbed by its positive
score gap. Consequently (LC4) gives
\[
             \tfrac12\int a(s)\|H(s)\|^2\,d\nu(s)
                         \le (1+Ch)\,2L .                  \tag{LC5}
\]
[Reference score gradient on the regular facet](POSITIVE_CORRELATION_CONTINUATION.md#the-full-interface-form-at-a-fixed-positive-correlation)

The full interaction in (LC2) equals the projected Mehler
interaction of \(H\) up to \(O(h)L\). On the compact strip the
two kernels differ by \(O(h)\) by smoothness. For the absolute
error, apply the bound \(2\min(q(x),q(y))\) and layer cake.
If \(B\subset E_h\) has Gaussian mass \(u\), its first
absolute normal moment is at least \(u^2/(4\nu(I))\), by
filling the centered normal interval with density \(\nu(I)\).
The transverse projection of \(1_B\,d\gamma_2\) onto the normal
coordinate has density at most \(\nu(I)\). Apply this bound
to each level set of \(q\) and integrate. It gives
\(\int_{E_h^2}|\langle w(x),w(y)\rangle|\,d\gamma_2^2
\le 8\nu(I)\int |t|q\,d\gamma_2=O(L)\), as claimed.

For a scalar balanced speed \(f\) on one simplex ray, the
planar gap (S37) gives
\[
 \int a f^2\,d\nu-2\iint M_\rho(s\tau,s'\tau)
                  f(s)f(s')\,d\nu(s)d\nu(s')
        \ge \tfrac23 \rho\kappa_3(\rho)\|f\|_{L^2(\nu)}^2 .
\]
The factor two comes from the squared norm of the one-pair
label incidence, and the displayed form is \(\rho\mathfrak q_\rho\)
in the normalization of (B6). The projection of \(f\) onto the single planar
rotation speed has squared norm at most \(\|f\|_2^2/3\),
because the rotation has identical squared norm on each of
the three rays. Thus the squared distance in (S37) bounds the
remaining two-thirds. Apply this scalar inequality to each
coordinate of \(H/\sqrt2\). The zero Gaussian arclength
mean of \(H\) gives exact balance in each coordinate. Hence
for some \(\xi>0\) the projected interaction is at most
\((1-\xi)\int a\|H\|^2/2\), using the uniform upper bound
on \(a\). The coordinate sum converges by boundedness of the
Mehler interaction on this facet. Inserting (LC5) and the
\(O(h)L\) projection error into (LC2) proves (LC1) for
sufficiently small \(h\). If \(L=0\), the positive score
gaps off the reference facet imply \(p_i=1_{T_i}\) almost
everywhere. The feature norm bound then forces \(V=e_{T(x)}\)
almost everywhere.
[The explicit planar simplex spectral gap, (S37)](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)

The argument uses the one-pair spectral gap. With several
facets, the feature fluxes can carry directions orthogonal
to all three label anchors, and their cross interactions
need not use the ordinary signed label incidence. An exact
several-facet or junction comparison for the augmented
relaxation remains open.
