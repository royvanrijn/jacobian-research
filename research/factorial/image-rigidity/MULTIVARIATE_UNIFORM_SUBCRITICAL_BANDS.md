# Uniform eventual coprimality on shared subcritical bands

For distinct ordinary exponent columns \(a,b\), remove equal
rows and suppose every remaining row has \(a_i,b_i>0\).
Put \(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\),
\(W=\sum_i d_i^2/m_i\), and
\(\rho(q)=q(1-q)\sum_i d_i^2/(a_i+d_iq)\).

**Theorem.** Fix \(0<W_-\le W_+<\infty\) and \(0<\kappa<1\).
There is a terminal threshold depending only on these bounds
such that, for every fully shared ordinary profile satisfying
\[
 W_-\le W\le W_+,\qquad \max_{0<q<1}\rho(q)\le1-\kappa,
\]
all later consecutive moment pairs are coprime over every
characteristic-zero field. The total degrees need not agree.
The threshold is unevaluated.

The additional points beyond the
[equal-degree band theorem](MULTIVARIATE_UNIFORM_SUBCRITICAL_EQUAL_DEGREE.md)
are transversality of the active stationary arc to the
physical order drift, and a uniform Airy zero band near
its simple fold. The latter can handle unbounded integer
degree differences.

## The active arc and vertical fibers

Use the centered stationary target
\(w=-\operatorname{Log}(n^\Delta t)-C_0\) and the limit
phase of the
[bounded-band contour arguments](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md).
Here \(C_0=\sum_i(b_i\log b_i-a_i\log a_i-d_i)\) and
\(\Delta=\sum_i d_i\). At the same physical \(t\), the next
order uses \(w'=w-\Delta\log(1+1/n)\). The centered scales
of any proposed common-root sequence stay in a fixed
annulus by the uniform radius estimate.

Pass along a proposed counterexample sequence with
\(n\to\infty\) to a retained and diffuse phase limit.
Write \(\tau\ge0\) for its diffuse mass and \(G\) for its
chord-normalized coefficient phase. Its centered stationary
map is \(X=G'+\log(1-q)-\log q\).
The connected proper two-sheet cover, its unique simple
fold \(q_c\) over \(w_c\), and the finite contour homotopies
of the equal-degree band proof apply independently of
\(\Delta\): this degree difference enters only in
the change of centered parameter between adjacent orders.
These homotopies identify a locally fixed actual coefficient
at every regular saddle. At least one coefficient is
nonzero by the same winding argument on a real gap or
on the negative real ray in the purely diffuse limit.

Let \(D(w)\) be the difference between the two upper
stationary values with a local choice of ordering. The
equal-modulus arc is active where both actual saddle
coefficients are nonzero. As in the equal-degree proof,
activity forces individual zeros to approach every regular
point of that arc, by the real and imaginary exponential
comparison with \(D'(w)\ne0\), or its local Rouché form.
The bounded normalized scales, the positive-real
single-saddle estimate, the unique transverse real
negative boundary point, and the covering-disk
maximum principle therefore identify the only active
arc as the one from that negative boundary to \(w_c\).
These statements concern the limit of the single-order
polynomials, so no adjacent parameter drift is involved.

We prove that no vertical fiber occurs along this active arc.
For a putative vertical pair write \(q_j=x+i y_j\),
\(0<y_1<y_2\). The
[ordinary exterior-quadrant fiber bound](MULTIVARIATE_REGULAR_SADDLE_PHASE_SEPARATION.md#a-bound-on-fibers-in-the-exterior-quadrants)
passes to each retained and diffuse limit by local Rouché
on the preimages, so \(0\le x\le1\).
For the same reason an upper ramification point has real
part in \([0,1]\), since it is a double preimage.

The stationary value of the Rice phase, up to a
common function of the target, is
\(P(q)=G(q)-qG'(q)-\log(1-q)\).
Put \(Q=P+X\), and let \(p_i=-a_i/d_i\) be the
retained poles and \(\ell_i=|x-p_i|\).
Write
\[
 M(\ell)=\frac12\log\frac{\ell^2+y_2^2}{\ell^2+y_1^2}.
\]
Equality of the real parts of the two stationary values,
and equality of their targets, give
\[
 \sum_{\rm retained} a_i M(\ell_i)
       +\frac{\tau}{2}(y_2^2-y_1^2)=M(1-x),\qquad
 \sum_{\rm retained} b_i M(\ell_i)
       +\frac{\tau}{2}(y_2^2-y_1^2)=M(x).                 \tag{1}
\]
For a retained row, the difference of its real part
at the same \(x\) equals the corresponding log-modulus
difference, because the linear terms cancel.
The diffuse chord term is \(\tau q(q-1)/2\) and gives
the same displayed quadratic contribution to both
differences. The identities follow also by taking
finite retained truncations and passing to the limit.
All summands on the left are nonnegative. Since each
\(a_i,b_i\ge1\), (1) forces each
\(\ell_i\ge\max(x,1-x)\).

Set \(f(y)=\Re X'(x+iy)\). It has the form
\[
 f(y)=\tau+\sum_{\rm retained}
       \frac{|d_i|\ell_i}{\ell_i^2+y^2}
       -\frac{x}{x^2+y^2}-\frac{1-x}{(1-x)^2+y^2}.
 \tag{2}
\]
The ratio of its positive and negative parts as written
in (2) is nondecreasing as a function of \(y^2\).
For each retained term against each denominator term,
the differentiated ratio has a nonnegative contribution
because the retained squared distance is larger.
The constant \(\tau\) gives a strictly increasing ratio
when it is positive. Uniform convergence of the sums
on the vertical segment justifies differentiation.
The identities \(X(q_1)=X(q_2)\) and \(P'=-qX'\)
give \(\int f=0\) and \(\int yf=0\) between \(y_1\) and
\(y_2\). A nonzero function with at most one sign
change cannot satisfy both equations, by subtracting
the sign-change location times the first equation
from the second. The only remaining case is an
identically zero \(f\). For \(0<x<1\), the ratio in
(2) is strict: if \(x\ne1/2\) the two nonzero
denominator distances differ, and if \(x=1/2\)
all retained poles have distance strictly larger
than \(1/2\). At \(x=0\) identically zero requires
\(\tau=0\), all retained distances one, and
\(\sum |d_i|=1\). This is just one retained row
\((1,2)\), with no diffuse mass. The case \(x=1\)
is its reversal.

In the exceptional \((1,2)\) limit, the vertical fibers
have \(q=i y\) and target
\(\log(y+y^{-1})-i\pi/2\) up to a real centering constant.
Their equal-modulus component is precisely a horizontal
ray from the fold to infinity. It cannot be active:
activity is locally constant on the regular component,
whereas the centered annulus bounds all limiting root
scales. Thus on the active arc \(\Re D'\) never vanishes.
The height is a strictly monotone coordinate there.
The same is true through the regular lower endpoint,
where the two real branch values are strictly ordered.
In particular each height contains at most one point
of the closed active arc, including its fold endpoint.

## The elementary Airy class and its transverse band

At the limiting fold the actual local cubic relative class
is nonzero and has fixed integer weights on its three
descending sectors, by the
[uniform cubic transport](MULTIVARIATE_UNIFORM_SUBCRITICAL_EQUAL_DEGREE.md#uniformity-at-the-interior-fold).
In fact exactly one sector has weight zero. To see this,
use the canonical phase \(u^3/3-\lambda u\).
On the equal-modulus ray \(\lambda<0\), scale the two
saddles to \(i,-i\) for the phase \(z^3/3+z\).
The upper descending paths lie on
\(\Im(z^3/3+z)=2/3\), with the branch through \(i\)
given by
\[
 x=(y-1)\sqrt{\frac{y+2}{3y}},\qquad z=x+i y,\quad y>0 .
\]
Indeed substituting in the level equation gives
\(x^2=(y-1)^2(y+2)/(3y)\), and at \(i\) these are
the two directions of decreasing real phase.
Their ends are the upper right cubic ray and the
negative real cubic ray. Conjugation gives the lower
right ray and the same negative ray for the lower
saddle. There are no other stationary points on these
descending paths. Truncating inside a fixed chart
and connecting to the descended exterior chain
therefore gives the two true local Gaussian packets
at each sufficiently small fixed \(\lambda<0\).
Their integer coefficients are, up to orientation,
the endpoint weights on the upper right and lower
right rays. Rotation gives the analogous description
on each of the other two equal-modulus rays, with a
different pair of endpoint weights. These statements
persist for the finite-profile cubic charts by local
contour transport and a strict phase gap on the
truncated tails. If all three endpoint weights were
nonzero, all three regular arcs at the fold would
be active, contradicting the unique active arc.
At least two weights are nonzero because their
sum is zero and the class is nonzero. Hence the
class is one elementary Airy contour, with a nonzero
integer multiplier. Choose the rotation of the cubic
coordinate so that its active ray is \(\lambda<0\).

The active ray leaves the fold in a direction
that is not horizontal in the \(w\)-plane.
To verify this without assuming a uniform angle
bound across all profiles, write \(q_c=x+i y\).
Here \(0\le x\le1\) by the quadrant count.
At a critical point \(P'=-qX'=0\) and
\(Q'=(1-q)X'=0\). The \(y\)-derivatives of their
real parts therefore imply
\[
 \tau+\sum_i\frac{a_i}{\ell_i^2+y^2}
     =\frac1{(1-x)^2+y^2},\qquad
 \tau+\sum_i\frac{b_i}{\ell_i^2+y^2}
     =\frac1{x^2+y^2}.
\]
As in (1), every retained distance is at least
\(\max(x,1-x)\). The ratio in (2) is then strictly
increasing at its zero, except for the same two
one-row limits. In the strict case
\(\partial_y\Re X'(q_c)>0\), so \(X''(q_c)\)
is not real. In the analytic cubic chart,
\((\lambda'(w_c))^3=2/X''(q_c)\), by comparing the
cubic coefficient in \(q-q_c\) and the mixed phase
derivative. Consequently \(\lambda'(w_c)\)
is not real in any rotation. In the exceptional
one-row limit its horizontal equal-modulus ray
is inactive by the explicit vertical fiber formula.
The three rays at a simple fold have directions
separated by \(2\pi/3\), so neither remaining ray
is horizontal there. The active rotation again gives
\(\Im\lambda'(w_c)\ne0\).

On a small fixed disk around the fold, use the
finite-profile cubic normal form
\(C_j(w)+u^3/3-\lambda_j(w)u\), where \(j\)
indexes the profile. The actual class transports
uniformly and remains the same elementary Airy
class. After normalizing out \(e^{nC_j(w)}\),
the actual moment has the expansion with the
joint Airy error from the
[local cubic proof](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#a-uniform-joint-norm-for-the-actual-cubic-class).
In particular at any zero in this disk,
\[
 |\mathrm{Ai}(\xi)|
       \le C n^{-1/3}|\mathrm{Ai}'(\xi)|,\qquad
 \xi=n^{2/3}\lambda_j(w).
 \tag{3}
\]
Indeed divide the analytic amplitude modulo
\(u^2-\lambda_j\). Its constant coefficient is
bounded away from zero, its linear coefficient
is bounded, and its integrated remainder satisfies
the joint bound (2) of that cubic proof after one
integration by parts. These bounds are uniform
over the approximating profiles. The remainder can
be absorbed on the \(\mathrm{Ai}\) side of (3).
The exterior chains have the same uniform strict
loss, so do not change the estimate.

Put \(\gamma=-\mathrm{Ai}(\xi)/\mathrm{Ai}'(\xi)\).
The denominator is nonzero by (3) and uniqueness
for the Airy equation. Shrinking the physical disk makes
\(|\gamma^2\xi|\le C^2|\lambda_j|\le1/4\), and
\(|\gamma|^3\le1/16\) for large \(n\).
Continue this zero backward along the family
\(\mathrm{Ai}(\xi(t))+t\gamma\mathrm{Ai}'(\xi(t))=0\).
The Airy equation \(\mathrm{Ai}''(\xi)=\xi\mathrm{Ai}(\xi)\)
and differentiation give
\[
 \xi'(t)=-\frac{\gamma}{1-t^2\gamma^2\xi(t)} .
\]
While \(|\xi(t)-\xi(1)|\le2|\gamma|\), the denominator
has modulus at least \(1-1/4-2|\gamma|^3\ge5/8\).
The derivative bound \(8|\gamma|/5\) makes the
bootstrap strict and permits continuation to \(t=0\).
There \(\xi(0)\) is an Airy zero, which is real
and negative. Hence \(|\Im\xi|\le2|\gamma|\).
[Airy equation and zero locus](https://dlmf.nist.gov/9.9#i)
Thus every actual zero in the disk lies in the band
\[
                   |\Im\lambda_j(w)|\le C'/n .            \tag{4}
\]
The constants are uniform for the given sequence near
the limiting fold.

## The adjacent order shift

A common root is nonreal, since the even-order moment
is positive at every real parameter. By conjugation use
the upper half-plane. Pass to limits of \(w,w'\) in the
closed strip. Each lies on the closed active arc:
off the equal-modulus set a single nonzero Gaussian
is leading, and on an inactive equal-modulus arc there
is also a single actual Gaussian. The positive real
boundary has a one-saddle exclusion everywhere; at the
negative boundary the only possible limit is its unique
equal-modulus point. The fold is allowed in this initial
localization. Since the two parameters have the same
height, the strict height-coordinate result forces
their limits to agree. In particular \(\Delta/n\to0\).

If their common limit is regular, including the lower
endpoint, use the actual two-saddle expansions.
Their finite-profile stationary difference \(D_j\)
and its derivatives converge on a small disk to
the limiting \(D\). Equation (3) of the
[uniform supercritical comparison](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md#uniform-comparison-of-consecutive-orders)
has the same form here, with both nonzero leading
coefficients and uniform analytic bounds. In particular
a zero gives \(\Re D_j(w)=O(n^{-1})\).
If \(\Delta\ne0\), dividing the two equations gives
a right-hand ratio tending to one since \(w'-w\to0\).
But the log modulus on the left is
\[
 -\Delta\,\Re D_j'(w)
      +O\bigl(n^{-1}+|\Delta|/n+\Delta^2/n\bigr).
\]
The limiting real derivative is nonzero on the active
arc, so the leading term is bounded away from zero
for bounded nonzero \(\Delta\), and dominates the
error when \(|\Delta|\to\infty\). This is impossible.
If \(\Delta=0\), the adjacent quotient tends on the
left to \(e^{D(w)}\ne1\) at a regular active point.
At the lower endpoint \(D=i\pi\) gives the same
contradiction.

If instead both parameters tend to the fold, there
are two cases. For bounded \(\Delta\) the full physical
simple-fold exclusion of the
[adjacent Airy matrix](MULTIVARIATE_SHARED_SIMPLE_FOLD_AND_REGULAR_CONTOURS.md#comparing-both-orders-on-the-same-cubic-contour)
is uniform along the profiles. After passage to a
subsequence \(\Delta\) is a fixed integer. When it is
nonzero, the matrix determinant at the fold is the
nonzero analytic multiplier derivative \(\Delta\,dq/du\)
times nonvanishing bounded amplitudes. When it is zero,
use the bounded-field and large-field estimates of the
equal-degree band proof. Both use the common local Airy
class, finite contour loss, and uniform analytic Stirling
estimates already verified.

For unbounded \(|\Delta|\), the band (4) applies at
both \(w,w'\) in the same profile coordinate. By
analytic convergence, \(\Im\lambda_j'(w)\) stays away
from zero with one sign on a sufficiently small fold
disk. Along the real segment from \(w\) to
\(w'=w-\Delta\log(1+1/n)\), integration therefore gives
\[
 |\Im(\lambda_j(w')-\lambda_j(w))|
       \ge c\,|\Delta|\log(1+1/n).
\]
This contradicts the sum of the two \(O(n^{-1})\)
zero-band widths. All proposed unbounded-terminal
sequences are excluded. Since the coefficients are
integers the resulting Bézout identity transfers to
every characteristic-zero field.
