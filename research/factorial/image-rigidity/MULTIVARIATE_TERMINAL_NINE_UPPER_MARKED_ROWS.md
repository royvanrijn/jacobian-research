# Terminal-nine angular separation with a marked unit row

Let the unequal rows of two distinct ordinary exponent columns
be fully shared. Let \(u\) be their order-eight outer geometric mean
after positive endpoint normalization and removal of binomial factors.
The order-eight and order-nine polynomials have disjoint root
arguments in either of the following cases.

- The profile contains a row with entries \(3,4\) and \(u\ge21/25\).
- The profile contains a row with entries \(4,5\) and \(u\ge22/25\).

The rows may have either orientation. The bounds allow all
remaining shared exponents and dimensions. The proof also
gives the respective sufficient residual conditions
\(v\ge19/20\) and \(v\ge121/125\) after one marked row is removed.

## Residual profiles and their coefficient enclosures

An order-eight unit row with smaller exponent \(a\) contributes
\(u_a=\sqrt{(8a+1)/(8a+8)}\) to the outer geometric mean.
The full mean is the product over the rows and the residual mean
is \(v=u/u_a\). The thresholds follow by squaring:
\[
 (21/25)^2> (25/32)(19/20)^2,\qquad
 (22/25)^2> (33/40)(121/125)^2.
\]
The specified lower bound on \(v\) in either case excludes another
copy of the same row. Swapping the two columns for both orders
sends each root to its reciprocal up to positive rescaling, so fix
the orientation in which the marked exponent increases.

The exact normalized marked coefficients at order \(j=8,9\) are
\[
   C_{j,k}=\binom jk\frac{(ja+k)!}
     {((ja)!)^{1-k/j}((j(a+1))!)^{k/j}} .
\]
Any residual order-eight reciprocal is at most
\(c=(v^{-2}-1)/7\), because \(v^{-2}\) is the product of the factors
\(1+7x\). Its matched order-nine reciprocal is at most \(8c/(9-c)\).
For nonflat residuals with adjacent outer means \(v,w\), use the
[linked logarithm ratio](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections)
\[\rho=\log w/\log v\]
and the same grouped coefficient enclosure.
At a residual interval midpoint \(v_0\) the reference arrays are
\(C_{8,k}v_0^{e_{8,k}}\) and \(C_{9,k}v_0^{(64/63)e_{9,k}}\),
where \(e_{j,k}=k(j-k)/(j-1)\). To these reference coefficients the
checker applies the three skew vectors, the two mean defect vectors
and the linked outer vector of that enclosure. It keeps each linear
term and quadratic product together as one interval scalar times
its full coefficient vector. The remaining per-coefficient errors
include both the truncated logarithm series and the cubic
exponential remainder. This bounds the flat residual \(v=1\) as
well by continuity on every retained closed interval ending at one.

The factorial roots in \(C_{j,k}\) are bracketed by integer powers.
Each reference power \(v_0^{p/q}\) is bracketed on the dyadic grid
\(2^{-52}\) by testing integer powers against \(v_0^p\). The
floating-point proposal affects only speed. This finer bracket
keeps the cancellations in the central Cayley coefficients visible
to the interval arithmetic.

## Root disks and angular tests

Use \(p=1/8\), \(\lambda=\sqrt p\), and the Cayley coordinate
\(t=(\lambda X-1)/(\lambda X+1)\).
The transformed polynomial is the sum of the normalized
physical coefficients times \((\lambda X-1)^k(\lambda X+1)^{j-k}\).
Apply this linear map to each grouped contribution separately.

For each rational residual subinterval the checker selects trial
centers from the fixed affine tables and chooses a rational radius
for each center. The interval Taylor coefficients \(D_h\) around
the center include its entire grouped model. If
\[
  |D_0|+\sum_{h=2}^j |D_h|r^h < |D_1|r
\]
holds with certified bounds on the two sides, Rouché's theorem
gives one root in the disk. The disks and their conjugates are
disjoint. They account for all \(j\) roots and fit within
\(|X|<2\), which maps into the physical open left half-plane.
Four disks for each order are strictly above the real axis.
The additional order-nine disk has real center and its single
root is real by conjugacy.

The seven levels \(b\), selected by affine interpolation in
the residual mean, mark physical arguments
\(\pi-2\arctan(\lambda b)\). For a point \(X=x+iy\) in an
upper disk the sign of
\[
                 (y-b)(1+pby)+pbx^2
\]
locates its side of the ray. Its derivative in \(b\) is
\(-1+p(x^2+y^2)-2pby<0\) throughout the upper disk.
The checker proves strict signs on the whole bounding rectangle
of each disk. It places the four nonic disks in alternating upper
sectors and the four octic disks in the four intervening sectors.
Conjugation and the real nonic root complete the separation.

The [exact checker](verify_multivariate_terminal_nine_upper_marked_rows.py)
certifies 36 octic and 48 nonic residual intervals for \((3,4)\),
and 46 at each order for \((4,5)\). Run it from the repository root.

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_upper_marked_rows.py
~~~
