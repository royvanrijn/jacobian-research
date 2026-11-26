# All-degree coprimality for ordinary tails in a short exponent window

Fix any finite ordinary factorial core and let \(J\ge4\) be an
integer at least its profile-independent
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For the nonempty appended ordinary tail, write
\[
 L=\min_i\min(a_i,b_i),\qquad
 d=\min_i|a_i-b_i|,\qquad
 R=\max_i\max(a_i,b_i)-L
\]
after removing equal rows. Thus every tail exponent belongs to
\([L,L+R]\) and every tail gap is at least \(d\ge1\).

**Theorem.** For every fixed \(A>0\) there is a computable lower
bound on \(L\), depending only on \(A,J\) and the core, above
which every such tail with \(d\ge2\) gives ordinary consecutive
coprimality in every degree under the condition
\[
          R\le A d+\frac{d\log(2+d)}{512J^3}.
\]
There is no restriction on the number of rows or their directions.
Every exponent within the window can vary independently. In
particular the result covers arbitrary independent placements
within a window of width \(Ad\), with unbounded gaps.
In addition to the fixed core, one can retain arbitrary further
rows whose larger exponents are at most \(JR\).

A finite-template consequence is that the
[count-simplex certificate](MULTIVARIATE_AFFINE_TAIL_COUNT_POLYTOPE.md)
passes for every fixed common-slope affine template with gaps at
least two. In particular, all tails of the types
\((cM+a_i,cM+b_i)\), with one fixed positive integer \(c\) and
arbitrary fixed integer offsets with \(|a_i-b_i|\ge2\), satisfy every-degree
ordinary consecutive coprimality for all large enough \(M\),
uniformly over their independent multiplicities. The bounds in
both statements use the unevaluated cutoff \(J\).

## An envelope with one or two unmarked indices

For a low pair of consecutive orders with terminal \(4\le N<J\),
let \(m\) be the even order and \(m'\) the odd order. Strip the
factorial of the smaller endpoint from each appended row. At
order \(m\) the \(i\)-th stripped row coefficient at index \(k\)
is a product of consecutive integers from \(m\ell_i+1\) to
\(L_i(k)\), where \(\ell_i=\min(a_i,b_i)\) and
\[
 L_i(k)=
 \begin{cases}
 m\ell_i+k d_i,& a_i<b_i,\
 m\ell_i+(m-k)d_i,& a_i>b_i
 \end{cases}
 \qquad (d_i=|a_i-b_i|).
\]
Put \(E(k)=\max_i L_i(k)\) and \(\mu=\min_{0\le k\le m}E(k)\)
on the integer grid. Let \(t_3\) be its third-smallest ordinate,
counting repetitions.

Each direction's envelope is strictly monotone, with change at
least \(d\) per unit distance in its direction. When both
directions occur, the two envelopes cross at most once, and
\(E\) decreases to their crossing and increases after it. The
same description holds with a monotone \(E\) if there is just
one direction or the crossing is outside the index interval.
In the crossing case, only the two integers adjacent to the
crossing can have \(E(k)<\mu+d\). Indeed the envelope changes
by at least \(d\) from either of those integers to the next
one farther from the crossing on that side. If the crossing
is at an integer, only that index can have \(E(k)<\mu+d\).
In the monotone case the endpoint minimizer is the only such
index. Consequently \(t_3\ge\mu+d\).

Every integer \(x\) with \(\mu<x\le t_3\) exceeds all
\(m\ell_i\), since \(E(k)\ge\max_i m\ell_i\) at every
index. Therefore \(x\) is a literal factor at index \(k\)
in some appended row exactly when \(E(k)\ge x\). The set
\(I=\{k:E(k)<x\}\) contains one or two adjacent indices.
The cardinality bound follows from \(x\le t_3\) and adjacency
from convexity.

If \(d\ge2\), at least one mark in this interval is a
nonmultiple of \(m\), since the interval contains \(d\)
consecutive integers. For large \(d\) it contains many
marks in the single residue class \(1\bmod m\).

## Common-slope affine templates

For a fixed template of common slope \(c\), apply the envelope
argument to the offset rows of positive mass at any point of
the count simplex. Negative offsets are allowed in the
envelope computations. The selected nonmultiple mark has
the form \(x(M)=mcM+b\), with \(m\nmid b\). At order \(m'\),
every available literal factor has leading coefficient
\(m'c\). Proportionality to \(x\) requires its integral
intercept to be \(m'b/m\), impossible because
\(\gcd(m,m')=1\). At order \(m\) the limiting literal
ordinates are zero on \(I\) and positive outside it.
Their lower convex polygon has total horizontal width at
most one at slope zero. The other polygon only has slope
zero. Thus the fixed count-simplex test passes at each
point and by compactness its uniform two-match separation
is strictly positive. Its exact prime-depth transfer gives
the asserted computable height bound for all independent
multiplicities.

## Unbounded windows and unbounded gaps

It suffices to prove a uniform bound on \(d\) at a possible
collision satisfying
\[
                  R\le \frac{d\log(2+d)}{256J^3}.              \tag{1}
\]
Indeed the displayed theorem hypothesis implies (1) for all
sufficiently large \(d\), with a threshold computable from
\(A,J\). Once \(d\) is bounded, the theorem hypothesis also
bounds \(R\). A single fixed common-slope menu of all types
with offsets in that bounded integer window gives the height
bound by the previous section. If retaining further rows
with upper exponents at most \(JR\), the curvature bound
below restricts their number to less than \(12JR\). Each
unequal retained row contributes at least \(1/(JR)\).
There are therefore finitely many possible retained cores
in this bounded-window case.

Suppose there is a physical common root. The
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
gives \(\sum_i d_i^2/\max(a_i,b_i)<12\), with the sum
over the full profile. One tail row already gives
\(L+R>d^2/12\). For large \(d\), (1) implies \(R\le d^2/24\),
so \(L>d^2/24>R\). There are at least
\(q=\lfloor d/m\rfloor\ge d/(2J)\) integers in the residue
class \(1\bmod m\) in \([\mu+1,\mu+d]\).
Take any \(q\) consecutive terms of that arithmetic progression.
They all belong to \([mL,2mL]\). The lower bound follows
from \(\mu\ge mL\), and at an index outside \(I\) some row
attains at least \(x\) with every exponent at most \(L+R<2L\).
Write \(B=J^2R\).

At least one of these marks has a prime divisor \(p>B\) for
all sufficiently large \(d\). To prove this, suppose all
prime factors are at most \(B\). No prime dividing \(m\)
divides any selected mark. At any other prime \(p\le B\),
remove the contribution of a term of maximal \(p\)-valuation.
For every \(e\ge1\) reached by some term, the terms divisible
by \(p^e\) lie in one index residue class modulo \(p^e\),
so after the removal at most \(\lfloor(q-1)/p^e\rfloor\)
remain. The product of the remaining prime powers over all
primes divides \((q-1)!\). The factor removed for each
prime is at most \(2mL\). It follows that
\[
                  (mL)^q\le (2mL)^{\pi(B)}(q-1)! .
\]
The elementary bound
\(\pi(B)\le\sqrt B+8(\log 2)B/\log B\) follows by bounding the
sum of \(\log p\) for \(2^{s-1}<p\le2^s\) with
\(\log\binom{2^s}{2^{s-1}}\le 2^s\log 2\), then separating
primes at \(\sqrt B\). As \(R\ge d\), inequality (1) implies
\[
   \frac{B}{\log B}\le\frac{d}{128J}
         \quad (d\ge2)
\]
using \(\log(2+d)\le2\log d\).
Also \(\sqrt B=o(d/J)\) effectively for fixed \(J\). Thus
\(\pi(B)\le d/(8J)\le q/4\) for all sufficiently large
\(d\), by direct use of these bounds. Since \(q\le d\),
the product inequality gives \((mL)^{3/4}\le2^{1/4}d\).
This contradicts \(mL>d^2/6\) for \(d\ge32\).

For the mark divisible by this prime, all other literal
tail factors in the marked order differ from \(x\) by a
nonzero integer of magnitude at most \(mR<B<p\).
In the other order, a tail factor \(y\) gives
\(my-m'x\) of magnitude at most \(mm'R<B<p\). This integer
is nonzero modulo \(m\) because \(m\nmid x\) and the orders
are coprime. Therefore all tail factors at the other order
are \(p\)-units. For large enough \(d\), \(p>B\ge J^2 d\)
also exceeds the fixed core's full stripped factors and the
orders, so the binomial coefficients are units at both orders.
Every further retained row with upper exponent at most \(JR\)
has all stripped factors at most \(J^2R\), hence contributes
only units.

At order \(m\) the native marked polygon has zero ordinates
exactly on \(I\); the other polygon is entirely horizontal.
Their common slopes thus support a rational common factor
of degree at most one, by the
[Newton polygon factor-degree comparison](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks).
A rational common factor must have even degree because the
even ordinary moment is strictly positive on the real axis.
This excludes a physical collision when \(d\) is large.
All bounds above can be taken simultaneously at the finitely
many pairs with \(4\le N<J\). The first two pairs use
[ordinary low-order coprimality](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md),
and the absolute angular cutoff handles \(N\ge J\).
Rational coprimality extends over every characteristic-zero field.


## A balanced unit-gap obstruction to the single-mark count test

There are common-slope templates with unit gaps for which every
single proportional mark leaves at least two shared polygon
slopes at a signed-balanced point of the count simplex. For
\(r\ge2\), consider the pair of orders \(m=2r\) and \(m+1\)
and the four types
\[
     (M,M+r),\quad(M+1,M),\quad
     (M+K+r,M+K),\quad(M+K,M+K+1),\qquad K=r+1 .
\]
Give each type the same count. Their signed gaps sum to zero.
The two normalized offset windows are disjoint. All factors
proportional to any one mark belong to just one of the two
windows. Reversing the columns of its two types reverses the
ordinates and changes the signs of both slope multisets.
It therefore suffices to analyze the first window.

A mark in degree \(j\in\{m,m+1\}\) has the form \(jM+h\),
with normalized level \(t=h/j\in(0,r]\). If \(t\) is
nonintegral, at most one of the two degrees has a factor
proportional to the mark, since the consecutive degrees are
coprime. Discard the common positive count multiplier. At
the marked degree the increasing row contributes for
\(k\ge\lceil h/r\rceil\), while the decreasing row contributes
for \(k\le j-h\) if \(h\le j\). If \(0<t<1\), then
\(\lceil h/r\rceil\le2\) and \(j-h\ge1\). Thus every ordinate
is at least one and both endpoints equal one, so the whole
lower polygon is horizontal. If \(t>1\), the decreasing row
contributes nothing and at least the first three ordinates
are zero, giving horizontal width at least two. For
nonintegral \(t\) the unmarked degree only has slope zero,
so both cases leave at least two shared slopes. When
\(t\ge2\) is an integer, the same initial-zero argument
applies in both degrees.

At \(t=1\), the ordinates are
\[
 (1,0,1,\ldots,1)\quad\text{and}\quad(1,0,0,1,\ldots,1)
\]
in degrees \(m,m+1\) respectively. Both polygons contain
slope \(-1\) once and \(1/(m-1)\) with multiplicity \(m-1\).
This covers every proportional mark at the balanced point.
These calculations establish failure of the individual-slope
certificate at a balanced point.
