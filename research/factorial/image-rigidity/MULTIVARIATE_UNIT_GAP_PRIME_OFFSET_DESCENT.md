# Prime-offset descent for same-orientation unit-gap factorial binomials

Suppose every unequal row in an ordinary binary profile is
\((a_i,a_i+1)\), \(a_i\ge0\), after reversing both columns
if necessary. Let \(r\) be the number of such rows,
\(A=\max_i a_i\), and \(B=A+1\). Equal rows give only scalar factors.
Use the constant-normalized moment polynomial
\[
    F_m(z)=\sum_{k=0}^{m}\binom mk
                   \prod_i (m a_i+1)_k\,z^k .
\]
For an integer \(t\ge2\) and a prime \(p>2tB+2\), consider the
actual adjacent orders \(p+t-1,p+t\).

**Theorem.** Suppose \(r\le2\), all the \(a_i\) agree, or \(t>rB\).
If \(F_{t-1}\) and \(F_t\) have coprime reductions modulo \(p\),
then \(F_{p+t-1},F_{p+t}\) are coprime over every characteristic-zero
field. For \(r\ge2\), under the stated alternatives every common
root at the large orders would be a \(p\)-adic unit and would
reduce to a common root of the smaller pair.

This provides an arithmetic descent at arbitrary positive offsets,
including offsets growing with \(p\). For fixed \(t\) the condition
on the smaller reductions excludes at most finitely many primes
whenever that ordinary pair is coprime. For example, one or two
same-orientation unit-gap rows, or any number with the same base,
are covered at every \(p+2\) satisfying the prime bound outside the
prime divisors of the order-\(1,2\) resultant.
Writing \(C=\prod_i(a_i+1)\) and \(U=\prod_i(2a_i+1)\),
that resultant is \(C[C+(2^r-2)U]\), by evaluating
\(F_2=1+2Uz+2^rUCz^2\) at the root of \(F_1=1+Cz\).
In particular a single row has no exceptional prime in this
range. The argument
extends the [prime and post-prime comparisons](MULTIVARIATE_UNIT_GAP_PRIME_TERMINALS.md)
to a comparison with a lower-order pair.

## Exact Newton polygons from a fixed tail

For \(m=p+u\), \(u\in\{t-1,t\}\), the \(p\)-adic valuation of
the coefficient at \(k\) in \(F_m\) is exactly
\[
          y_u(k)=\mathbf1_{\{u<k<p\}}
                       +\#\{i:k\ge p-u a_i\}.                   \tag{1}
\]
In fact \(m<p^2\), and the binomial valuation counts one carry
precisely in the displayed interval. Each factorial increment
starts at \(p a_i+u a_i\), below the first following multiple
\((a_i+1)p\), and its sampled interval never reaches
\((a_i+2)p\) because \(u(a_i+1)<p\).
The first multiple has valuation one since \(a_i+1<p\).
This proves (1) with no assumptions on the relative bases.
The standard unit and carry formulas used here are the prime
case of [Granville's factorial congruences](https://www.cecm.sfu.ca/organics/papers/granville/paper/binomial/html/node2.html).

Read polygon slopes from the constant coefficient to the
leading coefficient. The horizontal face runs through
\(0\le k\le u\), and its unit polynomial is exactly
\(F_u\bmod p\) by Lucas' theorem and
\((p+u)a_i+h\equiv u a_i+h\pmod p\).
Thus a common root of valuation zero must reduce to a root in
both of the stated smaller reductions.

For \(r=1\) the rest of the lower polygon is one edge from
\((u,0)\) to \((p+u,1)\). This follows because all intermediate
ordinates after \(u\) are at least one. Its slope is \(1/p\)
and its width is \(p\). A rational common factor with no unit
roots must use width a multiple of the reduced denominator
\(p\), hence would have degree \(p\).
But the even-order moment is strictly positive at all real
values of \(z\), as the integral of an even power of a nonzero
real binomial. Every rational common factor has even degree.
Since \(p\) is odd this gives the desired contradiction.

For \(r\ge2\), the first positive slope joins \((u,0)\) to
\((p-uA-1,1)\). Before that endpoint all intermediate
ordinates after \(u\) are one. After the endpoint they are
at least two, whereas the extrapolated edge at \(k=m\) has
height \(p/(p-uB-1)<2\).
The two first positive slopes are therefore
\(1/(p-(t-1)B-1)\) and \(1/(p-tB-1)\), and are different.

All remaining vertices before \(p\) are right endpoints of
plateaus, necessarily among \(k=p-u a_i-1\). From \(p\) to \(m\)
the ordinate is \(r\), so the final point there is the only
possible vertex in that interval. Subtracting the first tail
vertex's horizontal coordinate gives a single finite template
for the remaining polygon:
\[
 \left(u(A-a_i),\,1+\#\{h:a_h>a_i\}\right)
       \quad\text{for the distinct bases }a_i,\qquad
                           (uB+1,r).                            \tag{2}
\]
Indeed all prior polygon edges are strictly below the tail
slopes under the prime bound. To check this directly, candidate
tail slopes with positive height have denominator at most \(uB+1\)
and integer numerator at least one, so exceed either first
positive slope. The other sample points of (1) within constant
plateaus cannot be lower vertices after the first positive edge.

## Comparing the two small tail polygons

For \(r=2\), or when all bases agree, the tail is just the line
from \((0,1)\) to \((uB+1,r)\). In the first case every
intermediate point has ordinate at least two, and in the
second case there are no intermediate points in (2).
Its slope is different at the two orders. Thus only the
horizontal slope remains possible.

For general \(r\ge2\), each tail slope has the form
\[
         H/(u c+\epsilon),\qquad 1\le H\le r,
\]
with \(1\le c\le B-1\) and \(\epsilon=0\) for an interior
edge, or \(1\le c\le B\) and \(\epsilon=1\) for the final edge.
A coincidence of slopes at \(u=t\) and \(u=t-1\) requires
\[
          t(H_1 c_2-H_2 c_1)
                    =H_1(c_2-\epsilon_2)+H_2\epsilon_1.        \tag{3}
\]
Its right side is nonnegative and at most \(rB\). If \(t>rB\),
equality is possible only when that side is zero. This would
require \(c_2=\epsilon_2=1\) and \(\epsilon_1=0\).
The smaller-order final edge would then start at \(k=p-1\),
the point for a zero base. Write \(r_0\) for the number of zero
bases. Its positive slope would be \((r_0-1)/t\).

To compare that possible last vertex, divide the horizontal
coordinates in (2) by \(u\). The interior locations and heights
are independent of \(u\); the final point moves strictly to
the left from \(B+1/(t-1)\) to \(B+1/t\) with no height change.
The point for the zero base is the last interior point.
If it is a strict lower vertex at the smaller order, it remains
one at the larger order: the slopes coming from the fixed
preceding points are unchanged, and its outgoing chord slope
to the endpoint strictly increases in the divided coordinates.
If it is collinear at the smaller order, use the entire last
edge, starting earlier with \(c_2>1\).
In the strict case the final slope in the original larger-order
coordinates is \((r_0-1)/(t+1)\), smaller than \((r_0-1)/t\).
Every larger-order tail slope is at most its final slope, so
even this exceptional coincidence is impossible.

This proves disjointness of the positive slopes in all the
claimed cases for \(r\ge2\). Horizontal reduction then gives
the stated descent. The Newton polygon root correspondence is
the one used in the [prime-terminal proof](MULTIVARIATE_UNIT_GAP_PRIME_TERMINALS.md).
Rational Bézout transfers coprimality to characteristic zero.

## A residual test for smaller offsets

Even when \(t\le rB\), the possible nonzero common slopes
for \(r\ge2\) are confined to the two small templates in (2).
Their residue polynomials can also be formed without any
factorial arguments proportional to \(p\).
For an index \(k=p-w\) on a tail edge, with \(k<p\) and
\(b_i=u a_i\), the unit part of its rising factorial in row
\(i\) is the reduction of
\[
 U_i(u,w)=
 \begin{cases}
  (-1)^{w-b_i}/((w-b_i-1)!\,b_i!),& w>b_i,\\[2pt]
  -(a_i+1)(b_i-w)!/b_i!,& w\le b_i .
 \end{cases}
\]
To see this, count whether the interval reaches \((a_i+1)p\).
Stopping \(w-b_i\) places before that multiple gives the first
formula by Wilson's theorem. Crossing it gives the unit
\(a_i+1\) and the unit factorial block on each side, which give
the second formula.
The unit of the binomial coefficient at this index is
\((-1)^{w+1}u!(w-1)!/(u+w)!\), by stripping its unique factor
\(p\) and applying Wilson again.
At the final index \(k=p+u\) the row units instead are
\(-(a_i+1)(b_i+u)!/b_i!\), with binomial factor one.

For each edge of common slope, use these rational coefficients
at the edge's sample points, translated by the initial index
and with powers in steps of the reduced slope denominator.
All denominators are units under the prime bound. Removing the
initial monomial gives precisely the ordinary edge residue
polynomial modulo \(p\), by substituting a suitable power of a
root of \(p\) and reducing the exposed edge.
The two polynomials at a shared positive slope must have
a common nonzero residue root if that slope occurs in a
common factor. This yields a bounded-size extra modular test
alongside the smaller horizontal comparison. For example, the
bases \(7,2,2\) and offset \(t=3\) do have a common tail slope
\(1/10\), so the horizontal comparison alone does not follow
from Newton geometry for all offsets.

In that example the two tail residues at slope \(1/10\)
are linear in the residual variable. Up to nonzero scalars,
they have roots
\[
          -\frac1{742711167339724800000},
             \qquad
          -\frac1{7085086867792382341939200000}.
\]
The determinant of the unscaled coefficient pairs is
\(1122293/6816132755613941760000000\), with numerator
\(31\cdot41\cdot883\). The smaller pair has resultant
\(2^{14}3^5\cdot733886302043\).
Thus every prime \(p>50\) not dividing \(883\cdot733886302043\)
certifies this entire prime-offset family, despite the shared
Newton slope.

## Terminals \(mp+t\) with a higher prime quotient

Suppose \(r\ge2\), let \(m\ge2\), and require
\[
        (mp+t)B<p^2,\qquad p>\max(2tB+2,r(tB+1)).
\]
For terminals \(N=mp+t\) the same horizontal and tail tests apply,
with one additional possible common slope. Its denominator is \(p\),
its width is \((m-1)p\), and its residue polynomial in both orders
is a nonzero scalar times
\[
        \sum_{v=0}^{m-1}\binom{m-1}{v}
                  \prod_i(m a_i+1)_v\,((-1)^r z)^v .           \tag{4}
\]
Consequently, for \(m=2\) the descent theorem holds at \(N=2p+t\)
under the displayed prime bounds. For larger \(m\), one may also
exclude the new slope if (4) is irreducible of odd degree over
\(\mathbb F_p\). In general (4) gives the possible additional
factor widths by multiplying its divisor degrees by \(p\).

Indeed, at order \(mp+u\) write \(k=pv+w\) with \(0\le w<p\).
The valuation at \(k\) is \(rv+y_u(w)\), by the same first-level
carry calculation as (1). No second-level carry occurs since
all factorial arguments are less than \(p^2\).
Put \(q=p-uA-1\). The centered ordinate
\(y_u(w)-rw/p\) has its unique minimum at \(q\).
On \(0\le w\le u\) its lowest value is \(-ru/p\), greater than
\(1-rq/p\) since \(q-u>p/2\) and \(r\ge2\).
On \(u<w\le q\) its minimum is at \(q\) itself.
After \(q\) the ordinate is at least \(2-r\), strictly larger
because \(p>r(uA+1)\).
The polygon is therefore the old initial horizontal and first
positive edges, followed by the slope \(r/p\) between \(q\) and
\(q+(m-1)p\), followed by the old tail shifted to the last
block. The first positive slope is smaller than \(r/p\);
the tail slopes are larger by \(p>r(uB+1)\).
The first and last block assertions now follow directly from
(1) and (2). Thus the comparison of all slopes except \(r/p\)
is unchanged. The horizontal residues are still \(F_u\bmod p\).

On the central edge, let \(k=vp+q\).
At such an index the unit of the binomial coefficient divided
by its unit at \(k=q\) is \(\binom{m-1}{v}\).
This follows by stripping its one factor \(p\): in the
top-digit factorials the quotient is
\(m!/(v!(m-v-1)!)\), while all low-digit factorials and
the Wilson sign are independent of \(v\).
For row \(i\), adding \(vp\) to \(k=q\) raises its valuation by
\(v\), for it reaches the consecutive multiples
\((m a_i+1)p,\ldots,(m a_i+v)p\). These all have unit factors
\(m a_i+1,\ldots,m a_i+v\). The intervening complete residue
blocks add a factor \((-1)^v\).
Taking the product gives the relative unit in (4).
Since \(p>r\), the reduced slope denominator is \(p\), which
explains the displayed residue degrees and widths.
When all other slopes are excluded by the earlier tests,
a common factor must use only this central slope.
For \(m=2\) its only nonzero width is \(p\); for an irreducible
odd-degree residue its only nonzero width is \(p(m-1)\).
Both are odd, contrary to even-order positivity.
