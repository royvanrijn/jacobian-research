# Native prime separation in a widening factorial height band

Fix a finite ordinary two-monomial core \(C\). Append a nonempty
collection of unequal rows with smaller exponents in the integer
interval \([M,M+W]\) and gaps in \([1,D]\), with independently
chosen orientations and multiplicities. Here \(D\ge1\) is fixed,
whereas \(M\), \(W\), and the number of appended rows may vary.
The gap of a row is the absolute difference of its entries.

There is an explicit integer \(K(W,C,D,J)\) with the following
property, where \(J\ge16\) is any integer at least the
[effective absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).
If \(M>K\), every consecutive pair with terminal \(N\ge15\)
is coprime under any one of these conditions:

- \(N\) is even.
- The even order \(m\) in the pair exceeds \(2D\).
- Every appended gap is at least two, or every appended gap is one.

In particular, for \(D\le6\), or for either gap condition in the
third case, the original profile satisfies physical SFC(2) in
every degree. This uses the
[complete theorem through terminal fourteen](MULTIVARIATE_TERMINAL_FOURTEEN_ALL_PROFILES.md).
All statements hold over any characteristic-zero field. For fixed
\(C,D,J\) the height threshold satisfies
\[
                \log\log K=O_{C,D,J}((W+1)\log(W+2)).
                                                               \tag{1}
\]
Thus, for example, the result allows
\(W=o(\log\log M/\log\log\log M)\) as \(M\to\infty\).
It imposes no row-count cap and gives a uniform threshold for
a finite list of cores by taking the largest one. The cutoff
\(J\) itself is unevaluated.

The mark below uses the maximum of the factorial arguments at
each coefficient index. This choice handles arbitrary offsets
within the band, extending the setting of the
[common-upper-height theorem](MULTIVARIATE_COMMON_HEIGHT_BOUNDED_GAP_TAIL_CLOSURE.md).
The [affine count-simplex test](MULTIVARIATE_AFFINE_TAIL_COUNT_POLYTOPE.md)
gives a separate decision procedure for any fixed list of
types. Here the envelope gives a direct certificate throughout
the stated classes, with an explicit dependence on the band width.

## Finding a mark before considering primes

Take adjacent orders less than \(J\) and denote their even order
by \(m\ge4\). Work first with any nonempty set of row types from
the band. At order \(j\) a type \((a,b)\) has factorial
argument
\[
                         L_{j,i}(k)=(j-k)a+kb .
\]
Let \(f_j(k)=\max_i L_{j,i}(k)\) at the integers \(0\le k\le j\).
Since the maximum of affine functions is convex, every strict
sublevel set on this grid is an integer interval. Put
\(F=\min_{0\le k\le m} f_m(k)\). Take \(x=F+1\) if this
is nonzero modulo \(m\), and \(x=F+2\) otherwise. In either
case \(m\nmid x\). At a minimizer of \(f_m\), all smaller-end
arguments \(m\min(a,b)\) are at most \(F<x\). After those
scalar factorials have been removed, a selected type thus
includes the literal factor \(x\) precisely when
\(L_{m,i}(k)\ge x\). The indices where no selected type
contributes \(x\) form the nonempty interval
\[
                          I=\{k:f_m(k)<x\}.                   \tag{2}
\]

If \(x=F+1\), this interval has width at most one. In
fact, if three consecutive indices lay in \(I\), every line
at their middle index would be at most \(F-1\): use the right
neighbor for an increasing line and the left neighbor for
a decreasing line. This contradicts the definition of \(F\).
The same argument works when \(x=F+2\) and every gap is at
least two.

Suppose instead that \(x=F+2\) and \(I\) has width at
least two. At the center of any three consecutive indices
in \(I\), every line is at most \(F\), and some line
attains \(F\) and has absolute slope one. Since
\(F\equiv-1\pmod m\), this center is \(m-1\) for
a slope of \(+1\) and \(1\) for a slope of \(-1\).
Reversing both exponent columns if necessary, take the
center to be \(m-1\) and write \(F=mZ-1\).
A unit increasing row ends at exponent \(Z\). All increasing
upper exponents are at most \(Z\), because \(f_m(m)\le mZ\).
Every decreasing row ends at an exponent \(b\le Z-1\),
because its argument at \(m-1\) is \(mb+d\le mZ-1\),
where \(d\) is its gap. At \(m-2\), the increasing arguments
are at most \(mZ-2\), and the decreasing arguments \(mb+2d\)
are even. As \(mZ-1\le f_m(m-2)\le mZ\), a decreasing row
therefore attains \(mZ\). This gives
\[
                     2d=m(Z-b),\qquad d\ge m/2 .             \tag{3}
\]
Thus this exceptional case is impossible when \(D<m/2\),
including the all-unit-gap case.

When \(m=N\) is the upper order, there is a replacement
mark at order \(j=m-1\) even in the exceptional case. The
decreasing row from (3) has argument \(jb+2d=jZ+(Z-b)\)
at \(j-2\), so \(f_j(k)\ge jZ+1\) for all \(k\le j-2\).
For any decreasing row the bound at \(m-2\) gives
\(2d\le m(Z-b)\). Hence at \(j-1\) its argument is at most
\[
                 jb+d\le jZ-(m/2-1)(Z-b)\le jZ-1 .
\]
The increasing unit row has argument \(jZ-1\) exactly
there. At \(j\), the maximum is \(jZ\), attained by that
unit row. Thus the minimum of \(f_j\) is \(jZ-1\). Choose
\(x=jZ+1\), which is nonzero modulo \(j\). It gives
\(I=\{j-1,j\}\) and is above every smaller-end argument
by the same argument used at (2). The result in each
stated case is a marked order \(j\) and an integer \(x\)
such that \(j\nmid x\), \(I\) has width at most one, and
each index outside \(I\) samples \(x\) in at least one
selected type. In particular \(jM<x\le j(M+W+D)\).

## A count reserve and an exact divisor bound

Put
\[
 T=2D(W+1),\quad b_*=4JT,\quad Q=T b_*^T,\quad
 H=J^2(W+D),\quad L=\operatorname{lcm}(1,\ldots,H).
\]
There are at most \(T\) types. Let \(r\ge1\) be the
tail count. For \(s=0,\ldots,T\) put
\(t_s=(r/T)b_*^{-s}\). At least one of the \(T+1\)
disjoint intervals \([t_s/b_*,t_s)\) is free of positive
type multiplicities. For this \(s\), select all types
with multiplicity at least \(t=t_s\). At least one type
is selected since some multiplicity is at least \(r/T\).
Every selected type has multiplicity at least \(t\ge r/Q\)
and the total omitted count \(r_0\) is at most
\(Tt/b_*=t/(4J)\). Apply the envelope argument to this
selected set.

Write \(j\) for the marked order and \(j'\) for the other
one. Strip the smaller-end factorial of every row at each
order. Let \(B_\ell\) be \(\ell!\) times the product of the
full stripped ramps of the fixed core in order \(\ell\).
Each stripped core coefficient, including its binomial
coefficient, divides \(B_\ell\).

For a prime \(p\mid x\) put \(E=v_p(x)\) and \(L_p=v_p(L)\).
Assume \(E>L_p\). At the marked order all tail ramp
factors other than \(x\) differ from \(x\) by nonzero integers
of absolute value at most \(H\). Their \(p\)-valuations
are at most \(L_p\). In the other degree a tail ramp
factor \(G=j'M+v\) satisfies
\[
                 jG=j'x+(jv-j'u),\qquad x=jM+u .
\]
Here \(0<u\le j(W+D)\), \(0<v\le j'(W+D)\),
and \(jv-j'u\ne0\) since \(j\nmid x\) and
\(\gcd(j,j')=1\). The absolute value of this offset
is at most \(H\). The \(p\)-valuation of \(jG\) is
therefore at most \(L_p\) and likewise for \(G\).
This argument allows primes dividing \(j\).

At order \(j\), a coefficient valuation is \(E\) times
its literal count of \(x\) plus a remainder in
\([0,v_p(B_j)+rjD L_p]\). For any \(k_*\in I\)
that count is at most \(r_0\). At any \(k\notin I\)
the count is at least \(t\). All other-order ordinates
are in \([0,v_p(B_{j'})+rj'D L_p]\), so every slope \(q\)
of their lower Newton polygon has absolute value at most
the right endpoint of this interval. Subtracting \(qk\)
from the marked ordinates, an outside index then exceeds an
inside index by at least
\[
 (t-r_0)E-v_p(B_j B_{j'}^j)-r j(1+j')D L_p .
\]
Since \(j,j'<J\), \(j(1+j')\le J^2\) and
\(t-r_0\ge3r/(4Q)\). The difference is strictly positive
whenever
\[
                    E> 2Q\bigl(v_p(B_j B_{j'}^j)+J^2D L_p\bigr).
                                                               \tag{4}
\]
Any shared-slope face of the marked polygon then has
both endpoints in \(I\). Their cumulative width is at
most one. The Newton polygon product rule transfers
the slope multiset of a rational common factor to both
parents. The rule applies after clearing polynomial
denominators and Gauss normalization. See
[Filaseta and Rocks, Theorem 2.1](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf).
The even exact moment is strictly positive at every real
argument by the independent Gamma integral, since the
two distinct monomials cannot be proportional as functions.
Its nonconstant rational factors have even degree.
Consequently (4) certifies coprimality.

A common root therefore forces \(x\mid K\), where
\[
                  K=L^{\,2QJ^2D}
                       \left(\prod_{2\le \ell<J}B_\ell^J\right)^{2Q}.
\]
The primewise bound includes \(E\le L_p\) directly.
This is impossible for \(M>K\) because \(x>jM\).
The bound is independent of the tail multiplicities.
Using \(L\le H!\) gives (1) from the definitions of \(Q,H\).
All degrees at least \(J\) use the absolute theorem.
A rational Bézout identity extends each coprimality
assertion to any characteristic-zero field.

## A balanced obstruction to using just the literal count polygons

For every even \(m=2d\ge4\), consider orders \(m,m+1\) and
equal positive multiplicities \(r\) of the four types
\[
 (M,M+d),\qquad (M+1,M),\qquad
 (M+d,M+d+1),\qquad (M+2d,M+d).
                                                               \tag{5}
\]
The signed total gap is zero. The smaller exponents have band
width \(d\). For any \(\tau>0\) one can take
\(r=\lfloor \tau M/[2(d^2+1)]\rfloor\) as \(M\to\infty\)
to obtain total curvature tending to \(\tau\) and every
relative gap tending to zero.

In the literal affine count polygons of the two orders, every
mark has at least two shared slope units. This includes marks
with proportional affine factors at both orders, as permitted in
the [count-simplex certificate](MULTIVARIATE_AFFINE_TAIL_COUNT_POLYTOPE.md).
To see this, divide the multiplicities by \(r\), and index one
proportional class by the normalized intercept \(q\): its
possible factor at order \(j\) is \(j(M+q)\), which exists in
the integer factorial ramps only if \(jq\) is an integer.
All classes have \(0<q\le2d\). At \(0<q\le d\), only the
first two types in (5) count. At \(d<q\le2d\) the count is
the reversal of the same calculation with \(q-d\) in place
of \(q\). It suffices to treat the first interval.

For noninteger \(q\), at most one of \(j=m,m+1\) can have
\(jq\) integral, so the other polygon is horizontal.
Write \(t=jq\) at the order with the mark. If \(0<q<1\),
the increasing row counts for \(k\ge\lceil t/d\rceil\),
and the decreasing row counts for \(k\le j-t\). The
intervals cover the grid because \(\lceil t/d\rceil\le2\)
and \(j-t\ge1\). Both endpoint counts equal one, giving a
horizontal polygon of width \(j\). If \(1<q\le d\), only
the increasing row counts. Its first index is
\(\lceil t/d\rceil\ge3\), giving an initial horizontal edge
of width at least two.

For integer \(q=1\) the decreasing row counts only at zero.
The increasing row first counts at \(k=2\) for \(j=m\)
and at \(k=3\) for \(j=m+1\). The lower polygons have a
first slope \(-1\) of width one and a last slope
\(1/(m-1)\) of width \(m-1\) in both orders. The upper
order also has a horizontal edge of width one. For any
integer \(q=h\) with \(2\le h\le d\), the increasing row
starts at \(2h\) and \(2h+1\) in the two orders, respectively,
and the decreasing row contributes nothing. The final slope
is \(1/[2(d-h)+1]\) of width \(2(d-h)+1\) in both
orders. The upper order's horizontal width exceeds the
lower order's width by one, so all \(m\) slope units of
the lower order are shared. Reversing at \(q>d\) negates
the slopes in both orders and preserves the shared widths.
Multiplication by \(r\) also preserves them.

Thus the signed-gap-zero simplex can fail the literal
count test at an odd terminal, even for a fixed
four-type template and arbitrarily large heights.
The obstruction concerns the leading count polygons.
The higher valuations and residue coefficients of the
ordinary moments remain available for arithmetic separation.
