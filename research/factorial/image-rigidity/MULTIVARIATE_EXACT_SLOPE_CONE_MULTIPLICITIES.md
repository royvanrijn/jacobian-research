# Exact reduction of affine multiplicities by good prime factors

Fix a finite ordinary two-monomial core and \(t\) unequal row
types \((c_iM+a_i,c_iM+b_i)\), where \(c_i>0,a_i,b_i\) are
integers. The counts \(\boldsymbol m\in\mathbb Z_{\ge0}^t\)
can vary arbitrarily with the positive integer \(M\). At least
one moving row is present. Write \(d_i=b_i-a_i\).

At each consecutive pair \(n,N=n+1\), use the literal factors
of the [count-simplex test](MULTIVARIATE_AFFINE_TAIL_COUNT_POLYTOPE.md#the-rational-test).
A mark \(x_\mu(M)\) is a primitive affine integer form of
positive slope representing a proportionality class of stripped
moving factorial factors in the two orders. An occurrence of
that mark at a coefficient index contributes the row count.
Let \(K_\mu\subset\mathbb R_{\ge0}^t\) consist of count vectors
for which the two resulting literal Newton polygons share
at least two unit slopes, counting edge lengths.

**Theorem.** Each \(K_\mu\) is a finite union of closed rational
polyhedral cones. Beyond a computable \(M\) depending only on
the core, template and orders, every physical common root
forces
\[
        \boldsymbol m\in\bigcup_{\mu_0}
                               \bigcap_{\mu\ne\mu_0} K_\mu .     \tag{1}
\]
In particular there is a finite reduction for arbitrary counts.
At each \(4\le N<J\), with \(J\) the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
express (1) as a union of cones \(K\). Suppose every \(K\)
satisfies
\[
                 \dim(K\cap\{\boldsymbol d\cdot\boldsymbol m=0\})
                      \le 1.                                   \tag{2}
\]
Then all possible large-height failures occur in a computable
finite list of fixed count vectors or exact families
\[
          \boldsymbol m=r\boldsymbol c+\boldsymbol s,\qquad
          r\in\mathbb Z_{\ge0},\qquad
          \boldsymbol c\in\mathbb Z_{\ge0}^t\setminus\{0\},
          \qquad \boldsymbol d\cdot\boldsymbol c=0 .
\]
The [affine divisor-rank test](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md)
and zero-strength tests on these families then give an
effective all-degree test for the original arbitrary counts.

## All but one mark gives an exact slope constraint

Fix the orders. Choose a finite prime set \(\mathcal S\)
containing the prime factors of \(N!\), the core ramps in both
orders, the constant multipliers relating the moving factors
to their primitive marks, and every nonzero cross determinant
of two distinct marks. An outside prime \(p\notin\mathcal S\)
dividing \(x_\mu(M)\) cannot divide any other mark by the cross
determinants. The actual valuations at \(p\) are exactly
the literal arrays multiplied by \(v_p(x_\mu(M))\).
The shared slopes must accommodate the common factor by
the Newton polygon product rule, as in the count-simplex
transfer. Even positivity makes the common factor have degree
at least two. Thus \(\boldsymbol m\in K_\mu\) for this prime.

For distinct \(x_\mu=B_\mu M+A_\mu\) and \(x_\nu=B_\nu M+A_\nu\)
write \(\Delta=B_\nu A_\mu-B_\mu A_\nu\ne0\). If two positive
mark values are both supported on \(\mathcal S\), then
\[
            (B_\nu/\Delta)x_\mu-(B_\mu/\Delta)x_\nu=1 .
\]
This is a fixed equation in the two-fold product of the
finitely generated rational \(\mathcal S\)-unit group. The
effective height bound of
[Bérczes, Evertse and Győry, Theorem 2.1](https://math.leidenuniv.nl/reports/files/2008-16.pdf)
bounds both mark values from the pair and \(\mathcal S\).
There are finitely many pairs and all marks tend to infinity,
so beyond a computable threshold at most one mark lacks an
outside prime. This proves (1).

For each mark enumerate possible polygon vertex chains.
Having a convex chain below the ordinate array is given
by homogeneous rational linear inequalities in
\(\boldsymbol m\). The chain is then the convex minorant,
since it passes through its own vertices of the array.
Enumerate two distinct matched positions in the two unit
slope lists of each pair of chains. The matches impose
two homogeneous linear equalities in the counts. This
gives the finite cone description even with collinear
vertices or ties.

## Reduction on a fixed signed-gap slice

The [signed imbalance bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
at a collision and the fixed core give a finite list of
possible \(h=\boldsymbol d\cdot\boldsymbol m\).
The recession cone of the nonempty polyhedron
\(P=K\cap\{\boldsymbol d\cdot\boldsymbol m=h\}\) is exactly
\(K\cap\{\boldsymbol d\cdot\boldsymbol m=0\}\).
The nonnegative orthant makes \(P\) contain no lines.
If the recession cone is zero, \(P\) is bounded and its
integer vectors can be enumerated. Otherwise (2) gives
a single ray with primitive integer generator
\(\boldsymbol c\ge0\). Write \(\boldsymbol m=\boldsymbol v+
q\boldsymbol c\), with \(\boldsymbol v\) in the compact vertex
hull and \(q\ge0\). At an index with \(c_i>0\), set
\(r=\lfloor m_i/c_i\rfloor\). Then \(q-r\) belongs to
\([-v_i/c_i,1-v_i/c_i]\), so \(\boldsymbol m-r\boldsymbol c\)
belongs to a computable bounded box. Its integer values
give the required finite list of shifts.

The fixed count vectors are separated for large \(M\) by
the [all-directions affine ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md).
For an exact family apply the cited divisor-rank test at the
pair producing that family. At zero strength the reference
is the core's consecutive moments with upper argument
multiplied by \((N/n)^h\). A core with disjoint consecutive
arguments suffices for every such \(h\). For the empty core
and \(h\ne0\) the two endpoint binomials at zero strength
have distinct argument scales and no common root. For the
empty core and \(h=0\) the
[uniform weak-curvature cap](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
covers the zero-strength neighborhood. The absolute
curvature bound keeps \(r/M\) bounded at a collision, so
the compact positive-strength comparison covers the rest
once the template tests pass. The
[low ordinary orders](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)
and the absolute cutoff supply the remaining order pairs.

## All-degree closure for two opposed windows

Take four types
\[
        (U,U+D),\quad (U+e,U),\quad
        (V+D',V),\quad (V,V+e'),
\]
with \(U=cM+a,\ V=c'M+a'\), \(c,c'>0\) integers,
\(D>e\ge1,\ D'>e'\ge1\), and
\[
                       e/D+e'/D'\le1 .                          \tag{3}
\]
Assume \([a/c,(a+D)/c]\) and \([a'/c',(a'+D')/c']\) are
disjoint. The four counts \(x,y,z,w\) are independent.
For the empty core or any fixed ordinary core with disjoint
consecutive root arguments, these profiles are regular in
every degree above a computable height.

In each offset window take the first-factor mark in both
orders. At order \(j\), the first \(U\)-window factor is
\(jU+1\). Coprimality of \(n,N\) gives it no proportional
factor at the other order, since \(1/j\) is nonintegral
and cannot equal a fraction with the other denominator.
The disjoint normalized offset intervals prevent sharing
with the other window. Among the four first-factor marks
at least three have outside primes, so each window
supplies at least one. The literal arrays for a first
mark in the \(U\)-window are
\((y,x+y,\ldots,x+y,x)\) and a zero array at the other
order. When \(x\ne y\) the marked polygon has no horizontal
edge because the minimum is unique. This is incompatible
with a common factor, so \(x=y\).
The reflected argument in the \(V\)-window forces \(z=w\).

The bounded signed gap \((D-e)x-(D'-e')z\) now puts the
pair counts in exact shifted rays. The homogeneous
direction has positive entries proportional to
\((D'-e',D'-e',D-e,D-e)\). To verify full divisor rank for
this ray, use a relation on the logarithmic multipliers
and absorb the common endpoint normalization into the
raw weights \(\alpha_{j,k}\), as in the
[reflected-window divisor argument](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md#an-infinite-family-beyond-the-single-mark-test).
Put \(\alpha_{j,0}=-\sum_{k\ge1}\alpha_{j,k}\) and
\(A_{j,l}=\sum_{k\ge l}\alpha_{j,k}\).
The two divisor windows are disjoint. For every
\(l\ge\lfloor je/D\rfloor+1\), choose an offset \(h'\) of
the \(l\)-th \(D\)-block with \(h'>je\) and \(j\nmid h'\).
In the first such block take \(je+1\); later blocks
have at least two consecutive offsets. This divisor
comes only from the increasing \(D\)-row at order
\(j\). Its exponent in the relation is a positive
multiple of \(A_{j,l}\), so that tail sum vanishes.
In the \(V\)-window the divisors above \(je'\) come
only from the decreasing \(D'\)-row. Reflecting the
weights sends the tail sum at \(l\) to
\(-A_{j,j-l+1}\). Thus \(A_{j,k}=0\) also for
\(k\le j-\lfloor je'/D'\rfloor\). By (3) the two
ranges cover \(1,\ldots,j\). All weights in the
relation are zero, giving full rank at both orders.

The zero-strength tests pass for the given cores, so
the earlier reduction proves the claimed closure.
In particular the separated windows with gaps
\(+d,-1,-d,+1\) for any \(d\ge2\) admit arbitrary
independent multiplicities. Unequal slopes or gaps
in the two windows are permitted and can give a
nonzero limiting signed bias.
