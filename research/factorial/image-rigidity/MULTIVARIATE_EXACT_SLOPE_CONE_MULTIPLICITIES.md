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
and zero-strength tests on these families then give a
computable all-degree threshold for the original arbitrary
counts when the resulting tests pass.

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
First translate the height parameter by an integer so that
all smaller-endpoint intercepts are nonnegative.
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

## A circulation condition for unit gaps

Assume all absolute gaps are one and \(c_i\equiv1\bmod J!\).
For a paired order \(j\) put
\[
      \ell_i=\min(a_i,b_i),\qquad
                     q_i(j)=(j\ell_i+1)/(j c_i).
\]
At adjacent orders \(n,N\), make two vertex sets from the
distinct \(q_i(n)\) and \(q_i(N)\), respectively. Give type
\(i\) its own edge between the corresponding vertices.
Orient an increasing type from the lower to the upper
order and a decreasing type in reverse. Parallel types
are separate edges.

**Corollary.** Beyond a computable \(M\), a physical common
root requires \(\boldsymbol m\) to be a nonnegative integer
circulation on this directed graph. Thus a template whose
graphs have no directed cycle at every relevant pair is
coprime in every degree for large \(M\), uniformly in the
independent counts. A graph with just one directed simple
cycle, considered up to cyclic rotation and with distinct
edge identities, reduces the counts to one exact ray.

In fact the denominator of \(q_i(j)\) contains the part \(j\)
unreduced because \(j\ell_i+1\equiv1\bmod j\) and the slopes
are coprime to \(j\). A factor in the other order \(j'\)
cannot have that normalized offset, because \(j',c_i\)
are coprime to \(j\) for every type. Within order \(j\),
proportionality to the first mark forces a block offset
\(h'\) to satisfy \(h'c_i\equiv c_{i'}\equiv1\bmod j\),
where \(i'\) is the type contributing the factor. With
unit gaps \(1\le h'\le j\), this gives \(h'=1\).
Thus the good prime at each vertex forces the signed
sum of the counts there to vanish, by the endpoint
minimum argument in the two-window proof. At most one
vertex mark can lack a good prime. The other order then
has all its vertex sums zero, hence so does the full
signed gap, which for unit gaps is the signed sum of
all counts. Subtracting the other sums in the order
of the missing vertex forces its sum to vanish too.

Conservation at every vertex is the circulation condition.
Any positive edge of a finite nonnegative circulation
leads to a directed cycle by following positive outgoing
edges. Subtract the least count on that cycle and iterate.
Each step removes a positive edge, so this gives a
decomposition into directed simple cycles and proves both
claims, including parallel opposite edges.

There are arbitrarily long paths of overlapping first-mark
classes in the acyclic condition. Fix one adjacent pair
\(n,N\), begin with \(c_0=1,\ell_0=0\), and for each successive
type put
\[
 c_i=Q_i c_{i-1},\qquad
 \ell_i=Q_i\ell_{i-1}+(Q_i-1)/j_i,
\]
where \(Q_i>1\) is any integer congruent to one modulo
\(J!\) and \(j_i\) alternates between \(n\) and \(N\).
For neighboring types the offset difference at order \(j\)
is \((Q_i-1)/c_i\) times \(1/j_i-1/j\).
For non-neighboring types, coincidence would require
\(1/j\) strictly between \(1/N\) and \(1/n\), which no
integer order supplies. The incidence graph at the chosen
pair is a path. At other pairs it consists of isolated
edges and paths of length two. Every choice of orientations
therefore satisfies the acyclic test.

## Polynomial height extension

The exact cone test also applies to row types with polynomial
smaller exponents \(a_i(M)\in\mathbb Q[M]\) of positive degree and
positive leading coefficient, still with fixed nonzero signed
gaps. The heights and the assigned counts must be ordinary.
Use as marks the distinct primitive irreducible divisors of
all nonconstant literal ramp polynomials in the two orders.
For a mark \(F\), a ramp factor of order \(b\) along \(F\)
contributes \(b m_i\) at each coefficient index where it occurs.
This defines the literal polygons and cones \(K_F\).

Enlarge \(\mathcal S\) by the primes of all coefficient
denominators, constant multipliers and pairwise resultants of
the distinct primitive divisors. If \(p\notin\mathcal S\)
divides \(F(M)\), no other divisor vanishes there and the
literal arrays equal the native valuations divided by
\(v_p(F(M))\). The same finite cone calculation applies.
Every nonlinear divisor has at least two distinct complex
roots by irreducibility in characteristic zero. The effective
\(S\)-part inequality of
[Bugeaud, Evertse and Győry, Theorem 2.2](https://arxiv.org/pdf/1708.08290)
gives an outside prime for each nonlinear mark at every
sufficiently large \(M\). The preceding \(S\)-unit argument
still allows at most one omitted linear mark. Thus (1) applies
with the omission restricted to linear divisors. When none
are present one intersects all the cones.

Criterion (2) still reduces all potential failures to finitely
many fixed count vectors and exact balanced rays. For a fixed
count vector the curvature of the moving rows tends to zero
and its signed gap is fixed, so it is enough to check its
twisted zero-strength core pair. For a balanced ray use the
[polynomial-scale divisor-rank test](MULTIVARIATE_POLYNOMIAL_SCALE_DIFFUSE_EXCLUSION.md)
with constant \(C_i,S_i\). The zero-strength tests described
above suffice for an empty or angularly regular core. The
cutoff in degree then gives all-degree coprimality when the
tests pass.

For example, take two windows of types with signed gaps
\(+2,-1\) at the common smaller exponent \(M^2+1\) and
\(-2,+1\) at the common smaller exponent \(M^4+16\).
All multiplicities are arbitrary. The normalized ramp
polynomials in the first window are \(M^2+1+h'/j\) with
\(0<h'/j\le2\); their complex root moduli are at most
\(\sqrt3\). The second window's moduli are at least \(2\).
Thus irreducible divisors do not overlap between windows.
Within a window, different normalized offsets give
relatively prime polynomials. In either order \(j\), the
first offset \(1/j\) appears in both types of its window
and does not appear in the other order. All primitive
divisors are nonlinear because the ramp polynomials have
no real zeros. Apply the good prime condition to one
such divisor at a first offset to force equal counts in
that window. The signed-gap bound then gives equal main
counts in both windows up to one of finitely many shifts.

The divisor-rank proof for two opposed affine windows
applies with the same tail sums: every offset above the
short ramp that is not integral is exclusive to the long
row in its window and its order. Choose any primitive
divisor at that offset to read off its tail sum. Both
windows have long gap two and short gap one, so the two
reflected index ranges again cover every position. The
polynomial-scale test therefore yields all-degree closure
above a computable \(M\) for this nonempty template, alone
or with any fixed angularly regular core. At each tested
order every first-block factor is shared within its
window, since \(j\ge2\), so the
[separating-divisor test](MULTIVARIATE_POLYNOMIAL_HEIGHT_NATIVE_MARKS.md)
does not supply this conclusion.
