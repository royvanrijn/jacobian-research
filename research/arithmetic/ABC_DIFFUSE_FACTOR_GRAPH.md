# Conditional counts and diffuse-coordinate graphs at the balanced endpoint

At the normalized \(1/3\) endpoint and minimum nonatomic total
\(\Gamma=1/2\), fixing the outside factors on any two letters
allows a uniform bound on the number of possible completions.
The triple family then gives a three-partite graph of triangles
with near-quadratic size at the scale of the outside factor sets.
The input is the [bounded-root full-modulus count](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots).

## A uniform pair-fibre bound

Suppose a limiting profile at \(d=1/3\) has \(\Gamma=1/2\). Then its three nonatomic masses are all \(c_T=1/6\). All three letters are numerically positive and none is finite-scale, by the established alternatives for zero height and finite-scale letters. Hence every letter provides positive-weight full moduli of arbitrarily small normalized cost.

Fix \(\delta>0\). On each letter \(T\) choose a fixed finite list of positive atoms omitting at most \(\delta\) of its positive atomic mass. Choose the zero window and disjoint windows shrinking sufficiently slowly around the selected atoms, so that their masses tend to the corresponding atomic masses. At fixed truncation, choose one anatomic representation per primitive triple. Define \(D_T\) to be the vector of its original root coordinates outside the chosen windows on \(T\). The number of possible outside vectors satisfies
\[
 \log\#\operatorname{im}(D_T)\le (1/6+\delta+o(1))L\log X.
 \tag{EN1}
\]
This uses the ordinary root-coordinate-box count. The \(o(1)\) takes the height limit at each fixed truncation before the small-budget limit.

Fix the values of \(D_U,D_V\) for two letters. On the remaining letter \(T\) choose a positive-weight full modulus of cost at most \(\delta\), with fullness comparable to \(1/L\), and count its possible values by freezing the corresponding roots. The modulus block is a fixed set of degree positions in the box, so its height and fullness ranges are uniform in the frozen values. In each selected nonzero window on the two free letters, group its roots at the representative exponent \(\lfloor r/L\rfloor\). The zero windows and the exponent differences give a low rational multiplier of height \(X^{o(1)}\). The amplified (GW) count applies to the remaining positive-window groups and low multiplier, uniformly in the fixed roots and target. Choose its normalized loss at most \(\delta\), by choosing a fixed sufficiently large amplification order and taking the small-budget limit. Its exact represented ratio, with fixed factors reinstated, determines both coprime complete terms and then the third. At fixed truncation the exceptional-coefficient and representation losses may be taken arbitrarily small. Thus the fibre bound is
\[
 \log\#\{\text{triples with given }D_U,D_V\}
 \le (2\delta+o(1))L\log X .
 \tag{EN2}
\]
All three pair bounds hold for the same choices of the original outside windows. The modulus choices may differ across the three applications.

## Consequences of the saturated count

Sample a triple uniformly from each box. Shannon entropy is taken with natural logarithms. Write \(Z\) for the sampled triple and \(N_0=L\log X\). Along the saturating sequence \(H(Z)=(1/3+o(1))N_0\). Since each \(D_T\) is a function of \(Z\), (EN1)–(EN2) give
\[
 \begin{split}
 H(Z\mid D_U,D_V)&\le (2\delta+o(1))N_0,\\
 H(D_U,D_V)&\ge (1/3-2\delta+o(1))N_0,\\
 (1/6-3\delta+o(1))N_0
    &\le H(D_T)\le(1/6+\delta+o(1))N_0,\\
 I(D_U;D_V)&\le(4\delta+o(1))N_0.
 \end{split}
 \tag{EN3}
\]
For example, the lower bound on \(H(D_U)\) follows by subtracting the capacity of \(D_V\) from the lower bound on their joint entropy. The mutual-information bound is the corresponding joint deficit.

These estimates localize to arbitrary coordinate sub-blocks of the two outside vectors, including selections whose degree indices depend on the box. Write \(S\subset D_U\) and \(R\subset D_V\) for their coordinate vectors and \(b_S,b_R\) for their normalized root-box capacities. At each truncation the capacities are the sums of the corresponding upper root-size exponents divided by \(L\). The complementary coordinate vectors have capacities at most \(1/6+\delta-b_S+o(1)\) and \(1/6+\delta-b_R+o(1)\). The chain rule and (EN3) give
\[
\begin{split}
 (b_S+b_R-4\delta+o(1))N_0
   &\le H(S,R)\le (b_S+b_R+o(1))N_0,\\
 (b_S-4\delta+o(1))N_0
   &\le H(S)\le(b_S+o(1))N_0,\\
 I(S;R)&\le(4\delta+o(1))N_0.
\end{split}
\]
The first lower bound subtracts the two complementary capacities from the whole-pair joint entropy. The single-block lower bound follows either by subtracting the complementary capacity from \(H(D_U)\) or by subtracting \(H(R)\) from the first bound. All capacity errors here are uniform in coordinate selection, since the number of coordinates is fixed before the height limit. In particular a fixed exponential saving in the possible paired outside vectors, or in such a sub-block of them, excludes saturation.

There is also a combinatorial formulation. Make a three-partite hypergraph with a part for the possible values of each \(D_T\) and one edge per sampled triple. Each part has size at most \(X^{(1/6+\delta+o(1))L}\). The maximum codegree of a pair in different parts is at most \(X^{(2\delta+o(1))L}\) by (EN2). Greedily retain an edge and discard all edges sharing two of its vertices. Each step discards at most three times that codegree. The resulting three-partite linear hypergraph therefore retains at least \(X^{(1/3-2\delta+o(1))L}\) edges. Thus at the level of normalized \(X\)-exponents the endpoint produces a dense partial Latin square on the outside factors.

## Rectangles with controlled prime overlap

Assume \(d=1/3,\ \Gamma=1/2\) and put \(c=1/6\). Fix
\(m\ge2\) and \(0<\xi<c\). For all sufficiently small fixed \(\delta>0\),
the thinned graph contains \(m\) rows and at least
\[
 X^{(\xi-4m\delta+o(1))L}
 \tag{FG1}
\]
columns sharing all these rows. A cell is one actual primitive triple.
A row or column is incident to its own cells, and the \(C\)-label of
a cell is incident to that cell alone. Let \(P(D)=\prod_e x_e^e\)
denote the outside product represented by an outside vector \(D\).
For each outside label \(D\) and each complete term \(N\) of any
nonincident cell, the rectangle can be chosen to satisfy
\[
 \operatorname{rad}\gcd(P(D),N)<X^{\xi L}.
 \tag{FG2}
\]
Its outside labels are distinct within each letter. On any fixed
finite subrectangle one can strip the outside products to pairwise
coprime factors at arbitrarily small numerical cost, with each
remaining prime appearing only on its designated letter in its
incident cells.

### Counting shared prime divisors

For any fixed \(N\le X\), at most \(X^{(c+\delta-\xi+o(1))L}\)
outside vectors on a given letter fail (FG2). Indeed, choose a
squarefree divisor \(g\mid N\) of size at least \(X^{\xi L}\)
dividing \(P(D)\), with divisor-many choices. Assign its primes
to the root coordinates they divide. At fixed truncation there
are \(J^{\omega(g)}=X^{o_X(1)}\) choices with \(J\) the number
of coordinates. If the product of the primes assigned to
coordinate \(e\) is \(g_e\), its root is a multiple of \(g_e\).
Counting multiples up to each box endpoint divides the trivial
capacity by \(\prod_e g_e=g\), saving \(\xi\) in the normalized
exponent. This count applies to all outside vertices of each
letter and hence also after thinning. Height errors at a
fixed truncation are absorbed before taking the small-budget
limit.

### Selecting many columns

Let \(E\) be the edges of the thinned graph viewed from \(A,B\),
with the unique \(C\)-label and completed triple attached to
each edge. By (EN2) and thinning, \(|E|\ge X^{(2c-2\delta+o(1))L}\)
and the size of each part is at most \(X^{(c+\delta+o(1))L}\).
Any pair of vertices from different parts determines at most
one edge. A star here consists of one column and \(m\) ordered
row neighbors. The power mean applied to the column degrees
bounds the number of these stars below by
\[
 |E|^m/|B|^{m-1}\ge
 X^{((m+1)c-(3m-1)\delta+o(1))L}.
\]
Repeated rows occur in at most
\(O_m(X^{(m(c+\delta)+o(1))L})\) stars.

For the nonincidences within one star, first fix the edge of the
complete term \(N\) and hence \(N\) itself. If the tested outside
label is on another row, choose it last using the divisor count.
If the label is the \(C\)-vertex on another row, choose that
label by the divisor count and use its pair with the common
column to recover its row. The column in a single-column star
is incident to every edge. Thus stars with a failed internal
nonincidence number at most
\[
 O_m\!\left(X^{((m+1)(c+\delta)-\xi+o(1))L}\right).
\]
For sufficiently small \(\delta\), the remaining stars account for
the lower-bound exponent. Fix an \(m\)-tuple of rows with at least
\(X^{(c-(4m-1)\delta+o(1))L}\) such columns.

For a fixed column its \(m\) edges determine all their complete
terms. Put a directed conflict from this column to a different
column if any of the latter's \(B\) or \(C\) products fails (FG2)
with a complete term of the former. For a \(B\)-product apply the
divisor count directly. For a \(C\)-product on a fixed row apply
the divisor count to that label and recover the column from the
row-label pair. The outgoing degree is therefore at most
\(X^{(c+\delta-\xi+o(1))L}\). In the undirected graph of conflicts
the average degree is at most twice this. A graph on \(M\) vertices
with average degree \(\bar d\) has an independent set of size at
least \(M/(1+\bar d)\): in a random vertex order retain vertices
that precede all their neighbors and use convexity of the reciprocal.
This gives (FG1). The selected columns give all cross-column
nonincidences of \(B,C\) labels, and those of row labels were
already checked within each column.

### Isolating the outside prime supports

By the [first-moment and support bounds](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#consequences-at-the-minimum-total-nonatomic-mass),
take a fixed \(R\) above the limiting support and drop all
coordinates with \(eL>R\). Their numerical weight tends to zero.
In the anatomic decomposition the only possibly nonsquarefree base
among the positive-index roots is the packed root at
\(k=\lceil\varepsilon_0^{-1}\rceil\). Since the truncation loss
\(\varepsilon_0=o(L)\), this root is absent from \(eL\le R\).
The retained roots have disjoint prime supports on each label, so
their product is the radical of the retained outside product.
The normalized root mass retained on each outside vertex is at
least \(c-o(1)\), because the atom windows omit the entire
nonatomic part in the limit. Thus \(\xi<c\) makes a repeated
outside product a prohibited overlap. Repeated rows were already
excluded and distinct columns were chosen. [Bernert–Browning–Lichtman–Teräväinen, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

For a fixed \(m\times n\) subrectangle, remove from each retained
outside factor the full power of every prime appearing in any
complete term of a nonincident cell. The removal costs at most
\(3mnR\xi\) in numerical exponent for each factor, since a retained
prime has valuation at most \(R/L\) and (FG2) bounds the radical
overlap with each term. The prime supports left on distinct roles
are disjoint. In a common incident cell, primitivity already gives
disjointness on separate letters. The removed factors and the
negligible high-degree tail enter the complementary factors of the
original identities. For fixed \(m,n\), taking \(\xi\) sufficiently
small before \(\delta\) makes the numerical loss arbitrarily small.
The stripping statement uses only fixed-size subrectangles.

## Separating complete values in the rectangles

At \(d=1/3,\ \Gamma=1/2\) put \(c=1/6\). A column wedge
consists of two distinct edges with the same outside \(B\)-vertex.
For any fixed \(0<\tau<1\), there is \(\eta(\tau)>0\) such that,
for sufficiently small fixed \(\delta\), the number of ordered
column wedges whose complete \(C\)-terms satisfy
\[
 \left|\frac{C_1}{K C_2}-1\right|\ll X^{-\tau}
 \tag{FG3}
\]
is at most \(X^{(3c-\eta(\tau)+o(1))L}\), uniformly in the positive
target \(K\), which may depend on the box and height. The same
holds for the complete \(A\)-terms of a column wedge. By the power
mean, all wedges together have exponent at least \(3c-5\delta\).

### Counting distinct terms and ratios in strips

Take either one copy or two independent copies of a single letter
from the box. The combined normalized root mass is at most \(2/3\),
the first moments sum to at most \(2\), and the nonatomic mass
is \(c\) or \(2c\), respectively. The
[one-root grouping](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#one-grouped-root-on-the-free-pair)
applies to these copies with the same bounds: one copy contains a
compact positive interval of nonatomic mass at least \(1/12\).
For phase width \(\theta\) it retains extra mass \(t\ge\theta/60\)
there and a finite selection of atoms of the copies. The respective
costs of fixing roots outside the windows are at most
\(c-t+t/16\) or \(2c-t+t/16\).

After fixing those roots the complete term, or the ratio of two
complete terms, is a fixed factor times the grouped monomial.
When applying the [real strip count](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights),
keep one tuple for each distinct monomial value. This gives exact
fiber bound one irrespective of common prime factors of the two
terms. The calculation for [fixed ratio strips](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#counts-near-a-fixed-ratio)
with \(\theta\) a fixed absolute multiple of \(\tau\) gives an
absolute \(k_0>0\) with the following bounds, uniform in the
positive strip target. A strip of relative width \(O(X^{-\tau})\)
contains at most \(X^{(c-k_0\tau+o(1))L}\) distinct complete values
of a letter or \(X^{(2c-k_0\tau+o(1))L}\) distinct ratios of two
copies of that letter. The single-value count also bounds the
number of outside vectors visible in that strip, up to divisor-many
representations per complete term at fixed truncation.

### Recovering outside labels from a ratio

Call a wedge good at threshold \(\xi\) if both comparisons of one
\(C\) outside product with the opposite complete \(C\)-term
satisfy (FG2). For a good wedge with reduced ratio \(C_1/C_2=u/v\),
the number of possible ordered pairs of outside \(C\)-labels with
this ratio is at most \(X^{(2\xi+o(1))L}\). To see this, drop roots
at scaled degrees above a fixed \(R\) as in the preceding
section. These roots cost \(X^{o(L)}\) by the moment and support
bounds, and the remaining roots are squarefree and pairwise
coprime within each term. In the first label, the product of all
outside root primes dividing \(C_2\) is less than \(X^{\xi L}\).
Choose that product and allocate its primes to the prescribed
degree positions, with divisor-many possibilities per product.
Every other outside root prime is visible in \(u\), giving
divisor-many possibilities for these roots and their positions.
The second label is handled with \(v\). This estimate allows
common factors in the complementary parts of the complete terms.

For each ordered pair of \(C\)-labels there is at most one wedge
at any fixed column by pair uniqueness. Hence good wedges
satisfying (FG3) cost exponent at most
\(3c-k_0\tau+2\xi+\delta\). Wedges failing the overlap bound
cost at most \(3(c+\delta)-\xi\): fix the opposite edge and use
the divisor count for the other \(C\)-label, which recovers the
other row together with the common column. Taking \(\xi\) small
compared with \(k_0\tau\) and \(\delta\) small compared with
\(\xi\) proves (FG3).

### Separated rectangles

For fixed \(m\ge2\) and fixed \(0<\tau<1\), one can choose a
rectangle with \(m\) rows and \(X^{(\epsilon+o(1))L}\) columns
for some \(\epsilon>0\), satisfying (FG2) and with all complete
\(C\)-values pairwise separated by more than \(X^{1-\tau}\).
The complete \(B\)-values from any two distinct columns can
also be required to have this separation. In particular, the
count at this endpoint cannot come from an affine realization
with a single complete value for each row and column and a row
diameter at most \(X^{1-\tau}\), since every column would then
have close complete \(C\)-terms.

All complete values at the endpoint are \(X^{1+o(1)}\), so a
difference at most \(X^{1-\tau}\) implies the relative strip
with target one and exponent \(\tau/2\). In the star count,
discard stars containing such a pair of \(C\)-values. First
choose the offending two-row wedge by (FG3), then choose the
remaining \(m-2\) row labels. This costs at most exponent
\((m+1)c-\eta(\tau/2)+(m-2)\delta\), which has a strict
deficit from the star lower bound for small \(\delta\).
Fix the rows as before. Add a conflict from a fixed column
to any other column if any of their complete \(C\)-values,
or any of their complete \(B\)-values, are too close.
For each target value in the fixed column the single-value
strip estimate bounds the number of possible outside labels
in the other column. A \(B\)-label is the column itself, and
a \(C\)-label with its fixed row recovers the column by pair
uniqueness. Thus these new conflicts have outgoing degree
at most \(X^{(c-k_0\tau/2+o(1))L}\). Take the overlap threshold
\(\xi\) smaller than this saving \(k_0\tau/2\) and take
\(\delta\) sufficiently small. The independent-set argument
still gives a positive column exponent. Separation of
\(B\)-values is only asserted across distinct columns.

### Prime sizes in a growing rectangle

Every outside root coordinate has normalized mass at most
\(\delta+o(1)\), for the selections in (EN1). Otherwise a
coordinate of mass strictly greater than \(\delta\) along a subsequence
would have bounded scaled degree by the moment bound.
The limiting measure outside the shrinking atom and zero windows
could have an atom at their centers if the windows shrank too
rapidly. Choose the windows slowly enough to recover the full
mass of each selected atom, as in (EN1). The limit of the outside
measure then has at most \(\delta\) total atomic mass. A coordinate
of mass strictly greater than \(\delta\) would give an atom of at least
that mass, a contradiction. In particular every prime in an
outside product is at most \(X^{(\delta+o(1))L}\).

A fixed bounded degree cutoff \(R\) discards only \(o(1)\)
numerical weight. Below it the roots are squarefree, so a single
prime power or one entire degree coordinate contributes at most
\(R\delta+o(1)\) to the numerical weight. For any fixed endpoint
profile the outside numerical weight on each letter is at least
the first moment of its nonatomic measure, which is positive.
Consequently an outside factor of positive numerical weight
cannot be supplied by a bounded number of primes or degree
positions as \(\delta\downarrow0\). This rules out using a bounded
number of private prime powers per label to match the diffuse
weights, irrespective of how those prime powers are prescribed.

If the column exponent in (FG1) is larger than \(\delta\),
all but at most \(X^{(\delta+o(1))L}\) columns have no prime of
their outside \(B\)-product exclusive to that column among the
outside \(B\)-products. Indeed, assign one distinct prime to
each column having such a prime and count the available primes.
This is compatible with small pairwise radical overlaps. The
ability to isolate all primes on a fixed subrectangle thus
does not extend to simultaneous isolation on the entire
growing rectangle.

## Complete-value codegrees from common powers

The separated rectangles need not have a single complete value
per outside label. Common-power decompositions give a lower bound
on how much the complete values must vary over those labels.

Suppose each complete \(C\)-value in the box can be written
\(C=WZ^q\) with \(q\ge2\) fixed before the height limit and with
at most \(X^{(u+o(1))L}\) possible coefficients \(W\). For example,
in the anatomic factorization \(u\) can be the normalized root
mass of the degree positions whose indices are not divisible by
\(q\), since the exceptional coefficient can be frozen at the
truncation loss. The value of \(q\) may vary between successive
small-budget boxes.

For fixed \(D\ne0\), at most \(X^{(2u+o(1))L}\) ordered pairs
of complete \(C\)-values satisfy \(C_2-C_1=D\).
Indeed, fix the two coefficients and use the
[uniform binary-form bound](ABC_OMITTED_COORDINATE_GRAPH.md#uniform-binary-form-bound)
for \(W_2Z_2^q-W_1Z_1^q=D\). Its constant is uniform in the
coefficients and nonzero difference of polynomial height.
The same argument handles a fixed affine difference
\(rC_2-sC_1=D\ne0\) with positive integer \(r,s\) of polynomial
height, by absorbing them into the coefficients.
[Zeng, Lemma 3.1](https://hexagonmath.org/pdf/2611.00621v1)

### Multiplicity over outside labels

Make the bipartite graph on the actual complete \(A,B\)
values, with an edge for each primitive triple. It is simple
because \(C=A+B\). For distinct complete values \(A,A'\), each
common complete \(B\) yields a pair of complete \(C\)-values
with fixed nonzero difference \(A'-A\), and that pair determines
\(B\). Thus the codegree of two complete vertices is at most
\(X^{(2u+o(1))L}\), in either orientation.

At \(d=1/3,\ \Gamma=1/2\), write \(c=1/6\). Suppose an outside
\(A\)-vector has at most \(X^{(\alpha+o(1))L}\) complete values,
and an outside \(B\)-vector has at most \(X^{(\beta+o(1))L}\),
after taking the outside-window loss to zero. Then
\[
 2\alpha+\beta+2u\ge c,\qquad \alpha+2\beta+2u\ge c .
 \tag{FD1}
\]
Indeed, if the complete graph has \(M\) edges and \(N_A,N_B\)
vertices, Cauchy–Schwarz bounds
\[
 M^2/N_B\le\sum_b\deg(b)^2
    \le M+X^{(2u+o(1))L}N_A^2.
\]
With edge exponent \(2c\) and vertex exponents at most
\(c+\alpha,c+\beta\), the diagonal term can dominate the
lower bound only if \(\beta\ge c\). Otherwise the codegree term
gives the first inequality. The other orientation gives the second.
These comparisons hold with arbitrarily small fixed outside-window
losses, which are then sent to zero. In particular
\(\alpha+\beta\ge (2c-4u)/3\). A realization with only one
complete value over each outside row and column needs \(u\ge1/12\).

If two letters \(B,C\) have such decompositions with the same
exponent and dictionary exponents \(u_B,u_C\), fixing any
complete \(A\) gives at most \(X^{(u_B+u_C+o(1))L}\) incident
triples: fix the two coefficients in the nonzero difference
\(C-B=A\). Hence the complete \(A\)-values have cardinality at
least \(X^{(2c-u_B-u_C-o(1))L}\), and some outside \(A\)-vector
has at least \(X^{(c-u_B-u_C-o(1))L}\) complete values whenever
that exponent is positive.

### Differences along two fixed rows

The rectangles of (FG1) can be taken with two fixed rows and
arbitrarily close to \(c\) in their column exponent, taking
\(\xi\) near \(c\) and \(\delta\) small. Their complete \(C\)-values
are distinct: an equality to a complete \(C\) in a different
cell would violate (FG2), since the outside radical on either
cell is at least \(X^{(c-o(1))L}\).
For fixed rows, any ordered pair of complete \(C\)-values
determines at most divisor-many columns. To see this, choose
an outside representation of the first complete value and
use pair uniqueness with the first row. It follows from the
fixed-difference bound that such rectangles have at least
\(X^{(c-2u-o(1))L}\) distinct nonzero complete \(C\)-differences
between the two rows, taking the column exponent to \(c\).
If both complete row values are fixed and the complete
\(B\)-values agree within each column, those differences are
constant. This gives the same threshold \(u\ge1/12\) even
when the rectangle uses only the fixed two rows.

## Higher complete-value codegrees with mixed powers

For two fixed coefficients and a nonzero difference, the curve
\(U x^p-V y^q=D\) is absolutely irreducible. Use Eisenstein
in \(x\) at a simple linear factor of \(Vy^q+D\) over
\(\mathbb C[y]\). For integer \(x,y\) with respective height
bounds \(X^{h_x},X^{h_y}\), the rectangular curve estimate bounds
the count by
\[
 X^{\min(h_x/q,h_y/p)+\varepsilon}
\]
up to a constant depending on the degrees and \(\varepsilon\).
Indeed the largest monomial scale is
\(X^{\max(ph_x,qh_y)}\); the height bound for the homogenizing
coordinate is \(1\). For complete terms at most \(X\), the
extra exponent is at most \(1/(pq)\). If \(p,q\) have a
common divisor greater than one, the uniform binary power bound
gives \(X^\varepsilon\) by rewriting both powers at that
common exponent. [Browning–Heath-Brown, equation (3)](https://arxiv.org/pdf/math/0506497)

### The curve for several fixed complete rows

Fix \(m\ge2\) distinct complete row values \(A_i\). At a common
complete column \(t\), suppose
\[
 t+A_i=U_i x_i^{q_i}\le X,\qquad q_i\ge2,\quad U_i\ge1,
 \qquad 1\le i\le m.
 \tag{MP1}
\]
For fixed coefficient and degree choices the number of common
complete columns is \(O_{\mathbf q,\varepsilon}(X^{1/Q+\varepsilon})\),
where \(Q=\prod_i q_i\). The same proof works for affine
replacements \(r_i t+s_i\) with nonzero leading coefficients
and distinct zeros, when the complete values on the right
are positive and at most \(X\).

The curve given by the equations \(x_i^{q_i}=(t+A_i)/U_i\)
has function-field degree \(Q\) over \(\mathbb C(t)\).
Adjoin the roots successively. At \(t=-A_i\) the right sides
from the earlier rows are units with nonzero constant term.
All their roots are formal unit power series in \(t+A_i\),
so the preceding field embeds in \(\mathbb C((t+A_i))\) and
has an integer-valued valuation with \(t+A_i\) of order one.
The next equation is Eisenstein for that valuation. The
coordinate ring over \(\mathbb C[t]\) is free with basis the
products of powers \(x_i^{j_i}\), \(0\le j_i<q_i\), by
the monic equations in separate variables. It injects into the
generic field, so the affine curve is geometrically integral.

Eliminate \(t\) to view the curve in the \(m\) root coordinates.
Its projection to \(x_i\) has degree \(Q/q_i\), since
\(\mathbb C(x_i)\) contains \(t\) and has degree \(q_i\)
over \(\mathbb C(t)\). A generic hyperplane fiber therefore
gives degree at least \(Q/q_{\min}\), where \(q_{\min}=\min_i q_i\).
Eliminating \(t\) using the coordinate of degree \(q_{\min}\)
gives the matching degree upper bound by Bézout, using the
remaining \(m-1\) equations of degrees \(q_i\).
All root coordinates have height at most \(X^{1/q_{\min}}\).
The coefficient-uniform affine curve bound now gives the
exponent \((1/q_{\min})/(Q/q_{\min})=1/Q\).
[Binyamini–Cluckers–Kato, Corollary 8](https://arxiv.org/pdf/2401.03982v2#page=11)

### Degree dictionaries and the number of common columns

Allow a dictionary of coefficient and exponent choices
of normalized exponent \(u_i\) for row \(i\), with all exponents
in its dictionary at least \(q_{i,*}\). The count of common
complete columns has normalized exponent at most
\[
 \sum_{i=1}^m u_i + \frac{1}{L\prod_{i=1}^m q_{i,*}}+o(1).
 \tag{MP2}
\]
The constants above can depend on the degrees, held fixed
before the height limit. For two rows the rectangular bound
or the common-divisor norm bound may also be used.

In the complete-value graph at the minimum nonatomic endpoint,
let the outside-to-complete multiplicity exponents be
\(\alpha,\beta\) on the two sides as in (FD1). Suppose each
complete \(C\)-value has a coefficient dictionary of normalized
size \(u\) for exponents at least \(q_*\). Then
\[
 m\alpha+(m-1)\beta+mu+\limsup\frac1{Lq_*^m}\ \ge\ \frac16,
 \tag{MP3}
\]
and the same inequality holds with \(\alpha,\beta\) exchanged.
Indeed, count ordered stars of \(m\) distinct complete \(A\)-values
at a complete \(B\)-vertex. If \(\beta<1/6\), their mean column
degree tends to infinity with a positive exponent.
The power mean applied to \((\deg B-m+1)_+^m\) gives star
exponent at least \((m+1)/6-(m-1)\beta\). On the other hand,
the row-tuple count and (MP2) bound it by
\(m(1/6+\alpha)+mu+\limsup 1/(Lq_*^m)\). For \(\beta\ge1/6\)
the inequality already holds. The symmetric count gives
the other orientation.

In particular an additive rectangle with \(m\) fixed complete
rows and column exponent tending to \(1/6\) requires
\(mu+\limsup 1/(Lq_*^m)\ge1/6\). The residual term tends to
zero if \(q_*\gg L^{-\theta}\) and \(m\theta>1\), even when
the exponents vary and need have no common divisor.

## Fixed affine differences on a narrow degree band

The common-power representation can be replaced by closeness of
the valuation degrees over a narrow band. Consider the endpoint
\(d=1/3,\ \Gamma=1/2\) with the thinned outside-factor graph.
Suppose the \(C\)-letter has root mass at least \(1/3-1/100\)
on the scaled degree interval \([3-1/500,3+1/500]\).
Fix two distinct outside \(A\)-rows, and fix positive integers
\(r_1,r_2\) and a nonzero integer \(D\), all of polynomial height
in \(X\). Then the number of their common columns whose complete
\(C\)-terms satisfy
\[
 r_2 C_2-r_1 C_1=D
 \tag{AD1}
\]
is at most \(X^{(1/6-\eta+o(1))L}\) for an absolute
\(\eta>0\), uniformly in \(r_1,r_2,D\). The outside-window
loss \(\delta\) is taken sufficiently small. In particular, the
fixed complete row and column realization of a rectangle of
near-full column exponent is excluded for this \(C\)-profile
even if the degree indices have no common power.

Set \(\xi=1/1000\). First count the columns in which the outside
\(C_1\)-product \(P\) satisfies
\(\operatorname{rad}\gcd(P,r_2 D)\ge X^{\xi L}\).
The [divisor allocation estimate](#counting-shared-prime-divisors)
applies to the fixed nonzero integer \(r_2D\) of polynomial
height and bounds the number of \(C_1\)-labels by
\(X^{(1/6+\delta-\xi+o(1))L}\). The first fixed row with such a
label determines the column by pair uniqueness. For \(\delta\)
sufficiently small this gives a strict saving.

For the other columns choose on \(C_1\) a block of outside
coordinates in the enlarged band \([2.997,3.003]\) with root
mass between \(7/100\) and \(71/1000\), up to limiting errors.
At least \(1/6-1/100\) of its nonatomic mass lies in the original
band. The outside windows retain this mass, and each outside
coordinate has limiting mass at most the omission parameter
\(\delta\), so taking \(\delta<1/1000\) gives such a block by
successive addition. Fix the coordinates of this block at
cost at most \(71/1000\). From their full product remove the
full powers of primes dividing \(r_2D\); the remaining modulus
\(n\) satisfies \(r_2 C_2\equiv D\pmod n\) with a unit target.
This removal depends only on the fixed factors and on
\(r_2,D\). The original roots in this bounded scale band are
squarefree, and their prime supports are disjoint. Therefore
the removal loses at most \(\xi\) in normalized root mass
and at most \(3.003\xi\) in numerical weight.

On \(C_2\) freeze all coordinates outside the original narrow
band, enlarged by a vanishing amount to capture its limiting
mass. This costs at most \(1/100\). Group the remaining
coordinates at power \(q=\lfloor3/L\rfloor\), leaving the
degree corrections in a wider rational cofactor. This gives
\(h=1/3,\ E=1/1500\), up to limiting errors, in (QD1), with the frozen factors absorbed
into the unit target. Primitivity across the two terms of
(AD1) is unnecessary. Each prime of \(n\) is absent from
\(r_2D\), and thus from \(C_2\) by the congruence.

Use \(kL\to\kappa=2.997\), residue coefficient \(b\ge.069\),
and numerical weight bounded by
\[
 .07\cdot2.997-.001\cdot3.003\le w\le .071\cdot3.003 .
\]
In the [one-root count](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-quantitative-congruence-estimate-for-one-grouped-root)
take \(\alpha=19/40,\ M=40,\ r=2/25\), giving \(B=4/75\).
The Taylor order bound is \(779/500<\kappa\), the determinant
threshold is at least \(53751/10^6>B\), and
\(2B/(1-\alpha)=64/315<w\). The count for each fixing is at most
\(X^{(w/\kappa+o(1))L}\). Each counted pair determines the
complete \(C_2\); (AD1) then determines \(C_1\). Its
outside labels are determined up to divisor-many choices.
Including both freezes, the exponent is at most
\[
 .071+.01+\frac{.071\cdot3.003}{2.997}<.153<\frac16 .
\]
Together with the large-overlap branch this proves the assertion.
The decimals in the parameter choices are exact rationals.

## Combining congruences from several fixed rows

At \(d=1/3,\ \Gamma=1/2\), put \(c=1/6\). Suppose the
nonatomic measure of the \(C\)-letter gives positive mass to some
bounded interval \([s,S]\) with \(s>6\). There are a fixed integer
\(m\ge2\) and \(\eta>0\) such that the following holds in the
thinned outside graph for sufficiently small window loss \(\delta\).
Fix \(m\) outside rows. Count their common columns for which, with
one integer parameter \(t\) in an interval of length at most \(X\),
\[
 C_i=r_i t+s_i\qquad (1\le i\le m).
 \tag{CT1}
\]
The nonzero integer slopes \(r_i\) and integer offsets \(s_i\)
are fixed and of polynomial height, and the zeros \(-s_i/r_i\)
are distinct. The count is at most
\(X^{(c-\eta+o(1))L}\), uniformly in these choices. In
particular an additive realization with fixed complete rows and
a single complete value in each common column cannot have column
exponent approaching \(c\).

Choose \(z\) strictly between \(1/s\) and \(c\). By using
continuity endpoints and moving the interval slightly within
\((6,\infty)\) if needed, we may assume that its limiting mass
outside the atomic windows is bounded below by a fixed \(b>0\).
Take \(m\) large enough that \(z/m<b\). For each row use a
fixed block of outside degree positions of scaled degree in
\([s,S]\), of normalized root mass between
\(z/m+o(1)\) and \(z/m+\delta+o(1)\).
These blocks exist by adding positions in the band one at a
time. Each outside coordinate has limiting mass at most
\(\delta\), by the slow-window selection used above.
The total cost of freezing the block factors is at most
\(z+m\delta+o(1)\), and those factors have combined numerical
weight at least \(sz+o(1)\).

Let
\[
 M=\left|\prod_i r_i\prod_{i<j}(r_i s_j-r_j s_i)\right|.
\]
It is a nonzero fixed integer of polynomial height. For each
row discard the columns for which its \(C\) outside product
has radical overlap at least \(X^{\xi L}\) with \(M\).
The [divisor allocation bound](#counting-shared-prime-divisors)
gives at most \(X^{(c+\delta-\xi+o(1))L}\) such columns,
since the outside \(C\)-label and the fixed row determine the
column. A fixed complete row has only divisor-many outside
representations, so the same estimate works when rows are
fixed by their complete values in the thinned graph.

In every remaining column strip all primes dividing \(M\) from
the frozen block on each row. These primes cost at most
\(S\xi+o(1)\) in numerical weight per block, since the roots in
the fixed scaled interval are squarefree and have disjoint
prime supports. Call the remaining block moduli \(n_i\).
The primes of \(n_i\) are coprime to \(r_i\), and the moduli
are pairwise coprime. Any shared prime would divide both
\(r_i t+s_i\) and \(r_jt+s_j\) and therefore divides
\(r_i s_j-r_j s_i\). After fixing the blocks, the congruences
\(r_i t\equiv -s_i\pmod{n_i}\) specify at most one residue
class modulo \(\prod_i n_i\). Choose \(\xi>0\) with
\(mS\xi<(sz-1)/2\). Then the product exceeds \(X\), so
the residue class contains at most one integer in the stated
interval. That integer determines each complete \(C_i\), hence
the outside labels up to divisor-many representations, and
the first row and its \(C\)-label determine the column.
The good-column bound is \(X^{(z+m\delta+o(1))L}\).
Now take \(\delta\) sufficiently small relative to \(c-z\) and
\(\xi\) to obtain the stated positive saving.

For the additive complete-value graph on the thinned edges,
write \(\alpha,\beta\) for the outside-to-complete
multiplicity exponents on the two sides as in (FD1). In the
star count of (MP3), (CT1) gives a common-column exponent at
most \(c-\eta\), independently of any power dictionary, by
taking \(r_i=1\) and \(s_i\) the distinct complete row
values. In a complete star the outside column label can vary
from row to row. The bad-overlap estimate still counts all
possible values of \(t\), because a complete row has divisor-many
outside representations and each outside row and \(C\)-label
recover the actual edge. The good-overlap reconstruction
counts \(t\) itself. The star calculation now gives
\[
 m\alpha+(m-1)\beta\ge\eta,qquad
 m\beta+(m-1)\alpha\ge\eta,
 \tag{CT2}
\]
after sending the outside-window loss to zero. The diagonal
case, in which the corresponding mean column degree can
remain bounded, already has \(\beta\ge c\) or \(\alpha\ge c\).
The use of the thinned edges preserves exponent \(2c\)
in this limit.

## Fixed affine differences on separated degree bands

Several affine relations can supply high-degree factors for a
single congruence on a free row. A shallow residue record keeps
the cost of the congruence count small even when the combined
full modulus has large numerical weight.
[One-root quantitative congruence estimate (QD2)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-quantitative-congruence-estimate-for-one-grouped-root)

Throughout, work in the \(d=1/3,\ \Gamma=1/2\) limiting graph. Put \(c=1/6\). The complete \(C\)-letter has total root mass \(1/3\), nonatomic mass \(c\), first moment \(1\), and bounded scaled-degree support. Fix the outside windows with omitted atomic mass at most \(\delta\). Every outside coordinate then has normalized mass at most \(\delta+o(1)\), by the diffuse factor graph argument. At fixed truncation, all bounded scaled-degree roots are squarefree because the packed index \(k_0\) satisfies \(k_0L\to\infty\). A bounded-scope prime overlap therefore loses only its radical mass at the level of original root factors.

### Selecting moduli and shallow residue records

Fix \(m+1\) outside rows of the thinned graph, labeled \(0,\ldots,m\), with complete \(C\)-values \(C_j\) in each common column. Consider the columns satisfying fixed affine relations
\[
 a_j C_0-b_j C_j=D_j,\qquad 1\le j\le m,
 \tag{AL1}
\]
where \(a_j,b_j>0\) and \(D_j\ne0\) have uniformly polynomial height. Assume \(D_j/a_j\) are distinct for different \(j\). In a genuine additive rectangle with fixed distinct complete rows and with a common complete \(B\)-value in each column, take \(a_j=b_j=1\) and \(D_j=A_0-A_j\); the assumptions hold.

Define the fixed nonzero integer
\[
 F=\prod_{j=1}^m a_j D_j
   \prod_{1\le i<j\le m}(a_jD_i-a_iD_j).
\]
Since \(m\) is fixed, \(F\) has polynomial height. Discard columns where \(\operatorname{rad}\gcd(P(C_j),F)\ge X^{\xi L}\) for at least one \(j\ge1\), where \(P(C_j)\) is the numerical outside product and \(\xi>0\) is fixed. The divisor-allocation bound on that single outside vector gives at most \(X^{(c+\delta-\xi+o(1))L}\) possibilities, and the fixed row with its \(C\)-label determines the column by pair uniqueness. The union still has a strict saving for \(\delta<\xi\).

On each remaining column, lock a prescribed set of outside root positions on letter \(C_j\), of normalized mass in \([\mu_j,\mu_j+\delta]+o(1)\) in a bounded scaled-degree range \([\kappa,B_j]\). The chosen positions can be the same in all rows if desired. For a layer of positive nonatomic mass exceeding \(\mu_j\), such a set exists independently of the particular column by successively adding sufficiently small outside coordinate boxes. The normalized sizes of fixed box positions, and not their particular numerical values, determine this selection. Strip all primes dividing \(F\) from each locked root. The stripped root mass loss per row is at most \(\xi+o(1)\), because the original roots are squarefree and pairwise coprime. The full product \(n_j\) of the remaining locked roots to their original degrees divides \(C_j\). Any prime shared by \(C_i,C_j\) divides \(a_jD_i-a_iD_j\), by eliminating \(C_0\) from (AL1). Hence the \(n_j\) are pairwise coprime. Their product \(n\) is \(k\)-full for \(k=\lfloor\kappa/L\rfloor\), and (AL1) determines a unit residue \(C_0\bmod n\) by CRT.

Write \(\mu=\sum_j\mu_j\). The total freeze cost is at most \(\mu+m\delta+o(1)\), including the unstripped roots. Uniformly in their values, the numerical weight \(w\) of \(n\), and the full possible record mass \(b_n=\log Q_n/(L\log X)\), satisfy
\[
 w\ge\kappa(\mu-m\xi)+o(1),\qquad
 b_n\ge \mu-m\xi+o(1),\qquad
 Q_n=\prod_{p^e\parallel n}p^{\lfloor e/k\rfloor}.
 \tag{AL2}
\]
The bound on \(b_n\) follows prime by prime: every retained prime occurs to the original degree \(e\ge k\), so \(p\) occurs at least once in \(Q_n\). The prime power notation allows the same inequality using a base valuation greater than one, although the present bounded-scale roots are squarefree.

Let \(0<b_0<\mu-m\xi\) be fixed. Choose a divisor \(Q\mid Q_n\) by adding allowed copies of individual primes until its normalized logarithmic size first reaches \(b_0\). Each allowed copy has normalized contribution \(\log p/(L\log X)\le\delta+o(1)\), because its prime came from an outside coordinate. Repetition of a prime up to its allowed multiplicity does not increase that increment size. Thus such a shallow record has mass in \([b_0,b_0+\delta+o(1)]\). The construction is deterministic for the fixed modulus, so it produces no extra freeze cost. The residue partition therefore uses the selected mass in place of \(w/\kappa\).

On \(C_0\) fix all factors outside any chosen free set and express the residual value as \(v z^q\) at one common approximate base exponent. Let \(f\) be the normalized freeze cost on \(C_0\), let \(h\) bound its grouped-root exponent \(\log_X H(z)/L\), and let \(E\) bound the numerical exponent of \(v\). Suppose the remaining data satisfy the strict inequalities of (QD2) with record lower bound \(b_0\), fullness \(\kappa\), and weight lower bound \(\kappa(\mu-m\xi)\). Then the good columns are counted at normalized exponent at most
\[
 \mu+f+b_0+(m+1)\delta .
 \tag{AL3}
\]
Indeed (QD2) counts \(v,z\) for each fixing and hence counts the complete \(C_0\). Each \(C_0\) has divisor-many outside representations at fixed truncation, and the fixed outside row \(0\) together with one such label determines the column by pair uniqueness. If \(\mu+f+b_0<c\), taking \(\delta\) sufficiently small gives a strict saving for all columns, combining (AL3) and the bad branch.

### Two rows, unequal bands

Suppose the \(C\)-measure has root mass \(17/60\) in
\[
 [30/17-.003,\ 30/17+.003]
\]
and mass \(1/20\) in \([9.99,10.01]\), with more than \(6/125=.048\) of nonatomic mass in the high band. For example, symmetry about each displayed center gives total first moment one. The main band may have an atom of mass \(1/6\) at \(30/17\), with the rest of its mass nonatomic, so the example can have total nonatomic mass exactly \(1/6\).

Use one locking row (\(m=1\)), \(\mu=6/125\), \(\xi=1/1000\), \(\kappa=9.99\), and \(b_0=1/25\). The free \(C_0\) freezes its high band, of mass \(f=1/20\), and groups the main band at \(q=\lfloor(30/17)/L\rfloor\). This gives \(h=17/60\) and \(E\le(.003)h=17/20000\). Take in (QD2)
\[
 \alpha=1/2,\qquad M=59,\qquad r=1/6,\qquad B=hr+EM.
\]
The exact margins are
\[
 \alpha(M+1)r=5<9.99,\qquad
 B=17527/180000<1/10=b_0\alpha(M+1)r/2,
\]
\[
 4B=17527/45000<46953/100000
       =9.99(.048-.001),\qquad r<30/17.
\]
The good exponent is at most \(69/500+2\delta=0.138+2\delta<1/6\). Thus every fixed affine relation of the stated type loses a fixed exponent on this mixed-band profile. The two bands are centered at \(30/17\) and \(10\).

### Several rows with a very light high band

Set \(f=1/200\), \(h=197/600\), and \(a=540/197\). Consider a formal \(C\)-measure with mass \(f\), all nonatomic, symmetric about \(20\) in \([19.99,20.01]\). Put the remaining mass \(h\) in two equal clusters symmetric about \(a-1/100\) and \(a+1/100\), each within radius \(1/100000\) of its center. For example at each lower center put atomic mass \(1/12\), and distribute the remaining mass of that cluster symmetrically and continuously. The total nonatomic mass is then \(1/6\), and the first moment is \(ha+20f=1\).

Use twelve locking rows, each with lock mass \(1/250\) from the high band, and take
\[
 \mu=6/125,\quad \xi=10^{-6},\quad
 \kappa=19.99,\quad b_0=21/500.
\]
Each row has high nonatomic mass \(1/200>1/250\), so it supplies the required block. On \(C_0\), freeze its high band and group both main clusters at \(q=\lfloor a/L\rfloor\). The two signs of the correction each have limiting numerical exponent \(E=h/200=197/120000\): each main cluster has mass \(h/2\) at mean signed displacement \(1/100\) from \(a\). The height of their rational ratio is bounded by the larger of the two correction weights. In (QD2) take
\[
 \alpha=1/2,\qquad M=63,\qquad r=5/16.
\]
Then \(r<a\) and the exact margins are
\[
 \alpha(M+1)r=10<19.99,\qquad
 B=49447/240000<21/100=b_0\alpha(M+1)r/2,
\]
\[
 4B=49447/60000<23982003/25000000
       =19.99(.048-12\cdot10^{-6}).
\]
The cost is \(19/200+13\delta=0.095+13\delta<1/6\). Thus a thirteen-row complete additive rectangle cannot have the near-full column exponent for this profile. The selected high mass in any individual row is below \(0.005\); the larger modulus comes from genuinely simultaneous locking.

Here is also a comparison to the natural one-row adaptation of the QD method that groups the entire two main clusters and freezes its high band. For any choice of common center, the correction exponent is at least \(E=h/200\). Indeed the maximum of the positive and negative correction integrals is minimized at the common mean \(a\), where the two integrals agree. With one locking row, let \(\mu_1\) be its root-mass freeze, \(w\) its available weight, and \(b\) its record coefficient after dividing by a possible fixed amplification order. A saving needs \(b\le b_*-\mu_1\), with \(b_*=c-f=97/600\). Since the maximal main degree is below \(69/25=2.76\), and the entire high band has mass \(f\) and maximum degree \(20.01\), write \(l_+=69/25\), \(A=(20.01-l_+)f\), giving \(w\le A+l_+\mu_1\). Consequently
\[
 wb\le \frac{(A+l_+b_*)^2}{4l_+}
       =\frac{4930487}{192000000}.
\]

A necessary condition from the determinant and reconstruction inequalities of (QD2), even with fixed amplification, is
\[
 E<\frac{wb}{32(2h-b)}.
 \tag{AL4}
\]
For proof without amplification, put \(K=2h/(\alpha b)>1\) and \(W=w(1-\alpha)/2\). The determinant forces \((\alpha b(M+1)/2-h)r>EM\), while reconstruction forces \(hr+EM<W\). Eliminating \(r\) gives
\[
 E<W\,\frac{M+1-K}{M(M+1)}
   \le \frac{W}{4(K-1)}
   =\frac{wb\alpha(1-\alpha)}{8(2h-b\alpha)}.
\]
Here \(x=M+1>K\) and \((x-K)/x^2\le 1/(4K)\), while \(x-1\ge x(1-1/K)\). Using \(\alpha(1-\alpha)\le1/4\) proves (AL4). With amplification by \(\ell\), apply the same steps with \(\ell h,\ell E,\ell b\); the resulting condition on \(E\) is at least as restrictive. For the stated measure the right-hand side is at most
\[
 \frac{4930487}{3041280000}
   <\frac{197}{120000}=E.
\]
This separates the simultaneous-lock certificate from the one-lock certificate with the same main-group representation, including arbitrary choices of \(k\), shallow record, and fixed amplification.

### A phase-concentration criterion for one letter

The single-letter version of the variable-mesh grouping supplies a different class of affine codegree restrictions, without confinement to a fixed degree band. Suppose on \(C\) there is a sequence of spacings \(a\in(0,1]\), widths \(\theta\to0\), and continuity windows \(\|r/a\|\le\theta\) in a fixed compact positive degree range. Let \(t\) be the captured nonatomic mass. Assume \(t/\theta\to\infty\). For each chosen spacing, select finitely many positive atoms omitting at most \(t/16\) of the total atomic mass, with weighted discrepancy at most \(\theta t/2\) from the unshifted lattice. This condition concerns only the \(C\)-letter. 

Fix one sufficiently advanced choice of \(a,\theta,t\). On \(C_0\) keep the windows, shrinking neighborhoods of those selected atoms, and the zero window. Group the positive scales at their nearest multiple of \(a\), using the common degree \(\lfloor a/L\rfloor\). The signed numerical correction costs at most \(2a\theta t+o(1)\): the retained nonatomic contribution is at most \(a\theta t\), unlisted atom contributions on the windows are at most \(a\theta t/16\), and selected atoms contribute at most \(a\theta t/2+o(1)\). The root height is at most \(X^{(2/a+o(1))L}\) by the first-moment bound, and the freeze cost is at most \(c-t+t/16\).

Use just one locking row \(C_1\), with fixed difference \(a_1 C_0-b_1 C_1=D_1\). Choose a compact positive scaled-degree range \([b,B]\) on it containing positive nonatomic mass. It can supply any sufficiently small fixed target numerical weight \(w>0\) at root cost at most \(w/b+o(1)\): remove coordinates of fixed positive macroscopic mass, take that cutoff slowly to zero, and accumulate the remaining root positions to the target numerical weight. Their increments tend to zero, since all their degrees lie in the chosen compact range. Let \(c_{\mathrm p}>0\) be the dimension-one constant in (UW). Choose \(K>\max(1,256/c_{\mathrm p})\) fixed and set \(w=K\theta\). Discard bad columns by the degree-independent divisor test for \(a_1D_1\), choosing \(\xi\) much smaller than \(w/B\) and \(\delta<\xi\). After stripping the bad primes the surviving modulus weight \(w'\) lies in \([w/2,w+o(1)]\). Since \(t/\theta\to\infty\), choose the phase instance so that \(w/b<t/16\) and \(8w/t<b\).

Take \(\kappa=8w/t\), \(C_1=32/(at)\) as the short height cap, and use (UW) without amplification. The actual modulus weight \(w'\ge w/2\) gives the required short height bound, since \((2/a)\kappa/w'\le32/(at)\), with a strict margin if the strip loss is chosen smaller and the height cap enlarged slightly. The sufficient wider condition is \((2a\theta t/w') C_1\le128\theta/w<c_{\mathrm p}\). The chart cost is at most \(2w'/\kappa\le t/4+o(1)\). For fixed frozen factors the exact monomial gives the complete \(C_0\), then the affine relation fixes \(C_1\), and the first fixed row recovers the column up to divisor-many choices. Thus the good exponent is at most \(c-t+t/16+t/16+t/4=c-5t/8\), and the bad columns have their own fixed saving \(\xi-\delta\).

For example a one-letter measure with mass \(1/6\) at scale \(3\) and mass \(1/6\) in the symmetric middle-thirds Cantor distribution on \([2,4]\) meets this concentration test along the meshes \(a=3^{-k}\). Its atom aligns exactly, and the first few zero ternary digits supply captured mass of order at least \(\theta^{\log2/\log3}\), so \(t/\theta\) is unbounded. Only this \(C\)-measure enters the codegree test.

## Conditional entropy of the complete values

The codegree restrictions apply to the distribution of complete
values over outside labels. In particular, a small fraction of
labels with many completions cannot account for the required
multiplicity in (CT2).

Use the uniform distribution on the thinned edge set at
\(d=1/3,\ \Gamma=1/2\), with \(c=1/6\) and \(N_0=L\log X\).
For an edge, put
\[
 I_A=-N_0^{-1}\log\Pr(A\mid D_A),\qquad
 I_B=-N_0^{-1}\log\Pr(B\mid D_B),
 \tag{EI1}
\]
where probabilities are evaluated at the sampled values.
Suppose the codegree of any \(m\) distinct complete row
values is at most \(X^{(c-\eta+o(1))L}\) uniformly, for some
\(0<\eta<c\). The complete graph retains the attached outside
labels of each edge. For any fixed \(\gamma>0\), the proportion
of edges on which
\[
 mI_A+(m-1)I_B\le\eta-\gamma
 \tag{EI2}
\]
is exponentially small in \(N_0\) for sufficiently small
outside-window loss. If the codegree bound holds in both
orientations, the same assertion holds with \(A,B\) interchanged.
In particular, after the window losses tend to zero, the
limiting entropies satisfy
\[
 mH(A\mid D_A)+(m-1)H(B\mid D_B)\ge(\eta-o(1))N_0,
 \tag{EI3}
\]
as well as the transposed inequality under the same two-sided
hypothesis. This applies to the joint congruence bound (CT1)
and to any uniform complete-value codegree saving with these
hypotheses.

For the proof, restricting to \(I_A\le u\) leaves at most
\(X^{uL}\) choices of \(A\) for each label \(D_A\), since each
chosen value has conditional probability at least \(X^{-uL}\).
The analogous restriction on \(I_B\) with parameter \(v\)
leaves at most \(X^{vL}\) complete column values per label.
Write \(n_A,n_B\) for the sizes of the two complete vertex
classes of the restricted graph. If \(M\) is its edge count
and \(K\) its \(m\)-row codegree bound, the ordered star count
gives
\[
 M\ll_m n_B+(K n_A^m n_B^{m-1})^{1/m}.
 \tag{EI4}
\]
In fact the number of ordered stars is at least
\(n_B(M/n_B-m+1)_+^m\), by applying the power mean to the
nonnegative shifted degrees. The same number is at most
\(Kn_A^m\). With the two outside capacities \(c+\delta+o(1)\)
the two exponents in (EI4) are consequently bounded by
\[
 c+\delta+v,\qquad
 2c+(2-1/m)\delta+u+(1-1/m)v-\eta/m,
\]
up to limiting errors. Cover the region (EI2) by finitely
many rectangles with upper corners \(u,v\) chosen on a grid
fine enough that \(mu+(m-1)v\le\eta-\gamma/2\).
Then both exponents save a fixed amount from \(2c\) for
small enough \(\delta\). The full thinned set has exponent at
least \(2c-2\delta+o(1)\), proving the probability bound.
Taking expectations of the nonnegative information variables
and then letting \(\gamma\) decrease to zero gives (EI3).

The same argument can use a coefficient and exponent dictionary
without bounding its worst-case size. For each complete
\(C\)-value choose one representation \(C=WZ^q\) with
\(q\ge q_*\), and put \(I_D=-N_0^{-1}\log\Pr((W,q))\).
For fixed \(m\), let \(\rho=\limsup 1/(Lq_*^m)\). If
\(\rho<c\), the mixed-power curve bound (MP2) implies that
\[
 mI_A+(m-1)I_B+mI_D+\rho\ge c-\gamma
 \tag{EI5}
\]
outside an exponentially small fraction, for any fixed
\(0<\gamma<c-\rho\) and sufficiently small window loss.
Indeed a restriction \(I_D\le u_D\) leaves at most
\(X^{u_D L}\) dictionary choices. In (EI4) its codegree
exponent is at most \(m u_D+\rho+o(1)\), and the same finite
grid argument has a strict deficit on the complementary
event in (EI5). This also yields the transposed inequality
and replaces \(u\) in the worst-case version (MP3) by the
normalized entropy of the selected dictionary.

## Cross ratios along a shared outside column

Fix two outside rows on the \(A\)-letter. For a column on the
\(B\)-letter write \((A_i,B_i,C_i)\) for its two incident
primitive triples. The complete \(B_i\) may vary with the row.
Let \(P_B\) be the product of the common outside \(B\)-factors.
For the reduced cross ratio
\[
 \frac rs=\frac{A_1/C_1}{A_2/C_2}
 \quad\text{in lowest terms},
 \qquad P_B\mid r-s,\quad r\ne s.
 \tag{XR1}
\]
Both \(A_i\) and \(C_i\) are units modulo \(P_B\) and agree
there, so the reduced ratio still gives a unit congruence
\(r/s\equiv1\bmod P_B\). If \(r=s\), the two primitive
fractions \(A_i/C_i\) are equal, so both complete triples
and their chosen outside labels agree. This contradicts
distinct rows.

For any fixed \(r/s\), the nonzero integer \(r-s\) has
polynomial height. It has divisor-many possible outside
products \(P_B\) by (XR1). Each such product assigns its prime
powers to the fixed degree positions in divisor-many ways,
since every root divides \(P_B\) and the number of positions
is fixed before the height limit. Hence the number of
columns producing that ratio is \(X^{o(L)}\) at fixed
truncation, with arbitrarily small exponent losses.
Consequently any dictionary of reduced cross ratios of size
at most \(X^{(v+o(1))L}\) covers at most that many columns
up to such losses. A near-full rectangle needs normalized
dictionary exponent at least \(c\).

The numerical weight of the common outside column gives another
necessary condition. If \(P_B=X^{w+o(1)}\), then the reduced
height \(\max(r,s)\) is at least \(X^{w+o(1)}\).
The conclusions use the outside product and the primitive
fractions, so they still apply when that outside label supports
many different complete \(B\)-values. Permuting the letters
and using the signed additive equation gives the corresponding
restrictions on the other two kinds of columns.

## Target ratios in a graph with a common full divisor

The powerful-number graph at one target counts common neighbours
by recovering the omitted factors from congruences. When the
targets vary, an identity in the complete cross ratio allows
one to count common neighbours by the matched products or ratios
of their targets. The full modulus can be held constant over
many different outside target labels.
[Powerful-number graph](POWERFUL_SUM_GRAPHS.md#codegrees-for-any-omitted-exponent)

### Recovery from a matched target ratio

Suppose two primitive triples \(A_i+B_i=N_i\) have reduced
cross ratio \(r/s=A_1B_2/(A_2B_1)\ne1\). Then
\[
 \big((r-s)A_1-rN_1\big)\big((r-s)A_2+sN_2\big)
   =-rsN_1N_2 .
 \tag{TR1}
\]
This follows by expanding and using
\(sA_1(N_2-A_2)=rA_2(N_1-A_1)\).
The product on the right is nonzero and of polynomial height.
Given \(r/s\) and \(N_1N_2\), the product determines \(N_1,N_2\)
up to divisor-many choices. Its signed factor pairs in (TR1)
then recover \(A_1,A_2\) up to divisor-many choices.
The same conclusion follows given \(r/s\) and the reduced
target ratio \(N_1/N_2\). If \(D=\gcd(N_1,N_2)\), both primitive
ratios \(A_i/B_i\) are \(-1\bmod D\), hence \(D\mid r-s\).
The scale \(D\) therefore has divisor-many choices, giving
the individual targets from their reduced ratio. The case
\(r/s=1\) gives identical reduced primitive fractions
\(A_i/B_i\), so it cannot occur with distinct outside row labels.

### Counting complete cross ratios by a full modulus

Suppose all targets share a fixed \(k\)-full divisor
\(n=X^{w+o(1)}\), with \(w>0\) and \(k\asymp 1/L\). With
outside labels fixed, write the free terms
\[
 A=P\,V_A\prod_i x_i^{e_i},\qquad
 B=Q\,V_B\prod_j y_j^{f_j}.
 \tag{TR2}
\]
In the endpoint boxes the short slots are the roots of the
selected positive atomic windows, with representative powers
of order \(1/L\) and heights \(X^{O(L)}\). The corrections
and zero-window factors give the wide factors of height
\(X^{o(1)}\). The fixed outside products are \(P,Q\), including
the frozen exceptional coefficients at the small truncation
loss. All variables are units modulo \(n\).

Choose one coloring of the primes by the short slots on both
letters. In each triple the short slots have disjoint prime
supports by primitivity and the anatomic decomposition.
A random coloring respects a fraction \(X^{-o_X(1)}\) of the
triples on average, because each triple uses \(o(\log X)\)
primes and the slot count is fixed before the height limit.
Keep a coloring that retains at least this fraction of the
edges. Now ratios of two roots from the same slot have pairwise
disjoint prime supports across different slots.

For fixed outside labels the relation \(A/B\equiv-1\bmod n\)
fits the [bounded-root estimate (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots).
There are \(X^{o(L)}\) choices of the exact ratio and hence
of its primitive triple, by the amplification described below.
For two fixed left labels with outside products \(P_1,P_2\),
the right label \(Q\) is shared within each wedge. It cancels
from the complete cross ratio, which is \(1\bmod n\).
After removing the fixed factor \(P_1/P_2\), the remaining
monomial uses the slot ratios \(x_{i,1}/x_{i,2}\),
\(y_{j,2}/y_{j,1}\) and the wide factor
\((V_{A,1}/V_{A,2})(V_{B,2}/V_{B,1})\). Thus (GW) also
bounds the number of exact complete cross ratios by
\(X^{o(L)}\), uniformly in the two left labels.

For rational short slots of disjoint prime support, the
exact-fibre argument used in (GW) applies to their numerators
and denominators separately. For fixed value of the signed
monomial and wide height \(n^\rho\), orient every power to
be positive by inverting a signed slot if necessary. On
each prime in a narrow numerator or denominator, retain at
most the ceiling of the matching valuation in the fixed
value divided by the slot power. The retained parts have
divisor-many choices. The residual part to power \(|a_i|\)
divides the opposite part of the wide factor. This gives fibre cost at
most \(n^{2\rho\sum_i1/|a_i|+\varepsilon}\) by the divisor
bound. To amplify, color primes by sign in the short slots
and independently by sign in the wide factor. Then a fixed
\(\ell\)-fold product preserves the support separation, has
no cancellation in any slot, and has divisor-many preimages.
Apply (GW) to these products and take \(\ell\)-th roots.
The counting loss tends to zero after normalization by
\(L\log X\), by first choosing \(\ell\) large and then taking
the small-budget limit. In the absence of short slots,
rational reconstruction gives the same vanishing loss
directly, since the wide height is \(X^{o(1)}\) and the
full modulus has positive numerical weight.

For a fixed box with at least one short slot, suppose the short roots in
(TR2) have heights at most \(n^{C_0/k}\), the wide factors
have heights at most \(n^{\rho_0}\), and (GW) applies to the
amplified heights \(2\ell C_0/k,4\ell\rho_0\).
If \(j_0\) is the total number of short slots and \(a_i\)
their signed powers, both the edge multiplicity and the
cross-ratio count cost at most \(X^\varepsilon\) times
\[
 K_0=n^{(j_0^2+j_0)/(\ell k)+8\rho_0\sum_i1/|a_i|}.
\]
Take one witness per outside edge in the colored graph.
For each pair of row labels consider only the matched
products \(N_1N_2\) and the matched reduced ratios \(N_1/N_2\)
among witnesses at common columns. If the smaller of the
two set sizes is bounded by \(D_0\) in both orientations,
(TR1) bounds the codegrees by \(K_0D_0X^\varepsilon\).
Consequently, for outside part sizes \(V_A,V_B\), the
original family count is
\[
 \ll X^\varepsilon K_0
 \min\{V_A V_B,\,
       V_A+V_B+\sqrt{K_0D_0}
             \min(V_A\sqrt{V_B},V_B\sqrt{V_A})\}.
 \tag{TR3}
\]
The coloring loss and all divisor losses are absorbed by
smaller intermediate exponents. The full product and
quotient sets of a positive polynomial-height target
dictionary of size \(T\) both have size \(T^2X^{o_X(1)}\)
at the exponent scale. The product statement follows by
divisor counting, and the quotient statement then follows
by Cauchy–Schwarz and equality of their multiplicative
energies. The smaller matched sets in (TR3) use the graph
pairing.

### Distinct target ratios in near-full rectangles

For every fixed \(m\ge2\) and \(\epsilon>0\), at
\(d=1/3,\ \Gamma=1/2\) there are \(m\) fixed outside
\(A\)-rows and at least \(X^{(c-\epsilon+o(1))L}\) common
outside \(B\)-columns, with the following property. For
each pair of those rows, their reduced target ratios
\(C_i/C_j\) are all distinct across the selected columns.
The complete \(B\)-values in one column can vary between
the rows.

To prove this, first freeze a prescribed positive-weight
full block on \(C\) at arbitrarily small normalized cost
\(\mu>0\). Such a block exists in a compact positive degree
band with nonatomic mass. Choose a fixing that retains
at least \(X^{(2c-\mu-2\delta+o(1))L}\) of the thinned edges
and use the fixed slot coloring above. All targets now
share the same numerical modulus. For any two fixed rows
and any fixed reduced target ratio, the modular count and
(TR1) give only \(X^{o(L)}\) common columns in that fibre.
The actual triples recovered by (TR1) give the outside
columns because their representations were prescribed.

Counting ordered stars of \(m\) distinct rows over the
columns by the power mean gives at least
\(X^{((m+1)c-m\mu-O_m(\delta)+o(1))L}\) such stars.
Divide by the number of row tuples to find a fixed tuple
with \(X^{(c-m\mu-O_m(\delta)+o(1))L}\) common columns.
Greedily retain columns and discard every other column
sharing a target ratio for any of the finitely many row
pairs with a retained column. Each step discards only
\(X^{o(L)}\) columns by the uniform fibre bound. Choosing
\(\mu,\delta\) sufficiently small proves the assertion.
The argument interchanges the two addends.

If all complete targets on the endpoint outside graph have
a multiplicatively separated form \(C_{ij}=p_i q_j z_{ij}\)
with positive rational factors, and the \(z_{ij}\) on each
row range in a dictionary of normalized exponent at most
\(u\), all target ratios for two fixed rows have normalized
dictionary exponent at most \(2u\). The selected rectangles
therefore require \(u\ge c/2=1/12\). The same test applies
to a specified near-full rectangle within a fixed common
modulus family. In particular exact multiplicative
separation \(C_{ij}=p_i q_j\) cannot supply all the targets
of the endpoint graph.


## Multiplicative separation of the complete addends

At \(d=1/3,\ \Gamma=1/2\), write \(c=1/6\). For any fixed
\(\epsilon>0\), the thinned outside graph contains a subgraph
with at least \(X^{(2c-\epsilon+o(1))L}\) edges on which the
complete integers on both free letters separate multiplicatively
into row and column factors. Here row and column indices are
the outside labels on the free letters. For \(A+B=C\), this
means positive rational factors \(a_i,b_j,a'_i,b'_j\) with
\[
 A_{ij}=a_i b_j,\qquad B_{ij}=a'_i b'_j .
 \tag{MS1}
\]
The third letter need not separate. The edge subgraph need not
be a complete bipartite graph.

Freeze a positive-weight full block on \(C\) at arbitrarily
small root cost \(\mu>0\) as in the [target-ratio argument](#distinct-target-ratios-in-near-full-rectangles).
Retain one prime coloring with a negligible normalized loss.
This time require every prime of each complete integer to have
the color of its letter, and require the positive short slots
on the free letters to have separate subcolors. For each
triple these color requirements are compatible because the
complete values are pairwise coprime and the short slots on a
letter have disjoint prime supports. Each edge imposes
\(o(\log X)\) conditions, giving the same coloring count as
for (TR2).

Write \(R_{ij}=A_{ij}/B_{ij}\). For any four-cycle in this
graph,
\[
 S=\frac{R_{ij}R_{i'j'}}{R_{ij'}R_{i'j}}
 \tag{MS2}
\]
is \(1\) modulo the fixed full block. All outside products
cancel in (MS2). After grouping the atomic windows as in
(TR2), each narrow slot contains a ratio of up to four grouped
roots and the wide factor is a product and quotient of the
four wide corrections. The narrowed slots keep disjoint
supports by the coloring, their heights are \(X^{O(L)}\),
and the wide factor has height \(X^{o(1)}\). The amplified
(GW) count for rational narrow slots in the [proof of (TR3)](#counting-complete-cross-ratios-by-a-full-modulus)
therefore leaves only \(X^{o(L)}\) exact values of \(S\)
over all four-cycles.

If the graph has \(E\) edges and part sizes \(M_A,M_B\),
two applications of Cauchy–Schwarz give at least
\(E^4/(M_A^2M_B^2)\) ordered four-cycle homomorphisms. Indeed
the sum of common-neighbour counts over ordered row pairs is
at least \(E^2/M_B\), and one sums their squares. Those
with a repeated row or column number
\(O(E(M_A+M_B))\), which has smaller exponent for
\(\mu,\delta\) sufficiently small. Averaging over the first
edge gives a fixed anchor \((i_0,j_0)\) with
\[
 X^{(2c-3\mu-O(\delta)+o(1))L}
\]
opposite edges \((i,j)\) for which both cross edges
\((i,j_0),(i_0,j)\) exist. In this calculation the first
step has \(E\ge X^{(2c-\mu-2\delta+o(1))L}\) and
\(M_A,M_B\le X^{(c+\delta+o(1))L}\).
Pigeonhole one value of \(S\) in (MS2) with
\((i',j')=(i_0,j_0)\). The rational odds \(R_{ij}\) then
separate because \(R_{i,j_0}\) depends only on the row,
\(R_{i_0,j}\) only on the column, and \(S,R_{i_0,j_0}\)
are fixed. At every prime of color \(A\) the valuation of
the odds is the valuation of the complete \(A\)-value. At
every prime of color \(B\) its negative is the valuation of
the complete \(B\)-value. Splitting the identity by these
two colors proves (MS1). Taking \(\mu,\delta\) sufficiently
small proves the edge exponent asserted above. The argument
also applies to a signed relation of two positive free
values and a positive target, since their ratio then has a
fixed sign modulo the target block and (MS2) stays \(1\).

## A bound on the nonatomic degree support

Let \(U,V,T\) be any permutation of the letters at the
same minimum nonatomic endpoint. Put
\[
 w_U=\int_0^\infty r\,d\nu_{U,\mathrm{na}}(r),\qquad
 w_V=\int_0^\infty r\,d\nu_{V,\mathrm{na}}(r).
\]
Then
\[
 \operatorname{supp}\nu_{T,\mathrm{na}}
 \ \subseteq\ [0,\;6(2-\max\{w_U,w_V\})].
 \tag{MS3}
\]
The proof uses proportionality of the two common column
parameters in a fixed-row realization of (MS1). It does
not assume multiplicative separation of the third letter.

By symmetry take \(w_U=\max(w_U,w_V)\) and use \(U\) for
the rows. If the assertion fails, choose a compact band
\([s,S]\) carrying positive nonatomic \(T\)-mass \(h\)
with \(s>6(2-w_U)\). Take \(z\) strictly between
\((2-w_U)/s\) and \(c\), and choose \(m\ge2\) with
\(z/m<h\). The small frozen common full block on \(T\)
can use positions disjoint from the subsequent band locks:
if it is chosen inside that band, it can have arbitrarily
small mass and still leave more than \(z/m\) there on every
row. Choose the overlap threshold \(\xi>0\) so small that
\[
 w_U+sz-m(m-1)S\xi>2 .
 \tag{MS4}
\]
Use (MS1) on the two free letters, with its preliminary
\(\mu,\delta\) small compared with \(\xi/m\) and
\((c-z)/m\).

The resulting subgraph has a fixed \(m\)-tuple of rows
with common-column exponent \(c-3m\mu-O_m(\delta)\) and
the internal nonincidences (FG2). In detail, its ordered
star exponent is at least
\((m+1)c-3m\mu-O_m(\delta)\). As in the
[internal star count](#selecting-many-columns), those
with a failed nonincidence are at most
\(X^{((m+1)(c+\delta)-\xi+o(1))L}\), by the
divisor-allocation estimate and pair uniqueness. Dividing
the good stars by the number of row tuples gives the
claimed number of columns. Only the \(T\) outside
products and their overlaps with the other rows' complete
\(T\)-values are needed.

Take the two column parameters \(u_j,v_j\) to be the
complete \(U,V\) values on the first fixed row. By
(MS1), on these rows
\[
 T_{ij}=\varepsilon_U\alpha_i u_j+
        \varepsilon_V\beta_i v_j ,
 \tag{MS5}
\]
with the signs of the additive equation and fixed positive
rational \(\alpha_i,\beta_i\). The coefficients can be
read from the ratios of actual complete values on any one
column. Their numerators and denominators are units at
every \(T\)-colored prime. Each \(u_j\) is divisible by the
same outside product \(P_U\) on the first row, of numerical
weight at least \(w_U-o(1)\). This lower bound follows
from retaining the nonatomic degree mass outside the
shrinking atomic windows and using the first-moment
convergence at \(\Gamma=1/2\).

On each row freeze outside root positions in the band
until their normalized mass lies in
\([z/m,z/m+\delta]+o(1)\). Each position has mass at most
\(\delta+o(1)\), so the preceding choice of \(m\) ensures
these blocks exist. The total freeze costs at most
\(z+m\delta+o(1)\). The recorded blocks determine all
their prime factors. Remove from every block the full
powers of its primes occurring in another block. By (FG2)
the removed normalized root mass per row is at most
\((m-1)\xi+o(1)\). The degrees in the band are at most
\(S/L\), and the bounded-scale roots are squarefree with
disjoint prime supports within a letter, so the numerical
weight lost is at most \(S(m-1)\xi+o(1)\) per row.
The remaining block moduli \(n_i\) are pairwise coprime
and are coprime to \(P_U\) by the letter coloring.

Consider two columns with the same frozen blocks. For
each \(i\), congruence (MS5) holds with zero right-hand
value modulo \(n_i\) for both columns. Eliminate the
\(v\)-term using the two column values \(v_j,v_{j'}\).
Since \(\alpha_i\) is a unit, it follows that
\[
 n_i\mid u_jv_{j'}-u_{j'}v_j .
\]
The determinant is also divisible by \(P_U\), so it is
divisible by \(P_U\prod_i n_i\), of numerical weight
strictly larger than \(2\) by (MS4). Both complete column
parameters are positive and at most \(X\). The determinant
has absolute value at most \(2X^2\) and hence vanishes.
The two parameter pairs then have the same ratio. They
are primitive pairs from the first row, so they coincide,
determine the same triple, and have the same outside
column label.

Thus each fixing gives at most one column, with normalized
count exponent at most \(z+m\delta\). This is strictly
below \(c-3m\mu-O_m(\delta)\) for the stated parameter
choices, a contradiction. This proves (MS3).

### Trimming the heights of the two column parameters

The same determinant permits further locks on the first fixed
row. Write \(s_*=\sup\operatorname{supp}\nu_{T,\mathrm{na}}\).
On \(U\) allow locks on its inside coordinates, disjoint from
the fixed outside product \(P_U\). On \(V\) allow any prescribed
whole-coordinate locks. Let \(p<c\) be their combined limiting
root cost and \(w\) their numerical weight. Then
\[
 w_U+w+s_*(c-p)\le2 .
 \tag{MS9}
\]
These locks include finite sets of selected atomic windows on
either letter. Eligible quantiles on the free coordinates can
also be used whenever they are realized by whole-coordinate
locks with the stated limiting costs and weights.

If the inequality fails, take \(s<s_*\) with
\(w_U+w+s(c-p)>2\), a compact positive-mass \(T\)-band
above \(s\), and \(z<c-p\) with \(w_U+w+sz>2\).
The earlier construction gives \(m\) rows with nearly full
column exponent and pairwise coprime target locks, now of
total mass \(z+m\delta+o(1)\). Also freeze the prescribed
\(U,V\) factors on the first row at cost \(p+o(1)\).
For two columns with the same records, the product of
these latter factors divides the determinant of their
parameters. It is coprime to \(P_U\) and to the target
moduli, by disjoint original positions and the letter
coloring. The determinant now has a fixed divisor of
weight at least \(w_U+w+sz-m(m-1)S\xi+o(1)\), which
exceeds two for small enough \(\xi\). The locks again
determine the column. The cost \(p+z+m\delta+o(1)\)
contradicts the selected column exponent.

In particular, use high atomic windows on \(U\) and the
high tail on \(V\) above a threshold \(s>0\). For every
\(s\le s_*\),
\[
 sc\le
 \int\min(r,s)\,d\nu_{U,\mathrm{at}}(r)
 +\int\min(r,s)\,d\nu_V(r),
 \tag{MS10}
\]
with the analogous inequality on interchanging \(U,V\).
Indeed approximate the positive excess moments
\(\int(r-s)_+\,d\nu_{U,\mathrm{at}}\) and
\(\int(r-s)_+\,d\nu_V\) by a finite set of atoms and a
continuity tail cutoff, respectively. These give the lock
quantity \(w-sp\) in (MS9), with \(s\) in place of \(s_*\).
If a choice has \(p\ge c\), the same inequality for those
locks is automatic from \(w_U+w\le2\). The full moments
are one, and
\(w_U+\int r\,d\nu_{U,\mathrm{at}}+\int r\,d\nu_V=2\),
giving (MS10). The trimmed condition implies (MS3).

### A formal profile excluded by trimming

Give \(U\) an atom of mass \(1/6\) at \(6/5\), and mass
\(1/6\) uniform on \([4,28/5]\). Give \(V\) atoms of
masses \(1/100,47/300\) at \(16,1\), and mass \(1/6\)
uniform on \([3,26/5]\). Give \(T\) an atom of mass
\(1/6\) at \(2\), and mass \(1/6\) uniform on \([1,7]\).
Every total mass is \(1/3\), every nonatomic mass is
\(1/6\) and every first moment is one. The nonatomic
weights are \(4/5,41/60,2/3\), respectively. Thus all
three untrimmed bounds (MS3) hold.

The pooled mass at scaled degrees at least four exceeds
\(1/3\). Its numerical weight is
\[
 \frac4{25}+\frac45+\frac{23}{55}+\frac{11}{24}<2.
\]
Hence even the fractional top \(1/3\) has weight below two,
and the three-factor packing tests hold. But choosing the
\(V\)-atom at \(16\) in (MS9) gives
\[
 \frac45+\frac4{25}+7\left(\frac16-\frac1{100}\right)
   =\frac{617}{300}>2.
\]
The profile therefore violates the trimmed bound.

## Power dictionaries after separating the free letters

The multiplicative separation removes the outside-to-complete
multiplicity terms from some of the dictionary restrictions.
Take \(m\) fixed rows in the [separated subgraph](#multiplicative-separation-of-the-complete-addends)
with column exponent tending to \(c=1/6\), as provided by the
star count. Their complete targets are distinct linear forms
\(\ell_i(u,v)\) of the two primitive complete integers on the
first row. More precisely, the forms are pairwise nonproportional.
Proportional forms would give two rows with equal odds of the
free letters on every common column. On one column the
primitive fractions then determine identical triples, contrary
to distinct outside rows. Thus \(T_1/T_2\) determines \(u/v\),
and this primitive ratio determines the column on its own.

### Several rows with a common power

Suppose that each \(T_i\) in such a rectangle has a representation
\(T_i=W_i x_i^q\) with positive integer \(W_i,x_i\) and a common \(q\ge2\). For fixed coefficients
on \(m\ge3\) rows, the number of columns is
\[
 O_{m,q,\varepsilon}\!\left(X^{2/q^{m-1}+\varepsilon}\right).
 \tag{MS6}
\]
Indeed \(\ell_1,\ell_2\) span the pencil of linear forms.
Putting \(t=T_2/T_1\), the ratios \(T_i/T_1\) are affine
linear functions of \(t\) with distinct simple zeros for
\(i\ge2\). The projective point \([x_1:\cdots:x_m]\) lies on
the curve obtained by pulling this pencil back under
coordinatewise \(q\)-th powers. In the chart \(x_1\ne0\),
the coordinates \(x_i/x_1\) for \(i\ge2\) are Kummer roots
of the affine forms in \(t\), with the fixed \(W_i\)
included in the coefficients. Successive Eisenstein arguments
at their distinct zeros give a geometrically integral curve,
as in the [mixed-power fibre product](#the-curve-for-several-fixed-complete-rows).
The map to \([x_1:x_2]\) has degree \(q^{m-2}\).
Its centre \(x_1=x_2=0\) misses the curve, since the
remaining coordinates must all vanish there by the equations
\(W_i x_i^q=a_i W_1 x_1^q+b_i W_2 x_2^q\).
A generic hyperplane in this pencil has exactly \(q^{m-2}\)
distinct points, giving the same projective degree. The
actual points have projective height at most \(X^{1/q}\),
and each point determines the column by its target ratio.
The uniform projective curve estimate therefore gives (MS6).
[Binyamini–Cluckers–Kato, Theorem 1](https://arxiv.org/pdf/2401.03982v2#page=2)

In particular, if the \(W_i\) range in a dictionary of
normalized exponent \(u_0\), saturation of the rectangle
requires
\[
 m u_0+\limsup_{L\to0}\frac{2}{Lq^{m-1}}\ge\frac16 .
 \tag{MS7}
\]
All the degree-dependent constants are taken at fixed
truncation before the height limit. This applies to a
common-power dictionary of the original endpoint family,
because the separated subgraph need not change its
representations of complete target values. If the full family
initially allows varying \(q\ge Q\), one can first retain one
exponent at negligible normalized cost. A nontrivial root
gives \(q\le\log_2 X\), and for a root equal to one one may
choose the minimum eligible integer \(q\) without changing
the coefficient. Thus a global dictionary with \(q\ge Q\)
also needs \(m u_0+\limsup 2/(LQ^{m-1})\ge c\).
The heterogeneous estimates below apply to a fixed set of
row-wise power choices.

### A weighted projective bound for different powers

For fixed choices \(T_i=W_i x_i^{q_i}\) on \(m\ge3\) rows,
let
\[
 M=\operatorname{lcm}(q_1,\ldots,q_m),\qquad
 g=\gcd(q_1,\ldots,q_m),\qquad
 D=\frac{\prod_i q_i}{M}.
\]
Then the column count is \(O_{\mathbf q,\varepsilon}(X^{2/D+\varepsilon})\).
This recovers (MS6) when the powers are equal. It also gives
\(O_{Q,\varepsilon}(X^{1/(3Q^2)+\varepsilon})\) for
the three powers \(2Q,3Q,6Q\). For the four powers
\(abc,abd,acd,bcd\) with pairwise coprime \(a,b,c,d\),
the gcd is one and \(D=(abcd)^2\).

Set \(w_i=M/q_i\), so the weights \(w_i\) have gcd one, and
use as projective coordinates all monomials in the \(x_i\)
of weighted degree \(W=M/g\). Each is an integer at most
\(X^{W/M}=X^{1/g}\), and the pure coordinates
\(Y_i=x_i^{q_i/g}\) are among them. Work over \(\mathbb C\)
on the chart \(x_1\ne0\), scaling \(x_1\) to one by the
weighted action. Over \(t=T_2/T_1\) the resulting \(m-1\)
roots satisfy monic Kummer equations with distinct simple
zeros, exactly as in (MS6). Their field has degree
\(\prod_{i\ge2}q_i\) over \(\mathbb C(t)\).
The residual group \(\mu_{w_1}\) acts by \(x_i\mapsto
\zeta^{w_i}x_i\), and acts faithfully since \(\gcd_i w_i=1\).

The ratio field of the monomial image is the fixed field
of this residual action. To verify generation, the primitive
weight relation for any pair \(i,j\) is a difference of
monomials of weighted degree \(\operatorname{lcm}(w_i,w_j)\).
Pad the two monomials by the same power of \(x_i\) to
reach degree \(W\), which is a common multiple of all
the weights. The pair relations generate the lattice
\(\{(e_i)\in\mathbb Z^m:\sum_i w_i e_i=0\}\). At each
prime, choose an index with \(p\)-adic unit weight.
The pair relations from that index to all the others
span the lattice after localization at \(p\), since their
other coordinates are \(p\)-adic units. This gives full
index one over \(\mathbb Z\). Setting \(x_1=1\), the
resulting monomial ratios generate exactly the characters
fixed by \(\mu_{w_1}\). Together with \(t\), which is
a fixed multiple of \((Y_2/Y_1)^g\), they generate the
fixed field. Its degree over \(\mathbb C(t)\) is \(D\).

The map of the image curve to \([Y_1:Y_2]\) has degree
\(D/g\). This linear projection is defined on the projective
closure. If \(Y_1=Y_2=0\), all \(Y_i\) vanish by the
relations expressing \(W_iY_i^g\) in the pencil of the
first two. For any monomial coordinate \(Z_{\mathbf e}\)
there is a homogeneous relation
\(Z_{\mathbf e}^W=\prod_i Y_i^{e_i w_i}\), so all
coordinates would vanish, which is impossible. A generic
hyperplane of this pencil therefore has \(D/g\) points,
giving that degree for the geometrically integral image
curve. Binyamini–Cluckers–Kato's projective estimate at
height \(X^{1/g}\) yields the stated exponent. The pure
coordinates \(Y_1,Y_2\) give the target ratio \(t\), so
each projective point still determines the column. In
using dictionaries one can take the smaller of this
exponent and any bound for the same degree choices.

### Three rows with different large powers

Suppose instead that \(T_i=W_i x_i^{q_i}\) on three rows,
with positive integer \(W_i,x_i\) and arbitrary \(q_i\ge Q\ge2\). For fixed choices of
\((W_i,q_i)\) the count is
\[
 O_{\mathbf q,\varepsilon}
  \!\left(X^{3/Q^{3/2}+\varepsilon}\right).
 \tag{MS8}
\]
The three pairwise nonproportional linear forms have a
linear relation with three nonzero coefficients. First
remove \(g=\gcd(x_1,x_2,x_3)\). Each prime of \(g\) divides
both \(T_1,T_2\). It therefore divides the numerator of
the fixed nonzero determinant of \(\ell_1,\ell_2\),
because their denominators are units at target primes
and \(\gcd(u,v)=1\). With the letter coloring this follows
directly. One can also include the cleared denominators
of the forms in this fixed polynomial-height integer.

For any fixed polynomial-height nonzero integer \(F\),
the number of \(g\le X\) supported on primes dividing
\(F\) is \(X^\varepsilon\) for every \(\varepsilon>0\),
uniformly in \(F\). To see this use Rankin's bound
\(X^\eta\prod_{p\mid F}(1-p^{-\eta})^{-1}\) for fixed
\(0<\eta<1\). The logarithm of the product is
\(O_\eta(1+\omega(F)^{1-\eta})=o(\log X)\) by ordering
the primes. Thus fixing \(g\) has arbitrarily small
exponent cost. The divided \(x_i/g\) form a primitive
integer triple with heights at most \(X^{1/Q}\) and
satisfy a generalized Fermat equation with fixed
nonzero coefficients, including the powers of \(g\).
Browning–Verzobio's coefficient-uniform bound for
primitive triples has box exponent at most
\(1/\sqrt{q_1}+1/\sqrt{q_2}+1/\sqrt{q_3}\) in a cube.
Each recovered triple restores \(T_1/T_2\), so (MS8)
follows. [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

Consequently a dictionary of \((W_i,q_i)\) of normalized
exponent \(u_0\), with \(q_i\ge Q\) at fixed truncation,
requires \(3u_0+\limsup 3/(LQ^{3/2})\ge1/6\).
In particular, if \(LQ^{3/2}\to\infty\), the coefficient
and exponent choices must have exponent at least \(1/18\).
The complete rows and columns are allowed to vary over
their outside labels throughout.

## Complete values determined by single outside labels

At \(d=1/3,\ \Gamma=1/2\), for every fixed \(\epsilon>0\)
one can retain \(X^{(1/3-\epsilon+o(1))L}\) triples on which
each complete value is a function of the outside label on
its own letter. In particular the near-full outside graph
has a realization by complete additive values. This uses
(MS1) successively in all three orientations and a bound
for the multiplicative energy of the factor dictionaries.

### Freezing factors related on every edge

Suppose the three parts of a linear hypergraph have
size at most \(X^{(c+\delta+o(1))L}\), its edge count has
exponent at least \(2c-\eta\), and each vertex in part \(i\)
has a positive rational factor \(f_i\) of polynomial height.
Assume \(f_1 f_2=f_3\) on all edges. There is a subfamily of
exponent at least \(2c-5\eta-8\delta+o(1)\) on which all
three factors are constant.

For each part, classify labels by the size of the fibre of
their factor, using dyadic bins for that size. Pigeonhole
one bin on each part over edges, losing only a power of
\(\log X\) because the height of each part is polynomial.
If fibre sizes on part \(i\) are in \([K_i,2K_i)\), put
\[
 v_i=c+\delta-\frac{\log K_i}{L\log X}+o(1).
\]
Take the limiting upper-capacity error so that \(v_i\ge0\).
The factor dictionary on that part has exponent at most
\(v_i\). Pair uniqueness of the original hypergraph implies
that any fixed two factor values occur together on at most
\(4K_iK_j\) edges. The number of matched factor pairs in
these bins consequently has exponent at least
\(v_i+v_j-\eta-2\delta+o(1)\).

For each part separately color primes by numerator or
denominator sign in that part's factor. Averaging over
the three independent colorings retains a fraction
\(X^{-o_X(1)}\) of the edges, because each factor has
polynomial height and hence \(o(\log X)\) prime divisors.
Apply the preceding pair lower bound to the retained
edges. Within each dictionary \(S_i\), products have no
cancellation, so each fixed product has divisor-many
factor representations. Its multiplicative energy is
therefore at most \(|S_i|^2 X^{o_X(1)}\). For two different
dictionaries the energy of their products, or of their
quotients, is bounded by the geometric mean of the two
single-dictionary energies. This follows by writing the
mixed energy as a sum over equal quotients and using
Cauchy–Schwarz.

Under \(f_1 f_2=f_3\), a matched pair on any two parts
gives the value on the third by product or quotient.
Cauchy–Schwarz on the matched pairs now gives at least
\(X^{(v_i+v_j-2\eta-4\delta+o(1))L}\) distinct output values.
Thus \(v_k\ge v_i+v_j-2\eta-4\delta+o(1)\) for each
orientation. Adding the two inequalities with outputs
different from \(k\) gives \(v_k\le 2\eta+4\delta+o(1)\).
Pigeonhole values of the factors on any two parts, at
total cost at most \(4\eta+8\delta+o(1)\). The third factor
is determined by the relation, proving the stated
retained exponent.

### Applying separation on all three pairs

The proof of (MS1) works on any already selected subgraph
whose edge exponent is arbitrarily close to \(2c\). In
the anchor estimate it multiplies the edge deficit by
three and adds the common-block freeze and the negligible
coloring and congruence losses. The common full block is
available on any choice of the third letter, and it can
have arbitrarily small normalized cost. The three pair
orientations can therefore be applied successively, keeping
the edge exponent arbitrarily close to \(2c\). All
factorizations obtained earlier remain valid on the next
subgraph.

For example, let \(i,j,k\) denote outside labels of
\(A,B,C\). Separation in the \((A,B)\) and \((A,C)\)
orientations writes
\[
 A=a_i b_j=\widetilde a_i d_k .
\]
Its factors are positive rationals of polynomial height:
they can be taken from the actual anchor terms and from
the valuation projections of their cross ratios in (MS2).
The exact three-vertex relation is
\((\widetilde a_i/a_i)d_k=b_j\).
Apply the factor-freezing argument to make \(b_j\)
constant on a still near-full edge set. There \(A\)
depends only on \(i\). Repeat for \(B\) and \(C\) using
their two factorizations. The number of steps is fixed,
so the initial losses can be chosen to meet any prescribed
positive final loss. This proves the claim. The
factorization assignments need only hold on the retained
edges throughout.

### Consequence for all three nonatomic supports

In the retained family, a fixed complete value permits
only divisor-many outside labels at fixed truncation.
Thus the graph can also be thinned by a subpower factor
to make the retained label-to-complete-value maps
injective. The ordered star count then gives, for any
fixed \(m\), \(m\) outside rows with distinct fixed
complete row values and a common column exponent
arbitrarily close to \(c\). Each column has a single
complete value on the other free letter.

Apply the [fixed-row Chinese remainder bound](#combining-congruences-from-several-fixed-rows)
to any choice of target letter. Its affine slopes are
\(1\) or \(-1\), according to the signs of the equation,
and its zeros are distinct because the complete rows
are distinct. Positive diffuse mass on a compact band
above scaled degree six would give a strict column
saving, contradicting the preceding star count.
Consequently every letter at the minimum nonatomic
endpoint satisfies
\[
 \operatorname{supp}\nu_{T,\mathrm{na}}\subseteq[0,6].
 \tag{MS11}
\]
The trimmed inequalities still apply with any positive
band below the top of this support and may use additional
information about the individual degree profiles.

The one-parameter argument also permits high-degree locks on the
common free column. At the minimum endpoint, for any \(s>0\)
with positive nonatomic \(T\)-mass above \(s\), and any
opposite letter \(U\), one necessarily has
\[
 \frac{s}{6}\le\int\min(r,s)\,d\nu_U(r).
 \tag{MS12}
\]
To see this, work in the same complete additive rectangle,
with complete \(U\) as the common parameter. Freeze prescribed
whole \(U\)-factors of mass \(p\) and weight \(w\). Their complete
product \(F\) divides the parameter, and the remaining integer
parameter is at most \(X^{1-w+o(1)}\). On sufficiently many rows
freeze disjoint high-degree \(T\)-blocks of total mass \(z\)
as in (CT1), with \(p+z<c\). Their stripped congruences
determine the reduced parameter if \(sz>1-w\), using small
overlap and window losses. Division by \(F\) is allowed at
all retained \(T\) primes by primitivity, since \(F\) is the
same for every column sharing the freeze. The resulting
contradiction proves \(w+s(c-p)\le1\); when \(p\ge c\) the
same inequality follows from \(w\le1\). Take the \(U\) tail
above \(s\), approximating its excess moment by continuity
cutoffs, to get (MS12). In particular an atom above six
on \(U\) forces strictly smaller nonatomic support on every
other letter.

### Two-letter packing with the third value fixed

On the complete-label graph, prescribed blocks on any two
letters with combined root mass \(p<c\) have numerical weight
at most one. Indeed there are at most \(X^{(c+\delta+o(1))L}\)
complete values on the third letter. For each one, fix the
values of the two blocks. Their products are coprime on every
primitive witness. The congruence for one remaining complete
term modulo its own block is zero; modulo the other block
it is a fixed signed multiple of the third value. Their
product determines that term uniquely if its weight exceeds
one. This would give exponent \(c+\delta+p+o(1)\), strictly
below the retained graph's exponent. In particular one can
select full atomic windows or take eligible nonatomic
quantiles on these two letters, with total mass strictly below
\(c\). Continuity gives the respective boundary quantile
conditions at mass \(c\).

The lock choices may even depend on the edge if divisors of
the original roots are used. Fix \(\varepsilon>0\). Among the
complete-label edges, those with factors \(d_{W,e}\mid x_{W,e}\)
on a pair of letters satisfying
\[
 \prod_{W,e} d_{W,e}\le X^{(c-\varepsilon)L},
 \qquad
 \prod_{W,e} d_{W,e}^{e}>X^{1+\varepsilon}
 \tag{MS13}
\]
have exponent at most \(2c-\varepsilon+\delta+o(1)\).
To count them, sum over the dyadic sizes of all selected
divisors before enumerating their values. At fixed truncation
there are only a power of \(\log X\) many box choices, and
the product of upper divisor sizes is at most
\(X^{(c-\varepsilon+o(1))L}\). With the third complete
value fixed the same Chinese remainder argument gives at
most one edge for every record. This proof never requires
counting the completions of the individual recorded roots.

For example, consider the pooled weak measure of the two
letters restricted to a bounded band of scaled degrees at
most \(R\), with continuity endpoints. If its fractional
upper quantile at a fixed \(t<c\) is greater than \(1+\gamma\)
for some \(\gamma>0\), then all but an exponentially small
fraction of the near-full complete-label edges have a prime
in this band with normalized logarithmic size at least
\(\gamma/(4R)\). The bounded-scale roots are squarefree
outside the packed index, so the prelimit measure can be
split into its individual prime contributions without
changing the fractional quantile. Rank those contributions
by scaled degree and round down at the last prime. If all
prime increments are smaller than \(\gamma/(4R)\), the
numerical weight lost is less than \(\gamma/4\), and the
root mass is at most \(t\). The quantile convergence on
the bounded band makes this an instance of (MS13) with a
fixed positive margin. Thus a fractional packing violation
forces macroscopic primes even in degree windows whose
roots were not assumed prime.

## Quantitative trimming near the minimum nonatomic total

Consider any balanced limiting endpoint with \(d=1/3\), total root mass \(1/3\) on each letter, and
\[
 c=\frac16,\qquad \gamma=\Gamma-\frac12,\qquad 0\le\gamma\le c/20 .
\]
For distinct letters \(U,V,T\), let
\[
 \Delta=c_U+c_V-\frac13,\qquad
 w_U=\int r\,d\nu_{U,\mathrm{na}}(r).
\]
The endpoint pair bounds give \(0\le\Delta\le2\gamma\) and \(c-\gamma\le c_W\le c+\gamma\) on every letter. The signs of the three terms may be any fixed signs of a primitive abc relation. The normalization makes all complete magnitudes at most \(X^{1+o(1)}\).

Fix a compact positive degree band \([s,S]\) on \(T\), with nonatomic mass \(H>0\), and an integer \(m\ge2\). Allow prescribed locks on \(U\) disjoint from its outside product, and arbitrary locks on \(V\), of combined normalized root cost \(p\) and numerical weight \(w\). For any \(t>0\) such that
\[
 t<mH,\qquad p+t<c_V-3m\Delta ,
\]
one has
\[
 \boxed{\quad
   w_U+w+st\le
   2+m(m-1)S\,\max\{0,3m\Delta+c_T-c_U\}.
   \quad}
 \tag{NT1}
\]
In particular the right-side error is at most \(m(m-1)S(6m+2)\gamma\). The band mass is nonatomic so no large-coordinate rounding penalty is needed inside the band. At exact minimum mass this specializes to the finite-band version of the [trimmed support bound](#trimming-the-heights-of-the-two-column-parameters).

### Rectangles with a controlled loss away from the minimum

The [pair-fibre count](#a-uniform-pair-fibre-bound) works with these masses also. Fix enough positive atom windows plus the zero windows on each letter that the omitted atomic mass is at most \(\delta\). Then the outside part sizes have exponents at most \(c_W+\delta\). With any two outside labels fixed, their retained atomic factors are finitely many grouped narrow roots with small numerical cofactor. The remaining letter has positive nonatomic mass \(c_W\ge c-\gamma\), and so supplies a positive full modulus at arbitrarily small root cost. The grouped full-modulus count gives a pair-fibre exponent as small as desired. Thin to a linear outside graph, with count exponent \(d-o_\delta(1)\). This includes the outside representation and the finite-atom-window errors.

As in [multiplicative separation](#multiplicative-separation-of-the-complete-addends), freeze a common full block on \(T\) at positive but arbitrarily small mass cost \(\mu\). It may be chosen disjoint from the eventual high band locks within each target; if chosen from that band, take its mass small compared to the available band slack. Color primes globally by letter and by their narrow slots, at negligible normalized cost. The four-cycle cross ratios of the complete free ratios have only \(X^{o(L)}\) values by the anchored odds count, since all outside products cancel and all narrow slots and the wide cofactor have the same bounds. The ordered four-cycle homomorphism lower bound is \(E^4/(M_U^2M_V^2)\) for \(E\) edges and part sizes \(M_U,M_V\). Repeated rows or columns are negligible: their normalized exponent is at most \(3c+3\gamma+O(\delta)\), whereas the cycle lower bound is at least \(4c-4\gamma-4\mu-O(\delta)\). The hypothesis \(\gamma\le c/20\) leaves a fixed margin. Average over an anchored edge, losing one factor of \(E\), and fix one cross ratio. The resulting separated edge exponent is at least
\[
 3d-2(c_U+c_V)-3\mu-O(\delta)
    =d-2\Delta-3\mu-O(\delta).
 \tag{NT2}
\]
On these edges the complete free values separate multiplicatively into row and column rational factors, by the global letter coloring and the exact odds identity.

Count ordered \(m\)-stars above outside \(V\)-columns. They have exponent at least
\[
 m(d-2\Delta-3\mu)-(m-1)c_V-O_m(\delta).
\]
The mean column degree grows with positive exponent, since \(c_U-3\Delta\ge c-7\gamma>0\), and \(\mu,\delta\) can be small. Thus repeated rows are negligible. Ignore only stars in which a \(T\) outside label has radical overlap at least \(X^{\xi L}\) with the complete \(T\) on another row of the same column. By the [divisor allocation estimate](#counting-shared-prime-divisors), the \(T\) label has exponent at most \(c_T+\delta-\xi\) after the offending complete term is fixed. Its pair with the common column determines its row by thinning. Counting the known complete edge and the remaining \(m-2\) rows shows the bad stars have exponent at most
\[
 (m-1)c_U+c_V+c_T-\xi+O_m(\delta).
\]
Subtracting the lower star exponent shows that it is sufficient to take
\(\xi>3m\Delta+c_T-c_U\), with \(\mu,\delta\) then arbitrarily small. Divide the good stars by the capacity of the row tuples. There are \(m\) fixed rows with internal overlaps below \(X^{\xi L}\) on at least
\[
 X^{(c_V-3m\Delta-3m\mu-O_m(\delta)+o(1))L}
\]
columns, all with the separated complete odds. No control on prime overlap between different columns is used.

### Determinant of two columns

On row one, write its two complete free values \(u_j,v_j\) as the column parameters. The fixed-row targets have the form
\[
 T_{ij}=\epsilon_U a_i u_j+\epsilon_V b_i v_j
\]
with fixed positive rational coefficients, which are units at every \(T\)-colored prime. The fixed outside \(U\)-product \(P_U\) divides all \(u_j\), with numerical weight at least \(w_U-o(1)\), by retaining the entire nonatomic part outside the shrinking atomic windows. This lower bound needs only weak convergence, so any residual escaped moment is harmless.

Lock the prescribed \(U,V\) factors of cost \(p\) and weight \(w\) on row one. On each of the \(m\) rows, choose an outside \(T\)-block in \([s,S]\) of mass \(t/m+o(1)+O(\delta)\), rounding at the small outside-coordinate cap. These blocks exist since \(t/m<H\). For each block strip all its primes that occur in another fixed block. The internal overlap bound loses at most \((m-1)\xi+o(1)\) root mass per row. Bounded-degree outside roots are squarefree, so the numerical weight loss is at most \(m(m-1)S\xi+o(1)\) in total. The recorded blocks determine the intersections and hence the stripping. The remaining locks on \(T\) are pairwise coprime. The common target divisor used in (NT2) is not included in the determinant weight and need not be counted a second time.

The choices for all these locks have exponent at most \(p+t+O_m(\delta)\), strictly below the column exponent by hypothesis and a sufficiently small choice of \(\mu,\delta\). Thus at least two distinct columns share the records. The determinant \(u_jv_l-u_lv_j\) is nonzero, since equality of complete primitive ratios on row one would identify that triple and its prescribed outside column label. It is divisible by the fixed outside \(U\) factor and the fixed selected \(U,V\) factors, using disjoint original coordinates and the letter coloring. It is also divisible by each stripped target block, by eliminating one of \(u,v\) from the two target congruences for that row. These target blocks are coprime to the free-letter factors by coloring. Therefore the product of weight at least
\[
 w_U+w+st-m(m-1)S\xi+o(1)
\]
divides this nonzero determinant, whose size is \(O(X^{2+o(1)})\). Letting \(\xi\) decrease to \(\max(0,3m\Delta+c_T-c_U)\), and choosing \(\mu,\delta\) afterward for each fixed slack, gives (NT1).

### An explicit separation from \(\Gamma=1/2\) in a broad spectral sector

Take a formal profile of mass \(1/3\) on each letter, all heights one, and \(\Gamma=1/2\).

- \(U\): an atom of mass \(1/6\) at \(6/5\) and nonatomic mass \(1/6\) uniformly on \([4,28/5]\).
- \(V\): atoms of masses \(1/100,47/300\) at \(16,1\) and nonatomic mass \(1/6\) uniformly on \([3,26/5]\).
- \(T\): an atom of mass \(1/6\) at \(41/20\), nonatomic mass \(1/12\) uniformly on \([0,2]\), and nonatomic mass \(1/12\) uniformly on \([34/5,7]\).

Their nonatomic numerical weights are \(4/5,\ 41/60,\ 79/120\), respectively. All three bare support caps at minimum mass hold. Even the fractional top total mass \(1/3\) has weight below two: the high \(V\) atom and all degrees at least four already have mass exceeding \(1/3\) and their entire numerical weight is
\[
 \frac4{25}+\frac{23}{40}+\frac45+\frac{23}{55}<2 .
\]
Thus the full three-factor packing condition holds.

In (NT1) use \(m=2\), the \(T\) band \([34/5,7]\), \(t=31/200\), and the high \(V\) atom as the prescribed lock \(p=1/100,\ w=4/25\). Their left side is
\[
 \frac45+\frac4{25}+\frac{34}{5}\frac{31}{200}
 =2+\frac7{500}.
\]
As long as \(\gamma\le1/20000\), the budget right side \(c_V-6\Delta\) is at least \(1/6-13/20000>33/200=p+t\), while twice the band mass exceeds \(t\). The quantitative error on the right is at most \(196\gamma\le49/5000<7/500\). Thus the same strict spectral inequalities rule out balanced endpoints in this sector whenever \(\Gamma\le1/2+1/20000\). Small perturbations of the stated weights and interval masses preserve all the strict comparisons. The strict inequalities define the spectral sector for this separation.

### A count-exponent saving in the same sector

The proof also applies to limiting box profiles with \(d<1/3\). In that setting put
\[
 \Delta=c_U+c_V-d
\]
and assume \(d=1/3-\varepsilon\) and \(|c_W-1/6|\le\gamma_0\) on each of the three letters, where \(4\varepsilon+7\gamma_0<1/6\). The usual pair-fibre proof gives \(d\le c_U+c_V\), since the remaining letter has positive nonatomic mass and supplies arbitrarily cheap positive full moduli. The anchored exponent in (NT2) is still \(3d-2(c_U+c_V)=d-2\Delta\). The four-cycle lower bound has exponent at least \(4c-4\varepsilon-4\gamma_0-o(1)\), while repeated rows or columns cost exponent at most \(3c+3\gamma_0+o(1)\), giving the required positive margin. The selected graph's mean column degree has exponent at least \(c_U-3\Delta\ge c-3\varepsilon-7\gamma_0>0\). Thus the stated star count and overlap subtraction apply without change. In particular (NT1) holds if \(p+t<c_V-3m\Delta\) and \(t<mH\), with this new \(\Delta\).

For an explicit sector, require \(|c_W-1/6|\le10^{-5}\), a \(T\) nonatomic band \([34/5,7]\) of mass greater than \(31/400\), and prescribed first-row locks of root cost at most \(101/10000\) and combined numerical weight \(w\) such that \(w_U+w\ge959/1000\). Any such limiting profile satisfies
\[
 d<\frac13-\frac1{10000}.
\]
Indeed suppose \(d\ge1/3-1/10000\) and take \(m=2,\ t=31/200\). The case \(d>1/3\) is irrelevant by the global limiting small-budget bound. Now \(\Delta\le.00012\), so \(c_V-6\Delta\ge.1659366\ldots>.1651\ge p+t\). Also
\[
  m(m-1)S\max\{0,3m\Delta+c_T-c_U\}\le .01036 .
\]
But the left side of (NT1) is at least \(.959+(34/5)(31/200)=2.013\), a contradiction. This supplies a normalized count-exponent saving for the sector, including for profiles whose count exponent falls strictly below the balanced endpoint. The formal example lies strictly inside its spectral weight and mass hypotheses.

## Positive nonatomic numerical weight on every letter

There are absolute \(\gamma_0,\beta_0>0\) such that at a
balanced endpoint of exponent \(d=1/3\) with
\(0\le\Gamma-1/2\le\gamma_0\),
\[
 \int r\,d\nu_{T,\mathrm{na}}(r)>\beta_0
 \tag{DN1}
\]
on each letter \(T\). The count fixes a complete opposite
value and uses a narrow degree interval on \(T\).
The full modulus comes from the third letter, whose
nonatomic numerical weight is controlled by (VM7).

### The complete-label reduction close to minimum mass

Put \(\gamma=\Gamma-1/2\) and \(c=1/6\).
The pair inequalities give \(c-\gamma\le c_W\le c+\gamma\)
on each letter. In particular all three letters provide a
positive-weight full modulus at arbitrarily small root cost.
The [outside pair-fibre bound](#a-uniform-pair-fibre-bound)
then gives a thinned outside graph of exponent arbitrarily
close to \(2c\), with capacities \(c+\gamma+\delta\).
For any fixed tolerance \(\epsilon>0\), if \(\gamma\)
is sufficiently small in terms of \(\epsilon\), one can
retain a subgraph of exponent at least \(2c-\epsilon\)
whose complete values are functions of their own labels.
The retained complete part sizes have exponents at most
\(c+\gamma+\delta\).

Indeed an application of (MS1) to a subgraph of exponent
\(2c-\lambda\) has new edge deficit at most
\(3\lambda+4\gamma+3\mu+O(\delta)\), with \(\mu\)
the cost of the newly frozen common divisor. This follows
from the same four-cycle count with both part capacities
at most \(c+\gamma+\delta\). For the required small
parameters, the repeated-row and repeated-column counts
remain negligible. The subsequent factor-freezing step
in the [complete-label proof](#freezing-factors-related-on-every-edge)
works with \(\gamma+\delta\) in place of \(\delta\), giving
new deficit at most \(5\lambda+8(\gamma+\delta)\).
There are only three operations of each kind. Choose
\(\gamma,\delta,\mu\) sufficiently small for their final
deficit to be less than \(\epsilon\). This argument uses
the upper capacities and positive nonatomic masses, so
their exact equality is unnecessary.

### A one-root count for a concentrated numerical moment

Fix a letter \(T\). The [positive weight of both types on
the opposite pair](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#atomic-and-nonatomic-weight-near-the-minimum-mass)
gives one of its two opposite letters, say \(V\), with
uniformly positive nonatomic mass in a fixed positive
compact range \([a,b]\). Denote this mass lower bound by
\(\zeta>0\), and call the remaining letter \(U\).
Use the complete-label graph just proved, so \(U\) can
be fixed at cost \(c+\gamma+\delta\).

The [height and tail bounds](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#locking-a-high-degree-tail)
at near-minimum \(\Gamma\) ensure that for a fixed \(R\) and
sufficiently small \(\gamma\) the \(T\)-degree measure
on some interval \([1/2,R]\) has numerical weight at
least \(1/3\). In fact its full numerical height tends
to one as \(\gamma\to0\), the numerical weight above
a fixed cutoff from (VM2) is \(O(\gamma)\), including
escape, and the weight below \(1/2\) is at most \(1/6\).
Use a uniform small width \(\sigma>0\), to be chosen
below. Cover this fixed range by a bounded number of
intervals of width at most \(\sigma\) with continuity
endpoints, allowing slightly larger enclosing endpoints.
At least one interval \(I\), bounded away from zero,
has mass \(t\ge t_*>0\) for a fixed \(t_*\) depending
only on \(R,\sigma\).

On \(V\) freeze an outside block in \([a,b]\) with root
mass \(\mu\in[\lambda t,2\lambda t]+o(1)\), where
\(\lambda>0\) is a sufficiently small absolute constant.
This is possible by the diffuse mass \(\zeta\), taking
\(\lambda\le\zeta/4\) and outside-coordinate rounding
errors small compared with \(\lambda t_*\). Its numerical
weight \(w\) lies between \(a\lambda t\) and
\(2b\lambda t\) in the limit, and use fullness
\(kL\to\kappa=a/2\). Choose \(\lambda\) also so that
\[
 \mu\le t/16,\qquad 2w/\kappa\le t/16 .
\]
For example one can take
\(\lambda\le\min\{\zeta/4,1/64,a/(128b)\}\).
The selected \(V\)-modulus divides the other free
letter, so with complete \(U\) fixed, \(T\) is in a
fixed unit residue class.

On \(T\) leave the \(I\)-coordinates and the degrees
below a small cutoff \(\tau\) free. Freeze the rest.
Write \(B_T=\int r\,d\nu_{T,\mathrm{na}}\). The cost of
freezing is at most
\[
 c+\gamma+B_T/\tau-t,
\]
since the total atomic mass is at most \(c+\gamma\)
and the nonatomic mass above \(\tau\) is at most
\(B_T/\tau\). Choose \(\tau\) a continuity cutoff in
\([\tau_0/2,\tau_0]\), with
\(0<\tau_0<\min(1/4,\sigma t_*)\), staying below \(I\).
Group the \(I\)-roots at their midpoint exponent and
put the corrections and the low roots into the wider
factor. The narrow root has normalized height at most
\(t\), while the wider factor has numerical height at
most \(\sigma t+\tau/3<2\sigma t\) in the limit.
The single-root estimate (UW) applies with narrow
height constant \(C_1=1/\lambda\) relative to \(n^{1/k}\).
In fact \(t\kappa/w\le1/(2\lambda)\), and the product
of \(C_1\) with the wider exponent relative to \(n\)
is at most \(2\sigma/(a\lambda^2)\). Choose \(\sigma\)
so small that this is strictly below the dimension-one
constant in (UW). The choices of \(t_*,\tau_0\) are
made after \(\sigma\) is fixed. Taking one tuple for
each exact grouped value gives exact fibre bound one,
as the complete \(T\) and fixed \(U\) determine the
primitive triple. The modulus count therefore costs
at most \(2w/\kappa\) in the normalized exponent.

If \(B_T\le\beta_0=\tau_0 t_*/32\), then
\(B_T/\tau\le t/16\). Adding the fixed \(U\) and all
the locks and congruence classes gives exponent at
most \(2c+2\gamma+\delta-13t/16+o(1)\).
For sufficiently small \(\gamma,\delta\) this is
strictly below the exponent \(2c-t_*/8\) obtainable
from the complete-label graph. This contradiction
proves (DN1). The constants are uniform in the choice
of the three letters. In particular the lower bound
does not depend on how many positive atom scales
occur on the letter being counted.

## Limits of the graph condition

The pair-entropy identities alone are compatible with exact saturation. In a finite abelian group of size \(n\), let \(U,V\) be independent uniform elements and put \(W=U+V\). Each pair is independent and determines the third. Taking a direct product and doing addition coordinatewise also satisfies the corresponding sub-coordinate identities. With \(\log n\) at the single-letter capacity scale, these are the equalities suggested by (EN3). The arithmetic relation between the actual root vectors is needed to eliminate this possibility.

A graph in which every edge belongs to a unique triangle can still
have \(n^{2-o(1)}\) triangles on \(O(n)\) vertices. For the
Ruzsa–Szemerédi construction, take a set \(D\subseteq[1,n]\)
with no nontrivial three-term arithmetic progression and
\(|D|=n^{1-o(1)}\). On three parts take the triples
\[
 (x,x+a,x+2a),\qquad 1\le x\le n,\qquad a\in D .
\]
Each pair lies in at most one such triple. The three projected graph
edges of any triangle have difference parameters \(a,b,c\in D\)
with \(a+b=2c\), so \(a=b=c\). The graph has exactly \(n|D|\)
triangles. This shows that excluding additional triangles in the
projected factor graph would still require an arithmetic restriction
to obtain a fixed exponent saving. [Behrend's progression-free sets](https://pmc.ncbi.nlm.nih.gov/articles/instance/1078964/pdf/pnas01693-0039.pdf)

### Private prime powers and exact additivity

The progression example also admits primitive additive values with
private prime powers at every label. The following construction
applies to any finite additive pattern.

Let \(A,B,C\) be finite label sets with integer assignments \(u_i,v_j,w_k\) that are injective within their own parts. Let \(E\) be any finite nonempty set of triples of labels with
\[
 u_i+v_j=w_k\qquad ((i,j,k)\in E).
\]
For each label prescribe a positive exponent \(e\). There exist positive integer values \(A_i,B_j,C_k\) preserving these additive relations, such that each label has a distinct prime of exact prescribed valuation, and that prime divides no other value in the three parts. Every selected triple \(A_i+B_j=C_k\) is primitive.

To prove this, take \(H\) at least \(|E|\) and let \(M_0\) be the product of all primes at most \(H\). Consider
\[
 A_i=L+M_0u_i,\qquad
 B_j=L'+M_0v_j,\qquad
 C_k=L+L'+M_0w_k .
\]
The additive relations hold for every \(L,L'\). Prescribe \(L\equiv L'\equiv1\) modulo each prime at most \(H\), so \(A_i,B_j\) are units there. Choose one distinct prime \(p>H\) for each label, also greater than the total number of labels and than the absolute differences among assigned integers in each part. In particular, \(p\nmid M_0\). The following prescriptions can all be made modulo \(p^{e+1}\).

- At a prime assigned to \(A_i\), take \(L\equiv p^e-M_0u_i\pmod{p^{e+1}}\). The other \(A\)-values are units. Choose \(L'\bmod p\) to avoid all the \(B\)- and \(C\)-zeros, and lift it arbitrarily.
- At a prime assigned to \(B_j\), do the symmetric construction.
- At a prime assigned to \(C_k\), impose \(L+L'\equiv p^e-M_0w_k\pmod{p^{e+1}}\). The other \(C\)-values are units. Choose \(L\bmod p\) so that neither an \(A_i\) nor a \(B_j\) vanishes, and lift subject to the prescribed sum.

The number of forbidden residues is at most the total number of labels, which is less than \(p\). Applying CRT to the small primes and the assigned prime powers gives residues \(L_0,L'_0\) modulo a fixed integer \(M\) with all the asserted local properties.

Choose the integer representatives \(L_0,L'_0\) so that all the finitely many differences
\[
 \Delta_{ij}=L'_0-L_0+M_0(v_j-u_i)
 \qquad ((i,j,k)\in E)
\]
are nonzero. This is possible by changing \(L'_0\) by a multiple of \(M\). Vary the two shifts together as \(L=L_0+Mt,\ L'=L'_0+Mt\). Any common divisor of \(A_i,B_j\) on an edge divides the fixed integer \(\Delta_{ij}\). At a prime \(p\mid M\), the CRT conditions already prevent simultaneous divisibility of the two terms. A prime \(p\nmid M\) dividing one of the finitely many \(\Delta_{ij}\) is greater than \(H\). For every edge, at most one residue \(t\bmod p\) makes both terms divisible by \(p\), since \(M\) is invertible modulo \(p\). There are at most \(|E|\le H<p\) forbidden residues in total. Choose an allowed residue for each such prime, and apply CRT once more to \(t\). Taking \(t\) sufficiently large makes every term positive and ensures \(\gcd(A_i,B_j)=1\) on every edge. The exact assigned valuations and their privacy are preserved because \(t\) only changes \(L,L'\) by multiples of \(M\).

For the preceding progression pattern \(y=x+a,\ z=x+2a\), choose
\(u_x=x,\ v_y=-2y,\ w_z=-z\). The construction realizes
\(n^{2-o(1)}\) primitive triples whose projected graph has unique
triangles on each edge, with arbitrary prescribed private
valuations. The construction gives no lower bound on the numerical
weights of the private powers relative to the complete terms, and
no upper bound on the remaining radicals.
