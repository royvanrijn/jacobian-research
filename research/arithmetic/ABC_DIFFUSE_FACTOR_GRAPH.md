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

### Three-letter divisor packing

Suppose the three sets of complete values in a primitive additive
graph have sizes \(N_A,N_B,N_C\). Use fixed signs on the letters
to write the relation as \(a+b=c\), where signed terms are allowed
and \(|a|,|b|,|c|\le H\). Assign each vertex a divisor record with
\(K_W\) possible records on letter \(W\), and let \(F_W\) be the
full divisor of the corresponding complete value specified by
its record. Let \(F\) be a lower bound for \(F_AF_BF_C\) on the
edges. Then the number of edges is at most
\[
 \left((1+2H/F)\,K_AK_BK_C\,N_AN_BN_C\right)^{1/2}.
 \tag{MS14}
\]
Indeed, within a fixed triple of records that occurs on an
edge, the \(F_W\) are pairwise coprime by primitivity.
Let \(a_r,b_{-r}\) be the sizes of the complementary residue
classes of the first two fibres modulo \(F_C\), and let \(n_C\)
be the size of the third fibre. For a fixed sum and residue,
the difference of two possible first terms is divisible by
\(F_A\), by \(F_B\) since the second terms also differ by it up
to sign, and by \(F_C\). Thus at most \(1+2H/F\) pairs in the
residue class can give the same sum. The edges for that residue
are bounded by the smaller of \(a_r b_{-r}\) and
\((1+2H/F)n_C\). Taking geometric means and using Cauchy–Schwarz
over \(r\) gives a bound of \(((1+2H/F)n_An_Bn_C)^{1/2}\) for
the three fibres. Sum this bound over record triples and use
\(\sum_{\rho}\sqrt{n_{W,\rho}}\le\sqrt{K_WN_W}\) on each letter.

This applies to edge-dependent choices \(d_{W,e}\mid x_{W,e}\)
at any fixed truncation, with \(F_W=\prod_e d_{W,e}^e\).
Keep one triple of dyadic size lists for the selected divisors
at subpower cost. Within it, choose one record independently
for each vertex among all the selected records it is offered.
Every edge retains its chosen records with probability at
least \(X^{-o_X(1)}\), since the number of such records per
vertex is bounded by a fixed power of the divisor function
of its complete value. Some assignment therefore preserves
the edge count up to a subpower factor. If the total normalized
root cost is at most \(\mu\), the dyadic lists give
\(K_AK_BK_C\le X^{(\mu+o_X(1))L}\). For numerical weight at
least \(w\), use \(F\ge X^w\). With heights at most
\(X^{1+o_X(1)}\), any such family with complete part exponents
at most \(c_W\) has normalized edge exponent bounded by
\[
 \frac{c_A+c_B+c_C+\mu+(1-w)_+/L}{2}.
 \tag{MS15}
\]
Here \(w\) is measured on the ordinary logarithmic height scale
and the graph exponent is divided by \(L\). In particular,
for a fixed margin \(w>1\) the last term is zero in the limiting
graph. A fixed margin below one makes the bound vacuous in the
small-budget limit.

For a balanced endpoint with \(d=1/3\) and
\(\Gamma=1/2+\gamma\), the
[near-minimum complete-label reduction](#the-complete-label-reduction-close-to-minimum-mass)
gives exponent at least \(2c-68\gamma/3\), up to arbitrarily
small losses, where \(c=1/6\) and \(\sum_W c_W=3c+\gamma\).
Consequently prescribed whole blocks on all three letters
of total mass below \(c-139\gamma/3\) cannot have total
numerical weight above one. More generally, for fixed
\(\eta>0\), the edges admitting a choice of that weight by
a fixed positive margin and at cost at most
\(c-139\gamma/3-\eta\) have a normalized count saving of
at least \(\eta/2\) relative to the lower exponent of the
complete-label reduction, in the limit of vanishing window
and extraction losses.

For a bounded continuity band of degrees at most \(R\), this
also tests the fractional upper quantile of the pooled
**three** degree measures. If its weight at
\(t<c-139\gamma/3\) exceeds \(1+\zeta\), a fixed
\(\zeta>0\), all but an exponentially small fraction of the
retained edges have a prime in the band of normalized
logarithmic size at least \(\zeta/(4R)\). As in (MS13), the
squarefree prime contributions in the bounded band can be
ranked and rounded down at the last contribution. In the
absence of such a prime, rounding loses at most \(\zeta/4\)
in numerical weight, so the chosen divisors contradict the
three-letter count bound with a fixed margin.

For a formal profile illustrating the new restriction, use the
same degree measure on each letter,
\[
 \nu=\tfrac1{20}\delta_{36/5}+\tfrac7{60}\delta_{6/5}
      +\tfrac16\,\operatorname{Unif}[2,4],
\]
where the uniform law has mass one. The total root mass is
\(1/3\), the nonatomic mass is \(c\), and the numerical
height is one. The three nonatomic supports lie in \([0,6]\)
and the integral in (MS12) is at least \(s/6\) for \(0\le s\le4\)
by concavity and the endpoint values. For any two letters,
the fractional top mass \(c\) has weight at most
\(2(9/25)+4(c-1/10)=74/75\), using degree at most four
after the two high atoms. Yet the three high atomic windows
together have mass \(3/20<c\) and weight \(27/25>1\).
The count saving extends below the extremal exponent.
Let \(d=2c-\varepsilon\), \(0\le\varepsilon\le10^{-3}\),
and suppose \(|c_W-c|\le10^{-3}\) on each letter.
Set \(\gamma=\sum_W c_W-3c\). The
[quantitative separation argument](#the-complete-label-reduction-close-to-minimum-mass)
now gives a complete-label graph of normalized exponent at least
\(2c-35\varepsilon-68\gamma/3\), up to arbitrarily small losses.
Indeed the simultaneous anchor deficit is at most
\(7\varepsilon+4\gamma\): substitute the initial exponent
\(2c-\varepsilon\) in the seventh power of the edge count.
The nondegenerate cycle count still dominates because
\(3d>2(c_i+c_j)+\max(c_i,c_j)\) for every pair in this range.
The same vector inequalities and freeze add at most four
times the anchor deficit plus \(8\gamma/3\).
Combining the lower exponent with (MS15) for \(w>1\) yields
\[
 \varepsilon\ \ge\ \frac{c-\mu-139\gamma/3}{70}.
 \tag{MS16}
\]
For the displayed profile the full high windows have
\(\mu=3/20\) and \(\gamma=0\), so it gives
\(d\le1/3-1/4200\). The same conclusion holds when \(d\) is
already below the stated range.

### Divisor records for longer additive relations

For \(m\ge3\), take \(m\) finite sets of signed integers of
magnitude at most \(H\), with the additive relation
\(\sum_{i=1}^m x_i=0\). The \(i\)-th set has \(N_i\) values
and \(K_i\) divisor records. On each counted tuple, assume
the recorded full divisors of its entries are pairwise
coprime and have product greater than \(2H\). Then the
number of tuples is at most
\[
 \prod_{i=1}^m N_i^{(m-2)/(m-1)}K_i^{1/(m-1)}.
 \tag{MS17}
\]
For a fixed tuple of records with divisors \(F_1,\ldots,F_m\),
take the uniform distribution on its counted tuples. The
entry \(x_i\) is determined by its residues at the other
\(m-1\) moduli, since its own residue is zero and the full
product exceeds \(2H\). Write the residue array in columns
indexed by the moduli, omitting the known diagonal zeros.
The additive relation implies that in each column any
\(m-2\) entries determine the last. Conditional on all
earlier columns, the entropy of the column therefore
equals the joint entropy with any one entry omitted.
Sum the \(m-1\) resulting subadditivity bounds to bound
the column entropy by \((m-2)/(m-1)\) times the sum of
its entry entropies, with the same conditioning. On summing
over columns, each row contributes at most its joint
entropy, since conditioning on all earlier columns gives
at most the entropy conditioned on the earlier entries
of that row alone. The joint entropy of the tuple is
therefore at most \((m-2)/(m-1)\sum_i\log n_i\), where
\(n_i\) are the respective fibre sizes. Exponentiate and
sum over the record tuples. Hölder on each record
dictionary gives (MS17).

The argument also allows choices depending on the tuple.
For fixed \(m\) and fixed truncation the same dyadic size
and random assignment argument as in (MS15) applies, with
a separate copy of the record map on each position of
the relation. In particular the \(m\) sets may coincide.

## Joint record entropy with one progression fibre

Take a distribution on finitely many signed additive tuples
\(\boldsymbol x=(x_1,\ldots,x_m)\), \(m\ge3\), with
\(\sum_i x_i=0\). On each slot let \(K_i\) be a function of \(x_i\)
recording a divisor \(F_i\) of \(x_i\). Assume the divisors on each
tuple are pairwise coprime. Write \(\boldsymbol K=(K_1,\ldots,K_m)\).
For a fixed record vector, put \(Q=\prod_i F_i\), and let \(J\)
bound the maximum occupancy in a residue class modulo \(Q\) on
any one distinguished slot. The slot may depend on the record
vector. It suffices to use the least available occupancy bound.
Then
\[
 (m-1)H(\boldsymbol x)-(m-2)\sum_i H(x_i\mid\boldsymbol K)
 \le (m-1)H(\boldsymbol K)+\mathbb E\log J .
 \tag{JR1}
\]
No lower bound on \(Q\) relative to the numerical heights is needed.

To prove (JR1), condition on a record vector and write
\(\mathcal R\) for the residue array at the \(F_i\), omitting the
known diagonal zeros. Write \(\mathcal R_i\) for its \(i\)-th row
and \(a=(m-2)/(m-1)\). The column argument in (MS17) gives
\(H(\mathcal R)\le a\sum_i H(\mathcal R_i)\).
With \(\mathcal R\) and the distinguished value \(x_j\) fixed,
any \(m-2\) of the remaining \(m-1\) entries determine the last.
Average the \(m-1\) resulting subadditivity bounds, then use that
conditioning reduces entropy. This gives
\[
 \begin{aligned}
 H(\boldsymbol x\mid\mathcal R)
 &\le H(x_j\mid\mathcal R)+a\sum_{i\ne j}H(x_i\mid\mathcal R,x_j)\\
 &\le a\sum_i H(x_i\mid\mathcal R)+\frac{\log J}{m-1}.
 \end{aligned}
\]
The last step follows because \(\mathcal R\) gives a full residue
of \(x_j\bmod Q\). Also \(H(\mathcal R_i)+H(x_i\mid\mathcal R)
\le H(x_i)\), since the row is a function of its value and of the
fixed record. Add the inequalities and average over records.
For each record vector this also gives the finite count
\[
 \#\{\boldsymbol x\text{ with these records}\}
 \le J^{1/(m-1)}\prod_i N_i^{(m-2)/(m-1)},
 \tag{JR2}
\]
where \(N_i\) bound the respective record fibre sizes. This follows
by taking the uniform distribution on that fibre.

### Application to complete-label graphs

At \(d=1/3\), \(\Gamma=1/2\), use a regularized near-full primitive
graph of complete values, with \(c=1/6\) and \(n=X^L\).
On each letter choose a record function \(K_T\) of the complete
value specifying a divisor. The product of the three divisors
on the edges under consideration has numerical weight at least
\(w-o(1)\). On a letter \(T\) allow all those values to lie in an
interval of length at most \(X^{h_T+o(1)}\). The \(h_T\) are the
exponents of the interval lengths, so an interval need not begin
at zero. Set
\[
 s_T=\max_{U\ne T}\sup\operatorname{supp}\nu_{U,\mathrm{na}},
 \qquad
 \eta=\min_T\min\left(c,\frac{(h_T-w)_+}{s_T}\right).
\]
Here \(s_T>0\), since both opposite nonatomic measures have mass
\(c\). The progression bound (AC1), applied on an appropriate
distinguished slot, gives \(\log J\le(\eta+o(1))\log n\) up to
the arbitrarily small graph and approximation losses. It applies
uniformly to all record vectors because (AC1) is uniform in
the modulus and residue class, with a finite grid of interval
exponents. One can regularize all three parts at once by pruning
vertices below the chosen slightly reduced degree threshold.
The part caps bound the number of removed edges.

Use the uniform edge distribution and introduce the normalized
entropies \(a_T=H(K_T)/\log n\), \(a=H(\boldsymbol K)/\log n\).
Apply (JR1) for \(m=3\), with \(H(\boldsymbol x)\) close to
\(2c\log n\) and each \(H(x_T)\le(c+o(1))\log n\).
Since \(H(x_T\mid\boldsymbol K)\le H(x_T)-H(K_T)\), the
limiting inequality is
\[
 2a-\sum_T a_T+\eta\ge c .
 \tag{JR3}
\]
The stronger finite form (JR1) keeps the actual conditional
entropies \(H(x_T\mid\boldsymbol K)\). For records given by fixed
blocks, the total divisor weight is the sum of their weights;
for variable divisors use a lower bound uniform on the
considered edges.

There is also a pointwise version outside exponentially sparse
edges. With \(\mathcal I=-\log\Pr(\boldsymbol K)/\log n\) and
\(\mathcal I_T=-\log\Pr(K_T)/\log n\), one has
\(2\mathcal I-\sum_T\mathcal I_T+\eta\ge c\) up to any fixed
positive tolerance, after taking sufficiently small regularization
and profile losses. In fact, on each of finitely many information
boxes with \(\mathcal I\le u\) and \(\mathcal I_T\ge v_T\), there
are at most \(n^u\) joint records. The lower degree bound
\(n^{c-o(1)}\) for each retained complete vertex makes its
record fibre size at most \(n^{c-v_T+o(1)}\), by counting
the edges with that record on letter \(T\). Applying (JR2)
and summing gives exponent
\(u+(3c-\sum_T v_T+\eta)/2+o(1)\). A fixed violation gives
a strict deficit from \(2c\) on a sufficiently fine grid.

In particular, if the product of the dictionary sizes has
normalized exponent at most \(p<c\), then
\(2a-\sum_Ta_T\le\sum_Ta_T\le p\). Thus
\[
 w+(c-p)s_T\le h_T \qquad(T=A,B,C).
 \tag{JR4}
\]
This gives an interval-width consequence even when all three
complete terms have numerical height one. For the broad symmetric
profile, a near-full family confined on any letter to an interval
of length \(X^{3/4+o(1)}\) requires a dictionary exponent at
least \(5/48\) for a complete degree-three atom factor on any
letter. Indeed \(s_T=4,w=1/2\) and
\(1/6-(3/4-1/2)/4=5/48\).

### Excluding three small atom dictionaries with a deterministic relation

At the broad symmetric profile
\(c\delta_3+c\operatorname{Unif}[2,4]\), take all three records
to be the entire degree-three atom factors of their letters.
Their product has numerical weight \(3/2\), so \(J=1\) on
these fibres and \(\eta=0\). By (CA2), each \(a_T\ge c/2\).
Therefore (JR3) gives
\[
 a\ge\frac{c+\sum_T a_T}{2}\ge\frac{5c}{4}=\frac5{24}.
 \tag{JR5}
\]
The complete values on any pair of letters determine the edge.
Their pair entropy is \(2c+o(1)\) in normalized units and their
marginal entropies total at most \(2c+o(1)\). The record
functions on that pair consequently have mutual information
\(o(\log n)\). In particular
\[
 \liminf\sum_T\frac{H(K_T\mid(K_U)_{U\ne T})}{\log n}
 \ \ge\ \liminf\frac{3c-\sum_T a_T}{2}.
\]
Thus a family where every atom record is determined by the other
two up to vanishing normalized conditional entropy must have
\(a_T\to c\) for all three letters. If each marginal is at the
individual lower endpoint \(c/2\), the normalized conditional
entropy of each single atom factor given the other two is at
least \(c/4=1/24\), by (JR5) and the two marginal upper bounds
on their pair entropy.

## Repeated locks force a common full factor

At \(d=1/3,\ \Gamma=1/2\), put \(c=1/6\) and use the
near-full complete-label additive graph. Fix one target letter
\(T\). Suppose a prescribed block of its original root positions
has limiting normalized root cost \(\beta\) and numerical
weight \(w\), with scaled degrees in a fixed compact interval
\([s,S]\subset(0,\infty)\). All its roots are squarefree because
these are bounded scaled degrees below the packed index.
Let \(m\ge2\) be a fixed even integer with \(m\beta<c\).

If \(mw>1\), then for every \(\epsilon>0\) the graph can be
thinned at arbitrarily small normalized cost to one on which
all complete values of \(T\) have a common full factor carried
by this block, with numerical weight at least
\[
 W_*-\epsilon,\qquad W_*=\frac{mw-1}{m-1}.
 \tag{RP1}
\]
In the common factor each prime is assigned to the same original
degree position throughout the retained graph. This last property
is needed both for iterating the argument and for combining the
common modulus with other disjoint blocks.

### Additive energy and the cost of one iteration for \(m\ge4\)

Retain one outside label and one admissible original representation
per complete value in each part, at fixed-truncation divisor loss.
For notation in the count, put \(n=X^L\), take part caps
\(n^{c+\delta+o(1)}\), and suppose the current retained edge
exponent is \(2c-\eta+o(1)\). First discard target vertices
of graph degree less than \(n^{c-\lambda}\), with
\(\lambda>\eta+\delta\). Their total edge count has a strict
saving, so the surviving graph still has exponent
\(2c-\eta+o(1)\). Keep this regularization before selecting
any further target subset.

The number of nondegenerate four-cycles in the projection
onto the other two complete letters is at least
\(n^{4c-4\eta-4\delta+o(1)}\), for small \(\eta,\delta\).
Each produces an additive quadruple of complete target values
\(T_{11}+T_{22}=T_{12}+T_{21}\). A fixed ordered quadruple
comes from at most \(n^{c+\delta+o(1)}\) cycles, since the
first complete row value determines the other three row and
column values. Thus the target part has additive energy at
least \(n^{3c-4\eta-5\delta+o(1)}\).
The degree \(m\) Fourier moment of its indicator exponential
sum is the number of ordered equal sums of \(m/2\) values on
each side. By Hölder interpolation of the second and fourth
moments,
\[
 E_m\ge E_4^{(m-2)/2}|\mathcal T|^{-(m-4)/2}
 \ge n^{(m-1)c-A_m+o(1)},\qquad
 A_m=2(m-2)\eta+(3m-7)\delta .
\]
Repeated target values cost only \(O_m(|\mathcal T|^{m-2})\)
tuples: if the two equal entries cancel, the remaining
nontrivial relation determines one of the remaining values;
otherwise use their common value as the coefficient-two
variable. They can be discarded.

Suppose a common full factor \(F\) on the specified block
has already been selected, with every prime of \(F\) in its
fixed assigned degree position. Let \(W=\log_X F\).
Divide that fixed prime subset out of the block roots and
call the resulting block of each vertex its residual block.
Its full numerical weight is \(w-W+o(1)\), and its radical
is coprime to \(F\). Let \(P(t)\) denote this residual
squarefree radical. Choose a fixed \(\xi>0\) below
\(\epsilon/(2mS)\).
For a tuple in which every pair satisfies
\[
 \log_X\operatorname{rad}\gcd(P(t_i),P(t_j))<\xi L ,
\]
record all \(m\) residual blocks. In each block strip every
prime appearing in any other recorded block. In total this
loses at most \(m(m-1)S\xi\) of numerical weight: a single
root prime occurs to its original bounded scaled degree.
The stripped moduli are pairwise coprime and divide the
reduced values \(t_i/F\). If
\(W<W_*-\epsilon\), their total numerical weight plus \(W\)
is at least
\[
 W+m(w-W)-m(m-1)S\xi-o(1)>1 .
\]
Hence their product exceeds \(2X/F\) for large \(X\).

For each fixed record tuple, the fixed-record bound in (MS17) gives at most
\((\prod_i a_i)^{(m-2)/(m-1)}\) such additive tuples, where
\(a_i\) is the size of its residual-record fibre.
The stripping is deterministic from the record tuple.
Each residual-record dictionary has exponent at most
\(\beta\), since its roots divide by fixed factors already
assigned to their positions. Hölder over all records
therefore gives exponent at most
\[
 \frac{m(m-2)}{m-1}(c+\delta)+\frac{m\beta}{m-1}
   =(m-1)c-\frac{c-m\beta}{m-1}+O_m(\delta)
\]
for tuples with no large overlap. The strict \(m\beta<c\)
makes this smaller than the energy exponent if all graph
losses are sufficiently small.

Consequently the tuples with a pair overlap at least
\(\xi\) have exponent at least \((m-1)c-A_m\).
Fix an offending ordered pair of slot positions, losing only
a constant factor. For fixed values of that pair there are at
most \(|\mathcal T|^{m-3}\) completions of a signed additive
tuple, because all but one remaining value can be chosen freely.
Thus the pairs with such overlap have exponent at least
\(2c-A_m-(m-3)\delta\). Some fixed \(t_0\) has at least
\(n^{c-A_m-(m-2)\delta+o(1)}\) such neighbours.
For each choose one squarefree divisor
\(g\mid P(t_0)\) dividing its residual block, with
\(\log_X g\ge\xi L\). There are divisor-many choices of
\(g\), and at fixed truncation \(J^{\omega(g)}=X^{o_X(1)}\)
ways to assign its primes to the neighbour's \(J\) degree
positions. Pigeonhole both the divisor and those assignments.
The retained neighbours therefore share the same new full
factor \(G=\prod_e g_e^e\), with numerical weight at least
\(s\xi+o(1)\), coprime to \(F\).

Since each retained neighbour had graph degree at least
\(n^{c-\lambda}\), retaining all its edges gives the next
edge deficit at most \(A_m+(m-2)\delta+\lambda+o(1)\).
For example take \(\lambda=\eta+2\delta+\delta'\) with
\(\delta'>0\) arbitrarily small. The update then costs at
most
\[
 (2m-3)\eta+(4m-7)\delta+\delta'+o(1).
\]
Replace \(F\) by \(FG\), whose prime assignments remain
fixed, and repeat unless \(W\ge W_*-\epsilon\).
Each step increases \(W\) by at least \(s\xi/2\) for large
\(X\), and \(F\) divides every currently retained \(T\)-value
at most \(X\), so at most \(1+\lceil2/(s\xi)\rceil\) steps
are possible. All parameters are fixed before the height
limit. Choose \(\eta,\delta,\delta'\) initially small
enough that the finitely many possible updates leave all
the stated strict margins and the requested final graph
tolerance. This proves (RP1).

### Two locks on a fixed difference

For \(m=2\), use the same target energy bound
\(E_4\ge n^{3c-4\eta-5\delta+o(1)}\) and regularization.
Write \(r(D)\) for the number of ordered target pairs of
difference \(D\). Then \(\sum_D r(D)^2=E_4\) and
\(\sum_D r(D)=|\mathcal T|^2\). The zero-difference
term \(r(0)^2=|\mathcal T|^2\) has strictly smaller
exponent than the energy. Hence some fixed nonzero
\(D\) has \(r(D)\ge n^{c-4\eta-7\delta+o(1)}\).
For a pair with this difference and overlap less than
\(\xi L\) in logarithmic root mass, record the two
residual blocks and strip their shared primes. If
\(W<2w-1-\epsilon\), the resulting coprime full moduli
have product exceeding twice the height bound for the divided target
values, by the same weight calculation with \(m=2\).
One modulus divides the first value and the other fixes
it modulo the second modulus using \(D/F\).
There is at most one pair for these records. Their
exponent is at most \(2\beta<c\), so pairs with large
overlap retain the lower exponent above.

Each such pair determines a radical divisor of \(D\)
in its residual blocks of normalized root mass at least
\(\xi\). Assign those primes to original degree positions
on the first value and pigeonhole the divisor and
assignments as before. Each first value occurs in just
one fixed-difference pair. The resulting set of first
values retains the pair exponent, so the next edge
deficit is at most \(4\eta+7\delta+\lambda\) up to
arbitrarily small losses. Iteration proves (RP1)
also for \(m=2\).

### A support consequence for high atomic blocks

If the block in (RP1) lies entirely above degree six, it
is disjoint from any supplying nonatomic band on \(T\),
by the complete-label support bound. Write
\(s_V^*=\sup\operatorname{supp}\nu_{V,\rm na}\) for each
letter \(V\). Then
\[
   mw+(m-1)c\,\max_V s_V^*\le m .
 \tag{RP2}
\]
For \(mw\le1\) this follows from the support bound \(\max_V c s_V^*\le1\).
Otherwise select the common factor as in (RP1) and take
a complete additive star with fixed rows on one opposite
letter and common complete columns on the other.
All its \(T\)-terms are divisible by \(F\); the column
parameter already has a fixed residue modulo \(F\), using
any one fixed row. For each \(u<s_T^*\), choose a
continuity band of positive diffuse mass with lower degree
at least \(u\) and upper degree at most the high block's
lower degree. The chosen diffuse outside roots have no
prime in \(F\), since the prime assignment of \(F\) is
fixed in the disjoint high positions. Use the fixed-row
CRT count to lock this band across sufficiently many rows
at combined mass \(z<c\) and weight at least \(u z\).
Strip fixed row-difference primes exactly as in (CT1).
The surviving moduli are pairwise coprime and coprime
to \(F\). If \(W+uz>1\) with a fixed margin, the column
is determined by CRT, including the fixed congruence
modulo \(F\). The same bad-overlap estimate works since
the retained graph stays arbitrarily close to full.
This contradiction shows \(W_*+cs_T^*\le1\) on taking
\(\epsilon\to0\). For a supplying band on a different
letter \(V\), take the common columns on \(T\) itself.
The complete column parameter is divisible by \(F\),
and primitivity makes all \(V\)-block primes coprime
to \(F\). The same count gives \(W_*+cs_V^*\le1\).
This proves (RP2).

For a strict formal example, give both other letters an
atom of mass \(1/6\) at three and nonatomic mass \(1/6\)
uniform on \([2.9,3.1]\). Give \(T\) an atom of mass
\(1/25\) at \(15\) and one of mass \(19/150\) at
\(16503/15200\), together with nonatomic mass \(97/600\)
uniform on \([1.49,1.51]\) and mass \(1/200\) uniform on
\([3.99,4]\). These measures all have total mass \(1/3\),
nonatomic mass \(1/6\) and numerical moment one. The
fractional three-letter packing at mass \(1/6\) has
weight at most
\[
 0.6+4(0.005)+3.1(1/6-0.04-0.005)<1 .
\]
The single-parameter trimmed bound also holds: at \(s=3.1\)
the truncated moment of \(T\) is \(0.519525>3.1/6\), and
the ratio of this truncated moment to \(s\) is
nonincreasing, so the opposite diffuse supports are allowed.
For \(T\) itself its opposite letters' moments are untrimmed
at the end of its diffuse support and there is no restriction.
But take \(m=4,\ \beta=0.04,\ w=0.6,\ s_T^*=4\).
Then \(m\beta=0.16<c\) while
\(mw+(m-1)c s_T^*=4.4>4\). This violates (RP2).

### Common factors on several letters

Suppose a near-full complete-label graph has a common full
factor \(F_i\) on each letter \(i\), where empty factors are
allowed. The factors are pairwise coprime by primitivity.
For each \(i\), the opposite two complete parts have
complementary residues on all their edges modulo \(F_i\)
(with fixed signs if the equation is differently oriented).
A single residue pair retains an arbitrarily near-full
edge exponent. Indeed, if \(E_r\) is the edge count in a
matching residue pair and \(a_r,b_r\) the relevant two
part sizes, then
\[
 E=\sum_r E_r
 \le (\max_r E_r)^{1/2}\sum_r\sqrt{a_rb_r}
 \le(\max_r E_r)^{1/2}\sqrt{|\mathcal A||\mathcal B|}.
\]
Thus \(\max_r E_r\ge E^2/(|\mathcal A||\mathcal B|)\).
Apply this for all three \(F_i\). Each complete part now
lies in a single residue class modulo
\(Q=\prod_iF_i\), with numerical weight
\(W=\log_X Q=\sum_iW_i\). Each of these finitely many
thinnings multiplies the incoming edge deficit only by a
fixed constant and preserves an arbitrarily near-full graph.
As each part still has many complete values, \(Q\) is
at most their height bound \(X^{1+o_X(1)}\).

For each selected letter \(i\), let a prescribed block of
root cost \(\beta_i\), numerical weight \(w_i\), and compact
positive scaled range \([s_i,S_i]\) be given. Choose an
even integer \(m_i\ge2\) with \(m_i\beta_i<c\). Let \(F_i\)
be the common factor currently extracted from this block,
with the fixed prime assignments as in (RP1); put \(w_i=0\)
for a part where no block is used. Fix a small threshold
\(\epsilon_0>0\). If
\(W+m_i(w_i-W_i)>1+\epsilon_0\) on some selected part, apply
the preceding repeated-lock extraction on that part,
with the following change to the range argument.
Write each target integer as \(r_i+Q y\), where \(r_i\)
is the fixed residue on its part. Its reduced range has
length at most \(X/Q\). Balanced additive tuples, or
a fixed nonzero difference when \(m_i=2\), retain the
same additive count after this affine change.

Record the residual locks once the common assigned part on
letter \(i\) has been removed. Their primes are coprime to
\(Q\): they are disjoint from the assigned common factors
on their own letter, and coprime to the other letters'
common factors by primitivity. After stripping pair
intersections in a private tuple, each reduced variable
\(y\) has a fixed residue modulo its own residual modulus,
namely \(-r_i/Q\) modulo that modulus. This need not be
the zero residue, but the residue-entropy lemma is
unchanged for prescribed affine residues. In each
residue column the missing entry is now a known constant;
the remaining entries satisfy the corresponding known
linear sum. The private record product exceeds \(X/Q\)
as long as
\(W+m_i(w_i-W_i)-m_i(m_i-1)S_i\xi>1\), using
two-variable CRT interpretation for \(m_i=2\).
Choose a common small \(\xi>0\) to ensure this whenever
the extraction is triggered. The energy and overlap
arguments then give a new common assigned factor on
letter \(i\), coprime to the previous \(Q\), of weight at
least \((\min_i s_i)\xi+o(1)\).

After each extraction, synchronize residues for the new
factor and keep the previously synchronized residues.
Each step has the same finite-deficit amplification as
before, with the additional residue-class thinning
above. There are at most
\(1+\lceil 2/((\min_i s_i)\xi)\rceil\) extraction steps
unless a contradiction has occurred. By starting with
small enough losses, the graph remains arbitrarily
near-full. Thus the process terminates with
\[
 W+m_i(w_i-W_i)\le1+\epsilon_0
 \qquad\text{on every selected letter.}
\]
Taking \(\epsilon_0\) to zero along subsequences of the
bounded weights gives the limiting necessary inequalities
with right side one.

If all selected blocks lie above the nonatomic supports,
the resulting synchronized modulus \(Q\) can be included
in every fixed-row CRT reconstruction. For any target the
nonatomic row locks are prime-disjoint from its own
assigned common high factor and coprime to the opposite
common factors by primitivity. The column parameter
already has one fixed residue modulo the entire \(Q\).
Let \(s_*=\max_i\sup\operatorname{supp}\nu_{i,\mathrm{na}}\).
It follows that \(W+c s_*\le1\). Putting \(t=1-W\), the
residual inequalities yield
\[
 t\ge c s_*,
 \qquad
 t+\sum_{\text{selected }i}(w_i-t/m_i)_+\le1 .
\]
If \(\sum_i1/m_i\le1\), the function on the left is
nondecreasing in \(t\), and consequently
\[
 c s_*+\sum_{\text{selected }i}(w_i-c s_*/m_i)_+\le1 .
 \tag{RP3}
\]

On each letter put an atom of mass \(1/25\) at
\(29/4\), an atom of
mass \(1/6-1/25\) at \(q=423603/175400\), nonatomic
mass \(1/1000\) uniform on \([3.99,4]\), and the remaining
nonatomic mass \(1/6-1/1000\) uniform on
\([q-0.001,q+0.001]\). The three numerical weights per
letter are \(0.29\) on the high atom, \(0.003995\) on the
diffuse tail, and \(0.706005\) on the combined low mass.
Thus the moment is one. The three-letter fractional
top-\(c\) weight is at most
\(0.87+0.012+(1/6-0.123)\cdot2.417<0.988\),
and (MS12) holds because the truncated moment at four
is \(0.87>4/6\). The independent repeated locks with
\(m_i=4\) and \(\beta_i=0.04,w_i=0.29\) only force
total common weight \(0.16\), compatible with the support
bound. But the synchronized residual inequalities force
\(W\ge4(3\cdot0.29)-3=0.48\); the diffuse limit
\(s_*=4\) requires \(W\le1/3\). Therefore (RP3) supplies
a strict additional restriction in this sector.

## A trimmed packing bound using mixed additive moments

At the minimum nonatomic endpoint, set \(c=1/6\) and
\(s_*=\max_T\sup\operatorname{supp}\nu_{T,\mathrm{na}}\).
Choose arbitrary prescribed whole-coordinate blocks on the
three letters (including empty blocks), of total normalized
root cost \(p<c\) and total numerical weight \(w\). Then
\[
 w+s_*(c-p)\le1.
 \tag{MXP1}
\]
In particular the three distinguished blocks can all be atomic.
The weighted contribution \(s_*(c-p)\) is available from
nonatomic blocks in repeated copies of just one letter.

### Rich mixed additive relations

Write \(n=X^L\). Take the complete-label graph with edge
exponent \(2c-\eta\) and complete part sizes at most
\(n^{c+\delta}\), where \(\eta,\delta\) can be arbitrarily
small at the exact endpoint. Let \(S_A,S_B,S_C\) be
the exponential sums of the three complete vertex sets.
The graph's additive triples give
\[
 n^{2c-\eta+o(1)}
 \le \left|\int_0^1 S_A(\theta)S_B(\theta)
                    \overline{S_C(\theta)}\,d\theta\right|,
\]
with conjugations changed for a different signed orientation.
The integral may also count other triples in the part sets.
On the subset where some \(|S_i|<n^{c-\lambda}\), its
absolute contribution is at most
\(3n^{2c-\lambda+\delta+o(1)}\), by Cauchy--Schwarz
and the \(L^2\) norms of the other two sums. Fix
\(\lambda>\eta+\delta\). Thus the common large spectrum has
measure at least \(n^{-c-\eta-3\delta+o(1)}\), by the
trivial joint amplitude upper bound.

For any fixed nonnegative even counts \(m_i\), and
\(m=m_A+m_B+m_C\), this gives at least
\[
 n^{(m-1)c-m\lambda-\eta-3\delta+o(1)}
\]
ordered tuples in the balanced mixed additive relation with
\(m_i/2\) plus and \(m_i/2\) minus slots from part \(i\).
Indeed its count is \(\int\prod_i |S_i|^{m_i}\).

### One distinguished lock per letter

If (MXP1) fails, choose \(s<s_*\) and \(0<z<c-p\) with
\(w+sz>1\). A compact positive continuity band \([s,S]\)
above this lower cutoff has positive nonatomic mass on
some letter. With sufficiently many copies of that letter,
select outside root positions of total mass \(z+o(1)\)
in this band. The normalized outside coordinate cap can
be taken arbitrarily small, so there is arbitrarily small
rounding loss. On at most one further slot per letter
put its prescribed block. Add empty slots so that the
per-letter counts are even, \(m\ge4\), and each letter
supplying outside locks has a slot of a different letter.

Keep the prime coloring by letter, at subpower cost by
the fixed-truncation divisor bound. Fix one admissible
original representation of each complete value, retaining
the outside label that determines it. For an outside slot
\(h\) call its record bad if its radical shares root mass
at least \(\xi\) with the union of the primes of the other
locks on its letter. Fix all complete values except the
one in \(h\) and one slot of a different letter. The union
is then fixed and has polynomial height. The
[divisor allocation estimate](#counting-shared-prime-divisors)
gives at most \(n^{c+\delta-\xi+o(1)}\) outside labels
for \(h\). The remaining slot is determined by the
additive relation. The total bad exponent, for this or
any other outside slot, is at most
\[
 (m-1)(c+\delta)-\xi .
\]
This treats all possible overlaps of each outside slot in
one count. The distinguished blocks of different letters
are already coprime by the coloring.

For good tuples retain the distinguished block on each
letter and strip from every outside block the primes
that occur in any other block. Since the outside slots
were good and their degrees are at most \(S\), the full
numerical weight lost is at most \(S k\xi+o(1)\), where
\(k\) is the number of outside slots. The surviving
moduli are pairwise coprime and their product exceeds
twice the height bound if \(w+sz-Sk\xi>1\).
With the original records fixed, the stripping is
deterministic. The residue count (MS17) and Hölder over
the original record dictionaries give a good-tuple exponent
at most
\[
 \frac{m(m-2)}{m-1}(c+\delta)+\frac{p+z}{m-1}.
\]
At the exact endpoint first fix the positive numerical
and cost margins and then let the graph losses decrease.
Both exponents fall strictly below the mixed moment
lower bound, which proves (MXP1).

### Consequences for degree profiles

For any finite collection of atomic windows of combined
mass below \(c\), (MXP1) bounds the selected excess
weight above \(s_*\) by \(1-cs_*\). When the combined atomic mass
above \(s_*\) is strictly below \(c\), increasing finite
selections gives the bound for the whole upper tail.
In particular, if \(s_*=6\), an atom at \(r>6\)
cannot occur. Its mass is at most \(1/r<c\) by the
numerical moment, so selecting it would violate (MXP1).
Hence all three full degree measures are supported in
\([0,6]\) in this case.

The distinguished locks may also be chosen as divisors
depending on the edge. With a fixed positive violation
margin, bin the costs and weights of the witnessed records
on each letter before randomly assigning one admissible
record to each complete value. The fixed-truncation
divisor bound preserves the edge exponent, as in (MS15).
Each record then satisfies the same binned bounds even
when its mixed tuple uses a value from a different
edge. Reserve a fixed small cost slack below \(c-p\) for
the diffuse locks and apply the preceding argument.
In particular an adaptive sector with
\(p\le c-\epsilon\) and \(w+s(c-p)>1+\epsilon\) has a
fixed exponent deficit when \([s,S]\) supplies a fixed
positive nonatomic mass and the graph losses are small.

### A deficit inequality for balanced profiles

Put \(d=2c-\varepsilon\), \(\rho=\max_T|c_T-c|\),
and \(\gamma=\sum_T c_T-3c\). Assume
\(0\le\varepsilon\le10^{-3}\) and \(\rho\le10^{-3}\),
as in (MS16). Use at most one prescribed block per
letter, of total cost \(p\) and weight \(w\). For each of
\(k\ge1\) additional outside slots \(h\), choose a
positive continuity band \([s_h,S_h]\) with nonatomic
mass strictly exceeding a chosen lock mass \(z_h>0\).
The slots can use different bands and letters. Write
\[
 \mu=p+\sum_h z_h,\qquad
 W_{\rm lock}=w+\sum_h s_hz_h,\qquad
 S_{\rm sum}=\sum_h S_h .
\]
Assume \(\mu<c\), and pad with empty slots as in the
mixed moment argument, using \(m\) total slots. Then
\[
 (m+1)\left(35\varepsilon+\frac{68}{3}\gamma+2\rho\right)
 \ \ge\
 \min\left\{\frac{(W_{\rm lock}-1)_+}{S_{\rm sum}},
                \frac{c-\mu+\rho}{m-1}\right\}.
 \tag{MXP2}
\]
If \(W_{\rm lock}>1\), strip overlap of root mass less
than \(\xi\) on each outside slot, for any fixed
\(\xi<(W_{\rm lock}-1)/S_{\rm sum}\). The graph deficit
\(\eta\) is at most \(35\varepsilon+68\gamma/3\) by the
complete-label reduction used in (MS16), up to arbitrarily
small additional losses. With part caps \(c+\rho\),
the moment lower exponent is
\((m-1)c-(m+1)\eta-(m+3)\rho\).
Comparing it with the bad- and good-record exponents
shows that \((m+1)(\eta+2\rho)\) must be at least
\(\min\{\xi,(c-\mu+\rho)/(m-1)\}\).
The second term uses
\(m(m-2)/(m-1)=m-1-1/(m-1)\).
Let the auxiliary losses decrease and \(\xi\) increase
to its stated upper limit. When \(W_{\rm lock}\le1\),
the right side is zero and follows from the pair
capacity bound \(2\gamma\ge-3\varepsilon\).
Thus (MP2) holds for every admissible selection and
can be optimized over all its bands and lock sizes.

### A quantitative saving in a three-letter sector

Consider balanced boxes, with \(d=2c-\varepsilon\), and
assume \(|c_T-c|\le1/50000000\) on every letter. Suppose
each letter has a prescribed block of root cost at most
\(41/1000\) and numerical weight at least \(289/1000\).
On at least one letter require nonatomic mass above
\(9/10000\) in \([399/100,4]\). Then
\[
 d<\frac13-\frac1{50000000}.
 \tag{MXP3}
\]
This sector contains the symmetric high-atom profile
used to separate (RP3) from the independent tests.

Indeed let \(\rho=1/50000000\) and suppose
\(0\le\varepsilon\le\rho\). The
[complete-label reduction below \(1/3\)](#three-letter-divisor-packing)
gives deficit \(\eta\le35\varepsilon+68\rho\le103\rho\),
with part caps \(c+\rho\), all up to arbitrarily small
additional losses. Use 43 outside slots on the supplying
letter of combined root cost \(z=383/10000\).
Each slot needs less than \(9/10000\). Add its distinguished
slot and take two slots on each other letter, one with its
distinguished block. This makes \(m=48\), with total
record cost at most \(p+z=1613/10000\). Set
\(\xi=11/100000\). The stripped good-weight lower bound is
\[
 3\frac{289}{1000}+\frac{399}{100}\frac{383}{10000}
 -4\cdot43\frac{11}{100000}
 =1.000897>1 .
\]
The mixed moment lower exponent approaches
\((m-1)c-(m+1)\eta-(m+3)\rho\) as the auxiliary
Fourier cutoff loss tends to zero. For the bad tuples,
\(\xi>(m+1)(\eta+2\rho)\) with room. For the good tuples
the gap available from the record count is
\[
 \frac{c-1613/10000}{47}>
 49(103\rho)+\left(51+\frac{48\cdot46}{47}\right)\rho .
\]
Both branches contradict the mixed moment lower bound.
All the comparisons are strict under arbitrarily small
truncation and box losses, proving (MXP3).

### Atomic factor information in complete-label graphs

At the minimum nonatomic endpoint, take a graph whose edge
exponent is arbitrarily close to \(2c\) and complete part
caps are \(c+o(1)\). Choose any prescribed root block on
each of a fixed set \(I\) of letters, with numerical weights
\(w_i\). Let \(F_i\) denote the complete product of that
block, as a function of the chosen complete vertex, and
let \(w=\sum_{i\in I}w_i\).
For a uniformly sampled edge write \(N_0=L\log X\) and
\(\mathcal I_i=-N_0^{-1}\log\Pr(F_i)\).
With \(s_*=\max_j\sup\operatorname{supp}\nu_{j,\rm na}>0\),
the limiting weights satisfy
\[
 \liminf\sum_{i\in I}\frac{H(F_i)}{N_0}
 \ \ge\ \max\left\{0,\,
       c-\frac{(1-w)_+}{s_*}\right\}.
 \tag{MXP4}
\]
The same lower bound, minus any fixed positive tolerance,
holds for \(\sum_{i\in I}\mathcal I_i\) outside an
exponentially sparse edge subset, with sufficiently small
initial graph and atom-window losses.

Indeed, write \(u\) for the displayed candidate lower
bound. If \(u>0\), cover the event
\(\sum_{i\in I}\mathcal I_i\le u-\epsilon\)
by finitely many information boxes with nonnegative upper
corners \(u_i\) satisfying \(\sum_i u_i\le u-\epsilon/2\).
On each such box, which selects edges of the complete
graph, the block value on letter \(i\) has dictionary
size at most \(n^{u_i}\), by the probability definition.
The proof of (MXP1) uses only the dictionary
sizes of its distinguished slots, not the normalized raw
root costs that bounded them. For small fixed graph loss,
choose a nonatomic band lower degree \(s<s_*\) close to
\(s_*\), and filler mass just below \(c-\sum_i u_i\).
Then the blocked mixed additive tuple has total root
dictionary cost strictly below \(c\) and numerical weight
strictly above one. The same overlap and entropy argument
gives a fixed loss from exponent \(2c\) for the selected
edges. Summing the finitely many boxes proves sparsity,
and taking expectations of the nonnegative information
variables proves the stated lower bound.

For the symmetric profile \(c\delta_3+c\operatorname{Unif}[2,4]\)
use its full atomic windows as the three possible blocks.
Each has weight \(1/2\). The inequality gives normalized
information at least \(1/24\) for every individual atomic
block, at least \(1/6\) for every pair, and hence total
normalized entropy at least \(1/4\). The same lower
bounds hold for the grouped atomic radical roots: given
the grouped root, allocation of its prime factors to
the fixed atom-window positions has only divisor-many
possibilities, so the atomic full block function is
determined up to conditional entropy \(o(N_0)\).

## Converting high-degree blocks to common factors with diffuse locks

At the minimum nonatomic endpoint let \(s_*\) denote
the largest nonatomic support endpoint, and put \(c=1/6\).
Suppose a prescribed block on one letter has normalized
root cost \(\beta<c\), with degrees in a compact band
\([s,S]\) above \(s_*\). Allow previously assigned common
factors, from disjoint blocks above \(s_*\), whose
synchronized product has numerical weight at least \(W_0\)
on an arbitrarily near-full graph. The new block can be
made common up to arbitrarily small root mass if
\[
 W_0+cs-(s-s_*)\beta>1 .
 \tag{MC1}
\]
All previously assigned common factors persist.
In particular the numerical weight of the new common
factor can approach the full numerical weight of the block.

### Extraction using different kinds of locks

Use the residue-class selection for (RP3) to synchronize
the assigned factors. Write their product as \(Q\), its
weight as \(W\), and each complete term as \(r_i+Qy\)
with \(r_i\) fixed on its letter. On each part the \(y\)'s
lie in an interval of length at most \(X/Q\) up to the
height factor of the original box. The mixed moment
argument for (MXP1) goes through for the reduced sets.
Their three-term equation has a fixed additive offset,
which adds a complex exponential of modulus one to the
triple Fourier integrand. Balanced mixed relations have
the same counts in the \(y\)'s as in the complete values.
A lock coprime to \(Q\) prescribes a residue of \(y\),
and the counting proof of (MS17) works with each
diagonal entry of the residue array prescribed.

Fix a desired positive residual root tolerance \(b_0\).
Suppose the current residual part of the new block has
normalized root cost \(b\ge b_0\), after dividing out all
assigned prime powers. Its cost agrees with the residual
dyadic root cap up to vanishing error. Choose a fixed
\(u<s_*\) close to \(s_*\), with a bounded band of positive
nonatomic mass above \(u\). Also choose a fixed \(\kappa>0\)
sufficiently small that \(c-\kappa>\beta\) and
\[
 W_0+s(c-\kappa)-(s-u)\beta>1 .
\]
Take the largest integer \(k\ge1\) with \(kb<c-\kappa\)
and put \(z=c-\kappa-kb\). On \(k\) slots of the block's
letter lock its residual part. On outside slots lock
mass \(z\) from the supplying diffuse band. Add empty
slots to make every count even and to leave a slot on
a different letter for the overlap argument. The total
number of slots is bounded independently of \(b\ge b_0\).
All these locks are coprime to \(Q\). On their own letter
they are disjoint from the assigned positions or primes,
and across letters the prime coloring and primitivity
give coprimality. Their full total weight plus \(W\)
is at least
\[
 W+ksb+uz
 =W+u(c-\kappa)+(s-u)kb
 \ge W_0+s(c-\kappa)-(s-u)\beta>1 ,
\]
up to arbitrarily small dyadic losses. Here
\(kb\ge c-\kappa-b\), regardless of the parity of \(k\).
The dictionary cost is at most \(c-\kappa\) up to
rounding from the outside positions.

Count and strip the outside overlaps by the union
argument for (MXP1). A uniform small overlap threshold
makes their numerical loss less than one third of the
weight margin. For the other \(k\) slots use a fixed
threshold \(\xi>0\) on the numerical weight of each
full pairwise gcd. If all such weights are below \(\xi\),
retain only the largest power of each shared prime.
The full product becomes the least common multiple.
The loss in numerical weight is at most
\(\binom{k}{2}\xi\). Prime by prime this follows because
every occurrence beyond the largest is counted in at
least one of the pairwise minima of the exponents.
Choose \(\xi\) making this loss less than another third
of the weight margin, uniformly for the finitely many
possible \(k\). The product of the resulting coprime
moduli exceeds the reduced range, with a fixed strict
height margin.

The residue count (MS17), with the strictly smaller
record cost than \(c\), bounds these tuples below the
mixed moment exponent. The outside-overlap branch
also has a strict deficit by (MXP1). Consequently the
mixed moment count persists on tuples containing a
pair of residual-block slots with full gcd of weight at
least \(\xi\). If \(k=1\), there is no such pair and the
assumed near-full graph is already impossible. For
\(k\ge2\), fix an offending slot pair. For fixed values
of this pair all but one of the other entries determine
the last, so there are \(n^{2c-o(1)}\) offending value
pairs as the graph losses tend to zero. Fix a pivot
with \(n^{c-o(1)}\) such neighbours. Pigeonhole the full
gcd divisor of its residual block and the neighbour's
degree assignments for the primes of that divisor.
There are only subpower-many such records at fixed
truncation. On the neighbour set keep the full original
powers of those primes, which adds a new common factor
of numerical weight at least \(\xi\).

As in (RP1), first prune the target vertices below degree
\(n^{c-\lambda}\), for \(\lambda\) slightly exceeding the
current graph deficit plus the part-cap error. Thus the
selected neighbours retain a near-full graph. The mixed
moment loss, pair-completion count, degree pruning and
residue synchronization multiply graph deficits by
bounded constants for the bounded slot counts. Each
iteration increases common weight by at least \(\xi\).
There are therefore finitely many iterations, bounded
by \(1+\lceil 2/\xi\rceil\), since the assigned factors
divide the complete values. Choose the initial losses
small enough for this many iterations. Eventually
\(b<b_0\), where the missing numerical weight is at
most \(S b_0+o(1)\), proving (MC1).

### Support and further packing after conversion

The same synchronized affine coordinates in the proof of
(MXP1) give
\[
 W+w+s_*(c-p)\le1
 \tag{MC2}
\]
for additional prescribed blocks of cost \(p<c\) and
weight \(w\), disjoint from the common primes. The
numerical product is now tested against \(X/Q\) and
the dictionary cost remains \(p\). The outside supplying
band can be taken disjoint from the common positions
because their degrees are above \(s_*\).
Taking \(p=w=0\) recovers \(W+cs_*\le1\).

In particular every atom at \(r>12-s_*\) can be made
common up to arbitrary tolerance. Its mass \(\beta\)
is at most \(1/r<c\), and (MXP1) gives
\((r-s_*)\beta\le1-cs_*\). Therefore
\(cr-(r-s_*)\beta\ge c(r+s_*)-1>1\).
Use a sufficiently small continuity window about \(r\)
in (MC1). Iterate over finite lists of such atoms and
take increasing lists in (MC2), to get
\[
 cs_*+\sum_T\int_{(12-s_*,\infty)}r\,d\nu_{T,\mathrm{at}}(r)\le1 .
 \tag{MC3}
\]
The cutoff can be lowered as common weight increases.
In a graph already synchronized at weight \(W\), the
same argument with (MC2) converts any further atom
at \(r>s_*\) of mass below \(c\) if
\(r+s_*>12(1-W)\). The separate criterion (MC1) also
allows conversion below this cutoff when that atom has
sufficiently small mass.

For a formal example, give \(T\) an atom of mass \(9/200\)
at \(8\), an atom of mass \(73/600\) at \(84/73\), and
nonatomic mass \(c\) uniform on \([2,4]\). Give each
opposite letter an atom of mass \(c\) at three and
nonatomic mass \(c\) uniform on \([2,4]\). Every first
moment is one. Here \(s_*=4\), and the largest possible
left side of ordinary mixed packing is at most
\[
 4c+\frac9{200}(8-4)=\frac{127}{150}<1 .
\]
The only even \(m\) with \(m(9/200)<c\) is \(m=2\).
Its residual repeated weight \(2(9/25)\) lies below one,
so the original repeated-lock test supplies no common
factor. But (MC1) applies because
\(8c-(8-4)(9/200)=173/150>1\).
Its common weight \(9/25\) violates (MC2) at zero
additional lock cost, since \(9/25+4c=77/75>1\).

### An iterated inequality for the full atomic spectrum

The conversion and packing inequalities can be applied
without specifying a finite list in advance. Set
\(D_0=1-cs_*\). For an atom at \(r>s_*\) of mass
\(0<\beta_{T,r}<c\), write
\(u_{T,r}=(r-s_*)(c-\beta_{T,r})\). Define
\[
 D_{j+1}
 =D_0-\sum_{\substack{T,\ r>s_*\\
              0<\beta_{T,r}<c,\ u_{T,r}>D_j}}
                r\beta_{T,r}.
 \tag{MC4}
\]
Every iterate must be nonnegative. In particular the
iteration must converge from above to a nonnegative
fixed point \(D\).

To justify this even for infinitely many atoms, the
inductive assertion is that every finite subset selected
at a given stage can be made common up to arbitrary
tolerance. For any finite subset at the next stage, the
inequalities \(u_{T,r}>D_j\) have a positive margin.
Take a sufficiently large finite subset of the preceding
stage so that its common numerical weight approximates
\(D_0-D_j\) within that margin. The previous induction
applies to this finite subset. Then (MC1) converts the
remaining desired atoms one by one and the previous
weight can only increase. All windows in this construction
can be disjoint. Apply (MC2) to finite selections and
increase them to obtain \(D_{j+1}\ge0\).
The right-hand side of (MC4) is nondecreasing in \(D_j\),
and \(D_1\le D_0\), so the sequence decreases. Its selected
sets increase to the set from the strict inequality
\(u_{T,r}>D\), giving the fixed point by monotone convergence.

For every atom at \(r>6\) left outside these selected
sets, its root mass satisfies
\[
 c-\frac{D}{r-s_*}
 \ \le\ \beta_{T,r}\ \le\ \frac{D}{r-s_*}.
 \tag{MC5}
\]
The lower bound is the selection rule. For the upper
bound use (MC2) after arbitrarily close finite
approximations to the selected common weight, with the
remaining atom as the extra lock. The condition
\(\beta_{T,r}<c\) follows from its moment and \(r>6\).
In particular for \(\sigma>0\) each remaining atom
at degree at least \(6+\sigma\) has mass at least
\(c\sigma/(6-s_*+\sigma)\), since \(D\le c(6-s_*)\).
There can be only finitely many such remaining scales
on each letter, uniformly in the profile once
\(\sigma\) is fixed.

The iteration can exclude a profile after its first
selection passes the synchronized packing test. On
each letter take atoms of masses \(1/100,1/50,41/300\)
at \(12,4,115/41\), respectively, and diffuse mass \(c\)
uniform on \([2,3]\). The first moment is one and
\(s_*=3\). At \(D_0=1/2\) the degree-twelve atoms are
selected, with total weight \(9/25\). Even adding all
degree-four atoms as extra locks in (MC2) gives only
\(9/25+1/2+3/50=23/25\).
Now \(D_1=7/50\) and the degree-four atoms have
\(u_{T,4}=11/75>7/50\). Their added weight is \(6/25\),
so \(D_2=-1/10\). At the first common weight \(9/25\)
the older unsupplemented repeat test cannot begin
converting those atoms. The largest admissible even
number of degree-four locks is eight, and
\(9/25+8(2/25)=1\).

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

For an explicit loss, freeze cheap positive-weight full
blocks on all three letters in the same graph and perform
the prime coloring by letter and short slot once. Let
\(E\) be the resulting edge count with exponent at least
\(2c-\eta_0\), where \(\eta_0\) is arbitrarily small.
Write \(M_i\) for the three outside capacities. For a fixed
pair \(i,j\) of letters, each edge is in its usual
bipartite projection, with no multiplicity by pair
uniqueness. Call an edge good for this pair if it has
at least \(\tau_{ij}=E^3/(10^4 M_i^2M_j^2)\) nondegenerate
opposite four-cycle anchors in that projection.
If at least \(E/10\) edges were bad, their projected
subgraph would have more than \(E_{\rm bad}\tau_{ij}\)
four-cycles. Indeed the ordered cycle count is at
least \(E_{\rm bad}^4/(M_i^2M_j^2)\) before degeneracies,
and the degenerate count is \(O(E_{\rm bad}(M_i+M_j))\)
of strictly smaller exponent for small parameters. This
contradicts the defining upper bound on the bad edges.

Thus at least \(0.7E\) edges have the stated lower bound
on anchors in all three projections. Average over the
at most \(E^3\) triples of anchors, one in each
projection. Some fixed triple has at least a constant
times \(E^7/(M_A M_B M_C)^4\) compatible edges. This
has exponent at least \(2c-7\eta_0-4\gamma-O(\delta)\)
because \(c_A+c_B+c_C=3c+\gamma\). Pigeonhole the three
complete cross ratios simultaneously at negligible
normalized cost as in (MS2). This gives multiplicative
separation in all three orientations on a single graph,
now with deficit \(\eta\le4\gamma\) up to arbitrarily
small additional loss.

Apply the [factor-energy argument](#freezing-factors-related-on-every-edge)
once to the three product relations, viewed as a single
coordinatewise product relation in \((\mathbb Q_{>0})^3\).
In each coordinate invert the output factor, so the
product of the three part vectors is \((1,1,1)\).
Bin the vector fibre sizes and sign-color each coordinate
of each part independently. A fixed vector product still
has only divisor-many representations, by applying the
integer divisor bound to every rational coordinate.
The same energy calculation therefore applies. If
\(v_i=c_i-\log K_i/(L\log X)\) is the vector dictionary cap
up to window loss, with \(K_i\) its fibre-size bin, a
matched pair on \(i,j\) has exponent at least
\(v_i+v_j-\eta-\Delta_{ij}\), where
\(\Delta_{ij}=c_i+c_j-2c\). Thus the output cap satisfies
\(v_k\ge v_i+v_j-2\eta-2\Delta_{ij}\) up to window losses.
Adding the inequalities for the two outputs other than
\(i\) gives \(v_i\le2\eta+c_i-c+\gamma\).

Fix the factor vectors on the two parts with the
smallest \(c_i-c+\gamma\). Every complete letter has
an opposite factor fixed, which gives its single-label
realization. The two chosen quantities sum to at most
\(8\gamma/3\), because all three sum to \(4\gamma\).
The final deficit is therefore at most
\(5\eta+8\gamma/3\le68\gamma/3\), up to arbitrarily
small additional loss. This gives the claimed near-minimum
complete-label graph for every tolerance \(\epsilon>0\).

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

## An explicit one-letter mesh discrepancy

At \(d=1/3,\ \Gamma=1/2\), let \(U,T\) be distinct letters
and let \(0<a\le1\). For every compact positive band
\([s,S]\) with \(s>0\), \(S<2s\), and positive nonatomic
\(T\)-mass in its interior, one has
\[
 \int \|r/a\|\,d\nu_U(r)>\frac{s}{11500}.
 \tag{DN2}
\]
In particular the left side is at least
\(\sup\operatorname{supp}\nu_{T,\mathrm{na}}/11500\),
by testing bands approaching the positive support endpoint.
Only one letter's phase discrepancy enters this inequality.

Work in the complete-label graph with fixed complete rows
on the remaining letter. For any fixed \(m\), its ordered
star count produces \(m\) distinct rows with common \(U\)
columns of normalized exponent arbitrarily close to
\(c=1/6\). In a common column with complete \(U\)-value \(u\),
the corresponding \(T\)-terms have the form
\(T_i=\varepsilon u+t_i\) with a common \(\varepsilon=\pm1\)
and distinct fixed offsets \(t_i\). In particular all pairwise
differences \(t_i-t_j\) are nonzero. By the support bound,
we may take \(s\le6\). Make harmless positive changes to
the band edges inside the original interval to obtain
continuity endpoints, keeping its lower bound \(s\), upper
bound \(S<2s\), and positive nonatomic mass \(H\).

Set \(\mu=0.083\). Take \(m\) fixed large enough that
\(H>\mu/m\). On each \(T_i\) choose a prescribed outside
block of roots in the band of normalized mass in
\([\mu/m,\mu/m+\delta+o(1)]\). Here the upper rounding
error comes from the maximum outside coordinate mass,
and \(\delta\) can be arbitrarily small. Fix the selected
roots. Their total freeze cost is at most
\(\mu+m\delta+o(1)\). The selected block positions are fixed
before the column values or their prime factors are chosen.

Put \(F=\prod_{i<j}(t_i-t_j)\) and choose
\(0<\xi<1/(10000m)\). Discard columns for which some outside
\(T_i\)-product has radical overlap of mass \(\xi\) or more
with \(F\). The fixed-row divisor-allocation estimate and
pair uniqueness bound these columns by
\(X^{(c+\delta-\xi+o(1))L}\); the incoming star loss and
\(\delta\) are taken smaller after \(\xi\) is fixed.
Strip all primes of \(F\) from the recorded roots.
Their bases are squarefree in the fixed bounded degree
band. The total stripped base mass is at most
\(m\xi+o(1)<0.0001+o(1)\). The surviving full products
\(n_i\mid T_i\) are pairwise coprime, since every common
prime divisor of \(T_i,T_j\) divides \(t_i-t_j\).

Let \(n=\prod_i n_i\). The congruences
\(\varepsilon u+t_i\equiv0\bmod n_i\) define a fixed
unit residue of \(u\bmod n\), since each witness is
primitive. Choose \(kL\to s\), below all selected degrees.
All modulus valuations satisfy \(k\le e<2k\) for small
\(L\), by squarefreeness and \(S<2s\). Thus the full shallow
record is \(Q_n=\operatorname{rad}(n)\). Its normalized
lower and upper masses, and the numerical weight of \(n\),
can use
\[
 b=0.0829,\qquad b_Q\le0.0831,\qquad w\ge0.0829s.
\]
Here take \(m\delta<0.0001\), with the usual limiting
errors absorbed by the strict margins. All estimations
are uniform in the fixed root values. Equivalently pass
to a subsequence of the modulus weights if (QD2) is stated
with a limiting \(w\).

On \(U\), round every scaled degree to a nearest
nonnegative multiple \(ia\) of \(a\) and put the \(i\)-th
power of its original root into \(z\). Take \(q=\lfloor
a/L\rfloor\), leaving a wider rational multiplier \(v\)
from the signed rounding errors. Fix the exceptional
coefficient as part of the box and absorb it into the
unit congruence target. All primes of \(z,v\) come from
\(u\), so both are units modulo \(n\).
The first-moment and mass bounds give
\[
 h\le\frac1a+\frac16\le\frac{7}{6a},\qquad
 E\le a\int \|r/a\|\,d\nu_U(r).
\]
Indeed \(i\le r/a+1/2\), and rounding \(ia\) to \(iqL\)
has total numerical cost \(o(1)\), since \(L\sum i\beta_e=o(1)\)
for each fixed \(a\). The bounded continuous error
\(\|r/a\|\) converges by weak convergence, even at nearest
rounding ties.

Suppose the discrepancy is at most \(s/11500\).
Apply (QD2) with \(\alpha=1/2\),
\(M=\lceil113/a\rceil\), and its short exponent parameter
\(r_0=as/113\). Since \(s\le6\), \(r_0<qL\) for small
\(L\). The Taylor order bound is
\(\alpha(M+1)r_0\le115s/226<s\). Also
\[
 B=hr_0+EM
 \le\left(\frac7{678}+\frac{114}{11500}\right)s
 <0.020725\,s .
\]
The determinant threshold satisfies
\(b\alpha(M+1)r_0/2\ge (0.0829/4)s=0.020725s\).
Likewise \(4B<0.0829s\le w\), giving the reconstruction
condition. The number of grouped pairs has normalized
cost at most \(0.0831\), by the exact shallow record.
Each pair gives \(u\), hence its column. Counting the
recorded target roots as well gives normalized exponent
at most \(0.1662<c\). The bad-overlap columns also have
a fixed saving. This contradicts the near-full star.

The constant depends on the lower scale of the supplying band.
Combining it with (VM7) also gives a uniform positive
one-letter lower bound, independent of the profile and of
\(a\in(0,1]\): for every \(U\), at least one of the two
other letters has uniformly positive diffuse numerical
weight, and hence its nonatomic support endpoint is
uniformly bounded below after dividing by its mass \(c\).

### Explicit weight outside any finite atom list

Put
\[
 \theta_U=\frac1{11500}
   \max_{T\ne U}\sup\operatorname{supp}\nu_{T,\rm na},
 \qquad H=\left\lceil\frac{2c}{\theta_U}\right\rceil.
\]
For any list of \(m\ge1\) positive atoms of \(U\), let
\(W\) denote their combined numerical weight. Then
\[
 1-W\ \ge\ \frac{\theta_U}{4H^{m-1}} .
\]
If \(W<c\), the inequality is immediate. Otherwise at
least one selected scale \(r_*\) is at least one, since
the total atomic mass is \(c\). Put \(d=\lceil r_*\rceil\).
Simultaneous Dirichlet approximation of the \(m-1\)
numbers \(d r_i/r_*\) for the other selected scales gives
some \(1\le k\le H^{m-1}\) for which
\(\|kd r_i/r_*\|\le1/H\). This follows by putting the
\(H^{m-1}+1\) consecutive multiples of their fractional
vectors into \(H^{m-1}\) half-open cubes. When \(m=1\)
take \(k=1\). The mesh \(a=r_*/(dk)\le1\) aligns
\(r_*\) exactly. The selected phase contribution is at
most \(c/H\le\theta_U/2\), and the unselected contribution
is at most \((dk/r_*)(1-W)\le2H^{m-1}(1-W)\).
Apply the preceding lower bound \(\theta_U\).

For one selected scale \(r\), no approximation is needed:
choose \(a=r/\lceil r\rceil\) to get the sharper
\(1-W\ge (r/\lceil r\rceil)\theta_U\).

### A consequence for finite lists of atoms

Put \(\theta=\max_{T\ne U}\sup\operatorname{supp}
\nu_{T,\mathrm{na}}/11500>0\), and let \(H\) be an integer
at least \(2c/\theta\). Given a list of at most \(p\) positive
atoms on \(U\), let \(E\) be the numerical weight outside
it. If the selected weight is at least \(c\), some selected
scale \(r_*\) is at least one, since the selected mass is
at most \(c\). Put \(d_*=\lceil r_*\rceil\). Simultaneous
Dirichlet approximation to the other at most \(p-1\)
numbers \(d_*r_i/r_*\) gives \(1\le k\le H^{p-1}\)
for which the mesh \(a=r_*/(d_*k)\le1\) aligns \(r_*\)
and has discrepancy at most \(1/H\) at the other
selected scales. The omitted part contributes at most
\(E/a\le2 H^{p-1}E\). By (DN2),
\[
 E\ge\frac{\theta}{4 H^{p-1}} .
 \tag{DN3}
\]
If the selected weight is less than \(c\), the same
inequality follows from \(E>1-c\) and the support bound
\(\theta\le6/11500\). Thus this is an explicit uniform
distance from any such finite atomic truncation, given
the support of the opposite diffuse measures.

### A one-letter bound for several shifted meshes

Work at \(d=1/3,\ \Gamma=1/2\), and let \(c=1/6\). Let \(U\) be one letter and let
\[
 s_*=\max_{T\ne U}\sup\operatorname{supp}\nu_{T,\mathrm{na}}
\]
be the larger opposite nonatomic support endpoint. The [complete-label realization](#complete-values-determined-by-single-outside-labels) and the [unequal-height modulus count (UW)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height) give the following one-letter obstruction.

Write the complete \(U\)-value in each original box as
\[
 u=Y\prod_{i=1}^m Z_i^{q_i},\qquad
 H(Z_i)\le X^{(H_i+o(1))L},\qquad H(Y)\le X^{E+o(1)} .
\]
The grouping must be fixed at the level of the original degree coordinates: each \(Z_i\) is a signed product of integer powers of designated original roots, and \(Y\) contains the signed exponent differences. Powers and assignments may depend on the fixed box. In particular the same original root may appear in more than one grouped coordinate. Set
\[
 C_i=\max(1,64H_i),\qquad B_m=m^2+m,\qquad \ell=B_m .
\]
Write \(c_m^*>0\) for the (UW) wider-height constant in dimension \(m\). Then
\[
 \boxed{\quad E\prod_{i=1}^m C_i
           \ \ge\ \frac{3c_m^* s_*}{320\,\ell^{m+1}} .\quad}
 \tag{AF1}
\]
This controls the numerical height of the correction for a fixed number of grouped roots.

#### Count in the complete-label graph

Fix an opposite target letter \(T\), and take a compact positive continuity band \([s,S]\) on it with nonatomic mass \(H>0\) and \(S<2s\). Choose a fixed number \(n_0\) of rows with \(n_0H>1/40\). The complete-label graph supplies those fixed distinct complete rows with common \(U\)-column exponent arbitrarily close to \(c\), by the usual ordered-star count. On each of those rows the complete target has form
\[
 T_j=\epsilon u+t_j
\]
for fixed distinct integers \(t_j\) and the same sign \(\epsilon\in\{1,-1\}\). These statements use only exact single-label values and hold for each sign orientation of the primitive signed relation.

Select blocks from the original outside roots of \(T_j\) in \([s,S]\), making their aggregate root mass \(1/40+o(1)\) up to an arbitrarily small outside-window rounding error. This is possible because the band contributes nonatomic mass \(H\) to every row and the outside coordinate size cap tends to zero with the window loss. Keep the original freeze cost at most \(1/32\). Let \(D=\prod_{i<j}(t_i-t_j)\), a fixed nonzero integer of polynomial height. Discard columns in which an outside target product shares radical of mass at least \(\xi\) with \(D\). The [divisor allocation estimate](#counting-shared-prime-divisors) gives a strict exponent saving on those columns, taking the rectangle and window losses much smaller than the fixed \(\xi>0\). Strip primes dividing \(D\) from the \(n_0\) selected outside blocks, with total removed root mass at most \(1/160+o(1)\) by choosing \(\xi\) small enough. The bounded-degree outside roots are squarefree, and the stripped blocks are pairwise coprime because a common prime of two target terms divides their fixed difference. Their product \(n\) is a full modulus with \(kL\to\kappa=s\) and numerical weight \(w\) satisfying
\[
  \frac{3s}{160}\le w\le\frac{s}{16}.
\]
All constants have harmless strict margins if the band is chosen inside \((s,2s)\). The selected large factors cost at most \(1/32\), even after stripping, since the originally selected root values enumerate the stripped ones.

The target congruences impose a single residue for \(u\bmod n\). It is a unit, since each of the complete triples is primitive and \(n\) comes from the opposite target terms. The factor \(s\) in \(\kappa\) cancels with the lower weight bound in all the narrow caps: \(C_i\) strictly exceed the required ratios \(\kappa H_i/w\) if \(H_i>0\), and a constant cap also suffices for \(H_i=0\).

Choose one grouped tuple for every distinct complete \(u\).
A product of \(\ell\) grouped tuples satisfies the congruence
with the common residue raised to \(\ell\). Its exact monomial
is the product of the \(\ell\) complete positive integers.
For a fixed such product, the ordered factorizations into
\(\ell\) complete values number \(X^{o_X(1)}\) at fixed
truncation. This bounds both the exact fibres of (UW) after
amplification and the multiplicity of the product map, even
when the grouped roots have overlapping prime supports.

If the left side of (AF1) is at most \(3c_m^*s/(320\ell^{m+1})\), the amplified wider-height condition holds strictly, since
\[
 \frac{\ell E}{w}\prod_i(\ell C_i)
     =\frac{\ell^{m+1}E\prod_iC_i}{w}\le c_m^*/2.
\]
The chart exponent after taking the \(\ell\)-th root is at most \(B_mw/(\ell\kappa)\le1/16\). Including the locks gives count exponent at most \(1/32+1/16<c\). Each complete value gives its outside column up to divisor-many representations. This contradicts the near-full columns; the discarded bad-overlap columns have their own strict saving. Let \(s\) approach \(s_*\) from below to obtain (AF1).

#### Shifted meshes

For \(p\) phases \(\phi_1,\ldots,\phi_p\in[0,1)\) and \(0<a\le1\), define
\[
 D_p(a;\phi)=\int \min_j \|r/a-\phi_j\|\,d\nu_U(r).
\]
Assign each degree to the nearest point \(a(i+\phi_j)\), with a deterministic rule for ties. Use a first grouped rational root containing each root raised to \(i\), and \(p\) further grouped roots, one the product of roots assigned to each phase. Their exponents are \(\lfloor a/L\rfloor\) and \(\lfloor a\phi_j/L\rfloor\), respectively. Multipliers \(i\) may be negative, which is allowed for the first rational root. Its height coefficient is at most \(2/a\), since \(|i|\le r/a+3/2\), and the total moment and root mass are one and \(1/3\). Each phase root has coefficient at most \(1/3\). Use a dummy root \(1\) with exponent \(1\) for any phase whose assigned exponent is zero. The remainder has numerical weight at most \(aD_p(a;\phi)+o(1)\), since the minimum distance function is bounded and continuous and the floor error contributes \(o(1)\) by the first-moment bound. Thus (AF1) gives
\[
 D_p(a;\phi)\ \ge\
 \frac{3c_{p+1}^*}{320\cdot128\cdot32^p\,B_{p+1}^{p+2}}\,s_*.
 \tag{AF3}
\]
The constant is uniform in the phases and in the mesh spacing.

### Diffuse tails and mesh discrepancy near minimum mass

There are absolute \(\gamma_1,\alpha_0>0\) such that
each letter of a balanced \(d=1/3\) endpoint with
\(\Gamma-1/2\le\gamma_1\) has atomic numerical weight
greater than \(\alpha_0\). The argument uses stability
of the fixed-row count under the small graph deficit.

First, for every \(s>6\), the nonatomic mass of any letter
above \(s\) tends uniformly to zero as \(\Gamma\to1/2\).
Otherwise the (VM2) tail moment bound allows one to
enclose a fixed positive amount of this mass in a fixed
bounded interval with lower endpoint greater than six.
Use the [near-minimum complete-label reduction](#the-complete-label-reduction-close-to-minimum-mass)
and choose a fixed \(m\) large enough for the high-band
locks in (CT1). Their cost and weight margins depend
only on the positive amount of mass retained, the band,
and \(m\). For small enough \(\gamma,\delta\), the raw
outside capacities are at most \(c+\gamma+\delta\), so
the bad-overlap branch retains a strict deficit from
\(c\) as well. The fixed complete additive rows have
common-column exponent arbitrarily close to \(c\)
as \(\gamma\to0\), and give the same contradiction as
at minimum mass. More quantitatively, for fixed \(s>6\) and finite \(S\)
the nonatomic mass \(H\) in \([s,S]\) is
\(O_{s,S}(\sqrt{\gamma})\). The recursions in the
[near-minimum complete-label reduction](#the-complete-label-reduction-close-to-minimum-mass)
give complete-label edge deficit at most
\(68\gamma/3\), up to arbitrarily small additional
loss. Thus the common column exponent after selecting
\(m\) distinct rows is at least
\(c-((74/3)m-1)\gamma\) up to such loss.
Take \(z=(c+1/s)/2\), \(m=\max\{2,\lceil z/H\rceil+1\}\)
and overlap threshold \(\xi=25m\gamma\) in (CT1).
The bad-overlap branch has a strict column deficit.
The good branch also has one if
\(z<c-25m\gamma\) and \(sz>1+25m^2 S\gamma\).
Both inequalities hold for small \(\gamma\) if \(H\)
is larger than a sufficiently large multiple of
\(\sqrt{\gamma}\). The fixed number of rows for each
profile is chosen before taking the height limit.
Together with (VM2) this implies
\[
 \int(r-6)_+\,d\nu_{U,\mathrm{na}}(r)\longrightarrow0
 \quad\text{uniformly as }\Gamma\to1/2 .
 \tag{DN4}
\]

Second, the mesh proof of (DN2) works close to minimum
mass when the supplying band has \(a_0\le s<S\le7\),
\(S<2s\), and nonatomic mass at least a fixed \(H_0>0\).
Indeed choose \(m\) using \(H_0\) and then the fixed
\(\xi\) for the overlap pruning. Take the complete-label
deficit and \(\gamma,\delta\) small enough that the
\(m\)-star has common-column exponent greater than
both \(0.1662\) and \(c+\gamma+\delta-\xi\).
The one-letter height bounds still have \(l_U=1/3\)
and \(h_U\le1\). The determinant margins are unchanged,
and \(s\le7\) still gives \(r_0<a\).
For any letter \(U\), a supplying band on an opposite
letter exists with \(a_0,H_0\) uniform by (VM7).
To obtain it, intersect its uniformly positive
nonatomic mass in a compact positive range with
\([0,7]\), using the uniform mass bound above seven
just proved. Then select a band from a fixed finite
cover of that compact range by intervals of ratio
less than two. Hence \(D_U(1)=\int\|r\|\,d\nu_U(r)\)
is greater than \(a_0/11500\) uniformly for small
\(\gamma\).

On the atomic part, \(\|r\|\le r\), and on the
nonatomic part, \(\|r\|\le|r-6|\). Write \(e_U\)
for the escaped moment and \(\alpha_U\) for the
atomic numerical weight. We obtain
\[
 D_U(1)\le 2\alpha_U+6c_U-h_U+e_U+
                  2\int(r-6)_+\,d\nu_{U,\mathrm{na}}(r).
\]
The last terms beyond \(2\alpha_U\) tend uniformly
to zero by (VM1), (VM2), and (DN4), while \(D_U(1)\)
has a positive lower bound. This proves the asserted
atomic moment bound for each letter. Together with
(DN1) it bounds both types of numerical weight away
from zero near minimum mass.

## A count saving from a high-degree band on one letter

A limiting box profile with
\(\lvert c_W-1/6\rvert\le10^{-5}\) on all three letters
and nonatomic mass greater than \(2/25\) in
\([34/5,7]\) on one letter satisfies
\[
 d<\frac13-\frac1{10000}.
 \tag{DN5}
\]
There are no numerical-weight or lock hypotheses on
the other two letters.

Suppose \(d=2c-\varepsilon\) with
\(0\le\varepsilon\le10^{-4}\), and write
\(\rho=10^{-5}\). The pair-fibre argument, with
\(c_W>0\) on every letter, thins to a linear graph
of exponent \(d\) up to an arbitrarily small loss.
In the [simultaneous anchor selection](#the-complete-label-reduction-close-to-minimum-mass)
the outside part capacities sum to at most \(3c+3\rho\),
so the initial deficit from \(2c\) after separation is
at most \(7\varepsilon+12\rho\) up to the same losses.
The degenerate-cycle contribution stays strictly smaller,
since \(d\) is within \(10^{-4}\) of \(2c\) and the
part exponents are within \(\rho\) of \(c\).

Only the two letters opposite the high band need full
labels. Freeze just the factor vector on the high-band
letter in the [vector energy argument](#the-complete-label-reduction-close-to-minimum-mass).
For initial deficit \(\eta\) its vector dictionary exponent
is at most \(2\eta+4\rho\), since both pair excesses
\(c_i+c_j-2c\) are at most \(2\rho\).
The remaining graph has each opposite complete value
determined by its own label and count exponent at least
\(2c-21\varepsilon-40\rho\), up to arbitrarily small loss.

Each opposite complete value has divisor-many outside
representations at fixed truncation. Keeping one label
per complete value on each opposite part, chosen to
retain the largest edge count, loses only a negligible
normalized factor. Take two distinct fixed complete
rows in this thinned graph. The star count
gives common columns with exponent at least
\(c-42\varepsilon-83\rho\ge c-0.00503\) up to loss.
In each column the two high-band terms differ by a
fixed nonzero shift. Apply the fixed-row congruence
argument (CT1) with total lock mass \(z=4/25\), by
choosing \(z/2\) from the indicated band on each row.
After discarding columns with radical overlap of mass
at least \(\xi=0.0051\) with the fixed shift, their
exponent is at most \(c+\rho-\xi=c-0.00509\) up to
arbitrary window loss. Strip shift primes from the
two recorded blocks. The combined numerical weight
of their remaining coprime moduli is at least
\[
 \frac{34}{5}\frac4{25}-2\cdot7\cdot0.0051=1.0166>1 .
\]
Thus the common complete column parameter is uniquely
reconstructed for each fixing. Their exponent is at
most \(z=0.16\) up to the rounding losses, also strictly
below the star exponent. This disproves the supposed
range for \(d\), proving (DN5).

## An all-scale progression bound for complete letters

Work at \(d=1/3,\ \Gamma=1/2\), with \(c=1/6\) and \(n=X^L\). Fix a column letter \(U\) and put \(s_*=\max_{T\ne U}\sup\operatorname{supp}\nu_{T,\rm na}\). The complete-label realization gives simple additive graphs whose edge exponent is arbitrarily close to \(2c\), every part cap is arbitrarily close to \(c\), and the exact values of an edge are determined by its labels. Keep one outside representation per complete value, at fixed-truncation divisor loss. The complete \(U\)-vertices can be regularized to have degrees at least \(n^{c-\delta_0}\), where \(\delta_0>0\) is arbitrarily small. Remove the vertices below this cutoff after taking the input graph deficit and part-cap error sufficiently smaller than \(\delta_0\). The removed edges have a strict saving from the original count exponent by the \(U\)-part cap.

### Interval and residue conditions

For a fixed \(0<s<s_*\), a fixed \(\tau\ge0\), and any prescribed exponent tolerance \(\epsilon>0\), take the regularized graph with losses sufficiently small in terms of these parameters and the profile. Uniformly for every interval \(I\) and every integer progression \(r+Q\mathbb Z\) with \(Q\ge1\) and \(|I|/Q\le X^\tau\), the regularized complete column set obeys
\[
 \bigl|\mathcal U\cap I\cap(r+Q\mathbb Z)\bigr|
 \ \le\ X^{(\min(c,\tau/s)+\epsilon+o(1))L}.
 \tag{AC1}
\]
The \(o(1)\) term has the usual small-budget and fixed-truncation order. The case \(Q>2X\) is immediate on the ambient complete-value range, so the numerical height of \(Q\) can be assumed polynomial. The exponent bound is uniform in \(I,r,Q\), even when they depend on the height. In particular one can take a finite grid of \(\tau\)'s and use the same graph losses to obtain, with any additional prescribed tolerance, the statement for all \(0\le\tau\le1\) simultaneously. Taking \(s\uparrow s_*\) gives the limiting coefficient \(1/s_*\).

### Counting good stars over an arbitrarily small column subset

Let \(\mathcal S=\mathcal U\cap I\cap(r+Q\mathbb Z)\). A compact positive continuity band \([s,S]\) on an opposite target \(T\) has positive nonatomic mass \(H\), for some \(S>s\). To see this when \(s\) is in a gap of the support, use a narrower band farther above \(s\) and use \(s\) as its lower degree bound. Choose \(p>\tau/s\) with \(p<\min(c,\tau/s+\epsilon/3)\), when \(\tau/s<c\) and the total part-cap bound is insufficient. Choose a fixed \(m\) so that \(mH>p\). The rows lie on the remaining letter. Each fixed \(u\in\mathcal S\) belongs to at least \(n^{m(c-\delta_0)+o(1)}\) ordered stars of distinct rows. For that \(u\), as a row varies its complete target \(t_i=\varepsilon u+a_i\) determines the row, where \(\varepsilon\) is fixed and the complete row offsets \(a_i\) are distinct.

For each slot let \(P_i\) be the original outside
radical of its target. A star is bad at slot \(i\) if
\(P_i\) shares normalized radical mass at least \(\xi\)
with \(Q\prod_{j\ne i}P_j\). Take \(\xi>0\) fixed with
\(s(p-m\xi)>\tau\). For a fixed column \(u\), hold all
row slots except \(i\) fixed. The product with which
\(P_i\) is compared is then fixed and has polynomial
height. The [divisor allocation estimate](#counting-shared-prime-divisors) counts at most \(n^{c+\delta-\xi+o(1)}\) offending outside target labels, where \(\delta\) includes the outside-window cap error. A complete target and the fixed \(u\) determine its row, so the whole bad-star count for this \(u\) is at most \(n^{m(c+\delta)-\xi+o(1)}\), times a constant depending on \(m\). Take \(m(\delta_0+\delta)\ll\xi\). This comparison gives a saving relative to the good-star lower bound for each \(u\) separately, regardless of \(|\mathcal S|\).

Sum over \(u\in\mathcal S\) and average over the at most \(n^{m(c+\delta)+o(1)}\) row tuples. Some fixed tuple has at least
\[
 |\mathcal S|\, n^{-m(\delta_0+\delta)+o(1)}
\]
good columns from \(\mathcal S\). This order of counting is needed when \(\mathcal S\) may have cardinality of substantially smaller exponent than \(c\).

### CRT reconstruction within the given interval

On the fixed rows, record outside blocks from the supplying band, of combined normalized mass \(p\) up to arbitrarily small errors. Since each outside coordinate cap can be made arbitrarily small, the selected mass can reach \(p/m\) on each row at arbitrarily small overshoot. Their positions are prescribed from the box, and the original root values have at most \(n^{p+o(1)}\) records. The original roots in the bounded scale band are squarefree. Strip all primes dividing \(Q\) and all primes occurring in two or more recorded blocks. On the good columns the stripped root mass is at most \(m\xi\), using the single union tested at each slot. The surviving full divisors have pairwise disjoint primes and are coprime to \(Q\), and their product has numerical weight at least \(s(p-m\xi)\), up to errors.

Each such divisor of \(t_i\) gives a fixed affine residue for \(u\). Together with \(u\equiv r\bmod Q\) they give a progression with step \(Q X^{s(p-m\xi)+o(1)}>|I|\), so there is at most one column for each record. Combining with the good-star lower bound gives \(|\mathcal S|\le n^{p+m(\delta_0+\delta)+o(1)}\). Take the losses small compared with \(\epsilon\). The total-part cap handles \(\tau/s+\epsilon\ge c\). The comparison also works for \(\tau=0\), or directly an ordinary progression with \(|I|\le Q\) meets the interval in at most two values. This proves (AC1).

### Consequences for ordinary and generalized progressions

**Ordinary progressions.** A progression of any polynomial-height step and \(X^{O(L)}\) terms has intersection exponent \(o(1)\) with the regularized \(U\)-part, in the small-budget limit. More quantitatively (AC1) gives exponent at most \(\zeta/s_*+o(1)\) for a progression of \(X^\zeta\) terms, up to arbitrary extraction loss. At the symmetric \(c\delta_3+c\operatorname{Unif}[2,4]\) profile, a subset of fixed normalized exponent \(v>0\) therefore requires any containing progression to have at least \(X^{4v-o(1)}\) terms. In particular an interval containing such a subset has at least that numerical width.

**Generalized progressions.** A generalized arithmetic progression of fixed rank \(R\) and coefficient-box volume at most \(X^{\kappa L}\), with \(\kappa\) fixed, meets the regularized \(U\)-part in at most \(X^{((1-1/R)\kappa+o(1))L}\) values, with arbitrarily small input losses. Slice along a side of maximal coefficient exponent. Its exponent is at least the actual volume exponent divided by \(R\), and each slice is an ordinary progression of \(X^{O(L)}\) terms and contributes exponent \(o(1)\). If the selected step exceeds \(2X\), its slice contributes at most one complete value. This works irrespective of duplicate representations. More generally the same proof saves the exponent of the largest side for any presentation for which the number of other coefficient choices has the corresponding bound. For each fixed \(R\), a generalized progression of normalized volume \(c+\delta\) therefore excludes a near-full column subset whenever \((1-1/R)(c+\delta)<c\).

**Previously frozen column factors.** Suppose the selected columns in addition have a prescribed full factor \(F\) from \(U\), after freezing roots of normalized cost \(b\). Fix their values and apply the same proof with the known residue \(u\equiv0\bmod F\) in addition to \(u\equiv r\bmod Q\). The product of the opposite row locks is coprime to \(F\) by primitivity, so no new strip loss is needed for \(F\). If \(|I|/\operatorname{lcm}(F,Q)\le X^\tau\) for the frozen records, the resulting upper column exponent is \(b+\tau/s_*\), up to the same limiting errors. In particular, with \(Q=1\), a fixed reliable numerical weight \(w=\log_X F\) reduces the interval exponent in (AC1) from \(h=\log_X|I|\) to \((h-w)_+\), with only the standalone record cost \(b\).

### Modular progressions with a short parameter range

Let \(\mathcal S\) consist of the regularized complete \(U\)-values in an interval of numerical width at most \(X^h\), for which
\[
 u\equiv A+d j\pmod M,\qquad 0\le j<J\le X^\zeta .
\]
Here \(j\) is integral, \(M\) has polynomial height, and \(M\ge X^\beta\). For fixed exponents in the small-budget limit, the normalized exponent of \(\mathcal S\) is at most
\[
 \min\!\left(c,\;\frac{(h+\zeta-\beta)_++\zeta}{s_*}\right)
 \tag{AC2}
\]
up to arbitrary graph and extraction losses. This bound is uniform in \(A,d,M\), including noninvertible steps and repeated representations. In particular a modular progression with \(X^{O(L)}\) terms and modulus of numerical weight at least \(1-o(1)\) contains only \(X^{o(L)}\) regularized complete values in the limiting exponent sense.

To prove the claim, choose \(s<s_*\) and \(p>(h+\zeta-\beta)_+/s\), with an arbitrarily small fixed margin. If \(p\ge c\), the total part cap suffices. Choose finitely many opposite supplying rows and strip all mutual overlaps and primes of \(M\) by good-star averaging on \(\mathcal S\), as in (AC1). Each recorded lock tuple costs at most \(n^{p+o(1)}\), has a fixed CRT condition \(u\equiv K\bmod N\) with \(\gcd(M,N)=1\), and has weight \(\log_X N\ge sp-o(1)\) up to arbitrarily small stripping losses. We fix the ordinary column representation losses small enough that the selected rows retain \(\mathcal S\) up to a smaller fixed normalized exponent loss.

For one such record, all differences of admissible pairs \((u,j)\) lie in the lattice
\[
 \Lambda=\{(x,y)\in\mathbb Z^2:x\equiv0\pmod N,\ x-dy\equiv0\pmod M\}.
\]
It has index \(MN\), since any \(y\) permits one class of \(x\) modulo \(MN\) by CRT. But any determinant of two differences has absolute value at most \(2X^{h+\zeta}<MN\) for large height. Thus all admissible pairs lie on an affine line. If the line has constant \(j\), the same CRT modulus \(MN\) directly gives at most one \(u\) in the interval. If it has constant \(u\) the claim is immediate. Otherwise the primitive integral direction has nonzero integral \(j\)-component, and the full set of integral points on the line inside the \(j\)-range projects to an ordinary progression of at most \(J+1\) terms in \(u\). Apply (AC1) to this ordinary progression, at exponent \(\zeta/s_*\) up to arbitrary losses. That application is valid uniformly over records using graphs with the smaller of the tolerances needed by the two counts. Counting the records and allowing \(p\) to approach \((h+\zeta-\beta)_+/s\), then \(s\uparrow s_*\), proves (AC2).

## A uniform obstruction to bounded-rank additive progressions

For every fixed rank bound \(R\ge1\) there is
\(\delta_R>0\) with the following property at
\(d=1/3,\Gamma=1/2\). Within any one complete letter, a
generalized arithmetic progression of rank at most \(R\)
and coefficient-box volume at most
\(X^{(c+\delta_R)L}\) contains at most
\(X^{(c-\delta_R+o(1))L}\) distinct complete values,
where \(c=1/6\). The progression may depend on the
height and the box. The result applies to arbitrary
presentations, including ones that represent a value
more than once.

To choose parameters uniformly, use (DN1) and (MS11).
A fixed positive lower bound for the nonatomic moment,
together with nonatomic mass \(c\) and support in
\([0,6]\), supplies an absolute \(b_0>0\) and \(H>0\)
so that a compact continuity band within \([b_0,7]\)
has nonatomic mass at least \(H\) on every letter.
Indeed take a sufficiently small lower cutoff that its
complement carries at least half the moment, and then
divide by the upper degree bound to obtain a mass bound.
Let \(k_0\in(0,1]\) be the constant in the
[one-letter real strip count](#counting-distinct-terms-and-ratios-in-strips).
Fix \(a_0>0\) below \(\min(H/4,c/(4R))\), choose
\(B>64R/(k_0 b_0)\) with \(B\ge2\), and put
\(a_j=a_0B^{-j}\) for \(0\le j\le R+1\).

Write the progression as \(\{b+\sum_i d_i u_i:
0\le u_i<N_i\}\), and put \(n=X^L\). One of the
\(R+1\) intervals \([a_{j+1},a_j)\) contains no side
exponent \(\log N_i/\log n\). Thus each side has
exponent below \(a_{j+1}\) or at least \(a_j\).
Put \(p=a_j/4\), \(\tau=b_0p/4\), and
\(\xi=b_0p/28\).

If one side of exponent at least \(a_j\) has step
\(|d_i|\ge X^{1-\tau}\), hold the other coefficients
fixed. There are at most \(n^{c+\delta_R-a_j}\) such
slices. If \(|d_i|>2X\), each slice contains at most one
eligible value for large \(X\). Otherwise \(|d_i|\) has
polynomial height and every slice lies in one residue
class modulo it. Select outside roots of normalized
cost \(p\) in the supplying band. The divisor-allocation
bound discards at most \(n^{c-\xi+o(1)}\) values with
radical overlap of root mass at least \(\xi\) with
\(|d_i|\), up to arbitrarily small window losses.
Strip these fixed step primes on the other values.
There is weight at least \(b_0p-7\xi\) left in the
stripped full blocks. Their product with \(|d_i|\)
exceeds the complete height range since
\(1-\tau+b_0p-7\xi>1\). Each slice and block record
therefore contains at most one complete value.
The good exponent is at most \(c+\delta_R-a_j+p\).

If every side of exponent at least \(a_j\) has step
below \(X^{1-\tau}\), fix the sides of exponent below
\(a_{j+1}\). Their total cost is at most \(R a_{j+1}\).
Every remaining slice has diameter at most
\(R n^{c+\delta_R} X^{1-\tau}\). At sufficiently small
\(L\) this has relative width at most \(O(X^{-\tau/2})\)
among the complete values of the box, which have
numerical height tending to one. The strip estimate
bounds the exponent by
\[
 c-k_0\tau/2+R a_{j+1}<c-k_0\tau/4 .
\]
All margins are uniform once \(R\) is fixed, since
\(a_j\ge a_R>0\). Take \(\delta_R\) small relative to
these margins and absorb the graph and truncation
losses. This proves the assertion.

In particular a subset of normalized size \(c-o(1)\)
of any complete letter has additive doubling tending
to infinity. For a fixed doubling bound, the
Green–Ruzsa form of Freiman's theorem puts the subset
inside a proper progression of bounded rank and volume
a bounded multiple of its cardinality. Its finite
subgroup is trivial in the integers, and the progression
contradicts the count above. This conclusion applies
to low-degree profiles as well. The rank bound depends
on the doubling bound, so the argument permits
unbounded doubling with zero limiting exponent.
[Green–Ruzsa, Theorem 1.1](https://arxiv.org/pdf/math/0505198v2)

## Atomic common factors at the broad symmetric endpoint

Suppose all three limiting measures are
\[
 \nu_W=c\delta_3+c\,\operatorname{Unif}[2,4],\qquad c=1/6,
\]
at \(d=1/3,\Gamma=1/2\). The near-full complete-label graph excludes a common factor of positive numerical weight carried by the original roots of any one degree-three atom window. Consequently each full atom-window factor dictionary needs normalized exponent at least \(c/2\) on any arbitrarily near-full edge subset. The same lower bound holds for the normalized information of each atom factor outside an exponentially sparse family, and hence for each individual atom factor entropy.

### A modulus already common to the atom window

Suppose a factor \(F\) divides all complete \(A\)-values on an arbitrarily
near-full graph and is carried by their full atom-window factors. Assume
each nonzero valuation of \(F\) is \((3+o(1))/L\), uniformly across its
primes. For instance, this holds when \(F\) is the greatest common divisor
of two original full atom-window factors. Pass to a subsequence with
numerical weight \(W=\lim\log_XF>0\). For \(Q=\operatorname{rad}F\), the
normalized root mass is \(b=\lim\log_X Q/L=W/3>0\).

Since each \(A\)-value lies in residue zero modulo \(F\) in an interval of
numerical width at most one, (AC1) with any opposing supplying degree
\(s<4\) shows \(c\le(1-W)_+/s\). To apply it, regularize this near-full
graph on \(A\) after the common-factor selection, with arbitrarily small
additional losses. Its retained \(A\)-part has size exponent arbitrarily
close to \(c\), because the opposite part caps remain \(c\). Each opposite
profile has positive nonatomic mass above every \(s<4\). Its locks are
automatically coprime to \(F\) on every original edge by primitivity. In
the limit \(s\uparrow4\) this gives \(W\le1/3\), hence \(b\le1/9\).

On the other two letters retain the atom window and the diffuse interval \([3-6t,3+6t]\), where \(t=(101/200)b\) is its mass on each letter. Freeze all remaining original factors, at normalized cost \(2(c-t)\) up to arbitrarily small errors. Group the retained roots unweighted at degree \(q=\lfloor3/L\rfloor\) in the ratio of those two letters. Its narrow and wider heights have bounds
\[
 h=c+t\le401/1800,\qquad E=3t^2.
\]
For the wider bound one uses the equal positive and negative deviation moments from each identical retained profile; their sum on one letter is \(3t^2\), since the conditional diffuse radius is \(6t\). Negligible packed and endpoint errors can be frozen or carried in the correction. Restoring the fixed factors multiplies the unit target by a fixed unit, and the original ratio is congruent to a fixed sign modulo \(F\) on every edge. No modulus record charge is needed for this already common factor.

Take \(kL\to3\), below all selected valuations, and use the shallow residue classes modulo \(Q\), of normalized cost \(b\). To check the first one-root determinant, take \(M+1=\lceil4h/b\rceil\). Then
\[
 \frac{3}{M}\left(\frac b2-\frac{h}{M+1}\right)
 \ \ge\ \frac{3b^2}{16h}
 \ \ge\ \frac{5400}{6416}b^2
 \ >\ 3(101/200)^2 b^2 = E .
\]
In detail \(h/(M+1)\le b/4\) and \(M<4h/b\). All parameters are fixed after the positive limiting \(b\) is chosen. One can therefore take a rank fraction \(\alpha<1\) arbitrarily close to one and a normalized horizontal degree \(r_0\) close below \(3/(\alpha(M+1))\), so \(h r_0+ME< b\alpha(M+1)r_0/2\) and \(\alpha(M+1)r_0<3\). The bounded-degree Newton proof uses only rank drop, so no Hensel terminal bound is needed.

On each chart a nonzero auxiliary polynomial of fixed wider degree, with polynomial coefficient heights, is available. The Newton-breakpoint argument recovers the remaining outside factors of both free letters. Primitivity puts each counted prime on a single letter, and its grouped multiplier is one. A prime outside the support of the nonzero auxiliary coefficients can occur only at a valuation breakpoint. Each such exceptional position has arbitrarily small outside cap. The outside roots omitted from the interval are already frozen. Thus after fixing the chart and frozen records there are only subpower and arbitrarily small exponent choices of the two outside labels. Apply (EN2) to complete the triples. The resulting normalized upper bound is
\[
 2(c-t)+b=2c-\frac b{100}<2c .
 \tag{CA1}
\]
This contradicts an arbitrarily near-full edge count.
[The one-root determinant and valuation breakpoints](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#newton-breakpoints-in-the-residue-determinant)

### Sparsity of atom-window overlaps

Fix \(\xi>0\). Among regularized complete vertices on a letter,
ordered pairs whose atom-window factors have gcd of numerical
weight at least \(\xi\) have exponent at most \(2c-\epsilon_\xi\)
for some \(\epsilon_\xi>0\) when the graph and window losses are
sufficiently small. Otherwise there are \(n^{2c-o(1)}\) such
pairs, and a fixed first vertex has \(n^{c-o(1)}\) partners
by the part cap. Pigeonhole their gcd with the first vertex
at divisor-many cost. The resulting divisor is carried by
the exact atom-window factors of all selected partners, and
its prime valuations remain in the window at three.
Each partner has near-full graph degree, so the selected
edges have exponent \(2c-o(1)\), contradicting (CA1).

### Individual atomic information

Suppose some arbitrarily near-full edge family on the same profile has a
dictionary of full atom-window factors on \(A\) of exponent \(u<c/2\).
Thin to injective complete-label maps and prune low-degree \(A\)-vertices,
so every retained \(A\)-vertex still has degree exponent close to \(c\).
Prime color globally by letter: independently assign each prime a letter,
and retain triples for which all their primes have the corresponding
letter colors. This loses only \(X^{o_X(1)}\) at fixed budget, since every
triple has \(O(\log X/\log\log X)\) distinct primes. Factors from
different complete parts are then coprime. Regularize degrees after
coloring if necessary.

Use ordered pairs of distinct edges sharing a complete \(C\)-vertex. The count is at least \(n^{3c-o(1)}\), in normalized exponents with \(n=X^L\), by the degree second moment and subtraction of the diagonal. With signs absorbed in the entries this gives distinct pairs in
\[
 a_1+b_1=a_2+b_2 .
\]
The quadruple determines its \(C\)-vertex. On each of \(a_1,a_2\) mark the whole prescribed atom-window factor, of weight \(1/2+o(1)\) and dictionary cost at most \(u\). On both \(B\)-entries mark outside factors from one bounded positive diffuse band, with combined root mass \(z>0\), choosing \(z\) small enough that \(2u+z<c\). Their combined numerical weight is at least \(s z+o(1)>0\) for a fixed lower degree \(s>0\).

Quadruples in which the two \(B\)-outside records have radical overlap at
least \(\xi>0\) in normalized root mass are exponentially fewer: fixing
\(b_1\) bounds the possible \(b_2\) outside labels by \(n^{c-\xi+o(1)}\)
by divisor allocation, and after \(a_1\) is chosen the last entry is
determined. Choose all window losses much smaller than \(\xi\). In the
remaining quadruples strip shared \(B\) primes, losing at most
\(2S\xi+o(1)\) numerical weight, where \(S\) bounds their degrees. If the
full gcd of the two \(A\)-atom locks also has numerical weight below
\(\xi'\), keep just the larger of the two powers at each shared prime.
This loses at most \(\xi'\) numerical weight. The surviving full divisors
are pairwise coprime. Take \(\xi,\xi'\) so small that their total weight
is still strictly greater than one. The records counted before this
deterministic stripping have total exponent at most \(2u+z<c\).

Use the fixed-record residue packing inequality (MS17) for the four-term
equation, with its numerical height at most \(X^{1+o(1)}\). If the product
of the four coprime moduli exceeds twice the ambient height, its number of
tuples for fixed records is at most the product of the four record-fibre
sizes to power \(2/3\). Summing by Hölder gives normalized exponent at
most \((8c+2u+z)/3<3c\). The quadruple count at exponent \(3c-o(1)\) thus
persists where the \(A\)-atom gcd has numerical weight at least \(\xi'\).

For each ordered \(A\) pair there are at most \(|\mathcal B|\) completions, since one of the \(B\)-entries determines the other. The \(n^{3c-o(1)}\) count therefore gives \(n^{2c-o(1)}\) ordered
\(A\) pairs with atom gcd of weight at least \(\xi'\), contradicting
the overlap bound.

All margins in the extraction can be fixed before taking the limits: choose \(u<c/2\), then \(z,\xi,\xi'\), and then the initial near-full graph and outside-window tolerances. For the common-factor obstruction one can pass to a subsequence with \(b\) in the compact positive range supplied by \(\xi'\), and take the extraction tolerances smaller than the resulting fixed saving. Finally, if atom-factor pointwise information were below some fixed \(u<c/2\) on a subset of edges with exponent arbitrarily close to \(2c\), that subset would have a dictionary of exponent at most \(u\), by the definition of information. The above contradiction rules this out. For the original near-full graph with uniform edge measure, let
\(F_T\) denote the full atom-window factor on letter \(T\).
Taking expectations and then \(u\uparrow c/2\) gives
\[
 \liminf\frac{H(F_T)}{\log n}\ge\frac c2=\frac1{12}
 \qquad(T=A,B,C).
 \tag{CA2}
\]

## Repeated atomic records at the broad symmetric endpoint

Assume the minimum nonatomic endpoint with
\[
 d=\frac13,\qquad c=\frac16,\qquad
 \nu_W=c\delta_3+c\,\operatorname{Unif}[2,4]\quad (W=A,B,C).
\]
Use the near-full complete-label graphs, outside coordinate caps and
atom windows of the endpoint argument. Put \(n=X^L\). A prescribed
atomic record on one complete part is a product of chosen prime powers
from its shrinking degree-three window, where a chosen prime is kept
with its assigned full exponent. The chosen divisor can depend on the
complete vertex, whose original representation is fixed. Its
numerical weight is \(\log_X F\). We consider finitely many such
record types, and can put each record's weight into a fixed small
interval before applying the proof.

For type \(i\), let \(w_i\) be the limiting numerical weight and let
\(\alpha_i\ge0\) be an upper dictionary exponent on an arbitrarily
near-full edge family. Types can belong to the same complete part.
For every fixed nonnegative integral multiplicity vector
\((m_i)\), one has
\[
 \sum_i m_i\alpha_i\ \ge\
 \min\left\{c,\frac{\sum_i m_i w_i-1/3}{4}\right\}.
 \tag{RA1}
\]
For the pointwise version use
\(I_i=-\log\Pr(F_i=\text{actual})/\log n\), with probabilities
from the uniform edge distribution. Formula (RA1) holds with
\(\alpha_i\) replaced by \(I_i\), up to any fixed positive tolerance,
outside an exponentially sparse edge family. The corresponding
inequality for \(\sum_i m_i H(F_i)/\log n\) follows by expectation,
with a limit inferior at the exact endpoint.

Use the [sparsity of atom-window overlaps](#sparsity-of-atom-window-overlaps).
The same pair bound holds for subfactors and includes repeated
vertices, since diagonal pairs have normalized exponent at most \(c\).

### Counting with repeated records

Write \(P=\sum_i m_i\alpha_i\) and \(W=\sum_i m_iw_i\).
A strict failure gives \(P<c\) and \(W+4(c-P)>1\). Fix \(s<4\)
and \(0<z<c-P\) with \(W+sz>1\). Pick a fixed compact nonatomic
band above \(s\) close to four on a supplying letter. Use as many
separate copies of this letter as necessary to provide outside
lock cost \(z+o(1)\) in the band. The rounding loss tends to zero
with the outside-coordinate cap. Pick one separate slot per
occurrence of each prescribed atom type. Add empty slots so the
counts on all signed complete parts are even, the total count
\(m\ge4\), and there is an empty slot on a different part from
the supplying letter. Signs on the slots are chosen to obtain
the balanced mixed additive relation.

Globally color primes by their complete letter before taking the
mixed moments. This loses only a subpower factor at fixed budget.
Prune low-degree complete vertices after this selection and fix
one admissible original representation and outside label per
complete value. Deleting unused vertices is harmless.
The common large Fourier spectrum estimate in *Rich mixed additive
relations* yields \(n^{(m-1)c-o(1)}\) ordered tuples in this
balanced relation, with all graph losses made sufficiently small
after the fixed margins, multiplicities, and thresholds.

Exclude tuples where any pair of atomic slots of the same letter
has full atom gcd of numerical weight at least \(\xi>0\). Fix the
offending two complete values, then fix all but one of the other
slots. The last value is determined by the signed additive
relation. The preceding pair bound and the complete part caps
give an exponent strictly below \((m-1)c\). Thus all these
exclusions have a fixed positive deficit once graph losses
are sufficiently small.

For an outside slot, exclude the tuples where its radical shares
normalized root mass at least \(\xi'>0\) with the union of the
other locks on its letter. Fix all complete values except this
slot and an empty slot of a different letter. The union of
the other same-letter locks is fixed, and has polynomial height
since the number of slots is fixed. The outside divisor
allocation bound permits at most
\(n^{c-\xi'+o(1)}\) choices of outside label for the offending
slot, with the fixed-truncation caps made arbitrarily close to
their endpoint values. The label determines its complete value.
The empty slot is then determined by the additive relation. This
again yields a fixed deficit below \((m-1)c\). Prime coloring
already ensures coprimality between different letters.

On each remaining tuple, delete from every outside lock all its
primes shared with any other lock. The total numerical weight
lost is at most \(Sk\xi'+o(1)\), where \(k\) is the fixed number
of outside slots and \(S\) is the band's upper degree bound.
At each prime shared among the atomic locks retain the largest
power, breaking ties deterministically. The atomic weight lost
is at most the sum of the pairwise logarithmic gcds, hence at
most \(\binom{m}{2}\xi+o(1)\). At a single prime, this follows
by summing the minima with a slot attaining the largest
valuation. Take \(\xi,\xi'\) so small that the product of the
remaining pairwise coprime moduli has numerical weight strictly
above one. It therefore exceeds twice the absolute height of
all signed entries.

The total original record cost is \(P+z+o(1)<c\). Stripping is
deterministic once the original record tuple is fixed, so the
fixed-record bound (MS17), summed by Hölder over the original
dictionaries, gives the upper exponent
\[
 \frac{m(m-2)}{m-1}c+\frac{P+z}{m-1}+o(1)<(m-1)c .
\]
This contradicts the mixed moment lower bound after the sparse
exceptional tuples are removed. It proves (RA1).

For the information version, cover a strict failure event in
the finitely many nonnegative information coordinates with
positive multiplicity by finitely many boxes with small fixed
mesh and upper corners \(\alpha_i\). Each such box determines
an edge subfamily and restricts type \(i\)'s values to at most
\(n^{\alpha_i}\) possibilities. The strict failure margin persists
at the upper corners. Applying the dictionary version to a
near-full box would contradict (RA1), yielding the pointwise
sparsity assertion. Nonnegativity gives the entropy statement.

### Information of atomic subfactors

If a prescribed atom subfactor has normalized radical mass
\(\beta\in[0,c]\), its numerical weight is \(3\beta\). Repeating
it \(j\ge1\) times gives
\[
 \alpha\ge
 \max\left\{0,\sup_{j\ge1}
       \min\left(\frac{c}{j},
          \frac{3\beta}{4}-\frac{1}{12j}\right)\right\}.
 \tag{RA2}
\]
In particular every divisor retaining half of the degree-three
atom's root mass, \(\beta=1/12\), needs dictionary exponent
at least \(1/24\), by taking \(j=4\). This is stronger than the
earlier tests for a single distinguished slot per letter and for
the full atom alone. For example, as formal assignments to these
listed types, give the three full atoms information \(1/12\)
each and the half atom on \(A\) information \(1/30\).
The one-distinguished-per-letter test (MXP4) on the half atom
and either other full atom asks only for a combined \(5/48\),
whereas these information values sum to \(7/60>5/48\).
The half atom alone needs zero by that test, two full atoms
require \(1/6\) exactly, and adding a third listed type on the
remaining letter keeps the sum at least \(c\). All individual
full atom lower bounds also hold. Yet four copies of the half
atom have total weight one and total cost \(2/15<c\), violating
(RA1).

One may also take in (JR3) a positive-weight atom subfactor
\(K_A\) on \(A\) and the entire atom factors \(F_B,F_C\).
Their weight sum is strictly greater than one. Since
\(H(F_B,F_C)\le H(F_B)+H(F_C)\), we obtain
\[
 \frac{H(K_A\mid F_B,F_C)}{\log n}\ \ge\
 \frac12\left(c+\frac{H(K_A)-H(F_B)-H(F_C)}{\log n}\right)-o(1).
 \tag{RA3}
\]
If the two full atom entropies approach \(c/2\) each, then the
conditional entropy of any positive-weight subfactor on \(A\)
is at least half its own entropy in normalized units.
For a half-mass atom subfactor, (RA2) gives a lower bound \(1/48\).

The formal assignment of dictionary exponents
\(\alpha_i=\beta_i/2\) to every atomic subfactor of radical mass
\(\beta_i\) satisfies every instance of (RA1). Indeed \(W=6P\), and
\(W+4(c-P)=4c+2P\le1\) when \(P\le c\).

### Same-letter additive quadruples

Let \(\mathcal A\) be the complete part on any of the three
letters, and let \(K(x)\) be its full degree-three atom factor.
Take ordered solutions
\[
    x_1+x_2=x_3+x_4,\qquad x_i\in\mathcal A.
\]
The mixed additive Fourier lower bound with four slots all on
this part gives at least \(n^{3c-o(1)}\) quadruples. Discard all
quadruples where one of the six pairs has full atom gcd of
numerical weight at least a fixed \(\xi\in(0,1/12)\).
Each exceptional pair family has exponent strictly below \(2c\)
by the common-factor obstruction. With the pair fixed, at most
\(|\mathcal A|\) choices of one other value determine the last.
Thus the discarded quadruples have exponent strictly below
\(3c\), and at least \(n^{3c-o(1)}\) remain. The retained family
is invariant under exchanging the two positions of either side
and under exchanging sides. These permutations act transitively
on the four positions.

On a retained quadruple put \(K_i=K(x_i)\). Deterministically
retain only the largest valuation at each prime across the
four full atom factors. As in the proof of (RA1), the sum of the
six logarithmic gcds bounds the removed numerical weight.
The remaining pairwise coprime divisors have product of weight
at least \(2-6\xi-o(1)>1\). Fix a tuple
\(\boldsymbol K=(K_1,K_2,K_3,K_4)\) of original records.
The stripping depends only on these records. Thus the conditional
entropy proof of (MS17) applies in this fibre, to give
\[
  H((x_i)_i\mid\boldsymbol K)
       \le \frac23\sum_{i=1}^4 H(x_i\mid\boldsymbol K)
       \le \frac23\sum_{i=1}^4\bigl(H(x_i)-H(K_i)\bigr).
\]
Here all probabilities are for the uniform distribution on the
retained quadruples. The last inequality uses that \(K_i\) is a
function of \(x_i\). All four record marginals are equal by the
permutation symmetry. Write \(a=H(K_i)/\log n\) and
\(J=H(\boldsymbol K)/\log n\). The complete part cap gives
\(H(x_i)\le(c+o(1))\log n\), while the uniform quadruple entropy
is at least \((3c-o(1))\log n\). Since \(\boldsymbol K\) is a
function of the quadruple, adding \(H(\boldsymbol K)\) to the
conditional estimate and rearranging gives
\[
       3J-8a\ge c-o(1).                               \tag{RA4}
\]
In particular
\[
 \frac{H(K_4\mid K_1,K_2,K_3)}{\log n}
       \ge J-3a\ge \frac{c-a}{3}-o(1).                \tag{RA5}
\]
The same lower bound holds for any one record conditioned on the other three. If the dictionary exponent is at most
\(c/2+o(1)\), each such conditional entropy is at least
\(c/6-o(1)=1/36-o(1)\) in normalized units.

Consequently an atom dictionary of exponent strictly below
\(c\) cannot be given an injective identification with elements
of an abelian group that respects all same-letter additive
quadruples. Such a Freiman 2-homomorphism would give
\(K_1+K_2=K_3+K_4\) in that group and let any three records
determine the fourth, contradicting (RA5). More generally the
conditional entropy conclusion holds on any permutation-
symmetric quadruple subfamily with exponent approaching
\(3c\), since the same pair exclusions and fibre bound apply.

## A local-density sector for the atomic common-factor obstruction

Assume \(d=1/3,\Gamma=1/2\), put \(c=1/6\), and suppose
\[
   \nu_W=c\delta_3+\sigma_W\qquad(W=A,B,C),
\]
where each \(\sigma_W\) is nonatomic of mass \(c\). All
three numerical heights are one. Fix a target letter \(T\)
and let \(U,V\) be its opposing letters. Put
\[
  s_T^{\rm opp}
   =\max\{\sup\operatorname{supp}\sigma_U,
                       \sup\operatorname{supp}\sigma_V\},
  \qquad B_T=(1-cs_T^{\rm opp})/3 .
\]
If \(B_T\le0\), (AC1) already excludes a positive common
factor from \(T\)'s atom window. Otherwise suppose that
on each side of three in a fixed neighborhood the measure
\(\sigma_W\), for \(W=U,V\), dominates Lebesgue measure
times \(k_W>0\). Suppose also that this neighborhood on
each side has \(\sigma_W\)-mass strictly greater than
\(\tfrac12(\tfrac12+\epsilon)B_T\) for some \(\epsilon>0\).
Set \(k_{\rm harm}=2/(k_U^{-1}+k_V^{-1})\).
A sufficient condition for the common-factor obstruction
on \(T\) is
\[
 (1+2\epsilon)^2
      \left(c+(\tfrac12+\epsilon)B_T\right)<3k_{\rm harm}.
 \tag{LC1}
\]

### Counting with the common divisor

A common atom-window factor \(F\) of numerical weight
\(W_F>0\) on a near-full subgraph has root mass
\(b=W_F/3>0\). Apply (AC1) to the \(T\)-part in its fixed
residue zero modulo \(F\), using the best supplying
opposing letter and a cutoff increasing to
\(s_T^{\rm opp}\). This gives
\(c\le(1-W_F)/s_T^{\rm opp}\), so \(b\le B_T\).

Put \(t=(1/2+\epsilon)b\). In both free letters retain the
atom and the diffuse mass \(t/2\) nearest three on each
side. The nonatomic condition makes this an interval
selection with continuous endpoints. For a measure on
one side with density lower bound \(k\), the first
distance moment of the nearest mass \(m\) is at most
\(m^2/(2k)\): its distance quantile at mass \(u\le m\)
is at most \(u/k\), and one integrates in \(u\). By the
assumed supply every selected interval is within the
density neighborhood.

Freeze the rest of both free letters, at root cost
\(2(c-t)+o(1)\). Group the retained roots at degree
\(q=\lfloor3/L\rfloor\) in their signed ratio. Its
narrow height is at most \(h=c+t+o(1)\). The positive
error for that ratio is at most the positive degree
deviation from one letter plus the negative degree
deviation from the other. The negative error obeys
the reverse bound. Both wider heights are at most
\[
 E+o(1),\qquad
 E=\frac{t^2}{8}\left(\frac1{k_U}+\frac1{k_V}\right)
   =\frac{t^2}{4k_{\rm harm}}.
\]
No reflection symmetry between the diffuse profiles is
used, since equal masses were selected separately in
each one-sided neighborhood. (LC1) gives the strict
inequality \(E<3b^2/(16h)\).

The one-root Taylor and Newton proof for the broad
symmetric common-factor obstruction now applies without
change. In detail, \(M+1=\lceil4h/b\rceil\) gives
\[
 E<\frac 3M\left(\frac b2-\frac h{M+1}\right),
\]
since \(h/(M+1)\le b/4\) and \(M<4h/b\). Take rank
fraction close to one and horizontal degree just below
the Taylor precision bound. The shallow modulus record
has cost \(b\), the wider grouped ratio is a unit, and
the bounded-degree Newton argument recovers all retained
outside labels of both free letters up to arbitrarily
small cost. The fixed chart and freeze cost is therefore
\[
 2(c-t)+b=2c-2\epsilon b<2c.
\]
All extractions can be taken below this positive margin
after choosing a positive limiting \(b\). A quantitative
lower bound on \(b\) gives a positive uniform margin
for the corresponding same-letter atom-gcd sparsity
argument.

If (LC1) holds for all three choices of \(T\), the proof
of repeated atomic records applies. With
\(s_*=\max_W\sup\operatorname{supp}\sigma_W\), its
dictionary or pointwise information inequality is
\[
 \sum_i m_i\alpha_i
 \ \ge\
  c-\frac{(1-\sum_i m_i w_i)_+}{s_*},
 \tag{LC2}
\]
for any fixed nonnegative integral multiplicities,
with nonnegativity understood if the displayed bound
is negative. In particular full atoms still need
individual information at least \(c/2\), repeating
twice, and every half-root-mass subfactor needs
\(c/4\), repeating four times.

### Uniform neighborhoods of unequal widths

For arbitrary \(\ell_W\in[1,11/10]\), take
\[
 \sigma_W=c\,\operatorname{Unif}[3-\ell_W,3+\ell_W].
\]
Then each measure has mean contribution \(3c\) and
the complete numerical height is one. Each \(k_W\)
can be taken at least \(5/66\), each \(B_T\le1/9\),
and \(\epsilon=1/200\) gives
\[
 (1+2\epsilon)^2
 \left(c+(\tfrac12+\epsilon)B_T\right)
 \le \frac{4090601}{18000000}<\frac5{22}
 \le3k_{\rm harm}.
\]
There is sufficient one-sided mass inside any common
neighborhood with radius, for example, \(1/2\):
its mass is at least \(5/132>101/3600\).

For the symmetric half-width \(\ell=9/8\) one has
\(B_T=5/48\) and \(3k_{\rm harm}=2/9\). Choose
\(\epsilon=1/400\). The left side of (LC1) is
\[
 \left(\frac{201}{200}\right)^2\frac{4205}{19200}
   =\frac{169886205}{768000000}<\frac29 .
\]
The same radius \(1/2\) supplies the required one-sided
mass.

## Anchor sets for arbitrary nonatomic masses

Let \(d>0\) be a normalized limiting box-count exponent, with
nonatomic masses \(c_A,c_B,c_C\) and \(\Gamma=c_A+c_B+c_C\).
Write \(n=X^L\). Assume the indicated modulus letter has
a positive-weight full block of arbitrarily small normalized
root cost, with fullness at least \(\kappa/L\) for a fixed
\(\kappa>0\) after the block is chosen. A positive nonatomic
mass suffices. The support of its restriction to \((0,\infty)\)
contains a compact positive-degree band of positive mass.
After excluding finitely many atomic windows, all individual
outside coordinates have arbitrarily small mass. Greedily
choose a block of small positive mass from that band.
Its numerical weight is positive, and its minimum power is
at least a fixed positive constant times \(1/L\).

If \(0<\Gamma<2d\), all three nonatomic masses are positive.
Indeed, if \(c_A=0\) and \(c_B>0\), (GE1) with a cheap
positive modulus on \(B\) gives \(c_C\ge d\). The same
argument using \(C\) then gives \(c_B\ge d\), contradicting
\(\Gamma<2d\). This is also the relevant check for
finite-scale letters, whose nonatomic mass is zero. When
\(\Gamma\ge2d\), the statements requiring three cheap
moduli below assume them explicitly. A statement about
just one projection only needs a cheap modulus on its
complementary letter, together with the pair-fibre bound
for the projections being used.

Fix enough positive atomic windows on each letter, as
well as a small window at zero. The outside labels record
all the remaining original factors, with part caps
\(M_T\le n^{c_T+o(1)}\). Here the normalized tolerance
can be made arbitrarily small in the sequential limits.
The finite selected atomic windows give boundedly many
grouped roots of height \(n^{O(1)}\). The grouping and
zero-window errors have ordinary numerical height
\(X^{o(1)}\). One first fixes the finite atom list and all
amplification orders, then takes sufficiently small
window widths and box and height tolerances. This is the
same order of limits as in (GE1) and (TR2)–(TR3).

For a pair of fixed outside labels, recording a cheap
positive full complementary modulus and applying the
amplified grouped ratio bound gives a pair fibre of
arbitrarily small exponent. The primitive ratio
determines the triple. If all three pair projections
are needed, greedily thin to a linear hypergraph of
edge exponent arbitrarily close to \(d\), since deleting
the collisions for any chosen edge costs at most three
times the maximum pair fibre. If only one projection
is needed, thin that projection alone. Freeze all cheap
moduli used below at total root cost as small as desired.
This preserves pair uniqueness. Color primes once by
complete letter and grouped slot, at fixed-truncation
subpower cost, exactly as in (TR2). We now use \(E\) for
the remaining edge count. The multiplicative separation
arguments below apply with these actual finite \(E\)
and part sizes \(M_T\).

### A single anchor in one projection

Consider the \(U,V\) outside projection, whose target is
\(T\). Write \(i,j\) for its row and column labels and
\(R_{ij}=U_{ij}/V_{ij}\) for the ratio of the positive
complete magnitudes. The signed additive relation makes
all \(R_{ij}\) congruent to the same sign modulo the
fixed full block on \(T\). Thus for two rows \(i,i_0\)
and any shared column \(j\), the quotient
\(R_{ij}/R_{i_0j}\) is \(1\) modulo that block. The
outside product of \(V\) cancels, and the remaining
outside ratio is the fixed ratio of the \(U\)-labels
\(i,i_0\). The remaining grouped numerator and
denominator slots are exactly the two-edge wedge in
(TR2)–(TR3). Their signs are handled by the sign coloring
in the amplification there. In particular the exact
complete quotients have at most \(n^{o(1)}\) values
uniformly for a fixed pair \(i,i_0\). The uniformity
includes the congruence target, which is a unit and can
depend on the two outside labels. If no narrow slots
are present, rational reconstruction suffices as in
(TR2).

For a fixed anchor row \(i_0\), retain all edges
\((i,j)\) such that \((i_0,j)\) is an edge of the
unthinned-for-this-step graph. No distinct-row
restriction is imposed. Summing over the anchors gives
\(\sum_j \deg(j)^2\ge E^2/M_V\). Thus some anchor
retains at least \(E^2/(M_U M_V)\) edges. On each
retained row \(i\), keep one most frequent exact
quotient \(R_{ij}/R_{i_0j}\). The total loss is
at most the uniform \(n^{o(1)}\) factor, regardless
of how many rows there are. The resulting exact identity
is
\[
    R_{ij}=r_i R_{i_0j}.
\]
At each prime of color \(U\), the valuation of this
ratio is the valuation of \(U_{ij}\), and at each
prime of color \(V\) its negative is the valuation of
\(V_{ij}\). Therefore both complete magnitudes
separate into a rational row factor times a rational
column factor. The factors have polynomial height
since they can be taken from ratios of actual anchor
terms by valuation projections. The retained normalized
exponent is at least
\[
   2d-c_U-c_V.
   \tag{OA1}
\]
The expression is useful when positive. No four-cycle,
non-backtracking path, or minimum degree condition is
required.

### Freeze the joint vector record

The factor-freezing step can be sharpened independently
of the choice of anchors. Consider a graph of exponent
\(D\) with part caps \(a_i\), and put \(G=\sum_i a_i\).
Suppose the vector factor product relation is available.
Use the regularized dictionary exponents \(v_i\) as above.
In addition to the earlier inequalities, the TOTAL
number of matched vector factor triples has exponent
at most \((v_A+v_B+v_C)/2\). To see this, use the
same numerator–denominator sign coloring on each
dictionary and each vector coordinate. A same-part
product has divisor-many representations, so each
same-part multiplicative energy is at most the square
of its dictionary size times \(X^{o_X(1)}\). The mixed
product energy for dictionaries \(i,j\) is at most
the geometric mean of their two same-part energies
by Cauchy–Schwarz applied to their quotient energies.
A further Cauchy–Schwarz inequality over the output
dictionary \(k\) bounds the number of matched pairs
by
\[
  \bigl(|S_i||S_j||S_k|\bigr)^{1/2}X^{o_X(1)}.
\]
Each matched pair determines the third vector, so
this is also the bound for joint vector records.
All coloring and dyadic selections cost only \(X^{o_X(1)}\)
for fixed truncation. A triple of factors can be
frozen as a single joint record.

Sum the three previously established inequalities
\(v_i+v_j-v_k\le 2(a_i+a_j-D)\). This gives
\(\sum_i v_i\le4G-6D\), so the joint record freeze
costs at most \(2G-3D\). It yields an exponent at least
\[
   4D-2G.                                           \tag{OA5}
\]
This improves (OA3). Using \(D=4d-2\Gamma\) from
(OA2) gives \(16d-10\Gamma\), or
\(1/3-16\varepsilon-10\gamma\) in the near-minimum
notation. The original part caps \(c_T\) are retained.
For pure abstract product dictionaries the joint
record estimate is sharp: take three disjoint lists
of \(K\) primes \(p_i,q_j,r_k\) and records
\((p_i/q_j,q_j/r_k,r_k/p_i)\) as \(i,j,k\) vary.
The three dictionary sizes are \(K^2\) and the joint
dictionary has size \(K^3\).

### Refining the labels by small anchor sets

There is a better choice when part caps may increase.
Keep exponent \(d\) before separating, rather than
selecting one anchor per projection. In a projection
with row part \(U\) and column part \(V\), discard
columns of degree less than \(E/(10M_V)\); this
removes at most \(E/10\) edges. Pick a set of
\[
 K_{UV}\ll (M_U M_V/E)\log X+1
\]
anchor rows that meet every remaining column.
For example take \(K\) independent uniform row
samples. A column of degree at least \(E/(10M_V)\)
fails with probability at most
\(\exp(-KE/(10M_UM_V))\). Since the part sizes are
polynomial in \(X\), choose the constant multiplying
\(\log X\) to union bound over the columns. Repetition
of sampled rows does not affect the upper bound.
Give each good column \(j\) one of its sampled
neighbors \(i_0(j)\), called its assigned anchor.
The reference graph here is the graph just after
the common freezes and coloring. The choice of
reference edges need not survive subsequent
selections.

On an edge refine the row label \(i\) to
\((i,i_0(j))\). This enlarges its cap by only
\(\Delta_{UV}=c_U+c_V-d\) in the limiting exponent.
For each refined row the outside ratio in the
two-edge comparison is fixed, so choose a most
frequent quotient as before. The identity
\(R_{ij}=r_{(i,i_0(j))}R_{i_0(j),j}\) is an exact
separation into a function of the refined row and
a function of the original column. It continues
to hold if that column is itself refined.

Cyclically orient the three projections as above,
make the three anchor choices on the reference
graph, and discard the light columns in all three
projections. The discarded union has size at most
\(3E/10\). Each part acts as the row part of just
one of the three projections, so it gains just one
anchor index. Pair uniqueness persists after these
refinements, because an extended label determines
the original label. Select most frequent quotient
values on all refined rows successively in the
three projections, losing only \(n^{o(1)}\)
overall. Thus all three exact separations hold at
edge exponent \(d\) with enlarged total part cap
\[
    G=\Gamma+\sum_{\{U,V\}}\Delta_{UV}=3\Gamma-3d .
\]
Every complete value is now a function of
its own refined outside label on the selected subgraph by the joint
record freeze (OA5), with final exponent
\[
   4d-2G=10d-6\Gamma
       =\frac13-10\varepsilon-6\gamma .             \tag{OA6}
\]
In the displayed positive range all three cheap
moduli are automatic by \(0<\Gamma<2d\). At the
endpoint the conclusion applies through
\(\Gamma<5/9\). It gives no realization of positive
exponent when \(\Gamma\ge5/9\), and it does not
assert the original unsplit outside label cap.


For example three-letter divisor records of total
cost \(\mu\) and weight above the shortest numerical height have (MS15)
upper exponent \((G+\mu)/2\) on the refined
realization. Comparing with (OA6) gives
\(23d\le15\Gamma+\mu\), or
\[
  \varepsilon\ \ge\
   \frac{1/6-\mu-15\gamma}{23}.                     \tag{OA7}
\]
In (MS14), CRT reconstruction bounds the difference of two
possible values of the first term for a fixed third term and fixed
records. Use the shortest letter as the first term, or use any
smaller containing interval, to get the stated weight threshold. All comparisons
extend to a subfamily on which the lock exists
uniformly after choosing its dyadic size lists,
by beginning the reduction on that subfamily.


### Effective outside entropy in place of the spectral cap

On the uniform edge distribution of a current family,
let \(h_T=H(\text{outside label on }T)/\log n\) and
let \(h=\sum_T h_T\) denote its limiting total outside
entropy, taking a subsequence if needed. For any fixed
\(\tau>0\) put
\(I_T=-\log\Pr(\text{actual outside label on }T)/\log n\).
The nonnegative sum \(\sum_T I_T\) has mean \(h+o(1)\).
Markov's inequality keeps a positive fraction on
\(\sum_T I_T\le h+\tau+o(1)\), uniformly at fixed
\(\tau\). Partition this region with a fixed mesh
\(\tau\) in the three coordinates. One box retains a
positive fraction and has nonnegative upper
corners \(b_T\) with
\(\sum_T b_T\le h+4\tau+o(1)\).
Each marginal then takes at most \(n^{b_T}\)
values, by its defining probability lower threshold.
No information independence is being assumed.
Further thinning maintains these caps.

Use \(b_T\) as the actual upper cap of each original
part everywhere in the small-anchor proof. The
original exponent is still \(d\), the pair cap
obeys \(b_U+b_V\ge d-o(1)\) by the pair-fibre
bound, and the enlarged total cap is at most
\(3\sum_T b_T-3d+o(1)\). Consequently taking
\(\tau\) to zero gives the effective version of
(OA6), with exponent at least
\[
           10d-6h.                                  \tag{OA8}
\]
Likewise (OA7) becomes \(23d\le15h+\mu\) whenever
the three-letter lock comparison is applicable.
These versions apply to formal \(\Gamma\) anywhere
in its range. A positive exponent in (OA8) requires \(h<5d/3\).

If \(\Gamma>0\), the cheap modulus prerequisite
is automatic in the positive (OA8) range, even
when \(0<\Gamma<2d\) is not known. In fact if
\(c_A=0\), the argument from (GE1) gives two
other positive masses \(c_B,c_C\ge d\).
The pair-fibre bounds for \((A,B)\) and \((A,C)\)
are then available, because their complementary
letters have cheap moduli. Each such pair has
joint entropy \(d-o(1)\) on the original uniform
edge family, since any fibre has
\(n^{o(1)}\) edges at arbitrarily small normalized
loss. Also \(H(A)/\log n=o(1)\), by \(c_A=0\)
and the outside cap. Consequently each of
\(H(B)/\log n\) and \(H(C)/\log n\) is at least
\(d-o(1)\), giving \(h\ge2d\) in the limiting
statement. Thus a zero-mass modulus exception
cannot occur in the positive effective range.


### A matching obstruction for anchor selection

For \(m\) labels in each part, the edges \((i,i,i)\),
\(1\le i\le m\), satisfy exact pair uniqueness and all
projections are matchings. Every single anchor retains
just one edge. A cover by anchors needs \(m\) rows,
regardless of the orientation. The labels can even carry
primitive positive additive values with all three height
exponents one, by setting \(A_i=2H+2i+1\),
\(B_i=2H+2i+2\), \(C_i=4H+4i+3\), where
\(H\asymp X\) and \(m=o(X)\). The two addends are
consecutive. This model tests the incidence argument only.
The numerical values have no imposed small-radical budget.

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

### A primitive matching with saturated internal root entropy

There is a family of positive primitive integer triples \(A+B=C\)
which is a matching in the three complete-value sets and whose
two varying original root coordinates on each letter have
asymptotically full joint entropy and vanishing normalized
mutual information. The coordinate powers are distinct,
and all six roots in each triple are squarefree and pairwise
coprime. Thus complete matching, exact additivity and
internal coordinate independence are compatible even with
this arithmetic form of the root vectors.

The construction uses powers one and two and does not have a
small radical. The small-radical hypothesis is the additional condition
needed to rule out this mechanism.

**Construction and primitive density**

For positive integers \(x,y\in[N,2N]\), put
\[
 A=x^2(x+4y),\qquad
 B=y^2(5x+2y),\qquad
 C=(x+y)^2(x+2y).
\]
Expansion gives
\[
 C=x^3+4x^2y+5xy^2+2y^3=A+B.
\]
The six root factors are the linear forms
\[
 x,\quad x+4y,\quad y,\quad 5x+2y,\quad
 x+y,\quad x+2y .
\]
Each is of size comparable to \(N\). The coefficient
determinants of any two distinct forms are nonzero and have
prime factors only among \(2,3,5\). This follows directly
from the displayed coefficient vectors: apart from \((5,2)\)
and \((0,1)\), the vectors are \((1,s)\) with
\(s\in\{0,4,1,2\}\). Determinants with \((5,2)\) have
absolute values \(2,18,3,8\), and its determinant with
\((0,1)\) has absolute value five.

Require \(x\equiv y\equiv1\bmod900\). At \(p=2\) the
only even root is \(x+y\equiv2\bmod4\). At \(p=3\) the
only divisible root is \(x+2y\equiv3\bmod9\). At \(p=5\)
the only divisible root is \(x+4y\equiv5\bmod25\).
All three local conditions give exact valuation one.
For \(p\ge7\), two root factors are simultaneously
divisible by \(p\) only when \(x\equiv y\equiv0\bmod p\).
That bad event has density \(1/p^2\).
For each of the six linear forms, divisibility by
\(p^2\) also has density \(1/p^2\), because its
coefficients are coprime. The union of bad residues
therefore has density at most \(7/p^2<1\) at \(p\).

For any fixed cutoff \(P\ge7\), CRT and counting
points in residue classes give a limiting surviving
proportion at least
\[
    \frac{1}{900^2}\prod_{7\le p\le P}
                       (1-7/p^2),
\]
which is bounded below by a positive constant
independent of \(P\). To discard the primes \(p>P\),
the number of parameter pairs with some \(p>P\)
dividing both \(x,y\) is at most
\[
 O\!\left(\sum_{P<p\le 2N}(N/p+1)^2\right)
       =O(N^2/P+N\log N+N).
\]
For a square factor \(p^2\mid L(x,y)\), with
\(L\) one of the six forms, one has \(p\ll\sqrt N\).
At least one coefficient of \(L\) is a unit modulo
\(p^2\), so fixing the other parameter gives
\(O(N/p^2+1)\) choices. Summing over \(p>P\) gives
\(O(N^2/P+N^{3/2})\) for each form. Choose \(P\) so
large that the \(N^2/P\) error is smaller than half
the surviving lower proportion, then take \(N\) large.
This leaves \(\gg N^2\) pairs for which all six
roots are squarefree and pairwise coprime. The
resulting triples are primitive.

**Matching and coordinate entropy**

For each of the three letters the map from \((x,y)\)
to its two root values is an invertible linear map
over \(\mathbb Q\), with fixed nonsingular matrix.
For a fixed complete value \(T\), the root raised to
power two has support on \(\operatorname{rad}(T)\) because
its square divides \(T\). The fixed-radical integer bound
gives \(N^{o(1)}\) possibilities uniformly since \(T\ll N^3\).
The other root is then determined, and the root pair determines
\((x,y)\). Thus each complete vertex is incident to at most
\(N^{o(1)}\) selected pairs.
Greedy selection deleting conflicts on any of the
three complete values leaves a matching of size
\(N^{2-o(1)}\). In particular all three complete
marginal entropies approach \(2\log N\) in the
uniform edge distribution.

The two root coordinates on each letter jointly
determine the parameter pair, so their joint
entropy is \((2-o(1))\log N\) on the matching.
Each one has entropy at most \(\log N+O(1)\)
because it is an integer of size \(O(N)\).
By entropy subadditivity each one therefore has
entropy \((1-o(1))\log N\), and their mutual
information is \(o(\log N)\). This gives full
normalized coordinate entropy in all three
letters simultaneously.

Finally, on the retained pairs the radical of
\(ABC\) is the product of the six squarefree
coprime roots. It is comparable to \(N^6\),
whereas \(C\) is comparable to \(N^3\). Thus its radical has height exponent two. The small-budget
regime is still open for matchings of this type.

### Complete-letter curve covers of the cubic matching

The cubic family above has exact additive and multiplicative
records and still requires \(N^{1-o(1)}\) bounded-degree curves
to cover any of its complete-letter projections.
Indeed, the Jacobian of the complete map \((x,y)\to(A,B)\)
at \((1,1)\) is
\[
 \det\begin{pmatrix}11&4\\5&16\end{pmatrix}=156.
\]
Thus the map is dominant over characteristic zero.
The map to either other complete pair differs from it
by an invertible linear substitution. Any nonzero
polynomial in such a pair of total degree at most \(D\)
therefore pulls back to a nonzero polynomial of total
degree at most \(3D\). It has at most \(O(DN)\) zeros in
the parameter box by the grid polynomial zero bound.
In detail, write it in \(y\), exclude the at most \(3D\)
values of \(x\) at which its leading coefficient vanishes,
and then bound the number of \(y\)'s at a remaining \(x\)
by the degree in \(y\). This bound applies to arbitrary
coefficients of the proposed complete-letter curve.
It gives the claimed cover-size lower bound at fixed \(D\).
It even gives \(M/O(DN)\) for a budget-dependent \(D\).

There is a corresponding rank formulation. Let
\(r=\binom{D+2}{2}\), and sample \(r\) matching edges
independently and uniformly. The evaluation matrix of
the degree-at-most-\(D\) complete-pair monomials has full
rank except with probability \(O_D(N/M)=N^{-1+o(1)}\).
After each of the first \(r-1\) evaluations, choose
a nonzero polynomial in the current evaluation kernel.
Failure of the next evaluation to increase rank implies
that this polynomial vanishes at the next point. Its
conditional probability is at most \(O(DN)/M\).
Sum these estimates through the \(r\)-th evaluation.

The same cubic family lies in a single rational projective
image given by its three homogeneous cubic forms. This
distinguishes a projective parametrization from an affine
complete-letter curve cover. At small budgets the diffuse
projective case is covered by a different residue argument.
[Projective parametrizations](ABC_RATIONAL_PROJECTIVE_CARRIERS.md)


## Additive energy of wholly diffuse complete-value dictionaries

Let \(P,Q\) be two complete-value dictionaries extracted
from anatomic boxes in the sequential small-budget limits.
They may be the very same dictionary. Write \(n=X^L\).
Assume each letter has normalized root cap tending to the
same number \(l>0\), each weak scaled-degree measure
\(\nu_P,\nu_Q\) has mass \(l\) and no atoms, and their
full numerical heights tend to \(h_P,h_Q\), respectively.
Set
\[
  e_W=h_W-\int r\,d\nu_W(r)\quad(W=P,Q).
\]
Thus \(e_W\) is the escaped first moment. Fix a choice of sign
\(\sigma\in\{-1,1\}\). There is a positive normalized deficit
\(\delta\), depending on the limiting two-letter profile,
such that, uniformly for nonzero integers \(D\) allowed
to vary with \(X\),
\[
   \#\{(u,v)\in P\times Q:\ u+\sigma v=D\}
           \ \le\ n^{l-\delta+o(1)} .               \tag{DS1}
\]
No coprimality condition on \(u,v\) is required. One can
take a fixed positive \(\delta\) throughout any family
with a uniform strict quantile margin, bounded retained
degrees, and uniform small-coordinate errors as in the
proof below.

### Selecting enough numerical weight

Let \(\rho=\nu_P+\nu_Q\), of total mass \(2l\), and
write \(F_\rho(t)\) for its fractional upper moment
at mass \(t\). Since it is nonatomic and has positive
mass, its upper half has strictly larger mean than its
lower half. To see strictness, select a median \(r_0\)
separating the two halves, which is possible by nonatomicity.
On the upper half \(r>r_0\) almost everywhere and on the
lower half \(r<r_0\) almost everywhere. In particular
\[
 F_\rho(l)>\tfrac12\int r\,d\rho(r).
\]
Choose a letter \(W_0\) with \(e_{W_0}=\max(e_P,e_Q)=e_*\).
It follows that
\[
 F_\rho(l)+e_*
   >\tfrac12(h_P+h_Q)\ge h_{\min}:=\min(h_P,h_Q).
 \tag{DS2}
\]
By continuity of the fractional upper moment, one can use
a mass \(p<l\) with a fixed positive margin in (DS2).
Moreover, one can restrict the selected degrees to a
bounded interval: increase a continuity cutoff and
use integrability to make the lost finite moment as
small as desired. Equivalently one can pick continuity
bands in the two letters, or their prescribed fractions
near a continuity quantile, with combined mass at most
\(p\), numerical weight at least \(w\), and
\[
  w+e_*>h_{\min}.
\]
Choose a compact containment bound \(S\) for all selected
scaled degrees; an initial lower cutoff \(s>0\) can be
imposed with arbitrarily small loss of numerical weight.

These selections are realized by whole original root
coordinates, with arbitrarily small mass and moment
losses. Indeed, each individual coordinate in a
fixed compact scaled band has mass tending to zero.
A coordinate of nonvanishing mass would give a point
atom along a subsequence. Choose continuity intervals
for the retained pools and rank prelimit coordinates by
degree, with a final round down or round up at this
vanishing cost. The full moment on bounded continuity
windows converges.

Choose a continuity cutoff \(H>S\) arbitrarily large and
retain the original full tail on the selected endpoint
\(W_0\). Its normalized root cost is at most
\(\zeta+o(1)\), with \(\zeta\) arbitrarily small,
by the uniform bound on first moments and \(H\) large.
Its full numerical weight tends to
\(e_*+\int_{(H,\infty)}r\,d\nu_{W_0}(r)\ge e_*\).
The bounded selection and the tail on this endpoint
use disjoint degree positions. If \(e_*=0\), the tail
may be omitted. Choose a strict margin so that
\(p+\zeta<l\), and choose \(\xi>0\) satisfying
\[
     w+e_* -2S\xi>h_{\min}.
 \tag{DS3}
\]
All choices are fixed before taking sufficiently small
box and height tolerances.

### Uniform shift count

Use the roots of the bounded-degree selections as the
band records, and the root tuple determining the full
tail power on \(W_0\) as the extra record. At fixed
truncation the original bounded-degree anatomic roots
are squarefree and have disjoint prime supports
within each value. The single packed root of possible
non-squarefree type has scaled degree tending to
infinity and may be in the retained tail. Its possible
values are still bounded by its box root cap.

Any prime common to \(u\) and \(v\) divides \(D\).
Call a pair bad if the product of the selected
unpowered band roots on either endpoint has radical
overlap with \(\operatorname{rad}(D)\) at least \(n^\xi\).
Whenever a pair exists, \(|D|\le2X^{O(1)}\).
For each endpoint, the usual outside divisor allocation
count bounds the number of bad values by
\(n^{l-\xi+o(1)}\), uniformly in \(D\): choose a divisor
of \(\operatorname{rad}(D)\) of the stated mass, allocate
its primes to the selected positions, then count only
root multiples at those positions. The divisor and
assignment counts are \(X^{o_X(1)}\) at fixed truncation.
The assigned divisor saves \(\xi\) in the normalized
root exponent. A single endpoint determines its partner
through the fixed \(D\), so this also bounds bad pairs.

On each remaining pair strip all primes dividing \(D\)
from both bounded-band locks. Keep the one full tail
lock intact. The two remaining full endpoint divisors
are coprime, because a common prime would divide \(D\)
and hence cannot survive on the endpoint without the
tail. The strip loses numerical weight at most
\(2S\xi+o(1)\) by squarefreeness in the bands and their
upper degree bound. The product of the two coprime
divisors has weight strictly exceeding \(h_{\min}\)
by (DS3). Either endpoint is zero modulo its own
divisor and lies in a fixed residue determined by \(D\)
modulo the other divisor. Choose the shorter endpoint.
Its height is less than the product of the moduli
for sufficiently small losses, and so the fixed record
admits at most one pair.

Stripping is deterministic given the records and the fixed
\(D\). The unstripped band/root and tail records have
at most \(n^{p+\zeta+o(1)}\) possibilities. Together the
good and bad branches give exponent at most
\[
    \max(p+\zeta,l-\xi)+o(1)<l .
\]
This proves (DS1). The use of exactly one tail lock is
what avoids any need to control the numerical weight of
overlap at escaped scaled degrees.

### Energy and iterated sums

If \(P=Q=S\) and \(|S|=n^{l+o(1)}\), apply (DS1) with
\(\sigma=-1\) to every nonzero shift. Then
\[
 E_4(S)=\sum_D r_{S-S}(D)^2
 \le |S|^2+n^{l-\delta+o(1)}|S|^2
 \le n^{3l-\delta'+o(1)}
\]
for some \(\delta'>0\). With either sign, bounding the
largest nonzero representation count in this way
gives
\[
             |S\pm S|\ge n^{l+\delta+o(1)} ,
\]
up to decreasing \(\delta\). The difference-zero
contribution is lower order because
\(|S|^2-r_{S-S}(0)=|S|^2-|S|\). The same representation
argument gives a uniform fixed-target upper exponent
\((m-1)l-\delta\) for every fixed \(m\ge2\) sums of
\(m\) elements of \(S\), by fixing all but two and
using the plus-sign pair estimate. Young's inequality
also transfers the squared energy deficit to higher
equal-sum energies. The conclusions apply to saturated marginal dictionaries
with arbitrary edge density, including a matching-shaped family.


Coprime pairs from actual limiting boxes have an additional
uniform version through the [pooled degree inequality](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-pooled-degree-inequality-and-primitive-shift-fibres).
That version retains both escaped moments and can also use
prime subproducts inside large atomic coordinates. Here the
one-tail selection permits arbitrary common divisors.

## Primitive two-divisor fibres attaining the unsaved count

At budgets tending to \(3/2\), there are primitive triples
with balanced radical masses and two fixed full divisors
whose fibre of distinct targets attains the scale
\(n^{1/3-t-o(1)}\). Here \(t\) is the normalized root cost
of that divisor record.
For \(x\in[3N,4N]\), \(y\in[N,2N]\), use
\[
 U=x^3(x+2y),\qquad V=y^3(2x+y),\qquad
 T=(x+y)^3(x-y).
\]
All are positive of size comparable to \(N^4\), and
\(U=V+T\) by expansion. Fix \(2/3<\alpha<1\).
Take primes \(p,q\) in disjoint fixed-multiple intervals
at scale \(N^\alpha\), and consider \(x=p u,\ y=q v\)
in the indicated boxes. The moduli \(p,q\) are
larger than the cofactor ranges \(O(N^{1-\alpha})\),
so \(p\mid x\) and \(q\mid y\) hold to valuation one,
while \(p\nmid y\) and \(q\nmid x\).

Require \(x\equiv1\), \(y\equiv12\bmod72\) and
\(\gcd(u,v)=1\). Uniformly for each such \(p,q\)
there are \(\gg N^{2(1-\alpha)}\) pairs. Indeed the
two fixed congruence classes of \(u,v\bmod72\)
do not allow simultaneous divisibility by two or
three. At each other prime the joint forbidden
density is \(1/\ell^2\) for \(\ell\) that prime.
The elementary finite-prime sieve and the bound
\(O(H^2/P+H\log H+H)\) for common primes above \(P\)
in ranges \(O(H)\), \(H=N^{1-\alpha}\), prove the
claim uniformly. For every prime \(\ell\ge5\),
the six root linear forms
\[
 x,\quad x+2y,\quad y,\quad 2x+y,\quad x+y,\quad x-y
\]
have pairwise independent coefficient vectors over
\(\mathbb F_\ell\). Any common prime of two forms
must then divide both \(x,y\), which is excluded
by the parameter choices. At two, only the letter
\(V\) is divisible, with fixed valuation seven.
At three, only \(V\) is divisible. Thus these are
primitive triples, and \(F=p^3\mid U,\ G=q^3\mid V\)
are full prime-power divisors of the indicated terms.

One may secure matching lower radical sizes as well
as the immediate upper caps. Allow \(p,q\) to vary
over their respective prime intervals. The prime
number theorem gives \(\gg N^\alpha/\log N\) choices
in each interval. Each parameter pair \(x,y\)
occurs at most once in this count, since \(2\alpha>1\)
and either integer cannot contain two different
prime factors from its respective interval.
Thus there are \(\gg N^2/(\log N)^2\) distinct
qualified pairs. For any fixed \(\eta>0\), the
number of integers \(z\ll N\) with
\(\operatorname{rad}(z)\le N^{1-\eta}\) is at most
\(N^{1-\eta+o(1)}\) by the fixed-radical support
bound. Each possible integer value of any of the
six nonzero forms has \(O(N)\) parameter preimages.
Discarding parameters where any of their radicals
is too small therefore costs \(N^{2-\eta+o(1)}\),
a lower order term. Taking \(\eta\) slowly to zero
makes each linear-factor radical \(N^{1-o(1)}\) on the
retained pairs. Distinct root factors are coprime
apart from the fixed small-prime contributions on \(V\),
so all three complete radicals are \(N^{2-o(1)}\).

Pigeonhole \(p,q\) after this discard to obtain
\(N^{2(1-\alpha)-o(1)}\) parameter pairs for one
fixed record. Since \(FG\asymp N^{6\alpha}\gg N^4\),
each target in this fibre determines both free
complete values by CRT. A fixed free value \(U\)
has only \(N^{o(1)}\) parameter preimages, since
\(x^3\mid U\) and fixing \(x\) determines \(y\).
It follows that the fibre has
\(N^{2(1-\alpha)-o(1)}\) distinct targets.

Normalize by \(n=N^6\), so each of the three radical
masses is \(1/3+o(1)\) and the record root cost is
\(t=\alpha/3+o(1)\). The fibre has
\(n^{1/3-t-o(1)}\) distinct targets, meeting the
unsaved scale with the radical cap imposed on the
target itself. Its numerical heights are \(\asymp N^4\),
so this construction has budget tending to \(3/2\).
The two principal degrees per letter occur at scaled
positions \(3/2\) and \(9/2\). The construction shows
that the two full divisors and even the three radical caps
can meet the unsaved scale with low valuation orders.
