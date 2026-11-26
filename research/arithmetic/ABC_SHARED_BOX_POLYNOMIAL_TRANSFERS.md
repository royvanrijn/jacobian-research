# Polynomial transfers and fresh primes on shared-value abc boxes

A rational transfer built from a shared-value box has an exact
obstruction to reusing the primes already present in the box.
The monomial identities on the three-axis difference graph
reduce to original rows. A quantitative resultant argument
then bounds specialization onto small outside support,
including coefficient records of bounded degree that vary
with height.

## Three-term monomial identities on a difference graph

Let \(G\) be any finite simple graph. Assign a variable
\(x_v\) to each vertex and a linear form \(L_e=x_v-x_w\)
to each edge, with either orientation. Work over a field
of characteristic zero. Common translation can be divided
out by setting one variable to zero. This does not add
linear dependencies among the edge forms.

**GP1.** Suppose three nonzero signed scalar multiples of
monomials in the \(L_e\) sum to zero identically. Cancel their
common monomial. Either all terms are constants, or their
edge supports are exactly

- the three single edges of a triangle, one per term, or
- the three pairs of opposite edges of a four-clique, one
  pair per term.

In both nonconstant cases each used edge has exponent one.
Up to individual signs, term ordering and a common scalar,
the identities are the triangle relation and the Plücker
relation, respectively.

To prove it, the three reduced monomials are pairwise
coprime by unique factorization and the identity. They are
homogeneous, so all have the same degree. A degree appearing
only once cannot cancel. Thus if one is constant they are
all constant. Color the edges in the three nonconstant terms
red, blue and green, respectively.

For any red edge \(e\) and blue edge \(f\), the green term
vanishes identically on \(L_e=L_f=0\). Its restriction is a
product in the integral polynomial ring of that linear
subspace, so some green edge form lies in the linear span
of \(L_e,L_f\). If the edges are disjoint, no other edge
form is in that span. If the edges meet, the only
additional possible edge is the one closing their triangle.
The same reasoning applies to any two colors. This is
the graphical description of the rank-two flats: the
only multiple ones are the three forms on a triangle.

Choose a rainbow triangle with \(ab\) red, \(bc\) blue
and \(ac\) green. Any additional red edge has to meet both
\(bc\) and \(ac\). It is therefore \(cd\) for some vertex
\(d\) outside the triangle. The closing triangles require
\(ad\) blue and \(bd\) green. These six edges partition
exactly as the opposite pairs of the four-clique. In
this pattern a further red edge would have to be \(ce\)
with \(e\ne d\), by the same first two tests, but then it
would miss the blue edge \(ad\). The same argument in the
other colors rules out all further edges. This proves
the support assertion.

On the intersection belonging to a rainbow triangle,
the only vanishing edge forms are its three edges.
Use two independent normal coordinates \(u,v\) and work
over the rational function field of the intersection.
The three initial homogeneous terms are nonzero scalar
multiples of \(u^a,v^b,(u\pm v)^c\), with \(a,b,c>0\).
A least degree appearing only once cannot cancel. Two
initial powers of distinct linear forms are linearly
independent, so the least degree cannot appear just twice
either. Hence \(a=b=c\). If this common degree exceeds
one, the third power has a nonzero mixed coefficient,
which the first two cannot cancel. Thus the exponents are
all one. Every edge of the two listed support patterns
occurs in a rainbow triangle. Finally the first two
products in either pattern are linearly independent, and
the familiar triangle or Plücker relation furnishes the
unique coefficient ratio. This proves (GP1).
[Graphic arrangements and their triangle and four-clique multinets, Aprodu–Farkas–Raicu–Suciu, §9.8](https://arxiv.org/pdf/2512.10160v1)

Put
\[
                    \mathcal R_G=\mathbb Q[x_v:v\in V(G)]
                         [L_e^{-1}:e\in E(G)].
\]
After clearing Laurent powers, (GP1) classifies \(f\) with
\(f,1-f\in\mathcal R_G^\times\). For a tripartite graph,
a nonconstant such \(f\) is the ratio of two differences
on one triangle, with the coefficients fixed by its
triangle relation. In a complete shared-value box these
triangles are exactly its original rows.

## Two output marks require a new difference factor

**GP2.** On a \(K_4\)-free graph no two distinct rational
functions \(f,g\), with at least one nonconstant, can have all
five of \(f,1-f,g,1-g,f-g\) units of \(\mathcal R_G\).

For a proof, restrict to the configuration where all finite
potentials are distinct. The functions give two distinct
ordered output marks together with \(0,1,\infty\). By (GP1)
each function is constant or a triangle ratio, so both are
invariant under translating and scaling all input potentials.
Regard this as a map on spherical configurations with one
extra input mark. Send that extra mark to infinity by a
Möbius map before evaluating the finite potentials. The
result is independent of the normalization by affine
invariance, and local holomorphic normalizations make it
a holomorphic map of ordered spherical configuration spaces.
Chen and Salter's Theorem 4.4 says that such a map is, up to
one common projective transformation, constant or a selection
of input marks. Normalizing the first three output marks
to \(0,1,\infty\) eliminates the common transformation.
At least one output is nonconstant, so five input marks
are selected, with three used as normalizers.

Return to the chart where the extra input is infinity.
The extra mark has to be one of the three normalizers:
otherwise at least one output is a cross-ratio of four
finite input marks. Its value and complement use all six
differences on those marks, each in exactly one of their
three relatively prime polynomial terms, contrary to
\(K_4\)-freeness. Permute the three output boundaries so
that infinity maps to infinity. This preserves the unit
tests, since the difference of two transformed values is
a constant times their old difference divided by their
unit denominators. The two outputs now have the form
\((x_i-x_a)/(x_b-x_a)\) and \((x_j-x_a)/(x_b-x_a)\)
on four distinct finite marks. Their values, complements
and difference again require all six edges. This proves
(GP2).
[Chen–Salter, Theorem 4.4](https://link.springer.com/article/10.1007/s00209-025-03893-6)
· [The change-of-markings argument with the full input support](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#universal-arithmetic-support-forces-a-change-of-markings)

## Fresh primes for a polynomial factor on actual boxes

Use the setup of
[the shared-value box count (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem).
There are \(h=3m\) potential vertices, residue alphabet bounds
\(R_1,R_2,R_3\), their product \(R\), and a bound \(K\) for
the complete labels in each pair class. The lengths \(L_e\)
in the potential graph are the actual complete values.
All potentials are determined up to translation by the
tree lengths. Denote their union of prime supports by \(S_{\rm box}\).
For a positive integer \(Z\), write
\(\operatorname{rad}_{\rm out}(Z;S)\) for its radical outside
a set of primes \(S\).

Use translation-invariant polynomials on the potentials,
specified as polynomials in the coordinates obtained by
setting one potential to zero. Fix bounds on their total
degrees and polynomial coefficient heights. A finite
collection of polynomials may vary with height and the
specified packet tuple, but is fixed before selecting the
box. Let \(P\) be a nonconstant primitive irreducible
integer polynomial of this kind that is not a scalar
multiple of any crossed difference \(L_e\). Choose a
vertex \(v\) on which \(P\) depends and denote its axis by
\(i\). Reanchor at a different potential; translation
invariance ensures that \(P\) still depends on the
coordinate of \(v\). Let \(\mathcal H\) be a fixed bounded
collection of nonzero polynomials not divisible by \(P\)
over \(\mathbb Q\), with the same degree and height bounds.
Let \(S_0\) be a prescribed pool of polynomial prime product.

**GP3.** For every \(\varepsilon>0\), the number of ordered
boxes on which \(P\) or an element of \(\mathcal H\) vanishes,
or on which
\[
 \operatorname{rad}_{\rm out}\!\left(
      |P|;\ S_0\cup S_{\rm box}\cup
           \bigcup_{H\in\mathcal H}\{p:p\mid H\}\right)
                    \le Z
\]
is at most
\[
 O\!\left(
       R^m K^{h-2}
       + Z X^\varepsilon\, \frac{R^m K^{h-2}}{R_i}
     \right).                                                \tag{GP3}
\]
Constants are uniform in the chosen polynomials and pool,
with dependence on the fixed parameters and \(\varepsilon\).
The height parameter in the inner allowance is the original
\(X\). The bound does not require near-fullness.

To prove this, omit \(v\), choose the residues of all other
vertices, and choose the \(h-2\) labels on a spanning tree
of their connected tripartite graph. This costs at most
\(R^m K^{h-2}/R_i\) and fixes all the anchor potentials up to
translation. Use the tree coordinate differences, so at
fixed anchor residue lists the candidate labels belong to
a product of \(h-2\) sets of size at most \(K\).
Each coordinate change between those differences and
anchored potentials is invertible and linear. All potentials
in the chosen gauge have magnitude \(O_m(X)\).

Over the rational function field of the anchors, \(P\)
is irreducible as a polynomial in the remaining variable
\(x_v\). This is Gauss's lemma, since \(P\) was irreducible
and depends on \(x_v\). Every incident edge form has nonzero
resultant with \(P\) in \(x_v\): a zero resultant would
make the edge form divide \(P\) before specialization.
The same nonzero-resultant condition holds with every
member of \(\mathcal H\) that depends on \(x_v\).
For the other members use their nonzero anchor values.
The resulting anchor polynomials, together with the
leading coefficient of \(P\), have bounded degree. Each is
nonzero before specialization.

For each fixed anchor residue list, a nonzero anchor
polynomial of degree at most \(D\) has at most
\(D K^{h-3}\) zeros on the tree-label grid. Pad smaller
label sets to size \(K\), then induct on the number of
variables by separating the leading coefficient in the
last variable and using the univariate root bound.
For a bad anchor there are at most \(R_i K\) choices of
the omitted potential, by using any one of its neighbors.
The bad anchor contribution is therefore
\(O(R^m K^{h-2})\). Any vanishing of \(P\) or a member of
\(\mathcal H\) on a good anchor has at most a bounded number
of choices of \(x_v\) and fits in the second term.

For a good anchor all nonzero resultants, fixed forbidden
polynomial values and fixed edge labels have polynomial
height. Put their supports and \(S_0\) into one fixed
pool. At a prime outside it, if the evaluated \(P\) vanishes
then no forbidden factor or any original complete label
vanishes there. For the moving factors this is the
univariate resultant identity, and all other factors
are fixed on the anchors. Hence the outside radical of
the evaluated \(P\) off this fixed pool is at most the
displayed outside radical. Enumerate its squarefree
support at cost at most \(Z\), then enumerate the
polynomial-height signed values of \(P\) supported on
the augmented pool, at uniform height-subpower cost.
Each nonzero value has at most \(\deg P\) preimages in
\(x_v\). This proves the good anchor contribution and
(GP3).
[Uniform supported-value enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Related resultant method on independent primitive directions](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#fresh-radical-costs-of-extra-polynomial-factors)

## Consequences for actual rational transfers

Consider a rational function of the potential differences.
Write it with a primitive integral polynomial presentation
\((F_1,F_2,F_3)\) satisfying \(F_1+F_2=F_3\), with the
three polynomials pairwise coprime over \(\mathbb Q\), and
assume fixed degree and polynomial coefficient-height
bounds. Clear numerical scalars before evaluation.
The integer values, when nonzero, give a positive primitive
additive triple after signs and removal of their numerical
gcd. Suppose the rational function is neither constant nor
one of the ratios coming from an original triangle.
By (GP1), some irreducible nonconstant factor \(P\) of one
of the three polynomials lies outside the edge forms.
Those factors satisfy fixed polynomial height bounds. For
instance a Kronecker substitution with base larger than
the degrees preserves all coefficient data in the
product and its factors, and the fixed-degree univariate
factor height bound applies after substitution.

Use (GP3) with the other two output polynomials in
\(\mathcal H\). At its counted private primes \(P\)
divides the integral value of its output term, and neither
other term is divisible. This follows by taking \(P\)
primitive integral and applying Gauss's lemma to its
polynomial factorization. The primes therefore survive
primitive reduction. Thus (GP3) also bounds the number
of boxes on which the primitive output has outside
radical at most \(Z\) off the original box labels and \(S_0\).

The same assertion applies to a normalized configuration
of at least two generically distinct rational outputs,
with at least one nonconstant. For the support of
such a configuration use the three letters of every
individual normalized output and the primitive
determinant between each pair of outputs. If any
individual output uses a forbidden factor, the
preceding bound applies. Otherwise (GP2) gives a
forbidden numerator factor of a pair difference, whose
denominators are products of allowed edge forms up to
constants. Include these denominators, and if necessary
the input denominators of both chosen output coordinates,
among the forbidden factors. At the counted fresh
primes the two rational coordinates are integral and
coincide modulo the prime. Their normalized primitive
determinant consequently retains the prime. Factor
heights and degree bounds remain uniform after forming
the pair difference.

For a simultaneous packet stratum, let \(n=X^L\),
\(d=s_1+s_2+s_3\), \(t=1/6-d\). At fixed positive
margins and fixed \(m\), (LB4) gives \(R^m K^{h-1}\)
boxes up to arbitrarily small normalized power losses,
where \(R_i=n^{s_i+O(\eta+\kappa)}\) and
\(K=n^{t+O(\eta)}\). A catalogue of at most \(n^a\)
polynomial transfers of the indicated bounded degree
and height may be chosen separately for each packet tuple.
If
\[
           a<t,\qquad a+b<t+\min_i s_i
\]
both hold with fixed positive margins, all but a
power-small fraction of the boxes avoid a cheap
new output of any map in the catalogue. Here cheap
means an outside radical at most \(n^b\) off the full
original label support and the prescribed pool.
All arithmetic exceptions in (GP3), including exact
vanishing of any nonzero bounded-degree polynomial used
in the argument, are included in the bound.
[The forced box exponent (LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)

## Joint private support of several transfer factors

**GP4.** Let \(\mathcal P\) be a prescribed collection of at
most a fixed number of pairwise nonassociated irreducible
forbidden factors as in (GP3). For each \(P\in\mathcal P\),
let \(s(P)\) be the largest price of an axis containing a
vertex on which \(P\) depends. On all but a power-small
fraction of the actual near-full boxes, each evaluated
factor \(P\) contributes radical greater than
\(n^{t+s(P)-\sigma}\), for any fixed \(\sigma>0\), on primes
absent from the original box labels, all other factors in
\(\mathcal P\), and the prescribed pool. In particular
these contributions may be added over the factors. Choose
the pivot vertex attaining \(s(P)\), and use (GP3) with
the other factors as the forbidden collection. Take
a finite union of exceptions. If the factors arise from
normalized rational outputs, apply the argument with
the indicated denominator and other-term tests as
needed. Duplicate factors are counted once. The private
contributions then survive primitive reduction in a role
or an output secant where the factor occurs.

The original-label private-prime bound is also valid
with these additional forbidden factors in the pool.
For an edge \(e\), omit its endpoint on the larger-price
axis and write \(x_u\) for the other endpoint, a known
anchor. Add the support of \(P|_{x_v=x_u}\) for every
\(P\in\mathcal P\) to the fixed anchor pool. These are
nonzero on all but a power-small fraction of boxes by
the same tree grid test, since \(L_e\nmid P\) as
polynomials. At every remaining prime dividing \(L_e\),
evaluation of \(P\) agrees with its indicated anchor
value modulo the prime. Apply (PV1) with this enlarged
pool. Thus jointly the box edge labels carry the
private exponents of (PV2) and the new polynomial
factors carry the exponents in (GP4). These prime
supports are pairwise disjoint across both collections.
[Private primes of actual box edges (PV1–PV2)](ABC_SHARED_BOX_PRIVATE_PRIME_SUPPORT.md#a-bound-for-one-position)

For example, take all the differences between two
potentials in the same axis as \(\mathcal P\). These
linear factors are distinct and none is a cross edge.
Their private prime supports, outside every original
box label and the fixed pool, have combined radical
greater than
\[
                 n^{\binom m2(3t+d-3\sigma)} .
\]
This also excludes the case where those differences,
after their fixed packet factors are removed, have
uniformly much smaller radical. The fixed-label shift
argument separately improves this example by a factor
of two in the outside radical exponent.
[Same-axis shift bound (PA4)](ABC_PRIVATE_RADICALS_ON_ARRAY_DIFFERENCES.md#same-axis-shifts-in-one-actual-array)

For a single simple example, extend a triangle by a
fourth potential in the same axis as one of its vertices.
The Plücker ratio then uses the difference in that axis.
This missing edge is an immediate case of (GP3), so
its generally fresh primes cannot be charged to the
old cross edges. The same statement applies to a
nonlinear residual numerator or denominator factor.
A strict count for the original boxes still needs
an upper bound using the three original radical caps,
or a forced arithmetic transfer with a sufficiently
small outside support on the actual restricted boxes.
