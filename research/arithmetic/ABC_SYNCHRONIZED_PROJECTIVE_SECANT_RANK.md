# Joint power-class rank of two boxes and their primitive secants

Two actual boxes on the same residue lists give two arrays of
complete labels and an array of differences between corresponding
labels. Almost all synchronized
pairs have independent power classes on all three sets. The odds
of their primitive secant directions then have an exact joint
rank with the original data. The unknown normalization gcd
cancels in the odds.

## A finite count on synchronized pairs

Use the [shared-value box setup (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem),
with \(m\ge2\), \(E_0=3m^2\) indexed edges, residue alphabet bounds
\(R_i\), \(R=R_1R_2R_3\), and pair-class bound \(K\). A synchronized
pair is an ordered pair of boxes on the same ordered residue lists.
Their potential vectors are \(x^{(1)},x^{(2)}\), each considered
up to common translation within its own array. Let \(L_e^{(j)}\)
be the two complete positive labels in position \(e\). Define
\[
       D_e=L_e^{(2)}-L_e^{(1)}.
\]
Use the signed edge convention from the original additive
triangles, so the \(D_e\) are signed differences between the
potentials of \(x^{(2)}-x^{(1)}\). Distinct translations of
the two arrays do not change any \(D_e\).

For each array use either its original labels or their individual
least order-\(M\) multipliers \(Z_M(L_e^{(j)})\) as observables.
The third set of observables consists of the \(|D_e|\).
For original labels on both arrays let \(q\ge2\) be arbitrary
and set \(\Delta=1\). If either array uses multipliers, let
\(q\mid M\), \(q\ge2\), and set \(\Delta=X^{1/M}\). On a pair
with all \(D_e\ne0\) write \(z_f\) for the resulting \(3E_0\)
positive observables, including their position and array indices.
Use \({\cal D}_{q,P}\) as in (BX1), with a prescribed pool \(P\)
of prime product at most \(X^B\), \(B\) fixed.

**SQ1.** The number of synchronized pairs with some \(D_e=0\),
or with a relation
\[
      {\cal D}_{q,P}\!\left(\prod_f z_f^{b_f}\right)\le Z,
      \qquad \mathbf b\not\equiv0\bmod q ,
\]
for an allowed order is
\[
    \ll_{m,B,\varepsilon} X^\varepsilon R^m K^{6m-2}
                 \left(\frac{Z\Delta^2}{K}
                          +\frac{Z\Delta}{\sqrt K}\right),
       \qquad Z\ge1 .
 \tag{SQ1}
\]
The bound is uniform in \(M\), the pool and the original
box dictionaries. No original radical cap is needed for this
finite count.
[Full edge-rank count on a single box (BX1)](ABC_SHARED_BOX_FULL_POWER_RANK.md#counting-all-relations-including-matchings)

### Distinct roots of the original and shifted forms

Fix the common residue lists and the first array. There are
at most \(K^{3m-1}\) possibilities for that array by recording
a spanning tree. When one potential of the second array is
unknown, each original label incident to it is, up to sign, a
slope-one affine function with root at a known neighbor.
The signed shifted label at the same position is obtained by
subtracting the fixed first label in the same sign convention.
Its root is the second-array neighbor plus a known offset
depending only on the two corresponding first-array potentials.

We may discard pairs having a zero \(D_e\) or any collision
between the roots of two distinct such forms at a common
vertex. The same statement is used in both choices of fixed
array. A zero \(D_e\) fixes the difference of two indexed
potentials in the second array. A root collision with distinct
neighbors also fixes the difference of those two second-array
potentials, with the first array and the indices fixed. For
the same neighbor, the two distinct forms are the original
and shifted form, whose roots differ by the nonzero first
label. To count a fixed difference of two second-array
potentials, omit one of them and record a spanning tree
on the remaining connected graph. Its \(3m-2\) labels
cost at most \(K^{3m-2}\). The difference equation recovers
the omitted potential without any extra choices. There are
at most \(R^m\) common residue lists and only boundedly many
choices of the two positions and the forms. Thus all these
exceptions cost \(O_m(R^m K^{6m-3})\).

### Recovering active forms

For a proposed relation choose a prime \(p\mid q\) and divide
the weights by their common \(p\)-power as in (BX1). At least
one reduced weight is a unit modulo \(p\), and valuations
outside the original defect support satisfy the reduced
relation modulo \(p\). Keep one array fixed so that at
least one active observable still varies. If a shifted
observable is active, either orientation works. Otherwise
keep fixed the array opposite an active original position.
We thus consider active original and shifted forms on one
variable complete tripartite array, with at most two forms
on a given pair of vertices. All active observables from
the fixed array are fixed integers.

If two active forms meet at a vertex of the variable array,
fix every other potential in it using \(3m-2\) tree labels.
Pool all fixed observables and all differences of the roots
of incident original and shifted forms. These roots are
pairwise distinct by the preliminary pruning. The pool
depends only on the known potentials and has polynomial
prime product. A prime outside the pool divides at most one
incident form. Hence for the two selected active forms, each
valuation off the original defect support and this pool is
divisible by \(p\). For an original form tested through its
least multiplier, also remove its high-valuation tail of
radical at most \(\Delta\). Its remaining multiplier
valuations are the negatives of its original valuations
modulo \(p\). Enumerating the original defect support together
with these at most two tails costs at most \(Z\Delta^2\).
The two selected canonical coefficients are supported on
this record and the fixed pool. The fixed-value binary
norm count in the one-vertex proof of (BX1) recovers the
unknown potential in \(X^\varepsilon\) choices after the
support is fixed, uniformly over orders and weights. For a
signed shifted form first fix its sign, among two choices.
The affine roots remain distinct, so the eliminated equation
has a nonzero constant term. Combining with the known
arrays yields the first term in (SQ1).

If no vertex meets two active forms, the active forms
constitute a matching. Choose one between vertices \(a,b\)
and fix every other potential of the variable array using
\(3m-3\) tree labels. The other active forms are now fixed.
Include their prime supports in the pool. Off this known
pool the selected form has \(p\)-power defect at most
\(Z\Delta\), considering the possible tail of the one form
tested through its multiplier. Choose a known vertex in
the remaining axis. If the active form is original, use the
three labels on its original triangle. If it is shifted,
use the three shifted labels of that triangle. They satisfy
the same signed additive equation by subtraction. At fixed
common residue lists, each dictionary has at most \(K\)
possible values, also after translating by its fixed first
label. For the shifted triangle, split into at most eight
sign patterns and use the resulting positive absolute values.
The arbitrary-order one-grid bound (OG6) counts at most
\(O(X^\varepsilon(K+Z\Delta K^{3/2}))\) triangles, and
their two incident lengths and signs determine \(a,b\)
relative to the known vertex. This gives the second term
of (SQ1), because the unrestricted contribution of the
two omitted vertices is \(K^2\). The \(1/K\) term and the
preliminary exceptions fit in (SQ1) as \(Z,\Delta,K\ge1\).
[Arbitrary-order grid counts (OG6)](ABC_ONE_POWER_GRID_ADDITIVE_GRAPHS.md#an-arbitrary-order-defect-on-a-single-letter)

## Exact rank for primitive direction odds

In the synchronized packet stratum of (CS3), let \(n=X^L\),
\(d=s_1+s_2+s_3\), and \(t=1/6-d>0\). Use its fixed positive
price and \(t\) margins. The forced synchronized pair count
has exponent \(md+(6m-2)t\), up to the small fixed tolerance
losses. Let \(a=0\) for original labels and \(a=1/K_0\) if
multipliers are used with \(M\ge K_0/L\). For any fixed margin
\[
                  0<\rho<t/2-a,
 \tag{SQ2}
\]
take the packet tolerances and the inner allowance
\(\varepsilon/L\) sufficiently small. By (SQ1), outside a
power-small fraction of the pairs, the valuation defect of
every allowed nonzero weight vector on all \(3E_0\) observables
is greater than \(n^\rho\). Indeed the two savings in (SQ1)
use \(\rho+2a<t\) and \(\rho+a<t/2\). Both follow from
(SQ2). All parameters are fixed at the budget before taking
the height limit. [Synchronized pair lower bound (CS3)](ABC_CELL_PRIMITIVE_SECANT_BOUNDS.md#fresh-direction-radicals-on-synchronized-actual-boxes)

In any triangle of a pair, its three \(D_e\) have the signed
additive relation. The primitive secant divides all three
by the same positive integer. Therefore the odds of any two
nonzero primitive coordinates have absolute value \(|D_e|/|D_{e'}|\)
on that triangle. Expressed as weights on the \(E_0\) shifted
positions, the odds from all triangles generate the sum-zero
integer lattice. In fact the graph joining two edge positions
when they belong to the same triangle is connected. For edges
in one role, change their endpoints one at a time, using an
edge to the third axis to link edges sharing an endpoint.
Each other role joins the first through a triangle. A spanning
tree of this graph gives an integral basis of edge differences
for the sum-zero lattice.

It follows that the subgroup generated by all primitive
secant odds modulo \(q\)-th powers has order \(q^{E_0-1}\)
over \(\mathbb Q\) on the generic pairs. Together with both
original arrays of observables it has order \(q^{3E_0-1}\).
The order is exact for composite \(q\) as well: the sum-zero
lattice modulo \(q\) is the kernel of the split sum map
\((\mathbb Z/q)^{E_0}\to\mathbb Z/q\), and all nonzero weight
vectors on the direct sum are tested by (SQ1).

The same two orders hold over \(F(\mu_q)\) for any field
\(F\) under the root-discriminant bound of (BX4) with this
\(\rho\) and the same allowed \(q\). In particular, for a
fixed repair order \(M\) include the prime support of \(M\)
in the prescribed pool and use
\(\operatorname{rd}(F)\le n^{\rho/2}\). The degrees of the
Kummer extensions obtained by adjoining the corresponding
\(q\)-th roots are then \(q^{E_0-1}\) and \(q^{3E_0-1}\),
respectively. The proof of (BX4) uses only the defect bound on
every nonzero vector and hence applies to the three sets of
observables.
[Base-change stability (BX4)](ABC_SHARED_BOX_FULL_POWER_RANK.md#stability-under-base-change)

This rank theorem concerns directions computed from the actual
primitive quotients, including when the rows are stored by their
least repairs. It imposes no small-radical hypothesis on the
physical differences or on their primitive normalization gcds.
It rules out a low-rank description of these secants over a
small-root-discriminant base in the simultaneous packet regime.
It does not supply an upper bound on the rank forced by the
original three radical caps.

## Any fixed number of synchronized arrays

The same margin works for any fixed \(k\ge2\) arrays on the
same residue lists. Use the \(k\) sets of original observables
and the \(\binom{k}{2}\) sets of differences for all ordered-by-index
pairs. With at least one multiplier array retain the requirement
\(q\mid M\). The finite exception bound is
\[
 \ll_{k,m,B,\varepsilon} X^\varepsilon R^m K^{k(3m-1)}
        \left(\frac{Z\Delta^2}{K}
                 +\frac{Z\Delta}{\sqrt K}\right).
 \tag{SQ3}
\]
Indeed choose a variable array occurring in an active form
and fix all \(k-1\) others. The moving forms at an incident
edge can include its original form and its differences with
several fixed arrays. For the same neighbor, two different
shifted roots coincide only when the corresponding complete
labels of two fixed arrays agree. Count any such equality
by making one of those two arrays variable and omitting an
endpoint, as in the \(D_e=0\) test. For distinct neighbors,
any collision again fixes a difference of two potentials
of the variable array. All these exceptions cost
\(O_{k,m}(R^m K^{k(3m-1)-1})\). The remaining one-vertex and
matching recoveries are exactly the same, with every other
array costing \(K^{3m-1}\) choices before the restriction.

If \(B_{\mathcal Q}\) counts the actual boxes on ordered
lists \(\mathcal Q\), Hölder gives
\[
 \sum_{\mathcal Q}B_{\mathcal Q}^k
       \ge \frac{(\sum_{\mathcal Q} B_{\mathcal Q})^k}
                   {(\#\{\mathcal Q\})^{k-1}}.
\]
Using (LB4) and the residue alphabet bounds, the exponent
is at least \(md+k(3m-1)t\) up to the tolerance losses.
Thus (SQ3) again proves a power-small exception under (SQ2).
The integer lattice argument for all primitive direction odds
now gives a joint subgroup of order
\[
       q^{\,kE_0+\binom{k}{2}(E_0-1)}
\]
together with the \(k\) original sets of observables. The
same exponent gives the Kummer degree under (BX4) with
the stated \(\rho\). All tolerances may depend on fixed \(k,m\).
Neither is taken to grow with the height.
[Forced box count (LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)
