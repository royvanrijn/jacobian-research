# Fixed-colour coherent cubic norm graphs

For a triangle the
[multiplicative surface](ABC_COHERENT_CUBE_CLASS_MORDELL_LANG.md)
has no positive-dimensional abelian translates away
from its boundary. For a longer cycle such translates
can lie in the same geometric open set. For imaginary
quadratic cube norms, a restriction on the real periods
of the coordinate elliptic curves still rules out
positive-rank cosets containing prime witnesses.

Fix a connected simple graph \(G=(V,\mathcal E)\) with
\(n=|V|\ge2\), \(m=|\mathcal E|\), and a direction on each
edge. Write \(\tau(G)\) for its number of spanning
trees. Fix a positive squarefree \(d\) and a positive
cube-free integer \(w_e\) for every edge. Put
\[
 E_e:y^2=x^3-dw_e^2/4,\qquad
 a_e=w_e\sqrt{-d}/2,\qquad t_e=(y_e+a_e)/(y_e-a_e).
\]
A coherent witness assigns distinct rational primes
\(p_i\) to the vertices and cube norm generators
\(\alpha_{ij}=(u_{ij}+v_{ij}\sqrt{-d})/2\) to the
directed edges. As in the
[prime path construction](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#witness-selection),
\(N(\alpha_{ij})=(p_i p_j)^3\), and
\(\alpha_{ij}/\bar\alpha_{ij}\) is a quotient of vertex
values, so its oriented product on a cycle is one.
The fixed-colour condition is \(v_e=w_e b_e^3\) for
nonzero signed integers \(b_e\).

## Degree of the graph locus

**Proposition MG1.** In \(\prod_e E_e\), the closed graph
locus obtained by requiring the \(t_e\) to be quotients
of vertex values is geometrically irreducible of
dimension \(n-1\) and has polarized degree
\[
                           3^{m-n+1}(n-1)!\tau(G)
 \tag{1}
\]
for the product of the origin polarizations.

Over an algebraic closure, let \(Z_G\) be the closure
in \((\mathbf P^1)^m\) of the image of the torus of
vertex values modulo a common factor, using
\(t_{ij}=z_j/z_i\). A spanning tree recovers all the
vertex quotients, so this torus has dimension \(n-1\).
Each \(t_e:E_e\to\mathbf P^1\) is finite flat of degree
three. Its Kummer equation over the torus is
\[
                   x_e^3=\frac{4a_e^2 t_e}{(t_e-1)^2}.
\]
For each edge \(e=(i,j)\), the divisor \(z_i=z_j\)
occurs with valuation \(-2\) in its right-hand side,
and valuation zero for every other edge, since \(G\)
is simple. Thus the \(m\) Kummer classes are independent
over the function field of the torus. The product
inverse image of \(Z_G\) has a field of degree \(3^m\)
as its generic fibre. Finite flatness ensures that
all its components dominate \(Z_G\), so it is
geometrically irreducible and generically reduced.

A choice of \(n-1\) distinct coordinate hyperplane
classes on \((\mathbf P^1)^m\) intersects \(Z_G\) in
one point generically if the corresponding edges are
a spanning tree, and has intersection number zero
otherwise. Indeed on the tree-coordinate torus the
other coordinates are single-valued monomials, while
a non-tree choice of \(n-1\) edges contains a cycle
and has a relation among the chosen coordinates.
Therefore the degree of \(Z_G\) for the sum of the
coordinate hyperplane classes is \((n-1)!\tau(G)\).
The pullback of that polarization is three times the
specified origin polarization. Multiplying by the
degree \(3^m\) of the cover and dividing by \(3^{n-1}\)
gives (1). For a triangle this recovers degree \(18\).

## An obstruction for real translates

**Lemma MG2.** Let \(E:y^2=x^3+A\) with \(A<0\) real,
\(a=\sqrt A\), and \(t=(y+a)/(y-a)\). Suppose
\(B_i\in E(\mathbf R)\), \(n_i>0\) are integers, and
\(\varepsilon_i\in\{1,-1\}\). If a nonempty product
\[
                   \prod_i t([n_i]P+B_i)^{\varepsilon_i}
 \tag{2}
\]
is constant as a meromorphic function on \(E\), then
two terms have the same \(n_i,B_i\) and opposite signs.
The terms can therefore all be paired and canceled.

Uniformize \(E(\mathbf C)\) by \(\mathbf C/\Lambda\)
using its real invariant differential. The real locus
is connected, since \(x^3+A\) has only one real root,
so each \(B_i\) has a real lift. The imaginary periods
form \(b\mathbf Z\) for some \(b>0\). The zero and pole
of \(t\) are opposite nonreal points of order three,
each of multiplicity three. Their imaginary coordinates
modulo \(b\) are the two distinct nonzero thirds:
a third with imaginary coordinate zero would be
represented by a real point. For \(n=n_i\) the zeros
and poles of \(t([n]P+B_i)\) therefore have nearest
lifts to the real axis at heights \(b/(3n)\) and
\(-b/(3n)\), in a fixed order independent of \(B_i\).
Choose a term of maximal \(n\) and a closest pole of
its underlying function. At that point no term of
smaller \(n\) can contribute and no zero of another
underlying function with the same \(n\) can occur.
For (2) to have zero valuation there, another pole
with maximal \(n\) and opposite \(\varepsilon\) is
required. At a shared pole the two arguments of \(t\)
are its same pole on \(E\), which implies equality of
the translations \(B_i\). Cancellation and induction
give the conclusion. Negative \(n_i\) can be brought
to this form using \(t(-P)=t(P)^{-1}\).

The restriction to real translates has geometric
content. If \(T\) is the nonreal point with \(t(T)=0\),
then \(\prod_{j=-1}^1 t(P+[j]T)=-1\). The divisor
telescopes and the value is \(-1\) at a nonzero point
of order two, where the other two factors are inverse.
Combining two such triples using a second translate
whose coset of \(\langle T\rangle\) is different gives
a product of six varying functions equal to \(1\).
This supplies an elliptic translate in the open
geometric six-cycle locus, with distinct coordinate
functions. It uses nonreal translates and does not
provide rational prime witnesses.

## Uniform count for bridgeless graphs

**Theorem MG3.** Suppose every edge of \(G\) lies on a cycle.
With \(d,(w_e)\) fixed, the number of coherent prime
witnesses on \(G\) is at most
\[
             C_G^{\,1+\sum_{e\in\mathcal E}\operatorname{rank}E_e(\mathbf Q)}
 \tag{3}
\]
for a constant depending only on the graph and no
restriction on prime sizes. In particular the
arithmetic count holds for graphs without triangles.

A witness defines a rational point on the graph locus by
\(x_e=p_i p_j/b_e^2,\ y_e=u_e/(2b_e^3)\). All these
coordinates are finite with \(x_e\ne0\), and the square
classes of the \(x_e\) are distinct on distinct edges.
The arithmetic condition excludes any positive-rank
intersection of a translated abelian subvariety of the
graph locus with the rational Mordell–Weil group that
contains this witness, as follows.

For two different positive cube-free parts \(v,w\),
the curves \(E_v,E_w\) are not \(\mathbf Q\)-isogenous.
They are isomorphic over
\(\mathbf Q(\zeta_3,(v/w)^{1/3})\), using the two
coordinate scalings by powers two and three of the
cube root. Its order-three automorphism over
\(\mathbf Q(\zeta_3)\) changes the isomorphism by a
nontrivial cubic CM unit. All the CM endomorphisms
are defined over \(\mathbf Q(\zeta_3)\), so every
nonzero geometric homomorphism between these curves
also changes by that unit. This proves the assertion.
For equal cube-free parts,
\(\operatorname{End}_{\mathbf Q}(E_w)=\mathbf Z\).

Now take a coset contained in the graph locus, with a
rational witness as base point, and suppose its
intersection with the product of the rational groups
has positive rank. The identity component of the
Zariski closure of the relative rational subgroup is
a positive-dimensional abelian variety over \(\mathbf Q\).
It admits a nonzero morphism from a coordinate
elliptic curve over \(\mathbf Q\): choose a nonzero
projection to that curve, then use its dual and a
polarization to get a morphism back. Restrict the
coset along its image. By the isogeny restriction just
proved, all varying coordinates are on the same
\(E_w\) and consist of rational translates of integer
multiples of the parameter. Pick a cycle containing
one varying coordinate. On that cycle the product
of \(t_e\) with the orientation signs is constant.
The fixed coordinates contribute nonzero constants.
Lemma MG2 gives two distinct varying edges with
identical \(x\)-coordinates at the base point, also
when a sign of a multiplying integer has been
absorbed by inversion. This contradicts the distinct
prime square classes. Hence the stated coset has
rank-zero intersection.

Apply [Gao–Ge–Kühne, Theorem 1.1](https://pmihes.centre-mersenne.org/item/10.5802/pmihes.26.pdf)
using dimension \(m\) and degree (1). The graph locus
is defined over \(\mathbf Q\), because quadratic
conjugation inverts all \(t_e\) and preserves the
vertex-quotient relation. The theorem gives at most
\(C_G^{1+\sum_e\operatorname{rank} E_e(\mathbf Q)}\)
cosets. Each coset containing a prime witness
has at most a graph-dependent number of rational
points, because their differences are in the rational
torsion subgroup of the product. For an elementary
uniform torsion bound here, every \(E_e\) acquires
the model \(y^2=x^3+1\) over a local extension of
degree at most six at \(5\), and also at \(7\).
Prime-to-residue-characteristic torsion injects into
the good reduction over the corresponding finite
residue fields, so the full rational torsion is
bounded absolutely by combining the two places.
Finally the edge square classes determine their
unordered prime pairs, so connectedness gives at
most two choices for the vertex assignment once one
edge orientation is selected. The \(b_e\) are then
determined up to sign by \(x_e=p_i p_j/b_e^2\).
This bounded multiplicity gives (3).

The bound concerns fixed colour tuples and a fixed
quadratic field. The original congruence-lattice
count varies these parameters, so an estimate for
their total occupied weight is still required. The
real-period lemma uses the imaginary quadratic sign
through the nonreal \(3\)-torsion zero and pole; it
does not give the corresponding count for real
quadratic norms.

## The remaining orientation condition at cut vertices

**Proposition MG4.** Suppose \(d>3\) and the distinct vertex
primes are unramified split primes in
\(\mathbf Q(\sqrt{-d})\). For arbitrary primitive edge cube
norms of the stated sizes, the \(t_e\)-cycle conditions
are equivalent to consistency of the prime orientation
at each vertex separately on each connected component
of \(G-v\). A single coherent orientation for the whole
graph exists exactly when these signs at \(v\) agree
across the components of \(G-v\), at every vertex \(v\).
With these signs agreed, changing the signs of the
generators as needed makes them coherent across all
paths. In particular the cycle test already forces the
global prime orientations for a graph without cut
vertices.

To see this, choose a reference prime ideal
\(\mathfrak P_v\) above \(p_v\) for each vertex.
For an edge directed out of \(v\), the valuation of
\(t_e\) at \(\mathfrak P_v\) is either \(+3\) or \(-3\);
reverse its sign if the edge is directed into \(v\).
For nonincident edges the valuation is zero. Vanishing
of every signed cycle sum of the valuations at
\(\mathfrak P_v\) is equivalent to the valuation
function on the directed edges being a difference of
vertex values. This follows by integration along
paths from one vertex. On edges of \(G-v\) the
function is zero, so the potential is constant on
each component of \(G-v\). This gives exactly the
component-wise sign condition on the incidences at
\(v\), and conversely such a potential is obtained
by assigning the specified constant on each component.

If this holds at every vertex, traverse a simple
cycle, using the conjugate edge generator whenever
the traversal reverses the edge. At each vertex
the two incident generators use opposite prime
ideals to exponent three. Their product around
the cycle therefore generates the rational ideal
\((\prod p_v^3)\) over the vertices of the cycle.
The only units for \(d>3\) are \(\pm1\), so the
product is rational and the \(t_e\)-cycle relation
holds. In the other direction a \(t_e\)-cycle
relation implies all the stated valuation conditions.
A global orientation requires that the sign depend
only on the vertex. With these signs consistent,
the adjacent oriented prime ideal classes differ by
\(3\)-torsion, hence all vertices are in one coset
by connectedness. The common anchor construction
gives coherent generators with the same ideals
as the given ones. Their ratios are units \(\pm1\),
which accounts for the allowed sign changes.

There are unbounded genuine cube norm examples of
the distinction. In any fixed \(d>3\), choose five
distinct odd rational primes split completely in
the Hilbert class field. Both orientations over
each prime are principal. Take \(G\) to be two
triangles meeting at one vertex, choose the
orientation at their common vertex oppositely in
the two triangles, and construct the principal
cube edges separately in the triangles. Every
cycle test holds, and the generators are
prime-primitive. The common vertex has incompatible
orientations on the two components of its deleted
graph, so the same edges have no simultaneous
global prime orientation. Chebotarev gives such
examples with arbitrarily large distinct primes.

## Rational cancellation for both quadratic signs

**Lemma MG5.** Let \(E:y^2=x^3+A\) over \(\mathbf Q\),
where \(A\in\mathbf Q^\times\) is nonsquare, and let
\(T=(0,\sqrt A)\). For rational translations \(B_i\),
the conclusion of Lemma MG2 holds for either sign of
\(A\). In particular, if a product of the form (2)
is constant with \(B_i\in E(\mathbf Q)\), then its
signed factors pair off exactly. The proof uses the
order of a pole modulo rational torsion.

Write \(F=E(\mathbf Q)_{\rm tors}\). This group is
cyclic: for odd \(p\) a full rational \(E[p]\)
would give \(\mu_p\subset\mathbf Q\) by the Weil
pairing, and full rational \(2\)-torsion is impossible
because \(x^3+A\) has only one real root.
Also \(T\) has order three and \(T\notin F\).
Two underlying factors can share a zero or pole
only if, for some choice of the signs on \(T\),
\[
                   [n_i]P+B_i=\pm T,\qquad
                   [n_j]P+B_j=\pm T .
\]
Eliminating \(P\) gives \([n_j]B_i-[n_i]B_j\) equal
to an integral multiple of \(T\). A rational
integral multiple of \(T\) is zero. Hence the
factors belong to the same block for the common
value of \(B_i/n_i\) in \(E(\mathbf Q)\otimes\mathbf Q\).

In each block change the parameter so that every
translation belongs to \(F\). In detail, let \(g\)
be the gcd of the multiplying integers in the
block. A Bézout combination of the \(B_i\) gives
a rational point \(B_\star\) with free class equal
to \(g\) times the common value. Choose \(D\) with
\([g]D=B_\star\) over an algebraic closure and replace
\(P\) by \(P-D\). The translations
\(U_i=B_i-[n_i/g]B_\star\) are now rational torsion.

Choose maximal \(n\) in the block, with its translation
\(U\). There is a torsion point \(P\) with
\([n]P=T-U\) whose image in \(E_{\rm tors}/F\) has
order exactly \(3n\). To verify this, at each prime
\(p\) use
\(E[p^\infty]\simeq(\mathbf Q_p/\mathbf Z_p)^2\)
and place the cyclic subgroup \(F[p^\infty]\) along
the first coordinate. At \(p\ne3\) the right-hand
side \(T-U\) has second coordinate zero. Choose the
second coordinate of a solution with exact order
\(p^{v_p(n)}\), and solve independently in the
first coordinate. If \(p=3\) and \(F[3^\infty]\ne0\),
the second coordinate of \(T\) is nonzero:
otherwise \(T\) would belong to the order-three
line of \(F\). Choose a solution in that coordinate
of order \(3^{v_3(n)+1}\). If \(F[3^\infty]=0\),
choose a coordinate in which \(T\) is nonzero to
do the same. The obtained orders also give the
exact order in the quotient, since \([n]P=T-U\)
implies \([3n]P\in F\).

At this selected point any other zero or pole
\([m]P+U'=\pm T\) implies that \([m]\bar P\) has
order three in \(E_{\rm tors}/F\). Hence
\(\gcd(3n,m)=n\), and maximality forces \(m=n\).
The two arguments have to equal the same sign
of \(T\), since their difference is rational.
Then \(U'=U\). A zero total valuation of the
product requires an identical factor with the
opposite exponent. Undoing the common translation
and iterating proves cancellation.

**Corollary MG6.** The fixed-colour graph bound (3)
holds for real quadratic cube norm witnesses as well,
using the curves \(E_e:y^2=x^3+dw_e^2/4\) for
positive squarefree \(d>1\) and \(a_e=w_e\sqrt d/2\).
More generally it bounds the rational points of the
graph locus with finite nonzero \(x_e\) whose edge
\(x\)-coordinates have distinct square classes.
The bound for rational points uses the same
Mordell–Weil rank sum.

The divisor and polarization calculation of
Proposition MG1 is geometric and hence unchanged.
On a positive-rank rational coset through a
point in question, the isogeny argument in MG3
again forces all varying coordinates in a
one-dimensional rational slice to have the
same cube-free part. Their translations are now
rational, so Lemma MG5 applies to every cycle on
which at least one coordinate varies. Its
cancellation would equate two edge \(x\)-coordinates
at the original point up to the inversion that
preserves \(x\). Distinct square classes exclude
this. The same finite torsion and uniform
Mordell–Lang count therefore applies. The
positive-squarefree restriction \(d>1\) ensures
that the \(3\)-torsion pole is nonrational; the
geometric six-cycle example shows why a claim for
arbitrary translations would require a different
condition.
