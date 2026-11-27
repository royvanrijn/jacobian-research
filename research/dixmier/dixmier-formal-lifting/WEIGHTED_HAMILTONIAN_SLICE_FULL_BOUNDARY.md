# Zero-sum weights and the full Newton boundary for rational slices

The positive-sum [pole test](WEIGHTED_PLANE_HAMILTONIAN_SLICES.md#a-finite-test-at-arbitrary-signed-weights)
decides exactness for a fixed weighted homogeneous Hamiltonian.
Changing the sign of the weight preserves homogeneity of a
fixed polynomial, with the sign of its degree reversed. For a
full polynomial, changing the sign of the weight changes which
face is exposed. The following calculation supplies zero-sum
weights and distinguishes these two uses of the weight.

**Proposition ZeroSumWeightedHamiltonianSlice.**
Let \(k\) have characteristic zero and let \(p\in k[x,y]\)
be nonconstant and homogeneous for a nonzero weight of sum
zero. Then a rational \(G\) with \(J(p,G)=1\) exists exactly
when the weighted degree of \(p\) is nonzero. A rank-one Weyl
operator with a nonzero scalar commutator mate satisfies the
same necessary condition whenever its principal symbol in a
zero-sum weight is nonconstant.

It is enough to treat the primitive weight \((1,-1)\), since
the other orientation follows by exchange of variables.
Write \(t=xy\) and \(d=\deg_{(1,-1)}p\). Then \(J(p,t)=dp\),
so \(t/(dp)\) is a slice when \(d\ne0\). For \(d=0\) we have
\(p=P(t)\), and in \(k(t)(x)\) the equation is
\[
                  -P'(t)x\,\partial_xG=1 .
\]
The associated differential has nonzero residue at \(x=0\).
This decides the formerly missing weight sum for every
weighted homogeneous plane polynomial over \(k\).
The pole test applies to negative-sum homogeneous symbols
by reversing their weight.

For the Weyl assertion write \([y,x]=1\) and \(N=xy\).
The ordered monomials of weight zero are \(x^j y^j\), which
are the falling factorials in \(N\), by induction using
\(xN=(N-1)x\). Thus the degree-zero subalgebra is \(k[N]\).
The same subalgebra is described by weight zero in symmetric
ordering, which preserves the grading and the scalar subspace.
If a member of a scalar pair has degree zero and nonconstant
principal symbol, its highest operator is a nonconstant
\(H(N)\). Let the other member have degree \(b\) and highest
operator \(Q_b\). If \(b<0\), the scalar commutator is impossible
since all its terms have weight less than zero. If \(b=0\),
the highest two operators commute in \(k[N]\) and again all
other terms have strictly negative weight. If \(b>0\), then
\([N,Q_b]=bQ_b\), so
\[
       [H(N),Q_b]=Q_b\bigl(H(N+b)-H(N)\bigr)\ne0 .
\]
This weight-\(b\) term is unique from the two highest
operators. It is nonzero by characteristic zero and the
domain property of the Weyl algebra. This excludes the
zero-degree case in both orientations.

**Corollary PlaneKellerFullNewtonBoundarySliceTest.**
Let \(F,Q\in k[x,y]\) satisfy \(J(F,Q)=\kappa\in k^\times\).
For every nonzero integral weight, every nonconstant initial
form of either component has a rational Hamiltonian slice.
In particular this applies to each non-origin vertex and
each edge of the Newton polygon of the component, using
the convex hull of its support without adjoining the origin.

Over an algebraically closed field, a lattice polygon with
prescribed nonzero vertex coefficients admits a polynomial
passing these tests at every nonzero weight if and only if
every non-origin vertex is off the diagonal and every edge
admits a passing homogeneous polynomial with those endpoints.
This last criterion concerns the boundary tests.

To verify necessity for the plane Keller pair, take any
weight \((r,s)\) and put \(u=\deg F,\ v=\deg Q\),
\(n=u+v-r-s\). The rescaled polynomials
\(F_\epsilon=\epsilon^u F(\epsilon^{-r}x,\epsilon^{-s}y)\)
and \(Q_\epsilon=\epsilon^v Q(\epsilon^{-r}x,\epsilon^{-s}y)\)
satisfy \(J(F_\epsilon,Q_\epsilon)=\kappa\epsilon^n\).
Both are polynomials in \(\epsilon\), so \(n\ge0\).
When the leading form \(p\) of \(F\) is nonconstant and the
first surviving coefficient \(q_j\) of \(Q_\epsilon\) has
\(j<n\), it satisfies \(J(p,q_j)=0\). It is algebraic over
\(k(p)\) by the characteristic-zero differential criterion.
Lift it by the simple-root formal recursion in
\(k(x,y)[[\epsilon]]\), substituting \(F_\epsilon\) for
\(p\) in the coefficients of its separable minimal polynomial.
Its lift commutes with \(F_\epsilon\) for the Jacobian
bracket by differentiating that polynomial and cancelling its
invertible derivative. Subtract \(\epsilon^j\) times the lift.
After finitely many subtractions the coefficient of
\(\epsilon^n\) has bracket \(\kappa\) with \(p\).
Dividing by \(\kappa\) gives the claimed slice. This
argument uses ordinary polynomial multiplication for every
sign of \(r+s\). The existing
[Weyl theorem](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
uses positive weight sum, and the graded proof above supplies
the zero-sum case.

The diagonal is excluded at any non-origin vertex by the
monomial case of the pole test. For the gluing direction,
use the [endpoint normalization argument](WEIGHTED_PLANE_HAMILTONIAN_SLICES.md#gluing-the-tests-around-a-newton-polygon)
on every edge. Its scaling step depends on the primitive
lattice direction of the edge, so the coefficient ratio
can be matched at any weight sum. Different edges share
only endpoint coefficients. Vertex and edge tests therefore
give the stated equivalence, also for segments by testing
the two endpoints and the full segment.

For example \(x^2(x+y)^3+xy\) passes the positive-sum tests.
The edge from \((2,3)\) to \((5,0)\) is the only edge
with outward normal of positive sum. Its face polynomial
is \(x^2(x+y)^3\), a monomial in the polynomial coordinates
\(x,x+y\) with distinct exponents. The other two outward
normals are \((-2,1)\) and \((-1,-4)\), both of negative
sum. Their common vertex \((1,1)\) is exposed by
\((-1,-1)\) and contributes \(xy\), which fails the
rational-slice test. This separates the two boundary criteria.
