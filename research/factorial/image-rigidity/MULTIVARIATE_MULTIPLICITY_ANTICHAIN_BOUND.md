# Uniform antichain bounds for collision multiplicities

Fix an ordinary two-monomial core, possibly empty, and a collection
of ordinary unequal rows that may be appended. A multiplicity
profile uses a finite multiset from the collection.
Let \(J\ge4\) be an
[absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For each \(2\le n<N=n+1<J\), form the vector of ratios
\[
  X_{j,k}(a,b)=\frac{((j-k)a+kb)!}{(ja)!},
       \qquad j=n,N,\quad 1\le k\le j ,
\]
for a single row, with \(X_{j,0}=1\). Suppose these
vectors, for all the rows under consideration, generate
a group of rank at most \(R\) at each of the finitely
many pairs. Any collection of \(R\) named row types
satisfies the hypothesis. The rows can have arbitrary
orientations, heights, gaps and private endpoints.

**Theorem.** There is an explicit bound \(B(J,R)\),
independent of the core and of the rows' exponents, on
the number of different collision profiles along any
chain of strictly increasing multisets from the collection.
Here a collision profile has distinct exponent columns and a
common root at some adjacent moment pair. Equal rows are omitted.
For a finite list of \(r\ge1\) types the exceptional count
vectors are contained in a union of at most \(B(J,r)\)
rational affine hyperplanes with strictly positive
normals. More precisely, their set is a union of at most
that many lattice cosets intersected with the admissible
integer counts, each an antichain for componentwise order. Thus any multiplicity box of side \(L+1\)
contains at most \(B(J,r)(L+1)^{r-1}\) exceptional
vectors, uniformly in the exponent heights.

The rank bound concerns the multiplicative group of whole
coefficient vectors, rather than the number of primes
in the entries. The theorem gives a bound on the number
of exceptions on a chain. It does not exclude a single
exceptional profile or provide a theorem for chains whose
group rank is unbounded.

## Resultant equations in the row multiplier group

At fixed \(n,N\), let \(C_{j,k}\) be the constant-normalized
core coefficients including the binomial factors.
For a count vector or multiset \(m\), the full moment
coefficients after constant normalization are
\(C_{j,k} X_{j,k}(m)\), where \(X(m)\) is the
componentwise product of the row vectors with their
multiplicities. This identity follows directly from
the factorial moment formula. The physical collision is
equivalent to
\[
 P(X(m))=0,\qquad
 P(X)=\operatorname{Res}_z\!\left(
      \sum_{k=0}^n C_{n,k} X_{n,k}z^k,\
      \sum_{k=0}^N C_{N,k} X_{N,k}z^k\right).
\]
This uses the original common argument. The Sylvester
formula makes \(P\) a rational polynomial in \(D=n+N\)
variables of total degree at most \(D\), since each
matrix entry is either a constant, zero or a scalar
multiple of one variable. It is nonzero. For instance
one can specialize all intermediate coefficient variables
to zero and the two leading coefficients to one. This
gives \(1+z^n\) and \(1+z^N\), which have no common
root since a common root would have \(z^{N-n}=1\).
In particular the number of monomials in \(P\) is
at most \(\binom{2D}{D}\).

Use the nonzero terms \(c_\alpha X^\alpha\) of \(P\).
At any zero, partition the terms into minimal
vanishing blocks, choosing one anchor \(\alpha_0\)
per block. Each block has at least two members.
For that block the tuple of monomial ratios
\[
                     \bigl(X^{\alpha-\alpha_0}\bigr)
                               _{\alpha\ne\alpha_0}
\]
belongs to the image of the row group, which has rank at
most \(R\). The tuple solves a linear equation with right
side one with nonzero rational coefficients. Minimality
of the vanishing block makes the solution nondegenerate,
because a vanishing subsum of the terms not including the
anchor would be a smaller vanishing block.
By [Evertse, Schlickewei and Schmidt, Theorem 1.1](https://arxiv.org/pdf/math/0409604),
a block with \(s\) nonanchor terms has at most
\(\exp((6s)^{3s}(R+1))\) such tuples. When \(s=1\) the
tuple is uniquely determined. The entire block partition
and the corresponding tuple at every block constitute
finite data which can occur on only a bounded number
of classes. For example, let \(L=\binom{4J}{2J}\).
The number of partitions for a fixed pair is at most
\(L^L\), by labeling blocks with elements of the set
of terms. The number of joint tuple choices is at most
\(\lceil\exp(L(6L)^{3L}(R+1))\rceil\).
Their product times \(J\) is a valid \(B(J,R)\).

## Curvature forces every collision class to be an antichain

Suppose two strictly nested multisets \(m_0\subsetneq m_1\)
have the same data at one adjacent pair. Let \(q\) be
their difference multiset. Then every anchored ratio has
value one on \(X(q)=X(m_1)/X(m_0)\).
For every integer \(t\ge0\), the ordinary profile
\(m_0+tq\) has those same anchored ratios. Dividing each
resultant block by the anchor term therefore shows that
each block still vanishes. In particular all these exact
ordinary profiles would have a collision at the same pair.
But \(q\) contains an unequal row. The total curvature is
additive and grows without bound in \(t\), contradicting
[absolute curvature exclusion at \(S\ge12\)](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md).
Thus each class meets a strictly increasing chain at most
once. Terminal two separates by even real positivity,
and terminals at least \(J\) are excluded by the cutoff.
Summing the class bounds at the other terminals proves
the chain claim.

For finitely many types, the preimage of fixed anchored
ratios is either empty or a coset of the kernel of a
homomorphism from \(\mathbb Z^r\) into positive rational
tuples. The kernel is a saturated lattice because its
image has no torsion. Its real span \(V\) is rational.
For each coset meeting the admissible counts, the
argument above gives \(V\cap\mathbb R_{\ge0}^r=\{0\}\).
Indeed any nonzero point in that intersection would
give a nonzero rational point, by rational linear
feasibility after normalizing the sum of coordinates to
one. Clearing denominators gives a nonnegative integral
kernel vector, by saturation. Repeating it from any
nonnegative member of the coset contradicts curvature.
Every admissible member of a tuple class is an exact
collision, since its block sums vanish.

Finally \(V^\perp\) has a strictly positive rational
vector. For example project the compact standard simplex
orthogonally to \(V^\perp\). Its convex image does not
contain zero since \(V\) misses the simplex. The nearest
point of that image to zero has positive inner product
with the projection of every vertex, by differentiating
squared distance along line segments. Its coordinates
are therefore all strictly positive. Approximate it
by a rational vector of \(V^\perp\), keeping all
coordinates positive. Its value is constant on the
corresponding coset, yielding the stated hyperplane.
A fixed-coordinate line in the box is a monotone chain,
so counting its other \(r-1\) coordinates proves the box
bound.

In particular the result limits exact collisions uniformly
over growing exponents with finitely many multiplicity
generators. A general argument that removes the remaining
isolated classes is still required for SFC(2).


For a single multiset increment \(q\ne0\), there is also
an elementary bound along \(m_0+tq\). At fixed orders
the resultant evaluated on these counts is
\(\sum_\alpha a_\alpha \rho_\alpha^t\) with
\(\rho_\alpha>0\). Combine equal bases and remove zero
terms. The remaining exponential polynomial cannot
vanish identically, since integer repetitions eventually
reach the absolute curvature exclusion. An exponential
polynomial with \(s\) distinct real exponents has at
most \(s-1\) distinct real zeros: divide by the lowest
exponential and apply Rolle induction to its derivative,
which has one fewer term. Thus along this single
increment ray one may use \(J(L-1)\) exceptions, where
\(L\) is the monomial bound above.
