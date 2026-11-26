# Long Newton edges give angular separation through Galois pair orbits

Let \(P_n,P_N\) be the ordinary moment polynomials of any
two distinct exponent columns, with \(N=n+1\ge3\).
Suppose \(P_N\) is irreducible over \(\mathbb Q\).
A Newton polygon edge is called primitive here if its
horizontal width \(w\) and integer vertical rise \(h\)
satisfy \(\gcd(w,h)=1\). The rise may have either sign.

**Theorem.** Suppose that at a prime \(p\), the upper
polynomial has a primitive edge of prime width
\[
                               N/2<w<N,\qquad p\ne w .
\]
Then the two moments have disjoint root arguments under
independent positive rescalings.

More generally, allow a finite collection of primitive
upper edges, at possibly different primes, with widths
\(w_i\) satisfying
\[
    N/2<w_i<N,\quad p_i\nmid w_i,\quad
    \gcd(w_i,(N-w_i)!)=1,\quad \gcd(N,w_1,\ldots,w_r)=1 .      \tag{1}
\]
The same conclusion holds. In place of the factorial
gcd condition on \(w_i\), it suffices that the corresponding
local splitting automorphism can fix the complementary
roots while acting as a \(w_i\)-cycle on the edge, or that
its permutation of the complement has order coprime
to \(w_i\).

This is a test on the original factorial moments and their
actual native primes. It uses a long edge that need not
be the endpoint chord. When it or other exact separations
hold at \(3\le N<J\), the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
supplies every remaining terminal. The first pair is
covered by strict positivity of the quadratic on the
real line.

## A long cycle from a primitive edge

Fix an edge at \(p\) from index \(a\) to index \(b\),
with \(w=b-a\) and slope \(h/w\) in the
constant-to-leading convention. Set \(\pi^w=p\) and
substitute \(z=\pi^{-h}x\). Divide the polynomial by
the appropriate power of \(\pi\) to make its minimum
coefficient valuation zero. The reduction is
\[
                        x^a(c_a+c_bx^w)
\]
with \(c_a,c_b\ne0\) in \(\mathbb F_p\).
Indeed, the supporting line can meet no coefficient
point strictly between the endpoints, since its slope
has reduced denominator \(w\).
The edge is taken maximal, so there is no further
coefficient on its supporting line.

Work over a finite unramified field \(K\) containing
\(\mu_w\) and the residue roots of the displayed binomial.
Since \(p\nmid w\), the \(w\) nonzero residue roots are
simple. Each lifts by Hensel's lemma in \(K(\pi)\).
The physical root associated to each lift has
valuation \(-h/w\), with reduced denominator \(w\).
It therefore generates \(K(\pi)\) over \(K\), since its
degree is at least \(w\) by its valuation and
\([K(\pi):K]=w\). The extension \(K(\pi)/K\) is cyclic
because \(K\) contains \(\mu_w\).
Its generator gives one cycle on the \(w\) roots
corresponding to the edge.

Extend the generator to the local splitting field.
It permutes the complementary \(N-w\) roots.
If the order of that permutation is coprime to \(w\),
raising the automorphism to that order fixes the
complement and leaves a single \(w\)-cycle.
In particular this holds under (1).
When \(w\) is prime and \(w>N/2\), the complementary
permutation acts on fewer than \(w\) points and its
order is prime to \(w\).
The local permutation restricts to an element of the
global Galois group, since it permutes all roots and
fixes their rational coefficient field.

## Why the global group is transitive on unordered pairs

Write \(G\) for the global Galois group. The irreducibility
hypothesis makes \(G\) transitive on its \(N\) roots.
Any nontrivial system of blocks has a common block size
\(b\) dividing \(N\), with \(1<b\le N/2\).
For a cycle of width \(w_i>N/2\), a block containing a
fixed vertex is fixed setwise by the cycle.
Such a block cannot also contain a moving vertex, since
then it would contain the whole orbit of size \(w_i>b\).
Thus the moving vertices form a union of blocks, so
\(b\mid w_i\) for every available long cycle.
The final gcd condition in (1) rules out every nontrivial
block system. For a single prime width \(N/2<w<N\),
this gcd condition follows from \(N<2w\).

Now fix any orbit on unordered pairs and regard it as
the edge set of a graph on the \(N\) roots.
The graph is nonempty and regular by transitivity.
Every nonempty \(G\)-invariant graph is connected in a
primitive action: otherwise its connected components
are a nontrivial block system, since positive degree
and transitivity exclude isolated vertices.
At a vertex fixed by any one long cycle, its neighbors
in the long orbit are either all of that orbit or none.
If the graph's degree were below \(w_i\), every vertex
fixed by the cycle would have no neighbor in the long
orbit. This would disconnect the graph, since both
sets are nonempty. Thus any nonempty invariant graph
has degree at least \(w_i\).
If the orbit graph were proper, its nonempty complement
would also be invariant and have degree at least \(w_i\),
contradicting \(2w_i>N-1\). Hence there is only one
orbit on unordered pairs.
This also covers a long cycle with exactly one fixed
vertex.

The
[Galois quotient criterion](MULTIVARIATE_GALOIS_ORBIT_ANGULAR_SEPARATION.md)
now gives angular separation. Its upper unordered
quotient traces are all distinct and each has degree
\(\binom N2\), exceeding the corresponding bound
\(\binom{N-1}{2}\) at the lower order.

## Applying the test to a delayed factorial mark

A useful source of a long edge is a native delayed
mark at an offset \(h_0\) in one increasing row of gap
\(d\), with first sampled index \(j=\lceil h_0/d\rceil\).
At a prime \(p>N\), remove all row factorial scalars
\((N\ell_i)!\). Suppose the marked factor has valuation
\(E>0\) and all other sampled row factors are units.
Then the binomial coefficients are units. The coefficient
valuations before \(j\) are zero and thereafter equal
\(E\). The upper polygon has a horizontal edge of width
\(j-1\) and one edge of width \(w=N-j+1\) and rise \(E\).
The latter is primitive when \(\gcd(E,w)=1\).
An originally decreasing mark is handled by exchanging
the two columns.

For example, take the univariate ordinary profile
\((38,39)\) at \(N=9\). The rising factors are
\(343,\ldots,351\); at \(p=347\) there is precisely
one mark at offset five. The positive edge has width
five and rise one. Irreducibility of the upper polynomial
follows by reduction modulo \(11\): its monic reduction is
\[
   z^9+2z^8+7z^7+8z^6+8z^5+5z^4+8z^3+4z^2+4z-1.
\]
In the residue quotient, the third iterated
\(11\)-th power of \(z\) is
\[
   9z^8+3z^7+2z^6- z^5+2z^4+5z^3+2z^2+7z-2 .
\]
Its difference from \(z\) is coprime with the reduction,
and the ninth iterated power equals \(z\).
These are the degree-nine finite-field irreducibility
conditions. Thus the long-edge test proves angular
separation at this pair using an index-five mark.
The example illustrates the certificate at a single
pair and does not require a search through terminals.
