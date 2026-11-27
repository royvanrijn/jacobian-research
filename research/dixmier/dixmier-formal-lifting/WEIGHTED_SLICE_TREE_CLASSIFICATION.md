# Root multiplicities and plane trees for Hamiltonian slices

For positive-degree mixed weights and differential-order faces,
a rational primitive gives a Belyi function whose zero and pole
multiplicities determine a weighted plane tree. These trees
classify the passing root positions up to change of parameter.
For negative-degree mixed weights, a tree constructs passing
root positions exactly when the pole budget permits them.

Use the following consequence of
[Pakovich and Zvonkin, Propositions 2.8 and 2.11 and Theorems 2.24 and 3.1](https://arxiv.org/pdf/1306.4141v1).
Take two partitions of the same integer \(N>0\), with \(V\) parts
between them and greatest common divisor one. Rational functions
whose fibres over \(0,\infty\) have these partitions and whose
fibre over \(1\) has partition \((V-1,1^{N-V+1})\) correspond to
weighted bipartite plane trees. The parts are the vertex degrees,
where the degree is the sum of incident edge weights. The
correspondence uses \(S/(S-1)\) as the Belyi function with
vertex values \(0,1\) and face value \(\infty\). The edge weights
are positive integers. The cyclic order at each vertex is part
of the tree. Existence is equivalent to \(N\ge V-1\).
The functions are taken up to a Möbius change of source.
For \(V>2\), the face of degree \(V-1\) distinguishes its
preimage of \(1\). For \(V=2\), the function has one zero
and one pole, and its deck group acts transitively on the
preimages of \(1\), so choosing one gives the same isomorphism
class. These realizations can be taken over
\(\overline{\mathbb Q}\), as in Section 2.3 of that source.

## Differential-order faces

**Theorem LeadingFactorSliceTreeClassification.**
Let \(g\) be a positive integer and let \(m_1,\ldots,m_s\)
be positive integers, with \(s\ge2\). Over an algebraically
closed field of characteristic zero there is a polynomial \(h\)
with these distinct-root multiplicities such that \(dx/h^{1/g}\)
has a rational primitive in \(k(x,h^{1/g})\) exactly under the
following conditions.

First \(M=\sum m_i/g\) is an integer at least \(s\) and
no \(m_i=g\). Form the list of nonzero signed integers
\[
                  g-m_i\ (1\le i\le s),\qquad
                  M-s\text{ copies of }g.
\]
Let \(\Delta\) be their positive greatest common divisor and
\(N\) the sum of the positive entries divided by \(\Delta\).
The final condition is \(N\ge M-1\).
The admissible \(h\), modulo affine changes of \(x\) and
multiplication by a nonzero constant, are indexed over
\(\mathbb C\) by the weighted bipartite plane trees for
the positive and negative parts of the divided list.
Every class has an algebraic representative.

The [leading-factor residue calculation](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md)
gives the equivalence between exactness and the polynomial
equation \(B'-(h'/gh)B=1\). Indeed its trace and pole argument
uses only exactness for this equivalence. It also gives the
first conditions, a simple root of \(B\) at every root of
\(h\), and \(M-s\) additional simple roots. The
[rational-map calculation](POLYNOMIAL_WEYL_LEADING_FACTOR_RESIDUE.md#a-rational-map-from-the-residue-equation)
gives \(S^\Delta\) proportional to \(B^g/h\), with \(S(\infty)=1\).
The divided list is the finite signed divisor pattern of \(S\).
Outside its zeros and poles, only infinity can ramify, with
index \(M-1\). The resulting tree has \(V=M\).

Conversely use a tree with this degree list to realize \(S\).
Put the distinguished point of \(S^{-1}(1)\) at \(x=\infty\).
The logarithmic derivative \(S'/S\) has simple poles at the
finite zeros and poles of \(S\) and no other finite poles or
zeros. At infinity it has order \(M\) in \(1/x\). Thus
\[
                    B=\frac{g}{\Delta}\frac{S}{S'},\qquad
                    h=\frac{B^g}{S^\Delta}
\]
are polynomials. The finite order of \(h\) at a point of
signed multiplicity \(r/\Delta\) is \(g-r\), which is \(m_i\)
or zero according to the list. Logarithmic differentiation
gives \(B'-(h'/gh)B=1\).

For a fixed admissible \(h\), \(B\) is unique. A nonzero
difference of two solutions would have \(g\)-th power
proportional to \(h\). The perfect-power exclusion for
multiple roots in the leading-factor note rules this out.
Equivalently, in the perfectly \(g\)-th-power pattern all
\(s\) original entries are negative and the \(M-s\) positive
entries are \(g\). Then \(N=M-s<M-1\).
The normalization \(S(\infty)=1\) determines \(S\)
from \(h,B\), and two source coordinates placing the
distinguished point at infinity differ by an affine change.
The tree existence theorem and its algebraic realizations
give the sufficiency over any algebraically closed field of
characteristic zero. For necessity, any fixed candidate
and the coefficients of an exactness certificate lie in a
finitely generated field that embeds into \(\mathbb C\).

For a face \(p=y^g h(x)\) at the differential-order weight,
the relative differential on the generic level is a scalar
multiple of \(dx/h^{1/g}\) after adjoining a \(g\)-th root
of the generic level to the constants. Exactness descends
by normalized trace. For constant \(h\) the differential
is exact. With a single root of multiplicity \(m_1\),
direct integration gives exactness precisely when
\(m_1\ne g\). Thus the multiple-root theorem completes
the multiplicity test at the axis weight.

## Mixed weights of positive degree

**Theorem MixedPositiveWeightSliceTreeClassification.**
Use the notation of the
[mixed-weight interpolation test](WEIGHTED_PLANE_HAMILTONIAN_SLICES.md#strictly-mixed-weights-and-positive-degree).
Assume that the symbol has at least two monomials and
the weights, after swapping variables if necessary,
are \(r<0<s\), \(r+s>0\), and the degree is positive.
Fix the endpoints of the support segment and the
multiplicities \(m_1,\ldots,m_j\) of the nonzero roots of \(G\),
with \(\sum_i m_i=\deg G\) as prescribed by the endpoints.
Let \(m=\beta-\nu\rho\) and write \(\nu=a/b\) in lowest
terms, with \(a,b>0\).

Some distinct nonzero root positions admit a slice exactly
when all of the following hold. The number \(m\) is an integer
at least \(j\). All \(\ell_i=b-a m_i\) and
\(\ell_\infty=a\deg G-bm\) are nonzero.
In the list
\[
                \ell_1,\ldots,\ell_j,\quad
                m-j\text{ copies of }b,\quad \ell_\infty
\]
divide by the positive greatest common divisor \(\Delta\).
If \(N\) denotes the sum of the positive divided entries, the
remaining condition is \(N\ge m\).
The root positions that pass, up to scaling the invariant
parameter \(t\), are indexed over \(\mathbb C\) by the
corresponding weighted bipartite plane trees with a marked
vertex of signed degree \(\ell_\infty/\Delta\).
Every class has an algebraic representative.

By the interpolation test a slice is equivalent to a polynomial
\(A\) of degree at most \(m\) satisfying
\[
                            A'-\frac{a}{b}\frac{G'}{G}A=t^{m-1}.      \tag{1}
\]
The polynomial has simple zeros at the roots of \(G\);
the nonresonance condition there is \(\ell_i\ne0\).
If it vanished at \(0\), the equation there would force
vanishing order \(m\), leaving no degree for any root of \(G\).
Any further finite zero is simple, by (1) at a nonzero
ordinary point. Taking the logarithmic derivative of
\(A^b/G^a\) gives \(b t^{m-1}/A\). Since the logarithmic
derivative of a rational function is \(O(1/t)\) at infinity,
\(\deg A=m\), and the nonzero coefficient of \(1/t\)
forces \(\ell_\infty\ne0\).

Take the rational \(\Delta\)-th root \(S\) of a constant
times \(A^b/G^a\), normalized to \(S(0)=1\).
Its divisor has exactly the displayed pattern divided by
\(\Delta\), including the signed order at infinity.
Away from its zeros and poles the log derivative
\(S'/S=(b/\Delta)t^{m-1}/A\) vanishes only at \(t=0\),
to order \(m-1\). Thus the distinguished preimage of \(1\)
has index \(m\), and the vertex count is \(V=m+1\).
The tree criterion gives \(N\ge m\); the vertex at infinity
must be marked.

For the converse, take such a marked tree. Put its
distinguished preimage of \(1\) at \(0\) and the marked
vertex at infinity. Put
\[
                          A=\frac{b}{\Delta}t^{m-1}\frac{S}{S'} .
\]
This is regular and nonzero at \(0\). It has simple zeros
at the finite points of the divisor of \(S\), no other
finite poles or zeros, and degree \(m\) from the nonzero
signed order of \(S\) at infinity. At a finite point with
signed order \(\ell_i/\Delta\), put a root of multiplicity
\(m_i=(b-\ell_i)/a\) in \(G\). The finite entries equal to
\(b\) are exactly the additional roots of \(A\); actual root
entries satisfy \(\ell_i<b\). The full divisor identity
makes \(A^b/G^a\) proportional to \(S^\Delta\), proving (1).
The two marked points fix the source coordinate up to scaling.
When \(m=1\) the two-vertex deck symmetry described above
still respects the marked vertex.

To check uniqueness of the tree for fixed root positions,
a nonzero difference between solutions of (1) would be a
constant times \(G^{a/b}\). This requires every \(b\mid m_i\)
and \(k=a\deg G/b\le m\). Equality \(k=m\) would give
\(\ell_\infty=0\). When \(k<m\), every original root
entry is negative since \(a m_i/b\) is a positive integer
different from one. The only positive entries are the
\(m-j\) copies of \(b\), with \(m-j>0\) and \(\Delta=b\).
This gives \(N=m-j<m\), again impossible.
The field assertions follow as in the differential-order case.

For example, with \(m=2\), \(a/b=1/2\), and multiplicities
one and four, the divided entries are \(1,-2,1\), with
the last entry marked. The unique tree gives a single
configuration up to scaling. The representative
\(G=(t-1)(t-2)^4\) has \(A=-2(t-1)(t-2)\) in (1).

Several exposed edges can be realized at once with their prescribed
root multiplicities. Apply this theorem on each positive-degree
mixed edge and use the
[boundary gluing theorem](WEIGHTED_PLANE_HAMILTONIAN_SLICES.md#gluing-the-tests-around-a-newton-polygon)
on a polygon satisfying its other face and vertex conditions.
The vertex normalizations rescale the invariant parameters,
which preserve the root multiplicities and the slice property.

## Mixed weights of negative degree

**Theorem MixedNegativeWeightSliceMultiplicityExistence.**
Use the notation of the
[negative-degree interpolation test](WEIGHTED_PLANE_HAMILTONIAN_SLICES.md#strictly-mixed-weights-and-negative-degree).
Fix the endpoints and multiplicities of the nonzero roots
of \(G\), so \(\mu>0\) and the \(\alpha_i>0\) are fixed
positive rationals. Put
\[
  M=\mathbf 1_{\mu\in\mathbb Z}\mu+
                       \sum_{\alpha_i\in\mathbb Z}\alpha_i .
\]
There are distinct nonzero root positions admitting a
rational Hamiltonian slice if and only if
\[
                       M\ \ge\ \#\{i:\alpha_i\notin\mathbb Z\}.
\]
The positions can be chosen algebraic. This is an existence
criterion for the multiplicity pattern. The interpolation
test still decides given root positions.

The necessity is the nonintegral-root bound in the
interpolation theorem. When all the exponents are integral,
its direct integration works at arbitrary root positions.
Otherwise let \(b\ge2\) be their least common denominator,
let \(j\) be the number of nonintegral \(\alpha_i\), and
assume \(e=M-j\ge0\). Put \(N=b(\mu+\sum_i\alpha_i)\).
We construct a polynomial \(S\) of degree \(N\), branched
at most over \(0,1,\infty\), with the following finite
fibres.

Over \(0\), put one point of degree \(b\mu\) if \(\mu\) is
nonintegral, one of degree \(b(\alpha_i+1)\) for every
nonintegral \(\alpha_i\), and \(e\) points of degree \(b\).
The sum of these degrees is \(N\) by \(e=M-j\).
Over \(1\), reserve one point of degree \(\mu\) if \(\mu\)
is integral and one of degree \(\alpha_i+1\) for every
integral \(\alpha_i\). Fill the partition to total \(N\)
with points of degree one.

To justify that the latter prescription fits, let \(p\) be the
number of points over \(0\) and let \(q_0\) be the number
of reserved points over \(1\).
Their total ramification defect is
\[
              \mathbf 1_{\mu\in\mathbb Z}(\mu-1)
                         +\sum_{\alpha_i\in\mathbb Z}\alpha_i
                 = p-1 .
\]
Each reserved point contributes at least one to this defect,
except the possible reservation \(\mu=1\), so \(q_0\le p\).
On the zero side, all degrees other than the possible marked
degree \(b\mu\) are at least two, giving \(N\ge2p-1\).
Therefore the reserved degrees sum to \(p-1+q_0\le N\).
After filling with degree-one points the two partitions have
\(N+1\) parts in total. Their gcd is one, since a common
divisor at least two would bound their total length by \(N\).
The weighted tree criterion realizes the partitions and
gives a single pole of degree \(N\) over the third value.
Place that pole at infinity to obtain the polynomial \(S\).
The preceding tree correspondence is applied with the two
finite vertex values \(0,1\).

Mark as \(u=0\) the reserved point of degree \(\mu\) over
\(1\) if \(\mu\) is integral, or the prescribed point of
degree \(b\mu\) over \(0\) otherwise. Assign the roots
\(\lambda_i\) to the remaining reserved points over \(1\)
for the integral \(\alpha_i\), and to their prescribed
points over \(0\) for the nonintegral \(\alpha_i\).
These assignments use distinct points even when several
degrees coincide. Put
\[
                         H=\frac{bS}{uS'}.
\]
At an integral nonzero root, \(H\) has a pole of order
\(\alpha_i\). At a nonintegral root or an extra zero of
\(S\), it has a simple zero. At \(u=0\), \(H\) has a
pole of order \(\mu\) if \(\mu\) is integral and is a
unit otherwise. It is a unit at infinity. Consequently
\[
                S \ \propto\ u^{b\mu}
                           \prod_i(u-\lambda_i)^{b\alpha_i} H^b
\]
by comparison of divisors. Logarithmic differentiation,
using \(S'/S=b/(uH)\), gives the required equation
\[
        uH'+\left(\mu+\sum_i\frac{\alpha_i u}{u-\lambda_i}\right)H=1 .
\]
The source normalizations and assignments can be made over
\(\overline{\mathbb Q}\). They work over every algebraically
closed characteristic-zero field.

For instance, when \(\mu=2\) and \(\alpha_1=\alpha_2=1/2\)
the zero partition is \((3,3)\) and the reserved degree
over \(1\) is two. One may use \(S=(1-u^2)^3\) and
\(H=(u^2-1)/(3u^2)\), with roots at \(1,-1\).

## The order-twelve profile

**Corollary OrderTwelveLeadingFactorTreeClassification.**
In the order-\((12,42)\) Weyl reduction let \(H,L\) be the
numbers of multiplicity-fifteen and multiplicity-three
roots of \(h\). All profiles with \(H\ge1\), \(L\ge0\),
\(H\equiv L\bmod 2\), \(5H+L\ge8\), and \(L\le H+2\)
are realized by the leading-factor residue equation.
For each profile all complex configurations, up to affine
change, correspond to plane trees with trivalent internal
vertices and with \(L\) ordinary leaves and
\[
                              F=1+(H-L)/2
\]
leaves of a second color. This realizes only the
leading-factor exactness condition.

The [order-twelve reduction](POLYNOMIAL_WEYL_ORDER_TWELVE_RESIDUAL.md#the-remaining-differential-order-pair)
provides the multiplicities and the degree bound for the Weyl
application. In the residue equation, \(g=6\),
\(M=(5H+L)/2\) and \(\Delta=3\). When \(L>0\) the list
contains an entry three. When \(L=0\), the \(M-H>0\)
additional entries are six. The two partitions are \(3^H\)
for the poles
and \(1^L2^{(3H-L)/2}\) for the zeros, so \(N=3H\)
and \(N\ge M-1\) is exactly the last inequality.

In the unweighted bipartite dessin every pole vertex has
valence three and every zero vertex has valence one or two.
Suppress each bivalent vertex. The degree-one faces then
correspond to terminal loops at distinct pole vertices.
Their two original edges connect one bivalent vertex to
one trivalent vertex and use two of its three incidences.
Deleting these \(F\) loops leaves a connected graph with no
cycles, by Euler's formula. The pole vertices that lost a
loop are now the second color of leaf. The first color is
the simple-zero vertices, and the remaining pole vertices
have valence three. Conversely attach one loop at each leaf
of the second color and subdivide every edge between pole
vertices, including the loops, by a bivalent vertex.
The cyclic order of the loop edges has no extra choice up to
isomorphism since they bound a degree-one face. This gives
the inverse construction. When \(L=0\), this is the
cubes-and-squares tree construction of
[Pakovich and Zvonkin, Example 2.13](https://arxiv.org/pdf/1306.4141v1).

The same description applies locally to the known
[multiplicity-fifteen crossing](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md#the-two-roots-and-real-coefficients).
Its first symbol at weights \((-1,3)\) is \(y^2L(x^3y)^2\),
with \(L\) having nonzero roots of multiplicities four and
one. In \(t=(x^3y)^{-1}\), the mixed-face data are
\(\nu=1/3\), \(m=3\), and root multiplicities eight and two.
The signed entries are \(-5,1,3,1\), with the last entry
marked. The weighted star has one edge of weight three
and two of weight one. Marking one weight-one leaf gives
two cyclic orders up to marked isomorphism, exchanged by
reflection. This accounts for the two conjugate local
configurations found by the explicit quadratic calculation.

In particular, the smallest profile without a small root
has \(H=2,L=0\). Its single tree is an edge between two
leaves of the second color. One representative is
\[
            h=(x^2+15)^{15},\qquad B=x(x^2+15)(2x^2/675+1/15),
\]
for which direct differentiation gives \(B'-(h'/6h)B=1\).
The resultant's small-root conditions are vacuous in this
profile. Exclusion of an actual Weyl completion therefore
requires information beyond this residue equation.
