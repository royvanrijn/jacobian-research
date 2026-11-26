# An eigenvalue obstruction at split real-cubic primes

In the [independent three-pair case](CM_NEWTON_THREE_BALANCED_GLOBAL.md),
one simultaneous ordinary-reduction test would be a bounded eigenvalue
expression whose valuation is nonnegative at every place of every
one-failing-pair Frobenius field. The [trace polynomial obstruction](CM_NEWTON_SPLIT_POLYNOMIAL_OBSTRUCTION.md)
covers fixed polynomial combinations of relative traces. Here the
obstruction also covers Laurent expressions in individually labelled
eigenvalues, even without symmetry under exchanging roots.

Let \(F=\mathbb Q(\zeta_7+\zeta_7^{-1})\), and let \(p\) range over the
rational primes split in \(F(i)\) with any fixed auxiliary congruence
\(p\equiv1\pmod{\ell^{2N}}\) as in the trace polynomial obstruction.
Use its finite-field models: a trace \(a\in O_F\) is a unit at two primes
over \(p\) and has valuation one at the third, and all its real conjugates
have magnitude less than \(2\sqrt p\). Put \(\varpi=\sqrt p\) and label the
roots of the three conjugate polynomials as \(\alpha_i,\beta_i\), with
\(\alpha_i\beta_i=p\). Write \(t_i=\alpha_i/\varpi\). Normalizing any fixed
Laurent monomial in the six roots and \(\sqrt p\) by a fixed power of
\(\sqrt p\) gives a Laurent monomial in \(\varpi,t_1,t_2,t_3\), since
\(\beta_i=\varpi/t_i\).

**Proposition.** Fix a nonconstant Laurent polynomial
\[
P(\varpi,t_1,t_2,t_3)=\sum_{r\in\mathbb Z}\varpi^r P_r(t_1,t_2,t_3),
\qquad P_r\in\mathbb Q[t_1^{\pm1},t_2^{\pm1},t_3^{\pm1}],
\]
with finite support in \(r\). Over the one-failing-pair finite-field models
above, either its archimedean absolute values are unbounded or its value
has negative valuation at a place over \(p\) in models at arbitrarily
large \(p\). In the latter case one may keep a prescribed auxiliary
torsion level rational. This concerns expressions evaluated on varying
finite-field varieties, rather than a density assertion about a fixed
variety over a number field.

## Orders at the bad and ordinary pairs

Normalize the local valuation by \(v_p(p)=1\). If label \(i\) is bad,
the two eigenvalues there have slope \(1/2\), and \(t_i\) is a unit with
\[
 t_i+t_i^{-1}=\sigma_i(a)/\varpi.
\]
Its residue is one of the roots of \(X^2+1\). Let \(\zeta\) denote the
corresponding root in \(\mathbb Q_p\), available since \(p\) splits in
\(F(i)\). Then
\[
 t_i+t_i^{-1}=(t_i-\zeta)(1-(t_i\zeta)^{-1}),\qquad
 \frac{t_i-\zeta}{\varpi}\equiv\frac{u_i}{2}\pmod{\mathfrak m},
\]
where \(u_i\) is the residue of \(\sigma_i(a)/p\) and \(\mathfrak m\) is
the maximal ideal of the local eigenvalue field. For an ordinary label
\(j\), put \(s_j=2v_p(t_j)\in\{1,-1\}\). Then
\(t_j/\varpi^{s_j}\) has residue \(u_j^{-s_j}\), where \(u_j\) is the
unit trace residue. This follows by taking the unit root of the ordinary
quadratic to reduce to the trace, and using the root product \(p\).

For any one model and fixed place of \(F\) over \(p\), these three root
orderings can be exchanged independently by conjugating the place in
the eigenvalue field over \(F(\varpi)\). Indeed the three discriminants
\(\Delta_j=\sigma_j(a)^2-4p\) have valuation one at different primes of
\(F\) above \(p\), and zero at the other two primes. A square-class
relation in \(F(\varpi)\) would make a product of some \(\Delta_j\) equal
to a square in \(F\) or to \(p\) times a square. The valuations allow only
the empty product in the first case and the product of all three in the
second. The latter product is totally negative and \(p\) times a square
is totally positive. Thus the three square classes remain independent.
This supplies all signs \(s_j\) at the ordinary labels and either residue
root at the bad label. Varying the place of the cyclic cubic field
changes which label is bad, while permuting the three leading trace
residues.

For each bad label \(i\) write the other labels as \(j,k\), and expand
\[
 P_r(t_i,t_j,t_k)=\sum_{b,c}f_{r,b,c}(t_i)t_j^b t_k^c.
\]
For a nonzero coefficient let \(m_{r,b,c}\) be its order at \(\zeta\).
For primes outside a finite set, its first nonzero local coefficient
has unit residue. In the orientation \((i,\zeta,s_j,s_k)\), twice the
valuation of the associated term is
\[
 r+m_{r,b,c}+s_jb+s_kc.                                                \tag{1}
\]
After removing the minimum power of \(\varpi\), the leading sum has
distinct residue monomials of the form
\(u_i^{m_{r,b,c}}u_j^{-s_jb}u_k^{-s_kc}\). Indeed equal exponents force
equal \(b,c,m_{r,b,c}\), and the equality of (1) then forces equal \(r\).
The resulting nonzero Laurent residue polynomial has coefficients in
\(\mathbb Q(i)\) before reduction. After clearing coordinate powers and
multiplying by its complex conjugate, the bounded grid argument in the
[trace polynomial obstruction](CM_NEWTON_SPLIT_POLYNOMIAL_OBSTRUCTION.md#avoiding-fixed-cancellations-on-actual-frobenius-data)
can ensure it stays nonzero in actual models at arbitrarily large \(p\).
Hence each negative formal minimum in (1) is realized at a place in the
eigenvalue field of an actual model. The trace target can simultaneously
be any fixed interior point of the archimedean trace cube, with the
chosen auxiliary torsion level.

## The bounded and integral cases

For archimedean embeddings all \(|t_i|=1\). Arbitrary triples of points
of the unit circle away from \(\pm1\) can be approached by the actual
models: choose real trace targets \(t_i+t_i^{-1}\) and choose each root
ordering along its continuous branch. A nonzero Laurent polynomial is
nonzero somewhere on this product of circle arcs, by the univariate
identity theorem in each variable. If the highest supported \(r\) is
positive, choose such a target for \(P_r\). The corresponding evaluations
are unbounded. Thus in the archimedean bounded case every supported
\(r\) is nonpositive.

Suppose the expression is also integral at all places over \(p\) on all
the models. Equation (1), with both signs of each ordinary valuation,
forces
\[
 \operatorname{ord}_{\zeta_+}(f_{r,b,c}),\ \operatorname{ord}_{\zeta_-}(f_{r,b,c})
 \ \geq\ |b|+|c|-r.                                                    \tag{2}
\]
Here \(\zeta_+,\zeta_-\) are the two roots of \(X^2+1\).
For a fixed supported \(r\), let \(A_j\) be the maximum absolute exponent
in its \(j\)th variable. The power span in variable \(i\) is at most
\(2A_i\). Multiplying the coefficient in (2) by \(t_i^{A_i}\) produces
a polynomial of degree at most \(2A_i\), with both root multiplicities
unchanged. Thus \(A_i\geq |b|+|c|-r\) for every pair of exponents supported
on the other two labels. In particular \(A_i\geq A_j-r\) whenever
\(i\ne j\). Interchanging the indices rules out \(r<0\).

When \(r=0\), all three maxima are one common \(A\). If \(A>0\),
choose a term with \(|b|=A\) in one of the ordinary labels. The degree
bound at bad label \(i\) forces \(c=0\). The corresponding nonzero
coefficient must be proportional to
\(t_i^{-A}(t_i^2+1)^A\), because its polynomial of degree at most
\(2A\) has multiplicity \(A\) at each of the two roots. There is then a
term whose exponents in the first two labels both have magnitude \(A\).
Applying (2) with the third label bad would give \(A\geq2A\), a
contradiction. Therefore \(A=0\) and the only uniformly bounded and
integral expression is a rational constant.

The distinction from an exact multiplicative relation is important.
The same models have all three independent root swaps and maximal
expected Frobenius multiplicative rank.
[Frobenius ramification and independence](CM_NEWTON_FROBENIUS_RAMIFICATION_AND_DENOMINATORS.md#inertia-and-the-frobenius-weyl-group).

Noot's tensor-representation
argument uses integral sums of specially selected eigenvalue products
whose valuation depends on the tensor weight geometry. Its hypothesis is
the odd tensor power of the standard rank-two representations, whereas
here the Frobenius spectrum consists of a direct sum of rank-two pairs.
The proof above excludes the analogous bounded integral fixed Laurent
test in the individual roots of these pairs. It leaves open arguments
using global relations across residue primes or more general
denominators. [Noot, Theorem 2.2 and its proof, printed pp. 168–170](https://www.numdam.org/item/CM_1995__97_1-2_161_0.pdf).
