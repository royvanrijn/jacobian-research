# Exclusion of polynomial Weyl pairs of differential orders four and six

**Theorem PolynomialFourSixStringExclusion.** Let \(k\) be a
characteristic-zero field and let \(P,Q\in k[x][D]\), where
\([D,x]=1\). If \(P\) and \(Q\) have differential orders four and
six, then \([P,Q]\) is not a nonzero scalar.

The proof uses the retained
[constant-leading case](MONIC_POLYNOMIAL_FOUR_SIX_STRING_EXCLUSION.md)
and
[linear leading-factor case](LOW_DEGREE_LEADING_FACTOR_POLYNOMIAL_FOUR_SIX_EXCLUSION.md).
For higher degree the
[residue constraint](POLYNOMIAL_FOUR_SIX_LEADING_FACTOR_RESIDUE.md)
supplies a root of multiplicity at least three, and the
[higher-root local gate](POLYNOMIAL_FOUR_SIX_HIGHER_ROOT_LOCAL_GATE.md)
controls the coefficients there.

Suppose \([P,Q]=\kappa\in k^\times\). Write their leading
coefficients as \(\alpha h^2,\beta h^3\). The two initial cases
allow us to assume \(\deg h\ge2\). Extend scalars to the
algebraic closure and choose a root \(b\) of \(h\) of multiplicity
\(m\ge3\). Put \(\tau=x-b\), and write \(p_i,q_i\) for the
normal-ordered coefficients in \(D\). By the local gate,
\[
 \tau^i\mid p_i,q_i\quad (i\ge2),\qquad
 p_1(b)q_0'(b)-q_1(b)p_0'(b)=\kappa.
 \tag{1}
\]

## The weighted orders near the root

For a real parameter \(0\le s\le1\), give \(\tau^j D^i\) weight
\(i-sj\), and denote the resulting maximal weights of \(P,Q\)
by \(A(s),B(s)\), respectively. Both functions are continuous
since they are maxima over finite sets of affine functions.

For \(0<s<1\) a normal-ordering contraction loses exactly
\(1-s>0\), since it removes one \(D\) and one \(\tau\). The
commutator has weight at most \(A(s)+B(s)-(1-s)\). Its symbol
at that bound is the Poisson bracket of the leading symbols
with \(\{\xi,\tau\}=1\). Because \([P,Q]=\kappa\), one always has
\[
                         A(s)+B(s)\ge1-s.                 \tag{2}
\]
For all sufficiently small positive \(s\), the leading
monomials are proportional to \(\tau^{2m}\xi^4\) and
\(\tau^{3m}\xi^6\). These are the least \(\tau\)-exponents
in the coefficients of highest differential order. Hence
\(A(s),B(s)>0\) and \(3A(s)=2B(s)\) near zero.

Except at finitely many \(s\), each leading symbol is a
single monomial. At a parameter of this kind where both maximal
weights are positive, let \((i,j)\) and \((i',j')\) be the
exponents of \((\xi,\tau)\) in the leading monomials. Their
Poisson bracket, if nonzero, is a nonzero multiple of
\(\xi^{i+i'-1}\tau^{j+j'-1}\). It would have to be the leading
symbol of the nonzero scalar commutator. Thus its two
exponents would be zero, requiring the two original exponent
vectors to sum to \((1,1)\). The possibilities either include
\((0,0)\), of weight zero, or include \((0,1)\), of negative
weight. This contradicts the positivity of both maximal
weights. Consequently the bracket vanishes and the two
exponent vectors are proportional.

On each interval of unique leading monomials within the
region of positivity, the proportionality makes \(A(s)/B(s)\)
constant. Continuity across the finitely many transitions
then gives \(3A(s)=2B(s)\) on the full connected interval
of positivity starting near zero. If this interval ended at
some \(s_0<1\), both weights would tend to zero there by the
same equality and continuity. This contradicts (2).
The equality \(3A(s)=2B(s)\) therefore persists for
\(0<s<1\) and at \(s=1\).

At \(s=1\), (1) bounds both integral weights by one.
Terms of differential order at least two have weight at most
zero, as do terms of order zero. At least one of the first
order coefficients is nonzero at \(b\), by the second part
of (1), so at least one of \(A(1),B(1)\) equals one. This is
impossible for integral weights satisfying \(3A(1)=2B(1)\).

## Other differential orders

The [gcd-two theorem](POLYNOMIAL_GCD_TWO_AND_REDUCED_ORDER_BOUNDS.md)
extends the exclusion to every nondividing pair of polynomial
differential orders with gcd two.
