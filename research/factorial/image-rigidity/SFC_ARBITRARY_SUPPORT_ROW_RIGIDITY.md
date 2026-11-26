# Single-order factorial row rigidity for arbitrary supports

An exponent profile with \(s\) labeled columns is a finite
multiset \(A\) of nonnegative integer rows \(b\in\mathbb N^s\).
Its degree-\(m\) factorial moment form is
\[
       P_{A,m}(c)=
       \sum_{|\alpha|=m}\binom m\alpha c^\alpha
                            \prod_{b\in A}(b\cdot\alpha)! .
                                                               \tag{1}
\]
Appended rows that are constant across all columns multiply
each form by a scalar depending only on \(m\).

For a directed column pair \(i\to j\), define
\[
 G_{ij}^+(A)=\sum_{b\in A}\max(b_j-b_i,0),\qquad
 G_{ij}^-(A)=\sum_{b\in A}\max(b_i-b_j,0).
\]
For two profiles put
\[
 D_{ij}(A,A')=
 \max\{G_{ij}^+(A)+G_{ij}^-(A'),\
       G_{ij}^-(A)+G_{ij}^+(A')\}.
                                                               \tag{2}
\]

**Theorem.** Suppose \(P_{A,m}\) and \(P_{A',m}\) agree
after a nonzero overall rescaling and independent nonzero
rescalings of the coefficient variables. For each column
pair with \(m>D_{ij}(A,A')\), the full rows having
\(b_i\ne b_j\) agree, including every other coordinate
and all multiplicities. In particular, if the pairs
satisfying the inequality contain a spanning tree on
the column labels, the two multisets of nonconstant rows
are identical. Conversely identical nonconstant rows give
proportional moment forms at every order.

This gives a finite exact rigidity threshold for any
support size. It depends on gap sums along column pairs
and is independent of the common exponents of the rows.
In particular no distinct nonconstant row profile can be
obtained by an exact merger or splitting that preserves
all moment forms up to diagonal coefficient changes.
The theorem concerns equality of forms, allowing
arbitrary diagonal scalings at the tested order.

## Rational ratios on a lattice simplex

Write \(W_A(\alpha)=\prod_{b\in A}(b\cdot\alpha)!\) for
the coefficient in (1) with its multinomial factor
removed. For \(|t|=m-1\), the ratio
\(W_A(t+e_j)/W_A(t+e_i)\) is the evaluation of the
rational function
\[
   \mathcal R_{A,ij}(t)=
   \frac{\displaystyle
       \prod_{\substack{b\in A\\b_j>b_i}}
        \prod_{h=1}^{b_j-b_i}(b\cdot t+b_i+h)}
        {\displaystyle
       \prod_{\substack{b\in A\\b_i>b_j}}
        \prod_{h=1}^{b_i-b_j}(b\cdot t+b_j+h)} .
                                                               \tag{3}
\]
All its factors are positive on the nonnegative real
simplex. Equivalence of the forms implies that
\(\mathcal R_{A,ij}(t)/\mathcal R_{A',ij}(t)\) is the same
nonzero constant on every integer point of this
simplex. The constant is positive rational, as is seen
from any one point, even if the allowed coefficient
scales have been taken in a larger field.

Eliminate one coordinate using \(\sum t_\nu=m-1\).
The cross-multiplied polynomial has total degree at
most \(D_{ij}(A,A')\). A polynomial of degree at most
\(m-1\) in the remaining variables is determined by its
values at the nonnegative integer points of total
degree at most \(m-1\). Indeed expand it in the basis
\(\prod_\nu\binom{t_\nu}{\gamma_\nu}\). The value at
\(t=\gamma\) determines that coefficient after the
ones of smaller total degree, since every other basis
element of equal or greater degree vanishes there.
Consequently \(m>D_{ij}(A,A')\) upgrades the sample
equalities to a rational identity on the full affine
hyperplane \(\sum t_\nu=m-1\).

## Recovering complete rows from the factors

The derivative of a factor of (3) in the direction
\(e_j-e_i\) is nonzero, positive for a numerator factor
and negative for a denominator factor. Thus none is a
constant on the affine hyperplane. Two proportional
factors are proportional by a positive constant since
they are positive on the simplex. Their directional
derivatives must have the same sign. In particular
there is no common factor to cancel between the
numerator and denominator of (3).

Normalize each affine-linear factor by the absolute
value of that directional derivative. Unique
factorization now identifies the two multisets of
normalized numerator factors and the two multisets
of normalized denominator factors in the rational
identity. To see what one factor determines, let
\(u\) be the lower index of the row in the transfer
and put \(a=b_u,\ d=|b_j-b_i|>0\).
The normalized factor is
\[
                     \frac{b\cdot t+a+h}{d},
                         \qquad 1\le h\le d .
\]
Its derivatives in the directions \(e_\nu-e_u\)
recover the entire normalized vector
\(\mathbf v=(b-a\mathbf 1)/d\). Its value at
\((m-1)e_u\) is the slot \((ma+h)/d\).
The sign of the original derivative specifies
whether \(u=i\) or \(u=j\). Thus each normalized
factor records both the full relative row vector
and a slot of the same row's smaller exponent.

Fix a sign and a vector \(\mathbf v\), and form the
signed difference of row multiplicities between the
two profiles. If it does not vanish, choose the
largest gap \(d\) occurring with nonzero signed
multiplicity, and choose a smaller exponent \(a\)
with such a multiplicity. Among the \(d\) integers
\(ma+1,\ldots,ma+d\), at least one is coprime to
\(d\). Its slot has reduced denominator \(d\) and
cannot be contributed by a row of smaller gap.
Rows with the same gap but different smaller
exponents have disjoint slot intervals, because
\(m\ge d\), which follows from (2). Once the gap,
smaller exponent and relative vector are fixed, the
entire row is fixed as \(a\mathbf 1+d\mathbf v\).
The selected slot therefore gives a nonzero signed
factor count, contradicting the identity.

This recovers exactly the complete rows unequal
across \(i,j\). Any nonconstant row is unequal
across some edge of a spanning tree. Recovery along
the edges gives its multiplicity, proving the claim.
The same argument works with distinct or coincident
columns. When the columns are distinct, (1) is the
factorial moment of the polynomial with that support.

For example the three-column profiles with even and
odd parity rows in \(\{0,1\}^3\) have the same row
projections onto each pair of columns. At order three
their stripped coefficients at \(\alpha=(1,1,1)\)
are respectively \(8\) and \(6\), while their
coefficients on the two-column faces agree. Every
pair has \(D_{ij}=2\), so this is an instance where
the threshold detects the full rows but the edge
profiles alone do not.

## Why the gap bound matters

No order threshold independent of the row gaps
can distinguish arbitrary profiles. Fix \(s\ge2\)
and a finite bound \(B\). For each
\(\alpha\in\mathbb N^s\) with \(1\le|\alpha|\le B\),
choose distinct nonnegative integer vectors \(u_\alpha,
v_\alpha\) satisfying \(\alpha\cdot u_\alpha=
\alpha\cdot v_\alpha\), which is always possible
in at least two dimensions. The nonzero polynomial
\[
                  \prod_\alpha (X^{u_\alpha}-X^{v_\alpha})
\]
has signed coefficients describing a difference of
two row multisets. Under \(X_i\mapsto T^{\alpha_i}\)
it vanishes for each selected \(\alpha\), so the
histograms of \(b\cdot\alpha\) on the two multisets
are exactly equal. Therefore their raw forms (1)
agree at every order through \(B\). The difference
cannot be supported only on constant rows, since
projection to coordinate one would then distinguish
them. A common fixed core can be added to ensure
distinct columns, and a common translation of all
rows preserves the construction. This is the
finite-order switching construction in
[exact factorial row replacement](MULTIVARIATE_EXACT_FACTORIAL_ROW_REPLACEMENTS.md)
applied to the count directions of arbitrary support.

For two columns, the theorem specializes to the
[ordered slot rigidity criterion](MULTIVARIATE_SINGLE_ORDER_ROW_RIGIDITY.md).
In higher support, using only the two-column edge
restrictions of the moment form would recover pairwise
projections of the rows. The transfer ratios used here
vary the other column counts as well, and therefore
recover the full rows. The result imposes no
restriction on row transformations that preserve
consecutive coprimality or the
[strong factorial window](https://arxiv.org/abs/1304.3956)
without preserving the moment forms.
