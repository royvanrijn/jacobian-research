# Exact total-nonnegativity order of the fixed-shape auxiliary triangle

Let \(d\ge1\) be an integer and \(c\ge0\) be real. Define the auxiliary
coefficient triangle \(R(c,d)\) by
\[
B_0=1,\qquad B_{n+1}=(z+dn+c)B_n-zB_n',\qquad
R_{n,k}(c,d)=[z^k]B_n .
\]
In the factorial adjoint, \(c=\gamma+d-1\), with the shape \(\gamma\)
held fixed as the row varies. Write \(\mathrm{TN}_r\) for
nonnegativity of all minors of sizes at most \(r\), with row and
column indices in increasing order.

**Theorem.** If \(c\) is a nonnegative integer, \(R(c,d)\) is
\(\mathrm{TN}_r\) for every \(r\). If \(c\) is noninteger, the largest
such \(r\) is \(\lfloor c\rfloor+1\). For every \(r\ge1\), each
order-\(r\) minor is a polynomial with nonnegative rational
coefficients in \(x\) after substituting \(c=r-1+x\).
These statements allow row and column indices arbitrarily far beyond
the positivity order. For \(c>r-1\), every order-\(r\) minor with
ordered row and column lists satisfying \(i_a\ge j_a\) is strictly
positive. The earlier fixed-shape production argument
gave total nonnegativity of the leading block through index \(N\)
under \(c\ge N\). The exact criterion for that block, when
\(N\ge1\), is that \(c\) is an integer or \(c\ge N-1\).
[Factorial adjoint](HIGHER_GAP_ADJOINT_RANK.tex) ·
[Leading-block positivity](gap_three_hurwitz_structure.tex)

## Separate the gap from the shape

Put \(h(t)=1-(1-dt)^{1/d}\). The recurrence gives
\[
\sum_{n\ge0}B_n(z)\frac{t^n}{n!}
    =(1-dt)^{-c/d}e^{zh(t)}.
\]
For a formal series \(g\) with \(g(0)=1\) and \(f(t)=t+O(t^2)\),
write \(\mathcal R[g,f]_{n,k}=(n!/k!)[t^n]g(t)f(t)^k\).
Multiplication of column generating functions gives
\[
R(c,d)=\mathcal R[1,h]\,\mathcal R[(1-u)^{-c},u].                 \tag{1}
\]
The first factor depends only on the gap. The second factor is the
lower triangular Toeplitz matrix with coefficients
\[
H_j(c)=[u^j](1-u)^{-c},\qquad H_j(c)=0\quad(j<0).
\]
Its rows and columns are scaled by \(n!\) and \(1/k!\), respectively.

Here is a production-matrix proof that the first factor in (1) is
totally nonnegative. Since \(h'(t)=(1-h(t))^{1-d}\), differentiation
of its \(k\)-th column generating function gives production entries
\[
P_{i,k}=\frac{i!}{(k-1)!}[u^{\,i-k+1}](1-u)^{1-d}
\quad(k\ge1),\qquad P_{i,0}=0 .
\]
Equivalently, conjugate by \(\operatorname{diag}(0!,1!,\ldots)\)
the product of the lower triangular Toeplitz matrix of
\((1-u)^{1-d}\) and the upper shift with entry \(k\) at
\((k-1,k)\). The Toeplitz matrix is a product of \(d-1\) lower
triangular matrices of ones. Its finite leading sections are
products of nonnegative lower bidiagonal matrices. The shift
is also totally nonnegative, so Cauchy–Binet gives the same
property for \(P\). To pass to the coefficient triangle, its
leading section through \(N\) is the product of
\(\operatorname{diag}(1,\mathcal R[1,h]_{0..N-1,0..N-1})\) with
the matrix obtained by placing \((1,0,\ldots,0)\) above the
production block \(P_{0..N-1,0..N}\). Minors of the latter
matrix either lie in the production block, vanish, or expand
positively along the new first row. Induction proves the assertion.
[Production-matrix formula, Sokal, Theorem 2.19](https://arxiv.org/pdf/2105.05583)

## Shape positivity and sharpness

For row indices \(i_1<\cdots<i_r\) and column indices
\(j_1<\cdots<j_r\), the minor of the Toeplitz factor vanishes if
\(i_a<j_a\) for some \(a\). Otherwise reverse rows and columns
and use the partitions
\(\lambda_a=i_{r+1-a}+a\) and \(\mu_a=j_{r+1-a}+a\). The minor is
\[
\det[H_{\lambda_a-\mu_b-a+b}(c)]_{a,b=1}^r
     =s_{\lambda/\mu}(1^c).
\]
For positive integer \(c\), this is the skew Jacobi–Trudi formula
at \(c\) variables all equal to one. For real \(c\), the notation
denotes its polynomial continuation. The Littlewood–Richardson
expansion expresses it as a nonnegative integer combination of
\(s_\eta(1^c)\) with \(\ell(\eta)\le r\). Read each row of a
contributing tableau from right to left after the earlier rows.
The largest entry in row \(a\) is at most \(a\). A larger
entry would violate the lattice-word condition before the row
has supplied any smaller entries. The hook-content specialization is
\[
s_\eta(1^c)=
\prod_{(a,b)\in\eta}\frac{c+b-a}{h_\eta(a,b)},
\]
where \(h_\eta(a,b)\) is the hook length. The identity extends
polynomially from integer \(c\). If \(c=r-1+x\), every numerator
factor is a polynomial in \(x\) with nonnegative coefficients.
For nonnegative integer \(c\), skew tableau counting gives
nonnegativity directly at every size. Combining these properties
with (1) by finite Cauchy–Binet proves the positive assertions.
If \(c>r-1\) and \(i_a\ge j_a\), the skew shape admits a
semistandard tableau by filling each row with its row number,
so its Littlewood–Richardson expansion has a nonzero term.
All hook-content factors in that term are then strictly positive.
In Cauchy–Binet use the row index set itself as the intermediate
set. The first factor has determinant one on this principal
minor, while the second gives the strictly positive skew value.
[Jacobi–Trudi and Littlewood–Richardson formulas, Morales–Pak–Panova, §§10.1–10.2](https://www.math.ucla.edu/~pak/papers/excited.pdf) ·
[Hook-content formula, Krattenthaler, Theorem 1 at \(q\to1\)](https://arxiv.org/pdf/math/9807068)

For sharpness the inverse triangle has column generating series
\((1-u)^c\) in column zero. Indeed
the inverse of \(\mathcal R[g,h]\) is
\(\mathcal R[1/g(h^{-1}),h^{-1}]\) and \(g(t)=(1-h(t))^{-c}\).
Triangular cofactor inversion gives the exact obstruction
\[
\det R(c,d)[1,\ldots,N;\,0,\ldots,N-1]
   =(-1)^N R(c,d)^{-1}_{N,0}
   =N!\binom cN .
\]
The inverse entry is unchanged when computed in the leading
section through \(N\). At \(N=\lfloor c\rfloor+2\) for noninteger
\(c\), exactly one factor in the numerator of \(\binom cN\) is
negative. Thus the asserted positivity order is sharp for
every integer gap. The same minor gives the necessary condition
for the leading block whenever \(\lfloor c\rfloor+2\le N\).
If \(c\ge N-1\), all orders through \(N\) are nonnegative, and
the full block has determinant one.

At the physical parameter \(\gamma=n+b-1>0\), put
\(N=n+d-2\) and \(c=N+b\). For every \(d\ge2\), \(n\ge1\)
and \(b\ge0\), the leading block through \(N\) contains every
auxiliary row used in the cyclic Wronskian. Each minor of
that block with \(i_a\ge j_a\) is strictly positive. Indeed,
\(c>r-1\) for every \(r\le N\), and the only minor of order
\(N+1\) has determinant one. Its weak nonnegativity also
follows by specializing the prior leading-block theorem to
\(N=n+d-2\).
[Cyclic Wronskian rows and leading-block theorem](gap_three_hurwitz_structure.tex)

## An obstruction at the rotated adjacent rows

The gap-three adjoint uses the adjacent rows \(B,C\) of the
fixed-shape triangle through the rotated quotient
\[
W(z)=\frac{C(\alpha z)B(\overline\alpha z)
              -B(\alpha z)C(\overline\alpha z)}
              {(\alpha-\overline\alpha)z},
\qquad \alpha=e^{i\pi/3}.
\]
Ordered two-row minors control the same divided difference in
two independent real variables. Rotation introduces signed sine
weights, and total nonnegativity of the input triangle alone
does not give \(\mathrm{TN}_r\) of the resulting generalized
Hurwitz matrix. This persists when every coefficient of \(W\) is
strictly positive.

For an explicit control take the infinite Pascal triangle
\((\binom nk)_{n,k\ge0}\), with its \(k\)-th column scaled by
\(d_k\) and \(n\)-th row scaled by \(1/d_n\), where
\(d_4=7\) and all other \(d_j=1\).
Pascal's triangle is a planar binomial path matrix and is totally
nonnegative. The scalings preserve this property and give the
monic adjacent rows
\[
B=(1+z)^3,\qquad C=((1+z)^4+6z^4)/7 .
\]
With \(Q=1+z+z^2\), the quotient is
\[
W(z)=\frac{Q^3+6(-z^3+3z^5+z^6)}7
   =z^6+3z^5+\tfrac67z^4+\tfrac17z^3+
       \tfrac67z^2+\tfrac37z+\tfrac17.
\]
The expression follows by expanding the divided difference
of \(z^4(1+v)^3\) and using
\((\alpha^j-\overline\alpha^j)/(\alpha-\overline\alpha)
=\sin(j\pi/3)/\sin(\pi/3)\).
In the convention \(H_3(W)_{i,j}=[z^{6-(3j-i)}]W\), the minor
on rows \(2,3\) and columns \(1,2\) is \(3/7-6/7=-3/7\).
These two rows do not satisfy the factorial recurrence. The result
isolates the need for that recurrence in any transfer to the
actual cyclic Wronskian. Total nonnegativity of its generalized
Hurwitz matrix and unrestricted Strong Factorial remain open here.
[Cyclic Wronskian and factorial discriminant](gap_three_hurwitz_structure.tex)
