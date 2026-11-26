# A coefficient-uniform window for homogeneous binary factorial moments

**Theorem.** Fix \(d\ge1\). For every nonzero homogeneous
polynomial \(f(X,Y)\) of degree \(d\) over a
characteristic-zero field, every window of length
\[
                         2d^3+7d^2+8d
\]
at positive orders contains a nonzero factorial moment.
The bound is independent of the coefficients, including
their degree of algebraicity. If \(g(x)=f(x,1-x)\)
has degree zero or one, the respective lengths one and
two suffice. More generally the displayed bound can
use the actual positive degree of \(g\) in place of \(d\).

The argument bounds the size of a differential equation
for a normalized generating series and rules out every
nonzero polynomial solution of that equation. It allows
repeated roots of \(g\), including at integration endpoints.
The bound is in terms of degree and can exceed the
actual support size.

## Normalization and a rational differential system

Over \(\mathbb C\), integration of monomials against the
exponential product measure, followed by \(X=rx\),
\(Y=r(1-x)\), gives
\[
           m_n=\frac{\Lambda_2(f^n)}{(dn+1)!}
                  =\int_0^1 g(x)^n\,dx.
\]
The polynomial \(g\) is nonzero, since homogeneity
recovers \(f(X,Y)\) from \(g(X/(X+Y))\).
The case of constant \(g\) is immediate. If \(g(x)=a+bx\)
with \(b\ne0\), then
\(m_n=((a+b)^{n+1}-a^{n+1})/(b(n+1))\).
Two successive zero values would give two successive
equality-of-powers conditions for \(a+b\) and \(a\).
If one endpoint is zero the other is nonzero. If both
are nonzero the quotient would have to be one.
This proves the length-two claim.

Write \(e=\deg g\ge1\) for the rest of the proof. For
sufficiently small \(T\), set
\[
 F(T)=\sum_{n\ge0}m_nT^n=\int_0^1\frac{dx}{1-Tg(x)}.
\]
Put \(z=1/T\) and let \(r_1(z),\ldots,r_e(z)\) be the
local roots of \(g(r)=z\), away from the critical values.
For \(|z|\) large, choose the logarithm values so that
\[
 \ell_j(z)=\int_0^1\frac{dx}{x-r_j(z)}
              =\log(1-r_j(z))-\log(-r_j(z)).
\]
Partial fractions give
\[
                  F(1/z)=-\sum_{j=1}^e
                              \frac{z}{g'(r_j(z))}\,\ell_j(z). \tag{1}
\]

Let \(C(z)\) be the companion matrix in the convention
\[
        C(z)v_j=r_j v_j,\qquad
              v_j=(1,r_j,\ldots,r_j^{e-1})^{\mathsf T}.
\]
Its last row encodes reduction of \(r^e\) modulo
\(g(r)-z\). The first entry of that row depends affinely
on \(z\), and every other matrix entry is constant.
Define
\[
 D(z)=g'(C(z)),\qquad
 N(z)=D(z)^{-1} C(z)^{-1}(C(z)-I)^{-1}.
\]
These matrices are rational over \(\mathbb C(z)\).
Their inverses exist there because \(g(r)-z\) is
separable, and roots equal to \(0\) or \(1\) occur
only at the corresponding endpoint values.
If \(e_j^{\mathsf T}\) denotes the \(j\)-th coordinate
row, indexed from zero, let the \(k\)-th row of \(M(z)\)
be \(k e_{k-1}^{\mathsf T}D(z)^{-1}\) for \(k\ge1\),
and let its zeroth row be zero. Differentiation gives
\[
                 v_j'=M v_j,\qquad
                 \ell_j' v_j=N v_j ,
\]
using \(r_j'=1/g'(r_j)\) and
\(\ell_j'=r_j'/(r_j(r_j-1))\).

A rational system with state dimension \(2e\) is
\[
 U'=A(z)U,\qquad
 A(z)=\begin{pmatrix}M&N\\0&M\end{pmatrix},\qquad
 b(z)=\bigl(z e_0^{\mathsf T}D(z)^{-1},\,0\bigr).
\]
The \(2e\) solution columns
\[
                       \binom{v_j}{0},\qquad\binom{\ell_j v_j}{v_j}
\]
are a fundamental solution matrix, since the
Vandermonde matrix of the \(v_j\)'s is invertible.
The scalar outputs under the row \(b\) are
\(z/g'(r_j)\) and \(z\ell_j/g'(r_j)\).
In particular (1) is a scalar output.

Every scalar output tends to zero as \(z\to0\) on a
small sector. To see this, let \(r_*\) be the limiting
root of one inverse branch and \(u\ge1\) its multiplicity
as a zero of \(g\). Write
\(g(r)=(r-r_*)^u a(r)\) with \(a(r_*)\ne0\).
Local inversion in \(z^{1/u}\) yields
\(g'(r_j(z))\) of order \(z^{(u-1)/u}\).
Thus \(z/g'(r_j(z))=O(z^{1/u})\) in absolute value.
The two logarithms are bounded when their arguments
tend to nonzero values. If an argument tends to zero,
its logarithm has at most logarithmic growth on the
sector. This also covers \(r_*=0\) and \(r_*=1\).
Any continuation of the logarithms changes the
associated outputs only by constant multiples of
the first type of output.

## Scalar elimination and a uniform width

Set \(b_0=b\) and \(b_{j+1}=b_j'+b_j A\).
Let \(r\le2e\) be the first index for which the rows
\(b_0,\ldots,b_r\) are dependent over \(\mathbb C(z)\).
Then \(r\ge1\), and a rational relation with nonzero
last coefficient gives a scalar differential equation of
order \(r\) for all the outputs. Its solution space near
an ordinary point is exactly the output space.
Indeed the first \(r\) rows have rank \(r\) there, and
the fundamental solution matrix is invertible.
Thus the outputs give arbitrary initial jets of
orders zero through \(r-1\). The scalar existence
and uniqueness theorem gives the assertion.

In particular this equation has no nonzero solution
that is polynomial in \(T=1/z\). Such a solution, as
a function of \(z\), would be a scalar output near an
ordinary point. Continue the identity to a small
sector at \(z=0\), avoiding the finite algebraic
singularities along the path. All scalar outputs
decay there. A nonzero polynomial in \(1/z\)
cannot decay there.

For the degree bound, put
\[
                  B(z)=\det D(z)\det C(z)\det(C(z)-I).
\]
Every power \(C^j\) with \(0\le j\le e-1\) is affine
in \(z\). Indeed its rows come from reduction of
\(r^k\) with \(k\le2e-2\). In replacing a power \(r^e\)
using \(g(r)=z\), the resulting \(z\)-term has
\(r\)-degree at most \(e-2\) and needs no further reduction.
Thus \(D(z)\) is affine. The determinants and adjugates
of \(C\) and \(C-I\) are affine because their \(z\)-part
has rank one. It follows that \(B\), \(B A\) and \(B b\)
have polynomial entries of degrees at most \(e+2\).
For example the three adjugates in \(BN\) have total
degree at most \(e-1+1+1\), and the extra factor \(z\)
in \(Bb\) costs one degree.

Write \(q=e+2\). By induction \(b_j=L_j/B^{j+1}\)
with polynomial rows \(L_j\) of degrees at most
\((j+1)q\). The derivative identity giving the
induction has numerator
\[
            B L_j'-(j+1)B'L_j+L_j(BA).
\]
Take signed \(r\)-row minors of \(L_0,\ldots,L_r\)
in \(r\) columns where the first \(r\) rows are
independent. They give a relation
\(\sum_{j=0}^r c_j L_j=0\) with \(c_r\ne0\).
The relation lifts from these columns to the full rows
by the independence of the first \(r\) full rows
and uniqueness of the dependence. Thus the scalar
equation can be taken as
\[
               \sum_{j=0}^r c_j B^j\,\partial_z^{\,j} F(1/z)=0.
\]
For each \(j\),
\[
 \deg(c_j B^j)\le
             q\left(\frac{(r+1)(r+2)}2-1\right)=K.
\]
Indeed each minor omits the row of degree \((j+1)q\)
and the factor \(B^j\) adds at most \(jq\).

Now \(\partial_z=-T\theta\), where
\(\theta=T\partial_T\), and its \(j\)-th power is
\((-1)^jT^j\theta(\theta+1)\cdots(\theta+j-1)\).
Multiplying the converted equation by \(T^K\)
gives an annihilator of the form
\[
                     \mathscr L=\sum_{i=0}^W T^i P_i(\theta),
                     \qquad W=K+r .
\]
The conversion preserves the absence of nonzero
polynomial solutions. If \(m_n,\ldots,m_{n+W-1}\)
all vanished, let \(Q(T)=\sum_{j<n}m_jT^j\).
Then \(F-Q\) has order at least \(n+W\), and
\(\mathscr L\) does not lower that order.
But \(\mathscr L Q\) has degree at most \(n+W-1\).
The identity \(\mathscr L Q=-\mathscr L(F-Q)\)
forces \(\mathscr L Q=0\), contradicting \(Q(0)=1\).
Finally \(r\le2e\) gives
\(W\le(e+2)(2e^2+3e)+2e=2e^3+7e^2+8e\).

All identities concern the complex embedding of a field
generated by the finitely many coefficients of the phase.
Any finitely generated characteristic-zero field embeds
in \(\mathbb C\), so the same nonvanishing conclusion
holds over the original coefficient field. The factorial
normalizing factors are nonzero throughout.

## Why an abstract recurrence bound is insufficient

A coefficient-uniform conclusion uses more than
P-recursivity and eventual nonvanishing of each phase.
Even coefficient positivity and all-order positivity
on nonzero real pairs at even indices do not supply
the missing conclusion in the abstract recurrence setting.

For \(r\ge0\) set
\[
 V_r(X,Y)=\prod_{j=0}^{r-1}(X^2+jY^2)
              +\prod_{j=r}^{2r-1}(2X^2+jY^2),
\]
and define \(U_0=1\) and
\[
               U_{2r}=V_r+XYV_{r-1}\quad(r\ge1),
             \qquad
               U_{2r+1}=((r+1)X+Y)V_r .
\]
Each \(U_n\) is homogeneous of degree \(n\), with
strictly positive coefficients on all its degree-\(n\)
monomials. The second product in \(V_r\) gives the
positive coefficients on its even powers. The two
terms in \(U_n\) then fill the parity classes.
For nonzero real \(X,Y\) and \(r\ge2\),
\(V_r\ge(X^2+Y^2)V_{r-1}>0\). The first product
satisfies this by its appended factor. For the
second product the ratio is
\[
       \frac{(2X^2+(2r-2)Y^2)(2X^2+(2r-1)Y^2)}
                    {2X^2+(r-1)Y^2}
           \ge X^2+Y^2 .
\]
Since \(2|XY|\le X^2+Y^2\), every \(U_{2r}\) is
positive on nonzero real pairs. At \(r=1\) this is
the positive definite quadratic \(3X^2+2XY+Y^2\).

Both products in \(V_r\) satisfy first-order
rational recurrences in \(r\). Addition, shifts,
multiplication by \(r+1\), and interlacing the
even and odd subsequences preserve P-recursivity.
For interlacing, use the generating series with
\(T\) replaced by \(T^2\) and apply the chain rule.
Thus this is one symbolic P-recursive family.

For a fixed complex pair with \(Y\ne0\) put
\(x=X^2/Y^2\). The second product, scaled by
\(Y^{-2r}\), is \(\prod_{j=r}^{2r-1}(2x+j)\).
Its modulus is bounded below by a positive constant
times \((2r-1)!/(r-1)!\) for large \(r\), because
\(\sum_{j=r}^{2r-1}j^{-1}\) is bounded.
The first product is bounded above by a constant
times \((r-1)!r^{|x|}\).
Their quotient in modulus therefore tends to
infinity in favor of the second product.
Successive ratios satisfy \(V_r/V_{r-1}\sim4rY^2\),
by taking the corresponding ratio of the second products.
It follows that every fixed nonzero phase is eventually
nonvanishing in both subsequences. The case \(Y=0\)
follows directly from the products.

Nevertheless at \((X,Y)=(i\sqrt{k},1)\), \(k\ge1\),
the exact zero set of \(V_r\) is \(k+1\le r\le2k\).
Before this block both products have the same nonzero
sign, and after it the second product is nonzero.
At the stated phase the two terms of \(U_{2r}\)
are respectively real and purely imaginary. Thus
\(U_n\) has an exact zero block from \(2k+3\) to
\(4k+1\), of length \(2k-1\), with all other positive
indices nonzero. For every fixed odd prime \(p\),
choosing \(k=(p^s-1)/2\) puts that entire block
strictly between the native orders \(p^s-1\) and
\(p^{s+1}-1\). The abstract recurrence family
demonstrates the role of the polynomial-kernel
condition in the differential argument.
