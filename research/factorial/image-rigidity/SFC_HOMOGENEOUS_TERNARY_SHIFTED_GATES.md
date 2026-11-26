# Universal limit tests for shifted homogeneous ternary windows

Fix a strictly positive integer direction \(a\) and three distinct
ordinary exponent offsets \(b_1,b_2,b_3\) of the same total degree.
The ordinary exponents are \(Ta+b_j\). For each fixed window
\(m,m+1,m+2\), the following conditions are sufficient for its
regularity at every sufficiently large integer \(T\).

- If the offsets are noncollinear, require that the three forms
  \[
      [z^j]\exp\big(j(Az+Bz^2/2+Cz^3/6)\big),
                  \qquad j=m,m+1,m+2,
  \]
  have no common solution with \(C\ne0\).
- If the offsets are collinear, require that the three forms
  \[
      [z^j]\frac{\exp\big(j(Az+Uz^2/(2(1-2z)))\big)}
                     {\sqrt{1-2z}},
                  \qquad j=m,m+1,m+2,
  \]
  have no common solution in \(A,U\).

Both are ordinary finite polynomial tests with rational
coefficients. In the first test homogeneity permits \(C=1\).
The thresholds refer to the prescribed fixed window. In either
case, failure on an unbounded sequence of translations forces
a solution of the respective test. These implications leave
open the uniform family of tests as the start varies. For collinear
offsets, the [dyadic Gaussian comparison](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md)
excludes every even start and reduces odd starts to one valuation
pair.

## Cumulant normalization

Use the [tilted Gamma identity](SFC_TERNARY_COMMON_TRANSLATION_CUBIC_OBSTRUCTION.md#exact-tilted-gamma-normalization)
with \(t=1/T\), coefficients \(C_j\) normalized by \(\max_j|C_j|=1\),
and the common factor from the equal total degrees omitted.
At order \(m\), put \(Q=\sum_j C_j Y^{b_j}\), with
\(\delta_\rho=Y_\rho-1\). Let \(S,L,H\) be the value, gradient
and Hessian of \(Q\) at \(Y=\mathbf1\), and let
\(D_a=\operatorname{diag}(1/a_\rho)\), \(V=L^{\mathsf T}D_a L\).
The only order dependence of their expectations is through
\(\epsilon=t/m\). The cumulants of \(Q\) satisfy
\[
\begin{split}
 \kappa_1(Q)&=S+O(t),\
 \kappa_2(Q)&=t V/m+O(t^2),\
 \kappa_3(Q)&=t^2 K/m^2+O(t^3),\qquad
        K=2\sum_\rho L_\rho^3/a_\rho^2+3L^{\mathsf T}D_a HD_aL,\
 \kappa_q(Q)&=O(t^{q-1})\qquad(q\ge4).
\end{split}                                                   \tag{1}
\]
Only cumulants up to the largest fixed test order are used.

Here is the combinatorial reason for the general order. Expand
each copy of \(Q-S\) into monomials in \(\delta_\rho\) and group
the slots by the copies they came from. A cumulant of \(q\)
products is the sum of slot partitions that connect all \(q\)
groups, multiplying the elementary cumulants on each block.
To see connectivity, for an arbitrary slot partition group its
components on the \(q\) product labels. The moment-cumulant
formula on product labels then expresses the full moment sum
as the sum of products of the connected sums, proving the
assertion by its triangular recursion. Independence allows
a block only within one coordinate label. A block of size
\(h\ge2\) contributes \(O(t^{h-1})\) and a singleton contributes
\(O(t)\), by the tilted Gamma cumulants. Connecting \(q\)
groups costs at least \(q-1\) units of \(\sum_{\text{blocks}}(h-1)\).
At \(q=2\) the leading connected partition is one pair of
linear slots. At \(q=3\) the leading possibilities are one
triple of linear slots or two pairs with a central quadratic
factor. These give, respectively, the coordinate third
cumulants and the Hessian term in (1).

## Noncollinear offsets

Suppose the three moments vanish on an unbounded sequence and
pass to a projective coefficient limit. Scale by
\(R=\max(|S|,\sqrt{t|V|},t^{2/3})\).
If \(R/t^{2/3}\to\infty\) on a subsequence, (1) gives limits
of the cumulants of \(Q/R\) of the form \(A,B/m,0,\ldots\),
with \(A,B\) not both zero. In fact higher cumulants vanish
by \(t^{q-1}/R^q\to0\), and the first two errors are negligible.
If \(B=0\) the resulting moments are nonzero powers of \(A\).
If \(B\ne0\), the vanishing moments require the same number
\(A/\sqrt{-B}\) to be a normalized Hermite node at three
consecutive orders. At adjacent orders equality of two
nonzero nodes is impossible by
[scaled Hermite isolation](SFC_SCALED_HERMITE_ROOT_ISOLATION.md).
A zero node is impossible at the even order.
Thus \(R=O(t^{2/3})\). In particular \(S\to0\), \(V\to0\),
and the limiting gradient \(L\ne0\) by noncollinearity.

Take limits of \(S/t^{2/3}\) and \(V/t^{1/3}\), denoted
by \(A,B\), and normalize \(Q\) by \(t^{2/3}\).
Its cumulants tend to \(A,B/m,K_*/m^2,0,\ldots\), where
\(K_*\) is the limit of \(K\). At \(S=V=0\) its expression
reduces to
\[
           K_* = 3\sum_j C_j(b_j\cdot D_aL)^2
                          -\sum_\rho L_\rho^3/a_\rho^2 .
\]
This is nonzero by the
[circumradius bound](SFC_TERNARY_COMMON_TRANSLATION_CUBIC_OBSTRUCTION.md#independent-offsets-and-the-cubic-obstruction).
The generating series for these limiting moments is
\(\exp(Az+Bz^2/(2m)+K_*z^3/(6m^2))\). Rescaling \(z\) by \(m\)
gives the first universal test.

## Collinear offsets

Write \(b_j=b_0+\lambda_j u\) with distinct real \(\lambda_j\),
and \(J_1=\sum_j C_j\lambda_j\). Here \(L=S b_0+J_1u\).
Use \(R=\max(|S|,\sqrt t\,|J_1|,t)\).
If \(R/t\to\infty\), the Taylor expansion divided by \(R\)
has a leading scalar plus a Gaussian linear form, since all
terms of degree at least two in \(\delta\) are \(O(t/R)\)
in every fixed norm. The component \(S b_0\) in \(L\)
vanishes on this normalization. At least one of the leading
scalar and linear coefficients is nonzero by definition of
\(R\); the variance of \(u\cdot(\delta/\sqrt t)\) has
nonzero limit \(\sum_\rho u_\rho^2/(ma_\rho)\).
The same scaled-Hermite argument as above excludes this case.
Thus \(S=O(t)\), \(J_1=O(\sqrt t)\).

At the projective limit \(\sum_j C_j\lambda_j^2\ne0\) by
Vandermonde. The Hessian limit is
\((\sum_j C_j\lambda_j^2)u u^{\mathsf T}\), since the other
pieces of the Hessian expansion have factors \(S,J_1,L\).
Divide the Taylor phase by \(t\), and pass to limits of \(S/t\)
and \(J_1/\sqrt t\). The resulting moments are those of
\(A+B Z+Z^2\) with \(Z\) centered real Gaussian of variance
\(1/m\), after division by the nonzero quadratic coefficient
and a fixed rescaling of the real Gaussian linear form.
Its exponential moment series is
\[
  (1-2z/m)^{-1/2}
       \exp\!\left(Az + \frac{B^2 z^2}{2m(1-2z/m)}\right)
\]
by completing the square, coefficientwise. Rescaling \(z\)
by \(m\) gives the second universal test with \(U=B^2\).

## The window at orders 2, 3, 4

Both tests can be excluded in this entire window, so every
homogeneous ternary support on a positive common-translation
ray satisfies the window at orders \(2,3,4\) for sufficiently
large \(T\). For the cubic test, the moments at \(2,3\) give
\(B=-2A^2\) and \(C=9A^3\). Thus \(A\ne0\) and the fourth
moment is \(A^4\). For the quadratic test the second moment
gives \(U=-2A^2-2A-3/2\). After this substitution the third
and fourth moments, up to nonzero scalar factors, are
\[
  -18A^3-54A^2-45A-17,\qquad
  -320A^4-1280A^3-1680A^2-1056A-327 .
\]
They have leading coefficients prime to \(7\). The Euclidean
remainders modulo \(7\), starting with the quartic divided by
the cubic, are \(3A^2-2A-3,\ 1-2A,\ 2,\ 0\). The two polynomials
are therefore coprime over \(\mathbb Q\).
