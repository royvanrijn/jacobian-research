# The cubic conditional comparison fails in the six-cell \((2,2,2)\) block

Let \(V=V_{(2,2,2)}^-\) be the \((12)\)-minus space of the
six-cell interface. The simplest extension of the exterior-hook
construction to this shape chooses a homogeneous cubic equivariant
field on the four untied coordinates and asks for the analogue of
(S91f) at each fixed common cutoff in (S92a). That conditional
inequality fails in the positive Young basis of (S69).

Use the notation of [the compressed interface and Gaussian
reductions](SIMPLEX_INTERFACE_SYMMETRY.md).
The four untied means are denoted by \(m=(m_1,m_2,m_3,m_4)\),
with \(0>m_1>m_2>m_3>m_4\) after translating the tied mean to
zero. The remaining-label representation on \(V\) is the
\((2,1,1)\) of \(S_4\), which is the exterior square of its
standard representation. Let
\[
 H=\begin{pmatrix}1&-1&0&0\\1&1&-2&0\\1&1&1&-3\end{pmatrix},
 \qquad W_\wedge(m)=-Hm\wedge H(m^{\odot2})
\]
with exterior coordinates \(12,13,23\). The field is unchanged
by adding a common constant to the untied means.

In the rational positive off-diagonal Young convention of (S74),
order the minus tableaux by their successive untied boxes
\(ABCD,ACBD,CABD\), where \(A=(1,2)\), \(B=(2,2)\), \(C=(3,1)\)
and \(D=(3,2)\). In this basis the field is
\[
 W(m)=X W_\wedge(m),\qquad
 X=\begin{pmatrix}2&1&0\\-4/3&4/3&2/3\\1&-1&1\end{pmatrix}.
\]
Indeed the adjacent transpositions of the four untied labels
have Young matrices
\[
 \begin{pmatrix}-1&0&0\\0&-1/3&8/9\\0&1&1/3\end{pmatrix},
 \quad
 \begin{pmatrix}-1/2&3/4&0\\1&1/2&0\\0&0&-1\end{pmatrix},
 \quad
 \operatorname{diag}(1,-1,-1).
\]
For each corresponding permutation matrix \(P_i\) on the four
untied coordinates, direct multiplication gives
\(X\bigwedge^2(HP_iH^\top(HH^\top)^{-1})=R(P_i)X\).
These rational matrices are a positive diagonal rescaling of
the orthonormal Young basis. For \(b=m_1-m_2>0\),
\(c=m_2-m_3>0\) and \(d=m_3-m_4>0\), the three entries of \(W\) are
\[
 \begin{split}
 &3b(2bc+bd+2c^2+2cd+d^2),\\
 &2d(3b^2+6bc+3bd+2c^2+2cd),\\
 &6cd(c+d).
 \end{split}
\]
In particular the field is strictly positive in the interior
of the ordered chamber. Up to scale this is the only nonzero
equivariant polynomial field of degree at most three. The
first three-label sign projection \(A_3\) is one-dimensional,
and on it the group of the three other untied labels acts by
the sign, by vertical three-strip branching for \((2,2,2)\).
The projected polynomial must therefore be a multiple of the
Vandermonde of those three untied means. Equivariance determines
the other three projections, which detect the field by (S89).

For cutoff \(u\) write \(F_i=\Phi(u-m_i)\),
\(f_i=\phi(u-m_i)/F_i\), and \(f_0=\phi(u)/\Phi(u)\). Let \(B_i\)
be the self conditional first and second moments of the \(i\)-th
untied normal,
\[
 B_i=(m_i-f_i,\ m_i^2+1-(u+m_i)f_i).
\]
Use \(B^{(l)}\) for the analogous four pairs with \(m_l\)
replaced by zero. For an array \(B\) of these pairs, set
\(\mathcal W(B)=-X(HB_{\cdot1}\wedge HB_{\cdot2})\).
Conditional independence gives the exact field expectation
\(\mathcal W(B)\) because the terms with a repeated index in the
wedge vanish.

After division by the common positive self density and
multiplication by \(f_0\), the residual of the proposed
conditional comparison is
\[
 E_u(m)=(2f_0+\textstyle\sum_l f_l)W(m)
       -2f_0\mathcal W(B)-\sum_l f_l T_l\mathcal W(B^{(l)}).
\]
The adjacent scalar ratio \(f_l/f_0\) is the same Gaussian
identity as in (S92a), independent of the choice of shape.
Here \(T_1=\operatorname{diag}(1,1,-2)\) by (S85), since only
the tableau \(CABD\) has the first three labels in a column.

At the rational point \(m=(-1,-2,-3,-4)\), take \(u\to+\infty\).
The ratios \(f_l/f_1\) tend to zero for \(l>1\), since
\(\phi(u-m_l)/\phi(u-m_1)=
\exp(u(m_l-m_1)+(m_1^2-m_l^2)/2)\) and \(F_l,F_1\to1\).
The self field satisfies
\(\mathcal W(B)-W(m)=O((1+u)f_1)\) by the two displayed
moment formulas and \(H\mathbf1=0\). Thus the self residual
divided by \(f_1\) tends to zero because \((1+u)f_0\to0\).
Every adjacent field tends to the polynomial field with the
appropriate replaced mean. Consequently
\[
 \lim_{u\to+\infty}\frac{E_u(m)}{f_1}
   =W(-1,-2,-3,-4)-T_1W(0,-2,-3,-4)
   = (24,32,12)-(66,68,-24)=(-42,-36,36).
\]
The first two coordinates are strictly negative for sufficiently
large cutoffs. This breaks the conditional analogue of (S91f)
for the unique positive cubic field. The [integrated cubic proof](SIX_CELL_222_GAP.md) gives the
full block gap using the regrouping below.

## An equivalent integrated comparison

Write \(P(u)=\prod_i F_i\), let \(g_0\) be the tied-coordinate
\(N(0,1/2)\) density, and put \(G=g_0P/f_0\).
The integrated unit-variance comparison for this field has
residual \(\int_{\mathbb R}G E_u\,du\). An integration by parts
gives the identity
\[
 \begin{split}
 \int_{\mathbb R}G E_u\,du
 =\int_{\mathbb R}G\bigg[
 &(u+f_0)\bigl(W(m)-\mathcal W(B)\bigr)\\
 &+\sum_l f_l\bigl(\mathcal W(B^{l\to u})
                  -T_l\mathcal W(B^{(l)})\bigr)\bigg]du ,
 \end{split}
\]
where \(B^{l\to u}\) is \(B\) with pair \(l\) replaced by \((u,u^2)\).
Indeed \(G'/G=f_0-u+\sum_l f_l\) and
\(B'_{lr}=f_l(u^r-B_{lr})\) for \(r=1,2\). The field expectation
is affine in the pair at each fixed coordinate, since wedge
terms with repeated indices vanish. Therefore
\(\frac{d}{du}\mathcal W(B)=
\sum_lf_l(\mathcal W(B^{l\to u})-\mathcal W(B))\).
Substituting these identities and integrating by parts proves
the rewrite. The boundary term vanishes by Gaussian decay
of \(g_0/f_0\) and polynomial bounds for the truncated moments
at both endpoints. The [integrated proof](SIX_CELL_222_GAP.md)
shows this regrouped integral is nonnegative coordinatewise.
