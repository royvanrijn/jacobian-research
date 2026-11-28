# The six-cell $(2,2,1,1)$ interface gap

**Theorem (Quadratic six-cell gap).** For every fixed $0<\rho<1$ and $d\ge5$, the full isotypic block of
Young shape $(2,2,1,1)$ has the explicit gap
\[
 \mathfrak q_\rho(h)\ \ge\
 \frac{\ell}{\sqrt{1-\rho^2}}(1-\rho)\|h\|^2.
 \tag{QH1}
\]
Here $\ell=6\nu_{12}(F_{12})$ is the local coefficient scale in
(S77) of the [simplex interface analysis](SIMPLEX_INTERFACE_SYMMETRY.md).
The key is a positive quadratic field on the ordered facet.
Its direct conditional residual has negative coordinates at large
cutoff. Integration by parts gives a different conditional formula
for the full integral, with nonnegative prefix sums.

## A positive quadratic field in the Young cone

Index the untied coordinates by $1,2,3,4$, corresponding to labels
$3,4,5,6$. The $(12)$-minus Young space has six standard tableaux.
After the fixed first vertical pair, write the untied boxes as
$A=(1,2)$, $B=(2,2)$, $C=(3,1)$, $D=(4,1)$ and order the basis by
$ABCD,ACBD,ACDB,CABD,CADB,CDAB$. Use the rational positive
off-diagonal Young convention of (S74). The last vector $e_6$ has
sign on labels $1,2,3,4$ in the full representation and also sign on
the remaining labels $5,6$.

Let $R$ be this representation restricted to the four untied
labels and, for $i<j$ with complementary pair $a<b$, write
$\pi_{ij}=(i,j,a,b)$ in one-line notation. Define
\[
 W(m)=\sum_{i<j} R(\pi_{ij})e_6
                    (m_i-m_j)(m_a-m_b).
 \tag{QH2}
\]
An untied-label permutation changes the ordering of each pair by
the same signs in the vector and scalar factors. Thus (QH2) is
equivariant. It is translation invariant, homogeneous quadratic,
and multi-affine as a polynomial of the four original coordinates.
In particular, when the four coordinates are conditionally independent,
the field of their conditional means equals their conditional field
expectation.

Let $b=m_1-m_2$, $c=m_2-m_3$, $d=m_3-m_4$. The positive diagonal
matrix
\[
 D=\operatorname{diag}(5/9,10/27,10/9,10/9,5/12,5/3)
\]
puts the field in the form
\[
 D^{-1}W(m)=
 (bd,\ 2bc+3bd+2c^2+2cd,\ c(b+c+d),\ c(b+c+d),\
        4bc+3bd+4c^2+4cd,\ bd)^T.
 \tag{QH3}
\]
All entries are strictly positive on the ordered interior. The
field is unique up to scale among equivariant polynomial fields of
degree at most two. In (S90) the projection to the last basis
vector detects the field. Its scalar polynomial is alternating in
each pair $(1,2)$ and $(3,4)$ of untied labels, so it is divisible
by $(m_1-m_2)(m_3-m_4)$ and hence is determined up to scale.

Write $T_i=I-3A_i$ for the adjacent transport (S85) at untied index
$i$. In this basis $T_1=\operatorname{diag}(1,1,1,-2,-2,-2)$, since
the first three full labels lie in a column in the last three
tableaux. Then $T_i=R((1i))T_1R((1i))^{-1}$. We need the following
finite prefix comparison. For $t_0\ge t_1\ge t_2\ge t_3\ge t_4$
and $1\le r\le4$, let $t^{i\to t_0}$ replace only position $i$
by $t_0$. Then
\[
 P_r(t)=\sum_{i\le r}(I-T_i)W(t^{i\to t_0})\ \ge\ 0
       \quad\text{coordinatewise}.
 \tag{QH4}
\]
For reproducibility, put $a=t_0-t_1$, $b=t_1-t_2$, $c=t_2-t_3$,
$d=t_3-t_4$ and $s=a+b+c+d$. In the columns below the entries
are the six coordinates of $D^{-1}P_r$ in order. Set
$A=2bc+3bd+2c^2+2cd$, $B=4bc+3bd+4c^2+4cd$,
\[
 C_2=4ab+4ac+3ad+4b^2+8bc+4bd+4c^2+4cd,\qquad
 D_2=4ab+16ac+12ad+4b^2+20bc+13bd+16c^2+16cd,
\]
\[
 D_3=3ab+4ac+3b^2+19bc+12bd+16c^2+16cd.
\]
Then the four columns are
\[
\begin{array}{c|cccc}
 &1&2&3&4\\ \hline
1&0&0&4bs&4bd\\
2&0&8(b+c)s&4(3b+2c)s&4A\\
3&0&C_2&c(a+4b+4c+4d)&4c(b+c+d)\\
4&3cs&(b+4c)s&4cs&4c(b+c+d)\\
5&3(4ac+3ad+B)&D_2&D_3&4B\\
6&3d(a+b)&3bd&b(a+b+c+4d)&4bd
\end{array}
\]
These formulas follow by substituting (QH2) into the definition of
$P_r$ and applying the adjacent-generator matrices of (S74).
The exact script reproduces them and checks every monomial coefficient.
In particular, for $f_1\ge f_2\ge f_3\ge f_4\ge0$, Abel summation
of (QH4) gives
\[
 \sum_i f_i(I-T_i)W(t^{i\to t_0})
       =\sum_{r=1}^4(f_r-f_{r+1})P_r(t)\ge0,\qquad f_5=0.
 \tag{QH5}
\]

## The integrated Gaussian comparison

Apply the unit-variance scaling of (S91d) to the quadratic
field. Translation invariance allows the source means to have
tied mean zero and untied means $0>m_1>m_2>m_3>m_4$. Write $u$ for
the common cutoff, $F_i=\Phi(u-m_i)$,
$f_i=\phi(u-m_i)/F_i$, $f_0=\phi(u)/\Phi(u)$,
$g_0(u)=\pi^{-1/2}e^{-u^2}$, and
\[
 t_i=m_i-f_i,\qquad t_0=-f_0,\qquad
 G(u)=g_0(u)\prod_i F_i/f_0.
\]
The adjacent scalar density ratio $f_i/f_0$ is the same
completion of squares as in (S92a). Its conditional means
replace one coordinate by $t_0$ and leave the other three $t_j$
unchanged. After multiplication by the positive factor $f_0$,
the conditional unit-variance row comparison has residual
\[
 E(u)=\left(2f_0+\sum_i f_i\right)W(m)-2f_0 W(t)
                    -\sum_i f_i T_iW(t^{i\to t_0}).
 \tag{QH6}
\]
The full scaled comparison is therefore $\int G E\,du\ge0$
coordinatewise.

The key exact identities are $G'/G=f_0-u+\sum_i f_i$ and
$t_i'=f_i(u-t_i)$. They follow by differentiating $\Phi$ and
the upper truncated means. With $g=u+f_0$, integrate the
$G'$ contribution against $W(m)-W(t)$. The boundary terms
vanish by Gaussian decay of $g_0/f_0$ and polynomial growth of
the truncated means, exactly as in the
[cubic integrated comparison](NONHOOK_CUBIC_CONDITIONAL_OBSTRUCTION.md).
Since $W$ is affine in each single coordinate, this gives
\[
 \begin{split}
 \int G E\,du
 &=\int G\left[g(W(m)-W(t))
    +\sum_i f_i\bigl(W(t^{i\to u})-T_iW(t^{i\to t_0})\bigr)\right]du\\
 &=\int G\left[g\bigl(2\mathrm d W(t)[f]+W(f)\bigr)
    +\sum_i f_i(I-T_i)W(t^{i\to t_0})\right]du .
 \end{split}
 \tag{QH7}
\]
Here $\mathrm d W(t)[f]$ denotes the directional derivative
of the field. The second line uses $u-t_0=g$ and $m=t+f$.

The truncated mean map is strictly increasing in $m_i$, with
derivative equal to the conditional variance. The inverse
Mills value $f_i$ is also strictly increasing in $m_i$, since
its derivative is the Mills value times the positive residual
mean below $u$. Thus $t_0>t_1>\cdots>t_4$ and
$f_1>\cdots>f_4>0$. From (QH3) the directional derivative
$\mathrm d W(t)[f]$ is nonnegative, since the directions in the
three ordered gaps are all positive. The field $W(f)$ is
strictly positive in each coordinate. Also $g=u-t_0>0$.
The last sum in (QH7) is nonnegative by (QH5). The integrated
residual is therefore strictly positive in every coordinate.

For comparison, (QH6) itself fails coordinatewise at some
cutoffs. In the normalization $\widehat W=108W/5$, at
$m=(-1,-2,-3,-4)$ the corresponding residual satisfies
$(108/5)E(u)/f_1\to\widehat W(m)-T_1\widehat W(0,-2,-3,-4)$ as
$u\to+\infty$. Its value is
$(-12,-40,-24,264,531,180)^T$.
Indeed $f_i/f_1\to0$ for $i>1$, and the self residual divided
by $f_1$ is $O(f_0)$. The integrated identity (QH7) supplies the desired comparison.

## Operator consequence

The self and adjacent field integrals scale by $\sigma^2$ under
unit-variance normalization, where $\sigma=\sqrt{1-\rho^2}$.
The unit-variance source is $z=\rho x/\sigma$, and
$\sigma^2W(z)=\rho^2W(x)$. The row normalization in (S91d)
then converts (QH7) into
$\mathsf K_\lambda W(x)\le\rho a_\rho(x)W(x)$ on the ordered
facet. Positivity of the folded Young kernel (S69) makes
$\sqrt{a_\rho}W$ a positive supersolution of its normalized
self-adjoint compact operator at factor $\rho$. The strict integrated residual and the compact eigenvector
argument in (S91) give operator norm strictly below $\rho$ on
the full multiplicity space. The lower bound $a_\rho/\rho\ge\ell/\sigma$
and the unused-coordinate Hermite factors of (S94) give (QH1).


Run the [exact checker](check_six_cell_2211_quadratic_gap.py)
from the repository root to reconstruct the Young matrices
and verify all four columns of the prefix comparison.

~~~sh
.venv/bin/python research/plurality-stability/check_six_cell_2211_quadratic_gap.py
~~~

The [replay record](six_cell_2211_quadratic_gap_verification.json)
keeps the finite check and its scope. The Gaussian comparison
and the operator argument are written proofs.
