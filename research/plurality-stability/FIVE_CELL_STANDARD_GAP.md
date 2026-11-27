# The balanced five-cell interface gap through correlation 41/50

The standard sector of the five-cell Gaussian simplex interface has
gap at least \(\ell/1000\) at every \(0<\rho\le41/50\), where \(\ell\)
is the distance between two Gaussian cell-moment vectors. The other
five-cell sectors already have positive gaps by
[(S76)](SIMPLEX_INTERFACE_SYMMETRY.md#strict-supersolution-and-the-five-cell-reduction).
For every fixed correlation in the certified range and every ambient
dimension \(d\ge4\), (S42) gives a positive gap off geometric rotations
with a gap constant independent of \(d\). The
[reflection-chamber comparison](FIVE_CELL_A3_BRIDGE.md) covers
\(41/50\le\rho<1\), so the five-cell interface gap holds at
every positive correlation. The
[Gaussian Neumann limit](STANDARD_SECTOR_NEAR_ONE.md) gives the
exact first-order scale of the standard-sector gap near one.

## The three facet gaps

In the active space \(E=\mathbf 1^\perp\subset \mathbb R^5\), the
reference facet is \(F=\{x_1=x_2\ge x_3,x_4,x_5\}\). Its pair plane is
\(H=\{x\in E:x_1=x_2\}\). Write \(g_i=x_1-x_{i+2}\) for \(1\le i\le3\).
Then
\[
 x_1=x_2=\tfrac15(g_1+g_2+g_3),\qquad
 x_{i+2}=x_1-g_i.
\]
The corresponding three linear functionals on \(H\) are represented
by \((e_1+e_2)/2-e_{i+2}\). Under standard Gaussian measure on \(H\),
\(g\) has covariance \(C_0=I+\mathbf1\mathbf1^\top/2\). Thus
\(\mu=\nu_{12}/w\), with \(w=\nu_{12}(F)\), is the law of \(g\)
conditioned on its positive orthant. Write \(p\) for the probability
of that orthant. For a centered trivariate Gaussian of covariance
\(C\), sign symmetry gives
\[
 P_C=\frac18+\frac1{4\pi}\sum_{i<j}\arcsin\!\left(
             \frac{C_{ij}}{\sqrt{C_{ii}C_{jj}}}\right).
 \tag{FV1}
\]
In fact the expansion of the product of the three sign indicators
has no odd terms in expectation. For the two-sign terms, rotation
in the plane of the two normalized linear functionals gives
\(2\arcsin(C_{ij}/\sqrt{C_{ii}C_{jj}})/\pi\).
In particular \(p=1/8+3\arcsin(1/3)/(4\pi)\). The surface measure is
\(w=\phi_0p\) with \(\phi_0=(2\pi)^{-1/2}\), and (S77) gives
\(\ell=5w\).

Let \(\sigma_\rho=\sqrt{1-\rho^2}\). By (S49), the standard scalar
speeds are the permutation-invariant functions of \(g\), and balance
is exactly \(\int h\,d\mu=0\). On this balanced space put
\(D_\rho=\sigma_\rho a_\rho/(\rho\ell)\).
With the four-dimensional Mehler kernel \(M_\rho\), let \(A_\rho\)
have quadratic form
\[
 \frac1{5\rho}\iint h(x)h(y)
 \left[2M_\rho(x,y)+\sum_{l=3}^5 M_\rho(x,(2l)y)\right]d\mu(x)d\mu(y).
\]
The scaled standard form divided by \(w\ell\) is
\[
 \frac{\sigma_\rho q_{\rm std,\rho}(h)}{w\ell}
   =\int D_\rho h^2\,d\mu-\langle h,\sigma_\rho A_\rho h\rangle_\mu .
 \tag{FV2}
\]
Each coefficient \(C_m\) in \(A_\rho=\sum_{m\ge1}\rho^{m-1}C_m\)
is nonnegative by the ambient Hermite expansion (S11). The constant
coefficient vanishes by balance.

## A reciprocal bound with finitely many Gaussian tilts

Write \(s_i=\rho g_i/\sigma_\rho\) and \(q_0=\ell/(2\phi_0)=5p/2\).
In the cone integral for the score derivative from (B7), the normal
\(n=(e_1-e_2)/\sqrt2\) is positively tilted by \(z_n\). Equivalently
that normal coordinate has Rayleigh law \(N\), with density
\(n e^{-n^2/2}\) for \(n>0\), and is independent of the standard
Gaussian pair-plane coordinates. The three untied comparisons give
\[
 q_0D_\rho(g)=
 \Pr\{X_i+N/\sqrt2\ge -s_i\text{ for all }i\},
 \qquad X\sim N(0,C_0).
 \tag{FV3}
\]
In particular \(D_\rho\ge1\) by the value at \(s=0\), and \(D_\rho\)
is pointwise nondecreasing in \(\rho>0\). To bound the joint probability,
represent the \(X_i\) as \(S/\sqrt2+Z_i\) with independent standard
normals. Given \(S,N\), the three success events are independent with
probabilities that are nonnegative increasing functions of
\((S+N)/\sqrt2\). For two such functions the expected product is at
least the product of expectations, by the covariance formula using
two independent copies. Applying this twice gives the same product
bound for the three functions.

Each marginal failure probability \(p_{\rm fail}(s)\) is at most
\(Z(s)/4\), with
\[
 Z(s)=\tfrac13\big(e^{-s^2/3}+e^{-4s^2/9}+e^{-4s^2/3}\big).
\]
The one-dimensional marginal has exactly the variances used in
[(TW2)–(TW3)](TETRAHEDRAL_STANDARD_GAP.md#a-reciprocal-bound-from-the-two-facet-gaps).
In detail \(p_{\rm fail}(0)=1/4\) by the Rayleigh integration there.
The Gaussian upper-tail ratio after conditioning on \(N\) is bounded
by \(2\overline\Phi(\sqrt{2/3}\,s)\), whose polar integral is bounded by
\(Z(s)\) by an upper Riemann sum of three equal angular pieces.
Applying \((1-z/4)^{-1}\le1+z/3\) yields the reciprocal envelope
\[
 D_\rho^{-1}\le E_\rho
    =q_0\prod_{i=1}^3 (1+Z(s_i)/3).
 \tag{FV4}
\]
At the left endpoint \(a=0\), use \(E_0=1\), since the limiting
coefficient \(D_0\) is identically one.

## The planar Hermite tail and the finite coefficients

The argument for [(TW4)](TETRAHEDRAL_STANDARD_GAP.md#uniform-spectral-comparison-on-a-correlation-interval)
extends to the three-dimensional pair plane. At each ambient degree,
the absolute value of the quadratic term for a transported kernel is
at most the same-facet term by the orthogonal action on Hermite chaos.
The kernel weights add to five, cancelling the normalization of \(A_\rho\)
other than \(1/\rho\). On extension by zero from \(F\) to \(H\), the
same-facet bound is diagonal in the planar Hermite degree \(k\). Its
coefficient after ambient truncation at \(N\) and scaling by
\(\sigma_\rho/\rho\) is
\[
 \lambda_k(\rho)=\rho^{k-1}
 \left(1-\sigma_\rho\sum_{j<J}\frac{{2j\choose j}}{4^j}\rho^{2j}\right),
 \qquad J=\max\{0,\lfloor (N-k)/2\rfloor+1\}.
\]
The odd normal degrees vanish on \(H\). Reindex the even ones and
use the decreasing coefficients \({2j\choose j}/4^j\) to obtain
\(\lambda_k(\rho)\le \rho^N\). For \(k\ge1\), \(\lambda_k\) is
nondecreasing by the derivative calculation in (TW4).
The planar constant moment vanishes by ordinary balance.
The zero extension contributes \(p\) to the squared norm and \(p^2\)
to the product of two moment integrals. It follows on the slab
\([a,c]\) that the scalar tail cap is \(\tau=c^N/p\) and that one can
subtract from it the terms
\[
 \sum_{k=1}^N(c^N-\lambda_k(c))
 \sum_{|\beta|=k}\Big|\int h H_\beta^{\rm norm}\,d\mu\Big|^2 ,
 \tag{FV5}
\]
with \(H_\beta^{\rm norm}\) the orthonormal planar Hermite polynomials.
The \(p^2\) factors cancel in these subtracted terms.

Use the following orthogonal coordinate vectors in \(E\):
\[
 (3,3,-2,-2,-2),\quad (0,0,1,-1,0),\quad(0,0,1,1,-2),\quad
 (1,-1,0,0,0).
\]
Their squared norms are \((30,2,6,2)\).
The four raw coordinates on \(F\) are
\((2(g_1+g_2+g_3),-g_1+g_2,-g_1-g_2+2g_3,0)\).
After transport by \((23)\) they are
\((-3g_1+2g_2+2g_3,g_2,-g_2+2g_3,g_1)\).
Expansion in Hermite products with these variances gives the finite
kernel coefficients. For invariant \(h\), all three adjacent terms
have the same quadratic form as \((23)\), so it suffices to give that
transport weight three and the self term weight two. Average monomials
under permutation of the \(g_i\). Each monomial is replaced by the
normalized average of its orbit and its constant part can be
discarded. The denominator of a degree-\(m\) ambient term divides
\(5\cdot30^m m!\) after the \(1/5\) normalization. The denominator of
the planar term divides \(30^m m!\).

For \(m\le N\), take rational upper bounds on
\(\max_{\rho\in[a,c]}\sigma_\rho \rho^{m-1}\), calculated at the
endpoints and at the critical point \(\sqrt{(m-1)/m}\) if present.
Rounding the \(\lambda_k(c)\) upward gives rational lower bounds on
the planar subtractions. Combining this finite comparison and (FV5)
therefore gives
\(\sigma_\rho A_\rho\preceq K+\tau I\) for all \(\rho\in[a,c]\),
with \(K\) a finite self-adjoint symmetric polynomial kernel.

## Exact moments and the balanced weighted comparison

To reserve \(\gamma=1/1000\) in (FV2), let
\(e=E_a/(1-\tau-\gamma)\) provided the denominator is positive.
Monotonicity of \(D_\rho\) and \(D_a\ge1\) give
\(D_\rho-\tau-\gamma\ge e^{-1}\). The finite comparison is thus
sufficient if its norm on the balanced subspace of
\(L^2(e^{-1}\mu)\) is less than \(1\).
For the normalized orbit averages \(v_i\), the Gram matrix of their
mean-zero representing vectors is
\[
 G_{ij}=\int e v_i v_j\,d\mu-
        \frac{\int e v_i\,d\mu\,\int e v_j\,d\mu}{\int e\,d\mu}.
 \tag{FV6}
\]
Here the vectors represent \(h\mapsto\int h v_i\,d\mu\), so
ordinary balance is retained. If \(Q\) is the symmetric rational
coefficient matrix of \(K\), the nonzero eigenvalues of \(QG\) are
those of the finite weighted operator. For even \(j\), its operator
norm is bounded by \(\big[\operatorname{tr}((QG)^j)\big]^{1/j}\).

Each term in \(E_a\) is an exponential
\(\exp(-\sum_i r_i g_i^2)\) with
\(r_i\in\{0,a^2/(3\sigma_a^2),4a^2/(9\sigma_a^2),4a^2/(3\sigma_a^2)\}\).
Its nonzero factors in the product each carry weight \(1/9\).
A tilt replaces \(C_0\) by
\(C=(C_0^{-1}+2\operatorname{diag}(r_i))^{-1}\) and multiplies its
density by \(\sqrt{\det C/\det C_0}\).
The scalar \(q_0\) in (FV4) cancels the factor \(p^{-1}\) from \(\mu\),
leaving \(5/2\). Thus the required moments are finite combinations of
orthant Gaussian moments of dimension at most three. Since the
permutation orbit of each tilt is included, it is enough to evaluate
one representative covariance and average over permutations of the
moment exponents.

Write \(M_C(b)=\int_{\mathbb R_+^d}x^b\phi_C(x)\,dx\), where
\(\phi_C\) is the centered Gaussian density with covariance \(C\).
Use (FV1) for \(M_C(0)\), with the corresponding pair and half-line
versions in dimensions two and one. For \(b_i>0\), put \(v=b-e_i\).
The moment recursion is
\[
 M_C(b)=\sum_j C_{ij}\begin{cases}
   v_j M_C(v-e_j),&v_j>0,\\[2pt]
   (2\pi C_{jj})^{-1/2} M_{C\mid j}(v_{\widehat j}),&v_j=0,
 \end{cases}
 \tag{FV7}
\]
where \(C\mid j=C_{\widehat j,\widehat j}
-C_{\widehat j,j} C_{j,\widehat j}/C_{jj}\).
Indeed \(x_i\phi_C=-\sum_j C_{ij}\partial_j\phi_C\), so integration
by parts lowers the interior degree and possibly contributes
the face at \(x_j=0\). The Gaussian on that face has conditional
covariance \(C\mid j\). Every recursion lowers degree or dimension.
This gives ball enclosures without angular quadrature even for the tilted
trihedral cone. To compute \(\int e v_i v_j\,d\mu\), fix one representative
of \(v_i\) and average the exponents of \(v_j\), using the invariance of
\(e\,d\mu\).

## Certified range

The [checker](check_five_cell_weighted_gap.py) generates exact rational
coefficient matrices and computes the weighted Gram and even power
traces with Arb. The following are strict outward-rounded bounds.

| Slab | Ambient degree \(N\) | Trace power \(j\) | Upper bound on \(\operatorname{tr}((QG)^j)\) |
|---|---:|---:|---:|
| \([0,1/3]\) | 6 | 2 | 0.310292 |
| \([1/3,1/2]\) | 8 | 2 | 0.772691 |
| \([1/2,3/5]\) | 10 | 4 | 0.490624 |
| \([3/5,7/10]\) | 14 | 4 | 0.841668 |
| \([7/10,3/4]\) | 16 | 4 | 0.998126 |
| \([3/4,4/5]\) | 20 | 16 | 0.901935 |
| \([4/5,41/50]\) | 21 | 16 | 0.895836 |

The complete slab cover proves the stated gap. The replay is internal.
Run from the repository root with Python 3.12.3 and python-flint 0.9.0:

~~~sh
.venv/bin/python research/plurality-stability/check_five_cell_weighted_gap.py
~~~

The [replay record](five_cell_weighted_gap_verification.json) retains the
input degrees, powers, and Arb bounds. Independent mathematical review
and formal verification remain open.
