# Ordinary sextic Keller maps of weight \((1,-1,-2)\)

**Theorem SexticKellerGradedOneNegativeOneNegativeTwo.**
Let \(k\) be a field of characteristic zero. If a Keller map
\(F:k^3\to k^3\) of ordinary degree at most six is equivariant
for the diagonal action with weight vector \((1,-1,-2)\) in
linear source and target frames, then \(F\) is a polynomial
automorphism. In those frames its positive-weight component
is a scalar times the positive-weight source coordinate.
The ordinary degree minimum for an equivariant counterexample
of this weight is seven, over every such field.

## Reduction to four polynomial ideals over \(\mathbb Q\)

Use source coordinates \(x,y,z\) of weights \(1,-1,-2\), and
order the destination coordinates by the same weights.
Equivariance forces \(F(0)=0\) and its linear part to be diagonal.
The Keller condition makes all three diagonal entries nonzero.
Rescale the three outputs to make them one, so \(\det JF=1\).

Put \(u=xy,\ v=x^2z\). On \(x\ne0\), the components are
\[
            xL(u,v),\quad x^{-1} B(u,v),\quad x^{-2} A(u,v).
\]
For a component of weight \(e\), \(x^e u^iv^j\) has \(x\)-exponent
\(e+i+2j\) and ordinary degree \(e+2i+3j\). Both \(i,j\)
are nonnegative. The exponent condition and the cap of six give
all allowed terms:
\[
\begin{aligned}
 L&=1+l_1u+l_2 u^2+l_3v+l_4uv,\\
 B&=u+b_2u^2+b_3u^3+(b_4+b_5u+b_6u^2)v+b_7v^2,\\
 A&=v+a_1u^2+a_2u^3+a_3u^4+(a_4u+a_5u^2)v+(a_6+a_7u)v^2 .
\end{aligned}
\]
Each displayed coefficient is free before the determinant equations.
Distinct monomials in \(u,v\) lift to distinct source monomials.

The quotient substitution \((x,y,z)\mapsto(x,u,v)\) has
Jacobian \(x^3\). The quotient outputs are \(LB,L^2A\), and
the first output is \(xL\). The chain rule therefore gives
\[
 \{LB,L^2A\}_{u,v}=L^2
 \quad\Longleftrightarrow\quad
 L\{B,A\}+2A\{B,L\}+B\{L,A\}=1 .
\]
Here \(\{S,T\}=S_uT_v-S_vT_u\).
This is the quotient bracket of
[Kistner and Shaska, Lemma 2.1, Corollary 2.4 and equation (28)](sources/kistner-shaska-graded-keller-2608.02863v1.pdf).
Their conditional bounded classification in Theorem 8.5
uses a different box of degrees in \(u,v\). Our proof uses the
ordinary source degree and a separate set of coefficient ideals.

If \(L\ne1\), take the first nonzero coefficient in the order
\(l_4,l_2,l_3,l_1\). Over the algebraic closure, nonzero diagonal
rescalings of \(y,z\) set this coefficient to one, and inverse
rescalings of the corresponding outputs keep the linear part
the identity. These changes preserve the degree and formula
for the allowed supports. In the \(l_2\) case choose a square
root for the \(u\)-scaling. The earlier coefficients are zero,
so the four cases cover every geometric point with \(L\ne1\).

In each case take all coefficients in \(u,v\) of the residual
of the bracket equation. The
[exact checker](verify_sextic_graded_112.py)
generates these coefficient ideals from the displayed supports.
Singular's exact characteristic-zero Gröbner basis computation
returns the reduced basis \(\{1\}\) in each chart.
The four checks show that \(L=1\). The argument over the
algebraic closure shows the same conclusion for the original
coefficients. The positive-weight component is thus \(x\) in the
normalized frame, and the
[bounded coordinate theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
proves inversion. The [degree-seven example](sources/gao-2608.00222v2.pdf),
Theorem 3.3, has the same equivariance and a collision over
\(\mathbb Q\), hence over every characteristic-zero field.

The calculation uses the four charts over \(\mathbb Q\).
The checker verifies a triangular positive control and the
normalized degree-seven example, and detects the failure of
its determinant when the degree-seven monomial is removed.

Run `.venv/bin/python research/jacobian/monodromy-theory/verify_sextic_graded_112.py`
from the repository root with SymPy 1.14.0 and Singular 4.3.2.
The four exact unit ideals were replayed internally using those versions.
