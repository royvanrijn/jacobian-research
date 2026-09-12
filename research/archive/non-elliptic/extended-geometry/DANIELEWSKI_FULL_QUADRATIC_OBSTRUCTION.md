# The full quadratic Danielewski obstruction

The coupled Danielewski contraction equation has no solution whose three
target coordinates have ambient degree at most two.

This closes all 127 quadratic support charts at once, including the 28
large-support charts left by the earlier chart-by-chart frontier.

## 1. Normalized ansatz

Use

\[
x=bc,\qquad w=b((bc)^2+1),
\]

and let

\[
\begin{aligned}
\mathcal Q=\operatorname{span}_{\mathbb Q}\{&
1,a,c,x,w,\\
&ac,ax,aw,a^2,c^2,x^2,w^2,cw,cx,xw\}.
\end{aligned}                                      \tag{1}
\]

Only \(a,c,w\) have nonzero source-linear terms.  Therefore any triple
\((Y_1,Y_2,Y_3)\in\mathcal Q^3\) with nonzero constant Jacobian has an
invertible linear coefficient matrix on \((a,c,w)\).  A constant target
linear change normalizes it to

\[
(A,C,W)=(a,c,w)+v_0x+\sum_{j=1}^{10}v_jN_j,        \tag{2}
\]

where \(v_0,v_j\in\mathbb Q^3\) and

\[
(N_1,\ldots,N_{10})
=(ac,ax,aw,a^2,c^2,x^2,w^2,cw,cx,xw).             \tag{3}
\]

The determinant at the source origin is then \(-1\).  Thus the required
identity is

\[
\det\frac{\partial(A,C,W)}{\partial(a,b,c)}+1=0.   \tag{4}
\]

## 2. Exact coefficient ideal

The coefficient \(v_0\) must be retained: \(x=bc\) has zero source-linear
jet, so normalizing the coefficient matrix on \((a,c,w)\) does not remove
it.  Introduce 33 independent coefficients \(z_0,\ldots,z_{32}\) for the
three rows in (2).  Expanding (4) gives 102 distinct source monomials.

Three independent linear coefficient equations have rank three.  Solving
them exactly and removing scalar-duplicate nonzero equations leaves

\[
\boxed{
30\text{ coefficient variables and }89\text{ equations}.
}                                                    \tag{5}
\]

Let \(I_{\mathbb Q}\) be their ideal.

## 3. Rational unit certificate

Order the 30 remaining coefficient variables by residual monomial first
and target row second.  In other words, use column-major order on the
three-by-eleven coefficient matrix after the three linear eliminations.

Singular's exact modular-over-\(\mathbb Q\)
`modGB("slimgb",I,1)` algorithm reconstructs the rational basis from
prime-field calculations and gives

\[
\boxed{
\operatorname{GB}(I_{\mathbb Q})=\{1\}.
}                                                    \tag{6}
\]

The returned basis is over \(\mathbb Q\); it is not an inference from the
three exploratory finite-field reductions and does not rely on a
bounded-height coefficient search.  An independent `msolve` computation
over \(\mathbb Q\) also returns the empty variety.

Equation (6) proves that the normalized coefficient scheme is empty.
Since every nonzero constant-Jacobian triple admits the normalization
(2), no such triple exists.

## 4. Theorem

**Theorem.**  Let \(k\) be a characteristic-zero field.  For

\[
P=x(x^2+1),
\]

no three elements of the ambient-degree-two space

\[
k[a,c,x,w]_{\le2}/(cw-P(x))
\]

pull back under

\[
(a,b,c)\longmapsto(a,c,bc,b((bc)^2+1))
\]

to a polynomial map \(\mathbb A^3\to\mathbb A^3\) with nonzero constant
Jacobian.

**Proof.**  A solution over \(k\) is defined over a finitely generated
characteristic-zero subfield.  The linear normalization is valid after
base extension, and its coefficients satisfy the rational universal ideal
\(I_{\mathbb Q}\).  But (6) says this ideal is the unit ideal. \(\square\)

## 5. Consequences

The entire ambient-quadratic coupled layer is closed:

- no support restriction remains;
- mixing \(a\) into all three target coordinates does not help;
- simultaneous use of \(cx,xw\) and all five diagonal quadratic
  directions does not help;
- the first possible coupled contraction has ambient degree at least
  three.

Thus

\[
\boxed{
\text{every surviving target triple has ambient degree }\ge3.
}                                                    \tag{7}
\]

The de Rham compiler and leading-\(a\)-face theorem now become the primary
tools: degree-three and higher searches should be reduced by flux and
Poisson rank before any full coefficient elimination.

## 6. Reproduction

Run

```bash
.venv/bin/python scripts/verify_danielewski_full_quadratic_obstruction.py
```

The checker reconstructs the normalized target including all three
residual \(x\)-coefficients, expands all 102 source coefficients, performs
the three exact linear eliminations, verifies the
30-variable/89-equation ledger, invokes Singular's exact modular
reconstruction over \(\mathbb Q\) in column-major coefficient order, and
requires the returned basis to be exactly \(\{1\}\).
