# The Meng--Yang Schur-descent bridge

## Status and source

This note imports Lemmas A.1--A.2 and Proposition A.3 of Meng and Yang,
*A five-variable counterexample to the Hessian conjecture, and the
low-dimensional status of the Jacobian and Hessian conjectures*,
[arXiv:2607.22198v2](https://arxiv.org/abs/2607.22198v2).  The paper was
submitted on July 24, 2026 and revised on July 27, 2026.  The proof is
reproduced below. For an affine source coordinate of a Keller map, the
doubled construction also gives an equivalence of the entire descended
gradient map with a product of the original map.

The headline identities were also replayed from the authors'
[verification repository at commit
`7b23ada`](https://github.com/malyang/hc5-counterexample/commit/7b23ada77a7369a4e75ce3db70227dc716a09de6):
the dependency-free checker confirms degree \(14\), \(42\) monomials,
determinant \(128\), and the \(y=0\) collision for the v2 representative;
the symbolic Schur checker confirms
\(\det M=0\), bordered determinant \(-4\), and descended determinant \(128\).

The exact conclusion is slightly more conditional than the shorthand
"constant bordered Hessian descends": the unbordered Hessian pencil must
also be singular identically.  This second hypothesis is automatic for the
doubled Keller potentials used by Meng--Yang, but not for an arbitrary
polynomial affine-linear in one variable.

## 1. General Schur descent

Let \(K\) be a field of characteristic zero, let
\(w=(w_1,\ldots,w_m)\), and write

\[
 \Phi(t,w)=tA(w)+B(w).
\]

Put

\[
 g=\nabla_w A,\qquad
 M(s,w)=\operatorname{Hess}_w B+s\operatorname{Hess}_w A.
\]

Then

\[
 \operatorname{Hess}_{t,w}\Phi
 =
 \begin{pmatrix}
  0 & g^T\\
  g & M(t,w)
 \end{pmatrix}.
\]

> **Meng--Yang Schur-descent lemma.**  Suppose
>
> \[
>  \det\operatorname{Hess}_{t,w}\Phi=c\in K^\times,
>  \qquad
>  \det M(s,w)\equiv0
> \]
>
> identically in \(s,w\).  For every
> \(\lambda\in K^\times\) and \(\mu\in K\), set
>
> \[
>  \psi_{\lambda,\mu}(w)
>  =B(w)+\frac{\lambda}{2}A(w)^2+\mu A(w).
> \]
>
> Then
>
> \[
>  \det\operatorname{Hess}\psi_{\lambda,\mu}=-\lambda c.
> \tag{1.1}
> \]

Thus \(\mu\) does not affect the Hessian determinant, while \(\lambda\)
scales it.

### Proof

Quadratically repair the affine pivot:

\[
 \widehat\Phi(t,w)
 =\Phi(t,w)-\frac{(t-\mu)^2}{2\lambda}.
\]

Its critical equation in \(t\) has the polynomial solution

\[
 t^*=\mu+\lambda A(w),
\]

and its critical value is \(\psi_{\lambda,\mu}\).  The \(t,t\) pivot of
\(\operatorname{Hess}\widehat\Phi\) is \(-1/\lambda\), so the Schur
complement identity gives

\[
 \det\operatorname{Hess}_{t,w}\widehat\Phi
 =-\frac1\lambda
   \det\operatorname{Hess}\psi_{\lambda,\mu}.          \tag{1.2}
\]

For any scalar \(a\), vector \(g\), and square matrix \(M\),

\[
 \det
 \begin{pmatrix}
  a&g^T\\
  g&M
 \end{pmatrix}
 =a\det M-g^T\operatorname{adj}(M)g.                  \tag{1.3}
\]

Because \(\det M(s,w)\equiv0\), changing the bordered corner from \(0\)
to \(-1/\lambda\) does not change the determinant.  Hence
\(\det\operatorname{Hess}\widehat\Phi=c\), and (1.2) gives (1.1).

Equivalently, direct differentiation gives

\[
 \operatorname{Hess}\psi_{\lambda,\mu}
 =M(\mu+\lambda A,w)+\lambda gg^T,
\]

and the rank-one determinant identity gives the same conclusion.

## 2. Collision transfer

The constant-Hessian statement holds for every \(\mu\), but preservation
of a specified collision imposes a normalization.

Suppose

\[
 (t_0,w_+)\ne(t_0,w_-),\qquad
 \nabla\Phi(t_0,w_+)=\nabla\Phi(t_0,w_-).
\]

Because \(\partial_t\Phi=A\), the equality already implies

\[
 A(w_+)=A(w_-)=A_0.
\]

Choose

\[
 \mu=t_0-\lambda A_0.                                \tag{2.1}
\]

At both points,

\[
\begin{aligned}
 \nabla\psi_{\lambda,\mu}
 &=\nabla B+(\mu+\lambda A)\nabla A\\
 &=\nabla B+t_0\nabla A
 =\nabla_w\Phi(t_0,w),
\end{aligned}
\]

so

\[
 \nabla\psi_{\lambda,\mu}(w_+)
 =\nabla\psi_{\lambda,\mu}(w_-).
\]

Equation (2.1) keeps the specified pair at its original points \(w_\pm\).

## 3. The doubled-Keller bridge

Let \(F:K^n\to K^n\) be a Keller map that is affine-linear in the source
coordinate \(t=x_n\):

\[
 F(x',t)=F_0(x')+tF_1(x'),\qquad
 x'=(x_1,\ldots,x_{n-1}).
\]

For the doubled potential

\[
 \Phi(x',t,y)=\langle y,F(x',t)\rangle=tA+B
\]

one has

\[
 A=\langle y,F_1\rangle,\qquad B=\langle y,F_0\rangle,
\]

and the usual doubling identity gives

\[
 c=\det\operatorname{Hess}\Phi
   =(-1)^n(\operatorname{Jac}F)^2\in K^\times.        \tag{3.1}
\]

Here \(w=(x',y)\) contains \(n-1\) source variables and \(n\) dual
variables.  For every \(s\), the polynomial \(B+sA\) is linear in \(y\).
Consequently the \(n\) rows of \(M(s,w)\) indexed by the dual variables
are supported in only the \(n-1\) columns indexed by \(x'\).  They are
dependent, and therefore

\[
 \det M(s,w)\equiv0.                                  \tag{3.2}
\]

This row-count is the structural reason Schur descent applies to a doubled
Keller map; it replaces any homogeneous-face search for the same
determinant cancellation.

If

\[
 x'_+\ne x'_-,\qquad
 F(x'_+,t_0)=F(x'_-,t_0),
\]

then at the two points \((x'_\pm,y=0)\) the doubled gradients agree and
\(A=0\).  Taking \(\mu=t_0\) in (2.1) yields the
\((2n-1)\)-variable family

\[
 \psi_\lambda
 =B+\frac{\lambda}{2}A^2+t_0A
\]

with

\[
 \det\operatorname{Hess}\psi_\lambda
 =(-1)^{n+1}\lambda(\operatorname{Jac}F)^2             \tag{3.3}
\]

and a gradient collision at \((x'_\pm,0)\).  This is Meng--Yang's
odd-dimensional Schur-descent bridge.

### Polynomial equivalence of the descended map

The same construction with
\(\psi_{\lambda,\mu}=B+\lambda A^2/2+\mu A\) works for every \(\mu\).
Write \(t=\mu+\lambda A(x',y)\) and
\(\eta=D_{x'}F(x',t)^{\mathsf T}y\), where \(D_{x'}F\) is taken
with \(t\) fixed. The map
\[
 {\cal C}(x',y)=(x',t,\eta)
\]
is a polynomial isomorphism of affine spaces. Its inverse fixes \(x'\) and
recovers \(y\) by
\[
 y=DF(x',t)^{-\mathsf T}
   \begin{pmatrix}\eta\\(t-\mu)/\lambda\end{pmatrix}.
 \tag{3.4}
\]
Indeed, the last column of \(DF(x',t)\) is \(F_1(x')\), so multiplication
by \(DF(x',t)^{\mathsf T}\) gives precisely the \(\eta\) and \(A\)
equations defining \({\cal C}\). The inverse is polynomial by the
adjugate identity and the assumption \(\det DF\in K^\times\).
Differentiating \(B+\lambda A^2/2+\mu A\) gives
\[
 \nabla\psi_{\lambda,\mu}(x',y)=(\eta,F(x',t)).
 \tag{3.5}
\]
Thus \(\nabla\psi_{\lambda,\mu}\) is \(F\times
\mathrm{id}_{\mathbb A^{n-1}}\) after the change \({\cal C}\) in the source
and a block permutation in the target. This proof applies without a
collision hypothesis.

For each target \((\eta_0,b)\) over any commutative \(K\)-algebra, the
two formulas give inverse morphisms between its fibre and \(F^{-1}(b)\).
The isomorphism preserves the full fibre scheme and is natural in the
base algebra. In particular, over a number field it preserves whether a
given fibre has a point over the field or over any of its completions.
An arithmetic example can be transferred whenever its map is in the stated
affine-linear source form.

For the generic invariants, the product adjoins the same independent
\(\eta\)-variables to the source and target function fields. It preserves
the degree and the permutation action of the Galois group of the normal
closure, before and after extending the constant field to an algebraic
closure. If \(S\) is the finite normalization of the target in the source
function field for \(F\), the new normalizing algebra is \(S[\eta]\).
It is finite over the new target ring, normal, and has the required
fraction field. The source open immersion into the normalization is
transported by the same product. Finiteness over an open set of the target
persists under the product, and finiteness on a neighborhood of
\((\eta_0,b)\) implies finiteness on a neighborhood of \(b\) by pullback
along the constant section \(\eta=\eta_0\).

Consequently each pair \((\lambda,\mu)\), with \(\lambda\ne0\), gives the
same map-equivalence class for a fixed \(F\). If \(F\) is noninvertible,
every such descended map is noninvertible. A collision of \(F\) over any
target lifts by taking the same \(\eta_0\) at both source points and using
(3.4). For the common-\(t_0\) collision above, \(\mu=t_0\) places both
lifted points on \(y=0\).

### A jointly affine block of source coordinates

The proof also applies to \(1\le k\le n\) jointly affine source
coordinates \(t=(t_1,\ldots,t_k)\), with
\(F(x',t)=F_0(x')+F_1(x')t\) and \(F_1(x')\) an \(n\)-by-\(k\) matrix.
Set \(A=F_1(x')^{\mathsf T}y\) and \(B=y^{\mathsf T}F_0(x')\).
For any polynomial \(\chi\) on \(\mathbb A^k\) whose gradient is a
polynomial automorphism, let \(\psi=B+\chi(A)\). Now put
\(t=\nabla\chi(A)\) and \(\eta=D_{x'}F(x',t)^{\mathsf T}y\).
Equations (3.4) and (3.5) still hold with the last block in (3.4)
replaced by \((\nabla\chi)^{-1}(t)\). The resulting gradient is
equivalent to \(F\times\mathrm{id}_{\mathbb A^{n-k}}\) in dimension
\(2n-k\). In particular, one can take a nondegenerate quadratic \(\chi\)
plus a linear term. Application with \(k>1\) requires joint affine
dependence in the indicated source block.

If \(n\ge2\) and \(k=n-1\), the original Keller map in this construction
is an automorphism. Write the single remaining variable as \(s\), and let
\(V(s)\) be the span over \(K(s)\) of the columns \(v_i(s)\) of \(F_1(s)\).
It has dimension \(n-1\) since \(DF\) is invertible. The coefficient of
\(t_i\) in the Jacobian determinant is
\(\det(v_i',v_1,\ldots,v_{n-1})=0\), so every \(v_i'\) lies in \(V\).
Hence the derivative of \(v_1\wedge\cdots\wedge v_{n-1}\) is proportional
to the wedge itself. The ratios of its nonzero coordinates have derivative
zero and are in \(K\). Thus \(V\) is the extension of a fixed hyperplane
over \(K\). Choose target linear coordinates whose first row annihilates
that hyperplane. The map in these coordinates is
\((a(s),B_0(s)+B_1(s)t)\). Its Jacobian is \(a'(s)\det B_1(s)\), a unit of
\(K[s]\). Both factors are therefore units, and \(s\) and then \(t\) can be
recovered polynomially. In particular, a jointly affine two-coordinate
block in a three-dimensional Keller map supplies only invertible gradients
by this descent. The planar affine case is also covered over the reals by
[Theorem 1 of *Global injectivity of planar non-singular maps polynomial in one variable*](https://arxiv.org/abs/2302.05394).

## 4. The five-variable counterexample

For the three-variable Keller counterexample used by Meng--Yang,
\(\operatorname{Jac}F=-2\), the pivot is \(t=x_3\), and the two points

\[
 (1,-3/2,13/2),\qquad(-1,3/2,13/2)
\]

have the same image.  Thus \(c=-4\).  With
\(\lambda=1\), \(\mu=t_0=13/2\), and

\[
 \psi=B+\frac12A^2+\frac{13}{2}A,
\]

equation (1.1) gives

\[
 \det\operatorname{Hess}\psi=4.
\]

Clearing the denominator by setting

\[
 \Psi=2\psi=A^2+13A+2B
\]

multiplies the five-variable Hessian determinant by \(2^5\), so

\[
 \det\operatorname{Hess}\Psi=128.
\]

The polynomial has total degree \(14\), and its gradient agrees at

\[
 (1,-3/2,0,0,0),\qquad(-1,3/2,0,0,0).
\]

Therefore \(\mathrm{HC}_5\) is false.  Stabilization gives falsehood in
every dimension at least five; together with the known truth in dimensions
at most three, this leaves only \(\mathrm{HC}_4\) unresolved.

The v2 representative above is the fixed parent used by the retained
[normal-jet formulas](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md) and
[quintic graph exclusion](../jacobian/affine-complements/FULL_COMPLEX_QUINTIC_GRAPH_EXCLUSION_AND_CUBIC_NILPOTENCE.md).
The Schur theorem retains its explicit structural hypotheses; arbitrary
variable elimination is not justified by these formulas.

## References

- Guowu Meng and Liang Yang, [*A five-variable counterexample to the
  Hessian conjecture, and the low-dimensional status of the Jacobian and
  Hessian conjectures*, arXiv:2607.22198v2](https://arxiv.org/abs/2607.22198v2),
  especially Appendix A.
- Meng--Yang, [exact-arithmetic verification
  repository](https://github.com/malyang/hc5-counterexample).
