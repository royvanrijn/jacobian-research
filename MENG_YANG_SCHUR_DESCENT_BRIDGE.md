# The Meng--Yang Schur-descent bridge

## Status and source

This note imports Lemmas A.1--A.2 and Proposition A.3 of Meng and Yang,
*A five-variable counterexample to the Hessian conjecture, and the
low-dimensional status of the Jacobian and Hessian conjectures*,
[arXiv:2607.22198v2](https://arxiv.org/abs/2607.22198v2).  The paper was
submitted on July 24, 2026 and revised on July 27, 2026.  The proof is
reproduced below because the construction is a reusable determinant lemma,
not a coefficient-specific cancellation.

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

This distinction is important: arbitrary \(\mu\) preserves the determinant,
whereas (2.1) is what transfers a given equal-gradient pair.

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

The July 24 arXiv v1 displayed the equally constant-Hessian representative
\(A^2+11A+2B\), with its collision at nonzero dual coordinates.  The
current v2 uses \(A^2+13A+2B\) because the collision is visible at \(y=0\).
The repository's
[nonlinear toric descent audit](HC5_NONLINEAR_TORIC_DESCENT.md) is pinned
to the v1 representative; that convention does not change its exact
calculation.

## 5. Import rule and limits

The projective invariant transported by this bridge is presently the top
degree of the actual affine-gradient compactification, not the full
multidegree list and not the polar degree of the homogenized potential.
The doubled and descended potentials have the same generic degree three,
but their ambient dimensions and gradient degrees are respectively
\((6,7)\) and \((5,13)\).  Their aggregate Segre corrections are therefore
\(7^6-3\) and \(13^5-3\); no componentwise Segre transform under Schur
descent is claimed.  See
[`PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md`](PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md)
for the canonical all-dimensional convention and registry.

When a proposed one-variable descent has the form \(\Phi=tA+B\), the first
test should be the two structural hypotheses

\[
 \det\operatorname{Hess}_{t,w}\Phi\in K^\times,
 \qquad
 \det\operatorname{Hess}_w(B+sA)\equiv0.
\]

Once they hold, (1.1) is the determinant proof.  Homogeneous Hessian-face
calculations may still be needed to establish one of those hypotheses in a
different setting, but should not be used to rediscover the rank-one update
identity.

The lemma does **not** make a second descent automatic.  The term \(A^2\)
usually consumes the dual-linearity that forced (3.2), so the descended
potential need not have another affine pivot with a singular reduced
Hessian pencil.  In particular, the bridge explains the passage from six
to five variables but does not settle \(\mathrm{HC}_4\).

The [continuation note](SCHUR_DESCENT_CONTINUATIONS.md) records the exact
remainder when the reduced pencil is not identically singular, a simultaneous
multi-pivot extension, and the codimension-one block-affine obstruction that
closes the apparent pure-source one-shot route to \(\mathrm{HC}_4\).

## References

- Guowu Meng and Liang Yang, [*A five-variable counterexample to the
  Hessian conjecture, and the low-dimensional status of the Jacobian and
  Hessian conjectures*, arXiv:2607.22198v2](https://arxiv.org/abs/2607.22198v2),
  especially Appendix A.
- Meng--Yang, [exact-arithmetic verification
  repository](https://github.com/malyang/hc5-counterexample).
