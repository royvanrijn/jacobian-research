# Rank-one sextic reduction and the complete coordinate Meng chart

## Status

Let

\[
 \psi=q_2+h_3+h_4+h_6
\]

be a collision-normalized Meng potential in four variables.  Here \(q_2\)
is nondegenerate and \(h_d\) is an arbitrary homogeneous form of degree
\(d\).  There is no support restriction.

> **Theorem `HC4T11`.**  If
> \(\operatorname{Hess}(h_6)\) has generic rank one and
> \(\det\operatorname{Hess}(\psi)\) is a nonzero constant, then
> \(\nabla\psi\) cannot identify a nonzero antipodal pair.

Together with `HC4T31`, `HC4T21`, and `HC4CQ1`, this gives:

> **Corollary `HC4TC1`.**  No potential
> \(q_2+h_3+h_4+h_6\) in the four-variable coordinate Meng chart has both
> nonzero constant Hessian determinant and a nonzero antipodal gradient
> collision.

This closes the coordinate chart with homogeneous support
\(\{2,3,4,6\}\).  It does not remove quintic or higher homogeneous
layers, and it does not treat a non-coordinate coisotropic embedding.

We work after scalar extension to an algebraic closure.  This preserves
the constant-Hessian identity and any collision.

## 1. The sextic is a sixth power

The rank-one polynomial-Hessian normal form says that, up to a linear
term, a polynomial with Hessian rank one is composite-univariate.  Since
\(h_6\) is homogeneous of degree six, the linear remainder vanishes and
the univariate polynomial has only its degree-six term.  After a linear
change,

\[
 h_6=cx^6,\qquad c\ne0.                                \tag{1.1}
\]

Thus the sextic Hessian has the constant three-dimensional kernel

\[
 W=\langle u,v,w\rangle.                               \tag{1.2}
\]

In binary-root language, (1.1) is the highest-weight point \(L^6\).  It
is more special than the binary-sextic `SIC(2)` nullcone \(L^4Q\), but the
same projective-root mechanism will reappear on the lower faces.

Put

\[
 H_0=\operatorname{Hess}(q_2),\quad
 A=\operatorname{Hess}(h_3),\quad
 B=\operatorname{Hess}(h_4),\quad
 C=\operatorname{Hess}(h_6).
\]

Spatial scaling gives

\[
 \det(H_0+\lambda A+\lambda^2B+\lambda^4C)=\det H_0.    \tag{1.3}
\]

The degree-ten coefficient is

\[
 \bar C\,\det(B_W),                                    \tag{1.4}
\]

where \(\bar C\ne0\) is the one-dimensional sextic quotient block.
Consequently the ternary Hessian \(B_W\) is singular.

Over \(K(x)\), the three-variable singular-Hessian theorem supplies a
constant-in-\(W\) projective kernel for \(B_W\).  It is also constant in
\(x\).  Indeed, quartic homogeneity gives

\[
 B_W(\rho x,\rho W)=\rho^2B_W(x,W).                    \tag{1.5}
\]

After replacing \(W\) by \(\rho W\), the projective kernel scheme at
\(\rho x\) equals the one at \(x\).  The nonzero one-dimensional
\(x\)-base modulo scaling is a point.  Hence its kernel line, or its
kernel plane in rank one, can be chosen over the constant field.

We now split according to the generic rank of \(B_W\).

## 2. Quartic residual rank two

Let \(v\in W\) span the constant kernel of \(B_W\).  The degree-nine
coefficient in (1.3) is

\[
 \bar C\,
 \operatorname{tr}\!\left(
   \operatorname{adj}_2(B_W)A_W
 \right),                                               \tag{2.1}
\]

where \(\operatorname{adj}_2\) denotes the rank-two compound.  In a
congruence chart \(B_W=\operatorname{diag}(b_1,b_2,0)\), this is

\[
 \bar C\,b_1b_2 A_{vv}.                                 \tag{2.2}
\]

Therefore

\[
 D_vh_6=0,\qquad D_v^2h_4=0,\qquad D_v^2h_3=0.          \tag{2.3}
\]

The common-direction reduction of `HC4T31` excludes the collision.

## 3. Quartic residual rank one

Let \(K\subset W\) be the constant two-plane kernel of \(B_W\).  The
degree-nine face vanishes identically.  At degree eight there appear to
be two competing terms:

\[
 \det B
 \quad\text{and}\quad
 \bar C\,\operatorname{Mix}(B_W,A_W,A_W).               \tag{3.1}
\]

The first term is automatically zero.  A symmetric \(1+3\) block matrix
whose \(3\times3\) principal block has rank one has total rank at most
three: the two-plane kernel of the principal block meets the kernel of
the one cross functional nontrivially.

Terms of the form
\(\bar C\,\operatorname{Mix}(B_W,B_W,H_{0,W})\) also vanish because
\(\operatorname{rank}B_W=1\).  On a chart
\(B_W=\operatorname{diag}(b,0,0)\), the complete degree-eight face is
therefore

\[
 \bar C\,b\det(A_K)=0.                                  \tag{3.2}
\]

Choose constant coordinates on \(K\).  Because \(h_3\) is cubic,

\[
 A_K=
 \begin{pmatrix}a&b'\\ b'&d\end{pmatrix}
\]

has linear-form entries, and (3.2) says

\[
 ad-(b')^2=0.                                           \tag{3.3}
\]

This is exactly the binary \(\operatorname{Sym}^2\) nullcone used in the
`SIC(2)` root classification.  Unique factorization of linear forms
turns (3.3) into

\[
 A_K=\ell
 \begin{pmatrix}r^2&rs\\rs&s^2\end{pmatrix}             \tag{3.4}
\]

for constants \(r,s\) and a linear form \(\ell\), including the evident
degenerate cases.  The constant repeated-root direction
\((-s,r)\in K\) lies in the kernel of \(A_K\).  It consequently satisfies
(2.3), and `HC4T31` again applies.

This is the precise place where the SIC binary-root machinery strengthens
the Schur-face calculation: determinant zero is not left as a moving
rational cone; degree one forces the \(\operatorname{Sym}^2\) nullcone
root to be constant.

## 4. Quartic residual rank zero

If \(B_W=0\), then homogeneity gives

\[
 h_4=ax^4+x^3\ell(W).                                   \tag{4.1}
\]

Thus the full Hessian \(B\) has rank at most two.  The possible
degree-seven term with three \(B\)-columns and one \(A\)-column vanishes.
Every \(C\)-selected term containing \(B_W\) also vanishes.  The
degree-seven face of (1.3) is exactly

\[
 \bar C\,\det(A_W)=0.                                   \tag{4.2}
\]

Apply the ternary singular-Hessian theorem to \(h_3\) over \(K(x)\).
Cubic homogeneity,

\[
 A_W(\rho x,\rho W)=\rho A_W(x,W),                      \tag{4.3}
\]

makes its projective kernel scheme constant over the one-dimensional
base by the same argument as (1.5).  Hence there is a nonzero constant
\(v\in W\) with \(D_v^2h_3=0\).  Since \(B_W=0\) and \(h_6=cx^6\), this
direction again satisfies (2.3).  This finishes the rank-one sextic
theorem.

## 5. Exhaustion of the coordinate chart

The top coefficient of (1.3) first makes
\(\det\operatorname{Hess}(h_6)=0\), so the generic sextic Hessian rank is
at most three.  The four possibilities are now:

| sextic Hessian rank | result |
|---:|---|
| \(3\) | `HC4T31` |
| \(2\) | `HC4T21` |
| \(1\) | `HC4T11` above |
| \(0\) | \(h_6=0\), then `HC4CQ1` |

These exhaust \(q_2+h_3+h_4+h_6\) and prove `HC4TC1`.  This is a
complete coordinate chart only for the displayed homogeneous support:
quintic and higher homogeneous layers are not removed by this exhaustion.

## Reproduction

Run:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_one_reduction.py
```

The checker verifies every displayed determinant face, the disappearance
of both same-weight contaminants at degree eight, the
\(\operatorname{Sym}^2\) repeated-root parametrization, the homogeneous
rank-one and rank-zero quartic forms, and the one-base scaling identities
on all cubic and quartic monomials.  It also replays the common-direction
descent `HC4T31`.

The external structural inputs are de Bondt's
[small-rank polynomial-Hessian normal forms](https://arxiv.org/abs/1609.03904),
the low-dimensional
[Gordan--Noether/singular-Hessian classification](https://arxiv.org/abs/1501.05168),
the Hessian conjecture in dimensions at most three, and
[Moh's plane degree bound](https://www.math.purdue.edu/~ttm/jacobian.pdf).
