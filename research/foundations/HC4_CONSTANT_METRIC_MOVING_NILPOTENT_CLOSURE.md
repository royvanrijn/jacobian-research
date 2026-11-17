# Constant-metric closure for moving nilpotent `HC4` pencils

## Status and scope

This note continues the nilpotent relative-pencil formulation `HC4RSD17` and
the constant-frame classification `HC4RSD42`.

The point is deliberately different from `HC4RSD42`: the relative nilpotent
endomorphism is now allowed to vary polynomially.  Only the Hessian metric is
assumed constant.

> **Theorem HC4RSD54 — constant metric, arbitrarily moving nilpotent frame.**
> Let `K` have characteristic zero, let
>
> \[
> S=\operatorname{Hess}\psi\in\operatorname{Mat}_4(K),
> \qquad \det S\ne0,
> \]
>
> be constant, and let
>
> \[
> T(x)=\operatorname{Hess}A(x).
> \]
>
> Assume
>
> \[
> \det(S+sT(x))=\det S
> \qquad\text{for every }s.
> \tag{0.1}
> \]
>
> Equivalently, `HC4RSD17` says that
>
> \[
> N(x)=S^{-1}T(x)
> \]
>
> is a polynomially varying nilpotent matrix.  Then every gradient pencil
>
> \[
> \nabla(\psi+sA):\mathbb A^4\longrightarrow\mathbb A^4
> \]
>
> is a polynomial automorphism.  In particular no `HC4` counterexample can
> have constant Hessian metric, even when the Jordan flag of `N(x)` moves.

Thus the genuine moving-frame frontier requires **both** the metric `S(x)` and
the nilpotent relative frame `N(x)` to move.

## 1. Constant congruence removes the metric

After scalar extension if necessary, choose a constant
\(C\in GL_4(K)\) with

\[
C^{\mathsf T}SC=I.
\tag{1.1}
\]

Put

\[
\widetilde\psi(y)=\psi(Cy),
\qquad
\widetilde A(y)=A(Cy).
\]

Then

\[
\operatorname{Hess}\widetilde\psi=I,
\qquad
H(y):=\operatorname{Hess}\widetilde A(y)=C^{\mathsf T}T(Cy)C.
\tag{1.2}
\]

The matrix `H(y)` is symmetric because it is a Hessian.  It is also similar
to the relative nilpotent matrix:

\[
C^{-1}N(Cy)C
 =C^{-1}S^{-1}T(Cy)C
 =C^{\mathsf T}T(Cy)C
 =H(y),
\tag{1.3}
\]

where (1.1) gives \(S^{-1}=CC^{\mathsf T}\).  Hence `H(y)` is nilpotent for
every `y`.

Therefore

\[
F_s(y)=y+s\nabla\widetilde A(y)
\tag{1.4}
\]

has nilpotent symmetric nonlinear Jacobian

\[
J(s\nabla\widetilde A)=sH(y).
\]

## 2. Apply the dimension-four symmetric-nilpotent theorem

De Bondt and van den Essen proved that, in characteristic zero and dimension
at most four, every polynomial map

\[
x+H(x)
\]

whose nonlinear Jacobian `JH` is both symmetric and nilpotent is a polynomial
automorphism.  Their 2005 paper removes the homogeneity assumption present in
the earlier low-dimensional result.

Applying that theorem to (1.4) shows that `F_s` is a polynomial automorphism
for every scalar `s`.

Finally

\[
\nabla(\psi+sA)(Cy)
 =C^{-\mathsf T}\bigl(y+s\nabla\widetilde A(y)\bigr)
   +\text{constant},
\tag{2.1}
\]

so the original gradient pencil is obtained from `F_s` by invertible affine
changes on source and target.  It is therefore a polynomial automorphism.

This proves `HC4RSD54`.

## 3. Why this matters for the remaining frontier

The previous constant-Jordan theorem `HC4RSD42` fixed `N` and allowed `S` to
move.  `HC4RSD54` is transverse to it:

- `HC4RSD42`: **constant `N`**, moving Hessian metric allowed;
- `HC4RSD54`: **constant `S`**, arbitrarily moving nilpotent `N` allowed.

Thus a genuinely new `HC4` mechanism must live in the intersection of the two
moving phenomena:

\[
\boxed{\text{both }S(x)\text{ and the nilpotent flag of }N(x)\text{ must move}.}
\]

This also explains why trying to prove that every nilpotent flag is constant
is too strong: symmetric Hessian-nilpotent polynomial matrices do move.  In
dimension four their movement is nevertheless already harmless when measured
in a constant metric.

## 4. External theorem used

M. de Bondt and A. van den Essen,
*Nilpotent symmetric Jacobian matrices and the Jacobian Conjecture II*,
Journal of Pure and Applied Algebra **196** (2005), 135--148.
The main result states that the Jacobian conjecture holds for maps `x+H` with
`JH` nilpotent and symmetric when `n<=4`, without requiring `H` to be
homogeneous.
